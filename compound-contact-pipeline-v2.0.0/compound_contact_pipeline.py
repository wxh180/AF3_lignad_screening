#!/usr/bin/env python3
"""Analyze protein-ligand contacts across compound seed ensembles."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import statistics
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from dataclasses import dataclass
from itertools import combinations, repeat
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import Chain, MMCIFParser, Model, NeighborSearch, ShrakeRupley, Structure
from Bio.PDB.Polypeptide import is_aa


__version__ = "2.0.0"


class PipelineError(RuntimeError):
    """Raised for invalid pipeline inputs or unsafe operations."""


@dataclass(frozen=True)
class DatasetSpec:
    dataset_dir: Path
    dataset_rel: str
    compound: str
    ligand_resname: str | None
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None
    ligand_smiles: str | None = None
    ligand_file: Path | None = None


@dataclass(frozen=True)
class StructureJob:
    cif_path: Path
    structure_id: str
    dataset_rel: str
    compound: str
    ligand_resname: str | None
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None
    ligand_smiles: str | None = None
    ligand_file: Path | None = None


@dataclass(frozen=True)
class QCRecord:
    severity: str
    stage: str
    message: str
    compound: str = ""
    structure_id: str = ""


@dataclass(frozen=True, order=True)
class ResidueKey:
    chain: str
    resseq: int
    icode: str
    resname: str


@dataclass(frozen=True)
class ContactRecord:
    structure_id: str
    compound: str
    residue: ResidueKey
    is_contact: bool
    min_contact_distance_angstrom: float | None
    capped_distance_angstrom: float


@dataclass(frozen=True)
class StructureMetrics:
    ligand_resname: str | None
    n_protein_residues: int
    n_contact_residues: int
    minimum_pair_distance_angstrom: float | None
    clash_pair_count: int
    buried_sasa_total_A2: float | None = None
    interface_area_A2: float | None = None


@dataclass(frozen=True)
class StructureResult:
    job: StructureJob
    status: str
    contacts: tuple[ContactRecord, ...]
    metrics: StructureMetrics
    qc: tuple[QCRecord, ...]


@dataclass(frozen=True)
class MinimizationRecord:
    job: StructureJob
    status: str
    reasons: tuple[str, ...]
    staged_cif_path: Path | None
    output_relative_path: Path | None
    candidate_result: StructureResult | None
    initial_potential_kcal_mol: float | None
    final_potential_kcal_mol: float | None
    backbone_rmsd_A: float | None
    ligand_rmsd_A: float | None
    before_vdw_clash_pair_count: int | None
    after_vdw_clash_pair_count: int | None
    before_maximum_vdw_overlap_A: float | None
    after_maximum_vdw_overlap_A: float | None
    contact_set_jaccard: float | None
    config: "MinimizationConfig"
    package_versions: tuple[tuple[str, str], ...]


COMMON_MONATOMIC_IONS = {
    "AL",
    "BA",
    "BR",
    "CA",
    "CD",
    "CL",
    "CO",
    "CS",
    "CU",
    "F",
    "FE",
    "GA",
    "HG",
    "IOD",
    "K",
    "LI",
    "MG",
    "MN",
    "NA",
    "NI",
    "RB",
    "SR",
    "ZN",
}

R_KCAL_PER_MOL_K = 0.00198720425864083
UNIT_TO_MOLAR = {"M": 1.0, "mM": 1e-3, "uM": 1e-6, "nM": 1e-9, "pM": 1e-12}


@dataclass(frozen=True)
class EnergyValue:
    compound: str
    method: str
    value_kcal_mol: float
    structure_id: str = ""
    uncertainty_kcal_mol: float | None = None
    reference_compound: str = ""
    replicate: str = ""


@dataclass(frozen=True)
class AffinityValue:
    compound: str
    metric: str
    value: float
    unit: str
    concentration_M: float
    temperature_K: float
    replicate: str
    delta_g_standard_kcal_mol: float | None


def parse_chain_list(value: str) -> tuple[str, ...] | None:
    chains = tuple(part.strip() for part in value.split(",") if part.strip())
    return chains or None


def resolve_inside(root: Path, relative_value: str) -> Path:
    resolved_root = root.resolve()
    candidate = (resolved_root / relative_value).resolve()
    try:
        candidate.relative_to(resolved_root)
    except ValueError as exc:
        raise PipelineError(f"dataset path is outside ROOT: {relative_value}") from exc
    return candidate


def find_seed_dirs(base: Path) -> list[Path]:
    candidates = sorted(
        (path.resolve() for path in base.rglob("seed-*") if path.is_dir()),
        key=lambda path: path.as_posix(),
    )
    candidate_set = set(candidates)
    return [
        path
        for path in candidates
        if not any(parent in candidate_set for parent in path.parents)
    ]


def find_cif_files(base: Path) -> list[Path]:
    resolved_base = base.resolve()
    files = {
        cif.resolve()
        for seed_dir in find_seed_dirs(resolved_base)
        for cif in seed_dir.rglob("*.cif")
        if cif.is_file()
    }
    return sorted(files, key=lambda path: path.relative_to(resolved_base).as_posix())


def _parse_enabled(value: str, row_number: int) -> bool:
    normalized = value.strip().lower()
    if not normalized:
        return True
    if normalized == "true":
        return True
    if normalized == "false":
        return False
    raise PipelineError(
        f"manifest row {row_number}: enabled must be true or false, got {value!r}"
    )


def load_manifest(root: Path, manifest_path: Path) -> list[DatasetSpec]:
    resolved_root = root.resolve()
    if not manifest_path.is_file():
        raise PipelineError(f"manifest is not a file: {manifest_path}")

    required = {"dataset_dir", "compound", "ligand_resname"}
    specs: list[DatasetSpec] = []
    seen_directories: set[Path] = set()
    compound_settings: dict[
        str,
        tuple[
            str | None,
            tuple[str, ...] | None,
            str | None,
            str | None,
            Path | None,
        ],
    ] = {}

    with manifest_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        missing_columns = sorted(required - fieldnames)
        if missing_columns:
            raise PipelineError(
                "manifest is missing required columns: " + ", ".join(missing_columns)
            )

        for row_number, row in enumerate(reader, start=2):
            dataset_value = (row.get("dataset_dir") or "").strip()
            compound = (row.get("compound") or "").strip()
            ligand_resname = (row.get("ligand_resname") or "").strip()
            for column, value in (
                ("dataset_dir", dataset_value),
                ("compound", compound),
                ("ligand_resname", ligand_resname),
            ):
                if not value:
                    raise PipelineError(
                        f"manifest row {row_number}: {column} must not be empty"
                    )

            enabled = _parse_enabled(row.get("enabled") or "", row_number)
            dataset_dir = resolve_inside(resolved_root, dataset_value)
            if not enabled:
                continue
            if not dataset_dir.is_dir():
                raise PipelineError(
                    f"manifest row {row_number}: dataset directory does not exist: "
                    f"{dataset_value}"
                )
            if dataset_dir in seen_directories:
                raise PipelineError(
                    f"manifest row {row_number}: duplicate enabled dataset directory: "
                    f"{dataset_value}"
                )
            seen_directories.add(dataset_dir)

            protein_chains = parse_chain_list(row.get("protein_chains") or "")
            ligand_chain = (row.get("ligand_chain") or "").strip() or None
            ligand_smiles = (row.get("ligand_smiles") or "").strip() or None
            ligand_file_value = (row.get("ligand_file") or "").strip()
            ligand_file = (
                resolve_inside(resolved_root, ligand_file_value)
                if ligand_file_value
                else None
            )
            if ligand_smiles and ligand_file:
                raise PipelineError(
                    f"manifest row {row_number}: exactly one of ligand_smiles "
                    "and ligand_file is allowed"
                )
            settings = (
                ligand_resname,
                protein_chains,
                ligand_chain,
                ligand_smiles,
                ligand_file,
            )
            previous = compound_settings.get(compound)
            if previous is not None and previous != settings:
                raise PipelineError(
                    f"manifest row {row_number}: compound {compound!r} has "
                    "conflicting ligand or chain settings"
                )
            compound_settings[compound] = settings

            specs.append(
                DatasetSpec(
                    dataset_dir=dataset_dir,
                    dataset_rel=dataset_dir.relative_to(resolved_root).as_posix(),
                    compound=compound,
                    ligand_resname=ligand_resname,
                    protein_chains=protein_chains,
                    ligand_chain=ligand_chain,
                    ligand_smiles=ligand_smiles,
                    ligand_file=ligand_file,
                )
            )

    return sorted(specs, key=lambda spec: (spec.dataset_rel, spec.compound))


def _job_from_spec(root: Path, cif_path: Path, spec: DatasetSpec) -> StructureJob:
    return StructureJob(
        cif_path=cif_path,
        structure_id=cif_path.relative_to(root).as_posix(),
        dataset_rel=spec.dataset_rel,
        compound=spec.compound,
        ligand_resname=spec.ligand_resname,
        protein_chains=spec.protein_chains,
        ligand_chain=spec.ligand_chain,
        ligand_smiles=spec.ligand_smiles,
        ligand_file=spec.ligand_file,
    )


def build_structure_jobs(
    root: Path,
    manifest_path: Path | None,
    default_ligand: str | None,
    default_protein_chains: tuple[str, ...] | None,
    default_ligand_chain: str | None,
    default_ligand_smiles: str | None = None,
    default_ligand_file: Path | None = None,
) -> list[StructureJob]:
    resolved_root = root.resolve()
    if not resolved_root.is_dir():
        raise PipelineError(f"ROOT is not a directory: {root}")

    jobs: list[StructureJob] = []
    if manifest_path is not None:
        specs = load_manifest(resolved_root, manifest_path)
        candidate_files = {
            cif_path
            for spec in specs
            for cif_path in find_cif_files(spec.dataset_dir)
        }
        for cif_path in sorted(
            candidate_files,
            key=lambda path: path.relative_to(resolved_root).as_posix(),
        ):
            matches = [
                spec for spec in specs if cif_path.is_relative_to(spec.dataset_dir)
            ]
            if len(matches) > 1:
                matched = ", ".join(spec.dataset_rel for spec in matches)
                raise PipelineError(
                    f"CIF {cif_path.relative_to(resolved_root).as_posix()} matches "
                    f"more than one dataset: {matched}"
                )
            if len(matches) == 1:
                jobs.append(_job_from_spec(resolved_root, cif_path, matches[0]))
    else:
        for seed_dir in find_seed_dirs(resolved_root):
            dataset_dir = seed_dir.parent
            compound = dataset_dir.name
            dataset_rel = dataset_dir.relative_to(resolved_root).as_posix()
            for cif_path in sorted(
                (path.resolve() for path in seed_dir.rglob("*.cif") if path.is_file()),
                key=lambda path: path.relative_to(resolved_root).as_posix(),
            ):
                jobs.append(
                    StructureJob(
                        cif_path=cif_path,
                        structure_id=cif_path.relative_to(resolved_root).as_posix(),
                        dataset_rel=dataset_rel,
                        compound=compound,
                        ligand_resname=default_ligand,
                        protein_chains=default_protein_chains,
                        ligand_chain=default_ligand_chain,
                        ligand_smiles=default_ligand_smiles,
                        ligand_file=default_ligand_file,
                    )
                )

    unique_jobs = {job.cif_path: job for job in jobs}
    return sorted(unique_jobs.values(), key=lambda job: job.structure_id)


def residue_key(residue) -> ResidueKey:
    chain = residue.get_parent()
    _, resseq, icode = residue.get_id()
    normalized_icode = "" if icode in (None, " ", "?") else str(icode)
    return ResidueKey(
        chain=str(chain.id),
        resseq=int(resseq),
        icode=normalized_icode,
        resname=residue.get_resname().strip(),
    )


def detect_ligand_resname(model) -> str:
    candidates = {
        residue.get_resname().strip()
        for chain in model
        for residue in chain
        if not is_aa(residue, standard=False)
        and not residue.get_id()[0].startswith("W")
        and residue.get_resname().strip() not in COMMON_MONATOMIC_IONS
    }
    if not candidates:
        raise PipelineError("no ligand candidate found")
    if len(candidates) != 1:
        names = ", ".join(sorted(candidates))
        raise PipelineError(f"multiple ligand candidates found: {names}")
    return next(iter(candidates))


def _invalid_structure_result(
    job: StructureJob,
    stage: str,
    message: str,
    qc_prefix: tuple[QCRecord, ...] = (),
) -> StructureResult:
    qc = (
        *qc_prefix,
        QCRecord(
            severity="error",
            stage=stage,
            message=message,
            compound=job.compound,
            structure_id=job.structure_id,
        ),
    )
    return StructureResult(
        job=job,
        status="invalid",
        contacts=(),
        metrics=StructureMetrics(None, 0, 0, None, 0),
        qc=qc,
    )


def build_sasa_entity(residues, entity_id: str):
    entity = Structure.Structure(entity_id)
    model = Model.Model(0)
    entity.add(model)
    chains = {}
    for residue in residues:
        chain_id = str(residue.get_parent().id)
        if chain_id not in chains:
            chains[chain_id] = Chain.Chain(chain_id)
            model.add(chains[chain_id])
        chains[chain_id].add(deepcopy(residue))
    return entity


def entity_sasa(entity) -> float:
    calculator = ShrakeRupley(probe_radius=1.4)
    calculator.compute(entity, level="S")
    return float(entity.sasa)


def calculate_buried_sasa(
    protein_residues, ligand_residues
) -> tuple[float, float]:
    receptor_sasa = entity_sasa(build_sasa_entity(protein_residues, "receptor"))
    ligand_sasa = entity_sasa(build_sasa_entity(ligand_residues, "ligand"))
    complex_sasa = entity_sasa(
        build_sasa_entity([*protein_residues, *ligand_residues], "complex")
    )
    buried_total = max(0.0, receptor_sasa + ligand_sasa - complex_sasa)
    return buried_total, buried_total / 2.0


def analyze_structure(
    job: StructureJob, cutoff: float = 4.5, clash_cutoff: float = 2.0
) -> StructureResult:
    parser = MMCIFParser(QUIET=True)
    try:
        structure = parser.get_structure(job.cif_path.stem, str(job.cif_path))
        models = list(structure)
    except Exception as exc:
        return _invalid_structure_result(job, "parse", f"failed to parse CIF: {exc}")

    if not models:
        return _invalid_structure_result(job, "parse", "CIF contains no models")

    qc: list[QCRecord] = []
    if len(models) > 1:
        qc.append(
            QCRecord(
                severity="warning",
                stage="model_selection",
                message=f"analyzed first model and ignored {len(models) - 1} additional model(s)",
                compound=job.compound,
                structure_id=job.structure_id,
            )
        )
    model = models[0]

    try:
        ligand_resname = job.ligand_resname or detect_ligand_resname(model)
    except PipelineError as exc:
        return _invalid_structure_result(
            job, "ligand_detection", str(exc), tuple(qc)
        )

    ligand_residues = [
        residue
        for chain in model
        if job.ligand_chain is None or str(chain.id) == job.ligand_chain
        for residue in chain
        if residue.get_resname().strip() == ligand_resname
        and not residue.get_id()[0].startswith("W")
    ]
    if not ligand_residues:
        return _invalid_structure_result(
            job,
            "ligand_selection",
            f"ligand {ligand_resname} not found",
            tuple(qc),
        )
    ligand_atoms = [atom for residue in ligand_residues for atom in residue.get_atoms()]
    if not ligand_atoms:
        return _invalid_structure_result(
            job,
            "ligand_selection",
            f"ligand {ligand_resname} has no atoms",
            tuple(qc),
        )

    ligand_object_ids = {id(residue) for residue in ligand_residues}
    protein_residues = [
        residue
        for chain in model
        if job.protein_chains is None or str(chain.id) in job.protein_chains
        for residue in chain
        if id(residue) not in ligand_object_ids
        and not residue.get_id()[0].startswith("W")
        and is_aa(residue, standard=False)
    ]
    protein_atoms = [atom for residue in protein_residues for atom in residue.get_atoms()]
    if not protein_atoms:
        return _invalid_structure_result(
            job, "protein_selection", "no protein atoms found", tuple(qc)
        )

    neighbor_search = NeighborSearch(protein_atoms)
    minimum_by_residue: dict[ResidueKey, float] = {}
    for ligand_atom in ligand_atoms:
        for protein_atom in neighbor_search.search(ligand_atom.get_coord(), cutoff):
            key = residue_key(protein_atom.get_parent())
            distance = float(
                np.linalg.norm(ligand_atom.get_coord() - protein_atom.get_coord())
            )
            previous = minimum_by_residue.get(key)
            if previous is None or distance < previous:
                minimum_by_residue[key] = distance

    contacts = []
    for protein_residue in protein_residues:
        key = residue_key(protein_residue)
        minimum = minimum_by_residue.get(key)
        contacts.append(
            ContactRecord(
                structure_id=job.structure_id,
                compound=job.compound,
                residue=key,
                is_contact=minimum is not None,
                min_contact_distance_angstrom=minimum,
                capped_distance_angstrom=minimum if minimum is not None else cutoff,
            )
        )
    contacts.sort(key=lambda record: record.residue)

    minimum_pair_distance = float("inf")
    clash_pair_count = 0
    protein_coordinates = np.asarray(
        [atom.get_coord() for atom in protein_atoms], dtype=float
    )
    for ligand_atom in ligand_atoms:
        distances = np.linalg.norm(
            protein_coordinates - ligand_atom.get_coord(), axis=1
        )
        minimum_pair_distance = min(
            minimum_pair_distance, float(np.min(distances))
        )
        clash_pair_count += int(np.count_nonzero(distances <= clash_cutoff))

    buried_sasa = None
    interface_area = None
    try:
        buried_sasa, interface_area = calculate_buried_sasa(
            protein_residues, ligand_residues
        )
    except Exception as exc:
        qc.append(
            QCRecord(
                severity="warning",
                stage="sasa",
                message=f"failed to calculate buried SASA: {exc}",
                compound=job.compound,
                structure_id=job.structure_id,
            )
        )

    metrics = StructureMetrics(
        ligand_resname=ligand_resname,
        n_protein_residues=len(protein_residues),
        n_contact_residues=sum(record.is_contact for record in contacts),
        minimum_pair_distance_angstrom=minimum_pair_distance,
        clash_pair_count=clash_pair_count,
        buried_sasa_total_A2=buried_sasa,
        interface_area_A2=interface_area,
    )
    return StructureResult(
        job=job,
        status="valid",
        contacts=tuple(contacts),
        metrics=metrics,
        qc=tuple(qc),
    )


def summarize_values(values) -> dict[str, float | int | None]:
    numeric = [float(value) for value in values if value is not None]
    if not numeric:
        return {
            "count": 0,
            "mean": None,
            "std": None,
            "median": None,
            "min": None,
            "max": None,
        }
    return {
        "count": len(numeric),
        "mean": statistics.fmean(numeric),
        "std": statistics.pstdev(numeric),
        "median": statistics.median(numeric),
        "min": min(numeric),
        "max": max(numeric),
    }


def contact_set_jaccard(left: set, right: set) -> float:
    if not left and not right:
        return 1.0
    return len(left & right) / len(left | right)


def _residue_label(key: ResidueKey) -> str:
    return f"{key.resname} {key.chain}:{key.resseq}{key.icode}"


def _residue_summary_row(
    key: ResidueKey,
    records: list[ContactRecord],
    n_valid_structures: int,
    compound: str | None,
) -> dict[str, object]:
    values = [record.capped_distance_angstrom for record in records]
    distance_stats = summarize_values(values)
    n_contacts = sum(record.is_contact for record in records)
    row: dict[str, object] = {}
    if compound is not None:
        row["compound"] = compound
    row.update(
        {
            "chain": key.chain,
            "resseq": key.resseq,
            "icode": key.icode,
            "resname": key.resname,
            "label": _residue_label(key),
            "n_valid_structures": n_valid_structures,
            "n_present": len(records),
            "n_contacts": n_contacts,
            "contact_frequency": n_contacts / len(records),
            "mean_capped_distance": distance_stats["mean"],
            "std_capped_distance": distance_stats["std"],
            "median_capped_distance": distance_stats["median"],
            "min_capped_distance": distance_stats["min"],
            "max_capped_distance": distance_stats["max"],
        }
    )
    return row


def aggregate_residues(
    results: list[StructureResult] | tuple[StructureResult, ...],
) -> list[dict[str, object]]:
    valid_counts: dict[str, int] = defaultdict(int)
    groups: dict[tuple[str, ResidueKey], list[ContactRecord]] = defaultdict(list)
    for result in results:
        if result.status != "valid":
            continue
        valid_counts[result.job.compound] += 1
        for record in result.contacts:
            groups[(result.job.compound, record.residue)].append(record)

    rows = [
        _residue_summary_row(
            key=residue,
            records=records,
            n_valid_structures=valid_counts[compound],
            compound=compound,
        )
        for (compound, residue), records in groups.items()
    ]
    return sorted(
        rows,
        key=lambda row: (
            str(row["compound"]),
            str(row["chain"]),
            int(row["resseq"]),
            str(row["icode"]),
            str(row["resname"]),
        ),
    )


def aggregate_global_residues(
    results: list[StructureResult] | tuple[StructureResult, ...],
) -> list[dict[str, object]]:
    valid_results = [result for result in results if result.status == "valid"]
    groups: dict[ResidueKey, list[ContactRecord]] = defaultdict(list)
    for result in valid_results:
        for record in result.contacts:
            groups[record.residue].append(record)

    rows = [
        _residue_summary_row(
            key=residue,
            records=records,
            n_valid_structures=len(valid_results),
            compound=None,
        )
        for residue, records in groups.items()
    ]
    return sorted(
        rows,
        key=lambda row: (
            str(row["chain"]),
            int(row["resseq"]),
            str(row["icode"]),
            str(row["resname"]),
        ),
    )


def aggregate_compounds(
    results: list[StructureResult] | tuple[StructureResult, ...],
) -> list[dict[str, object]]:
    grouped: dict[str, list[StructureResult]] = defaultdict(list)
    for result in results:
        grouped[result.job.compound].append(result)

    rows = []
    for compound in sorted(grouped):
        compound_results = grouped[compound]
        valid = [result for result in compound_results if result.status == "valid"]
        contact_stats = summarize_values(
            [result.metrics.n_contact_residues for result in valid]
        )
        minimum_stats = summarize_values(
            [result.metrics.minimum_pair_distance_angstrom for result in valid]
        )
        clash_stats = summarize_values(
            [result.metrics.clash_pair_count for result in valid]
        )
        buried_stats = summarize_values(
            [result.metrics.buried_sasa_total_A2 for result in valid]
        )
        interface_stats = summarize_values(
            [result.metrics.interface_area_A2 for result in valid]
        )
        contact_sets = [
            {record.residue for record in result.contacts if record.is_contact}
            for result in valid
        ]
        pairwise = [
            contact_set_jaccard(left, right)
            for left, right in combinations(contact_sets, 2)
        ]
        rows.append(
            {
                "compound": compound,
                "n_discovered": len(compound_results),
                "n_valid": len(valid),
                "n_invalid": len(compound_results) - len(valid),
                "n_warnings": sum(
                    record.severity == "warning"
                    for result in compound_results
                    for record in result.qc
                ),
                "mean_contact_residues": contact_stats["mean"],
                "std_contact_residues": contact_stats["std"],
                "mean_minimum_pair_distance": minimum_stats["mean"],
                "std_minimum_pair_distance": minimum_stats["std"],
                "mean_clash_pair_count": clash_stats["mean"],
                "std_clash_pair_count": clash_stats["std"],
                "mean_buried_sasa_total_A2": buried_stats["mean"],
                "std_buried_sasa_total_A2": buried_stats["std"],
                "mean_interface_area_A2": interface_stats["mean"],
                "std_interface_area_A2": interface_stats["std"],
                "mean_pairwise_contact_jaccard": (
                    statistics.fmean(pairwise) if pairwise else None
                ),
            }
        )
    return rows


def _read_tsv_rows(
    path: Path, required_columns: set[str]
) -> list[tuple[int, dict[str, str]]]:
    if not path.is_file():
        raise PipelineError(f"TSV input is not a file: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or ())
        missing = sorted(required_columns - fieldnames)
        if missing:
            raise PipelineError(
                f"{path.name} is missing required columns: {', '.join(missing)}"
            )
        rows = [
            (
                row_number,
                {key: (value or "").strip() for key, value in row.items()},
            )
            for row_number, row in enumerate(reader, start=2)
        ]
    if not rows:
        raise PipelineError(f"{path.name} contains no data rows")
    return rows


def _finite_float(value: str, label: str, row_number: int) -> float:
    try:
        number = float(value)
    except ValueError as exc:
        raise PipelineError(
            f"row {row_number}: {label} must be numeric, got {value!r}"
        ) from exc
    if not math.isfinite(number):
        raise PipelineError(f"row {row_number}: {label} must be finite")
    return number


def bootstrap_mean_ci(values: list[float]) -> tuple[float | None, float | None]:
    if len(values) < 2:
        return None, None
    generator = np.random.default_rng(20260901)
    array = np.asarray(values, dtype=float)
    indices = generator.integers(0, len(array), size=(10_000, len(array)))
    means = array[indices].mean(axis=1)
    low, high = np.percentile(means, [2.5, 97.5])
    return float(low), float(high)


def load_energy_results(
    path: Path,
    known_compounds: set[str],
    known_structure_ids: set[str],
) -> list[EnergyValue]:
    rows = _read_tsv_rows(path, {"compound", "method", "value_kcal_mol"})
    values: list[EnergyValue] = []
    identities: set[tuple[str, str, str, str, str]] = set()
    rbfe_references: dict[tuple[str, str], set[str]] = defaultdict(set)

    for row_number, row in rows:
        compound = row["compound"]
        method = row["method"]
        if not compound:
            raise PipelineError(f"row {row_number}: compound must not be empty")
        if compound not in known_compounds:
            raise PipelineError(f"row {row_number}: unknown compound {compound!r}")
        if not method:
            raise PipelineError(f"row {row_number}: method must not be empty")
        value = _finite_float(row["value_kcal_mol"], "value_kcal_mol", row_number)
        structure_id = row.get("structure_id", "")
        if structure_id and structure_id not in known_structure_ids:
            raise PipelineError(
                f"row {row_number}: unknown structure_id {structure_id!r}"
            )

        uncertainty_text = row.get("uncertainty_kcal_mol", "")
        uncertainty = (
            _finite_float(uncertainty_text, "uncertainty_kcal_mol", row_number)
            if uncertainty_text
            else None
        )
        if uncertainty is not None and uncertainty < 0:
            raise PipelineError(
                f"row {row_number}: uncertainty_kcal_mol must not be negative"
            )
        reference = row.get("reference_compound", "")
        if reference and reference not in known_compounds:
            raise PipelineError(
                f"row {row_number}: unknown reference compound {reference!r}"
            )
        replicate = row.get("replicate", "")
        identity = (compound, method, structure_id, reference, replicate)
        if identity in identities:
            raise PipelineError(
                f"row {row_number}: duplicate energy row identity for {identity!r}"
            )
        identities.add(identity)

        if method.lower() == "rbfe":
            if not reference:
                raise PipelineError(
                    f"row {row_number}: rbfe requires reference_compound"
                )
            rbfe_references[(compound, method)].add(reference)

        values.append(
            EnergyValue(
                compound=compound,
                method=method,
                value_kcal_mol=value,
                structure_id=structure_id,
                uncertainty_kcal_mol=uncertainty,
                reference_compound=reference,
                replicate=replicate,
            )
        )

    for (compound, method), references in rbfe_references.items():
        if len(references) > 1:
            raise PipelineError(
                f"{compound}/{method} has mixed reference compounds: "
                + ", ".join(sorted(references))
            )
    return sorted(
        values,
        key=lambda item: (
            item.compound,
            item.method,
            item.reference_compound,
            item.structure_id,
            item.replicate,
        ),
    )


def aggregate_energy_values(
    values: list[EnergyValue] | tuple[EnergyValue, ...],
) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str], list[EnergyValue]] = defaultdict(list)
    for value in values:
        groups[(value.compound, value.method, value.reference_compound)].append(
            value
        )

    rows = []
    for (compound, method, reference), group in sorted(groups.items()):
        numeric = [value.value_kcal_mol for value in group]
        stats = summarize_values(numeric)
        ci_low, ci_high = bootstrap_mean_ci(numeric)
        uncertainty_stats = summarize_values(
            [value.uncertainty_kcal_mol for value in group]
        )
        rows.append(
            {
                "compound": compound,
                "method": method,
                "reference_compound": reference,
                "n_values": len(group),
                "mean": stats["mean"],
                "std": stats["std"],
                "median": stats["median"],
                "min": stats["min"],
                "max": stats["max"],
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "mean_reported_uncertainty": uncertainty_stats["mean"],
            }
        )
    return rows


def load_affinity_results(
    path: Path, known_compounds: set[str], convert_ki: bool = False
) -> list[AffinityValue]:
    rows = _read_tsv_rows(path, {"compound", "metric", "value", "unit"})
    metric_names = {"kd": "Kd", "ki": "Ki", "ic50": "IC50"}
    values: list[AffinityValue] = []
    for row_number, row in rows:
        compound = row["compound"]
        if compound not in known_compounds:
            raise PipelineError(f"row {row_number}: unknown compound {compound!r}")
        metric_key = row["metric"].lower()
        if metric_key not in metric_names:
            raise PipelineError(
                f"row {row_number}: metric must be Kd, Ki, or IC50"
            )
        metric = metric_names[metric_key]
        number = _finite_float(row["value"], "value", row_number)
        if number <= 0:
            raise PipelineError(f"row {row_number}: affinity value must be positive")
        unit = row["unit"]
        if unit not in UNIT_TO_MOLAR:
            raise PipelineError(
                f"row {row_number}: unsupported affinity unit {unit!r}"
            )
        temperature_text = row.get("temperature_K", "")
        temperature = (
            _finite_float(temperature_text, "temperature_K", row_number)
            if temperature_text
            else 298.15
        )
        if temperature <= 0:
            raise PipelineError(
                f"row {row_number}: temperature_K must be positive"
            )
        concentration = number * UNIT_TO_MOLAR[unit]
        convert = metric == "Kd" or (metric == "Ki" and convert_ki)
        delta_g = (
            R_KCAL_PER_MOL_K * temperature * math.log(concentration)
            if convert
            else None
        )
        values.append(
            AffinityValue(
                compound=compound,
                metric=metric,
                value=number,
                unit=unit,
                concentration_M=concentration,
                temperature_K=temperature,
                replicate=row.get("replicate", ""),
                delta_g_standard_kcal_mol=delta_g,
            )
        )
    return values


def aggregate_affinity_values(
    values: list[AffinityValue] | tuple[AffinityValue, ...],
) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, float], list[AffinityValue]] = defaultdict(list)
    for value in values:
        groups[(value.compound, value.metric, value.temperature_K)].append(value)

    rows = []
    for (compound, metric, temperature), group in sorted(groups.items()):
        concentration_stats = summarize_values(
            [value.concentration_M for value in group]
        )
        delta_stats = summarize_values(
            [value.delta_g_standard_kcal_mol for value in group]
        )
        delta_numeric = [
            value.delta_g_standard_kcal_mol
            for value in group
            if value.delta_g_standard_kcal_mol is not None
        ]
        ci_low, ci_high = bootstrap_mean_ci(delta_numeric)
        rows.append(
            {
                "compound": compound,
                "metric": metric,
                "temperature_K": temperature,
                "n_values": len(group),
                "mean_concentration_M": concentration_stats["mean"],
                "std_concentration_M": concentration_stats["std"],
                "median_concentration_M": concentration_stats["median"],
                "min_concentration_M": concentration_stats["min"],
                "max_concentration_M": concentration_stats["max"],
                "mean_delta_g_standard_kcal_mol": delta_stats["mean"],
                "std_delta_g_standard_kcal_mol": delta_stats["std"],
                "delta_g_ci95_low": ci_low,
                "delta_g_ci95_high": ci_high,
            }
        )
    return rows


def format_tsv_value(value: object) -> str:
    """Return a deterministic, spreadsheet-friendly TSV cell value."""
    if value is None:
        return ""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (float, np.floating)):
        numeric = float(value)
        if not math.isfinite(numeric):
            return ""
        return f"{numeric:.6f}"
    return str(value)


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    """Write rows as a deterministic tab-separated table."""
    if fieldnames is None and not rows:
        raise PipelineError(f"cannot infer TSV columns for empty table: {path.name}")

    columns = list(rows[0]) if fieldnames is None else list(fieldnames)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {name: format_tsv_value(row.get(name)) for name in columns}
            )


def _validated_inventory_paths(
    output_dir: Path, relative_paths: list[str]
) -> list[Path]:
    """Resolve inventory entries while ensuring each remains below output_dir."""
    resolved_output = output_dir.resolve()
    validated: list[Path] = []
    for relative_path in relative_paths:
        relative = Path(relative_path)
        if (
            not relative_path
            or relative.is_absolute()
            or ".." in relative.parts
            or not relative.parts
        ):
            raise PipelineError(
                f"inventory path is outside output directory: {relative_path}"
            )
        unresolved_candidate = resolved_output / relative
        resolved_parent = unresolved_candidate.parent.resolve()
        try:
            resolved_parent.relative_to(resolved_output)
        except ValueError as exc:
            raise PipelineError(
                f"inventory path is outside output directory: {relative_path}"
            ) from exc
        candidate = resolved_parent / unresolved_candidate.name
        if candidate == resolved_output:
            raise PipelineError(
                f"inventory path is outside output directory: {relative_path}"
            )
        validated.append(candidate)
    return validated


def prepare_output_directory(output_dir: Path, overwrite: bool) -> None:
    """Create output_dir or safely remove only files from its prior inventory."""
    output_dir = Path(output_dir)
    if not output_dir.exists():
        output_dir.mkdir(parents=True)
        return
    if output_dir.is_symlink() or not output_dir.is_dir():
        raise PipelineError(f"output path is not a directory: {output_dir}")

    if not any(output_dir.iterdir()):
        return
    if not overwrite:
        raise PipelineError(
            f"output directory is not empty; pass --overwrite to reuse it: {output_dir}"
        )

    inventory = output_dir / "generated_files.json"
    if not inventory.exists():
        return
    if inventory.is_symlink() or not inventory.is_file():
        raise PipelineError(f"invalid generated-file inventory: {inventory}")

    try:
        payload = json.loads(inventory.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise PipelineError(
            f"could not read generated-file inventory: {inventory}"
        ) from exc
    relative_paths = payload.get("files") if isinstance(payload, dict) else None
    if not isinstance(relative_paths, list) or not all(
        isinstance(value, str) for value in relative_paths
    ):
        raise PipelineError(
            f"generated-file inventory must contain a string list named 'files': {inventory}"
        )

    candidates = _validated_inventory_paths(output_dir, relative_paths)

    # Validate the full inventory before deleting any entry.
    for candidate in candidates:
        if (
            candidate.exists()
            and candidate.is_dir()
            and not candidate.is_symlink()
        ):
            raise PipelineError(
                f"generated-file inventory cannot list a directory: {candidate}"
            )

    for candidate in candidates:
        if candidate.exists() or candidate.is_symlink():
            candidate.unlink()


def claim_output_path(path: Path, output_dir: Path) -> Path:
    """Validate that a new output is contained and will not replace a file."""
    path = Path(path)
    resolved_output = Path(output_dir).resolve()
    if path.is_symlink():
        raise PipelineError(
            f"output path was not listed in the preceding inventory: {path}"
        )
    candidate = path.resolve()
    try:
        candidate.relative_to(resolved_output)
    except ValueError as exc:
        raise PipelineError(f"output path is outside output directory: {path}") from exc
    if candidate == resolved_output:
        raise PipelineError(f"output path is outside output directory: {path}")
    if candidate.exists():
        raise PipelineError(
            f"output path was not listed in the preceding inventory: {path}"
        )
    return candidate


def write_generated_inventory(
    output_dir: Path, generated_paths: list[Path]
) -> Path:
    """Write the sorted set of generated paths used by a future safe overwrite."""
    resolved_output = Path(output_dir).resolve()
    relative_paths: set[str] = {"generated_files.json"}
    for path in generated_paths:
        candidate = Path(path).resolve()
        try:
            relative = candidate.relative_to(resolved_output)
        except ValueError as exc:
            raise PipelineError(
                f"generated path is outside output directory: {path}"
            ) from exc
        if candidate == resolved_output:
            raise PipelineError(
                f"generated path is outside output directory: {path}"
            )
        relative_paths.add(relative.as_posix())

    inventory = claim_output_path(
        resolved_output / "generated_files.json", resolved_output
    )
    inventory.write_text(
        json.dumps({"files": sorted(relative_paths)}, indent=2) + "\n",
        encoding="utf-8",
    )
    return inventory


def jobs_to_rows(
    jobs: list[StructureJob],
    minimization_records: list[MinimizationRecord] | None = None,
    analysis_source: str = "original",
) -> list[dict[str, object]]:
    """Serialize resolved structure assignments using a fixed schema."""
    record_by_id = {
        record.job.structure_id: record
        for record in (minimization_records or [])
    }
    return [
        {
            "structure_id": job.structure_id,
            "dataset_dir": job.dataset_rel,
            "compound": job.compound,
            "ligand_resname": job.ligand_resname,
            "protein_chains": (
                ",".join(job.protein_chains) if job.protein_chains else None
            ),
            "ligand_chain": job.ligand_chain,
            "cif_path": job.cif_path.as_posix(),
            "analysis_source": analysis_source,
            "original_cif_path": job.cif_path.as_posix(),
            "analysis_cif_path": (
                record_by_id[job.structure_id].output_relative_path.as_posix()
                if analysis_source == "minimized"
                and job.structure_id in record_by_id
                and record_by_id[job.structure_id].status == "accepted"
                and record_by_id[job.structure_id].output_relative_path is not None
                else job.cif_path.as_posix()
                if analysis_source == "original"
                else None
            ),
        }
        for job in jobs
    ]


def contacts_to_rows(
    results: list[StructureResult] | tuple[StructureResult, ...],
) -> list[dict[str, object]]:
    """Flatten per-residue structure contacts using a fixed schema."""
    rows: list[dict[str, object]] = []
    for result in results:
        for contact in result.contacts:
            residue = contact.residue
            rows.append(
                {
                    "structure_id": contact.structure_id,
                    "compound": contact.compound,
                    "chain": residue.chain,
                    "resseq": residue.resseq,
                    "icode": residue.icode,
                    "resname": residue.resname,
                    "is_contact": contact.is_contact,
                    "min_contact_distance_angstrom": (
                        contact.min_contact_distance_angstrom
                    ),
                    "capped_distance_angstrom": contact.capped_distance_angstrom,
                }
            )
    return rows


def structure_results_to_rows(
    results: list[StructureResult] | tuple[StructureResult, ...],
    analysis_source: str = "original",
) -> list[dict[str, object]]:
    """Serialize one structure-level QC/geometry row per analyzed CIF."""
    rows: list[dict[str, object]] = []
    for result in results:
        metrics = result.metrics
        rows.append(
            {
                "structure_id": result.job.structure_id,
                "compound": result.job.compound,
                "status": result.status,
                "analysis_source": analysis_source,
                "ligand_resname": metrics.ligand_resname,
                "n_protein_residues": metrics.n_protein_residues,
                "n_contact_residues": metrics.n_contact_residues,
                "minimum_pair_distance_angstrom": (
                    metrics.minimum_pair_distance_angstrom
                ),
                "clash_pair_count": metrics.clash_pair_count,
                "buried_sasa_total_A2": metrics.buried_sasa_total_A2,
                "interface_area_A2": metrics.interface_area_A2,
            }
        )
    return rows


def minimization_records_to_rows(
    records: list[MinimizationRecord] | tuple[MinimizationRecord, ...],
    original_results: list[StructureResult] | tuple[StructureResult, ...],
) -> list[dict[str, object]]:
    """Serialize minimization protocol, outcome, and same-system QC values."""
    original_by_id = {
        result.job.structure_id: result for result in original_results
    }
    rows: list[dict[str, object]] = []
    for record in records:
        original = original_by_id[record.job.structure_id]
        candidate = record.candidate_result
        config = record.config
        initial = record.initial_potential_kcal_mol
        final = record.final_potential_kcal_mol
        rows.append(
            {
                "structure_id": record.job.structure_id,
                "compound": record.job.compound,
                "status": record.status,
                "reason": "; ".join(record.reasons),
                "protocol": config.protocol,
                "protein_forcefield": config.protein_forcefield,
                "ligand_forcefield": config.ligand_forcefield,
                "platform": config.platform,
                "device_index": config.device_index,
                "ph": config.ph,
                "pocket_radius_A": config.pocket_radius_A,
                "tolerance_kj_mol_nm": config.tolerance_kj_mol_nm,
                "output_cif": (
                    record.output_relative_path.as_posix()
                    if record.output_relative_path is not None
                    else None
                ),
                "initial_potential_kcal_mol": initial,
                "final_potential_kcal_mol": final,
                "delta_potential_kcal_mol": (
                    final - initial
                    if initial is not None and final is not None
                    else None
                ),
                "backbone_rmsd_A": record.backbone_rmsd_A,
                "ligand_rmsd_A": record.ligand_rmsd_A,
                "before_minimum_pair_distance_A": (
                    original.metrics.minimum_pair_distance_angstrom
                ),
                "after_minimum_pair_distance_A": (
                    candidate.metrics.minimum_pair_distance_angstrom
                    if candidate is not None
                    else None
                ),
                "before_fixed_clash_pair_count": (
                    original.metrics.clash_pair_count
                ),
                "after_fixed_clash_pair_count": (
                    candidate.metrics.clash_pair_count
                    if candidate is not None
                    else None
                ),
                "before_vdw_clash_pair_count": record.before_vdw_clash_pair_count,
                "after_vdw_clash_pair_count": record.after_vdw_clash_pair_count,
                "before_maximum_vdw_overlap_A": (
                    record.before_maximum_vdw_overlap_A
                ),
                "after_maximum_vdw_overlap_A": record.after_maximum_vdw_overlap_A,
                "before_contact_residue_count": (
                    original.metrics.n_contact_residues
                ),
                "after_contact_residue_count": (
                    candidate.metrics.n_contact_residues
                    if candidate is not None
                    else None
                ),
                "contact_set_jaccard": record.contact_set_jaccard,
                "package_versions": ";".join(
                    f"{name}={version}"
                    for name, version in sorted(record.package_versions)
                ),
            }
        )
    return rows


def qc_to_rows(
    records: list[QCRecord]
    | tuple[QCRecord, ...]
    | list[StructureResult]
    | tuple[StructureResult, ...],
) -> list[dict[str, object]]:
    """Serialize QC records, accepting either records or structure results."""
    flattened: list[QCRecord] = []
    for value in records:
        if isinstance(value, StructureResult):
            flattened.extend(value.qc)
        else:
            flattened.append(value)
    return [
        {
            "severity": record.severity,
            "stage": record.stage,
            "compound": record.compound,
            "structure_id": record.structure_id,
            "message": record.message,
        }
        for record in flattened
    ]


def safe_filename_map(values) -> dict[str, str]:
    """Map arbitrary labels to unique, deterministic filename components."""
    result: dict[str, str] = {}
    used: set[str] = set()
    for value in sorted({str(value) for value in values}):
        base = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._") or "item"
        candidate = base
        if candidate in used:
            suffix = hashlib.sha1(value.encode("utf-8")).hexdigest()[:8]
            candidate = f"{base}_{suffix}"
            counter = 2
            while candidate in used:
                candidate = f"{base}_{suffix}_{counter}"
                counter += 1
        used.add(candidate)
        result[value] = candidate
    return result


def build_heatmap_matrix(
    rows: list[dict[str, object]], value_key: str
) -> tuple[np.ndarray, list[str], list[ResidueKey]]:
    """Build a compound-by-residue matrix, retaining absence as NaN."""
    compounds = sorted({str(row["compound"]) for row in rows})
    residues = sorted(
        {
            ResidueKey(
                str(row["chain"]),
                int(row["resseq"]),
                str(row["icode"]),
                str(row["resname"]),
            )
            for row in rows
        }
    )
    compound_index = {value: index for index, value in enumerate(compounds)}
    residue_index = {value: index for index, value in enumerate(residues)}
    matrix = np.full((len(compounds), len(residues)), np.nan, dtype=float)
    for row in rows:
        residue = ResidueKey(
            str(row["chain"]),
            int(row["resseq"]),
            str(row["icode"]),
            str(row["resname"]),
        )
        value = row.get(value_key)
        if value is not None:
            matrix[
                compound_index[str(row["compound"])], residue_index[residue]
            ] = float(value)
    return matrix, compounds, residues


def _save_figure(
    figure: plt.Figure, path: Path, output_root: Path, *, dpi: int = 200
) -> Path:
    claimed = claim_output_path(path, output_root)
    claimed.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(claimed, dpi=dpi, bbox_inches="tight")
    return claimed


def _plot_root(plots_dir: Path, output_root: Path | None) -> Path:
    return Path(plots_dir) if output_root is None else Path(output_root)


def _staggered_contact_labels(
    axis: plt.Axes,
    rows: list[dict[str, object]],
    label_frequency: float,
) -> None:
    candidates = [
        row
        for row in rows
        if float(row["contact_frequency"]) >= label_frequency
    ]
    for index, row in enumerate(candidates):
        label = f"{row['resname']}{row['resseq']}{row['icode']}"
        axis.annotate(
            label,
            (int(row["resseq"]), float(row["contact_frequency"])),
            textcoords="offset points",
            xytext=(0, 6 + 10 * (index % 4)),
            ha="center",
            va="bottom",
            rotation=45,
            fontsize=6.5,
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": "white",
                "edgecolor": "none",
                "alpha": 0.75,
            },
        )


def plot_compound_chain_contacts(
    rows: list[dict[str, object]],
    plots_dir: Path,
    cutoff: float,
    label_frequency: float,
    output_root: Path | None = None,
) -> list[Path]:
    """Plot capped distance and contact frequency for every compound/chain."""
    grouped: dict[tuple[str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["compound"]), str(row["chain"]))].append(row)
    logical_stems = {
        key: f"{key[0]}_chain_{key[1]}_contacts" for key in grouped
    }
    filename_map = safe_filename_map(logical_stems.values())
    root = _plot_root(plots_dir, output_root)
    generated: list[Path] = []

    for key in sorted(grouped):
        compound, chain = key
        chain_rows = sorted(
            grouped[key],
            key=lambda row: (
                int(row["resseq"]),
                str(row["icode"]),
                str(row["resname"]),
            ),
        )
        x = np.asarray([int(row["resseq"]) for row in chain_rows])
        means = np.asarray(
            [float(row["mean_capped_distance"]) for row in chain_rows]
        )
        standard_deviations = np.asarray(
            [float(row["std_capped_distance"]) for row in chain_rows]
        )
        frequencies = np.asarray(
            [float(row["contact_frequency"]) for row in chain_rows]
        )

        figure, (distance_axis, frequency_axis) = plt.subplots(
            2,
            1,
            figsize=(12, 8),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 1.5]},
        )
        try:
            distance_axis.errorbar(
                x,
                means,
                yerr=standard_deviations,
                fmt="o-",
                markersize=3,
                linewidth=1,
                capsize=2,
                color="tab:blue",
                ecolor="0.55",
                label="Mean capped distance ± population SD",
            )
            distance_axis.axhline(
                cutoff,
                color="tab:red",
                linestyle="--",
                linewidth=1,
                alpha=0.7,
                label=f"Contact cutoff ({cutoff:g} Å)",
            )
            distance_axis.set_ylabel("Capped ligand distance (Å)")
            distance_axis.set_title(
                f"{compound}: ligand-contact geometry, chain {chain}"
            )
            distance_axis.grid(alpha=0.25)
            distance_axis.legend(frameon=False)

            frequency_axis.plot(
                x,
                frequencies,
                "s-",
                markersize=3,
                linewidth=1,
                color="tab:orange",
            )
            frequency_axis.set_xlabel("Residue number")
            frequency_axis.set_ylabel("Contact frequency")
            frequency_axis.set_ylim(-0.02, 1.08)
            frequency_axis.grid(alpha=0.25)
            _staggered_contact_labels(
                frequency_axis, chain_rows, label_frequency
            )
            figure.tight_layout()
            filename = filename_map[logical_stems[key]] + ".png"
            generated.append(
                _save_figure(figure, Path(plots_dir) / filename, root)
            )
        finally:
            plt.close(figure)
    return generated


def _heatmap_ticks(
    axis: plt.Axes, residues: list[ResidueKey], maximum: int = 40
) -> None:
    if not residues:
        return
    step = max(1, math.ceil(len(residues) / maximum))
    positions = list(range(0, len(residues), step))
    labels = [
        f"{residues[index].resname}{residues[index].resseq}"
        f"{residues[index].icode}"
        for index in positions
    ]
    axis.set_xticks(positions, labels, rotation=60, ha="right", fontsize=7)


def plot_chain_heatmaps(
    rows: list[dict[str, object]],
    plots_dir: Path,
    cutoff: float,
    output_root: Path | None = None,
) -> list[Path]:
    """Plot masked compound-by-residue frequency and distance heatmaps."""
    chains = sorted({str(row["chain"]) for row in rows})
    logical_stems = [
        f"{metric}_heatmap_chain_{chain}"
        for chain in chains
        for metric in ("contact_frequency", "capped_distance")
    ]
    filename_map = safe_filename_map(logical_stems)
    root = _plot_root(plots_dir, output_root)
    generated: list[Path] = []
    configurations = (
        (
            "contact_frequency",
            "contact_frequency",
            "Contact frequency",
            "Oranges",
            0.0,
            1.0,
        ),
        (
            "capped_distance",
            "mean_capped_distance",
            "Mean capped distance (Å)",
            "Blues_r",
            0.0,
            cutoff,
        ),
    )

    for chain in chains:
        chain_rows = [row for row in rows if str(row["chain"]) == chain]
        for stem_prefix, value_key, colorbar_label, color_map, minimum, maximum in configurations:
            matrix, compounds, residues = build_heatmap_matrix(
                chain_rows, value_key
            )
            width = max(7.0, min(20.0, 0.24 * len(residues) + 3.0))
            height = max(3.0, min(12.0, 0.4 * len(compounds) + 2.2))
            figure, axis = plt.subplots(figsize=(width, height))
            try:
                cmap = plt.get_cmap(color_map).with_extremes(bad="lightgray")
                image = axis.imshow(
                    np.ma.masked_invalid(matrix),
                    aspect="auto",
                    interpolation="nearest",
                    cmap=cmap,
                    vmin=minimum,
                    vmax=maximum,
                )
                axis.set_yticks(range(len(compounds)), compounds)
                _heatmap_ticks(axis, residues)
                axis.set_xlabel("Residue")
                axis.set_ylabel("Compound")
                axis.set_title(f"{colorbar_label}, chain {chain}")
                figure.colorbar(image, ax=axis, label=colorbar_label)
                figure.tight_layout()
                logical_stem = f"{stem_prefix}_heatmap_chain_{chain}"
                filename = filename_map[logical_stem] + ".png"
                generated.append(
                    _save_figure(figure, Path(plots_dir) / filename, root)
                )
            finally:
                plt.close(figure)
    return generated


def _numeric_or_nan(value: object) -> float:
    return np.nan if value is None else float(value)


def plot_structural_summary(
    rows: list[dict[str, object]],
    plots_dir: Path,
    output_root: Path | None = None,
) -> Path:
    """Plot structural evidence without presenting it as an affinity ranking."""
    ordered = sorted(rows, key=lambda row: str(row["compound"]))
    compounds = [str(row["compound"]) for row in ordered]
    x = np.arange(len(compounds))
    figure, axes = plt.subplots(
        4,
        1,
        figsize=(max(8.0, 0.65 * len(compounds) + 3.5), 12),
        sharex=True,
    )
    try:
        axes[0].bar(x, [int(row["n_valid"]) for row in ordered], color="tab:blue")
        axes[0].set_ylabel("Valid seeds")

        axes[1].bar(
            x,
            [
                _numeric_or_nan(row.get("mean_pairwise_contact_jaccard"))
                for row in ordered
            ],
            color="tab:green",
        )
        axes[1].set_ylabel("Contact-set\nJaccard")
        axes[1].set_ylim(0, 1.05)

        axes[2].bar(
            x,
            [
                _numeric_or_nan(row.get("mean_buried_sasa_total_A2"))
                for row in ordered
            ],
            yerr=[
                0.0
                if row.get("std_buried_sasa_total_A2") is None
                else float(row["std_buried_sasa_total_A2"])
                for row in ordered
            ],
            capsize=3,
            color="tab:purple",
        )
        axes[2].set_ylabel("Buried SASA\n(Å²)")

        axes[3].bar(
            x,
            [
                _numeric_or_nan(row.get("mean_clash_pair_count"))
                for row in ordered
            ],
            yerr=[
                0.0
                if row.get("std_clash_pair_count") is None
                else float(row["std_clash_pair_count"])
                for row in ordered
            ],
            capsize=3,
            color="tab:red",
        )
        axes[3].set_ylabel("Clash pairs")
        axes[3].set_xticks(x, compounds, rotation=45, ha="right")
        axes[3].set_xlabel("Compound")
        for axis in axes:
            axis.grid(axis="y", alpha=0.25)
        figure.suptitle(
            "Structural pose evidence (not a binding-affinity estimate)", y=1.01
        )
        figure.tight_layout()
        return _save_figure(
            figure,
            Path(plots_dir) / "compound_structural_summary.png",
            _plot_root(plots_dir, output_root),
        )
    finally:
        plt.close(figure)


def plot_minimization_before_after(
    records: list[MinimizationRecord] | tuple[MinimizationRecord, ...],
    original_results: list[StructureResult] | tuple[StructureResult, ...],
    plots_dir: Path,
    output_root: Path | None = None,
) -> Path:
    """Plot paired geometry QC without presenting potential energy as affinity."""
    completed = [
        record
        for record in records
        if record.candidate_result is not None
        and record.initial_potential_kcal_mol is not None
        and record.final_potential_kcal_mol is not None
    ]
    original_by_id = {
        result.job.structure_id: result for result in original_results
    }
    figure, axes = plt.subplots(
        4,
        1,
        figsize=(max(9.0, 0.55 * len(completed) + 4.0), 13),
        sharex=True,
    )
    try:
        if not completed:
            for axis in axes:
                axis.text(
                    0.5,
                    0.5,
                    "No completed minimizations",
                    ha="center",
                    va="center",
                    transform=axis.transAxes,
                )
        else:
            x = np.arange(len(completed), dtype=float)
            labels = [record.job.structure_id for record in completed]
            colors = [
                "tab:green" if record.status == "accepted" else "tab:orange"
                for record in completed
            ]
            for index, record in enumerate(completed):
                original = original_by_id[record.job.structure_id]
                candidate = record.candidate_result
                axes[0].plot(
                    [index - 0.12, index + 0.12],
                    [
                        original.metrics.clash_pair_count,
                        candidate.metrics.clash_pair_count,
                    ],
                    "o-",
                    color=colors[index],
                )
                axes[1].plot(
                    [index - 0.12, index + 0.12],
                    [
                        record.before_vdw_clash_pair_count,
                        record.after_vdw_clash_pair_count,
                    ],
                    "o-",
                    color=colors[index],
                )

            axes[0].set_ylabel("Fixed-cutoff\nclash pairs")
            axes[1].set_ylabel("van der Waals\nclash pairs")
            width = 0.34
            axes[2].bar(
                x - width / 2,
                [record.backbone_rmsd_A for record in completed],
                width,
                label="Protein backbone",
                color="tab:blue",
            )
            axes[2].bar(
                x + width / 2,
                [record.ligand_rmsd_A for record in completed],
                width,
                label="Ligand heavy atoms",
                color="tab:purple",
            )
            axes[2].axhline(
                completed[0].config.max_backbone_rmsd_A,
                color="tab:blue",
                linestyle="--",
                linewidth=1,
            )
            axes[2].axhline(
                completed[0].config.max_ligand_rmsd_A,
                color="tab:purple",
                linestyle="--",
                linewidth=1,
            )
            axes[2].set_ylabel("Aligned RMSD (Å)")
            axes[2].legend(frameon=False, ncol=2)
            axes[3].bar(
                x,
                [
                    record.final_potential_kcal_mol
                    - record.initial_potential_kcal_mol
                    for record in completed
                ],
                color=colors,
            )
            axes[3].axhline(0.0, color="0.3", linewidth=0.8)
            axes[3].set_ylabel(
                "Same-system potential-energy\nchange (kcal/mol)"
            )
            axes[3].set_xticks(x, labels, rotation=45, ha="right")
            axes[3].set_xlabel("Structure")
        for axis in axes:
            axis.grid(axis="y", alpha=0.25)
        figure.suptitle(
            "Restrained minimization geometry QC (not binding affinity)", y=1.01
        )
        figure.tight_layout()
        return _save_figure(
            figure,
            Path(plots_dir) / "minimization_before_after.png",
            _plot_root(plots_dir, output_root),
        )
    finally:
        plt.close(figure)


def _energy_axis_label(method: str, reference: str) -> str:
    normalized = method.lower()
    if normalized == "rbfe":
        qualifier = f" vs {reference}" if reference else ""
        return f"Relative binding free energy ΔΔG{qualifier} (kcal/mol)"
    if normalized in {"vina", "autodock_vina"}:
        return "Docking score (kcal/mol)"
    if normalized in {"mmgbsa", "mm/gbsa", "mmpbsa", "mm/pbsa"}:
        return f"Estimated {method} energy (kcal/mol)"
    return f"{method} energy-like value (kcal/mol)"


def plot_energy_methods(
    values: list[EnergyValue] | tuple[EnergyValue, ...],
    summary_rows: list[dict[str, object]],
    plots_dir: Path,
    output_root: Path | None = None,
) -> list[Path]:
    """Plot each computational method/reference series without mixing methods."""
    grouped: dict[tuple[str, str], list[EnergyValue]] = defaultdict(list)
    for value in values:
        grouped[(value.method, value.reference_compound)].append(value)
    logical_stems = {
        key: (
            f"energy_{key[0]}_vs_{key[1]}" if key[1] else f"energy_{key[0]}"
        )
        for key in grouped
    }
    filename_map = safe_filename_map(logical_stems.values())
    summaries = {
        (
            str(row["compound"]),
            str(row["method"]),
            str(row["reference_compound"]),
        ): row
        for row in summary_rows
    }
    root = _plot_root(plots_dir, output_root)
    generated: list[Path] = []

    for key in sorted(grouped):
        method, reference = key
        method_values = grouped[key]
        compounds = sorted({value.compound for value in method_values})
        figure, axis = plt.subplots(
            figsize=(max(7.0, 0.75 * len(compounds) + 3.5), 5.5)
        )
        try:
            for position, compound in enumerate(compounds):
                replicates = [
                    value
                    for value in method_values
                    if value.compound == compound
                ]
                jitter = (
                    np.linspace(-0.12, 0.12, len(replicates))
                    if len(replicates) > 1
                    else np.asarray([0.0])
                )
                for offset, value in zip(jitter, replicates):
                    axis.errorbar(
                        position + offset,
                        value.value_kcal_mol,
                        yerr=value.uncertainty_kcal_mol,
                        fmt="o",
                        color="tab:blue",
                        ecolor="0.55",
                        capsize=2,
                        alpha=0.75,
                    )
                summary = summaries[(compound, method, reference)]
                mean = float(summary["mean"])
                standard_deviation = summary.get("std")
                axis.errorbar(
                    position,
                    mean,
                    yerr=(
                        None
                        if standard_deviation is None
                        else float(standard_deviation)
                    ),
                    fmt="D",
                    color="black",
                    capsize=4,
                    markersize=5,
                )
                low = summary.get("ci95_low")
                high = summary.get("ci95_high")
                if low is not None and high is not None:
                    axis.vlines(
                        position,
                        float(low),
                        float(high),
                        color="tab:red",
                        linewidth=3,
                        alpha=0.7,
                    )
            axis.set_xticks(range(len(compounds)), compounds, rotation=45, ha="right")
            axis.set_ylabel(_energy_axis_label(method, reference))
            title = f"Computational energy evidence: {method}"
            if reference:
                title += f" (reference {reference})"
            axis.set_title(title)
            axis.grid(axis="y", alpha=0.25)
            figure.tight_layout()
            filename = filename_map[logical_stems[key]] + ".png"
            generated.append(
                _save_figure(figure, Path(plots_dir) / filename, root)
            )
        finally:
            plt.close(figure)
    return generated


def plot_experimental_affinity(
    values: list[AffinityValue] | tuple[AffinityValue, ...],
    summary_rows: list[dict[str, object]],
    plots_dir: Path,
    output_root: Path | None = None,
) -> Path:
    """Plot measured concentrations and separately derived standard ΔG values."""
    keys = sorted(
        {
            (value.compound, value.metric, value.temperature_K)
            for value in values
        }
    )
    labels = [
        f"{compound}\n{metric}, {temperature:g} K"
        for compound, metric, temperature in keys
    ]
    summaries = {
        (
            str(row["compound"]),
            str(row["metric"]),
            float(row["temperature_K"]),
        ): row
        for row in summary_rows
    }
    figure, (concentration_axis, delta_axis) = plt.subplots(
        2,
        1,
        figsize=(max(8.0, 0.8 * len(keys) + 3.5), 9),
        sharex=True,
    )
    try:
        has_delta_g = False
        for position, key in enumerate(keys):
            group = [
                value
                for value in values
                if (value.compound, value.metric, value.temperature_K) == key
            ]
            jitter = (
                np.linspace(-0.12, 0.12, len(group))
                if len(group) > 1
                else np.asarray([0.0])
            )
            concentration_axis.scatter(
                position + jitter,
                [value.concentration_M for value in group],
                color="tab:blue",
                alpha=0.75,
            )
            summary = summaries[key]
            mean_concentration = float(summary["mean_concentration_M"])
            concentration_std = summary.get("std_concentration_M")
            if concentration_std is None:
                concentration_error = None
            else:
                standard_deviation = float(concentration_std)
                concentration_error = np.asarray(
                    [
                        [min(standard_deviation, mean_concentration * 0.999999)],
                        [standard_deviation],
                    ]
                )
            concentration_axis.errorbar(
                position,
                mean_concentration,
                yerr=concentration_error,
                fmt="D",
                color="black",
                capsize=4,
            )

            delta_values = [
                value.delta_g_standard_kcal_mol
                for value in group
                if value.delta_g_standard_kcal_mol is not None
            ]
            if delta_values:
                has_delta_g = True
                delta_axis.scatter(
                    position + jitter[: len(delta_values)],
                    delta_values,
                    color="tab:purple",
                    alpha=0.75,
                )
                delta_mean = summary.get("mean_delta_g_standard_kcal_mol")
                delta_std = summary.get("std_delta_g_standard_kcal_mol")
                delta_axis.errorbar(
                    position,
                    float(delta_mean),
                    yerr=None if delta_std is None else float(delta_std),
                    fmt="D",
                    color="black",
                    capsize=4,
                )
                low = summary.get("delta_g_ci95_low")
                high = summary.get("delta_g_ci95_high")
                if low is not None and high is not None:
                    delta_axis.vlines(
                        position,
                        float(low),
                        float(high),
                        color="tab:red",
                        linewidth=3,
                        alpha=0.7,
                    )

        concentration_axis.set_yscale("log")
        concentration_axis.set_ylabel("Measured concentration (M; log scale)")
        concentration_axis.set_title("Experimental affinity measurements")
        concentration_axis.grid(axis="y", alpha=0.25, which="both")
        delta_axis.set_ylabel("Derived standard ΔG° (kcal/mol)")
        delta_axis.set_xlabel("Compound and experimental metric")
        delta_axis.set_xticks(range(len(keys)), labels, rotation=45, ha="right")
        delta_axis.grid(axis="y", alpha=0.25)
        if not has_delta_g:
            delta_axis.text(
                0.5,
                0.5,
                "No Kd (or explicitly converted Ki) values",
                transform=delta_axis.transAxes,
                ha="center",
                va="center",
            )
        figure.tight_layout()
        return _save_figure(
            figure,
            Path(plots_dir) / "experimental_affinity.png",
            _plot_root(plots_dir, output_root),
        )
    finally:
        plt.close(figure)


def energy_values_to_rows(
    values: list[EnergyValue] | tuple[EnergyValue, ...],
) -> list[dict[str, object]]:
    """Serialize validated computational-energy evidence."""
    return [
        {
            "compound": value.compound,
            "method": value.method,
            "value_kcal_mol": value.value_kcal_mol,
            "structure_id": value.structure_id,
            "uncertainty_kcal_mol": value.uncertainty_kcal_mol,
            "reference_compound": value.reference_compound,
            "replicate": value.replicate,
        }
        for value in values
    ]


def affinity_values_to_rows(
    values: list[AffinityValue] | tuple[AffinityValue, ...],
) -> list[dict[str, object]]:
    """Serialize measured affinities with normalized molar and ΔG fields."""
    return [
        {
            "compound": value.compound,
            "metric": value.metric,
            "value": value.value,
            "unit": value.unit,
            "concentration_M": value.concentration_M,
            "temperature_K": value.temperature_K,
            "replicate": value.replicate,
            "delta_g_standard_kcal_mol": (
                value.delta_g_standard_kcal_mol
            ),
        }
        for value in values
    ]


def build_parser() -> argparse.ArgumentParser:
    """Create the command-line parser."""
    parser = argparse.ArgumentParser(
        description="Analyze compound contacts across seed-* CIF ensembles."
    )
    parser.add_argument("root", type=Path)
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--ligand")
    parser.add_argument("--protein-chain", action="append", default=[])
    parser.add_argument("--ligand-chain")
    parser.add_argument("--minimize", choices=("openmm",))
    parser.add_argument(
        "--analysis-source", choices=("original", "minimized"), default="original"
    )
    parser.add_argument("--ligand-smiles")
    parser.add_argument("--ligand-file", type=Path)
    parser.add_argument(
        "--openmm-platform",
        choices=("CPU", "CUDA", "OpenCL", "Reference"),
        default="CPU",
    )
    parser.add_argument("--openmm-device-index")
    parser.add_argument("--minimize-ph", type=float, default=7.4)
    parser.add_argument("--pocket-radius", type=float, default=6.0)
    parser.add_argument("--max-backbone-rmsd", type=float, default=0.5)
    parser.add_argument("--max-ligand-rmsd", type=float, default=1.5)
    parser.add_argument("--minimization-tolerance", type=float, default=10.0)
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--clash-cutoff", type=float, default=2.0)
    parser.add_argument("--energy-results", type=Path)
    parser.add_argument("--affinity-results", type=Path)
    parser.add_argument("--convert-ki", action="store_true")
    parser.add_argument("--label-frequency", type=float, default=0.5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument(
        "--out-dir", type=Path, default=Path("compound_contact_results")
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate cross-option constraints before reading structures."""
    if not args.root.is_dir():
        raise PipelineError(f"ROOT is not a directory: {args.root}")
    if (
        not math.isfinite(args.cutoff)
        or not math.isfinite(args.clash_cutoff)
        or args.cutoff <= 0
        or args.clash_cutoff <= 0
    ):
        raise PipelineError("distance cutoffs must be finite and positive")
    if not 0.0 <= args.label_frequency <= 1.0:
        raise PipelineError("--label-frequency must be between 0 and 1")
    if args.workers < 1:
        raise PipelineError("--workers must be at least 1")
    if args.fail_fast and args.workers != 1:
        raise PipelineError("--fail-fast requires --workers 1")
    minimization_numbers = (
        args.minimize_ph,
        args.pocket_radius,
        args.max_backbone_rmsd,
        args.max_ligand_rmsd,
        args.minimization_tolerance,
    )
    if any(not math.isfinite(value) or value <= 0 for value in minimization_numbers):
        raise PipelineError("minimization numeric values must be finite and positive")
    if args.analysis_source == "minimized" and not args.minimize:
        raise PipelineError("--analysis-source minimized requires --minimize")
    if args.minimize and args.workers != 1:
        raise PipelineError("--minimize requires --workers 1")
    if args.openmm_device_index and args.openmm_platform not in {"CUDA", "OpenCL"}:
        raise PipelineError("--openmm-device-index requires a GPU platform")
    if args.ligand_smiles and args.ligand_file:
        raise PipelineError("provide exactly one of --ligand-smiles and --ligand-file")
    if args.manifest and (
        args.ligand
        or args.protein_chain
        or args.ligand_chain
        or args.ligand_smiles
        or args.ligand_file
    ):
        raise PipelineError(
            "manifest chain and ligand settings are authoritative"
        )


def validate_minimization_configuration(
    args: argparse.Namespace, jobs: list[StructureJob]
) -> None:
    """Validate the ligand chemistry required for a requested minimization."""
    if not args.minimize:
        return
    if any(bool(job.ligand_smiles) == bool(job.ligand_file) for job in jobs):
        raise PipelineError(
            "minimization requires exactly one ligand chemistry source per job"
        )
    if args.manifest is None and len({job.compound for job in jobs}) != 1:
        raise PipelineError(
            "global ligand chemistry can only be applied to one discovered compound"
        )


def _protein_chain_tuple(values: list[str]) -> tuple[str, ...] | None:
    chains: list[str] = []
    for value in values:
        for chain in parse_chain_list(value) or ():
            if chain not in chains:
                chains.append(chain)
    return tuple(chains) or None


def _analyze_jobs(
    jobs: list[StructureJob], args: argparse.Namespace
) -> list[StructureResult]:
    if args.workers == 1:
        results: list[StructureResult] = []
        for job in jobs:
            result = analyze_structure(job, args.cutoff, args.clash_cutoff)
            if args.fail_fast and result.status != "valid":
                message = result.qc[-1].message if result.qc else "invalid structure"
                raise PipelineError(
                    f"structure failed under --fail-fast: "
                    f"{job.structure_id}: {message}"
                )
            results.append(result)
        return results

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        return list(
            executor.map(
                analyze_structure,
                jobs,
                repeat(args.cutoff),
                repeat(args.clash_cutoff),
            )
        )


def _minimization_config(args: argparse.Namespace):
    from restrained_openmm import MinimizationConfig

    return MinimizationConfig(
        ph=args.minimize_ph,
        pocket_radius_A=args.pocket_radius,
        tolerance_kj_mol_nm=args.minimization_tolerance,
        platform=args.openmm_platform,
        device_index=args.openmm_device_index,
        max_backbone_rmsd_A=args.max_backbone_rmsd,
        max_ligand_rmsd_A=args.max_ligand_rmsd,
    )


def _contact_residue_set(result: StructureResult) -> set[ResidueKey]:
    return {record.residue for record in result.contacts if record.is_contact}


def run_minimization_jobs(
    jobs: list[StructureJob],
    original_results: list[StructureResult],
    args: argparse.Namespace,
    temporary_root: Path,
    minimize_one=None,
) -> list[MinimizationRecord]:
    """Run requested minimizations sequentially and retain every outcome."""
    from restrained_openmm import (
        MinimizationRequest,
        candidate_relative_path,
        evaluate_acceptance,
        minimize_structure,
    )

    if len(jobs) != len(original_results):
        raise PipelineError("minimization jobs and original results must align")
    minimize_one = minimize_one or minimize_structure
    config = _minimization_config(args)
    temporary_root = Path(temporary_root)
    cache_path = temporary_root / "system_generator_cache.json"
    records: list[MinimizationRecord] = []

    for job, original in zip(jobs, original_results):
        if job.structure_id != original.job.structure_id:
            raise PipelineError("minimization jobs and original results must align")
        if original.status != "valid":
            records.append(
                MinimizationRecord(
                    job, "not_attempted", ("original structure is invalid",),
                    None, None, None, None, None, None, None, None, None,
                    None, None, None, config, (),
                )
            )
            continue

        staged_path = temporary_root / candidate_relative_path(
            job.structure_id, accepted=True
        )
        staged_path.parent.mkdir(parents=True, exist_ok=True)
        ligand_resname = original.metrics.ligand_resname or job.ligand_resname
        if ligand_resname is None:
            raise PipelineError(
                f"valid structure has no selected ligand: {job.structure_id}"
            )
        request = MinimizationRequest(
            cif_path=job.cif_path,
            output_path=staged_path,
            structure_id=job.structure_id,
            compound=job.compound,
            ligand_resname=ligand_resname,
            protein_chains=job.protein_chains,
            ligand_chain=job.ligand_chain,
            ligand_smiles=job.ligand_smiles,
            ligand_file=job.ligand_file,
        )
        try:
            backend = minimize_one(request, config, cache_path)
            candidate_job = StructureJob(
                cif_path=backend.output_path,
                structure_id=job.structure_id,
                dataset_rel=job.dataset_rel,
                compound=job.compound,
                ligand_resname=ligand_resname,
                protein_chains=job.protein_chains,
                ligand_chain=job.ligand_chain,
                ligand_smiles=job.ligand_smiles,
                ligand_file=job.ligand_file,
            )
            candidate = analyze_structure(
                candidate_job, args.cutoff, args.clash_cutoff
            )
            if candidate.status != "valid":
                message = (
                    candidate.qc[-1].message
                    if candidate.qc
                    else "minimized candidate is invalid"
                )
                raise PipelineError(message)
            decision = evaluate_acceptance(
                backbone_rmsd_A=backend.backbone_rmsd_A,
                ligand_rmsd_A=backend.ligand_rmsd_A,
                before_fixed_clashes=original.metrics.clash_pair_count,
                after_fixed_clashes=candidate.metrics.clash_pair_count,
                before_vdw_clashes=backend.before_vdw_clash_pair_count,
                after_vdw_clashes=backend.after_vdw_clash_pair_count,
                max_backbone_rmsd_A=config.max_backbone_rmsd_A,
                max_ligand_rmsd_A=config.max_ligand_rmsd_A,
            )
        except Exception as exc:
            if args.fail_fast:
                raise PipelineError(
                    f"minimization failed under --fail-fast: "
                    f"{job.structure_id}: {exc}"
                ) from exc
            records.append(
                MinimizationRecord(
                    job, "failed", (str(exc),), None, None, None,
                    None, None, None, None, None, None, None, None, None,
                    config, (),
                )
            )
            continue

        status = "accepted" if decision.accepted else "rejected"
        records.append(
            MinimizationRecord(
                job=job,
                status=status,
                reasons=decision.reasons,
                staged_cif_path=backend.output_path,
                output_relative_path=candidate_relative_path(
                    job.structure_id, accepted=decision.accepted
                ),
                candidate_result=candidate,
                initial_potential_kcal_mol=backend.initial_potential_kcal_mol,
                final_potential_kcal_mol=backend.final_potential_kcal_mol,
                backbone_rmsd_A=backend.backbone_rmsd_A,
                ligand_rmsd_A=backend.ligand_rmsd_A,
                before_vdw_clash_pair_count=backend.before_vdw_clash_pair_count,
                after_vdw_clash_pair_count=backend.after_vdw_clash_pair_count,
                before_maximum_vdw_overlap_A=backend.before_maximum_vdw_overlap_A,
                after_maximum_vdw_overlap_A=backend.after_maximum_vdw_overlap_A,
                contact_set_jaccard=contact_set_jaccard(
                    _contact_residue_set(original),
                    _contact_residue_set(candidate),
                ),
                config=config,
                package_versions=backend.package_versions,
            )
        )
    return records


def select_analysis_results(
    jobs: list[StructureJob],
    original_results: list[StructureResult],
    minimization_records: list[MinimizationRecord],
    args: argparse.Namespace,
) -> list[StructureResult]:
    """Select original or strictly accepted minimized coordinates for analysis."""
    if args.analysis_source == "original":
        return original_results
    if len(jobs) != len(original_results) or len(jobs) != len(minimization_records):
        raise PipelineError("minimization records must align with structure jobs")

    selected: list[StructureResult] = []
    for job, original, record in zip(jobs, original_results, minimization_records):
        if record.status == "accepted" and record.candidate_result is not None:
            selected.append(record.candidate_result)
            continue
        stage = (
            "minimization_acceptance"
            if record.status == "rejected"
            else "minimization"
        )
        reason = "; ".join(record.reasons) or f"minimization {record.status}"
        selected.append(
            _invalid_structure_result(job, stage, reason, original.qc)
        )
    return selected


def _validate_auto_detected_ligands(results: list[StructureResult]) -> None:
    detected: dict[str, set[str]] = defaultdict(set)
    for result in results:
        if (
            result.status == "valid"
            and result.job.ligand_resname is None
            and result.metrics.ligand_resname is not None
        ):
            detected[result.job.compound].add(result.metrics.ligand_resname)
    for compound, ligand_names in sorted(detected.items()):
        if len(ligand_names) > 1:
            raise PipelineError(
                f"compound {compound!r} has different auto-detected ligands: "
                + ", ".join(sorted(ligand_names))
            )


def _core_table_specs(
    jobs: list[StructureJob],
    results: list[StructureResult],
    compound_residue_rows: list[dict[str, object]],
    global_residue_rows: list[dict[str, object]],
    compound_rows: list[dict[str, object]],
    minimization_records: list[MinimizationRecord] | None = None,
    original_results: list[StructureResult] | None = None,
    analysis_source: str = "original",
) -> list[tuple[str, list[dict[str, object]], list[str]]]:
    tables = [
        (
            "run_manifest.tsv",
            jobs_to_rows(jobs, minimization_records, analysis_source),
            [
                "structure_id",
                "dataset_dir",
                "compound",
                "ligand_resname",
                "protein_chains",
                "ligand_chain",
                "cif_path",
                "analysis_source",
                "original_cif_path",
                "analysis_cif_path",
            ],
        ),
        (
            "per_structure_contacts.tsv",
            contacts_to_rows(results),
            [
                "structure_id",
                "compound",
                "chain",
                "resseq",
                "icode",
                "resname",
                "is_contact",
                "min_contact_distance_angstrom",
                "capped_distance_angstrom",
            ],
        ),
        (
            "structure_summary.tsv",
            structure_results_to_rows(results, analysis_source),
            [
                "structure_id",
                "compound",
                "status",
                "analysis_source",
                "ligand_resname",
                "n_protein_residues",
                "n_contact_residues",
                "minimum_pair_distance_angstrom",
                "clash_pair_count",
                "buried_sasa_total_A2",
                "interface_area_A2",
            ],
        ),
        (
            "compound_residue_summary.tsv",
            compound_residue_rows,
            [
                "compound",
                "chain",
                "resseq",
                "icode",
                "resname",
                "label",
                "n_valid_structures",
                "n_present",
                "n_contacts",
                "contact_frequency",
                "mean_capped_distance",
                "std_capped_distance",
                "median_capped_distance",
                "min_capped_distance",
                "max_capped_distance",
            ],
        ),
        (
            "global_residue_summary.tsv",
            global_residue_rows,
            [
                "chain",
                "resseq",
                "icode",
                "resname",
                "label",
                "n_valid_structures",
                "n_present",
                "n_contacts",
                "contact_frequency",
                "mean_capped_distance",
                "std_capped_distance",
                "median_capped_distance",
                "min_capped_distance",
                "max_capped_distance",
            ],
        ),
        (
            "compound_summary.tsv",
            compound_rows,
            [
                "compound",
                "n_discovered",
                "n_valid",
                "n_invalid",
                "n_warnings",
                "mean_contact_residues",
                "std_contact_residues",
                "mean_minimum_pair_distance",
                "std_minimum_pair_distance",
                "mean_clash_pair_count",
                "std_clash_pair_count",
                "mean_buried_sasa_total_A2",
                "std_buried_sasa_total_A2",
                "mean_interface_area_A2",
                "std_interface_area_A2",
                "mean_pairwise_contact_jaccard",
            ],
        ),
    ]
    if minimization_records is not None:
        if original_results is None:
            raise PipelineError(
                "original results are required for minimization reporting"
            )
        tables.append(
            (
                "minimization_summary.tsv",
                minimization_records_to_rows(
                    minimization_records, original_results
                ),
                [
                    "structure_id",
                    "compound",
                    "status",
                    "reason",
                    "protocol",
                    "protein_forcefield",
                    "ligand_forcefield",
                    "platform",
                    "device_index",
                    "ph",
                    "pocket_radius_A",
                    "tolerance_kj_mol_nm",
                    "output_cif",
                    "initial_potential_kcal_mol",
                    "final_potential_kcal_mol",
                    "delta_potential_kcal_mol",
                    "backbone_rmsd_A",
                    "ligand_rmsd_A",
                    "before_minimum_pair_distance_A",
                    "after_minimum_pair_distance_A",
                    "before_fixed_clash_pair_count",
                    "after_fixed_clash_pair_count",
                    "before_vdw_clash_pair_count",
                    "after_vdw_clash_pair_count",
                    "before_maximum_vdw_overlap_A",
                    "after_maximum_vdw_overlap_A",
                    "before_contact_residue_count",
                    "after_contact_residue_count",
                    "contact_set_jaccard",
                    "package_versions",
                ],
            )
        )
    return tables


def _energy_table_specs(
    values: list[EnergyValue], summary_rows: list[dict[str, object]]
) -> list[tuple[str, list[dict[str, object]], list[str]]]:
    if not values:
        return []
    return [
        (
            "energy_values_normalized.tsv",
            energy_values_to_rows(values),
            [
                "compound",
                "method",
                "value_kcal_mol",
                "structure_id",
                "uncertainty_kcal_mol",
                "reference_compound",
                "replicate",
            ],
        ),
        (
            "energy_summary.tsv",
            summary_rows,
            [
                "compound",
                "method",
                "reference_compound",
                "n_values",
                "mean",
                "std",
                "median",
                "min",
                "max",
                "ci95_low",
                "ci95_high",
                "mean_reported_uncertainty",
            ],
        ),
    ]


def _affinity_table_specs(
    values: list[AffinityValue], summary_rows: list[dict[str, object]]
) -> list[tuple[str, list[dict[str, object]], list[str]]]:
    if not values:
        return []
    return [
        (
            "affinity_values_normalized.tsv",
            affinity_values_to_rows(values),
            [
                "compound",
                "metric",
                "value",
                "unit",
                "concentration_M",
                "temperature_K",
                "replicate",
                "delta_g_standard_kcal_mol",
            ],
        ),
        (
            "affinity_summary.tsv",
            summary_rows,
            [
                "compound",
                "metric",
                "temperature_K",
                "n_values",
                "mean_concentration_M",
                "std_concentration_M",
                "median_concentration_M",
                "min_concentration_M",
                "max_concentration_M",
                "mean_delta_g_standard_kcal_mol",
                "std_delta_g_standard_kcal_mol",
                "delta_g_ci95_low",
                "delta_g_ci95_high",
            ],
        ),
    ]


def _run_pipeline(
    args: argparse.Namespace,
    temporary_root: Path,
    minimize_one=None,
) -> int:
    """Execute one validated analysis run inside a run-scoped staging area."""
    validate_args(args)
    jobs = build_structure_jobs(
        root=args.root,
        manifest_path=args.manifest,
        default_ligand=(args.ligand.strip() if args.ligand else None),
        default_protein_chains=_protein_chain_tuple(args.protein_chain),
        default_ligand_chain=(
            args.ligand_chain.strip() if args.ligand_chain else None
        ),
        default_ligand_smiles=(
            args.ligand_smiles.strip() if args.ligand_smiles else None
        ),
        default_ligand_file=(
            resolve_inside(args.root, str(args.ligand_file))
            if args.ligand_file
            else None
        ),
    )
    validate_minimization_configuration(args, jobs)
    if not jobs:
        raise PipelineError("no CIF files found under seed-* directories")

    print(f"Analyzing {len(jobs)} CIF structure(s)...", file=sys.stderr)
    original_results = _analyze_jobs(jobs, args)
    _validate_auto_detected_ligands(original_results)
    valid_original_results = [
        result for result in original_results if result.status == "valid"
    ]
    if not valid_original_results:
        raise PipelineError("no valid structures")

    minimization_records = None
    if args.minimize:
        minimization_records = run_minimization_jobs(
            jobs,
            original_results,
            args,
            temporary_root,
            minimize_one=minimize_one,
        )
        results = select_analysis_results(
            jobs, original_results, minimization_records, args
        )
    else:
        results = original_results
    valid_results = [result for result in results if result.status == "valid"]
    if args.analysis_source == "minimized" and not valid_results:
        raise PipelineError("no valid minimized structures")

    known_compounds = {job.compound for job in jobs}
    known_structure_ids = {job.structure_id for job in jobs}
    energy_values = (
        load_energy_results(
            args.energy_results, known_compounds, known_structure_ids
        )
        if args.energy_results
        else []
    )
    affinity_values = (
        load_affinity_results(
            args.affinity_results,
            known_compounds,
            convert_ki=args.convert_ki,
        )
        if args.affinity_results
        else []
    )

    compound_residue_rows = aggregate_residues(results)
    global_residue_rows = aggregate_global_residues(results)
    compound_rows = aggregate_compounds(results)
    energy_summary_rows = aggregate_energy_values(energy_values)
    affinity_summary_rows = aggregate_affinity_values(affinity_values)

    table_specs = _core_table_specs(
        jobs,
        results,
        compound_residue_rows,
        global_residue_rows,
        compound_rows,
        minimization_records=minimization_records,
        original_results=original_results,
        analysis_source=args.analysis_source,
    )
    table_specs.extend(_energy_table_specs(energy_values, energy_summary_rows))
    table_specs.extend(
        _affinity_table_specs(affinity_values, affinity_summary_rows)
    )

    output_dir = Path(args.out_dir)
    prepare_output_directory(output_dir, args.overwrite)
    claimed_tables = {
        name: claim_output_path(output_dir / name, output_dir)
        for name, _, _ in table_specs
    }
    qc_path = claim_output_path(output_dir / "qc.tsv", output_dir)

    generated: list[Path] = []
    if minimization_records is not None:
        for record in minimization_records:
            if (
                record.staged_cif_path is None
                or record.output_relative_path is None
            ):
                continue
            target = claim_output_path(
                output_dir / record.output_relative_path, output_dir
            )
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(record.staged_cif_path, target)
            generated.append(target)
    for name, rows, fieldnames in table_specs:
        target = claimed_tables[name]
        write_tsv(target, rows, fieldnames)
        generated.append(target)

    qc_records = [record for result in results for record in result.qc]
    if minimization_records is not None and args.analysis_source == "original":
        for record in minimization_records:
            if record.status == "accepted":
                continue
            qc_records.append(
                QCRecord(
                    severity="error",
                    stage=(
                        "minimization_acceptance"
                        if record.status == "rejected"
                        else "minimization"
                    ),
                    message="; ".join(record.reasons)
                    or f"minimization {record.status}",
                    compound=record.job.compound,
                    structure_id=record.job.structure_id,
                )
            )
    plots_dir = output_dir / "plots"
    plot_operations = [
        (
            "compound contact panels",
            lambda: plot_compound_chain_contacts(
                compound_residue_rows,
                plots_dir,
                args.cutoff,
                args.label_frequency,
                output_root=output_dir,
            ),
        ),
        (
            "chain heatmaps",
            lambda: plot_chain_heatmaps(
                compound_residue_rows,
                plots_dir,
                args.cutoff,
                output_root=output_dir,
            ),
        ),
        (
            "structural summary",
            lambda: [
                plot_structural_summary(
                    compound_rows, plots_dir, output_root=output_dir
                )
            ],
        ),
    ]
    if energy_values:
        plot_operations.append(
            (
                "computational energy",
                lambda: plot_energy_methods(
                    energy_values,
                    energy_summary_rows,
                    plots_dir,
                    output_root=output_dir,
                ),
            )
        )
    if affinity_values:
        plot_operations.append(
            (
                "experimental affinity",
                lambda: [
                    plot_experimental_affinity(
                        affinity_values,
                        affinity_summary_rows,
                        plots_dir,
                        output_root=output_dir,
                    )
                ],
            )
        )
    if minimization_records is not None:
        plot_operations.append(
            (
                "minimization before/after",
                lambda: [
                    plot_minimization_before_after(
                        minimization_records,
                        original_results,
                        plots_dir,
                        output_root=output_dir,
                    )
                ],
            )
        )

    plot_failed = False
    for label, operation in plot_operations:
        try:
            generated.extend(operation())
        except Exception as exc:
            plot_failed = True
            message = f"{label} plot failed: {exc}"
            print(f"WARNING: {message}", file=sys.stderr)
            qc_records.append(
                QCRecord(
                    severity="error",
                    stage="plotting",
                    message=message,
                )
            )

    write_tsv(
        qc_path,
        qc_to_rows(qc_records),
        ["severity", "stage", "compound", "structure_id", "message"],
    )
    generated.append(qc_path)

    if plot_failed:
        return 1

    inventory = write_generated_inventory(output_dir, generated)
    generated.append(inventory)
    for path in sorted(generated, key=lambda value: value.as_posix()):
        print(path)
    minimization_incomplete = (
        minimization_records is not None
        and any(record.status != "accepted" for record in minimization_records)
    )
    return 1 if minimization_incomplete else 0


def run_pipeline(args: argparse.Namespace, minimize_one=None) -> int:
    """Execute one analysis run and return its process exit code."""
    with tempfile.TemporaryDirectory(prefix="compound-contact-") as temp_dir:
        return _run_pipeline(args, Path(temp_dir), minimize_one=minimize_one)


def main(argv: list[str] | None = None) -> int:
    """Command-line entry point."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return run_pipeline(args)
    except PipelineError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
