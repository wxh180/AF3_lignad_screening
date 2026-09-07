# Compound Contact and Energy Aggregation Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a portable Python CLI that recursively analyzes protein–ligand contacts across compound/seed CIF ensembles, aggregates structural metrics, imports method-labeled energy or affinity data, and creates reproducible tables and figures.

**Architecture:** Keep deployment to one executable `compound_contact_pipeline.py`, but organize it internally as typed data models followed by discovery, structure analysis, aggregation, auxiliary-data, output, plotting, and orchestration functions. Tests exercise each boundary through public functions and finish with a real CLI run on synthetic mmCIF files. Energy calculations remain external; this tool validates, aggregates, and visualizes their results.

**Tech Stack:** Python 3.10+, standard library, Biopython, NumPy, Matplotlib, pytest.

**Spec:** `docs/superpowers/specs/2026-09-01-compound-contact-energy-pipeline-design.md`

## Global Constraints

- Support Python 3.10 and newer.
- Runtime dependencies are limited to Biopython, NumPy, and Matplotlib.
- Analyze only the first model in each CIF and warn when additional models exist.
- Contact means an atom-pair distance less than or equal to the configured cutoff.
- Missing residues do not receive cutoff-filled values and do not enter residue denominators.
- Protein-chain filtering never filters ligand discovery.
- Never label a contact, SASA, clash, Vina, or MM/GBSA metric as measured binding affinity.
- Never aggregate values across different energy methods or RBFE reference compounds.
- Do not invoke shell commands, docking software, molecular-dynamics software, or network services from the pipeline.
- A nonempty output directory requires `--overwrite`; overwrite cleanup is restricted to paths in the preceding generated-file inventory.
- All table and figure orderings must be deterministic.

## File Structure

- Create `compound_contact_pipeline.py`: complete executable and importable pipeline.
- Create `tests/test_compound_contact_pipeline.py`: unit, integration, plot-smoke, and CLI tests.
- Create `datasets.example.tsv`: valid manifest example.
- Create `energy_results.example.tsv`: valid method-labeled energy example.
- Create `affinity_results.example.tsv`: valid experimental-affinity example.
- Create `README.md`: installation, inputs, commands, outputs, interpretation, and physical-energy workflow.

---

### Task 1: Dataset Models, Recursive Discovery, and Manifest Resolution

**Files:**
- Create: `compound_contact_pipeline.py`
- Create: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Produces: `PipelineError`, `DatasetSpec`, `StructureJob`, `QCRecord`, `find_seed_dirs()`, `find_cif_files()`, `load_manifest()`, and `build_structure_jobs()`.
- Consumes: only standard-library paths, CSV parsing, and dataclasses.

- [ ] **Step 1: Write failing discovery and manifest tests**

Add tests that create nested empty `.cif` files and assert deterministic discovery, outermost-seed deduplication, duplicate-basename preservation, path containment, overlap rejection, and manifest-free compound assignment:

```python
from pathlib import Path

import pytest

import compound_contact_pipeline as ccp


def touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def test_find_cif_files_recurses_and_deduplicates_nested_seed_dirs(tmp_path):
    first = touch(tmp_path / "C19" / "seed-1" / "model.cif")
    nested = touch(tmp_path / "C19" / "seed-1" / "seed-extra" / "nested.cif")
    second = touch(tmp_path / "GSH" / "seed-2" / "model.cif")

    found = ccp.find_cif_files(tmp_path)

    assert found == [first.resolve(), nested.resolve(), second.resolve()]


def test_build_structure_jobs_preserves_duplicate_basenames(tmp_path):
    touch(tmp_path / "C19" / "seed-1" / "model.cif")
    touch(tmp_path / "C19" / "seed-2" / "model.cif")

    jobs = ccp.build_structure_jobs(
        root=tmp_path,
        manifest_path=None,
        default_ligand="WXI",
        default_protein_chains=("A",),
        default_ligand_chain="B",
    )

    assert [job.structure_id for job in jobs] == [
        "C19/seed-1/model.cif",
        "C19/seed-2/model.cif",
    ]
    assert {job.compound for job in jobs} == {"C19"}


def test_manifest_rejects_dataset_path_outside_root(tmp_path):
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\n../outside\tC19\tWXI\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="outside ROOT"):
        ccp.load_manifest(tmp_path, manifest)


def test_manifest_rejects_overlapping_enabled_datasets(tmp_path):
    touch(tmp_path / "screen" / "C19" / "seed-1" / "model.cif")
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\n"
        "screen\tALL\tWXI\n"
        "screen/C19\tC19\tWXI\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="more than one dataset"):
        ccp.build_structure_jobs(tmp_path, manifest, None, None, None)
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'find_cif or structure_jobs or manifest'
```

Expected: collection fails because `compound_contact_pipeline` does not exist.

- [ ] **Step 3: Implement models, recursive discovery, and manifest validation**

Create these exact data models and function signatures:

```python
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path


class PipelineError(RuntimeError):
    pass


@dataclass(frozen=True)
class DatasetSpec:
    dataset_dir: Path
    dataset_rel: str
    compound: str
    ligand_resname: str | None
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None


@dataclass(frozen=True)
class StructureJob:
    cif_path: Path
    structure_id: str
    dataset_rel: str
    compound: str
    ligand_resname: str | None
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None


@dataclass(frozen=True)
class QCRecord:
    severity: str
    stage: str
    message: str
    compound: str = ""
    structure_id: str = ""


def parse_chain_list(value: str) -> tuple[str, ...] | None:
    chains = tuple(part.strip() for part in value.split(",") if part.strip())
    return chains or None


def resolve_inside(root: Path, relative_value: str) -> Path:
    root = root.resolve()
    candidate = (root / relative_value).resolve()
    try:
        candidate.relative_to(root)
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
        path for path in candidates
        if not any(parent in candidate_set for parent in path.parents)
    ]


def find_cif_files(base: Path) -> list[Path]:
    files = {
        cif.resolve()
        for seed_dir in find_seed_dirs(base)
        for cif in seed_dir.rglob("*.cif")
        if cif.is_file()
    }
    root = base.resolve()
    return sorted(files, key=lambda path: path.relative_to(root).as_posix())
```

Implement `load_manifest(root, manifest_path)` by opening the file as UTF-8
with `newline=""` and constructing `csv.DictReader(handle, delimiter="\t")`.
Apply exact required-column validation, restrict `enabled` to `true`/`false`,
check that directories exist, reject conflicting compound settings, and include
row numbers in `PipelineError` messages. Implement
`build_structure_jobs(root, manifest_path, default_ligand,
default_protein_chains, default_ligand_chain)` so manifest rows are
authoritative, each CIF matches exactly one row, and manifest-free jobs derive
`compound` from each outermost seed directory's parent.

- [ ] **Step 4: Run focused tests and verify GREEN**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'find_cif or structure_jobs or manifest'
```

Expected: all focused tests pass.

- [ ] **Step 5: Commit Task 1**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: discover compound seed datasets"
```

---

### Task 2: CIF Parsing, Ligand Detection, Contacts, and Geometric Clashes

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: `StructureJob`, contact cutoff, and clash cutoff.
- Produces: `ResidueKey`, `ContactRecord`, `StructureMetrics`, `StructureResult`, `detect_ligand_resname()`, and `analyze_structure()`.

- [ ] **Step 1: Add a synthetic mmCIF writer and failing contact tests**

Build synthetic structures with Biopython rather than embedding proprietary data:

```python
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure
import numpy as np


def write_test_cif(path: Path, residue_specs, model_count: int = 1) -> Path:
    structure = Structure.Structure("test")
    serial = 1
    for model_index in range(model_count):
        model = Model.Model(model_index)
        structure.add(model)
        chains = {}
        for chain_id, hetflag, resseq, icode, resname, atoms in residue_specs:
            chain = chains.get(chain_id)
            if chain is None:
                chain = Chain.Chain(chain_id)
                chains[chain_id] = chain
                model.add(chain)
            residue = Residue.Residue((hetflag, resseq, icode), resname, "")
            chain.add(residue)
            for atom_name, element, coord in atoms:
                atom = Atom.Atom(
                    atom_name,
                    np.asarray(coord, dtype=float),
                    1.0,
                    1.0,
                    " ",
                    f" {atom_name:<3}"[:4],
                    serial,
                    element=element,
                )
                serial += 1
                residue.add(atom)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = MMCIFIO()
    writer.set_structure(structure)
    writer.save(str(path))
    return path


def test_analyze_structure_separates_protein_and_ligand_chain_filters(tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("C", " ", 20, " ", "GLY", [("CA", "C", (20.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (4.5, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(cif, "C19/seed-1/model.cif", "C19", "C19", "WXI", ("A",), "B")

    result = ccp.analyze_structure(job, cutoff=4.5, clash_cutoff=2.0)

    assert result.status == "valid"
    assert len(result.contacts) == 1
    assert result.contacts[0].is_contact is True
    assert result.contacts[0].min_contact_distance_angstrom == pytest.approx(4.5)


def test_noncontact_is_capped_but_raw_contact_distance_is_blank(tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "far.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (10.0, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(cif, "C19/seed-1/far.cif", "C19", "C19", "WXI", ("A",), "B")

    contact = ccp.analyze_structure(job, 4.5, 2.0).contacts[0]

    assert contact.is_contact is False
    assert contact.min_contact_distance_angstrom is None
    assert contact.capped_distance_angstrom == 4.5


def test_auto_detection_refuses_ambiguous_nonprotein_residues(tmp_path):
    cif = write_test_cif(
        tmp_path / "mix" / "seed-1" / "model.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_L01", 1, " ", "L01", [("C1", "C", (3.0, 0.0, 0.0))]),
            ("B", "H_L02", 2, " ", "L02", [("C1", "C", (4.0, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(cif, "mix/seed-1/model.cif", "mix", "mix", None, ("A",), None)

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "invalid"
    assert "multiple ligand candidates" in result.qc[0].message


def test_insertion_code_clashes_and_extra_model_warning_are_preserved(tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "models.cif",
        [
            ("A", " ", 10, "B", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (2.0, 0.0, 0.0))]),
        ],
        model_count=2,
    )
    job = ccp.StructureJob(cif, "C19/seed-1/models.cif", "C19", "C19", "WXI", ("A",), "B")

    result = ccp.analyze_structure(job, cutoff=4.5, clash_cutoff=2.0)

    assert result.contacts[0].residue.icode == "B"
    assert result.metrics.clash_pair_count == 1
    assert any(record.stage == "model_selection" for record in result.qc)


@pytest.mark.parametrize(
    ("residue_specs", "message"),
    [
        (
            [("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))])],
            "ligand WXI not found",
        ),
        (
            [("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (0.0, 0.0, 0.0))])],
            "no protein atoms",
        ),
    ],
)
def test_missing_required_structure_components_are_invalid(tmp_path, residue_specs, message):
    cif = write_test_cif(tmp_path / "C19" / "seed-1" / "missing.cif", residue_specs)
    job = ccp.StructureJob(cif, "C19/seed-1/missing.cif", "C19", "C19", "WXI", ("A",), "B")

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "invalid"
    assert any(message in record.message for record in result.qc)


def test_malformed_cif_is_an_invalid_result(tmp_path):
    cif = touch(tmp_path / "C19" / "seed-1" / "bad.cif")
    cif.write_text("not an mmCIF", encoding="utf-8")
    job = ccp.StructureJob(cif, "C19/seed-1/bad.cif", "C19", "C19", "WXI", None, None)

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "invalid"
    assert result.qc[0].stage == "parse"
```

- [ ] **Step 2: Run contact tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'analyze_structure or noncontact or auto_detection'
```

Expected: tests fail because the analysis types and functions are absent.

- [ ] **Step 3: Implement structure parsing and geometry**

Add exact data models:

```python
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
```

Implement ligand detection with this explicit ion set and rule:

```python
COMMON_MONATOMIC_IONS = {
    "AL", "BA", "BR", "CA", "CD", "CL", "CO", "CS", "CU", "F",
    "FE", "GA", "HG", "IOD", "K", "LI", "MG", "MN", "NA", "NI",
    "RB", "SR", "ZN",
}


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
```

In `analyze_structure(job, cutoff, clash_cutoff)`, parse with `MMCIFParser(QUIET=True)`, select `models[0]`, identify ligand residues independently of the protein-chain filter, and create one record for every present protein residue. Use `NeighborSearch(protein_atoms)` for cutoff contacts and NumPy coordinate differences for the exact minimum distance and inclusive clash-pair count:

```python
minimum_pair_distance = float("inf")
clash_pair_count = 0
protein_coordinates = np.asarray([atom.get_coord() for atom in protein_atoms])
for ligand_atom in ligand_atoms:
    distances = np.linalg.norm(protein_coordinates - ligand_atom.get_coord(), axis=1)
    minimum_pair_distance = min(minimum_pair_distance, float(distances.min()))
    clash_pair_count += int(np.count_nonzero(distances <= clash_cutoff))
```

Catch parsing and selection failures and return an invalid `StructureResult` with a stage-specific `QCRecord`. Add a warning when `len(models) > 1`. Sort contact records by `ResidueKey`.
An invalid result uses `StructureMetrics(None, 0, 0, None, 0)` so downstream
code never has to special-case a missing metrics object.

- [ ] **Step 4: Run contact tests and verify GREEN**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'analyze_structure or noncontact or auto_detection'
```

Expected: focused tests pass.

- [ ] **Step 5: Commit Task 2**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: analyze ligand contacts in cif structures"
```

---

### Task 3: Buried SASA and Structure-Level Metrics

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: selected protein and ligand residues inside `analyze_structure()`.
- Produces: `build_sasa_entity()` and `calculate_buried_sasa()`; populates the SASA fields in `StructureMetrics` without invalidating contacts on SASA failure.

- [ ] **Step 1: Write failing SASA tests**

```python
def test_buried_sasa_is_positive_for_touching_receptor_and_ligand(tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "touching.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(cif, "C19/seed-1/touching.cif", "C19", "C19", "WXI", ("A",), "B")

    metrics = ccp.analyze_structure(job, 4.5, 2.0).metrics

    assert metrics.buried_sasa_total_A2 is not None
    assert metrics.buried_sasa_total_A2 > 0.0
    assert metrics.interface_area_A2 == pytest.approx(metrics.buried_sasa_total_A2 / 2.0)


def test_sasa_failure_becomes_warning_without_invalidating_contacts(monkeypatch, tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(cif, "C19/seed-1/model.cif", "C19", "C19", "WXI", ("A",), "B")
    monkeypatch.setattr(ccp, "calculate_buried_sasa", lambda protein, ligand: (_ for _ in ()).throw(ValueError("bad SASA")))

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "valid"
    assert result.metrics.buried_sasa_total_A2 is None
    assert any(record.stage == "sasa" for record in result.qc)
```

- [ ] **Step 2: Run SASA tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k sasa
```

Expected: tests fail because SASA is not calculated.

- [ ] **Step 3: Implement isolated and complex SASA calculations**

Use copied residues so SASA calculations cannot mutate the parsed structure used for contacts:

```python
from copy import deepcopy

from Bio.PDB import Chain, Model, ShrakeRupley, Structure


def build_sasa_entity(residues, entity_id: str):
    entity = Structure.Structure(entity_id)
    model = Model.Model(0)
    entity.add(model)
    chains = {}
    for residue in residues:
        chain_id = residue.get_parent().id
        if chain_id not in chains:
            chains[chain_id] = Chain.Chain(chain_id)
            model.add(chains[chain_id])
        chains[chain_id].add(deepcopy(residue))
    return entity


def entity_sasa(entity) -> float:
    calculator = ShrakeRupley(probe_radius=1.4)
    calculator.compute(entity, level="S")
    return float(entity.sasa)


def calculate_buried_sasa(protein_residues, ligand_residues) -> tuple[float, float]:
    receptor_sasa = entity_sasa(build_sasa_entity(protein_residues, "receptor"))
    ligand_sasa = entity_sasa(build_sasa_entity(ligand_residues, "ligand"))
    complex_sasa = entity_sasa(build_sasa_entity([*protein_residues, *ligand_residues], "complex"))
    buried_total = max(0.0, receptor_sasa + ligand_sasa - complex_sasa)
    return buried_total, buried_total / 2.0
```

Call this function after contact geometry succeeds. Catch only the SASA exception at that boundary, preserve valid contact output, set both SASA fields to `None`, and append a warning QC record.

- [ ] **Step 4: Run SASA and prior tests**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'sasa or analyze_structure or noncontact'
```

Expected: all selected tests pass.

- [ ] **Step 5: Commit Task 3**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: calculate buried ligand interface sasa"
```

---

### Task 4: Residue, Compound, and Global Aggregation

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: ordered `StructureResult` values.
- Produces: `summarize_values()`, `contact_set_jaccard()`, `aggregate_residues()`, `aggregate_compounds()`, and `aggregate_global_residues()`, each returning ordered `list[dict[str, object]]` rows suitable for TSV output.

- [ ] **Step 1: Write failing denominator and Jaccard tests**

```python
def make_contact(structure_id, compound, residue, is_contact, capped):
    return ccp.ContactRecord(
        structure_id,
        compound,
        residue,
        is_contact,
        capped if is_contact else None,
        capped,
    )


def make_valid_result(job, contacts, buried_sasa=100.0):
    metrics = ccp.StructureMetrics(
        ligand_resname=job.ligand_resname,
        n_protein_residues=len(contacts),
        n_contact_residues=sum(contact.is_contact for contact in contacts),
        minimum_pair_distance_angstrom=3.0,
        clash_pair_count=0,
        buried_sasa_total_A2=buried_sasa,
        interface_area_A2=buried_sasa / 2.0,
    )
    return ccp.StructureResult(job, "valid", tuple(contacts), metrics, ())


def test_residue_frequency_uses_only_structures_where_residue_is_present(tmp_path):
    residue = ccp.ResidueKey("A", 10, "", "ALA")
    job1 = ccp.StructureJob(tmp_path / "one.cif", "C/seed-1/one.cif", "C", "C", "LIG", None, None)
    job2 = ccp.StructureJob(tmp_path / "two.cif", "C/seed-2/two.cif", "C", "C", "LIG", None, None)
    first = make_valid_result(job1, [make_contact(job1.structure_id, "C", residue, True, 3.0)])
    second = make_valid_result(job2, [])

    row = ccp.aggregate_residues([first, second])[0]

    assert row["n_valid_structures"] == 2
    assert row["n_present"] == 1
    assert row["n_contacts"] == 1
    assert row["contact_frequency"] == 1.0


def test_contact_set_jaccard_has_explicit_empty_set_behavior():
    assert ccp.contact_set_jaccard(set(), set()) == 1.0
    assert ccp.contact_set_jaccard(set(), {"A:1"}) == 0.0
    assert ccp.contact_set_jaccard({"A:1", "A:2"}, {"A:2", "A:3"}) == pytest.approx(1 / 3)


def test_global_residue_summary_uses_present_structures_across_compounds(tmp_path):
    residue = ccp.ResidueKey("A", 10, "", "ALA")
    job1 = ccp.StructureJob(tmp_path / "one.cif", "C1/seed-1/one.cif", "C1", "C1", "LIG", None, None)
    job2 = ccp.StructureJob(tmp_path / "two.cif", "C2/seed-1/two.cif", "C2", "C2", "LIG", None, None)
    first = make_valid_result(job1, [make_contact(job1.structure_id, "C1", residue, True, 3.0)])
    second = make_valid_result(job2, [make_contact(job2.structure_id, "C2", residue, False, 4.5)])

    row = ccp.aggregate_global_residues([first, second])[0]

    assert row["n_present"] == 2
    assert row["n_contacts"] == 1
    assert row["contact_frequency"] == 0.5
```

- [ ] **Step 2: Run aggregation tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'frequency_uses or jaccard'
```

Expected: tests fail because aggregation functions are absent.

- [ ] **Step 3: Implement deterministic statistics and grouping**

Use population standard deviation and blank statistics for empty input:

```python
import statistics
from collections import defaultdict
from itertools import combinations


def summarize_values(values) -> dict[str, float | int | None]:
    numeric = [float(value) for value in values if value is not None]
    if not numeric:
        return {"count": 0, "mean": None, "std": None, "median": None, "min": None, "max": None}
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
```

In `aggregate_residues()`, create a compound-to-valid-structure count, group
records by `(compound, ResidueKey)`, and calculate frequency with `n_present`.
Return each row with these fields in this order:

```python
[
    "compound", "chain", "resseq", "icode", "resname", "label",
    "n_valid_structures", "n_present", "n_contacts", "contact_frequency",
    "mean_capped_distance", "std_capped_distance",
    "median_capped_distance", "min_capped_distance", "max_capped_distance",
]
```

In `aggregate_compounds()`, include invalid structures in discovered/invalid
counts but calculate numeric metrics only from valid structures. Its stable
columns are `compound`, `n_discovered`, `n_valid`, `n_invalid`, `n_warnings`,
`mean_contact_residues`, `std_contact_residues`,
`mean_minimum_pair_distance`, `std_minimum_pair_distance`,
`mean_clash_pair_count`, `std_clash_pair_count`, `mean_buried_sasa_total_A2`,
`std_buried_sasa_total_A2`, `mean_interface_area_A2`,
`std_interface_area_A2`, and `mean_pairwise_contact_jaccard`. Calculate mean
pairwise Jaccard over `combinations(contact_sets, 2)` and leave it blank for one
valid seed. Implement `aggregate_global_residues()` with the same
missing-residue denominator rule but without a compound key. Sort all rows by
compound, chain, residue number, insertion code, and residue name.

- [ ] **Step 4: Run aggregation and full current tests**

```bash
pytest -q tests/test_compound_contact_pipeline.py
```

Expected: all current tests pass.

- [ ] **Step 5: Commit Task 4**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: aggregate contacts across compounds and seeds"
```

---

### Task 5: Computed-Energy and Experimental-Affinity Inputs

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: TSV paths, known compounds, known structure IDs, and `convert_ki`.
- Produces: `EnergyValue`, `AffinityValue`, `load_energy_results()`, `aggregate_energy_values()`, `load_affinity_results()`, `aggregate_affinity_values()`, and `bootstrap_mean_ci()`.

- [ ] **Step 1: Write failing validation, isolation, and conversion tests**

```python
def test_energy_aggregation_never_mixes_methods(tmp_path):
    table = tmp_path / "energy.tsv"
    table.write_text(
        "compound\tmethod\tvalue_kcal_mol\treplicate\n"
        "C19\tvina_score\t-8.0\t1\n"
        "C19\tmmgbsa\t-30.0\t1\n",
        encoding="utf-8",
    )

    values = ccp.load_energy_results(table, {"C19"}, set())
    rows = ccp.aggregate_energy_values(values)

    assert [row["method"] for row in rows] == ["mmgbsa", "vina_score"]
    assert [row["mean"] for row in rows] == [-30.0, -8.0]


def test_energy_table_rejects_unknown_nonblank_structure_id(tmp_path):
    table = tmp_path / "energy.tsv"
    table.write_text(
        "compound\tmethod\tvalue_kcal_mol\tstructure_id\n"
        "C19\tmmgbsa\t-30.0\tmissing.cif\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="unknown structure_id"):
        ccp.load_energy_results(table, {"C19"}, {"C19/seed-1/model.cif"})


def test_kd_converts_to_standard_free_energy_but_ic50_does_not(tmp_path):
    table = tmp_path / "affinity.tsv"
    table.write_text(
        "compound\tmetric\tvalue\tunit\ttemperature_K\treplicate\n"
        "C19\tKd\t1\tnM\t298.15\t1\n"
        "C19\tIC50\t2\tnM\t298.15\t1\n",
        encoding="utf-8",
    )

    values = ccp.load_affinity_results(table, {"C19"}, convert_ki=False)

    assert values[0].delta_g_standard_kcal_mol == pytest.approx(-12.278, abs=0.01)
    assert values[1].delta_g_standard_kcal_mol is None


def test_bootstrap_mean_ci_is_reproducible():
    assert ccp.bootstrap_mean_ci([1.0, 2.0, 3.0]) == ccp.bootstrap_mean_ci([1.0, 2.0, 3.0])


def test_rbfe_series_rejects_mixed_reference_compounds(tmp_path):
    table = tmp_path / "energy.tsv"
    table.write_text(
        "compound\tmethod\tvalue_kcal_mol\treference_compound\treplicate\n"
        "C19\trbfe\t1.0\tREF1\t1\n"
        "C19\trbfe\t1.2\tREF2\t2\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="mixed reference compounds"):
        ccp.load_energy_results(table, {"C19", "REF1", "REF2"}, set())


def test_ki_conversion_requires_explicit_opt_in(tmp_path):
    table = tmp_path / "affinity.tsv"
    table.write_text(
        "compound\tmetric\tvalue\tunit\nC19\tKi\t5\tnM\n",
        encoding="utf-8",
    )

    without_conversion = ccp.load_affinity_results(table, {"C19"}, convert_ki=False)[0]
    with_conversion = ccp.load_affinity_results(table, {"C19"}, convert_ki=True)[0]

    assert without_conversion.delta_g_standard_kcal_mol is None
    assert with_conversion.delta_g_standard_kcal_mol is not None
```

- [ ] **Step 2: Run auxiliary-data tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'energy or standard_free or bootstrap'
```

Expected: tests fail because auxiliary data types and loaders are absent.

- [ ] **Step 3: Implement strict TSV normalization and method-specific grouping**

Create exact models and constants:

```python
import math

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


def bootstrap_mean_ci(values: list[float]) -> tuple[float | None, float | None]:
    if len(values) < 2:
        return None, None
    generator = np.random.default_rng(20260901)
    array = np.asarray(values, dtype=float)
    indices = generator.integers(0, len(array), size=(10_000, len(array)))
    means = array[indices].mean(axis=1)
    low, high = np.percentile(means, [2.5, 97.5])
    return float(low), float(high)
```

Parse numeric fields with `math.isfinite`, reject nonpositive affinity/temperature, reject negative uncertainty, and define the energy duplicate key as `(compound, method, structure_id, reference_compound, replicate)`. A blank energy `structure_id` is compound-level. Validate that each `rbfe` method series has one nonblank reference compound. Group energy rows only by `(compound, method, reference_compound)`. Normalize affinities to molar units, calculate `R*T*log(concentration_M)` for `Kd` and opt-in `Ki`, never for `IC50`, and summarize by `(compound, metric, temperature_K)`.

- [ ] **Step 4: Run auxiliary and full tests**

```bash
pytest -q tests/test_compound_contact_pipeline.py
```

Expected: all current tests pass.

- [ ] **Step 5: Commit Task 5**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: integrate energy and affinity evidence"
```

---

### Task 6: Safe Output Directory, TSV Serialization, and Generated-File Inventory

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: output path, jobs, structure results, QC records, and ordered aggregate row dictionaries.
- Produces: `jobs_to_rows()`, `contacts_to_rows()`, `structure_results_to_rows()`, `qc_to_rows()`, `format_tsv_value()`, `write_tsv()`, `prepare_output_directory()`, `claim_output_path()`, and `write_generated_inventory()`.

- [ ] **Step 1: Write failing formatting and containment tests**

```python
import json


def test_write_tsv_formats_floats_and_blanks_missing_values(tmp_path):
    output = tmp_path / "values.tsv"

    ccp.write_tsv(output, [{"name": "A", "value": 1.23456789, "missing": None}])

    assert output.read_text(encoding="utf-8") == "name\tvalue\tmissing\nA\t1.234568\t\n"


def test_nonempty_output_directory_requires_overwrite(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    (output / "existing.txt").write_text("keep", encoding="utf-8")

    with pytest.raises(ccp.PipelineError, match="--overwrite"):
        ccp.prepare_output_directory(output, overwrite=False)


def test_overwrite_rejects_inventory_path_escape_before_deleting_files(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    kept = output / "kept.tsv"
    kept.write_text("keep", encoding="utf-8")
    (output / "generated_files.json").write_text(
        json.dumps({"files": ["kept.tsv", "../victim.txt", "generated_files.json"]}),
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="outside output directory"):
        ccp.prepare_output_directory(output, overwrite=True)

    assert kept.exists()


def test_unlisted_existing_output_cannot_be_replaced(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    unrelated = output / "compound_summary.tsv"
    unrelated.write_text("unrelated", encoding="utf-8")

    ccp.prepare_output_directory(output, overwrite=True)

    with pytest.raises(ccp.PipelineError, match="not listed in the preceding inventory"):
        ccp.claim_output_path(unrelated, output)
    assert unrelated.read_text(encoding="utf-8") == "unrelated"
```

- [ ] **Step 2: Run output tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'write_tsv or output_directory or inventory'
```

Expected: tests fail because output helpers are absent.

- [ ] **Step 3: Implement narrow overwrite behavior and deterministic TSVs**

```python
import json


def format_tsv_value(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.6f}"
    if isinstance(value, bool):
        return "1" if value else "0"
    return str(value)


def write_tsv(
    path: Path,
    rows: list[dict[str, object]],
    fieldnames: list[str] | None = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None and not rows:
        raise PipelineError(f"cannot infer TSV columns for empty table: {path.name}")
    columns = list(rows[0]) if fieldnames is None else list(fieldnames)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({name: format_tsv_value(row.get(name)) for name in columns})


def _validated_inventory_paths(output_dir: Path, relative_paths: list[str]) -> list[Path]:
    resolved_output = output_dir.resolve()
    validated = []
    for relative_path in relative_paths:
        candidate = (resolved_output / relative_path).resolve()
        try:
            candidate.relative_to(resolved_output)
        except ValueError as exc:
            raise PipelineError(f"inventory path is outside output directory: {relative_path}") from exc
        validated.append(candidate)
    return validated


def claim_output_path(path: Path, output_dir: Path) -> Path:
    resolved_output = output_dir.resolve()
    if path.is_symlink():
        raise PipelineError(f"output path was not listed in the preceding inventory: {path}")
    candidate = path.resolve()
    try:
        candidate.relative_to(resolved_output)
    except ValueError as exc:
        raise PipelineError(f"output path is outside output directory: {path}") from exc
    if candidate.exists():
        raise PipelineError(f"output path was not listed in the preceding inventory: {path}")
    return candidate
```

In `prepare_output_directory(output_dir, overwrite)`, reject a nonempty
directory without overwrite. With overwrite, read `generated_files.json` when
present, validate every listed path before unlinking any path, then unlink only
existing regular files or symlinks in that validated list. Do not recursively
delete directories. Before every table or figure write, call
`claim_output_path()`; this prevents a new run from replacing an unrelated file
that was not listed in the prior inventory. `write_generated_inventory()`
writes sorted POSIX-relative paths, including `generated_files.json` itself,
only after all requested tables and figures succeed.

Add row adapters with fixed schemas. `jobs_to_rows()` emits `structure_id`,
`dataset_dir`, `compound`, `ligand_resname`, `protein_chains`, `ligand_chain`,
and `cif_path`. `contacts_to_rows()` emits structure/compound plus residue-key,
contact, and capped-distance fields. `structure_results_to_rows()` emits status,
resolved ligand, protein-residue count, contact-residue count, minimum distance,
clash count, and both SASA metrics. `qc_to_rows()` emits `severity`, `stage`,
`compound`, `structure_id`, and `message`. Pass these exact fields to
`write_tsv()` when a table such as `qc.tsv` is empty so it still receives a
header.

- [ ] **Step 4: Run output and full tests**

```bash
pytest -q tests/test_compound_contact_pipeline.py
```

Expected: all current tests pass.

- [ ] **Step 5: Commit Task 6**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: write deterministic protected outputs"
```

---

### Task 7: Contact Heatmaps, Per-Chain Panels, and Evidence-Specific Summary Plots

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: aggregate row dictionaries and an output plots directory.
- Produces: `safe_filename_map()`, `build_heatmap_matrix()`, `plot_compound_chain_contacts()`, `plot_chain_heatmaps()`, `plot_structural_summary()`, `plot_energy_methods()`, and `plot_experimental_affinity()`; plotting functions return generated `Path` objects.

- [ ] **Step 1: Write failing filename, mask, and plot-smoke tests**

```python
def test_safe_filename_map_resolves_sanitization_collisions_deterministically():
    mapping = ccp.safe_filename_map(["A/B", "A B"])

    assert mapping["A/B"] != mapping["A B"]
    assert mapping == ccp.safe_filename_map(["A/B", "A B"])


def test_contact_plot_functions_create_nonempty_pngs(tmp_path):
    rows = [
        {
            "compound": "C19",
            "chain": "A",
            "resseq": 10,
            "icode": "",
            "resname": "ALA",
            "label": "ALA A:10",
            "n_valid_structures": 2,
            "n_present": 2,
            "n_contacts": 1,
            "contact_frequency": 0.5,
            "mean_capped_distance": 3.75,
            "std_capped_distance": 0.75,
        }
    ]

    paths = ccp.plot_compound_chain_contacts(rows, tmp_path, cutoff=4.5, label_frequency=0.5)
    paths.extend(ccp.plot_chain_heatmaps(rows, tmp_path, cutoff=4.5))

    assert paths
    assert all(path.exists() and path.stat().st_size > 0 for path in paths)


def test_heatmap_matrix_uses_nan_for_absent_compound_residue_pairs():
    rows = [
        {"compound": "C1", "chain": "A", "resseq": 10, "icode": "", "resname": "ALA", "contact_frequency": 1.0},
        {"compound": "C2", "chain": "A", "resseq": 11, "icode": "", "resname": "GLY", "contact_frequency": 0.5},
    ]

    matrix, compounds, residues = ccp.build_heatmap_matrix(rows, "contact_frequency")

    assert compounds == ["C1", "C2"]
    assert [residue.resseq for residue in residues] == [10, 11]
    assert np.isnan(matrix[0, 1])
    assert np.isnan(matrix[1, 0])
```

- [ ] **Step 2: Run plot tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'safe_filename or plot_functions'
```

Expected: tests fail because plotting functions are absent.

- [ ] **Step 3: Implement headless, masked, method-separated plotting**

Set the backend before importing `pyplot`:

```python
import hashlib
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def safe_filename_map(values) -> dict[str, str]:
    result = {}
    used = set()
    for value in sorted(set(values)):
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


def build_heatmap_matrix(rows, value_key: str):
    compounds = sorted({str(row["compound"]) for row in rows})
    residues = sorted({
        ResidueKey(str(row["chain"]), int(row["resseq"]), str(row["icode"]), str(row["resname"]))
        for row in rows
    })
    compound_index = {value: index for index, value in enumerate(compounds)}
    residue_index = {value: index for index, value in enumerate(residues)}
    matrix = np.full((len(compounds), len(residues)), np.nan, dtype=float)
    for row in rows:
        residue = ResidueKey(str(row["chain"]), int(row["resseq"]), str(row["icode"]), str(row["resname"]))
        matrix[compound_index[str(row["compound"])], residue_index[residue]] = float(row[value_key])
    return matrix, compounds, residues
```

Implement two-panel compound/chain figures using residue number on the x-axis, blue mean capped distance ± population SD above, orange contact frequency below, a dashed cutoff line, light grids, and staggered labels at or above `label_frequency`. For heatmaps, create the union of residue keys per chain, fill absent compound–residue combinations with `np.nan`, plot `np.ma.masked_invalid(matrix)`, and use a copied colormap whose bad color is `lightgray`. Make one frequency and one capped-distance heatmap per chain. Each plotting function receives the run's output root separately from its plots directory and calls `claim_output_path()` immediately before `savefig()`.

Implement structural summary as four vertical panels for valid-seed count, Jaccard consistency, buried SASA, and clash count. Implement one energy figure per `(method, reference_compound)` with individual replicate points and mean/SD/CI, preserving method-specific y-axis labels. Plot experimental concentration on a log scale and standard ΔG in a separate panel. Close every figure in `finally` after saving at 200 dpi.

- [ ] **Step 4: Run plot and full tests**

```bash
pytest -q tests/test_compound_contact_pipeline.py
```

Expected: all tests pass without GUI requirements or Matplotlib warnings.

- [ ] **Step 5: Commit Task 7**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: visualize compound contact evidence"
```

---

### Task 8: CLI Orchestration, Parallel Execution, Examples, and Documentation

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`
- Create: `datasets.example.tsv`
- Create: `energy_results.example.tsv`
- Create: `affinity_results.example.tsv`
- Create: `README.md`

**Interfaces:**
- Consumes: every interface from Tasks 1–7.
- Produces: `build_parser()`, `validate_args()`, `run_pipeline()`, `main()`, executable CLI behavior, example inputs, and user documentation.

- [ ] **Step 1: Write failing CLI and end-to-end tests**

```python
import subprocess
import sys


def test_cli_rejects_fail_fast_with_multiple_workers(tmp_path):
    parser = ccp.build_parser()
    args = parser.parse_args([str(tmp_path), "--fail-fast", "--workers", "2"])

    with pytest.raises(ccp.PipelineError, match="--workers 1"):
        ccp.validate_args(args)


def test_cli_end_to_end_writes_core_tables_and_contact_plots(tmp_path):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    output = tmp_path / "results"

    completed = subprocess.run(
        [
            sys.executable,
            str(Path(ccp.__file__).resolve()),
            str(root),
            "--ligand",
            "WXI",
            "--protein-chain",
            "A",
            "--ligand-chain",
            "B",
            "--out-dir",
            str(output),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    for name in [
        "run_manifest.tsv",
        "per_structure_contacts.tsv",
        "structure_summary.tsv",
        "compound_residue_summary.tsv",
        "global_residue_summary.tsv",
        "compound_summary.tsv",
        "qc.tsv",
        "generated_files.json",
    ]:
        assert (output / name).exists()
    assert list((output / "plots").glob("*.png"))


def test_cli_rejects_conflicting_auto_detected_ligands_within_compound(tmp_path):
    root = tmp_path / "inputs"
    for seed, ligand in [("seed-1", "L01"), ("seed-2", "L02")]:
        write_test_cif(
            root / "mix" / seed / "model.cif",
            [
                ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
                ("B", f"H_{ligand}", 1, " ", ligand, [("C1", "C", (3.0, 0.0, 0.0))]),
            ],
        )

    completed = subprocess.run(
        [
            sys.executable,
            str(Path(ccp.__file__).resolve()),
            str(root),
            "--out-dir",
            str(tmp_path / "results"),
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "different auto-detected ligands" in completed.stderr


def test_plot_failure_preserves_tables_and_returns_one(monkeypatch, tmp_path):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    output = tmp_path / "results"
    args = ccp.build_parser().parse_args([
        str(root), "--ligand", "WXI", "--out-dir", str(output)
    ])

    def fail_plot(*args, **kwargs):
        raise RuntimeError("plot failed")

    monkeypatch.setattr(ccp, "plot_compound_chain_contacts", fail_plot)

    assert ccp.run_pipeline(args) == 1
    assert (output / "compound_summary.tsv").exists()
    assert "plot failed" in (output / "qc.tsv").read_text(encoding="utf-8")
    assert not (output / "generated_files.json").exists()


def test_cli_continues_past_malformed_cif_when_one_structure_is_valid(tmp_path):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "good.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    bad = touch(root / "C19" / "seed-2" / "bad.cif")
    bad.write_text("broken", encoding="utf-8")
    output = tmp_path / "results"

    completed = subprocess.run(
        [sys.executable, str(Path(ccp.__file__).resolve()), str(root), "--ligand", "WXI", "--out-dir", str(output)],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0
    qc_text = (output / "qc.tsv").read_text(encoding="utf-8")
    assert "C19/seed-2/bad.cif" in qc_text
    assert "parse" in qc_text


def test_cli_returns_nonzero_when_no_structure_is_valid(tmp_path):
    root = tmp_path / "inputs"
    bad = touch(root / "C19" / "seed-1" / "bad.cif")
    bad.write_text("broken", encoding="utf-8")

    completed = subprocess.run(
        [sys.executable, str(Path(ccp.__file__).resolve()), str(root), "--ligand", "WXI"],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "no valid structures" in completed.stderr
```

- [ ] **Step 2: Run CLI tests and verify RED**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k cli
```

Expected: tests fail because the CLI orchestration is absent.

- [ ] **Step 3: Implement parser, validation, orchestration, and exit semantics**

Create exact options:

```python
import argparse
from concurrent.futures import ProcessPoolExecutor


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Analyze compound contacts across seed-* CIF ensembles.")
    parser.add_argument("root", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--ligand")
    parser.add_argument("--protein-chain", action="append", default=[])
    parser.add_argument("--ligand-chain")
    parser.add_argument("--cutoff", type=float, default=4.5)
    parser.add_argument("--clash-cutoff", type=float, default=2.0)
    parser.add_argument("--energy-results", type=Path)
    parser.add_argument("--affinity-results", type=Path)
    parser.add_argument("--convert-ki", action="store_true")
    parser.add_argument("--label-frequency", type=float, default=0.5)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--fail-fast", action="store_true")
    parser.add_argument("--out-dir", type=Path, default=Path("compound_contact_results"))
    parser.add_argument("--overwrite", action="store_true")
    return parser


def validate_args(args) -> None:
    if not args.root.is_dir():
        raise PipelineError(f"ROOT is not a directory: {args.root}")
    if args.cutoff <= 0 or args.clash_cutoff <= 0:
        raise PipelineError("distance cutoffs must be positive")
    if not 0.0 <= args.label_frequency <= 1.0:
        raise PipelineError("--label-frequency must be between 0 and 1")
    if args.workers < 1:
        raise PipelineError("--workers must be at least 1")
    if args.fail_fast and args.workers != 1:
        raise PipelineError("--fail-fast requires --workers 1")
    if args.manifest and (args.ligand or args.protein_chain or args.ligand_chain):
        raise PipelineError("manifest chain and ligand settings are authoritative")
```

`run_pipeline(args)` performs these operations in order:

1. Validate arguments and build ordered jobs.
2. Analyze sequentially when `workers == 1`; otherwise use `ProcessPoolExecutor(max_workers=args.workers).map()` so result order matches job order.
3. Under fail-fast, stop immediately when a sequential result is invalid.
4. Validate one auto-detected ligand residue name per compound; conflicting names are fatal.
5. Require at least one valid structure.
6. Load and validate optional energy and affinity tables using the discovered compound and structure-ID sets.
7. Build all aggregate rows in memory.
8. Prepare the output directory only after configuration and input validation succeeds.
9. Claim each target path, then write core tables, optional normalized/summary tables, and figures while collecting generated paths.
10. Write `qc.tsv`, then write the generated-file inventory only if every requested plot succeeds.
11. Return `1` when any plot fails after tables are written; otherwise return `0`.

Wrap `main()` so `PipelineError` writes `ERROR: <message>` to standard error and returns `2`; unexpected exceptions are not swallowed. Print generated paths to standard output only after success.

- [ ] **Step 4: Create exact example TSV files**

`datasets.example.tsv`:

```tsv
dataset_dir	compound	ligand_resname	protein_chains	ligand_chain	enabled
screen/C19	C19	WXI	A	B	true
screen/GSH	GSH	GSH	A,C	L	true
```

`energy_results.example.tsv`:

```tsv
compound	method	value_kcal_mol	structure_id	uncertainty_kcal_mol	reference_compound	replicate
C19	mmgbsa	-31.2	screen/C19/seed-1/model.cif	1.8		1
C19	mmgbsa	-29.9	screen/C19/seed-2/model.cif	1.6		2
GSH	rbfe	1.4		0.4	C19	1
```

`affinity_results.example.tsv`:

```tsv
compound	metric	value	unit	temperature_K	replicate
C19	Kd	25	nM	298.15	1
C19	Kd	31	nM	298.15	2
GSH	IC50	2.1	uM	298.15	1
```

- [ ] **Step 5: Write README usage and interpretation guidance**

Document:

- Python dependencies and installation command.
- Manifest and manifest-free layouts.
- The two primary commands from the specification.
- Every output table and plot.
- Contact cutoff, capped-distance, missing-residue denominator, clash, SASA, and Jaccard definitions.
- Why structural metrics and Vina are not binding affinity.
- How to import ensemble MM/GBSA, OpenFE RBFE, or experimental Kd results.
- The formula `ΔG° = RT ln(Kd / 1 M)` and the opt-in Ki/never-automatic IC50 rules.
- An AF3-specific warning that predicted pose consistency is useful prioritization evidence but does not replace energetic sampling or experiment.

Include this quick-start command verbatim:

```bash
python compound_contact_pipeline.py /path/to/project \
  --manifest datasets.tsv \
  --cutoff 4.5 \
  --out-dir compound_contact_results
```

- [ ] **Step 6: Run CLI tests and the complete suite**

```bash
pytest -q tests/test_compound_contact_pipeline.py -k cli
pytest -q
python -m py_compile compound_contact_pipeline.py
python compound_contact_pipeline.py --help
```

Expected: all tests pass, compilation exits zero, and help lists every specified option.

- [ ] **Step 7: Commit Task 8**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py README.md datasets.example.tsv energy_results.example.tsv affinity_results.example.tsv
git commit -m "feat: complete compound contact analysis pipeline"
```

---

## Final Verification

- [ ] Run the complete verification set from a clean working tree:

```bash
pytest -q
python -m py_compile compound_contact_pipeline.py
python compound_contact_pipeline.py --help >/dev/null
git diff --check
git status --short
```

- [ ] Confirm every acceptance criterion in the specification maps to evidence:

| Acceptance criterion | Evidence |
| --- | --- |
| Arbitrary nested discovery | Task 1 recursive/nested discovery test |
| Unique IDs for repeated filenames | Task 1 duplicate-basename test |
| Correct missing-residue and invalid-structure denominators | Task 4 denominator test plus Task 8 malformed-CIF test |
| Core tables and figures | Task 8 end-to-end test |
| Energy methods remain separate | Task 5 method-isolation test |
| Kd converts and IC50 does not | Task 5 affinity test |
| Structural metrics are not affinity | README text review and method-specific plot labels |
| No docking or MD dependency | End-to-end test runs using only declared Python dependencies |

- [ ] Inspect generated PNG files from the end-to-end temporary run for readable labels, gray missing-data cells, uncropped axes, and method-specific units.

- [ ] Confirm `git status --short` is empty before completion is reported.
