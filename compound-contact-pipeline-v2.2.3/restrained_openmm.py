"""NumPy geometry/QC and a lazily loaded restrained OpenMM backend.

The optional OpenMM stack is intentionally loaded only on demand. Keeping this
module importable with NumPy alone lets the contact-only pipeline remain light.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time
from typing import Callable, Sequence

import numpy as np


class MinimizationDependencyError(RuntimeError):
    """The optional physical-minimization environment is unavailable."""


class MinimizationBackendError(RuntimeError):
    """A selected structure could not be safely minimized."""


def load_openmm_stack(*, import_module=None):
    """Load optional packages only when physical minimization is requested."""
    if import_module is None:
        from importlib import import_module
    from types import SimpleNamespace

    try:
        return SimpleNamespace(
            openmm=import_module("openmm"),
            app=import_module("openmm.app"),
            unit=import_module("openmm.unit"),
            generators=import_module("openmmforcefields.generators"),
            toolkit=import_module("openff.toolkit"),
        )
    except ImportError as exc:
        raise MinimizationDependencyError(
            "OpenMM minimization requires the optional packages; install and activate "
            "environment-openmm.yml (OpenMM, OpenFF Toolkit, openmmforcefields)."
        ) from exc


@dataclass(frozen=True)
class MinimizationRequest:
    cif_path: Path
    output_path: Path
    structure_id: str
    compound: str
    ligand_resname: str
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None
    ligand_smiles: str | None
    ligand_file: Path | None


@dataclass(frozen=True)
class BackendMinimizationResult:
    output_path: Path
    initial_potential_kcal_mol: float
    final_potential_kcal_mol: float
    backbone_rmsd_A: float
    ligand_rmsd_A: float
    before_vdw_clash_pair_count: int
    after_vdw_clash_pair_count: int
    before_maximum_vdw_overlap_A: float
    after_maximum_vdw_overlap_A: float
    package_versions: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class MinimizationProgressEvent:
    operation: str
    state: str
    elapsed_seconds: float | None = None


ProgressCallback = Callable[[MinimizationProgressEvent], None]


def _emit_progress(
    callback: ProgressCallback | None,
    event: MinimizationProgressEvent,
) -> None:
    if callback is None:
        return
    try:
        callback(event)
    except Exception:
        return


@dataclass(frozen=True)
class MinimizationStage:
    """One bounded restrained-minimization stage."""

    name: str
    max_iterations: int
    backbone_k_kcal_A2: float
    distant_k_kcal_A2: float
    pocket_sidechain_k_kcal_A2: float
    ligand_k_kcal_A2: float


DEFAULT_STAGES = (
    MinimizationStage("hydrogen_relaxation", 500, 10.0, 10.0, 10.0, 10.0),
    MinimizationStage("pocket_relaxation", 1000, 10.0, 10.0, 0.0, 1.0),
    MinimizationStage("gentle_relaxation", 1000, 2.0, 2.0, 0.0, 0.2),
)


@dataclass(frozen=True)
class MinimizationConfig:
    """Stable, versioned parameters for restrained OpenMM minimization."""

    protocol: str = "openmm_restrained_v1"
    protein_forcefield: str = "amber/protein.ff14SB.xml"
    ligand_forcefield: str = "gaff-2.2.20"
    ph: float = 7.4
    pocket_radius_A: float = 6.0
    tolerance_kj_mol_nm: float = 10.0
    platform: str = "CPU"
    device_index: str | None = None
    max_backbone_rmsd_A: float = 0.5
    max_ligand_rmsd_A: float = 1.5
    stages: tuple[MinimizationStage, ...] = DEFAULT_STAGES


def kcal_per_A2_to_kj_per_nm2(value: float) -> float:
    """Convert a force constant from kcal mol^-1 A^-2 to kJ mol^-1 nm^-2."""
    return value * 4.184 * 100.0


@dataclass(frozen=True)
class OverlapMetrics:
    """Element-aware protein-ligand van der Waals overlap measurements."""

    clash_pair_count: int
    maximum_overlap_A: float
    minimum_pair_distance_A: float | None


@dataclass(frozen=True)
class AcceptanceDecision:
    """Immutable result of the conservative minimized-pose acceptance gates."""

    accepted: bool
    reasons: tuple[str, ...]


VDW_RADII_A = {
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
}


def _coordinates(values: object, name: str) -> np.ndarray:
    coordinates = np.asarray(values, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1:] != (3,):
        raise ValueError(f"{name} must have shape (n, 3)")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"{name} must be finite")
    return coordinates


def _element_radius(element: str) -> float:
    normalized = str(element).strip().capitalize()
    try:
        return VDW_RADII_A[normalized]
    except KeyError as exc:
        raise ValueError(f"unsupported element for van der Waals QC: {element!r}") from exc


def vdw_overlap_metrics(
    protein_coordinates_A: object,
    protein_elements: Sequence[str],
    ligand_coordinates_A: object,
    ligand_elements: Sequence[str],
    threshold_A: float = 0.4,
) -> OverlapMetrics:
    """Measure element-aware protein-ligand clashes without changing inputs."""
    protein = _coordinates(protein_coordinates_A, "protein_coordinates_A")
    ligand = _coordinates(ligand_coordinates_A, "ligand_coordinates_A")
    if len(protein_elements) != len(protein):
        raise ValueError("protein_elements must match protein_coordinates_A")
    if len(ligand_elements) != len(ligand):
        raise ValueError("ligand_elements must match ligand_coordinates_A")
    if not np.isfinite(threshold_A):
        raise ValueError("threshold_A must be finite")

    protein_radii = np.asarray([_element_radius(element) for element in protein_elements])
    ligand_radii = np.asarray([_element_radius(element) for element in ligand_elements])
    if len(protein) == 0 or len(ligand) == 0:
        return OverlapMetrics(0, 0.0, None)

    distances = np.linalg.norm(protein[:, None, :] - ligand[None, :, :], axis=2)
    overlaps = protein_radii[:, None] + ligand_radii[None, :] - distances
    inclusive_threshold = threshold_A - np.finfo(float).eps * max(1.0, abs(threshold_A))
    return OverlapMetrics(
        clash_pair_count=int(np.count_nonzero(overlaps >= inclusive_threshold)),
        maximum_overlap_A=float(np.max(overlaps)),
        minimum_pair_distance_A=float(np.min(distances)),
    )


def _rmsd(reference: np.ndarray, mobile: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.sum((reference - mobile) ** 2, axis=1))))


def align_on_backbone(
    reference_backbone_A: object,
    mobile_backbone_A: object,
    reference_ligand_A: object,
    mobile_ligand_A: object,
) -> tuple[float, float]:
    """Return backbone and protein-aligned ligand RMSD using row-vector Kabsch."""
    reference_backbone = _coordinates(reference_backbone_A, "reference_backbone_A")
    mobile_backbone = _coordinates(mobile_backbone_A, "mobile_backbone_A")
    reference_ligand = _coordinates(reference_ligand_A, "reference_ligand_A")
    mobile_ligand = _coordinates(mobile_ligand_A, "mobile_ligand_A")
    if len(reference_backbone) != len(mobile_backbone):
        raise ValueError("backbone coordinate sets must have the same length")
    if len(reference_ligand) != len(mobile_ligand):
        raise ValueError("ligand coordinate sets must have the same length")
    if len(reference_backbone) == 0:
        raise ValueError("backbone coordinate sets must not be empty")
    if len(reference_ligand) == 0:
        raise ValueError("ligand coordinate sets must not be empty")

    reference_centroid = reference_backbone.mean(axis=0)
    mobile_centroid = mobile_backbone.mean(axis=0)
    covariance = (mobile_backbone - mobile_centroid).T @ (
        reference_backbone - reference_centroid
    )
    left, _, right_transpose = np.linalg.svd(covariance)
    rotation = left @ right_transpose
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_transpose
    translation = reference_centroid - mobile_centroid @ rotation
    aligned_backbone = mobile_backbone @ rotation + translation
    aligned_ligand = mobile_ligand @ rotation + translation
    return _rmsd(reference_backbone, aligned_backbone), _rmsd(
        reference_ligand, aligned_ligand
    )


def select_pocket_residues(
    protein_coordinates_A: object,
    protein_residue_ids: Sequence[str],
    ligand_coordinates_A: object,
    radius_A: float,
) -> frozenset[str]:
    """Select protein residues with any atom at or inside the ligand radius."""
    protein = _coordinates(protein_coordinates_A, "protein_coordinates_A")
    ligand = _coordinates(ligand_coordinates_A, "ligand_coordinates_A")
    if len(protein_residue_ids) != len(protein):
        raise ValueError("protein_residue_ids must match protein_coordinates_A")
    if not np.isfinite(radius_A) or radius_A <= 0:
        raise ValueError("radius_A must be finite and positive")
    if len(protein) == 0 or len(ligand) == 0:
        return frozenset()

    distances = np.linalg.norm(protein[:, None, :] - ligand[None, :, :], axis=2)
    close_atoms = np.any(distances <= radius_A, axis=1)
    return frozenset(
        str(residue_id)
        for residue_id, is_close in zip(protein_residue_ids, close_atoms)
        if is_close
    )


def evaluate_acceptance(
    *,
    backbone_rmsd_A: float,
    ligand_rmsd_A: float,
    before_fixed_clashes: int,
    after_fixed_clashes: int,
    before_vdw_clashes: int,
    after_vdw_clashes: int,
    max_backbone_rmsd_A: float,
    max_ligand_rmsd_A: float,
) -> AcceptanceDecision:
    """Apply stable, conservative pose-preservation and clash-count gates."""
    numeric_inputs = {
        "backbone_rmsd_A": backbone_rmsd_A,
        "ligand_rmsd_A": ligand_rmsd_A,
        "before_fixed_clashes": before_fixed_clashes,
        "after_fixed_clashes": after_fixed_clashes,
        "before_vdw_clashes": before_vdw_clashes,
        "after_vdw_clashes": after_vdw_clashes,
        "max_backbone_rmsd_A": max_backbone_rmsd_A,
        "max_ligand_rmsd_A": max_ligand_rmsd_A,
    }
    for input_name, value in numeric_inputs.items():
        if not np.isfinite(value):
            raise ValueError(f"{input_name} must be finite")

    reasons: list[str] = []
    if backbone_rmsd_A > max_backbone_rmsd_A:
        reasons.append(
            f"backbone RMSD {backbone_rmsd_A:.3f} A exceeds "
            f"{max_backbone_rmsd_A:.3f} A"
        )
    if ligand_rmsd_A > max_ligand_rmsd_A:
        reasons.append(
            f"ligand RMSD {ligand_rmsd_A:.3f} A exceeds {max_ligand_rmsd_A:.3f} A"
        )
    if after_fixed_clashes > before_fixed_clashes:
        reasons.append(
            f"fixed clash count increased from {before_fixed_clashes} to {after_fixed_clashes}"
        )
    if after_vdw_clashes > before_vdw_clashes:
        reasons.append(
            "van der Waals clash count increased from "
            f"{before_vdw_clashes} to {after_vdw_clashes}"
        )
    return AcceptanceDecision(accepted=not reasons, reasons=tuple(reasons))


def candidate_relative_path(structure_id: str, accepted: bool) -> Path:
    """Return the deterministic, contained publication path for a candidate CIF."""
    structure_path = Path(structure_id)
    if (
        not structure_id
        or structure_path.is_absolute()
        or any(part in {"", ".", ".."} for part in structure_path.parts)
    ):
        raise ValueError(f"unsafe structure_id: {structure_id!r}")
    if structure_path.suffix != ".cif":
        raise ValueError(f"structure_id must end in .cif: {structure_id!r}")
    bucket = "accepted" if accepted else "rejected"
    output_name = f"{structure_path.stem}_minimized.cif"
    return Path("minimized_structures") / bucket / structure_path.with_name(output_name)


BACKBONE_NAMES = frozenset({"N", "CA", "C", "O"})


def _residue_identity(residue) -> str:
    return repr((residue.chain.id, residue.id, residue.insertionCode, residue.name))


def _atom_identity(atom) -> tuple[str, str, str]:
    return (_residue_identity(atom.residue), atom.name, atom.element.symbol)


def _read_input_atom_records(cif_path):
    """Capture first-model input identities before OpenMM normalizes them."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    data = MMCIF2Dict(str(cif_path))
    columns = {key.removeprefix("_atom_site."): value for key, value in data.items()
               if key.startswith("_atom_site.")}
    records = {}
    models = columns.get("pdbx_PDB_model_num", ["1"] * len(columns["id"]))
    for i, serial in enumerate(columns["id"]):
        if models[i] != models[0]:
            continue
        if serial in records:
            raise MinimizationBackendError(f"duplicate first-model atom_site.id: {serial}")
        records[serial] = {key: values[i] for key, values in columns.items()}
    return records


def _bind_input_atom_records(topology, records):
    """Associate by preserved atom_site.id, not changed names or chain IDs."""
    chain_labels = {}
    for atom in topology.atoms():
        if atom.id not in records:
            raise MinimizationBackendError(f"OpenMM atom has no input identity: {atom.id}")
        row = records[atom.id]
        if atom.element is None or atom.element.symbol.upper() != row["type_symbol"].upper():
            raise MinimizationBackendError(f"input/OpenMM element mismatch for atom {atom.id}")
        chain = atom.residue.chain
        label = row["label_asym_id"]
        if id(chain) in chain_labels and chain_labels[id(chain)] != label:
            raise MinimizationBackendError("OpenMM chain merges distinct input label chains")
        chain_labels[id(chain)] = label
        # Keep distinct label-chain segments separate throughout preparation.
        # User-facing author-chain selection is performed against raw records.
        chain.id = label


def _select_residues(topology, request, input_records=None):
    # Match the contact pipeline's extended amino-acid classification, lazily.
    from Bio.PDB.Polypeptide import is_aa

    residues = list(topology.residues())
    def selection_identity(residue):
        if input_records is None:
            return residue.chain.id, residue.name
        rows = [input_records[a.id] for a in residue.atoms()]
        identities = {(row["auth_asym_id"], row["label_comp_id"]) for row in rows}
        if len(identities) != 1:
            raise MinimizationBackendError("OpenMM residue merges distinct input residue identities")
        return next(iter(identities))

    identities = {id(r): selection_identity(r) for r in residues}
    ligands = [r for r in residues if identities[id(r)][1] == request.ligand_resname and
               (request.ligand_chain is None or identities[id(r)][0] == request.ligand_chain)]
    if len(ligands) != 1:
        raise MinimizationBackendError(
            f"{request.structure_id}: exactly one ligand residue is required; found "
            f"{len(ligands)}. Select a ligand chain or remove extra ligand copies."
        )
    ligand = ligands[0]
    protein = [r for r in residues if r is not ligand and is_aa(identities[id(r)][1], standard=False)
               and (request.protein_chains is None or identities[id(r)][0] in request.protein_chains)]
    if not protein:
        raise MinimizationBackendError("no selected protein residues")
    selected_ids = {id(r) for r in protein} | {id(ligand)}
    return [r for r in residues if id(r) in selected_ids], ligand


def _write_original_cif(topology, positions_nm, input_records, output_path, app, unit):
    """Publish raw label/author namespaces without changing prepared identities.

    OpenMM Topology carries only one chain/name namespace, and its writer emits
    that namespace into both label and author columns. Write the label topology
    (including correct bond references), then restore the author columns from
    the original atom-site records using a CIF codec, never text substitution.
    """
    from copy import deepcopy
    from io import StringIO
    from Bio.PDB import MMCIFIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    publication = deepcopy(topology)
    rows = [input_records[a.id] for a in publication.atoms()]
    for atom, row in zip(publication.atoms(), rows):
        atom.name = row["label_atom_id"]
        atom.residue.name = row["label_comp_id"]
        atom.residue.id = row["label_seq_id"]
        atom.residue.chain.id = row["label_asym_id"]
        insertion = row.get("pdbx_PDB_ins_code", "")
        atom.residue.insertionCode = "" if insertion in {".", "?"} else insertion
    encoded = StringIO()
    app.PDBxFile.writeFile(publication, unit.Quantity(positions_nm, unit.nanometer), encoded, keepIds=True)
    encoded.seek(0)
    data = MMCIF2Dict(encoded)
    for field in ("id", "group_PDB", "label_atom_id", "label_comp_id", "label_asym_id",
                  "label_seq_id", "auth_atom_id", "auth_comp_id", "auth_asym_id",
                  "auth_seq_id", "pdbx_PDB_ins_code"):
        if all(field in row for row in rows):
            data[f"_atom_site.{field}"] = [row[field] for row in rows]
    writer = MMCIFIO()
    writer.set_dict(data)
    writer.save(str(output_path))


def _identity_indices(atoms):
    indices = {}
    for index, atom in enumerate(atoms):
        identity = _atom_identity(atom)
        if identity in indices:
            raise MinimizationBackendError(f"duplicate atom identity: {identity}")
        indices[identity] = index
    return indices


def _original_positions(original_atoms, prepared_atoms, positions_nm):
    """Restore original atom order, never positional-index order after H addition."""
    positions = _coordinates(positions_nm, "minimized positions")
    indices = _identity_indices(prepared_atoms)
    _identity_indices(original_atoms)
    if len(positions) != len(prepared_atoms):
        raise MinimizationBackendError("position count does not match prepared topology")
    try:
        return positions[[indices[_atom_identity(atom)] for atom in original_atoms]]
    except KeyError as exc:
        raise MinimizationBackendError(f"missing original atom: {exc.args[0]}") from exc


def _restraint_constants(atoms, ligand_id, pocket, stage):
    constants = []
    for atom in atoms:
        residue_id = _residue_identity(atom.residue)
        if atom.element.atomic_number == 1:
            value = 0.0
        elif residue_id == ligand_id:
            value = stage.ligand_k_kcal_A2
        elif atom.name in BACKBONE_NAMES:
            value = stage.backbone_k_kcal_A2
        elif residue_id not in pocket:
            value = stage.distant_k_kcal_A2
        else:
            value = stage.pocket_sidechain_k_kcal_A2
        constants.append(kcal_per_A2_to_kj_per_nm2(value))
    return constants


def _ensure_ligand_connectivity(topology, ligand, molecule, positions_nm):
    """Recover AF3 bondless ligand connectivity only after chemistry validation."""
    ligand_atoms = list(ligand.atoms())
    heavy = [atom for atom in ligand_atoms if atom.element.atomic_number != 1]
    ligand_ids = {id(atom) for atom in ligand_atoms}
    heavy_ids = {id(atom) for atom in heavy}
    existing_heavy_bonds = []
    for atom1, atom2 in topology.bonds():
        in1 = id(atom1) in ligand_ids
        in2 = id(atom2) in ligand_ids
        if in1 != in2:
            raise MinimizationBackendError("covalently attached ligands are unsupported")
        if in1 and id(atom1) in heavy_ids and id(atom2) in heavy_ids:
            existing_heavy_bonds.append((atom1, atom2))
    if existing_heavy_bonds:
        return

    molecular_atoms = list(molecule.atoms)
    molecular_heavy = [
        index for index, atom in enumerate(molecular_atoms)
        if atom.atomic_number != 1
    ]
    message = "coordinate-derived connectivity does not match supplied chemistry"
    if len(heavy) != len(molecular_heavy) or not heavy:
        raise MinimizationBackendError(message)

    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError as exc:
        raise MinimizationDependencyError(
            "AF3 bondless-ligand recovery requires RDKit from environment-openmm.yml"
        ) from exc

    all_atoms = list(topology.atoms())
    atom_positions = _coordinates(positions_nm, "input positions")
    if len(all_atoms) != len(atom_positions):
        raise MinimizationBackendError("input position count does not match topology")
    position_index = {id(atom): index for index, atom in enumerate(all_atoms)}

    editable = Chem.RWMol()
    conformer = Chem.Conformer(len(heavy))
    for index, atom in enumerate(heavy):
        editable.AddAtom(Chem.Atom(int(atom.element.atomic_number)))
        x, y, z = atom_positions[position_index[id(atom)]] * 10.0
        conformer.SetAtomPosition(index, (float(x), float(y), float(z)))
    inferred = editable.GetMol()
    inferred.AddConformer(conformer)
    try:
        rdDetermineBonds.DetermineConnectivity(inferred)
    except Exception as exc:
        raise MinimizationBackendError(
            f"could not infer ligand connectivity from coordinates: {exc}"
        ) from exc

    inferred_adjacency = {index: set() for index in range(len(heavy))}
    for bond in inferred.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        inferred_adjacency[i].add(j)
        inferred_adjacency[j].add(i)

    molecular_adjacency = {index: set() for index in molecular_heavy}
    for bond in molecule.bonds:
        i, j = bond.atom1_index, bond.atom2_index
        if i in molecular_adjacency and j in molecular_adjacency:
            molecular_adjacency[i].add(j)
            molecular_adjacency[j].add(i)

    candidates = {
        i: [
            j for j in molecular_heavy
            if heavy[i].element.atomic_number == molecular_atoms[j].atomic_number
            and len(inferred_adjacency[i]) == len(molecular_adjacency[j])
        ]
        for i in inferred_adjacency
    }
    order = sorted(
        inferred_adjacency,
        key=lambda i: (len(candidates[i]), -len(inferred_adjacency[i])),
    )
    mapping = {}

    def match(depth):
        if depth == len(order):
            return True
        i = order[depth]
        for j in candidates[i]:
            if j in mapping.values():
                continue
            if all(
                (other in inferred_adjacency[i])
                == (mapped in molecular_adjacency[j])
                for other, mapped in mapping.items()
            ):
                mapping[i] = j
                if match(depth + 1):
                    return True
                del mapping[i]
        return False

    if not match(0):
        raise MinimizationBackendError(message)

    for bond in inferred.GetBonds():
        topology.addBond(
            heavy[bond.GetBeginAtomIdx()],
            heavy[bond.GetEndAtomIdx()],
        )


def _ligand_hydrogens(topology, ligand, molecule):
    """Match explicit heavy connectivity and supply custom Modeller H variants.

    CIF need not carry bond orders: exact chemistry comes from the supplied
    molecule. No bonds are guessed from atom names or interatomic distances.
    """
    ligand_atoms = list(ligand.atoms())
    heavy = [a for a in ligand_atoms if a.element.atomic_number != 1]
    molecular_atoms = list(molecule.atoms)
    molecular_heavy = [i for i, a in enumerate(molecular_atoms) if a.atomic_number != 1]
    message = ("ligand element/bond graph does not match supplied chemistry "
               "(or missing heavy atoms)")
    if len(heavy) != len(molecular_heavy) or not heavy:
        raise MinimizationBackendError(message)
    heavy_indices = {id(a): i for i, a in enumerate(heavy)}
    ligand_indices = {id(a) for a in ligand_atoms}
    adjacency = {i: set() for i in range(len(heavy))}
    existing_h = {i: [] for i in adjacency}
    for a, b in topology.bonds():
        if (id(a) in ligand_indices) != (id(b) in ligand_indices):
            raise MinimizationBackendError("covalently attached ligands are unsupported")
        if id(a) not in ligand_indices:
            continue
        if id(a) in heavy_indices and id(b) in heavy_indices:
            i, j = heavy_indices[id(a)], heavy_indices[id(b)]
            adjacency[i].add(j)
            adjacency[j].add(i)
        elif a.element.atomic_number == 1 and id(b) in heavy_indices:
            existing_h[heavy_indices[id(b)]].append(a.name)
        elif b.element.atomic_number == 1 and id(a) in heavy_indices:
            existing_h[heavy_indices[id(a)]].append(b.name)
    if sum(map(len, existing_h.values())) != len(ligand_atoms) - len(heavy):
        raise MinimizationBackendError(message)
    molecular_adjacency = {i: set() for i in molecular_heavy}
    hydrogen_counts = {i: 0 for i in molecular_heavy}
    for bond in molecule.bonds:
        i, j = bond.atom1_index, bond.atom2_index
        if i in molecular_adjacency and j in molecular_adjacency:
            molecular_adjacency[i].add(j)
            molecular_adjacency[j].add(i)
        elif i in molecular_adjacency and molecular_atoms[j].atomic_number == 1:
            hydrogen_counts[i] += 1
        elif j in molecular_adjacency and molecular_atoms[i].atomic_number == 1:
            hydrogen_counts[j] += 1
    candidates = {
        i: [j for j in molecular_heavy
            if heavy[i].element.atomic_number == molecular_atoms[j].atomic_number
            and len(adjacency[i]) == len(molecular_adjacency[j])
            and len(existing_h[i]) <= hydrogen_counts[j]]
        for i in adjacency
    }
    order = sorted(adjacency, key=lambda i: (len(candidates[i]), -len(adjacency[i])))
    mapping = {}

    def match(depth):
        if depth == len(order):
            return True
        i = order[depth]
        for j in candidates[i]:
            if j not in mapping.values() and all(
                (other in adjacency[i]) == (mapped in molecular_adjacency[j])
                for other, mapped in mapping.items()
            ):
                mapping[i] = j
                if match(depth + 1):
                    return True
                del mapping[i]
        return False

    if not match(0):
        raise MinimizationBackendError(message)
    names = {a.name for a in ligand_atoms}
    definitions = []
    counter = 1
    for i, atom in enumerate(heavy):
        definitions.extend((name, atom.name) for name in existing_h[i])
        for _ in range(hydrogen_counts[mapping[i]] - len(existing_h[i])):
            while f"H{counter}" in names:
                counter += 1
            name = f"H{counter}"
            names.add(name)
            definitions.append((name, atom.name))
    return definitions


def _load_ligand_chemistry(request, toolkit):
    if (request.ligand_smiles is not None) == (request.ligand_file is not None):
        raise MinimizationBackendError("exactly one ligand chemistry source is required")
    if request.ligand_smiles is not None:
        return toolkit.Molecule.from_smiles(request.ligand_smiles, allow_undefined_stereo=False)
    molecules = toolkit.Molecule.from_file(
        str(request.ligand_file), file_format="SDF", allow_undefined_stereo=False,
    )
    if isinstance(molecules, (list, tuple)):
        if len(molecules) != 1:
            raise MinimizationBackendError("ligand SDF must contain exactly one molecule")
        molecules = molecules[0]
    for conformer in molecules.conformers or ():
        _coordinates(conformer.magnitude, "SDF conformer positions")
    return molecules


def _physical_state(context, unit):
    state = context.getState(getEnergy=True, getPositions=True, groups={0})
    energy = float(state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole))
    if not np.isfinite(energy):
        raise MinimizationBackendError("physical potential energy must be finite")
    positions = _coordinates(
        state.getPositions(asNumpy=True).value_in_unit(unit.nanometer), "minimized positions",
    )
    return energy, positions


def minimize_structure(
    request: MinimizationRequest,
    config: MinimizationConfig,
    cache_path: Path,
    progress_callback: ProgressCallback | None = None,
) -> BackendMinimizationResult:
    """Minimize a selected noncovalent complex without changing original atom IDs."""
    stack = load_openmm_stack()
    source = (f"SMILES {request.ligand_smiles!r}" if request.ligand_smiles is not None
              else f"SDF {request.ligand_file}")
    try:
        return _minimize_structure(
            request,
            config,
            cache_path,
            stack,
            progress_callback=progress_callback,
        )
    except Exception as exc:
        raise MinimizationBackendError(
            f"{request.structure_id} ({source}): {exc}"
        ) from exc


def _minimize_structure(
    request,
    config,
    cache_path,
    stack,
    progress_callback: ProgressCallback | None = None,
):
    preparation_started = time.perf_counter()
    _emit_progress(
        progress_callback,
        MinimizationProgressEvent("preparation", "started"),
    )
    mm, app, unit = stack.openmm, stack.app, stack.unit
    molecule = _load_ligand_chemistry(request, stack.toolkit)
    input_records = _read_input_atom_records(request.cif_path)
    cif = app.PDBxFile(str(request.cif_path))
    _bind_input_atom_records(cif.topology, input_records)
    retained, ligand = _select_residues(cif.topology, request, input_records)
    ligand_id = _residue_identity(ligand)
    input_positions_nm = _coordinates(
        cif.positions.value_in_unit(unit.nanometer), "input positions"
    )
    _ensure_ligand_connectivity(
        cif.topology, ligand, molecule, input_positions_nm
    )
    ligand_hydrogens = _ligand_hydrogens(cif.topology, ligand, molecule)
    modeller = app.Modeller(cif.topology, cif.positions)
    retained_ids = {_residue_identity(r) for r in retained}
    modeller.delete([r for r in modeller.topology.residues()
                     if _residue_identity(r) not in retained_ids])
    # Crystallographic cell metadata must not select SystemGenerator's PME path.
    modeller.topology.setPeriodicBoxVectors(None)
    original_topology = modeller.topology
    original_atoms = list(original_topology.atoms())
    _identity_indices(original_atoms)
    original_nm = _coordinates(modeller.positions.value_in_unit(unit.nanometer), "original positions")
    protein_indices = [i for i, a in enumerate(original_atoms)
                       if a.element.atomic_number != 1 and _residue_identity(a.residue) != ligand_id]
    ligand_indices = [i for i, a in enumerate(original_atoms)
                      if a.element.atomic_number != 1 and _residue_identity(a.residue) == ligand_id]
    backbone_indices = [i for i in protein_indices if original_atoms[i].name in BACKBONE_NAMES]
    protein_elements = [original_atoms[i].element.symbol for i in protein_indices]
    ligand_elements = [original_atoms[i].element.symbol for i in ligand_indices]
    original_A = original_nm * 10.0
    pocket = select_pocket_residues(
        original_A[protein_indices], [_residue_identity(original_atoms[i].residue) for i in protein_indices],
        original_A[ligand_indices], config.pocket_radius_A,
    )
    before = vdw_overlap_metrics(original_A[protein_indices], protein_elements,
                                 original_A[ligand_indices], ligand_elements)
    variants = [(ligand_hydrogens
                 if _residue_identity(r) == ligand_id else None)
                for r in original_topology.residues()]
    generator = stack.generators.SystemGenerator(
        forcefields=[config.protein_forcefield], small_molecule_forcefield=config.ligand_forcefield,
        molecules=[molecule], cache=str(cache_path), forcefield_kwargs={"constraints": app.HBonds},
        nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
    )
    platform = mm.Platform.getPlatformByName(config.platform)
    # A separate Modeller keeps the output topology entirely pre-hydrogenation.
    prepared = app.Modeller(original_topology, modeller.positions)
    prepared.addHydrogens(generator.forcefield, pH=config.ph, variants=variants, platform=platform)
    prepared_atoms = list(prepared.topology.atoms())
    prepared_indices = _identity_indices(prepared_atoms)
    _original_positions(original_atoms, prepared_atoms,
                        prepared.positions.value_in_unit(unit.nanometer))
    system = generator.create_system(prepared.topology)
    force = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    for parameter in ("k", "x0", "y0", "z0"):
        force.addPerParticleParameter(parameter)
    force.setForceGroup(31)
    heavy_indices = [i for i, a in enumerate(original_atoms) if a.element.atomic_number != 1]
    for i in heavy_indices:
        force.addParticle(prepared_indices[_atom_identity(original_atoms[i])], [0., *original_nm[i]])
    system.addForce(force)
    integrator = mm.VerletIntegrator(0.001 * unit.picosecond)
    properties = {} if config.device_index is None else {"DeviceIndex": config.device_index}
    context = mm.Context(system, integrator, platform, properties)
    try:
        context.setPositions(prepared.positions)
        initial_energy, _ = _physical_state(context, unit)
        _emit_progress(
            progress_callback,
            MinimizationProgressEvent(
                "preparation",
                "completed",
                time.perf_counter() - preparation_started,
            ),
        )
        for stage in config.stages:
            constants = _restraint_constants(original_atoms, ligand_id, pocket, stage)
            for slot, i in enumerate(heavy_indices):
                force.setParticleParameters(slot, prepared_indices[_atom_identity(original_atoms[i])],
                                            [constants[i], *original_nm[i]])
            force.updateParametersInContext(context)
            stage_started = time.perf_counter()
            _emit_progress(
                progress_callback,
                MinimizationProgressEvent(stage.name, "started"),
            )
            mm.LocalEnergyMinimizer.minimize(
                context, config.tolerance_kj_mol_nm * unit.kilojoule_per_mole / unit.nanometer,
                stage.max_iterations,
            )
            final_energy, final_nm = _physical_state(context, unit)
            _emit_progress(
                progress_callback,
                MinimizationProgressEvent(
                    stage.name,
                    "completed",
                    time.perf_counter() - stage_started,
                ),
            )
    finally:
        del context
        del integrator
    finalization_started = time.perf_counter()
    _emit_progress(
        progress_callback,
        MinimizationProgressEvent("finalization", "started"),
    )
    output_nm = _original_positions(original_atoms, prepared_atoms, final_nm)
    final_A = output_nm * 10.0
    backbone_rmsd, ligand_rmsd = align_on_backbone(
        original_A[backbone_indices], final_A[backbone_indices],
        original_A[ligand_indices], final_A[ligand_indices],
    )
    after = vdw_overlap_metrics(final_A[protein_indices], protein_elements,
                                final_A[ligand_indices], ligand_elements)
    _write_original_cif(original_topology, output_nm, input_records, request.output_path, app, unit)
    from importlib.metadata import PackageNotFoundError, version
    versions = []
    for name, module in (("openmm", mm), ("openff-toolkit", stack.toolkit),
                         ("openmmforcefields", stack.generators)):
        module_version = getattr(module, "__version__", None)
        if module_version is None:
            try:
                module_version = version(name)
            except PackageNotFoundError:
                module_version = "unknown"
        versions.append((name, str(module_version)))
    _emit_progress(
        progress_callback,
        MinimizationProgressEvent(
            "finalization",
            "completed",
            time.perf_counter() - finalization_started,
        ),
    )
    return BackendMinimizationResult(
        request.output_path, initial_energy, final_energy, backbone_rmsd, ligand_rmsd,
        before.clash_pair_count, after.clash_pair_count,
        before.maximum_overlap_A, after.maximum_overlap_A, tuple(versions),
    )
