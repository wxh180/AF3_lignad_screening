from dataclasses import FrozenInstanceError
import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest
from types import SimpleNamespace

import restrained_openmm as rom


OPENMM_STACK_AVAILABLE = all(
    importlib.util.find_spec(name) is not None
    for name in ("openmm", "openmmforcefields", "openff.toolkit")
)


def test_openmm_environment_contains_protocol_dependencies():
    """Catch a published optional environment that cannot run the protocol."""
    environment = Path(__file__).parents[1] / "environment-openmm.yml"
    text = environment.read_text(encoding="utf-8")

    for dependency in (
        "openmm",
        "openmmforcefields",
        "openff-toolkit",
        "rdkit",
        "ambertools",
        "biopython",
        "numpy",
        "matplotlib",
        "pytest",
    ):
        assert dependency in text


@pytest.mark.skipif(
    not OPENMM_STACK_AVAILABLE,
    reason="optional OpenMM/OpenFF stack is not installed",
)
def test_real_openmm_backend_minimizes_complete_tiny_complex(tmp_path):
    """Smoke-test the physical backend on CPU when its optional stack exists."""
    from Bio.PDB import MMCIFIO

    cif_path = tmp_path / "tiny.cif"
    atom_names = ["N", "CA", "C", "O", "OXT", "CB", "C1", "O1"]
    residue_names = ["ALA"] * 6 + ["LIG"] * 2
    chain_ids = ["A"] * 6 + ["L"] * 2
    elements = ["N", "C", "C", "O", "O", "C", "C", "O"]
    coordinates_A = np.array(
        [
            [-1.20, 0.00, 0.00],
            [0.00, 0.00, 0.00],
            [1.50, 0.00, 0.00],
            [2.15, 1.00, 0.00],
            [2.15, -1.00, 0.00],
            [0.00, 1.50, 0.00],
            [4.50, 0.00, 0.00],
            [5.70, 0.00, 0.00],
        ]
    )
    data = {
        "data_": "tiny",
        "_atom_site.group_PDB": ["ATOM"] * 6 + ["HETATM"] * 2,
        "_atom_site.id": [str(index) for index in range(1, 9)],
        "_atom_site.type_symbol": elements,
        "_atom_site.label_atom_id": atom_names,
        "_atom_site.label_alt_id": ["."] * 8,
        "_atom_site.label_comp_id": residue_names,
        "_atom_site.label_asym_id": chain_ids,
        "_atom_site.label_entity_id": ["1"] * 6 + ["2"] * 2,
        "_atom_site.label_seq_id": ["1"] * 6 + ["."] * 2,
        "_atom_site.pdbx_PDB_ins_code": ["?"] * 8,
        "_atom_site.Cartn_x": [str(value) for value in coordinates_A[:, 0]],
        "_atom_site.Cartn_y": [str(value) for value in coordinates_A[:, 1]],
        "_atom_site.Cartn_z": [str(value) for value in coordinates_A[:, 2]],
        "_atom_site.occupancy": ["1.0"] * 8,
        "_atom_site.B_iso_or_equiv": ["0.0"] * 8,
        "_atom_site.auth_seq_id": ["1"] * 8,
        "_atom_site.auth_comp_id": residue_names,
        "_atom_site.auth_asym_id": chain_ids,
        "_atom_site.auth_atom_id": atom_names,
        "_atom_site.pdbx_PDB_model_num": ["1"] * 8,
        "_struct_conn.id": ["ligand_bond"],
        "_struct_conn.conn_type_id": ["covale"],
        "_struct_conn.ptnr1_label_asym_id": ["L"],
        "_struct_conn.ptnr1_label_comp_id": ["LIG"],
        "_struct_conn.ptnr1_label_seq_id": ["."],
        "_struct_conn.ptnr1_label_atom_id": ["C1"],
        "_struct_conn.ptnr1_symmetry": ["1_555"],
        "_struct_conn.ptnr2_label_asym_id": ["L"],
        "_struct_conn.ptnr2_label_comp_id": ["LIG"],
        "_struct_conn.ptnr2_label_seq_id": ["."],
        "_struct_conn.ptnr2_label_atom_id": ["O1"],
        "_struct_conn.ptnr2_symmetry": ["1_555"],
        "_struct_conn.pdbx_dist_value": ["1.20"],
    }
    writer = MMCIFIO()
    writer.set_dict(data)
    writer.save(str(cif_path))
    output_path = tmp_path / "tiny_minimized.cif"
    request = rom.MinimizationRequest(
        cif_path=cif_path,
        output_path=output_path,
        structure_id="tiny/seed-1/tiny.cif",
        compound="tiny",
        ligand_resname="LIG",
        protein_chains=("A",),
        ligand_chain="L",
        ligand_smiles="CO",
        ligand_file=None,
    )

    result = rom.minimize_structure(
        request, rom.MinimizationConfig(platform="CPU"), tmp_path / "cache.json"
    )

    assert result.output_path.stat().st_size > 0
    assert math.isfinite(result.initial_potential_kcal_mol)
    assert math.isfinite(result.final_potential_kcal_mol)
    assert result.after_vdw_clash_pair_count <= result.before_vdw_clash_pair_count


def test_dependency_loader_names_the_optional_environment():
    """Catch missing dependencies leaking an unactionable ImportError."""
    def missing(name):
        raise ModuleNotFoundError(name)

    with pytest.raises(rom.MinimizationDependencyError, match="environment-openmm.yml"):
        rom.load_openmm_stack(import_module=missing)


def test_dependency_module_import_does_not_load_optional_packages():
    """Catch eager imports breaking contact-only installations."""
    import subprocess
    import sys

    completed = subprocess.run(
        [sys.executable, "-c", "import sys; import restrained_openmm; "
         "assert not any(n.split('.')[0] in {'openmm', 'openff', 'openmmforcefields'} "
         "for n in sys.modules)"], capture_output=True, text=True,
    )
    assert completed.returncode == 0, completed.stderr


def backend_request(tmp_path, **changes):
    from dataclasses import replace
    request = rom.MinimizationRequest(
        tmp_path / "input.cif", tmp_path / "result.cif", "C/model.cif", "C",
        "LIG", ("A",), "L", "CO", None,
    )
    request = replace(request, **changes)
    if not request.cif_path.exists():
        write_backend_cif(request.cif_path, backend_topology())
    return request


def write_backend_cif(path, topology, positions_nm=None):
    """Faithful atom_site codec: both namespaces start equal, like PDBxFile.writeFile."""
    from Bio.PDB import MMCIFIO
    atoms = list(topology.atoms())
    if positions_nm is None:
        positions_nm = np.zeros((len(atoms), 3))
    data = {"data_": "fixture"}
    fields = {
        "id": [str(i+1) for i in range(len(atoms))],
        "group_PDB": ["ATOM" if a.residue.name in {"ALA", "ILE", "HIS"} else "HETATM" for a in atoms],
        "type_symbol": [a.element.symbol for a in atoms],
        "label_atom_id": [a.name for a in atoms],
        "label_comp_id": [a.residue.name for a in atoms],
        "label_asym_id": [a.residue.chain.id for a in atoms],
        "label_seq_id": [a.residue.id for a in atoms],
        "pdbx_PDB_ins_code": [a.residue.insertionCode or "?" for a in atoms],
        "label_alt_id": ["."]*len(atoms), "occupancy": ["1"]*len(atoms),
        "B_iso_or_equiv": ["0"]*len(atoms), "pdbx_PDB_model_num": ["1"]*len(atoms),
    }
    for name in ("atom_id", "comp_id", "asym_id", "seq_id"):
        fields[f"auth_{name}"] = list(fields[f"label_{name}"])
    for i, axis in enumerate("xyz"):
        fields[f"Cartn_{axis}"] = [str(float(p[i])*10) for p in positions_nm]
    data.update({f"_atom_site.{name}": values for name, values in fields.items()})
    writer = MMCIFIO()
    writer.set_dict(data)
    writer.save(str(path) if isinstance(path, Path) else path)


def read_backend_output(path):
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    data = MMCIF2Dict(str(path))
    return {"names": data["_atom_site.label_atom_id"],
            "positions_nm": np.array([data[f"_atom_site.Cartn_{a}"] for a in "xyz"], dtype=float).T*.1,
            "data": data}


class BackendResidue(SimpleNamespace):
    def atoms(self):
        return iter(self._atoms)


class BackendTopology:
    """Only the unavailable OpenMM topology container, not selection logic."""
    def __init__(self):
        self._residues, self._atoms, self._bonds = [], [], []
        self.box_vectors = None

    def setPeriodicBoxVectors(self, vectors):
        self.box_vectors = vectors

    def getPeriodicBoxVectors(self):
        return self.box_vectors

    def residue(self, chain, name, number, atoms, insertion=""):
        residue = BackendResidue(chain=SimpleNamespace(id=chain), name=name,
                                  id=number, insertionCode=insertion,
                                  index=len(self._residues), _atoms=[])
        self._residues.append(residue)
        for name, symbol in atoms:
            atom = SimpleNamespace(name=name, element=SimpleNamespace(
                symbol=symbol, atomic_number={"C": 6, "N": 7, "O": 8, "H": 1, "Na": 11}[symbol]),
                residue=residue, index=len(self._atoms), id=str(len(self._atoms)+1))
            self._atoms.append(atom)
            residue._atoms.append(atom)
        return residue

    def atoms(self):
        return iter(self._atoms)

    def residues(self):
        return iter(self._residues)

    def bonds(self):
        return iter(self._bonds)


def backend_topology():
    topology = BackendTopology()
    topology.residue("A", "ALA", "10", [("N", "N"), ("CA", "C"),
                                              ("C", "C"), ("O", "O"), ("CB", "C")])
    topology.residue("A", "ALA", "20", [("CB", "C")], insertion="B")
    ligand = topology.residue("L", "LIG", "9", [("O1", "O"), ("C1", "C")])
    topology._bonds.append(tuple(ligand.atoms()))
    topology.residue("B", "ALA", "1", [("CA", "C")])
    topology.residue("A", "HOH", "2", [("O", "O")])
    topology.residue("A", "NA", "3", [("NA", "Na")])
    topology.residue("A", "XYZ", "4", [("C1", "C")])
    return topology


def backend_molecule():
    # Reversed heavy-atom order relative to CIF, plus explicit graph hydrogens.
    atoms = [SimpleNamespace(atomic_number=n, molecule_atom_index=i)
             for i, n in enumerate([6, 8, 1, 1, 1, 1])]
    bonds = [SimpleNamespace(atom1_index=i, atom2_index=j)
             for i, j in [(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)]]
    return SimpleNamespace(atoms=atoms, bonds=bonds, conformers=None)


def test_backend_selection_retains_only_selected_protein_and_one_ligand(tmp_path):
    """Catch ligand selection incorrectly inheriting the protein chain filter."""
    topology = backend_topology()
    retained, ligand = rom._select_residues(topology, backend_request(tmp_path))
    assert [(r.chain.id, r.name, r.id) for r in retained] == [
        ("A", "ALA", "10"), ("A", "ALA", "20"), ("L", "LIG", "9")]
    assert ligand.id == "9"


@pytest.mark.parametrize("copies", [0, 2])
def test_backend_requires_exactly_one_ligand(tmp_path, copies):
    """Catch accepting absent or ambiguous ligand selections."""
    topology = backend_topology()
    topology._residues = [r for r in topology.residues() if r.name != "LIG"]
    for i in range(copies):
        topology.residue("L", "LIG", str(i), [("C1", "C")])
    with pytest.raises(rom.MinimizationBackendError, match="exactly one ligand residue"):
        rom._select_residues(topology, backend_request(tmp_path))


def test_backend_identity_mapping_survives_reordering_and_omits_transient_atoms():
    """Catch positional-index mapping silently assigning another atom's position."""
    original = backend_topology()
    prepared = backend_topology()
    prepared.residue("A", "ALA", "30", [("H", "H")])
    atoms = list(prepared.atoms())[::-1]
    positions = np.arange(len(atoms)*3).reshape(-1, 3)
    mapped = rom._original_positions(list(original.atoms()), atoms, positions)
    np.testing.assert_array_equal(mapped[0], [36, 37, 38])
    np.testing.assert_array_equal(mapped[-1], [3, 4, 5])
    assert mapped.shape == (12, 3)


def test_backend_identity_mapping_rejects_lost_or_duplicate_identity():
    """Catch a removed original atom or duplicate identity being silently accepted."""
    atoms = list(backend_topology().atoms())
    with pytest.raises(rom.MinimizationBackendError, match="missing original atom"):
        rom._original_positions(atoms, atoms[1:], np.zeros((11, 3)))
    with pytest.raises(rom.MinimizationBackendError, match="duplicate atom identity"):
        rom._original_positions(atoms, [*atoms, atoms[0]], np.zeros((13, 3)))


def test_backend_restraint_assignment_uses_original_pocket_and_heavy_atoms():
    """Catch free backbone, restrained pocket side chains, or restrained hydrogens."""
    topology = backend_topology()
    ligand = list(topology.residues())[2]
    topology.residue("A", "ALA", "10", [("H", "H")])
    atoms = [*list(topology.atoms())[:8], list(topology.atoms())[-1]]
    pocket = frozenset({rom._residue_identity(list(topology.residues())[0])})
    expected = [
        [4184, 4184, 4184, 4184, 4184, 4184, 4184, 4184, 0],
        [4184, 4184, 4184, 4184, 0, 4184, 418.4, 418.4, 0],
        [836.8, 836.8, 836.8, 836.8, 0, 836.8, 83.68, 83.68, 0],
    ]
    for stage, want in zip(rom.DEFAULT_STAGES, expected):
        assert rom._restraint_constants(atoms, rom._residue_identity(ligand), pocket, stage) == pytest.approx(want)


def test_backend_ligand_hydrogen_definitions_follow_graph_not_atom_order():
    """Catch attaching ligand hydrogens by atom order instead of chemical graph."""
    topology = backend_topology()
    ligand = list(topology.residues())[2]
    definitions = rom._ligand_hydrogens(topology, ligand, backend_molecule())
    assert sorted(parent for _, parent in definitions) == ["C1", "C1", "C1", "O1"]
    assert len({name for name, _ in definitions}) == 4


def test_backend_ligand_graph_mismatch_requires_explicit_connectivity():
    """Catch guessing bonds from coordinates or accepting missing heavy atoms."""
    topology = backend_topology()
    topology._bonds.clear()
    with pytest.raises(rom.MinimizationBackendError, match="_struct_conn"):
        rom._ligand_hydrogens(topology, list(topology.residues())[2], backend_molecule())


class BackendQuantity:
    def __init__(self, value):
        self.value = np.asarray(value)

    def value_in_unit(self, unit):
        return self.value / unit


@pytest.fixture
def backend_stack(monkeypatch):
    """Substitute only third-party chemistry/OpenMM execution and file codec."""
    import copy
    stack = SimpleNamespace(energy=float(-12), nonfinite_positions=False,
                            sdf_molecules=[backend_molecule()], fail_template=False,
                            fail_after_stage=0, covalent_to_unselected=False, periodic_input=False)
    unit = SimpleNamespace(nanometer=1., angstrom=.1, picosecond=1.,
                           kilojoule_per_mole=1., kilocalorie_per_mole=4.184,
                           Quantity=lambda values, scale: BackendQuantity(np.asarray(values)*scale))

    class Molecule:
        @staticmethod
        def from_smiles(smiles, *, allow_undefined_stereo):
            if allow_undefined_stereo:
                raise ValueError("stereochemistry must be defined")
            return backend_molecule()

        @staticmethod
        def from_file(path, *, file_format, allow_undefined_stereo):
            if file_format != "SDF" or allow_undefined_stereo:
                raise ValueError("exact SDF chemistry required")
            return stack.sdf_molecules

    class PDBxFile:
        def __init__(self, path):
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            self.topology = backend_topology()
            raw = MMCIF2Dict(str(path))
            chain_field = ("label_asym_id" if len(set(raw["_atom_site.label_asym_id"])) >
                           len(set(raw["_atom_site.auth_asym_id"])) else "auth_asym_id")
            for i, atom in enumerate(self.topology.atoms()):
                atom.residue.chain.id = raw[f"_atom_site.{chain_field}"][i]
                atom.residue.name = raw["_atom_site.auth_comp_id"][i]
                atom.name = raw["_atom_site.auth_atom_id"][i]
                if atom.residue.name in {"HID", "HIE"}:
                    atom.residue.name = "HIS"
                if atom.residue.name == "ILE" and atom.name == "CD":
                    atom.name = "CD1"
            if stack.periodic_input:
                self.topology.setPeriodicBoxVectors(np.eye(3))
            if stack.covalent_to_unselected:
                atoms = list(self.topology.atoms())
                self.topology._bonds.append((atoms[7], atoms[8]))
            self.positions = BackendQuantity(np.array([
                [0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3], [1, 1, 0],
                [30, 0, 0], [3, 0, 0], [4, 0, 0], [50, 0, 0],
                [50, 0, 0], [50, 0, 0], [50, 0, 0]], dtype=float)*.1)

        @staticmethod
        def writeFile(topology, positions, file, keepIds=False):
            if not keepIds:
                raise ValueError("identity preservation required")
            write_backend_cif(file, topology, positions.value_in_unit(1.))

    class Modeller:
        def __init__(self, topology, positions):
            self.topology, self.positions = copy.deepcopy(topology), positions

        def delete(self, residues):
            removed = {rom._residue_identity(r) for r in residues}
            kept = [a for a in self.topology.atoms()
                    if rom._residue_identity(a.residue) not in removed]
            indices = [a.index for a in kept]
            self.topology._residues = [r for r in self.topology.residues()
                                      if rom._residue_identity(r) not in removed]
            self.topology._atoms = kept
            kept_ids = {id(a) for a in kept}
            self.topology._bonds = [(a, b) for a, b in self.topology.bonds()
                                    if id(a) in kept_ids and id(b) in kept_ids]
            self.positions = BackendQuantity(self.positions.value[indices])
            for i, atom in enumerate(kept):
                atom.index = i
            for i, residue in enumerate(self.topology.residues()):
                residue.index = i

        def addHydrogens(self, forcefield, *, pH, variants, platform):
            if pH != 7.4 or len(variants[2]) != 4:
                raise ValueError("hydrogen preparation contract")
            # Insert a transient hydrogen first, shifting every original index.
            residue = list(self.topology.residues())[0]
            hydrogen = SimpleNamespace(name="H", element=SimpleNamespace(symbol="H", atomic_number=1),
                                       residue=residue, id="99", index=0)
            self.topology._atoms.insert(0, hydrogen)
            self.positions = BackendQuantity(np.vstack([[0, 0, .1], self.positions.value]))
            for i, atom in enumerate(self.topology.atoms()):
                atom.index = i

    class SystemGenerator:
        def __init__(self, *, forcefields, small_molecule_forcefield, molecules,
                     cache, forcefield_kwargs, nonperiodic_forcefield_kwargs):
            if (forcefields != ["amber/protein.ff14SB.xml"] or
                    small_molecule_forcefield != "gaff-2.2.20" or
                    forcefield_kwargs != {"constraints": "HBonds"} or
                    nonperiodic_forcefield_kwargs != {"nonbondedMethod": "NoCutoff"}):
                raise ValueError("physical model changed")
            self.forcefield = object()

        def create_system(self, topology):
            if topology.getPeriodicBoxVectors() is not None:
                raise ValueError("periodic input incorrectly selects PME")
            if stack.fail_template:
                raise ValueError("No template found: missing heavy atoms")
            return SimpleNamespace(addForce=lambda force: setattr(stack, "force", force))

    class Force:
        def __init__(self, expression):
            self.expression, self.particles, self.parameters = expression, [], []

        def addPerParticleParameter(self, name):
            self.parameters.append(name)

        def setForceGroup(self, group):
            self.group = group

        def addParticle(self, index, parameters):
            self.particles.append((index, list(parameters)))

        def setParticleParameters(self, slot, index, parameters):
            self.particles[slot] = (index, list(parameters))

        def updateParametersInContext(self, context):
            context.restraints = copy.deepcopy(self.particles)

    class Context:
        def __init__(self, system, integrator, platform, properties):
            self.stage, self.restraints = 0, []

        def setPositions(self, positions):
            self.positions = positions

        def getState(self, *, getEnergy=False, getPositions=False, groups=-1):
            # Including group 31 contaminates the user-visible physical energy.
            energy = (stack.energy if self.stage >= stack.fail_after_stage else -12)
            if groups != {0} and groups != 1:
                energy = 9999
            values = self.positions.value.copy()
            if self.stage:
                values += [0, .02, 0]
            if stack.nonfinite_positions and self.stage >= stack.fail_after_stage:
                values[0, 0] = np.nan
            return SimpleNamespace(getPotentialEnergy=lambda: BackendQuantity(energy),
                                   getPositions=lambda asNumpy: BackendQuantity(values))

    class Minimizer:
        @staticmethod
        def minimize(context, tolerance, maxIterations):
            # Consume assigned restraints; rejects wrong units/groups/stage state.
            expected = [
                [4184]*8,
                [4184, 4184, 4184, 4184, 0, 4184, 418.4, 418.4],
                [836.8, 836.8, 836.8, 836.8, 0, 836.8, 83.68, 83.68],
            ][context.stage]
            actual = [p[0] for _, p in context.restraints]
            if (actual != pytest.approx(expected) or stack.force.group != 31 or
                    stack.force.parameters != ["k", "x0", "y0", "z0"] or
                    stack.force.expression != "k*((x-x0)^2+(y-y0)^2+(z-z0)^2)" or
                    tolerance != 10 or maxIterations != [500, 1000, 1000][context.stage]):
                raise ValueError("incorrect restraint protocol")
            if [i for i, _ in context.restraints] != list(range(1, 9)):
                raise ValueError("restrained transient H or stale identity mapping")
            reference_nm = [[0, 0, 0], [.1, 0, 0], [0, .2, 0], [0, 0, .3],
                            [.1, .1, 0], [3, 0, 0], [.3, 0, 0], [.4, 0, 0]]
            if not np.allclose([p[1:] for _, p in context.restraints], reference_nm):
                raise ValueError("restraint targets must remain the original positions in nm")
            context.stage += 1

    stack.unit = unit
    stack.toolkit = SimpleNamespace(Molecule=Molecule, __version__="test-openff")
    stack.app = SimpleNamespace(PDBxFile=PDBxFile, Modeller=Modeller, HBonds="HBonds", NoCutoff="NoCutoff")
    stack.generators = SimpleNamespace(SystemGenerator=SystemGenerator, __version__="test-forcefields")
    stack.openmm = SimpleNamespace(CustomExternalForce=Force, Context=Context,
                                  LocalEnergyMinimizer=Minimizer, VerletIntegrator=lambda dt: dt,
                                  Platform=SimpleNamespace(getPlatformByName=lambda name: name),
                                  __version__="test-openmm")
    monkeypatch.setattr(rom, "load_openmm_stack", lambda: stack)
    return stack


def test_backend_minimizes_with_physical_energy_and_serializes_original_atoms(tmp_path, backend_stack):
    """Catch restrained energy reporting, stale positions, or transient H publication."""
    request = backend_request(tmp_path)
    result = rom.minimize_structure(request, rom.MinimizationConfig(), tmp_path / "cache.json")
    output = read_backend_output(request.output_path)
    assert output["names"] == ["N", "CA", "C", "O", "CB", "CB", "O1", "C1"]
    assert output["positions_nm"][0] == pytest.approx([0, .02, 0])
    assert output["positions_nm"][-1] == pytest.approx([.4, .02, 0])
    assert result.initial_potential_kcal_mol == pytest.approx(-12 / 4.184)
    assert result.final_potential_kcal_mol == pytest.approx(-12 / 4.184)
    assert result.backbone_rmsd_A == pytest.approx(0, abs=1e-12)
    assert result.ligand_rmsd_A == pytest.approx(0, abs=1e-12)
    assert result.before_vdw_clash_pair_count == result.after_vdw_clash_pair_count == 3
    assert result.before_maximum_vdw_overlap_A == pytest.approx(1.22)
    assert result.after_maximum_vdw_overlap_A == pytest.approx(1.22)
    assert dict(result.package_versions)["openmm"] == "test-openmm"


def test_backend_emits_progress_around_existing_three_stages(tmp_path, backend_stack):
    request = backend_request(tmp_path)
    events = []

    result = rom._minimize_structure(
        request,
        rom.MinimizationConfig(),
        tmp_path / "cache.json",
        backend_stack,
        progress_callback=events.append,
    )

    assert result.output_path == request.output_path
    assert [(event.operation, event.state) for event in events] == [
        ("preparation", "started"),
        ("preparation", "completed"),
        ("hydrogen_relaxation", "started"),
        ("hydrogen_relaxation", "completed"),
        ("pocket_relaxation", "started"),
        ("pocket_relaxation", "completed"),
        ("gentle_relaxation", "started"),
        ("gentle_relaxation", "completed"),
        ("finalization", "started"),
        ("finalization", "completed"),
    ]
    completed = [event for event in events if event.state == "completed"]
    assert all(event.elapsed_seconds is not None for event in completed)
    assert all(event.elapsed_seconds >= 0 for event in completed)


def test_backend_progress_callback_failure_cannot_fail_minimization(
    tmp_path, backend_stack
):
    request = backend_request(tmp_path)

    def broken_callback(event):
        raise RuntimeError("display failure")

    result = rom._minimize_structure(
        request,
        rom.MinimizationConfig(),
        tmp_path / "cache.json",
        backend_stack,
        progress_callback=broken_callback,
    )

    assert result.output_path.exists()


@pytest.mark.parametrize("smiles, file", [(None, None), ("CO", Path("ligand.sdf"))])
def test_backend_requires_exactly_one_chemistry_source(tmp_path, backend_stack, smiles, file):
    """Catch silently selecting one conflicting chemistry source."""
    with pytest.raises(rom.MinimizationBackendError, match="exactly one.*chemistry"):
        rom.minimize_structure(backend_request(tmp_path, ligand_smiles=smiles, ligand_file=file),
                               rom.MinimizationConfig(), tmp_path / "cache.json")


@pytest.mark.parametrize("count", [0, 2])
def test_backend_sdf_requires_one_molecule(tmp_path, backend_stack, count):
    """Catch dropping molecules from multi-molecule SDFs or accepting empty files."""
    backend_stack.sdf_molecules = [backend_molecule() for _ in range(count)]
    with pytest.raises(rom.MinimizationBackendError, match="one molecule"):
        rom.minimize_structure(backend_request(tmp_path, ligand_smiles=None, ligand_file=Path("ligand.sdf")),
                               rom.MinimizationConfig(), tmp_path / "cache.json")


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), float("-inf")])
def test_backend_nonfinite_energy_is_an_error(tmp_path, backend_stack, bad):
    """Catch non-finite potential energy reaching a successful QC result."""
    backend_stack.energy = bad
    with pytest.raises(rom.MinimizationBackendError, match="finite"):
        rom.minimize_structure(backend_request(tmp_path), rom.MinimizationConfig(), tmp_path / "cache.json")
    assert not (tmp_path / "result.cif").exists()


def test_backend_nonfinite_positions_are_an_error(tmp_path, backend_stack):
    """Catch non-finite transient coordinates being hidden by output filtering."""
    backend_stack.nonfinite_positions = True
    with pytest.raises(rom.MinimizationBackendError, match="finite"):
        rom.minimize_structure(backend_request(tmp_path), rom.MinimizationConfig(), tmp_path / "cache.json")
    assert not (tmp_path / "result.cif").exists()


def test_backend_template_error_identifies_structure_and_chemistry(tmp_path, backend_stack):
    """Catch unhelpful raw template/missing-heavy-atom errors escaping the boundary."""
    backend_stack.fail_template = True
    with pytest.raises(rom.MinimizationBackendError, match="C/model.cif.*SMILES.*missing heavy atoms"):
        rom.minimize_structure(backend_request(tmp_path), rom.MinimizationConfig(), tmp_path / "cache.json")


def test_backend_rejects_nonfinite_sdf_conformer(tmp_path, backend_stack):
    """Catch ignoring invalid supplied SDF coordinates even when CIF coordinates are finite."""
    backend_stack.sdf_molecules[0].conformers = [SimpleNamespace(magnitude=np.array([[np.nan, 0, 0]]))]
    with pytest.raises(rom.MinimizationBackendError, match="SDF.*finite"):
        rom.minimize_structure(backend_request(tmp_path, ligand_smiles=None, ligand_file=Path("ligand.sdf")),
                               rom.MinimizationConfig(), tmp_path / "cache.json")


def test_backend_rejects_covalent_ligand_before_unselected_chain_deletion(tmp_path, backend_stack):
    """Catch topology filtering hiding an unsupported covalent ligand bond."""
    backend_stack.covalent_to_unselected = True
    with pytest.raises(rom.MinimizationBackendError, match="covalently attached"):
        rom.minimize_structure(backend_request(tmp_path), rom.MinimizationConfig(), tmp_path / "cache.json")


@pytest.mark.parametrize("field, bad", [("energy", float("nan")), ("nonfinite_positions", True)])
def test_backend_stops_on_nonfinite_values_in_intermediate_stage(tmp_path, backend_stack, field, bad):
    """Catch checking only initial/final states and overlooking a failed stage."""
    backend_stack.fail_after_stage = 2
    setattr(backend_stack, field, bad)
    with pytest.raises(rom.MinimizationBackendError, match="finite"):
        rom.minimize_structure(backend_request(tmp_path), rom.MinimizationConfig(), tmp_path / "cache.json")
    assert not (tmp_path / "result.cif").exists()


def test_backend_crystal_cell_does_not_enable_periodic_electrostatics(tmp_path, backend_stack):
    """Catch crystal mmCIF box vectors silently selecting SystemGenerator's PME default."""
    backend_stack.periodic_input = True
    result = rom.minimize_structure(backend_request(tmp_path), rom.MinimizationConfig(), tmp_path / "cache.json")
    assert result.final_potential_kcal_mol == pytest.approx(-12 / 4.184)


def rewrite_backend_input(path, change):
    from Bio.PDB import MMCIFIO
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    data = MMCIF2Dict(str(path))
    change(data)
    writer = MMCIFIO()
    writer.set_dict(data)
    writer.save(str(path))


def test_backend_preserves_author_chain_with_distinct_label_chains(tmp_path, backend_stack):
    """Catch OpenMM label-chain promotion breaking author-chain selection/publication."""
    request = backend_request(tmp_path, protein_chains=("X",), ligand_chain="X")
    def change(data):
        data["_atom_site.auth_asym_id"] = ["X" if c in {"A", "L"} else c
                                           for c in data["_atom_site.label_asym_id"]]
    rewrite_backend_input(request.cif_path, change)
    result = rom.minimize_structure(request, rom.MinimizationConfig(), tmp_path / "cache.json")
    data = read_backend_output(result.output_path)["data"]
    assert data["_atom_site.auth_asym_id"] == ["X"]*8
    assert data["_atom_site.label_asym_id"] == ["A"]*6 + ["L"]*2
    assert data["_atom_site.auth_seq_id"] == ["10"]*5 + ["20", "9", "9"]


def test_backend_restores_original_atom_name_after_openmm_normalization(tmp_path, backend_stack):
    """Catch publishing normalized ILE/CD1 instead of the selected input's ILE/CD."""
    request = backend_request(tmp_path)
    def change(data):
        for namespace in ("label", "auth"):
            data[f"_atom_site.{namespace}_comp_id"][:5] = ["ILE"]*5
            data[f"_atom_site.{namespace}_atom_id"][4] = "CD"
    rewrite_backend_input(request.cif_path, change)
    result = rom.minimize_structure(request, rom.MinimizationConfig(), tmp_path / "cache.json")
    output = read_backend_output(result.output_path)
    assert output["names"][:5] == ["N", "CA", "C", "O", "CD"]
    assert output["data"]["_atom_site.label_comp_id"][:5] == ["ILE"]*5
    assert output["positions_nm"][4] == pytest.approx([.1, .12, 0])


@pytest.mark.parametrize("raw_name", ["HID", "HIE"])
def test_backend_restores_selected_residue_name_after_normalization(tmp_path, backend_stack, raw_name):
    """Catch selecting or publishing by OpenMM's HIS replacement instead of raw name."""
    # These names are not protein residues in the contact pipeline's is_aa
    # classification, so exercise them as an explicitly selected residue.
    request = backend_request(tmp_path, ligand_resname=raw_name)
    def change(data):
        for namespace in ("label", "auth"):
            data[f"_atom_site.{namespace}_comp_id"][6:8] = [raw_name]*2
    rewrite_backend_input(request.cif_path, change)
    result = rom.minimize_structure(request, rom.MinimizationConfig(), tmp_path / "cache.json")
    data = read_backend_output(result.output_path)["data"]
    assert data["_atom_site.label_comp_id"][6:] == [raw_name]*2
    assert data["_atom_site.auth_comp_id"][6:] == [raw_name]*2


def test_backend_raw_identity_round_trip_with_real_openmm_codec(tmp_path):
    """Exercise actual parser normalization and the publication writer when installed."""
    app = pytest.importorskip("openmm.app")
    from openmm import unit
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    request = backend_request(tmp_path, protein_chains=("X",), ligand_chain="X")
    def change(data):
        for namespace in ("label", "auth"):
            data[f"_atom_site.{namespace}_comp_id"][:5] = ["ILE"]*5
            data[f"_atom_site.{namespace}_atom_id"][4] = "CD"
            data[f"_atom_site.{namespace}_comp_id"][5] = "HID"
        data["_atom_site.auth_asym_id"] = ["X" if c in {"A", "L"} else c
                                           for c in data["_atom_site.label_asym_id"]]
    rewrite_backend_input(request.cif_path, change)
    records = rom._read_input_atom_records(request.cif_path)
    parsed = app.PDBxFile(str(request.cif_path))
    assert list(parsed.topology.atoms())[4].name == "CD1"
    assert list(parsed.topology.residues())[1].name == "HIS"
    rom._bind_input_atom_records(parsed.topology, records)
    retained, ligand = rom._select_residues(parsed.topology, request, records)
    assert [(r.chain.id, r.name) for r in retained] == [("A", "ILE"), ("L", "LIG")]
    rom._write_original_cif(parsed.topology, parsed.positions.value_in_unit(unit.nanometer),
                            records, request.output_path, app, unit)
    raw = MMCIF2Dict(str(request.cif_path))
    written = MMCIF2Dict(str(request.output_path))
    for field in ("label_atom_id", "label_comp_id", "label_asym_id", "label_seq_id",
                  "auth_atom_id", "auth_comp_id", "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code"):
        assert written[f"_atom_site.{field}"] == raw[f"_atom_site.{field}"]
    # Restoration must not mutate the preparation topology used for mapping.
    assert list(parsed.topology.atoms())[4].name == "CD1"
    assert list(parsed.topology.residues())[1].name == "HIS"


def test_kcal_per_A2_converts_to_kj_per_nm2():
    """Catch a missed Angstrom-to-nanometre squared conversion."""
    assert rom.kcal_per_A2_to_kj_per_nm2(1.0) == pytest.approx(418.4)


def test_protocol_has_three_bounded_stages():
    """Catch drift in the stable restrained-minimization protocol contract."""
    config = rom.MinimizationConfig()

    assert config.protocol == "openmm_restrained_v1"
    assert config.protein_forcefield == "amber/protein.ff14SB.xml"
    assert config.ligand_forcefield == "gaff-2.2.20"
    assert config.ph == 7.4
    assert config.pocket_radius_A == 6.0
    assert config.tolerance_kj_mol_nm == 10.0
    assert config.platform == "CPU"
    assert config.device_index is None
    assert config.max_backbone_rmsd_A == 0.5
    assert config.max_ligand_rmsd_A == 1.5
    assert rom.DEFAULT_STAGES == (
        rom.MinimizationStage("hydrogen_relaxation", 500, 10.0, 10.0, 10.0, 10.0),
        rom.MinimizationStage("pocket_relaxation", 1000, 10.0, 10.0, 0.0, 1.0),
        rom.MinimizationStage("gentle_relaxation", 1000, 2.0, 2.0, 0.0, 0.2),
    )
    assert config.stages == rom.DEFAULT_STAGES


def test_vdw_overlap_boundary_is_inclusive():
    """Catch an overlap comparison that excludes the documented boundary."""
    result = rom.vdw_overlap_metrics(
        np.array([[0.0, 0.0, 0.0]]), ["C"],
        np.array([[3.0, 0.0, 0.0]]), ["C"],
    )

    assert result.clash_pair_count == 1
    assert result.maximum_overlap_A == pytest.approx(0.4)
    assert result.minimum_pair_distance_A == pytest.approx(3.0)


def test_vdw_overlap_rejects_unsupported_element():
    """Catch silently assigning a generic radius to unsupported chemistry."""
    with pytest.raises(ValueError, match="unsupported element"):
        rom.vdw_overlap_metrics(
            np.array([[0.0, 0.0, 0.0]]), ["Zn"],
            np.array([[3.0, 0.0, 0.0]]), ["C"],
        )


def test_vdw_overlap_normalizes_standard_element_case():
    """Catch treating standard Biopython-style element case as unknown."""
    result = rom.vdw_overlap_metrics(
        np.array([[0.0, 0.0, 0.0]]), ["CL"],
        np.array([[3.0, 0.0, 0.0]]), ["cl"],
    )

    assert result.clash_pair_count == 1
    assert result.maximum_overlap_A == pytest.approx(0.5)


def test_alignment_removes_rigid_backbone_motion_without_mutating_inputs():
    """Catch a Kabsch transform that does not map row-vector coordinates."""
    reference_backbone = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    reference_ligand = np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]])
    rotation = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([4.0, -3.0, 2.0])
    mobile_backbone = reference_backbone @ rotation + translation
    mobile_ligand = reference_ligand @ rotation + translation
    original_mobile_backbone = mobile_backbone.copy()
    original_mobile_ligand = mobile_ligand.copy()

    backbone_rmsd, ligand_rmsd = rom.align_on_backbone(
        reference_backbone, mobile_backbone, reference_ligand, mobile_ligand
    )

    assert backbone_rmsd == pytest.approx(0.0, abs=1e-12)
    assert ligand_rmsd == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_array_equal(mobile_backbone, original_mobile_backbone)
    np.testing.assert_array_equal(mobile_ligand, original_mobile_ligand)


def test_alignment_reports_ligand_motion_after_protein_alignment():
    """Catch calculating ligand RMSD before applying the protein alignment."""
    reference_backbone = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]]
    )
    reference_ligand = np.array([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]])
    rotation = np.array([[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([4.0, -3.0, 2.0])
    mobile_backbone = reference_backbone @ rotation + translation
    mobile_ligand = reference_ligand @ rotation + translation
    mobile_ligand[0] += np.array([0.5, 0.0, 0.0])

    backbone_rmsd, ligand_rmsd = rom.align_on_backbone(
        reference_backbone, mobile_backbone, reference_ligand, mobile_ligand
    )

    assert backbone_rmsd == pytest.approx(0.0, abs=1e-12)
    assert ligand_rmsd == pytest.approx(np.sqrt(0.125), abs=1e-12)


def test_select_pocket_residues_includes_radius_boundary():
    """Catch excluding a protein residue exactly at the pocket radius."""
    residues = rom.select_pocket_residues(
        np.array([[0.0, 0.0, 0.0], [20.0, 0.0, 0.0]]),
        ["A:10:ALA", "A:20:GLY"],
        np.array([[6.0, 0.0, 0.0]]),
        radius_A=6.0,
    )

    assert residues == frozenset({"A:10:ALA"})


@pytest.mark.parametrize("accepted, bucket", [(True, "accepted"), (False, "rejected")])
def test_candidate_relative_path_is_deterministic_and_contained(accepted, bucket):
    """Catch output paths that omit status routing or alter the structure ID."""
    relative_path = rom.candidate_relative_path("C19/seed-1/model.cif", accepted)

    assert relative_path == Path(
        f"minimized_structures/{bucket}/C19/seed-1/model_minimized.cif"
    )


@pytest.mark.parametrize("structure_id", ["/tmp/model.cif", "C19/../model.cif"])
def test_candidate_relative_path_rejects_unsafe_structure_id(structure_id):
    """Catch an absolute or traversing ID escaping the output directory."""
    with pytest.raises(ValueError, match="structure_id"):
        rom.candidate_relative_path(structure_id, accepted=True)


def test_acceptance_allows_metrics_at_all_gates():
    """Catch rejecting a candidate that exactly meets every acceptance limit."""
    decision = rom.evaluate_acceptance(
        backbone_rmsd_A=0.5,
        ligand_rmsd_A=1.5,
        before_fixed_clashes=2,
        after_fixed_clashes=2,
        before_vdw_clashes=3,
        after_vdw_clashes=3,
        max_backbone_rmsd_A=0.5,
        max_ligand_rmsd_A=1.5,
    )

    assert decision.accepted is True
    assert decision.reasons == ()


def test_acceptance_rejects_ligand_drift():
    """Catch failing to reject ligand RMSD above its configured gate."""
    decision = rom.evaluate_acceptance(
        backbone_rmsd_A=0.1,
        ligand_rmsd_A=1.6,
        before_fixed_clashes=2,
        after_fixed_clashes=0,
        before_vdw_clashes=3,
        after_vdw_clashes=0,
        max_backbone_rmsd_A=0.5,
        max_ligand_rmsd_A=1.5,
    )

    assert not decision.accepted
    assert decision.reasons == ("ligand RMSD 1.600 A exceeds 1.500 A",)


def test_acceptance_reports_all_rejections_in_stable_order():
    """Catch missing acceptance gates or unstable QC reason ordering."""
    decision = rom.evaluate_acceptance(
        backbone_rmsd_A=0.6,
        ligand_rmsd_A=1.6,
        before_fixed_clashes=2,
        after_fixed_clashes=3,
        before_vdw_clashes=3,
        after_vdw_clashes=4,
        max_backbone_rmsd_A=0.5,
        max_ligand_rmsd_A=1.5,
    )

    assert decision.accepted is False
    assert decision.reasons == (
        "backbone RMSD 0.600 A exceeds 0.500 A",
        "ligand RMSD 1.600 A exceeds 1.500 A",
        "fixed clash count increased from 2 to 3",
        "van der Waals clash count increased from 3 to 4",
    )
    with pytest.raises(FrozenInstanceError):
        decision.accepted = True


@pytest.mark.parametrize("nonfinite", [float("nan"), float("inf"), float("-inf")])
@pytest.mark.parametrize(
    "input_name",
    [
        "backbone_rmsd_A",
        "ligand_rmsd_A",
        "before_fixed_clashes",
        "after_fixed_clashes",
        "before_vdw_clashes",
        "after_vdw_clashes",
        "max_backbone_rmsd_A",
        "max_ligand_rmsd_A",
    ],
)
def test_acceptance_rejects_nonfinite_metrics_and_limits(input_name, nonfinite):
    """Catch NaN or infinity silently passing the conservative QC gates."""
    inputs = {
        "backbone_rmsd_A": 0.1,
        "ligand_rmsd_A": 0.2,
        "before_fixed_clashes": 2,
        "after_fixed_clashes": 1,
        "before_vdw_clashes": 3,
        "after_vdw_clashes": 2,
        "max_backbone_rmsd_A": 0.5,
        "max_ligand_rmsd_A": 1.5,
    }
    inputs[input_name] = nonfinite

    with pytest.raises(ValueError, match=rf"{input_name} must be finite"):
        rom.evaluate_acceptance(**inputs)
