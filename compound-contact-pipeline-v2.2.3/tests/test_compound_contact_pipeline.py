import csv
import json
import subprocess
from types import SimpleNamespace
import sys
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

import compound_contact_pipeline as ccp
import restrained_openmm as rom


_PRODUCTION_LOAD_FREESASA = ccp._load_freesasa
_PRODUCTION_CALCULATE_BURIED_SASA = ccp.calculate_buried_sasa


@pytest.fixture(autouse=True)
def _reference_sasa_backend_when_freesasa_is_unavailable(monkeypatch):
    fake_path = Path(__file__).resolve().parent / "fakes"
    existing_pythonpath = __import__("os").environ.get("PYTHONPATH", "")
    monkeypatch.setenv(
        "PYTHONPATH",
        str(fake_path) + ((":" + existing_pythonpath) if existing_pythonpath else ""),
    )
    try:
        _PRODUCTION_LOAD_FREESASA()
    except ccp.PipelineError:
        monkeypatch.setattr(
            ccp,
            "_load_freesasa",
            lambda: SimpleNamespace(__version__="test-biopython-reference"),
        )
        monkeypatch.setattr(
            ccp,
            "calculate_buried_sasa",
            ccp.calculate_buried_sasa_biopython_reference,
        )


def test_version_is_2_2_3():
    assert ccp.__version__ == "2.2.3"




def test_environment_requires_freesasa_python():
    environment = Path("environment-openmm.yml").read_text(encoding="utf-8")
    assert "freesasa-python>=2.2" in environment


def test_sasa_backend_description_reports_version_and_threading(monkeypatch):
    monkeypatch.setattr(
        ccp, "_load_freesasa", lambda: SimpleNamespace(__version__="2.2.1")
    )
    assert ccp._sasa_backend_description() == "FreeSASA 2.2.1 Lee-Richards; 20 slices; 1 thread/worker"

def touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return path


def write_test_cif(path: Path, residue_specs, model_count: int = 1) -> Path:
    structure = Structure.Structure("test")
    serial = 1
    for model_index in range(model_count):
        model = Model.Model(model_index, serial_num=model_index + 1)
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


def test_load_smiles_table_accepts_user_facing_headers(tmp_path):
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\n"
        "EOAI10089186\tCC1=CC=CC=C1\n"
        "EOAI10403337\tCCO\n",
        encoding="utf-8",
    )

    entries = ccp.load_smiles_table(table)

    assert entries == {
        "EOAI10089186": "CC1=CC=CC=C1",
        "EOAI10403337": "CCO",
    }


def test_smiles_table_resolves_compound_id_embedded_in_folder_name(tmp_path):
    root = tmp_path / "B56a_EOAI10089186"
    touch(root / "seed-1" / "model.cif")
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\nEOAI10089186\tCC1=CC=CC=C1\n",
        encoding="utf-8",
    )
    jobs = ccp.build_structure_jobs(root, None, "LIG", ("A",), "F")

    resolved = ccp.apply_smiles_table(jobs, table)

    assert {job.ligand_smiles for job in resolved} == {"CC1=CC=CC=C1"}
    assert {job.ligand_compound_id for job in resolved} == {"EOAI10089186"}
    assert {job.ligand_chemistry_source for job in resolved} == {
        f"smiles_table:{table.resolve().as_posix()}"
    }


def test_explicit_ligand_smiles_takes_priority_over_smiles_table(tmp_path):
    root = tmp_path / "B56a_EOAI10089186"
    touch(root / "seed-1" / "model.cif")
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\nEOAI10089186\tCCO\n",
        encoding="utf-8",
    )
    jobs = ccp.build_structure_jobs(
        root, None, "LIG", ("A",), "F", "CC1=CC=CC=C1", None
    )

    resolved = ccp.apply_smiles_table(jobs, table)

    assert {job.ligand_smiles for job in resolved} == {"CC1=CC=CC=C1"}
    assert {job.ligand_compound_id for job in resolved} == {None}
    assert {job.ligand_chemistry_source for job in resolved} == {
        "command_line_smiles"
    }


def test_smiles_table_rejects_duplicate_compound_ids(tmp_path):
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\nEOAI10089186\tCCO\nEOAI10089186\tCCC\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="duplicate Compound ID"):
        ccp.load_smiles_table(table)


def test_smiles_table_requires_one_unambiguous_match_per_unresolved_job(tmp_path):
    root = tmp_path / "B56a_EOAI10089186_EOAI10403337"
    touch(root / "seed-1" / "model.cif")
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\nEOAI10089186\tCCO\nEOAI10403337\tCCC\n",
        encoding="utf-8",
    )
    jobs = ccp.build_structure_jobs(root, None, "LIG", ("A",), "F")

    with pytest.raises(ccp.PipelineError, match="multiple compound IDs"):
        ccp.apply_smiles_table(jobs, table)


def test_smiles_table_reports_missing_match_before_minimization(tmp_path):
    root = tmp_path / "B56a_UNKNOWN"
    touch(root / "seed-1" / "model.cif")
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\nEOAI10089186\tCCO\n",
        encoding="utf-8",
    )
    jobs = ccp.build_structure_jobs(root, None, "LIG", ("A",), "F")

    with pytest.raises(ccp.PipelineError, match="no Compound ID from"):
        ccp.apply_smiles_table(jobs, table)


def test_manifest_rejects_dataset_path_outside_root(tmp_path):
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\n../outside\tC19\tWXI\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="outside ROOT"):
        ccp.load_manifest(tmp_path, manifest)


def test_minimization_requires_exact_ligand_chemistry(tmp_path):
    root = tmp_path / "root"
    touch(root / "C19" / "seed-1" / "model.cif")
    jobs = ccp.build_structure_jobs(root, None, "WXI", None, "B", None, None)
    args = ccp.build_parser().parse_args(
        [str(root), "--ligand", "WXI", "--minimize", "openmm"]
    )

    with pytest.raises(ccp.PipelineError, match="ligand chemistry"):
        ccp.validate_minimization_configuration(args, jobs)


def test_minimization_rejects_global_chemistry_for_multiple_compounds(tmp_path):
    root = tmp_path / "root"
    touch(root / "C19" / "seed-1" / "model.cif")
    touch(root / "GSH" / "seed-1" / "model.cif")
    jobs = ccp.build_structure_jobs(root, None, "WXI", None, "B", "CC", None)
    args = ccp.build_parser().parse_args(
        [str(root), "--ligand", "WXI", "--ligand-smiles", "CC", "--minimize", "openmm"]
    )

    with pytest.raises(ccp.PipelineError, match="one discovered compound"):
        ccp.validate_minimization_configuration(args, jobs)


def test_cli_version_outputs_pipeline_version(capsys):
    with pytest.raises(SystemExit) as exc_info:
        ccp.build_parser().parse_args(["--version"])

    assert exc_info.value.code == 0
    assert capsys.readouterr().out == "2.2.3\n"


def test_manifest_rejects_smiles_and_ligand_file_together(tmp_path):
    root = tmp_path / "root"
    (root / "C19").mkdir(parents=True)
    sdf = root / "C19.sdf"
    sdf.write_text("fixture", encoding="utf-8")
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\tligand_smiles\tligand_file\n"
        "C19\tC19\tWXI\tCC\tC19.sdf\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="exactly one"):
        ccp.load_manifest(root, manifest)


def test_manifest_rejects_ligand_file_outside_root(tmp_path):
    root = tmp_path / "root"
    (root / "C19").mkdir(parents=True)
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\tligand_file\n"
        "C19\tC19\tWXI\t../outside.sdf\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="outside ROOT"):
        ccp.load_manifest(root, manifest)


def test_manifest_rejects_conflicting_chemistry_for_one_compound(tmp_path):
    root = tmp_path / "root"
    (root / "C19-a").mkdir(parents=True)
    (root / "C19-b").mkdir(parents=True)
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\tligand_smiles\n"
        "C19-a\tC19\tWXI\tCC\n"
        "C19-b\tC19\tWXI\tCCC\n",
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match="conflicting ligand or chain"):
        ccp.load_manifest(root, manifest)


def test_contact_only_manifest_chemistry_is_optional(tmp_path):
    root = tmp_path / "root"
    (root / "C19").mkdir(parents=True)
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\nC19\tC19\tWXI\n",
        encoding="utf-8",
    )

    spec = ccp.load_manifest(root, manifest)[0]

    assert spec.ligand_smiles is None
    assert spec.ligand_file is None


def test_example_manifest_has_versioned_chemistry_columns_and_parses(tmp_path):
    root = tmp_path / "root"
    (root / "screen" / "C19").mkdir(parents=True)
    sdf = root / "screen" / "GSH" / "GSH.sdf"
    sdf.parent.mkdir(parents=True)
    sdf.write_text("fixture", encoding="utf-8")
    manifest = Path(__file__).parents[1] / "datasets.example.tsv"

    specs = ccp.load_manifest(root, manifest)

    assert manifest.read_text(encoding="utf-8").splitlines()[0].split("\t") == [
        "dataset_dir",
        "compound",
        "ligand_resname",
        "protein_chains",
        "ligand_chain",
        "ligand_smiles",
        "ligand_file",
        "enabled",
    ]
    assert specs[0].ligand_smiles == "CC(C)C1=CC=C(C=C1)C(C)C(=O)O"
    assert specs[0].ligand_file is None
    assert specs[1].ligand_smiles is None
    assert specs[1].ligand_file == sdf.resolve()


def test_example_smiles_table_contains_supplied_compound_and_omits_invalid_entry():
    table = Path(__file__).parents[1] / "compound_smiles.example.tsv"

    entries = ccp.load_smiles_table(table)

    assert entries["EOAI10089186"] == (
        "CC1=CC=C(C)C(=C1)C1=NN2C=CN3C(=O)N(CC(=O)NCC4=CC=CC=C4)N=C3C2=C1"
    )
    assert "EOAI10404813" not in entries


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


def test_analyze_structure_separates_protein_and_ligand_chain_filters(tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("C", " ", 20, " ", "GLY", [("CA", "C", (20.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (4.5, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(
        cif, "C19/seed-1/model.cif", "C19", "C19", "WXI", ("A",), "B"
    )

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
    job = ccp.StructureJob(
        cif, "C19/seed-1/far.cif", "C19", "C19", "WXI", ("A",), "B"
    )

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
    job = ccp.StructureJob(
        cif, "mix/seed-1/model.cif", "mix", "mix", None, ("A",), None
    )

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
    job = ccp.StructureJob(
        cif, "C19/seed-1/models.cif", "C19", "C19", "WXI", ("A",), "B"
    )

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
def test_missing_required_structure_components_are_invalid(
    tmp_path, residue_specs, message
):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "missing.cif", residue_specs
    )
    job = ccp.StructureJob(
        cif, "C19/seed-1/missing.cif", "C19", "C19", "WXI", ("A",), "B"
    )

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "invalid"
    assert any(message in record.message for record in result.qc)


def test_malformed_cif_is_an_invalid_result(tmp_path):
    cif = touch(tmp_path / "C19" / "seed-1" / "bad.cif")
    cif.write_text("not an mmCIF", encoding="utf-8")
    job = ccp.StructureJob(
        cif, "C19/seed-1/bad.cif", "C19", "C19", "WXI", None, None
    )

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "invalid"
    assert result.qc[0].stage == "parse"


def analysis_job(tmp_path):
    cif_path = write_test_cif(
        tmp_path / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    return ccp.StructureJob(
        cif_path=cif_path,
        structure_id="C19/seed-1/model.cif",
        dataset_rel="C19",
        compound="C19",
        ligand_resname="WXI",
        protein_chains=("A",),
        ligand_chain="B",
    )


def test_analyze_structure_fast_mode_never_calls_sasa(monkeypatch, tmp_path):
    job = analysis_job(tmp_path)

    def forbidden_sasa(protein_residues, ligand_residues):
        raise AssertionError("SASA must not run in fast mode")

    monkeypatch.setattr(ccp, "calculate_buried_sasa", forbidden_sasa)
    result = ccp.analyze_structure(job, compute_sasa=False)

    assert result.status == "valid"
    assert result.metrics.buried_sasa_total_A2 is None
    assert result.metrics.interface_area_A2 is None
    assert all(record.stage != "sasa" for record in result.qc)


def test_analyze_structure_default_retains_full_sasa(monkeypatch, tmp_path):
    job = analysis_job(tmp_path)
    calls = []

    def fake_sasa(protein_residues, ligand_residues):
        calls.append((len(protein_residues), len(ligand_residues)))
        return 42.0, 21.0

    monkeypatch.setattr(ccp, "calculate_buried_sasa", fake_sasa)
    result = ccp.analyze_structure(job)

    assert calls == [(1, 1)]
    assert result.metrics.buried_sasa_total_A2 == pytest.approx(42.0)
    assert result.metrics.interface_area_A2 == pytest.approx(21.0)


def test_buried_sasa_is_positive_for_touching_receptor_and_ligand(tmp_path):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "touching.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(
        cif, "C19/seed-1/touching.cif", "C19", "C19", "WXI", ("A",), "B"
    )

    metrics = ccp.analyze_structure(job, 4.5, 2.0).metrics

    assert metrics.buried_sasa_total_A2 is not None
    assert metrics.buried_sasa_total_A2 > 0.0
    assert metrics.interface_area_A2 == pytest.approx(
        metrics.buried_sasa_total_A2 / 2.0
    )




def _single_atom_residue(chain_id, residue_name, element, coord):
    structure = Structure.Structure("sasa-test")
    model = Model.Model(0)
    structure.add(model)
    chain = Chain.Chain(chain_id)
    model.add(chain)
    residue = Residue.Residue((" ", 1, " "), residue_name, "")
    chain.add(residue)
    atom = Atom.Atom(
        element,
        np.asarray(coord, dtype=float),
        1.0,
        1.0,
        " ",
        f" {element:<3}"[:4],
        1,
        element=element,
    )
    residue.add(atom)
    return residue


def test_freesasa_backend_uses_safe_default_lee_richards_path(monkeypatch):
    calls = []
    init_values = []
    setter_calls = []

    class FakeResult:
        def __init__(self, area):
            self._area = area

        def totalArea(self):
            return self._area

    class FakeParameters:
        def __init__(self, values=None):
            init_values.append(values)

        def setProbeRadius(self, value):
            setter_calls.append(("probe-radius", value))

        def setNSlices(self, value):
            setter_calls.append(("n-slices", value))

        def setNThreads(self, value):
            setter_calls.append(("n-threads", value))

        def algorithm(self):
            return "LeeRichards"

    class FakeFreeSASA:
        LeeRichards = "LeeRichards"
        ShrakeRupley = "ShrakeRupley"
        __version__ = "test-freesasa"
        Parameters = FakeParameters

        @staticmethod
        def calcCoord(coords, radii, parameters):
            calls.append((list(coords), list(radii), parameters))
            return FakeResult([100.0, 50.0, 120.0][len(calls) - 1])

    monkeypatch.setattr(
        ccp, "calculate_buried_sasa", _PRODUCTION_CALCULATE_BURIED_SASA
    )
    monkeypatch.setattr(ccp, "_load_freesasa", lambda: FakeFreeSASA)
    protein = [_single_atom_residue("A", "ALA", "C", (0.0, 0.0, 0.0))]
    ligand = [_single_atom_residue("B", "LIG", "O", (3.0, 0.0, 0.0))]

    buried, interface = ccp.calculate_buried_sasa(protein, ligand)

    assert buried == pytest.approx(30.0)
    assert interface == pytest.approx(15.0)
    assert len(calls) == 3
    assert calls[0][1] == [1.7]
    assert calls[1][1] == [1.52]
    assert calls[2][1] == [1.7, 1.52]
    assert init_values == [None]
    assert setter_calls == [
        ("probe-radius", 1.4),
        ("n-slices", 20),
        ("n-threads", 1),
    ]


def test_real_freesasa_tracks_reference_when_installed(monkeypatch):
    freesasa = pytest.importorskip("freesasa")
    monkeypatch.setattr(ccp, "_load_freesasa", lambda: freesasa)
    monkeypatch.setattr(
        ccp, "calculate_buried_sasa", _PRODUCTION_CALCULATE_BURIED_SASA
    )
    protein = [_single_atom_residue("A", "ALA", "C", (0.0, 0.0, 0.0))]
    ligand = [_single_atom_residue("B", "LIG", "O", (3.0, 0.0, 0.0))]

    fast_buried, fast_interface = ccp.calculate_buried_sasa(protein, ligand)
    ref_buried, ref_interface = ccp.calculate_buried_sasa_biopython_reference(
        protein, ligand
    )

    assert fast_buried > 0.0
    assert fast_interface == pytest.approx(fast_buried / 2.0)
    assert fast_buried == pytest.approx(ref_buried, rel=0.20)
    assert ref_interface == pytest.approx(ref_buried / 2.0)


def test_freesasa_dependency_error_has_install_hint(monkeypatch):
    def missing_import(name):
        assert name == "freesasa"
        raise ModuleNotFoundError("No module named 'freesasa'")

    monkeypatch.setattr(ccp.importlib, "import_module", missing_import)

    with pytest.raises(ccp.PipelineError, match="freesasa-python"):
        _PRODUCTION_LOAD_FREESASA()


def test_pipeline_preflights_freesasa_before_discovering_structures(monkeypatch, tmp_path):
    args = ccp.build_parser().parse_args([str(tmp_path)])

    monkeypatch.setattr(
        ccp,
        "_load_freesasa",
        lambda: (_ for _ in ()).throw(ccp.PipelineError("FreeSASA unavailable")),
    )

    def forbidden_discovery(*args, **kwargs):
        raise AssertionError("structure discovery must not start before SASA dependency check")

    monkeypatch.setattr(ccp, "build_structure_jobs", forbidden_discovery)

    with pytest.raises(ccp.PipelineError, match="FreeSASA unavailable"):
        ccp.run_pipeline(args)

def test_sasa_failure_becomes_warning_without_invalidating_contacts(
    monkeypatch, tmp_path
):
    cif = write_test_cif(
        tmp_path / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 1, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    job = ccp.StructureJob(
        cif, "C19/seed-1/model.cif", "C19", "C19", "WXI", ("A",), "B"
    )

    def raise_sasa_error(protein_residues, ligand_residues):
        raise ValueError("bad SASA")

    monkeypatch.setattr(ccp, "calculate_buried_sasa", raise_sasa_error)

    result = ccp.analyze_structure(job, 4.5, 2.0)

    assert result.status == "valid"
    assert result.metrics.buried_sasa_total_A2 is None
    assert any(record.stage == "sasa" for record in result.qc)


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
    job1 = ccp.StructureJob(
        tmp_path / "one.cif", "C/seed-1/one.cif", "C", "C", "LIG", None, None
    )
    job2 = ccp.StructureJob(
        tmp_path / "two.cif", "C/seed-2/two.cif", "C", "C", "LIG", None, None
    )
    first = make_valid_result(
        job1, [make_contact(job1.structure_id, "C", residue, True, 3.0)]
    )
    second = make_valid_result(job2, [])

    row = ccp.aggregate_residues([first, second])[0]

    assert row["n_valid_structures"] == 2
    assert row["n_present"] == 1
    assert row["n_contacts"] == 1
    assert row["contact_frequency"] == 1.0


def test_contact_set_jaccard_has_explicit_empty_set_behavior():
    assert ccp.contact_set_jaccard(set(), set()) == 1.0
    assert ccp.contact_set_jaccard(set(), {"A:1"}) == 0.0
    assert ccp.contact_set_jaccard(
        {"A:1", "A:2"}, {"A:2", "A:3"}
    ) == pytest.approx(1 / 3)


def test_global_residue_summary_uses_present_structures_across_compounds(tmp_path):
    residue = ccp.ResidueKey("A", 10, "", "ALA")
    job1 = ccp.StructureJob(
        tmp_path / "one.cif", "C1/seed-1/one.cif", "C1", "C1", "LIG", None, None
    )
    job2 = ccp.StructureJob(
        tmp_path / "two.cif", "C2/seed-1/two.cif", "C2", "C2", "LIG", None, None
    )
    first = make_valid_result(
        job1, [make_contact(job1.structure_id, "C1", residue, True, 3.0)]
    )
    second = make_valid_result(
        job2, [make_contact(job2.structure_id, "C2", residue, False, 4.5)]
    )

    row = ccp.aggregate_global_residues([first, second])[0]

    assert row["n_present"] == 2
    assert row["n_contacts"] == 1
    assert row["contact_frequency"] == 0.5


def test_compound_summary_counts_invalid_structures_and_compares_contact_sets(
    tmp_path,
):
    residues = [
        ccp.ResidueKey("A", 1, "", "ALA"),
        ccp.ResidueKey("A", 2, "", "GLY"),
        ccp.ResidueKey("A", 3, "", "SER"),
    ]
    jobs = [
        ccp.StructureJob(
            tmp_path / f"{index}.cif",
            f"C/seed-{index}/{index}.cif",
            "C",
            "C",
            "LIG",
            None,
            None,
        )
        for index in range(1, 4)
    ]
    first = make_valid_result(
        jobs[0],
        [make_contact(jobs[0].structure_id, "C", residue, True, 3.0) for residue in residues[:2]],
    )
    second = make_valid_result(
        jobs[1],
        [make_contact(jobs[1].structure_id, "C", residue, True, 3.0) for residue in residues[1:]],
    )
    invalid = ccp.StructureResult(
        jobs[2],
        "invalid",
        (),
        ccp.StructureMetrics(None, 0, 0, None, 0),
        (ccp.QCRecord("error", "parse", "bad", "C", jobs[2].structure_id),),
    )

    row = ccp.aggregate_compounds([first, second, invalid])[0]

    assert row["n_discovered"] == 3
    assert row["n_valid"] == 2
    assert row["n_invalid"] == 1
    assert row["mean_pairwise_contact_jaccard"] == pytest.approx(1 / 3)


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
    assert ccp.bootstrap_mean_ci([1.0, 2.0, 3.0]) == ccp.bootstrap_mean_ci(
        [1.0, 2.0, 3.0]
    )


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
        "compound\tmetric\tvalue\tunit\nC19\tKi\t5\tnM\n", encoding="utf-8"
    )

    without_conversion = ccp.load_affinity_results(
        table, {"C19"}, convert_ki=False
    )[0]
    with_conversion = ccp.load_affinity_results(table, {"C19"}, convert_ki=True)[
        0
    ]

    assert without_conversion.delta_g_standard_kcal_mol is None
    assert with_conversion.delta_g_standard_kcal_mol is not None


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        (
            "C19\tmmgbsa\t-30.0\t1.0\t1\nC19\tmmgbsa\t-29.0\t-0.1\t2\n",
            "uncertainty",
        ),
        (
            "C19\tmmgbsa\t-30.0\t1.0\t1\nC19\tmmgbsa\t-30.0\t1.0\t1\n",
            "duplicate",
        ),
    ],
)
def test_energy_table_rejects_negative_uncertainty_and_duplicate_rows(
    tmp_path, rows, message
):
    table = tmp_path / "energy.tsv"
    table.write_text(
        "compound\tmethod\tvalue_kcal_mol\tuncertainty_kcal_mol\treplicate\n"
        + rows,
        encoding="utf-8",
    )

    with pytest.raises(ccp.PipelineError, match=message):
        ccp.load_energy_results(table, {"C19"}, set())


def test_write_tsv_formats_floats_and_blanks_missing_values(tmp_path):
    output = tmp_path / "values.tsv"

    ccp.write_tsv(
        output, [{"name": "A", "value": 1.23456789, "missing": None}]
    )

    assert (
        output.read_text(encoding="utf-8")
        == "name\tvalue\tmissing\nA\t1.234568\t\n"
    )


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
        json.dumps(
            {"files": ["kept.tsv", "../victim.txt", "generated_files.json"]}
        ),
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

    with pytest.raises(
        ccp.PipelineError, match="not listed in the preceding inventory"
    ):
        ccp.claim_output_path(unrelated, output)
    assert unrelated.read_text(encoding="utf-8") == "unrelated"


def test_overwrite_removes_only_files_from_valid_prior_inventory(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    generated = output / "old.tsv"
    generated.write_text("old", encoding="utf-8")
    unrelated = output / "notes.txt"
    unrelated.write_text("keep", encoding="utf-8")
    inventory = output / "generated_files.json"
    inventory.write_text(
        json.dumps({"files": ["generated_files.json", "old.tsv"]}),
        encoding="utf-8",
    )

    ccp.prepare_output_directory(output, overwrite=True)

    assert not generated.exists()
    assert not inventory.exists()
    assert unrelated.read_text(encoding="utf-8") == "keep"


def test_overwrite_unlinks_inventory_symlink_without_deleting_target(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    target = output / "target.tsv"
    target.write_text("keep", encoding="utf-8")
    link = output / "old.tsv"
    link.symlink_to(target.name)
    inventory = output / "generated_files.json"
    inventory.write_text(
        json.dumps({"files": ["generated_files.json", "old.tsv"]}),
        encoding="utf-8",
    )

    ccp.prepare_output_directory(output, overwrite=True)

    assert not link.exists()
    assert not link.is_symlink()
    assert target.read_text(encoding="utf-8") == "keep"


def test_generated_inventory_is_sorted_and_includes_itself(tmp_path):
    output = tmp_path / "results"
    output.mkdir()
    table = output / "z.tsv"
    plot = output / "plots" / "a.png"
    table.touch()
    plot.parent.mkdir()
    plot.touch()

    inventory = ccp.write_generated_inventory(output, [table, plot])

    payload = json.loads(inventory.read_text(encoding="utf-8"))
    assert payload == {
        "files": ["generated_files.json", "plots/a.png", "z.tsv"]
    }


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

    paths = ccp.plot_compound_chain_contacts(
        rows, tmp_path, cutoff=4.5, label_frequency=0.5
    )
    paths.extend(ccp.plot_chain_heatmaps(rows, tmp_path, cutoff=4.5))

    assert paths
    assert all(path.exists() and path.stat().st_size > 0 for path in paths)


def test_heatmap_matrix_uses_nan_for_absent_compound_residue_pairs():
    rows = [
        {
            "compound": "C1",
            "chain": "A",
            "resseq": 10,
            "icode": "",
            "resname": "ALA",
            "contact_frequency": 1.0,
        },
        {
            "compound": "C2",
            "chain": "A",
            "resseq": 11,
            "icode": "",
            "resname": "GLY",
            "contact_frequency": 0.5,
        },
    ]

    matrix, compounds, residues = ccp.build_heatmap_matrix(
        rows, "contact_frequency"
    )

    assert compounds == ["C1", "C2"]
    assert [residue.resseq for residue in residues] == [10, 11]
    assert np.isnan(matrix[0, 1])
    assert np.isnan(matrix[1, 0])


def test_evidence_summary_plot_functions_create_nonempty_pngs(tmp_path):
    compound_rows = [
        {
            "compound": "C19",
            "n_valid": 2,
            "mean_pairwise_contact_jaccard": 0.75,
            "mean_buried_sasa_total_A2": 220.0,
            "std_buried_sasa_total_A2": 10.0,
            "mean_clash_pair_count": 1.0,
            "std_clash_pair_count": 0.5,
        }
    ]
    energy_values = [
        ccp.EnergyValue("C19", "vina", -7.0, replicate="1"),
        ccp.EnergyValue("C19", "vina", -8.0, replicate="2"),
    ]
    affinity_values = [
        ccp.AffinityValue(
            "C19", "Kd", 10.0, "nM", 1e-8, 298.15, "1", -10.91
        )
    ]

    paths = [ccp.plot_structural_summary(compound_rows, tmp_path)]
    paths.extend(
        ccp.plot_energy_methods(
            energy_values,
            ccp.aggregate_energy_values(energy_values),
            tmp_path,
        )
    )
    paths.append(
        ccp.plot_experimental_affinity(
            affinity_values,
            ccp.aggregate_affinity_values(affinity_values),
            tmp_path,
        )
    )

    assert all(path.exists() and path.stat().st_size > 0 for path in paths)


def test_cli_rejects_fail_fast_with_multiple_workers(tmp_path):
    parser = ccp.build_parser()
    args = parser.parse_args(
        [str(tmp_path), "--fail-fast", "--workers", "2"]
    )

    with pytest.raises(ccp.PipelineError, match="--workers 1"):
        ccp.validate_args(args)


@pytest.mark.parametrize(
    ("option", "value"),
    [("--cutoff", "nan"), ("--clash-cutoff", "inf")],
)
def test_cli_rejects_nonfinite_distance_cutoffs(tmp_path, option, value):
    args = ccp.build_parser().parse_args([str(tmp_path), option, value])

    with pytest.raises(ccp.PipelineError, match="finite and positive"):
        ccp.validate_args(args)


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--minimize-ph", "nan"),
        ("--pocket-radius", "0"),
        ("--max-backbone-rmsd", "inf"),
        ("--max-ligand-rmsd", "-1"),
        ("--minimization-tolerance", "0"),
    ],
)
def test_cli_rejects_invalid_minimization_numeric_values(
    tmp_path, option, value
):
    args = ccp.build_parser().parse_args([str(tmp_path), option, value])

    with pytest.raises(ccp.PipelineError, match="finite and positive"):
        ccp.validate_args(args)


def test_cli_rejects_minimized_analysis_without_minimization(tmp_path):
    args = ccp.build_parser().parse_args(
        [str(tmp_path), "--analysis-source", "minimized"]
    )

    with pytest.raises(ccp.PipelineError, match="--minimize"):
        ccp.validate_args(args)


def test_cli_allows_minimization_with_multiple_workers(tmp_path):
    args = ccp.build_parser().parse_args(
        [str(tmp_path), "--minimize", "openmm", "--workers", "16"]
    )

    ccp.validate_args(args)
    assert args.workers == 16


def test_cli_quiet_defaults_false_and_can_be_enabled(tmp_path):
    parser = ccp.build_parser()
    assert parser.parse_args([str(tmp_path)]).quiet is False
    assert parser.parse_args([str(tmp_path), "--quiet"]).quiet is True


def test_progress_reporter_quiet_suppresses_info(capsys):
    reporter = ccp.ProgressReporter(quiet=True)
    reporter.info("hidden progress")
    assert capsys.readouterr().err == ""


def test_progress_reporter_failure_is_side_effect_free():
    class BrokenStream:
        def write(self, value):
            raise OSError("stream unavailable")

        def flush(self):
            raise OSError("stream unavailable")

    reporter = ccp.ProgressReporter(stream=BrokenStream())
    reporter.info("best effort only")


def test_quiet_does_not_suppress_cli_errors(tmp_path, capsys):
    code = ccp.main([str(tmp_path), "--quiet", "--cutoff", "0"])
    captured = capsys.readouterr()
    assert code == 2
    assert "ERROR:" in captured.err


def test_cli_rejects_device_index_on_cpu(tmp_path):
    args = ccp.build_parser().parse_args(
        [str(tmp_path), "--openmm-device-index", "0"]
    )

    with pytest.raises(ccp.PipelineError, match="GPU"):
        ccp.validate_args(args)


def test_cli_rejects_manifest_with_global_ligand_chemistry(tmp_path):
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text("", encoding="utf-8")
    args = ccp.build_parser().parse_args(
        [str(tmp_path), "--manifest", str(manifest), "--ligand-smiles", "CC"]
    )

    with pytest.raises(ccp.PipelineError, match="authoritative"):
        ccp.validate_args(args)


def test_cli_rejects_global_smiles_and_ligand_file_together(tmp_path):
    args = ccp.build_parser().parse_args(
        [
            str(tmp_path),
            "--ligand-smiles",
            "CC",
            "--ligand-file",
            "C19.sdf",
        ]
    )

    with pytest.raises(ccp.PipelineError, match="exactly one"):
        ccp.validate_args(args)


def minimization_fixture(tmp_path, *, ligand_x=1.5):
    root = tmp_path / "inputs"
    cif_path = write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            (
                "A",
                " ",
                10,
                " ",
                "ALA",
                [("CA", "C", (0.0, 0.0, 0.0))],
            ),
            (
                "B",
                "H_WXI",
                1,
                " ",
                "WXI",
                [("C1", "C", (ligand_x, 0.0, 0.0))],
            ),
        ],
    )
    args = ccp.build_parser().parse_args(
        [
            str(root),
            "--ligand",
            "WXI",
            "--protein-chain",
            "A",
            "--ligand-chain",
            "B",
            "--ligand-smiles",
            "CC",
            "--minimize",
            "openmm",
        ]
    )
    jobs = ccp.build_structure_jobs(
        root, None, "WXI", ("A",), "B", "CC", None
    )
    original_results = [ccp.analyze_structure(jobs[0])]
    return cif_path, args, jobs, original_results


def fake_minimization_backend(*, ligand_x=2.5, ligand_rmsd_A=0.5, fail=None):
    def minimize(request, config, cache_path, progress_callback=None):
        if fail is not None:
            raise fail
        protein_chain = (request.protein_chains or ("A",))[0]
        ligand_chain = request.ligand_chain or "B"
        write_test_cif(
            request.output_path,
            [
                (
                    protein_chain,
                    " ",
                    10,
                    " ",
                    "ALA",
                    [("CA", "C", (0.0, 0.0, 0.0))],
                ),
                (
                    ligand_chain,
                    f"H_{request.ligand_resname}",
                    1,
                    " ",
                    request.ligand_resname,
                    [("C1", "C", (ligand_x, 0.0, 0.0))],
                ),
            ],
        )
        return rom.BackendMinimizationResult(
            output_path=request.output_path,
            initial_potential_kcal_mol=12.0,
            final_potential_kcal_mol=4.0,
            backbone_rmsd_A=0.1,
            ligand_rmsd_A=ligand_rmsd_A,
            before_vdw_clash_pair_count=1,
            after_vdw_clash_pair_count=0,
            before_maximum_vdw_overlap_A=1.9,
            after_maximum_vdw_overlap_A=0.1,
            package_versions=(("openmm", "test"),),
        )

    return minimize


def test_minimization_candidate_qc_skips_sasa(monkeypatch, tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    calls = []
    real_analyze = ccp.analyze_structure

    def tracking_analyze(job, cutoff=4.5, clash_cutoff=2.0, compute_sasa=True):
        calls.append((job.cif_path, compute_sasa))
        return real_analyze(
            job,
            cutoff,
            clash_cutoff,
            compute_sasa=compute_sasa,
        )

    monkeypatch.setattr(ccp, "analyze_structure", tracking_analyze)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )

    assert records[0].status == "accepted"
    assert calls == [(records[0].staged_cif_path, False)]
    assert records[0].candidate_result.metrics.buried_sasa_total_A2 is None


def test_minimization_remains_sequential_with_multiple_cpu_workers(tmp_path):
    root = tmp_path / "inputs"
    for seed in ("seed-1", "seed-2"):
        write_test_cif(
            root / "C19" / seed / "model.cif",
            [
                ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
                ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (1.5, 0.0, 0.0))]),
            ],
        )
    args = ccp.build_parser().parse_args(
        [
            str(root),
            "--ligand", "WXI",
            "--protein-chain", "A",
            "--ligand-chain", "B",
            "--ligand-smiles", "CC",
            "--minimize", "openmm",
            "--workers", "4",
        ]
    )
    jobs = ccp.build_structure_jobs(root, None, "WXI", ("A",), "B", "CC", None)
    original_results = [ccp.analyze_structure(job) for job in jobs]
    entered = []
    active = False
    backend = fake_minimization_backend()

    def tracking_backend(request, config, cache_path, progress_callback=None):
        nonlocal active
        assert active is False
        active = True
        entered.append(request.structure_id)
        try:
            return backend(
                request,
                config,
                cache_path,
                progress_callback=progress_callback,
            )
        finally:
            active = False

    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=tracking_backend,
    )

    assert [record.status for record in records] == ["accepted", "accepted"]
    assert entered == [job.structure_id for job in jobs]


def _tracking_analysis_modes(monkeypatch):
    real_analyze = ccp.analyze_structure
    modes = []

    def tracking_analyze(job, cutoff=4.5, clash_cutoff=2.0, compute_sasa=True):
        modes.append(compute_sasa)
        return real_analyze(
            job,
            cutoff,
            clash_cutoff,
            compute_sasa=compute_sasa,
        )

    monkeypatch.setattr(ccp, "analyze_structure", tracking_analyze)
    return modes


def test_minimized_analysis_runs_fast_original_fast_candidate_then_full_accepted(
    monkeypatch, tmp_path
):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"
    modes = _tracking_analysis_modes(monkeypatch)

    assert ccp.run_pipeline(args, minimize_one=fake_minimization_backend()) == 0
    assert modes == [False, False, True]


def test_minimized_analysis_rejected_candidate_variant_skips_final_sasa(
    monkeypatch, tmp_path
):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"
    modes = _tracking_analysis_modes(monkeypatch)

    with pytest.raises(ccp.PipelineError, match="no valid minimized structures"):
        ccp.run_pipeline(
            args,
            minimize_one=fake_minimization_backend(ligand_rmsd_A=1.6),
        )
    assert modes == [False, False]


def test_analysis_source_original_sasa_runs_full_original_fast_candidate(
    monkeypatch, tmp_path
):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "original"
    args.out_dir = tmp_path / "results"
    modes = _tracking_analysis_modes(monkeypatch)

    assert ccp.run_pipeline(args, minimize_one=fake_minimization_backend()) == 0
    assert modes == [True, False]


def test_no_minimization_full_sasa_only(monkeypatch, tmp_path):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    args = ccp.build_parser().parse_args(
        [
            str(root),
            "--ligand", "WXI",
            "--protein-chain", "A",
            "--ligand-chain", "B",
            "--out-dir", str(tmp_path / "results"),
        ]
    )
    modes = _tracking_analysis_modes(monkeypatch)

    assert ccp.run_pipeline(args) == 0
    assert modes == [True]


def test_pipeline_uses_smiles_table_for_minimization_and_reports_provenance(
    tmp_path,
):
    root = tmp_path / "B56a_EOAI10089186"
    write_test_cif(
        root / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("F", "H_LIG", 1, " ", "LIG", [("C1", "C", (1.5, 0.0, 0.0))]),
        ],
    )
    table = tmp_path / "compounds.tsv"
    table.write_text(
        "Compound ID\tSmiles\nEOAI10089186\tCC\n", encoding="utf-8"
    )
    output = tmp_path / "results"
    args = ccp.build_parser().parse_args(
        [
            str(root),
            "--ligand",
            "LIG",
            "--protein-chain",
            "A",
            "--ligand-chain",
            "F",
            "--smiles-table",
            str(table),
            "--minimize",
            "openmm",
            "--analysis-source",
            "minimized",
            "--out-dir",
            str(output),
        ]
    )

    assert ccp.run_pipeline(
        args, minimize_one=fake_minimization_backend()
    ) == 0
    with (output / "run_manifest.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["ligand_compound_id"] == "EOAI10089186"
    assert row["ligand_chemistry_source"] == (
        f"smiles_table:{table.resolve().as_posix()}"
    )


def test_accepted_minimization_selects_candidate_and_preserves_input(tmp_path):
    cif_path, args, jobs, original_results = minimization_fixture(tmp_path)
    original_bytes = cif_path.read_bytes()

    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )
    args.analysis_source = "minimized"
    full_minimized = [ccp.analyze_structure(records[0].candidate_result.job)]
    selected = ccp.select_analysis_results(
        jobs,
        original_results,
        records,
        args,
        full_minimized_results=full_minimized,
    )

    assert records[0].status == "accepted"
    assert selected[0].status == "valid"
    assert selected[0].job.structure_id == jobs[0].structure_id
    assert selected[0].metrics.clash_pair_count == 0
    assert cif_path.read_bytes() == original_bytes


def test_rejected_minimization_is_excluded_from_minimized_analysis(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)

    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(ligand_rmsd_A=1.6),
    )
    args.analysis_source = "minimized"
    selected = ccp.select_analysis_results(
        jobs,
        original_results,
        records,
        args,
        full_minimized_results=[None],
    )

    assert records[0].status == "rejected"
    assert selected[0].status == "invalid"
    assert selected[0].qc[-1].stage == "minimization_acceptance"
    assert "ligand RMSD" in selected[0].qc[-1].message


def test_final_full_analysis_failure_never_falls_back_to_fast_candidate(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )
    args.analysis_source = "minimized"
    fast_candidate = records[0].candidate_result
    final_invalid = ccp._invalid_structure_result(
        fast_candidate.job,
        "sasa",
        "final scientific analysis failed",
    )

    selected = ccp.select_analysis_results(
        jobs,
        original_results,
        records,
        args,
        full_minimized_results=[final_invalid],
    )

    assert selected[0].status == "invalid"
    assert selected[0].qc[-1].stage == "sasa"
    assert "final scientific analysis failed" in selected[0].qc[-1].message
    assert records[0].candidate_result is fast_candidate
    assert fast_candidate.status == "valid"
    assert fast_candidate.metrics.buried_sasa_total_A2 is None


def test_failed_minimization_becomes_qc_record(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)

    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(
            fail=rom.MinimizationBackendError("parameterization failed")
        ),
    )

    assert records[0].status == "failed"
    assert records[0].candidate_result is None
    assert records[0].reasons == ("parameterization failed",)


def test_invalid_original_structure_is_not_minimized(tmp_path):
    _, args, jobs, _ = minimization_fixture(tmp_path)
    invalid = ccp._invalid_structure_result(jobs[0], "parse", "bad CIF")
    called = False

    def unexpected_backend(request, config, cache_path):
        nonlocal called
        called = True
        raise AssertionError("backend must not be called")

    records = ccp.run_minimization_jobs(
        jobs, [invalid], args, tmp_path / "staging", unexpected_backend
    )

    assert not called
    assert records[0].status == "not_attempted"
    assert records[0].reasons == ("original structure is invalid",)


def test_original_analysis_ignores_minimization_outcome(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(ligand_rmsd_A=1.6),
    )

    selected = ccp.select_analysis_results(
        jobs, original_results, records, args
    )

    assert selected is original_results


def test_minimization_summary_has_stable_schema_and_blank_failed_metrics(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(
            fail=rom.MinimizationBackendError("parameterization failed")
        ),
    )

    rows = ccp.minimization_records_to_rows(records, original_results)
    expected_columns = [
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
    ]

    assert list(rows[0]) == expected_columns
    assert rows[0]["reason"] == "parameterization failed"
    assert rows[0]["initial_potential_kcal_mol"] is None
    summary = tmp_path / "minimization_summary.tsv"
    ccp.write_tsv(summary, rows, expected_columns)
    fields = summary.read_text(encoding="utf-8").splitlines()[1].split("\t")
    assert fields[13] == ""
    assert fields[16] == ""


def test_minimization_summary_formats_completed_metrics_and_versions(monkeypatch, tmp_path):
    monkeypatch.setattr(
        ccp, "_load_freesasa", lambda: SimpleNamespace(__version__="2.2.1")
    )
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )

    row = ccp.minimization_records_to_rows(records, original_results)[0]

    assert row["delta_potential_kcal_mol"] == -8.0
    assert row["before_fixed_clash_pair_count"] == 1
    assert row["after_fixed_clash_pair_count"] == 0
    assert row["output_cif"] == (
        "minimized_structures/accepted/C19/seed-1/model_minimized.cif"
    )
    assert row["package_versions"] == "freesasa=2.2.1;openmm=test"


def test_source_aware_rows_record_original_and_minimized_paths(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )
    args.analysis_source = "minimized"
    full_minimized = [ccp.analyze_structure(records[0].candidate_result.job)]
    selected = ccp.select_analysis_results(
        jobs,
        original_results,
        records,
        args,
        full_minimized_results=full_minimized,
    )

    manifest_row = ccp.jobs_to_rows(
        jobs, records, analysis_source="minimized"
    )[0]
    structure_row = ccp.structure_results_to_rows(
        selected, analysis_source="minimized"
    )[0]

    assert manifest_row["analysis_source"] == "minimized"
    assert manifest_row["original_cif_path"] == jobs[0].cif_path.as_posix()
    assert manifest_row["analysis_cif_path"] == records[0].output_relative_path.as_posix()
    assert structure_row["analysis_source"] == "minimized"


def test_minimization_before_after_plot_is_nonempty(tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )

    output = ccp.plot_minimization_before_after(
        records, original_results, tmp_path
    )

    assert output.exists()
    assert output.stat().st_size > 0


def test_run_pipeline_publishes_minimization_outputs_and_inventory(tmp_path):
    cif_path, args, _, _ = minimization_fixture(tmp_path)
    original_bytes = cif_path.read_bytes()
    args.out_dir = tmp_path / "results"

    status = ccp.run_pipeline(args, minimize_one=fake_minimization_backend())

    assert status == 0
    output_cif = (
        args.out_dir
        / "minimized_structures/accepted/C19/seed-1/model_minimized.cif"
    )
    assert output_cif.is_file()
    assert (args.out_dir / "minimization_summary.tsv").is_file()
    assert (args.out_dir / "plots/minimization_before_after.png").is_file()
    inventory = json.loads(
        (args.out_dir / "generated_files.json").read_text(encoding="utf-8")
    )
    assert output_cif.relative_to(args.out_dir).as_posix() in inventory["files"]
    assert cif_path.read_bytes() == original_bytes


def test_original_analysis_returns_one_but_keeps_tables_for_rejected_candidate(tmp_path):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.out_dir = tmp_path / "results"

    status = ccp.run_pipeline(
        args,
        minimize_one=fake_minimization_backend(ligand_rmsd_A=1.6),
    )

    assert status == 1
    assert (args.out_dir / "compound_summary.tsv").is_file()
    assert (
        args.out_dir
        / "minimized_structures/rejected/C19/seed-1/model_minimized.cif"
    ).is_file()
    assert "minimization_acceptance" in (
        args.out_dir / "qc.tsv"
    ).read_text(encoding="utf-8")


def test_minimized_analysis_with_no_accepted_structure_fails_cleanly(tmp_path):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"

    with pytest.raises(ccp.PipelineError, match="no valid minimized structures"):
        ccp.run_pipeline(
            args,
            minimize_one=fake_minimization_backend(ligand_rmsd_A=1.6),
        )


def test_minimized_analysis_no_accepted_publishes_rejection_summary_and_reason(
    tmp_path, capsys
):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"

    with pytest.raises(ccp.PipelineError, match="no valid minimized structures"):
        ccp.run_pipeline(
            args,
            minimize_one=fake_minimization_backend(ligand_rmsd_A=1.6),
        )

    captured = capsys.readouterr()
    assert "REJECTED" in captured.err
    assert "ligand RMSD" in captured.err
    summary = args.out_dir / "minimization_summary.tsv"
    assert summary.is_file()
    summary_text = summary.read_text(encoding="utf-8")
    assert "rejected" in summary_text
    assert "ligand RMSD" in summary_text


def test_minimized_analysis_no_accepted_publishes_failure_summary_and_reason(
    tmp_path, capsys
):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"

    with pytest.raises(ccp.PipelineError, match="no valid minimized structures"):
        ccp.run_pipeline(
            args,
            minimize_one=fake_minimization_backend(
                fail=rom.MinimizationBackendError("parameterization failed")
            ),
        )

    captured = capsys.readouterr()
    assert "FAILED" in captured.err
    assert "parameterization failed" in captured.err
    summary = args.out_dir / "minimization_summary.tsv"
    assert summary.is_file()
    summary_text = summary.read_text(encoding="utf-8")
    assert "failed" in summary_text
    assert "parameterization failed" in summary_text


def test_minimized_pipeline_reports_phase_and_candidate_progress(tmp_path, capsys):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"

    code = ccp.run_pipeline(args, minimize_one=fake_minimization_backend())
    captured = capsys.readouterr()

    assert code == 0
    assert "Phase 1/3: Pre-minimization analysis" in captured.err
    assert "Phase 2/3: Restrained OpenMM minimization" in captured.err
    assert "candidate geometry QC completed" in captured.err
    assert "Phase 3/3: Final minimized scientific analysis" in captured.err
    assert "Batch: 1/1" in captured.err
    assert "generated_files.json" in captured.out


def test_quiet_minimized_pipeline_suppresses_normal_progress(tmp_path, capsys):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.quiet = True
    args.out_dir = tmp_path / "results"

    code = ccp.run_pipeline(args, minimize_one=fake_minimization_backend())
    captured = capsys.readouterr()

    assert code == 0
    assert "Phase 1/3" not in captured.err
    assert "candidate geometry QC" not in captured.err
    assert "Batch: 1/1" not in captured.err
    assert "generated_files.json" in captured.out


def test_cli_end_to_end_writes_core_tables_and_contact_plots(tmp_path):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            (
                "A",
                " ",
                10,
                " ",
                "ALA",
                [("CA", "C", (0.0, 0.0, 0.0))],
            ),
            (
                "B",
                "H_WXI",
                1,
                " ",
                "WXI",
                [("C1", "C", (3.0, 0.0, 0.0))],
            ),
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


def test_cli_end_to_end_integrates_energy_and_affinity_tables(tmp_path):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            (
                "A",
                " ",
                10,
                " ",
                "ALA",
                [("CA", "C", (0.0, 0.0, 0.0))],
            ),
            (
                "B",
                "H_WXI",
                1,
                " ",
                "WXI",
                [("C1", "C", (3.0, 0.0, 0.0))],
            ),
        ],
    )
    energy = tmp_path / "energy.tsv"
    energy.write_text(
        "compound\tmethod\tvalue_kcal_mol\treplicate\n"
        "C19\tmmgbsa\t-30.0\t1\n"
        "C19\tmmgbsa\t-31.0\t2\n",
        encoding="utf-8",
    )
    affinity = tmp_path / "affinity.tsv"
    affinity.write_text(
        "compound\tmetric\tvalue\tunit\ttemperature_K\treplicate\n"
        "C19\tKd\t25\tnM\t298.15\t1\n"
        "C19\tKd\t30\tnM\t298.15\t2\n",
        encoding="utf-8",
    )
    output = tmp_path / "results"
    args = ccp.build_parser().parse_args(
        [
            str(root),
            "--ligand",
            "WXI",
            "--energy-results",
            str(energy),
            "--affinity-results",
            str(affinity),
            "--out-dir",
            str(output),
        ]
    )

    assert ccp.run_pipeline(args) == 0
    for name in [
        "energy_values_normalized.tsv",
        "energy_summary.tsv",
        "affinity_values_normalized.tsv",
        "affinity_summary.tsv",
    ]:
        assert (output / name).exists()
    assert list((output / "plots").glob("energy_*.png"))
    assert (output / "plots" / "experimental_affinity.png").exists()


def test_cli_parallel_analysis_preserves_sorted_structure_order(tmp_path):
    root = tmp_path / "inputs"
    for seed in ["seed-2", "seed-1"]:
        write_test_cif(
            root / "C19" / seed / "model.cif",
            [
                (
                    "A",
                    " ",
                    10,
                    " ",
                    "ALA",
                    [("CA", "C", (0.0, 0.0, 0.0))],
                ),
                (
                    "B",
                    "H_WXI",
                    1,
                    " ",
                    "WXI",
                    [("C1", "C", (3.0, 0.0, 0.0))],
                ),
            ],
        )
    output = tmp_path / "results"
    args = ccp.build_parser().parse_args(
        [
            str(root),
            "--ligand",
            "WXI",
            "--workers",
            "2",
            "--out-dir",
            str(output),
        ]
    )

    assert ccp.run_pipeline(args) == 0
    rows = output.joinpath("structure_summary.tsv").read_text(
        encoding="utf-8"
    ).splitlines()
    structure_ids = [row.split("\t", 1)[0] for row in rows[1:]]
    assert structure_ids == sorted(structure_ids)


def test_parallel_fast_analysis_preserves_input_order_and_skips_sasa(
    tmp_path, capsys
):
    root = tmp_path / "inputs"
    for seed, x in [("seed-2", 3.2), ("seed-1", 3.0)]:
        write_test_cif(
            root / "C19" / seed / "model.cif",
            [
                ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
                ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (x, 0.0, 0.0))]),
            ],
        )
    jobs = ccp.build_structure_jobs(root, None, "WXI", ("A",), "B")
    args = ccp.build_parser().parse_args([str(root), "--workers", "2"])
    reporter = ccp.ProgressReporter()

    results = ccp._analyze_jobs(
        jobs,
        args,
        compute_sasa=False,
        reporter=reporter,
        phase_label="Pre-minimization analysis",
        phase_number=1,
        phase_total=3,
    )

    assert [result.job.structure_id for result in results] == [
        job.structure_id for job in jobs
    ]
    assert all(result.metrics.buried_sasa_total_A2 is None for result in results)
    stderr = capsys.readouterr().err
    assert "Phase 1/3: Pre-minimization analysis" in stderr
    assert "2 structures | 2 CPU workers | SASA disabled" in stderr
    assert "Phase completed:" in stderr


def test_cli_rejects_conflicting_auto_detected_ligands_within_compound(
    tmp_path,
):
    root = tmp_path / "inputs"
    for seed, ligand in [("seed-1", "L01"), ("seed-2", "L02")]:
        write_test_cif(
            root / "mix" / seed / "model.cif",
            [
                (
                    "A",
                    " ",
                    1,
                    " ",
                    "ALA",
                    [("CA", "C", (0.0, 0.0, 0.0))],
                ),
                (
                    "B",
                    f"H_{ligand}",
                    1,
                    " ",
                    ligand,
                    [("C1", "C", (3.0, 0.0, 0.0))],
                ),
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
            (
                "A",
                " ",
                1,
                " ",
                "ALA",
                [("CA", "C", (0.0, 0.0, 0.0))],
            ),
            (
                "B",
                "H_WXI",
                1,
                " ",
                "WXI",
                [("C1", "C", (3.0, 0.0, 0.0))],
            ),
        ],
    )
    output = tmp_path / "results"
    args = ccp.build_parser().parse_args(
        [str(root), "--ligand", "WXI", "--out-dir", str(output)]
    )

    def fail_plot(*args, **kwargs):
        raise RuntimeError("plot failed")

    monkeypatch.setattr(ccp, "plot_compound_chain_contacts", fail_plot)

    assert ccp.run_pipeline(args) == 1
    assert (output / "compound_summary.tsv").exists()
    assert "plot failed" in (output / "qc.tsv").read_text(encoding="utf-8")
    assert not (output / "generated_files.json").exists()


def test_cli_continues_past_malformed_cif_when_one_structure_is_valid(
    tmp_path,
):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "good.cif",
        [
            (
                "A",
                " ",
                1,
                " ",
                "ALA",
                [("CA", "C", (0.0, 0.0, 0.0))],
            ),
            (
                "B",
                "H_WXI",
                1,
                " ",
                "WXI",
                [("C1", "C", (3.0, 0.0, 0.0))],
            ),
        ],
    )
    bad = touch(root / "C19" / "seed-2" / "bad.cif")
    bad.write_text("broken", encoding="utf-8")
    output = tmp_path / "results"

    completed = subprocess.run(
        [
            sys.executable,
            str(Path(ccp.__file__).resolve()),
            str(root),
            "--ligand",
            "WXI",
            "--out-dir",
            str(output),
        ],
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
        [
            sys.executable,
            str(Path(ccp.__file__).resolve()),
            str(root),
            "--ligand",
            "WXI",
        ],
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "no valid structures" in completed.stderr
