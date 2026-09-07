from pathlib import Path

import numpy as np
from Bio.PDB import Atom, Chain, MMCIFIO, Model, Residue, Structure

import benchmark_sasa as bs


def test_build_benchmark_result_reports_speedup_and_percent_difference():
    result = bs.build_benchmark_result(
        fast_buried=90.0,
        fast_interface=45.0,
        fast_seconds=1.0,
        reference_buried=100.0,
        reference_interface=50.0,
        reference_seconds=10.0,
    )

    assert result["speedup_vs_biopython"] == 10.0
    assert result["buried_sasa_percent_difference"] == -10.0
    assert result["freesasa_lee_richards_buried_sasa_A2"] == 90.0
    assert result["biopython_shrake_rupley_buried_sasa_A2"] == 100.0


def test_load_residues_selects_requested_protein_and_ligand_chains(tmp_path: Path):
    structure = Structure.Structure("bench")
    model = Model.Model(0)
    structure.add(model)

    protein_chain = Chain.Chain("A")
    model.add(protein_chain)
    ala = Residue.Residue((" ", 1, " "), "ALA", "")
    protein_chain.add(ala)
    ala.add(
        Atom.Atom(
            "CA",
            np.array([0.0, 0.0, 0.0]),
            1.0,
            1.0,
            " ",
            " CA ",
            1,
            element="C",
        )
    )

    ligand_chain = Chain.Chain("F")
    model.add(ligand_chain)
    lig = Residue.Residue(("H_LIG_F", 1, " "), "LIG_F", "")
    ligand_chain.add(lig)
    lig.add(
        Atom.Atom(
            "O1",
            np.array([3.0, 0.0, 0.0]),
            1.0,
            1.0,
            " ",
            " O1 ",
            2,
            element="O",
        )
    )

    cif_path = tmp_path / "model.cif"
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(str(cif_path))

    protein, ligand = bs.load_residues(
        cif_path,
        ligand_resname="LIG_F",
        protein_chains={"A"},
        ligand_chain="F",
    )

    assert len(protein) == 1
    assert protein[0].get_resname() == "ALA"
    assert len(ligand) == 1
    assert ligand[0].get_resname() == "LIG_F"


def test_benchmark_residues_uses_requested_repeats_and_reports_speedup():
    ticks = iter([0.0, 1.0, 1.0, 2.0, 2.0, 12.0, 12.0, 22.0])
    calls = {"fast": 0, "reference": 0}

    def timer():
        return next(ticks)

    def fast_fn(protein, ligand):
        calls["fast"] += 1
        return 90.0, 45.0

    def reference_fn(protein, ligand):
        calls["reference"] += 1
        return 100.0, 50.0

    result = bs.benchmark_residues(
        ["protein"],
        ["ligand"],
        repeats=2,
        fast_fn=fast_fn,
        reference_fn=reference_fn,
        timer=timer,
    )

    assert calls == {"fast": 2, "reference": 2}
    assert result["freesasa_lee_richards_seconds"] == 1.0
    assert result["biopython_shrake_rupley_seconds"] == 10.0
    assert result["speedup_vs_biopython"] == 10.0


def test_main_prints_machine_readable_summary(monkeypatch, tmp_path, capsys):
    cif_path = tmp_path / "model.cif"
    cif_path.write_text("placeholder", encoding="utf-8")

    monkeypatch.setattr(bs, "load_residues", lambda *args, **kwargs: (["p"], ["l"]))
    monkeypatch.setattr(
        bs,
        "benchmark_residues",
        lambda *args, **kwargs: {
            "freesasa_lee_richards_buried_sasa_A2": 90.0,
            "freesasa_lee_richards_interface_area_A2": 45.0,
            "freesasa_lee_richards_seconds": 1.0,
            "biopython_shrake_rupley_buried_sasa_A2": 100.0,
            "biopython_shrake_rupley_interface_area_A2": 50.0,
            "biopython_shrake_rupley_seconds": 10.0,
            "buried_sasa_percent_difference": -10.0,
            "speedup_vs_biopython": 10.0,
        },
    )

    exit_code = bs.main(
        [
            str(cif_path),
            "--ligand",
            "LIG_F",
            "--protein-chain",
            "A,B,C",
            "--ligand-chain",
            "F",
            "--repeats",
            "2",
        ]
    )

    assert exit_code == 0
    stdout = capsys.readouterr().out
    assert "metric\tvalue" in stdout
    assert "speedup_vs_biopython\t10" in stdout
    assert "buried_sasa_percent_difference\t-10" in stdout
