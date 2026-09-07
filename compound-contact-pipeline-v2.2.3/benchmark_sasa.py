#!/usr/bin/env python3
"""Benchmark v2.2.3 FreeSASA Lee-Richards against the Biopython reference."""

import argparse
from pathlib import Path
from statistics import median
from time import perf_counter

from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

import compound_contact_pipeline as ccp


def build_benchmark_result(
    *,
    fast_buried: float,
    fast_interface: float,
    fast_seconds: float,
    reference_buried: float,
    reference_interface: float,
    reference_seconds: float,
) -> dict[str, float]:
    percent_difference = 0.0
    if reference_buried != 0.0:
        percent_difference = (
            (fast_buried - reference_buried) / reference_buried * 100.0
        )
    speedup = (
        float("inf") if fast_seconds == 0.0 else reference_seconds / fast_seconds
    )
    return {
        "freesasa_lee_richards_buried_sasa_A2": fast_buried,
        "freesasa_lee_richards_interface_area_A2": fast_interface,
        "freesasa_lee_richards_seconds": fast_seconds,
        "biopython_shrake_rupley_buried_sasa_A2": reference_buried,
        "biopython_shrake_rupley_interface_area_A2": reference_interface,
        "biopython_shrake_rupley_seconds": reference_seconds,
        "buried_sasa_percent_difference": percent_difference,
        "speedup_vs_biopython": speedup,
    }


def load_residues(
    cif_path: Path,
    *,
    ligand_resname: str,
    protein_chains: set[str] | None,
    ligand_chain: str | None,
):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))
    models = list(structure)
    if not models:
        raise ValueError(f"no models found in {cif_path}")
    model = models[0]

    ligand_residues = [
        residue
        for chain in model
        if ligand_chain is None or str(chain.id) == ligand_chain
        for residue in chain
        if residue.get_resname().strip() == ligand_resname
        and not residue.get_id()[0].startswith("W")
    ]
    ligand_ids = {id(residue) for residue in ligand_residues}
    protein_residues = [
        residue
        for chain in model
        if protein_chains is None or str(chain.id) in protein_chains
        for residue in chain
        if id(residue) not in ligand_ids
        and not residue.get_id()[0].startswith("W")
        and is_aa(residue, standard=False)
    ]
    if not protein_residues:
        raise ValueError("no protein residues selected")
    if not ligand_residues:
        raise ValueError(f"ligand {ligand_resname} not found")
    return protein_residues, ligand_residues


def benchmark_residues(
    protein_residues,
    ligand_residues,
    *,
    repeats: int = 1,
    fast_fn=None,
    reference_fn=None,
    timer=perf_counter,
) -> dict[str, float]:
    if repeats < 1:
        raise ValueError("repeats must be >= 1")
    fast_fn = fast_fn or ccp.calculate_buried_sasa
    reference_fn = reference_fn or ccp.calculate_buried_sasa_biopython_reference

    fast_times = []
    reference_times = []
    fast_result = None
    reference_result = None

    for _ in range(repeats):
        start = timer()
        fast_result = fast_fn(protein_residues, ligand_residues)
        fast_times.append(timer() - start)

    for _ in range(repeats):
        start = timer()
        reference_result = reference_fn(protein_residues, ligand_residues)
        reference_times.append(timer() - start)

    fast_buried, fast_interface = fast_result
    reference_buried, reference_interface = reference_result
    return build_benchmark_result(
        fast_buried=fast_buried,
        fast_interface=fast_interface,
        fast_seconds=median(fast_times),
        reference_buried=reference_buried,
        reference_interface=reference_interface,
        reference_seconds=median(reference_times),
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compare v2.2.3 FreeSASA Lee-Richards buried SASA against the "
            "Biopython Shrake-Rupley reference on one mmCIF structure."
        )
    )
    parser.add_argument("cif", type=Path)
    parser.add_argument("--ligand", required=True, dest="ligand_resname")
    parser.add_argument("--protein-chain", default=None)
    parser.add_argument("--ligand-chain", default=None)
    parser.add_argument("--repeats", type=int, default=1)
    return parser


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    protein_chains = None
    if args.protein_chain:
        protein_chains = {
            part.strip() for part in args.protein_chain.split(",") if part.strip()
        }
    protein_residues, ligand_residues = load_residues(
        args.cif,
        ligand_resname=args.ligand_resname,
        protein_chains=protein_chains,
        ligand_chain=args.ligand_chain,
    )
    result = benchmark_residues(
        protein_residues,
        ligand_residues,
        repeats=args.repeats,
    )
    print("metric\tvalue")
    for key, value in result.items():
        print(f"{key}\t{value:.6g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
