#!/usr/bin/env python3
#python contacts_cif.py path/to/structure.cif --ligand {ligand_ID} --cutoff 4.5
import argparse
import csv
import math
from collections import defaultdict

from Bio.PDB import MMCIFParser, NeighborSearch, Selection
from Bio.PDB.Polypeptide import is_aa


def dist(a, b):
    ax, ay, az = a.get_coord()
    bx, by, bz = b.get_coord()
    dx, dy, dz = ax - bx, ay - by, az - bz
    return math.sqrt(dx*dx + dy*dy + dz*dz)


def residue_key(residue):
    """Stable key: (resname, chain_id, resseq, icode)"""
    chain = residue.get_parent()
    resname = residue.get_resname()
    resseq, icode = residue.get_id()[1], residue.get_id()[2]
    return (resname, chain.id, resseq, icode or "")


def format_residue_key(key):
    resname, chain, resseq, icode = key
    icode_str = icode if icode and icode != " " else ""
    return f"{resname} Chain {chain} {resseq}{icode_str}"


def find_contacts_cif(
    cif_path,
    ligand_name="GSH",
    cutoff=4.5,
    include_waters=False,
    standard_aa_only=True,
):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struc", cif_path)

    # Collect ligand and protein atoms
    ligand_residues = []
    protein_atoms = []

    for model in structure:
        for chain in model:
            for residue in chain:
                resname = residue.get_resname().strip()
                hetflag = residue.get_id()[0]  # ' ' for standard, 'H_' for het, 'W' for water in many files

                # ligand residues
                if resname == ligand_name:
                    ligand_residues.append(residue)
                    continue

                # skip waters unless requested
                if hetflag.startswith("W"):
                    if include_waters:
                        # Usually not needed for "protein_atoms", but kept for completeness
                        pass
                    continue

                # protein residues (amino acids)
                if is_aa(residue, standard=standard_aa_only):
                    protein_atoms.extend(residue.get_atoms())

    if not ligand_residues:
        print(f"❌ Ligand {ligand_name} not found in {cif_path}")
        return []

    if not protein_atoms:
        print("❌ No protein atoms found (polymer not detected).")
        return []

    ns = NeighborSearch(list(protein_atoms))

    all_reports = []

    for lig in ligand_residues:
        # Identify this ligand instance (chain, resseq, icode)
        chain = lig.get_parent()
        lig_resseq, lig_icode = lig.get_id()[1], lig.get_id()[2]
        lig_label = f"{lig.get_resname()} Chain {chain.id} {lig_resseq}{'' if (not lig_icode or lig_icode==' ') else lig_icode}"

        # Map contacting residue -> minimal distance
        contacts_min = defaultdict(lambda: float("inf"))

        lig_atoms = list(lig.get_atoms())
        for latom in lig_atoms:
            # candidate neighboring protein atoms
            neighbors = ns.search(latom.get_coord(), cutoff)
            for patom in neighbors:
                pres = patom.get_parent()
                # only amino-acid residues (avoid accidentally picking hets)
                if not is_aa(pres, standard=False):
                    continue
                # compute exact distance and keep minimum per residue
                d = dist(latom, patom)
                key = residue_key(pres)
                if d < contacts_min[key]:
                    contacts_min[key] = d

        # Prepare sorted report for this ligand instance
        rows = sorted(contacts_min.items(), key=lambda kv: (kv[0][1], kv[0][2], kv[1]))  # by chain, resnum, then distance
        report = {
            "ligand": lig_label,
            "contacts": [
                {
                    "resname": k[0],
                    "chain": k[1],
                    "resnum": k[2],
                    "icode": k[3],
                    "min_distance": round(v, 3),
                }
                for k, v in rows
            ],
        }
        all_reports.append(report)

    # Print a human-readable summary
    for rep in all_reports:
        print(f"\nResidues within {cutoff} Å of {rep['ligand']}:")
        if not rep["contacts"]:
            print("  (none)")
            continue
        for c in rep["contacts"]:
            icode = c["icode"] if c["icode"] and c["icode"] != " " else ""
            print(
                f"  {c['resname']:>3s}  Chain {c['chain']}  {c['resnum']}{icode:1s}   min d = {c['min_distance']:.3f} Å"
            )
        print(f"Total: {len(rep['contacts'])}")

    return all_reports


def main():
    ap = argparse.ArgumentParser(description="Find protein residues within cutoff of a ligand in an mmCIF file.")
    ap.add_argument("cif", help="Input mmCIF file")
    ap.add_argument("--ligand", default="GSH", help="3-letter ligand code (default: GSH)")
    ap.add_argument("--cutoff", type=float, default=4.5, help="Distance cutoff in Å (default: 4.5)")
    ap.add_argument("--csv", default=None, help="Optional CSV output filename")
    ap.add_argument("--include-waters", action="store_true", help="Include waters as potential neighbors (off by default)")
    ap.add_argument(
        "--strict-aa",
        action="store_true",
        help="Only treat standard amino acids as protein (default off = include all amino acids).",
    )
    args = ap.parse_args()

    reports = find_contacts_cif(
        args.cif,
        ligand_name=args.ligand.upper(),
        cutoff=args.cutoff,
        include_waters=args.include_waters,
        standard_aa_only=args.strict-aa if hasattr(args, "strict-aa") else args.strict_aa
        if hasattr(args, "strict_aa") else False,
    )

    if args.csv and reports:
        # Flatten to CSV: one row per contact
        with open(args.csv, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["ligand", "resname", "chain", "resnum", "icode", "min_distance"])
            for rep in reports:
                for c in rep["contacts"]:
                    writer.writerow([rep["ligand"], c["resname"], c["chain"], c["resnum"], c["icode"], c["min_distance"]])
        print(f"\nCSV written to: {args.csv}")


if __name__ == "__main__":
    main()
