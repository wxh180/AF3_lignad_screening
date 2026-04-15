#!/usr/bin/env python3
"""
Analyze ligand contacts across all CIF files under directories matching seed-*.

Features:
1. Scans all CIF files under seed-* directories
2. Computes per-residue minimum ligand distance per structure
3. Fills non-contact residues with cutoff value (default 4.5 Å)
4. Aggregates mean, std, and contact frequency
5. Exports:
   - summary TSV
   - per-structure TSV
   - overall mean-distance plot
   - overall contact-frequency plot
   - one mean-distance plot per chain
   - one contact-frequency plot per chain

Examples:
python contacts_seed_batch.py . --ligand GSH --cutoff 4.5
python contacts_seed_batch.py /path/to/project --ligand ATP --cutoff 4.5 --chain A
"""

import argparse
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from Bio.PDB import MMCIFParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa


def atom_dist(a, b):
    ax, ay, az = a.get_coord()
    bx, by, bz = b.get_coord()
    dx, dy, dz = ax - bx, ay - by, az - bz
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def residue_key(residue):
    """
    Stable residue key:
    (chain_id, resseq, icode, resname)
    """
    chain = residue.get_parent()
    resname = residue.get_resname().strip()
    resseq, icode = residue.get_id()[1], residue.get_id()[2]
    icode = "" if icode in (None, " ") else str(icode)
    return (chain.id, int(resseq), icode, resname)


def residue_label(key):
    chain_id, resseq, icode, resname = key
    return f"{resname} {chain_id}:{resseq}{icode}"


def find_cif_files(root):
    root = Path(root)
    cifs = []
    for seed_dir in root.glob("seed-*"):
        if seed_dir.is_dir():
            cifs.extend(seed_dir.rglob("*.cif"))
    return sorted(cifs)


def extract_structure_contacts(cif_path, ligand_name="GSH", cutoff=4.5, chain_filter=None):
    """
    Returns:
      residue_min_dist: dict[(chain, resseq, icode, resname)] = min_distance_to_ligand
      all_protein_residues: set of residue keys found in the structure
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(cif_path.stem, str(cif_path))

    ligand_residues = []
    protein_atoms = []
    all_protein_residues = set()

    for model in structure:
        for chain in model:
            if chain_filter and chain.id != chain_filter:
                continue

            for residue in chain:
                resname = residue.get_resname().strip()
                hetflag = residue.get_id()[0]

                if resname == ligand_name:
                    ligand_residues.append(residue)
                    continue

                if hetflag.startswith("W"):
                    continue

                if is_aa(residue, standard=False):
                    key = residue_key(residue)
                    all_protein_residues.add(key)
                    protein_atoms.extend(list(residue.get_atoms()))

    if not ligand_residues:
        print(f"WARNING: ligand {ligand_name} not found in {cif_path}")
        return {}, all_protein_residues

    if not protein_atoms:
        print(f"WARNING: no protein atoms found in {cif_path}")
        return {}, all_protein_residues

    ns = NeighborSearch(protein_atoms)
    residue_min_dist = defaultdict(lambda: float("inf"))

    for lig in ligand_residues:
        for latom in lig.get_atoms():
            neighbors = ns.search(latom.get_coord(), cutoff)
            for patom in neighbors:
                pres = patom.get_parent()
                if not is_aa(pres, standard=False):
                    continue
                key = residue_key(pres)
                d = atom_dist(latom, patom)
                if d < residue_min_dist[key]:
                    residue_min_dist[key] = d

    return dict(residue_min_dist), all_protein_residues


def aggregate_contacts(cif_files, ligand_name="GSH", cutoff=4.5, chain_filter=None, fill_value=4.5):
    """
    For each CIF:
      - compute min distance per residue if contacting ligand within cutoff
      - for all protein residues in that structure with no contact, assign fill_value

    Returns:
      summary: list of dicts
      per_structure_data: dict[cif_name][residue_key] = distance
    """
    per_structure_data = {}
    global_residue_keys = set()

    for cif in cif_files:
        residue_min_dist, structure_residues = extract_structure_contacts(
            cif, ligand_name=ligand_name, cutoff=cutoff, chain_filter=chain_filter
        )

        structure_distances = {}
        for rkey in structure_residues:
            if rkey in residue_min_dist:
                structure_distances[rkey] = residue_min_dist[rkey]
            else:
                structure_distances[rkey] = fill_value

        per_structure_data[cif.name] = structure_distances
        global_residue_keys.update(structure_residues)

    summary = []
    for rkey in sorted(global_residue_keys, key=lambda x: (x[0], x[1], x[2], x[3])):
        vals = []
        n_contact = 0
        for _, rmap in per_structure_data.items():
            if rkey in rmap:
                vals.append(rmap[rkey])
                if rmap[rkey] < fill_value:
                    n_contact += 1

        if not vals:
            continue

        chain_id, resseq, icode, resname = rkey
        summary.append({
            "chain": chain_id,
            "resseq": resseq,
            "icode": icode,
            "resname": resname,
            "label": residue_label(rkey),
            "n_structures": len(vals),
            "n_contacts": n_contact,
            "contact_frequency": n_contact / len(vals),
            "mean_distance": float(np.mean(vals)),
            "std_distance": float(np.std(vals, ddof=0)),
        })

    return summary, per_structure_data


def save_summary_tsv(summary, out_tsv):
    with open(out_tsv, "w") as fh:
        fh.write(
            "chain\tresseq\ticode\tresname\tlabel\tn_structures\tn_contacts\tcontact_frequency\tmean_distance\tstd_distance\n"
        )
        for row in summary:
            fh.write(
                f"{row['chain']}\t{row['resseq']}\t{row['icode']}\t{row['resname']}\t"
                f"{row['label']}\t{row['n_structures']}\t{row['n_contacts']}\t"
                f"{row['contact_frequency']:.4f}\t{row['mean_distance']:.4f}\t{row['std_distance']:.4f}\n"
            )


def save_per_structure_tsv(per_structure_data, out_tsv):
    all_keys = sorted(
        {k for rmap in per_structure_data.values() for k in rmap.keys()},
        key=lambda x: (x[0], x[1], x[2], x[3])
    )

    with open(out_tsv, "w") as fh:
        header = ["structure", "chain", "resseq", "icode", "resname", "label", "distance"]
        fh.write("\t".join(header) + "\n")

        for structure_name, rmap in per_structure_data.items():
            for key in all_keys:
                if key not in rmap:
                    continue
                chain, resseq, icode, resname = key
                fh.write(
                    f"{structure_name}\t{chain}\t{resseq}\t{icode}\t{resname}\t{residue_label(key)}\t{rmap[key]:.4f}\n"
                )


def annotate_frequency_labels(ax, rows, threshold=0.5):
    """Annotate residues on the frequency panel with simple staggered offsets."""
    candidates = [row for row in rows if row["contact_frequency"] >= threshold]
    candidates.sort(key=lambda row: (row["resseq"], row["resname"], row["chain"]))

    offset_cycle = [
        (0, 6),
        (0, 16),
        (0, 26),
        (8, 10),
        (-8, 18),
        (8, 28),
        (-8, 36),
    ]

    last_resseq = None
    cluster_index = 0

    for row in candidates:
        resseq = row["resseq"]
        if last_resseq is None or abs(resseq - last_resseq) > 2:
            cluster_index = 0
        else:
            cluster_index += 1

        xytext = offset_cycle[cluster_index % len(offset_cycle)]
        label = f"{row['resname']}{row['resseq']}"
        ax.annotate(
            label,
            (row["resseq"], row["contact_frequency"]),
            textcoords="offset points",
            xytext=xytext,
            ha="center",
            va="bottom",
            fontsize=6.5,
            rotation=45,
            color="black",
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": "none", "alpha": 0.7},
        )
        last_resseq = resseq



def plot_summary(summary, out_png, title=None, fill_value=4.5):
    if not summary:
        print("No data to plot.")
        return

    x = [row["resseq"] for row in summary]
    y = [row["mean_distance"] for row in summary]
    yerr = [row["std_distance"] for row in summary]
    freq = [row["contact_frequency"] for row in summary]

    fig, (ax1, ax2) = plt.subplots(
        2, 1, figsize=(12, 8), dpi=200, sharex=True,
        gridspec_kw={"height_ratios": [3, 1.4]}
    )

    ax1.errorbar(
        x, y, yerr=yerr,
        fmt="o-", markersize=3, linewidth=1,
        ecolor="gray", elinewidth=0.8, capsize=2,
        color="tab:blue", label="Mean distance ± SD"
    )
    ax1.axhline(fill_value, color="red", linestyle="--", linewidth=1, alpha=0.6,
                label=f"no-contact fill = {fill_value:.1f} Å")
    ax1.set_ylabel("Mean ligand-contact distance (Å)")
    ax1.grid(alpha=0.25)
    ax1.set_title(title if title else "Residue contact distance and frequency across seed-* CIF files")
    ax1.legend(frameon=False, loc='best')

    ax2.plot(x, freq, "s-", markersize=3, linewidth=1, color="tab:orange")
    ax2.set_xlabel("Residue index")
    ax2.set_ylabel("Contact frequency")
    ax2.set_ylim(-0.02, 1.02)
    ax2.grid(alpha=0.25)

    annotate_frequency_labels(ax2, summary, threshold=0.5)

    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight")
    plt.close(fig)


def plot_per_chain(summary, out_prefix, fill_value=4.5):
    chains = sorted(set(row["chain"] for row in summary))
    for chain in chains:
        sub = [row for row in summary if row["chain"] == chain]
        if not sub:
            continue

        x = [row["resseq"] for row in sub]
        y = [row["mean_distance"] for row in sub]
        yerr = [row["std_distance"] for row in sub]

        freq = [row["contact_frequency"] for row in sub]

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(12, 8), dpi=200, sharex=True,
            gridspec_kw={"height_ratios": [3, 1.4]}
        )
        ax1.errorbar(
            x, y, yerr=yerr,
            fmt="o-", markersize=3, linewidth=1,
            ecolor="gray", elinewidth=0.8, capsize=2,
            color="tab:blue", label="Mean distance ± SD"
        )
        ax1.axhline(fill_value, color="red", linestyle="--", linewidth=1, alpha=0.6,
                    label=f"no-contact fill = {fill_value:.1f} Å")
        ax1.set_ylabel("Mean ligand-contact distance (Å)")
        ax1.grid(alpha=0.25)
        ax1.set_title(f"Residue contact distance and frequency for chain {chain}")
        ax1.legend(frameon=False, loc='best')

        ax2.plot(x, freq, "s-", markersize=3, linewidth=1, color="tab:orange")
        ax2.set_xlabel("Residue index")
        ax2.set_ylabel("Contact frequency")
        ax2.set_ylim(-0.02, 1.02)
        ax2.grid(alpha=0.25)

        annotate_frequency_labels(ax2, sub, threshold=0.5)

        fig.tight_layout()
        out_png = f"{out_prefix}_chain_{chain}.png"
        fig.savefig(out_png, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved per-chain plot to {out_png}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", help="Root directory containing seed-* subdirectories")
    ap.add_argument("--ligand", required=True, help="Ligand residue name, e.g. GSH")
    ap.add_argument("--cutoff", type=float, default=4.5, help="Distance cutoff in Å (default: 4.5)")
    ap.add_argument("--chain", default=None, help="Optional protein chain ID to restrict analysis")
    ap.add_argument("--out-prefix", default="contacts_seed_batch", help="Output file prefix")
    args = ap.parse_args()

    cif_files = find_cif_files(args.root)
    if not cif_files:
        print("No CIF files found under seed-* directories.")
        return

    print(f"Found {len(cif_files)} CIF files.")

    summary, per_structure_data = aggregate_contacts(
        cif_files,
        ligand_name=args.ligand,
        cutoff=args.cutoff,
        chain_filter=args.chain,
        fill_value=args.cutoff,
    )

    out_prefix = Path(args.out_prefix)

    save_summary_tsv(summary, str(out_prefix) + ".tsv")
    save_per_structure_tsv(per_structure_data, str(out_prefix) + "_per_structure.tsv")

    plot_title = f"{args.ligand} contacts across seed-* CIF files"
    if args.chain:
        plot_title += f" (chain {args.chain})"

    plot_summary(summary, str(out_prefix) + ".png", title=plot_title, fill_value=args.cutoff)
    plot_per_chain(summary, str(out_prefix), fill_value=args.cutoff)

    print(f"Saved summary to {out_prefix}.tsv")
    print(f"Saved per-structure distances to {out_prefix}_per_structure.tsv")
    print(f"Saved combined distance/contact-frequency plot to {out_prefix}.png")


if __name__ == "__main__":
    main()
