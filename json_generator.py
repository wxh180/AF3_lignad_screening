#!/usr/bin/env python3
import requests
import json
import os

# ---- USER INPUT ----
uniprot_ids = [
    "Q16772", "O15217", "P0CG30", "P09211", "P78417",
    "Q9Y2Q3", "P28161", "Q9H4Y5", "P08263", "P30711",
    "O43813", "O14880", "P46439", "Q99735", "P09488",
    "P09210", "P21266", "P10620", "Q7RTV2", "Q03013",
    "P0CG29", "V9HWE9"
]

# Folder to save json files
OUTPUT_DIR = "AF3_json_inputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def fetch_sequence(uniprot_id):
    """Fetch protein sequence from UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    r = requests.get(url)

    if not r.ok:
        print(f" Error fetching {uniprot_id}")
        return None

    fasta = r.text.strip().split("\n")
    seq = "".join(fasta[1:])  # skip header
    return seq


def write_json(uniprot_id, sequence):
    """Write alphafold3 JSON format for one protein + GSH ligand."""
    data = {
        "name": f"AF3_{uniprot_id}",
        "modelSeeds": [1],
        "sequences": [
            {
                "protein": {
                    "id": ["A", "B"],
                    "sequence": sequence
                }
            },
            {
                "ligand": {
                    "id": ["C", "D"],
                    "ccdCodes": ["GSH"]
                }
            }
        ],
        "dialect": "alphafold3",
        "version": 3
    }

    out_path = os.path.join(OUTPUT_DIR, f"{uniprot_id}.json")
    with open(out_path, "w") as f:
        json.dump(data, f, indent=4)

    print(f"✅ Generated: {out_path}")


# ---- MAIN LOOP ----
for uid in uniprot_ids:
    seq = fetch_sequence(uid)
    if seq:
        write_json(uid, seq)
