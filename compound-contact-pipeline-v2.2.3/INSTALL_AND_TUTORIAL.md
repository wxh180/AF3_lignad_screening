# Installation and tutorials

This guide explains how to install and run version 2.2.3 of the compound
contact, restrained-minimization, and binding-evidence pipeline on a Linux
workstation or server.

The project is currently a script-based package. It is run as:

```bash
python compound_contact_pipeline.py ROOT [options]
```

It is not installed from PyPI and does not currently provide a system-wide
`compound-contact-pipeline` command.

## 1. Choose an installation

Use the lightweight installation if you need contact analysis, plots, buried
SASA, clash detection, or imported energy/affinity summaries. Use the complete
Conda environment if you also want restrained OpenMM minimization.

| Installation | Contact analysis | Plots | Energy/affinity import | OpenMM minimization |
|---|---:|---:|---:|---:|
| Lightweight Python environment | yes | yes | yes | no |
| `environment-openmm.yml` | yes | yes | yes | yes |

## 2. Obtain the project files

Keep these files together in one directory:

```text
compound-contact-pipeline/
├── compound_contact_pipeline.py
├── benchmark_sasa.py
├── restrained_openmm.py
├── environment-openmm.yml
├── datasets.example.tsv
├── energy_results.example.tsv
├── affinity_results.example.tsv
└── tests/
```

In the commands below, replace `/path/to/compound-contact-pipeline` with that
directory:

```bash
cd /path/to/compound-contact-pipeline
```

## 3. Lightweight installation

Python 3.10 or newer is required. The following approach does not modify the
system Python installation:

```bash
cd /path/to/compound-contact-pipeline
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install "numpy>=1.24" "matplotlib>=3.7" "biopython>=1.83" "freesasa>=2.2"
```

Verify the installation:

```bash
python compound_contact_pipeline.py --version
python compound_contact_pipeline.py --help
```

The version command should print:

```text
2.2.3
```

To run the automated tests, install Pytest and execute the suite:

```bash
python -m pip install pytest
python -m pytest -q
```

Reactivate this environment in a later terminal with:

```bash
cd /path/to/compound-contact-pipeline
source .venv/bin/activate
```

## 4. Complete OpenMM installation

OpenMM minimization additionally requires OpenMM, OpenFF Toolkit,
openmmforcefields, RDKit, and AmberTools. Conda or Mamba is recommended because
these packages include compiled dependencies.

If Conda is available:

```bash
cd /path/to/compound-contact-pipeline
conda env create -f environment-openmm.yml
conda activate compound-contact-openmm
python compound_contact_pipeline.py --version
```

Mamba can be substituted for the environment-creation command:

```bash
mamba env create -f environment-openmm.yml
conda activate compound-contact-openmm
```

### FreeSASA requirement and v2.2.3 SASA acceleration

v2.2.3 requires the FreeSASA Python bindings for all full scientific analyses.
The supplied `environment-openmm.yml` installs `freesasa-python>=2.2` from
conda-forge. To update an existing environment, run:

```bash
conda env update -n compound-contact-openmm -f environment-openmm.yml --prune
conda activate compound-contact-openmm
python - <<'PY'
import freesasa
print("FreeSASA import: OK")
print(getattr(freesasa, "__version__", "version not exposed"))
PY
```

The production calculation uses FreeSASA Lee-Richards with a 1.4 A probe, 20
slices per atom, the same element-based atomic radii as the retained Biopython
reference, and one FreeSASA thread per worker process. With `--workers 16`, up
to 16 structures are analyzed concurrently without nested FreeSASA threading.
Treat v2.2.3 SASA values as a new quantitative baseline and recompute all
compounds before structural-SAR comparisons.

Before a large production run on a new server, benchmark one representative
minimized structure against the retained Biopython Shrake-Rupley reference:

```bash
python benchmark_sasa.py /path/to/minimized_model.cif \
  --ligand LIG_F \
  --protein-chain A,B,C \
  --ligand-chain F \
  --repeats 1
```

The output is a two-column table containing both buried-SASA values, both
runtimes, the percentage difference in buried SASA, and
`speedup_vs_biopython`. For large complexes, start with `--repeats 1` because
the Biopython reference path is intentionally slow.

### AlphaFold 3 custom ligands without bond records

AlphaFold 3 may write custom-ligand coordinates without `_struct_conn` or
`_chem_comp_bond` records. v2.2.1 can recover this case automatically when
`--ligand-smiles` or `--ligand-file` supplies the exact chemistry. RDKit first
proposes heavy-atom connectivity from the ligand coordinates; the pipeline
accepts that proposal only when its element-labeled graph matches the supplied
chemistry. Chemical bond orders and GAFF parameters are still taken from the
supplied SMILES/SDF. If the coordinate-derived graph disagrees, minimization
stops with a chemistry-validation error.

When every minimization fails or is rejected, inspect
`minimization_summary.tsv`; v2.2.1 writes this diagnostic table before exiting
and also prints each failure/rejection reason in the progress log.

If an environment with this name already exists, update it instead:

```bash
conda env update -n compound-contact-openmm -f environment-openmm.yml --prune
conda activate compound-contact-openmm
```

Confirm that OpenMM can see its computational platforms:

```bash
python - <<'PY'
from openmm import Platform

for index in range(Platform.getNumPlatforms()):
    platform = Platform.getPlatform(index)
    print(index, platform.getName())
PY
```

On a GPU server, `CUDA` should appear if the OpenMM build and NVIDIA runtime
are compatible. CPU minimization remains available when CUDA is not listed.

Run the tests after installation:

```bash
python -m pytest -q
```

## 5. Organize the input structures

The pipeline recursively finds CIF files below directories named `seed-*`.
For example:

```text
/data/ligand_screen/
├── C19/
│   ├── seed-1/model.cif
│   ├── seed-2/model.cif
│   └── seed-3/model.cif
└── GSH/
    ├── seed-1/model.cif
    └── seed-2/model.cif
```

In manifest-free mode, the parent of each outermost `seed-*` directory is
treated as a compound dataset. In this example, the inferred compound names
are `C19` and `GSH`.

Use a manifest for multi-compound work, whenever chain assignments differ, or
when minimization is requested. It makes compound and ligand identities
explicit and reproducible.

## Tutorial 1: Analyze one compound without minimization

Assume the input directory is:

```text
/data/WorkDir/Wei/AF3_local/NSD2_C19/
├── seed-1/model.cif
├── seed-2/model.cif
└── seed-3/model.cif
```

Assume the ligand residue name is `WXI`, the protein is chain `A`, and the
ligand is chain `B`:

```bash
cd /path/to/compound-contact-pipeline
source .venv/bin/activate

python compound_contact_pipeline.py \
  /data/WorkDir/Wei/AF3_local/NSD2_C19 \
  --ligand WXI \
  --protein-chain A \
  --ligand-chain B \
  --cutoff 4.5 \
  --clash-cutoff 2.0 \
  --out-dir /data/WorkDir/Wei/AF3_local/NSD2_C19/contact_results
```

The analysis uses original CIF coordinates. It does not modify any input
structure.

Inspect these outputs first:

- `compound_summary.tsv`: dataset-level structural summary.
- `structure_summary.tsv`: status, clashes, and buried SASA for each CIF.
- `compound_residue_summary.tsv`: residue contact frequency and distance.
- `qc.tsv`: warnings and invalid structures.
- `plots/`: per-chain plots and cross-compound summaries.

## Tutorial 2: Analyze several compounds with a manifest

Create `/data/ligand_screen/datasets.tsv` as a tab-separated file:

```text
dataset_dir	compound	ligand_resname	protein_chains	ligand_chain	ligand_smiles	ligand_file	enabled
C19	C19	WXI	A	B		C19/C19.sdf	true
GSH	GSH	GSH	A,C	L		GSH/GSH.sdf	true
```

Important rules:

- `dataset_dir` and `ligand_file` are relative to the root directory.
- `compound` is the identifier used in tables and figure labels.
- `protein_chains` accepts comma-separated chain IDs.
- Use exactly one of `ligand_smiles` or `ligand_file` for minimization.
- A row with `enabled=false` is skipped.
- Do not combine a manifest with the command-line `--ligand`,
  `--protein-chain`, or `--ligand-chain` options.

Run the analysis:

```bash
conda activate compound-contact-openmm

python /path/to/compound-contact-pipeline/compound_contact_pipeline.py \
  /data/ligand_screen \
  --manifest /data/ligand_screen/datasets.tsv \
  --cutoff 4.5 \
  --workers 4 \
  --out-dir /data/ligand_screen/contact_results
```

`--workers 4` parallelizes contact-only structure analysis. The output is
deterministic across worker counts.

The most useful multi-compound figures are:

- `plots/contact_frequency_heatmap_chain_<chain>.png`
- `plots/capped_distance_heatmap_chain_<chain>.png`
- `plots/compound_structural_summary.png`

Gray heatmap cells mean the residue was absent for that compound; they are not
zero-frequency contacts.

## Tutorial 3: Run restrained minimization on the CPU

Minimization requires exact ligand chemistry. Prefer a curated, one-molecule
SDF containing bond orders, formal charge, and stereochemistry. The ligand's
heavy-atom identity and bond graph must match the selected CIF residue.

For a single manifest-free dataset:

```bash
conda activate compound-contact-openmm

python /path/to/compound-contact-pipeline/compound_contact_pipeline.py \
  /data/WorkDir/Wei/AF3_local/NSD2_C19 \
  --ligand WXI \
  --protein-chain A \
  --ligand-chain B \
  --ligand-file C19.sdf \
  --minimize openmm \
  --analysis-source minimized \
  --openmm-platform CPU \
  --workers 4 \
  --out-dir /data/WorkDir/Wei/AF3_local/NSD2_C19/minimized_results
```

`--ligand-file C19.sdf` is resolved relative to the root directory, so the
example expects:

```text
/data/WorkDir/Wei/AF3_local/NSD2_C19/C19.sdf
```

The protocol performs three restrained L-BFGS stages. Protein backbone and
distant protein atoms remain restrained while pocket side chains and the
ligand relax conservatively.

Original CIF files are never overwritten. Completed candidates are written to:

```text
minimized_results/minimized_structures/accepted/
minimized_results/minimized_structures/rejected/
```

Review:

- `minimization_summary.tsv`
- `structure_summary.tsv`
- `qc.tsv`
- `plots/minimization_before_after.png`

An accepted structure must satisfy all configured geometry gates. A rejected
candidate is preserved for inspection but is not used when
`--analysis-source minimized` is selected.

### Use a compound-ID/SMILES table

When a folder name contains a compound ID such as `B56a_EOAI10089186`, the
pipeline can retrieve the exact SMILES automatically. Copy or edit
`compound_smiles.example.tsv`; its required headers are:

```text
Compound ID	Smiles
```

Then run from any directory using an explicit table path:

```bash
conda activate compound-contact-openmm

python /path/to/compound-contact-pipeline/compound_contact_pipeline.py \
  /EM2_NAS/AF3_workdir/AB56aC_EV1057/af_output/B56a_EOAI10089186/B56a_EOAI10089186 \
  --ligand LIG \
  --protein-chain A,B,C \
  --ligand-chain F \
  --smiles-table /path/to/compound-contact-pipeline/compound_smiles.example.tsv \
  --minimize openmm \
  --analysis-source minimized \
  --openmm-platform CUDA \
  --openmm-device-index 0 \
  --workers 16 \
  --out-dir minimized_cuda_results
```

Replace `LIG` with the exact ligand residue name used in the CIF. The lookup
does not guess residue identity; it supplies ligand bond orders and chemistry
after the ligand residue is selected. `run_manifest.tsv` records the matched
compound ID and the table used. The supplied source list contained `#NAME?`
for `EOAI10404813`, so that entry is intentionally absent until a valid SMILES
is available.

## Tutorial 4: Run minimization on an NVIDIA GPU

First verify that OpenMM lists the `CUDA` platform as shown in the installation
section. Then select a GPU explicitly:

```bash
conda activate compound-contact-openmm

python /path/to/compound-contact-pipeline/compound_contact_pipeline.py \
  /data/WorkDir/Wei/AF3_local/NSD2_C19 \
  --ligand WXI \
  --protein-chain A \
  --ligand-chain B \
  --ligand-file C19.sdf \
  --minimize openmm \
  --analysis-source minimized \
  --openmm-platform CUDA \
  --openmm-device-index 0 \
  --workers 16 \
  --out-dir /data/WorkDir/Wei/AF3_local/NSD2_C19/minimized_cuda_results
```

Use `--openmm-device-index 1` for the second visible GPU, and so forth.
`--workers 16` uses up to 16 CPU processes for the analysis phases. OpenMM
minimization itself remains sequential on the selected GPU, so this does not
launch 16 GPU minimizations at once.

SASA scheduling is automatic in v2.2. With `--analysis-source minimized`, the
original structures and minimized-candidate QC skip SASA; only accepted
minimized structures receive the final full SASA calculation. In v2.2.3 that
final calculation uses compiled FreeSASA Lee-Richards. There is no `--sasa-mode` option.

Progress and timing messages are printed to stderr by default. Use `--quiet`
to suppress normal progress while retaining warnings and errors.

Interpret the main progress labels as follows:

```text
OpenMM preparation...       CPU-side topology/chemistry/system setup
hydrogen_relaxation [CUDA]  LocalEnergyMinimizer stage on selected CUDA device
pocket_relaxation [CUDA]    LocalEnergyMinimizer stage on selected CUDA device
gentle_relaxation [CUDA]    LocalEnergyMinimizer stage on selected CUDA device
candidate geometry QC       CPU-side post-minimization geometry/contact analysis
Final minimized scientific analysis  CPU-side full analysis including FreeSASA
```

High CPU utilization with little or no GPU utilization during `OpenMM
preparation...` is expected. GPU activity should be visible during the named
CUDA minimization stages.

If CUDA initialization fails, repeat the same run with
`--openmm-platform CPU`. This distinguishes a GPU/runtime problem from a ligand
or protein-parameterization problem.

## Tutorial 5: Compare original and minimized analyses

The safest comparison is to run minimization twice into two output
directories, changing only `--analysis-source`.

Set the data root once:

```bash
DATASET_ROOT=/data/ligand_screen
```

Analyze original coordinates while still collecting minimization QC:

```bash
python compound_contact_pipeline.py "$DATASET_ROOT" \
  --manifest "$DATASET_ROOT/datasets.tsv" \
  --minimize openmm \
  --analysis-source original \
  --out-dir "$DATASET_ROOT/results_original"
```

Analyze only accepted minimized coordinates:

```bash
python compound_contact_pipeline.py "$DATASET_ROOT" \
  --manifest "$DATASET_ROOT/datasets.tsv" \
  --minimize openmm \
  --analysis-source minimized \
  --out-dir "$DATASET_ROOT/results_minimized"
```

Compare:

- `structure_summary.tsv` for clash and SASA changes.
- `compound_residue_summary.tsv` for altered contacts.
- `plots/contact_frequency_heatmap_chain_<chain>.png` for contact patterns.
- `minimization_summary.tsv` for ligand/backbone RMSD and acceptance.

Do not interpret the reported OpenMM potential-energy decrease as binding
energy. It compares the same parameterized complex before and after local
relaxation and is only a minimization quality-control value.

## Tutorial 6: Adjust minimization gates carefully

The defaults are:

```text
pocket radius                 6.0 Å
maximum protein backbone RMSD 0.5 Å
maximum ligand RMSD           1.5 Å
minimization tolerance        10 kJ mol⁻¹ nm⁻¹
```

Override them only for a documented structural reason:

```bash
DATASET_ROOT=/data/ligand_screen

python compound_contact_pipeline.py "$DATASET_ROOT" \
  --manifest "$DATASET_ROOT/datasets.tsv" \
  --minimize openmm \
  --analysis-source minimized \
  --pocket-radius 7.0 \
  --max-backbone-rmsd 0.5 \
  --max-ligand-rmsd 1.5 \
  --minimization-tolerance 10 \
  --out-dir "$DATASET_ROOT/results_minimized"
```

Do not loosen the RMSD thresholds simply to convert a rejected pose into an
accepted one. A large ligand displacement or severe ligand–backbone overlap
usually calls for pose rebuilding or redocking.

## Tutorial 7: Import calculated energy values

The pipeline does not calculate binding free energy from contacts or from the
minimized complex potential energy. It imports results from a separate,
method-specific calculation.

Create a tab-separated `energy_results.tsv`:

```text
compound	method	value_kcal_mol	structure_id	uncertainty_kcal_mol	reference_compound	replicate
C19	mmgbsa	-31.2	C19/seed-1/model.cif	1.8		1
C19	mmgbsa	-29.9	C19/seed-2/model.cif	1.6		2
GSH	rbfe	1.4		0.4	C19	1
```

Then run:

```bash
python compound_contact_pipeline.py /data/ligand_screen \
  --manifest /data/ligand_screen/datasets.tsv \
  --energy-results /data/ligand_screen/energy_results.tsv \
  --out-dir /data/ligand_screen/results_with_energy
```

Use consistent method labels such as `vina`, `mmgbsa`, or `rbfe`. Different
methods are summarized and plotted separately. For RBFE, every row must name
the same reference compound for a given comparison series.

Inspect:

- `energy_values_normalized.tsv`
- `energy_summary.tsv`
- `plots/energy_<method>.png`

Vina scores are docking scores; MM/GBSA values are endpoint estimates; RBFE
values are relative free energies. Do not pool them into a single numerical
ranking.

## Tutorial 8: Import experimental affinity measurements

Create a tab-separated `affinity_results.tsv`:

```text
compound	metric	value	unit	temperature_K	replicate
C19	Kd	25	nM	298.15	1
C19	Kd	31	nM	298.15	2
GSH	IC50	2.1	uM	298.15	1
```

Run:

```bash
python compound_contact_pipeline.py /data/ligand_screen \
  --manifest /data/ligand_screen/datasets.tsv \
  --affinity-results /data/ligand_screen/affinity_results.tsv \
  --out-dir /data/ligand_screen/results_with_affinity
```

Supported metrics are `Kd`, `Ki`, and `IC50`; supported units are `M`, `mM`,
`uM`, `nM`, and `pM`. Kd is automatically converted to standard binding free
energy. Ki is converted only when `--convert-ki` is supplied. IC50 is never
automatically converted because it is assay-dependent.

Inspect:

- `affinity_values_normalized.tsv`
- `affinity_summary.tsv`
- `plots/experimental_affinity.png`

## Tutorial 9: Rerun safely after changing settings

The pipeline refuses to write into a nonempty output directory unless
`--overwrite` is supplied:

```bash
DATASET_ROOT=/data/ligand_screen

python compound_contact_pipeline.py "$DATASET_ROOT" \
  --manifest "$DATASET_ROOT/datasets.tsv" \
  --out-dir "$DATASET_ROOT/contact_results" \
  --overwrite
```

Overwrite is inventory-controlled. The program removes only files recorded in
the preceding `generated_files.json`; it does not recursively delete the
output directory or unrelated files.

For an important comparison, prefer a new output directory rather than
overwriting the earlier analysis.

## 6. Reading exit codes and QC results

Check the shell exit status immediately after a run:

```bash
echo $?
```

| Exit status | Meaning |
|---:|---|
| `0` | Requested analysis completed successfully |
| `1` | Partial result, such as a plotting failure or rejected/failed minimization |
| `2` | Configuration error or no usable structure for the requested analysis source |

A status of `1` can still leave valid tables and original-coordinate results.
Always inspect `qc.tsv` and `structure_summary.tsv` before drawing conclusions.

## 7. Common problems

### No CIF files found

Confirm that at least one directory name begins with `seed-` and that a CIF
file occurs somewhere below it:

```bash
DATASET_ROOT=/data/ligand_screen
find "$DATASET_ROOT" -type f -path '*/seed-*/*.cif' -print | head
```

### Ligand not found or several ligands detected

Use an explicit manifest or supply `--ligand RESNAME` and, when necessary,
`--ligand-chain CHAIN`. Confirm the actual CIF residue and chain names in
ChimeraX, Coot, or the mmCIF atom-site records.

### Ligand chemistry does not match the CIF

The minimization chemistry source must represent the same heavy atoms,
elements, connectivity, formal charge, and stereochemistry as the selected CIF
ligand. Regenerate or curate the SDF rather than allowing atom identities to be
guessed from coordinates.

### No or multiple compound IDs match the SMILES table

The script matches complete IDs embedded in the inferred compound name,
dataset path, or structure path. Rename the dataset to include exactly one
table ID (for example `B56a_EOAI10089186`), correct the table, or provide an
explicit `--ligand-smiles`/`--ligand-file` override.

### Undefined stereochemistry

Use an isomeric SMILES or a curated SDF with defined stereocenters. The
pipeline rejects ambiguous stereochemistry rather than choosing an arbitrary
isomer.

### Protein force-field template failure

Inspect the model for missing atoms, incomplete residues, nonstandard residues,
covalent modifications, metal coordination, or ambiguous termini. OpenMM
minimization is not appropriate until the chemical model is explicitly
prepared.

### Candidate rejected after minimization

Inspect `minimization_summary.tsv` for:

- fixed 2 Å clashes before and after;
- element-aware van der Waals clashes;
- protein-backbone RMSD;
- ligand-heavy-atom RMSD after protein alignment; and
- the acceptance reason.

Do not weaken the acceptance gates automatically. Severe overlaps usually
need rebuilding or redocking.

### Cryo-EM structure has a map

Use OpenMM only for conservative, map-free geometry preparation. Perform the
final model refinement in Phenix with the experimental map and a validated
ligand restraint CIF.

## 8. Recommended routine workflow

For an AF3 or docking ensemble:

1. Run contact-only analysis on all original structures.
2. Inspect `qc.tsv`, clash counts, contact consistency, and predicted-model
   confidence information from the original modeling workflow.
3. Prepare exact ligand chemistry as a curated SDF.
4. Run restrained minimization while retaining original-coordinate analysis.
5. Inspect every rejected candidate and the before/after QC figure.
6. Run the accepted-minimized analysis into a separate output directory.
7. Use Vina only for rapid triage; use consistently sampled MM/GBSA for an
   approximate energetic comparison or RBFE for a suitable congeneric series.
8. Treat experimental Kd under a controlled assay protocol as the decision
   standard.

Contact frequency, buried SASA, pose consistency, and minimized potential
energy are useful structural evidence. None of them independently establishes
binding affinity.
