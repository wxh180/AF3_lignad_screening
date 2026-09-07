# Compound contact, restrained-minimization, and binding-evidence pipeline

`compound_contact_pipeline.py` recursively analyzes protein–ligand contacts in
all CIF files below `seed-*` directories, compares compounds across structural
ensembles, and combines the structural results with optional computed-energy
or experimental-affinity tables.

The script deliberately keeps four evidence classes separate:

1. original or conservatively minimized pose geometry: contacts, distances,
   clashes, buried solvent-accessible
   surface area (SASA), and contact-set consistency;
2. same-system OpenMM potential-energy change used only as minimization QC;
3. method-labeled computational results such as Vina, MM/GBSA, or RBFE; and
4. experimental Kd, Ki, or IC50 measurements.

Structural pose metrics are useful for quality control and prioritization, but
they are not binding-affinity measurements and are never reported as such.

For step-by-step environment setup and worked examples, see
[Installation and tutorials](INSTALL_AND_TUTORIAL.md).

## Requirements

- Python 3.10 or newer
- NumPy
- Matplotlib
- Biopython

Install the runtime dependencies in a virtual environment:

```bash
python -m pip install "numpy>=1.24" "matplotlib>=3.7" "biopython>=1.83"
```

For development and tests, also install `pytest`:

```bash
python -m pip install pytest
python -m pytest -q
```

Restrained minimization has a separate optional environment because OpenMM,
ligand parameterization, and chemistry toolkits are not needed for contact-only
runs:

```bash
conda env create -f environment-openmm.yml
conda activate compound-contact-openmm
python compound_contact_pipeline.py --version
```

The expected version for this release is `2.2.1`.

## v2.2 performance behavior

SASA scheduling is automatic; there is no `--sasa-mode` option. With
`--analysis-source minimized`, the original structures and minimized-candidate
QC use fast geometry analysis without SASA. Only accepted minimized structures
receive the final full scientific analysis with buried SASA. With
`--analysis-source original`, the original structures receive full SASA and
minimized candidates remain geometry-QC only.

`--workers N` controls CPU analysis concurrency only. Restrained OpenMM
minimization remains sequential on the platform/device selected by
`--openmm-platform` and `--openmm-device-index`, so increasing `--workers` does
not launch multiple minimizations on one GPU. Progress and timing messages go
to stderr by default. `--quiet` suppresses those normal progress messages
without suppressing warnings or errors.

For a 16-core workstation using CUDA device 0:

```bash
python compound_contact_pipeline.py ./ \
  --ligand LIG_F \
  --protein-chain A,B,C \
  --ligand-chain F \
  --ligand-smiles "$LIGAND_SMILES" \
  --minimize openmm \
  --analysis-source minimized \
  --openmm-platform CUDA \
  --openmm-device-index 0 \
  --workers 16 \
  --out-dir minimized_cuda_results
```

## v2.2.1 AF3 ligand compatibility

AlphaFold 3 custom-ligand mmCIF files may contain ligand coordinates without
`_struct_conn` or `_chem_comp_bond` records. During OpenMM preparation, v2.2.1
handles this conservatively: when the selected ligand has zero explicit
heavy-atom bonds, RDKit proposes connectivity from the 3D coordinates and the
resulting element-labeled graph must be isomorphic to the supplied SMILES/SDF
chemistry before any bonds are added to the OpenMM topology. Bond orders,
aromaticity, formal charges, and GAFF parameterization continue to come from
the supplied chemistry, not from coordinate guessing. Partially bonded or
inconsistent ligand topologies still fail validation rather than being
silently repaired.

If no minimized structure is accepted, `minimization_summary.tsv` is still
written before the pipeline reports `no valid minimized structures`, and live
progress now prints the specific failure/rejection reason.

## Quick start

With an explicit dataset manifest:

```bash
python compound_contact_pipeline.py /path/to/project \
  --manifest datasets.tsv \
  --cutoff 4.5 \
  --out-dir compound_contact_results
```

Without a manifest, the parent directory of each outermost `seed-*` directory
is treated as one dataset, and that parent's basename becomes the compound
name:

```bash
python compound_contact_pipeline.py /path/to/project \
  --ligand GSH \
  --protein-chain A \
  --ligand-chain L \
  --cutoff 4.5 \
  --out-dir compound_contact_results
```

The script processes only the first model in each CIF. A warning is written to
`qc.tsv` when additional models are present.

## Input layouts

### Manifest-free discovery

Discovery is recursive. For this layout, the inferred compounds are `C19` and
`GSH`:

```text
project/
├── screen/
│   ├── C19/
│   │   ├── seed-1/model.cif
│   │   └── seed-2/model.cif
│   └── GSH/
│       ├── seed-1/model.cif
│       └── seed-2/model.cif
└── compound_contact_pipeline.py
```

Only outermost matching seed directories are selected, so a nested `seed-*`
directory cannot cause a CIF to be counted twice. Structure IDs are unique,
root-relative POSIX paths such as `screen/C19/seed-1/model.cif`.

If `--ligand` is omitted, each structure must contain exactly one ligand
candidate after amino acids, water, and common monatomic ions are excluded.
All valid structures assigned to one compound must resolve to the same ligand
residue name.

Repeat `--protein-chain` to select more than one protein chain:

```bash
--protein-chain A --protein-chain C
```

Comma-separated values such as `--protein-chain A,C` are also accepted.
Protein-chain filtering and `--ligand-chain` filtering are independent.

### Manifest-controlled discovery

Copy `datasets.example.tsv` and edit one row per dataset directory. Paths in
`dataset_dir` are relative to `ROOT` and must stay inside it.

| Column | Required | Meaning |
|---|---:|---|
| `dataset_dir` | yes | Directory containing one or more `seed-*` trees |
| `compound` | yes | Stable compound identifier used in all outputs |
| `ligand_resname` | yes | CIF ligand residue name |
| `protein_chains` | no | Comma-separated protein chain IDs |
| `ligand_chain` | no | Ligand chain ID |
| `ligand_smiles` | minimization only | Exact isomeric SMILES; mutually exclusive with `ligand_file` |
| `ligand_file` | minimization only | Root-relative, one-molecule SDF; mutually exclusive with `ligand_smiles` |
| `enabled` | no | `true` by default; use `false` to skip the row |

Manifest ligand and chain settings are authoritative. Do not combine
`--manifest` with `--ligand`, `--protein-chain`, or `--ligand-chain`.

## Command-line options

```text
ROOT                         root of the recursive dataset search
--manifest PATH              authoritative dataset TSV
--ligand RESNAME             default ligand residue name
--protein-chain CHAIN        repeatable protein-chain filter
--ligand-chain CHAIN         independent ligand-chain filter
--minimize openmm            opt in to restrained OpenMM minimization
--analysis-source SOURCE     original (default) or accepted minimized poses
--ligand-smiles SMILES       exact chemistry for one manifest-free compound
--ligand-file PATH           root-relative one-molecule SDF alternative
--smiles-table PATH          TSV/CSV mapping Compound ID to Smiles
--openmm-platform PLATFORM   CPU (default), CUDA, OpenCL, or Reference
--openmm-device-index INDEX  CUDA/OpenCL device index
--minimize-ph VALUE          hydrogen-addition pH; default 7.4
--pocket-radius ANGSTROM     mobile-pocket radius; default 6.0
--max-backbone-rmsd VALUE    acceptance gate in Å; default 0.5
--max-ligand-rmsd VALUE      acceptance gate in Å; default 1.5
--minimization-tolerance V   L-BFGS tolerance in kJ mol⁻¹ nm⁻¹; default 10
--cutoff ANGSTROM            inclusive contact cutoff; default 4.5
--clash-cutoff ANGSTROM      inclusive atom-pair clash cutoff; default 2.0
--energy-results PATH        method-labeled computed-energy TSV
--affinity-results PATH      experimental-affinity TSV
--convert-ki                 explicitly permit Ki-to-ΔG° conversion
--label-frequency FRACTION   contact-label threshold from 0 to 1; default 0.5
--workers N                  CPU analysis process count; default 1
--fail-fast                  stop at the first invalid CIF; requires one worker
--quiet                      suppress normal progress and timing messages
--out-dir PATH               output directory
--overwrite                  replace only files in the preceding inventory
```

Without `--fail-fast`, a malformed or otherwise invalid CIF is recorded and
the run continues. The command succeeds when at least one structure is valid.
Configuration errors and runs with no valid structures return status 2. A plot
failure preserves completed tables and returns status 1. Requested
minimization with any rejected, failed, or unattempted structure also returns
status 1 after writing the available original-analysis outputs. Configuration
errors and minimized analysis with no accepted valid structure return status 2.

### Automatic SMILES lookup

Use `--smiles-table` when compound IDs are embedded in dataset names, for
example `B56a_EOAI10089186`. The table must be TSV or CSV and contain columns
named `Compound ID` and `Smiles` (capitalization and spaces are ignored).
`compound_smiles.example.tsv` contains the valid entries from the supplied
EOAI list.

```bash
python compound_contact_pipeline.py \
  /path/to/B56a_EOAI10089186 \
  --ligand LIG \
  --protein-chain A,B,C \
  --ligand-chain F \
  --smiles-table compound_smiles.example.tsv \
  --minimize openmm \
  --analysis-source minimized \
  --openmm-platform CUDA \
  --openmm-device-index 0 \
  --workers 16 \
  --out-dir minimized_cuda_results
```

Each unresolved structure must match exactly one table ID using a complete
non-alphanumeric-delimited token. Missing, ambiguous, partially blank, and
duplicate IDs are errors. Explicit `--ligand-smiles`, `--ligand-file`, or
manifest chemistry takes priority and is not replaced. The matched ID and
chemistry source are recorded in `run_manifest.tsv`.

## Structural calculations

### Residue contacts and capped distances

A residue is a contact when any selected protein atom–ligand atom distance is
less than or equal to `--cutoff`. For every protein residue present in a valid
structure:

- a contacting residue receives its minimum atom-pair distance;
- a present non-contact residue receives the cutoff as its
  `capped_distance_angstrom`; and
- `min_contact_distance_angstrom` remains blank for a non-contact residue.

Contact frequency uses `n_contacts / n_present`. A residue absent from one
structure is not silently treated as a non-contact and does not enter that
residue's denominator. Cross-compound heatmaps show an absent
compound–residue pair in gray rather than as zero.

Summary standard deviations are population standard deviations (`ddof=0`).

### Clashes

`clash_pair_count` is the number of selected protein atom–ligand atom pairs at
or below `--clash-cutoff`. It is a geometry warning, not an energy term.

### Buried SASA

Biopython's Shrake–Rupley calculation uses a 1.4 Å probe. The script reports:

```text
buried_sasa_total_A2 = SASA(receptor) + SASA(ligand) - SASA(complex)
interface_area_A2    = buried_sasa_total_A2 / 2
```

SASA failure is nonfatal and produces blank SASA fields plus a QC warning.
Buried surface area can describe interface size but is not a binding free
energy.

### Pose consistency

For each compound, the script calculates the Jaccard similarity between every
pair of valid structures' contacting-residue sets and reports their mean. Two
empty contact sets have Jaccard similarity 1.0. A single valid structure has no
pairwise comparison, so its Jaccard summary is blank.

## Restrained steric-clash minimization

Minimization is opt-in and never overwrites an input CIF. Each structure needs
exact ligand chemistry from one—and only one—manifest `ligand_smiles` or
`ligand_file` field. The CIF must contain the ligand's heavy atoms and explicit
heavy-atom connectivity; bond orders, formal charge, stereochemistry, and
protonation are supplied by the SMILES or one-molecule SDF rather than guessed
from coordinates.

For a manifest-free C19/NSD2 dataset containing one compound, a typical CPU run
is:

```bash
python compound_contact_pipeline.py /data/WorkDir/Wei/AF3_local/NSD2_C19 \
  --ligand WXI \
  --protein-chain A \
  --ligand-chain B \
  --ligand-file C19.sdf \
  --minimize openmm \
  --analysis-source minimized \
  --out-dir C19_contact_minimized_results
```

Use `--openmm-platform CUDA --openmm-device-index 0` to select one CUDA GPU.
OpenMM minimization is deliberately sequential on that selected device.
`--workers N` is still useful because it parallelizes the CPU analysis phases
before and after minimization without creating concurrent GPU jobs.

The named `openmm_restrained_v1` protocol uses Amber ff14SB for protein,
GAFF 2.2.20 for ligand, nonperiodic `NoCutoff`, hydrogen-bond constraints, and
hydrogens added at pH 7.4. It runs three bounded L-BFGS stages:

| Stage | Iterations | Backbone | Distant protein | Pocket side chains | Ligand |
|---|---:|---:|---:|---:|---:|
| Hydrogen relaxation | 500 | 10 | 10 | 10 | 10 |
| Pocket relaxation | 1,000 | 10 | 10 | 0 | 1 |
| Gentle relaxation | 1,000 | 2 | 2 | 0 | 0.2 |

Force constants are in kcal mol⁻¹ Å⁻². Pocket residues have any protein atom
within 6 Å of a ligand heavy atom. A candidate is accepted only when protein
backbone RMSD is at most 0.5 Å, ligand-heavy-atom RMSD after protein alignment
is at most 1.5 Å, and neither the fixed 2 Å clash count nor the element-aware
van der Waals clash count increases. Change the two RMSD gates or pocket radius
only with an explicit scientific reason; a rejected severe overlap generally
needs pose rebuilding or redocking, not weaker gates.

Accepted and rejected completed candidates are written separately under
`minimized_structures/accepted/` and `minimized_structures/rejected/`.
`--analysis-source minimized` includes accepted candidates only. Rejected,
failed, and unattempted candidates become explicit invalid rows; the pipeline
never silently substitutes original coordinates. `--analysis-source original`
keeps the original contact analysis while still producing minimization QC.

The reported initial/final potential energies compare the same parameterized
system before and after relaxation. They are useful for detecting a failed or
pathological minimization, but are not binding energies and must not rank
different compounds. The restraint force is excluded from these physical
energy queries.

Implementation references: [OpenMM L-BFGS minimizer](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html),
[OpenMM positional restraints](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Restraining%20Atom%20Positions.html), and
[openmmforcefields SystemGenerator](https://github.com/openmm/openmmforcefields/blob/main/openmmforcefields/generators/system_generators.py).

### Minimization troubleshooting

The v2.2 progress labels show which resource is expected to be active:

```text
OpenMM preparation...       CPU-side topology/chemistry/system setup
hydrogen_relaxation [CUDA]  LocalEnergyMinimizer stage on selected CUDA device
pocket_relaxation [CUDA]    LocalEnergyMinimizer stage on selected CUDA device
gentle_relaxation [CUDA]    LocalEnergyMinimizer stage on selected CUDA device
candidate geometry QC       CPU-side post-minimization geometry/contact analysis
Final minimized scientific analysis  CPU-side full analysis including SASA
```

High CPU use with little or no GPU activity during `OpenMM preparation...` is
expected. GPU utilization should become visible during the named CUDA
minimization stages.

- An `environment-openmm.yml` dependency message means the base environment is
  active; create and activate the optional conda environment shown above.
- A ligand template or element/bond-graph mismatch usually means the SMILES/SDF
  does not exactly match the selected CIF residue, the CIF lacks explicit
  ligand heavy-atom connectivity, or the residue has missing heavy atoms.
- Undefined stereochemistry is rejected. Supply an isomeric SMILES or curated
  SDF rather than allowing the toolkit to choose a stereoisomer.
- A missing protein force-field template usually indicates incomplete residues,
  nonstandard covalent chemistry, or terminal-state ambiguity. Prepare the
  model explicitly; do not reinterpret the error as weak binding.
- CUDA/OpenCL device indices are accepted only with the corresponding GPU
  platform. Re-run on CPU to separate platform setup problems from topology or
  chemistry problems.

For a map-supported cryo-EM C19/NSD2 model, use this map-free OpenMM stage only
for conservative geometry preparation. Finalize the experimental model with
Phenix real-space refinement using the cryo-EM map and a validated ligand
restraint CIF; OpenMM minimization does not replace map-aware refinement.

## Binding energy and physical affinity

This pipeline does not invent binding energy from contact distances. It imports
results from a separately executed, explicitly named method, validates their
identity, and visualizes each method independently. A defensible evaluation
ladder is:

1. **Docking score (fast triage).** Redock or locally optimize consistently,
   then import Vina results with `method=vina`. Vina uses an empirically tuned
   scoring function and receptor/ligand approximations; its score is not a
   measured affinity and should not be mixed numerically with MM/GBSA or RBFE.
2. **Ensemble MM/GBSA (intermediate cost).** Parameterize every complex using
   the same force-field and protonation protocol, run explicit-solvent MD with
   independent replicates, inspect equilibration, and evaluate comparable
   equilibrated snapshots. Import replicate estimates with `method=mmgbsa`.
   Report convergence and sensitivity; endpoint estimates depend strongly on
   sampling, dielectric choices, force fields, and entropy treatment.
3. **Relative binding free energy (higher rigor for a congeneric series).** Use
   an alchemical workflow such as OpenFE, preserve a single reference compound
   per series, use independent repeats, and inspect overlap, uncertainty, and
   cycle closure. Import relative results with `method=rbfe` and populate
   `reference_compound`. RBFE values are ΔΔG values, not absolute ΔG values.
4. **Experimental affinity (decision standard).** Measure Kd when possible
   under a controlled assay protocol and import the replicates with
   `--affinity-results`. Experimental values remain separate from computed
   evidence in both tables and figures.

Useful method references are the original
[AutoDock Vina paper](https://doi.org/10.1002/jcc.21334), an
[MM/PBSA and MM/GBSA methods review](https://pmc.ncbi.nlm.nih.gov/articles/PMC4487606/),
and the current
[OpenFE RBFE protocol documentation](https://docs.openfree.energy/en/v1.10.0/tutorials/rbfe_cli_tutorial.html).

### Importing computed energies

Copy `energy_results.example.tsv`. Required columns are `compound`, `method`,
and `value_kcal_mol`.

| Column | Meaning |
|---|---|
| `compound` | Must match a discovered compound |
| `method` | Method label; methods are never pooled |
| `value_kcal_mol` | Finite method result in kcal/mol |
| `structure_id` | Optional; when present, must match `run_manifest.tsv` |
| `uncertainty_kcal_mol` | Optional nonnegative uncertainty |
| `reference_compound` | Required for `method=rbfe` |
| `replicate` | Replicate identifier used in duplicate validation |

An energy row may be compound-level by leaving `structure_id` blank. Duplicate
row identities are rejected. Every RBFE series must use one nonblank reference
compound. The summary provides population SD and a reproducible 95% bootstrap
confidence interval when at least two values are available.

### Importing experimental affinities

Copy `affinity_results.example.tsv`. Accepted metrics are `Kd`, `Ki`, and
`IC50`; accepted units are `M`, `mM`, `uM`, `nM`, and `pM`. Values and
temperatures must be positive. A blank temperature defaults to 298.15 K.

For a dissociation constant, the standard binding free energy is:

```text
ΔG° = R T ln(Kd / 1 M)
```

The implementation first converts Kd to molar concentration and uses
`R = 0.00198720425864083 kcal mol⁻¹ K⁻¹`. Kd is converted automatically. Ki is
converted only with `--convert-ki`, because its thermodynamic interpretation
depends on the inhibition model and assay. IC50 is never converted
automatically: it is assay-dependent and is not generally a dissociation
constant.

Example run with both evidence tables:

```bash
python compound_contact_pipeline.py /path/to/project \
  --manifest datasets.tsv \
  --energy-results energy_results.tsv \
  --affinity-results affinity_results.tsv \
  --out-dir compound_contact_results
```

## Outputs

All floating-point TSV cells use six decimal places. Missing numeric values are
blank instead of `nan`.

| Output | Contents |
|---|---|
| `run_manifest.tsv` | Resolved assignments, ligand-chemistry provenance, and original/analysis source paths |
| `per_structure_contacts.tsv` | One row per present protein residue in every valid structure |
| `structure_summary.tsv` | Selected analysis source, status, and structural metrics for every CIF |
| `compound_residue_summary.tsv` | Per-compound residue contact and distance statistics |
| `global_residue_summary.tsv` | Exploratory residue summary across compounds |
| `compound_summary.tsv` | Valid/invalid counts, pose consistency, SASA, and clash summaries |
| `minimization_summary.tsv` | Protocol, status/reason, output path, before/after geometry, and same-system potential energy |
| `energy_values_normalized.tsv` | Validated optional computed-energy values |
| `energy_summary.tsv` | Method- and reference-specific energy statistics |
| `affinity_values_normalized.tsv` | Measured values, molar normalization, and allowed ΔG° conversions |
| `affinity_summary.tsv` | Metric- and temperature-specific affinity statistics |
| `qc.tsv` | Structure and plotting errors/warnings with stage and identity |
| `generated_files.json` | Exact generated paths for narrowly scoped overwrite |

Figures below `plots/` are:

- `<compound>_chain_<chain>_contacts.png`: mean capped distance ± population
  SD and contact frequency;
- `contact_frequency_heatmap_chain_<chain>.png`: compound-by-residue contact
  frequency;
- `capped_distance_heatmap_chain_<chain>.png`: compound-by-residue capped
  distance;
- `compound_structural_summary.png`: valid seeds, contact-set consistency,
  buried SASA, and clashes, explicitly labeled as non-affinity evidence;
- `minimization_before_after.png`: paired fixed/van der Waals clash counts,
  backbone and ligand RMSD gates, and same-system potential-energy change;
- `energy_<method>.png` (or a reference-qualified variant): individual and
  summarized computed values, kept method-specific; and
- `experimental_affinity.png`: measured concentrations and separately derived
  standard ΔG° values.

Filenames are sanitized deterministically; sanitization collisions receive a
stable short hash.

## Safe overwrite behavior

A nonempty output directory is rejected unless `--overwrite` is supplied.
Every successful run records its generated files in `generated_files.json`.
On overwrite, the script validates the full preceding inventory and removes
only listed regular files or symlinks. It never recursively deletes directories
or removes unrelated files. A new output also cannot replace an unlisted file.

## AlphaFold 3 and predicted-pose caution

When CIF ensembles come from AlphaFold 3 or a related structure predictor,
repeatable contacts can be useful pose-prioritization evidence. They do not
demonstrate equilibrium occupancy, include the thermodynamic sampling needed
for binding free energy, or replace assay validation. The
[AlphaFold 3 paper](https://www.nature.com/articles/s41586-024-07487-w)
evaluates biomolecular structure prediction; this pipeline therefore treats
predicted-pose consistency as structural evidence only.
