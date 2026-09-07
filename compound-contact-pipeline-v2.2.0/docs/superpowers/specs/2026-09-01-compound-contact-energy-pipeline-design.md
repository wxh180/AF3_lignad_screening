# Compound Contact and Energy Aggregation Pipeline Design

**Date:** 2026-09-01

## Purpose

Build a portable command-line pipeline that analyzes protein–ligand contacts in
AlphaFold 3 or other mmCIF ensembles stored below arbitrarily nested compound
directories. The pipeline will aggregate results across seed structures,
compare compounds visually, and combine those structural results with
method-labeled energy or experimental-affinity data without representing a
contact-distance statistic as binding free energy.

## Scientific interpretation

The pipeline will maintain three distinct evidence classes:

1. **Structural evidence:** contacts, capped distances, buried solvent-accessible
   surface area, steric clashes, and contact-set consistency across seeds.
2. **Computed energy evidence:** Vina scores, MM/GBSA estimates, or alchemical
   relative binding free energies imported from the program that calculated
   them.
3. **Experimental evidence:** measured dissociation constants and their optional
   conversion to standard binding free energy.

Results from different energy methods will never be averaged together or
collapsed into a composite affinity score. Lower-is-better directionality and
units will be displayed for each method, and every plot and table will preserve
the method name.

The contact pipeline itself will not run molecular dynamics, prepare force-field
topologies, assign ligand protonation states, or parameterize compounds. Those
steps require chemical information that is not reliably recoverable from an
arbitrary prediction CIF. MM/GBSA or OpenFE calculations will be run separately
with a consistent preparation protocol, and their results will be ingested by
this pipeline.

## Deliverables

- `compound_contact_pipeline.py`: executable Python CLI with importable,
  testable functions.
- `datasets.example.tsv`: documented example dataset manifest.
- `energy_results.example.tsv`: documented example computed-energy input.
- `affinity_results.example.tsv`: documented example experimental-affinity
  input.
- `README.md`: installation, directory-layout, command, output, and scientific
  interpretation guidance.
- `tests/test_compound_contact_pipeline.py`: unit and end-to-end tests using
  synthetic mmCIF fixtures.

The single executable script is retained for easy transfer to a compute server;
its internal functions will remain separated by responsibility.

## Runtime dependencies

- Python 3.10 or newer
- Biopython
- NumPy
- Matplotlib
- pytest for the test suite only

No pandas, SciPy, molecular-dynamics engine, docking engine, or shell command is
required for contact analysis and aggregation.

## Input discovery

### Recursive discovery

The positional `ROOT` directory is searched recursively for directories whose
names match `seed-*`. Every `.cif` file below such a directory is a candidate
structure. Resolved file paths are deduplicated, and processing order is the
lexicographic order of paths relative to `ROOT` so repeated runs are
deterministic.

Only a directory's outermost matching `seed-*` ancestor contributes a file. This
prevents a CIF below nested `seed-*` directories from being discovered twice.

Each structure receives a collision-resistant `structure_id` equal to its POSIX
path relative to `ROOT`, including the `.cif` suffix. Two structures with the
same filename in different directories therefore remain distinct.

### Dataset manifest

When `--manifest datasets.tsv` is supplied, the manifest is authoritative. It
is tab-delimited and contains one row per dataset directory.

Required columns:

| Column | Meaning |
| --- | --- |
| `dataset_dir` | Directory relative to `ROOT`; only `seed-*` directories below it are included for this row. |
| `compound` | Human-readable compound identifier used in tables and plots. |
| `ligand_resname` | Exact, case-sensitive mmCIF ligand residue name. |

Optional columns:

| Column | Meaning | Default |
| --- | --- | --- |
| `protein_chains` | Comma-separated protein chain IDs. | All protein chains |
| `ligand_chain` | Ligand chain ID. | All chains containing the named ligand |
| `enabled` | `true` or `false`. | `true` |

Example:

```tsv
dataset_dir	compound	ligand_resname	protein_chains	ligand_chain	enabled
screen/C19	C19	WXI	A	B	true
screen/GSH	GSH	GSH	A,C	L	true
```

Manifest validation fails before CIF processing if a required value is empty,
a dataset directory is outside `ROOT`, an enabled dataset directory does not
exist, a compound identifier is duplicated with conflicting ligand settings,
or a discovered CIF matches more than one enabled dataset row.

### Manifest-free fallback

Without `--manifest`:

- The directory immediately above each `seed-*` directory supplies the compound
  name.
- `--ligand RESNAME` applies one ligand residue name to every compound.
- If `--ligand` is absent, ligand auto-detection is attempted independently for
  each CIF.

Auto-detection excludes amino-acid residues, water, and a built-in list of
common monatomic ions. It succeeds only when exactly one remaining residue name
is present. Zero or multiple candidates produce a QC failure; the pipeline does
not guess. All successfully auto-detected CIFs assigned to one compound must
yield the same residue name. A disagreement invalidates that compound and makes
the run exit nonzero rather than mixing different ligands.

## Structure parsing and selection

- Biopython `MMCIFParser(QUIET=True)` parses each file.
- Only the first model is analyzed. A CIF with additional models receives a QC
  warning that records the number ignored.
- Protein residues are residues accepted by `is_aa(..., standard=False)` on the
  selected protein chains.
- Water is always excluded.
- Ligand residues match both `ligand_resname` and, when supplied,
  `ligand_chain`.
- Protein-chain selection never filters ligand discovery.
- Multiple residues matching the same ligand definition are treated as one
  ligand group; minimum distances are taken over all matching ligand atoms.
- Disordered atoms use Biopython's selected conformer. The selected alternate
  location and occupancy behavior is recorded in the README rather than
  silently mixing conformers.

A structure is invalid for quantitative aggregation when parsing fails, the
ligand is absent, no ligand atoms remain, or no protein atoms remain. It is
retained in `qc.tsv` but excluded from valid-structure denominators.

## Contact calculations

For every valid structure, a `NeighborSearch` over protein atoms identifies
protein atoms within `--cutoff` of every ligand atom. The default cutoff is
4.5 Å and contact is defined inclusively as distance less than or equal to the
cutoff.

The stable residue key is:

```text
(chain_id, residue_number, insertion_code, residue_name)
```

For each protein residue present in a structure, output:

- `is_contact`: `1` when at least one atom pair is within the cutoff, otherwise
  `0`.
- `min_contact_distance_angstrom`: the minimum atom-pair distance when contact
  exists; blank for a non-contact residue.
- `capped_distance_angstrom`: the minimum contact distance when contact exists;
  otherwise the cutoff.

The capped value preserves the requested original visualization behavior while
making clear that a value equal to the cutoff is not the residue's measured
minimum distance.

### Steric clashes

The minimum ligand–protein atom distance and the number of atom pairs at or
below `--clash-cutoff` are recorded per structure. The default clash cutoff is
2.0 Å. This is a geometric warning, not a force-field energy.

### Buried surface area

Biopython `ShrakeRupley` calculates solvent-accessible surface area using a
1.4 Å probe for the isolated receptor, isolated ligand, and complex in the same
coordinates. The complex contains only the selected protein chains and selected
ligand group, matching the receptor and ligand used in the isolated
calculations. The following values are reported:

```text
buried_sasa_total_A2 = SASA_receptor + SASA_ligand - SASA_complex
interface_area_A2 = buried_sasa_total_A2 / 2
```

If SASA calculation fails, contact analysis remains valid, the SASA fields are
blank, and the failure is recorded as a QC warning.

## Aggregation rules

### Residue-level compound summary

For each `(compound, chain, residue key)`:

- `n_valid_structures`: structures with a valid protein–ligand analysis for the
  compound.
- `n_present`: valid structures in which the residue exists.
- `n_contacts`: present structures in which the residue contacts the ligand.
- `contact_frequency = n_contacts / n_present`.
- Mean, population SD, median, minimum, and maximum capped distance over the
  `n_present` structures.

A residue absent from one seed does not receive an artificial cutoff value and
does not enter that residue's denominator.

### Compound-level summary

For each compound, report:

- CIF files discovered, valid structures, invalid structures, and warning count.
- Mean and SD for the number of contacting residues per valid structure.
- Mean and SD for minimum protein–ligand distance, clash-pair count,
  `buried_sasa_total_A2`, and `interface_area_A2`.
- Contact-set consistency across seeds, calculated as the mean pairwise Jaccard
  similarity of residue contact sets. A compound with one valid seed has a blank
  Jaccard value and an explanatory QC note. For a seed pair, two empty contact
  sets have Jaccard similarity `1.0`; one empty and one nonempty set have
  similarity `0.0`.

No structural metric is renamed or interpreted as affinity.

### Global residue summary

A global table repeats residue statistics across all compounds for exploratory
inspection. Its frequencies are descriptive only because compounds may have
unequal numbers of valid seed structures.

## Energy-result integration

`--energy-results energy_results.tsv` accepts tab-delimited computed values.

Required columns:

| Column | Meaning |
| --- | --- |
| `compound` | Compound identifier matching the analysis dataset. |
| `method` | Method label such as `vina_score`, `mmgbsa`, or `rbfe`. |
| `value_kcal_mol` | Numeric score or energy in kcal/mol. |

Optional columns:

| Column | Meaning |
| --- | --- |
| `structure_id` | Relative CIF path when the value belongs to one seed or pose. |
| `uncertainty_kcal_mol` | Reported uncertainty for this value. |
| `reference_compound` | Reference compound for a relative free-energy value. |
| `replicate` | Replicate identifier. |

Validation rejects non-finite values, unknown compounds, unknown structure IDs,
negative uncertainty, and mixed reference compounds within one RBFE comparison
series. Duplicate rows with identical identifying fields are rejected rather
than silently double-counted. The identifying fields are `compound`, `method`,
`structure_id`, `reference_compound`, and `replicate`; blank optional fields
remain part of that identity. A blank `structure_id` denotes a compound-level
value and therefore is not checked against discovered structures.

Aggregation occurs within `(compound, method, reference_compound)` only. The
pipeline reports count, mean, population SD, median, minimum, maximum, and a
95% bootstrap confidence interval when at least two independent values exist.
Imported per-value uncertainty is retained and plotted but is not combined with
between-replicate variability. Bootstrap intervals use 10,000 resamples and a
fixed random seed of `20260901` so repeated runs are identical.

Vina values are labeled `Vina score (kcal/mol)`, not binding free energy.
MM/GBSA values are labeled `MM/GBSA estimate (kcal/mol)`. RBFE values are labeled
`relative binding free energy ΔΔG (kcal/mol)` and display the reference compound.

## Experimental-affinity integration

`--affinity-results affinity_results.tsv` accepts:

| Column | Requirement |
| --- | --- |
| `compound` | Required; must match a dataset compound. |
| `metric` | Required; `Kd`, `Ki`, or `IC50`. |
| `value` | Required positive finite number. |
| `unit` | Required; `M`, `mM`, `uM`, `nM`, or `pM`. |
| `temperature_K` | Optional; defaults to 298.15 K. |
| `replicate` | Optional replicate identifier. |

For `Kd`, the pipeline calculates:

```text
delta_G_standard_kcal_mol = R * T * ln(Kd / 1 M)
```

using `R = 0.00198720425864083 kcal mol^-1 K^-1`. `Ki` is reported and can be
converted only when `--convert-ki` is explicitly supplied. `IC50` is never
converted automatically because the conversion depends on assay conditions and
mechanism.

Concentrations are normalized to molar units. Experimental data are summarized
within `(compound, metric, temperature_K)` and plotted separately from computed
values. Free energies are calculated per replicate before summary statistics
are computed.

## Outputs

All outputs are written below `--out-dir`, which defaults to
`compound_contact_results`. Without `--overwrite`, a nonempty output directory
is rejected. Each successful run writes a machine-readable inventory of the
files it generated. With `--overwrite`, only paths listed in the preceding
inventory are replaced or removed after containment checks; unrelated files in
the output directory are never removed.

### Tables

- `run_manifest.tsv`: resolved dataset and CIF assignments.
- `per_structure_contacts.tsv`: one row per present protein residue per valid
  structure.
- `structure_summary.tsv`: one row per discovered CIF, including structural
  metrics and processing status.
- `compound_residue_summary.tsv`: residue statistics within each compound.
- `global_residue_summary.tsv`: exploratory residue statistics across compounds.
- `compound_summary.tsv`: compound structural metrics and pose consistency.
- `energy_values_normalized.tsv` and `energy_summary.tsv` when computed energies
  are supplied.
- `affinity_values_normalized.tsv` and `affinity_summary.tsv` when experimental
  affinities are supplied.
- `qc.tsv`: errors and warnings with dataset, compound, structure ID, stage, and
  message.
- `generated_files.json`: exact relative paths generated by the completed run,
  used to make a later `--overwrite` operation narrowly scoped.

Floating-point TSV values use six decimal places; missing numeric values are
empty fields rather than strings such as `nan`.

### Figures

- `plots/<compound>_chain_<chain>_contacts.png`: mean capped distance ± SD and
  contact frequency, with labels for frequencies at or above
  `--label-frequency` (default 0.5).
- `plots/contact_frequency_heatmap_chain_<chain>.png`: compound-by-residue
  contact-frequency heatmap; a residue absent from all valid structures for one
  compound is gray rather than displayed as zero frequency.
- `plots/capped_distance_heatmap_chain_<chain>.png`: compound-by-residue mean
  capped-distance heatmap; absent compound–residue combinations are gray.
- `plots/compound_structural_summary.png`: valid-seed count, contact-set Jaccard,
  buried SASA, and clash summary without an affinity ranking.
- `plots/energy_<method>.png`: individual values plus method-specific mean and
  SD/CI when computed energies are supplied.
- `plots/experimental_affinity.png`: measured affinity and converted standard
  binding free energy when experimental data are supplied.

Filenames are sanitized deterministically. A filename collision after
sanitization receives a short stable hash suffix.

## Command-line interface

Primary usage:

```bash
python compound_contact_pipeline.py ROOT \
  --manifest datasets.tsv \
  --cutoff 4.5 \
  --clash-cutoff 2.0 \
  --energy-results energy_results.tsv \
  --affinity-results affinity_results.tsv \
  --out-dir compound_contact_results
```

Manifest-free usage:

```bash
python compound_contact_pipeline.py ROOT \
  --ligand WXI \
  --protein-chain A \
  --out-dir compound_contact_results
```

Additional options:

- `--protein-chain` may be repeated and applies only in manifest-free mode.
- `--ligand-chain` applies only in manifest-free mode.
- `--workers N` uses deterministic process-based parallel analysis; default `1`.
- `--label-frequency FLOAT` controls residue annotations; default `0.5`.
- `--convert-ki` enables the explicitly qualified Ki-to-ΔG conversion.
- `--fail-fast` stops at the first per-structure error; default behavior records
  the error and continues.
- `--overwrite` permits replacement of files in an existing output directory.

The CLI exits nonzero for invalid configuration, invalid auxiliary tables, no
CIF files, no valid structures, or an existing non-empty output directory
without `--overwrite`. Individual CIF failures do not produce a nonzero exit if
at least one valid structure remains, unless `--fail-fast` is active.

## Error handling and logging

- User-facing progress and warnings go to standard error.
- Standard output contains a concise final list of generated tables and figures.
- Every recoverable structure-level problem is also written to `qc.tsv`.
- Configuration and schema errors identify the file, row, column, and invalid
  value when applicable.
- Plotting failure for one chain or optional energy method is recorded and does
  not discard completed tables; any plotting failure makes the overall exit
  status nonzero so automation can detect incomplete output.

## Testing strategy

Tests will be written before implementation and will cover:

1. Recursive `seed-*` discovery, nested-seed deduplication, deterministic order,
   and duplicate CIF basenames.
2. Manifest precedence, path containment, overlapping dataset rejection, chain
   parsing, and manifest-free compound assignment.
3. Ligand auto-detection success and ambiguous-candidate failure.
4. Separation of protein-chain and ligand-chain filters.
5. Inclusive cutoff behavior, minimum contact distance, capped non-contact
   distance, insertion codes, and missing-residue denominators.
6. Missing ligand, missing protein, malformed CIF, multiple models, and
   continue-versus-fail-fast behavior.
7. Clash counting, buried-SASA identities, and SASA warning fallback.
8. Compound summaries and pairwise contact-set Jaccard similarity.
9. Energy schema validation, method isolation, structure-ID validation,
   reference-compound validation, and bootstrap reproducibility with a fixed
   random seed.
10. Unit conversion, Kd-to-ΔG calculation, opt-in Ki conversion, and refusal to
    convert IC50.
11. TSV formatting, sanitized filenames, plot smoke tests, overwrite
    protection, and an end-to-end CLI run.

Synthetic CIF fixtures keep the test suite independent of proprietary compound
structures and external programs.

## Recommended physical-energy workflow outside this script

### Medium-cost relative ranking: ensemble MM/GBSA

Use one consistent receptor preparation, ligand protonation protocol, force
field, charge model, solvent model, minimization/MD schedule, and snapshot count
for every compound. Compute single-trajectory MM/GBSA over multiple sampled
snapshots and independent repeats, then import replicate estimates and
uncertainties. Interpret differences primarily as relative ranking and validate
against known experimental data.

Amber provides an official `MMPBSA.py` tutorial:
<https://ambermd.org/tutorials/advanced/tutorial3/py_script/section1.htm>.

### Higher-rigor comparison for related compounds: RBFE

For congeneric compounds with a defensible atom mapping, use an explicit-solvent
relative binding free-energy protocol with independent repeats and a connected
transformation network. Import each ΔΔG with its reference compound and
uncertainty.

OpenFE documents its OpenMM relative free-energy protocol and repeat-based
uncertainty model:
<https://docs.openfree.energy/en/stable/reference/api/openmm_rfe.html>.

### Rapid empirical pre-screen: Vina score-only

Vina can score pre-positioned, properly protonated PDBQT poses quickly. Use the
same receptor preparation and scoring box for every compound. Treat these
numbers as empirical pose scores and do not relabel them as experimental or
alchemical binding free energies.

AutoDock Vina documentation:
<https://autodock-vina.readthedocs.io/en/stable/vina.html>.

## Acceptance criteria

The implementation is acceptable when all tests pass and a synthetic example
run demonstrates that it:

- discovers arbitrary nested compound datasets;
- preserves unique structure identity across repeated filenames;
- produces correct residue and compound denominators in the presence of missing
  residues and invalid structures;
- writes every specified core table and contact figure;
- validates and separates energy methods;
- converts Kd correctly without converting IC50;
- reports structural metrics without calling them affinity; and
- completes without external docking or MD software when no energy or affinity
  tables are supplied.
