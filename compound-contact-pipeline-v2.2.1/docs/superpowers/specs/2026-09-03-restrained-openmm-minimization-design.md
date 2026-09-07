# Restrained OpenMM Minimization Design

**Date:** 2026-09-03  
**Target release:** 2.0.0

## Purpose

Extend the compound contact pipeline with an optional, reproducible OpenMM
minimization stage that relieves local protein-ligand steric clashes without
silently changing the predicted binding mode. The extension must preserve all
input CIF files, retain the current analysis as the default behavior, and make
every minimized structure auditable through before/after geometry and energy
quality-control measurements.

This is pose preparation, not affinity prediction. Potential energy from a
single minimized complex is not a binding free energy and will never be placed
in the pipeline's imported binding-energy tables or used to rank compounds.

## Approaches considered

### 1. Restrained OpenMM minimization (selected)

Use an Amber protein force field and a small-molecule template generator,
followed by staged L-BFGS minimization. This is portable across AF3 and docking
ensembles, can add hydrogens consistently, provides explicit energies and
coordinates, and supports quantitative movement limits.

The cost is a separate, heavier optional environment and a requirement for the
ligand's exact chemical identity. These constraints are preferable to guessing
bond orders, protonation, charge, or stereochemistry from Cartesian coordinates.

### 2. Phenix geometry or real-space refinement

Phenix is preferable when a cryo-EM map and ligand restraint CIF are available,
because map-aware refinement can retain experimental support while improving
geometry. It is not selected as the batch default because most AF3 or docking
ensembles have no map, Phenix is separately licensed, and its inputs and
acceptance criteria differ from force-field minimization. The README will keep
this as the recommended cryo-EM-specific follow-up.

### 3. Generic external-command hook

A command template could support arbitrary minimizers, but it would make input
contracts, error handling, coordinate identity, and QC difficult to guarantee.
It is deferred until a concrete second backend is required.

## Scope

Version 2.0 will add:

- optional `--minimize openmm` execution;
- exact ligand chemistry fields in manifest and manifest-free modes;
- restrained, staged local minimization;
- accepted and rejected minimized mmCIF outputs without modifying originals;
- element-aware clash metrics in addition to the existing fixed-distance
  severe-clash count;
- before/after potential energy, displacement, clash, and contact QC;
- an explicit choice of original or accepted minimized coordinates for the
  existing aggregation pipeline;
- a before/after minimization summary plot; and
- an optional OpenMM environment specification and setup documentation.

Version 2.0 will not add molecular dynamics, solvent equilibration, MM/GBSA,
RBFE, docking, automatic protonation-state enumeration, missing-heavy-atom
modeling, covalent-ligand handling, map refinement, or affinity ranking.

## Backward compatibility and activation

Without `--minimize`, all existing commands, required dependencies, output
names, and analysis semantics remain available. OpenMM-related packages are
imported lazily only after minimization is requested. Hydrogen-free inputs,
including normal AF3 outputs, retain the existing numerical definitions.

The new command pattern is:

```bash
python compound_contact_pipeline.py ROOT \
  --manifest datasets.tsv \
  --minimize openmm \
  --analysis-source minimized \
  --out-dir compound_contact_results
```

`--analysis-source` accepts `original` or `minimized` and defaults to
`original`. Selecting `minimized` requires `--minimize openmm`. This explicit
choice prevents an optional preparation step from silently changing the
coordinates underlying established contact summaries.

The release will expose `--version` and report `2.0.0`.

## Ligand chemistry contract

OpenMM parameterization must not infer chemical identity from ligand atom
positions. Each enabled compound must provide exactly one of:

- `ligand_smiles`: an isomeric SMILES describing the exact protonation,
  tautomeric state, formal charge, and stereochemistry; or
- `ligand_file`: a root-relative one-molecule SDF containing the exact ligand.

These become optional manifest columns but are required when minimization is
requested. Manifest-free minimization uses mutually exclusive
`--ligand-smiles` and `--ligand-file`; a single global chemistry definition is
permitted only when all discovered jobs belong to one compound.

Ligand files must resolve inside `ROOT`. A file containing zero or multiple
molecules, undefined stereochemistry, non-finite coordinates, or chemistry
that does not match the selected mmCIF ligand topology is a configuration or
structure-level error as appropriate. The pipeline will not repair or guess a
mismatch.

The OpenFF molecule must contain all hydrogens in its chemical graph. The input
mmCIF may omit hydrogens; OpenMM `Modeller.addHydrogens` will add them after the
protein and ligand templates are registered. The selected ligand topology must
match the supplied molecule by element and bond pattern; atom ordering need not
match. A ligand whose mmCIF lacks the bonds needed for graph matching fails
with guidance to regenerate the mmCIF with explicit `_struct_conn` records; no
distance-based bond inference is attempted. Version 2.0 supports one selected
ligand residue per structure. Multiple copies are rejected with an actionable
QC message instead of being combined.

## Force field and physical model

The pinned default protocol is:

- protein: `amber/protein.ff14SB.xml`;
- ligand: `gaff-2.2.20` through `openmmforcefields`;
- ligand charges: user-supplied SDF partial charges when present, otherwise
  the template generator's documented AM1-BCC route;
- nonbonded treatment: nonperiodic `NoCutoff`;
- constraints: bonds involving hydrogen;
- pH for standard protein residue hydrogen placement: 7.4;
- minimizer: OpenMM L-BFGS through `LocalEnergyMinimizer`;
- convergence tolerance: 10 kJ mol^-1 nm^-1; and
- platform: CPU by default, with explicit `CPU`, `CUDA`, `OpenCL`, or
  `Reference` selection.

Force-field names, package versions, platform, device index, pH, tolerance,
iteration limits, and restraint constants are written to the minimization
summary. Defaults remain pinned within a protocol version so results do not
change merely because a newer force field becomes installed.

The complex used for minimization contains only the selected protein chains
and the single selected ligand. Waters, ions, unselected chains, and other
heterogens are excluded, matching the scope of the contact calculation. A
structure with missing protein heavy atoms or any unparameterized retained
residue fails clearly; the pipeline does not build missing heavy atoms.

## Restraint protocol

The pocket is the set of protein residues with any heavy atom within 6.0 A of
any ligand heavy atom in the original pose. Protein backbone atoms are `N`,
`CA`, `C`, and `O`. Harmonic positional restraints use each atom's original
position and OpenMM `CustomExternalForce`. Hydrogens are never position
restrained.

All force constants below use the energy expression `k*r^2`:

1. **Hydrogen relaxation, 500 iterations.** Restrain all original heavy atoms
   at 10 kcal mol^-1 A^-2. This relaxes newly added hydrogens before heavy atoms
   can respond.
2. **Pocket relaxation, 1,000 iterations.** Restrain backbone and non-pocket
   protein heavy atoms at 10 kcal mol^-1 A^-2; leave pocket side-chain heavy
   atoms free; restrain ligand heavy atoms at 1 kcal mol^-1 A^-2.
3. **Gentle relaxation, 1,000 iterations.** Reduce backbone and non-pocket
   protein restraints to 2 kcal mol^-1 A^-2 and ligand restraints to
   0.2 kcal mol^-1 A^-2; keep pocket side chains free.

Each stage has a bounded iteration count. A failure or non-finite coordinate or
energy stops that structure's minimization and becomes a QC error. Version 2.0
runs minimizations sequentially; `--minimize openmm` therefore requires
`--workers 1`. This avoids parameter-cache races and accidental competition
between multiple GPU contexts. Contact analysis retains its existing parallel
mode when minimization is off.

## Geometry and movement QC

Hydrogens added internally for parameterization and minimization are transient:
the published minimized CIF contains only the selected protein and ligand atoms
that existed in the input, with their original atom names and residue identity.
This keeps the existing contact and inclusive `--clash-cutoff` calculations
comparable before and after minimization. RMSD and the new element-aware clash
metric use heavy atoms only, so transient proton placement cannot change them.

An additional element-aware overlap metric uses a documented Bondi-style van
der Waals radius table. For every noncovalent protein-ligand heavy-atom pair:

```text
overlap_A = radius_protein + radius_ligand - distance_A
```

Pairs with overlap at least 0.4 A are counted as van der Waals clashes. The
maximum overlap is also reported. Unsupported elements cause a QC error rather
than silently receiving a generic radius. Covalently attached ligands are out
of scope.

Protein backbone atoms are least-squares aligned between original and final
coordinates. The pipeline reports:

- aligned protein-backbone RMSD;
- ligand-heavy-atom RMSD after applying that same protein alignment;
- original and final fixed-cutoff clash-pair counts;
- original and final van der Waals clash-pair counts;
- original and final maximum van der Waals overlap;
- original and final minimum protein-ligand heavy-atom distance;
- original and final contacting-residue counts and contact-set Jaccard;
- initial and final force-field potential energy; and
- potential-energy change for this same parameterized system.

Restraint energy is assigned to a separate OpenMM force group and excluded from
the reported physical-system potential energy. Energy differences remain
within-structure optimization diagnostics only.

## Acceptance and rejection

A completed result is accepted only when all of the following hold:

- coordinates and reported energies are finite;
- atom identities required for QC are preserved;
- aligned protein-backbone RMSD is at most 0.5 A;
- aligned ligand-heavy-atom RMSD is at most 1.5 A;
- the fixed-cutoff clash-pair count does not increase; and
- the van der Waals clash-pair count does not increase.

The RMSD thresholds are configurable with `--max-backbone-rmsd` and
`--max-ligand-rmsd`. The pocket radius is configurable with
`--pocket-radius`. Force constants and stage iteration counts remain part of
the named protocol instead of becoming a large collection of loosely tested
CLI knobs.

Acceptance does not assert that a binding mode is correct. A pose may pass
because it required little movement yet still be unsupported biologically. A
large original ligand-backbone penetration should be redocked or rebuilt;
minimization rejection is a safeguard, not a request to weaken the gates.

## Data flow

1. Discover jobs and validate manifest/chemistry settings.
2. Analyze every original CIF using heavy-atom contact and clash definitions.
3. For each originally valid job, build the selected OpenMM complex, add
   hydrogens, parameterize, minimize through three stages, and calculate QC.
4. Write every completed candidate to a run-scoped temporary directory.
5. Classify candidates as accepted or rejected from objective QC gates.
6. If `--analysis-source original`, aggregate the original results. If it is
   `minimized`, reanalyze only accepted candidates under their original
   `structure_id`; failed or rejected jobs become invalid post-minimization
   results so denominators remain explicit.
7. Prepare the protected output directory, then publish tables, figures, and
   minimized CIFs under the generated-file inventory.

Temporary files are removed after publication. Input files are opened
read-only and never replaced.

## Outputs

When minimization is requested, add:

- `minimization_summary.tsv`: one row per discovered structure with protocol,
  status, reason, before/after QC, energies, and relative output path;
- `plots/minimization_before_after.png`: paired clash, overlap, RMSD, and
  potential-energy-change panels, grouped by compound;
- `minimized_structures/accepted/<original-parent>/model_minimized.cif`; and
- `minimized_structures/rejected/<original-parent>/model_minimized.cif` for a
  completed result that fails an acceptance gate.

No CIF is written for setup or parameterization failure. The summary and
`qc.tsv` retain the failure reason. Accepted and rejected paths are included in
`generated_files.json` and pass the existing output-containment checks.

Core `run_manifest.tsv` gains `analysis_source`, `original_cif_path`, and
`analysis_cif_path`. Existing `structure_id` values remain the original
root-relative identifiers even when accepted minimized coordinates are used.

`structure_summary.tsv` records the selected `analysis_source`. Contact,
residue, compound, energy, and affinity tables otherwise retain their current
roles. OpenMM potential energies never enter `energy_values_normalized.tsv` or
`energy_summary.tsv`.

## Error handling

Configuration errors fail before minimization, including missing optional
dependencies, invalid force-field or platform names, missing ligand chemistry,
and conflicting CLI options. Per-structure topology, parameterization,
minimization, acceptance, and write failures are recorded in `qc.tsv` and
`minimization_summary.tsv` and processing continues unless `--fail-fast` is
active.

If `--analysis-source minimized` yields no accepted valid structures, the run
returns status 2. If original analysis is selected, minimization failures do
not erase valid original structural analysis, but they remain visible and make
the run return status 1 to signal incomplete optional preparation.

## Optional dependencies and reproducibility

Base contact analysis keeps its current lightweight requirements. A new
`environment-openmm.yml` documents a conda-forge environment containing Python,
OpenMM, openmmforcefields, OpenFF Toolkit, RDKit, AmberTools, NumPy,
Matplotlib, Biopython, and pytest.

The code checks optional imports only when `--minimize openmm` is selected and
emits a single installation command on failure. The minimization table records
resolved package versions. A parameter cache is shared across compounds only
inside the run-scoped temporary area and is deleted at the end of the run.

## Testing strategy

Tests are written before production changes and cover:

1. ligand chemistry option and manifest validation, path containment, and
   compound consistency;
2. lazy optional-dependency failure without affecting the base CLI;
3. heavy-atom selection and element-aware overlap calculations at boundary
   values;
4. Kabsch backbone alignment and ligand RMSD under translation/rotation;
5. pocket and restraint-group assignment;
6. staged protocol constants and kcal/A^2 to kJ/nm^2 conversion;
7. objective acceptance and every rejection reason;
8. deterministic minimized output paths and preservation of structure IDs;
9. original-versus-minimized analysis-source behavior and denominator rules;
10. summary serialization, QC records, plot creation, and inventory safety;
11. continuation and fail-fast behavior after a per-structure minimization
    error; and
12. an optional integration test, skipped when the OpenMM environment is not
    installed, that minimizes a tiny complete protein-ligand fixture and
    verifies finite energy plus reduced overlap.

Tests that do not require OpenMM use a dependency-injected fake backend and
real coordinate/QC calculations. The optional integration test exercises the
actual third-party stack when available.

## Documentation

The README will include:

- a base-install and separate OpenMM-environment install path;
- manifest and manifest-free minimization examples;
- exact ligand chemistry and protonation requirements;
- protocol constants and acceptance criteria;
- interpretation of every before/after field;
- a warning that minimization and complex potential energy are not binding
  affinity calculations; and
- a cryo-EM note recommending Phenix real-space refinement with the map and
  ligand restraints when experimental density is available.

## Implementation references

- [OpenMM `LocalEnergyMinimizer`](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.LocalEnergyMinimizer.html)
  documents L-BFGS minimization, convergence tolerance, and bounded iterations.
- The [OpenMM positional-restraint cookbook](https://openmm.github.io/openmm-cookbook/latest/notebooks/cookbook/Restraining%20Atom%20Positions.html)
  demonstrates harmonic restraints with `CustomExternalForce`.
- The [openmmforcefields documentation](https://github.com/openmm/openmmforcefields/blob/main/README.md)
  documents ff14SB, GAFF 2.2.20, exact small-molecule identity requirements,
  template matching, charges, and parameter caching.
- [Phenix real-space refinement](https://phenix-online.org/documentation/reference/real_space_refine)
  remains the map-aware route for experimental cryo-EM models.

## Acceptance criteria

The release is complete when:

- all prior tests still pass with minimization disabled;
- the new test suite passes with warnings treated as errors;
- `--help` and `--version` describe the 2.0 interface;
- a base environment can still run contact analysis without importing OpenMM;
- an OpenMM environment can produce accepted/rejected minimized CIFs and a
  complete before/after summary from a synthetic example;
- originals remain byte-for-byte unchanged;
- minimized analysis never includes a rejected or failed candidate; and
- documentation never labels complex potential energy, clash reduction, or
  minimization as binding affinity.
