# Compound Contact Pipeline v2.2.0 Performance and Progress Design

Date: 2026-09-07
Baseline: compound-contact-pipeline v2.1.0
Status: Approved design

## 1. Purpose

Version 2.2.0 will reduce unnecessary CPU work during restrained OpenMM workflows and make long-running jobs visibly traceable without changing the scientific definitions of contacts, clash metrics, buried SASA, interface area, or the existing three-stage minimization protocol.

The main performance problem in v2.1.0 is that `analyze_structure()` always computes buried SASA. When minimization is requested, every original structure is fully analyzed before OpenMM begins, and every minimized candidate is fully analyzed again before acceptance. The acceptance decision does not require SASA. For `--analysis-source minimized`, this means expensive Shrake-Rupley calculations are performed on original structures that are not used for the final scientific summaries and on minimized candidates that may later be rejected.

A second limitation is that v2.1.0 rejects `--minimize` when `--workers` is greater than 1, even though CPU-side structure analysis can be parallelized independently from sequential single-GPU minimization.

A third limitation is observability: the command can remain at `Analyzing N CIF structure(s)...` for a long time with no indication of which CPU or GPU stage is active.

## 2. Goals

v2.2.0 will:

1. Automatically skip SASA when a structure is being analyzed only for minimization eligibility or acceptance QC.
2. Compute full SASA only for structures that contribute to the selected scientific analysis source.
3. Allow `--workers N` to parallelize CPU analysis while keeping OpenMM minimization sequential on the selected GPU.
4. Provide default progress and stage timing output to stderr.
5. Add `--quiet` to suppress normal progress while preserving warnings and errors.
6. Preserve deterministic result ordering and the existing output table schemas.
7. Preserve the v2.1.0 restrained OpenMM scientific protocol and acceptance criteria.

## 3. Non-goals

v2.2.0 will not add:

- concurrent minimization jobs on one GPU;
- multi-GPU scheduling;
- asynchronous CPU/GPU streaming;
- resume/checkpoint support;
- a persistent performance/profiling results table;
- new scientific metrics;
- changes to force fields, minimization stages, restraint strengths, tolerances, pH handling, clash definitions, RMSD gates, or acceptance criteria.

These are intentionally deferred so v2.2.0 remains a focused performance and observability release.

## 4. Existing v2.1.0 Behavior

The current `analyze_structure()` performs parsing, ligand/protein selection, contact analysis, fixed-cutoff clash analysis, and then `calculate_buried_sasa()`. `calculate_buried_sasa()` builds receptor, ligand, and complex entities and invokes Biopython `ShrakeRupley` three times.

The current `_run_pipeline()` always calls `_analyze_jobs()` on all original CIFs before minimization. `run_minimization_jobs()` then minimizes each structure sequentially and calls `analyze_structure()` again on each minimized candidate before acceptance evaluation.

The current argument validation rejects `--minimize` unless `--workers 1` is used.

The existing OpenMM backend creates a `SystemGenerator`, adds hydrogens, creates the force-field system and `Context`, and then performs the three minimization stages:

1. `hydrogen_relaxation`
2. `pocket_relaxation`
3. `gentle_relaxation`

Those stages and their current parameters remain unchanged in v2.2.0.

## 5. Chosen Architecture

v2.2.0 will use a phase-based architecture rather than streaming or concurrent GPU execution.

### Phase 1: CPU analysis

Analyze original structures in parallel using up to `--workers N` processes. Whether SASA is enabled depends on the selected analysis source.

### Phase 2: sequential OpenMM minimization

Minimize structures one at a time on the selected OpenMM platform/device. Candidate structural QC is performed without SASA.

### Phase 3: final CPU full analysis

When `--analysis-source minimized` is selected, perform full analysis including SASA only for accepted minimized structures, using up to `--workers N` processes.

Aggregation, plotting, and table generation then proceed from the selected full scientific results.

This architecture is preferred because it is deterministic, easy to debug, compatible with the current control flow, and avoids CPU contention with OpenMM setup during GPU minimization.

## 6. Internal Analysis Modes

`analyze_structure()` will gain an internal parameter:

```python
analyze_structure(job, cutoff=4.5, clash_cutoff=2.0, compute_sasa=True)
```

`compute_sasa=True` remains the default so existing direct callers preserve v2.1.0 behavior.

### Fast geometry analysis (`compute_sasa=False`)

The function still performs:

- mmCIF parsing and first-model selection;
- ligand detection and ligand-chain selection;
- protein-chain selection;
- contact residue identification;
- per-residue minimum ligand distance;
- global minimum protein-ligand atom distance;
- fixed-cutoff clash pair count;
- all existing QC associated with parsing and selection.

It does not call `calculate_buried_sasa()`. `buried_sasa_total_A2` and `interface_area_A2` remain `None` for this intermediate result.

### Full scientific analysis (`compute_sasa=True`)

The function performs all fast geometry operations plus the existing buried-SASA/interface-area calculation with unchanged scientific definitions.

No public `--sasa-mode` option will be added. SASA scheduling is automatic.

## 7. Exact Execution Rules

### 7.1 No minimization

For a normal analysis run without `--minimize`:

```text
original CIFs -> FULL analysis including SASA -> aggregation/reports
```

CPU analysis uses up to `--workers N` processes.

### 7.2 `--minimize openmm --analysis-source original`

Execution order:

```text
original CIFs
    -> FULL analysis including SASA, parallel CPU
    -> sequential OpenMM minimization
    -> FAST minimized-candidate QC, no SASA
    -> acceptance/rejection
    -> reports from original FULL results
```

The minimized candidate exists only for minimization QC and therefore does not receive SASA.

### 7.3 `--minimize openmm --analysis-source minimized`

Execution order:

```text
original CIFs
    -> FAST analysis, no SASA, parallel CPU
    -> sequential OpenMM minimization
    -> FAST candidate QC, no SASA
    -> acceptance/rejection
    -> accepted minimized CIFs only
    -> FULL analysis including SASA, parallel CPU
    -> reports from accepted minimized FULL results
```

Rejected and failed minimized candidates never receive SASA.

The final result list must preserve one position per original input job. Accepted structures are represented by their final full minimized result. Rejected, failed, or not-attempted structures remain invalid results with the same minimization-stage semantics used by v2.1.0.

## 8. Worker Semantics

`--workers N` will mean CPU analysis concurrency only.

For example:

```bash
--minimize openmm --openmm-platform CUDA --openmm-device-index 0 --workers 16
```

means:

```text
CPU analysis workers:      up to 16
OpenMM minimization jobs:  1 at a time
GPU concurrency:           1
```

The validation rule that currently rejects `--minimize` with `--workers != 1` will be removed.

The existing rule that `--fail-fast` requires `--workers 1` will remain. This avoids introducing nondeterministic fail-fast behavior in multiprocessing.

CPU multiprocessing will retain deterministic result ordering. Completion/progress messages may appear in completion order, but returned `StructureResult` objects must align with input job order.

## 9. Minimization and Acceptance Data Flow

`run_minimization_jobs()` will continue to operate sequentially. For each valid original job it will:

1. construct the existing `MinimizationRequest`;
2. run the existing OpenMM backend;
3. build a candidate `StructureJob` from the staged minimized CIF;
4. call fast candidate analysis with `compute_sasa=False`;
5. evaluate acceptance using the same v2.1.0 inputs:
   - backbone RMSD;
   - ligand RMSD;
   - before/after fixed clash counts;
   - before/after VDW clash counts;
   - current configured thresholds;
6. compute contact-set Jaccard exactly as before;
7. store the fast candidate result in the minimization record for QC/provenance.

For `--analysis-source minimized`, accepted candidate paths are then passed through a separate full-analysis phase. The final scientific result should not simply reuse the fast `candidate_result`, because its SASA fields are intentionally unset.

## 10. Progress and Timing Output

Progress is enabled by default and written to stderr so stdout can remain usable for generated-file paths and machine-readable workflows.

### 10.1 CPU phase messages

A phase header will state the number of structures, worker count, and whether SASA is enabled. Example:

```text
=== Phase 1/3: Pre-minimization analysis ===
25 structures | 16 CPU workers | SASA disabled

[ 1/25] structure_A completed    0.72 s
[ 2/25] structure_B completed    0.75 s
...
Phase completed: 2.8 s
```

When multiprocessing is used, per-structure completion messages may reflect completion order. This does not alter result ordering.

### 10.2 OpenMM messages

Before each potentially long operation, a message is emitted. The OpenMM backend will expose lightweight progress/timing events without changing the scientific calculation.

Representative output:

```text
=== Phase 2/3: Restrained OpenMM minimization ===
25 structures | CUDA device 0 | sequential

[1/25] structure_A
    OpenMM preparation...
    preparation completed ........ 5.84 s
    hydrogen_relaxation [CUDA] ...
    hydrogen_relaxation .......... 0.91 s
    pocket_relaxation [CUDA] .....
    pocket_relaxation ............ 1.63 s
    gentle_relaxation [CUDA] .....
    gentle_relaxation ............ 1.11 s
    candidate geometry QC ........ 0.18 s
    ACCEPTED
    total ........................ 9.93 s

Batch: 1/25 | elapsed 00:00:10 | ETA 00:04:00
```

The exact labels should make it possible to distinguish CPU preparation from active CUDA minimization. In particular, `OpenMM preparation...` is printed before system preparation, while each minimization-stage label is printed immediately before `LocalEnergyMinimizer.minimize()`.

### 10.3 ETA

Batch ETA is based on completed structures and elapsed wall time. It is informational only and must never affect control flow. ETA may be omitted until at least one item has completed.

### 10.4 Quiet mode

A new `--quiet` flag suppresses normal phase, progress, and timing messages. Warnings, validation errors, minimization failures, and other diagnostic errors remain visible on stderr.

Progress reporting must be best-effort and side-effect free. A formatting or timing failure must not alter scientific calculations or acceptance outcomes.

## 11. OpenMM Instrumentation Boundary

The OpenMM backend will gain a minimal optional progress callback or equivalent reporting hook that can emit events around:

- ligand/topology preparation;
- hydrogen addition and force-field/system construction;
- Context creation as part of preparation;
- `hydrogen_relaxation`;
- `pocket_relaxation`;
- `gentle_relaxation`;
- final geometry/output completion.

The callback must not be required for backend use and must have a no-op default so existing tests and direct callers remain valid.

Instrumentation will not change:

- `MinimizationConfig` scientific defaults;
- force fields;
- `SystemGenerator` chemistry behavior;
- restraint force construction;
- stage iteration limits;
- minimization tolerance;
- pH handling;
- RMSD calculations;
- VDW overlap calculations;
- output coordinate handling.

## 12. Output Compatibility

Existing scientific TSV columns and plot definitions will remain unchanged.

Intermediate fast results may contain `None` for buried SASA and interface area, but these intermediate results must not be passed into final scientific aggregation when the selected analysis source requires full results.

No runtime columns will be added to the scientific summary tables in v2.2.0. Timing information is operational metadata and will remain console-only in this release.

The existing minimization manifest/provenance behavior remains intact.

## 13. Error Handling

Existing invalid-structure and minimization failure semantics remain in place.

Additional requirements:

- A fast analysis failure is handled exactly like the corresponding v2.1.0 structural analysis failure except SASA-specific warnings cannot occur when SASA is disabled.
- Failure during final full analysis of an accepted minimized structure produces an invalid final scientific result and an appropriate QC record; it must not silently fall back to the fast result.
- OpenMM progress callbacks must not be allowed to turn an otherwise successful minimization into a scientific failure.
- `--quiet` affects only normal informational output.
- `--fail-fast` continues to require `--workers 1` and retains current failure semantics.

## 14. Determinism

The v2.2.0 optimization changes when calculations occur, not their scientific definitions.

For equivalent valid structures, full analysis in v2.2.0 must use the same contact, clash, and SASA code paths as v2.1.0. Parallel CPU execution must return results aligned with the original sorted job list.

The sequential OpenMM ordering remains the original job ordering so the shared `system_generator_cache.json` behavior is not made concurrent in this release.

## 15. Testing Strategy

The release must include automated tests demonstrating:

1. `compute_sasa=False` never invokes `calculate_buried_sasa()`.
2. `compute_sasa=True` retains existing SASA behavior.
3. Direct `analyze_structure()` callers retain full analysis by default.
4. `--analysis-source minimized` uses fast analysis for originals.
5. Minimized candidate acceptance QC uses fast analysis.
6. Rejected candidates never receive full/SASA analysis.
7. Accepted minimized candidates receive exactly one final full/SASA analysis.
8. `--analysis-source original` performs full original analysis but no minimized-candidate SASA.
9. A non-minimized run continues to perform full analysis.
10. `--minimize --workers 16` passes argument validation.
11. `--fail-fast --workers >1` remains rejected.
12. Parallel CPU analysis returns deterministic input-aligned results.
13. OpenMM minimization remains sequential even when `--workers >1`.
14. `--quiet` suppresses normal progress output but not warnings/errors.
15. Progress callbacks are optional and cannot change minimization results.
16. Stage timing/progress events occur around the existing three stage names.
17. Existing v2.1.0 OpenMM protocol-contract tests continue to pass without scientific parameter changes.
18. Existing output schemas remain unchanged.

Where possible, orchestration tests should use mocks/stubs rather than requiring OpenMM/CUDA. Existing backend-specific tests should continue to validate the real protocol contract independently.

## 16. Scientific Equivalence Requirement

For a fixed dataset and unchanged configuration, v2.2.0 full scientific analysis should match v2.1.0 within the same numerical behavior already present in the underlying libraries.

The intended invariant is:

```text
v2.1.0 scientific metrics ~= v2.2.0 scientific metrics
```

The expected difference is reduced wall-clock CPU time and improved runtime visibility, not altered structure interpretation.

## 17. Performance Expectation

The main guaranteed reduction is elimination of unnecessary Shrake-Rupley work:

- with `--analysis-source minimized`, original structures do not receive SASA;
- minimized candidates used for QC do not receive SASA;
- rejected candidates never receive SASA;
- only accepted minimized structures receive final SASA.

Additional speedup comes from parallelizing CPU full/fast analysis using `--workers N`.

No fixed speedup factor is specified because runtime depends on protein size, ligand size, number of structures, SASA cost, storage performance, ligand parameterization, and OpenMM setup/minimization behavior.

## 18. Files Expected to Change

Implementation is expected to focus on:

- `compound_contact_pipeline.py`
  - internal analysis mode;
  - phase orchestration;
  - worker validation semantics;
  - final full-analysis phase;
  - progress/ETA/quiet behavior;
- `restrained_openmm.py`
  - optional progress/timing instrumentation only;
- `tests/test_compound_contact_pipeline.py`
  - SASA scheduling, workers, ordering, progress, and orchestration tests;
- `tests/test_restrained_openmm.py`
  - optional progress hook and protocol-preservation tests;
- `README.md` and `INSTALL_AND_TUTORIAL.md`
  - updated worker semantics, progress output, `--quiet`, and v2.2 usage examples.

No unrelated refactoring is planned.

## 19. Release Criteria

v2.2.0 is ready when:

- all existing tests pass;
- all new tests above pass;
- `--minimize openmm --workers 16` is supported;
- OpenMM remains sequential on a single selected device;
- unnecessary SASA is demonstrably skipped according to the execution rules;
- final scientific output schemas are unchanged;
- default progress identifies CPU preparation versus CUDA minimization stages;
- `--quiet` suppresses informational progress;
- documentation describes the new worker semantics accurately;
- no v2.1.0 scientific minimization parameters have changed.
