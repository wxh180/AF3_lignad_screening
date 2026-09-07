# Compound Contact Pipeline v2.2.0 Performance and Progress Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build v2.2.0 so expensive SASA work is performed only for the selected scientific analysis source, CPU analysis uses `--workers N`, OpenMM remains sequential on one selected device, and long runs emit safe stage-aware progress/timing output.

**Architecture:** Keep the existing phase-oriented single-process orchestration in `compound_contact_pipeline.py`, but split structure analysis into fast geometry and full scientific modes with `compute_sasa`. CPU analysis phases use `ProcessPoolExecutor` with completion-order reporting and input-order result reconstruction; `run_minimization_jobs()` stays serial and passes a no-op-by-default progress callback into `restrained_openmm.py`. For `--analysis-source minimized`, accepted candidate CIFs are analyzed a second time in one final full/SASA CPU phase and only those full results enter scientific aggregation.

**Tech Stack:** Python 3, `concurrent.futures`, Biopython `MMCIFParser`/`ShrakeRupley`, NumPy, matplotlib, OpenMM, OpenFF Toolkit, openmmforcefields, pytest.

**Spec:** `docs/superpowers/specs/2026-09-07-compound-contact-v2.2-performance-design.md`

## Global Constraints

- Baseline is compound-contact-pipeline v2.1.0; release version is `2.2.0`.
- Do not change contact definitions, fixed clash definitions, buried SASA/interface-area definitions, plotting definitions, or scientific TSV schemas.
- Keep `MinimizationConfig.protocol == "openmm_restrained_v1"`, protein force field `amber/protein.ff14SB.xml`, ligand force field `gaff-2.2.20`, pH handling, tolerance, restraint strengths, RMSD gates, VDW gates, and the three existing minimization stages unchanged.
- OpenMM minimization remains sequential and uses one selected platform/device; do not add concurrent single-GPU jobs or multi-GPU scheduling.
- `--workers N` controls CPU analysis concurrency only.
- `--fail-fast --workers >1` remains invalid.
- Progress is informational, writes to stderr, and must never change scientific control flow or minimization success/failure.
- `--quiet` suppresses normal progress/timing output only; warnings and errors remain visible.
- Parallel analysis may report completion out of order, but returned `StructureResult` values must remain aligned with sorted input jobs.
- Do not add runtime columns to scientific output tables.

---

### Task 1: Add fast versus full structure analysis without changing default scientific behavior

**Files:**
- Modify: `compound_contact_pipeline.py:613-760`
- Test: `tests/test_compound_contact_pipeline.py` near the existing structure-analysis tests

**Interfaces:**
- Consumes: existing `StructureJob`, `StructureResult`, `StructureMetrics`, `calculate_buried_sasa()`.
- Produces: `analyze_structure(job, cutoff=4.5, clash_cutoff=2.0, compute_sasa=True) -> StructureResult`.
- Contract: `compute_sasa=False` leaves `buried_sasa_total_A2` and `interface_area_A2` as `None` and does not emit SASA warnings; omitted/default `compute_sasa` preserves v2.1.0 behavior.

- [ ] **Step 1: Write tests proving SASA can be skipped and remains enabled by default**

Add these tests using the existing `write_test_cif()` helper:

```python
def analysis_job(tmp_path):
    cif_path = write_test_cif(
        tmp_path / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (3.0, 0.0, 0.0))]),
        ],
    )
    return ccp.StructureJob(
        cif_path=cif_path,
        structure_id="C19/seed-1/model.cif",
        dataset_rel="C19",
        compound="C19",
        ligand_resname="WXI",
        protein_chains=("A",),
        ligand_chain="B",
    )


def test_analyze_structure_fast_mode_never_calls_sasa(monkeypatch, tmp_path):
    job = analysis_job(tmp_path)

    def forbidden_sasa(protein_residues, ligand_residues):
        raise AssertionError("SASA must not run in fast mode")

    monkeypatch.setattr(ccp, "calculate_buried_sasa", forbidden_sasa)
    result = ccp.analyze_structure(job, compute_sasa=False)

    assert result.status == "valid"
    assert result.metrics.buried_sasa_total_A2 is None
    assert result.metrics.interface_area_A2 is None
    assert all(record.stage != "sasa" for record in result.qc)


def test_analyze_structure_default_retains_full_sasa(monkeypatch, tmp_path):
    job = analysis_job(tmp_path)
    calls = []

    def fake_sasa(protein_residues, ligand_residues):
        calls.append((len(protein_residues), len(ligand_residues)))
        return 42.0, 21.0

    monkeypatch.setattr(ccp, "calculate_buried_sasa", fake_sasa)
    result = ccp.analyze_structure(job)

    assert calls == [(1, 1)]
    assert result.metrics.buried_sasa_total_A2 == pytest.approx(42.0)
    assert result.metrics.interface_area_A2 == pytest.approx(21.0)
```

- [ ] **Step 2: Run the new tests and verify they fail for the intended reason**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'fast_mode_never_calls_sasa or default_retains_full_sasa'
```

Expected: the fast-mode test fails because `analyze_structure()` does not yet accept `compute_sasa`.

- [ ] **Step 3: Add `compute_sasa=True` and guard only the existing SASA block**

Change the signature to:

```python
def analyze_structure(
    job: StructureJob,
    cutoff: float = 4.5,
    clash_cutoff: float = 2.0,
    compute_sasa: bool = True,
) -> StructureResult:
```

Replace the unconditional SASA section with:

```python
    buried_sasa = None
    interface_area = None
    if compute_sasa:
        try:
            buried_sasa, interface_area = calculate_buried_sasa(
                protein_residues, ligand_residues
            )
        except Exception as exc:
            qc.append(
                QCRecord(
                    severity="warning",
                    stage="sasa",
                    message=f"failed to calculate buried SASA: {exc}",
                    compound=job.compound,
                    structure_id=job.structure_id,
                )
            )
```

Do not alter parsing, chain selection, contact generation, fixed clash counting, or `StructureMetrics` field names.

- [ ] **Step 4: Run structure-analysis tests**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'analyze_structure or sasa or contact'
```

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "perf: allow structure analysis to skip sasa"
```

---

### Task 2: Add safe progress primitives, `--quiet`, and new worker validation semantics

**Files:**
- Modify: `compound_contact_pipeline.py:1-20, 2377-2460`
- Test: `tests/test_compound_contact_pipeline.py:947-1010` and new progress tests nearby

**Interfaces:**
- Produces: `ProgressReporter(quiet=False, stream=None)` with `info(message: str) -> None`.
- Produces: `_format_duration(seconds: float) -> str` for human-readable elapsed/ETA values.
- CLI adds `--quiet` as a boolean flag.
- Validation permits `--minimize openmm --workers N` for any `N >= 1`; `--fail-fast --workers >1` remains rejected.

- [ ] **Step 1: Replace the obsolete minimization-worker validation test and add quiet/progress safety tests**

Replace `test_cli_rejects_minimization_with_multiple_workers` with:

```python
def test_cli_allows_minimization_with_multiple_workers(tmp_path):
    args = ccp.build_parser().parse_args(
        [str(tmp_path), "--minimize", "openmm", "--workers", "16"]
    )

    ccp.validate_args(args)
    assert args.workers == 16
```

Add:

```python
def test_cli_quiet_defaults_false_and_can_be_enabled(tmp_path):
    parser = ccp.build_parser()
    assert parser.parse_args([str(tmp_path)]).quiet is False
    assert parser.parse_args([str(tmp_path), "--quiet"]).quiet is True


def test_progress_reporter_quiet_suppresses_info(capsys):
    reporter = ccp.ProgressReporter(quiet=True)
    reporter.info("hidden progress")
    assert capsys.readouterr().err == ""


def test_progress_reporter_failure_is_side_effect_free():
    class BrokenStream:
        def write(self, value):
            raise OSError("stream unavailable")

        def flush(self):
            raise OSError("stream unavailable")

    reporter = ccp.ProgressReporter(stream=BrokenStream())
    reporter.info("best effort only")


def test_quiet_does_not_suppress_cli_errors(tmp_path, capsys):
    code = ccp.main([str(tmp_path), "--quiet", "--cutoff", "0"])
    captured = capsys.readouterr()
    assert code == 2
    assert "ERROR:" in captured.err
```

- [ ] **Step 2: Run the focused CLI/progress tests and verify failure**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'multiple_workers or quiet or progress_reporter'
```

Expected: failures for missing `--quiet`, missing `ProgressReporter`, and the old validation rule.

- [ ] **Step 3: Implement the safe reporter and duration formatter**

Add `import time` and define near the top-level utility helpers:

```python
class ProgressReporter:
    def __init__(self, quiet: bool = False, stream=None):
        self.quiet = quiet
        self.stream = sys.stderr if stream is None else stream

    def info(self, message: str) -> None:
        if self.quiet:
            return
        try:
            print(message, file=self.stream, flush=True)
        except Exception:
            return


def _format_duration(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"
```

The broad exception handler is intentional only for best-effort progress output; do not use this pattern for scientific work.

- [ ] **Step 4: Add `--quiet` and remove only the minimization/worker restriction**

Add to `build_parser()`:

```python
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="suppress normal progress and timing messages",
    )
```

Delete only:

```python
    if args.minimize and args.workers != 1:
        raise PipelineError("--minimize requires --workers 1")
```

Keep:

```python
    if args.fail_fast and args.workers != 1:
        raise PipelineError("--fail-fast requires --workers 1")
```

- [ ] **Step 5: Run focused tests**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'multiple_workers or quiet or progress_reporter or fail_fast'
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: add safe progress controls and cpu worker semantics"
```

---

### Task 3: Make CPU analysis phase-aware, timed, parallel, and deterministically ordered

**Files:**
- Modify: `compound_contact_pipeline.py:1-20, 2499-2527`
- Test: `tests/test_compound_contact_pipeline.py` near `test_cli_parallel_analysis_preserves_sorted_structure_order`

**Interfaces:**
- Consumes: `ProgressReporter`, `analyze_structure(..., compute_sasa=...)`.
- Produces: `_timed_analyze_job(index, job, cutoff, clash_cutoff, compute_sasa) -> tuple[int, StructureResult, float]`.
- Produces: `_analyze_jobs(jobs, args, *, compute_sasa=True, reporter=None, phase_label="Structure analysis", phase_number=None, phase_total=None) -> list[StructureResult]`.
- Contract: parent process emits progress; worker processes never write progress directly; returned results are aligned with input jobs.

- [ ] **Step 1: Add direct tests for fast parallel analysis, deterministic ordering, and parent-side progress**

Add the following tests:

```python
def test_parallel_fast_analysis_preserves_input_order_and_skips_sasa(tmp_path, capsys):
    root = tmp_path / "inputs"
    for seed, x in [("seed-2", 3.2), ("seed-1", 3.0)]:
        write_test_cif(
            root / "C19" / seed / "model.cif",
            [
                ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
                ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (x, 0.0, 0.0))]),
            ],
        )
    jobs = ccp.build_structure_jobs(root, None, "WXI", ("A",), "B")
    args = ccp.build_parser().parse_args([str(root), "--workers", "2"])
    reporter = ccp.ProgressReporter()

    results = ccp._analyze_jobs(
        jobs,
        args,
        compute_sasa=False,
        reporter=reporter,
        phase_label="Pre-minimization analysis",
        phase_number=1,
        phase_total=3,
    )

    assert [result.job.structure_id for result in results] == [
        job.structure_id for job in jobs
    ]
    assert all(result.metrics.buried_sasa_total_A2 is None for result in results)
    stderr = capsys.readouterr().err
    assert "Phase 1/3: Pre-minimization analysis" in stderr
    assert "2 structures | 2 CPU workers | SASA disabled" in stderr
    assert "Phase completed:" in stderr
```

Keep the existing CLI-level deterministic-order test as a regression test.

- [ ] **Step 2: Run the new direct test and verify failure**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'parallel_fast_analysis_preserves_input_order'
```

Expected: `_analyze_jobs()` does not yet accept the phase/compute/reporting keyword arguments.

- [ ] **Step 3: Add the timed top-level worker helper**

Import `as_completed`:

```python
from concurrent.futures import ProcessPoolExecutor, as_completed
```

Add:

```python
def _timed_analyze_job(
    index: int,
    job: StructureJob,
    cutoff: float,
    clash_cutoff: float,
    compute_sasa: bool,
) -> tuple[int, StructureResult, float]:
    started = time.perf_counter()
    result = analyze_structure(
        job,
        cutoff,
        clash_cutoff,
        compute_sasa=compute_sasa,
    )
    return index, result, time.perf_counter() - started
```

Keep it top-level so it is pickleable by `ProcessPoolExecutor`.

- [ ] **Step 4: Refactor `_analyze_jobs()` to report completion order but return input order**

Use this signature:

```python
def _analyze_jobs(
    jobs: list[StructureJob],
    args: argparse.Namespace,
    *,
    compute_sasa: bool = True,
    reporter: ProgressReporter | None = None,
    phase_label: str = "Structure analysis",
    phase_number: int | None = None,
    phase_total: int | None = None,
) -> list[StructureResult]:
```

At entry, set `reporter = reporter or ProgressReporter(quiet=True)`, capture a phase start time, and emit:

```python
    prefix = (
        f"Phase {phase_number}/{phase_total}: "
        if phase_number is not None and phase_total is not None
        else ""
    )
    reporter.info(f"=== {prefix}{phase_label} ===")
    reporter.info(
        f"{len(jobs)} structures | {args.workers} CPU workers | "
        f"SASA {'enabled' if compute_sasa else 'disabled'}"
    )
```

For one worker, call `_timed_analyze_job()` in a loop, preserve the current fail-fast check, and emit one line after each completion. For multiple workers, use `executor.submit()` plus `as_completed()`, assign each result into a pre-sized list by returned index, and emit progress in future-completion order. The result reconstruction must be equivalent to:

```python
    ordered: list[StructureResult | None] = [None] * len(jobs)
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [
            executor.submit(
                _timed_analyze_job,
                index,
                job,
                args.cutoff,
                args.clash_cutoff,
                compute_sasa,
            )
            for index, job in enumerate(jobs)
        ]
        completed = 0
        for future in as_completed(futures):
            index, result, elapsed = future.result()
            ordered[index] = result
            completed += 1
            reporter.info(
                f"[{completed:>{len(str(len(jobs)))}}/{len(jobs)}] "
                f"{result.job.structure_id} completed    {elapsed:.2f} s"
            )
    results = [result for result in ordered if result is not None]
    if len(results) != len(jobs):
        raise PipelineError("parallel analysis did not return one result per job")
```

Finish both serial and parallel paths with:

```python
    reporter.info(
        f"Phase completed: {time.perf_counter() - phase_started:.1f} s"
    )
```

- [ ] **Step 5: Run direct and existing parallel-order tests**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'parallel_fast_analysis or parallel_analysis_preserves_sorted_structure_order or fail_fast'
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "perf: parallelize timed cpu analysis phases"
```

---

### Task 4: Add optional, failure-safe OpenMM preparation and stage progress events

**Files:**
- Modify: `restrained_openmm.py:1-20, 58-105, 592-711`
- Test: `tests/test_restrained_openmm.py` near the backend execution tests and protocol-contract test

**Interfaces:**
- Produces: `MinimizationProgressEvent(operation: str, state: str, elapsed_seconds: float | None = None)`.
- Produces: optional `progress_callback: Callable[[MinimizationProgressEvent], None] | None = None` parameter on `minimize_structure()`.
- Internal `_emit_progress()` catches callback exceptions.
- Event operations are exactly `preparation`, each `MinimizationStage.name`, and `finalization`; states are `started` and `completed`.

- [ ] **Step 1: Add backend tests for event sequence and callback failure isolation**

Add `Callable` support in implementation later; add tests using the existing `backend_stack` fixture:

```python
def test_backend_emits_progress_around_existing_three_stages(tmp_path, backend_stack):
    request = backend_request(tmp_path)
    events = []

    result = rom._minimize_structure(
        request,
        rom.MinimizationConfig(),
        tmp_path / "cache.json",
        backend_stack,
        progress_callback=events.append,
    )

    assert result.output_path == request.output_path
    assert [(event.operation, event.state) for event in events] == [
        ("preparation", "started"),
        ("preparation", "completed"),
        ("hydrogen_relaxation", "started"),
        ("hydrogen_relaxation", "completed"),
        ("pocket_relaxation", "started"),
        ("pocket_relaxation", "completed"),
        ("gentle_relaxation", "started"),
        ("gentle_relaxation", "completed"),
        ("finalization", "started"),
        ("finalization", "completed"),
    ]
    completed = [event for event in events if event.state == "completed"]
    assert all(event.elapsed_seconds is not None for event in completed)
    assert all(event.elapsed_seconds >= 0 for event in completed)


def test_backend_progress_callback_failure_cannot_fail_minimization(tmp_path, backend_stack):
    request = backend_request(tmp_path)

    def broken_callback(event):
        raise RuntimeError("display failure")

    result = rom._minimize_structure(
        request,
        rom.MinimizationConfig(),
        tmp_path / "cache.json",
        backend_stack,
        progress_callback=broken_callback,
    )

    assert result.output_path.exists()
```

Use the existing `backend_request(tmp_path, **changes)` helper defined in `tests/test_restrained_openmm.py:147`; do not create a duplicate request fixture.

- [ ] **Step 2: Run the progress tests and verify failure**

Run:

```bash
pytest -q tests/test_restrained_openmm.py -k 'emits_progress or callback_failure'
```

Expected: `_minimize_structure()` does not yet accept `progress_callback`.

- [ ] **Step 3: Add progress event type and safe emitter**

Update imports:

```python
from typing import Callable, Sequence
import time
```

Add after `BackendMinimizationResult`:

```python
@dataclass(frozen=True)
class MinimizationProgressEvent:
    operation: str
    state: str
    elapsed_seconds: float | None = None


ProgressCallback = Callable[[MinimizationProgressEvent], None]


def _emit_progress(
    callback: ProgressCallback | None,
    event: MinimizationProgressEvent,
) -> None:
    if callback is None:
        return
    try:
        callback(event)
    except Exception:
        return
```

- [ ] **Step 4: Thread the optional callback through `minimize_structure()`**

Change the public signature to:

```python
def minimize_structure(
    request: MinimizationRequest,
    config: MinimizationConfig,
    cache_path: Path,
    progress_callback: ProgressCallback | None = None,
) -> BackendMinimizationResult:
```

Pass it through:

```python
        return _minimize_structure(
            request,
            config,
            cache_path,
            stack,
            progress_callback=progress_callback,
        )
```

Change the internal signature to:

```python
def _minimize_structure(
    request,
    config,
    cache_path,
    stack,
    progress_callback: ProgressCallback | None = None,
):
```

- [ ] **Step 5: Instrument preparation, each existing minimization stage, and finalization**

At the beginning of `_minimize_structure()`:

```python
    preparation_started = time.perf_counter()
    _emit_progress(
        progress_callback,
        MinimizationProgressEvent("preparation", "started"),
    )
```

Leave the existing chemistry/topology/SystemGenerator/hydrogen/system/context code scientifically unchanged. Immediately after `initial_energy, _ = _physical_state(context, unit)` emit:

```python
        _emit_progress(
            progress_callback,
            MinimizationProgressEvent(
                "preparation",
                "completed",
                time.perf_counter() - preparation_started,
            ),
        )
```

Inside the unchanged `for stage in config.stages:` loop, emit `started` immediately before `LocalEnergyMinimizer.minimize()` and `completed` immediately after `_physical_state()`:

```python
            stage_started = time.perf_counter()
            _emit_progress(
                progress_callback,
                MinimizationProgressEvent(stage.name, "started"),
            )
            mm.LocalEnergyMinimizer.minimize(
                context,
                config.tolerance_kj_mol_nm
                * unit.kilojoule_per_mole
                / unit.nanometer,
                stage.max_iterations,
            )
            final_energy, final_nm = _physical_state(context, unit)
            _emit_progress(
                progress_callback,
                MinimizationProgressEvent(
                    stage.name,
                    "completed",
                    time.perf_counter() - stage_started,
                ),
            )
```

After releasing the context and before output-coordinate/QC work, emit `finalization started`; after `_write_original_cif()` and package-version collection, emit `finalization completed`. Do not move or change any scientific calculation.

- [ ] **Step 6: Run backend progress and protocol-contract tests**

Run:

```bash
pytest -q tests/test_restrained_openmm.py -k 'progress or protocol_has_three_bounded_stages or backend_minimizes_with_physical_energy'
```

Expected: PASS, including the unchanged protocol-contract assertions.

- [ ] **Step 7: Commit**

```bash
git add restrained_openmm.py tests/test_restrained_openmm.py
git commit -m "feat: expose safe openmm stage progress events"
```

---

### Task 5: Keep minimization sequential, force fast candidate QC, and report OpenMM/ETA progress

**Files:**
- Modify: `compound_contact_pipeline.py:2544-2683`
- Test: `tests/test_compound_contact_pipeline.py:1040-1274` plus new sequencing/progress tests

**Interfaces:**
- Consumes: `MinimizationProgressEvent`, `ProgressReporter`.
- `run_minimization_jobs(..., reporter: ProgressReporter | None = None)` remains a serial loop over input jobs.
- Candidate calls are exactly `analyze_structure(candidate_job, args.cutoff, args.clash_cutoff, compute_sasa=False)`.
- Injected test minimizers adopt the optional keyword `progress_callback=None`.

- [ ] **Step 1: Update the fake minimizer signature and add tests proving candidate SASA is skipped and job execution stays serial with multiple CPU workers**

Change the existing fake helper signature to:

```python
    def minimize(request, config, cache_path, progress_callback=None):
```

When `progress_callback` is provided, the fake may emit no events; orchestration tests should not require backend instrumentation.

Add a candidate-SASA test:

```python
def test_minimization_candidate_qc_skips_sasa(monkeypatch, tmp_path):
    _, args, jobs, original_results = minimization_fixture(tmp_path)
    calls = []
    real_analyze = ccp.analyze_structure

    def tracking_analyze(job, cutoff=4.5, clash_cutoff=2.0, compute_sasa=True):
        calls.append((job.cif_path, compute_sasa))
        return real_analyze(
            job,
            cutoff,
            clash_cutoff,
            compute_sasa=compute_sasa,
        )

    monkeypatch.setattr(ccp, "analyze_structure", tracking_analyze)
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_minimization_backend(),
    )

    assert records[0].status == "accepted"
    assert calls == [(records[0].staged_cif_path, False)]
    assert records[0].candidate_result.metrics.buried_sasa_total_A2 is None
```

Add a serial-order test with at least two jobs and `args.workers = 4`; the fake backend appends `request.structure_id` at entry and asserts the previous call has returned before the next begins. The final assertion must equal the original `jobs` order.

- [ ] **Step 2: Run the new tests and verify the SASA test fails**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'candidate_qc_skips_sasa or minimization_remains_sequential'
```

Expected: candidate analysis currently uses full SASA.

- [ ] **Step 3: Add an OpenMM-event formatter in the pipeline layer**

Add a helper that receives the structure reporter and configured platform:

```python
def _openmm_progress_callback(
    reporter: ProgressReporter,
    platform: str,
):
    def callback(event) -> None:
        suffix = f" [{platform}]" if event.operation in {
            "hydrogen_relaxation",
            "pocket_relaxation",
            "gentle_relaxation",
        } else ""
        if event.state == "started":
            label = (
                "OpenMM preparation"
                if event.operation == "preparation"
                else event.operation
            )
            reporter.info(f"    {label}{suffix} ...")
            return
        label = (
            "preparation"
            if event.operation == "preparation"
            else event.operation
        )
        elapsed = 0.0 if event.elapsed_seconds is None else event.elapsed_seconds
        reporter.info(f"    {label} completed ........ {elapsed:.2f} s")

    return callback
```

Because `ProgressReporter.info()` is already failure-safe, this callback is also operationally safe.

- [ ] **Step 4: Refactor `run_minimization_jobs()` progress without introducing concurrency**

Add the optional reporter parameter:

```python
def run_minimization_jobs(
    jobs: list[StructureJob],
    original_results: list[StructureResult],
    args: argparse.Namespace,
    temporary_root: Path,
    minimize_one=None,
    reporter: ProgressReporter | None = None,
) -> list[MinimizationRecord]:
```

Set `reporter = reporter or ProgressReporter(quiet=True)`, record `batch_started = time.perf_counter()`, and emit the phase header from the caller rather than creating worker threads here.

At the beginning of each valid job, emit:

```python
        item_started = time.perf_counter()
        reporter.info(f"[{len(records) + 1}/{len(jobs)}] {job.structure_id}")
```

Call the backend with:

```python
            backend = minimize_one(
                request,
                config,
                cache_path,
                progress_callback=_openmm_progress_callback(
                    reporter, config.platform
                ),
            )
```

Then time candidate QC and call fast analysis only:

```python
            qc_started = time.perf_counter()
            candidate = analyze_structure(
                candidate_job,
                args.cutoff,
                args.clash_cutoff,
                compute_sasa=False,
            )
            reporter.info(
                "    candidate geometry QC completed "
                f"........ {time.perf_counter() - qc_started:.2f} s"
            )
```

After appending each completed/failed/not-attempted record, emit status, item total, elapsed batch time, and ETA. Use completed record count for ETA:

```python
        completed_count = len(records)
        batch_elapsed = time.perf_counter() - batch_started
        eta_seconds = (
            batch_elapsed / completed_count * (len(jobs) - completed_count)
            if completed_count
            else 0.0
        )
        reporter.info(
            f"Batch: {completed_count}/{len(jobs)} | "
            f"elapsed {_format_duration(batch_elapsed)} | "
            f"ETA {_format_duration(eta_seconds)}"
        )
```

Do not use `args.workers` to create OpenMM concurrency.

- [ ] **Step 5: Keep failure semantics and contact-set Jaccard unchanged**

Verify the acceptance call still uses:

```python
            decision = evaluate_acceptance(
                backbone_rmsd_A=backend.backbone_rmsd_A,
                ligand_rmsd_A=backend.ligand_rmsd_A,
                before_fixed_clashes=original.metrics.clash_pair_count,
                after_fixed_clashes=candidate.metrics.clash_pair_count,
                before_vdw_clashes=backend.before_vdw_clash_pair_count,
                after_vdw_clashes=backend.after_vdw_clash_pair_count,
                max_backbone_rmsd_A=config.max_backbone_rmsd_A,
                max_ligand_rmsd_A=config.max_ligand_rmsd_A,
            )
```

and Jaccard still compares `_contact_residue_set(original)` with `_contact_residue_set(candidate)`.

- [ ] **Step 6: Run minimization orchestration regressions**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'minimization or accepted_minimization or rejected_minimization or failed_minimization'
```

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "perf: use fast serial minimization qc with progress"
```

---

### Task 6: Add the final accepted-minimized full/SASA phase and wire exact analysis-source scheduling

**Files:**
- Modify: `compound_contact_pipeline.py:2650-2715, 2979-3040`
- Test: `tests/test_compound_contact_pipeline.py:1177-1274, 1395-1450` plus new SASA scheduling tests

**Interfaces:**
- Produces: `_analyze_accepted_minimized(jobs, minimization_records, args, reporter, phase_number, phase_total) -> list[StructureResult | None]` aligned one-to-one with input jobs.
- `select_analysis_results(..., full_minimized_results: list[StructureResult | None] | None = None)` uses full minimized results only for accepted entries.
- `_run_pipeline()` chooses `compute_sasa` automatically from minimization and `analysis_source`.

- [ ] **Step 1: Add orchestration tests for the three SASA scheduling modes**

Add a tracking helper around `analyze_structure()` that records `(structure_id, compute_sasa, cif_path)` and delegates to the real function. Add these tests:

```python
def test_minimized_analysis_runs_fast_original_fast_candidate_then_full_accepted(
    monkeypatch, tmp_path
):
    root = tmp_path / "inputs"
    write_test_cif(
        root / "C19" / "seed-1" / "model.cif",
        [
            ("A", " ", 10, " ", "ALA", [("CA", "C", (0.0, 0.0, 0.0))]),
            ("B", "H_WXI", 1, " ", "WXI", [("C1", "C", (1.5, 0.0, 0.0))]),
        ],
    )
    output = tmp_path / "results"
    args = ccp.build_parser().parse_args(
        [
            str(root), "--ligand", "WXI", "--protein-chain", "A",
            "--ligand-chain", "B", "--ligand-smiles", "CC",
            "--minimize", "openmm", "--analysis-source", "minimized",
            "--out-dir", str(output),
        ]
    )
    real_analyze = ccp.analyze_structure
    modes = []

    def tracking_analyze(job, cutoff=4.5, clash_cutoff=2.0, compute_sasa=True):
        modes.append(compute_sasa)
        return real_analyze(
            job, cutoff, clash_cutoff, compute_sasa=compute_sasa
        )

    monkeypatch.setattr(ccp, "analyze_structure", tracking_analyze)
    assert ccp.run_pipeline(args, minimize_one=fake_minimization_backend()) == 0
    assert modes == [False, False, True]
```

Add a rejected-candidate variant using `ligand_rmsd_A=1.6` and assert modes are `[False, False]` with no final full call.

Add an `--analysis-source original` variant and assert modes are `[True, False]`.

Add a no-minimization variant and assert the only analysis mode is `[True]`.

- [ ] **Step 2: Run the new scheduling tests and verify failure**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'runs_fast_original or rejected_candidate_variant or analysis_source_original_sasa or no_minimization_full_sasa'
```

Expected: current pipeline runs full analysis before minimization and reuses the candidate fast/full result instead of a dedicated final phase.

- [ ] **Step 3: Implement aligned full analysis for accepted minimized structures**

Add:

```python
def _analyze_accepted_minimized(
    jobs: list[StructureJob],
    minimization_records: list[MinimizationRecord],
    args: argparse.Namespace,
    reporter: ProgressReporter,
    phase_number: int,
    phase_total: int,
) -> list[StructureResult | None]:
    accepted_indices: list[int] = []
    accepted_jobs: list[StructureJob] = []
    for index, record in enumerate(minimization_records):
        if record.status != "accepted" or record.candidate_result is None:
            continue
        accepted_indices.append(index)
        accepted_jobs.append(record.candidate_result.job)

    aligned: list[StructureResult | None] = [None] * len(jobs)
    if not accepted_jobs:
        return aligned

    full_results = _analyze_jobs(
        accepted_jobs,
        args,
        compute_sasa=True,
        reporter=reporter,
        phase_label="Final minimized scientific analysis",
        phase_number=phase_number,
        phase_total=phase_total,
    )
    for index, result in zip(accepted_indices, full_results):
        aligned[index] = result
    return aligned
```

- [ ] **Step 4: Make `select_analysis_results()` require full accepted results for minimized scientific analysis**

Change the signature to include:

```python
    full_minimized_results: list[StructureResult | None] | None = None,
```

For `analysis_source == "original"`, return `original_results` exactly as before. For minimized analysis, require `full_minimized_results` to be aligned in length, and for each accepted record append the aligned full result. If a full accepted analysis returns an invalid result, preserve that invalid result rather than falling back to `candidate_result`. Rejected/failed/not-attempted entries continue to use `_invalid_structure_result()` with the existing minimization stages/reasons.

- [ ] **Step 5: Refactor `_run_pipeline()` into exact phase scheduling**

Create one reporter after validation:

```python
    reporter = ProgressReporter(quiet=args.quiet)
```

Remove the old unconditional:

```python
    print(f"Analyzing {len(jobs)} CIF structure(s)...", file=sys.stderr)
    original_results = _analyze_jobs(jobs, args)
```

Use these rules:

```python
    if not args.minimize:
        original_results = _analyze_jobs(
            jobs,
            args,
            compute_sasa=True,
            reporter=reporter,
            phase_label="Structure analysis",
            phase_number=1,
            phase_total=1,
        )
    else:
        total_phases = 3 if args.analysis_source == "minimized" else 2
        original_results = _analyze_jobs(
            jobs,
            args,
            compute_sasa=(args.analysis_source == "original"),
            reporter=reporter,
            phase_label="Pre-minimization analysis",
            phase_number=1,
            phase_total=total_phases,
        )
```

Before `run_minimization_jobs()`, emit:

```python
        reporter.info(
            f"=== Phase 2/{total_phases}: Restrained OpenMM minimization ==="
        )
        device = (
            f" device {args.openmm_device_index}"
            if args.openmm_device_index is not None
            else ""
        )
        reporter.info(
            f"{len(jobs)} structures | {args.openmm_platform}{device} | sequential"
        )
```

Call minimization with `reporter=reporter`.

For `analysis_source == "minimized"`, call `_analyze_accepted_minimized(..., phase_number=3, phase_total=3)` and pass its aligned list into `select_analysis_results()`. For `analysis_source == "original"`, do not run a final minimized SASA phase.

- [ ] **Step 6: Ensure final full-analysis failure never falls back to fast candidate data**

Add a test that makes the final `compute_sasa=True` analysis return an invalid `StructureResult` for an accepted candidate and assert the selected result is invalid with its final-analysis QC. The test must also assert the fast `candidate_result` remains present only in the `MinimizationRecord` for QC/provenance.

- [ ] **Step 7: Run all orchestration and schema regressions**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k 'minimized_analysis or original_analysis or run_pipeline or stable_schema or source_aware or parallel_analysis'
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "perf: defer full sasa to selected analysis source"
```

---

### Task 7: Verify quiet/progress behavior end-to-end and preserve all scientific output contracts

**Files:**
- Modify: `tests/test_compound_contact_pipeline.py`
- Modify only if a test exposes a defect: `compound_contact_pipeline.py`, `restrained_openmm.py`

**Interfaces:**
- No new public interfaces.
- Confirms operational output is stderr-only and generated paths remain stdout-only.

- [ ] **Step 1: Add end-to-end progress/quiet tests around a fake minimization backend**

Add one normal-progress test that runs a single accepted minimized structure and asserts stderr contains all orchestration labels that the fake backend can legitimately trigger:

```python
def test_minimized_pipeline_reports_phase_and_candidate_progress(tmp_path, capsys):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.out_dir = tmp_path / "results"

    code = ccp.run_pipeline(args, minimize_one=fake_minimization_backend())
    captured = capsys.readouterr()

    assert code == 0
    assert "Phase 1/3: Pre-minimization analysis" in captured.err
    assert "Phase 2/3: Restrained OpenMM minimization" in captured.err
    assert "candidate geometry QC completed" in captured.err
    assert "Phase 3/3: Final minimized scientific analysis" in captured.err
    assert "Batch: 1/1" in captured.err
    assert "generated_files.tsv" in captured.out
```

Add a quiet version:

```python
def test_quiet_minimized_pipeline_suppresses_normal_progress(tmp_path, capsys):
    _, args, _, _ = minimization_fixture(tmp_path)
    args.analysis_source = "minimized"
    args.quiet = True
    args.out_dir = tmp_path / "results"

    code = ccp.run_pipeline(args, minimize_one=fake_minimization_backend())
    captured = capsys.readouterr()

    assert code == 0
    assert "Phase 1/3" not in captured.err
    assert "candidate geometry QC" not in captured.err
    assert "Batch: 1/1" not in captured.err
    assert "generated_files.tsv" in captured.out
```

- [ ] **Step 2: Run the complete pipeline test module**

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py
```

Expected: PASS.

- [ ] **Step 3: Run the complete OpenMM backend test module**

Run:

```bash
pytest -q tests/test_restrained_openmm.py
```

Expected: PASS; the existing `test_protocol_has_three_bounded_stages` must remain unchanged and passing.

- [ ] **Step 4: Commit any test-only additions or narrowly required corrections**

```bash
git add compound_contact_pipeline.py restrained_openmm.py tests/test_compound_contact_pipeline.py tests/test_restrained_openmm.py
git commit -m "test: cover v2.2 progress and scientific compatibility"
```

---

### Task 8: Update version and user documentation for v2.2.0

**Files:**
- Modify: `compound_contact_pipeline.py:20`
- Modify: `README.md`
- Modify: `INSTALL_AND_TUTORIAL.md`

**Interfaces:**
- CLI reports version `2.2.0`.
- Documentation explains that `--workers` affects CPU analysis only and that OpenMM stays serial per selected device.
- Documentation shows `--workers 16` for the user's 16-core CUDA workflow and documents `--quiet`.

- [ ] **Step 1: Add/adjust the version assertion**

Add to `tests/test_compound_contact_pipeline.py`:

```python
def test_version_is_2_2_0():
    assert ccp.__version__ == "2.2.0"
```

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k version_is_2_2_0
```

Expected: FAIL because the baseline still reports `2.1.0`.

- [ ] **Step 2: Update the version constant**

Change:

```python
__version__ = "2.2.0"
```

Run:

```bash
pytest -q tests/test_compound_contact_pipeline.py -k version_is_2_2_0
```

Expected: PASS.

- [ ] **Step 3: Update README usage and performance notes**

Add a v2.2 section that states all of the following explicitly:

```text
- SASA scheduling is automatic; there is no --sasa-mode option.
- With --analysis-source minimized, original and candidate-QC analyses skip SASA; only accepted minimized structures receive full SASA.
- --workers N controls CPU analysis concurrency only.
- OpenMM minimization remains sequential on --openmm-platform/--openmm-device-index.
- Progress and timings go to stderr by default.
- --quiet suppresses normal progress without suppressing warnings/errors.
```

Include this exact practical example:

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

- [ ] **Step 4: Update INSTALL_AND_TUTORIAL troubleshooting**

Document how to interpret progress labels:

```text
OpenMM preparation...       CPU-side topology/chemistry/system setup
hydrogen_relaxation [CUDA]  LocalEnergyMinimizer stage on selected CUDA device
pocket_relaxation [CUDA]    LocalEnergyMinimizer stage on selected CUDA device
gentle_relaxation [CUDA]    LocalEnergyMinimizer stage on selected CUDA device
candidate geometry QC       CPU-side post-minimization geometry/contact analysis
Final minimized scientific analysis  CPU-side full analysis including SASA
```

Also document that high CPU with no GPU use during `OpenMM preparation` is expected, whereas GPU activity should be visible during the named CUDA stages.

- [ ] **Step 5: Run documentation-sensitive CLI smoke checks**

Run:

```bash
python compound_contact_pipeline.py --version
python compound_contact_pipeline.py --help | grep -E -- '--workers|--quiet|--analysis-source|--openmm-platform'
```

Expected: version prints `2.2.0`, and all four options are present in help.

- [ ] **Step 6: Commit**

```bash
git add compound_contact_pipeline.py README.md INSTALL_AND_TUTORIAL.md tests/test_compound_contact_pipeline.py
git commit -m "docs: release compound contact pipeline v2.2.0"
```

---

### Task 9: Full regression verification and release archive

**Files:**
- Verify: `compound_contact_pipeline.py`
- Verify: `restrained_openmm.py`
- Verify: `tests/test_compound_contact_pipeline.py`
- Verify: `tests/test_restrained_openmm.py`
- Verify: `README.md`
- Verify: `INSTALL_AND_TUTORIAL.md`
- Create outside git tree: `compound-contact-pipeline-v2.2.0.zip`

**Interfaces:**
- Final release artifact contains the source, tests, environment file, examples, README, and installation/tutorial documentation.

- [ ] **Step 1: Run the entire automated suite from the repository root**

Run:

```bash
pytest -q
```

Expected: all tests pass with no changes to the v2.1.0 protocol-contract expectations other than the intended new v2.2 tests.

- [ ] **Step 2: Verify the source contains no obsolete minimization-worker restriction**

Run:

```bash
grep -n -- '--minimize requires --workers 1' compound_contact_pipeline.py && exit 1 || true
grep -n '__version__ = "2.2.0"' compound_contact_pipeline.py
grep -n 'MinimizationStage("hydrogen_relaxation", 500, 10.0, 10.0, 10.0, 10.0)' restrained_openmm.py
grep -n 'MinimizationStage("pocket_relaxation", 1000, 10.0, 10.0, 0.0, 1.0)' restrained_openmm.py
grep -n 'MinimizationStage("gentle_relaxation", 1000, 2.0, 2.0, 0.0, 0.2)' restrained_openmm.py
```

Expected: no obsolete worker restriction, version is 2.2.0, and all three protocol stages remain byte-for-byte equivalent in parameter values.

- [ ] **Step 3: Review git diff against the imported v2.1.0 baseline for scientific drift**

Run:

```bash
git diff 9aa752d -- compound_contact_pipeline.py restrained_openmm.py
```

Review specifically that modifications are limited to SASA scheduling, CPU orchestration, progress/timing hooks, worker validation, and versioning. Reject any unintended force-field, restraint, cutoff, RMSD, acceptance, contact, or SASA-formula change.

- [ ] **Step 4: Confirm repository is clean and commits are present**

Run:

```bash
git status --short
git log --oneline --decorate -12
```

Expected: empty `git status --short` and a sequence of focused v2.2 commits after design/plan commits.

- [ ] **Step 5: Build a clean release archive from tracked files**

From the parent directory of the repository, run:

```bash
rm -f compound-contact-pipeline-v2.2.0.zip
git -C compound-contact-pipeline-v2.2.0-dev archive \
  --format=zip \
  --prefix=compound-contact-pipeline-v2.2.0/ \
  -o ../compound-contact-pipeline-v2.2.0.zip \
  HEAD
sha256sum compound-contact-pipeline-v2.2.0.zip
```

Expected: one reproducible release ZIP plus a SHA-256 checksum for handoff.

- [ ] **Step 6: Final release verification**

Extract the archive into a temporary directory and run:

```bash
rm -rf /tmp/compound-contact-pipeline-v2.2.0-check
mkdir -p /tmp/compound-contact-pipeline-v2.2.0-check
unzip -q compound-contact-pipeline-v2.2.0.zip -d /tmp/compound-contact-pipeline-v2.2.0-check
cd /tmp/compound-contact-pipeline-v2.2.0-check/compound-contact-pipeline-v2.2.0
python compound_contact_pipeline.py --version
pytest -q
```

Expected: version `2.2.0` and full test suite PASS from the packaged release.
