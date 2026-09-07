# Restrained OpenMM Minimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Release version 2.0.0 of the compound-contact pipeline with opt-in staged OpenMM minimization, conservative pose-acceptance gates, minimized CIF outputs, and auditable before/after QC.

**Architecture:** Keep the existing CLI and aggregation flow in `compound_contact_pipeline.py`, and isolate optional molecular-mechanics imports and coordinate work in a new `restrained_openmm.py` module. The main pipeline first analyzes original structures, calls a dependency-injected minimizer when requested, accepts or rejects each candidate from pure QC functions, and explicitly selects original or accepted minimized coordinates for aggregation.

**Tech Stack:** Python 3.10+, NumPy, Biopython, Matplotlib, pytest; optional OpenMM, openmmforcefields, OpenFF Toolkit, RDKit, and AmberTools.

**Spec:** `docs/superpowers/specs/2026-09-03-restrained-openmm-minimization-design.md`

## Global Constraints

- Existing behavior remains the default when `--minimize` is absent.
- Input CIFs are read-only and must remain byte-for-byte unchanged.
- OpenMM-related imports occur only when `--minimize openmm` is requested.
- Protein parameters default to `amber/protein.ff14SB.xml`; ligand parameters default to `gaff-2.2.20`.
- The protocol is nonperiodic `NoCutoff`, constrains bonds involving hydrogen, adds hydrogens at pH 7.4, and defaults to the CPU platform.
- Stage limits are 500, 1,000, and 1,000 L-BFGS iterations with the force constants defined in the specification.
- Default acceptance limits are protein-backbone RMSD <= 0.5 A and ligand-heavy-atom RMSD <= 1.5 A, with no increase in either clash count.
- Minimized complex potential energy is same-system QC only and must never enter binding-energy or affinity summaries.
- Minimized structures preserve original `structure_id` values; rejected and failed candidates are never included in minimized-coordinate aggregation.
- Minimized output mmCIF files contain only selected protein/ligand atoms that existed in the input; transient added hydrogens are not published.
- Version 2.0 supports one selected ligand residue per structure and noncovalent ligands only.
- Every production-code change follows a witnessed RED-GREEN-REFACTOR cycle.

## File Structure

- Modify `compound_contact_pipeline.py`: version, manifest/CLI schema, minimization orchestration, analysis-source selection, summary serialization, and plot registration.
- Create `restrained_openmm.py`: optional dependency loader, protocol constants, geometry/QC primitives, OpenMM topology preparation, staged minimization, and candidate mmCIF writing.
- Modify `tests/test_compound_contact_pipeline.py`: configuration, orchestration, output, backward-compatibility, and CLI tests.
- Create `tests/test_restrained_openmm.py`: pure geometry, acceptance, dependency-boundary, and optional real-backend tests.
- Modify `datasets.example.tsv`: add optional `ligand_smiles` and `ligand_file` columns.
- Create `environment-openmm.yml`: reproducible optional conda-forge environment.
- Modify `README.md`: v2.0 usage, protocol, interpretation, troubleshooting, and cryo-EM/Phenix guidance.

---

### Task 1: Versioned Chemistry and CLI Configuration

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`
- Modify: `datasets.example.tsv`

**Interfaces:**
- Consumes: existing `DatasetSpec`, `StructureJob`, `load_manifest`, `build_structure_jobs`, `build_parser`, and `validate_args`.
- Produces: `__version__ = "2.0.0"`; `DatasetSpec.ligand_smiles`; `DatasetSpec.ligand_file`; `StructureJob.ligand_smiles`; `StructureJob.ligand_file`; `validate_minimization_configuration(args, jobs) -> None`.

- [ ] **Step 1: Add failing manifest chemistry tests**

Add focused tests proving that chemistry is optional for contact-only runs,
that one of SMILES/SDF is mandatory for minimization, that both are rejected,
that SDF paths cannot escape `ROOT`, and that one compound cannot carry
conflicting chemistry:

```python
def test_minimization_requires_exact_ligand_chemistry(tmp_path):
    root = tmp_path / "root"
    touch(root / "C19" / "seed-1" / "model.cif")
    jobs = ccp.build_structure_jobs(root, None, "WXI", None, "B", None, None)
    args = ccp.build_parser().parse_args(
        [str(root), "--ligand", "WXI", "--minimize", "openmm"]
    )
    with pytest.raises(ccp.PipelineError, match="ligand chemistry"):
        ccp.validate_minimization_configuration(args, jobs)


def test_manifest_rejects_smiles_and_ligand_file_together(tmp_path):
    root = tmp_path / "root"
    (root / "C19").mkdir(parents=True)
    sdf = root / "C19.sdf"
    sdf.write_text("fixture", encoding="utf-8")
    manifest = tmp_path / "datasets.tsv"
    manifest.write_text(
        "dataset_dir\tcompound\tligand_resname\tligand_smiles\tligand_file\n"
        "C19\tC19\tWXI\tCC\tC19.sdf\n",
        encoding="utf-8",
    )
    with pytest.raises(ccp.PipelineError, match="exactly one"):
        ccp.load_manifest(root, manifest)
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run:

```bash
python -m pytest -q \
  tests/test_compound_contact_pipeline.py::test_minimization_requires_exact_ligand_chemistry \
  tests/test_compound_contact_pipeline.py::test_manifest_rejects_smiles_and_ligand_file_together
```

Expected: FAIL because the extended signature/fields and validation function do
not exist.

- [ ] **Step 3: Extend immutable dataset and job models**

Add nullable fields and propagate them through `_job_from_spec` and both job
discovery modes:

```python
@dataclass(frozen=True)
class DatasetSpec:
    dataset_dir: Path
    dataset_rel: str
    compound: str
    ligand_resname: str | None
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None
    ligand_smiles: str | None = None
    ligand_file: Path | None = None


@dataclass(frozen=True)
class StructureJob:
    cif_path: Path
    structure_id: str
    dataset_rel: str
    compound: str
    ligand_resname: str | None
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None
    ligand_smiles: str | None = None
    ligand_file: Path | None = None
```

Extend `build_structure_jobs` with `default_ligand_smiles` and
`default_ligand_file`. Resolve manifest ligand files with `resolve_inside` and
compare `(ligand_resname, protein_chains, ligand_chain, ligand_smiles,
ligand_file)` for per-compound consistency.

- [ ] **Step 4: Add CLI flags and cross-option validation**

Add:

```python
__version__ = "2.0.0"

parser.add_argument("--version", action="version", version=__version__)
parser.add_argument("--minimize", choices=("openmm",))
parser.add_argument(
    "--analysis-source", choices=("original", "minimized"), default="original"
)
parser.add_argument("--ligand-smiles")
parser.add_argument("--ligand-file", type=Path)
parser.add_argument(
    "--openmm-platform",
    choices=("CPU", "CUDA", "OpenCL", "Reference"),
    default="CPU",
)
parser.add_argument("--openmm-device-index")
parser.add_argument("--minimize-ph", type=float, default=7.4)
parser.add_argument("--pocket-radius", type=float, default=6.0)
parser.add_argument("--max-backbone-rmsd", type=float, default=0.5)
parser.add_argument("--max-ligand-rmsd", type=float, default=1.5)
parser.add_argument("--minimization-tolerance", type=float, default=10.0)
```

`validate_args` rejects non-finite/nonpositive numeric values, minimized
analysis without minimization, minimization with workers other than one, device
index on non-GPU platforms, manifest plus global ligand chemistry, and global
SMILES plus global ligand file. `validate_minimization_configuration` checks
that every job has exactly one chemistry source and that manifest-free global
chemistry is applied to only one discovered compound.

- [ ] **Step 5: Run all configuration tests and verify GREEN**

Run:

```bash
python -m pytest -q tests/test_compound_contact_pipeline.py -k \
  "manifest or minimization or cli_rejects or version"
```

Expected: PASS with no warnings.

- [ ] **Step 6: Update the example manifest and test its header**

Set the header to:

```tsv
dataset_dir	compound	ligand_resname	protein_chains	ligand_chain	ligand_smiles	ligand_file	enabled
```

Use a SMILES example in one row and a root-relative SDF example in another,
leaving the mutually exclusive field blank. Add a test that parses both rows
after creating the referenced fixture SDF.

- [ ] **Step 7: Commit the configuration slice**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py datasets.example.tsv
git commit -m "feat: validate minimization chemistry inputs"
```

---

### Task 2: Dependency-Free Geometry and Acceptance QC

**Files:**
- Create: `restrained_openmm.py`
- Create: `tests/test_restrained_openmm.py`

**Interfaces:**
- Consumes: NumPy only at import time.
- Produces: `MinimizationConfig`; `BackendMinimizationResult`; `OverlapMetrics`; `AcceptanceDecision`; `kcal_per_A2_to_kj_per_nm2(value: float) -> float`; `vdw_overlap_metrics(protein_coordinates_A, protein_elements, ligand_coordinates_A, ligand_elements, threshold_A=0.4) -> OverlapMetrics`; `align_on_backbone(reference_backbone_A, mobile_backbone_A, reference_ligand_A, mobile_ligand_A) -> tuple[float, float]`; `select_pocket_residues(protein_coordinates_A, protein_residue_ids, ligand_coordinates_A, radius_A) -> frozenset[str]`; `evaluate_acceptance(...) -> AcceptanceDecision`; and `candidate_relative_path(structure_id: str, accepted: bool) -> Path`.

- [ ] **Step 1: Write failing force conversion and protocol tests**

```python
def test_kcal_per_A2_converts_to_kj_per_nm2():
    assert rom.kcal_per_A2_to_kj_per_nm2(1.0) == pytest.approx(418.4)


def test_protocol_has_three_bounded_stages():
    config = rom.MinimizationConfig()
    assert [stage.max_iterations for stage in config.stages] == [500, 1000, 1000]
    assert config.protein_forcefield == "amber/protein.ff14SB.xml"
    assert config.ligand_forcefield == "gaff-2.2.20"
```

- [ ] **Step 2: Run those tests and verify RED**

Run: `python -m pytest -q tests/test_restrained_openmm.py -k "converts or protocol"`

Expected: FAIL because `restrained_openmm` does not exist.

- [ ] **Step 3: Add protocol dataclasses and constants**

Implement immutable types with these stable fields:

```python
@dataclass(frozen=True)
class MinimizationStage:
    name: str
    max_iterations: int
    backbone_k_kcal_A2: float
    distant_k_kcal_A2: float
    pocket_sidechain_k_kcal_A2: float
    ligand_k_kcal_A2: float


@dataclass(frozen=True)
class MinimizationConfig:
    protocol: str = "openmm_restrained_v1"
    protein_forcefield: str = "amber/protein.ff14SB.xml"
    ligand_forcefield: str = "gaff-2.2.20"
    ph: float = 7.4
    pocket_radius_A: float = 6.0
    tolerance_kj_mol_nm: float = 10.0
    platform: str = "CPU"
    device_index: str | None = None
    max_backbone_rmsd_A: float = 0.5
    max_ligand_rmsd_A: float = 1.5
    stages: tuple[MinimizationStage, ...] = DEFAULT_STAGES
```

Define `DEFAULT_STAGES` exactly as:

```python
DEFAULT_STAGES = (
    MinimizationStage("hydrogen_relaxation", 500, 10.0, 10.0, 10.0, 10.0),
    MinimizationStage("pocket_relaxation", 1000, 10.0, 10.0, 0.0, 1.0),
    MinimizationStage("gentle_relaxation", 1000, 2.0, 2.0, 0.0, 0.2),
)
```

Use `kcal_per_A2_to_kj_per_nm2(value) = value * 4.184 * 100.0`.

- [ ] **Step 4: Write failing overlap, alignment, pocket, path, and acceptance tests**

Cover the exact 0.4 A overlap boundary, unsupported elements, rigidly rotated
and translated coordinates, ligand RMSD after protein alignment, a residue at
the inclusive 6.0 A pocket boundary, safe deterministic output paths, accepted
metrics, and each rejection reason:

```python
def test_vdw_overlap_boundary_is_inclusive():
    result = rom.vdw_overlap_metrics(
        np.array([[0.0, 0.0, 0.0]]), ["C"],
        np.array([[3.0, 0.0, 0.0]]), ["C"],
    )
    assert result.clash_pair_count == 1
    assert result.maximum_overlap_A == pytest.approx(0.4)


def test_acceptance_rejects_ligand_drift():
    decision = rom.evaluate_acceptance(
        backbone_rmsd_A=0.1,
        ligand_rmsd_A=1.6,
        before_fixed_clashes=2,
        after_fixed_clashes=0,
        before_vdw_clashes=3,
        after_vdw_clashes=0,
        max_backbone_rmsd_A=0.5,
        max_ligand_rmsd_A=1.5,
    )
    assert not decision.accepted
    assert decision.reasons == ("ligand RMSD 1.600 A exceeds 1.500 A",)
```

- [ ] **Step 5: Run the geometry tests and verify RED**

Run: `python -m pytest -q tests/test_restrained_openmm.py -k "overlap or alignment or pocket or path or acceptance"`

Expected: FAIL on missing functions.

- [ ] **Step 6: Implement pure QC functions**

Use an explicit heavy-element radius table (`C=1.70`, `N=1.55`, `O=1.52`,
`F=1.47`, `P=1.80`, `S=1.80`, `Cl=1.75`, `Br=1.85`, `I=1.98`) and reject an
unknown element. Implement Kabsch alignment for row-vector coordinates and
apply the protein-derived rotation/translation to the ligand. Return immutable
results and stable, ordered rejection messages. `candidate_relative_path`
must replace `.cif` with `_minimized.cif`, reject absolute/parent-traversal
structure IDs, and prefix `minimized_structures/accepted` or
`minimized_structures/rejected`.

- [ ] **Step 7: Run all pure QC tests and verify GREEN**

Run: `python -m pytest -q tests/test_restrained_openmm.py`

Expected: PASS except the explicitly skipped optional OpenMM integration test.

- [ ] **Step 8: Commit the QC slice**

```bash
git add restrained_openmm.py tests/test_restrained_openmm.py
git commit -m "feat: add minimization geometry quality control"
```

---

### Task 3: Lazy OpenMM Minimization Backend

**Files:**
- Modify: `restrained_openmm.py`
- Modify: `tests/test_restrained_openmm.py`

**Interfaces:**
- Consumes: `MinimizationConfig` and the exact structure/ligand selections in `MinimizationRequest`.
- Produces: `MinimizationRequest`, `load_openmm_stack`, `minimize_structure(request, config, cache_path) -> BackendMinimizationResult`.

- [ ] **Step 1: Write a failing lazy-dependency test**

Expose an importer seam solely for deterministic dependency testing:

```python
def test_dependency_loader_names_the_optional_environment():
    def missing(name):
        raise ModuleNotFoundError(name)

    with pytest.raises(rom.MinimizationDependencyError, match="environment-openmm.yml"):
        rom.load_openmm_stack(import_module=missing)
```

Also test that importing `restrained_openmm` itself does not add `openmm`,
`openff`, or `openmmforcefields` to `sys.modules`.

- [ ] **Step 2: Run the dependency tests and verify RED**

Run: `python -m pytest -q tests/test_restrained_openmm.py -k dependency`

Expected: FAIL because the loader and exception do not exist.

- [ ] **Step 3: Implement lazy imports and request/result contracts**

Define:

```python
@dataclass(frozen=True)
class MinimizationRequest:
    cif_path: Path
    output_path: Path
    structure_id: str
    compound: str
    ligand_resname: str
    protein_chains: tuple[str, ...] | None
    ligand_chain: str | None
    ligand_smiles: str | None
    ligand_file: Path | None


@dataclass(frozen=True)
class BackendMinimizationResult:
    output_path: Path
    initial_potential_kcal_mol: float
    final_potential_kcal_mol: float
    backbone_rmsd_A: float
    ligand_rmsd_A: float
    before_vdw_clash_pair_count: int
    after_vdw_clash_pair_count: int
    before_maximum_vdw_overlap_A: float
    after_maximum_vdw_overlap_A: float
    package_versions: tuple[tuple[str, str], ...]
```

`load_openmm_stack` imports `openmm`, `openmm.app`, `openmm.unit`,
`openmmforcefields.generators`, and `openff.toolkit` only inside the function,
and translates missing imports into one installation message.

- [ ] **Step 4: Write failing backend boundary tests**

Using a small fake OpenMM stack at the loader boundary, test:

- exactly one ligand residue is required;
- only selected protein chains and the ligand survive topology filtering;
- ligand chemistry is loaded from exactly one source;
- all original heavy atoms are restrained in stage 1;
- stage 2/3 restraint group constants match the protocol;
- restraint forces use force group 31;
- physical energy queries exclude force group 31;
- transient hydrogens are omitted from output; and
- non-finite energy/positions become `MinimizationBackendError`.

The fake stack must execute the module's real selection, restraint-assignment,
identity-mapping, finite-value, and output-position functions; only third-party
OpenMM/OpenFF objects are substituted.

- [ ] **Step 5: Run backend boundary tests and verify RED**

Run: `python -m pytest -q tests/test_restrained_openmm.py -k backend`

Expected: FAIL because backend helpers and `minimize_structure` are missing.

- [ ] **Step 6: Implement chemistry loading and topology preparation**

Load isomeric SMILES with undefined stereochemistry disallowed, or require one
molecule from SDF. Parse with `PDBxFile`, delete unselected residues, verify one
ligand residue, register the molecule with `SystemGenerator`, add hydrogens at
the configured pH, then create a nonperiodic system with hydrogen-bond
constraints. Translate template mismatch and missing-heavy-atom errors into a
message that identifies the structure and chemistry source.

- [ ] **Step 7: Implement the three minimization stages**

Create one `CustomExternalForce` with per-particle `k`, `x0`, `y0`, and `z0`,
using `k*((x-x0)^2+(y-y0)^2+(z-z0)^2)`. Assign group 31. Update per-particle
constants before each stage, call `updateParametersInContext`, then call
`LocalEnergyMinimizer.minimize` with the configured tolerance and stage limit.
Use `VerletIntegrator(0.001 ps)` and the selected platform/device properties.

Query initial/final energy from group 0 only and convert to kcal/mol. Verify all
values and positions are finite. Map final positions back to the original
selected atom identities, omit newly added hydrogens, and write with
`PDBxFile.writeFile(..., keepIds=True)`.

- [ ] **Step 8: Calculate backend geometry QC and verify GREEN**

Calculate pocket membership from original heavy atoms, Kabsch-align on matched
backbone heavy atoms, calculate ligand RMSD with that transform, and calculate
before/after van der Waals overlaps. Run:

```bash
python -m pytest -q tests/test_restrained_openmm.py -k \
  "dependency or backend or restraint or energy or transient"
```

Expected: PASS with no optional packages installed.

- [ ] **Step 9: Commit the backend slice**

```bash
git add restrained_openmm.py tests/test_restrained_openmm.py
git commit -m "feat: add restrained OpenMM backend"
```

---

### Task 4: Minimization Orchestration and Analysis-Source Selection

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: `restrained_openmm.MinimizationRequest`, `MinimizationConfig`, `BackendMinimizationResult`, `evaluate_acceptance`, and `candidate_relative_path`.
- Produces: `MinimizationRecord`; `run_minimization_jobs(jobs, original_results, args, temporary_root, minimize_one=None) -> list[MinimizationRecord]`; `select_analysis_results(jobs, original_results, minimization_records, args) -> list[StructureResult]`.

- [ ] **Step 1: Write failing orchestration tests with an injected backend**

Create an injected backend that writes a real synthetic candidate CIF and
returns deterministic energy/RMSD/overlap values. Test accepted, rejected,
backend-failed, and originally-invalid jobs separately. Assert that original
input bytes are unchanged.

```python
def test_rejected_minimization_is_excluded_from_minimized_analysis(tmp_path):
    records = ccp.run_minimization_jobs(
        jobs,
        original_results,
        args,
        tmp_path / "staging",
        minimize_one=fake_drifting_backend,
    )
    selected = ccp.select_analysis_results(jobs, original_results, records, args)
    assert selected[0].status == "invalid"
    assert selected[0].qc[-1].stage == "minimization_acceptance"
```

- [ ] **Step 2: Run orchestration tests and verify RED**

Run: `python -m pytest -q tests/test_compound_contact_pipeline.py -k "minimization_is or minimized_analysis"`

Expected: FAIL on missing record/functions.

- [ ] **Step 3: Implement records, config construction, and sequential loop**

Define the orchestration record with nullable fields for unattempted/failed
jobs:

```python
@dataclass(frozen=True)
class MinimizationRecord:
    job: StructureJob
    status: str
    reasons: tuple[str, ...]
    staged_cif_path: Path | None
    output_relative_path: Path | None
    candidate_result: StructureResult | None
    initial_potential_kcal_mol: float | None
    final_potential_kcal_mol: float | None
    backbone_rmsd_A: float | None
    ligand_rmsd_A: float | None
    before_vdw_clash_pair_count: int | None
    after_vdw_clash_pair_count: int | None
    before_maximum_vdw_overlap_A: float | None
    after_maximum_vdw_overlap_A: float | None
    contact_set_jaccard: float | None
    config: MinimizationConfig
    package_versions: tuple[tuple[str, str], ...]
```

The serializer derives original/final fixed clash, minimum-distance, and
contact-count values from the original and candidate `StructureResult` objects.

Only originally valid jobs are attempted. Convert all backend exceptions into
per-structure records and `QCRecord(stage="minimization", severity="error")`;
re-raise under `--fail-fast`. Analyze each completed candidate with a
`StructureJob` whose `cif_path` is temporary but whose `structure_id` is
unchanged, then call `evaluate_acceptance`.

- [ ] **Step 4: Implement explicit source selection**

For `original`, return original results untouched. For `minimized`, return the
candidate analysis only for accepted records. Convert rejected, failed, and
not-attempted records into invalid `StructureResult` objects carrying the
original job and a precise minimization QC stage. Never fall back silently to
original coordinates.

- [ ] **Step 5: Verify orchestration GREEN and old behavior unchanged**

Run:

```bash
python -m pytest -q tests/test_compound_contact_pipeline.py -k \
  "minimization or analysis_source or end_to_end"
python -m pytest -q tests/test_compound_contact_pipeline.py
```

Expected: all existing and new tests PASS.

- [ ] **Step 6: Commit the orchestration slice**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: orchestrate minimized coordinate analysis"
```

---

### Task 5: Minimized Files, Tables, and Before/After Plot

**Files:**
- Modify: `compound_contact_pipeline.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: `MinimizationRecord` and the existing protected output inventory.
- Produces: `minimization_records_to_rows`; `plot_minimization_before_after`; updated `jobs_to_rows`, `structure_results_to_rows`, `_core_table_specs`, and `run_pipeline`.

- [ ] **Step 1: Write failing serialization and plotting tests**

Assert a stable `minimization_summary.tsv` schema, blank numeric fields for
failed jobs, six-decimal formatting, original/final source paths in
`run_manifest.tsv`, `analysis_source` in `structure_summary.tsv`, and a nonempty
PNG from accepted/rejected test records.

Required minimization columns are:

```text
structure_id, compound, status, reason, protocol, protein_forcefield,
ligand_forcefield, platform, device_index, ph, pocket_radius_A,
tolerance_kj_mol_nm, output_cif, initial_potential_kcal_mol,
final_potential_kcal_mol, delta_potential_kcal_mol, backbone_rmsd_A,
ligand_rmsd_A, before_minimum_pair_distance_A,
after_minimum_pair_distance_A, before_fixed_clash_pair_count,
after_fixed_clash_pair_count, before_vdw_clash_pair_count,
after_vdw_clash_pair_count, before_maximum_vdw_overlap_A,
after_maximum_vdw_overlap_A, before_contact_residue_count,
after_contact_residue_count, contact_set_jaccard, package_versions
```

- [ ] **Step 2: Run output tests and verify RED**

Run: `python -m pytest -q tests/test_compound_contact_pipeline.py -k "minimization_summary or before_after or analysis_path"`

Expected: FAIL because serializers and plot do not exist.

- [ ] **Step 3: Implement serialization and plot**

Serialize package versions as sorted `name=version` pairs separated by `;`.
Build a four-panel figure showing paired fixed-clash counts, paired van der
Waals clash counts, backbone/ligand RMSD with threshold lines, and same-system
potential-energy change. Use compound colors, distinguish rejected candidates,
and label energy as `Same-system potential-energy change (kcal/mol)`.

- [ ] **Step 4: Integrate staging publication with inventory protection**

Create one run-scoped `TemporaryDirectory`, run minimization before output
publication, then call `prepare_output_directory`. For each completed record,
derive an accepted/rejected relative path with `candidate_relative_path`, claim
it with `claim_output_path`, create parents, and copy the staged candidate with
`shutil.copy2`. Add the file to `generated_files.json`.

Append minimization QC to `qc.tsv`. Add `minimization_summary.tsv` and the plot
only when minimization is requested. When original analysis is selected, any
failed/rejected requested minimization returns status 1 after writing complete
original tables; minimized analysis with no accepted valid structures raises
`PipelineError("no valid minimized structures")`.

- [ ] **Step 5: Write and pass an end-to-end injected-backend test**

Call `run_pipeline(args, minimize_one=fake_backend)` directly. Verify core
tables, minimization table/plot, accepted and rejected CIF paths, inventory
membership, unchanged originals, original/minimized source semantics, and exit
codes 0/1/2.

Run:

```bash
python -m pytest -q tests/test_compound_contact_pipeline.py -k \
  "minimization or inventory or end_to_end"
```

Expected: PASS.

- [ ] **Step 6: Commit the reporting slice**

```bash
git add compound_contact_pipeline.py tests/test_compound_contact_pipeline.py
git commit -m "feat: report minimization outcomes"
```

---

### Task 6: Environment, Documentation, Optional Integration, and Full Verification

**Files:**
- Create: `environment-openmm.yml`
- Modify: `README.md`
- Modify: `tests/test_restrained_openmm.py`
- Modify: `tests/test_compound_contact_pipeline.py`

**Interfaces:**
- Consumes: completed v2.0 CLI and backend.
- Produces: reproducible installation instructions, optional real OpenMM smoke coverage, and final verified release artifacts.

- [ ] **Step 1: Write a failing environment-contract test**

```python
def test_openmm_environment_contains_protocol_dependencies():
    text = Path("environment-openmm.yml").read_text(encoding="utf-8")
    for dependency in (
        "openmm", "openmmforcefields", "openff-toolkit", "rdkit",
        "ambertools", "biopython", "numpy", "matplotlib", "pytest",
    ):
        assert dependency in text
```

- [ ] **Step 2: Run it and verify RED**

Run: `python -m pytest -q tests/test_restrained_openmm.py::test_openmm_environment_contains_protocol_dependencies`

Expected: FAIL because `environment-openmm.yml` does not exist.

- [ ] **Step 3: Create the optional conda environment**

Create a conda-forge-only YAML with Python 3.11, the dependencies in Step 1,
and no unpinned alternative channels. Keep protocol force-field identifiers in
the program/summary even if package patch versions later change.

- [ ] **Step 4: Add the optional real-backend integration test**

Mark the test with `pytest.mark.skipif` unless all optional modules are
available. Use a tiny complete protein-ligand mmCIF fixture with explicit ligand
bonds and exact chemistry, run `minimize_structure` on CPU, and assert:

```python
assert result.output_path.stat().st_size > 0
assert math.isfinite(result.initial_potential_kcal_mol)
assert math.isfinite(result.final_potential_kcal_mol)
assert result.after_vdw_clash_pair_count <= result.before_vdw_clash_pair_count
```

If the optional stack cannot be installed in the current execution environment,
the skip is expected; the dependency boundary and fake-stack backend tests
remain mandatory and cannot be skipped.

- [ ] **Step 5: Rewrite README usage and scientific interpretation**

Document base versus OpenMM installs, version check, manifest fields, one exact
C19-style command, all protocol constants, accepted/rejected paths, table
columns, analysis-source behavior, CPU/GPU platform selection, dependency and
topology-mismatch troubleshooting, and the distinction between minimization
energy and binding free energy. Include the official OpenMM minimizer,
positional-restraint, and openmmforcefields links from the specification.

Add a cryo-EM note: for the user's C19/NSD2 or similar map-supported models,
use Phenix real-space refinement with the experimental map and ligand restraint
CIF after geometry preparation; do not substitute map-free OpenMM minimization
for final experimental refinement.

- [ ] **Step 6: Run focused documentation/configuration checks**

Run:

```bash
python -m pytest -q tests/test_restrained_openmm.py \
  tests/test_compound_contact_pipeline.py -W error
python compound_contact_pipeline.py --version
python compound_contact_pipeline.py --help
python -m py_compile compound_contact_pipeline.py restrained_openmm.py
git diff --check
```

Expected: all mandatory tests PASS, optional real OpenMM test either PASS or is
explicitly SKIPPED, version prints `2.0.0`, help lists new flags, compilation
succeeds, and the diff has no whitespace errors.

- [ ] **Step 7: Verify the lazy-import boundary in a base environment**

Run a Python subprocess that imports `compound_contact_pipeline`, asserts
`openmm`, `openff`, and `openmmforcefields` are absent from `sys.modules`, and
runs a contact-only synthetic CLI example. Expected: success without optional
packages.

- [ ] **Step 8: Review scientific labels and inspect generated plot**

Search all source/docs for `binding energy`, `affinity`, and `potential energy`.
Confirm no minimization metric is routed into imported energy summaries. Render
and inspect `minimization_before_after.png` for legibility, clipped labels,
threshold lines, and correct accepted/rejected encoding.

- [ ] **Step 9: Run the complete regression suite and commit**

```bash
python -m pytest -q -W error
git diff --check
git status --short
git add README.md environment-openmm.yml tests compound_contact_pipeline.py restrained_openmm.py datasets.example.tsv
git commit -m "feat: release restrained minimization workflow"
```

Expected: every mandatory test passes with warnings treated as errors and only
the optional real-backend test may be skipped when dependencies are absent.
