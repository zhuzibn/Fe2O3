# Fe2O3 deterministic regression and focused checks

This suite protects the existing Class 4 fixed-length LLG/LLGS implementation.
It does not validate hematite parameters, crystallographic mapping, the complete
sixth-order energy, or the accuracy of the legacy solvers.

The user accepted **cq-20261120-v1** on 2026-11-20. The tracked
`baseline_manifest.json` pins the unchanged repeat-1 suite under ignored
`baseline/cq-20261120-v1/`. The root API now uses validated configuration,
a versioned result and protected output. Full-duration/physical validation
remains outside this short suite.

## Run focused checks

From the repository root, in MATLAB:

```matlab
addpath('tests');
run_cpu_checks
run_gpu_checks
run_interface_checks
check_benchmark_contract
```

CPU checks require MATLAB. GPU/interface checks additionally require Parallel Computing
Toolbox and a supported GPU. The tested environment is MATLAB R2024b Update 1
on Windows. No Octave compatibility is claimed. Hash/path helpers use the
JVM supplied with MATLAB; no additional package installation is needed.

The CPU checks cover all six masks on three finite lattices, independent
direct-index neighbor matrices, reciprocity, field derivatives, guarded
anisotropy arithmetic, the retained dipole tensor and implicit-Gilbert torque
checks. GPU checks exercise RK4 stages and analytical convergence, actual
the canonical validator's time construction, sparse-versus-dense samples, chunk continuity,
spin/DMI alignment and the existing dipole cache. Heun/PC4 checks characterize
norm/finite/chunk behavior only; they do not establish formal solver order.

GPU instrumentation is isolated under ignored `.cq-checks/<unique>/`, with
`checks.log` and a disposable field wrapper. Source scripts are not modified.

Interface checks exercise complete validation, numeric classes, explicit
results/provenance, old/new restarts and plotting, output collisions including
late writers and Windows junctions, explicit replacement and shadowed sources.
They create disposable files under `candidates/` and run only tiny GPU cases.

## Run accepted regression

```matlab
addpath('benchmark');
candidate = run_regression; % Or run_regression('D:/chosen/fresh/directory')
compare_regression(candidate);
```

This verifies the pinned suite before reading fixed restart inputs, creates
six new `result.mat` files and legacy-schema comparison records from actual
new-run values, and compares all cases. Missing baseline data fail; the command
never generates or promotes a reference. Source identity includes untracked
production helpers and dirty Git status.

## Historical reference preparation

The following command was used before migration. It deliberately fails the
original-source hash guard in the migrated checkout; do not update its pins
to manufacture a replacement baseline. Its two original suites remain evidence.

```matlab
addpath('benchmark');
proposal = capture_legacy_reference;
% Optional: capture_legacy_reference('D:/chosen/fresh/directory')
metrics = compare_regression(fullfile(proposal,'repeat-2'), ...
                             fullfile(proposal,'repeat-1'));
addpath('tests');
check_benchmark_contract(proposal)
```

Through the installed WSL-to-Windows launcher:

```bash
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('benchmark'); capture_legacy_reference"
```

Preparation pins all 13 original MATLAB files to
`legacy_source_manifest.json`. It refuses changed production source. Each run
uses an isolated copy of original `main.m`; `substitutions.json` records exact
control changes and metadata hooks. The original `clear all` is isolated in a
separate function scope. All numerical kernels remain byte-identical. The
ordinary default simulation is never launched by this command.

All cases use `tstep=1e-16 s`. Other physical defaults remain unchanged unless
the case below explicitly overrides them.

| Case | Core | Chunks × steps; stride | Purpose |
|---|---|---|---|
| small-relax | 5×9×6 | 3×8; 2 | Default alternating spins, zero current |
| driven-restart | 5×9×6 | 3×8; 1 | Fixed nonuniform input, SOT 6e12 and STT 1e9 A/m², STT DL/FL enabled |
| default-core-short | 20×20×6 | 2×8; 2 | Default geometry; short numerical regression, not a 25 ps baseline |
| heun-legacy | 5×9×6 | 3×8; 1 | Same nonuniform input, Heun, zero current |
| pc4-legacy | 5×9×6 | 3×8; 1 | Same input, PC4, bootstrap and chunk transitions |
| dipole-cache | 5×9×6 | 2×52; 1 | Same input, RK4, dipoles enabled, chunk-local refresh |

The restart vector at occupied core `(r,c,k)` is the normalized
`[(-1)^(k+1),0.07*sin(r+2*c+k),0.05*cos(2*r-c+3*k)]`.
Actual input arrays are saved. Normal initialization keeps the source's exact
trigonometric residuals. Production RNG branch behavior is recorded without
claiming stochastic simulation; no thermal field is implemented.

## Files and ownership

Each fresh `benchmark/candidates/<unique>/reference-proposal/` contains:

- `status.json`: incomplete until both repeats and comparison succeed; then
  explicitly `reference-proposal-awaiting-review`, with manifest hashes.
- `repeat-1/` and `repeat-2/`: independent complete suites and manifests.
- Each suite's `original-source/`: unchanged source snapshots.
- Each case: isolated source, substitution record, original
  `legacy_workspace.mat`, explicit `record.mat`, `initial_input.mat`, and
  `Relax_0.mat` only when the case uses restart input.

The regression record schema is `1.0`: case specification, effective captured
configuration, occupancy/six masks and array order, actual initial spin
arrays, all spin/DMI samples and times, units/sampling convention and
source/runtime/RNG provenance. This benchmark record is distinct from the
production `fe2o3_result` schema described in the root README. Original legacy workspaces can contain GPU arrays;
comparison uses host arrays in `record.mat` and hashes other files as bytes.

Candidate destinations must be fresh and outside `benchmark/baseline/`.
Canonicalization resolves existing parent aliases/junctions and Windows case.
Ordinary capture/comparison never accepts data or updates baseline manifests.
`compare_regression(candidate)` requires an accepted `baseline_manifest.json`
and fails when it is absent. An explicit second path permits proposal/repeat
comparison; it does not promote that path.

Only explicit review authorizes a new accepted suite. No MAT data enter Git.
A clone needs separate transfer of the complete `baseline/cq-20261120-v1/`
directory, including supporting files. Keep its internal `manifest.json`
unchanged (its historical proposal status is intentional); acceptance is
recorded by the external tracked `baseline_manifest.json`. Verify the transfer
with `compare_regression(fullfile(pwd,'benchmark','baseline','cq-20261120-v1'))`
after `addpath('benchmark')`. Missing or modified artifacts fail; do not regenerate.
No promotion command is provided. Old local `.review-20261013/`
files and incomplete proposals are retained and are not accepted references.

## Comparison contract

Comparison requires all six cases in declared order, matching configuration,
initial arrays, mapping, units and sampling, expected shapes, and every
timestamp and saved spin/DMI sample. It checks artifact hashes and recomputes
norm/finite/empty-site invariants from raw reference and candidate arrays.
Failures raise MATLAB errors, including for missing references.

Accepted same-environment limits are `2e-12` for spin components and occupied
norm errors, `2e-10 T` for DMI samples, exact zeros at empty sites, and exact
time/input/mapping equality. These derive from the existing focused double
precision checks and accepted repeat measurements; they are not
cross-platform or full-duration tolerances. Do not widen them automatically.

`check_benchmark_contract` uses the accepted suite; an optional proposal path
uses its two historical repeats. It corrupts disposable copies, updates their
hashes where necessary to reach semantic checks, and requires the expected
error identifiers. It also checks protected paths and rechecks original hashes.
No benchmark timing claim is made: initialization/compilation/cache effects
are not separated as a performance experiment.

See [the current verification report](../docs/cq01-05-verification.md) for
actual commands, artifacts, measurements, limitations and acceptance status.
