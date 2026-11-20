# Archived CQ01–05 approved implementation plan

Completed 2026-11-20. This is the reviewed staged plan and its pre-migration
observations, retained as historical design/approval evidence. References to
pending gates and unchanged production files below describe those earlier
stages. The user subsequently accepted repeat-1 and authorized stages 5–7 by
replying “continue”. Current outcomes and limitations are in
[cq01-05-verification.md](cq01-05-verification.md).

# CQ-01–CQ-05: preserve behavior while making runs reproducible and inspectable

Revision 1 — 2026-11-20 — **APPROVED; baseline accepted and stages 5–7 authorized by the user's “continue”.**
Repository: `D:\Git_control\spin_dynamics\Fe2O3`.
Inspected `main` at `523500cf85d32462361776499937eb4ce25e7fd7`.
Stages 0–4 have passed, including explicit acceptance of repeat-1 and its
tolerances. Execution evidence is in
[docs/cq01-05-verification.md](cq01-05-verification.md). The accepted
suite is `benchmark/baseline/cq-20261120-v1/` (manifest SHA-256
`94e2b7a2ed194a0f9b2d2e5ec62774a890fa46665f0a872a5e731aec7fb13bb2`)
recorded in `benchmark/baseline_manifest.json`. Continue stages 5–7 without
another approval unless a material design/model change becomes necessary.
No commit, push, branch/worktree creation, dependency installation, or network use.

The previous root plan described the completed 2026-10-15 import. Its full
technical content and approval context already survive in
[the archived import plan](code-20261013-replacement-plan.md); the only
difference is that the archive labels completion and later PR #18 fixes.
Current limitations remain in [the methodology](model-methodology.md)
and [diff.md](../diff.md). No unfinished implementation is discarded here.

## Goal and Success Criteria

Improve all five catalog features [1] for the maintained `main.m` path and a
small callable runner, without changing the model or numerical algorithms.
The smallest coherent change is a function boundary around the existing
scripts, one configuration definition/validator, explicit result assembly,
protected persistence, and focused regression/check commands. Keep the script
kernels; do not replace them with a new simulation framework.

| Feature | Observable acceptance outcome |
|---|---|
| CQ-01 | A separately accepted, hashed, deterministic reference suite; comparison enforces every expected case, every saved spin/DMI sample, time values, input compatibility and raw-array invariants; missing references fail. |
| CQ-02 | `cfg=fe2o3_default_config(); result=fe2o3_run(cfg,output_file,overwrite)` is the documented computational entry; `main` calls it; defaults occur once; invalid or unknown inputs fail before integration. |
| CQ-03 | A versioned `fe2o3_result` contains effective inputs, derived quantities, initial state/mapping, trajectories and clocks, source/runtime provenance and known limitations; save/load round trip passes. |
| CQ-04 | Keep the script's `Relax.mat` default filename, but refuse an existing destination unless replacement is explicit; reject writes into accepted-reference storage even with overwrite requested. |
| CQ-05 | Runnable CPU and small GPU checks recover the useful historical evidence, test the new interface/output responsibilities, and distinguish numerical agreement from physical validation. |

Non-goals: new physics, physical-source reconciliation, topology changes,
independently gated torque redesign, Heun/PC4 repair, CPU production backend,
thermal noise, final-endpoint saving, new observables, GUI, full 25 ps runs,
performance optimization, or size-convergence/material claims. Existing
results and scientific defaults stay unchanged. Unknown external consumers
of arbitrary saved workspace temporaries cannot be promised compatibility;
the documented plot/restart consumers will be adapted explicitly.

**Approval stages:** after approval of this plan, implement only the focused
harness/reference preparation first. Present its actual artifacts, hashes,
errors and environment for explicit baseline acceptance. Do not refactor the
production interface until that gate passes. Plan approval is not acceptance
of data that have not yet been generated.

## Model Classification and Equation Family

Class 4 material-specific atomistic model; fixed-length LLG/LLGS [2, 3].
This is the classification of the existing implementation, not certification
of its hematite parameterization. One normalized classical Fe-labeled spin
per occupied site, one common moment, six software site masks; odd/even
layers initially point approximately along +x/-x. Those masks are not six
chemical species or a verified crystallographic sublattice mapping.
Electronic structure and transport are outside this classical model.

Inputs are explicit numerical/material/drive controls and optional restart
spins. Outputs remain sampled spin components, sampled DMI components and
time, plus configuration/provenance metadata. Invariants: finite values,
unit occupied-spin norms after updates, zero empty sites, matching sample
and clock counts. No stochastic field or variable-length dynamics applies.
Field/torque amplitudes are in tesla with gamma in rad/(s T), following the
retained implementation and shared equation notation [4].
Physical material/site/order provenance remains **REFERENCE MISSING** for
this exact implementation; no such definition is replaced in this task.

## Existing Behavior and Reference Implementations

Preimplementation source inspection, 2026-11-20, confirmed the Fe2O3 gaps in the matrix
(snapshot 2026-11-19) [1]. The matrix was read, not assumed current or edited.

| Feature | Actual source/artifact evidence | Remaining gap |
|---|---|---|
| CQ-01: Missing | No maintained benchmark directory/manifest or accepted trajectory inventory. Ignored `.review-20261013/` contains historical `*-pass.mat`, `dmi-operator.mat`, and `smoke/final.mat`. | Historical artifacts are not an accepted reference for current code. No automatic promotion or overwrite. |
| CQ-02: Partial | `main.m:3–108` clears caller state and defines controls/coefficients. Selected checks at `main.m:13–15,99–101`, `integrate_llg.m:5–9,116–118`. | No complete input API; restart shapes/masks unchecked; derived and editable values intermixed. |
| CQ-03: Missing | `main.m:134` calls unqualified `save('Relax.mat')` after cleanup. | Whole workspace, no explicit version/schema or comprehensive code/runtime identity. |
| CQ-04: Missing | Same save overwrites; restart reads fixed `Relax_0.mat`. | No destination preflight or baseline ownership guard. |
| CQ-05: Partial | Fields and RHS already separated into scripts/function; PR #18/#20 reports embed focused checks. | No maintained suite; ignored harnesses are stale, path/fixture dependent. |

Specific current evidence and conflicts:

- `git diff --exit-code cf19a99 HEAD -- '*.m'` returned 1: the only difference
  is seven explanatory DMI comments. Historical statements about exact source
  equality are dated evidence; they must not be described as a fresh pass.
- `.review-20261013/verify_structure.m` expects removed `ato_s`;
  `verify_shell_dmi_visualization.m` extracts shift/DMI blocks from
  `integrate_llg.m` that now reside in `field_calc.m`. Do not run these old
  scripts unmodified or overwrite their artifacts.
- The PR #20 report [5] preserves a recoverable RK4 harness, but hard-codes an
  old checkout and requires a disposable field-logging wrapper. Port its
  harness to maintained checks, not its historical source into production.
- `BFSTT=chi*BDSTT` uses the SOT ratio; DLT flags also gate corresponding FLT.
  Preserve these known behaviors, including metadata that describes them.
- `symbols.tex` [4] gives `ele=1.602176565e-19 C`, whereas this repository uses
  `1.602e-19 C` and recomputes gamma. Preserve repository arithmetic/constants.

ASTRO context [6]: inspected its `AGENTS.md`, `README.md`,
`astro_default_config.m`, `astro_build_production_result.m`, and
`benchmark/run_baseline_benchmark.m`. Its current candidate runner already
has fresh-path/mode protection beyond the matrix's dated description. Its
configuration/result organization is useful; its species meanings, model,
solver restrictions, default values, tolerances and baseline ownership do
not transfer to Fe2O3. No macrospin reference implementation was needed.

Preserve unrelated work: README's existing local-checkout paragraph and
untracked `Fe2O3-DMI-02-conversation-2026-10-21.md`. Only `plan.md` is changed
in this planning turn. It is already ignored by `.gitignore:2`; do not alter
that rule or stage it. No `CHANGELOG.md` currently exists.

## State Variables, Geometry, and Boundaries

Preserve component arrays `(W+4,L+4,H)` and histories `(W+4,L+4,H,Nsample)`;
MATLAB column-major linear indices are software ordering. Store requested
core dimensions separately from padded dimensions because `distrib.m`
mutates `natomW/natomL`. Retain occupancy and six masks/labels in results.
Keep the historical reversal of layer-index direction versus physical z
visible; it does not supply crystallographic coordinates.

| Physical detail | Adopted value or rule | Builder or interaction mapping | Qualifying reference and exact locator | Evidence status |
|---|---|---|---|---|
| Layer/lattice type and species | Existing three-layer repeat, common Fe-labeled moment, six masks | `distrib.m`; `systemgeneration.m` | Exact physical mapping not established by inspected local records [7, 8] | REFERENCE MISSING; freeze existing bookkeeping, no reconstruction |
| Lattice vectors, basis, spin coordinates and stacking | No verified complete coordinate table | Existing array offsets only | None established | REFERENCE MISSING; do not invent coordinates or claim a validated crystal |
| Lattice scales/thickness | `a=5.067e-10 m`, `c=13.882e-10 m`, `tz=c/6*H` | `main.m:39–42` | Exact value/normalization locators not available in inspected evidence | REFERENCE MISSING; preserve coefficients and operation order |
| J1 | `1.3740e-21 J/link`; code interior count 1; interlayer selection | All six `exchangej_.m` branches | No verified source-shell mapping | REFERENCE MISSING for material assignment |
| J2 | `-1.1696e-21 J/link`; code interior count 3; same-layer selection | Same | Same | REFERENCE MISSING for material assignment |
| J3 | `4.0413e-21 J/link`; code interior count 3; adjacent-layer selection | Same | Same | REFERENCE MISSING for material assignment |
| J4 | `2.8041e-21 J/link`; code interior count 6; adjacent-layer selections | Same | Same | REFERENCE MISSING for material assignment |
| Exchange sign, normalization and pair counting | Positive J antialigns under the code-consistent energy; unit spins, reciprocal directed pairs and energy factor 1/2 | `field_calc.m:80–82`; independent operator check planned | [4] supplies notation only; exact material convention unresolved | Preserve code; verify algebra, not literature assignment |
| DMI bonds | Actual layer parity times J1/J3/−J4 sums; D1 fixed zero | `dmi.m`, `field_calc.m:108–116` | Physical vectors/axis approximation and shell assignments not established | REFERENCE MISSING; no nonzero-D1 extension |
| Boundaries | Two empty lateral padding cells; omit out-of-range layer bonds | `distrib.m`, `exchangej_.m` | No verified physical surface termination | REFERENCE MISSING for physical interpretation; preserve finite operator |

Topology checks enumerate every site/mask, reciprocal link, interior count,
and truncated edge/surface on small lattices. These are software checks of
retained behavior. The target and verification geometry remain the existing
finite padded geometry. No periodic-first conversion is appropriate for this
migration, and no new periodic/mixed boundary selector will be introduced.
Height must be a positive multiple of three (at least three): the present
last-layer guards only terminate motif case 2. Other heights can read beyond
the last layer. This is an early input guard, not a topology repair.

Missing physical construction evidence blocks new construction/interpretation,
not preservation of the existing operators during this software-only task.

## Hamiltonian

Retain the partial energy convention already documented in [3], with
notation cross-checked against [4, `eq:atomistic-hamiltonian`]:

$$
E_{\mathrm{ex}+2+Z}=\frac12\sum_{i\ne j}J_{ij}\mathbf s_i\cdot\mathbf s_j
-K_1\sum_i s_{iz}^2-\mu_s\sum_i\mathbf B_{\mathrm{ext}}\cdot\mathbf s_i.
$$

This is total energy in joules, using both directed occurrences of a
reciprocal pair; J is energy/link and K is energy/site. For configured D1=0,
the retained code energy is $E_D=-\mathbf x^TW\mathbf y$ with
$W=\operatorname{diag}((-1)^k)(D_3A_3-D_4A_4)$ on occupied sites [3].
Check $W=-W^T$ before using it as an algebraic energy reference.

The complete sixth-order energy/axis mapping remains **UNRESOLVED / REFERENCE
MISSING**. Do not invent an energy observable or full-energy conservation
test. The optional dipole pair energy can be checked only for the existing
reciprocal retained tensor operator; physical truncation/normalization remain
unresolved. Its scalar test energy is $-\tfrac12\mu_s\sum_i\mathbf s_i\cdot
\mathbf B_{\mathrm{dip},i}$ when those algebraic conditions hold, in joules.
SOT/STT remain outside conservative energy. No thermal term exists.

## Parameters and Evidence

All values below are **current code settings** [2], not newly sourced
hematite parameters. Material provenance is **REFERENCE MISSING** unless a
qualifying locator is explicitly given; [7, 8] do not establish these exact
values. Keep expressions and their evaluation order to avoid migration drift.

| Inputs/constants | Current value and unit | Proposed ownership / dependency |
|---|---|---|
| `natomW,natomL,natomH` | `20,20,6` core counts | Editable configuration; padding derived internally |
| `J1,J2,J3,J4` | Values in the structure table, J/link | Editable configuration; same field sign |
| `Ksim1,hk6` | `-8e-22,1.76e-27 J/site` | Editable configuration; preserve existing sixth-order formula |
| `g,mus,alp` | `2,4.2313*mub J/T,0.001` | Editable configuration; positive g/moment, nonnegative damping |
| `a,c` | `5.067e-10,13.882e-10 m` | Editable configuration; positive; affect existing density/thickness expressions |
| `DMIenable,Dsim3,Dsim4` | `1,1.6e-23 J,1.28e-23 J` | Nominal inputs; effective weights zero when disabled; Dsim1 fixed at zero |
| `dipolemode` | `0` | Existing 0/1 selector; fixed 50-step chunk-local cache schedule |
| Dipole geometry | `dist_=0.5067e-9 m`; vector components `±0.5067e-9`, `±0.25335e-9`, `±0.4388e-9 m` | Existing `dipole_.m` literals, preserve rounded vectors independently of editable `a`; record in model metadata |
| `Hext` | `[0,0,0] T` | Finite three-component vector; historical name retained |
| `jcSOT,jcSTT` | `0,0 A/m²` | Finite signed currents |
| `SOT_DLT,SOT_FLT,thetaSH` | `1,1,0.2` | Binary selectors and finite signed efficiency |
| SOT angles / ratio | `psj_theta=90/180*pi`, `psj_phi=210/180*pi` rad; nominal `chi=5` | Preserve trigonometric residuals; effective chi zero if SOT_FLT=0 |
| `STT_DLT,STT_FLT,etaSTT,psjSTT` | `0,0,0.8,[0,0,1]` | Existing selectors/efficiency/unit polarization |
| `chiSTT` | Effective `0`, or `0.2` if STT_FLT=1 | Retained derived legacy value; still unused in actual BFSTT, no new ineffective editable knob |
| `rk4` | `1` | Keep 0/1/2 selections; 0/2 explicitly experimental, zero torque only |
| `tstep,gpusave,gpurun_number,savetstep` | `1e-16 s,1e-12 s,25,5000` | Positive numerical controls; enforce integral chunk/step ratio and stride divisibility |
| Distribution flags | `load_fixed_atom_distrib=0`, `save_fixed_atom_distrib=0` | Retain mutually exclusive controls; save flag documented legacy no-op, no new behavior |
| Restart filename | `Relax_0.mat` | Explicit `restart_file` input; validate content before GPU integration |
| RNG/noise | `rng('shuffle')` only in normal flag branch; no random model draws/no thermal field | Preserve existing branch behavior and record states; no new production seed/temperature option |
| Constants | `mub=9.274e-24 J/T`, `ele=1.602e-19 C`, `hbar=6.58211951440e-16 eV s`, `mu_0=1.25663706143592e-6 T m/A`, `kb=1.38064852e-23 J/K` | `constantfile.m` authoritative; do not expose editable duplicate constants; kb unused |
| Derived values | `gam=g*mub/(hbar*ele)`, `tz=c/6*H`, `ms=12*mus/c/a^2/sqrt(3)*2` | Recompute once after validation; save actual values; constants file's initial gam is superseded as today |
| Derived drives/time | `psjSHE`, component aliases, RE/TM-labelled equal amplitudes, effective chi, `gpusteps`, `runtime`, sample count and t | Internal only; reject user-supplied derived aliases |
| Output controls | `output_file=''`, `overwrite=false` for callable API | Empty means return without saving; `main` explicitly requests `Relax.mat` |

All required configuration fields must be present; no silent merge of partial
structs. Reject unknown fields, complex/nonfinite values, incorrect shapes,
nonbinary flags, invalid integer controls and unsupported selector values.
Validate all derived values for finite, consistent outputs. Ratios within
`32*eps(max(1,abs(ratio)))` of an integer are floating-point integral;
otherwise reject rather than silently round timing. This proposed guard keeps
default arithmetic and excludes previously unvalidated nonintegral timing.
No arbitrary literature-derived limits on signed J, K, drive or efficiency.
Polarization must be unit to `1e-12`; never silently normalize user inputs.

## Effective Fields and Units Check

Preserve $\mathbf B_i=-\mu_s^{-1}\partial E/\partial\mathbf s_i$ [3; 4,
`eq:atomistic-effective-field`, `eq:atomistic-pair-fields`]. Differentiating
both reciprocal pair occurrences cancels 1/2, giving
$\mathbf B_{\mathrm{ex},i}=-\mu_s^{-1}\sum_jJ_{ij}\mathbf s_j$;
$\mathbf B_{2,i}=2K_1s_{iz}\hat{\mathbf z}/\mu_s$ and the Zeeman field is
$\mathbf B_{\mathrm{ext}}$. Joules divided by J/T gives tesla; no extra
$\mu_0$ belongs here. These are algebraic checks of the preserved convention.

For DMI, preserve $\mathbf B_D=(S_{D,y},-S_{D,x},0)/\mu_s$; use the
antisymmetric occupied-site W above for finite-difference checks. The entire
sixth-order field is unchanged from `field_calc.m:83–107` and [3, Anisotropy]:
normalized stage angles, guarded poles/empty sites, and unnormalized stage z
in the second-order term. Link that full expression rather than rederive an
unsupported energy. Arithmetic/pole/regression tests apply to hk6; no claim
of complete energy validation.

Retain the dipole tensor and its cache from [3, Dipole approximation], with
$\mu_0\mu_s/(4\pi r^3)$ in tesla and the existing rounded vectors/distances.
Do not replace their scalar distance with a recomputed vector norm. Record
this as a legacy implementation choice, not a verified crystal metric.

## Dynamical Equation and Nonconservative Torques

The unchanged RK4 derivative is [2, `atomgpurk4.m`; 4, `eq:ll-with-flt`]:

$$
\dot{\mathbf s}=\frac{\gamma}{1+\alpha^2}\left[-\mathbf s\times\mathbf B
-\alpha\mathbf s\times(\mathbf s\times\mathbf B)
+\sum_{a=\mathrm{SOT,STT}}\left\{(\alpha B_{D,a}+B_{F,a})\mathbf s\times\boldsymbol\sigma_a
+(-B_{D,a}+\alpha B_{F,a})\mathbf s\times(\mathbf s\times\boldsymbol\sigma_a)\right\}\right].
$$

Gamma is rad/(s T); drive amplitudes are T. Preserve all damping cross terms,
signs and the uniform application to occupied sites. SOT/STT here label two
polarization-torque channels, not a spatial transport model.
With `hbar` numerically in eV s, retain the existing amplitudes
`DLT*hbar/2*efficiency*jc/(ms*tz)`, and both actual FL amplitudes equal effective
`chi` times the associated DL amplitude. Preserve `chiSTT`'s unused status.
Nonzero effective amplitudes continue to require RK4, regardless of enabled
flags. No torque equation conversion, microscopic claim or LLB departure.
General material applicability of these coefficients remains REFERENCE MISSING.

## Initial Conditions, Solver, and Outputs

Keep `systemgeneration.m` initialization including floating-point sine/cosine
residuals, six masks and empty zeros. Restart continues to use the last saved
frame, not the unsaved final integration state. Accept legacy
`mmx/mmy/mmz` MAT input and the new explicit result schema; require matching
component histories, padded dimensions, finite values, occupied unit norms
within `2e-12`, and zero empty entries. Reject incompatible input without
repairing or silently rescaling it. New-schema input additionally checks
stored mapping. Legacy input is labelled as lacking that provenance.

RK4 recomputes exchange/anisotropy/DMI at all four substages; dipoles are
cached. Heun/PC4 retain fixed fields, and PC4's first three updates use RK4.
No kernel, normalization, cache schedule or temporal sample change.
Defaults still perform 250,000 updates over 25 ps and save 50 pre-update
frames at 0:0.5:24.5 ps. Store requested/effective duration and sample policy
separately so the unsaved endpoint is explicit.

Proposed schema `fe2o3_result`, version `1.0`:

- `configuration`: exact validated inputs; `derived`: actual constants,
  coefficients, padded dimensions, update/sample counts and solver/cache settings.
- `geometry`: occupancy, all six site masks/labels, MATLAB ordering and boundary
  description, with unresolved physical coordinate mapping labelled.
- `initial_state`: actual x/y/z arrays; initialization mode and restart path/hash
  when used. Store arrays even for deterministic initialization.
- `trajectory`: existing `mmx/mmy/mmz`, `hdmi_allx/y/z`, matching `t`, units,
  `(W,L,H,time)` order and pre-update/endpoint policy. CPU arrays only.
- `provenance`: UTC creation time, exact entry/arguments, Git HEAD/branch/status,
  SHA-256 of all executed project M-files (including untracked/new files),
  MATLAB/platform/toolbox/GPU/driver details, elapsed time, RNG before/after
  initialization/run and algorithm. Record real lookup errors/unavailable
  fields; do not hard-code “no Git”. Verify project function resolution so a
  shadowed kernel cannot be attributed to the wrong source. If files change
  during a run, reject reproducibility claims and report it explicitly.
- `limitations`: retained model-source gaps, experimental solvers when selected,
  torque gating, cached dipoles, and unsaved endpoint. No fabricated energy data.

Save only this named struct. Make `main.m` a thin caller and adapt
`One_data_curve.m` to unpack the required fields when that struct is present,
while retaining its old workspace-input path. Keep its plots/site indices
unchanged. `systemgeneration.m` gets only explicit restart-path/schema
handling, not an initialization rewrite. Document this intentional workspace
format migration and the two supported consumer paths.

Destination validation precedes initialization/GPU allocation when feasible;
recheck immediately before saving. Use an exclusive destination reservation
or equivalent non-clobbering finalization, so another file appearing during
integration is not overwritten. Existing files/directories, input/output
identity and canonical paths inside accepted-baseline storage are rejected;
resolve `..` and path aliases, with Windows case handling. `overwrite=true`
applies only to ordinary outputs, never accepted references or restart input.
Use distinct files for every test/candidate. No append output exists here.

## Implementation Mapping

The planning turn wrote only this plan. Stages 0–4 are now implemented and
verified; stages 5–7 remain pending baseline acceptance. Use the existing
architecture in this order, with no stage proceeding until its gate passes.
The verification-first stages justify the variation from new-model construction:
structure/topology/fields/RHS/integrator are retained and tested before any
production-interface migration.

| Stage | Implementation | Verification gate |
|---|---|---|
| 0. Protect and preflight | Record source hashes, working-tree changes, runtime; prepare fresh ignored scratch/candidate paths | No existing artifact changed; MATLAB/PCT/GPU suitability established; source revision pinned |
| 1. Structure and topology checks | Recover CPU mask/basis-probe checks for current scripts; independent expected offset table from the existing guide | All six masks, parity, shell reciprocity/counts and finite boundary omissions agree; physical mapping stays unresolved |
| 2. Static fields and RHS checks | Recover term/pole/DMI checks and vector LLGS tests | Independent analytical/finite-difference/implicit-Gilbert checks pass within stated scope; no physics repairs hidden in test recovery |
| 3. Integrator checks | Recover PR #20 stage and analytic-convergence harness; sampling/cache/legacy characterization | Current RK4 passes independent checks; legacy behaviors characterized, not certified; no full default simulation |
| 4. Freeze pre-migration behavior | Capture six small cases using isolated copies of original code; comparator and mutation tests; prepare proposal manifest | Repeatability and all comparison/error tests pass; user accepts concrete artifacts and tolerances before stage 5 |
| 5. Configuration boundary | One defaults function, validator and callable runner; `main` becomes caller; minimal restart/plot adapters | Config/invalid-input tests; original-vs-new effective inputs, masks, initial state and all six accepted trajectories match |
| 6. Results and persistence | Explicit schema/provenance and preflight/non-clobbering save in runner | Save/load consistency, source identity, bad destination/overwrite/baseline guard tests pass; all accepted artifacts unchanged |
| 7. Documentation and closure | README, methodology, completion report, changelog and directly affected source references/viewer excerpts | Commands/limits match actual results; plan closes to report/backlog; unrelated work preserved; diff checks pass |

Proposed affected files (helpers can be local functions; no general framework):

| File(s) | Smallest responsibility change |
|---|---|
| `main.m` | Replace duplicated input/derivation/save block with explicit default-config runner call; retain script entry and filename |
| New `fe2o3_default_config.m`, `fe2o3_validate_config.m`, `fe2o3_run.m` | Defaults; validation/derived data; local-scope orchestration/result/provenance/output helpers around existing scripts |
| `systemgeneration.m`, `One_data_curve.m` | Narrow restart and plotting compatibility adapters described above |
| New `tests/run_cpu_checks.m`, `tests/run_gpu_checks.m`, `tests/check_benchmark_contract.m` | Focused checks; local helpers/fixtures only when needed; report prerequisites and nonzero failures |
| New `benchmark/capture_legacy_reference.m`, `benchmark/run_regression.m`, `benchmark/compare_regression.m`, `benchmark/README.md` | Explicit pre-migration preparation, fresh candidate generation, strict comparator, ownership/transfer instructions |
| `benchmark/baseline_manifest.json` | Add only after user accepts actual data; cases, configuration/input identities, hashes, environment, tolerances, acceptance record |
| `.gitignore` | Ignore only new generated `benchmark/candidates/`, raw `benchmark/baseline/` and test scratch; existing cache patterns already covered |
| `README.md`, `docs/model-methodology.md`, new `docs/cq01-05-verification.md` | API/schema/commands and actual validation evidence; preserve pre-existing README addition |
| New root `CHANGELOG.md` | At implementation, use canonical dotfiles template; one grouped meaningful CQ feature entry |
| `AGENTS.md` | At implementation only: update durable entry/output/test/baseline rules contradicted by the approved interface; keep physics approval rules |
| Existing implementation views/guides | Update only source excerpts/links/line references affected by moving main parameters; validate all existing masks/interactions/boundaries, no redesign |

Expected unchanged numerical files: `constantfile.m`, `distrib.m`,
`exchangej_.m`, `dmi.m`, `dipole_.m`, `field_calc.m`, `integrate_llg.m`,
`atomgpurk4.m`, `atomgpu.m`, `atomgpupc4.m`. Do not edit them to make a failed
check pass without a revised, approved plan. Script workspace variables are
bound locally by the runner with explicit field mapping, not arbitrary caller
state. Boundary: configuration/initialization → retained fields → RHS with
verified field/torque inputs → retained integration → result/save → plotting.
No geometry or material derivation moves into the dynamical RHS.

No new unit-cell or exchange illustration is needed: neither geometry nor
interaction topology changes. Existing `docs/structure-exchange.js` and its
guide may contain moved source references; audit and patch only those affected.
Preserve historical PR/notebook snapshots and `docs/unit-cell.*` schematic labels.

The parent workspace adoption matrix is outside this repository. Leave it
unchanged in this plan-only task; provide evidence for its Fe2O3 entry at
completion. Writing that separate repository requires scope authorization.

## Verification and Acceptance Criteria

### Fresh checks actually performed during planning

- Read source, docs, ignored local verification filenames and selected harnesses;
  read catalog/matrix and ASTRO examples. No production/test implementation made.
- `git status --short --branch`, `git rev-parse HEAD`,
  `git diff -- README.md`, and comparison of old root plan to its archive.
- `git diff --exit-code cf19a99 HEAD -- '*.m'`: exit 1, comment-only DMI
  difference described above, not exact source equality.
- `git check-ignore -v plan.md`: `.gitignore:2` owns the existing rule.
- Planning-document checks: all 13 required template headings occur in order,
  local Markdown links resolve, Typora math delimiters and whitespace checks
  pass; the previous plan body is already preserved in its archive.
  `git diff --check` passed; `git diff --stat` shows only the pre-existing
  three-line README change because `plan.md` is ignored. README and conversation
  file SHA-256 hashes remained unchanged.
- `nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv,noheader`:
  exit 255, OS blocked NVML. This does not establish Windows MATLAB GPU absence.
- MATLAB readiness command (no integration or MAT output):

```bash
matlab -batch "fprintf('MATLAB %s; %s\n',version,computer); fprintf('Parallel license=%d\n',license('test','Distrib_Computing_Toolbox')); disp(ver('parallel')); try, fprintf('Available GPU count=%d\n',gpuDeviceCount('available')); catch ME, fprintf('GPU preflight unavailable: %s: %s\n',ME.identifier,ME.message); end"
```

Returned exit 0: MATLAB `24.2.0.2740171 (R2024b) Update 1`, `PCWIN64`,
Parallel Computing Toolbox 24.2, license test 1, available GPU count 1.
GPU name/driver/device allocation and kernel execution remain unchecked this
turn. Do not relabel historical RTX 3060 results as this turn's measurements.

### Protected baseline and case inventory

There is **no accepted baseline yet**. Missing references are a blocker for
production refactoring, not permission to call a fresh run “accepted”.
After plan approval, preparation uses copies pinned to the inspected revision,
with original source hashes and a strict logged allowlist of main-script
case-control/output substitutions. All physics/kernel files are byte-identical.
Do not regenerate the old algorithm inside a new runner and compare it to itself.
Keep the original workspace outputs plus a minimal extracted comparison record.

Use fresh `benchmark/candidates/<unique-id>/reference-proposal/` paths. Upon
separate data acceptance, store raw artifacts under ignored
`benchmark/baseline/<accepted-id>/` and record hashes/paths in the tracked
manifest. No MAT artifacts enter Git. A clone lacks these raw files: document
manual transfer and integrity verification. Routine run/compare commands never
create, regenerate or promote references and cannot write beneath baseline storage.

All small cases retain `tstep=1e-16 s`, all unspecified model defaults and
exact initial inputs. Configuration changes below are verification fixtures,
not production defaults:

| Required case | Core / updates / stride | Inputs and purpose |
|---|---|---|
| `small-relax` | 5×9×6; 3 chunks × 8; stride 2 | Current alternating state, default zero drives and DMI; sampling/chunk migration |
| `driven-restart` | 5×9×6; 3×8; stride 1 | Fixed nonuniform restart; `jcSOT=6e12`, `jcSTT=1e9 A/m²`, STT_DLT=STT_FLT=1; preserve actual chi-based torques |
| `default-core-short` | 20×20×6; 2×8; stride 2 | Current geometry and physical defaults, shortened duration only; not a 25 ps accepted production benchmark |
| `heun-legacy` | 5×9×6; 3×8; stride 1 | Same nonuniform restart, rk4=0, zero currents; characterize legacy behavior |
| `pc4-legacy` | 5×9×6; 3×8; stride 1 | Same restart, rk4=2, zero currents; crosses bootstrap and chunk boundaries |
| `dipole-cache` | 5×9×6; 2×52; stride 1 | Same restart, RK4, dipolemode=1, zero currents; exercise refresh and chunk-local reset |

Nonuniform fixture, core indices `(r,c,k)`: normalize the vector
`[(-1)^(k+1),0.07*sin(r+2*c+k),0.05*cos(2*r-c+3*k)]` at occupied sites and
set empty sites to zero. Store actual arrays and hashes, not just this recipe.
Default initialization is captured exactly, without replacing its trig residuals.
No RNG seed alone is used as a substitute for stored input arrays.

Generate each case twice in separate fresh destinations before acceptance.
Compare all `mmx/mmy/mmz`, `hdmi_allx/y/z` and every `t` value; require exact
case inventory, shapes, mapping, initial arrays, units, algorithm, sample
policy, effective configuration/coefficients. Source version is expected to
change across migration and is reported, not required identical. Output path,
wall time and unused RNG state are documented nonphysical differences.
Recompute finite/norm/empty invariants from both reference and candidate data.

Provisional same-environment limits: maximum spin-component error `2e-12`,
occupied norm error `2e-12`, DMI field error `2e-10 T`; exact empty zeros,
mask/initial arrays and timestamps (same unchanged construction). These are
based on the repository's focused double-precision PR #20 gates, not certified
full-run/cross-hardware tolerances. Establish them against the two repeats and
independent checks, report measured errors, and obtain acceptance with the
artifacts. Failure does not authorize widening limits. Full-duration or
cross-platform baselines are outside the initial suite.

Comparator tests mutate separate fixtures: intermediate-only spin change,
intermediate DMI change, wrong time value with unchanged shape, missing case,
wrong config/mapping/initial state, NaN/Inf, norm corruption with forged summary,
invalid dimensions/schema and baseline hash mismatch. Each must raise an error
usable by `matlab -batch`; a compatible repeat must pass. Mutating a baseline
fixture must not mutate the actual accepted artifact.

### Focused independent checks and interface gates

1. Structure/topology: current scripts on 5×9×6, 4×8×3 and 6×10×9 cores;
   verify all masks and padding, direct expected neighbor offsets rather than
   a second circshift implementation, reciprocal shell matrices, boundaries
   and DMI antisymmetry. Recheck historical 140-site counts: J1–J4 directed
   entries `114,328,276,544` and DMI `820`. Exact integer comparisons; any
   discrepancy stops for diagnosis. No material coordination claim follows.
2. Static/RHS: isolated exchange and K1/Zeeman algebra; scaled finite
   differences (step `1e-6`, error `<1e-8` in normalized test units), DMI energy
   after antisymmetry; hk6 poles/empty/nonunit stage arithmetic; direct retained
   dipole tensor checks. Compare scalar RK4 RHS with independent vector and
   implicit-Gilbert forms on unit test spins, including isolated/combined DL/FL
   channels. Scale by gamma for numerical error checks, tolerance `2e-12`.
   Check analytical precession and damping directions. No hk6 total-energy claim.
3. GPU RK4: port [5]'s field wrapper in an isolated temporary path, restoring
   path/cwd on exit. Gate stage/spin errors `<2e-12`, field errors `<2e-10 T`;
   each enabled state-dependent exchange/anisotropy/DMI field must vary at a
   chosen nonuniform stage. Reuse the synthetic nonlinear anisotropy test
   `gamma=mus=1`, `alpha=0.3`, `K1=0.75`, z0=0.4, duration 1, 8/16/32/64
   steps: decreasing errors and order between 3.8 and 4.2. It advances an extra
   update solely to expose the pre-update frame at t=1; production unchanged.
4. Time/cache: strides 1/2/8 for 3×8 updates; chunked vs single-chunk identity
   only with dipoles off; fresh DMI evaluation from every stored spin; exact
   sample times. With dipoles on, verify existing refresh/cache behavior across
   52-step chunks without asserting chunk-size invariance or coupled convergence.
5. Config: defaults/derived values match original arithmetic; unknown/missing
   fields, arrays where scalars required, invalid flags/solver, nonfinite data,
   bad dimensions/height, nonintegral timing, bad stride, both distribution flags,
   incompatible/nonunit restart and nonzero legacy torques fail before integration.
   Check no dependency on poison variables in caller workspace and no accidental
   loss of caller variables. Zero-current legacy selections remain available.
6. Save/load: returned and loaded result fields/arrays agree exactly; effective
   configuration matches actual run, correct padded dimensions and DMI/time
   alignment; both restart formats and both plotting-input formats work. Verify
   Git/hash identity including a modified-source fixture, with honest unavailable
   metadata behavior. Source fixture tests run outside accepted/source folders.
7. Output protection: existing output hash unchanged on rejection, no integrator
   call before failure; fresh save succeeds, explicit ordinary replacement works;
   relative/absolute/`..`/alias baseline paths and input/output collision fail;
   a destination created after preflight still cannot be clobbered. Accepted
   artifact hashes match before/after the whole suite.

Commands below are **proposed interfaces, not existing or executed commands**.
Use MATLAB-visible D: paths through the installed Windows MATLAB launcher.
`tempname` gives fresh test directories; no external Python package required.

```bash
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_cpu_checks"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_gpu_checks"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('benchmark'); capture_legacy_reference"
# Pause for acceptance of the actual reference proposal before migration.
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('benchmark'); run_regression"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); check_benchmark_contract"
git diff --stat
git diff --check
```

`run_regression` allocates a fresh candidate directory, runs exactly the
manifest's expected suite through the production API and compares; any missing
reference or failed check raises. `capture_legacy_reference` is a deliberate
one-time preparation command; it creates proposals only. Record exact expanded
paths, commands, code hashes, GPU details and results in the completion report.
No performance, full coupled-field convergence, stochastic, or physical
finite-size validation claim follows from this suite. Failures of legacy
cases must be reported for scope review, never silently repaired or omitted.

## References

[1] [Shared feature catalog](/home/zhuzibn/.agents/skills/spin-dynamics-development/references/code-quality.md),
CQ-01 through CQ-05 and adoption recording; [workspace matrix](../../docs/repository-quality-matrix.md),
2026-11-19 snapshot, Fe2O3 entry. Source-inspection evidence only.

[2] Current repository MATLAB sources at the revision above: `main.m`,
`constantfile.m`, `systemgeneration.m`, `distrib.m`, `exchangej_.m`, `dmi.m`,
`field_calc.m`, `dipole_.m`, `integrate_llg.m`, all three kernels and
`One_data_curve.m`. Implementation evidence, not qualifying material literature.

[3] [Model methodology](model-methodology.md), especially state/geometry,
energy/fields, parameters, integration, known limitations and historical
verification; [structure/exchange guide](structure-exchange-guide.md);
[issue record](../diff.md). No scientific evidence gap is closed by citing these.

[4] Local equation repository
`/mnt/c/Users/zzf-m/OneDrive/presentations/equations`: read `AGENTS.md`,
`README.md`, `spintronic_methods.tex`, `symbols.tex`,
`sections/hamiltonian.tex`, `sections/llg_equation.tex`, and
`sections/atomistic_boundaries.tex`. Locators: `eq:atomistic-hamiltonian`,
`eq:atomistic-effective-field`, `eq:atomistic-pair-fields`, `eq:llg-gilbert`,
`eq:ll-with-flt`, `eq:atomistic-boundary-workflow`. These govern consistency
checks and expose the constant/boundary distinctions above; they do not
provide hematite construction provenance. No equation files modified/compiled.

[5] [PR #20 report](pr20-change-report.html#evidence): embedded harness,
commands/log, reviewed `cf19a99`, stage checks and nonlinear analytical
convergence. [PR #18 report](pr18-time-01-change-report.html#evidence)
and [methodology summary](model-methodology.md#merged-pr-18-verification):
prior pole/time/DMI checks. Historical results only, not fresh execution.

[6] ASTRO files listed under Existing Behavior, rooted at
`/mnt/d/Git_control/spin_dynamics/atomistic_v2`. Read-only software examples.

[7] Local paper knowledge-base `papers.jsonl` and `paper-catalog.md`, searched
for hematite, Fe2O3, Dannegger and exact cited DOIs. Reviewed records found:
`AFM/2023-Magnetic properties of hematite revealed by an ab initio-PRB-Tobias Dannegger`
(title: *Magnetic properties of hematite revealed by an ab initio parameterized
spin model*, DOI `10.1103/PhysRevB.107.184426`) and
`AFM/2020-Electrical Switching of Tristate Antiferromagnetic Néel Order in-PRL-Yang Cheng`
(DOI `10.1103/PhysRevLett.124.027202`). **Metadata/discovery only:** their
summaries/extracted texts/PDFs are indexed under the separately permissioned
OneDrive papers tree, not inside the authorized processed repository. No
matching summary/full-text copy was found inside this KB root. Do not call
these papers absent or use metadata as physical evidence. Source artifacts
were not opened; full journal references/author lists are not reconstructed
from incomplete metadata. No new paper-derived claim is adopted.

[8] Reviewed presentation record
`/mnt/c/Users/zzf-m/OneDrive/presentation-knowledge/notes/ee237/week15a-afm/summary.md`
(read in full) and `slides.md`, slides 25–28, especially slide 28: hematite SOT
switching, citing PRL 129, 017203 (2022). This record supports discovery/context
only, not the code's values, basis or shell mapping. Underlying publication
not used as scientific evidence or retrieved. Catalog/index search had no
matching title; the targeted note-text search located this record.

## Assumptions, Missing References, Risks, and Approval Questions

- **Baseline prerequisite:** none accepted for current source. The six-case
  small suite preserves migration behavior only. Existing local results and
  historical MAT files are protected, not automatically designated baselines.
  Baseline proposal acceptance with measured tolerances is required before
  production refactoring; routine commands cannot bootstrap it.
- **Runtime prerequisite:** MATLAB/PCT/available GPU count confirmed; actual
  GPU allocation/kernels and device identity still need stage-0/3 checks.
  A GPU failure blocks GPU verification; no silent CPU substitution.
- **Reference gaps:** structure/basis/stacking, exact J-shell assignments and
  signs/normalization provenance, a/c/moment/damping/density/thickness,
  anisotropy energy/axes, DMI vectors/coefficients and dipole metric/truncation
  remain REFERENCE MISSING. Search scope: [7] KB metadata and local topic
  notes; [8] presentation catalog/index and targeted text/summary; existing
  methodology. Relevant source artifacts require separate access. No network
  or C2DB search was made because this task adds no physical construction or
  parameter claim. Gaps block new model assertions, not the retained-code
  software migration; defer source reconciliation rather than guessing.
- **Numerical gaps:** full coupled-field convergence, physical size/edge
  convergence, complete legacy solver validation, nonzero D1 and full 25 ps
  reference remain unverified. Keep these visible in result metadata/reports.
- **Compatibility:** approve the new single-struct output plus narrow
  restart/plot adapters; arbitrary leftover-workspace variables cease to be
  public output. Default filename, trajectory samples and physical values
  remain the same. Existing files now cause an early error by default.
- **Validation changes:** early rejection of malformed configurations,
  incomplete layer motifs, nonintegral time controls and invalid restart
  spins is intentional. No silently “fixed” physical input or solver selection.
- **Scope control:** if a retained kernel fails a focused check, document the
  counterexample and seek a separate approved numerical/physics change; do not
  widen this task. Preserve legacy limitations while testing their migration.
- **Documentation:** stages 0–4 updated README/methodology, benchmark commands,
  durable reference-handling rules, the grouped changelog and verification
  report. Stages 5–7 must update the API/output documentation and affected
  viewer source references. Parent
  matrix changes remain separately scoped. No global instruction changes.

**Approved decisions:** the user's “proceed” approves revision 1's staged scope, six small reference
cases, local ignored baseline storage with tracked manifest, callable API and
single-struct output/consumer adapters, and default non-overwrite policy.
Implementation starts with harness/reference preparation; approval does not
accept future reference artifacts. The concrete baseline proposal is the next
review boundary before production-interface edits.

**Concrete review now pending:** accept the repeat-1 suite linked at the top
as `benchmark/baseline/cq-20261120-v1/`, together with `2e-12` spin/norm and
`2e-10 T` DMI limits and exact times/inputs. Both repeats and all 22 corruption/
path checks passed; all 13 production sources remain unchanged. Approval
would authorize copying the accepted artifacts, writing the tracked acceptance
manifest and continuing the already scoped stages 5–7. Nothing has been
promoted, committed or pushed.
