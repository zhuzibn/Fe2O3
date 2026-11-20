# CQ-01–CQ-05 verification and completion

2026-11-20. **All seven approved stages completed.** Original source revision:
`523500cf85d32462361776499937eb4ce25e7fd7` (`main`). The user's “proceed” approved
plan revision 1; the subsequent “continue” accepted the concrete repeat-1 suite
and tolerances and authorized the configuration/output migration. The
[archived plan](cq01-05-plan.md) preserves the model specification, exact source
search boundaries, claim/source map and approval stages.

## Outcome and preserved scope

| Feature | Verified outcome |
|---|---|
| CQ-01 | Accepted six-case baseline, pinned artifact hashes, exact migration agreement, strict raw-sample/invariant comparator and missing-reference failures |
| CQ-02 | One defaults function, complete validation, callable local-scope runner and thin script entry |
| CQ-03 | Versioned `fe2o3_result`, actual inputs/geometry/initial state/trajectories, source/runtime/RNG provenance, old/new restart and plotting compatibility |
| CQ-04 | Existing-output refusal by default, protected references/restart inputs including aliases, exclusive publication and explicit atomic replacement |
| CQ-05 | Maintained CPU/GPU, interface and benchmark-contract checks with explicit evidence limits |

The ten original numerical/constant sources remain byte-identical to
[legacy_source_manifest.json](../benchmark/legacy_source_manifest.json).
Only `main.m`, `systemgeneration.m` and `One_data_curve.m` were adapted among
the original 13 files. This remains a Class 4 fixed-length LLG/LLGS code;
physics, defaults, boundary rules, legacy torque gating, solver formulas,
sampling and the unsaved final endpoint are unchanged. All six migration
cases match the accepted reference with **zero spin and DMI differences**;
all effective inputs, masks, initial arrays and timestamps match exactly.

Existing README local-checkout text and the untracked DMI conversation are
preserved. The conversation SHA-256 remains
`3a937c1cbd36708a2d4fae05b6d2d3aca8617d7a26b060f3736212613e7f670d`.
No commit, push, branch/worktree change, dependency installation, network use
or original-result overwrite occurred. The parent adoption matrix remains
unchanged; this report supplies evidence for a separately authorized update.

## Accepted reference and preparation evidence

Accepted raw directory: `benchmark/baseline/cq-20261120-v1/`.
[baseline_manifest.json](../benchmark/baseline_manifest.json) records acceptance,
identity, cases, environment and tolerances. Repeat-1 was copied byte-for-byte;
its internal historical proposal manifest was not rewritten. Every supporting
artifact was rehashed after copying and during comparisons.

The original two-repeat preparation remains under:

```text
D:\Git_control\spin_dynamics\Fe2O3\benchmark\candidates\tpefe2942b_e5df_4b04_9815_b3d4584c3544\reference-proposal
```

[Historical preparation status and measurements](../benchmark/candidates/tpefe2942b_e5df_4b04_9815_b3d4584c3544/reference-proposal/status.json).
The following preparation results predate production migration; the later
migration/interface results are recorded separately below.

| Suite | Manifest SHA-256 |
|---|---|
| [repeat-1](../benchmark/candidates/tpefe2942b_e5df_4b04_9815_b3d4584c3544/reference-proposal/repeat-1/manifest.json) | `94e2b7a2ed194a0f9b2d2e5ec62774a890fa46665f0a872a5e731aec7fb13bb2` |
| [repeat-2](../benchmark/candidates/tpefe2942b_e5df_4b04_9815_b3d4584c3544/reference-proposal/repeat-2/manifest.json) | `3514e8a627713dd453d0ca587ffe9ac9aa0834f34bd008317993079c40a4e185` |

Every case includes source copies, logged substitutions, fixed initial arrays,
original workspace output and a host-array comparison record. Each suite has
119 unique hashed artifacts (184 manifest references including reused source
entries); complete suite sizes are 3,569,436 and 3,569,421 bytes respectively.
The raw files are ignored and require separate transfer to another computer.

All cases use `tstep=1e-16 s`. Full overrides and preparation recipe are in
[benchmark/README.md](../benchmark/README.md).

| Case | Core | Updates | Saved frames | Max spin difference | Max DMI difference (T) | Max occupied norm error, either repeat |
|---|---|---:|---:|---:|---:|---:|
| small-relax | 5×9×6 | 24 | 12 | 0 | 0 | 1.11022302463e-16 |
| driven-restart | 5×9×6 | 24 | 24 | 0 | 0 | 2.22044604925e-16 |
| default-core-short | 20×20×6 | 16 | 8 | 0 | 0 | 1.11022302463e-16 |
| heun-legacy | 5×9×6 | 24 | 24 | 0 | 0 | 2.22044604925e-16 |
| pc4-legacy | 5×9×6 | 24 | 24 | 0 | 0 | 2.22044604925e-16 |
| dipole-cache | 5×9×6 | 104 | 104 | 0 | 0 | 2.22044604925e-16 |

All timestamps, configurations, mappings and initial arrays match exactly.
Reference and candidate raw arrays are finite and empty entries are exactly
zero. Proposed acceptance tolerances remain `2e-12` for spin components/norms,
`2e-10 T` for DMI, and exact times/inputs/masks. These are same-environment
short-case limits based on focused repository checks; they do not certify
cross-hardware agreement or a 25 ps production trajectory.

## Environment and commands actually executed

MATLAB `24.2.0.2740171 (R2024b) Update 1`, PCWIN64;
Parallel Computing Toolbox 24.2; NVIDIA GeForce RTX 3060 Laptop GPU,
compute capability 8.6, driver 610.60. Additional device/toolkit and RNG
metadata are embedded in each record and suite manifest. Planning-time
`nvidia-smi` failed because OS access to NVML was blocked; Windows MATLAB
successfully allocated the GPU and executed these checks.

Successful commands from WSL, via the installed MATLAB launcher:

```bash
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_cpu_checks"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_gpu_checks"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('benchmark'); capture_legacy_reference"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); check_benchmark_contract('D:/Git_control/spin_dynamics/Fe2O3/benchmark/candidates/tpefe2942b_e5df_4b04_9815_b3d4584c3544/reference-proposal')"
git diff --stat
git diff --check
```

All four final MATLAB commands exited 0. Initial preparation attempts exposed
two new harness bookkeeping errors when growing untyped empty structure
arrays; explicitly typed empty records fixed them. Those incomplete proposals
remain ignored and unaccepted. One launcher attempt failed with a WSL vsock
startup error before MATLAB; a direct launcher retry succeeded. These were
not numerical model failures. No source kernel was edited to obtain a pass.

## Focused evidence

CPU topology checks execute actual scripts against direct-index tables from
the implementation guide, without reimplementing their circshift arithmetic:

| Core | Occupied sites | J1/J2/J3/J4 directed entries | DMI entries |
|---|---:|---|---:|
| 5×9×6 | 140 | 114 / 328 / 276 / 544 | 820 |
| 4×8×3 | 48 | 32 / 106 / 72 / 140 | 212 |
| 6×10×9 | 270 | 240 / 654 / 584 / 1160 | 1744 |

All six masks, every basis probe/component, reciprocal shell matrix,
finite edge/surface omission and DMI antisymmetry check passed. Isolated
exchange, K1, external-field and D3/D4/off cases passed analytical and
Cartesian finite-difference tests (step `1e-6`, maximum error below `1e-8`
in normalized test units). hk6 passed arithmetic checks for poles, empty,
off-pole and nonunit stage inputs; no hk6 energy was invented. The retained
rounded-vector dipole tensor passed direct indexing checks. Unit-spin RHS
tests matched an independent implicit-Gilbert linear solve, including each
isolated DL/FL channel, combined torques, precession and damping signs.

GPU RK4 independent reconstruction:

| Quantity | Measured maximum error |
|---|---:|
| Stage spin | 2.77555756156e-17 |
| Fresh field | 1.7763568394e-14 T |
| Completed one-step spin | 2.22044604925e-16 |

The nonlinear analytic anisotropy case used synthetic `gamma=mus=1`,
`alpha=0.3`, `K1=0.75`, initial z=0.4, duration 1. Errors for 8/16/32/64
updates to the sampled endpoint were `1.80185888987e-7`,
`1.16962723951e-8`, `7.44221573257e-10`, `4.69210226228e-11`;
observed orders `3.94536531, 3.97417279, 3.98742585` passed the 3.8–4.2 gate.
An extra update exposes the pre-update frame at time 1; production sampling
has not been changed.

The pre-migration harness executed the actual `main.m` timing block. The
maintained migrated harness now uses the canonical validator clock. Strides 1/2/8
match the dense spin/DMI trajectory and exact times, and chunked/single-chunk
results match with dipoles disabled. Every sampled DMI field matches a fresh
evaluation. Dipoles refresh at global pre-update indices 0,50,52,102 for
two 52-step chunks and remain cached otherwise; fresh stage-2 dipoles differ
by `0.000264275883243 T`, confirming the retained approximation. Heun/PC4
finite/norm/empty and chunk-continuity checks pass, including PC4 bootstrap.
Raw invariant checks also cover the analytical and dipole runs.

Final GPU log:
`.cq-checks/tp36871d1a_8f32_4db2_b483_8fdf8a853122/checks.log`.
The earlier GPU pass was rerun after independent review strengthened the
clock/sparse-sampling tests; only the final run establishes those extra gates.

## Pre-migration comparator, protection and independent audit

All **22** rejection checks passed with the required `Fe2O3:*` identifier:
intermediate-only spin/DMI corruption, wrong time/config/mapping/initial state,
NaN, Inf, invalid raw norms with a forged summary, wrong shape/schema, missing
case, corrupt baseline hash/nonfinite baseline, missing accepted manifest,
existing directory/file, absolute/relative/`..`/Windows case baseline paths,
and a directory-junction alias. The junction points to the benchmark parent;
no accepted baseline directory or sentinel was created. Original suite hashes
and the full repeat comparison passed again afterward.

Mutation scratch:
`benchmark/candidates/tpe49b34a2_1463_489e_bce8_86019707d738/`.
The test-owned junction and corrupted fixtures remain in that ignored directory.
They are not scientific outputs or references.

A separate read-only verifier reviewed the checks and comparator, independently
reconstructed neighbor counts, hashed all 238 unique proposal artifacts and
all 13 production sources, checked both manifest hashes, and reconstructed
all 12 executed `main.m` copies exactly from source plus logged substitutions.
No unlogged modifications were found. The verifier did not independently
rerun MATLAB; numerical evidence comes from the actual commands above.

## Evidence boundaries

The local equation reference consulted for this task is
`/mnt/c/Users/zzf-m/OneDrive/presentations/equations/spintronic_methods.tex`,
with `symbols.tex`, `sections/hamiltonian.tex`, `sections/llg_equation.tex`
(`eq:llg-gilbert`, `eq:ll-with-flt`) and `sections/atomistic_boundaries.tex`.
It supplies notation/convention checks, not material parameter provenance.
No equation files were modified or compiled. The paper KB metadata records
for Dannegger 2023 (`10.1103/PhysRevB.107.184426`) and Cheng 2020
(`10.1103/PhysRevLett.124.027202`) were discovery evidence only; their separate
OneDrive source artifacts were not accessed. The reviewed EE237 Week 15a AFM
summary and slide 28 were also consulted for context, not code parameters.
Exact paths/search limitations remain in the [archived plan](cq01-05-plan.md).

Crystallographic construction/site mapping, exchange-value provenance,
physical DMI vectors and shell assignments, material normalization, complete
sixth-order energy, dipole truncation and several parameter locators remain
**REFERENCE MISSING / UNRESOLVED**. No full 25 ps run, full coupled-field or
physical size convergence, stochastic ensemble, or complete legacy-solver
accuracy validation is claimed. The implementation viewer keeps identical
masks/neighbor/DMI behavior; its moved configuration excerpts and links were
updated to the new defaults/validator locations.

## Migration, interface and persistence checks

Commands executed from WSL through the existing Windows MATLAB launcher:

```bash
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('benchmark'); run_regression;"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_cpu_checks; run_gpu_checks; run_interface_checks;"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_interface_checks; check_benchmark_contract;"
matlab -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_interface_checks;"
```

All exited 0. The last interface run covers review corrections (empty MATLAB
string rejection, clean-or-dirty Git status comparison and Windows-only path
case folding) plus a deliberately modified-source fixture. No numerical kernel
changed to obtain a pass.

| Evidence | Artifact / outcome |
|---|---|
| Six-case migrated suite | [candidate manifest](../benchmark/candidates/tp5b9c0742_c3e7_4080_94fb_a8d27b183562/manifest.json); spin error **0**, DMI error **0 T** for every case; exact configuration/mapping/initial/time equality |
| Migrated CPU checks | Three topology fixtures and all term/RHS checks pass with the preparation measurements above |
| Migrated GPU checks | [checks.log](../.cq-checks/tpf5e7ebc2_3c9a_45f1_878b_6ad6667b90db/checks.log); unchanged stage errors and fourth-order analytic convergence measurements above; sparse clocks/cache/legacy characterization pass |
| Final interface checks | [interface.log](../benchmark/candidates/tpd25c1858_354b_4c4a_80f6_2e23fbcc578a/interface.log); canonical defaults and numeric classes; invalid inputs; memory/save equivalence; old/new restart and plot consumers; source identity checks pass |
| Comparator/path mutations | Fixtures in `benchmark/candidates/tpea696300_2b5b_4987_8c7b_a05f50f6ad40/`; **23 expected rejections**, including independently missing acceptance and raw suites; original accepted artifacts reverified |

The interface suite checks all restart history frames (including corruption in
an earlier frame), dimensions, finite real double values, occupied norms and
empty zeros. It rejects mismatches in each of the six site masks and occupancy,
unsupported schema and output/restart collisions. New and legacy inputs select
exactly the last saved frame, without renormalization.

Output tests cover existing files/directories, missing parents, invalid overwrite
values, absolute/relative/`..`/case/junction aliases into reference storage,
a writer appearing after preflight, exclusive hard-link collision behavior,
MAT round trip with only `fe2o3_result`, and explicit atomic replacement on NTFS.
No unsafe overwrite fallback is implemented. The in-memory path creates no output.

Provenance checks compare recorded hashes with actual source bytes (including
untracked helpers), compare Git status with the actual clean-or-dirty status,
reject a shadowed field script, and execute a test-owned source copy with an
extra comment. Its changed hash is recorded while its numerical trajectory
remains identical. Input state, effective coefficients, actual restart hash,
runtime/GPU identity and RNG states are explicit fields. Source/Git identity is
checked before and after integration; it is not inferred from commit ID alone.

Independent read-only review checked all 184 accepted manifest entries and all
12 migrated result/record artifacts. It found no physics/default, full-history
restart, baseline-pinning or NTFS publication defect. The three interface issues
it identified were corrected and the final interface suite passed.

## Visualization and documentation closure

Updated only the moved DMI parameter excerpt/links and related descriptions.
Embedded `exchangej_.m`/`dmi.m` snapshots match current source bytes as text;
the DMI field excerpt matches `field_calc.m`. Browser verification exercised
**240 selections** (six masks × five interactions × two actual layers × two
lateral positions × two padding displays) and **1,056 neighbor highlights**.
Targets/neighbors, row counts and nonempty source excerpts passed throughout.
Illustrative crystal geometry and historical notebook/PR reports are unchanged.

Documentation: updated root README, benchmark README, methodology, this report,
implementation guide/viewer, project AGENTS, grouped CHANGELOG entry and plan.
The detailed reviewed plan is archived; root `plan.md` now links the evidence
and lists only separately scoped scientific/solver/output backlog items.

Final file review checked all local documentation links and rehashed all 119
distinct accepted artifacts. The ten numerical/constant sources and unrelated
conversation remain unchanged. `git diff --stat` shows only scoped tracked
edits; new helper/test/report files were inspected separately. Plain
`git diff --check` flags the preserved CRLF endings on ten added plotting-adapter
lines; `git -c core.whitespace=cr-at-eol diff --check` passes. The rest of that
plotting file is byte-identical, with its original CRLF endings retained.
