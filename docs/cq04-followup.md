# CQ-04 acceptance-manifest and relative-path repair — 2026-11-20

The user approved revision 1 of the root plan with **“proceed”**. This follow-up
protects `benchmark/baseline_manifest.json`, including direct, real-parent,
case and hard-link aliases, regardless of overwrite. Production and candidate
path helpers also anchor relative paths to MATLAB `pwd` before Java resolution.
Fresh saves and explicit ordinary-result replacement retain their contracts.
No commit, push, branch change, accepted-data change or dependency installation.

## Preserved model and boundaries

Unchanged Class 4 intended hematite model, fixed-length LLG/LLGS and unit occupied
spins. The full parameters, Hamiltonian/field/RHS definitions, six-mask geometry,
25 ps defaults, sampling and restart policy remain in the
[approved migration specification](cq01-05-plan.md) and
[methodology](model-methodology.md). No numerical/constant source changed.
The only production source edit is `fe2o3_output_path.m`; the benchmark candidate
path helper and focused tests implement the same software ownership contract.

Crystal/J-shell/DMI mapping, hk6 energy and material-source limitations remain
**REFERENCE MISSING** or unresolved. The plan consulted the reviewed EE237
week15a summary (`presentation-knowledge/notes/ee237/week15a-afm/summary.md`,
slides 25–28) for context only; no material value was adopted. The shared
`spintronic_methods.tex`, `symbols.tex` and `sections/llg_equation.tex`
(`eq:llg-gilbert`, `eq:ll-with-flt`) informed unchanged conventions, not material
evidence. No paper PDF, equation source or visualization was changed.

## Reproduction and implementation gates

On the preserved source copy, requesting its disposable sibling manifest with
`overwrite=true` returned an allowed destination. The sentinel bytes were checked;
no save was attempted. The replacement guard rejects the same ownership class.
`fe2o3_save_result` already invokes this guard before serialization and again
before exclusive publication; those mechanisms remain unchanged.

The first complete check exposed a second path defect: the relative benchmark
destination was resolved against Java's startup directory instead of MATLAB's
changed directory. The exact existing `protected-path-5` test failed. Anchoring
relative paths to `pwd` in both helpers repaired this contract; no tolerance or
expected rejection was weakened. All final checks below passed.

## Commands and results

Snapshot: `main`, HEAD `b90daac5d5407462db49850284c3c1f58b21ba29`, plus the
recorded pre-existing README changes and this patch. Runtime: Windows MATLAB
`24.2.0.2740171 (R2024b) Update 1`, JVM, existing PCT/NVIDIA GPU.

Final command, from WSL:

```bash
matlab -wait -logfile 'D:/Git_control/spin_dynamics/.sisyphus/cq-followup-20261120-v1/fe-final.log' -batch "cd('D:/Git_control/spin_dynamics/Fe2O3'); addpath('tests'); run_interface_checks; check_benchmark_contract;"
```

| Check | Final evidence |
|---|---|
| Disposable acceptance manifest | 12 rejections: direct, MATLAB-relative, `..`, hard link, case alias and actual Windows junction, each with overwrite false/true; sentinel bytes preserved |
| Actual tiny GPU interface | Defaults/types/invalid inputs, in-memory/saved equality, dirty-source hashes, late writer, explicit ordinary replacement, restart corruption/masks, plotting and shadow rejection pass |
| Benchmark contract | All 23 deliberate corruption/path cases reject, including the previously failing relative path and actual junction |
| Existing accepted suite | All six stored cases validate with zero comparison differences; this is reference/copy validation, not six new simulation trajectories |
| Source preservation | Nineteen of twenty captured root MATLAB files remain byte-identical; only the production path helper differs |

Final interface artifacts are under
`benchmark/candidates/tp397fdca9_5bb8_49c8_980b_b9e3cb9f9282/`; final mutation
fixtures are under `benchmark/candidates/tp2ed86dcf_f40b_4d5d_9ca4_3cc35f7f3bdc/`.
The before source/sentinel and logs are in the workspace's ignored
`.sisyphus/cq-followup-20261120-v1/`. Early launcher/harness setup failures and
the initial relative-path failure are not counted as passes. The independent
reviewer inspected code/ownership; runtime evidence comes from the commands above.

No full-duration simulation was needed for this path-only repair. Short numerical
agreement and provenance do not resolve the model's physical evidence gaps.
