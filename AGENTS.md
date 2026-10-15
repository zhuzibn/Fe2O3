# Project Agent Instructions

## Scope

These instructions apply to this repository. The project is a MATLAB/GPU
research implementation of a material-specific atomistic spin-dynamics model
for $\alpha$-Fe$_2$O$_3$.

Use the installed `spin-dynamics-development` skill for any physics, solver,
parameter, topology, boundary, torque, or validation work. Follow the global
instructions in addition to this file.

## Read before changing the model

Before proposing or making a relevant change, read:

1. `README.md`;
2. `docs/model-methodology.md`;
3. the MATLAB files directly involved;
4. the applicable sections of
   `/mnt/c/Users/zzf-m/OneDrive/presentations/equations/spintronic_methods.tex`;
5. qualifying reviewed-paper or presentation records needed to support the
   physical claim.

`Alpha-Fe2O3.pdf` is an annotated design notebook. Use it to interpret the
existing diagrams and implementation history, but do not treat it as sufficient
scientific evidence for a crystal structure, interaction topology, material
parameter, or benchmark.

## Model contract

- Classification: Class 4 material-specific crystallographic atomistic model.
- State: one normalized classical spin per occupied site.
- Equation family: fixed-length LLG/LLGS.
- Field convention: effective flux densities and field-equivalent torque
  amplitudes are in tesla; `gamma` is in $\mathrm{rad\,s^{-1}\,T^{-1}}$.
- Exchange convention: the implemented field is
  $\mathbf{B}_{\mathrm{ex},i}=-(1/\mu_s)\sum_jJ_{ij}\mathbf{s}_j$; positive
  $J$ favors antiparallel alignment under the corresponding energy convention.
- Anisotropy convention: the implemented second-order energy is $-K_1s_z^2$;
  the additional `hk6` field is preserved from `code_20261013`. Its complete
  energy/field mapping remains unresolved; use `docs/model-methodology.md`
  for the implemented formula and known pole defect.
- Conservative fields, nonconservative torques, integration, and simulation
  orchestration must remain distinguishable even though the current scripts
  are not yet architecturally separated.

Do not change a sign, unit, prefactor, pair-counting rule, layer/index
orientation, site mask, neighbor offset, boundary rule, stochastic convention,
or default parameter without an explicit derivation, qualifying evidence, and
user approval of the resulting behavior change.

## Sources of truth and unresolved items

Use this precedence:

1. Checked-in code and accepted baselines define current behavior.
2. `docs/model-methodology.md` records its intended interpretation and known
   validation gaps.
3. Original peer-reviewed sources define physical claims when cited at an exact
   locator.
4. `Alpha-Fe2O3.pdf` records historical reasoning and unresolved annotations.

Do not silently reconcile conflicts. The exchange-value provenance, complete
crystallographic mapping, DMI directed bonds, dipole truncation justification,
and several parameter locators are currently unresolved. Preserve the visible
`UNRESOLVED` or `REFERENCE MISSING` status until qualifying evidence and a
specific validation result justify removing it.

Repository code, the notebook, and `spintronic_methods.tex` may document
behavior or conventions, but they do not replace a qualifying source for a
material-specific structure or parameter.

## Planning requirement

Before any spin-dynamics code modification, create or update root `plan.md`
using the `spin-dynamics-development` skill. Include every required heading,
the compact model specification, claim-to-source mapping, staged verification
gates, and exact unresolved inputs. Ask the user to approve that revision
before editing implementation files.

Documentation-only corrections that do not alter or newly assert the model do
not require `plan.md`. A documentation change that selects a new physical
interpretation, parameter provenance, or supported behavior does require the
same evidence and approval as a code change.

## Implementation boundaries

- Keep changes surgical and preserve the existing script-based architecture
  unless a refactor is explicitly requested.
- Treat `distrib.m` and `exchangej_.m` together when changing site masks or
  neighbor topology. Verify every layer/site type rather than extrapolating
  from one branch.
- Remember that `dmi.m` documents a reversal between MATLAB layer-index
  direction and physical $z$.
- Do not call the padded `circshift` construction periodic without a verified
  neighbor audit. The occupied core currently sees empty padding laterally.
- No thermal field is currently implemented; the imported code has no
  `thermalenable` control.
- Only `atomgpurk4.m` accepts the SOT/STT-labeled torque arguments. Do not claim
  torque support for the Heun or predictor-corrector kernels.
- Do not claim that the predictor-corrector selection is runnable until its RK4
  bootstrap signature mismatch is fixed and the complete path is tested.
- The current kernels hold effective fields fixed within a nominal integration
  step. Preserve that behavior unless the approved task explicitly changes and
  validates the integrator.
- Do not add unrequested physics, dependencies, configurability, or output
  formats.

## Run and output safety

The default run is GPU-dependent and potentially expensive. Do not launch the
full default simulation merely as a syntax check. Prefer small deterministic
verification cases defined in an approved plan.

`main.m` writes `final.mat` without prompting. Tests and candidate runs must
use a separate, non-destructive output path and must not overwrite a user's
existing result or accepted baseline. Do not add generated `.mat` files to Git
unless the user explicitly approves a documented baseline policy.

Do not install MATLAB toolboxes, GPU software, or other dependencies without
explicit permission. If MATLAB or GPU support is unavailable, report which
checks remain unverified.

## Verification gates

For relevant changes, run the smallest applicable checks in this order:

1. Structure and mask counts on a small lattice.
2. Explicit neighbor lists, reciprocity, shell multiplicities, and boundary
   coordination for all six site masks.
3. One-term-at-a-time analytical or finite-difference field checks.
4. LLG/LLGS precession, damping, torque-sign, unit-norm, and finite-value checks.
5. Solver timestep convergence and, when fields are state dependent, substage
   field-evaluation checks.
6. Boundary- and size-convergence checks for spatial claims.
7. Comparison with a qualifying analytical, published, or accepted repository
   reference.

Do not use trajectory plausibility or agreement between two implementations as
a substitute for analytical and topology checks. Record exact commands,
parameters, hardware-relevant details, and outcomes. Do not invent results.

Before finalizing a file-changing task, run `git diff --stat` and
`git diff --check` when practical, inspect unexpected large diffs, and perform
the documentation-impact check required by the global instructions.

## Documentation maintenance

- Update `README.md` when requirements, entry points, configuration, outputs,
  or repository layout change.
- Update `docs/model-methodology.md` when implemented equations, units, signs,
  topology, parameters, boundaries, initial conditions, solver behavior,
  evidence status, or validation results change.
- Update this `AGENTS.md` only for durable workflow, ownership, source-of-truth,
  or verification rules.
- Use `$...$` and `$$...$$` math delimiters for Typora-compatible Markdown.
- Keep implemented behavior, published methodology, interpretation, and future
  proposals visibly separate.
