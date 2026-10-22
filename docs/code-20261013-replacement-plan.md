# Import code_20261013 with the approved DMI state correction

Approved revision: 2026-10-15, revision 1, with the user's explicit instruction
to remove the old plan rather than archive it. All ten difference-review
items were accepted in `diff.md`. This is the completed import plan retained
as a historical approval record; relative file paths below are relative to
the repository root.

**Status after merged PR #18 (`864d6b0`):** the import-time deferrals below
are historical. ANI-01 and OUT-01 are corrected; TIME-01 update counts and
saved-time labels are corrected, while saving the final endpoint remains
open. The later J3/J4 DMI operator passes the recorded antisymmetry checks.
See [current issue status](../diff.md#current-issue-status-2026-10-22) and
[current methodology](model-methodology.md). These later fixes do not alter
the original import decisions recorded below.

## Goal and Success Criteria

Replace the MATLAB implementation in the active Fe2O3 repository with the
user-supplied `code_20261013`, with one approved correction: in `dmi.m`, use
the current spin arrays for the added unshifted adjacent-layer terms.

Source: `C:\Users\zzf-m\OneDrive\projects\2026-Fe2O3 switching-xuezhong\code_20261013`.
Destination: `C:\Users\zzf-m\OneDrive\code_softwares\spin_dynamics\Fe2O3`.
Reviewed local commit: `d0ab252d55ef70328be2f45ea1cba51201288522` on `main`.

Success means all 12 source MATLAB files and `Readme.txt` are present;
their text matches the source except the approved DMI substitutions;
existing destination line endings are preserved; meaningful checks below
pass; and current documentation accurately describes the imported code and
deferred defects. Do not commit, push, modify the source folder, overwrite
simulation results, or run the full default simulation.

This is an explicitly requested transfer of an existing research implementation,
not construction or independent physical validation of a new hematite model.
The accepted parameters and formulas remain provisional where evidence is
missing. No additional physics, solver repair, optimization, or output-format
change is included.

## Model Classification and Equation Family

Class 4 material-specific atomistic model; fixed-length LLG/LLGS. Each occupied
array entry carries one classical Fe spin $\mathbf{s}_i$, normalized after
each solver update. Six color masks are implementation site labels, not a
verified crystallographic sublattice assignment. The supplied initial state
alternates spin direction between layers. All sites use one moment magnitude.
Electronic bands and electronic-structure parameter generation are outside
this implementation [1, distrib.m; systemgeneration.m; atomgpurk4.m].

Inputs are script variables, fields and torque amplitudes in tesla, and
time in seconds. Outputs are saved spin/field histories and the remaining
MATLAB workspace. Intended numerical invariants are finite occupied spins
and unit norm; the accepted anisotropy pole defect prevents a general
finite-value guarantee.

## Existing Behavior and Reference Implementations

The local version uses a 100 x 100 x 6 core, domain-wall initialization,
enabled dipole coupling, zero-current torques, axial second/fourth-order
anisotropy, nominal 200 ps duration, and `finalnew.mat` [2]. The accepted
replacement uses the settings below [1; 3].

Read: repository `AGENTS.md`, `README.md`, `docs/model-methodology.md`,
directly affected MATLAB sources, the source readme and plotting helper,
and the relevant equations files [4]. The ASTRO `AGENTS.md` and `README.md`
were consulted for reference context [6]; its production implementation,
parameters, baseline, and tests are not imported. There is no cross-model
comparison in this task.

The old root `plan.md` is a different audit, including descriptions of files
and behavior no longer present. The user explicitly approved removing that
audit and replacing it with this plan, without an archive. Existing interactive documentation embeds the pre-import DMI source
and signs; label it as historical rather than silently treating it as the
new implementation.

## State Variables, Geometry, and Boundaries

Core dimensions are 20 x 20 x 6. The identical occupancy builder adds two
empty rows/columns on each lateral side, giving 24 x 24 x 6 arrays. Occupied
entries have `atomtype_=1`; empty entries are masked. The three-layer motif
and all six site masks are unchanged. Lateral finite-range interactions see
empty padding. Layer-end branches omit unavailable neighboring layers for
the accepted six-layer geometry; do not claim arbitrary-height support.

| Physical detail | Adopted value or rule | Builder or interaction mapping | Qualifying reference and exact locator | Evidence status |
|---|---|---|---|---|
| Material/site mapping | Existing three-layer motif and six color masks | `distrib.m` | None established for the exact mapping | REFERENCE MISSING; preserved code, no reconstruction |
| Basis, lattice vectors, stacking | No complete crystallographic coordinates supplied | Array/mask representation only | None established | REFERENCE MISSING; no new geometric interpretation |
| Lattice constants | `a=5.067e-10 m`, `c=13.882e-10 m` | Density/thickness definitions in `main.m` | Catalog record [5] identifies a relevant paper but does not establish these exact code assignments | REFERENCE MISSING |
| Magnetic volume and thickness | `ms=24*mus/(sqrt(3)*a^2*c)`, `tz=c*natomH/6` | Torque coefficients in `main.m` | None established for this normalization | REFERENCE MISSING; user-selected definitions |
| Exchange topology | Existing J1-J4 masks, offsets, boundary branches; bulk counts described locally as 1,3,3,6 | `exchangej_.m`, with verified variable renames | No exact published shell mapping | REFERENCE MISSING; code equivalence is established, physical shell assignment is not |
| Exchange values and sign | Values in parameter table; positive J in the implemented field favors antiparallel alignment under the stated energy | `main.m`, `integrate_llg.m` | General convention [4]; material-specific provenance remains unresolved | REFERENCE MISSING for material assignment; unchanged coefficients |
| DMI directed bonds | Source signed neighbor sums plus added adjacent-layer terms, using current state | `dmi.m` | No independently verified material-specific bond table | REFERENCE MISSING; requested source behavior, with a software indexing correction |
| Boundaries and finite geometry | Existing padding and layer guards | Builder and interaction scripts | No physical surface-termination validation | REFERENCE MISSING; preserve existing rule |

Array dimensions, variable names, padding implementation, and storage layout
are software descriptions, not claims that crystallography has been verified.
Missing evidence blocks any new reconstruction or assertion of validated
physical geometry, not the user's explicit transfer of the reviewed code.

## Hamiltonian

The convention consistent with the retained exchange, second-order anisotropy,
and external field is [2; 4, eq:atomistic-effective-field]

$$
\mathcal H_{\mathrm{ex}+2+Z}
=\frac12\sum_{i\ne j}J_{ij}\mathbf s_i\cdot\mathbf s_j
-K_1\sum_i s_{iz}^2-\mu_s\sum_i\mathbf B_{\mathrm{ext}}\cdot\mathbf s_i.
$$

The half factor assumes both directed occurrences of each reciprocal pair.
Parameters are energies per site or link as used in the code, not energy
densities. The old $-K_2s_z^4$ contribution is removed.

A verified Hamiltonian/derivative mapping for the imported sixth-order
anisotropy and complete DMI bond table is UNRESOLVED. Preserve the source
fields without inventing a validated energy. Under the general convention
$\mathcal H_D=-\frac12\sum_{i\ne j}\mathbf D_{ij}\cdot
(\mathbf s_i\times\mathbf s_j)$, antisymmetric directed coefficients would
give $\mathbf B_{D,i}=\mu_s^{-1}\sum_j\mathbf s_j\times\mathbf D_{ij}$
[4, eq:atomistic-pair-fields]. Whether the imported physical bond assignment
is justified remains a separate scientific question. The optional truncated
dipole term is preserved but disabled. SOT and STT are nonconservative and
are not included in this Hamiltonian.

## Parameters and Evidence

All adopted numerical values below are transcribed from source `main.m` and
`constantfile.m` [1]. They are accepted code settings; exact material-specific
publication locators remain REFERENCE MISSING unless explicitly stated.

| Parameter | Adopted value/unit | Use and status |
|---|---|---|
| `J1,J2,J3,J4` | `1.3740e-21,-1.1696e-21,4.0413e-21,2.8041e-21 J` | Unchanged exchange-link coefficients; existing provenance conflict retained |
| `Ksim1` | `-8e-22 J` | On-site second-order term |
| `hk6` | `1.76e-27 J` | Source sixth-order field; energy mapping unresolved |
| `Dsim` | `5.6e-24 J` | DMI enabled; zero if `DMIenable=0` |
| `mus` | `4.2313*mub = 3.92410762e-23 J/T` | Same moment for every occupied spin |
| `a,c` | `5.067e-10,1.3882e-9 m` | Density and thickness definitions above |
| `ms,tz` | Approximately `1.525589445e6 A/m`, `1.3882e-9 m` for six layers | Torque normalization |
| `g,alp` | `2,0.001` | `gam=g*mub/(hbar*ele)`; damping unchanged |
| `jcSOT,jcSTT` | `6e12,0 A/m^2` | Constant drive settings, no pulse added |
| SOT flags | DLT=1, FLT=1 | Both enabled |
| SOT polarization | theta=90 degrees, phi=210 degrees; approximately `[-sqrt(3)/2,-1/2,0]` | Same polarization supplied to all occupied sites |
| `thetaSH,chi` | `0.2,5` | `BDSOT` approximately 0.18648 T; `BFSOT` approximately 0.93239 T |
| STT flags, polarization, efficiency | DLT=0, FLT=0; `[0,0,1]`; `etaSTT=0.8` | Zero amplitude; preserve source implementation |
| `Hext` | `[0,0,0] T` | Zero external field |
| `dipolemode,dipole_tstep` | `0,50` | Disabled; unchanged refresh rule if enabled |
| Dipole distance | `0.5067e-9 m` | Existing three-neighbor in-plane truncation; justification unresolved |
| `rk4,tstep` | `1,1e-16 s` | RK4 default, fixed fields during each step |
| `gpusave,gpurun_number,savetstep` | `1e-12 s,150,5000` | Nominal 150 ps; 300 saved frames; accepted timing defects |
| Distribution flags | Both zero | Deterministic occupancy; optional restart loads `Relax_0.mat`; save flag is ineffective |
| Thermal/stochastic field | None | Unused `thermalenable` and `T` removed; existing `rng('shuffle')` does not make the occupancy random |
| Constants | `mub=9.274e-24 J/T`, `ele=1.602e-19 C`, `hbar=6.58211951440e-16 eV s`, `mu_0=1.25663706143592e-6 T m/A`, `kb=1.38064852e-23 J/K` | Source constants unchanged; `kb` remains unused |

## Effective Fields and Units Check

Use effective flux density $\mathbf B_i=-\mu_s^{-1}\partial\mathcal H/
\partial\mathbf s_i$ in tesla [4]. The retained exchange and second-order
fields are $\mathbf B_{\mathrm{ex},i}=-\mu_s^{-1}\sum_jJ_{ij}\mathbf s_j$
and $\mathbf B_{2,i}=(0,0,2K_1s_{iz}/\mu_s)$.

The imported anisotropy field is exactly source `integrate_llg.m`, lines
162-176. Define $q=\sin(\arccos s_z)$, $u=s_x/q$, $v=s_y/q$; the source
sixth-order numerators are

$$
A_x=12h_6q^5uv^2(2v+1)(2v-1)(-3+4v^2),\quad
A_y=-12h_6q^5u^2v(2v+1)(2v-1)(-3+4v^2),
$$

$$
A_z=-6h_6s_zq^4(-1+2u^2)(16u^4-16u^2+1).
$$

The source multiplies these by `mask=~(abs(mz)==1)` and divides by `mus`,
then adds $2K_1s_z/\mu_s$ to z. These are a transcription of code, not a
new energy derivation. Retain the confirmed +z NaN defect (ANI-01).

Writing the signed source sums as $\mathbf S_{\mathrm{nex}}$ and
$\mathbf S_{\mathrm{pre}}$, the imported DMI field is

$$
\mathbf B_D=\frac{D}{\mu_s}
(S_{\mathrm{nex},y}+S_{\mathrm{pre},y},
-S_{\mathrm{nex},x}-S_{\mathrm{pre},x},0).
$$

Correct only the unshifted adjacent-layer references from trajectory buffers
`mmx/mmy/mmz` to `mmxtmp/mmytmp/mmztmp`. Keep the source masks, signs,
coefficients, offsets, unconditional DMI call, and zero-amplitude disabling
behavior. The optional dipole formula remains its existing tensor sum with
$\mu_0/(4\pi)$ and the source distances.

For energy-valued coefficients, $\mathrm J/(\mathrm{J/T})=\mathrm T$.
The `hbar` constant is numerically $\hbar/e$ in V s for torque amplitudes;
the source gamma formula multiplies it by `ele` to obtain J s. Preserve
these conventions, including damping-dependent torque cross terms. Material
normalization and energy/field consistency gaps remain visible.

## Dynamical Equation and Nonconservative Torques

For each SOT/STT channel $a$ with polarization $\boldsymbol\sigma_a$ and
field-equivalent amplitudes $B_{D,a},B_{F,a}$, the retained kernel implements
[1, atomgpurk4.m; 4, eq:ll-with-flt]

$$
\dot{\mathbf s}=\frac{\gamma}{1+\alpha^2}\left[
-\mathbf s\times\mathbf B-\alpha\mathbf s\times(\mathbf s\times\mathbf B)
+\sum_a\left\{(\alpha B_{D,a}+B_{F,a})\mathbf s\times\boldsymbol\sigma_a
+(-B_{D,a}+\alpha B_{F,a})\mathbf s\times(\mathbf s\times\boldsymbol\sigma_a)
\right\}\right].
$$

SOT amplitudes are the source $B_D=(\hbar/e)\theta_{SH}j_c/(2M_st_z)$ and
$B_F=5B_D$. STT is zero under the accepted defaults. Preserve the existing
source `BFSTT=chi*BDSTT` expression; repairing its use of the SOT ratio rather
than `chiSTT` is outside this import and has no effect with zero STT.
No noise, stochastic interpretation, or variable-length dynamics is added.

## Initial Conditions, Solver, and Outputs

Use source initialization: occupied odd layers along +x and even layers
along -x, up to floating-point trigonometric residuals; no lateral domain
wall. Optional `load_fixed_atom_distrib=1` loads the last spin frame from
`Relax_0.mat`, which is not supplied. Default loading is off.

Default RK4 keeps effective fields fixed across substages and normalizes
the result. Numerical equations in all three kernels remain unchanged.
Heun/PC torque limitations and the PC bootstrap signature mismatch remain.

Nominal duration is 150 ps at 0.1 fs; actual accumulated updates span
149.985 ps because TIME-01 is intentionally deferred. Sample stride 5000
gives 300 frames. Preserve the source `t` construction, initial/final sample
policy, and DMI-history off-by-one/zero slots (OUT-01).

Production output remains the source complete-workspace `save('final.mat')`
after its cleanup, including `mmx/mmy/mmz`, `t`, and `hdmi_allx/y/z`.
The helper `One_data_curve.m` runs manually on loaded workspace variables;
its fixed site indices and six-panel arrangement remain unchanged.
Verification must use a separate scratch directory and must never overwrite
root `final.mat`, `finalnew.mat`, or another existing result.

## Implementation Mapping

Replace the old root `plan.md` with this approved plan, without archiving
the old audit, as explicitly requested by the user. Its pre-replacement
SHA-256 was `27088352d0dbb78e10aa101a5eaca4366596b868f8260e5558a9c9db128200b4`.

Import the 12 `.m` files and `Readme.txt` from the manifest below. Preserve
local LF endings in existing MATLAB files to avoid unrelated whole-file
formatting churn; retain source endings in newly added files. The source
directory stays read-only. No local code file needs deletion.

| Stage | Implementation | Verification gate |
|---|---|---|
| 0. Protect inputs | Recheck branch/status, source hashes and destination changes; replace old plan as requested | Approved plan installed, no unexplained source/destination changes |
| 1. Structure | Import `distrib.m` and `systemgeneration.m` | On a 5 x 9 x 6 core, compare occupancy/padding/all six masks with pre-import code; check source initialization and zeros outside occupancy |
| 2. Interaction topology | Import renamed shifts, exchange/dipole scripts, and source DMI | All 30 shift definitions retain their offsets; exchange/dipole text is exactly equivalent after renaming; enumerate DMI coefficients across all masks and surfaces and report antisymmetry separately from physical validity |
| 3. Fields | Import source field formulas; apply only approved DMI current-state correction | DMI uses current-state arrays for every added adjacent-layer term; controlled stale-history test proves independence from trajectory history; anisotropy source equality and documented pole failure preserved |
| 4. RHS/torques | Import coefficients and source kernels | Kernel equations unchanged; one-spin CPU checks against explicit vector RHS, precession/damping direction and finite unit norm away from known singular states |
| 5. Integrator/output | Import timing, explicit gathers, cleanup, output, helper/readme | Source text equality except DMI correction; index trace reproduces deferred OUT-01/TIME-01; small GPU smoke if available uses isolated output |
| 6. Documentation | Update current behavior descriptions and mark older explorer snapshot historical | Accepted defaults/issues consistent across docs; links work; no claim that imported DMI or anisotropy is independently validated |

Complete each transfer/correction gate before advancing. A failed import or
approved indexing-correction check blocks completion. Scientific checks that
expose a pre-existing/source defect must be reported; do not silently repair
it. Explicitly deferred defects are expected known limitations, not pass
criteria for scientific correctness.

Preserve the script architecture and data-flow separation already present:
occupancy/initialization -> shifts and conservative fields -> torque inputs
and dynamical kernel -> integration and output. No geometry or coefficient
derivation moves into the dynamical kernel.

Documentation changes: update `README.md`, `docs/model-methodology.md`, and
`diff.md`. Correct only the stale model-contract/output statements in the
project `AGENTS.md`; keep its workflow rules and all global instructions.
Mark `docs/structure-exchange.html`, its guide, and the notebook atlas link
as pre-import historical snapshots, including their old shift names, source
line references and DMI signs. Do not redesign or reinterpret the explorer.
Keep `LICENSE`, notebook assets, and user results. Add ignores for the
verification scratch directory and the standard generated/cache directories
required by global instructions; do not change unrelated ignore rules.

## Verification and Acceptance Criteria

Already completed: normalized source comparison; all 30 shift-name mappings;
exchange/dipole equivalence; source constant calculations; a MATLAB
three-index/4-D-buffer probe; exact anisotropy-block probe reproducing ANI-01;
and the output index trace. These are not a full scientific validation.

After plan approval, use a disposable `.review-20261013/` harness with old
code snapshots and isolated results, ignored by Git. Do not add external
dependencies. The installed `matlab` launcher invokes Windows MATLAB, so
file arguments use Windows-visible paths. Planned invocation:

```bash
matlab -batch "run('C:/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/.review-20261013/verify_import.m')"
git diff --stat
git diff --check
```

The harness will perform the staged checks above without invoking production
`main.m`. Use exact equality for masks, renamed operator text, and imported
source text after line-ending normalization. For the DMI correction, assign
different current and first-history states; independently calculate the
expected contribution of each changed adjacent-layer term on all masks and
compare in tesla with tolerance $10^{-12}\max(1,\max|B|)$, then change only
history values and require no output change. Include layer surfaces and
lateral edges. This checks state selection and implemented algebra, not
material provenance.

Use a normalized one-spin state away from the poles, fixed fields, and
isolated torque channels to compare the kernel with the explicit RHS at a
small timestep; verify finite values and norm error below $10^{-12}$ in
double precision. Document MATLAB/GPU availability. If GPU support works,
run only a small deterministic multi-chunk smoke with the accepted source
formulas and a distinct scratch output path, retaining the original
production defaults in repository files. Otherwise report GPU execution as
unverified. Retain small-run artifacts for review; no accepted baseline is
created or replaced.

No new claims about published ground states, physical bond directions,
anisotropy energy, timestep convergence, or boundary/size convergence are
required for this import. Those scientific gates remain unresolved and are
not inferred from matching trajectories or plausible switching.

Source-byte SHA-256 manifest at review:

```text
One_data_curve.m  9cc42af0b3d74537ed000a76509714d1849e66328d275f6274edec52452f2de4
Readme.txt        4ae5f4aca77456ec927fb990f9400ddc29736deb67f3008f06f7803db7cea41e
atomgpu.m         6e760d47eaa34989a69def00bc97c56887652b710ac9eed74b8e837768f08816
atomgpupc4.m      a7fee072388b81853dc21e4c86f341e04b77115d1cd823144f4ec947c063a655
atomgpurk4.m      52c3d6dba363859ac2d222b614701e57107371886925b48a8c083622f7105a0f
constantfile.m    17b95d335c43a399f9a124d2e9fde2270b1ad0402cfc49ec4994f641cebd0f9f
dipole_.m         5b4437e33199c6a38a41adc8ff296e01a1563bc60afae4aacf24306b7d570502
distrib.m         ddecda1398d9eeeef68c01cc2eb1e65454492454b8e3307813e51e75431941ca
dmi.m             63cac3ee2c9ba8b53cc904edf81378a6a3a604e6ec3720508cf1d1a4c00caf44
exchangej_.m      2a873e3a8ded20270581c71115f385d903d900ffc03f5738d6119a5061f21073
integrate_llg.m   384d13ec6205c821be8f33912bd36471d09849fe3c0ea792a9bd4a7f33acd51a
main.m            40fc30a9bf230a5139e80bd107f809716d6415f67bb24a0c248e4e78d7ebd38d
systemgeneration.m 334c1ea1bfc25665de1aac171099cabf4ac85c57d965db558f2bffe1ffae03da
```

## References

1. User-supplied `code_20261013`, exact folder and hashes above. Read directly:
   all compared MATLAB files, `Readme.txt`, and `One_data_curve.m`. Evidence
   of executable source behavior only, not a qualifying physical reference.
2. Active repository at commit `d0ab252d55ef70328be2f45ea1cba51201288522`:
   `README.md`, `docs/model-methodology.md`, `AGENTS.md`, and MATLAB sources.
   Existing implementation/validation record only.
3. This conversation's ten explicitly accepted review items, recorded in
   `diff.md`, including the sole DMI code correction and deferred ANI-01,
   OUT-01, and TIME-01.
4. Local equation reference
   `C:\Users\zzf-m\OneDrive\presentations\equations\spintronic_methods.tex`;
   `sections/hamiltonian.tex` (eq:atomistic-hamiltonian,
   eq:atomistic-effective-field, eq:atomistic-pair-fields) and
   `sections/llg_equation.tex` (eq:llg-gilbert, eq:ll-with-flt).
   Conventions only; not material-specific evidence.
5. Processed paper catalog `_knowledge/papers.jsonl`, record ID
   `AFM/2023-Magnetic properties of hematite revealed by an ab initio-PRB-Tobias Dannegger`,
   title: *Magnetic properties of hematite revealed by an ab initio
   parameterized spin model*, DOI `10.1103/PhysRevB.107.184426`.
   Catalog metadata was inspected. The indexed summary and source PDF lie
   outside the pre-authorized `_knowledge` subtree and were not opened.
   This is a discovery record, not an inspected journal reference supporting
   any adopted coefficient or structure claim. Full author metadata and exact
   material locators are not supplied by this record and are not invented.
6. ASTRO reference repository
   `C:\Users\zzf-m\OneDrive\code_softwares\spin_dynamics\atomistic_model\2D\atomistic\atomistic_v2`:
   `AGENTS.md` and `README.md`, inspected for workflow/context only. Its
   ferrimagnetic parameters, algorithms, tests and accepted baseline do not
   validate the present hematite import.

## Assumptions, Missing References, Risks, and Approval Questions

The accepted decision ledger is authoritative for this transfer. No
material-specific claim is promoted from UNRESOLVED to validated. Evidence
search covered active/source documentation, the equations reference, and
the relevant processed-paper catalog entry. No new online search, source-PDF
access, or publication reconstruction was performed. Every missing-reference
row above remains an explicit limitation; further physical reconstruction,
parameter justification and energy validation require a separate task.

Approved decisions: all ten source-difference items; current-state indexing
correction only in DMI; ANI-01, OUT-01 and TIME-01 explicitly deferred.
The optional PC path, absent restart file, fixed-field substages, truncated
dipole approximation, and existing nondefault STT flag/ratio limitations
remain documented. No additional fix is implied by a successful import.

Approval received for revision 1 and all ten decisions. The user's final
instruction was: "Approve this plan, remove the old plan.md". Accordingly,
replace the old audit with this active plan without an archive. Proceed with
the import, sole DMI correction, scoped documentation updates, and small
checks. No further approval is needed unless the implementation scope changes.

## Execution Record

Implementation completed on 2026-10-15 under the approved decisions. The old
root audit was removed without an archive at the user's explicit request.
Source text equality passed for all 13 imported files, allowing line-ending
normalization and exactly 18 approved DMI current-state substitutions.

Structure, state-selection, RK4 vector/analytical, and small GPU/output checks
passed. The independent topology diagnostic found 544 directed DMI coefficient
antisymmetry mismatches on the 140-site lattice; this is recorded as DMI-02
without changing the accepted source signs. ANI-01, OUT-01 and TIME-01 remain
deferred. No full production run or scientific convergence claim is made.

See docs/model-methodology.md, Import verification, and diff.md for exact
commands, parameters, hardware, results, and remaining limitations. The
verification harness and outputs are retained only in the ignored local
.review-20261013 directory. No commit or push was performed.
