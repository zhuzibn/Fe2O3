# Issue status and historical import review

## Current issue status (2026-10-29)

PR #20 is merged at `d79ef2a`; its reviewed MATLAB source (`cf19a99`)
matches this checkout. PR #18 (`864d6b0`) supplied the earlier issue fixes.
The import decisions and reproduced failures below are preserved as history.

| Issue | Current status | Change and evidence |
|---|---|---|
| ANI-01: anisotropy pole NaN | Resolved as an arithmetic defect | PR #18 guarded the pole division. PR #20 evaluates sixth-order angles from normalized stage components with `hypot` and guarded denominators in `field_calc.m`. The earlier 50,000-vector bitwise comparison applies to PR #18 only. Sixth-order energy and material provenance remain unresolved. |
| OUT-01: DMI-history zeros and lag | Recording fix retained | PR #20 saves each pre-update spin and its first DMI field evaluation at the same `sample_index`. The zero-difference output checks were recorded for PR #18; no new end-to-end output check was run here. |
| TIME-01: chunk duration and time labels | Update counts and timestamps corrected; final endpoint still unsaved | Each chunk performs `gpusteps` updates. Current defaults evolve for 25 ps and save 50 aligned frames at 0–24.5 ps. The computed 25 ps endpoint is cleared before output. |
| DMI-02: coefficient antisymmetry | Resolved on tested J3/J4 operators | PR #20 retains the PR #18 J3/J4 operator and adds a J1 term with zero configured weight. Historical finite-lattice and shell-weight checks cover the zero-J1 operator. Nonzero J1, physical bond mapping, and parameter provenance remain unvalidated. |

The default saved times are `0, 0.5, ..., 24.5 ps`. Restarting from the last
saved spin therefore uses the 24.5 ps state. Output is now `Relax.mat`,
while restart still reads `Relax_0.mat`. Saving the final endpoint with a
matching timestamp and DMI field remains an open output decision.

PR #20 also sets the default SOT current to zero, moves fields into
`field_calc.m`, reevaluates exchange/anisotropy/DMI at RK4 substages, and
keeps current spins on the GPU. `atomgpurk4.m` now returns a derivative.
The PC4 bootstrap shares the RK4 path for three updates, removing the old
signature mismatch. Both SOT and STT nonzero amplitudes are rejected for
Heun/PC4. Their legacy fixed-field formulas and the dipole cache are retained.

Evidence and exact scope are collected in
[PR #20 verification](docs/model-methodology.md#merged-pr-20-verification),
the [PR #20 snapshot report](docs/pr20-change-report.html#evidence),
[historical PR #18 verification](docs/model-methodology.md#merged-pr-18-verification),
and the [DMI visualization record](docs/structure-exchange-guide.md#current-dmi-visualization-check-2026-10-22).
The reports retain their original pre-merge wording. This documentation pass
checked source equality with `cf19a99` and refreshed source excerpts;
no MATLAB simulation or material-model benchmark was run.

Open items include crystal-to-array and exchange-shell mappings, parameter
provenance, sixth-order anisotropy energy, dipole truncation/cache accuracy,
full coupled-field convergence, and the complete legacy Heun/PC4 paths. See
[maintained limitations](docs/model-methodology.md#known-limitations-and-open-decisions).

## Historical comparison: 35bd00b and code_20261013

The following comparison and 2026-10-15 import review describe the versions
named in their headings. Their source paths, deferred decisions, and failed
checks are retained as an audit trail. Current status is summarized above.

The main difference is that `code_20261013` represents a substantially different simulation setup, even though the solver and exchange/dipole kernels are largely unchanged.

| Area                       | `35bd00b`                     | `code_20261013`                          |
| -------------------------- | ----------------------------- | ---------------------------------------- |
| Grid before padding        | `100×100×6`                   | `20×20×6`                                |
| Default initial state      | Domain wall enabled           | Uniform alternating in-plane state       |
| Dipole interaction         | Enabled                       | Disabled                                 |
| SOT                        | Disabled, `jcSOT=0`           | DLT and FLT enabled, `jcSOT=6×10¹² A/m²` |
| STT                        | Flag enabled but current zero | Disabled                                 |
| Runtime                    | 200 ps                        | 150 ps                                   |
| Saved magnetization frames | 2000                          | 300                                      |
| Output file                | `finalnew.mat`                | `final.mat`                              |

### Physics changes

- Anisotropy changed substantially in [main.m](</mnt/c/Users/zzf-m/OneDrive/projects/2026-Fe2O3 switching-xuezhong/code_20261013/main.m:29>) and [integrate_llg.m](</mnt/c/Users/zzf-m/OneDrive/projects/2026-Fe2O3 switching-xuezhong/code_20261013/integrate_llg.m:159>):

  - `35bd00b`: axial second- and fourth-order terms,
    `Ksim1=1.8073e-25`, `Ksim2=1.7624e-27`, with only `hani_z`.
  - External version: `Ksim1=-8e-22` plus a new sixth-order term `hk6=1.76e-27`, producing nonzero `hani_x`, `hani_y`, and `hani_z`.

- The magnetization and thickness definitions changed:

  - `ms = mus/d³`, `tz = 6c`
  - becomes `ms = 24 mus/(c a²√3)`, `tz = c`

  These changes also alter the SOT amplitude. With the external version’s defaults, its formulas give approximately `BDSOT=0.186 T` and `BFSOT=0.932 T`.

- The exchange constants `J1–J4`, moment `mus`, damping, gyromagnetic ratio, timestep, and zero external field are unchanged.

### Initialization and distribution

[systemgeneration.m](</mnt/c/Users/zzf-m/OneDrive/projects/2026-Fe2O3 switching-xuezhong/code_20261013/systemgeneration.m:5>) is substantially different:

- `35bd00b` defaults to a domain-wall configuration.
- The external version defaults to alternating `+x/-x` layers.
- Setting `load_fixed_atom_distrib=1` now loads magnetization from `Relax_0.mat`, instead of loading the atomic distribution from `atomtypee.mat`.
- The old atomic-distribution save/load code was removed from [distrib.m](</mnt/c/Users/zzf-m/OneDrive/projects/2026-Fe2O3 switching-xuezhong/code_20261013/distrib.m:136>).

Consequently, `save_fixed_atom_distrib` remains exposed in `main.m` but no longer saves anything, and `load_fixed_atom_distrib` no longer does what its name says.

### DMI changes

DMI was rewritten, not merely reparameterized:

- `Dsim` changes from `-4.8065e-23 J` to `+5.6e-24 J`.
- Neighbor sums now contain signed directional differences and additional self terms.
- The final field changes from “next minus previous” to signed sums:
  - old: `nex - pre`
  - new: `nex + pre`, with a separate minus sign for `hdmi_y`
- DMI calculation is now called unconditionally; disabling it only makes `Dsim=0`.
- Complete DMI-field histories `hdmi_allx/y/z` are saved.

There is a likely correctness issue in [dmi.m](</mnt/c/Users/zzf-m/OneDrive/projects/2026-Fe2O3 switching-xuezhong/code_20261013/dmi.m:14>): several self terms use `mmx(:,:,i±1)`, `mmy(:,:,i±1)`, and `mmz(:,:,i±1)`. Those are four-dimensional trajectory buffers, whereas the current state is `mmxtmp/mmytmp/mmztmp`. This appears to mix the current neighbor state with the first stored timestep of each integration block.

Also, the first saved DMI-field sample in every block remains zero because `hdmi_tmp*(:,:,:,1)` is never populated before subsampling.

### Changes that are mostly cosmetic

- `exchangej_.m` and `dipole_.m` are exactly equivalent after mapping renamed directions such as:
  - `0p1 → r1`
  - `0n1 → l1`
  - `n10 → u1`
  - `p10 → d1`
- `atomgpu.m`, `atomgpupc4.m`, and `constantfile.m` are identical.
- `atomgpurk4.m` only corrects its header from “Heun” to “rk4” and wraps the function declaration; its numerical equations are unchanged.
- CRLF versus LF line endings were ignored.

### Additional files

The external folder adds:

- [One_data_curve.m](</mnt/c/Users/zzf-m/OneDrive/projects/2026-Fe2O3 switching-xuezhong/code_20261013/One_data_curve.m:1>) for plotting selected spins and layer maps.
- `Readme.txt`.

It omits the repository `LICENSE`.

One further numerical concern: the new anisotropy calculation divides by `sin(acos(mz))`. At `|mz|=1` this is singular, and multiplying the result by a mask afterward will not remove MATLAB `NaN` values.

Model classification: Class 4 material-specific atomistic, fixed-length LLGS. I performed a static source comparison and rename-normalization check; no MATLAB simulation or physical-validity benchmark was run. I consulted the Hamiltonian, atomistic-boundary, and LLGS equation references to distinguish physics changes from bookkeeping changes. No files were modified.

Documentation: no update needed — read-only comparison.

## Interactive replacement review (2026-10-15)

This section records the subsequent interactive review. Statements above about
no file changes or no MATLAB checks describe the earlier comparison only.
All ten differences were accepted. The approved replacement has now been
implemented with the sole DMI current-state correction. The old audit in
plan.md was removed without an archive, as explicitly requested.

Source folder: `C:\Users\zzf-m\OneDrive\projects\2026-Fe2O3 switching-xuezhong\code_20261013`.

### Accepted decisions

| Item | User decision | Scope |
|---|---|---|
| 1. Atom distribution and exchange | Accepted | Adopt the source version, including removal of the optional `atomtypee.mat` occupancy save/load block. Generation, padding, and all six masks match for identical dimensions. Exchange matches after verifying all 30 renamed shift definitions. Default dimensions and spin initialization remain separate review items. |
| 2. DMI | Accepted with one correction | Adopt the source DMI formulas and `Dsim=5.6e-24 J`. Replace the added adjacent-layer references `mmx/mmy/mmz(:,:,i+/-1)` with `mmxtmp/mmytmp/mmztmp(:,:,i+/-1)` so they use the current state. These shorthand expressions describe both `i-1` and `i+1`, not literal MATLAB syntax. DMI-history sampling remains a separate output-review item. |
| 3. Anisotropy | Accepted exactly as written for now | Adopt `Ksim1=-8e-22 J`, `hk6=1.76e-27 J`, removal of the active fourth-order term, and the source anisotropy implementation. Record the pole issue below; do not apply the proposed pole-handling correction in this replacement. |
| 4. Lattice size and initial state | Accepted exactly as in source | Adopt the 20 x 20 x 6 core (24 x 24 x 6 padded), odd layers initialized along +x and even layers along -x, removal of the domain-wall initialization option, and the optional last-spin-frame restart from Relax_0.mat (off by default). |
| 5. Magnetization density and thickness | Accepted exactly as in source | Adopt ms=24*mus/(sqrt(3)*a^2*c), a=5.067e-10 m, c=13.882e-10 m, and tz=c*natomH/6. For six layers, ms is approximately 1.5256e6 A/m and tz=1.3882 nm. At equal current and torque efficiencies, the inverse ms*tz factor increases by approximately 3.60 relative to the local definitions. Physical provenance remains unverified. |
| 6. Applied current and spin torques | Accepted exactly as in source | Adopt jcSOT=6e12 A/m^2, enabled SOT damping-like and field-like switches, polarization angles theta=90 degrees and phi=210 degrees, and chi=5. Disable STT_DLT; jcSTT remains zero. thetaSH=0.2 and zero external field are unchanged. Source defaults give BDSOT approximately 0.1865 T and BFSOT approximately 0.9324 T. |
| 7. Dipole and thermal settings | Accepted exactly as in source | Set dipolemode=0 by default and allocate dipole temporary arrays only when enabled. The dipole formulas are equivalent after verified shift-variable renaming; the update interval remains 50 steps when enabled. Remove unused thermalenable and T settings; neither version implements a thermal field. |
| 8. Duration, sampling, and output | Accepted exactly as in source; both issues deferred | Adopt the nominal 150 ps run, unchanged 0.1 fs timestep, stride 5000 (nominal 0.5 ps), 300 saved spin frames, final.mat workspace output, and added hdmi_allx/y/z histories. Preserve the source recording and time bookkeeping; defer OUT-01 and TIME-01 below. The DMI current-state correction approved at Item 2 still applies. |
| 9. Solvers and GPU transfers | Accepted exactly as in source | Adopt the RK4 header correction and wrapped declaration, plus explicit gather(sxx/syy/szz) when recording CPU spin histories. Heun, predictor-corrector, and constants match after line-ending normalization; RK4 numerical equations are unchanged. Preserve the default RK4 selection, fixed fields within each nominal step, and the existing predictor-corrector bootstrap signature mismatch. |
| 10. Added files and workspace cleanup | Accepted | Import One_data_curve.m and Readme.txt; adopt source clear statements before final.mat is saved. The plotting helper uses hard-coded sites and a six-layer layout and runs manually with results already loaded. Relax_0.mat is absent from the source. Retain repository documentation, license, and instructions; update documentation affected by the replacement. |

These decisions select the requested code behavior. They do not establish
published provenance or independent physical validation of the parameters
and interaction formulas.

### Deferred issue ANI-01: anisotropy produces NaN at the +z pole

- Current status: resolved by `3e20e90`, merged in PR #18; see the status table above.
- Import-time status: recorded; correction explicitly deferred by the user at Item 3.
- Location: source `integrate_llg.m`, lines 162-176.
- Symptom: for `(mx,my,mz)=(0,0,1)`, all three anisotropy-field components
  become `NaN`.
- Cause: `cosphi=mx/sin(acos(mz))` and `sinphi=my/sin(acos(mz))`
  evaluate `0/0`; multiplying the resulting sixth-order terms by the later
  zero mask does not remove `NaN`.
- Verification: executed the exact source block from `cosphi=` through
  `hani_z=` in CPU MATLAB via `matlab -batch`, using
  `mmxtmp=[0,0,1]`, `mmytmp=[0,0,0]`, `mmztmp=[1,-1,0]`,
  `hk6=1.76e-27`, `Ksim1=-8e-22`, and `muigpu=1`.
  Rows below correspond to +z, -z, and +x; columns are the three field
  components. The unit moment isolates arithmetic and is not a material test.

```text
NaN  NaN  NaN
0    0    1.6e-21
0    0    0
```

- The exact +z failure was reproduced and the +x control was finite. The -z
  case was finite in this execution; do not report a reproduced NaN there.
- Proposed future correction, not approved for this replacement: avoid the
  angular division at the poles and explicitly set the sixth-order
  contribution there to zero, preserving the second-order field. Verify
  finite pole values and unchanged formulas away from the poles before use.
- No full simulation or physical-validity benchmark was run.

### Deferred issue OUT-01: DMI-history zeros and one-step misalignment

- Current status: resolved by `71db4fd`, merged in PR #18; see the status table above.
- Import-time status: recorded; correction explicitly deferred by the user at Item 8.
- Location: source `integrate_llg.m`, lines 77-79, 230-232, and 254-259.
- Cause: each field evaluated from spin frame `ct1` is stored at `ct1+1`.
  Field frame 1 remains at its allocation value of zero. Spin and field
  histories are then sampled using the same `1:savetstep:end` indices.
- Consequence under source defaults: saved indices within each chunk are
  `[1,5001]`. The field at index 1 is an unwritten zero; the field at index
  5001 was evaluated from spin frame 5000. Across 150 chunks, 150 of the
  300 saved field frames are unwritten zeros and the other 150 lag their
  corresponding spin states by one integration step (0.1 fs).
- Verification: inspected the source assignments and reproduced their index
  trace in Python using 10000 slots, stride 5000, and `history[ct1+1]=ct1`.
  The sampled evaluation indices were `[0,5000]`, with 0 denoting an
  unwritten slot. This was an indexing check, not a dynamics simulation.
- Proposed future correction, not approved for this replacement: store fields
  at `ct1`, aligned with the state used to evaluate them. Under the current
  stride, the unsampled final slot does not affect the saved field frames;
  any future sampling change must check final-slot coverage explicitly.
- This issue concerns recording only; it is separate from the approved DMI
  current-state indexing correction, which affects the calculated field.

### Deferred issue TIME-01: inherited chunk and saved-time mismatch

- Current status: update counts and timestamps corrected by `ab7d82e`, merged
  in PR #18; saving the final endpoint remains open. See the status table above.
- Import-time status: recorded; correction explicitly deferred by the user at Item 8.
- Present in both local and source versions.
- Location: `main.m` time-vector construction and `integrate_llg.m`
  initialization of frame 1, `while ct1<gpusteps`, chunk carryover, and
  `1:savetstep:end` sampling.
- Cause: each nominal 10000-step chunk stores its initial state at index 1
  and executes only 9999 updates. The time vector nevertheless assumes
  10000 updates per chunk and begins at `tstep`, although the first saved
  spin state is the initial state at elapsed time zero.
- Consequence: 150 chunks at `tstep=1e-16 s` evolve for 149.985 ps, not the
  nominal 150 ps. Saved time labels do not exactly match the spin states;
  the requested final-time state is not saved. The analogous local run
  evolves for 199.98 ps instead of its nominal 200 ps.
- Verification: source inspection and integer update counts:
  `150 * (10000-1) = 1499850` and `200 * (10000-1) = 1999800`.
  No full simulation was run.
- A future correction must define chunk update counts, initial/final sample
  policy, and actual elapsed-time labels together and verify them on a small
  deterministic multi-chunk case. No timing correction is authorized in this
  replacement.

### Import execution and verification (2026-10-15)

- Imported all 12 MATLAB files and Readme.txt from the reviewed source.
  Existing MATLAB line endings were preserved. After normalization, every
  imported file matches source text except 18 approved substitutions in dmi.m.
- Corrected all six added adjacent-layer references per component from
  mmx/mmy/mmz to mmxtmp/mmytmp/mmztmp. All signs, coefficients, offsets,
  source anisotropy, output recording, and timing behavior remain as accepted.
- Installed the approved replacement plan at root plan.md. The previous audit
  was removed without an archive at the user's explicit request. The approved
  plan record is also retained in docs/code-20261013-replacement-plan.md.
- No commit or push was made. No full default simulation was run and no
  existing simulation result was overwritten.
- Structure gate: MATLAB verified the 5 x 9 x 6 core (9 x 13 x 6 padded),
  140 occupied sites, identical pre/post-import occupancy and six masks,
  and alternating +x/-x occupied-layer initialization.
- Static source gate: all 30 shift definitions retain their offsets; entire
  exchange and dipole scripts match after the verified shift renames.
- DMI correction gate: all six masks, surfaces, and edges agree with an
  independently calculated correction to the stale-history source terms;
  changing trajectory history alone produces exactly unchanged corrected
  DMI sums. Field tolerance was 1e-12 times max(1, maximum field magnitude).
- RK4 gate: analytical precession and damping directions passed. The scalar
  kernel agreed with independent vector-form fixed-field RK4 in zero-torque,
  isolated SOT/STT DL/FL, and combined-channel cases; maximum component
  difference was 0 in this execution, with finite unit spins.
- Deferred ANI-01, OUT-01, and TIME-01 were reproduced/preserved; none was fixed.
- Small GPU smoke: 140 occupied sites, 2 chunks of 7 updates at 1e-16 s,
  8 saved frames, finite spins, zero empty sites, maximum norm error
  1.1102230246251565e-16. Expected source workspace cleanup and DMI-history
  zero slots were preserved.
- Hardware: MATLAB R2024b Update 1 (24.2.0.2740171), PCWIN64;
  NVIDIA GeForce RTX 3060 Laptop GPU, compute capability 8.6,
  reported total memory 6441926656 bytes.
- Commands and retained scratch artifacts are documented in
  docs/model-methodology.md, section Import verification. The GPU result is
  .review-20261013/smoke/final.mat, not a root production result.
- git diff --check reports six trailing-whitespace locations copied verbatim
  from source main.m and integrate_llg.m. Those source lines were retained
  to honor the exact-source import; documentation whitespace was checked
  separately. No line-ending-only changes were introduced in existing files.

### Historical DMI-02: antisymmetry defect corrected in PR #18

- Current status: the merged J3/J4 operator passes the recorded finite-lattice
  antisymmetry checks. The failure below describes the superseded import.
  Physical mapping and material provenance remain unresolved.
- Discovered during the import topology audit. At that time, the accepted
  source formulas were retained; no additional sign correction was authorized
  or applied within the import.
- On the 5 x 9 x 6 core, the complete signed scalar coefficient matrix C
  had 140 occupied sites, 934 directed nonzero entries, and zero diagonal.
  There were 544 directed entries for which Cij was not equal to -Cji
  (272 unordered pairs). Both missing-bond and signed-value differences
  are included in this exact matrix comparison.
- A concrete example in padded MATLAB indices is the pair (4,3,2) and
  (3,3,1), for which the imported operator gave Cij=Cji=+1. These shifted
  terms were already present in the supplied source; the import correction
  to unshifted current-state reads did not change their signs.
- This prevented treating the imported operator as a verified ordinary pairwise
  z-directed DMI field under the declared antisymmetric-bond energy
  convention. Physical bond assignment and energy consistency were unresolved
  in that audit; matching source code is not scientific validation.
- Evidence: MATLAB basis-spin probes over all 140 occupied sites, saved as
  .review-20261013/dmi-operator.mat. The test exercises all site masks,
  layer surfaces, and lateral edges. It does not establish crystallography.
