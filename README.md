# Alpha-Fe2O3 atomistic spin dynamics

This repository contains a MATLAB/GPU research implementation of an atomistic
spin-dynamics model for $\alpha$-Fe$_2$O$_3$ (hematite). It represents
fixed-length spins on a material-specific, three-layer repeating site pattern
and evolves them with LLG/LLGS-style dynamics.

The code is a research prototype, not yet a validated simulation package. In
particular, the crystal-to-array mapping, interaction topology, parameter
provenance, boundary behavior, and numerical convergence still require focused
tests. See [docs/model-methodology.md](docs/model-methodology.md) for the
implemented equations, assumptions, evidence status, and open issues.

The implementation was imported from `code_20261013` after a ten-item review.
The sole correction during that initial import was that added adjacent-layer DMI terms
read the current spin state rather than the first frame of each chunk.
Merged PR #18 (`864d6b0`) subsequently corrected the anisotropy pole division,
DMI-history alignment, chunk update counts, and saved-time labels, and replaced
DMI with separate J3/J4 weights and alternating layer signs. See
[diff.md](diff.md) for current issue status and the historical import decisions.
Merged PR #20 (`d79ef2a`, reviewed head `cf19a99`) moves field evaluation into
`field_calc.m`, reevaluates exchange, anisotropy, and DMI at RK4 substages,
and keeps the evolving spin state on the GPU. It also changes the default
to a 25 ps zero-current run saved as `Relax.mat`.

## Requirements

- MATLAB with Parallel Computing Toolbox.
- A MATLAB-supported NVIDIA GPU.
- Enough CPU and GPU memory for the selected lattice and saved trajectory.

Focused import checks use MATLAB R2024b Update 1. Their scope and hardware
results are recorded in [the methodology](docs/model-methodology.md#import-verification).
The accepted short regression suite records MATLAB/toolbox/GPU versions and source
hashes; a full-duration production baseline is still unavailable.

## Running the model

Run from the repository root so that project scripts are on the MATLAB path.
`main.m` runs the unchanged 25 ps defaults and requests `Relax.mat`; an existing
file causes an error before integration:

```matlab
cd('/path/to/Fe2O3')
main
```

For a customized run, edit a complete configuration returned by the defaults
function. This short example returns the result without saving:

```matlab
cfg = fe2o3_default_config();
cfg.natomW = 5; cfg.natomL = 9; cfg.natomH = 6;
cfg.gpusave = 8e-16; cfg.gpurun_number = 3; cfg.savetstep = 2;
result = fe2o3_run(cfg);
% Choose an existing output directory and a fresh filename to save instead:
% result = fe2o3_run(cfg, 'my-short-run.mat');
```

`fe2o3_validate_config` rejects missing/unknown fields, malformed or nonfinite
inputs, incomplete three-layer motifs, nonintegral step counts and unsupported
torques. Numeric inputs are converted to double before coefficient arithmetic.
`gpusave/tstep` must be integral to floating-point tolerance and divisible by
`savetstep`; `psjSTT` must be a unit row vector within `1e-12`.
Constants remain in `constantfile.m`; derived values are not editable inputs.
The most important controls are:

| Setting | Meaning | Current default |
|---|---|---:|
| `rk4` | `1`: RK4, `0`: Heun, `2`: fourth-order predictor-corrector | `1` |
| `DMIenable` | Enable the implemented DMI field | `1` |
| `dipolemode` | Enable the truncated in-plane dipole field | `0` |
| `natomW`, `natomL`, `natomH` | Unpadded array dimensions | `20, 20, 6` |
| `tstep` | Integration step in seconds | `1e-16` |
| `gpusave` | Length of one GPU integration chunk in seconds | `1e-12` |
| `gpurun_number` | Number of chunks | `25` |
| `savetstep` | Save every Nth in-chunk state | `5000` |
| `jcSOT`, `jcSTT` | Drive current densities in A/m² | `0, 0` |
| `load_fixed_atom_distrib` | Load last spin frame from `Relax_0.mat` | `0` |

No thermal field or domain-wall initialization option is implemented in this
version. The atom distribution is deterministic; the retained normal branch still
calls `rng('shuffle')` without drawing random model inputs. The legacy `save_fixed_atom_distrib` flag
has no save action. `Relax_0.mat` is not supplied; enabling the restart option
requires `cfg.restart_file` (default `Relax_0.mat`) with either legacy
`mmx/mmy/mmz` double histories or the new `fe2o3_result` struct. All frames must
have compatible padded dimensions, finite real values, unit occupied spins
within `2e-12`, and zero empty sites. New results must also match occupancy,
all six masks and array ordering. The last saved frame is selected without
normalization; this is a spin restart, not a full solver/cache checkpoint.
RK4 is the default and has focused verification. Heun and predictor-corrector
remain legacy, unvalidated alternatives with fixed-field approximations.
The predictor-corrector bootstrap now shares the first three RK4 updates;
the former argument mismatch is removed. Nonzero SOT or STT amplitudes are
rejected unless `rk4=1`; enabled switches with zero current are allowed.

The callable API is `fe2o3_run(cfg, output_file, overwrite)`, with defaults
`output_file=''` (in memory) and `overwrite=false`. Explicit replacement is
`fe2o3_run(cfg, 'chosen-file.mat', true)`. The destination's parent directory
must already exist. Accepted-reference storage and the active restart input
are protected even with replacement requested, including resolved path aliases.
A completed temporary MAT file is published by exclusive hard-link creation;
explicit replacement uses an atomic move. A filesystem lacking those operations
fails without a fallback that could clobber an existing file. Windows NTFS was tested.

Saved MAT files now contain only `fe2o3_result`, schema `1.0`. Arbitrary legacy
workspace temporaries are no longer public outputs:

| Field | Contents |
|---|---|
| `configuration`, `derived` | Validated nominal inputs; effective constants, coefficients, dimensions, counts and clock |
| `geometry`, `initial_state` | Occupancy, all six masks, array order and actual starting spins |
| `trajectory` | Host arrays `mmx/mmy/mmz`, `hdmi_allx/y/z` and `t` |
| `units`, `sampling`, `model`, `limitations` | Tesla/seconds conventions, sample policy, solver/cache behavior and unresolved model issues |
| `provenance` | UTC time, output arguments, Git state, SHA-256 of root MATLAB sources including untracked helpers, runtime/GPU identity, RNG states, restart identity and elapsed time |

The runner rejects shadowed project kernels and source/Git changes during a run.
Git must be available for provenance. No accepted raw MAT files are committed;
transfer the separate baseline suite when preparing another computer.

To plot a loaded result, run `load('my-short-run.mat'); One_data_curve`.
Both the new struct and historical workspace files are supported. It uses the spins at array
indices `(8,4,2)` and `(8,4,3)`, the latter's time trace, and a fixed six-panel
layout for final y-component maps. It is not called automatically by `main`.

## Default configuration

The defaults enable DMI. Both SOT switches remain enabled, with polarization
approximately `[-0.866,-0.5,0]` and field-like/damping-like ratio 5, but
`jcSOT=0` makes both torque amplitudes zero. STT and dipole coupling are
disabled, and `Hext=[0,0,0] T`. Occupied odd/even layers
start along +x/-x. Anisotropy uses `Ksim1=-8e-22 J` and `hk6=1.76e-27 J`.
The default run performs 250,000 updates, evolving for 25 ps. Its 50 saved
spin and DMI frames are aligned at `t = 0, 0.5, ..., 24.5 ps`; the 25 ps
endpoint is computed but excluded from the saved trajectory.

The lattice generator adds two empty array rows and columns on every lateral
side. Thus the default arrays become $24\times24\times6$ after generation.
Only entries selected by `atomtype_` represent magnetic sites.

**ANI-01** and **OUT-01** are resolved in the merged code. **TIME-01** update
counts and timestamps are corrected; saving the final-time state remains an
open output decision. Restarting from the last saved default frame uses the
state at 24.5 ps. The **DMI-02** antisymmetry defect is resolved for the tested
finite lattices, while physical bond assignments and parameter provenance
remain unresolved. See [the issue record](diff.md) for evidence and scope.

## Repository layout

| Path | Responsibility |
|---|---|
| `main.m` | Thin default-run script requesting protected `Relax.mat` |
| `fe2o3_default_config.m`, `fe2o3_validate_config.m` | Canonical inputs and complete validation/derived coefficients |
| `fe2o3_run.m` | Local-scope orchestration, explicit result assembly and provenance |
| `fe2o3_load_restart.m`, `fe2o3_output_path.m`, `fe2o3_save_result.m`, `fe2o3_sha256.m` | Restart validation, output ownership, protected publication and byte identity |
| `tests/`, `benchmark/` | Maintained checks, accepted-reference manifest and fresh candidate comparison |
| `constantfile.m` | Physical constants used by the scripts |
| `systemgeneration.m` | Alternating initial spin state or optional restart |
| `distrib.m` | Occupancy and layer/site masks plus lateral padding |
| `integrate_llg.m` | GPU state, RK4 stages, legacy solver dispatch, normalization, and sampling |
| `field_calc.m` | Neighbor shifts, exchange/anisotropy/DMI fields, cached dipole field, and total field |
| `exchangej_.m` | Site-dependent $J_1$--$J_4$ neighbor sums |
| `dmi.m` | Alternating layer signs for J1/J3/J4 DMI sums; J1 weight is fixed at zero in the validator |
| `dipole_.m` | Truncated nearest-neighbor in-plane dipole field |
| `atomgpurk4.m` | LLGS derivative with SOT/STT-labeled torques; stages are assembled by `integrate_llg.m` |
| `atomgpu.m` | Fixed-field Heun update without torque arguments |
| `atomgpupc4.m` | Legacy fixed-field predictor-corrector update; order unvalidated |
| `One_data_curve.m` | Manual plots of selected spins and final layer maps |
| `Readme.txt` | Original source-folder file descriptions |
| `diff.md` | Current issue status, historical import decisions, and verification record |
| `Alpha-Fe2O3.pdf` | Original annotated methodology notebook |
| `docs/Alpha-Fe2O3-notebook-transcription.md` | Page-addressable digital transcription of the notebook |
| `docs/Alpha-Fe2O3-notebook.html` | Interactive color-preserving notebook atlas with annotation fields |
| `docs/Alpha-Fe2O3-pages/` | Rendered notebook pages used by the interactive atlas |
| `docs/structure-exchange-guide.md` | Current illustrated mask-to-array explanation and all 24 color/exchange mappings |
| `docs/structure-exchange.html`, `docs/structure-exchange.js` | Offline explorer with current source excerpts and signed DMI operator |
| `docs/unit-cell.html`, `docs/unit-cell.js` | Rotatable notebook-based 3D schematic with cell repetition; illustrative coordinates, not simulation geometry |
| `docs/model-methodology.md` | Maintained description of implemented behavior and evidence gaps |
| `AGENTS.md` | Project-specific instructions for coding agents |

## Methodology and evidence

For an introduction to `distrib.m` and `exchangej_.m`, start with the
[illustrated implementation guide](docs/structure-exchange-guide.md) or open
[the interactive explorer](docs/structure-exchange.html) in a browser. It also
appears in the notebook atlas. The explorer and guide follow the current
MATLAB masks, exchange topology,
shift-variable names, source excerpts, and the merged J3/J4 DMI operator.
Select **DMI** and either actual layer index for each color to see the
alternating signs, `Dsim3=1.6e-23 J` and `Dsim4=1.28e-23 J` weights, and
surface/edge omissions. J1 has zero configured DMI weight; J2 DMI is absent.
Selecting a neighbor traces
its shift through `exchangej_.m` and its weight in `dmi.m`.
DMI labels show algebraic weights rather than physical bond-vector arrows;
physical shell assignments and parameter provenance remain unresolved.
Notebook images and transcription remain historical records.
No server, build step, or external dependencies are needed; keep the companion
JavaScript file and notebook page images beside the HTML when sharing it.
The explorer illustrates current implementation indexing and historical
exchange/notebook agreement; physical crystal geometry remains unverified.

The [3D unit-cell schematic](docs/unit-cell.html) redraws the Fe-column figure
on notebook page 1. Drag to rotate, use the camera controls to zoom or change
view, and repeat the prism 1–4 times along each of the three directions.
The author requested a schematic: spacing and aspect ratio are illustrative,
and the displayed labels are not assigned to MATLAB sites. This offline view
also appears in the atlas and does not modify or run the simulation.

[Alpha-Fe2O3.pdf](Alpha-Fe2O3.pdf) is retained as the original visual record.
The PDF is ignored by Git; its rendered page images are included under `docs/`.
The atlas and explorers work with those images, while links to the original PDF
require a local copy named `Alpha-Fe2O3.pdf` in the repository root.
It contains useful crystal diagrams, layer-to-array sketches, derivations, and
neighbor maps, but also tentative annotations and abbreviated references. It
is not a substitute for the original publications or for executable
verification.

The maintained methodology therefore distinguishes among:

- behavior directly observed in the MATLAB implementation;
- interpretations recorded in the notebook;
- claims supported by identified publications; and
- unresolved claims that still need an exact source locator or a numerical
  test.

Do not silently change physics to make the code agree with the notebook or a
paper. Report the conflict, define the intended convention, and validate the
change first.

## Current validation status

A separate [PR #18 / TIME-01 change report](docs/pr18-time-01-change-report.html)
compares `4ba0e42` with `77dece1`, including checked timing behavior and the
remaining final-state output question. That reviewed MATLAB source is now
superseded in part by PR #20. The report retains its original pre-merge context;
[historical verification status](docs/model-methodology.md#merged-pr-18-verification)
summarizes which checks support the resolved issues.

A separate [PR #20 change report](docs/pr20-change-report.html) compares
`7392314` with the review snapshot `cf19a99`, now merged at `d79ef2a`.
It shows all six changed files, explains the default-run and solver changes, and embeds the focused
RK4 checks and the cached-dipole limitation. Current numerical statements match
that source; `dmi.m` additionally contains explanatory comments. The report
retains its original pre-merge wording. See
[PR #20 verification scope](docs/model-methodology.md#merged-pr-20-verification)
for the recorded results and remaining solver limitations.

Maintained focused CPU/GPU checks and a protected reference-preparation suite
are now available; see [benchmark/README.md](benchmark/README.md) for commands
and [the current verification report](docs/cq01-05-verification.md) for evidence.
The accepted `cq-20261120-v1` suite covers six short cases. The new interface
matches every saved spin/DMI sample and all effective inputs exactly on the
tested environment. CPU, GPU, input/restart/output checks and 23 comparator/path
rejection checks pass. Full-duration and physical size convergence remain open.
Before using results quantitatively, at minimum verify lattice and
site counts, neighbor multiplicities, field signs and units, spin-norm
preservation, boundary coordination, time/output alignment, and timestep and
size convergence. The detailed checklist is in
[docs/model-methodology.md](docs/model-methodology.md#validation-requirements).
