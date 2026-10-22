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

## Requirements

- MATLAB with Parallel Computing Toolbox.
- A MATLAB-supported NVIDIA GPU.
- Enough CPU and GPU memory for the selected lattice and saved trajectory.

Focused import checks use MATLAB R2024b Update 1. Their scope and hardware
results are recorded in [the methodology](docs/model-methodology.md#import-verification).
No reproducible environment manifest or full production benchmark is provided.

## Running the model

The entry point is `main.m`. Run it from the repository root so that all helper
scripts are on the MATLAB path:

```matlab
cd('/path/to/Fe2O3')
main
```

Edit the configuration block near the beginning of `main.m` before running.
The most important controls are:

| Setting | Meaning | Current default |
|---|---|---:|
| `rk4` | `1`: RK4, `0`: Heun, `2`: fourth-order predictor-corrector | `1` |
| `DMIenable` | Enable the implemented DMI field | `1` |
| `dipolemode` | Enable the truncated in-plane dipole field | `0` |
| `natomW`, `natomL`, `natomH` | Unpadded array dimensions | `20, 20, 6` |
| `tstep` | Integration step in seconds | `1e-16` |
| `gpusave` | Length of one GPU integration chunk in seconds | `1e-12` |
| `gpurun_number` | Number of chunks | `150` |
| `savetstep` | Save every Nth in-chunk state | `5000` |
| `load_fixed_atom_distrib` | Load last spin frame from `Relax_0.mat` | `0` |

No thermal field or domain-wall initialization option is implemented in this
version. The atom distribution is deterministic despite the console message
referring to a random distribution. The legacy `save_fixed_atom_distrib` flag
has no save action. `Relax_0.mat` is not supplied; enabling the restart option
requires a compatible file with `mmx`, `mmy`, and `mmz` histories.
The predictor-corrector selection currently has an argument mismatch in its
RK4 bootstrap call and must not be treated as runnable until repaired and
tested.

By default, the run writes `final.mat` in the repository root and overwrites
an existing file with that name without prompting. Use a separate working copy
or change the output path before running cases whose existing output must be
preserved. The file is a saved MATLAB workspace, so it contains more than the
trajectory arrays alone, including `mmx/mmy/mmz`, `t`, and DMI histories
`hdmi_allx/y/z`. The source cleanup removes the six site-specific masks and
several temporary arrays before saving.

To plot a loaded result, run `One_data_curve`. It uses the spins at array
indices `(8,4,2)` and `(8,4,3)`, the latter's time trace, and a fixed six-panel
layout for final y-component maps. It is not called automatically by `main`.

## Default configuration

The defaults enable DMI and both SOT components, with `jcSOT=6e12 A/m^2`,
polarization approximately `[-0.866,-0.5,0]`, and field-like/damping-like
ratio 5. Their field-equivalent amplitudes are approximately 0.1865 T and
0.9324 T. STT and dipole coupling are disabled. Occupied odd/even layers
start along +x/-x. Anisotropy uses `Ksim1=-8e-22 J` and `hk6=1.76e-27 J`.
The default run performs 1,500,000 updates, evolving for 150 ps. Its 300 saved
spin and DMI frames are aligned at `t = 0, 0.5, ..., 149.5 ps`; the 150 ps
endpoint is computed but excluded from the saved trajectory.

The lattice generator adds two empty array rows and columns on every lateral
side. Thus the default arrays become $24\times24\times6$ after generation.
Only entries selected by `atomtype_` represent magnetic sites.

**ANI-01** and **OUT-01** are resolved in the merged code. **TIME-01** update
counts and timestamps are corrected; saving the final-time state remains an
open output decision. Restarting from the last saved default frame uses the
state at 149.5 ps. The **DMI-02** antisymmetry defect is resolved for the tested
finite lattices, while physical bond assignments and parameter provenance
remain unresolved. See [the issue record](diff.md) for evidence and scope.

## Repository layout

| Path | Responsibility |
|---|---|
| `main.m` | Configuration, orchestration, cleanup, and workspace output |
| `constantfile.m` | Physical constants used by the scripts |
| `systemgeneration.m` | Alternating initial spin state or optional restart |
| `distrib.m` | Occupancy and layer/site masks plus lateral padding |
| `integrate_llg.m` | Neighbor shifts, effective fields, GPU loop, and sampling |
| `exchangej_.m` | Site-dependent $J_1$--$J_4$ neighbor sums |
| `dmi.m` | J3/J4 DMI weights with alternating target-layer signs |
| `dipole_.m` | Truncated nearest-neighbor in-plane dipole field |
| `atomgpurk4.m` | Fixed-field RK4 update with SOT/STT-labeled torques |
| `atomgpu.m` | Fixed-field Heun update without torque arguments |
| `atomgpupc4.m` | Fixed-field fourth-order predictor-corrector update |
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
surface/edge omissions. DMI on J1/J2 is absent. Selecting a neighbor traces
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
merged in this checkout. The report retains its original pre-merge context;
[current verification status](docs/model-methodology.md#merged-pr-18-verification)
summarizes which checks support the resolved issues.

Focused import checks are recorded in the methodology; there is no accepted
production benchmark or maintained automated test suite. Before using results
quantitatively, at minimum verify lattice and
site counts, neighbor multiplicities, field signs and units, spin-norm
preservation, boundary coordination, time/output alignment, and timestep and
size convergence. The detailed checklist is in
[docs/model-methodology.md](docs/model-methodology.md#validation-requirements).
