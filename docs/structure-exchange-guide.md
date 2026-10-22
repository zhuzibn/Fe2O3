# From colored sites to exchange sums

This guide follows the current checked-in `distrib.m`, `exchangej_.m`, and
`dmi.m`. Notebook pages 10–11 provide the historical exchange sketches;
page 8 records the historical J3/J4 DMI idea. The merged PR #18 code now
reuses those two neighbor sums with separate weights and alternating layer signs.
This documentation does not change the Class 4, fixed-length LLG/LLGS model
or establish its physical crystallography, shell assignments, or parameters.

Open the [interactive implementation explorer](structure-exchange.html), also
embedded in the [notebook atlas](Alpha-Fe2O3-notebook.html#implementation).
Select a target color, interaction, actual layer index, and lateral position. Select a neighbor
to connect its array position, shift variable, and current source lines.

## 1. Read the structure before reading the shifts

![Notebook page 10: colored structure, square-array encoding, and exchange annotations](Alpha-Fe2O3-pages/page-10.jpg)

The left-hand sketches identify the colored sites. The square boxes encode
their array positions. The right-hand annotations list neighbors of a target
atom identified by the annotation's color. Page 11 completes the red/black
entries; the explorer switches the source image with the selected target.

Three kinds of labels must remain distinct:

| Notation | Meaning |
|---|---|
| Notebook box labels `0` and `1`, also prefixed to neighbor lists | Labels referring to the square-box representation, as confirmed by the author; not occupancy values or multipliers of the exchange sum |
| `atomtype_(r,c,k)` | `1` for an occupied site; `0` for an empty entry |
| A color mask, e.g. `atomtype_layer1gr(r,c)` | `1` at that color's lateral positions; `0` elsewhere |

An atom drawn with a notebook label `0` still has code occupancy `1`.

Use `(r,c,k)` for **row, column, layer**. The notebook's `(A,i,layer)`
corresponds to this notation. In `exchangej_.m`, the loop variable `i` is
instead the **layer index**. The explorer uses unpadded core coordinates;
MATLAB array rows and columns are each larger by two after padding.

## 2. How `distrib.m` builds the arrays

The following rules are software index rules read from `distrib.m`, not
fractional crystallographic coordinates. Positions advance by four columns.

| Color | Motif layer | Odd core rows: columns | Even core rows: columns | Code mask |
|---|---:|---|---|---|
| Green | 1 | 1, 5, 9, … | 3, 7, 11, … | `atomtype_layer1gr` |
| Orange | 1 | 3, 7, 11, … | 1, 5, 9, … | `atomtype_layer1or` |
| Blue | 2 | 3, 7, 11, … | 1, 5, 9, … | `atomtype_layer2blue` |
| Purple | 2 | 4, 8, 12, … | 2, 6, 10, … | `atomtype_layer2p` |
| Red | 3 | 1, 5, 9, … | 3, 7, 11, … | `atomtype_layer3red` |
| Black | 3 | 4, 8, 12, … | 2, 6, 10, … | `atomtype_layer3black` |

The motif repeats with `mod(k-1,3)`. Blue reuses orange's lateral mask;
red reuses green's; black reuses purple's. These are position relationships,
not identities between atoms in different layers.

For a small core, the first three rows begin as follows (`.` means empty):

```text
             columns 123456789
Layer 1, row 1        G.O.G.O.G
         row 2        O.G.O.G.O
         row 3        G.O.G.O.G

Layer 2, row 1        ..BP..BP.
         row 2        BP..BP..B
         row 3        ..BP..BP.

Layer 3, row 1        R..KR..KR
         row 2        .KR..KR..
         row 3        R..KR..KR
```

`distrib.m` inserts the core into `3:end-2,3:end-2` of a larger array:

```text
requested W × L core → (W+4) × (L+4) array
core (r,c,k)         → padded (r+2,c+2,k)
two empty rows/columns surround every lateral side
```

There is no corresponding layer padding. Layer-boundary omissions are handled
explicitly in `exchangej_.m`.

## 3. Worked example: green target, J3

A green target at `(r,c,k)` couples to three purple sites in layer `k+1`.
For notebook layer 1, these are the three layer-2 coordinates in the green
`J3` annotation.

| Neighbor | Offset `(Δr,Δc,Δk)` | Shifted x array read at target |
|---|---|---|
| `(r,c-1,k+1)` | `(0,-1,+1)` | `mmx_r1(:,:,i+1)` |
| `(r+1,c+1,k+1)` | `(+1,+1,+1)` | `mmx_u1l1(:,:,i+1)` |
| `(r-1,c+1,k+1)` | `(-1,+1,+1)` | `mmx_d1l1(:,:,i+1)` |

The exact x-component implementation is at `exchangej_.m:38–39`:

```matlab
shift_matx_gr=atomtype_layer1gr.*(mmx_r1(:,:,i+1)+mmx_u1l1(:,:,i+1)+...
    mmx_d1l1(:,:,i+1));
```

The shift gathers the neighbor value into the target's array position. The
mask selects green **targets**, not purple sources. The source's color follows
from its position and layer occupancy. Lines 40–43 repeat the calculation for
the y and z components. Lines 57–59 combine green and orange target sums.

### Why shift signs look reversed

```matlab
shifted = circshift(spins,[dr,dc]);
% At (r,c,k), shifted reads spins(r-dr,c-dc,k), with array wrapping.
```

Thus `circshift(...,[0,+1])` gathers the column-`c-1` neighbor. The suffix
encodes the **array shift**, not the neighbor displacement. Definitions live
in `integrate_llg.m:95–133`.

## 4. Complete neighbor reference

To keep all 24 combinations readable, define four lateral offset triplets:

| Symbol | Three `(Δr,Δc)` offsets |
|---|---|
| P | `(0,+1), (+1,-1), (-1,-1)` |
| M | `(0,-1), (+1,+1), (-1,+1)` |
| P2 | `(0,+2), (+1,0), (-1,0)` |
| M2 | `(0,-2), (+1,0), (-1,0)` |
| O | `(0,0)` — one neighbor |

`P @ +1`, for example, means all three P offsets in layer `k+1`.
Rows separated by a semicolon in a J4 entry are its two three-neighbor groups.
Line numbers refer to `exchangej_.m`; surface branches are listed separately
below. These are the notebook's neighbor coordinates translated into the
code's repeated motif, with the handwritten destination-color annotations.

| Target | Term | Neighbor color and offset group | Code lines |
|---|---|---|---|
| Green | J1 | Red: O @ -1 | 12–14 |
| Green | J2 | Orange: P2 @ 0 | 18–23 |
| Green | J3 | Purple: M @ +1 | 38–43 |
| Green | J4 | Blue: P2 @ +1; Black: M @ -1 | 78–86 |
| Orange | J1 | Blue: O @ +1 | 8–14 |
| Orange | J2 | Green: M2 @ 0 | 25–30 |
| Orange | J3 | Black: P @ -1 | 50–55 |
| Orange | J4 | Purple: P @ +1; Red: M2 @ -1 | 88–96 |
| Blue | J1 | Orange: O @ -1 | 104–106 |
| Blue | J2 | Purple: P @ 0 | 110–115 |
| Blue | J3 | Red: M2 @ +1 | 131–136 |
| Blue | J4 | Black: P @ +1; Green: M2 @ -1 | 151–159 |
| Purple | J1 | Black: O @ +1 | 104–106 |
| Purple | J2 | Blue: M @ 0 | 118–123 |
| Purple | J3 | Green: P @ -1 | 138–143 |
| Purple | J4 | Red: P @ +1; Orange: M @ -1 | 161–169 |
| Red | J1 | Green: O @ +1 | 182–184 |
| Red | J2 | Black: M @ 0 | 195–200 |
| Red | J3 | Blue: P2 @ -1 | 208–213 |
| Red | J4 | Purple: M @ -1; Orange: P2 @ +1 | 258–266 |
| Black | J1 | Purple: O @ -1 | 178–184 |
| Black | J2 | Red: P @ 0 | 188–193 |
| Black | J3 | Orange: M @ +1 | 220–225 |
| Black | J4 | Blue: M @ -1; Green: P @ +1 | 248–256 |

Layer 1's layer-3 neighbors belong to the preceding motif block; layer 3's
layer-1 neighbors belong to the following block. This is not a wrap from the
last physical layer to the first.

The blue J4 y/z temporaries are named `shift_maty_black` and
`shift_matz_black` (154–159). They use the blue mask and are consumed correctly
at 172–173; the variable names do not change the target color.

## 5. Surfaces and lateral edges

For a stack starting at motif layer 1 and ending at motif layer 3:

| Surface target | J1 count | J2 count | J3 count | J4 count | Boundary code |
|---|---:|---:|---:|---:|---|
| First-layer green | 0 | 3 | 3 | 3 | 7–10, 63–69 |
| First-layer orange | 1 | 3 | 0 | 3 | 45–48, 71–76 |
| Last-layer red | 0 | 3 | 3 | 3 | 177–180, 241–246 |
| Last-layer black | 1 | 3 | 0 | 3 | 215–218, 234–239 |
| Interior, any color | 1 | 3 | 3 | 6 | Interior branches above |

These counts assume the target is away from lateral edges. At lateral edges,
some listed positions lie in the empty padding and contribute zero. Although
`circshift` wraps the allocated array, the two-cell border and maximum lateral
offset of two prevent an occupied core site from coupling to the opposite
occupied edge. The explorer fades absent neighbors instead of drawing them
as real atoms. Blue and purple are never the terminating layers in this stack.
Other stack terminations are not supported by this illustration.

## 6. What the comparison establishes

The audit compared the original PDF pages 10–11 with all six color branches,
all four interactions, and all three component sums. A MATLAB R2024b Update 1
CPU diagnostic used a `5 × 9 × 6` core (140 occupied sites), the actual
`distrib.m`, the actual shift definitions in `integrate_llg.m`, and the actual
`exchangej_.m`. Each occupied site was activated separately with x/y/z probe
amplitudes `1,2,-3`, and every output was compared exactly with direct indexing
of an independently transcribed notebook neighbor table. All comparisons
passed. Directed-bond counts were `114,328,276,544` for J1–J4, with zero
reciprocity mismatches in that finite lattice. This was an operator diagnostic,
not a spin trajectory or GPU benchmark.

This establishes code/notebook agreement on the audited lattice, including
both row parities, both surfaces, interior layers, and lateral edges. Physical
shell distances, crystallographic mapping, and parameter provenance remain
unresolved as described in [model methodology](model-methodology.md).
The exchange Hamiltonian convention in the shared
`presentations/equations/sections/hamiltonian.tex` informed the reciprocity
check; it was not used as evidence for the material-specific offsets.

### Historical visualization synchronization check (2026-10-15)

The updated explorer passed 100 browser combinations across site, interaction,
surface, and lateral-edge controls, including every neighbor selection, exact
embedded-source equality, mobile overflow, and absence of JavaScript errors.
An independent CPU basis-probe comparison of the displayed tables with the
actual MATLAB scripts passed all 140 occupied-site probes on a `5 × 9 × 6`
core, for J1–J4 and signed DMI in x/y/z, all six masks, both row parities,
surfaces, and lateral edges. MATLAB was R2024b Update 1
(`24.2.0.2740171`). No trajectories or MAT outputs were generated.
This verifies visualization/code agreement, not the physical model or DMI
antisymmetry. DMI-02 remained unresolved for that pre-PR-18 operator; these
results are historical, not a check of the current DMI implementation.

Commands used (local disposable harnesses, not a maintained test suite):

```bash
node --check docs/structure-exchange.js
node --check docs/unit-cell.js
/home/zhuzibn/.venvs/webpage-capture/bin/python /tmp/check_current_fe2o3_visualizations.py
node /tmp/export_fe2o3_visualization.cjs
matlab -batch "run('C:/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/.review-20261013/verify_visualizations.m')"
git diff --check
```

## 7. 3D notebook schematic

Open [the rotatable 3D view](unit-cell.html) to explore the first structure
figure on notebook page 1. It supports dragging or arrow-key rotation, zoom,
near-top and side views, spin-arrow and cell-outline toggles, and independent
1–4 repetitions along the two basal directions and the vertical direction.
Select a sphere to inspect its illustration label and the prisms sharing it.

The cell-outline toggle shows one dashed hexagonal outline per illustrated
layer, following the visual style of the notebook's first figure. Nearby
staggered marker heights are grouped into six illustrated layers per prism.
Where a layer includes both perimeter column groups, its outline joins their
respective heights rather than drawing two parallel outlines. Perimeter atoms
therefore coincide with guide vertices; the central column remains inside the
hexagon. Guides repeat with the prisms and share coincident edges. This is a
graphical grouping of the schematic markers, not an assignment of physical
crystallographic planes or refined atomic heights.

The author explicitly chose a **notebook-based schematic**, rather than a
crystallographic reconstruction requiring a complete coordinate table.
The drawing uses a regular hexagonal prism, a central column, alternating
perimeter column groups, and four A–D markers per column, inspired by the
notebook. Its vertical offsets, four-site spacings, and aspect ratio are
illustrative graphics choices, not adopted material parameters. Coincident
markers and shared prism edges are deduplicated as the drawing repeats.
The sphere count includes exterior columns and must not be interpreted as
the number of Fe atoms in a bulk unit cell or a modeled surface termination.

Within each column sequence, A–D are ordered from top to bottom to match the
notebook's first figure, with the $c$ axis pointing upward. The dashed guides
use those same marker heights. Repeated prisms retain this label order.

The visual reference is Yang Cheng, Sisheng Yu, Menglin Zhu, Jinwoo Hwang,
and Fengyuan Yang, *Electrical Switching of Tristate Antiferromagnetic Néel
Order in α-Fe₂O₃ Epitaxial Films*, Physical Review Letters **124**, 027202
(2020), [DOI: 10.1103/PhysRevLett.124.027202](https://doi.org/10.1103/PhysRevLett.124.027202),
Fig. 1(a), reproduced in the notebook. The original PDF's caption confirms a
Fe-only schematic with in-plane moments; oxygen is omitted in this view too.
The opposite A,D versus B,C arrow convention is described in Fig. 1 and
Sec. I of [Dannegger and coauthors, Physical Review B **107**, 184426 (2023)](https://doi.org/10.1103/PhysRevB.107.184426).
Both indexed local paper summaries and the relevant original PDF pages were
consulted with the author's permission. Neither supplied a full coordinate
table used by this drawing. No online structure search was performed.

Column-group coloring is an illustration aid, separate from the six
MATLAB site-mask colors. No one-to-one mapping from these drawn markers to
`(r,c,k)` is asserted, and no exchange bonds are inferred from drawn distances.
The MATLAB model, its parameters, and its unresolved physical mapping are
unchanged. `unit-cell.js` is documentation graphics only.

## 8. DMI links in the explorer

Select **DMI** to inspect the J3/J4 neighbors reused by `dmi.m` at merged
PR #18 (`864d6b0`). Blue links show J3; brown links show J4. Atom colors still
identify the six site masks. Labels show the signed energy weight of each
read, not a reconstructed physical DMI vector. The former `nex`/`pre` sums
and unshifted J1 reads are absent from the current code.

Using the offset symbols in section 4, the reused sums are:

| Target | J3 offsets, weight $q_kD_3$ | J4 offsets, weight $-q_kD_4$ |
|---|---|---|
| Green | M @ +1 | P2 @ +1; M @ -1 |
| Orange | P @ -1 | P @ +1; M2 @ -1 |
| Blue | M2 @ +1 | P @ +1; M2 @ -1 |
| Purple | P @ -1 | P @ +1; M @ -1 |
| Red | P2 @ -1 | M @ -1; P2 @ +1 |
| Black | M @ +1 | M @ -1; P @ +1 |

Here $q_k=1-2\operatorname{mod}(k,2)=(-1)^k$. An odd target layer uses
$-D_3$ on J3 and $+D_4$ on J4; an even target layer uses $+D_3$ and $-D_4$.
The sign is determined by the **actual target layer**, not its motif label:
green/orange occupy layers 1 and 4, blue/purple 2 and 5, red/black 3 and 6.
The explorer exposes both occurrences, including the two interior choices
for blue/purple. Increasing array layer index still points opposite to
physical $z$, as documented in `dmi.m`.

An interior target has nine reads (three J3, six J4). DMI on J1/J2 is omitted;
their exchange contributions remain unchanged. The same surface guards and
lateral padding as exchange supply the omissions. Away from lateral edges,
first-layer green and last-layer red retain six DMI neighbors; first-layer
orange and last-layer black retain three. Interior targets retain nine.
The explorer uses a small `9 × 17 × 6` core for teaching the index rules;
`main.m` defaults to `20 × 20 × 6`, padded to `24 × 24 × 6`.

`dmi.m:9–14` constructs energy-weighted sums from the current-state exchange
sums, without a uniform-sublattice approximation:

$$
\mathbf S_{D,i}=q_{k_i}(D_3\mathbf S_{J3,i}-D_4\mathbf S_{J4,i}).
$$

`integrate_llg.m:174–182` masks these sums and assembles the flux density:

$$
\mathbf B_{D,i}=(S_{D,i,y},-S_{D,i,x},0)/\mu_s.
$$

With $A_3,A_4$ denoting the implemented neighbor matrices on occupied sites,
the code comments specify $W=\operatorname{diag}(q_k)(D_3A_3-D_4A_4)$ and
$E_D=-\mathbf x^T W\mathbf y$, where $\mathbf x,\mathbf y$ contain the spin
components. $W$ and $\mathbf S_D$ are in joules, $\mu_s$ in J/T, and
$\mathbf B_D$ in tesla. This transcribes the accepted implementation; it
does not establish its material-specific physical interpretation.

`main.m:43–51` sets `Dsim3=1.6e-23 J`, `Dsim4=1.28e-23 J`, or both zero
when disabled. The explorer shows the enabled defaults and each signed weight.
The source labels the absolute values user-selected and the ratio 1.25 a
Fig. 3(b) estimate. That comment's publication locator and estimate are not
verified here: physical parameter provenance remains **UNRESOLVED**.

Selecting a neighbor highlights its J3/J4 sum in the complete `dmi.m` snapshot
and its actual shift in the relevant `exchangej_.m` branch. An omitted
surface contribution may have no corresponding shift in that branch; its
zero assignment or reduced neighbor group is shown instead. Parameter and
field snapshots are available in the expandable DMI source panel.

**Historical DMI-02:** the imported operator had 544 antisymmetry mismatches
on the 140-site audit. It was superseded by PR #18; retain that result as an
import record, not a description of the current operator. Physical directed
bond vectors, crystal-to-array mapping, shell assignments, and the applicability
of the z-only approximation remain **UNRESOLVED**. Notebook page 8 is still a
historical comparison and supplies no independent evidence for the current
signs or magnitudes. See [model methodology](model-methodology.md#dzyaloshinskii--moriya-interaction).

### Current DMI visualization check (2026-10-22)

The displayed tables were exported from `structure-exchange.js` and compared
against basis probes through the actual `distrib.m`, shift definitions,
`exchangej_.m`, `dmi.m`, and DMI field-assembly lines in `integrate_llg.m`.
MATLAB R2024b Update 1 (`24.2.0.2740171`, PCWIN64) ran these checks on CPU:

| Core | Occupied-site probes | J3 directed entries | J4 directed entries | DMI directed entries | $W+W^T$ nonzeros |
|---|---:|---:|---:|---:|---:|
| 5 × 9 × 6 | 140 | 276 | 544 | 820 | 0 |
| 9 × 17 × 6 | 468 | 1032 | 2048 | 3080 | 0 |

All six masks, both row and layer parities, both terminating surfaces, lateral
edges and corners passed. J1–J4 x/y/z sums matched exactly. The DMI sums and
assembled x/y/z fields matched for enabled defaults, D3 only, D4 only, and
both zero (absolute tolerances $10^{-37}$ J and $10^{-13}$ T). The complete
9 × 17 × 6 displayed occupancy and color masks also matched `distrib.m`.
The zero antisymmetry mismatch is an implementation result on these two
lattices; it does not resolve the physical bond-vector or parameter evidence
gaps or establish energy-derivative, dynamics, or size-convergence validation.

Chromium passed 240 combinations of color, actual layer, interaction, lateral
position, and padding display, including 1,056 neighbor selections. Checks
covered linked source highlights, coefficient signs, active-neighbor colors,
the atlas embed, local images, no JavaScript errors, and no page overflow at
390, 768, and 1440 pixels. Desktop and narrow screenshots were inspected.
The DMI comparison now uses the actual zero-padded `page-08.jpg` filename.
The two complete source snapshots (`dmi.m`, `exchangej_.m`), parameter/field
excerpts, and both numerical defaults matched the checked-in files.

Commands used (local disposable harnesses, not a maintained test suite):

```bash
node --check docs/structure-exchange.js
node /tmp/export_fe2o3_shell_dmi.cjs
matlab -batch "run('C:/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/.review-20261013/verify_shell_dmi_visualization.m')"
/home/zhuzibn/.venvs/webpage-capture/bin/python /tmp/check_fe2o3_shell_dmi.py
git diff --check
```

No MATLAB implementation files, trajectories, or MAT results were modified.
The shared `presentations/equations/spintronic_methods.tex` reference and its
`sections/hamiltonian.tex` (effective-field and pair-field conventions) and
`sections/atomistic_boundaries.tex` (retained-bond antisymmetry and open-edge
omissions) informed the checks. No new publication claims were adopted;
the merged code supplies implementation provenance, not scientific evidence.

## Maintaining this explanation

The MATLAB files remain the executable source of truth. Refresh the explorer's
neighbor tables, coefficients, masks, shift names, boundary omissions, source
snapshots, line ranges, and explanatory text together whenever relevant code
changes. Verify all six site masks and every interaction and boundary option.
`dmiBonds(color,k)` must reuse the J3/J4 offset lists with the actual layer
sign; do not restore independent per-color DMI lists. Refresh `dmiSource`
and `exchangeSource` as JSON-escaped complete MATLAB file text, and
`dmiParameterSource`/`dmiFieldSource` from the displayed source ranges. Check
those strings against the files after normalizing CRLF to LF, and update
the visible line ranges and `dmiDefaults` with any source changes.
The notebook atlas preserves the historical notebook; its embedded code
explorer must describe current code. The 3D notebook schematic remains an
explicit illustration without a verified MATLAB site mapping. Do not infer
physical coordinates or remove evidence gaps from visual agreement alone.
