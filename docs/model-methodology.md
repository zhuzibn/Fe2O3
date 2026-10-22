# Alpha-Fe2O3 model methodology

## Purpose and status

This document describes the behavior currently implemented by the MATLAB files
in this repository. It is intended to make the model inspectable and to expose
where scientific evidence or numerical validation is still missing.

This is **not yet a validated reproduction of a published hematite model**.
Where the implementation, the annotated notebook, and an identified
publication have not been reconciled, the affected item is marked
`UNRESOLVED`. Such an item must not be presented as an established material
property or silently changed in the code.

The MATLAB files were imported from `code_20261013` after the ten-item review
in [../diff.md](../diff.md). The sole initial import correction replaced 18 adjacent-layer
DMI reads from trajectory buffers with reads from the current spin state.
Merged PR #18 (`864d6b0`) corrected ANI-01 and OUT-01, corrected TIME-01
update counts and saved-time labels, and replaced the DMI operator with J3/J4
shell weights that pass the tested antisymmetry checks. The final-time state
is still excluded from the saved trajectory. See
[merged PR verification](#merged-pr-18-verification) for resolution evidence
and [import verification](#import-verification) for the historical defects.

Evidence labels used below are:

- **Implemented**: directly observed in the checked-in MATLAB code.
- **Notebook**: recorded in `Alpha-Fe2O3.pdf`, but not independently verified.
- **Cataloged source**: an applicable peer-reviewed paper is identified in the
  local paper catalog, but its exact supporting passage was not inspected for
  this documentation pass.
- **UNRESOLVED**: evidence, convention, or validation is incomplete.

The source-of-truth hierarchy is:

1. The MATLAB code defines current executable behavior.
2. This document records the intended interpretation and known gaps.
3. [`Alpha-Fe2O3.pdf`](../Alpha-Fe2O3.pdf) preserves the original design notes
   and spatial diagrams.
4. Original peer-reviewed publications must support physical structure,
   material parameters, and comparison claims.

A publication does not by itself authorize changing existing behavior. A
conflict must be reported and resolved through an approved implementation plan.

## Model classification

The intended model is a **Class 4 material-specific crystallographic atomistic
spin model** for $\alpha$-Fe$_2$O$_3$.

- Each occupied array site carries a classical normalized spin
  $\mathbf{s}_i=(s_{ix},s_{iy},s_{iz})$.
- The equation family is fixed-length LLG/LLGS. Every implemented integrator
  explicitly renormalizes the updated spin.
- The magnetic order is represented through the occupied-site pattern and
  alternating initial spin directions, rather than through two coherent
  macrospins.
- The code contains only a classical spin model. Electronic bands, electronic
  altermagnetic splitting, and charge transport are outside its scope.
- The code uses one magnetic-moment magnitude for every occupied site.

The exact crystallographic mapping of every array site, including lattice
vectors and fractional coordinates, remains `UNRESOLVED`. The notebook provides
visual mappings, while `distrib.m` supplies the executable occupancy pattern.

The [3D notebook view](unit-cell.html) is an author-requested schematic of the
page-1 Fe-column figure. Its graphical spacing and repetition controls are not
physical structure inputs. It does not assign the illustrated A–D column
labels to MATLAB masks, infer exchange bonds, or resolve the crystallographic
mapping. The simulation geometry and parameters remain unchanged.

## State and array representation

The spin components are stored in arrays of shape
`natomW x natomL x natomH`. The first two indices are the lateral array indices;
the third is the layer index. `distrib.m` expands the requested lateral
dimensions by four, placing the requested region inside a two-cell-wide empty
border on every side.

An occupancy mask `atomtype_` uses `1` for an occupied magnetic site and `0`
for an empty array entry. Component arrays are multiplied by this mask before
neighbor fields are constructed. The complementary mask

$$
\texttt{ato\_s}=|\texttt{atomtype\_}-1|
$$

sets all three components of an empty entry to one only while the scalar GPU
integrator is evaluated. This avoids division by zero during normalization;
the empty entries are masked out again afterward.

### Three-layer motif and site masks

The layer type repeats according to `mod(i-1,3)`:

| Motif layer | MATLAB case | Site masks | Reused lateral mask |
|---|---:|---|---|
| Layer 1 | `0` | green, orange | none |
| Layer 2 | `1` | blue, purple | blue = layer-1 orange |
| Layer 3 | `2` | red, black | red = layer-1 green; black = layer-2 purple |

The color names are software labels inherited from the diagrams; they are not
chemical species. The explicit parity rules that create these masks are in
`distrib.m`. The notebook pages 7--11 show the intended honeycomb-to-array
encoding and the same color labels.

The corrected site-distribution diagram on notebook page 10 and `distrib.m`
are identical at the array-mask level. A direct MATLAB check with
`natomW=3`, `natomL=13`, and `natomH=3` verified every displayed green,
orange, blue, purple, red, black, and empty position in the unpadded core. In
particular, the Layer-2 middle row ends `BP..B`: column 12 is empty and column
13 is blue. The page annotation states “empty” and “the code correctly
implement this.” This closes the diagram-to-implementation discrepancy only;
it does not resolve the qualifying-evidence gap for the physical
crystal-to-array mapping.

The default `natomH=6` contains two repetitions of the three-layer motif. The
code assumes that adjacent layer indices exist where referenced; arbitrary
values of `natomH`, especially values that truncate the motif differently,
have not been validated.

### Distribution controls

`distrib.m` always generates occupancy and all six masks. The normal path is
deterministic even though `main.m` calls `rng('shuffle')` and prints “use random
atom distribution.” The old `atomtypee.mat` save/load block has been removed.

Despite its legacy name, `load_fixed_atom_distrib=1` now loads the final
`mmx/mmy/mmz` spin frame from `Relax_0.mat` in `systemgeneration.m`.
That file is not supplied. No compatibility check against the generated
dimensions and masks is performed. The default is zero. The remaining
`save_fixed_atom_distrib` flag has no save action.

## Geometry and boundary behavior

The implementation uses array shifts rather than an explicit bond table.
`integrate_llg.m` constructs ten lateral shifted copies of each spin component
using `circshift` with offsets

$$
(0,\pm1),\quad (0,\pm2),\quad (\pm1,0),\quad
(\pm1,\pm1).
$$

Although `circshift` is periodic as an array operation, the two-cell empty
border is masked to zero before fields are calculated. With the currently used
maximum lateral shift of two, the occupied core therefore sees missing/zero
neighbors at its lateral boundary rather than occupied sites from the opposite
side. This is best described as an implicit open-boundary construction using
empty padding, not as a verified periodic bulk model.

The first and last physical layers are also treated as open in the $c$
direction by conditional omission of out-of-range interlayer contributions in
`exchangej_.m` and `dmi.m`.

The following remain `UNRESOLVED`:

- generalization of the finite-lattice boundary audit below to other supported
  sizes and a separate physical justification of surface coordination;
- whether zero-valued padded neighbors are exactly equivalent to omission for
  every future field term;
- the physical mapping from lateral array offsets to crystallographic lattice
  vectors; and
- finite-size convergence and edge sensitivity of spatial results.

## Interaction topology

The notebook states that the intended exchange coordination numbers are
$1,3,3,6$ for $J_1,J_2,J_3,J_4$, respectively. The implementation builds
site- and layer-dependent neighbor sums in `exchangej_.m`:

- $J_1$ is an interlayer contribution with one selected neighbor per interior
  site; a terminating surface can omit it.
- $J_2$ uses three same-layer lateral neighbors.
- $J_3$ uses three neighbors in one adjacent layer, selected by site mask.
- $J_4$ combines two three-neighbor groups in adjacent layers for an interior
  site; a surface can retain only one group.

The detailed offset selection differs among green, orange, blue, purple, red,
and black sites. `exchangej_.m` is currently the executable source of truth for
those selections; notebook pages 10--11 are the visual design record.

The maintained [implementation guide](structure-exchange-guide.md) supplies the complete
24-entry color/interaction mapping and a worked explanation of the reversed
`circshift` signs. The [interactive explorer](structure-exchange.html) connects
the notebook, array masks, neighbor positions, boundary omissions, and exact
current source excerpts. Its DMI view reuses the J3/J4 neighbor tables with
separate shell weights and signs determined by the actual target layer.
Labels show algebraic energy weights, not physical DMI vectors. Both motif
occurrences are selectable so the sign reversal is visible.
Its coordinates are array indices, not reconstructed
crystallographic positions. The author's confirmed notebook `0`/`1` box labels
are distinct from the code's occupied/empty mask values.

A CPU MATLAB R2024b Update 1 basis-probe audit on a $5\times9\times6$ core
(140 occupied sites) compared all three components of all four implemented
neighbor sums with an independently transcribed notebook table. Every output
matched exactly, including both row parities, layer surfaces, interior layers,
and lateral edges. Directed-bond counts were 114, 328, 276, and 544 for
$J_1$--$J_4$, with zero reciprocity mismatches on this finite lattice.
This establishes notebook/code agreement for that case; it is not a
crystallographic validation, energy-derivative test, or dynamics benchmark.

The physical assignment of every offset group to a crystallographic shell,
including distance, multiplicity, bond direction, and source-table row, is
`UNRESOLVED`. It must be verified before the topology is refactored or used as
evidence for a material-specific claim.

## Energy and effective-flux-density convention

The effective quantities named `h...` in the MATLAB code are added directly to
the external input and passed to a gyromagnetic ratio in
$\mathrm{rad\,s^{-1}\,T^{-1}}$. They are therefore interpreted here as
effective flux densities in tesla:

$$
\mathbf{B}_{\mathrm{eff},i}
=-\frac{1}{\mu_s}\frac{\partial \mathcal{H}}{
\partial \mathbf{s}_i}.
$$

This interpretation agrees with the field convention documented in the local
`spintronic_methods.tex` reference. The variable name `Hext` is historical;
its components must be supplied in tesla under the current implementation.

The conservative energy consistent with the implemented exchange and
second-order anisotropy fields is

$$
\mathcal{H}_{\mathrm{ex+2}}
=\frac{1}{2}\sum_{i\ne j}J_{ij}\mathbf{s}_i\cdot\mathbf{s}_j
-\sum_i K_1s_{iz}^2,
$$

when the pair sum contains both directed occurrences of each bond. In this
sign convention, positive $J$ favors antiparallel alignment and negative $J$
favors parallel alignment.

The fourth-order term is no longer active. A consistent complete energy for
the imported sixth-order field remains `UNRESOLVED`; the current code's DMI
energy expression is transcribed below. The expression above is not the
complete Hamiltonian.

### Exchange

For the four implemented neighbor sums $\mathbf{S}_{n,i}$,

$$
\mathbf{B}_{\mathrm{ex},i}
=-\frac{1}{\mu_s}\sum_{n=1}^{4}J_n\mathbf{S}_{n,i}.
$$

This is implemented component by component in `integrate_llg.m`. No extra
factor of $1/2$ appears in the field. Pair reciprocity passed the finite-lattice
notebook comparison described above. Energy/field finite-difference agreement
has not yet been tested.

### Anisotropy

The current anisotropy block is `integrate_llg.m:159–173`. It retains the
imported numerators with the ANI-01 denominator guard added in `3e20e90`.
Define $q=\sin(\arccos s_z)$ and $\bar q=q+\mathbf{1}_{q=0}$, where the
indicator is one only when the computed $q$ is zero. With
$u=s_x/\bar q$, $v=s_y/\bar q$, the sixth-order numerators are

$$
A_x=12h_6q^5uv^2(2v+1)(2v-1)(-3+4v^2),\quad
A_y=-12h_6q^5u^2v(2v+1)(2v-1)(-3+4v^2),
$$

$$
A_z=-6h_6s_zq^4(-1+2u^2)(16u^4-16u^2+1).
$$

The code multiplies these by `mask=~(abs(mz)==1)` and divides by `mus`,
then adds $2K_1s_z/\mu_s$ to z. The guarded division avoids the former +z
`0/0`, and the mask removes the sixth-order contribution at both exact poles.
Thus the pole fields are finite: $(0,0,2K_1s_z/\mu_s)$ for $s_z=\pm1$.
**ANI-01 is resolved** as an arithmetic defect; the off-pole formula is
unchanged. These expressions transcribe the code, not a new energy derivation.

The second-order term uses `Ksim1=-8e-22 J`, and the sixth-order coefficient
is `hk6=1.76e-27 J`. Their material provenance, the sixth-order energy
derivation, and the physical axis mapping remain `UNRESOLVED`.

### Dzyaloshinskii--Moriya interaction

**Implemented at merged PR #18.** `exchangej_.m` first computes the current
J3/J4 spin-neighbor sums. `dmi.m:9–14` reuses them with weights in joules:

$$
q_k=1-2\operatorname{mod}(k,2)=(-1)^k,\qquad
\mathbf S_{D,i}=q_{k_i}(D_3\mathbf S_{J3,i}-D_4\mathbf S_{J4,i}).
$$

The actual layer index $k$ determines the sign, not the three-layer motif
label or target color. The same color at $k+3$ has the opposite sign.
Array layer index increases opposite to physical $z$, as recorded in the
code comments. Surface and lateral-edge omissions are inherited from the
exchange sums. There are three J3 and six J4 neighbors per interior target;
there is no J1/J2 DMI or unshifted adjacent-layer read. No uniform-sublattice
approximation replaces individual neighboring spins.

`integrate_llg.m:174–182` masks the sums and forms the field in tesla:

$$
\mathbf B_{D,i}=(S_{D,i,y},-S_{D,i,x},0)/\mu_s.
$$

The code comments specify

$$
W=\operatorname{diag}(q_k)(D_3A_3-D_4A_4),\qquad
E_D=-\mathbf x^T W\mathbf y,
$$

where $A_3,A_4$ are the implemented neighbor matrices on occupied sites and
$\mathbf x,\mathbf y$ are the spin-component vectors. This is a transcription
of the accepted code's energy/field convention, not a new physical bond-vector
assignment or a literature validation.

`main.m:43–51` sets `Dsim3=1.6e-23 J` and `Dsim4=1.28e-23 J` when enabled,
and both zero when disabled. Its comments describe the absolute values as
user-selected and their ratio 1.25 as a Fig. 3(b) estimate. The exact
publication locator and estimate have not been verified in this update;
physical parameter provenance remains `UNRESOLVED`.

**DMI-02 resolved for the tested operators:** the imported operator had 934
directed nonzero coefficients and 544 antisymmetry mismatches on the 140-site lattice. It is
superseded by the merged operator above, which has zero antisymmetry mismatches
on the checked lattices. The historical import check below
is retained as evidence of that earlier code. Physical shell assignments,
crystal-to-array mapping, directed bond vectors, and the applicability of
the z-only approximation remain `UNRESOLVED`.

**Visualization/code verification (2026-10-22):** CPU basis probes through
the actual scripts matched the displayed sums and field assembly on
5 × 9 × 6 (140 sites) and 9 × 17 × 6 (468 sites) cores, including all masks,
layer parities, surfaces, and lateral edges. The current operator has 820
and 3080 directed entries respectively, with zero $W+W^T$ nonzeros on both
lattices. Enabled, isolated-shell, and disabled weights passed. See the
[visualization check record](structure-exchange-guide.md#current-dmi-visualization-check-2026-10-22)
for commands, tolerances, and browser checks. The separate
[PR review](pr18-time-01-change-report.html#evidence) also records a
normalized-coefficient energy finite-difference check, summarized below.
These implementation checks do not establish material-specific physical validity.

### Dipole approximation

The implemented dipole field retains three same-layer neighbors per occupied
site at a fixed distance $a=0.5067\ \mathrm{nm}$. For each retained neighbor
$j$, `dipole_.m` evaluates

$$
\mathbf{B}_{ij}
=\frac{\mu_0}{4\pi}
\left[
\frac{3\mathbf{r}_{ij}
(\boldsymbol{\mu}_j\cdot\mathbf{r}_{ij})}{|\mathbf{r}_{ij}|^5}
-\frac{\boldsymbol{\mu}_j}{|\mathbf{r}_{ij}|^3}
\right],
\qquad
\boldsymbol{\mu}_j=\mu_s\mathbf{s}_j.
$$

This is a nearest-neighbor in-plane truncation, not a full long-range dipolar
sum. Interlayer dipole terms and more distant lateral sites are omitted. The
field is recomputed every 50 integration steps and held at its previous value
between updates when enabled. The imported default is `dipolemode=0`.

### External field and total field

The conservative field passed to the integrator is

$$
\mathbf{B}_{\mathrm{eff}}
=\mathbf{B}_{\mathrm{ex}}
+\mathbf{B}_{\mathrm{ani}}
+\mathbf{B}_{\mathrm{DMI}}
+\mathbf{B}_{\mathrm{ext}}
+\mathbf{B}_{\mathrm{dip}}.
$$

Despite its name, `Hext` is interpreted as $\mathbf{B}_{\mathrm{ext}}$ in
tesla because no multiplication by $\mu_0$ occurs before it is added.

## LLG/LLGS dynamics and torques

With torque amplitudes set to zero, the scalar kernels implement the
Landau--Lifshitz form equivalent to fixed-length Gilbert dynamics:

$$
\frac{d\mathbf{s}}{dt}
=-\frac{\gamma}{1+\alpha^2}
\left[
\mathbf{s}\times\mathbf{B}_{\mathrm{eff}}
+\alpha\mathbf{s}\times
(\mathbf{s}\times\mathbf{B}_{\mathrm{eff}})
\right].
$$

`atomgpurk4.m` additionally applies two pairs of polarization torques. One pair
is labeled SOT and uses `psjSHE`; the other is labeled STT and uses `psjSTT`.
Both use the same damping-like and field-like vector structures. In the RK4
kernel, a damping-like amplitude $B_D$ and field-like amplitude $B_F$ contribute

$$
\frac{\gamma}{1+\alpha^2}
\left[
-B_D\,\mathbf{s}\times(\mathbf{s}\times\boldsymbol{\sigma})
+\alpha B_D\,\mathbf{s}\times\boldsymbol{\sigma}
+\alpha B_F\,\mathbf{s}\times
(\mathbf{s}\times\boldsymbol{\sigma})
+B_F\,\mathbf{s}\times\boldsymbol{\sigma}
\right].
$$

The two labels do not establish two different microscopic mechanisms in the
current kernel. In particular, the STT-labeled term is a polarization torque,
not a spatial-gradient adiabatic/nonadiabatic transport term.

The code forms the amplitudes as

$$
B_{D,\mathrm{SOT}}
=\frac{\hbar\theta_{\mathrm{SH}}j_{c,\mathrm{SOT}}}
{2M_st_z},
\qquad
B_{D,\mathrm{STT}}
=\frac{\hbar\eta_{\mathrm{STT}}j_{c,\mathrm{STT}}}
{2M_st_z},
$$

where `hbar` is stored in $\mathrm{eV\,s}$ so that the elementary-charge
conversion implicit in the usual SI expression is canceled numerically. This
unit interpretation should be checked explicitly before nonzero-current
results are used.

The default SOT current is `6e12 A/m^2`, with both components enabled and
polarization approximately `[-0.866,-0.5,0]`. With
`ms=24*mus/(sqrt(3)*a^2*c)` and `tz=c*natomH/6`, the amplitudes are
`BDSOT=0.18648 T` and `BFSOT=5*BDSOT=0.93239 T`. STT is disabled and its
current is zero. The inherited `BFSTT=chi*BDSTT` uses the SOT ratio rather
than `chiSTT`; this nondefault-path limitation is retained, not repaired.

## Numerical integration

Three update kernels are selectable:

| `rk4` | Kernel | Implemented behavior |
|---:|---|---|
| `1` | `atomgpurk4.m` | Four-stage RK4 update; includes the SOT/STT-labeled torque arguments |
| `0` | `atomgpu.m` | Two-stage Heun update; no torque arguments |
| `2` | `atomgpupc4.m` | Four-step Adams--Bashforth predictor and Adams--Moulton corrector; no torque arguments |

The predictor-corrector selection attempts to bootstrap its first three steps
with `atomgpurk4.m`, but the call in `integrate_llg.m` supplies only the
field-only argument list while `atomgpurk4.m` requires the longer torque-aware
signature. The selected `rk4=2` path is therefore expected to fail before it
reaches `atomgpupc4.m`; it has not been executed in this documentation pass.

Each update kernel normalizes the final spin to unit length. The effective field
is evaluated once before the scalar kernel call and is held fixed during the
RK4/Heun substages or predictor-corrector evaluations. Consequently, `rk4=1`
is not a fully coupled RK4 evaluation of a state-dependent effective field.

The current guard rejects non-RK4 integration only when an SOT switch is set.
It does not guard nonzero STT. Choosing Heun or predictor-corrector with a
nonzero STT amplitude therefore omits that torque silently. This is a known
implementation limitation.

The trajectory is processed in GPU chunks. Each kernel result is explicitly
gathered into CPU trajectory buffers. Every `savetstep`-th in-chunk state is
copied to `mmx_`, `mmy_`, and `mmz_`; after the loop these are masked and
renamed `mmx`, `mmy`, and `mmz`.

**TIME-01 timing corrected (`ab7d82e`):** `while ct1<gpusteps+1` performs
`gpusteps` updates per chunk. Frame 1 holds the chunk's initial state, and
the last update grows the spin buffer to slot `gpusteps+1`, which supplies
the next chunk's initial state. At the defaults, 150 × 10,000 updates at
`1e-16 s` evolve for 150 ps. Both spin and DMI histories use
`1:savetstep:gpusteps`, and `main.m:104` labels the saved frames with
`t=(0:gpurun_number*final_m_savestep-1)*savetstep*tstep`. For the checked
integral chunk/step ratios, labels match elapsed time and start at zero.

**Final endpoint still unsaved:** sampling includes each chunk's start and
excludes its end. The default 300 frames therefore cover
`0, 0.5, ..., 149.5 ps`. The computed 150 ps endpoint is cleared before
`final.mat` is written; a restart from the last saved frame uses 149.5 ps.
Whether to save the endpoint with a matching timestamp and DMI field remains
an open output decision. Non-integral `gpusave/tstep` behavior and the memory
cost of dynamically growing the extra slot have not been validated.

**OUT-01 resolved (`71db4fd`):** fields evaluated from state `ct1` are now
recorded at `ct1`, while its updated spin is stored at `ct1+1`. Field slot 1
is populated, and every saved field is paired with its evaluation state.
The small GPU checks recomputed each DMI field from its saved spin frame
with zero difference. This recording fix is separate from the initial
import's DMI current-state correction and the later DMI operator change.

## Initial conditions

The default has no domain wall. All occupied sites use polar angle
$90^\circ$; odd layers use azimuth $0^\circ$ and even layers use
$180^\circ$, giving alternating +x/-x directions up to trigonometric
roundoff. Empty entries are zero. The old `dwcalc` branch is removed.
The optional restart uses the last saved spin frame from `Relax_0.mat`.

The correspondence between these array angles and experimentally relevant
hematite magnetic domains is `UNRESOLVED`.

## Checked-in parameter set

The following values are transcribed from `main.m`. Values in meV are derived
here using the exact SI electron volt and are not separate inputs.

| Quantity | MATLAB value | Derived/readable value | Evidence status |
|---|---:|---:|---|
| Lattice scale `c` | $1.3882\times10^{-9}\ \mathrm{m}$ | $1.3882\ \mathrm{nm}$ | Imported; exact source locator unresolved |
| Lattice scale `a` | $5.067\times10^{-10}\ \mathrm{m}$ | $0.5067\ \mathrm{nm}$ | Imported; exact source locator unresolved |
| Effective thickness `tz` | `natomH*c/6` | $1.3882\ \mathrm{nm}$ for `natomH=6` | Imported; physical normalization unresolved |
| $J_1$ | $1.3740\times10^{-21}\ \mathrm{J}$ | $8.5758\ \mathrm{meV}$ | Implemented; provenance conflict |
| $J_2$ | $-1.1696\times10^{-21}\ \mathrm{J}$ | $-7.3001\ \mathrm{meV}$ | Implemented; provenance conflict |
| $J_3$ | $4.0413\times10^{-21}\ \mathrm{J}$ | $25.2238\ \mathrm{meV}$ | Implemented; provenance conflict |
| $J_4$ | $2.8041\times10^{-21}\ \mathrm{J}$ | $17.5018\ \mathrm{meV}$ | Implemented; provenance conflict |
| $K_1$ | $-8\times10^{-22}\ \mathrm{J}$ | $-4.9932\ \mathrm{meV}$ | Imported; source locator unresolved |
| `hk6` | $1.76\times10^{-27}\ \mathrm{J}$ | $0.000010985\ \mathrm{meV}$ | Imported; energy mapping and source locator unresolved |
| $D_3$ | $1.6\times10^{-23}\ \mathrm{J}$ | `Dsim3`; J3 neighbor weights | PR #18; absolute value user-selected according to code; physical provenance unresolved |
| $D_4$ | $1.28\times10^{-23}\ \mathrm{J}$ | `Dsim4`; $D_3/D_4=1.25$ | PR #18; Fig. 3(b) estimate cited in code not independently verified |
| $g$ | $2$ | $\gamma\approx1.7590\times10^{11}\ \mathrm{rad\,s^{-1}\,T^{-1}}$ | Implemented |
| $\mu_s$ | $4.2313\mu_B$ | $3.9241\times10^{-23}\ \mathrm{J\,T^{-1}}$ | Implemented; source locator unresolved |
| Magnetization density `ms` | `24*mus/(sqrt(3)*a^2*c)` | $1.5256\times10^6\ \mathrm{A\,m^{-1}}$ | Imported; physical normalization unresolved |
| Dipole distance | $0.5067\ \mathrm{nm}$ | three retained in-plane neighbors | Implemented; source locator unresolved |
| Gilbert damping $\alpha$ | $0.001$ | dimensionless | Implemented; provenance unresolved |
| Time step | $10^{-16}\ \mathrm{s}$ | $0.1\ \mathrm{fs}$ | Implemented; convergence unresolved |

The four exchange values convert to the LDA+$U$ numbers copied in the left-hand
table on notebook page 3, while the `main.m` comment cites the 2023 hematite
paper. This is a provenance conflict, not proof that either source is wrong.
The exact source table and the mapping between its shell labels and
$J_1$--$J_4$ must be resolved before the comment or values are changed.

## Claim-to-implementation and evidence map

| Physical detail | Current implementation | Location | Scientific evidence status |
|---|---|---|---|
| Hematite crystal/site motif | Three-layer repeating occupancy and six color-labeled site masks | `distrib.m`; notebook pp. 1, 7--11 | Page-10 diagram/implementation agreement verified for a $3\times13\times3$ core; `UNRESOLVED`: no lattice-vector/basis table with exact paper locator |
| Exchange shells | Site-dependent $J_1$--$J_4$ offset sums; intended bulk multiplicities $1,3,3,6$ | `exchangej_.m`; notebook pp. 2--4, 9--11 | Notebook/code agreement and reciprocity passed the finite-lattice audit above; `UNRESOLVED`: physical shell distances and source mapping |
| Exchange values | Four values listed above | `main.m`; notebook pp. 2--3 | `UNRESOLVED`: code comment and copied table provenance conflict |
| Anisotropy | Second-order axial term plus guarded sixth-order field | `integrate_llg.m:159–173` | ANI-01 arithmetic defect resolved; `UNRESOLVED`: sixth-order energy, axis and parameter provenance |
| DMI | $q_k(D_3S_{J3}-D_4S_{J4})$ sums; code energy $-\mathbf x^T W\mathbf y$ | `dmi.m:1–15`; `main.m:43–51`; `integrate_llg.m:174–182` | DMI-02 antisymmetry defect resolved on tested lattices; `UNRESOLVED`: physical bonds, shell mapping, and parameter provenance |
| Dipole | Three same-layer neighbors at $0.5067\ \mathrm{nm}$, refreshed every 50 steps | `dipole_.m`; notebook p. 6 | `UNRESOLVED`: truncation justification and source locator missing |
| LLG/LLGS | Fixed-length LL form in tesla convention | GPU kernels; local equation reference | Focused RK4 precession, damping, and torque checks passed; general convergence unresolved |
| Boundary | Empty two-cell lateral padding; omitted out-of-range layer bonds | `distrib.m`, `exchangej_.m`, `dmi.m` | Finite-lattice notebook/code agreement audited above; broader coordination and finite-size tests missing |
| Initial state | Uniform within each layer, odd/even +x/-x | `systemgeneration.m` | Imported initial condition; physical domain mapping unresolved |

## Validation requirements

No maintained automated test suite or accepted production benchmark is
provided. Focused merged-PR and historical import checks are recorded below.
An earlier manual structure check passed: the corrected page-10 diagram agrees with
all six `distrib.m` site masks for a $3\times13\times3$ unpadded core. Before
using the model for quantitative claims, add or perform the remaining checks
in order:

1. **Structure:** extend the completed page-10 mask check to occupied-site
   counts, three-layer repetition, padding, and additional deterministic array
   sizes.
2. **Topology:** enumerate the actual bond table implied by every offset branch;
   verify reciprocity, interior multiplicities $1,3,3,6$, and surface
   coordination.
3. **Static fields:** compare exchange, anisotropy, DMI, dipole, and external
   contributions independently against hand calculations. Use finite
   differences of a declared energy where applicable.
4. **Dynamics:** verify precession direction, damping direction, torque signs,
   unit-spin preservation, and finite values on one- and few-spin problems.
5. **Integrator:** recompute state-dependent fields at substages or explicitly
   validate the present fixed-field approximation; then perform timestep
   convergence for each supported solver.
6. **Boundaries:** compare the implicit padded-open behavior with an explicit
   neighbor table and test edge coordination and DMI truncation.
7. **Simulation outputs:** verify chunk-boundary time alignment, sampling,
   spatial observables, and non-destructive output naming.
8. **Scientific comparison:** establish size convergence and compare an
   appropriate observable with a source-supported analytical, experimental, or
   published numerical reference.

Passing one stage does not validate later stages. In particular, an apparently
reasonable trajectory cannot validate the crystallographic mapping or bond
table.

## Known limitations and open decisions

- No thermal field is implemented; unused `thermalenable` and `T` settings
  were removed.
- The default distribution is deterministic; restart spin arrays are not
  checked for compatibility with generated dimensions and site masks.
- No explicit lattice-vector, coordinate, or bond-table data structure exists.
- The exchange provenance and shell mapping are unresolved.
- Physical DMI bond orientation, shell mapping, parameter provenance, and
  material-specific energy interpretation remain unresolved. DMI-02 describes
  the superseded import operator, not the current J3/J4 implementation.
- `UNRESOLVED` — notebook symmetry reasoning: $C_3$ symmetry cancels the
  in-plane components in the symmetry-related sum; it does not force every
  individual bond's DMI vector to lie along $z$, as Dannegger et al., Fig. 4,
  shows ([author preprint, Sec. III D and Fig. 4](https://arxiv.org/pdf/2301.13094)).
  The individual bond vectors and the applicability of the implemented
  $z$-only DMI approximation still need to be established.
- Dipole coupling is truncated to three in-plane neighbors and updated only
  every 50 steps.
- The field is held fixed within each nominal RK4/Heun/PC step.
- The predictor-corrector bootstrap call does not match the current RK4
  function signature and is not presently a verified runnable path.
- Heun and predictor-corrector omit torque arguments; the runtime guard does
  not protect against nonzero STT in those modes.
- The final-time state is computed but not saved: default output ends at
  149.5 ps although evolution reaches 150 ps. ANI-01, OUT-01, and the TIME-01
  update-count/time-label defects are corrected; endpoint saving remains open.
- `final.mat` overwrites an existing output and stores the remaining workspace.
- No accepted production baseline or scientific convergence validation is
  provided by the import checks.

## Merged PR #18 verification

Documentation reconciliation on 2026-10-22 checked the merged source at
`864d6b0` against the reviewed head `77dece1`:

```bash
git diff --exit-code 77dece1 HEAD -- '*.m'
```

The comparison passed with no MATLAB differences. The
[PR #18 report](pr18-time-01-change-report.html#evidence) retains its
pre-merge metadata, exact review command, harness, and full log. Its recorded
checks therefore apply to this merged source; they were not rerun in full
for this documentation update. They used MATLAB R2024b Update 1
(`24.2.0.2740171`, PCWIN64) and an NVIDIA GeForce RTX 3060 Laptop GPU
(compute capability 8.6, 6,441,926,656 bytes reported memory).

| Issue/check | Recorded result and scope |
|---|---|
| ANI-01 | Old +z NaNs reproduced; merged ±z, +x, and empty-entry fields finite. All 50,000 seeded off-pole unit-vector fields matched the old formula exactly. A small +z-initialized GPU run remained finite. |
| OUT-01 | Every saved DMI field matched a fresh evaluation from its saved spin, with maximum difference 0, for strides 1, 2, and 8. |
| TIME-01 | Three chunks of eight updates saved 24, 12, and 3 frames at strides 1, 2, and 8, with correct elapsed-time labels. Three × eight versus one × 24 updates at stride 2 gave identical saved spins. |
| DMI-02 | Basis probes on 5 × 9 × 6, 6 × 10 × 9, and 4 × 8 × 3 cores gave 820, 1744, and 212 directed entries, respectively, with zero antisymmetry mismatches. A 140-site finite-difference check of the code energy, using normalized coefficients and perturbation $10^{-6}$, had maximum absolute error $3.98\times10^{-9}$ against a $10^{-8}$ tolerance. |
| Spin invariants | Small GPU cases had finite outputs, zero empty sites, and maximum occupied-spin norm error $1.11\times10^{-16}$. |

The GPU cases used a 5 × 9 × 6 core (140 occupied sites), default RK4 and
torques, `tstep=1e-16 s`, and disabled dipole coupling. The separate
[current DMI visualization check](structure-exchange-guide.md#current-dmi-visualization-check-2026-10-22)
also covers all masks, layer parities, surfaces, edges, isolated D3/D4
weights, and disabled DMI on 140- and 468-site lattices.

This documentation update additionally executed the exact current anisotropy
block on CPU via `matlab -batch`, with `(mx,my,mz)` equal to +z, -z, +x,
and an empty entry, `Ksim1=-8e-22`, `hk6=1.76e-27`, and unit moment.
It returned `(0,0,-1.6e-21)`, `(0,0,1.6e-21)`, and two zero vectors,
respectively, matching the expected arithmetic. MATLAB reported the same
version and platform as above. No trajectory or MAT output was produced.

Exact CPU command executed from the repository root:

```bash
matlab -batch "s=fileread('integrate_llg.m'); first=strfind(s,'%% calcuation for anisotropy field'); last=strfind(s,'%% calcuation for DMI field'); mmxtmp=[0,0,1,0]; mmytmp=[0,0,0,0]; mmztmp=[1,-1,0,0]; hk6=1.76e-27; Ksim1=-8e-22; muigpu=1; eval(s(first(1):last(1)-1)); actual=[hani_x;hani_y;hani_z]'; assert(isequal(actual,[0,0,-1.6e-21;0,0,1.6e-21;0,0,0;0,0,0])); fprintf('PASS current ANI-01: +z, -z, +x and empty entries are finite and match expected fields with unit moment.\n'); fprintf('MATLAB %s; %s; CPU only; no trajectory or MAT output.\n',version,computer);"
```

The default 150 ps duration is a source-derived update count, not a full-run
benchmark. Endpoint saving, non-integral chunk/step ratios, alternative
solver paths, production memory/runtime cost, timestep/size convergence,
and material-specific physical validation remain outside these checks.
The local Hamiltonian, LLGS, and atomistic-boundary equation references
listed below informed the convention and validation scope; no new
publication-derived material claim was adopted.

## Import verification

This section records the 2026-10-15 import, before the fixes above. Its
reproduced defects and source-equality checks describe that historical version.

Focused checks on 2026-10-15 used MATLAB R2024b Update 1 (24.2.0.2740171,
PCWIN64) and an NVIDIA GeForce RTX 3060 Laptop GPU, compute capability 8.6,
6,441,926,656 bytes reported total memory. They establish implementation
properties for these cases, not material-model validity.

| Check | Result |
|---|---|
| Imported source text | All 13 files match `code_20261013` after line-ending normalization except the 18 approved DMI substitutions |
| Structure and initialization | 5 x 9 x 6 core, 9 x 13 x 6 padded, 140 occupied sites; all six masks match pre-import code; occupied odd/even layers initialize along +x/-x |
| Shift and exchange/dipole equivalence | All 30 shifts retain the same offsets; complete exchange/dipole text matches after verified renames |
| DMI current-state correction | All six masks, surfaces, and edges agree with the independently calculated change from stale source reads; changing history alone leaves the corrected result exactly unchanged |
| DMI signed operator | 934 directed nonzeros, zero diagonal; **544 antisymmetry mismatches**, recorded as DMI-02 |
| RK4 kernel | Scalar kernel agrees with independent vector-form fixed-field RK4 for zero torque, isolated SOT/STT DL/FL, and combined channels; maximum component error 0 in this run; analytical precession and damping direction passed |
| Accepted source defects | +z anisotropy NaN reproduced; source DMI-history index trace and nominal/actual update-count mismatch preserved |
| Small GPU/output smoke | 140 occupied sites, 2 chunks of 7 updates, 8 saved frames; finite spins, zero empty sites, maximum occupied-spin norm error $1.11\times10^{-16}$; expected cleanup and unfilled DMI frames preserved |

In that imported version, the DMI coefficient pair between padded sites `(4,3,2)` and
`(3,3,1)` has $C_{ij}=C_{ji}=+1$, rather than opposite signs. This finding
belongs to the superseded operator; it must not be applied to the merged
J3/J4 operator without a new check.

The disposable local harness and old MATLAB snapshots are in ignored
`.review-20261013/`. Commands from the WSL workspace were:

```bash
matlab -batch "run('C:/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/.review-20261013/verify_structure.m')"
matlab -batch "run('C:/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/.review-20261013/verify_import.m')"
matlab -batch "run('C:/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/.review-20261013/verify_gpu.m')"
```

The GPU run used a scratch copy of `main.m` with core dimensions 5 x 9 x 6,
`gpusave=8e-16 s`, `gpurun_number=2`, and `savetstep=2`; other source
defaults, including `tstep=1e-16 s`, were retained. Its result is isolated at
`.review-20261013/smoke/final.mat`. Production defaults and existing result
files were not changed by the smoke run. No full default simulation,
timestep/size convergence study, or published-physics benchmark was run.

## References and local records

1. `Alpha-Fe2O3.pdf`, repository root, annotated methodology notebook, especially
   pp. 1--12. This is a design record, not independent scientific evidence.
2. *Magnetic properties of hematite revealed by an ab initio parameterized spin
   model*, Phys. Rev. B **107**, 184426 (2023),
   DOI: [10.1103/PhysRevB.107.184426](https://doi.org/10.1103/PhysRevB.107.184426).
   Local catalog ID:
   `AFM/2023-Magnetic properties of hematite revealed by an ab initio-PRB-Tobias Dannegger`.
   The exact passages supporting the implemented construction and parameters
   were not inspected in this documentation pass.
3. Yang Cheng and coauthors, *Electrical Switching of Tristate
   Antiferromagnetic Néel Order in $\alpha$-Fe$_2$O$_3$ Epitaxial Films*, Phys.
   Rev. Lett. **124**, 027202 (2020),
   DOI: [10.1103/PhysRevLett.124.027202](https://doi.org/10.1103/PhysRevLett.124.027202).
   Local catalog ID:
   `AFM/2020-Electrical Switching of Tristate Antiferromagnetic Néel Order in-PRL-Yang Cheng`.
4. V. V. Mazurenko and coauthors, Phys. Rev. B paper cited incompletely on
   notebook page 12, where the year is transcribed as 2005. The publication
   year, full title, DOI, exact table, and locator are
   `REFERENCE MISSING`.
5. Local equation reference:
   `/mnt/c/Users/zzf-m/OneDrive/presentations/equations/spintronic_methods.tex`,
   particularly `sections/hamiltonian.tex`, `sections/llg_equation.tex`, and
   `sections/atomistic_boundaries.tex`. The applicable equation labels are
   `eq:atomistic-hamiltonian`, `eq:atomistic-effective-field`,
   `eq:atomistic-pair-fields`, `eq:llg-gilbert`, and `eq:ll-with-flt`. These
   files define shared notation and conventions; they are not
   material-specific evidence for hematite.
