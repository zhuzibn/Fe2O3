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

`load_fixed_atom_distrib` and `save_fixed_atom_distrib` read or write
`atomtypee.mat`. The normal path is deterministic; no random occupation is
generated even though `main.m` calls `rng('shuffle')` and prints “use random
atom distribution.”

Loading a saved `atomtype_` happens after the color/site masks and `ato_s` have
been built. Neither the masks nor `ato_s` are regenerated after loading, and
compatibility with the loaded occupancy array is not checked. This is an
`UNRESOLVED` consistency risk.

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
- finite-size convergence and edge sensitivity of domain-wall results.

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

The [implementation guide](structure-exchange-guide.md) supplies the complete
24-entry color/interaction mapping and a worked explanation of the reversed
`circshift` signs. The [interactive explorer](structure-exchange.html) connects
the notebook, array masks, neighbor positions, boundary omissions, and exact
source excerpts. Its coordinates are array indices, not reconstructed
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
anisotropy fields is

$$
\mathcal{H}_{\mathrm{ex+ani}}
=\frac{1}{2}\sum_{i\ne j}J_{ij}\mathbf{s}_i\cdot\mathbf{s}_j
-\sum_i\left(K_1s_{iz}^2+K_2s_{iz}^4\right),
$$

when the pair sum contains both directed occurrences of each bond. In this
sign convention, positive $J$ favors antiparallel alignment and negative $J$
favors parallel alignment.

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

The implemented anisotropy is aligned with array $z$:

$$
\mathbf{B}_{\mathrm{ani},i}
=\frac{1}{\mu_s}
\left(0,0,2K_1s_{iz}+4K_2s_{iz}^3\right).
$$

It corresponds to the on-site energy
$-K_1s_{iz}^2-K_2s_{iz}^4$. The relationship between array $z$ and the exact
crystallographic axis needs an explicit source locator.

### Dzyaloshinskii--Moriya interaction

`dmi.m` constructs two site-dependent interlayer sums named `nex` and `pre`.
The file warns that increasing the MATLAB layer index corresponds to decreasing
physical $z$, and that the selected neighbors mix the exchange $J_3$ and $J_4$
topologies.

Writing the two constructed sums as
$\mathbf{S}_{\mathrm{nex},i}$ and $\mathbf{S}_{\mathrm{pre},i}$, the implemented
field is

$$
\mathbf{B}_{\mathrm{DMI},i}
=\frac{D}{\mu_s}
\begin{bmatrix}
S_{\mathrm{nex},y}-S_{\mathrm{pre},y}\\
-S_{\mathrm{nex},x}+S_{\mathrm{pre},x}\\
0
\end{bmatrix}.
$$

The notebook derives a related component form on page 5. The complete directed
bond list, $\mathbf{D}_{ij}$ orientation, antisymmetry check, and energy
derivative remain `UNRESOLVED`. Do not infer them only from the color diagram.

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
between updates.

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

The trajectory is processed in GPU chunks. Every `savetstep`-th in-chunk state
is gathered to CPU arrays `mmx_`, `mmy_`, and `mmz_`; after the loop, these are
masked and renamed `mmx`, `mmy`, and `mmz`. Exact saved-time alignment at chunk
boundaries has not been verified.

## Initial conditions

When `dwcalc=1`, the lattice is divided at the midpoint of the second lateral
index. For odd layers, the two halves use in-plane azimuths $185^\circ$ and
$5^\circ$ with polar angle $90^\circ$. For even layers, the azimuths are
$5^\circ$ and $185^\circ$ with polar angle $85^\circ$. This creates the
implemented two-domain initial state with a small even-layer out-of-plane
component.

When `dwcalc=0`, all occupied sites use polar angle $90^\circ$; odd and even
layers use azimuths $210^\circ$ and $30^\circ$, respectively.

The correspondence between these array angles and the experimentally relevant
hematite magnetic domains is `UNRESOLVED`.

## Checked-in parameter set

The following values are transcribed from `main.m`. Values in meV are derived
here using the exact SI electron volt and are not separate inputs.

| Quantity | MATLAB value | Derived/readable value | Evidence status |
|---|---:|---:|---|
| Interlayer scale `c` | $1.3882\times10^{-9}\ \mathrm{m}$ | $1.3882\ \mathrm{nm}$ | Implemented; comment cites Ref. [2], exact locator unresolved |
| Effective thickness `tz` | `natomH*c` | $8.3292\ \mathrm{nm}$ for `natomH=6` | Implemented modeling choice |
| $J_1$ | $1.3740\times10^{-21}\ \mathrm{J}$ | $8.5758\ \mathrm{meV}$ | Implemented; provenance conflict |
| $J_2$ | $-1.1696\times10^{-21}\ \mathrm{J}$ | $-7.3001\ \mathrm{meV}$ | Implemented; provenance conflict |
| $J_3$ | $4.0413\times10^{-21}\ \mathrm{J}$ | $25.2238\ \mathrm{meV}$ | Implemented; provenance conflict |
| $J_4$ | $2.8041\times10^{-21}\ \mathrm{J}$ | $17.5018\ \mathrm{meV}$ | Implemented; provenance conflict |
| $K_1$ | $1.8073\times10^{-25}\ \mathrm{J}$ | $0.001128\ \mathrm{meV}$ | Implemented; source locator unresolved |
| $K_2$ | $1.7624\times10^{-27}\ \mathrm{J}$ | $0.0000110\ \mathrm{meV}$ | Implemented; source locator unresolved |
| $D$ | $-4.8065\times10^{-23}\ \mathrm{J}$ | $-0.3000\ \mathrm{meV}$ | Implemented; sign/topology source unresolved |
| $g$ | $2$ | $\gamma\approx1.7590\times10^{11}\ \mathrm{rad\,s^{-1}\,T^{-1}}$ | Implemented |
| $\mu_s$ | $4.2313\mu_B$ | $3.9241\times10^{-23}\ \mathrm{J\,T^{-1}}$ | Implemented; source locator unresolved |
| Volume length `d` | $0.35\ \mathrm{nm}$ | $M_s=\mu_s/d^3\approx9.1524\times10^5\ \mathrm{A\,m^{-1}}$ | Implemented modeling choice; physical meaning unresolved |
| Dipole distance | $0.5067\ \mathrm{nm}$ | three retained in-plane neighbors | Implemented; source locator unresolved |
| Gilbert damping $\alpha$ | $0.001$ | dimensionless | Implemented; provenance unresolved |
| Temperature `T` | $100\ \mathrm{K}$ | unused | Implemented variable only |
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
| Anisotropy | $-K_1s_z^2-K_2s_z^4$ | `integrate_llg.m` | `UNRESOLVED`: axis and parameter source locators missing |
| DMI | Interlayer component field using mixed $J_3/J_4$ neighbor sets | `dmi.m`; `integrate_llg.m`; notebook pp. 4--5, 8 | `UNRESOLVED`: directed bonds and antisymmetry not verified |
| Dipole | Three same-layer neighbors at $0.5067\ \mathrm{nm}$, refreshed every 50 steps | `dipole_.m`; notebook p. 6 | `UNRESOLVED`: truncation justification and source locator missing |
| LLG/LLGS | Fixed-length LL form in tesla convention | GPU kernels; local equation reference | Implemented convention; kernel-level analytical tests missing |
| Boundary | Empty two-cell lateral padding; omitted out-of-range layer bonds | `distrib.m`, `exchangej_.m`, `dmi.m` | Finite-lattice notebook/code agreement audited above; broader coordination and finite-size tests missing |
| Domain-wall initial state | Two lateral domains with layer-dependent angles | `systemgeneration.m` | Implemented; physical domain mapping unresolved |

## Validation requirements

No automated tests or accepted benchmark outputs were found. One focused
manual structure check has passed: the corrected page-10 diagram agrees with
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
   domain-wall observables, and non-destructive output naming.
8. **Scientific comparison:** establish size convergence and compare an
   appropriate observable with a source-supported analytical, experimental, or
   published numerical reference.

Passing one stage does not validate later stages. In particular, an apparently
reasonable trajectory cannot validate the crystallographic mapping or bond
table.

## Known limitations and open decisions

- `thermalenable` and `T` do not generate thermal noise.
- The default distribution is deterministic; fixed-distribution loading does
  not validate compatibility with the site masks.
- No explicit lattice-vector, coordinate, or bond-table data structure exists.
- The exchange provenance and shell mapping are unresolved.
- DMI bond orientation and antisymmetry are not independently tested.
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
- Exact saved-time alignment across GPU chunks is unverified.
- `finalnew.mat` overwrites an existing output and stores the full workspace.
- No regression tests, convergence baselines, or documented validated run are
  present.

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
