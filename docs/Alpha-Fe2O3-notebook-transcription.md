# Alpha-Fe2O3 annotated notebook: digital transcription

## Status and use

This document is a page-addressable transcription of
[`Alpha-Fe2O3.pdf`](../Alpha-Fe2O3.pdf). It preserves the notebook's historical
notes, equations, diagram labels, questions, and tentative assignments in an
editable form.

It is **not** an independently verified model description or a scientific
source. In particular, the structure, interaction topology, parameter values,
DMI directions, and dipole approximation transcribed here remain subject to
the evidence and validation requirements in
[`model-methodology.md`](model-methodology.md).

Transcription conventions:

- Page numbers refer to the PDF page number, starting at 1.
- Stable item IDs such as `P05-DMI-DERIVATION` are provided for later discussion.
- “Literal note” records notebook wording without endorsing it.
- `[unclear]` marks handwriting or diagram content that cannot be read with
  sufficient confidence.
- `UNRESOLVED` means the notebook itself asks a question, conflicts internally,
  or lacks enough information to select a physical interpretation.
- Color names are diagram labels, not chemical species.
- Coordinate tuples and compact neighbor notation are retained as written.
  Their index semantics will be added after the author explains them.

## Page index

| PDF page | Stable ID | Notebook section |
|---:|---|---|
| 1 | `P01-CRYSTAL-STRUCTURE` | Crystal-structure comparison and repeating unit |
| 2 | `P02-EXCHANGE-SOURCES` | Exchange-shell numbering and two source tables |
| 3 | `P03-EXCHANGE-MAPPING` | $J_1$--$J_4$ diagram and source-table relabeling |
| 4 | `P04-DMI-SOURCES` | DMI shell assignment from source [2] |
| 5 | `P05-DMI-DERIVATION` | DMI directions, energy, and field derivation |
| 6 | `P06-DIPOLE-DERIVATION` | Dipole energy, field, and in-plane displacement vectors |
| 7--8 | `P07-P08-ATOMISTIC-ZHIFENG` | Six-layer lattice-to-array sketches |
| 9 | `P09-ATOMISTIC-ZHANG-XUE` | Typed notes about the code and domain-wall motion |
| 10 | `P10-NEIGHBOR-MAPS-A` | Clean layer masks and neighbor lists, first set |
| 11 | `P11-NEIGHBOR-MAPS-B` | Clean layer masks and neighbor lists, second set |
| 12 | `P12-REFERENCES` | Abbreviated references |

## Page 1 — Crystal structure

Source: [PDF page 1](../Alpha-Fe2O3.pdf#page=1)

### `P01-CRYSTAL-STRUCTURE-COMPARISON`

The page places two published crystal-structure illustrations side by side:

- Source `[1]`: a hexagonal cell containing Fe sites and alternating magnetic
  arrows. Handwritten labels divide the sites into three vertical sequences:
  black $A_1,B_1,C_1,D_1$, blue $A_2,B_2,C_2,D_2$, and orange
  $A_3,B_3,C_3,D_3$.
- Source `[2]`: a perspective cell whose Fe sites are similarly labeled
  $A_1$--$D_1$, $A_2$--$D_2$, and $A_3$--$D_3$. Adjacent panels show AFM and
  weak-ferromagnetic (WF) spin configurations. Handwritten numbers on the AFM
  chain mark $A=2$, $B=1$, $C=4$, and $D=3$.

Literal note:

> Comparing ① ② ③, these two structures are identical.

### `P01-CRYSTAL-STRUCTURE-REPEAT`

The lower diagram relabels a four-site vertical sequence as $A,B,C,D$ and
draws three outlined vertical groups through the cell.

Literal note:

> repeating one unit cell can produce the whole structure

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 2 — Exchange-interaction sources

Source: [PDF page 2](../Alpha-Fe2O3.pdf#page=2)

### `P02-EXCHANGE-NUMBERING`

The upper crystal diagram is annotated with site numbers $1$, $2$, $3$, $4'$,
and $4''$. A black outline marks a central vertical sequence. A line reading
“follow this numbering” points from the source-[3] coordination-shell diagram
to the source-[2] structural diagram.

Open handwritten question:

> why $r_{2-1'} \ne r_{2-4''}$ ?

Status: `UNRESOLVED`.

### `P02-SOURCE-3-TABLE`

The source-[3] figure centers atom 2 and labels coordination shells $1$, $3'$,
$1'$, $4'$, $4''$, and $3$. Its table is transcribed below in meV.

| Interaction | Samuelsen *et al.* (Ref. 30) | LSDA | LDA+$U$ |
|---|---:|---:|---:|
| $J^{(1)}$ | $-3.076$ | $9.905$ | $8.576$ |
| $J^{(3')}$ | $-0.528$ | $-5.71$ | $-7.3$ |
| $J^{(1')}$ | $20.313$ | $25.957$ | $25.224$ |
| $J^{(4')}$ and $J^{(4'')}$ | $12.554$ | $13.488$ | $17.502$ |
| $J^{(3)}$ | $1.056$ | $-0.497$ | $-0.073$ |

The printed caption says negative couplings are ferromagnetic under that
source's convention.

### `P02-SOURCE-2-TABLE`

The source-[2] table lists five isotropic exchange interactions in meV.
Handwritten green labels connect each sublattice pair to the numbered diagram.

| Interaction | Sublattices | Handwritten pair | This paper | Ref. [9] | Ref. [8] |
|---|---|---|---:|---:|---:|
| $J_1$ | A--B | 2--1 | $-3.21$ | $-3.5$ | $-8.58$ |
| $J_2$ | A--D | 2--3 | $-3.84$ | $-3.2$ | $7.3$ |
| $J_3$ | A--B | 2--1 | $-26.10$ | $-13.9$ | $-25.22$ |
| $J_4$ | A--C | 2--4 | $-15.71$ | $-9.8$ | $-17.5$ |
| $J_5$ | A--D | 2--3 | $-0.17$ | — | $0.07$ |

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 3 — Exchange mapping

Source: [PDF page 3](../Alpha-Fe2O3.pdf#page=3)

### `P03-EXCHANGE-ARROWS`

Two copies of the crystal diagram show the neighbor groups of two differently
positioned central atoms. The arrows are color-coded:

- green: $J_1$;
- black: $J_2$;
- orange: $J_3$;
- blue: $J_4$.

Literal note:

> 同理，可得 atom 3, 4 的 exchange interaction

Translation: “Similarly, the exchange interactions for atoms 3 and 4 can be
obtained.”

### `P03-EXCHANGE-RELABELING`

The handwritten relabeling beside the source-[3] table is:

$$
J_1 = J^{(1)},\qquad
J_2 = J^{(3')},\qquad
J_3 = J^{(1')},\qquad
J_4 = J^{(4')}\ \text{and}\ J^{(4'')}.
$$

The LDA+$U$ entries $25.224$ and $17.502$ are circled. Literal note:

> $J_3,J_4$ 最大

Translation: “$J_3$ and $J_4$ are the largest.”

Status: the mapping is a notebook interpretation and remains `UNRESOLVED` as
material-specific evidence.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 4 — DMI source diagrams

Source: [PDF page 4](../Alpha-Fe2O3.pdf#page=4)

### `P04-DMI-SHELLS`

The left structure repeats the numbered Fe sites used on page 2. Source `[2]`
shows DMI vectors around a central site, with bond groups labeled A--B, A--C,
and A--D. Handwritten overlays associate:

- the third shell with pair 2--$1'$ and three sites;
- the fourth shell with pairs 2--$4'$ and 2--$4''$ and six sites.

The clipped source sentence reads:

> The dominant contributions come from the third and fourth shells surrounding
> the central atom, with three and six sites on the shells, respectively.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 5 — DMI construction and derivation

Source: [PDF page 5](../Alpha-Fe2O3.pdf#page=5)

### `P05-DMI-DIRECTIONS`

Two crystal diagrams show orange $D_1$ arrows and blue $D_2$ arrows from two
different central sites. Figure 1 contains red question marks and the notes:

> there is DMI pointing down? `[word uncertain]`

> should point up? `[word uncertain]`

Literal note:

> 同理，可得 atom 3, 4 的 DMI

Translation: “Similarly, the DMI for atoms 3 and 4 can be obtained.”

Status: vector directions and the “pinning” annotations are `UNRESOLVED`.

### `P05-DMI-ENERGY`

The notebook writes the pair energy as

$$
E_{\mathrm{DMI}}=
\mathbf{D}_{ij}\cdot(\mathbf{S}_i\times\mathbf{S}_j).
$$

Literal note:

> 由于 $C_3$ 对称，$\mathbf{D}$ 只有 $z$ 分量，$D_z$ 向上与向下方向相反。

Translation: “Because of $C_3$ symmetry, $\mathbf{D}$ has only a $z$
component, and $D_z$ has opposite signs in the upward and downward
directions.”

The example signs drawn for Figure 1 are:

- $D_{1,z}<0$;
- for $D_{2,z}$, the three upward arrows are $>0$ and the three downward
  arrows are $<0$.

### `P05-CROSS-PRODUCT`

$$
\mathbf{S}_i\times\mathbf{S}_j=
\begin{vmatrix}
\mathbf{x}&\mathbf{y}&\mathbf{z}\\
S_{ix}&S_{iy}&S_{iz}\\
S_{jx}&S_{jy}&S_{jz}
\end{vmatrix}
=
\begin{bmatrix}
S_{iy}S_{jz}-S_{iz}S_{jy},
S_{iz}S_{jx}-S_{ix}S_{jz},
S_{ix}S_{jy}-S_{iy}S_{jx}
\end{bmatrix}.
$$

### `P05-DMI-DERIVATION`

For $\mathbf{D}$ along the $-z$ direction, the notebook derives

$$
\begin{aligned}
E={}&-\mathbf{D}\cdot(\mathbf{S}_{i-1}\times\mathbf{S}_i)
-\mathbf{D}\cdot(\mathbf{S}_i\times\mathbf{S}_{i+1})\\
={}&D(S_{i-1,x}S_{i,y}-S_{i-1,y}S_{i,x})
+D(S_{i,x}S_{i+1,y}-S_{i,y}S_{i+1,x})\\
={}&D(S_{i-1,x}S_{i,y}-S_{i-1,y}S_{i,x}
+S_{i,x}S_{i+1,y}-S_{i,y}S_{i+1,x}),
\end{aligned}
$$

followed by

$$
\mathbf{H}_{\mathrm{DMI}}
=-\frac{1}{\mu}\frac{\partial E_{\mathrm{DMI}}}{\partial\mathbf{S}_i}
=-\frac{D}{\mu}
\begin{bmatrix}
-S_{i-1,y}+S_{i+1,y},
S_{i-1,x}-S_{i+1,x},
0
\end{bmatrix}.
$$

This is a literal transcription; the sign, directed-bond convention, and
meaning of “up/down” have not been reconciled here.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 6 — Dipole energy and field

Source: [PDF page 6](../Alpha-Fe2O3.pdf#page=6)

### `P06-DIPOLE-SPIN-FORM`

The first energy expression is

$$
E_{\mathrm{dipole}}
=-\frac{\mu_s^2\mu_0}{8\pi}
\sum_{i\ne j}
\left[
\frac{3(\mathbf{S}_i\cdot\mathbf{r}_{ij})
(\mathbf{r}_{ij}\cdot\mathbf{S}_j)}{|\mathbf{r}_{ij}|^5}
-\frac{\mathbf{S}_i\cdot\mathbf{S}_j}{|\mathbf{r}_{ij}|^3}
\right].
$$

The following field is written as

$$
\mathbf{H}_{\mathrm{dipole}}
=-\frac{1}{\mu_s}\frac{\partial E}{\partial\mathbf{S}_i}
=\frac{\mu_s\mu_0}{8\pi}
\sum_{i\ne j}
\left[
\frac{3\mathbf{r}_{ij}(\mathbf{r}_{ij}\cdot\mathbf{S}_j)}
{|\mathbf{r}_{ij}|^5}
-\frac{\mathbf{S}_j}{|\mathbf{r}_{ij}|^3}
\right].
$$

### `P06-DIPOLE-MOMENT-FORM`

The notebook then rewrites the pair energy using magnetic moments:

$$
E_{\mathrm{dipole}}
=-\frac{\mu_0}{4\pi}\sum_{j\ne i}
\left[
\frac{3(\mathbf{R}_{ij}\cdot\boldsymbol{\mu}_i)
(\mathbf{R}_{ij}\cdot\boldsymbol{\mu}_j)}{R_{ij}^5}
-\frac{\boldsymbol{\mu}_i\cdot\boldsymbol{\mu}_j}{R_{ij}^3}
\right],
$$

and substitutes $\boldsymbol{\mu}_i=\mu_s\mathbf{S}_i$ and
$\boldsymbol{\mu}_j=\mu_s\mathbf{S}_j$. The final field line is

$$
\mathbf{H}_{\mathrm{dipole}}
=-\frac{1}{\mu_i}\frac{\partial E}{\partial\mathbf{S}_i}
=\frac{\mu_0}{4\pi}
\left[
\frac{3\mathbf{R}_{ij}
(\mathbf{R}_{ij}\cdot\mu_s\mathbf{S}_j)}{R_{ij}^5}
-\frac{\mu_s\mathbf{S}_j}{R_{ij}^3}
\right].
$$

The last line is transcribed as written; an explicit summation sign is not
visible there.

### `P06-IN-PLANE-VECTORS`

For “layer 1 green atoms,” a honeycomb sketch identifies three in-plane
displacements:

$$
\begin{aligned}
\mathbf{r}_1&=[\mathrm{dist},0,0],\\
\mathbf{r}_2&=\mathrm{dist}
\left[-\sin\frac{\pi}{6},\cos\frac{\pi}{6},0\right],\\
\mathbf{r}_3&=\mathrm{dist}
\left[-\sin\frac{\pi}{6},-\cos\frac{\pi}{6},0\right].
\end{aligned}
$$

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Pages 7--8 — Atomistic model: Zhifeng

Sources: [PDF page 7](../Alpha-Fe2O3.pdf#page=7) and
[PDF page 8](../Alpha-Fe2O3.pdf#page=8)

### `P07-P08-SIX-LAYER-DISTRIBUTION`

The section is titled “Atomistic model: Zhifeng.” Two structure diagrams at
the top of page 7 repeat the $J_1$--$J_4$ arrow construction for representative
sites. The remainder sketches six physical layers and their encoding into
rectangular array cells:

- Layer 1: two alternating hollow-site colors, green and orange.
- Layer 2: blue perimeter sites and purple interior sites.
- Layer 3: red perimeter sites and black interior sites.
- Layer 4: a hexagonal ring with six hollow perimeter markers; a short
  annotation under the layer label is `[unclear]`.
- Layer 5: a hollow center plus three hollow perimeter markers on alternating
  vertices.
- Layer 6: a hollow center plus three hollow perimeter markers in the
  complementary orientation.

Pages 10--11 contain cleaner versions of the three-layer array masks and are
transcribed below.

Literal notes on page 8:

> exchange interaction: See two pages later

> DMI: the DMI of each atom is identical to the corresponding $J_3$, $J_4$

The word “identical” is retained literally; whether it means identical neighbor
indices, shell membership, or another correspondence is `UNRESOLVED`.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 9 — Atomistic model: Zhang Xue

Source: [PDF page 9](../Alpha-Fe2O3.pdf#page=9)

### `P09-CODE-NOTES`

Typed text, lightly normalized only for spacing:

> The atom distribution I used in the code. (p3–p4)
>
> The exchange number is 1 3 3 6 for each exchange interaction.
>
> Since the boundary for every layer and atom is different, I try to add a
> shell in the code to simplify the calculation.

### `P09-SWITCHING-NOTES`

> In addition, I have read Dr Fan’s work [Phys. Rev. App. 11 034018 (2019)]. I
> feel this work is suddenly stopped. Recently, I have been thinking about why
> I have this feeling, and I think the rotation of the field is just like the
> type-x SOT switching (Maybe).
>
> As for Zhenghang’s work, I have read some papers about DW-induced SOT
> switching. This is a very good way to characterize this non-uniform behavior
> with DW velocity. However, our DW is curved; how to define this velocity
> becomes important. As shown in Prof. Yang Hyunsoo’s work in Communication
> Physics (2018), the nonuniform switching DW velocity is defined as

$$
v_{\mathrm{DW,fwd}}=\frac{d}{2t_{c,\mathrm{fwd}}},
$$

> where $t$ represents the time corresponding to $P_{\mathrm{sw}}=50\%$ during
> the first forward switching.
>
> Maybe we can consider the moving speed of DW on the diagonal line (if DW is
> from the lower left corner to the upper right corner). I mentioned to
> Zhenhang that we need to consider this matter.

A square sketch shows a diagonal from the lower-left to the upper-right corner
and a left-pointing arrow at the right edge.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 10 — Layer masks and exchange-neighbor maps, first set

Source: [PDF page 10](../Alpha-Fe2O3.pdf#page=10)

### `P10-LAYER-MASKS`

The heading is “空心 up” (“hollow sites up”). Clean diagrams show the Layer 1,
Layer 2, and Layer 3 honeycomb patterns, a compact $3\times3$ cell encoding for
each, and its horizontal repetition. Layer 2 and Layer 3 are shown as a merged
encoding.

The updated page contains a corrective annotation on the expanded Layer-2
grid. “empty” points to column 12 of the middle row, and the adjacent note says
“the code correctly implement this.” For the displayed 13-column unpadded
core, column 13 is blue; the next purple site would be column 14 and therefore
falls outside the displayed grid.

Code correspondence was checked directly against `distrib.m` with
`natomW=3`, `natomL=13`, and `natomH=3`. Using `.` for an empty cell, the six
color/site masks are:

```text
Layer 1: G.O.G.O.G.O.G
         O.G.O.G.O.G.O
         G.O.G.O.G.O.G

Layer 2: ..BP..BP..BP.
         BP..BP..BP..B
         ..BP..BP..BP.

Layer 3: R..BR..BR..BR
         .BR..BR..BR..
         R..BR..BR..BR
```

The corrected page-10 site-distribution diagram and the executable masks are
therefore identical for this representation. This establishes
notebook-to-code consistency only; it does not independently validate the
mapping as the crystallographic structure of hematite.

### `P10-LAYER1-NEIGHBORS`

The following compact coordinate lists are transcribed as written. A leading
`0` or `1`, trailing `2p`, `gr`, or `or`, and text colors are diagram keys;
their exact semantics are to be supplied.

```text
J1: (A, i, 2)  (A, i, 3)

J2: (A, i+2) (A-1, i) (A+1, i)
    (A, i-2) (A-1, i) (A+1, i)

J3: 0 (A, i-1, 2) (A+1, i+1, 2) (A-1, i+1, 2) 2p
    0 (A, i+1, 3) (A+1, i-1, 3) (A-1, i-1, 3) 2p

J4: 0 (A, i+2, 2) (A+1, i, 2) (A-1, i, 2) or
    0 (A, i-1, 3) (A+1, i+1, 3) (A-1, i+1, 3) 2p
    0 (A, i+1, 2) (A-1, i-1, 2) (A+1, i-1, 2) 2p
    0 (A, i-2, 3) (A-1, i, 3) (A+1, i, 3) gr
```

Handwritten color annotations refer to Layer-2 blue/purple and Layer-3
red/black sites.

### `P10-LAYER2-NEIGHBORS`

```text
J1: (A, i, 3)  (A, i, 1)

J2: 1 (A, i+1) (A-1, i-1) (A+1, i-1)
    1 (A, i-1) (A-1, i+1) (A+1, i+1)

J3: 0 (A, i+1, 1) (A+1, i-1, 1) (A-1, i-1, 1) gr
    0 (A, i-2, 3) (A+1, i, 3) (A-1, i, 3) gr

J4: 0 (A, i+1, 3) (A+1, i-1, 3) (A-1, i-1, 3) 2p
    0 (A, i-2, 1) (A+1, i, 1) (A-1, i, 1) gr
    0 (A, i+1, 3) (A+1, i-1, 3) (A-1, i-1, 3) gr
    0 (A, i-1, 1) (A+1, i+1, 1) (A-1, i+1, 1) or
```

The apparent repeated $J_4$ triplet is retained exactly as visible and is
`UNRESOLVED`.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 11 — Layer masks and exchange-neighbor maps, second set

Source: [PDF page 11](../Alpha-Fe2O3.pdf#page=11)

### `P11-LAYER-MASKS`

The page repeats the Layer 1, Layer 2, and Layer 3 diagrams. The horizontally
repeated middle grid overlays the Layer-2 and Layer-3 encodings. Handwritten
labels identify Layer-2 purple and Layer-1 green, and subsequent notes name
same-layer red/black plus Layer-1 or Layer-2 color targets.

### `P11-NEIGHBOR-LISTS`

```text
J1: (A, i, 2)  (A, i, 1)

J2: 1 (A, i+1) (A-1, i-1) (A+1, i-1)
    1 (A, i-1) (A-1, i+1) (A+1, i+1)

J3: 1 (A, i-1, 1) (A+1, i+1, 1) (A-1, i+1, 1) or
    1 (A, i+2, 2) (A+1, i, 2) (A-1, i, 2) or

J4: 1 (A, i-1, 2) (A+1, i+1, 2) (A-1, i+1, 2) or
    1 (A, i+1, 1) (A+1, i-1, 1) (A-1, i-1, 1) gr
    1 (A, i-1, 2) (A+1, i+1, 2) (A-1, i+1, 2) 2p
    1 (A-1, i, 1) (A+1, i, 1) (A, i+2, 1) or
```

The two visibly identical $J_4$ coordinate triplets have different colored
leading markers and different trailing keys. They are therefore retained as
separate rows rather than deduplicated.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Page 12 — Abbreviated references

Source: [PDF page 12](../Alpha-Fe2O3.pdf#page=12)

### `P12-REFERENCES`

The notebook provides only the following abbreviated references:

1. `[1] 2020, PRL, Yang Cheng`
2. `[2] 2023, PRB, Tobias Dannegger`
3. `[3] 2005, PRB, V. V. Mazurenko`

Exact titles, author lists, volumes, pages/article numbers, DOIs, and the
specific figure/table locators have not been supplied in the notebook and are
therefore `REFERENCE MISSING` in this transcription.

User meaning: _To be supplied._

Code correspondence: _To be supplied._

## Follow-up mapping worksheet

Use this table in the next discussion to record the author's explanations
without changing the literal transcription above.

| Stable item ID | Author-supplied meaning | Code file/function/lines | Evidence or validation note |
|---|---|---|---|
| `P01-CRYSTAL-STRUCTURE-COMPARISON` |  |  |  |
| `P01-CRYSTAL-STRUCTURE-REPEAT` |  |  |  |
| `P02-EXCHANGE-NUMBERING` |  |  |  |
| `P02-SOURCE-3-TABLE` |  |  |  |
| `P02-SOURCE-2-TABLE` |  |  |  |
| `P03-EXCHANGE-ARROWS` |  |  |  |
| `P03-EXCHANGE-RELABELING` |  |  |  |
| `P04-DMI-SHELLS` |  |  |  |
| `P05-DMI-DIRECTIONS` |  |  |  |
| `P05-DMI-ENERGY` |  |  |  |
| `P05-CROSS-PRODUCT` |  |  |  |
| `P05-DMI-DERIVATION` |  |  |  |
| `P06-DIPOLE-SPIN-FORM` |  |  |  |
| `P06-DIPOLE-MOMENT-FORM` |  |  |  |
| `P06-IN-PLANE-VECTORS` |  |  |  |
| `P07-P08-SIX-LAYER-DISTRIBUTION` |  |  |  |
| `P09-CODE-NOTES` |  |  |  |
| `P09-SWITCHING-NOTES` |  |  |  |
| `P10-LAYER-MASKS` |  |  |  |
| `P10-LAYER1-NEIGHBORS` |  |  |  |
| `P10-LAYER2-NEIGHBORS` |  |  |  |
| `P11-LAYER-MASKS` |  |  |  |
| `P11-NEIGHBOR-LISTS` |  |  |  |
| `P12-REFERENCES` |  |  |  |
