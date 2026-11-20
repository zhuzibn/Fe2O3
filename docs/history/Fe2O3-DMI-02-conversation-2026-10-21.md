# Fe2O3: Historical Commit Comparison and the DMI-02 Antisymmetry Problem

**Exported:** 2026-10-21

> **Historical status — archived 2026-11-20:** This conversation preserves the
> discussion and reasoning from 2026-10-21. References to “current code,” local
> file paths, line numbers, and test results describe the versions discussed
> at that time and may no longer apply.
>
> For current DMI-02 status, implemented behavior, and remaining evidence gaps,
> see the [maintained model methodology](../model-methodology.md#dzyaloshinskii--moriya-interaction).
> Use that maintained record when assessing the present implementation.
>
> The original conversation below is preserved verbatim. Its historical test
> results and literature claims were not revalidated when this file was archived.

---

## User

does these issue "Deferred issue ANI-01, Deferred issue OUT-01, Deferred issue TIME-01, Unresolved source issue DMI-02" also in commit "35bd00b"?

## Assistant

**Only TIME-01 is present in `35bd00b` among these four issues.**

| Issue | In `35bd00b`? | Evidence |
|---|---|---|
| **ANI-01** | **No** | Uses second- and fourth-order polynomial anisotropy, without the problematic angular division. The +z MATLAB probe returned finite values. |
| **OUT-01** | **No / not applicable** | Does not record the `hdmi_allx/y/z` histories that have the zero-frame and alignment defects. |
| **TIME-01** | **Yes** | Each 10,000-slot chunk performs only 9,999 updates. Its default 200 chunks evolve for **199.98 ps**, rather than 200 ps; saved time labels are also misaligned. |
| **DMI-02** | **Not reproduced** | Testing its actual DMI operator on the same 140-site lattice found **0 antisymmetry mismatches**, compared with **544** in the current code. |

I inspected the commit directly and ran CPU MATLAB probes without running the full simulation. Passing the DMI antisymmetry check does **not** establish that its physical topology or parameters are correct.

Documentation: no update needed — read-only historical comparison; repository files unchanged.

## User

last time you mention the current code correctly reflect the journal paper, can you explain the DMI-02 problem (antisymmetry) more clearly

## Assistant

**I should correct my earlier statement: the current code has not been verified as a faithful implementation of the journal paper.** Saying that it “correctly reflects the paper” was too strong. We verified several implementation details and agreement with the supplied `code_20261013` source, but **DMI-02 reveals an inconsistency in the implemented DMI field**.

The paper describes DMI through a pair-interaction energy—the antisymmetric part of the exchange tensor. See [Dannegger et al., *Physical Review B* **107**, 184426 (2023), author manuscript, Sec. III D, Eqs. (2)–(4)](https://arxiv.org/html/2301.13094#S3.SS4).

**What antisymmetry means here**

Consider two spins, $i$ and $j$. Using the repository’s energy-sign convention, their DMI energy is

$$
E_{ij}=-\mathbf D_{ij}\cdot(\mathbf s_i\times\mathbf s_j).
$$

Describing the **same bond** with its labels reversed must give the same energy. But the cross product changes sign:

$$
\mathbf s_j\times\mathbf s_i
=-\mathbf s_i\times\mathbf s_j.
$$

Therefore its directed DMI coefficient must also change sign:

$$
\boxed{\mathbf D_{ji}=-\mathbf D_{ij}.}
$$

This is a consistency requirement for describing one physical bond from its two ends. It does **not** require the actual fields on arbitrary spins to be equal and opposite. An alternative overall sign convention for the energy does not change this requirement.

**How this applies to your code**

The current [field assembly]( /mnt/c/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/integrate_llg.m:186) can be written as

$$
\mathbf B_i^{D}
=\frac{D}{\mu_s}\sum_j C_{ij}
(s_{jy},-s_{jx},0).
$$

Here $C_{ij}$ is the signed coefficient with which spin $j$ contributes to the field at site $i$, after combining all masks, shifts, and `nex + pre` terms.

For this to represent ordinary pairwise $z$-directed DMI, we need

$$
\boxed{C_{ji}=-C_{ij}.}
$$

For the actual padded MATLAB sites

- $i=(4,3,2)$, a blue site;
- $j=(3,3,1)$, a green site;

the code produces:

| Contribution | Current code | Consistent pair, keeping $C_{ij}=+1$ |
|---|---:|---:|
| Spin $j$ contributing to site $i$: $C_{ij}$ | $+1$ | $+1$ |
| Spin $i$ contributing to site $j$: $C_{ji}$ | $+1$ | $-1$ |

You can trace these directly to the positive `mmx_d1` term in the blue-site branch at [dmi.m:69](/mnt/c/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/dmi.m:69) and the positive `mmx_u1` term in the green-site branch at [dmi.m:29](/mnt/c/Users/zzf-m/OneDrive/code_softwares/spin_dynamics/Fe2O3/dmi.m:29).

**Why the same signs are a problem**

For a bond with $\mathbf D_{ij}=D\hat{\mathbf z}$, differentiating its energy gives both fields:

$$
\mathbf B_i=\frac{D}{\mu_s}(s_{jy},-s_{jx},0),
\qquad
\mathbf B_j=\frac{D}{\mu_s}(-s_{iy},s_{ix},0).
$$

The current example instead gives the **negative of the required second field**, while retaining the first field.

Thus the two sites are not receiving fields derived consistently from the same DMI bond energy. This is a field-calculation problem, separate from the saved-history problem OUT-01. It can affect trajectories and equilibrium predictions; the size of that effect has not been quantified.

The previous CPU MATLAB check found:

| Version | Occupied sites | Antisymmetry violations |
|---|---:|---:|
| `35bd00b` | 140 | 0 |
| Current code | 140 | 544 matrix entries, representing 272 unordered pairs |

These are exact integer-coefficient mismatches, not numerical roundoff. The old version passes this check, but that alone does not establish that its physical bonds or parameters match the paper.

**What remains to establish**

The paper’s Fig. 4 discusses three-dimensional bond DMI vectors and their symmetry-related contributions. Its effective $z$-component result does not supply a verified mapping to every signed MATLAB neighbor term. [Author manuscript, Sec. III D and Fig. 4](https://arxiv.org/html/2301.13094#S3.SS4).

A justified correction therefore requires mapping the physical bonds to array sites, assigning each directed $\mathbf D_{ij}$, and deriving both endpoint fields from one energy. Merely changing the global sign of `Dsim` would turn $(+1,+1)$ into $(-1,-1)$ and leave DMI-02 unresolved.

Documentation: no update needed — explanation only; no files changed.

## User

$export-conversation
