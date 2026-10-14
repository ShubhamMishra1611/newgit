# T3: Λ selection, head to head. (A) Discretuum plus anthropics vs (B) swampland quintessence

**Labels.**
- **ESTABLISHED**: cited by arXiv ID and equation number. The abstract page or full text was fetched in this session; source texts are saved in `sources/`.
- **CONJECTURE**: a claim made in the literature and not proven.
- **OUR ATTEMPT**: computed here; the script is named.

All code and outputs are in `phase3/targets/T3/`. Units: M_P = (8πG)^(-1/2) = 1 for the slopes.

## 0. Headline numbers

| Quantity | Result | Label / script |
|---|---|---|
| BP: N(Λ ∈ [0, ε]) ∝ ε^p, J = 4–12 | p = 0.994 to 1.008. Local density = our face-corrected smooth formula to ≤ 0.2% | OUR ATTEMPT, `T3_A1_bp_discretuum.py` |
| BP: nearest-neighbour spacings near Λ = 0 (6×10^5 levels per J) | Poisson: var s = 0.996–1.005, P(s < 0.1) = 0.095. KS p = 0.06–0.68 | OUR ATTEMPT |
| BP: min abs Λ | Exponential, mean 1/(2ρ(0)), with ρ(0) ∝ R^(J−2). [[A1C_POOLED]] | OUR ATTEMPT, `T3_A1c_min_statistics.py` |
| Rigid-CY factor of 2 | A convention. Ashok–Douglas (AD) and Denef–Douglas (DD) count N and −N separately. Phase 1 counted modulo the full SL(2,ℤ) | OUR ATTEMPT + ESTABLISHED (text), `T3_A2_rigid_cy_factor2.py` |
| MSW P(ρ_V ≤ ρ_obs), Planck 2018, R_G = 1 Mpc, s = 1 | **0.086** (range 0.071–0.140 over s). Median ρ_V/ρ_obs = 13.0 | OUR ATTEMPT, `T3_A3_msw_anthropic.py` |
| The same over filter scale R_G = 0.5–3 Mpc | P = 0.043–0.29. Median 30 down to 2.5 | OUR ATTEMPT |
| Minimal late-time slope, all-positive geometric type II sources (dilaton–volume plane) | IIB √(18/7) = 1.604 (R₆ + D9). IIA √(50/19) = 1.622 (R₆ + F₀/D8). √2 if positive non-geometric Q flux is allowed | OUR ATTEMPT, `T3_B1_swampland_slopes.py` |
| Formal mixed-sign minimum | √(2/3) (TCC value), only along strong coupling or zero volume | OUR ATTEMPT |
| DESI DR2 exponential fit λ = 0.698 (+0.173/−0.202) | 4.1σ below √2; 5.2σ below 1.604 (one-sided, using the upper error) | ESTABLISHED (2504.04226 abstract) + OUR ATTEMPT |

## 1. Baseline reproductions (checks of the machinery)

- **BP**, hep-th/0004134:
  - λ = λ_bare + ½Σn_i²q_i² (2.21); shell (2.22)–(2.24); spacing (2.26). The degeneracy is D = 2^J from n_i → −n_i (text after (2.29)).
  - The counting code matches brute force exactly for J = 3–5 (scratch test).
  - With D = 2^J, (2.26) matches the measured density at J = 4 (ratio 1.0002). At small n_max it underestimates, by 2%, 15%, 50% and 116% at J = 6, 8, 10, 12. The faces n_i = 0 are missing from (2.26) (§2a).
- **AD**, hep-th/0307049, §3.2 (T² toy), reproduced:
  - N_vac(L) = 2σ(L);
  - Σσ ~ π²N²/12;
  - the continuum volume (3.11), which gives π²L/3 = d/dL[2·π²L²/12];
  - vol(F) = ∫_F d²τ/(4y²) = π/12 (2.23).
- **DD**, hep-th/0404116, §5.1 (rigid CY, read from the arXiv text this session), reproduced:
  - "240 vacua at τ = 2i" (L = 150). We count 240 fluxes, which is 120 SL(2,ℤ) orbits.
  - Fig. 7: counts in discs around τ = 2i against the estimate 2πAL², with A in the metric 1/(4y²). At R = 0.15/0.20/0.25/0.30 the estimate is 627/1119/1756/2541. The counts with N and −N separate are 632/1202/1734/2592; the orbit counts are half of those.
- **MSW**, astro-ph/9701099 (ApJ text):
  - ⟨ρ_V⟩/σ³ρ̄ = 2.579/3.010/3.224 for s = 0/1/∞ (eq. 25): reproduced exactly.
  - The s = 1 quantiles 0.0590/1.49/11.1 σ³ρ̄ (text after eq. 30): we get 0.05897/1.485/11.11.
  - The closed form (29) agrees with direct integration of (9) to 10⁻⁶.
  - Table 1, all 18 entries: our σ_rec from their eqs. (39)–(52) agrees to the 2 printed digits (e.g. λ₀ = 0.7, R_G = 1 Mpc: 0.00333 vs 0.0033).
  - Our P(≤ρ*) agrees within the rounding of their σ (≤ 6% relative). Examples: λ₀ = 0.7 gives 0.053/0.121 against their 0.054/0.12; λ₀ = 0.9 at R_G = 1 Mpc, using their rounded σ, gives 0.878 against 0.90.
- **HKTT**, 0711.2512:
  - ρ = Vol^(1/3), τ = e^(−φ)√Vol (2.3). Canonical ρ̂ = √(3/2) ln ρ, τ̂ = √2 ln τ (2.5). We verified (2.5) exactly with a sympy minisuperspace reduction of the 10d string-frame action (difference 0).
  - Scalings reproduced from a single rule, V_E ∝ ρ^(w−3χ/2) τ^(χ−4): (3.3), (3.4), (4.1), (4.2), (4.4).
- **STT**, 2303.03418, eq. (II.3) (dilaton coupling² ≥ 4/(d−2)): in our basis every term has a τ̂-component (4−χ)/√2 ≥ √2.
- **Planck 2018**, 1807.06209, Table 2, column TT,TE,EE+lowE+lensing: Ω_bh² = 0.02237, Ω_ch² = 0.1200, H₀ = 67.36, n_s = 0.9649, ln(10¹⁰A_s) = 3.044, τ = 0.0544, Ω_m = 0.3153, Ω_Λ = 0.6847, σ₈ = 0.8111. CAMB 2.0.4 with one 0.06 eV neutrino gives σ₈ = 0.8109.

## 2. (A1) Bousso–Polchinski discretuum near Λ = 0

**Setup (OUR ATTEMPT, `bp_lib.py`, `T3_A1_bp_discretuum.py`).**
- Λ = (Q − R²)/2, with Q = Σq_i²n_i² and R² = 2|λ_bare|.
- We count distinct Λ, i.e. the orthant n_i ≥ 0, which removes BP's D = 2^J.
- Charges q_i ∈ [1, 2) are random. As float64 numbers they are rationally independent at these sizes.
- Exact enumeration of the shell uses meet-in-the-middle (MITM). Sums are in 80-bit long double; the smallest abs Λ values are re-checked in exact rational arithmetic.
- "Flux cutoff" means n_max = R/q_min, the largest flux on the shell.

**Smooth density with faces (OUR ATTEMPT, exact identity).**
- Exactly, #{n ∈ ℤ^J_≥0 : Q ≤ X} = 2^(−J) Σ_{U⊆[J]} #(ℤ^U ∩ E_X). Each point with support T is counted 2^|T| times on the right, inside 2^(J−|T|) subsets U.
- Replacing each count by its volume gives ρ(Λ) = 2^(−J) Σ_k ω_{k−1} e_k(1/q) X^((k−2)/2), where X = R² + 2Λ and e_k is the k-th elementary symmetric polynomial.
- This matches all measured densities to 0.02–0.2%. BP's (2.26) keeps only k = J.

**(a) Flatness.** Columns: fitted exponent p of N(Λ ∈ [0, ε]) ∝ ε^p over windows with N ≥ 10 and ε ≤ 10⁻³R²; then the pull test on disjoint shells (`T3_A1_summary.py`).

| J | R | n_max | p | χ²/dof on disjoint shells [[A1S]] |
|---|---|---|---|---|
| 4 | 1500 | 1458 | 1.0011 | [[A1S4]] |
| 6 | 330 | 280 | 0.9937 | [[A1S6]] |
| 8 | 90 | 90 | 1.0009 | [[A1S8]] |
| 10 | 46 | 34 | 1.0083 | [[A1S10]] |
| 12 | 30 | 29 | 0.9955 | [[A1S12]] |

- In each case ε runs from about 3 mean spacings up to 0.02R², which is 8–10 decades.
- We find no gaps and no power-law deviation from ε¹ down to the single-level scale. The count tracks the smooth formula, so the "prior" is flat up to its O(Λ/R²) curvature.

**(b) Spacing statistics** (about 6×10^5 consecutive levels centred on Λ = 0, unfolded with the smooth density).

| J | var(s) | P(s < 0.1) | P(s < 0.01) | KS vs e^(−s): p | Σ²(100) | Σ²(1000) |
|---|---|---|---|---|---|---|
| 4 | 1.000 | 0.0948 | 0.0101 | 0.65 | 98 | 932 |
| 6 | 1.002 | 0.0949 | 0.0100 | 0.68 | 96 | 944 |
| 8 | 1.005 | 0.0952 | 0.0099 | 0.65 | 102 | 1030 |
| 10 | 0.996 | 0.0954 | 0.0099 | 0.61 | 98 | 966 |
| 12 | 1.004 | 0.0954 | 0.0102 | 0.06 | 98 | 961 |

- Reference values: Poisson gives 1 / 0.0952 / 0.00995; GOE gives var 0.286 and P(s<0.1) = 0.0078.
- These are Poisson statistics. There is no level repulsion and no rigidity: Σ²(L) ≈ L up to L ~ 10³. At L ≳ 10⁴ the estimator noise (about 20 independent windows) prevents a test for Berry saturation.
- This is the Berry–Tabor behaviour of flat-torus spectra. Rigorous Poisson pair correlation is known in related cases: Marklof, math/0210198 (inhomogeneous forms), and the 2D flat-torus results summarized in 2305.18261. We did not find a theorem for the exact BP case: homogeneous diagonal forms with J ≥ 3 and the sign identification.

**(c) Minimum abs Λ.**
- Prediction for a Poisson process: min abs Λ ~ Exp with mean 1/(2ρ(0)). Here ρ(0) ≈ ω_{J−1}R^(J−2)/(2^J Π q_i) plus face terms. At fixed J this means min abs Λ ∝ n_max^(−(J−2)), up to the face corrections, which matter at small n_max.
- Test, `T3_A1c_min_statistics.py`: 40 shell radii per (J, R), with J = 4…24.
- [[A1C_TABLE]]
- At fixed flux cutoff R = 11 (`T3_A1_bp_discretuum.out`, part c′), [[A1CP]].

**(d) Sensitivity to the charges.**
- [[A1D1]]
- [[A1D2]]

**What A1 does and does not say.**
- *Does:* for a fixed positive-definite quadratic form with rationally independent charges, the Λ values near 0 form a locally Poisson, gap-free sequence. Their density is known analytically, including faces. So "flat prior near 0" is exact at every scale we can resolve.
- *Does:* flatness is a property of incommensurability, not of the number of vacua. Repeated or near-commensurate charges produce combs and gaps, and reduce the effective J.
- *Does not:* real flux landscapes are not BP quadratic forms. In them:
  - Λ depends on stabilized moduli and on back-reaction.
  - Fluxes are bounded by a tadpole, not by a shell.
  - Duality groups add degeneracies of the kind BP warned of (D up to E₇₍₇₎).
  - Finding the small-Λ vacuum is NP-hard (hep-th/0602072, not re-read here).
- None of these was tested here.

## 3. (A2) The rigid-CY factor of 2. Verdict: convention, documented by AD themselves

- AD §3.2 (T² with fixed dilaton, the same toy as DD §5.1), quoted: "We can furthermore take a₁ > 0, taking into account a₁ < 0 by multiplying this result by 2. The result is … N_vac(L) = 2σ(L)."
- They then match this to the continuum volume (3.11): "π²L/3, which agrees with the L derivative of the previous computation."
- DD §5.1 follows the same convention: their "240 vacua at τ = 2i" is #{A : 0 < 2abs(A)² ≤ 150} = 240, which counts A and −A. Their fig. 7 agrees with 2πAL² only in that convention (§1).
- The source of the factor is AD's introduction: "restricting z to a fundamental region … after which any two choices N ≠ N′ will lead to distinct flux vacua."
- **OUR ATTEMPT:** this sentence overlooks the central element −1 ∈ SL(2,ℤ). It fixes every τ and maps N → −N. Brute force (L = 10–80): #{fluxes with τ ∈ F} = 2Σσ + 2·#(orbits at i) + 4·#(orbits at ρ) exactly.
- So:
  - Phase 1's Σσ ~ π²L²/12 is the number of orbits of the full duality group.
  - The continuum π²L²/6 counts (τ ∈ F, N) pairs.
  - Neither is an arithmetic error.
- *Physical count.* If −1 is a gauge identification, the physical number is half the continuum formula. In AD's T² reading, −1 is the large diffeomorphism z → −z. In the IIB reading it is (−1)^(F_L)Ω, believed gauged (**CONJECTURE**, standard lore, not checked here).
- Generalization (**OUR ATTEMPT**, not verified on a CY): the continuum density must be divided by the order of the subgroup of the duality group that acts trivially on moduli space. That order is at least 2 for any IIB flux count.

## 4. (A3) Martel–Shapiro–Weinberg with Planck 2018

**Method** (MSW eqs. 4, 7, 9, 12–13, 21–22, 29, 33, 36; s = V/U = 1 unless stated).
- The only cosmological input is σ³ρ̄. In a matter-dominated era it is time-independent, and it equals ρ_m0·A³ with A = σ₀(R_G)/g₀.
- σ₀ is the linear rms today, from CAMB with the Planck 2018 parameters above (CDM + baryons, Gaussian filter as in MSW eq. 36).
- g₀ = D(1)/a_early = 0.78799 for Ω_m = 0.3153 (MSW's f = 1/g₀, eq. 48).
- ρ_obs/ρ_m0 = Ω_Λ/Ω_m = 2.1716.

**Results (OUR ATTEMPT).**

| R_G [Mpc] | M_f = (2π)^1.5 ρ_m0 R_G³ [M_⊙] | σ₀ | x* = ρ_obs/σ³ρ̄ | P(≤ρ_obs), s = 1 (s = 0, ∞) | median ρ_V/ρ_obs |
|---|---|---|---|---|---|
| 0.5 | 7.8e10 | 2.792 | 0.049 | 0.043 (0.082, 0.034) | 30.4 |
| 1.0 | 6.3e11 | 2.106 | 0.114 | **0.086** (0.140, 0.071) | 13.0 |
| 2.0 | 5.0e12 | 1.510 | 0.309 | 0.184 (0.254, 0.161) | 4.8 |
| 3.0 | 1.7e13 | 1.208 | 0.603 | 0.293 (0.368, 0.264) | 2.5 |

- MSW's own L*-galaxy filter formula, R_G ≈ h^(−1/3)(Ω_Bh²/0.015)^(−1/3) Mpc, gives 0.998 Mpc with Planck values. The 1 Mpc row is therefore the like-for-like update of their 5.4% (λ₀ = 0.7, h = 0.7, COBE + BBKS).
- The increase to 8.6% comes from a smaller σ: we find σ_rec(z = 1000) ≡ A/1001 = 0.00267, against MSW's 0.0033 (n_s < 1, lower normalization).
- Replacing CAMB by MSW's BBKS transfer (eq. 41) normalized to σ₈ changes P by ≤ 7% relative.
- The filter scale dominates the uncertainty: a factor 6 in R_G moves P by a factor 7.

**With the A1 prior.**
1. *Power law near 0*, p(ρ) ∝ ρ^α. The A1 fits bound α = p − 1 to ±0.01, which gives P = 0.083–0.089 (P(α = −0.05/+0.05) = 0.102/0.072). The anthropic result needs flatness only at the |α| ≲ 0.05 level, and BP gives far better than that.
2. *Discrete Poisson prior* with N_w vacua in [0, 10σ³ρ̄]:

   | N_w | mean P | 16–84% range | P(no vacuum below ρ_obs) |
   |---|---|---|---|
   | 100 | 0.081 | 0–0.157 | 33% |
   | 1000 | 0.086 | 0.063–0.108 | 0 |
   | 10⁵ | 0.0856 | 0.083–0.088 | 0 |

   Discreteness is irrelevant once N_w ≳ 10³.
3. *Comb prior* (near-commensurate charges, A1d) with spacing Δ: for Δ = 0.3 (1.0) σ³ρ̄, P ranges over 0–0.31 (0–0.63) depending on the phase. A gappy landscape makes the anthropic prediction depend on where the comb happens to fall.

## 5. (B1) Swampland quintessence: exponent vectors and minimal slopes

**Exponent vectors.**
- A 10d string-frame term weighted by e^(−χφ) (χ = 2 sphere, 1 disk/RP², 0 torus or RR) and scaling internally as ρ^w gives V_E ∝ ρ^(w−3χ/2) τ^(χ−4). Hence k = (−a√(2/3), −b/√2) for V ∝ ρ^a τ^b.
- Scaling sources: HKTT (3.3), (3.4), (4.1), (4.2), (4.4) as reproduced in §1. Casimir (ρ^(−2)τ^(−4)) and the 10d one-loop Λ (ρ³τ^(−4)) are our derivation (χ = 0).

| term | k = (k_ρ̂, k_τ̂) | abs k | sign |
|---|---|---|---|
| R₆ (geometric flux, KK monopole) | (0.816, 1.414) | 1.633 | + if internal curvature < 0 |
| H₃ | (2.449, 1.414) | 2.828 | + |
| NS5 | (1.633, 1.414) | 2.160 | + |
| F_p | ((p−3)·0.816, 2.828) | 2.83–3.74 | + |
| D_p / O_p | ((6−p)·0.408, 2.121) | 2.12–2.45 | + / − |
| Q, R non-geometric | (−0.816, 1.414), (−2.449, 1.414) | 1.633, 2.828 | ± |
| Casimir | (1.633, 2.828) | 3.266 | ± |

**All-positive potentials (OUR ATTEMPT).** −∇V/V = Σp_i k_i with p_i ≥ 0 and Σp_i = 1, so abs(∇V)/V ≥ dist(0, conv{k_i}) at every field point, not just asymptotically. We checked this on 2×10^5 random points.

| source set | minimal slope | realized by | ε = γ²/2 | along the attractor |
|---|---|---|---|---|
| IIB geometric (R₆, H₃, F₁,₃,₅, D3–D9, NS5) | √(18/7) = 1.6036 | R₆ + D9 | 9/7 | vol ↑, g_s ↓ (controlled) |
| IIA geometric (R₆, H₃, F₀–₆, D4–D8, NS5) | √(50/19) = 1.6222 | R₆ + F₀ (+D8, collinear) | 25/19 | vol ↑, g_s ↓ (controlled) |
| + positive Casimir / 10d one-loop Λ | unchanged | | | |
| + positive Q flux | √2 | R₆/H₃/NS5 + Q | 1 (w = −1/3) | fixed vol, g_s ↓ |

ODE integration of the STT autonomous system (2303.03418 eqs. A.1–A.2) reproduces ε = 1.3158, 1.2857 and 1.0000 from random initial data.

**Physical realizability.**
- IIA R₆ + F₀: Romans mass on a negatively curved six-manifold. With H = 0 there is no D6 tadpole, so this is plausibly realizable as a runaway.
- IIB R₆ + D9: positive net D9 tension needs anti-D9s. With O9s and tadpole cancellation, BPS D9 + O9 sum to zero. So it is formal in type IIB/I.
- The √2 case needs non-geometric Q flux: formal.

**Mixed signs (OUR ATTEMPT).**
- *Two-term solutions.* None of 5258 combinations has a stable two-term (+, −) scaling attractor with V > 0. Analytically, the potential transverse to the trajectory is Ae^(−as) − Be^(−bs), and wherever V > 0 its extremum is a maximum. This is STT's "at least three terms, not all positive".
- *Collinear (+, −, +) triples* stabilize the transverse (dilaton-type) direction:
  - The minimum, γ* = √(2/3) (ε = 1/3, accelerating), comes from R₆(<0) + O4 + F4, rolling to **strong 10d coupling** (d ln g_s/ds = +1.22), and from Q + O8 + F2, rolling to **zero volume**. Both leave the controlled regime.
  - Three triples lie on lines through the origin (R₆ + O3 + Casimir or F₅; F₁ + O9 + Q). They give formal dS critical points with an exactly flat direction. Here V depends only on M_KK/M_P ∝ (ρτ²)^(−1/2). But O3/O9 tadpoles cancelled by BPS branes cancel the tension too, and flux cancellation adds off-line terms. Formal only.
  - Among mixed triples moving to large volume and weak coupling, the minimum is √2 (non-geometric), and 1.622 with geometric sources (R₆ + O8 + F₀).

**Against the bounds and the data.**

| | value | distance from λ_DESI = 0.698 (+0.173, upper error) |
|---|---|---|
| TCC average, 2/√6 (1909.11063) | 0.816 | 0.7σ |
| strong asymptotic √2 (2101.11617) | 1.414 | 4.1σ |
| our all-positive IIB min | 1.604 | 5.2σ |
| our all-positive IIA min | 1.622 | 5.3σ |
| formal mixed min (uncontrolled) | 0.816 | 0.7σ |

- With Ω_K free, λ = 0.722 (+0.182), and the √2 distance is 3.8σ.
- The Gaussian reading of an asymmetric posterior is rough. The DESI λ is the slope traversed at z ≲ few. For a thawing field starting at rest it is ≈ abs(∇V)/V, so the pointwise bound above applies.

## 6. (B2) Literature: what is already known

**Anthropic Λ with modern inputs.**
- Sudoh et al. 1607.00180 (semi-analytic galaxy formation): P(<Λ_obs) = 6.7%, rising to 9.0–9.7% with metallicity weighting.
- Barnes et al. 1801.08781 (EAGLE): whether Λ is "predicted" "depends crucially on the measure".
- Totani et al. 1804.10395: supernova lethality pulls the mean toward Λ_obs.
- **Sorini–Peacock–Lombriser 2411.07301 (2024):** star-formation-efficiency weighting with a flat prior gives median 539 Λ_obs and **P(Λ ≤ Λ_obs) = 0.5%**.
- We found no 2024–2026 paper that redoes exactly the MSW spherical model with Planck 2018. Our 4–29% is therefore a like-for-like update, not the state of the art. The modern spread (0.5% to 10%) is driven by the observer weighting, not by the cosmological parameters.

**Swampland vs DESI.**
- 2504.04226: λ = 0.698. 2405.17396 (DR1 + supernovae): λ = 0.48–0.77. 2405.18747 (DR1): no preference over ΛCDM.
- 2405.09323: acceleration today needs λ ≲ √3; for λ > √2 the acceleration is transient.
- 2309.03938: an open universe gives a stable accelerating fixed point for γ > 2/√(d−2).
- 2209.11821: F-theory complex-structure limits with γ < √2 along the gradient flow, pending Kähler moduli; it also gives the convex-hull formulation.
- 2303.03418, 2306.07327: rolling-dilaton no-go and convex-hull criterion.
- 2208.08989: evidence for asymptotically decelerating cosmologies.
- 2509.15091: a momentum coupling to dark matter re-allows λ ≥ √2.
- 2510.21627: two-field exponentials emulate a shallower slope.
- 2507.03090: dark-dimension radion with a coupling to dark matter.
- Our B1 numbers are explicit instances of the STT/CRV statements, not a new principle.

## 7. What this says about a Λ-selection principle

1. **(A) rests on two premises; the computations settle one and expose the other.**
   - *Prior (settled within the toy).* A flat prior near 0 holds exactly for incommensurate BP charges at every resolvable scale. The anthropic number needs flatness only to |α| ≲ 0.05.
   - *Prior (what it hinges on).* Flatness fails as soon as charges repeat or are nearly commensurate. So (A) hinges on incommensurability, which no one has computed for a real stabilized landscape.
   - *Anthropic weight (the weak link).* In MSW's model with Planck 2018, ρ_obs sits at the 4–29% quantile, depending on the filter mass. Richer observer models give 0.5% (2411.07301).
   - So (A) is consistent with ρ_obs but predicts it only to about an order of magnitude, and the leading uncertainty is the observer model.
2. **(B), restricted to the universal dilaton–volume sector with standard sources, does not accelerate at all.**
   - All-positive potentials have abs(∇V)/V ≥ 1.60 everywhere, or ≥ √2 with non-geometric flux.
   - Mixed-sign constructions reach the TCC value √(2/3) only by running to strong coupling or zero volume.
   - The DESI-preferred λ ≈ 0.7 is 4–5σ from every controlled value we found.
   - (B) also leaves V₀, the size of Λ, unexplained. Within this sector it is not a Λ-selection principle. It is a no-go that pushes quintessence model-building into other moduli or couplings (2209.11821, 2507.03090, 2509.15091).
3. Neither principle, as computed here, selects ρ_obs. (A) survives as a statistical explanation conditional on an unverified incommensurability and an observer measure. (B) is, in its minimal form, in tension with the very data that motivated it.

## 8. Failures and gaps

- **A1:** this is a quadratic-form toy with J ≤ 24 and n_max ≤ 4800; there are no moduli and no tadpole. The near-commensurate study uses one J. Σ²(L) for L ≳ 10⁴ is noise-limited.
- **A2:** the "at least 2" generalization to real CY monodromy groups is not checked. Whether −1 is gauged is lore, not derived here.
- **A3:**
  - MSW's spherical model; ρ_V < 0 is ignored.
  - The choice of the s parameter and the filter dominate the result.
  - Radiation is neglected in g₀ (affects only late-time growth).
  - There is no baryonic physics.
  - The P values move by a factor 7 over plausible R_G.
- **B1:**
  - Only the dilaton–volume truncation; no axions, no α′ or g_s corrections. Tadpoles are treated qualitatively.
  - The existence of compact negatively curved manifolds with the needed cycles and O-planes is assumed.
  - The Casimir scaling is derived, not cited.
  - Collinear triples with three different k's are the only three-term case examined. Non-collinear four-or-more-term mixed sets are not.
- **B2:** Exa-based search, which may miss papers. Authors of 2405.17396 and 2305.18261 were not read.

## 9. Novelty assessment

- *Expected results confirmed rather than discoveries:*
  - The BP Poisson and flatness result (Berry–Tabor for diagonal forms).
  - The MSW update.
  - The B1 hull minima, which are explicit instances of STT/CRV.
- *Small new pieces, not found in our search:*
  - the exact face-corrected BP density (BP's D = 2^J estimate is off by up to a factor 2.2 at n_max ~ 30);
  - the explicit verdict on the AD/DD factor of 2;
  - the closed-form minima √(18/7) and √(50/19), with the pointwise (not only asymptotic) bound;
  - the observation that the only TCC-saturating mixed-sign attractors in this catalogue run to strong coupling or zero volume.
- None of these changes the status of either selection principle.
