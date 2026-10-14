# T3: Λ selection, head to head. (A) Discretuum plus anthropics vs (B) swampland quintessence

**Labels.**
- **ESTABLISHED**: cited by arXiv ID and equation number. The abstract page or full text was fetched in this session; source texts are in `sources/`.
- **CONJECTURE**: a claim made in the literature and not proven.
- **OUR ATTEMPT**: computed here; the script is named.

All files are in `phase3/targets/T3/`. Slopes are in reduced Planck units.

## 0. Headline numbers

| Quantity | Result | Label / script |
|---|---|---|
| BP: N(Λ ∈ [0, ε]) ∝ ε^p, J = 4–12, about 10 decades of ε | p = 0.994–1.008. χ²/dof on 42 disjoint shells: 0.89–0.95 (J = 6–12), 1.40 (J = 4) | OUR ATTEMPT, `T3_A1_bp_discretuum.py`, `T3_A1_summary.py` |
| BP: spacings near Λ = 0 (6×10^5 levels per J) | Poisson: var s = 0.996–1.005, P(s<0.1) = 0.095, Σ²(L) ≈ L up to L ~ 10³ | OUR ATTEMPT |
| BP: min abs Λ, J = 4–24, 960 shells | z = 2ρ(0)·min abs Λ: mean 0.943 ± 0.030; KS p = 0.28 against Exp(1) | OUR ATTEMPT, `T3_A1c_min_statistics.py` |
| Rigid-CY factor of 2 | Convention, stated by AD: they count N and −N separately. Phase 1 counted modulo the full SL(2,ℤ) | ESTABLISHED (text) + OUR ATTEMPT, `T3_A2_rigid_cy_factor2.py` |
| MSW P(ρ_V ≤ ρ_obs), Planck 2018, R_G = 1 Mpc, s = 1 | **0.086** (0.071–0.140 over s). Median ρ_V/ρ_obs = 13 | OUR ATTEMPT, `T3_A3_msw_anthropic.py` |
| The same, R_G = 0.5–3 Mpc | P = 0.043–0.29; median 30 down to 2.5 | OUR ATTEMPT |
| Minimal slope, all-positive geometric type II sources (dilaton–volume plane) | IIB √(18/7) = 1.604, IIA √(50/19) = 1.622. √2 with positive non-geometric Q flux. Holds at every field point | OUR ATTEMPT, `T3_B1_swampland_slopes.py` |
| Minimal slope, mixed signs (formal) | √(2/3), but only while running to strong 10d coupling or zero volume | OUR ATTEMPT |
| DESI DR2 λ = 0.698 (+0.173/−0.202) | 4.1σ below √2; 5.2σ below 1.604 | ESTABLISHED (2504.04226 abstract) + OUR ATTEMPT |

## 1. Baseline reproductions

- **BP**, hep-th/0004134.
  - λ = λ_bare + ½Σn_i²q_i² (2.21); shell (2.22)–(2.24); spacing (2.26). The degeneracy is D = 2^J from n_i → −n_i (text after (2.29)).
  - Our MITM code equals brute force exactly for J = 3–5 (`T3_A1_mitm_selftest.py`).
  - (2.26) with D = 2^J reproduces the measured density at J = 4 (ratio 1.0002). It underestimates by 2%, 15%, 50% and 116% at J = 6, 8, 10, 12 (n_max = 280, 90, 34, 29).
  - Phase 1's q = 10^−1.73 for J = 100 (D = 2^J) is reproduced: 10^0.3 vacua with abs λ < 10⁻¹²⁰ (part e).
- **AD**, hep-th/0307049, §3.2 (T² toy). Reproduced: N_vac(L) = 2σ(L); Σσ ~ π²N²/12; the volume (3.11) (π²L/3); vol(F) = π/12 (2.23).
- **DD**, hep-th/0404116, §5.1 (arXiv text read this session).
  - "240 vacua at τ = 2i" (L = 150): we count 240 fluxes, i.e. 120 orbits.
  - Fig. 7 estimate 2πAL² (A measured with g = 1/(4y²)): at R = 0.15/0.20/0.25/0.30 it gives 627/1119/1756/2541. Our counts with N and −N separate are 632/1202/1734/2592.
- **MSW**, astro-ph/9701099 (ApJ text).
  - ⟨ρ_V⟩/σ³ρ̄ = 2.579/3.010/3.224 (eq. 25): exact.
  - The s = 1 quantiles 0.0590/1.49/11.1 σ³ρ̄ (text after eq. 30): we get 0.05897/1.485/11.11.
  - Eq. (29) equals direct integration of (9) to 10⁻⁶.
  - Table 1:
    - σ_rec from their eqs. (39)–(52): all 18 entries match to the 2 printed digits.
    - P(≤ρ*) agrees within the rounding of their σ, ≤ 6% (λ₀ = 0.7: 0.053/0.121 against 0.054/0.12).
- **HKTT**, 0711.2512.
  - Canonical ρ̂ = √(3/2) ln ρ, τ̂ = √2 ln τ (2.5). We verified this exactly by a sympy minisuperspace reduction of the 10d action.
  - Scalings (3.3), (3.4), (4.1), (4.2), (4.4) are reproduced by one rule, V_E ∝ ρ^(w−3χ/2) τ^(χ−4).
- **Planck 2018**, 1807.06209, Table 2, column TT,TE,EE+lowE+lensing.
  - Ω_bh² = 0.02237, Ω_ch² = 0.1200, H₀ = 67.36, n_s = 0.9649, ln(10¹⁰A_s) = 3.044, τ = 0.0544, Ω_m = 0.3153, Ω_Λ = 0.6847, σ₈ = 0.8111.
  - CAMB 2.0.4 with one 0.06 eV neutrino gives σ₈ = 0.8109.

## 2. (A1) Bousso–Polchinski discretuum near Λ = 0

**Setup (OUR ATTEMPT, `bp_lib.py`).**
- Λ = (Q − R²)/2, with Q = Σq_i²n_i² and R² = 2|λ_bare|.
- We count distinct Λ values, i.e. orthant points, which removes BP's D = 2^J.
- The charges are random q_i ∈ [1, 2). Taken as exact float64 numbers, they are rationally independent at these sizes.
- The shell is enumerated exactly with MITM in 80-bit arithmetic. Minima are re-verified in exact rational arithmetic.
- Flux cutoff means n_max = R/q_min.

**Face-corrected density (OUR ATTEMPT).**
- Exact identity: #{n ∈ ℤ^J_≥0 : Q ≤ X} = 2^(−J) Σ_{U⊆[J]} #(ℤ^U ∩ E_X). A point with support T is counted 2^|T| × 2^(J−|T|) times on the right.
- Volumes then give ρ(Λ) = 2^(−J) Σ_k ω_{k−1} e_k(1/q) X^((k−2)/2), with X = R² + 2Λ. BP's (2.26) is the k = J term alone.
- This matches every measured window density to ≤ 0.2%.

**(a) Is the prior flat?** We counted N(Λ ∈ [0, ε]) and N(Λ ∈ [−ε, 0]) for ε from about 3 mean spacings to 0.02R². That is 10–11 decades; for example, at J = 12 ε runs from 10⁻¹⁰ to 18, with R² = 900.

| J | R | n_max | fitted p (N ≥ 10, ε ≤ 10⁻³R²) | χ²/dof on 42 disjoint shells (p-value) | max abs pull |
|---|---|---|---|---|---|
| 4 | 1500 | 1458 | 1.0011 | 1.40 (0.045) | 2.8 |
| 6 | 330 | 280 | 0.9937 | 0.94 (0.59) | 1.9 |
| 8 | 90 | 90 | 1.0009 | 0.95 (0.57) | 2.2 |
| 10 | 46 | 34 | 1.0083 | 0.89 (0.67) | 2.5 |
| 12 | 30 | 29 | 0.9955 | 0.95 (0.56) | 2.4 |

- There are no gaps and no power-law deviation down to the single-level scale. Counts are Poisson around the smooth formula.
- The mild excess at J = 4 (|n/n_p − 1| ≤ 1.1%) is consistent with the larger lattice-point discrepancy expected in low dimension. We did not analyse it further.

**(b) Spacing statistics** (about 6×10^5 consecutive levels centred on Λ = 0, unfolded with ρ).

| J | var(s) | P(s<0.1) | P(s<0.01) | KS p vs e^(−s) | Σ²(100) | Σ²(1000) |
|---|---|---|---|---|---|---|
| 4 | 1.000 | 0.0948 | 0.0101 | 0.65 | 98 | 932 |
| 6 | 1.002 | 0.0949 | 0.0100 | 0.68 | 96 | 944 |
| 8 | 1.005 | 0.0952 | 0.0099 | 0.65 | 102 | 1030 |
| 10 | 0.996 | 0.0954 | 0.0099 | 0.61 | 98 | 966 |
| 12 | 1.004 | 0.0954 | 0.0102 | 0.06 | 98 | 961 |

- Reference values: Poisson gives 1 / 0.0952 / 0.00995; GOE gives var 0.286 and P(s<0.1) = 0.0078.
- These are Poisson statistics: no level repulsion, no rigidity. Σ²(L) for L ≳ 10⁴ is too noisy to test Berry saturation.
- **ESTABLISHED context:**
  - Poisson pair correlation is proven for inhomogeneous forms under Diophantine conditions (Marklof, math/0210198).
  - It is proven for 2D flat tori; see the summary in 2305.18261.
  - We found no theorem for the BP case: homogeneous diagonal forms, J ≥ 3, with sign identification.

**(c) Minimum abs Λ** (`T3_A1c_min_statistics.py`). Poisson prediction: z = 2ρ(0)·min abs Λ ~ Exp(1).

| J | n_max | mean z per R (s.e. 0.10–0.18) | R-exponent of ⟨min abs Λ⟩: measured / smooth ρ / −(J−2) |
|---|---|---|---|
| 4 | 277–2492 | 0.92, 0.85, 1.00 | −1.96 / −2.00 / −2 |
| 6 | 35–138 | 0.80, 1.34, 1.08 | −3.70 / −3.92 / −4 |
| 8 | 25–64 | 0.86, 0.78, 1.02 | −5.53 / −5.71 / −6 |
| 10 | 16–36 | 0.98, 1.01, 0.80 | −7.58 / −7.35 / −8 |
| 12 | 11–20 | 0.90, 0.96, 1.06 | −8.49 / −8.77 / −10 |
| 16 | 7–12 | 0.95, 0.73, 0.94 | −11.2 / −11.2 / −14 |
| 20 | 7–10 | 1.04, 0.93, 0.93 | −13.9 / −13.8 / −18 |
| 24 | 6–9 | 0.90, 0.90, 0.98 | −16.0 / −16.3 / −22 |

- Pooled over all J = 4–24 (960 shells): mean z = 0.943 ± 0.030, var 0.85, KS p = 0.28.
- So min abs Λ ≈ 1/(2ρ(0)) at every J. The BP exponent −(J−2) is reached only for n_max ≫ J. For J = 12–24 at n_max ~ 10, the effective exponent is 0.7–0.85 times (J−2), because faces dominate.
- At fixed flux cutoff n_max ≤ 11 (part c′), exact min abs Λ falls from 2.0×10⁻³ (J = 4) to 1.15×10⁻¹³ (J = 24). Prediction 1/(2ρ) = 1.15×10⁻¹³.
- For BP's own J = 100 example the face terms raise the distinct-vacuum density by 10^2.2–10^4.2. For J = 300, q = 10^−0.87 the gain is 10^60.7, since there n_max = 10.5 < √J, the breakdown regime BP flagged (`T3_A1e_bp_faces.py`). This is an extrapolation beyond the tested J ≤ 24.

**(d) Sensitivity to the charges.**
- *Near-commensurate charges*, q_i² = 1 + δu_i (J = 8, R² = 400.37):
  - δR² = 120 and 12: Poisson (P(s<0.1) = 0.095, largest gap 16 spacings ≈ ln N).
  - δR² = 1.2: the level count near 0 is 2% of the flat expectation, and the largest gap is 2×10⁶ mean spacings. min abs Λ = 1.5×10⁻⁶, against a Poisson expectation of 4×10⁻⁸.
  - δR² ≤ 0.12: min abs Λ ≈ 0.14–0.19, i.e. dist(R², ℤ)/2.
  - δ = 0: 4 distinct values in abs Λ < 1.
- *Repeated charges* (J = 12, three distinct charges used 4 times each):
  - The number of distinct Λ values follows J_eff = 6 (prediction X²/Π q_a²): 315/1059/3578 against a predicted 283/1007/3475 at R = 8/11/15.
  - The mean degeneracy is 860–13 000.
  - At R = 15 the distinct density is 4×10³ per unit Λ, against 4.8×10⁷ for generic J = 12.

**What A1 does and does not say.**
- *Does:* for a fixed positive-definite form with rationally independent charges, Λ near 0 is a locally Poisson, gap-free sequence. Its density is known analytically including faces. "Flat prior" holds at every resolvable scale.
- *Does:* flatness is a consequence of incommensurability, not of the vacuum count. Degenerate or near-commensurate charges produce combs, gaps and a reduced effective J.
- *Does not:* real landscapes are not fixed quadratic forms. There Λ depends on stabilized moduli; tadpoles, not a shell, bound the fluxes; dualities add degeneracies (BP mention E₇₍₇₎).
- None of these effects was tested here.

## 3. (A2) The rigid-CY factor of 2. Verdict: convention, documented by AD

- AD §3.2 (the same toy as DD §5.1), quoted: "We can furthermore take a₁ > 0, taking into account a₁ < 0 by multiplying this result by 2. The result is … N_vac(L) = 2σ(L)." They match this to the continuum (3.11): "π²L/3, which agrees with the L derivative of the previous computation."
- DD follow the same convention. Their 240 vacua at τ = 2i is #{A : 0 < 2abs(A)² ≤ 150}, which counts A and −A. Their fig. 7 agrees with 2πAL² only in that convention.
- The root is AD's intro: "restricting z to a fundamental region …, after which any two choices N ≠ N′ will lead to distinct flux vacua." That ignores the central element −1 ∈ SL(2,ℤ), which fixes every τ and maps N → −N.
- **OUR ATTEMPT (brute force, L = 10–80):** #{fluxes with τ ∈ F} = 2Σσ + 2·#(orbits at i) + 4·#(orbits at ρ), exactly.
- So:
  - Phase 1's Σσ ~ π²L²/12 counts orbits of the full duality group.
  - π²L²/6 counts (τ ∈ F, N) pairs.
  - Neither contains an error.
- *Physical count.* If −1 is gauged, the physical count is half the continuum formula. On the T² it is the diffeomorphism z → −z; in IIB it is (−1)^(F_L)Ω (**CONJECTURE**, standard lore, not checked).
- Generalization (**OUR ATTEMPT**, unverified for real CY monodromy): divide the continuum density by the order of the subgroup of dualities that acts trivially on moduli.

## 4. (A3) Martel–Shapiro–Weinberg with Planck 2018

**Method** (MSW eqs. 4, 7, 9, 12–13, 21–22, 29, 33, 36; s = 1 unless stated).
- σ³ρ̄ = ρ_m0(σ₀/g₀)³ (MSW 52 with f = 1/g₀, eq. 48).
- σ₀ is the CAMB linear rms today: CDM + baryons, with MSW's Gaussian filter R_G.
- g₀ = 0.78799 and ρ_obs/ρ_m0 = 2.1716.

| R_G [Mpc] | M_f [M_⊙] | σ₀ | ρ_obs/σ³ρ̄ | P(≤ρ_obs), s = 1 (s = 0, ∞) | median ρ_V/ρ_obs |
|---|---|---|---|---|---|
| 0.5 | 7.8e10 | 2.792 | 0.049 | 0.043 (0.082, 0.034) | 30.4 |
| 1.0 | 6.3e11 | 2.106 | 0.114 | **0.086** (0.140, 0.071) | 13.0 |
| 2.0 | 5.0e12 | 1.510 | 0.309 | 0.184 (0.254, 0.161) | 4.8 |
| 3.0 | 1.7e13 | 1.208 | 0.603 | 0.293 (0.368, 0.264) | 2.5 |

**Reading the table (OUR ATTEMPT).**
- MSW's L*-galaxy filter formula, evaluated with Planck values, gives R_G = 0.998 Mpc. The 1 Mpc row is therefore the like-for-like update of their 5.4%.
- The rise to 8.6% comes from a smaller σ_rec = 0.00267, against their 0.0033.
- MSW's BBKS transfer (eq. 41) normalized to σ₈ changes P by ≤ 7%.
- The filter scale dominates: a factor 6 in R_G moves P by a factor 7.

**Prior from A1.**
1. *Power law*, ρ^α with α = p − 1 within ±0.01 from A1a: P = 0.083–0.089. For α = −0.05/+0.05, P = 0.102/0.072. Weinberg-type reasoning needs |α| ≲ 0.05.
2. *Poisson-discrete prior* with N_w vacua in [0, 10σ³ρ̄]:

   | N_w | mean P | 16–84% range | P(no vacuum below ρ_obs) |
   |---|---|---|---|
   | 100 | 0.081 | 0–0.157 | 33% |
   | 1000 | 0.086 | 0.063–0.108 | 0 |
   | 10⁵ | 0.0856 | 0.083–0.088 | 0 |

3. *Comb prior* (A1d) with spacing Δ = 0.3 (1.0) σ³ρ̄: P ranges over 0–0.31 (0–0.63) depending on where the comb falls.

## 5. (B1) Swampland quintessence: exponent vectors and minimal slopes

**Exponent vectors.**
- A string-frame term weighted by e^(−χφ) (χ = 2 sphere, 1 disk/RP², 0 torus or RR) and scaling internally as ρ^w gives V_E ∝ ρ^(w−3χ/2) τ^(χ−4). Hence k = (−a√(2/3), −b/√2) for V ∝ ρ^a τ^b.
- Every τ̂-component is (4−χ)/√2 ≥ √2. This is STT 2303.03418 eq. (II.3) at d = 4.

| term (source of scaling) | k = (k_ρ̂, k_τ̂) | abs k | sign |
|---|---|---|---|
| R₆ / geometric flux (HKTT 4.2) | (0.816, 1.414) | 1.633 | + if internal curvature < 0 |
| H₃ (3.3) | (2.449, 1.414) | 2.828 | + |
| NS5 (4.4) | (1.633, 1.414) | 2.160 | + |
| F_p (3.3) | ((p−3)·0.816, 2.828) | 2.83–3.74 | + |
| D_p / O_p (4.1) | ((6−p)·0.408, 2.121) | 2.12–2.45 | + / − |
| Q, R non-geometric (4.2) | (−0.816, 1.414), (−2.449, 1.414) | 1.633, 2.828 | ± |
| Casimir (our χ = 0, L⁻¹⁰ derivation) | (1.633, 2.828) | 3.266 | ± |

**All-positive potentials (OUR ATTEMPT).** −∇V/V = Σp_i k_i lies in conv{k_i}. Hence abs(∇V)/V ≥ dist(0, conv{k_i}) at every point; we checked this on 2×10^5 random points.

| sources | min slope | edge | ε = γ²/2 | attractor |
|---|---|---|---|---|
| IIB geometric (R₆, H₃, F₁,₃,₅, D3–D9, NS5) | √(18/7) = 1.6036 | R₆ + D9 | 9/7 | vol ↑, g_s ↓ |
| IIA geometric (R₆, H₃, F₀–₆, D4–D8, NS5) | √(50/19) = 1.6222 | R₆ + F₀ (D8 collinear) | 25/19 | vol ↑, g_s ↓ |
| + positive Casimir or 10d one-loop Λ | unchanged | | | |
| + positive Q flux | √2 | R₆/H₃/NS5 + Q | 1 (w = −1/3) | fixed vol, g_s ↓ |

ODE integration of STT eqs. (A.1)–(A.2) from random data gives ε = 1.3158, 1.2857 and 1.0000.

**Realizability.**
- IIA R₆ + F₀ (Romans mass on a negatively curved space, H = 0, so no D6 tadpole) is a plausible runaway.
- Net positive D9 tension needs anti-D9s; BPS D9 + O9 cancel. So it is formal.
- Q flux is formal.

**Mixed signs (OUR ATTEMPT).**
- *Two-term solutions.* No stable two-term (+, −) scaling attractor with V > 0 exists among 5258 combinations. The transverse potential Ae^(−as) − Be^(−bs) has a maximum wherever V > 0, as STT's three-term requirement implies.
- *Collinear (+, −, +) triples* stabilize the transverse direction:
  - γ* = √(2/3) (ε = 1/3, accelerating) occurs for R₆(<0) + O4 + F4, which runs to **strong 10d coupling** (d ln g_s/ds = +1.22), and for Q + O8 + F2, which runs to **zero volume**.
  - Lines through the origin (R₆ + O3 + Casimir or F₅; F₁ + O9 + Q) give formal dS points with an exactly flat direction. There V depends only on ρτ² ∝ (M_P/M_KK)². But BPS tadpole cancellation removes the O-plane tension, and flux cancellation adds off-line terms.
  - Restricted to large volume and weak coupling, the mixed minimum is √2 (non-geometric) and 1.622 with geometric sources.

**Against the bounds and the data** (one-sided, using the DESI upper error 0.173):

| | slope | distance from λ = 0.698 |
|---|---|---|
| TCC average 2/√6 (1909.11063) | 0.816 | 0.7σ |
| asymptotic 2/√(d−2) = √2 (2101.11617) | 1.414 | 4.1σ |
| our all-positive minima IIB / IIA | 1.604 / 1.622 | 5.2σ / 5.3σ |
| formal mixed minimum (uncontrolled) | 0.816 | 0.7σ |

- With Ω_K free (λ = 0.722 +0.182), the √2 distance is 3.8σ.
- A thawing field starts at rest, so its effective λ ≈ abs(∇V)/V. The pointwise bound therefore applies to the DESI-era slope as well as to the asymptotic one.

## 6. (B2) What is already in the literature

**Anthropic Λ.**
- Sudoh et al. 1607.00180 (semi-analytic galaxy formation): P(<Λ_obs) = 6.7–9.7%.
- Barnes et al. 1801.08781 (EAGLE): the outcome "depends crucially on the measure".
- Totani et al. 1804.10395: supernova lethality.
- **Sorini–Peacock–Lombriser 2411.07301:** star-formation weighting with a flat prior gives median 539 Λ_obs and **P = 0.5%**.
- We found no 2024–2026 redo of the exact MSW model with Planck 2018. Ours is a like-for-like update. The 0.5–30% spread across models comes from the observer weighting, not from the cosmological parameters.

**DESI vs swampland.**
- Fits: 2504.04226 (λ = 0.698); 2405.17396 (DR1 + supernovae, λ = 0.48–0.77); 2405.18747 (DR1, no preference).
- 2405.09323: λ ≲ √3 is needed for acceleration today.
- 2309.03938: an open universe gives a stable accelerating fixed point for γ > 2/√(d−2).
- 2209.11821: F-theory limits with γ < √2 along the gradient flow, pending Kähler moduli; convex-hull dS conjecture.
- 2303.03418 and 2306.07327: no-go and convex-hull criterion.
- 2208.08989: asymptotically decelerating cosmologies.
- 2509.15091: a momentum coupling to dark matter re-allows λ ≥ √2.
- 2510.21627: two fields emulate a shallower slope.
- 2507.03090: dark-dimension radion coupled to dark matter.
- Our B1 minima are explicit instances of STT/CRV, not a new principle.

## 7. What this says about a Λ-selection principle

1. **(A) hinges on two premises. We settled one in a toy and found the other dominant.**
   - *Flat prior.* It holds exactly in BP with incommensurate charges, at every resolvable scale, and only |α| ≲ 0.05 is needed.
   - *But flatness is inherited from incommensurability.* Repeated or near-commensurate charges destroy it. Whether real stabilized landscapes are "incommensurate enough" has not been computed by anyone, including us.
   - *Anthropic weight.* With MSW and Planck 2018, ρ_obs sits at the 4–29% quantile depending on the filter mass. Modern observer models give 0.5%.
   - So (A) is *consistent* with ρ_obs to about an order of magnitude and is *not predictive* beyond that. The observer model, not the landscape, is the leading uncertainty.
2. **(B), in the universal dilaton–volume sector with standard type II sources, cannot accelerate.**
   - All-positive potentials have abs(∇V)/V ≥ 1.60 everywhere, or ≥ √2 with non-geometric flux.
   - The mixed-sign route to √(2/3) exits perturbative control.
   - The DESI-preferred λ ≈ 0.7 is 4–5σ from every controlled value found.
   - (B) also leaves V₀ (the size of Λ) free. Here it is a no-go that pushes model-building to other moduli or couplings, not a selection principle.
3. Neither principle, as computed, selects ρ_obs. (A) is a conditional statistical explanation. Minimal (B) is in tension with the data that motivated it.

## 8. Failures and gaps

- **A1:**
  - A quadratic-form toy only: J ≤ 24, n_max ≤ 4800, no moduli or tadpoles.
  - The J = 100/300 face corrections are extrapolations.
  - A first run was OOM-killed at J = 24, R = 11.5 (`*_firstrun_OOMkilled.out`). The rerun used R ≤ 10.75 and reproduced all earlier numbers exactly.
  - The ε windows of part (c), with only 6 shifts each, are noisy; A1c supersedes them.
- **A2:** the "≥ 2" generalization to real monodromy groups was not checked.
- **A3:**
  - The spherical model; ρ_V < 0 is ignored.
  - The s and R_G choices dominate the result.
  - Radiation is ignored in g₀; there is no baryonic physics.
- **B1:**
  - A two-field truncation; no axions, no α′/g_s corrections; tadpoles treated qualitatively.
  - Negatively curved compact spaces with the needed cycles and O-planes are assumed.
  - The Casimir scaling is our derivation.
  - Only collinear three-term mixed sets were analysed.
- **B2:** Exa-based search, which may miss papers.

## 9. Novelty

- *Expected results confirmed:* BP Poisson flatness (Berry–Tabor-like), the MSW update, and the STT/CRV-type hull minima.
- *Small new pieces, not found in our search:*
  - the exact face-corrected BP density (BP's D = 2^J estimate is off by up to 2.2× at n_max ~ 30, and by 10^2–10^60 in their large-J examples);
  - the documented resolution of the AD/DD factor of 2;
  - the closed forms √(18/7) and √(50/19), with the pointwise bound;
  - TCC-saturating mixed attractors that necessarily leave control.
- None changes the status of either principle.
