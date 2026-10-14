# Phase 2 — Targets

Inputs: the five Phase 1 maps in `../phase1/` (A–E) and the sanity checks in `../phase3/sanity/`.

## Selection criteria
A target qualifies only if all four hold:
1. It bears directly on the mission: (1) a non-perturbative or consistency-based definition outside AdS, or (2) vacuum selection / Λ.
2. A published number exists that we can reproduce first.
3. The open part reduces to an explicit computation (sympy/numpy/scipy, hours not months) with a yes/no or numerical answer.
4. Any outcome is informative, including a negative one.

## T1 — Is the critical string the extremal corner of tree-level unitarity? (mission 1, Front E)
**Question.** Within natural families of dual-resonant amplitudes, does requiring every residue to be positive in partial waves single out the critical string (intercept and D) as a corner of the allowed region?

**Seed (OUR ATTEMPT, `phase3/targets/T1_exploration/`).**
- Family: Veneziano with intercept a, α(s) = a + s, and external scalars on the leading trajectory (m² = −a).
- Level 2 gives c₁ ∝ −(a−1)(3a+2) and c₀ ∝ D(a²−2a) + 8a² + 14a + 4.
- Numerically, up to level 20, the region is a ≤ 1 with D ≤ D_max(a), where D_max is increasing in a. Its unique maximal-D point is (a, D) = (1, 26).

**What to do.**
- (i) Make the (a, D) result rigorous or high-precision to large level, with large-n asymptotics.
- (ii) Run the same analysis for closed strings: the Virasoro–Shapiro λ-family of 2408.03362, and VS with a general intercept.
- (iii) Test the claim in 2102.08951 (footnote 14) that the VS residues are Gegenbauer-positive for D ≤ 24.148…, and find the exact critical D.
- (iv) Run the superstring analogue, where the expected critical D is 10.

**Why this target.**
- Front E concluded that every uniqueness theorem smuggles in one UV-shaped assumption.
- A statement of the form "unitarity alone puts the critical string at a corner" is checkable, and it is either new or traceable to a known source; either result is informative.
- Known baselines exist to reproduce: D ≤ 26 (2201.11575; 2608.25183), our sanity script S2, and the intercept fixed by 7-point factorization (2603.04485).

## T2 — Stress-test the explicit KKLT de Sitter candidates (mission 2, Front D)
**Question.** Are the five candidate dS vacua of 2406.13751 self-consistent by their own control criteria?

**What to do.**
- (a) **Heavy conifold.** Compare the conifold-modulus mass (BDGL 1809.06861 eq. 2.24) with the Kähler-modulus mass at the published parameters, where g_sM² ≈ 12–17 sits well below the 46.1 bound.
- (b) **KPV at g_sM ≈ 1.** Assess the anti-D3 (KPV) metastability with the known α′² corrections at g_sM ≈ 1, outside the fitted range of existing analyses.
- (c) **dS window vs flux discreteness.** Compare the dS window with the granularity of flux discreteness, and give an expected yield to set against 5/396.

**Baseline to reproduce first.**
- The KKLT numbers: σ_cr = 113.589, V_AdS = −1.997×10⁻¹⁵.
- The KPV (p/M)_crit = 0.0805.
- √g_s·M_min = 6.80.

**Why this target.**
- These are the only explicit dS candidates in string theory, and their parameters are published.
- The disputed points are numbers, not philosophy. A clear inconsistency, or a clear pass, is informative either way.

## T3 — Λ selection head-to-head: discretuum + anthropics vs swampland quintessence (mission 2, Front C)
**Question.** For each of the two leading "selection principles" for Λ, compute the quantity it actually hinges on.

**(A) Anthropic + landscape.**
- Weinberg-type arguments need the prior on Λ to be flat near 0.
- Measure the density and spacing statistics of Bousso–Polchinski vacua near Λ = 0 with incommensurate charges, and how they scale with J.
- Resolve the factor-of-2 between the exact and continuum flux counts found in Phase 1 for the rigid-CY toy.
- Compute the Martel–Shapiro–Weinberg probability P(ρ_Λ ≤ ρ_obs) with Planck-2018 inputs.

**(B) Swampland.**
- Compute the minimal late-time exponential slope obtainable from the standard type II source classes (convex hull of dilaton–volume exponent vectors).
- Compare it with the √2 and √(2/3) bounds and with the DESI DR2 fit λ ≈ 0.70.

**Why this target.**
- The mission asks for a selection principle. These two are the live candidates, and each rests on one checkable quantitative premise: a flat prior plus the anthropic weight, or a small slope.

## Considered and not chosen
- **A-H1** (three-loop D0-brane v⁸/r²¹): high value, but a three-loop matrix-QM computation is out of scope for this session.
- **A-H4** (small-N IKKT): mostly reproduction at N ≤ 8, with weak leverage on the Lorentzian question.
- **B-H1** (O(λ) DSSYK vs dS3): needs a one-loop sine-dilaton bulk propagator whose conventions are not fixed in the literature. Risk of a convention mismatch masquerading as a result.
- **B-H3** (celestial string discrepancy): sent to Phase 4 as a claim to verify, not a target.
- **E-H1/H3** (LP boundaries): good, but they overlap T1 and need a large linear-programming engineering effort.
