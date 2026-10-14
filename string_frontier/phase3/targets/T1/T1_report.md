# T1 — Is the critical string the extremal corner of tree-level unitarity?

Phase 3 compute report. All code and outputs are in `phase3/targets/T1/`. Every number below was either computed by a script listed here or read from a fetched paper. Labels: **ESTABLISHED** (with arXiv ID and eq./fig. where read), **OUR ATTEMPT** (with the script), **CONJECTURE**.

## 0. Setup and conventions

**Families (OUR notation).**
- Open "intercept family": A_{a,c} = Γ(−α_s)Γ(−α_t)/Γ(c−α_s−α_t), with α(x) = a + x and α′ = 1.
  - c = 0: a = 1 is the bosonic string.
  - c = 1: a = 0 is the type‑I stripped function.
- Residues at α_s = n (s_n = n − a): R_n(t) = (α_t+1)_{n−c}/n!. The overall sign (−1)^c is dropped.
- External mass: m² = −a + δ. Setting δ = 0 means the externals are the level‑0 state ("on trajectory").
- Literature dictionary: this is the r = 0 slice of the hypergeometric amplitude in 2409.09561 eqs. (2.2) and (2.11), with m²_there = −a. It is also the q → 1 limit of the Coon amplitude in 2201.12331 eq. (4).

**Partial waves.**
- Variables: t = κ(x−1), with κ = (s_n − 4m²)/2 and x = cos θ.
- Expansion: R_n = Σ_j c_j C_j^{(D−3)/2}(x).
- Unitarity condition: **ρ_j ≡ c_j/κ^j ≥ 0**.
  - Above threshold (κ > 0) this is the usual c_j ≥ 0.
  - Below threshold, (|p||p′|)^j C_j is the polynomial that multiplies |g_j|².
  - Check: a massive vector exchanged between equal‑mass scalars has residue ∝ g²(t+κ) ∝ g²κx. So c_1 ∝ κ, while ρ_1 ∝ g² > 0.
- **Normalization check:** our coefficients reproduce 2608.25183 eq. (5) and (S18), and 2201.11575 Table 1, exactly. For example, B_{5,0} at D = 7 is 307/92160 in both (`T1_0_baseline_output.txt`).
- **Monotonicity in D (ESTABLISHED, 2608.25183 eq. (23), Gegenbauer connection formula):** positivity at D implies positivity at every D′ ∈ (3, D]. So every allowed set has the form D ≤ D_crit.

**Method (`pwcore.py`).**
- Exact rational arithmetic (python‑flint).
- For each (n, j), sign(ρ_j) = sign P_{n,j}(ν), where ν = (D−3)/2 and P_{n,j} is an explicit polynomial in ν.
- Real roots are isolated rigorously with arb.
- D_crit(n) is the first ν > 0 at which some P_{n,j} changes sign.

## 1. Baselines (ESTABLISHED results, reproduced)

**Bosonic Veneziano, a = 1.** Our level n is level n−1 in 2201.11575 and 2608.25183 (`T1_0_baseline.py`).
- All spins of all levels n ≤ 150 are ≥ 0 at D = 26, checked in exact arithmetic.
- The only exact zeros are (n, j) = (2, 0) and (3, 1). These are B_{1,0} and B_{2,1} of 2608.25183, which proves positivity for all levels in 3 < D ≤ 26 (eq. (4)).
- At D = 26 + 10⁻³, the only negatives are at levels 2 and 3.

**Type‑I stripped, c = 1, a = 0.**
- All levels n ≤ 150 are ≥ 0 at D = 10. The only zero is (3, 0), matching 2201.11575 eq. (1.12): B_{3,0} ∝ 10 − D.
- 2201.11575 proves the result for D ≤ 6. 2502.20372 proves it for D ≤ 10.

Per‑level D_crit(n) (rigorous root isolation; spin in brackets):

| n | 2 | 3 | 4 | 5 | 6 | 8 | 10 | 20 | 40 | 100 | 150 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| bosonic | 26 (0) | 26 (1) | 26.4 (2) | 27 (1) | 26.589 | 27.156 | 27.108 | 27.588 | 27.862 | 28.869 | 29.755 |
| type I | ∞ | 10 (0) | 11 (1) | 10.2357 (0) | 10.904 | 11 | 11.133 | 11.770 | 12.634 | 13.994 | 14.650 |

The binding constraint is always a low level. D_crit(n) grows slowly with n, consistent with the fixed‑j large‑n positivity of 2201.11575 eqs. (4.17) and (4.22) (ESTABLISHED).

## 2. The (a, D) region of the open c = 0 family, on trajectory

### Exact low levels (OUR ATTEMPT, `T1_1_lowlevels_symbolic.py`)
- Level 1: ρ_0 = (1−a)/2 and ρ_1 = 1/(D−3).
- Level 2: ρ_0 = [D(a²−2a) + 8a² + 14a + 4]/(8(D−1)), ρ_1 = (1−a)/(2(D−3)), ρ_2 > 0.
- Level 3: ρ_0 = (1−a)(1+a)[(26a+30) − D(3−a)]/(48(D−1)) and ρ_1 = [D(3a²−6a−1) + 30a² + 48a + 26]/(24(D−3)(D+1)).
- Levels 4 and 5 are in the output file.

### Region (OUR ATTEMPT)
Scripts: `T1_1_region.py`, `T1_1_boundary_exact.py`, `T1_1_highlevel_check.py`, `T1_1_gapcheck.py`.

The allowed set is {(a, D): −21/29 < a ≤ 1, 3 < D ≤ D_max(a)}. D_max is continuous and strictly increasing, and it is the spin‑0 boundary of four pieces:

| a‑interval | binding level, spin | D_max(a) |
|---|---|---|
| (−21/29, a₁] | 3, 0 | (26a+30)/(3−a) |
| [a₁, a₂] | 5, 0 | smallest root D > 3 of D²(−a⁴+4a³+14a²−36a−45) + D(−90a⁴−120a³+800a²+2000a+1250) − (1304a⁴+8224a³+19464a²+20464a+8080) |
| [a₂, a₃] | 4, 0 | smallest root of D²(a⁴−4a³−4a²+16a) + D(54a⁴+36a³−228a²−288a−64) + 296a⁴+1336a³+2368a²+2000a+704 |
| [a₃, 1] | 2, 0 | (8a²+14a+4)/(a(2−a)) |

Transition points (exact resultants):
- **a₁ = 0.3883686782…**, D = 15.35346. a₁ is a root of 45a⁵ + 78a⁴ + 20a³ − 64a² − 61a + 30.
- **a₂ = 0.5157634056…**, D = 17.40519. a₂ is a root of a degree‑15 resultant factor.
- **a₃ = 0.5221547735…**, D = 17.48343. a₃ is a root of 18a⁶ + 51a⁵ − 25a⁴ − 104a³ + 32a² + 72a − 32.

Consequences:
- a > 1 is excluded at level 1.
- D_max = 4 at a = −3/5, and D_max → 3 at a = −21/29.
- **The maximum of D over the region is 26, attained only at (a, D) = (1, 26).**

At the corner, the level‑1 scalar decouples (ρ_0 = 0) and level 1 becomes a massless vector. D = 26 is the level‑2 scalar bound.

**Verification.**
- **All levels n ≤ 100:** rigorous per‑level roots plus an exact check that every spin is ≥ 0 at D = D_max(a) (or a rational just below it). Done at 12 values of a: −7/10, −3/5, −1/2, −1/4, 0, 1/4, 2/5, 1/2, 3/4, 9/10, 99/100, 1.
- **Grids:**
  - step 1/100 with n ≤ 30;
  - step 1/1000 on [0.3, 0.6] with n ≤ 25;
  - step 1/5000 with n ≤ 8.
- **Margin of the higher levels:** min over 6 ≤ n ≤ 12 of D_crit(n; a), minus D_max(a), is at least 0.0209. This was checked on a 1/1000 grid and a 1/20000 grid near a = 0.516.
- **This is not a proof for all levels.**

| a | −0.7 | −0.5 | −0.25 | 0 | 0.25 | 0.4 | 0.5 | 0.52 | 0.75 | 0.9 | 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D_max | 3.189 | 4.857 | 7.231 | **10** | 13.273 | 15.532 | 17.138 | 17.457 | 20.267 | 23.313 | **26** |
| min_{n≥6} D_crit (n ≤ 100) | 4.018 | 5.607 | 7.873 | 10.489 | 13.551 | 15.674 | 17.215 | – | 20.956 | 23.987 | 26.589 |
| D_crit(n=100) | 5.409 | 7.763 | 10.813 | 13.994 | 17.308 | 19.339 | 20.691 | – | 24.074 | 26.542 | 28.869 |

### Large n (OUR ATTEMPT)
- For every sampled a, D_crit(n; a) trends upward with n beyond n ≈ 7, up to n = 100. The growth is not monotone: at a = 1, D_crit is 27.156 at n = 8 and 26.935 at n = 9.
- For every daughter trajectory Δ ≤ 12, the minimum over n ≥ 6 occurs at n ≤ 14 and lies above D_max.
- So D_max(a) is fixed by levels 2–5, and no higher level cuts the region further in the data.
- Heuristic for fixed j as n → ∞:
  - the forward peak near x = 1 scales as n^a/Γ(a+1);
  - the backward peak near x = −1 scales as (−1)^{n+j} n^{2a−1}/Γ(2a);
  - so for a < 1 the forward peak dominates and fixed‑j coefficients are positive in every D.
- This is not proven. 2409.09561 attributes fixed‑j asymptotics at r = 0 for general m² to 2406.04410, which I did not read.

### Sub‑threshold correction
- The seed's statement "no region for a < −1/3" was a **sign‑convention artifact**. For a < −1/3 the level‑1 vector lies below 4m², and c_1 < 0 while ρ_1 > 0.
- 2201.12331 and 2302.12263 stop at m² = 1/3 (a = −1/3) for exactly this reason: 2201.12331 eq. (13) gives c_{1,1} ∝ 1 − 3m².
- 2409.09561 instead imposes positivity only for n ≥ 3m² (its footnote 2). Under the ρ rule, sub‑threshold levels never bind for c = 0, so that convention gives the same region as ours.
- At m² = 1/3 we get exactly D_max = 32/5 = 6.4, set by level 3, spin 0. This was checked exactly for all levels ≤ 60 (`T1_1_highlevel_check_c0_N60_output.txt`). 2201.12331 quotes "≈ 6.3" with a stated 4% numerical accuracy, so the two are consistent.

**Exact lemma (OUR ATTEMPT).** If R(x) has definite parity, then (1+x)R(x) is Gegenbauer‑positive in D if and only if R is.
- Proof: the parity sectors separate. The coefficients of x·R are nonnegative combinations of R's coefficients, by the three‑term recurrence.
- At s = n with massless externals, Γ(−s)Γ(−t)/Γ(−s−t) = u·F_I with u ∝ (1+x).
- **Corollary:** the a = 0 member of the c = 0 family and the type‑I stripped function have identical D_crit(n) at every level. Numerically, both give 10.489178 at n = 7 and 13.9938 at n = 100.
- This explains the "m² = 0 ↔ Neveu–Schwarz, d = 10" identification in 2201.12331.

## 3. External mass off the trajectory: the corner is not robust

These are numerical scans (OUR ATTEMPT): `T1_1_massscan*.py` and `T1_1_offtraj_a1.py`, with levels ≤ 30–50.

**a = 1, m² = −1 + δ.**
- D_max(δ) dips below 26 (25.000 at δ = 1/20) and returns to exactly 26 at δ = 1/4.
- It exceeds 26 first at **δ\* = 0.4022459 (m² ≈ −0.5978)**, where the binding constraint is level 4, spin 0.
- It touches 26 again at δ = 1/2 and then grows: 26.85 at m² = 0 and 27.18 at m² = +1/2 (level 11 binding, n ≤ 50).
- **So (1, 26) is not a maximum of D over the (a, m²) plane, even at fixed a = 1.**

**Level‑1 edge, m² = −(1+3a)/4 (ρ_0 at level 1 = 0).**
- For a ≤ 1, D_max = (2a+3)² + 1. This comes from level 2: the residue is ∝ x² − 1/(2a+3)². For a = 1 it gives 26.
- For a > 1 the edge continues upward (44.86 at a = 2), but there the level‑1 vector is a spinning tachyon (mass² 1 − a < 0).

**Externals no heavier than the level‑0 state (m² ≤ −a).**
- On a grid (a step 1/16, δ step 1/32, levels ≤ 30), the maximum is **26, attained only at (1, −1)**.

**CONJECTURE.** (1, 26) is the unique maximum of D over {(a, m²): m² ≤ −a}.

## 4. Superstring analogue: the c = 1 family (OUR ATTEMPT)

Scripts: `T1_1_region.py` (c = 1) and `T1_1_highlevel_check.py` (c = 1).

**Region.**
- Level 2 spin 0 gives ρ_0 = −a/4, so a ≤ 0.
- On a ∈ (−7/11, 0], D_max(a) = (8a+10)/(1−a). It is set by level 3 spin 0 everywhere on the a‑grid. Verified for n ≤ 100 at a = −3/5, −1/2, −1/4, −1/10, 0, with min_{n≥6} D_crit above D_max by at least 0.48.
- **Corner at (a, D) = (0, 10)**, the type‑I point.

**Comparison with the literature.**
- The a ≤ 0 bound is the r = 0 case of Theorem 1 of 2409.09561 (ESTABLISHED).
- Why the stripped function gives 10: its level‑3 residue is ∝ 9x² − 1, so D − 1 ≤ 9.
- The full gluon amplitude gives 10 independently at level 1, from the F⁴ factor: the coefficient 9/(D−1) − 1 in 2201.11575 eq. (1.6) (ESTABLISHED).

**CHR planar λ‑family (2408.03362 eq. (18)).**
- D_crit(λ) = (λ+2)² + 1 (level 3, spin 0) for λ ≤ 3. Larger λ is set by higher levels, and the λ ≥ 8 values are upper bounds only because they are not converged at n ≤ 30.
- D_crit(λ) is strictly increasing: 5.25, 7.25, 10, 13.25, 17, 26 at λ = 1/16, 1/2, 1, 3/2, 2, 3.
- **No corner at λ = 1.**
- The residues of this family coincide with those of the c = 1 family at a = (λ−1)/2 with externals of mass m² = −3a/4. On this line the level‑2 scalar decouples identically (ρ_0 = 0).
- The equivalence is ESTABLISHED in 2409.09561 eq. (3.9), in r ≠ 0 variables.
- So the λ‑deformation runs along the decoupling edge and leaves the on‑trajectory slice. This is the same phenomenon as in §3.

**Positivity at D = 4.** Positivity of A and M at D = 4 for every λ we tested (levels ≤ 30) reproduces the corresponding statement of 2408.03362.

## 5. Closed strings: no corner (OUR ATTEMPT)

Scripts: `T1_2_closed.py` (n ≤ 100), `T1_2_lambda.py`, `T1_2_closed_intercept.py`.

| amplitude | D_crit | set by (level, spin) | D_crit at n = 100 |
|---|---|---|---|
| Virasoro–Shapiro stripped (2408.03362 eq. (1), (3)) | **23** exactly (ν = 10) | (4, 4) | 43.28 |
| four‑dilaton type II (2102.08951 eq. (3.21)) | **24.14788648641709770…** | (4, 4) | 43.38 |
| closed bosonic tachyons (residue = open²; 2201.11575 eq. (3.22)) | **57** | (5, 8) | 87.57 |

**Footnote‑14 claim of 2102.08951 ("experimentally positive for D ≤ 24.148…").**
- **Confirmed and made exact:** D_crit is the unique real root of **22D³ + 109D² − 10602D − 117329 = 0**.
- It comes from level 4, spin 4. All levels ≤ 100 are positive at D = 24, and the inequality fails at D = 24.148.
- The two lowest per‑level values are 24.1479 (n = 4) and 24.9299 (n = 5).

**VS‑ratio intercept families.**
- Definition: M = Γ(−α_s)Γ(−α_t)Γ(−α_u)/Π Γ(1+α), with S = α_s + α_t + α_u ∈ ℤ.
- Residues: ±[(α_t+1)_{n−1−S}]²/(n!)². S = −1 contains the closed bosonic string at a = 1; S = 0 contains VS at a = 0.
- In both families D_crit(a) is smooth and increasing through the string point: 9, 23, 39.4, 57, 75.3 at a = −1/2, 0, 1/2, 1, 3/2 for S = −1. The S = 0 family gives the same values.
- Positivity never bounds a from above, because the level‑1 residue is a perfect square.
- "No spinning tachyon" forces a ≤ 1, and the maximum is then 57, not 26.

**CHR closed λ‑family M (eq. (7), (13)).**
- D_crit(λ) = 3(λ+2)² − 3 (level 3, spin 2) for small λ. It is 9.76, 15.75, 23, 30.97, 39.4, 57 at λ = 1/16, 1/2, 1, 3/2, 2, 3.
- It is strictly increasing, so **VS is not extremal**.
- These values coincide exactly with the S = 0 family at a = (λ−1)/2. This is the closed analogue of 2409.09561 eq. (3.9).

**Conclusion.** Four‑point positivity of closed strings sees neither D = 10 nor D = 26, and has no intercept corner. The corner phenomenon is an open‑string, on‑trajectory statement.

## 6. Novelty assessment (honest)

Papers checked (abstract pages fetched; full text read where marked): 2201.11575 (full), 2608.25183 (full), 2502.20372 (full), 2102.08951 (full), 2408.03362 (full), 2302.12263 (full), 2409.09561 (full), 2201.12331 (full), 2110.01578, 2406.04410, 2212.00764, 2603.04485, hep‑th/9712207 (search highlights only for the OCQ statement).

**Already known.**
- **Fock‑space version (classic).** The no‑ghost theorem for old covariant quantization: ghost‑free iff D = 26 and α(0) = 1, or D < 26 and α(0) ≤ 1. This is quoted in hep‑th/9712207, citing Brower 1972 and Goddard–Thorn 1972; I read it there only, eq. not checked. So (1, 26) is the corner of the OCQ ghost‑free region. That region is a different, larger object than ours: it allows D ≤ 25 at every a < 1.
- **Four‑point version on the on‑trajectory slice: substantially known numerically.**
  - 2201.12331: q → 1 limit of its Fig. 2/Fig. 3, for m² ∈ [−1, 1/3], giving 26, 10 and ≈6.3. It states that intercepts above 1 are non‑unitary, and it identifies m² = 0 with D = 10.
  - 2302.12263 Fig. 1 (left panel), the r = 0 line.
  - 2409.09561 Fig. 1 and Theorem 1. The theorem gives a ≤ 1 for c = 0 and a ≤ 0 for c = 1 at r = 0.
  - 2302.12263 and 2409.09561 also describe a "cusp" critical point at (r, m², D) = (−1/2, 0, 5).
- **Verdict: the qualitative statement "the critical string is the max‑D corner of the intercept slice" is partially known.** It is implicit in those figures, not stated as a theorem.

**Apparently new (my reading of the papers above).**
1. The exact four‑piece boundary D_max(a), with algebraic transition points, including the narrow level‑4 window.
2. Rigorous exact verification to level 100 (150 at a = 1).
3. Extension below threshold to a > −21/29 using the ρ rule, and the exact value 32/5 at m² = 1/3.
4. The parity lemma.
5. The off‑trajectory failure of the corner (δ\*).
6. The c = 1 boundary formula.
7. The closed‑string numbers 23 and 57.
8. The exact algebraic D_crit = 24.1478864864… for 2102.08951 eq. (3.21), with its level and spin.
9. The closed λ ↔ S = 0 equivalence.

All of items 1–9 are refinements; none changes the conceptual picture.

## 7. Failures, dead ends, caveats
- The seed's left boundary (a ≥ −1/3) was wrong because of the sub‑threshold sign. The seed also did not identify the level‑5 and level‑4 pieces; its value at a = 0.4 was already level 5.
- My first run of `T1_1_highlevel_check.py` had a tuple‑unpacking bug. It was fixed and rerun; only the rerun outputs are kept.
- My `T1_1_transitions.py` first reported an L5 → L2 transition at a = 0.519. That crossing is real but **not on the boundary**, because the level‑4 window intervenes. The output file is annotated.
- There is no all‑level proof for any a ≠ 0, 1. The (a, m²) statements are grid‑based. Results for λ ≥ 8 are not converged at n ≤ 30.
- §3 shows the corner needs a physical input: externals on the trajectory, or no heavier than the level‑0 state. Unitarity of the four‑point amplitude alone, with a free external mass, does not single out (1, 26). This mirrors Front E's "one UV‑shaped assumption".
- For a ≠ 1, the intercept‑shifted Veneziano formula is not the OCQ string amplitude. Four‑point positivity is necessary, not sufficient (cf. the seven‑point fixing of the intercept in 2603.04485).
