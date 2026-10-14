# RT1 — Adversarial review of phase3/targets/T1 (partial-wave positivity vs intercept a and dimension D)

Reviewer: RT1. Under review: `phase3/targets/T1/T1_report.md` (+ scripts, `pwcore.py`).
All my code and outputs are in `phase4/RT1_code/`. I wrote my own code before reading `pwcore.py`.

**How my code differs from theirs:**
- They project monomials onto Gegenbauers with python‑flint `fmpq_poly`, then isolate roots with arb.
- I multiply the residue's linear factors (b + κx) one at a time into the Gegenbauer series, using the three‑term recurrence x·C_j = [(j+1)C_{j+1} + (j+2ν−1)C_{j−1}]/(2(j+ν)). The arithmetic is exact rational (flint `fmpq` is used only as a fast number type). D_crit per level comes from an exact bisection, which relies on monotonicity in D.
- The low levels are cross‑checked symbolically in sympy with a third method: the monomial formula, with D symbolic.

Conventions I derived independently:
- t = κ(x−1), with κ = (s_n − 4m²)/2 and x = 1 + 2t/(s_n − 4m²).
- Residue convention: A ≈ R/(M² − s).

## Verdict table

| # | Claim | Verdict | One‑line reason | My script(s) |
|---|---|---|---|---|
| 1 | Open c=0 on‑trajectory region: −21/29 < a ≤ 1, D_max(a) built from 4 pieces (L3, L5, L4, L2 spin 0), transitions 0.3883687 / 0.5157634 / 0.5221548, unique max (1,26), D_max(0)=10, D_max(−1/3)=32/5 | **VERIFIED** (numerical, not a proof) | Reproduced exactly: every formula, all 3 transition polynomials, monotonicity, and no higher level cuts (917 a‑values to level 60; 15 a‑values to level 100). The lower endpoint holds down to a = −20. *Their own scans never went below a = −0.75*, although the level‑3 bound changes sign for a < −1; I closed that gap. The r=0 slice is partly known (see Novelty). | `rt1_lowlevels_symbolic.py`, `rt1_transitions_exact.py`, `rt1_region.py`, `rt1_margin.py`, `rt1_region_negative_a.py`, `rt1_empty_check.py`, `rt1_empty_cover.py`, `rt1_far_negative.py` |
| 2 | Below threshold the correct condition is ρ_j = c_j/κ^j ≥ 0; the m² ≤ 1/3 (a ≥ −1/3) cut in 2201.12331 / 2302.12263 is a convention artifact | **VERIFIED** | Explicit Feynman‑rule check: a Hermitian, real‑coupled Proca field below threshold gives c_1 ∝ e²(M²−4m²) < 0. The Veneziano level‑1 vector residue *is* that Proca residue with e² = 1/2. 2201.12331 eqs. (13)–(14) and 2302.12263 do impose a_{n,ℓ} ≥ 0 below threshold. | `rt1_subthreshold.py` |
| 3 | Parity lemma: for R of definite parity, (1+x)R is Gegenbauer‑positive iff R is; hence the a=0 c=0 amplitude and type‑I share D_crit(n) | **VERIFIED** | Short proof below (it is elementary). Random test agrees 11840/11840. D_crit(n) is identical for n = 2..60. | `rt1_parity_and_closed_explain.py` |
| 4 | a=1, m²=−1+δ: D_max first exceeds 26 at δ* = 0.4022459 (level 4, spin 0); edge D = (2a+3)²+1; CONJECTURE: (1,26) is the unique max for m² ≤ −a | **VERIFIED** (numbers); novelty overstated | δ* = 0.40224589837224 is a root of the exact level‑4 quartic. The line values and the edge formula are reproduced, and a 924‑point grid supports the conjecture. **Missed dictionary:** the free‑mass (a, m²) plane is *exactly* the hypergeometric (r, m²) plane of 2302.12263/2409.09561. The edge is 2409.09561's Theorem‑1 line 2r+m²+1 = 0, and "m² ≤ −a" means r ≤ 0. | `rt1_massscan.py`, `rt1_hypergeo_map.py` |
| 5 | c=1 family: a ≤ 0; D_max = (8a+10)/(1−a) from level 3 spin 0; corner (0,10) | **VERIFIED** | Reproduced on 128 a‑values to level 40 and 7 a‑values to level 100. Lower endpoint −7/11. The region is empty for a ∈ [−20, −7/11). The a ≤ 0 bound is the r=0 case of 2409.09561 Thm 1, as the report says. | `rt1_c1.py`, `rt1_far_negative.py` |
| 6 | Closed: VS stripped D_crit = 23 (4,4); four‑dilaton (2102.08951) D_crit = root of 22D³+109D²−10602D−117329 = 24.14788648641709770 (4,4), confirming fn. 14; closed bosonic tachyon 57 (5,8); no corner for general intercept | **VERIFIED** (numbers); "surprising" and "new" are overclaimed | All numbers reproduced exactly. Footnote 14 is confirmed as quoted. 2201.11575 already states that closed‑string 4‑pt "critical dimensions are generically higher than their open string counterparts". The S=0 and S=−1 families are the same residues shifted by one level, so they do not check each other. | `rt1_closed.py`, `rt1_parity_and_closed_explain.py` |
| 7 | CHR λ‑families: planar D_crit = (λ+2)²+1; closed 3(λ+2)²−3 at small λ; no corner at λ=1; planar is in 2409.09561 eq. 3.9; closed claimed new | **VERIFIED‑KNOWN** | Values reproduced, and the eq. 3.9 mapping is consistent. But 2408.03362 **Fig. 1** already maps the (λ, D) positivity region of *both* M (closed) and A (planar), through level 30. The closed closed‑form is one coefficient (level 3, spin 2). The planar result is literally the same statement as claim 4's edge formula; at λ=3 the planar family *is* the bosonic Veneziano amplitude, shifted one level. | `rt1_closed.py` |

**Label check:**
- The ESTABLISHED / OUR ATTEMPT / CONJECTURE labels are honest. Nothing labelled CONJECTURE is presented as established, and "not a proof for all levels" is stated.
- Two calibration problems:
  - (i) §6 lists items 5, 7 and 9 as "apparently new" while the same content sits in published figures (see Novelty).
  - (ii) The region claim "−21/29 < a ≤ 1" was stated for all a, but the authors' scans covered only a ≥ −3/4.

---

## Claim 1 — details

**Exact low levels.** Rederived with the monomial→Gegenbauer formula and D symbolic (`rt1_lowlevels_symbolic_output.txt`).
- Level 2: ρ_0 = [D(a²−2a)+8a²+14a+4]/(8(D−1)) and ρ_1 = (1−a)/(2(D−3)).
- Level 3: ρ_0 = (1−a)(1+a)[(26a+30) − D(3−a)]/(48(D−1)) and ρ_1 = [D(3a²−6a−1)+30a²+48a+26]/(24(D−3)(D+1)).
- Levels 4 and 5 give spin‑0 numerators identical to the report's (up to an overall sign that does not change the roots).
- At a = 1, my level‑2 ρ_0 = (26−D)/(8(D−1)). This is 2608.25183 eq. (5) exactly, which also checks normalisation.

**Transitions** (`rt1_transitions_exact_output.txt`):
- L3/L5: −8(a−3)(45a⁵+78a⁴+20a³−64a²−61a+30), root a₁ = 0.388368678214867.
- L2/L4: −4(a−2)(18a⁶+51a⁵−25a⁴−104a³+32a²+72a−32), root a₃ = 0.522154773496928.
- L4/L5: the resultant has a degree‑15 factor with root a₂ = 0.515763405625103. The same factor's other root, 0.39183, is the off‑boundary crossing.
- D3(a) = 3 exactly at a = −21/29. D_max(−1/3) = 32/5, D_max(0) = 10, D_max(−3/5) = 4.

**Region scan** (`rt1_region_grid_N60_output.txt`):
- Grid: 917 a‑values (step 1/500 on (−0.724, 1], plus step 1/20000 near a₁, a₂, a₃).
- Levels ≤ 8 by exact bisection (tolerance 1e‑10), then levels 9..60 checked exactly at D_max: all positive.
- max |D_max − piecewise formula| = 1.0e‑10, which is the bisection tolerance.
- Binding pieces on the grid:
  - L3 on [−0.724, 0.3882];
  - L5 on [0.3884, 0.5156];
  - L4 on [0.5158, 0.5220];
  - L2 on [0.5222, 1].
- D_max is strictly increasing on the grid.

**Level 100.** At 15 values a ∈ {−0.72, −0.7, −0.6, −1/2, −1/3, −1/4, 0, 1/4, 2/5, 1/2, 0.52, 3/4, 9/10, 0.99, 1}, all levels 13..100 are exactly positive at D_max (`rt1_region_selected_output.txt`).

**Higher‑level margin** (`rt1_margin_N40_output.txt`): min_{6≤n≤40} D_crit(n) − D_max is
- 0.839 at a = −0.724,
- 0.489 at a = 0,
- 0.0209 at a = 0.5158 (level 6 spin 0),

which reproduces their 0.0209. It is narrow but positive.

**Lower endpoint: the gap in their verification.**
- Their region runs started at a = −18/25 and −3/4 (headers of `T1_1_region_c0_N30/N40_output.txt`).
- For a < −1 the level‑3 spin‑0 factor (1+a)[…] changes sign and becomes **positive**, so emptiness below −1 is not automatic. Explicit example: at a = −2 levels 1–4 are all positive.
- What I checked:
  - every a ∈ [−6, −21/29) on a 1/1000 grid with levels ≤ 30, at D = 3+10⁻⁶: no survivors;
  - a 1/2000 grid with levels ≤ 14, where adjacent grid points always share a failing level, so no window can hide between them;
  - spot checks at a = −7 … −61/3 with levels ≤ 90. The first failing level is about 2|a|+1.
- Conclusion: no allowed window exists. The claim stands, but it was not verified by the authors.

**Caveat (shared with the report).** Everything here is levels ≤ 60–100, not all levels, except at a = 1, where 2608.25183 proves all levels.

## Claim 2 — the sub‑threshold sign (most consequential claim)

**Derivation.**
- Factorisation at s = M² gives Res = Σ_j g_j² N_j (s−4m²)^j C_j^ν(x), with N_j > 0 and real cubic couplings.
- The projector contraction depends only on invariants. Since q·k = 0, −q·q′ = t−u = (s−4m²)x and −q² = s−4m².
- So for spin 1, −q·q′ = (s−4m²)x. For spin 2, (q·P̃·q′)² − (q·P̃·q)(q′·P̃·q′)/(D−1) = (s−4m²)²[x² − 1/(D−1)] ∝ (s−4m²)² C_2.
- Therefore g_j² ∝ c_j/κ^j, and the necessary condition is ρ_j ≥ 0.

**Explicit Lagrangian check** (`rt1_subthreshold.py`, D = 5, real momenta continued to imaginary |p|).
- Setup: scalar QED with a Proca field of mass M, A ≈ Res/(M²−s). The scalar exchange residue is +g₀².
- Result: Res_V = e²(t−u) = e²(M²−4m²)x. The polarisation‑sum form and the formula agree numerically above and below threshold.
- Consequences:
  - a perfectly healthy vector below threshold has c_1 < 0;
  - the rule "c_j ≥ 0 for all poles" would declare it a ghost.
- The Veneziano level‑1 residue is (1−a)/2 + κx, and its vector part κx = ½(M²−4m²)x is exactly this Proca residue with e² = 1/2.

**What the papers actually do** (fetched and read):
- **2201.12331** (Figueroa–Tourkine), eq. (13): c_{1,1} = (1−3m²)q/(2(d−3)). They write: "for m² > 1/3, the Coon amplitude has ghosts in all dimensions … this stems from the relation between cos θ and t … which becomes singular at the first resonance". They then restrict to −1 ≤ m² ≤ 1/3 (eq. 14). **This is the artifact.**
- **2302.12263** (Cheung–Remmen), Unitarity Bounds section: "we analytically find m₀² ≤ 1/3 from a_{n,n} ≥ 0". Same artifact.
- **2403.00906** (hypergeometric Coon) repeats it at its eqs. (26)–(27): "−1 ≤ m² ≤ 1/3 … introduced by c_{1,1}".
- **2409.09561**, footnote 2: imposes B_{n,j} ≥ 0 only for n ≥ 3m², i.e. it *drops* sub‑threshold levels. This is weaker, not wrong. On the c=0 slice it gives the same region, because sub‑threshold levels never bind under the ρ rule: level 1 has ρ = ((1−a)/2, 1/(D−3)) and level‑2 ρ_0 > 0 for a ∈ (−0.724, −2/3).
- I found no paper that corrects the m² ≤ 1/3 bound. Two Exa searches turned up only papers that repeat the bound: 2201.12331, 2302.12263, 2403.00906, 2401.13031 and 2303.02149. My confidence that none exists is moderate.

**Verdict: VERIFIED**, including the literature attribution.

**Note.** The report's exact D_max(−1/3) = 32/5 = 6.4 agrees with 2201.12331's "d ≃ 6.3" only within that paper's stated 4% numerical accuracy. The report says the same; there is no discrepancy.

## Claim 3 — parity lemma

**Proof.**
- For ν > 0, x·C_j has nonnegative coefficients: (j+1)/(2(j+ν)) and (j+2ν−1)/(2(j+ν)). So R ≥ 0 implies xR ≥ 0, and hence (1+x)R ≥ 0.
- Conversely, xR lives in the opposite parity sector. The coefficients of (1+x)R in R's sector are therefore exactly R's.

**Application.**
- At s = n with m = 0: (t+1)_n = (t+1)_{n−1}·(n(1+x)/2).
- (t+1)_{n−1} has parity (−1)^{n−1} under t ↔ u.
- So the lemma applies, and D_crit(n) is equal for all n. I verified this for n = 2..60.

**Caveats.**
- The lemma holds for the ρ‑rule only when κ > 0. Below threshold, (1+x) effectively becomes (1−x) and the "if" direction fails. This is irrelevant here (m = 0).
- It is elementary. The m²=0 → d=10 fact itself is already in 2201.12331 Fig. 3; the lemma only explains it.

## Claim 4 — off‑trajectory mass

**Exact checks** (`rt1_massscan_exact_output.txt`, at a = 1, D = 26):
- level 2 spin 0 vanishes at δ = 5/52;
- level 3 at δ = 1/14 and 1/4;
- level 4 at δ = 1/4 and **0.40224589837224051**;
- level 5 at δ = 1/8 and 1/2 (double roots).

**Line scan** (levels ≤ 40):

| δ | D_max | binding level |
|---|---|---|
| 1/20 | 25.000 | 2 |
| 1/4 | 26 | 3 |
| 0.4022 | 25.9994 | 4 |
| 0.4023 | 26.0007 | 4 |
| 1/2 | 26 | 5 |
| m² = 0 | 26.849 | 8 |
| m² = +1/2 | 27.180 | 11 |

All match the report.

**Edge m² = −(1+3a)/4.**
- Level‑2 spin 0 gives D ≤ 4a²+12a+10 = (2a+3)²+1.
- D_max equals this for a ∈ {−1/2, …, 1} with levels ≤ 30.
- For a > 1 higher levels bind (44.857 at a = 2), as reported.

**Conjecture grid.**
- Grid: a step 1/16, m² ∈ [−(1+3a)/4, −a] in 33 steps, levels ≤ 20; 924 points.
- Max = 26, only at (1, −1). The conjecture is supported, still grid‑only.

**Missed dictionary (novelty problem; the numbers are correct).**
- Level by level, T1's c=0 residue at (a, m_ext²) equals the bosonic‑type hypergeometric residue of 2409.09561 eq. (2.4) up to the positive factor 1/(r+1)_n.
- The map is r = (4/3)(m_ext²+a), m² = (4/3)m_ext² + a/3. I verified identical D_crit(n) at 5 points (`rt1_hypergeo_map_output.txt`, max difference 0).
- Consequences:
  - on‑trajectory ↔ r = 0;
  - the δ‑line ↔ r = 4δ/3, m² = −1 + 4δ/3;
  - the level‑1 edge ↔ 2r + m² + 1 = 0, which is exactly the boundary in 2409.09561 Theorem 1, from B_{n,n−1} ∝ 2r+m²+a;
  - "m² ≤ −a" ↔ r ≤ 0.
- So "the corner is not robust off‑trajectory" is a statement about the already‑mapped (r, m²) plane: 2302.12263 Fig. 1 (left) and 2409.09561 Fig. 1. The report's dictionary in §0 mentions only the r = 0 slice.
- The specific δ* and the edge formula are not stated in those papers (as far as I could read; I could not view the figure images).

## Claim 5 — c=1 family

- Level 2: ρ_0 = −a/4. Level 3: ρ_0 = (a+1)(D(a−1)+8a+10)/(24(D−1)), so D_max = (8a+10)/(1−a).
- Binding (3,0) at all 128 grid points in (−7/11, 0]. Higher levels are positive to level 40, and to level 100 at 7 points.
- The margin to the next level is ≥ 0.24 (level 5).
- Lower end: −7/11. The region is empty on [−6, −7/11) (step 1/1000, levels ≤ 14, overlap check) and at a = −7 … −61/3.
- The a ≤ 0 bound is 2409.09561 Theorem 1 at r = 0 (read).
- By 2409.09561 eq. (3.4) this slice is the bosonic hypergeometric line r = −1/3, so its D‑contours are in their Fig. 1. The closed form is new but is a one‑coefficient computation.

## Claim 6 — closed strings

**Reproduced** (`rt1_closed_N40_output.txt`, levels ≤ 40):

| Amplitude | D_crit | set by (level, spin) |
|---|---|---|
| VS stripped | 23.0000000000 | (4,4) |
| four‑dilaton, P = (st)²+(tu)²+(su)² × VS | 24.147886486327 | (4,4); bisection lower bound, matches the cubic's real root 24.14788648641709770 to 1e‑10 |
| closed bosonic tachyon | 57 | (5,8) |

- Next levels: dilaton (5, 24.9299), closed bosonic (6, 57.3714).
- S = −1 family: 9, 23, 39.4, 57, 75.2857 at a = −1/2 … 3/2. It is increasing, so there is no corner.

**Literature checks.**
- 2102.08951 (Caron‑Huot, Mazáč, Rastelli, Simmons‑Duffin, "Sharp boundaries for the swampland"): the four‑dilaton amplitude is in §3.6. In the ar5iv rendering it is eq. (58). The report's "(3.21)" could not be confirmed. Footnote 14 reads "…experimentally admit a positive expansion … 2 ≤ D ≤ 24.148…", so it is confirmed.

**Why the 4‑tachyon amplitude need not see D = 26.**
- The closed residue is the square of the open residue (2201.11575 §1, KLT). So its spin‑0 coefficient is ⟨R²⟩ > 0 *identically*.
- The open‑string D>26 obstruction is a spin‑0 coefficient (level 2), so it can never appear in the closed scalar channel.
- Explicitly, at D = 27:
  - the open level 2 has ρ = (−0.0096, 0, 0.0032);
  - its square has ρ = (0.1033, 0, 0.0024, 0, 4.6e‑5), all positive.
- State‑space reading:
  - ghost⊗ghost states have positive norm;
  - ghost⊗good states sit in spin‑j channels, where good⊗good states with the same quantum numbers outweigh them;
  - the 4‑point test only sees sums over states with a given (mass, SO(D−1) irrep).
- This is not surprising. 2201.11575 states explicitly that closed‑string 4‑point "critical dimensions" are generically higher than the open ones.

**Families.** The S = 0 family at level n+1 has *exactly* the residue of the S = −1 family at level n, at the same a. So "the S=0 family gives the same values" is a tautology, not a cross‑check.

## Claim 7 — Cheung–Hillman–Remmen λ‑families (2408.03362, read)

**Values reproduced.**
- Planar: 5.2539, 7.25, 10, 13.25, 17, 26 at λ = 1/16, 1/2, 1, 3/2, 2, 3, all from (3,0).
- Closed: 9.7617 and 15.75 at λ = 1/16 and 1/2, from (3,2).
- Closed at larger λ: 23 at λ=1 from (4,4); 30.969 at 3/2 from (4,4); 39.4 at 2 from (5,6); 57 at 3 from (6,8).

**Derivations.**
- Closed level 3 is (κ²x²−1/4)², with κ = (λ+2)/2. Its spin‑2 coefficient gives D ≤ 12κ² − 3 = 3(λ+2)² − 3.
- As polynomials in x, CHR level n equals the S=0 (closed) or c=1 (planar) family at a = (λ−1)/2 with m² = −3a/4, as the report says.
- Further, c=1 at (a, −3a/4), level n+1, equals c=0 on the level‑1 edge (a, −(1+3a)/4), level n.
  - So **"planar D_crit = (λ+2)²+1" and claim 4's "edge D = (2a+3)²+1" are the same statement.**
  - At λ = 3 the planar CHR amplitude has exactly the bosonic Veneziano partial waves, shifted one level: next levels 26.0, 26.4, 26.589, the bosonic sequence.
  - This is consistent with 2409.09561 eqs. (3.4)+(3.9): B^planar(λ) ∝ B¹((λ−3)/6, −λ/3), and λ = 3 gives (0, −1).

**Novelty.**
- 2408.03362 Fig. 1 ("Values of λ and D consistent with positivity … through level n = 30. Orange: M and A both positive. Red: additional region where M alone is positive") already maps both D_crit(λ) curves numerically. The text also states "When D ≥ 9 positivity bounds λ from below" (closed, consistent with 3·4−3 = 9) and "When D ≥ 5 …" (planar, consistent with 5).
- The "closed is new" claim therefore reduces to writing a closed form for the small‑λ branch of a curve CHR already plotted.

## Novelty statements in §6, checked

**Correct as written:**
- 2201.12331: q→1 gives 26, 10 and ≈6.3; "intercepts greater than 1" are non‑unitary.
- 2302.12263 Fig. 1: q=1 panel, with the "critical case" r = −1/2, m² = 0, D = 5.
- 2409.09561: Thm 1, eq. 3.9, footnote 2, Fig. 1 with the D = 5 critical point.
- 2608.25183 (Chen–Yin, Aug 2026) exists. It proves bosonic positivity for 3 < D ≤ 26 at all levels; its eq. (5) and eq. (23) (connection formula) are as cited.
- hep‑th/9712207 exists (Daszkiewicz–Hasiewicz–Jaskólski, "Free strings in non‑critical dimensions"); the OCQ statement is the standard Brower / Goddard–Thorn result. The report honestly says "eq. not checked".
- 2603.04485 exists (Basile–Remmen–Staudt; seven‑point fixing of the intercept).
- 2201.11575 and 2502.20372 match the report's descriptions.

**Overclaimed:**
- Item 5 (off‑trajectory) is the hypergeometric (r, m²) plane.
- Item 7 (closed 23, 57): qualitatively in 2201.11575; VS = λ=1 of CHR Fig. 1.
- Item 9 (closed λ ↔ S=0): the numerical region is CHR Fig. 1; the equivalence itself is a correct and simple observation.
- Items 1–3 and 8 are fair as "exact refinements".
- Item 3 (below‑threshold extension) is genuinely useful: the published m² ≤ 1/3 cut is an artifact. However, the region with D ≥ 4 is probably already visible in 2409.09561 Fig. 1 under its drop‑sub‑threshold rule. I could not view the figure (arXiv/OSTI/Springer PDFs are blocked here), so I mark that as unverified.

## Bottom line
- No numerical claim of T1 failed. Every number I re‑derived independently agrees to the bisection tolerance (≤ 1e‑10).
- The main issues:
  - (a) an unverified gap for a < −3/4 in the region claim, which I closed: the region is empty down to a ≈ −20;
  - (b) novelty overclaims, because the off‑trajectory plane is the known hypergeometric plane, and the closed‑string and λ results are already in 2201.11575 and 2408.03362 Fig. 1;
  - (c) two results are presented separately but are the same statement: the claim‑4 edge and the claim‑7 planar formula.
- The physics claim most worth keeping is claim 2: the m² ≤ 1/3 bound used in 2201.12331, 2302.12263 and 2403.00906 comes from applying c_j ≥ 0 to sub‑threshold poles. A healthy real‑coupled vector violates that rule.
