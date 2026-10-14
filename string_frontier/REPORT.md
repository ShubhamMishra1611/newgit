# String Frontier: what one day of agentic research could and could not do

**Mission.**
1. A non-perturbative definition of string/M-theory valid outside AdS.
2. A principle that selects our vacuum: the Standard Model with Λ ≈ 10⁻¹²² in Planck units.

**Process.**
- Phase 1: five parallel frontier-mapping agents.
- Phase 2: target selection.
- Phase 3: sanity reproductions, then three parallel compute agents.
- Phase 4: seven adversarial reviewers, each a fresh agent that did not produce the work it reviews, plus one second referee on the most consequential claim.

**Labels used throughout.**
- **ESTABLISHED:** a primary source, cited by arXiv ID and equation number where read.
- **CONJECTURE.**
- **OUR ATTEMPT:** a computation run here.

**Review verdicts.** VERIFIED, VERIFIED-KNOWN, PLAUSIBLE, FAILED, or CITATION PROBLEM, as defined in `phase4/REVIEW_PROTOCOL.md`.

---

## 0. Bottom line (read this first)

**Neither mission problem was solved.**
- No non-perturbative definition outside AdS was produced.
- No vacuum-selection principle was found.
- No step taken here brings either within reach in a way we can defend.

The honest yield is smaller and real. These items survived independent adversarial recomputation:

1. **A genuine error in a published celestial-holography result.**
   - The pole-sum formula for the celestial four-gluon open-superstring amplitude in Stieberger–Taylor, arXiv:1806.05688, eqs. (3.30)/(3.31), is wrong: it drops the u-channel poles.
   - Their small-angle claim (corrections O(r^{β−1})) is false for Re β < 1, which includes the physical line β ∈ iℝ. The same claim is repeated in arXiv:2405.01643 eq. (3.7).
   - A corrected closed form matches direct integration to ≤10⁻¹⁷.
   - Three separate agents converged on this: the Front-B mapper found it, reviewer R1 diagnosed and corrected it, and a second referee RB4 re-derived everything with different methods. No published correction was found. (§5.5, `phase4/R1_verdicts.md`, `phase4/RB4_verdicts.md`)

2. **A sign-convention error in the dual-amplitude bootstrap literature.**
   - Below the two-particle threshold, Gegenbauer positivity must be imposed on ρ_j = c_j/κ^j, not on c_j.
   - A healthy real-coupled massive vector below threshold has c₁ < 0.
   - Three published papers impose c_j ≥ 0 (2201.12331 eqs. 13–14, 2302.12263, 2403.00906). Their cut at m² ≤ 1/3 is therefore an artifact, and the unitarity region extends further.
   - Reviewer RT1 VERIFIED this and found no published correction. (§5.2)

3. **Exact unitarity boundaries for Veneziano-type amplitudes.**
   - Open strings with intercept a and external states on the leading trajectory: the allowed region is −21/29 < a ≤ 1, 3 < D ≤ D_max(a), with an exact four-piece algebraic boundary. Its unique maximal-D point is the critical bosonic string (a, D) = (1, 26).
   - The superstring-type family has the analogous corner at (0, 10).
   - These are numerical results to level 100, not a proof for all levels. The qualitative corner was already visible in published figures. The exact boundary, the sub-threshold extension and the parity lemma are new as far as we found.
   - Closed strings and freely varying external masses show **no** such corner. So four-point unitarity alone does **not** select the critical string. (§5.2)

4. **An exact form of a published numerical observation.**
   - Footnote 14 of 2102.08951 says the type-II four-dilaton residues are "experimentally" Gegenbauer-positive for D ≤ 24.148…
   - The exact value is the real root of 22D³ + 109D² − 10602D − 117329 = 0, i.e. D = 24.14788648641709770…, set by level 4, spin 4. VERIFIED by RT1. (§5.2)

5. **Explicit control failures of the 2024 KKLT de Sitter candidates (2406.13751), from the authors' public data.**
   - In the warped-throat EFT, the conifold modulus of all five representative vacua has no minimum once the anti-D3 is added. The BDGL bound g_sM² ≳ 46 is missed by factors 2.75–3.9.
   - This survives a full convention audit, but it is conditional on an off-shell warp-factor assumption that is itself debated.
   - The uplift strongly backreacts on the dilaton and complex-structure sector.
   - Other T2 claims were downgraded to known (the α′² KPV loss of control at g_sM ≈ 1, the conifold mass ≈ 0.5 m_wKK) or FAILED (a "yield/efficiency" comparison). (§5.3)

6. **Λ selection (T3).** [filled in §5.4 after review]

**Reproductions** (sanity checks and Phase 1 checks), all re-verified by reviewers:
- D = 26 four independent ways, including BRST: Q² ∝ Σ_m[(c−26)(m³−m)/12 − 2m(A+1)] c₋ₘcₘ.
- Veneziano poles, duality, Regge behaviour, and the D ≤ 26 partial-wave bound with c₀ = −(D−26)/(8(D−1)) at level 2.
- One-loop σ-model β functions, matched to Polyakov's O(N) result and exact backgrounds.
- KKLT/KPV/BDGL numbers.
- Bootstrap LP bounds.
- α_min of the 10D graviton bootstrap.
- D0-brane phase shifts.
- ABJM Fermi gas.
- The dS one-loop logarithms.

**Failures reported honestly.** Several of our own claims FAILED review (§6):
- Front E's "SUSY excludes the λ-family": a missing parity assumption, and not new.
- KKLT dS-window edges: grid values, not exact.
- A sign mislabel in the D0 v¹⁰ coefficient.
- T2's yield comparison: ill-posed.
- "4.1σ" DESI tension: PLAUSIBLE only; quote ≳4σ, dataset-dependent.

---

## 1. Scorecard

| # | Claim (short) | Label | Review verdict | Where |
|---|---|---|---|---|
| S1 | D=26 from Virasoro: Fock c=D, c_gh=−26, a_gh=−1, Kac level-2 null state (h=−1, γ=3/2, c=26), BRST Q²=0 ⇔ c=26, A=−1, light-cone ζ(−1) | ESTABLISHED, reproduced | VERIFIED-KNOWN | `phase3/sanity/s1_*` |
| S2 | Veneziano: residues, pole-sum = amplitude, t-pole from s-tower tail, Regge, c₀(level 2) = −(D−26)/(8(D−1)) ⇒ D ≤ 26 | ESTABLISHED, reproduced | VERIFIED-KNOWN (2201.11575 eq. 1.3; 2608.25183) | `phase3/sanity/s2_*` |
| S3 | One-loop β functions: RNC −1/3, lattice tadpole −1/(4π), O(N) Polyakov, cigar, S³ WZW, EOM ↔ β | ESTABLISHED, reproduced | VERIFIED-KNOWN | `phase3/sanity/s3_*` |
| B4 | Stieberger–Taylor 1806.05688 (3.30)/(3.31) wrong (missing u-channel poles); corrected formula; small-angle claim false for Re β<1; repeated in 2405.01643 (3.7) | OUR ATTEMPT | **VERIFIED** ×2 (R1, RB4), no prior correction found | `phase1/B_checks`, `phase4/R1_*`, `phase4/RB4_*` |
| T1.2 | Sub-threshold positivity must use ρ_j = c_j/κ^j; published m² ≤ 1/3 cut is an artifact | OUR ATTEMPT | **VERIFIED** | `phase3/targets/T1`, `phase4/RT1_*` |
| T1.1 | c=0 on-trajectory region −21/29 < a ≤ 1, four-piece D_max(a), unique corner (1,26) | OUR ATTEMPT (numerical to level 100) | **VERIFIED** (not a proof; qualitative corner partially known) | same |
| T1.3 | Parity lemma: definite-parity R ⇒ (1+x)R positive iff R is | OUR ATTEMPT | **VERIFIED** | same |
| T1.4 | Corner not robust off-trajectory: δ* = 0.4022459; edge D = (2a+3)²+1 | OUR ATTEMPT | VERIFIED (numbers); novelty overstated (= 2302.12263/2409.09561 (r, m²) plane) | same |
| T1.5 | c=1 family: D_max = (8a+10)/(1−a), corner (0,10) | OUR ATTEMPT | **VERIFIED** | same |
| T1.6 | Closed strings: VS-stripped D_crit = 23 (4,4); 4-dilaton D_crit = root of 22D³+109D²−10602D−117329 = 24.1478864864…; closed bosonic tachyons 57; no corner | OUR ATTEMPT | VERIFIED (numbers); "surprising" overclaimed | same |
| T1.7 | CHR λ-families: planar (λ+2)²+1, closed 3(λ+2)²−3 | OUR ATTEMPT | VERIFIED-KNOWN (2408.03362 Fig. 1) | same |
| T2.3 | No conifold minimum with anti-D3 in warped EFT for all five 2406.13751 vacua (×2.75–3.88) | OUR ATTEMPT | **VERIFIED (conditional** on the off-shell warp-factor model) | `phase3/targets/T2`, `phase4/RT2_*` |
| T2.1 | z_cf stays on its F-flat relation; moves only through g_s | OUR ATTEMPT | VERIFIED | same |
| T2.2 | Uplift backreaction R_cs = 9.4…1.3; V_cs/V_D̄3 = 1.52…0.14 | OUR ATTEMPT | VERIFIED; mechanism already in 2512.17095 §7 | same |
| T2.4 | Conifold mass ≈ 0.5 m_wKK, inside the Kähler spectrum | OUR ATTEMPT | VERIFIED-KNOWN (2406.13751 fn. 44, BKS 1902.07724) | same |
| T2.5 | Lüst–Randall extrapolation S/S0 ≈ 0.8 | OUR ATTEMPT (labelled crude) | PLAUSIBLE | same |
| T2.6 | α′²-KPV: control lost at g_sM ≈ 1 | OUR ATTEMPT | VERIFIED-KNOWN; robust number is the tension correction 5.924/(g_sM)² = 5.4–10.7 | same |
| T2.7 | dS window vs flux steps; "74–114 of 396 expected vs 5 found" | OUR ATTEMPT | **FAILED** (ill-posed comparison) | same |
| T3.* | Λ selection: see §5.4 | OUR ATTEMPT | see §5.4 | `phase3/targets/T3`, `phase4/RT3_*` |
| E1 | Planar λ-family excluded by SUSY nonlinear relations except at λ=1 | OUR ATTEMPT | **FAILED** (overclaim: needs the 6-scalar parity assumption of 2601.11705v2; corollary of 2610.02124) | `phase1/E_checks`, `phase4/R2_*` |
| E2–E5 | LP g̃₃ ≥ −10.612487; α_min^IIB = 0.138949, α_min^IIA = 0.140260; E_{3/2}(ρ) = 8.892745; VS expansion; λ-family (g̃₃, g̃₄) | ESTABLISHED reproduced / OUR ATTEMPT | VERIFIED | same |
| A1 | D0 one-loop series 15/16, 0, 9009/4096, ∓182325/16384 | reproduction | **FAILED (minor)**: v¹⁰ sign mislabelled (Euclidean vs Minkowski) | `phase1/A_code`, `phase4/R1_*` |
| A2 | Three-loop sugra target (16875/1024)N₁³N₂v⁸/(R⁷M²⁷r²¹) | OUR ATTEMPT | VERIFIED-KNOWN (next term of BBPT eq. 12); proposed SU(2) test ill-posed | same |
| A3–A5 | BMN fuzzy-sphere spectrum; Firat vertex 1.62187e8 → 6561/4096; ABJM Z(1)=1/4, Z(2)=1/(16π) | reproduction | VERIFIED (A4: Firat's printed eq. 5 cannot reproduce his eq. 8; our branch is correct) | same |
| A6 | IKKT status citations | — | VERIFIED; one label should be CONJECTURE | same |
| B1–B3 | dS log coefficients −571/45, −31/45; w₁₊∞ Jacobi; DSSYK spectrum + Verlinde eq. (3) typo | reproduction | VERIFIED | `phase1/B_checks`, `phase4/R1_*` |
| C1 | BP charge q ≈ 10^−1.43…10^−1.73 at J=100 | reproduction | VERIFIED (BP's 10^−1.6 exact in BP's own normalization) | `phase1/C_checks`, `phase4/R2_*` |
| C2 | DESI λ=0.698: TCC √(2/3) at 0.7σ; √2 at 4.1σ | OUR ATTEMPT | PLAUSIBLE ("≳4σ, dataset-dependent") | same |
| C3 | Dark-dimension l = 0.085–8.5 μm, M̂ = 5e8–2.4e9 GeV | OUR ATTEMPT | VERIFIED (gap to paper = M_P convention × (1/l vs 2π/l)) | same |
| D1 | KKLT σ_cr = 113.589, V_AdS = −1.997e−15; dS windows | reproduction | σ, V VERIFIED; windows **FAILED (minor)**: exact [2.9732, 7.8612]e−9 (D/σ³), [2.5937, 9.7195]e−11 (D/σ²) | `phase1/D_dS_checks`, `phase4/R2_*` |
| D2 | KPV 0.0805; BDGL 6.80, 46.25 | reproduction | VERIFIED | same |
| D3 | The five representative 2406.13751 vacua have g_sM² = 16.82, 14.61, 14.33, 16.16, 11.94 | data | VERIFIED (other 25 vacua not checked) | same |
| E6 | 2025–26 citations | — | 11/12 VERIFIED; **CITATION PROBLEM**: 2601.11705 cited in superseded v1 form | `phase4/R2_*` |

---

## 2. Frontier map (Phase 1)

FRONTIER_MAP_PLACEHOLDER

---

## 3. Targets chosen, and why (Phase 2)

Full justification is in `phase2/TARGETS.md`. A target qualified only if it bore directly on the mission, had a published number to reproduce first, reduced to an explicit computation, and would be informative whatever the outcome.

- **T1: Is the critical string the extremal corner of tree-level unitarity?** (mission 1: consistency-defined flat-space S-matrix)
  - Front E found that every uniqueness theorem relies on one UV-shaped input.
  - A seed computation (ours) suggested that four-point partial-wave positivity alone puts (a, D) = (1, 26) at a corner. That is a checkable "consistency selects the string" statement.
- **T2: Stress-test the only explicit dS candidates (2406.13751).** (mission 2: existence of Λ > 0 vacua)
  - Their parameters are public, and the disputes (conifold, KPV, flux discreteness) are numerical.
- **T3: Λ selection head-to-head.** (mission 2)
  - Anthropic selection on a discretuum hinges on a flat prior near Λ = 0.
  - Swampland quintessence hinges on the achievable exponential slope versus DESI.
  - Each premise is computable.
- **Not chosen:**
  - three-loop D0 scattering (too large);
  - small-N IKKT (weak leverage);
  - O(λ) DSSYK vs dS₃ (unfixed conventions risk a fake result);
  - LP boundary engineering (overlaps T1).
  - The celestial discrepancy went to Phase 4 as a claim, where it became our strongest result.

---

## 4. Sanity reproductions (Phase 3, first step)

All code and outputs are in `phase3/sanity/`. R1 and R2 independently re-derived each of these: VERIFIED-KNOWN.

**S1, D = 26 from the Virasoro algebra** (`s1_virasoro_D26.py`). The central charges and intercept were extracted from explicit operator algebra, not input:
- **Fock space.** (L₂L₋₂ − L₋₂L₂ − 4L₀)|0⟩ = (c/2)|0⟩ gives c_matter = D for D = 1…5.
- **bc ghosts.**
  - c_ghost = −26, independent of cutoff 8 vs 12.
  - The general-λ check gives c = −2, −26, −74 for λ = 1, 2, 3, matching −3(2λ−1)² + 1.
  - L₀ c₁|0⟩ = −c₁|0⟩.
- **Verma module (symbolic c, h).**
  - Level-2 Kac determinant: 2h(16h² + 2h(c−5) + c).
  - (L₋₂ + γL₋₁²)|h⟩ is physical and null at h = −1 only for γ = 3/2, c = 26.
- **BRST.**
  - On (matter Verma) ⊗ (ghost Fock), Q² ∝ Σ_m[(c−26)(m³−m)/12 − 2m(A+1)] c₋ₘcₘ.
  - Explicit coefficients: −2(A+1) at m = 1, (c−8A−34)/2 at m = 2, 2(c−3A−29) at m = 3.
  - So Q² = 0 ⇔ c = 26, A = −1.
- **Light cone.** (D−2)ζ(−1)/2 = −1 ⇒ D = 26. ζ(−1) = −1/12 was cross-checked with a smooth cutoff, giving −0.083333329 at ε = 10⁻³.

**S2, Veneziano** (`s2_veneziano.py`), with α′ = 1, α(s) = 1 + s, tachyon externals:
- **Residues.** R_n(t) = (α(t)+1)_n/n!, degree n in t, exact for n ≤ 5.
- **Duality.** The s-channel pole sum equals Γ(−α_s)Γ(−α_t)/Γ(−α_s−α_t) to 15 digits.
- **t-channel pole from the s-channel tail.** (−α_t)·Σ_s → 1 as α_t → 0⁻ (1.0000643 at α_t = −10⁻⁴), matching the exact amplitude to ~10⁻⁶. The pole comes from the ζ(1+d) ~ 1/d tail of the s-channel tower.
- **Regge.** A/[Γ(−α_t)(−α_s)^{α_t}] − 1 falls like 1/|s|: 1.2×10⁻³, then 1.2×10⁻⁵, then 1.2×10⁻⁷ at |s| = 10², 10⁴, 10⁶.
- **Partial waves (exact in D).**
  - Level 2: c₂ = 25/(4(D−3)(D−1)) and c₀ = −(D−26)/(8(D−1)).
  - Level 3: c₁ = −(D−26)/(2(D−3)(D+1)).
  - Tightest bound over levels 1–8 is D ≤ 26. Per-level bounds: 132/5 (n = 4), 27 (n = 5), 26.589 (n = 6), 26.773 (n = 7), 27.156 (n = 8).
- **Low energy.** −ζ(2) is the first stringy coefficient.

**S3, one-loop σ-model β functions** (`s3_beta_functions.py`):
- RNC expansion g = δ − ⅓R xx, exact on S² and H².
- Lattice tadpole slope −0.07978 vs −1/(4π) = −0.07958, constant → ln 32/(4π).
- β^G = α′R_μν reproduces Polyakov's O(N) β(g) = −(N−2)g²/(2π) for N = 3, 4, 5.
- Witten cigar: R + 2∇∇Φ = 0, with β^Φ constant ⇒ k = α′/4 at D = 2.
- S³ with H flux: h = 2R².
- Linear dilaton: D = 26 − 6α′V².
- The dilaton Euler–Lagrange equation of the spacetime action is exactly −2√G e^{−2Φ}(G^{mn}β^G_mn − 4β^Φ)/α′.

---

## 5. Target computations (Phase 3) and their review

### 5.2 T1: unitarity corners (`phase3/targets/T1/T1_report.md`; review `phase4/RT1_verdicts.md`)

**Setup.**
- A = Γ(−α_s)Γ(−α_t)/Γ(−α_s−α_t − c), with α(s) = a + s.
- Residues are expanded in Gegenbauer C_j^{(D−3)/2}(cos θ), using exact rational arithmetic (python-flint) and rigorous root isolation.

**Results (OUR ATTEMPT, VERIFIED by RT1 with independent code):**
- **Open c = 0, on trajectory (m² = −a).**
  - Allowed region: −21/29 < a ≤ 1, 3 < D ≤ D_max(a), with D_max strictly increasing.
  - Exact low levels:
    - level 2: ρ₀ = [D(a²−2a) + 8a² + 14a + 4]/(8(D−1)), ρ₁ = (1−a)/(2(D−3));
    - level 3: ρ₀ = (1−a)(1+a)[(26a+30) − D(3−a)]/(48(D−1)).
  - Pieces of the boundary:
    - level 3, spin 0, giving D_max = (26a+30)/(3−a), up to a₁ = 0.3883687 (a root of 45a⁵ + 78a⁴ + 20a³ − 64a² − 61a + 30);
    - level 5, up to a₂ = 0.5157634;
    - a narrow level-4 window, up to a₃ = 0.5221548 (a root of 18a⁶ + 51a⁵ − 25a⁴ − 104a³ + 32a² + 72a − 32);
    - level 2, giving D_max = (8a² + 14a + 4)/(a(2−a)), up to a = 1.
  - Sample values: D_max(0) = 10; D_max = 32/5 at m² = 1/3. **The unique maximum is (a, D) = (1, 26).**
  - Verified to level 100 at 12–15 values of a and on grids of up to 917 points to level 60. Higher levels clear the boundary by ≥ 0.0209. This is not an all-level proof.
- **Sub-threshold sign rule (T1.2): VERIFIED, and a correction to the literature.**
  - For a pole below 4m², the CM momentum is imaginary. The physical condition is ρ_j = c_j/κ^j ≥ 0.
  - Published Coon/hypergeometric studies impose c_j ≥ 0 there, which wrongly excludes a healthy massive vector with c₁ ∝ e²(M² − 4m²) < 0.
- **Parity lemma (T1.3):** VERIFIED, with an elementary proof. It explains the a = 0 ↔ type-I identification.
- **Superstring-type c = 1 family:** a ≤ 0, D_max = (8a+10)/(1−a), corner (0, 10). VERIFIED.
- **Negative results.**
  - Off trajectory, D_max exceeds 26 for δ > δ* = 0.4022459 at a = 1. This is already the (r, m²) plane of 2302.12263/2409.09561.
  - Closed strings have no corner: VS stripped gives 23, the 4-dilaton amplitude gives 24.1478864864…, and closed bosonic tachyons give 57.
    - The closed tachyon residue is the square of the open one, so the spin-0 ghost that kills the open string above D = 26 cannot appear.
  - The CHR λ-families have no corner (VERIFIED-KNOWN, 2408.03362 Fig. 1).

**What T1 means for mission (1).** Four-point tree-level unitarity puts the critical open string at a corner *only* with an additional physical input (externals on the leading trajectory). It does not do so for closed strings. This is the same pattern as Front E's finding that every uniqueness proof has one UV-shaped input. **The "unitarity alone selects the string" hope FAILED** in the form we could test.

### 5.3 T2: stress test of the 2406.13751 KKLT dS candidates (`phase3/targets/T2/T2_report.md`; review `phase4/RT2_verdicts.md`)

**Data and baselines.**
- Data: paper Table 2, plus the authors' public repository output. The repository pickle was decoded with a no-execution stub reader into `data_extracted_data_dS.json`; the raw file is not redistributed.
- Baselines reproduced: KKLT σ_cr = 113.5891, V_AdS = −1.9966e−15, KPV (p/M)_crit = 0.08049, √g_sM_min = 6.801, g_sM² ≥ 46.25.

| Criterion | Ex1 | Ex2 | Ex3 | Ex4 | Ex5 | Review |
|---|---|---|---|---|---|---|
| g_sM² (bound ≳ 46) | 16.82 | 14.61 | 14.33 | 16.16 | 11.94 | VERIFIED (data) |
| Conifold minimum with anti-D3, warped EFT (c′=1.18, c″=1.7546) | none (×2.75) | none (×3.17) | none (×3.23) | none (×2.86) | none (×3.88) | **VERIFIED (conditional)**; holds even with the bound halved |
| Complex-structure sector inert under uplift (R_cs; V_cs/V_D̄3) | 9.4; 1.52 | 2.6; 0.56 | 1.5; 0.23 | 1.3; 0.14 | 2.1; 0.47 | VERIFIED (mechanism known: 2512.17095 §7) |
| KPV with α′² at g_sM ≈ 1 | uncontrolled | uncontrolled | uncontrolled | uncontrolled | uncontrolled | VERIFIED-KNOWN; robust: tension correction 5.4–10.7× |

**What T2 means for mission (2).**
- The best explicit dS candidates pass only in the authors' leading-order EFT.
- Under the warped-throat description of the conifold, which is the one debated in the literature, they fail. The verdict therefore hinges on an unresolved modelling question (the off-shell warp factor), not on numerics.
- This sharpens the open question. It does not close it.

### 5.4 T3: Λ selection head-to-head (`phase3/targets/T3/T3_report.md`; review `phase4/RT3_verdicts.md`)

T3_PLACEHOLDER

### 5.5 The celestial-amplitude error (Phase 1 discovery → Phase 4 confirmation)

**Setup** (1806.05688).
- F_I(s,u) = Γ(1−s)Γ(1−u)/Γ(1−s−u) (eq. 3.10).
- The celestial amplitude reduces to I(r, β) = ½∫₀^∞ w^{−β−1} F_I(rw, −w) dw (eq. 3.15), with s → s + i0.

**Findings.**
- **The pole-sum formula is wrong.** The paper's (3.30)/(3.31) reproduces only the s-channel discontinuity, I₊ − I₋ = −iπ r^β S(β), which holds to 10⁻³⁰. It misses the u-channel poles at w = −m. F grows exponentially in a wedge around the negative w axis, with half-angle arctan(K/π), K = r ln r − (r−1) ln(r−1), so the Hankel contour cannot be closed.
  - Example at r = 2, β = −1.3: direct I₊ = −0.0358979 + 0.1982757i; formula +0.1440557 + 0.1982757i.
  - The mismatch persists at the physical β = 0.7i, on every branch, iε choice and sign reading of the printed formula.
- **Corrected formula.** I_± = π/(2 sin πβ)·[(r e^{∓iπ})^β S(β) − U_r(β)], where U_r is the analytically continued alternating u-channel sum.
  - It matches direct integration to 10⁻¹⁵–10⁻³¹, using two different continuations (R1: a Lindelöf integral; RB4: a generating function).
  - It has the exact residue −½ at β = 0 for all r > 1, and no spurious poles at the integers β ≤ 1 other than 0, where the paper's formula has them.
- **Small-angle claim false.** For Re β < 1, I ~ ½Γ(−β)(ln r)^β at large r. This is O(1) on the physical line: at β = 0.7i, |I| ≈ 0.5 out to r = 10³², while the claimed r^{β−1} term is 10⁻³². The same claim appears verbatim in 2405.01643 eq. (3.7).
- **Minor slips.** The residue sign in (3.32) is wrong, and a sinh argument is missing a factor of π.
- **Novelty.** No erratum or correction was found. The ingredients appear in 2012.04208 (exponentially growing cross-channel residues) and 2106.11948 (the imaginary part equals the s-channel sum).

---

## 6. Adversarial review: process and all verdicts

| Reviewer | Scope | File |
|---|---|---|
| R1 | Fronts A, B; sanity S1 | `phase4/R1_verdicts.md` |
| R2 | Fronts C, D, E; sanity S2, S3 | `phase4/R2_verdicts.md` |
| RT1 | T1 | `phase4/RT1_verdicts.md` |
| RT2 | T2 | `phase4/RT2_verdicts.md` |
| RT3 | T3 | `phase4/RT3_verdicts.md` |
| RB4 | Independent second referee on the celestial error | `phase4/RB4_verdicts.md` |

Every reviewer wrote independent code before reading ours (`phase4/*_code/`). The verdicts are in the §1 scorecard.

**FAILED, with reasons:**
- **E1.** The SUSY nonlinear relations of 2601.11705 need a 6-scalar parity assumption (v2, footnote 3). The exclusion is also a corollary of 2610.02124. We cited the superseded v1.
- **D1.** The dS-window edges were log-grid points, not exact edges.
- **A1.** The v¹⁰/r¹⁹ coefficient is +182325/16384 in the stated (Minkowski) convention.
- **T2.7.** A flux step in K′ alone breaks the tadpole. The window was calibrated on the survivors. The Ξ distribution of the 396 is unknown. One flux choice (Ex4) contains 22 vacua. So the "yield/efficiency" comparison carries no information.
- **T3 failures:** see §5.4.

**Downgraded for novelty (VERIFIED-KNOWN):** A2, T1.7, T2.4, T2.6 (and T2.2's mechanism); T1.4 and T1.6 are overclaimed.

---

## 7. What is still open (relative to the mission)

**Mission (1): non-perturbative definition outside AdS.** Untouched at the level of the definition.
- What we have:
  - (i) Four-point unitarity does not by itself single out the string (T1). Any S-matrix-based definition needs extra input: higher-point factorization, a spectrum assumption, or monodromy.
  - (ii) One published celestial-holography formula is corrected (§5.5).
- Concrete next computations, from the Phase 1 handles:
  - the three-loop D0 v⁸/r²¹ test, with general N₁, N₂;
  - a Lorentzian IKKT "Pauli index" diagnostic;
  - the closed-SFT tachyon–dilaton coupling against −27/32;
  - O(λ) DSSYK vs dS₃;
  - celestial Virasoro–Shapiro, using the corrected method of §5.5.

**Mission (2): vacuum selection.**
- Existence of controlled dS is still open. T2 shows the best candidates fail a debated-but-standard warped-throat criterion.
- For Λ: see §5.4.
- No selection principle has been derived. The live options and what each hinges on are mapped in Front C and T3.

---

## 8. Reproducibility and process notes

**Running the code.**
- Every number labelled OUR ATTEMPT or reproduced comes from a script in this tree, with its output saved next to it as `*_output.txt` or `*.out`.
- Python 3.13 with sympy, numpy, scipy and mpmath; T1 also uses python-flint.

**Handling of fetched material.**
- arXiv was reachable only through the Exa fetch tool.
- Paper full texts fetched for reading are kept local and **are not committed**; see `.gitignore`. Earlier WIP commits on this branch briefly contained some of them; see the note in the session summary.
- The 2406.13751 data repository's pickle was treated as untrusted. It was never unpickled with executable loaders; a stub reader extracted plain numbers.

**Known limits of the review.**
- Reviewers could not view figure images; PDFs were blocked, and figure-based novelty judgements rest on captions and text.
- Several 2026 papers were checked only through abstracts or partial full text.
- Nothing here has been peer reviewed by humans.
