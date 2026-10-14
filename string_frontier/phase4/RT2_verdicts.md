# RT2 — Adversarial review of T2 (stress test of the 2406.13751 candidate de Sitter vacua)

**Reviewer:** RT2. I had not seen this work before.

**Under review:** `phase3/targets/T2/T2_report.md`, with its `code/` folder and outputs.

**Method.**
- I wrote my own code from the primary-source formulas, in `phase4/RT2_code/`:
  - `rt2_throat_checks.py` covers the data and claims 1–5.
  - `rt2_kpv.py` and `rt2_hsv_tip.py` cover claim 6.
  - `rt2_window.py` covers claim 7.
  - `*_output.txt` holds the outputs.
  - `rt2_literature_notes.txt` holds the literature and citation evidence.
- I read only the decoded JSON. The pickle was never loaded.
- I fetched the authors' `code/kahler_stabilisation.py` as text through Exa, to check how the uplift and the complex-structure (cs) F-term are normalised.
- I fetched the INSPIRE citation list of 2406.13751 (79 records), using `refersto:recid:2800646` through Exa.

## Verdict table

| # | Claim (T2) | Verdict | One-line reason | My script |
|---|---|---|---|---|
| D | ESTABLISHED data: Table 2 and g_sM² = 16.82, 14.61, 14.33, 16.15, 11.93 | **VERIFIED** (one unflagged discrepancy) | Repo dS values match Table 2 to <0.1%. g_sM² reproduced. The Ex3 SUSY point in the repo differs from paper eq. 5.39: \|z\| is 1.171e-7 vs 1.121e-7 (4.4%), and W0 is 0.02154 vs 0.0217. | rt2_throat_checks.py §A |
| 1 | z_cf stays on its F-flat relation. The residual moves ≤0.025. Δln g_s = −0.108…−0.016 | **VERIFIED** | Reproduced exactly: residual changes +0.002, −0.002, +0.008, +0.025, −0.017. It is the expected behaviour of the authors' unwarped z_cf metric. | §B |
| 2 | Uplift backreaction R_cs = 9.4, 2.6, 1.5, 1.3, 2.1; V_cs/V_D̄3 = 1.52, 0.56, 0.23, 0.14, 0.47 | **VERIFIED** (novelty overstated) | Reproduced to 3 digits. The normalisation was checked against the authors' code. The mechanism is already stated by the authors in 2512.17095 §7. Ex3 R_cs is 1.43 if the paper's eq. 5.39 is used. | §C |
| 3 | No conifold minimum in the warped EFT: bound missed ×2.75–3.88; c′ ≤ 0.30–0.43 needed | **VERIFIED** (conditional on the BDGL/DST off-shell warp-factor assumption, which T2 labels) | I re-derived β = 4πc′c″/(g_sM²) from BDGL eqs. 2.17 and 2.32 and scanned V(L) numerically: no critical point for any vacuum. The convention checks pass, including the n_cf = 2 throats, which are orientifold images. No reading of the conventions flips the verdict. | §D |
| 4 | m_ζ ≈ 5369, 12503, 5598, 612, 576 H ≈ 0.5 m_wKK; 23–853× the lightest Kähler mass; KK ambiguity ×23 | **VERIFIED-KNOWN** | Every number reproduced (the S³ KK estimate is my own derivation from eq. B.3). The qualitative point, that z_cf sits at the warped KK scale with no parametric separation, is in 2406.13751 fn. 44 together with Fig. 5, and in BKS 1902.07724 eq. 4.15. | §E |
| 5 | Lüst–Randall S/S0 ≈ 0.79–0.84 (crude) | **PLAUSIBLE** (the label is honest) | Arithmetic reproduced: 0.94^(46.25/g_sM²) gives 0.844–0.787. It cannot be checked beyond linear response; Lüst–Randall themselves say the potential far from S0 is undetermined. T2 marks it "crude/undetermined". | §F |
| 6 | KPV with α′²: factor 41–167, bracket up to 2e3, no metastable minimum, g_sM ≥ 21–33 for control; SV/HSV validation | **VERIFIED-KNOWN** | My independent HSV/SV implementation reproduces every number, including T2's failure to reproduce the HSV shark-fin tip (I get 16.3, HSV quote 20.36). "Control lost" is trivially true at g_sM ≈ 1 and already stated by the authors and by Schreyer 2402.13311 (not cited). 81–92% of the 41–167 comes from cot-terms of a sub-stringy NS5 for p = 1. | rt2_kpv.py, rt2_hsv_tip.py |
| 7 | dS window 1.6–4.0 vs K′ step ×54–280; window = 12–25% of a step; 74–114 of 396 expected vs 5 found ⇒ 4–7% "efficiency" | **FAILED** (the yield/efficiency part is ill-posed). The arithmetic is VERIFIED. | Every number is reproduced. But the K′ ladder breaks the tadpole, and the toy window is calibrated on the survivors and already fails for Ex1 and Ex2. The Ξ distribution and the number of the 396 actually attempted are unknown. Ex4 alone gives 22 dS vacua from one flux choice. The "4–7% efficiency" therefore carries no information, yet §9(v) lists it as a new result. | rt2_window.py |

## Details

### D. Data extracted (labelled ESTABLISHED in T2)

**dS point.** The repo dS-point values agree with 2406.13751 Table 2 and eqs. 5.13, 5.29, 5.42, 5.55, 5.69 to <0.13% in every entry.

**g_sM².** g_sM² = 16.818, 14.614, 14.329, 16.150, 11.928, so T2's table (16.15, 11.93) is correct. The values 16.16 and 11.94 in the review brief are rounding. At the SUSY point the values are 18.74, 15.24, 14.64, 16.41, 12.44.

**Unflagged discrepancy (Ex3 SUSY point).** Paper eq. 5.39 gives g_s = 0.0450, z = 1.121e-7, W0 = 0.0217, g_sM = 0.811. The repo has 0.04517, 1.1706e-7, 0.02154, 0.8131. All other SUSY points agree to ≤0.1%. The effect is small: Ex3's R_cs becomes 1.43 instead of 1.51, and its F-flat residual is −0.062 instead of −0.073. T2 should note which source it uses.

**Counts.**
- "396" appears only in the Fig. 3 caption (and again in the TASI Fig. 8.2 caption).
- "30 dS vacua in five compactifications" is correct: 1+1+1+22+5 (§5.4: 29 precursors give 22 dS vacua; §5.5: 59 KKLT points give 5).
- c_code = g_s²·c(3.32) is reproduced (ratio 1.00000).
- The authors' code (`potential_km`, `gradient_potential_km`) adds `c_uplift/cy_volume**(4/3) + gs**3*CSFterm/cy_volume**2` with the string-frame volume. So V_cs = CSFterm/V_E², as T2 says.

**Minor overstatement (§7).** T2 says "we re-derived their bookkeeping V_dS = V_F,Kähler + V_D̄3 + V_cs". This is circular: V_F,Kähler is computed as the remainder (`t2_heavy_conifold.py`, and my §C). It is a decomposition, not a re-derivation.

### 1. F-flat relation

**Reproduction.** My code gives residuals ln(2π|z|)+πK′/(g_sM) = (−0.006, 0.023, −0.073, 1.900, 0.121) at the SUSY point and (−0.004, 0.021, −0.066, 1.925, 0.104) at the dS point. The changes are +0.0022, −0.0015, +0.0075, +0.0251, −0.0171, as claimed.

**Points T2 should add.**
- **Ex4 does not satisfy the relation itself.** For Ex4 the relation (3.11) with Q_throat ≈ MK′ misses by e^{1.9} ≈ 6.7 in |z|, even though it is stable under the uplift. T2's text says z_cf "is slaved to (3.11)", which is only true up to an example-dependent constant.
- **Typo in the paper.** Paper eqs. 4.15 and 3.18 omit n_cf. The data follow eq. 3.11, and TASI 2512.17095 eq. 6.80 with n_cf = 2. T2 used the correct form.
- **Reason for the behaviour.** The behaviour follows from the authors' unwarped Kähler potential (2.34). This is a property of their EFT, not a defect in it.

### 2. Backreaction of the uplift on the dilaton and cs sector

**Reproduction.** My values:
- R_cs = 9.380, 2.636, 1.513, 1.307, 2.093.
- V_cs/V_D̄3 = 1.519, 0.558, 0.225, 0.138, 0.472.
- \|F_cs\|/\|W\| = 1.41, 1.15, 0.86, 0.78, 1.51.
- W0·V_E^{2/3} = 12.6, 6.1, 30.7, 30.4, 12.4.

The physics is simple. c ∝ g_s·exp(−(4π/3)K′/(g_sM)), so d ln c/d ln g_s ≈ 17–23, and a 2–11% drop of g_s cuts the uplift by 1.3–9.4×.

**Novelty.** The qualitative claim that "integrate out cs supersymmetrically, then uplift" fails is already in 2512.17095 (McAllister and Schachner, TASI lectures, Dec 2025, §7 near eqs. 7.55–7.56): "one complex structure direction remains almost massless … highly susceptible to backreaction … the usual two-step procedure … is not valid here". The paper itself also says the Ex1 vacuum breaks SUSY "in both the complex structure and Kähler sectors" (§5.1).
- T2's literature search missed this reference, although it appears in the INSPIRE citation list.
- The per-vacuum numbers appear to be new, and they are correct.
- T2's novelty item (i) must cite 2512.17095 §7.

### 3. BDGL conifold criterion (the convention audit requested)

**Re-derivation.** I re-derived the potential from BDGL eq. 2.17 (warped regime 3.1) and eq. 2.32:
- V ∝ |S|^{4/3}[L² + 4πc′c″/(g_sM²)], with L = ln(S/s_KS).
- Critical points sit at L = −3/4 ± √(9/16 − β), which is BDGL eq. 3.2.
- My numerical sanity check at g_sM² = 160 reproduces eq. 3.2 to 5 digits.

**Inputs.**
- c″ = 2^{1/3}/I(0) = 1.75464 (eq. 2.33), with I(0) = 0.71805.
- c′ = 1.18 (eq. 2.18).
- The bound is 46.255, matching Junghans eq. 2.9 (46.1 with c″ = 1.75).

**Results.**
- β = 1.547, 1.780, 1.816, 1.611, 2.181, against the requirement β ≤ 0.5625.
- The bound is missed by 2.750, 3.165, 3.228, 2.864, 3.878.
- A minimum would need c′_max = 0.429, 0.373, 0.366, 0.412, 0.304.
- A direct scan of V(L) on [−40, 5] finds 0 critical points for all five vacua. V rises monotonically in L, so the system runs to ζ → 0.
- The strong-warping condition holds: e^{−4A_IR}/t = 3.0e8, 1.1e9, 2.3e7, 9.5e4, 5.3e5.

**Convention traps checked. None flips the verdict.**
- **n_cf = 2.** With two independent throats and one anti-D3, the flux potential doubles relative to the anti-D3 term, the bound halves to 23.1, and the factors become 1.37–1.94. But 2406.13751 §4.3 requires that "the orientifold symmetry … exchange the two conifold singularities". So there is one physical throat, with the anti-D3 and its image in the covering space; eq. 2.29 counts an anti-D3 as 2 units. The single-throat bound applies unchanged. Even the halved bound would be violated.
- **g_s at the SUSY point instead of the dS point.** g_sM² = 12.4–18.7, still violated.
- **β normalisation.** BDGL's g_s(α′M)² and Im ρ factors cancel between eqs. 2.17 and 2.32, so β is convention-free within BDGL.
- **The anti-D3 "×N" convention.** N = 1.

**Label issue.** T2 labels the K and V formulas ESTABLISHED. The physical conclusion, however, rests on the off-shell e^{−4A} ∝ |Z|^{−4/3} assumption. Junghans himself says "it is not clear whether this assumption is justified" (2201.03572, below eq. 2.4), and Lüst–Randall dispute it. T2 does list this input in Table 4 and gives an "undetermined" Lüst–Randall row, which is acceptable. The bottom-line sentence in §0 should keep the qualifier "in the BDGL/DST warped EFT".

**Novelty.** T2 does not claim novelty here. No paper found applies the bound to these vacua; the 2512.17095 lectures cite BDGL only for the W0 bound.

### 4. Conifold mass vs the Kähler spectrum

**Reproduction.**
- From Junghans's A at c″ → 0, matched to the repo's V_D̄3: m_ζ/H = 5369, 12503, 5598, 612, 576.
- From B.18: m_wKK/H = 10734, 25888, 11647, 1237, 1254. The ratio m_ζ/m_wKK = 0.500, 0.483, 0.481, 0.495, 0.459.
- I checked analytically why the ratio is about 0.5: m_ζ²/m_wKK² = 0.0610·√(g_sM²).
- The fraction of Kähler masses above m_ζ is 0.58, 0.47, 0.67, 0.83, 0.68.
- m_ζ/m_Kähler,min = 623, 853, 349, 23, 60.

**My own KK estimate.** I derived an S³ (l = 1) KK estimate from the metric of eq. B.3, with R from eqs. B.6 and B.7: m² = 3e^{2A_IR}/(4πtR²), with R = 0.29–0.31 ℓ_s in Einstein frame. It is 21–23× the B.18 value. This confirms T2's factor of ~23, which scales as (g_sM²)^{1/4}.

**Stronger statement that needs no c′ or γ0.** By the authors' own B.18, 40, 15, 49, 70, 47% of the Kähler moduli are heavier than m_wKK. This is visible in their Fig. 5, where Ex1's Kähler masses reach 2.7e4 H against m_wKK ≈ 1e4 H.

**Why VERIFIED-KNOWN.** The headline, that the conifold modulus sits inside the Kähler spectrum at about the warped KK scale with no parametric separation, follows from 2406.13751 fn. 44 together with Fig. 5, and from BKS 1902.07724 eq. 4.15 ("edge of reliability"). Only the per-vacuum numbers are new.

### 5. Lüst–Randall extrapolation

Eq. 7.5 of 2206.04708 (S/S0 ≈ 0.94 at g_sM² ≈ 46) and eq. 5.26 (exponent 0.44) are quoted correctly. My extrapolations:
- ln-linear: 0.844, 0.822, 0.819, 0.838, 0.787.
- S-linear: 0.835–0.767.

Lüst–Randall state that the potential far from S0 is not determined, because the third constraint (7.3) is not imposed. T2 labels the extrapolation "crude", "not reliable beyond linear order" and "undetermined" in Table 4. This is honest.

### 6. KPV with α′² corrections

**My implementation.**
- KPV potential: HSV eq. 2.6.
- HSV correction: eqs. 2.11 and 2.12, with c1 = 8.825 and c2 = 1.891.
- SV correction: eq. 3.5. I checked the LaTeX on ar5iv; the F₂ terms contain (πp/M − (ψ − sin2ψ/2)), so T2's reading is right. Constants are from SV Table 1.

**My results.**
- (p/M)_crit for KPV = 0.08049.
- SV zero-energy line: p/M = 0.00847 at g_sM = 10 and 0.00407 at g_sM = 20, against SV eq. 4.1's 0.00835 and 0.00392.
- The first V_up > 0 metastable minimum appears at g_sM = 3.6, at p/M ≈ 0.024. SV quote 3.6 and 0.025.
- HSV Fig. 4 asymptotes: 5.367 and 2.104 (HSV: 5.238 and 2.086).
- HSV shark-fin tip from a 2D scan: (0.0815, 16.3), against HSV's (0.0826, 20.363). The lower edge sits about 0.005 below HSV eq. 2.13.

So T2's unresolved HSV mismatch is **independently reproduced from the published formula**. It is most likely a difference in HSV's minimum criterion, not a T2 bug. It does not affect the g_sM ≈ 1 statements.

**At the five vacua, p = 1.**
- ψ_KPV = 0.490, 0.490, 0.426, 0.378, 0.490.
- Curvature factor f/(g_sM)² = 41.1, 54.4, 113.5, 167.2, 81.7. The constant c1/(g_sM)² alone is 8.0–15.9.
- The SV DBI bracket is 92–1496.
- HSV gives no minimum below π/2 and a maximum at ψ ≈ 0.74–0.75. SV's only minimum is at ψ ≈ 1.69, with V ≈ −72 to −157 V_D̄3.
- Control needs g_sM ≥ 21.3 (p/M = 1/16) or ≥ 33.0 (1/20). Everything matches T2.

**Is the "loss of control" an artifact? Answer: it is trivially true, and the "no minimum" result carries no information.**
- **Expansion parameter is not small.** At g_sM ≈ 0.75–1.05, R_S³² = b0²g_sM α′ = 0.70–0.98 α′. The simplest published α′² coefficient is the anti-D3 tension correction 1 − 5.924/(g_sM)² (HSV eq. 2.15, and the authors' own 2512.17095 eq. 6.122). Its correction term is 5.4–10.7, so the truncated tension changes sign below g_sM = 2.43.
- **Already known.** That the expansion is not controlled at g_sM ~ 1 is stated in:
  - 2406.13751 §3.3 ("the series of α′ corrections is not controlled when the leading correction starts to matter");
  - 2512.17095 below eq. 6.124 ("at small g_sM the known corrections are found to be significant");
  - HSV (g_sM > 20);
  - Schreyer 2402.13311 (g_sM² ~ O(100), and (g_sM)² > 4.91 for control). T2 does not cite this paper.
- **The 41–167 factor is mostly an NS5 artefact.** 81–92% of it comes from the cot²ψ and cot⁴ψ extrinsic-curvature terms of an NS5 whose radius² is 0.10–0.20 α′. For p = 1 the NS5 description itself is not valid: HSV §2.3 says a single anti-D3 "would not polarize into an NS5-brane", and 2512.17095 §6.5.1 says p = 1 is not captured by KPV's Myers analysis. The robust headline number is 5.4–10.7 (c/(g_sM)²), not 41–167.
- **What T2 gets right.** T2 does say that the absence of a truncated minimum "is not evidence of instability", which is correct.
- **What T2 should drop or fix.** The §0 bottom line ("Control is lost, and neither truncated potential has a metastable minimum") and Table 3's "V ≈ −72 V_D̄3" values report features of a divergent truncation. They should be dropped or explicitly marked meaningless. Novelty item (iii) is VERIFIED-KNOWN: the numbers are new, the content is not.

### 7. dS window vs flux steps, and expected yield — FAILED (ill-posed and overclaimed)

**Reproduction.** My code reproduces everything:
- KKLT baseline σ = 113.5891. For n = 2, Ξ ∈ [1.007, 3.773]; for n = 3, Ξ ∈ [1.016, 2.687], i.e. D ∈ [2.97, 7.86]e-9.
- Calibrated x = 2.49–12.45. Windows [1.099, 1.762] … [1.006, 4.040], i.e. factors 1.60–4.02.
- Step factor exp((4π/3)/g_sM) = 54, 98, 193, 179, 276. w/step = 0.119–0.247.
- Ξ_prec from eq. B.16 = 3.93, 1.89, 1.52, 1.76, 3.16; the exact values are 4.65, 2.18, 1.73, 2.00, 3.65.
- Expected counts: 74 (mean w) and 114 (x = 11.36). The per-vacuum w spans 41–120, so "74–114" is a selected sub-range.
- T2's use of n_cf = 2 in the step is correct (see claim 1).

**Why the comparison is not meaningful.**
1. **The K′ ladder is not physical.** It holds the bulk fluxes fixed. But eq. 4.16 (Q_bulk = Q_O + 2 − MK′) means a unit K′ step changes the tadpole by M = 16–20, so the ladder is not a family of consistent vacua. T2 concedes this, which makes "window = 12–25% of a step" irrelevant to whether the window can be hit.
2. **The toy is calibrated on the survivors and fails on them.** The window comes from a single-modulus KKLT toy matched to the lightest Kähler mass of the five vacua that survived. Kähler data do not exist for the other 391. The toy already puts Ex1 and Ex2 outside their own windows, although they are dS vacua; T2 attributes this to cs relaxation, which the toy lacks. A model falsified on 2 of its 5 calibration points cannot predict the yield.
3. **The base rate is unknown.** The ln Ξ distribution of the 396 is not published; log-uniform is an assumption. How many of the 396 had a KKLT point, converged, or passed the Hessian and control cuts is also unknown. Ex1–3 each have a single KKLT point, so this filter is evidently strong.
4. **Wrong counting unit.** The window applies per AdS precursor, not per flux choice: Ex4 yields 22 dS vacua from 29 precursors of one flux choice. Comparing "5 flux choices" with "74–114 window hits" mixes units.

**Conclusion.** The "implied efficiency of all other criteria, 4–7%" is just 5 divided by an unsupported expectation. It carries no information about the dS window or the vacua. T2 hedges in §5 but promotes it to a bottom-line bullet (§0) and a novelty item (§9(v)). Both should be removed, keeping only the reproducible arithmetic (window factor, step factor, Ξ_prec outside the paper's "1 ≲ Ξ ≲ 3" for Ex1 at 3.9 and Ex5 at 3.2).

## Novelty and citation audit

**Search.** I fetched the INSPIRE citing list (79 records) using `refersto:recid:2800646` through Exa. T2's `refersto:arxiv:` query returns unrelated records, which explains its "unusable" result. I also ran Exa searches for critiques by Bena, Graña, Van Riet, Junghans, Lüst–Randall and Schreyer.

**What I found.**
- No paper evaluates the BDGL bound, the conifold mass, the α′²-KPV numbers, or the window for the five vacua specifically.
- Three relevant works are missing from T2:
  - **2512.17095** (McAllister and Schachner, TASI lectures, Dec 2025). Its §7 already says the two-step procedure fails because of the light PFV direction (claim 2), and its eq. 6.122 and following text say the α′² corrections are significant at small g_sM (claim 6).
  - **2402.13311** (Schreyer, JHEP 07 (2024) 075): g_sM² ~ O(100) and (g_sM)² > 4.91 for control (claim 6).
  - The repo-code normalisation (V_cs = CSFterm/V_E²) is confirmed from the authors' `kahler_stabilisation.py`.
- **Novelty items, item by item.**
  - (i) Numbers new; mechanism known from 2512.17095 §7.
  - (ii) Known from 2406.13751 fn. 44 together with Fig. 5, and from BKS 1902.07724.
  - (iii) Numbers new; conclusion known.
  - (iv) A failure of the toy, not a finding about the vacua.
  - (v) Ill-posed.
- **ID spot-checks** (abstract pages fetched): 1902.07724 BKS ✓, 2201.03572 Junghans ✓, 2512.17995 HLSS ✓, 2402.13311 ✓. 2208.02826, 2212.07437, 2206.04708 and 2505.00149 are consistent with search results and ar5iv. BKS eq. 4.15 and Lüst–Randall eqs. 5.26 and 7.5 say what T2 claims. **No CITATION PROBLEM found.**

## Label check

- The ESTABLISHED items are correctly sourced (BDGL, Junghans, B.x equations, Table 2).
- The physical conclusion of §3b depends on an assumption its own source flags as unclear (Junghans), so it is conditional. T2's Table 4 carries this, but §0 states it flatly.
- §7 "re-derived their bookkeeping" overstates a remainder decomposition.
- §0 and §9(v) present the yield comparison as a result rather than as an assumption-laden aside.
- The Lüst–Randall row is honestly labelled.
- No CONJECTURE is presented as ESTABLISHED, apart from the mild §3b/§0 issue above.

## Requested corrections to T2

1. Drop the "74–114 vs 5 / 4–7% efficiency" item from §0 and §9(v), or label it explicitly as uninformative.
2. In claim 6, lead with the anti-D3 tension correction 5.924/(g_sM)² = 5.4–10.7, which changes sign below g_sM = 2.43. Explain that the 41–167 comes from an NS5 description that is invalid for p = 1. Remove the face-value "V ≈ −72 V_D̄3" values or mark them as meaningless. Cite 2402.13311 and 2512.17095 eq. 6.122.
3. Cite 2512.17095 §7 for the failure of the two-step procedure (novelty item i).
4. Note the Ex3 SUSY-point mismatch with eq. 5.39, and that Ex4 misses the (3.11)/(MK′) relation by e^{1.9}.
5. Add the orientifold-image argument (2406.13751 §4.3) showing that the single-throat BDGL bound applies for n_cf = 2.
6. Add the simpler, input-free version of claim 4: by the authors' own B.18, 15–70% of the Kähler moduli lie above m_wKK.
