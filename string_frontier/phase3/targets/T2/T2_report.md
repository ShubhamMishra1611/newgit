# T2 — Stress test of the KKLT de Sitter candidates of arXiv:2406.13751

Labels:
- **ESTABLISHED**: stated in a primary source we fetched this session, cited by arXiv ID and equation number.
- **CONJECTURE**: claimed or disputed in the literature.
- **OUR ATTEMPT**: computed here; the script is named.

All files are under `phase3/targets/T2/`: code in `code/`, fetched texts in `sources/`, outputs `*_output.txt`.

IDs fetched and checked this session (title and authors match):
- 2406.13751 (v2) and 2505.00149
- 1809.06861, 2201.03572, 2206.04708
- 2208.02826, 2212.07437
- 1902.07724, 2512.17995

hep-th/0301240 and hep-th/0112197 were verified in Phase 1 and not re-fetched here.

## 0. Bottom line

The five vacua are genuine critical points of the authors' leading-order EFT. Running their own repository data through the outputs below gives these results.

- **Conifold modulus.** It is pinned only because their EFT uses the *unwarped* metric for z_cf. In the warped-throat EFT (BDGL/Junghans) it has **no minimum** once the anti-D3 is added (bound missed by factors 2.75–3.9), and without the anti-D3 its mass sits **inside the Kähler spectrum**.
- **Uplift backreaction.** The dilaton relaxes and cuts the would-be uplift by factors 1.3–9.4; the complex-structure F-term energy is 14–152% of the anti-D3 energy.
- **Anti-D3 (KPV).** At g_sM ≈ 0.75–1.05 the known α′² corrections are 40–170× the leading term at the would-be KPV minimum. Control is lost, and neither truncated potential has a metastable minimum.
- **dS window.** One K′ step moves the uplift coefficient by a factor 54–280; the window spans 1.6–4. With log-uniform Ξ (assumed), 74–114 of the 396 aligned PFVs would land in it; 5 dS flux choices are reported.

## 1. Data extracted (sources)

**Paper sources (2406.13751 v2, HTML).**
- Table 2.
- §5, eqs. 5.6–5.73: fluxes, SUSY-point and dS-point g_s, W0, z_cf, g_sM, V_E, V_dS, m_min/H_dS.
- Table 1: the counts 202,073 → … → 33,371 anti-D3 PFVs → 30 dS vacua.
- Fig. 3 caption: 396 points with 0.1 ≤ Ξ ≤ 10, g_s < 0.4, W0 < 0.1. This is the only place "396" appears.
- eq. 3.30 (M/p ≳ 12), eq. 3.31/3.32 (c = η z^{4/3}/(g_sM²Ṽ^{2/3}), η ≈ 2.6727), eq. 3.34/B.16 (Ξ, ζ ≈ 114.037).
- B.9 (γ ≈ 39.46), B.18 (m_wKK), (3.11) (|z| = (1/2π)exp(−2πQ_throat/(g_sM²n_cf)), n_cf = 2).
- Q_O = 2 + h¹¹ + h²¹ = 160 (Ex1–3) or 100 (Ex4–5), with Q_flux = Q_O + 2 and p = 1.

**Repository sources** (AndreasSchachner/kklt_de_sitter_vacua, git tree 83c2eb2…).
- `data/dS_examples/data.p` was decoded with `code/safe_read_pickle.py`, a stub unpickler that executes no code from the file. Output: `data_extracted_data_dS.json`, `extract_repo_data_output.txt`.
- Not in the paper text but in the repo:
  - mirror volume Ṽ (`volTilde`);
  - uplift coefficient (`c uplift`);
  - complex-structure F-term term (`CSFterm`);
  - V_AdS;
  - the full Kähler and complex-structure mass spectra in Hubble and AdS units.
- The uplift as coded (`code/kahler_stabilisation.py`) is V = V_F + c_code/𝒱_s^{4/3} + g_s³·CSFterm/𝒱_s², with 𝒱_s = g_s^{3/2}V_E.
- We checked that c_code = g_s²·c(3.32) to 10⁻⁴ for all five.

**Published nowhere:**
- the Ξ values (or any data) for the 396 aligned PFVs;
- how many of the 396 were taken through Kähler stabilisation.

**Table 1.** dS-point data; values marked † are repo-only.

| Ex | h¹¹,h²¹ | M | K′ | g_s | W0 | \|z_cf\| | g_sM | g_sM² | V_E | Ṽ † | V_dS | V_AdS † | m_min/H |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 150,8 | 16 | 26/5 | 0.0657 | 0.0115 | 2.82e-8 | 1.051 | 16.82 | 36459 | 435.1 | 1.937e-19 | −4.527e-19 | 8.62 |
| 2 | 150,8 | 16 | 93/19 | 0.0571 | 0.00490 | 7.93e-9 | 0.913 | 14.61 | 44312 | 334.0 | 1.692e-20 | −4.899e-20 | 14.66 |
| 3 | 150,8 | 18 | 40/11 | 0.0442 | 0.0222 | 8.73e-8 | 0.796 | 14.33 | 51265 | 155.8 | 4.983e-19 | −1.185e-18 | 16.03 |
| 4 | 93,5 | 20 | 17/5 | 0.0404 | 0.0539 | 1.97e-6 | 0.808 | 16.15 | 13401 | 4.54 | 2.341e-15 | −3.148e-15 | 26.16 |
| 5 | 93,5 | 16 | 29/10 | 0.0466 | 0.0304 | 8.70e-7 | 0.746 | 11.93 | 8266 | 11.43 | 2.113e-15 | −1.359e-15 | 9.54 |

Throat warp factor per throat, e^{−8π(K′/2)/(3g_sM)}, with the effective K per throat equal to K′/n_cf: 1.0e-9, 1.8e-10, 4.9e-9, 2.2e-8, 8.4e-8. OUR ATTEMPT; this is not published.

## 2. Baselines (`code/t2_baselines.py`, `t2_baselines_output.txt`)

All reproduced. ESTABLISHED formulas; OUR ATTEMPT for the numbers.
- **KKLT** (hep-th/0301240 eqs. 13–14): σ_cr = 113.5891, V_AdS = −1.9966×10⁻¹⁵.
- **KPV** (hep-th/0112197 eqs. 40, 45, b0² = 0.93266): (p/M)_crit = 0.08049.
- **BDGL** (1809.06861 eqs. 2.10, 2.18, 2.33, 3.3):
  - I(0) = 0.71805 and c″ = 2^{1/3}/I(0) = 1.7546;
  - √g_s M_min = 6.801;
  - g_sM² ≥ 64πc′c″/9 = 46.25. This becomes 46.13 with c″ = 1.75, which is 2201.03572 eq. 2.9's 46.1.
- c″ = 2/Σ with Σ = 2^{2/3}I(0) = 1.13983, the same Σ as 2406.13751 eq. B.8. So both papers use the same KS normalisation.

## 3. Task 2 — heavy-conifold self-consistency (`code/t2_heavy_conifold.py`)

### 3a. In the authors' own EFT

**OUR ATTEMPT (repo data).** z_cf is slaved to its F-flat relation (3.11).
- The residual ln(2π|z|) + πK′/(g_sM) changes between the SUSY point and the dS point by only +0.002, −0.002, +0.008, +0.025, −0.017.
- ln|z| itself drops by 1.59, 0.70, 0.29, 0.19, 0.52.
- So z_cf moves only because g_s moves: Δln g_s = −0.108, −0.042, −0.021, −0.016, −0.042.
- Reason: their Kähler potential (2.34) is −log(−i∫Ω∧Ω̄), with no warping term. With that metric the formal conifold mass is exponentially large (cf. BDGL eq. 2.25, which scales as 1/s²).
- The paper's only statement about this mode is footnote 44 (§5.1): z_cf "gets replaced by the redshifted tower … with masses of order m_wKK". We found no mention of the BDGL bound in the fetched text, which covers the main body and Appendices A–C.

**OUR ATTEMPT (repo data).** The uplift backreacts strongly on the dilaton and complex-structure sector.
- With |z| slaved and Ṽ ∝ g_s⁻³ (Ṽg_s³ is constant to <1%, checked), V_D̄3 ∝ g_s·exp(−(4π/3)K′/(g_sM)).
- Hence d ln V_D̄3/d ln g_s = 21.7, 23.5, 20.1, 18.6, 17.3.
- Factor by which the uplift at the SUSY complex-structure point exceeds the uplift at the dS complex-structure point (V_E fixed): R_cs = **9.38, 2.64, 1.51, 1.31, 2.09**.
- The coded complex-structure F-term energy CSFterm/V_E² divided by V_D̄3: **1.52, 0.56, 0.23, 0.14, 0.47**.
- In Ex1 this complex-structure F-term energy (2.99e-19) exceeds the anti-D3 energy (1.97e-19).
- The vacua are therefore partly F-term uplifted: |F_cs|/|W| ≈ 1.41, 1.15, 0.86, 0.78, 1.51, defined via V_cs/|V_AdS| = |F|²/3|W|².
- The authors minimise the full potential, so their solution is not in error. But "integrate out the complex-structure sector supersymmetrically, then uplift" fails here at O(1)–O(10).
- This links the vacua to 2512.17995 §5.5. Its eq. 5.67 asks for W0𝒱^{2/3} ≫ 1, citing these examples as the closest; we find W0V_E^{2/3} = 12.6, 6.1, 30.7, 30.4, 12.4 (frame of 𝒱 not checked).

### 3b. Warped-throat EFT

**ESTABLISHED.** Sources: BDGL 1809.06861 eqs. 2.17, 2.22, 2.24, 2.32, 3.2, 3.3 and Junghans 2201.03572 eqs. 2.1, 2.2, 2.5, 2.8, 2.9, 2.13, 2.16.
- K ⊃ c′ξ′|Z|^{2/3}/𝒱^{2/3}, with ξ′ = 9g_sM².
- V_ζ = ζ^{4/3}/(2π²c′𝒱^{4/3})·[πc′c″/(g_sM²) + (1/4)ln²(ζ/ζ0)].
- A critical point exists only if 9 − 64πc′c″/(g_sM²) ≥ 0.

**OUR ATTEMPT.** Applied to the five vacua.
- β = 4πc′c″/(g_sM²) = 1.55, 1.78, 1.82, 1.61, 2.18. A minimum needs β ≤ 0.5625.
- The bound is missed by factors 46.25/(g_sM²) = **2.75, 3.17, 3.23, 2.86, 3.88**.
- With the anti-D3, V(L) ∝ e^{4L/3}(L² + β) has no critical point: the gradient at the SUSY point is (4/3)β > 0, and the potential falls monotonically towards ζ → 0. **The mass with the anti-D3 is undefined; there is no minimum.**
- A minimum would need c′ ≤ 0.43, 0.37, 0.37, 0.41, 0.30, against the DST value 1.18 (BDGL eq. 2.18).
- The regime of BDGL eq. 3.1 holds strongly. e^{−4A_IR}/t from 2406.13751 eqs. B.3 and B.9 is 3e8, 1.1e9, 2.3e7, 9.5e4, 5.3e5, so the warped term dominates the unwarped log term by many orders of magnitude.

**OUR ATTEMPT — masses without the anti-D3.**
- Junghans' A at c″ → 0 is m_ζ² = ζ^{2/3}/(8π²c′²g_sM²𝒱^{2/3}).
- We fixed ζ/𝒱 by matching his V_D̄3 = c″ζ^{4/3}/(2πg_sM²𝒱^{4/3}) to the repo's anti-D3 energy. This gives m_ζ² = √(2πV_D̄3/(c″g_sM²))/(8π²c′²), independent of Ṽ and Λ0.

**Table 2.** Mass scales in units of H_dS.

| Ex | m_ζ/H | m_wKK/H (B.18) | m_KK(S³, l=1)/H | Kähler m/H: min, median, max | share of Kähler modes heavier than m_ζ | m_ζ / m_Kähler,min | lightest bulk cs+τ m/H |
|---|---|---|---|---|---|---|---|
| 1 | 5369 | 10734 | 2.4e5 | 8.6, 7362, 26894 | 0.58 | 623 | 21.6 |
| 2 | 12503 | 25888 | 5.7e5 | 14.7, 11076, 40577 | 0.47 | 853 | 48.5 |
| 3 | 5598 | 11647 | 2.6e5 | 16.0, 11339, 42646 | 0.67 | 349 | 57.2 |
| 4 | 612 | 1237 | 2.8e4 | 26.2, 2566, 18890 | 0.83 | 23 | 48.0 |
| 5 | 576 | 1254 | 2.6e4 | 9.5, 1009, 6547 | 0.68 | 60 | 23.4 |

Points about Table 2:
- **Check:** our m_wKK(B.18) for Ex1 is 1.07×10⁴ H. Fig. 5's caption gives "≈10⁴ H_dS".
- m_ζ/m_wKK(B.18) = 0.46–0.50. With our S³-radius KK estimate the ratio is 0.022.
- The two KK definitions differ by a factor of about 23 (B.18 implies a KK length √g_s M ℓ_s; the tip S³ radius is 0.29–0.31 ℓ_s, Einstein frame). Unresolved.
- Either way the separation is not parametric: BKS 1902.07724 eq. 4.15 (ESTABLISHED) gives m_KK²/m_Z² = c·f_n²·y_UV².
- **Conclusion:** in the warped EFT the conifold modulus lies inside the Kähler spectrum: 47–83% of the Kähler modes are heavier. It is 23–850× heavier than the lightest Kähler mode.
- Dependence on O(1) inputs: m_ζ ∝ c′⁻¹c″^{−1/4}γ0^{1/4}, where γ0 = 1 is Junghans' gluing constant. β ∝ c′c″.

**CONJECTURE (Lüst–Randall 2206.04708).**
- Imposing the constraint equations gives e^{−4A}(0,S) ∝ (S0/S)^{0.44} (eq. 5.26), not ∝ |S|^{−4/3}.
- At g_sM² ≈ 46 the shift is S/S0 ≈ 0.94 (eq. 7.5).
- Fig. 5's caption says the corrected potential "always has a minimum, irrespective of the value of g_sM²". But the text says the potential far from S0 at g_sM² ≪ 46 is not determined, because a third constraint is not imposed.

**OUR ATTEMPT (crude).** A linear-response extrapolation of eq. 7.5, using V_D̄3/V_flux ∝ 1/(g_sM²), gives S/S0 ≈ 0.84, 0.82, 0.82, 0.84, 0.79. It is not reliable beyond linear order.

## 4. Task 3 — KPV with α′² corrections at g_sM ≈ 1 (`code/t2_kpv_alpha.py`)

**Potentials implemented**, with F = πp/M − ψ + ½sin2ψ:
- KPV (HSV 2208.02826 eq. 2.6).
- HSV curvature correction (eqs. 2.11–2.12): −S(c1 + 2c2cot²ψ + c2cot⁴ψ)/(g_sM)².
- Schreyer–Venken full known α′² set (2212.07437 eq. 3.5 with A.16–A.23, c1…c8 from their Table 1). We read the text-extracted "πp/M − ψ − sin2ψ/2" as πp/M − (ψ − ½sin2ψ).

**Validation (OUR ATTEMPT).**
- SV zero-energy line: we get p/M = 0.00847 at g_sM = 10 and 0.00407 at 20. SV's fit (eq. 4.1) gives 0.00835 and 0.00392. Agreement is 1–8% for g_sM = 4–20.
- First positive-energy minimum: g_sM = 3.6 at p/M ≈ 0.024. SV quote 3.6 and 0.025.
- HSV Fig. 4 asymptotes: we get 5.37 and 2.10; HSV state 5.238 and 2.086.
- **Mismatch:** our HSV "shark-fin" tip is at (p/M, g_sM) = (0.0815, 16.2), against HSV's (0.0826, 20.36). Our lower edge also sits about 0.005 lower in p/M than their eq. 2.13 fit. Unresolved; we suspect a different criterion for very shallow minima. This does not affect the g_sM ≈ 1 conclusions below.

**OUR ATTEMPT — at the candidate parameters** (p = 1, M = 16–20, g_sM ∈ [0.7, 1.1]):
- The tip S³ radius² is b0²g_sM = 0.70–0.98 α′.
- At the leading-order KPV minimum (ψ_KPV = 0.49, 0.43, 0.38 for M = 16, 18, 20), the NS5 radius² is 0.10–0.20 α′.

**Table 3.** Corrections at ψ_KPV for the five vacua.

| Ex | curvature factor (c1+2c2cot²+c2cot⁴)/(g_sM)² | c1/(g_sM)² alone | HSV truncation | SV truncation |
|---|---|---|---|---|
| 1 | 41 | 8.0 | no minimum; one maximum at ψ = 0.75; runaway to ψ→0 | only minimum at ψ = 1.69 > π/2 (V ≈ −72 V_D̄3) |
| 2 | 55 | 10.6 | no minimum | ψ = 1.69 (−96 V_D̄3) |
| 3 | 114 | 13.9 | no minimum | ψ = 1.69 (−144 V_D̄3) |
| 4 | 167 | 13.5 | no minimum | ψ = 1.69 (−157 V_D̄3) |
| 5 | 82 | 15.9 | no minimum | ψ = 1.69 (−144 V_D̄3) |

Across the whole grid, the full SV DBI bracket at ψ_KPV is 83–1990 and |ΔV/V_KPV| is 80–1300.

- **Result:** the corrections are O(10²–10³), not small. Perturbative control of the anti-brane sector is lost.
- The truncated potentials taken at face value contain no metastable anti-D3/NS5 state below the equator. That is not evidence of instability, because they are outside their validity range. They also give no support for metastability.
- Keeping the curvature factor ≤ 0.1 at the KPV minimum needs g_sM ≥ 21 (p/M = 1/16) or ≥ 33 (1/20). This agrees with HSV's "g_sM > 20" (ESTABLISHED, 2208.02826 abstract).
- 2406.13751 §3.3 itself says these results "cannot rule out" metastability at g_sM ~ 1, and their §6 calls this "a particularly acute risk".

## 5. Task 4 — dS window vs flux discreteness (`code/t2_dS_window.py`)

**OUR ATTEMPT — KKLT single-modulus toy** with uplift D/σⁿ.
- Baseline (x = aσ = 11.36):
  - n = 2: Ξ ∈ [1.007, 3.773].
  - n = 3: Ξ ∈ [1.016, 2.687], i.e. D ∈ [2.97, 7.86]×10⁻⁹ (Phase-1 grid: [2.99, 7.71]×10⁻⁹).
- For n = 2 the window width w = ln(Ξ_max/Ξ_min) grows with x: 0.47 at x = 2.5, 1.32 at 11.36, 2.13 at 30.
- The toy gives m·ℓ_AdS ≈ 2x + 2.4.
- **Calibration** to each vacuum's lightest Kähler mass in AdS units (repo: 7.33, 9.78, 11.50, 18.43, 27.36) gives these windows:
  - [1.10, 1.76], [1.05, 1.98], [1.04, 2.16], [1.01, 2.96], [1.01, 4.04].
  - So w = 0.47–1.39.

**Precursor alignment Ξ.**
- Exact, from repo V_AdS and Ṽ: 4.65, 2.18, 1.73, 2.00, 3.65.
- Via the paper's eq. B.16: 3.93, 1.89, 1.52, 1.76, 3.16.
- Ex1 and Ex2 lie **outside** the calibrated toy window. Ex1, and Ex5 marginally, also exceed the paper's stated range "1 ≲ Ξ ≲ 3" (App. B). Ex1 and Ex2 exist only because of the relaxation in §3a (Ξ/R_cs = 0.50, 0.83, 1.14, 1.53, 1.74).

**Discreteness.**
- One unit of K′ (n_cf = 2, eq. 3.11) changes ln Ξ by (4π/3)/(g_sM) = 3.99, 4.59, 5.26, 5.19, 5.62.
- So c changes by a factor 54–280 per step, against a window factor of 1.6–4.
- The window covers 12–25% of one K′ step (w/step). Along a K′ ladder with the bulk sector fixed, at most one value fits, with probability 0.12–0.25.
- K′ steps also violate the tadpole unless the bulk fluxes compensate.
- Landing in the window therefore relies on the bulk-flux degrees of freedom (W0, g_s, Ṽ). Across the ensemble, 396 points over 4.6 e-folds give a mean spacing of 0.012 in ln Ξ, far below the window.

**Yield.**
- Assume ln Ξ is uniform over [0.1, 10] for the 396 points. This is unverified, since their Ξ values are unpublished.
- Then P(window) = w/ln 100 = 0.19–0.29, so 74–114 are expected and 5 dS flux choices are reported. The implied efficiency of all the other criteria is **4–7%**.
- The ratio 5/396 = 1.26% is not an established success rate, because the paper does not say how many of the 396 were tried.
- Multiplicity: Ex4 alone gives 22 dS vacua from 29 AdS precursors.

## 6. Task 5 — literature (2024–2026)

We ran three Exa searches (critiques of the candidate vacua; conifold mass in explicit KKLT; anti-D3 at g_sM ~ 1). The INSPIRE citation API gave unusable results through our tools.
- **2512.17995** (Hebecker–Lüst–Schachner–Schreyer):
  - Warping and mixed non-ISD/quantum corrections make KKLT F-term uplifts uncontrollable (abstract; eqs. 5.61, 5.65–5.67).
  - Eq. 5.67's condition W0𝒱^{2/3} ≫ 1 is "similar to the examples of" 2406.13751.
  - It does not evaluate the conifold, KPV or window questions for the five vacua.
- **2505.00149** (Schachner): proceedings review, with no new tests.
- **1902.07724** (BKS): KK modes and the conifold modulus have the same scaling (eq. 4.15); the KKLT EFT is "at the edge of reliability".
- **2206.04708** (Lüst–Randall), **2208.02826** (HSV), **2212.07437** (SV): general analyses, none at g_sM ~ 1 with the five parameter sets.

**We found no paper that computes tasks (2)–(4) for these five vacua.** This is a search result, not proof of absence.

## 7. Verdict per vacuum

**Table 4.** Verdict per vacuum and criterion.

| Criterion (inputs it depends on) | Ex1 | Ex2 | Ex3 | Ex4 | Ex5 |
|---|---|---|---|---|---|
| A. Conifold minimum with anti-D3, warped EFT (c′ = 1.18, c″ = 1.7546, off-shell e^{−4A} ∝ \|Z\|^{−4/3}) | FAIL ×2.75 | FAIL ×3.17 | FAIL ×3.23 | FAIL ×2.86 | FAIL ×3.88 |
| A′. Same with the Lüst–Randall potential (eq. 7.5; extrapolated) | undetermined, ~16% shift | undetermined | undetermined | undetermined | undetermined, ~21% shift |
| B. Conifold heavier than Kähler sector, warped EFT (c′, γ0) | no: 58% of Kähler modes heavier | no: 47% | no: 67% | no: 83% | no: 68% |
| C. Complex-structure sector inert under uplift (repo CSFterm; authors' EFT) | strongly violated (R_cs = 9.4; V_cs = 1.5 V_D̄3) | violated (2.6; 0.56) | mild (1.5; 0.23) | mild (1.3; 0.14) | violated (2.1; 0.47) |
| D. KPV with α′² (HSV/SV coefficients, NS5↔D5 map) | uncontrolled (×41) | uncontrolled (×55) | uncontrolled (×114) | uncontrolled (×167) | uncontrolled (×82) |
| E. Precursor Ξ in toy window (calibration heuristic) | outside | outside | inside | inside | inside |

The authors' own EFT (unwarped z_cf metric, leading-order KPV) passes everything for all five; we re-derived their bookkeeping V_dS = V_F,Kähler + V_D̄3 + V_cs from the repo. Flux quantization (odd fluxes) is flagged by the authors and was not tested here.

## 8. Failures and gaps

- We did not rebuild their multi-moduli potential (needs CYTools and GV data); we used repo outputs only.
- No Lüst–Randall computation with the anti-D3; our number there is a linear-response guess.
- The warped KK scale is ambiguous by a factor of about 23.
- The HSV shark-fin tip is not reproduced (16.2 vs 20.36).
- The SV coefficients rest on text extraction, validated against two quoted SV results only.
- The yield estimate assumes log-uniform Ξ.

## 9. Novelty

**Not new:**
- the BDGL-bound violation, which follows from public numbers and was noted in Phase 1;
- "g_sM > 20 needed" (HSV).

**New as far as our search shows (all OUR ATTEMPT, from public repository data):**
- (i) z_cf pinned to its F-flat relation, while the uplift relaxes through the dilaton with the factors and complex-structure F-term shares above.
- (ii) In the warped EFT the conifold mass lies inside the Kähler spectrum.
- (iii) Explicit α′²-corrected KPV numbers at the five parameter points.
- (iv) A calibrated KKLT toy fails for Ex1 and Ex2.
- (v) The window/step ratio of 0.12–0.25 and the implied 4–7% efficiency.
