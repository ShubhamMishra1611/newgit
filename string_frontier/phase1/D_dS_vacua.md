# Front D — de Sitter vacua in type IIB: KKLT, LVS, and the case against (Phase 1 map)

Labels: **ESTABLISHED** = derived in the cited source within its stated approximations; **CONJECTURE** = claimed, disputed, or depending on unverified assumptions; **OUR CHECK** = recomputed by us in `phase1/D_dS_checks/quick_checks.py` (output: `quick_checks_output.txt`). "Verified" means the arXiv abstract page was fetched in this session and the title and authors matched. Equation numbers are given only where we read them in fetched text.

## 1. State of the art

**KKLT building blocks**
- **ESTABLISHED (GKP flux background, hep-th/0105097; formulas as quoted in KKLT hep-th/0301240):**
  - The tadpole condition is χ(X)/24 = N_D3 + (1/2κ²T3)∫H3∧F3 (KKLT eq. 1).
  - The superpotential is W = ∫G3∧Ω (eq. 2).
  - The Kähler potential is K = −3ln[−i(ρ−ρ̄)] − ln[−i(τ−τ̄)] − ln[−i∫Ω∧Ω̄] (eq. 3).
  - V = e^K(Σ g^{ab̄}D_aW D_b̄W̄ − 3|W|²) reduces to the positive semi-definite no-scale form e^K Σ_{i≠ρ}|D_iW|² (eq. 4). ρ is therefore unfixed at tree level.
  - The throat warping is e^{A_min} ~ exp[−2πK/(3g_sM)] (eq. 7).
- **ESTABLISHED as 4d supergravity (KKLT hep-th/0301240):**
  - K = −3ln[−i(ρ−ρ̄)] (eq. 11) and W = W0 + A e^{iaρ} (eq. 12).
  - With ρ = iσ, the SUSY AdS condition is W0 = −A e^{−aσ}(1 + 2aσ/3) (eq. 13), and V_AdS = −a²A²e^{−2aσ}/(6σ) (eq. 14).
  - The anti-D3 energy is δV = 2a0⁴T3/(g_s⁴(Im ρ)³) (eq. 15), modelled as D/(Im ρ)³ (eq. 16). The text says a factor 8 was added "for convenience", but our PDF extraction of that factor is ambiguous.
  - Full uplifted potential (unnumbered, Sec. III B): V = (aAe^{−aσ}/2σ²)(σaAe^{−aσ}/3 + W0 + Ae^{−aσ}) + D/σ³.
  - Worked example: W0 = −10⁻⁴, A = 1, a = 0.1, giving σ_cr ≈ 113; D = 3×10⁻⁹ gives a dS minimum. The canonical field is φ = √(3/2) ln σ (Sec. IV).
  - **OUR CHECK:**
    - σ_cr = 113.589 and V_AdS = −1.997×10⁻¹⁵.
    - With D = 3×10⁻⁹ the dS minimum sits at σ = 114.89 with V = +1.77×10⁻¹⁷. The barrier maximum is at σ = 139.3, with V''/V = −26.7 in canonical Planck units.
    - A dS minimum exists for D ∈ [2.99, 7.71]×10⁻⁹ with uplift D/σ³, and for D ∈ [2.64, 9.60]×10⁻¹¹ with D/σ². The D/σ² case is the warped scaling V ∝ 𝒱_E^{−4/3}.
- **ESTABLISHED (KPV hep-th/0112197):**
  - p anti-D3s at the KS tip puff up into an NS5 on S² ⊂ S³.
  - V_eff(ψ) = A0[V2(ψ) − (2ψ − sin2ψ)/2π] (eq. 45), with V2 = (1/π)√(b0⁴sin⁴ψ + (πp/M − ψ + ½sin2ψ)²) (eq. 40), A0 = μ3M/g_s (eq. 43) and b0² ≈ 0.93266 (eq. 9).
  - The configuration is metastable for p/M ≲ 8% (text after eq. 45).
  - The decay rate is ≈ exp(−27b0¹²g_sM⁶/(512πp³)) (eq. 65).
  - **OUR CHECK:** critical p/M = 0.0805. At p/M = 0.03, ψ_min = 0.2198, against the eq.-46 estimate of 0.2167.
- **ESTABLISHED (supergravity statement):** the anti-D3 sits in a nilpotent goldstino multiplet. Kallosh–Wrase 1411.1121 (eq. 1.3) give W = W0 + Ae^{−aρ} − M²S and K = −3ln(ρ+ρ̄) + SS̄, with S² = 0. At S = 0 this gives V ⊃ M⁴/(ρ+ρ̄)³, the KKLT uplift (eq. 1.2). The worldvolume derivation is in 1411.1121 and 1502.07627.
- **CONJECTURE (contested):** whether the 4d picture survives in 10d.
  - 1707.08678 (Moritz–Retolaza–Westphal) says the single-modulus gaugino-condensate KKLT cannot be uplifted.
  - 1908.04788 (Kachru–Kim–McAllister–Zimet) finds exact 10d/4d agreement.
  - 1902.01410 (Hamada–Hebecker–Shiu–Soler) supports the uplift via a perfect-square quartic gaugino term.
  - 1902.01415 (Gautason–Van Hemelryck–Van Riet–Venken) finds that the "local" 10d procedure does not reproduce the EFT.
  - **Current balance:** largely answered, but how to treat the gaugino bilinear remains open.

**LVS**
- **ESTABLISHED (BBCQ hep-th/0502058; CQS hep-th/0505076):**
  - With K = K_cs − 2ln(𝒱 + ξ/(2g_s^{3/2})) and ξ = −χζ(3)/(2(2π)³) (0505076 eq. 24), and W = W0 + Σ A_i e^{−a_iT_i} (eq. 14), the large-volume potential is V = λ√τ_s e^{−2a_sτ_s}/𝒱 − μτ_s e^{−a_sτ_s}/𝒱² + ν/𝒱³ (0502058 eq. 48; 0505076 eqs. 47, 138), with λ ~ a²|A|², μ ~ a|AW0| and ν ~ ξ|W0|².
  - It has a non-SUSY AdS minimum at τ_s = (4νλ/μ²)^{2/3} with 𝒱 ∝ W0 e^{aτ_s} (0502058 eqs. 52–53). For P⁴[1,1,1,6,9], τ4 ≈ 4.11 (eq. 54).
  - **Uplift example:** an ε/𝒱^{4/3} uplift (0505076 eq. 137) with a = 2π/10, W0 = 10, g_s = 0.1 and ξ = 1.31 gives τ4 = 4.17 and 𝒱0 = 4.5×10¹⁰. Adding ε = 4.9×10⁻¹⁹ yields dS at τ4 = 4.213, 𝒱 = 5.9×10¹⁰.
  - **Neglecting α′ corrections:** this is consistent only for tiny W0. With instantons, A = 1 and g_s = 0.1, it requires W0 ≲ 10⁻⁷⁵ (0505076 around eq. 27).
- **ESTABLISHED (EFT with uplift; Junghans 2201.03572):**
  - K is given in eq. 2.1, with conifold modulus Z and nilpotent Y. W = W0 + A_se^{−a_sT_s} − (M/2πi)Z(ln(Z/Λ0³) − 1) + (iK/g_s)Z + (i√c″/(√π g_sM))Z^{2/3}Y (eq. 2.2).
  - The constants are c′ ≈ 1.18, c″ ≈ 1.75 and ξ′ = 9g_sM² (eq. 2.4). The potential is eq. 2.5.
  - **Vacuum solution:** τ_s = ξ^{2/3}/((2κ_s)^{2/3}g_s) + 1/(3a_s) + 4α/(15a_s) + O(g_s) (eq. 2.7), so τ_s ∝ 1/g_s. The vacuum energy is V0 ∝ (α − 1) (eq. 2.12).
  - A dS minimum requires α ∈ (1, 9/4) (eq. 2.18).
  - The conifold bound is g_sM² ≥ 64πc′c″/9 ≈ 46.1 (eq. 2.9).
- **CONJECTURE (critique; 2201.03572 and 2205.02856):**
  - Curvature, warping and g_s corrections are parametrically larger on-shell than off-shell, because τ_s ∝ 1/g_s and a cancellation that Junghans calls the "non-perturbative no-scale structure". Some corrections are unsuppressed.
  - In the explicit model he studies, at least one correction is O(1) wherever V0 > 0, assuming O(1) coefficients.
  - The minimal D3 tadpole needed for control ranges from O(500) to O(10⁶).
  - Other uplift proposals (D-terms, T-branes) were not examined here.

**Explicit constructions**
- **ESTABLISHED (small W0 and AdS vacua):**
  - "Perturbatively flat vacua" have W_eff(τ) = Σ_q A_q(M·q)/(2πi)² e^{2πiτ p·q} (1912.10047 eq. 10). The example with M = (−16, 50), K = (3, −4) and Q_flux = 124 gives τ = 6.856i and |W0| = 2.037×10⁻⁸ (eqs. 18–20).
  - 2107.09064 builds SUSY AdS vacua with |V| < 10⁻¹²³ and controlled α′ effects, using Gopakumar–Vafa invariants. 2107.09065 is "A Cosmological Constant That is Too Small".
- **ESTABLISHED at leading order, CONJECTURE as string vacua (McAllister–Moritz–Nally–Schachner 2406.13751):**
  - They study 33,371 anti-D3 PFVs. Of these, 396 are "aligned", meaning Ξ ≈ (|z_cf|^{4/3}/|W0|²)(𝒱_E^{2/3}Ṽ^{1/3}/(g_sM)²)ζ ~ 1 with ζ ≈ 114 (v2 eq. 3.34). They find 5 compactifications carrying 30 dS vacua.
  - The uplift is V_D̄3 = c/𝒱_E^{4/3} with c = ηz_cf^{4/3}/(g_sM²Ṽ^{2/3}) and η ≈ 2.6727 (2505.00149 eq. 33). KPV metastability needs M/p ≳ 12 (2406.13751 v2 eq. 3.30).
  - The parameters (2505.00149 Table 1) are g_sM = 0.75–1.05, M = 16–20, p = 1 and g_s = 0.04–0.066.
  - **Caveats the authors flag themselves:**
    - α′ corrections to the anti-D3 sector are "a particularly acute risk"; they expect future searches to reach g_sM ≳ 5 but not ≳ 15.
    - Odd flux quanta may not be allowed. Since Q_O ≡ 0 mod 4, even fluxes force p ∈ 2ℤ.
    - The conifold modulus is integrated out supersymmetrically, treated as heavy.
- **ESTABLISHED (explicit LVS dS models in leading-order EFT):**
  - Crinò–Quevedo–Valandro 2010.15903 and AbdusSalam–Crinò–Shukla 2206.12889 (h¹¹ = 4, K3-fibred). Finding an orientifold with both a large D3 tadpole and suitable O3-planes is the bottleneck.
  - 2507.00615 numerically scans about 80,000 Calabi–Yaus with h¹¹ ≤ 6. It finds KKLT-like vacua for 10⁻¹⁵ ≤ |W0| ≤ 10⁻², LVS vacua for 0.1 ≤ |W0| ≤ 100, and coexisting minima.

**The case against**
- **CONJECTURE (conifold destabilization; Bena–Dudas–Graña–Lüst 1809.06861):**
  - Anti-D3s destabilize the light conifold modulus S unless √g_s M > M_min = (8/3)√(πc′c″) ≈ 6.8√N (eq. 3.3), with c′ ≈ 1.18 (eq. 2.18) and c″ = 2^{1/3}/I(0) ≈ 1.75 (eq. 2.33).
  - With the tadpole |Q| = 34 this gives h < 4.6, so Λ_IR/Λ_UV > 0.2 (eqs. 3.5–3.6).
  - **OUR CHECK:** I(0) = 0.71805, c″ = 1.7546, M_min = 6.801 and g_sM² ≥ 46.25.
  - Lüst–Randall (2206.04708) argue the bound is "too strong", because neglected constraint equations from the 10d metric remove the instability.
- **CONJECTURE (singular bulk; Gao–Hebecker–Junghans 2009.03914):**
  - KKLT implies |∂h|/h ≳ g_sM² near the instanton divisor (eqs. 1.1, 3.7), so the warp factor turns negative over O(1) regions.
  - An escape requires τ_Σ/𝒱^{2/3} ≲ 1/(g_sM²) (eq. 4.2), which in turn needs h¹¹ ≳ (g_sM²)^{3/2} ≈ 42 (text).
- **CONJECTURE (throat size and tadpole):**
  - Carta–Moritz–Westphal 1902.01412: it is hard to fit a throat long enough to avoid decompactification inside a generic Calabi–Yau.
  - Bena–Blåbäck–Graña–Lüst 2010.10519: the tadpole conjecture for many moduli.
- **CONJECTURE (α′ corrections to KPV):**
  - Hebecker–Schreyer–Venken 2208.02826: leading-order KPV is recovered only for g_sM > 20, which in LVS implies a tadpole of order 10³–10⁴. Their fit (eq. 2.13) covers only g_sM ∈ [20.4, 100].
  - Schreyer 2402.13311 states the current bound as g_sM² > 144, or O(100) from the anti-D3-stack side.
  - Kim 2409.19048 starts a string-field-theory treatment of the KS throat.
- **CONJECTURE (2025):** Hebecker–Lüst–Schachner–Schreyer 2512.17995 find that warping and mixed non-ISD/quantum corrections dominate in KKLT, "rendering the scenario effectively uncontrollable with current methods". In LVS-like models the same corrections are suppressed by 𝒱^{1/2} or 𝒱^{1/6}.
- **CONJECTURE (dynamical objections):** Sethi 1709.03554 argues flux SUSY-breaking backgrounds generically roll. Danielsson–Van Riet 1804.01120 review the case that string theory has no dS vacua.
- **CONJECTURE (swampland):**
  - Obied–Ooguri–Spodyneiko–Vafa 1806.08362 propose |∇V| ≥ cV.
  - The refined version (1810.05506) is "|∇V| ≥ cV, or min∇∇V ≤ −c′V" (equation numbers not checked).
  - Any metastable dS minimum violates both forms. The candidate vacua and the dS conjecture therefore cannot both be right.
  - **OUR CHECK:** at the KKLT-toy barrier top, V''/V = −26.7, consistent with the refined form.

**OUR CHECK (control numbers of the five candidate vacua)**
- g_sM² = 16.8, 14.6, 14.3, 16.2 and 11.9. This is 0.26–0.36 of the 46.1 bound and about 0.1 of the 144 bound.
- p/M = 0.05–0.0625. That is below the leading-order KPV limit of 0.0805, but in the regime g_sM ~ 1, where no α′-corrected analysis exists.

**Where the field stands (2024–26):**
- Leading-order EFT existence is ESTABLISHED for KKLT (2406.13751) and LVS (2010.15903, 2206.12889).
- Survival beyond leading order is the open fight, on four fronts: the anti-D3 sector at g_sM ~ 1; the light conifold modulus; flux quantization; and warping, argued fatal to KKLT control (2512.17995) and milder for LVS, where Junghans' on-shell enhancement remains.

## 2. The three sharpest open questions

1. **Anti-D3 and conifold sector at g_sM ≈ 0.75–1.05, M = 16–20, p = 1 (the actual parameters of the 2406.13751 vacua).**
   - **Question:** does the fully corrected NS5/anti-D3 potential V_eff(ψ) have a local minimum with positive energy, and does the conifold modulus keep a minimum?
   - **Answer form:** two numbers. (a) The true critical (p/M)_crit(g_sM) at g_sM ≈ 1. (b) The true threshold (g_sM²)_crit for the conifold, which is 46.1 if BDGL/Junghans hold and lower if Lüst–Randall hold.
   - **How it is answered:** the question is closed once both thresholds are computed with α′ and backreaction corrections under control, for example by string field theory (2409.19048) or by a supersymmetric localization of the KS tip.
2. **Flux quantization in Calabi–Yau orientifolds.**
   - **Question:** are odd integer flux quanta allowed on the 3-cycles used in 2406.13751, given that odd fluxes on toroidal orientifolds need exotic O3-planes?
   - **If not:** does an even-flux candidate with p = 2 anti-D3s exist, satisfying M ≥ 24, Ξ ~ 1 and W0 ≪ 1?
   - **Answer form:** yes or no, decided by topology (half-cycles and the torsion/Freed–Witten structure of the orientifold) plus an explicit search.
3. **Existence of a controlled corner for LVS dS.**
   - **Question:** do explicit topological data (χ, κ_s, a_s, an orientifold with D3 tadpole Q) exist such that all of the following hold at once: α ∈ (1, 9/4), g_sM² ≥ 46.1 (or ≥ 144), g_sM ≳ 20 with MK ≤ Q, the warping corrections of 2512.17995 below 10%, and Junghans' on-shell-enhanced corrections below 10% with O(1) coefficients?
   - **Answer form:** a single geometry meeting all the inequalities, or a proof that Q_min(χ, κ_s, a_s) exceeds every known tadpole bound.

## 3. Key papers

| arXiv ID | Authors | Year | Title (short) | Verified? | Why it matters |
|---|---|---|---|---|---|
| hep-th/0105097 | Giddings, Kachru, Polchinski | 2001 | Hierarchies from fluxes | yes | ISD flux background, warped throats |
| hep-th/0112197 | Kachru, Pearson, Verlinde | 2001 | Brane/flux annihilation | yes | Anti-D3 metastability, p/M ≲ 0.08 |
| hep-th/0301240 | Kachru, Kallosh, Linde, Trivedi | 2003 | de Sitter vacua in string theory | yes | KKLT |
| hep-th/0502058 | Balasubramanian, Berglund, Conlon, Quevedo | 2005 | Systematics of moduli stabilisation | yes | LVS potential |
| hep-th/0505076 | Conlon, Quevedo, Suruliz | 2005 | Large-volume flux compactifications | yes | LVS spectrum, uplift example, α′ necessity |
| 1411.1121 | Kallosh, Wrase | 2014 | Emergence of spontaneously broken SUSY on anti-D3 | yes | Nilpotent goldstino uplift |
| 1707.08678 | Moritz, Retolaza, Westphal | 2017 | Towards de Sitter from 10D | yes | 10d obstruction to single-modulus KKLT |
| 1709.03554 | Sethi | 2017 | Supersymmetry breaking by fluxes | yes | Rolling, not static |
| 1804.01120 | Danielsson, Van Riet | 2018 | What if string theory has no dS vacua? | yes | Review of the skeptical case |
| 1806.08362 | Obied, Ooguri, Spodyneiko, Vafa | 2018 | De Sitter space and the swampland | yes | dS conjecture |
| 1810.05506 | Ooguri, Palti, Shiu, Vafa | 2018 | Distance and dS conjectures | yes | Refined conjecture |
| 1809.06861 | Bena, Dudas, Graña, Lüst | 2018 | Uplifting runaways | yes | Conifold destabilization, √g_s M > 6.8 |
| 1902.01410 | Hamada, Hebecker, Shiu, Soler | 2019 | Understanding KKLT from 10d | yes | Defends the uplift in 10d |
| 1902.01412 | Carta, Moritz, Westphal | 2019 | Gaugino condensation and small uplifts | yes | Throat-fitting problem |
| 1908.04788 | Kachru, Kim, McAllister, Zimet | 2019 | dS vacua from ten dimensions | yes | 10d = 4d agreement |
| 1912.10047 | Demirtas, Kim, McAllister, Moritz | 2019 | Vacua with small flux superpotential | yes | PFVs, W0 ≈ 2×10⁻⁸ |
| 2009.03914 | Gao, Hebecker, Junghans | 2020 | Control issues of KKLT | yes | Singular-bulk problem |
| 2010.15903 | Crinò, Quevedo, Valandro | 2020 | dS from anti-D3 in LVS | yes | Explicit LVS dS |
| 2107.09064 | Demirtas, Kim, McAllister, Moritz, Rios-Tascon | 2021 | Small cosmological constants | yes | Explicit KKLT AdS, |V| < 10⁻¹²³ |
| 2201.03572 | Junghans | 2022 | LVS dS vacua are probably in the swampland | yes | On-shell enhancement of corrections |
| 2205.02856 | Junghans | 2022 | Topological constraints in LVS | yes | Tadpole O(500)–O(10⁶) |
| 2206.04708 | Lüst, Randall | 2022 | Effective theory of warped compactifications | yes | Disputes the conifold bound |
| 2208.02826 | Hebecker, Schreyer, Venken | 2022 | Curvature corrections to KPV | yes | g_sM > 20 |
| 2406.13751 | McAllister, Moritz, Nally, Schachner | 2024 | Candidate de Sitter vacua | yes | First explicit leading-order KKLT dS |
| 2512.17995 | Hebecker, Lüst, Schachner, Schreyer | 2025 | Effective potentials, warping, F-term uplifting | yes | Warping: KKLT uncontrolled, LVS suppressed |

Also verified: 2010.10519, 2204.06009, 2206.12889, 2212.07437, 2402.13311, 2409.19048, 2501.03984, 2505.00149, 2505.15883, 2507.00615, 2512.17095 (TASI lectures, McAllister–Schachner), 1902.01415, 1912.09948, 1502.07627, 2107.09065.

## 4. Computable handles (Python, a few hours each)

1. **KKLT toy, dS window and flux discreteness.**
   - **Compute:** the D-window for an uplift D/σⁿ (n = 2, 3) as a function of aσ_cr, the CDL or Hawking–Moss exponent, and the swampland ratios at every critical point.
   - **Compare to:** KKLT eqs. 13–14 and the Sec. III B example; OUR CHECK numbers above.
   - **Observation:** at g_sM = 1, shifting K → K+1 rescales the warped uplift by e^{−8π/3} ≈ 2.3×10⁻⁴, while the dS window spans only a factor of about 2.6. Landing in the window therefore requires near-continuous tuning of W0 and g_s.
   - **What would count as new:** a predicted dS yield per aligned vacuum, tested against 2406.13751's 5 compactifications out of 396 aligned (data on GitHub: AndreasSchachner/kklt_de_sitter_vacua).
2. **KPV with α′² terms at the candidate parameters.**
   - **Compute:** implement eqs. 40/45 (reproduced: 0.0805), add the α′² NS5 and anti-D3 terms of 2208.02826 (eqs. 2.12–2.15) and 2212.07437, and evaluate at g_sM ∈ [0.7, 1.1], p/M = 1/16–1/20.
   - **Compare to:** the HSV eq. 2.13 fit, which is valid only for g_sM ≥ 20.4.
   - **What would count as new:** the size of the known corrections at ψ_min for Examples 1–5, which we have not seen reported, and whether the truncated potential keeps a minimum there.
3. **Conifold modulus versus the "heavy conifold" assumption.**
   - **Compute:** implement BDGL eqs. 2.17 and 2.32 (reproduced: 6.80 and 46.25) and Junghans eq. 2.5, then evaluate m_S² (BDGL eq. 2.24) against the Kähler-modulus mass for the five candidates (g_sM² = 11.9–16.8).
   - **Compare to:** eq. 3.3 and the Lüst–Randall corrected potential.
   - **What would count as new:** a direct test of whether integrating out z_cf supersymmetrically, as 2406.13751 does, is self-consistent in their own vacua.
4. **LVS dS window and minimal tadpole.**
   - **Compute:** implement 2201.03572 eq. 2.5, reproduce eqs. 2.6–2.9 and 2.12 and the window α ∈ (1, 9/4), and recover the 0505076 example (τ4 = 4.17, 𝒱 = 4.5×10¹⁰, then 4.213 and 5.9×10¹⁰). Fixing conventions is part of the task.
   - **Scan:** (g_s, M, K, W0, χ, κ_s, a_s) under the constraints of Q3. Output Q_min.
   - **Compare to:** 2205.02856 (O(500)–O(10⁶)) and the h¹¹ = 4 model of 2206.12889.
   - **What would count as new:** a closed-form Q_min surface.
5. **PFV small W0 and even fluxes.**
   - **Compute:** reproduce |W0| = 2.037×10⁻⁸ and τ = 6.856i from 1912.10047 eqs. 4–7 and 11–18. Then enumerate all (M, K) satisfying eq. 15 with Q ≤ 138, splitting the results into even and odd fluxes.
   - **Compare to:** the Denef–Douglas estimate W0,min ~ 6×10⁻⁷ (text near eq. 21).
   - **What would count as new:** the smallest W0 reachable with even fluxes only. This is a direct, cheap probe of Q2.
