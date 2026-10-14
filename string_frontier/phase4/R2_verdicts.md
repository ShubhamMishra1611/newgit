# Phase 4 — Reviewer R2 verdicts

Reviewer R2 (adversarial). I re-derived every claim with my own code before reading the phase-1/phase-3 scripts. Code and outputs are in `phase4/R2_code/`.

**How sources were checked.** Sources were fetched in this session with `mcp__Exa__web_fetch_exa`: arXiv abstract pages, HTML or PDF full texts, and explicit version URLs (`/abs/<id>v2`, `/html/<id>v1`). Equation numbers are given only where I read them.

**Exa staleness warning.** The unversioned abstract page for 2601.11705 returned by Exa is stale: it shows only v1. Fetching `/abs/2601.11705v2` shows that a v2 exists and that it changes the physics (see E1).

## Verdict table

| Claim | Verdict | One-line reason | Script |
|---|---|---|---|
| **E1** λ-family fails a₂₁/a₂₀ = 1/4, so it is "excluded by SUSY" except at λ=1 | **FAILED (overclaim). The numbers are VERIFIED; the result is not new.** | The ratios reproduce to 8 digits. However, the current version 2601.11705**v2** says the nonlinear relations need an extra *6-scalar parity condition*, and that without it "SUSY Ward identities can be solved without constraints" (v2 footnote 3). The ratio 1/4 also already follows from string monodromy (v2 eq. 4.14). The exclusion is an immediate corollary of 2610.02124 (eqs. 47, 76–80). | `E1_planar_lambda_family.py`, `E1_master_equation_test.py` |
| **E2** LP bound g̃₃ ≥ −10.612487 | **VERIFIED** | Exact per-spin minimax gives α* = 10.6124872188 and β* = 0.0671870. This equals the closed form 9/2 + (7/4)√(61/5) of 2011.02957 eq. 3.12. | `E2_g3_lowerbound.py` |
| **E3** α_min^IIB = 0.138949, α_min^IIA = 0.140260, E₃/₂(ρ) = 8.892745, E₅/₂(ρ) = 4.719492 | **VERIFIED** | Three independent evaluations agree: 6ζ(s)L(s,χ₋₃), the Fourier/Bessel expansion, and a lattice sum. The normalization E_s = Σ′τ₂^s/\|m+nτ\|^{2s} is the one that reproduces 2102.02847 eq. 5. A scan of the full fundamental domain puts the minimum at ρ. | `E3_eisenstein.py` |
| **E4** M_VS = 1/(stu) + 2ζ₃ + ζ₅Σs² + 2ζ₃² stu | **VERIFIED (known)** | Exact sympy series. The next term is (ζ₇/2)(Σs²)². This is textbook, and 2102.08951 eq. 3.42. | `E4_E5_VS_family.py` |
| **E5** λ-family (g̃₃, g̃₄) table | **VERIFIED** | My satellite-sum expansion reproduces all 8 entries, cross-checked against the hypergeometrics of 2408.03362 v1 eq. 10. The λ→0 limit (3, ½) and λ→∞ limit (2, ½) were checked analytically. | `E4_E5_VS_family.py` |
| **E6** 12 citations from 2025–26 | **CITATION PROBLEM (1 of 12)** | All 12 IDs exist, and the titles and authors match. 2601.11705 is cited in its superseded v1 form: the parity assumption is missing, and v1 eq. 1.6 is v2 eq. 1.11. The other 11 say what is claimed, with minor omissions noted below. | — |
| **C1** BP q ≈ 10^−1.43…10^−1.73 at J=100; lattice counts within 13% | **VERIFIED** (with a convention note) | I get −1.4288 / −1.7298. BP's own 10^−1.6 is reproduced exactly (−1.576) with \|2λ_bare\| = 1 and D = 1. The doc's −1.43 comes from choosing λ_bare = 1, not from degeneracy. My brute-force shell counts for J = 3, 4, 5 sit within 0.93–1.10 of eq. 2.24. | `C_checks_R2.py` |
| **C2** TCC √(2/3) at 0.7σ, √2 at 4.1σ from λ = 0.698 (+0.173/−0.202) | **PLAUSIBLE** | The arithmetic is right (0.68σ and 4.14σ). "4.1σ" is a Gaussian extrapolation of a 68% interval to the far tail, which the paper does not support. A likelihood that is Gaussian in λ² gives 5.6σ instead. Quote it as "≳4σ, dataset-dependent (DESY5)". | `C_checks_R2.py` |
| **C3** dark dimension l = 0.085–8.5 μm, M̂ = 5e8–2.4e9 GeV | **VERIFIED** (arithmetic) | Reproduced. The gap to the paper's 10⁹–10¹⁰ GeV is *both* the M_P convention (×2.9) and m = 1/l vs 2π/l (×1.85), not the M_P convention alone. Λ^−1/4 = 85.4 μm, against the paper's "≈88 μm". | `C_checks_R2.py` |
| **C4** 1903.01148 is not Palti's review; 1903.06239 is | **VERIFIED** | 1903.01148 is Xie–Luo, "Asymmetric Single Magnitude Four Error Correcting Codes". 1903.06239 is Palti, "The Swampland: Introduction and Review". | — |
| **D1** KKLT: σ_cr, V_AdS, dS minimum, windows | **FAILED (minor numerics)**; the rest is VERIFIED | σ_cr = 113.589, V_AdS = −1.9966e−15, the minimum at 114.894 with V = +1.7665e−17, the barrier at 139.28 and V″/V = −26.74 all verified. The exact window edges are D/σ³ ∈ [2.9732, 7.8612]e−9 (not [2.99, 7.71]) and D/σ² ∈ [2.5937, 9.7195]e−11 (not [2.64, 9.60]). The phase-1 values are logspace grid points (step 2.3%), not edges. The ratio is 2.64, not 2.58. | `D_checks_R2.py`, `D1_window_exact.py` |
| **D2** KPV 0.0805; I(0), c″, √g_s·M_min, g_sM² ≥ 46.25 | **VERIFIED** | (p/M)_crit = 0.080488, ψ_min(0.03) = 0.2198, I(0) = 0.71805, c″ = 1.75464, 6.8011, 46.255. c′ = 1.18 was taken as input and not re-derived. | `D_checks_R2.py` |
| **D3** five candidate vacua of 2406.13751 have g_sM² ∈ [11.9, 16.8]; content of 2512.17995 | **VERIFIED** (with a caveat) | Table 2 (M, g_sM) gives 16.82, 14.61, 14.33, 16.16 and 11.94. These are the five *representative* vacua (one per compactification) out of 30; the other 25 were not checked. The 2512.17995 abstract matches, including "effectively uncontrollable" and the 𝒱^{1/2} / 𝒱^{1/6} suppression. | `D_checks_R2.py` |
| **S2** Veneziano partial waves and per-level D bounds | **VERIFIED (known)** | Independent Gegenbauer projection, with a verified monomial→C_j^λ formula and t = −(s+4)(1−x)/2, s = n−1. It reproduces −(D−26)/(8(D−1)), −(D−26)/(2(D−3)(D+1)), and the per-level bounds 26.4, 27, 26.589019, 26.772677, 27.155578. Positivity holds at D = 26 for n ≤ 20 and fails at D = 26.01. The level-2 result is 2201.11575 eq. 1.3; sharpness at 26 is 2608.25183. | `S2_veneziano_pw.py` |
| **S3** one-loop β-function checks | **VERIFIED (textbook)** | All of the following reproduce: RNC −1/3; lattice tadpole slope −0.079590 vs −1/(4π), constant ln32/(4π); O(N) β = −(N−2)g²/(2π); cigar 1/K, giving K = α′/4 at D = 2; S³ h = 2R², H² = 24/R²; dilaton EL equation = −2 × [G^{mn}β^G − 4β^Φ]·√G e^{−2Φ}/α′. | `S3_beta_checks.py` |

---

## E1 in detail (the key claim)

### (a) Does 2601.11705 say this, and does it apply?

**The paper exists and v1 says it.** 2601.11705 is Elvang, Herderschee and Morales, "String Theory from Maximal Supersymmetry".
- **v1** (`/html/2601.11705v1`): the relations g²a₂,₀ = (2/5)a₀,₀², g²a₂,₁ = a₀,₀²/10 and g²a₃,₁ = 2g²a₃,₀ − a₁,₀a₀,₀ are **eq. (1.6)**. The stated inputs are only "N=4 SUSY in 4d, SU(4) R-symmetry, and tree-level factorization".

**v2 adds a parity assumption.** **v2** (`/abs/2601.11705v2`, `/html/2601.11705`) changes the abstract to "When imposing a parity condition on the 6-scalar amplitudes, we find…". The relations become **eq. (1.11)**. Footnote 3 states: *"Including the parity-odd terms in the ansatz, we find that SUSY Ward identities can be solved without constraints among the 4-point Wilson coefficients."* v2 also notes that the Coulomb branch of N=4 SYM violates the relations, so it must contain the parity-odd (Levi-Civita) contact terms. Parity-odd terms therefore occur in a physical N=4 theory.

**The phase-1 conclusion overstates this.** The phase-1 text says "SUSY requires exactly 1/4" and "no λ≠1 member can be the 4-point amplitude of a planar N=4 theory with tree-level 6-point factorization". That is **false as stated**: it is true only with the additional 6-scalar parity condition. 2610.02124 likewise needs scalar parity plus a discrete R-symmetry; its abstract says "a particular scalar-parity condition at six points".

**The definitions are compatible.** In 2408.03362 v1 (I confirmed the numbering in the v1 PDF):
- The planar family is eq. (18): A = −ρ/(st) + λΓ(λ−λs)Γ(λ−λt)/Γ(2λ−λs−λt)·₃F₂[…].
- Its satellite form is eq. (23).
- At λ = ρ = 1 it equals −Γ(−s)Γ(−t)/Γ(1−s−t) (eq. 20). CHR explicitly treat it as the stripped amplitude to be dressed with F⁴ for gluons (text after eq. 21).

The EHM definition is F(s,u) = −g²/(su) + Σa_{k,q}s^{k−q}u^q (eq. 1.3), with adjacent planar channels s and u, and the Veneziano case is F = −g²Γ(−s)Γ(−u)/Γ(1−s−u). So the map u_EHM ↔ t_CHR, g² ↔ ρ is exact. I checked numerically that the λ=1 satellite sum equals the Veneziano F to 15 digits. The ratio a₂₁/a₂₀ is independent of ρ, so the test is well posed.

### (b) My recomputation

I used two independent methods:
- **(M1)** the level sum Σ_n R̂(n,t)/μ(n)^{p+1}, with digamma t-derivatives;
- **(M2)** a bivariate polygamma Taylor expansion of the satellite sum (eq. 23).

M2 also matches the hypergeometric closed forms for a₀₀…a₃₀ exactly, and it is exactly crossing-symmetric: a₁₀ = a₁₁, a₂₀ = a₂₂, a₃₁ = a₃₂.

**a₂₁/a₂₀ by λ:**

| λ | a₂₁/a₂₀ |
|---|---|
| 0.25 | 0.0240979 |
| 0.5 | 0.0924938 |
| 0.9 | 0.2200934 |
| 1 | 0.25000000 |
| 1.1 | 0.2785314 |
| 2 | 0.4745940 |
| 4 | 0.6831190 |

These agree with phase 1. The ρ-free third-relation test, (2/5)a₀₀(2a₃₀ − a₃₁)/(a₁₀a₂₀), also fails for λ≠1: 1.366, 1.204, 1.032, 0.971, 0.798 and 0.639 at λ = 0.25, 0.5, 0.9, 1.1, 2 and 4.

**Stronger, all-orders test.** I plugged A_λ into the 2610.02124 master equation (45), (x₁x₃ − x₂x₄)B B + x₂x₄ B B = x₁x₃ B B with B = −suA/ρ.
- The residual is exactly quadratic in 1/ρ.
- For λ=1, its nontrivial root is ρ* = 1.000 at every kinematic point.
- For λ = 0.5, 0.9 and 2, the root is point-dependent: 1.458/1.278/1.236, 1.067/1.043/1.038 and 0.647/0.750/0.767 respectively.
- So no ρ solves the master equation for λ≠1.

### (c) Novelty

The result is **not new**:
1. **Monodromy already gives the ratio.** 2601.11705 v2 eq. (4.14) shows that *string monodromy alone* (with a₀,₀ = g²α′²π²/6) gives a₂,₀ = g²α′⁴π⁴/90 and a₂,₁ = g²α′⁴π⁴/360. So a₂₁/a₂₀ = 1/4 is the leading monodromy constraint, known since Huang et al. (2008.02293 era).
2. **2610.02124 already implies the exclusion.** Its eq. (47) solves the master equation to all orders. Eqs. (74)–(78) show that the forward measure must be Σδ(s′ − n/α′) under its analytic and positive dispersive assumptions. The CHR planar family has masses μ(n) = 1 + (n−1)/λ with positive forward residues, so for λ≠1 it is excluded immediately.
3. **But the specific case is not in the literature.** 2610.02124 cites 2408.03362 only in its §1 literature list. A grep found no explicit test of the CHR λ-family there or in 2601.11705 (whose §3.4/§5 exclude the IST, single massive exchange, the 1-loop amplitude and the sum of two Venezianos).

So the specific statement "the CHR planar λ-family fails the EHM relations" is not written anywhere I found, but it is a one-line corollary of results already published. Phase 1 itself says "probably a corollary of 2610.02124". That is honest, but the claim should be relabelled VERIFIED-KNOWN and restated **with the parity assumption**.

**Label problems in E_bootstrap.md:**
- §1 item 5 ("4d planar N=4 + SU(4)_R + tree factorization … give the nonlinear relations", labelled ESTABLISHED) must add "+ 6-scalar parity condition (2601.11705 v2; 2610.02124)".
- The "Bottom line" list of UV-shaped inputs should read "maximal SUSY + multi-point factorization + scalar parity".

## E6 citation spot-check (all abstract pages fetched in this session)

| ID | Title / authors as fetched | Status |
|---|---|---|
| 2601.11705 | Elvang, Herderschee, Morales, "String Theory from Maximal Supersymmetry" | **Exists, but cited in its superseded v1 form.** v2 requires the 6-scalar parity condition. Eq. 1.6 (v1) is eq. 1.11 (v2), and eqs. 1.7–1.8 (v1) are eqs. 1.12–1.13 (v2). |
| 2607.14230 | Berman, Caron-Huot, Chandra, Elvang, Herderschee, Lin, Morales, "Gravitational EFTs with Maximal SUSY and a Peculiar Parity" | OK: VS and IST corners; finite states near the first level isolate VS. |
| 2610.02124 | Huang, Wan, Wang, Zhou, "Infrared Consistency and the Uniqueness of String Amplitudes" | OK, but it also uses scalar parity plus a discrete R-symmetry. |
| 2605.11084 | Wan, Zhou, "Analytic Bootstrap of the Veneziano Amplitude" | OK: monodromy or splitting. |
| 2608.25183 | Chen, Yin, "Exact All-Level Positivity of the Bosonic Veneziano Amplitude" | OK: 3 < D ≤ 26. |
| 2603.04485 | Basile, Remmen, Staudt, "Higher-Spin and Higher-Point Constraints on Stringy Amplitudes" | OK: 7-point intercept; Gross satellite excluded. |
| 2508.09246 | Cheung, Remmen, Sciotti, Tarquini, "Strings from Almost Nothing" | OK. |
| 2502.20372 | Mansfield, "Positivity of the Veneziano Amplitude in Ten Dimensions" | OK (v2, Feb 2026). |
| 2606.19283 | Xu, "A Dispersive Bootstrap for the Virasoro–Shapiro Amplitude" | OK: island with Virasoro-inspired nonlinear ansatz. |
| 2604.15235 | Peng, Rodina, Tokareva, Xu, "Sampling the Graviton Pole and Deprojecting the Swampland" | OK: M/M_P ≲ 7.8 in D=5. |
| 2512.17828 | Eckner, Figueroa, Metayer, Tourkine, "Regge trajectories for UV completions of graviton scattering from polynomial boundedness" | OK. |
| 2610.01101 | Cao, "Analytic Dual Resonance and the String Bootstrap" | OK, with a minor omission: the result also needs a "joint crossing-spin bound" and "complex Regge decay", not just ADR plus an integer spectrum. |

## E2–E5 notes

**E2.** My method is independent of the authors' LP. For each even J, the condition is max_{x∈(0,1]}[(4J²/(d−2) − 3)x − βJ²(2J² − (5d−4))x²]; then I minimize the convex α(β). The result matches the analytic formula (3.12) of 2011.02957 in d = 4, 5 and 10. The paper prints β = 0.0671875; the correct value is 0.0671870, a trivial typo in the paper.

**E3.** The normalization is pinned by 2102.02847 eq. 5: α^IIA = ζ(3)/(32g_s^{3/2}) + π²g_s^{1/2}/96. This is exactly (2ζ₃τ₂^{3/2} + 4ζ₂τ₂^{−1/2})/64, i.e. the perturbative part of E₃/₂ with E_s = Σ′τ₂^s/\|m+nτ\|^{2s}.
- The IIA minimum is at g_s = 3√ζ₃/π = 1.04697 and equals π^{3/2}ζ₃^{1/4}/(24√3) = 0.140260117.
- E₅/₂(ρ) = 4.71949 in the same normalization, whose leading term is 2ζ₅τ₂^{5/2}.
- Matching to tree-level ζ₅ therefore needs c = ½ in H4.
- The minimum at ρ for all s > 0 is the Rankin/Cassels/Ennola/Diananda theorem on Epstein zeta functions. My scan of the fundamental domain agrees.
- Eq. 30 of 2102.02847 (supplement) was not read.

## C1–C4 notes

**C1.** BP hep-th/0004134 eqs. 2.23–2.26 were read. The text before eq. 2.25 says (√2κ₄)⁴λ_bare ~ 1, J = 100, D small ⇒ q ~ 10^−1.6.
- With J = 100, D = 1 and Δλ = 10⁻¹²⁰, eq. 2.24 gives log₁₀q = −1.576 for \|2λ_bare\| = 1, −1.429 for λ_bare = 1, and −1.724 for λ_bare = ¼.
- The phase-1 range "−1.43 … −1.73" is therefore mostly a units choice (λ_bare = 1) plus D = 2^J. It is not wrong, but BP's number is reproduced exactly only in BP's normalization.
- **Brute force:** I counted shell points with a meet-in-the-middle sum for incommensurate q, J = 3, 4, 5, with counts of 1600–9000. The ratios to the eq. 2.24 volume estimate are 0.93–1.10.

**C2.** The 2504.04226 abstract and text were read: λ = 0.698 (+0.173, −0.202) at 68.3%, "~3.5σ preference for nonzero λ" (99.95% CL from the posterior), V = V₀e^{−λφ}, data DESI DR2 + CMB distance priors + DESY5.
- **Arithmetic:** (√(2/3) − 0.698)/0.173 = 0.68 and (√2 − 0.698)/0.173 = 4.14.
- **The lower tail is Gaussian-like.** The Gaussian-in-λ extrapolation on the lower side gives 3.46σ for λ = 0, consistent with the paper's posterior-based 3.5σ. A likelihood Gaussian in λ² would give only 2.0σ there.
- **The upper tail is untested.** The upper tail beyond 1σ is not checked by anything in the paper. Quoting "4.1σ" for √2 is a Gaussian extrapolation; a λ²-Gaussian likelihood gives 5.6σ.
- The TCC number (≈0.7σ) is robust to this choice: 0.66–0.68σ.
- **Further caveats:** the result depends on the dataset (DESY5 SNe drive the preference). The TCC √(2/3) and the strong √2 are asymptotic or averaged bounds, while the fit probes Δφ ≪ 1, and the doc acknowledges this only for √2.
- **Verdict:** PLAUSIBLE. Replace "4.1σ" with "far outside the 68% interval (≈4–6σ under Gaussian-type tails; posterior tail not published)".

**C3.** 2205.12293 was read: eq. 2.1 gives M̂ = m^{n/(n+2)}M_pl^{2/(n+2)}; the text after eq. 3.2 gives Λ^{1/4} = 2.31 meV and "Λ^−1/4 ≈ 88 μm"; eq. 3.8 gives l ~ 0.1–10 μm; and λ ~ 10⁻¹–10⁻³.
- ħc/2.31 meV = 85.4 μm, so the paper's 88 μm is itself slightly off.
- With reduced M_P and m = 1/l: M̂ = 5.15e8–2.39e9 GeV. With non-reduced M_P: 1.5e9–7.0e9 GeV. With non-reduced M_P and m = 2π/l: 2.8e9–1.3e10 GeV. The paper's range needs both factors.

## D1–D3 notes

**D1, derivation.** From K = −3ln(2σ) and W = W₀ + Ae^{−aσ}, I derived V_F = a²A²e^{−2aσ}/(6σ) + aAe^{−aσ}W/(2σ²). This is identical to the KKLT Sec. III B formula, and V(σ_cr) equals KKLT eq. 14.

**D1, the eq. 16 ambiguity.** KKLT eq. 16 reads δV = **8D**/(Imρ)³ ("the factor of 8 is added for later convenience"). But Sec. III B then adds "**D/σ³**" and uses D = 3×10⁻⁹ with that normalization (Fig. 2). The phase-1 windows use the Sec. III B convention, so they are comparable to KKLT's example; the factor-8 inconsistency is internal to KKLT.

**D1, exact window edges.** I solved {V′ = 0, V = 0} for the lower edge and {V′ = 0, V″ = 0} for the upper edge:

| Uplift | Lower edge | Upper edge |
|---|---|---|
| D/σ³ | D = 2.97321e−9 (σ = 114.88) | D = 7.86124e−9 (σ = 121.58) |
| D/σ² | D = 2.59368e−11 | D = 9.71947e−11 |

- The phase-1 output itself labels its values "(grid)". They come from `np.logspace(-11,-7,400)`, so the edges are wrong by 0.6–2%.
- Phase 1 is also internally inconsistent: a lower edge at 2.99e−9 would give V_min ≈ 7e−18 at D = 3e−9, not the stated 1.77e−17.
- The qualitative point (window ratio ≈ 2.6, which needs near-continuous tuning) is unaffected: the exact ratio is 2.644.

**D2.** I(0) = ∫₀^∞(x coth x − 1)/sinh²x·(sinh 2x − 2x)^{1/3}dx = 0.718050. A side remark, which I did not check: the D doc writes "M_min … ≈ 6.8√N" for BDGL eq. 3.3. Verify whether that √N should be √p (number of anti-branes).

**D3.** 2406.13751 Table 2 was read. The entries are (M, g_sM) = (16, 1.051), (16, 0.913), (18, 0.796), (20, 0.808), (16, 0.746), so p/M = 1/16–1/20 and Ξ, KPV M > 12 (eq. 3.30), with η ≈ 2.6727 (eq. 3.32). The range 11.9–16.8 covers these 5 representatives only, not all 30 dS vacua; the full set is on GitHub and I did not check it.

## S2–S3 notes

**S2 kinematics.** For tachyon externals (m² = −1, α′ = 1), s + t + u = −4 and t = −(s+4)(1−x)/2. This matches 2201.11575 eq. 1.2 (their level index is shifted by one: their pole s = n is my α(s) = n+1).

**S2 Gegenbauer basis.** The standard normalization C_j^{(D−3)/2} is used, with C₁ = (D−3)x, so R₁ = 2x = [2/(D−3)]C₁. Positivity does not depend on any positive rescaling.

**S2 novelty.** All S2 statements are sanity reproductions; none claims novelty. The per-level bounds above 26 at n ≥ 4 are consistent with the sharp D ≤ 26 of 2608.25183, which is set by levels 2 and 3.

**S3 caveat.** The cigar check at D = 2 gives K = α′/4, i.e. string-scale curvature. The one-loop β function is not reliable there; the exact SL(2)/U(1) result is k = 9/4. This is fine as an algebraic check of the β-function formulas, but it should not be presented as a controlled background.
