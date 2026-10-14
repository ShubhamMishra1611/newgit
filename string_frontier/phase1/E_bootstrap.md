# FRONT E — Is the string S-matrix fixed by consistency? (Phase 1 map)

Conventions: "verified" means the arXiv abstract page was fetched this session and title/authors matched. Equation numbers are given only where I read them. **OUR CHECK** means I recomputed it; code and output are in `phase1/E_checks/`.
Corrections to the brief's guesses: 2308.03833 is *Bespoke Dual Resonance*. The Cheung–Hillman–Remmen paper on Virasoro–Shapiro (VS) uniqueness is **2408.03362** (2406.02665 is the open-string paper). The closest "Bootstrapping string theory" is 2310.10710 (Chiang–Huang–Weng). GPV 2102.02847 reports a *minimum of α* near the string value, not a kink.

## 1. State of the art

**Positivity and dispersion: the EFT side**
- **ESTABLISHED.** Analyticity plus causality forces the leading irrelevant couplings, such as (∂φ)⁴, to be positive (hep-th/0602178; eq. not checked). This was later sharpened by the EFT-hedron, which phrases the constraints as total positivity / Hankel conditions (2012.15849), by moment bounds (2011.00037), and by full-crossing nonlinear bounds (2011.02400) (eqs. not checked).
- **ESTABLISHED (2011.02957).**
  - Low-energy amplitude: M_low = −g²(1/s+1/t+1/u) − λ + g₂(s²+t²+u²) + g₃ stu + g₄(s²+t²+u²)² + … (eq. 2.3).
  - Twice-subtracted B_k sum rules (2.18–2.19) give heavy averages: g₂=⟨1/m⁴⟩, g₃=⟨(3−4𝒥²/(d−2))/m⁶⟩, g₄=⟨1/(2m⁸)⟩ (2.22), with 𝒥²=J(J+d−3).
  - First null constraint: ⟨𝒥²(2𝒥²−(5d−4))/m⁸⟩=0 (2.23).
  - Resulting bounds: −2(5d−4)/(d−2) ≤ g̃₃ ≡ g₃M²/g₂ ≤ 3 (3.6), and the optimal 3-functional bound in d=4 is g̃₃ ≥ −10.61249 (3.10–3.11).
  - **OUR CHECK:** a linear program over (m ≥ M, even J ≤ 398) gives α*=10.612487 and β*=0.067187, matching (3.10).
- **ESTABLISHED (2102.08951), gravity.** The graviton pole is handled by measuring couplings at impact parameter b ~ 1/M instead of in the forward limit (§3).
  - 10D maximal SUSY: M_susy = 8πG/(stu) + g₀ + g₂(s²+t²+u²) + … (3.34), with s²M_susy → 0 (3.35).
  - Sum rules: C₋₂: 8πG/(−u) = ⟨m²(2m²+u)P_J(1+2u/m²)⟩ (3.36), and the improved g₀ measurement (3.40).
  - Result: **0 ≤ g₀M⁶/(8πG) ≤ 3.000** (3.41). Type II gives 2ζ(3) ≈ 2.40 (3.42), so string theory is interior. "It would be interesting to find a model which saturates the bound."
  - Twice-subtracted gravitational dispersion relations are justified in 2202.08280. Four-dimensional two-sided graviton bounds, up to an IR logarithm, are in 2201.06602.
- **ESTABLISHED, numerical, non-perturbative 10D (2102.02847).**
  - Definition: T/(8πG) = s⁴(1/(stu) + α ℓ_P⁶ + …) (eq. 3), with α ≥ 0 from the sum rule (7)–(8).
  - String values: α^IIB = E_{3/2}(τ,τ̄)/2⁶ ≥ 0.1389 (4, 27, 30) and α^IIA ≥ 0.1403 (5).
  - Primal bootstrap: **α_min^Boot ≈ 0.13 ± 0.02** (10).
  - **OUR CHECK:** 0.138949 and 0.140260.
  - 2212.00151 (9, 10, 11D): the bootstrap minimum is "very close but not equal" to the string/M-theory minimum (abstract; eq. not checked).
  - **CONJECTURE (GPV):** string theory fills the entire allowed α range.

**The string amplitudes themselves**
- **ESTABLISHED.**
  - Veneziano: A_V = Γ(−s)Γ(−t)/Γ(−s−t) (2406.02665 eq. 4). Level truncation: A_V(s,−k) = Σ_{n<k}(1−k)_n/n!·1/(n−s) (eq. 5); **OUR CHECK** numerically for k=1,2,3,5.
  - SUSY-stripped form: F(s,u) = −g²Γ(−s)Γ(−u)/Γ(1−s−u), with a₀₀=g²ζ₂, a₁₀=g²ζ₃, a₂₀=g²ζ₄, a₂₁=g²ζ₄/4, a₃₀=g²ζ₅, a₃₁=g²(2ζ₅−ζ₂ζ₃) (2601.11705 eqs. 1.7–1.8).
  - Virasoro–Shapiro: M_VS = −Γ(−s)Γ(−t)Γ(−u)/[Γ(1+s)Γ(1+t)Γ(1+u)] (2408.03362 eq. 1). Residues: R(n,t) = Π_{k=1}^{n−1}(t+k)²/(n!)² (eq. 3), i.e. *double* zeros.
  - **OUR CHECK:** M_VS = 1/(stu) + 2ζ₃ + ζ₅(s²+t²+u²) + 2ζ₃² stu + …, consistent with 2102.08951 eq. 3.42.
  - Direct proofs of partial-wave positivity: superstrings in D ≤ 6 (2201.11575), Veneziano in D ≤ 10 (2502.20372), bosonic Veneziano for 3 < D ≤ 26 (2608.25183).
- **ESTABLISHED, with technical assumptions.** For s,t ≫ 1, weakly coupled massive higher-spin theories coincide asymptotically with Veneziano and have infinitely many trajectories (1607.04253). Caveat: Häring–Zhiboedov build a unitary 3-parameter exponentially soft deformation of Veneziano (2311.13631), which violates those assumptions (as noted in 2408.03362).
  - No meromorphic amplitude can have finitely many Regge trajectories (2405.21057), and the same holds with external gravitons (2512.17828).
- **ESTABLISHED: 2→2 consistency alone does NOT fix the string.** Unitary deformations exist: infinite-parameter dual-resonant families (2210.12163, 2302.12263), user-chosen spectra (2308.03833), and Regge-pit numerics (2401.08736).

**"String = unique solution of X+Y+Z" statements (precise inputs)**
1. **Planar 4-point (2406.02665), ESTABLISHED.**
   - Inputs: a tree-level dual-resonant sum A = Σ_n R(n,t)/(μ(n)−s) (eq. 3); polynomial residues; level truncation R(n,t) = c(n)Π_{k≤n}(t−ξ(k)) (eq. 7); planar crossing.
   - These leave a 3-parameter (q,r,w) family with μ(n) = [n]_q (eqs. 6, 15, 17–18).
   - Positivity forces q=1 (App. A, eqs. 25–26). Superpolynomial softness (Regge eq. 19) then forces r=w=0, which is Veneziano.
2. **Closed 4-point (2408.03362), ESTABLISHED within the stated scope.**
   - Inputs: level truncation of both M and ∂_tM, which gives R = c(n)Π(t−ξ(k))² (eq. 4).
   - Solution: μ(n) = (n+λ−1)/λ, ξ(n) = −(2n+λ−1)/(2λ), c(n) = λ^{2n−2}/((1+3λ)/2)²_{n−1} (eq. 5), with amplitude (7) and coefficients (10).
   - Regge limit: M ~ s^{2J} + s⁻²[−ρ²/t + (3λ−1)²/(4λ(1−λ+λt))] (11). Superpolynomial softness forces λ=ρ=1, which is VS.
   - Caveat: uniqueness is shown only for solutions continuously connected to VS, or if odd partial waves vanish (text after eq. 5).
3. **2508.09246 (PRL), claimed proof.**
   - Residues *must* vanish on "Regge zeros" α(t)+r=0 (eq. 6).
   - Adding (i) ultrasoftness and (ii) "minimal zeros" (every residue zero is a Regge zero) gives μ(n) ~ J(n) and α(t) ~ t (eq. 13), and then Veneziano or VS.
   - The authors call (ii) their least conservative assumption.
4. **Dispersive sum rules + unitarity + (string monodromy OR splitting/hidden-zero condition) gives Veneziano**, an analytic proof (2605.11084). Earlier monodromy numerics placed Type I within 10⁻⁴ and recovered D_crit = 10 (2310.10710; see also 2310.10729).
5. **4d planar N=4 + SU(4)_R + tree factorization at 4, 5, 6 points.**
   - These give the nonlinear relations **g²a₂₀ = (2/5)a₀₀², g²a₂₁ = a₀₀²/10, g²a₃₁ = 2g²a₃₀ − a₁₀a₀₀** (2601.11705 eq. 1.6).
   - Adding positivity up to k=7, the bounds sit within ~5·10⁻⁵ of Veneziano. Exact convergence is a **CONJECTURE**.
   - 2610.02124 proves the relations to all orders. Under analytic and positive dispersive assumptions they uniquely select Veneziano.
6. **Gravity: here uniqueness FAILS without an extra input.**
   - In 4d N=8 with SU(4)×SU(4) R-symmetry, "peculiar parity", 4/5/6-point factorization and positivity, the allowed region is non-convex with two corners: VS and an infinite-spin-tower (IST) amplitude with every spin at one mass (2607.14230).
   - "Finite number of states near the first mass level" is needed to isolate VS (2607.14230; also 2610.02124: "finite-spin support at the lowest massive pole is … an additional assumption").
7. Other routes:
   - Analytic dual resonance plus an integer spectrum fixes the scalar string amplitude for n ≥ 5 (2610.01101).
   - Seven-point factorization fixes the bosonic intercept among Koba–Nielsen-type amplitudes with minimal degeneracy, and the simplest factorizable Gross satellite deformation is excluded by unitarity (2603.04485).
   - 2312.07652 excludes bespoke and binary-geometry deformations.
- **OUR CHECK (conditional on 2601.11705 eq. 1.6).** The *planar* λ-family A_λ of 2408.03362 (eq. 18) violates the ρ-independent requirement a₂₁/a₂₀ = 1/4 except at λ=1. For λ = 0.25, 0.5, 0.9, 1, 1.1, 2, 4 the ratio is 0.0241, 0.0925, 0.2201, **0.25000000**, 0.2785, 0.4746, 0.6831 (two independent summations agree). So no λ≠1 member can be the 4-point amplitude of a planar N=4 theory with tree-level 6-point factorization. This is probably a corollary of 2610.02124.

**Numerical "spectrum-input" islands (ESTABLISHED numerically; interpretation CONJECTURE)**
- 2406.03543:
  - Inputs: a scalar at the gap and the next state at 2M². The corner sits at ā₁₀ ≈ 0.7336, ā₂₀ ≈ 0.6598 (eq. 1.10), versus Veneziano ζ₃/ζ₂ = 0.73076 and ζ₄/ζ₂ = 0.65797 (1.11; OUR CHECK).
  - Fixing |ḡ₀|² = 1/ζ₂ (3.4) gives an island.
  - There is a 1-parameter family of corners, M²_ℓ = (μ₂−1)ℓ+1 (1.16), so Veneziano is **not** singled out.
- 2412.13368 conjectures a convex-hull description and finds islands around Dirac–Born–Infeld (DBI) and VS.
- 2406.12959: Type II is *strictly inside* the allowed region. Its extremal "single trajectory" solutions are shown spurious in 2512.17828.
- 2510.07991: maximizing higher-spin couplings selects linear trajectories; with gravity, the trajectory passes through the graviton pole.
- 2606.19283: a VS island appears once a Virasoro-inspired nonlinear ansatz is imposed.
- 2604.15235: non-projective gravitational bounds, e.g. M/M_P ≲ 7.8 in D=5.

**Bottom line.** Every uniqueness theorem has one "UV-shaped" input that is not derived from IR consistency: level truncation, ultrasoftness or minimal zeros, monodromy or splitting, maximal SUSY plus multi-point factorization, and (for gravity) finite spin at the first level. The non-perturbative 10D bootstrap only shows that string theory is *compatible* with the boundary.

## 2. Three sharpest open questions

1. **The IST loophole in gravity.** Take 4d N=8 (or 10D IIB) tree-level EFTs with 4/5/6-point factorization, peculiar parity and positivity (2607.14230). Does some further *consistency* condition exclude the infinite-spin-tower corner, without assuming finite spin at the first massive level? Candidates are 7–8-point factorization, multipositivity of 5-point residues, and locality of the IST's 5-point extension. *Answered by:* a proof that the IST fails such a condition, or an explicit consistent IST 5- and 6-point amplitude.
2. **Is α_min^Boot = E_{3/2}(e^{iπ/3})/64 = 0.138949 exactly in 10D?** GPV gives 0.13 ± 0.02 (2102.02847 eq. 10); 2212.00151 says "close but not equal". *Answered by:* rigorous primal and dual brackets with width < 10⁻³ around 0.1389. A primal amplitude below 0.1389 would mean non-string theories exist; convergence to 0.1389 would mean string theory saturates the bound.
3. **Which UV input is necessary (non-SUSY, planar)?** Do the surviving 4-point deformations admit factorizing, positive 5- and 6-point extensions?
   - The survivors: the HZ unitary exponentially soft 3-parameter deformation (2311.13631); the CHR (r,w) family (2406.02665 eq. 18); and the 1, 5, 11, 18, 27 linearized partial-level-truncation solutions for λ_pt = 0..4 (text after 2406.02665 eq. 21).
   - *Answered by:* an extension, or an obstruction for each. "None survive" would make *higher-point factorization + positivity* the minimal principle that selects the string.

## 3. Key papers (all verified this session)

| arXiv | Authors | Yr | Title | Ver. | Why |
|---|---|---|---|---|---|
| hep-th/0602178 | Adams, Arkani-Hamed, Dubovsky, Nicolis, Rattazzi | 2006 | Causality, Analyticity and an IR Obstruction to UV Completion | Y | Origin of positivity |
| 2011.02957 | Caron-Huot, Van Duong | 2020 | Extremal Effective Field Theories | Y | Sum rules, null constraints, LP (reproduced) |
| 2012.15849 | Arkani-Hamed, T.-C. Huang, Y.-t. Huang | 2020 | The EFT-Hedron | Y | Positive geometry of coefficients |
| 2102.08951 | Caron-Huot, Mazáč, Rastelli, Simmons-Duffin | 2021 | Sharp Boundaries for the Swampland | Y | Gravity bounds; g₀ ≤ 3.000 vs 2ζ₃ |
| 2102.02847 | Guerrieri, Penedones, Vieira | 2021 | Where is String Theory? | Y | Non-perturbative α_min |
| 2212.00151 | Guerrieri, Murali, Penedones, Vieira | 2022 | Where is M-theory in the space of scattering amplitudes? | Y | 9/10/11D; small gap to string |
| 2202.08280 | Häring, Zhiboedov | 2022 | Gravitational Regge bounds | Y | Justifies twice-subtracted dispersion |
| 1607.04253 | Caron-Huot, Komargodski, Sever, Zhiboedov | 2016 | Strings from Massive Higher Spins… | Y | Asymptotic uniqueness |
| 2210.12163 | Cheung, Remmen | 2022 | Veneziano Variations… | Y | 4-point deformations exist |
| 2302.12263 | Cheung, Remmen | 2023 | Stringy Dynamics from an Amplitudes Bootstrap | Y | Infinite-parameter unitary family |
| 2311.13631 | Häring, Zhiboedov | 2023 | The Stringy S-matrix Bootstrap: Maximal Spin and Superpolynomial Softness | Y | Softness leaves no IR imprint; deformation |
| 2312.07652 | Arkani-Hamed, Cheung, Figueiredo, Remmen | 2023 | Multiparticle Factorization and the Rigidity of String Theory | Y | Higher points kill deformations |
| 2310.10710 | Chiang, Huang, Weng | 2023 | Bootstrapping string theory EFT | Y | Monodromy + unitarity → Type I, D=10 |
| 2310.10729 | Berman, Elvang, Herderschee | 2023 | Flattening of the EFT-Hedron… | Y | SUSY flattening; monodromy |
| 2406.02665 | Cheung, Hillman, Remmen | 2024 | Bootstrap Principle for the Spectrum and Scattering of Strings | Y | Veneziano unique (level truncation) |
| 2408.03362 | Cheung, Hillman, Remmen | 2024 | Uniqueness Criteria for the Virasoro-Shapiro Amplitude | Y | VS unique; λ-family |
| 2406.03543 | Berman, Elvang | 2024 | Corners and Islands in the S-matrix Bootstrap of the Open Superstring | Y | Spectrum-input islands |
| 2406.12959 | Albert, Knop, Rastelli | 2024 | Where is tree-level string theory? | Y | Type II interior; extremal trajectory |
| 2405.21057 | Eckner, Figueroa, Tourkine | 2024 | On the number of Regge trajectories for dual amplitudes | Y | No finite-trajectory amplitudes |
| 2410.21499 | Häring, Zhiboedov | 2024 | What is the graviton pole made of? | Y | Eikonal vs stringy; Tauberian theorem |
| 2508.09246 | Cheung, Remmen, Sciotti, Tarquini | 2025 | Strings from Almost Nothing | Y | Regge zeros + minimal zeros ⇒ string |
| 2601.11705 | Elvang, Herderschee, Morales | 2026 | String Theory from Maximal Supersymmetry | Y | Nonlinear relations (1.6) ⇒ Veneziano |
| 2607.14230 | Berman, Caron-Huot, Chandra, Elvang, Herderschee, Lin, Morales | 2026 | Gravitational EFTs with Maximal SUSY and a Peculiar Parity | Y | VS vs IST corners |
| 2610.02124 | Y.-t. Huang, Wan, Wang, Zhou | 2026 | Infrared Consistency and the Uniqueness of String Amplitudes | Y | All-order proof; gravity gap |
| 2605.11084 | Wan, Zhou | 2026 | Analytic Bootstrap of the Veneziano Amplitude | Y | Analytic uniqueness (monodromy/splitting) |

Also verified and cited above: 2011.00037, 2011.02400, 2201.06602, 2203.11950, 2201.11575, 2502.20372, 2608.25183, 2412.13368, 2510.07991, 2603.04485, 2606.19283, 2604.15235, 2512.17828, 2401.08736, 2610.01101, 2607.27300, 2308.03833, 2312.16282, 2008.02293, 2208.12233.

## 4. Computable handles (Python, a few hours each)

**H1 — Scalar-EFT linear program, with the VS-deformation curve overlaid.**
- *Done:* (3.10–3.11) of 2011.02957 reproduced (`checks.py`). Also computed are the 2408.03362 eq. (10) coefficients for the λ-family at ρ=0, cross-checked against the direct dual-resonant sum, with crossing verified (`chr_vs_family.py`). Units: M = μ(1) = 1.

  | λ | g̃₃ | g̃₄ |
  |---|---|---|
  | 0.5 | 2.92211 | 0.49296 |
  | 1 (VS minus graviton pole) | 2.78697 | 0.48622 |
  | 2 | 2.58218 | 0.48243 |
  | 5 | 2.31538 | 0.48574 |

  The limits from eqs. 14–15 are (3, ½) as λ→0 and (2, ½) as λ→∞ (OUR expansion).
- *To do:* generate null constraints from (2.19) up to k=8, compute the (g̃₃, g̃₄) region in d=4 and d=10, and overlay the curve.
- *New if* the curve tracks the numerical boundary over an interval. 2408.03362 says that boundary "has thus far eluded analytic understanding".

**H2 — Nonlinear SUSY tests of candidate deformations.**
- *Done:* the planar λ-family fails the 2601.11705 eq. (1.6) ratio a₂₁/a₂₀ = 1/4 unless λ=1 (`ehm_test_planar_family.py`, `ehm_ratio_crosscheck.py`).
- *To do:* the same test on Coon, the Häring–Zhiboedov deformation and the Cheung–Remmen hypergeometric amplitudes. Then read the gravity relations of 2607.14230 and test the closed-string family M_λ (eq. 7) and the IST (eq. 15).
- *New if:* an explicit exclusion catalogue, or a surviving non-string point.

**H3 — Gravity g₀ linear program with spectral input.**
- Reproduce 0 ≤ g₀M⁶/8πG ≤ 3.000 (2102.08951 eq. 3.41) using the 14 functionals described there: C₋₂ (3.36) times p²…p¹⁰, (3.40) times four powers, and the null constraint X₂.
- Then add VS-like inputs: integer levels m² = nM² and J ≤ 2n−2. Check whether 2ζ₃ = 2.404 (3.42) becomes extremal or a kink in (g₀, g₂), with VS at (2ζ₃, ζ₅, 2ζ₃²) (OUR CHECK).
- *New if:* a 10D gravitational analogue of the Veneziano corner without coupling input.

**H4 — Predicted string region for a two-coefficient GPV bootstrap.**
- Map the fundamental domain under τ ↦ (E_{3/2}/64, c·E_{5/2}), with c fixed by matching the tree-level VS coefficients (2ζ₃, ζ₅). The E_{5/2} normalization is not checked in a primary source and must be re-derived.
- *Done:* on the sampled boundary of the fundamental domain, both functions are minimized at τ = e^{iπ/3}: E_{3/2} = 8.892745 and E_{5/2} = 4.719492 (`gpv_eisenstein.py`).
- *New:* an explicit boundary curve and corner, giving a falsifiable target for question 2.

**H5 — Positivity windows.**
- Use the closed-form partial-wave formulas of 2408.03362v1 (footnotes to eq. 13) to push the λ-family positivity scan from n=30 to n ~ 300, with an asymptotic large-n analysis to find D_crit(λ).
- Separately, test the claim in 2102.08951 (footnote 14) that the residues of (3.21) are Gegenbauer-positive for D ≤ 24.148… by adapting the Bessel-recurrence method of 2608.25183.
- *New:* a proof, or the exact critical D.
