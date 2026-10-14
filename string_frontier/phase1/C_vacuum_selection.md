# FRONT C — Vacuum selection: landscape counting and swampland conjectures (Phase 1 map)

Labels: **ESTABLISHED** (proved, or a controlled calculation under stated assumptions), **CONJECTURE**, **OUR ATTEMPT** (computed here; code `phase1/C_checks/C_vacuum_selection_checks.py`, output `phase1/C_checks/C_vacuum_selection_checks.out`).
"Verified" means I fetched the arXiv abstract page in this session and the title and authors matched. Equation numbers are given only where I read them in the fetched text.

> **ID correction for the brief:** arXiv:1903.01148 is *not* Palti's review. It is a coding-theory paper ("Asymmetric Single Magnitude Four Error Correcting Codes"). Palti's "The Swampland: Introduction and Review" is **arXiv:1903.06239** (verified).

## 1. State of the art

### 1a. Landscape counting
- **ESTABLISHED (heuristic model).** Bousso–Polchinski (hep-th/0004134): λ = λ_bare + ½Σ n_i² q_i² (eq. 2.21). An acceptable λ needs a lattice point in a thin shell (2.22–2.24), so the spacing is Δλ_min = D·Πq_i / (ω_{J−1}|2λ_bare|^{J/2−1}) (2.26). One flux leaves far too large a gap (2.19–2.20); BP say J≈100, q≈10^−1.6 works.
  - **OUR ATTEMPT:** (2.26) with Δλ = 10^−120, λ_bare = 1, J = 100 gives q = 10^−1.43 (D = 1) to 10^−1.73 (D = 2^J). Brute-force lattice counts for J = 3–5 match (2.24) within 13%.
- **ESTABLISHED (asymptotic, continuous-flux approximation).** Ashok–Douglas (hep-th/0307049).
  - I_vac(L ≤ L_max) = (2πL_max)^K/(π^{n+1}K!) ∫_{F×H} det(−R−ω), K = b₃ = 2n+2 (eq. 1.5). The tadpole is η(N_RR,N_NS) = L (3.3), and η > 0 at supersymmetric vacua (3.9–3.10).
  - Reliable only for L ≫ K; AD suspect failure for 1 ≪ L ≪ K.
  - **OUR ATTEMPT:** the "10^500" is the prefactor: log₁₀[(2πL)^K/K!] = 525 at K = 300, L = 1000. This leaves out the moduli-space integral, which can be tiny (see the Siegel volumes, 2.25).
- **ESTABLISHED (continuous-flux approximation).** Denef–Douglas (hep-th/0404116).
  - N(L ≤ L*) = (2πL*)^{2m}/(2m)!·|det η|^{−1/2}∫d^{2m}z det g ρ(z) (3.15–3.16).
  - For n = 1, ρ = π^{−2}(2−|F|²+2|F|³/√(4+|F|²)) (3.19). The index density is π^{−m}det(R+ω) (3.23). The tadpole bound is L* = χ/24 (2.4–2.5).
- **OUR ATTEMPT (new check).** Rigid-CY toy of DD §5.1: W = Aτ+B, A,B ∈ ℤ[i], τ̄ = −B/A. I read this from a verbatim reproduction, not the arXiv PDF §5.
  - The exact SL(2,ℤ)-orbit count is N(L) = Σ_{n≤L}σ₁(n) → π²L²/12 (ratio 1.0004 at L = 3000).
  - AD (1.5) and DD (3.15), evaluated by me, give π²L²/6: exactly twice the exact count. Likely cause: the flux integral counts N and −N separately, although the central element −1 ∈ SL(2,ℤ) identifies them at fixed τ.
  - τ is uniform in the hyperbolic measure within 2%; 120 of 18604 vacua sit at τ = i (L = 150).
- **ESTABLISHED (numerical, 2025).** Chauhan–Cicoli–Krippendorf–Maharana–Piantadosi–Schachner (2501.03984) systematically enumerated vacua on a two-modulus CY. They fit N_vac ∝ N_max^{5.80±0.05}, against the continuum prediction (2πN_max)^6/6! · ∫… (their eq. 42). They report "local deviations from statistical expectations" (abstract).
- **ESTABLISHED (theorem).** Bakker–Grimm–Schnell–Tsimerman (2112.06995): with bounded tadpole, F-theory self-dual flux vacua (connected components of the vacuum locus) are finite. The proof uses definability of period maps in ℝ_an,exp (2112.08383 eq. 4.7). It covers complex structure and fluxes on one fixed CY₄, not Kähler moduli or the number of geometries.
- **CONJECTURE.** Grimm's tameness conjecture (2112.08383): consistent EFTs are definable in an o-minimal structure, ℝ_an,exp in the strong form. Grimm–Monnee (2311.09295) state three precise counting, complexity and dimension conjectures for the vacuum locus, using ∫G₄∧⋆G₄ ≤ L (eq. 3.2).
- **ESTABLISHED (explicit).** Grimm–van de Heisteeg (2404.12422): the first exact landscape of W = 0 vacua on a genuine CY₄. Along symmetry loci the vacua become dense if the tadpole is unbounded. With the tadpole bound imposed, only 10 vacua remain (abstract and text).
- **ESTABLISHED (estimate, not a count of controlled vacua).** Taylor–Wang (1511.03209): AD-type estimate of 𝒪(10^272,000) flux vacua on a single elliptic CY₄ (base B_max). Gauge group E₈⁹×F₄⁸×(G₂×SU(2))¹⁶; the Standard Model "appears to arise … only from a broken E₈" (abstract).
- **ESTABLISHED (constructive).** Halverson–Long–Sung (1706.02299) construct (4/3)×2.96×10^755 F-theory geometries.
- **ESTABLISHED (constructive).** Cvetič–Halverson–Lin–Liu–Tian (1903.00009): about 10^15 globally consistent F-theory models with the exact chiral Standard Model spectrum and gauge-coupling unification. Moduli stabilization and Λ are not addressed.
- **ESTABLISHED (theorem).** Denef–Douglas (hep-th/0602072): finding a Bousso–Polchinski vacuum that matches a target Λ is NP-complete.
- **CONJECTURE, with strong evidence.** Tadpole conjecture, Bena–Blåbäck–Graña–Lüst (2010.10519): n_stab < βN_flux. Since N_flux ≤ χ/24 ≈ h^{3,1}/4, any β < 4 forbids full F-theory stabilization; the refined form says β < 3 (2407.16758 eqs. 2.10–2.11).
  - Support: asymptotic Hodge theory fixes one modulus per sl(2) block (2204.05331).
  - Against the refined form: the non-geometric 2⁶ Gepner model reaches β = 3.4 and has fully massive vacua within the bound N_flux ≤ 40 (2407.16758, eq. 3.20).
- **ESTABLISHED (explicit, with stated mild assumptions).** Demirtas–Kim–McAllister–Moritz–Rios-Tascon (2107.09064): supersymmetric AdS₄ vacua with |Λ| < 10^−123.
- **CONJECTURE.** McAllister–Moritz–Nally–Schachner (2406.13751) build KKLT de Sitter vacua at leading order in α′ and g_s. Showing they survive subleading corrections is "an important open problem" (abstract). KKLT itself: hep-th/0301240.

### 1b. Anthropics and measures
- **ESTABLISHED (conditional calculation).** Weinberg, PRL 59 (1987) 2607 (not on arXiv), refined by Martel–Shapiro–Weinberg (astro-ph/9701099). With a flat prior on ρ_V, 1–3 times the matter density is "not very unlikely" (abstract). The result depends only on σ³ρ̄ at recombination. Equation numbers not checked.
- **CONJECTURE (proposals).** Causal-patch measure (hep-th/0605263); causal entropic principle (hep-th/0702115), whose Λ prediction is in "excellent agreement". Linde/Vilenkin volume weighting and the Hartle–Hawking and tunneling wavefunctions (no arXiv IDs) remain alternatives. No regulator-independent measure is established.
- **CONJECTURE (alternative to selection).** Kaloper–Padilla sequestering (1309.6562) decouples Standard Model vacuum energy from gravity. It predicts that w ≈ −1 is transient and that the universe will collapse (abstract).

### 1c. Swampland: which statements are theorems
- **ESTABLISHED (theorem within AdS/CFT).** Harlow–Ooguri (1810.05337, 1810.05338): no bulk global symmetries; every irreducible representation of a gauge group is realized (completeness); gauge groups are compact. Outside AdS these are **CONJECTURE**.
- **ESTABLISHED.** Adams–DeWolfe–Taylor (1006.1352): the 10d N=1 supergravities with U(1)^496 and E₈×U(1)^248 are inconsistent. So every 10d supergravity without a known inconsistency is realized in string theory.
- **CONDITIONAL, resting on CONJECTURE.** Montero–Vafa (2008.11729) assume the cobordism conjecture (McNamara–Vafa, 1909.10355, itself derived from no global symmetries). With 16 supercharges and d > 6, this bounds gauge rank, and in d = 8, 9 the only solutions are the string-theoretic ones. Kim–Shiu–Vafa (1905.08261) use probe-string unitarity to exclude some 10d and 6d theories.
- **CONJECTURE: Weak Gravity Conjecture** (hep-th/0601001): Λ_UV ≲ g·M_Pl, and there is a state with m/q ≤ (m/q) of an extremal black hole.
- **CONJECTURE: Distance Conjecture** (hep-th/0605264): towers with m ~ e^{−aΔφ} (OPSV eq. 4).
  - Sharpened version: λ ≥ 1/√(d−2) (2206.04063).
  - Emergent String Conjecture (1910.01135): every infinite-distance limit is a decompactification or a weakly coupled critical string limit. It implies the taxonomy rules ζ_a·ζ_b = 1/(d−2) + δ_ab/n_a and ζ·Z_QG = 1/(d−2) (2405.20332 eq. 2.10).
- **CONJECTURE: de Sitter conjectures.**
  - OOSV (1806.08362): |∇V| ≥ cV (eq. 1.1). For an exponential potential the late-time attractor has w = −1 + λ²/3 (1.3), so acceleration needs λ < √2.
  - Refined form (1810.05506): either |∇V| ≥ cV (2) or min ∇∇V ≤ −c′V (3).
  - TCC (1909.11063): V(φ₂) ≤ A·exp(−2(φ₂−φ₁)/√((d−1)(d−2))). In d = 4 this means a slope ≥ √(2/3) ≈ 0.816, averaged over large field ranges, and forbids long-lived metastable dS.
  - Strong asymptotic form: |∇V|/V ≥ 2/√(d−2) = √2 in d = 4 (2101.11617). It "precludes accelerated expansion".
- **CONJECTURE: AdS Distance Conjecture and scale separation** (1906.05225): towers with m ~ |Λ|^α, α = ½ for supersymmetric AdS.
- **CONJECTURE: Festina Lente** (1910.01648): charged black holes in dS must evaporate, which gives a lower bound on the masses of charged particles. The bound is satisfied by the Standard Model. The formula was not read.
- **CONJECTURE (scenario built on conjectures plus data): Dark dimension** (2205.12293).
  - Starting point: m ~ |Λ|^α with 1/d ≤ α ≤ ½ (eqs. 2.2, 2.8). Torsion-balance bounds give m ≳ 6.6 meV (3.2). Together these force α = ¼, so Λ ~ m⁴ (3.3–3.4).
  - Result: one extra dimension, l ~ λΛ^{−1/4} with λ ~ 10^−1–10^−3, and species scale M̂ = m^{1/3}M_Pl^{2/3} (eq. 2.1 with n = 1).
  - **OUR ATTEMPT:** with the reduced Planck mass, l = 0.085–8.5 μm and M̂ = 5×10^8–2.4×10^9 GeV (paper: 10^9–10^10 GeV; Planck-mass convention). Λ^{−1/4} = 85 μm (paper: ≈ 88 μm).
- **2025 data versus swampland.**
  - DESI DR2 (2503.14738) strengthened the case for dynamical dark energy, as summarized in the 2504.04226 abstract; I did not read the DESI text.
  - Akrami–Alestas–Nesseris (2504.04226) fit exponential quintessence with λ = 0.698 (+0.173, −0.202), a ~3.5σ preference for λ ≠ 0.
  - **OUR ATTEMPT:** using the upper error bar, the TCC average bound √(2/3) is 0.7σ away. The strong bound √2 is 4.1σ away, so it cannot hold for today's potential. It may still hold only asymptotically, as 2101.11617 itself allows.
  - Bedroya–Obied–Vafa–Wu (2507.03090) fit the dark-dimension radion with a coupling to dark matter: c′ ≈ 0.05 ± 0.01, with |c| close to √(2/3). Their reading of the apparent w < −1 is evidence for a dark-matter–dark-energy coupling.

## 2. The three sharpest open questions

1. **Does controlled metastable de Sitter exist?** Is there a 4d string compactification with every modulus stabilized, V > 0, a lifetime longer than the TCC bound τ ≲ H^{−1} log(1/H) (as quoted in 2507.03090), and parametric control over α′ and g_s corrections?
   - An explicit vacuum with computed subleading corrections, for example completing 2406.13751, refutes TCC and OOSV.
   - A no-go theorem covering KKLT/LVS-type uplifts settles it the other way.
2. **What is the asymptotic slope bound?** Does any controlled string compactification have an infinite-distance region with |∇V|/V < √2 in d = 4, which would allow eternal acceleration?
   - The answer is either a theorem over all asymptotic limits (for example using the taxonomy polytopes of 2405.20332) or an explicit example.
   - It matters for observation because DESI-era fits give λ ≈ 0.7.
3. **Is the discretuum really uniform near Λ = 0?** For a fixed CY orientifold with all moduli stabilized and N_flux ≤ L*, does #{vacua with |Λ| < ε} / #{vacua} tend to c·ε, with c finite and nonzero, down to ε ~ 10^−122?
   - This is the flat prior Weinberg-type arguments need. It must hold where h^{3,1} ≫ L*, exactly where AD (1.5) and DD (3.15) are untrustworthy and the tadpole conjecture may forbid full stabilization.
   - Settled by an equidistribution theorem for self-dual lattice points, or an explicit arithmetic obstruction (gaps).

## 3. Key papers

| arXiv ID | Authors | Year | Title (short) | Verified? | Why it matters |
|---|---|---|---|---|---|
| hep-th/0004134 | Bousso, Polchinski | 2000 | Quantization of four-form fluxes… | yes | Discretuum; Λ spacing (2.26) |
| hep-th/0307049 | Ashok, Douglas | 2003 | Counting flux vacua | yes | Index formula (1.5) |
| hep-th/0404116 | Denef, Douglas | 2004 | Distributions of flux vacua | yes | Vacuum density (3.15–3.23); rigid-CY toy |
| hep-th/0602072 | Denef, Douglas | 2006 | Computational complexity of the landscape I | yes | NP-completeness of finding a target Λ |
| 1511.03209 | Taylor, Wang | 2015 | F-theory geometry with most flux vacua | yes | 10^272,000 estimate |
| 1706.02299 | Halverson, Long, Sung | 2017 | Algorithmic universality in F-theory | yes | ~4×10^755 geometries |
| 1903.00009 | Cvetič, Halverson, Lin, Liu, Tian | 2019 | A quadrillion Standard Models from F-theory | yes | Counts of Standard-Model-like vacua |
| 2010.10519 | Bena, Blåbäck, Graña, Lüst | 2020 | The tadpole problem | yes | Tadpole conjecture |
| 2204.05331 | Graña, Grimm, van de Heisteeg, Herraez, Plauschinn | 2022 | Tadpole conjecture in asymptotic limits | yes | Mechanism for linear scaling |
| 2407.16758 | Becker, Brady, Graña, Morros, Sengupta, You | 2024 | Tadpole conjecture in non-geometric backgrounds | yes | β = 3.4 > 3; fully massive vacua |
| 2112.06995 | Bakker, Grimm, Schnell, Tsimerman | 2021 | Finiteness for self-dual classes… | yes | Finiteness theorem |
| 2112.08383 | Grimm | 2021 | Taming the landscape of effective theories | yes | Tameness conjecture |
| 2311.09295 | Grimm, Monnee | 2023 | Finiteness theorems and counting conjectures | yes | Precise counting conjectures |
| 2404.12422 | Grimm, van de Heisteeg | 2024 | Exact flux vacua, symmetries… | yes | Exact landscape; density on symmetric loci |
| 2501.03984 | Chauhan, Cicoli, Krippendorf, Maharana, Piantadosi, Schachner | 2025 | Deep observations of the IIB flux landscape | yes | Numerical test of DD scaling |
| 2107.09064 | Demirtas, Kim, McAllister, Moritz, Rios-Tascon | 2021 | Small cosmological constants in string theory | yes | AdS with |Λ| < 10^−123 |
| 2406.13751 | McAllister, Moritz, Nally, Schachner | 2024 | Candidate de Sitter vacua | yes | Explicit KKLT candidates |
| astro-ph/9701099 | Martel, Shapiro, Weinberg | 1997 | Likely values of the cosmological constant | yes (PDF and ADS) | Anthropic distribution |
| hep-th/0605263 | Bousso | 2006 | Holographic probabilities in eternal inflation | yes | Causal-patch measure |
| 1903.06239 | Palti | 2019 | The Swampland: introduction and review | yes | Review (the brief's 1903.01148 is wrong) |
| 1810.05337 | Harlow, Ooguri | 2018 | Constraints on symmetry from holography | yes | Theorem: no global symmetries |
| 1806.08362 / 1810.05506 | OOSV / OPSV | 2018 | dS conjecture / refined dS conjecture | yes | dS conjectures (1.1), (2)–(3) |
| 1909.11063 | Bedroya, Vafa | 2019 | Trans-Planckian censorship and the swampland | yes | TCC slope bound |
| 1910.01135 / 2405.20332 | Lee, Lerche, Weigand / Etheredge et al. | 2019/24 | Emergent strings; taxonomy of limits | yes | Emergent String Conjecture and its rules |
| 2205.12293 | Montero, Vafa, Valenzuela | 2022 | The dark dimension and the swampland | yes | Λ linked to a micron dimension |
| 2504.04226 / 2507.03090 | Akrami et al. / Bedroya et al. | 2025 | DESI exponential quintessence / evolving dark sector | yes | λ ≈ 0.7 against swampland bounds |

Also verified but not in the table: 2102.01111, 2212.06187, hep-th/0601001, hep-th/0605264, 1909.10355, 1810.05338, 2008.11729, 1006.1352, 1905.08261, 1906.05225, 2206.04063, 1910.01648, 2101.11617, hep-th/0301240, hep-th/0702115, 1309.6562, 2110.05511, 2503.14738.

## 4. Computable handles (Python, a few hours each)

1. **Bousso–Polchinski discretuum: level statistics and search hardness.**
   - Compute: incommensurate q_i for J = 4–12, and J ≈ 30–100 with M-theory charges (2.37). Enumerate the shell (2.22) with a meet-in-the-middle sum of squares. Measure the count against (2.24)–(2.26), the nearest-neighbour spacing distribution of λ, and how runtime scales when finding |λ| < ε.
   - Compare: BP's q ≈ 10^−1.6 at J = 100 (we get 10^−1.43 to 10^−1.73); NP-completeness (hep-th/0602072).
   - New: a measured spacing law at λ → 0, a direct test of the flat prior (open question 3), and an empirical hardness curve.
2. **Exact versus continuum flux counting, and the Λ distribution.**
   - Compute: extend our rigid-CY enumeration (Σσ₁, continuum overcounts by exactly 2) to DD's symmetric T⁶/ℤ₂ toy with L* = 16 (DD §3.1.3; density 1/24, eq. 3.21). Enumerate all vacua and Λ̂ = e^K|W|² (DD 1.5). Histogram Λ̂ near 0 and the moduli at enhanced-symmetry points.
   - Compare: DD (3.15), (3.19)–(3.21), (3.27); the N_max^{5.80} fit of 2501.03984.
   - New: exact finite-L corrections, the fraction of vacua at symmetric points, and whether W₀ → 0 is approached uniformly or with arithmetic gaps.
3. **Anthropic probability (Weinberg / Martel–Shapiro–Weinberg) with modern inputs.**
   - Compute: the collapsed fraction (Press–Schechter or MSW spherical collapse) with a Planck-2018 σ(M). Get dP/dρ_V for a flat prior, then P(ρ_V ≤ ρ_obs). Redo it with the measured discretuum density from handles 1–2 as the prior, and with quintessence at λ ≈ 0.7 instead of Λ.
   - Compare: MSW's 1–3× ρ_m being "not very unlikely" (astro-ph/9701099); the causal-entropic prediction (hep-th/0702115).
   - New: whether a swampland-compatible rolling potential changes the anthropic likelihood of today's dark-energy density.
4. **Swampland slope versus data: convex-hull test.**
   - Compute: (a) map exponential-quintessence λ to (w₀, w_a). (b) For V = Σ A_i e^{−k_i·φ}, the late-time slope is the distance from the origin to the convex hull of {k_i}. Find the minimum over the 10d→4d dilaton–volume exponent vectors of type II sources (curvature, H, F_p, O-planes, Casimir), each derived by hand.
   - Compare: OOSV (1.3), TCC √(2/3), strong bound √2 (2101.11617), DESI fit λ = 0.698 (our tension numbers: 0.7σ and 4.1σ).
   - New: a source-class minimal slope set against the DESI region (open question 2).
5. **Tadpole versus moduli scaling in solvable toys.**
   - Compute: for torus products or K3-type sl(2) blocks (2204.05331), use integer programming or lattice reduction to find the minimal N_flux giving a full-rank mass matrix. Extract β(n) for n ≲ 50, at generic points and on symmetric loci (2404.12422).
   - Compare: 44% in K3×K3 (2010.10519); β = 3.4 and 2.74 (2407.16758); critical β = 4 (2.11).
   - New: any family with β ≥ 4 refutes the F-theory consequence of the tadpole conjecture.
6. **Dark-dimension consistency.**
   - Compute: the 5d one-loop Casimir energy on S¹ as a function of field content, and the λ⁴ = Λ/m⁴ it implies, against eq. (3.5): m^{1/2} ≤ λ⁴ ≲ 1.
   - New: which fermion excess is required, to set against 2205.12293 and 2507.03090.
