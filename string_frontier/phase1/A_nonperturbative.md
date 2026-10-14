# Front A — Non-perturbative formulations of string/M-theory (Phase 1 map)

Date 2026-10-14. "v" = abstract page fetched this session, title and authors match. [UNVERIFIED-ID] = seen only through search. Equation numbers appear only where we read them. Labels: ESTABLISHED / CONJECTURE / OUR ATTEMPT. Every OUR ATTEMPT item was run, and the scripts and outputs are in `phase1/A_code/`.

## 1. State of the art

**BFSS / DLCQ**
- CONJECTURE: Uncompactified 11d M-theory equals the N→∞ limit of D0-brane quantum mechanics (hep-th/9610043 v). There is no proof. Its finite-N, 11d-Lorentz-invariant content is the open core.
- ESTABLISHED (derivation under stated assumptions): Finite N equals DLCQ M-theory with N units of p₋ (Susskind hep-th/9704080 v). Seiberg derives this as the limit of a boosted, vanishing spacelike circle, for T^p with p≤5 (hep-th/9710009 v). Sen does the same and finds new subtleties at p=6 (hep-th/9709220 v). Hellerman–Polchinski show that the lightlike limit diverges order by order because of longitudinal zero modes. It "generally exists nonperturbatively but is more complicated" (hep-th/9711037 v).
- ESTABLISHED: Eikonal tests at q₁₁=0. DLCQ supergravity gives L′ = N₂v²/2R + (15/16)N₁N₂v⁴/(R³M⁹r⁷) + (225/64)N₁²N₂v⁶/(R⁵M¹⁸r¹⁴) + … (BBPT hep-th/9706072 v, eq. 13, using h₋₋ from eq. 8). The one-loop matrix result matches the first term. The two-loop SU(2) value 225/32 (eq. 14, Table 1) times (N₁N₂²+N₁²N₂)/2 (eq. 15) matches the second. The one-loop c₁₂ (v⁶/r¹¹) vanishes (BBPT text). The v⁴/r¹⁰ term is not renormalized at two loops (hep-th/9705091 v). Two-loop three-graviton terms agree with supergravity (Okawa–Yoneya hep-th/9806108 v; Taylor–Van Raamsdonk hep-th/9806066 v).
- OUR ATTEMPT: We reproduced 15/16 and c₁₂=0 from the one-loop mode sum f(x)=6+2cosh2x−8cosh x. Higher one-loop coefficients come out as 9009/4096 (v⁸/r¹⁵) and −182325/16384 (v¹⁰/r¹⁹). We also derived the three-loop DLCQ-supergravity target from eq. 13: **(16875/1024) N₁³N₂ v⁸/(R⁷M²⁷r²¹)**. See `check_d0_phase_shift.py`.
- ESTABLISHED (tension): At three loops, matrix theory gives v⁶ terms in four-graviton scattering with no apparent supergravity counterpart (Dine–Echols–Gray hep-th/9810021 v). Compactifications on curved spaces fail to reproduce the supergravity v⁴ term (Douglas–Ooguri–Shenker hep-th/9702203 v). Finite-N matrix theory therefore matches supergravity only for protected terms. The N→∞ statement has not been tested beyond them.
- ESTABLISHED, with conditions: Herderschee–Maldacena compute the three-graviton amplitude with nonzero p₋ for all legs via a supersymmetric index. It agrees with 11d at finite N, and v3 matches the overall coefficient (2312.12592 v). They then assume a smooth large-N limit plus factorization, and soft theorems follow, implying Lorentz invariance of four-point amplitudes. Five or more points are covered only for restricted kinematics (2312.15111 v). Follow-ups: 2303.14200 v, 2510.15488 v, and 2602.11667 v, where localization reproduces the three-point momentum dependence.
- ESTABLISHED (numerical; 't Hooft/IIA regime): Lattice Monte Carlo reproduces E/N²=7.41T^{14/5} at 0.4≤T≤1.0 after large-N and continuum extrapolation (Berkowitz et al. 1606.04951 v). Earlier: 0707.4454 v, 1311.5607 v. A BMN-regulated precision test goes down to T=0.25λ^{1/3} (2210.04881 v). The confined phase of BFSS has been observed and proposed as an M-theory signature (2110.01312 v). None of this probes flat 11d scattering.
- ESTABLISHED: Bootstrap gives the first nontrivial lower bound on ⟨Tr X²⟩, within a factor ≈2 of Monte Carlo (Lin 2302.04416 v). Review: 2508.20970 v.
- ESTABLISHED: BMN mass deformation for DLCQ M-theory on the maximally supersymmetric pp-wave (hep-th/0202021 v). Dasgupta–Sheikh-Jabbari–Van Raamsdonk give the exact μ→∞ spectrum around fuzzy-sphere vacua (hep-th/0205185 v, eq. 1.3, Tables 1–2). The single-membrane coupling is g=(R/μN)^{3/2}.
- OUR ATTEMPT: We reproduced DSV Table 1 for N=2,3,4: bosons and fermions, zero modes N²−1, and Σm_B = Σm_F = (2/9)N(8N²+1) (`check_bmn_spectrum.py`). The fermion term coefficient is (R/2)·iμ/4R, not iμ/4.

**IKKT / type IIB matrix model**
- CONJECTURE: The large-N reduced 10d SYM is a constructive definition of type IIB string theory (hep-th/9612115 v).
- ESTABLISHED: Bosonic Yang–Mills integrals converge for D=3 (N>3), D=4 (N>2), and D≥5 (N≥2) (hep-th/9804199 v). The SU(N) supersymmetric index bulk term is Σ_{d|N}1/d² for D=10 and 1/N² for D=4,6 (hep-th/9803265 v; Monte Carlo check hep-th/9803117 v). The bosonic model shows no spontaneous breaking of Lorentz symmetry (hep-th/9811220 v).
- ESTABLISHED (semi-analytic and numerical, Euclidean): SO(10) breaks to SO(3). The Gaussian-expansion free energy is minimal at d=3 (1108.1293 v), and complex Langevin confirms it (2002.07410 v).
- Lorentzian history: KNT reported (3+1)-dimensional expansion (1108.1540 v). That configuration turned out to be Pauli-matrix dominated (2 large eigenvalues), an artifact of replacing e^{iS_b} by e^{βS_b} (1904.05914 v). Complex Langevin departs from the Pauli structure (1904.05919 v).
- ESTABLISHED (claim of the KEK group): In the conventional large-N limit the undeformed Lorentzian model is equivalent to the Euclidean one under contour deformation (2307.01681 v; 2604.19836 v, conclusions). A Lorentz-invariant mass term γ is needed to see Lorentzian signature. The bosonic (1+1)-dimensional expansion is a Lorentz-boost artifact (2407.03491 v). IR cutoffs "classicalize" the model artificially, and a fix is to gauge-fix the Lorentz symmetry (2404.14045 v).
- CONJECTURE (newest): With N≤128, an anisotropic mass ξ, and a real fermion mass m_f inspired by polarized IKKT (PTEP version, eq. 3.12), there is a phase with smooth, real (3+1)-dimensional expansion. The authors say the deformation "has to be carefully examined" (2604.19836 v, Apr 2026).
- ESTABLISHED (polarized IKKT): The supersymmetric mass deformation is dual to a Euclidean IIB background, with a D1 in a flux-supported cavity, and localizes to a moduli integral (Hartnoll–Liu 2409.18706 v). Supergravity solutions correspond one-to-one to fuzzy-sphere vacua (2410.18173 v). The original IKKT model differs from the vanishing-mass limit of polarized IKKT (2411.18678 v).
- CONJECTURE (Steinacker): On M^{3,1}×K with fuzzy K, the one-loop IKKT action contains a 4d Einstein–Hilbert term with G_N set by the KK scale of K, and vacuum energy does not act as a cosmological constant (2110.03936 v; 2303.08012 v). For the k=−1 FLRW bounce background, the one-loop analysis stabilizes early times. The 2025 revision (v2) says late times approach a(t)≈const with a decreasing dilaton, replacing the journal claim of acceleration followed by coasting a∝t (2310.11126 v). Dominance in the matrix integral is not shown.

**Closed string field theory**
- ESTABLISHED: Zwiebach built the full quantum closed bosonic SFT. Tree level is an L∞ algebra, and higher genus implies the BV master equation, with minimal-area vertices (hep-th/9206084 v). There are now consistent SFTs for all string theories, with unitarity given a momentum-contour prescription and UV finiteness (Sen–Zwiebach review 2405.19421 v; 1607.08244 v).
- ESTABLISHED, but local only: Background independence holds between nearby CFTs, classically (hep-th/9307088 v), quantum mechanically (hep-th/9311009 v), and for superstrings (1711.08468 v).
- ESTABLISHED: SFT fixes the normalization of type IIB D-instanton amplitudes, matching Green–Gutperle (2104.11109 v).
- ESTABLISHED: Hyperbolic vertices exactly solve the geometric master equation for L≤2 arcsinh 1 and tend to minimal-area vertices as L→∞ (1909.00033 v). The three-vertex local coordinates are known in closed form (Firat 2102.03936 v, eqs. 1.3, 1.5), and the tachyon cubic term is κ²V=−t²+(1/3)r⁶t³ (eqs. 1.7–1.8). A Mirzakhani-type recursion fixes all vertices from the cubic one, so the classical equation of motion becomes a quadratic integral equation (2409.02982 v; 2302.12843 v).
- OUR ATTEMPT: We reproduced Firat's 1.62187×10⁸ at λ=arcsinh(1)/π, and the λ→∞ limit 6561/4096 (`check_hyperbolic_vertex.py`). The reading that works uses continuous loggamma branches and ln r=(v−π/2)/λ.
- CONJECTURE: Closed-string tachyon vacuum. The cubic potential −t²+(6561/4096)t³ is spoiled by the quartic −3.0172t⁴. The ghost-dilaton coupling restores a critical point: d=2.2944t² gives V=4.4422t⁵ (Yang–Zwiebach hep-th/0506077 v, eqs. 1.1–1.4). This is proposed as a zero-action "nothing" state. Also: an AI-assisted algebraic grade expansion with unproven convergence (2603.29926 [UNVERIFIED-ID]).
- What SFT lacks: it expands in g_s around one CFT. Its string-field path integral has no non-perturbative definition (measure and contour are unknown), its background independence is only local, and it has no resurgent completion.

**M-theory**
- ESTABLISHED: ABJM is dual to M-theory on AdS₄×S⁷/Z_k (0806.1218 v). Its S³ partition function is an Airy function to all orders in 1/N (1106.4631 v; 1110.4066 v). From ABJM one obtains the R⁴ coefficient (1804.00949 v), the absence of D⁴R⁴ (1808.10554 v), a D⁶R⁴ match, and a D⁸R⁴ estimate (2412.14094 v). For the (2,0) theory, the A₁ theory saturates the minimal central charge c=25 (1507.05637 v).
- OUR ATTEMPT: From the Fermi gas we got Z(1)=1/4 and Z(2)=1/(16π) to 10⁻¹². The ratio Z/Z_Airy−1 falls from 5×10⁻⁶ at N=1 to 2.5×10⁻¹¹ at N=8; these are exponentially small corrections (`check_abjm_fermigas.py`).
- Missing: every complete definition is AdS or DLCQ. Nothing defines flat 11d M-theory, or any background without a boundary or null circle.

## 2. The three sharpest open questions

**Q1 — The BFSS flat-space limit at nonzero longitudinal transfer.** Take the matrix-model 2→2 amplitude between threshold bound states, A_N(N_i/R, p⊥_i).
- Does the limit N_i,R→∞ with N_i/R fixed exist?
- Does it equal the 11d supergravity four-graviton amplitude, ∝κ₁₁² × (R⁴ polarization structure)/(stu), for generic kinematics with N₁≠N₃ (p₋ exchanged)?
- Do the finite-N mismatches (Dine–Echols–Gray v⁶ at three loops) vanish in that limit?

An explicit large-N computation, or a proof, would settle it. Eikonal checks at q₁₁=0 and the index-based three-point result do not.

**Q2 — Does Lorentzian IKKT have a regulator-free (3+1)-dimensional phase?** Define T_ij(t) from the band-diagonal blocks of the A₀ eigenbasis. Is there an ordering of limits (N→∞, then γ, ξ, m_f→0, or Lorentz gauge fixing as in 2404.14045) in which all of the following hold?
- (a) The model does not flow to the Euclidean saddle.
- (b) Exactly three eigenvalues of T_ij(t) grow at late t.
- (c) The 3×3 spatial blocks are not Pauli-like: the Σ_{i≤3}A_i² block has O(n) comparable eigenvalues, not 2.

A clean "no" is possible. The massless limit of polarized IKKT is already known not to be IKKT (2411.18678).

**Q3 — Does closed SFT have a finite-distance vacuum?** Is there a Lorentz-scalar, translation-invariant solution Ψ* of hyperbolic closed bosonic SFT, Q_BΨ+Σ(1/n!)ℓ_n(Ψ^n)=0, satisfying both of the following?
- S[Ψ*]=0.
- The kinetic operator around Ψ* has empty cohomology at ghost number 2, meaning no closed-string states.

The test: the truncated critical action converges as the level or grade increases, using the Firat–Valdes-Meller integral equation. A "yes" would be the first non-local background change computed inside SFT, which is the missing evidence that SFT carries non-perturbative content.

## 3. Key papers

| arXiv | Authors | Year | Title (short) | Verified? | Why it matters |
|---|---|---|---|---|---|
| hep-th/9610043 | Banks, Fischler, Shenker, Susskind | 1996 | M theory as a matrix model | v | The flat 11d conjecture |
| hep-th/9612115 | Ishibashi, Kawai, Kitazawa, Tsuchiya | 1996 | Large-N reduced model as superstring | v | IIB matrix integral |
| hep-th/9704080 | Susskind | 1997 | Another conjecture about M(atrix) theory | v | Finite N = DLCQ |
| hep-th/9710009 | Seiberg | 1997 | Why is the matrix model correct? | v | DLCQ derivation; T^{p≤5} |
| hep-th/9711037 | Hellerman, Polchinski | 1997 | Compactification in the lightlike limit | v | Zero-mode subtleties |
| hep-th/9706072 | Becker, Becker, Polchinski, Tseytlin | 1997 | Higher order graviton scattering | v | 15/16 and 225/64 (eq. 13) |
| hep-th/9810021 | Dine, Echols, Gray | 1998 | Tree level supergravity and the matrix model | v | Three-loop v⁶ mismatch |
| hep-th/0202021 | Berenstein, Maldacena, Nastase | 2002 | Strings in flat space and pp waves | v | BMN matrix model |
| hep-th/0205185 | Dasgupta, Sheikh-Jabbari, Van Raamsdonk | 2002 | Matrix perturbation theory on a pp-wave | v | Exact fuzzy-sphere spectrum |
| 1606.04951 | Berkowitz, Rinaldi, Hanada, Ishiki, Shimasaki, Vranas | 2016 | Precision lattice test at large N | v | E/N²=7.41T^{14/5} |
| 2302.04416 | Lin | 2023 | Bootstrap bounds on D0-brane QM | v | Rigorous bounds |
| 2312.12592 | Herderschee, Maldacena | 2023 | Three point amplitudes in matrix theory | v | Nonzero-p₋ amplitude |
| 2312.15111 | Herderschee, Maldacena | 2023 | Soft theorems in matrix theory | v | Lorentz invariance, conditionally |
| 1904.05914 | Aoki, Hirasawa, Ito, Nishimura, Tsuchiya | 2019 | Structure of emergent 3d space | v | Pauli-matrix artifact |
| 2307.01681 | Hirasawa et al. | 2023 | Lorentzian IIB with novel regularization | v | Euclidean equivalence; γ term |
| 2604.19836 | Anagnostopoulos, Azuma, Hirasawa, Nishimura, Papadoudis, Tsuchiya | 2026 | (3+1)D from deformed Lorentzian IIB | v | Smooth (3+1)D with deformations |
| 2409.18706 | Hartnoll, Liu | 2024 | The polarised IKKT matrix model | v | Holography for a matrix integral |
| 2411.18678 | Komatsu, Martina, Penedones, Vuignier, Zhao | 2024 | Einstein gravity from a matrix integral II | v | IKKT ≠ massless limit |
| 2303.08012 | Steinacker | 2023 | One-loop emergent gravity in IKKT | v | 4d Einstein–Hilbert mechanism |
| hep-th/9206084 | Zwiebach | 1992 | Closed SFT: quantum action and BV | v | Foundation of SFT |
| hep-th/9311009 | Sen, Zwiebach | 1993 | Quantum background independence of CSFT | v | Local only |
| 1909.00033 | Costello, Zwiebach | 2019 | Hyperbolic string vertices | v | Exact vertices |
| 2409.02982 | Firat, Valdes-Meller | 2024 | Topological recursion for hyperbolic SFT | v | Cubic data fixes everything |
| hep-th/0506077 | Yang, Zwiebach | 2005 | A closed string tachyon vacuum? | v | Q3 baseline numbers |
| 2104.11109 | Sen | 2021 | Normalization of IIB D-instanton amplitudes | v | Non-perturbative data from SFT |
| 1110.4066 | Mariño, Putrov | 2011 | ABJM theory as a Fermi gas | v | Exact M-theory benchmark |

## 4. Computable handles (single researcher, hours, numpy/sympy/scipy)

**H1. D0-brane eikonal series. Partly done.**
- Compute: the one-loop proper-time mode sum (done), then the two-loop SU(2) v⁶/r¹⁴ coefficient, using the BBPT/Becker–Becker propagators in the background X¹=vt, X²=b.
- Compare to: 225/32 (BBPT eq. 14).
- New: the three-loop v⁸/r²¹ coefficient checked against our DLCQ-supergravity target 16875/1024·N₁³N₂/(R⁷M²⁷). We found no published three-loop two-body check this session. The protected/unprotected boundary is exactly what Q1 asks about.

**H2. BMN spectra. Base done.**
- Extend to reducible vacua (DSV Table 2), then to the second-order O((R/μN)³) energy shift of a non-BPS state, using the same Hessian code.
- Compare to: DSV §6 effective coupling λ=g²K.
- New: the one-loop interaction between two separated fuzzy spheres, a membrane analogue of v⁴/r⁷ that can be compared with probe supergravity in the pp-wave.

**H3. Hyperbolic SFT couplings. Cubic t³ done.**
- Use Firat's explicit w(z) (eq. 1.3, including the O(z²) term) to compute the tachyon–ghost-dilaton coupling t·d². The ghost dilaton is not primary, so the coupling depends on the second coefficient of the map, not just on r.
- Compare to: the minimal-area value −27/32 (Yang–Zwiebach eq. 1.3) in the limit λ→∞.
- New: these couplings as a function of λ, a first step toward testing whether the Yang–Zwiebach critical point survives in hyperbolic SFT (Q3).

**H4. Small-N IKKT.**
- (a) Euclidean bosonic Metropolis Monte Carlo, D=10, N=4–16. Exact check: ⟨S⟩=D(N²−1)/4, which follows from scaling A→(1+ε)A (our derivation). Measure the eigenvalues of T_μν=Tr A_μA_ν/N and confirm there is no spontaneous symmetry breaking (hep-th/9811220).
- (b) Lorentzian bosonic classical solutions with a Lorentz-invariant mass term, found by gradient descent from random starts (method of 1911.08132 v), at N=16–64. Measure a Pauli index, (λ₁+λ₂)/Σλ of the Σ_{i≤3}A_i(t)² blocks.
- New: a classical map of Pauli versus smooth structure as a function of γ and the 2604.19836 anisotropy ξ. It tests whether their deformation already biases the result toward smooth (3+1)D at tree level (Q2).

**H5. ABJM exact partition function. Base done.**
- Fit δ(N)=Z/Z_Airy−1 for k=1–4 and N up to 12 to the form Σ c_a μ^p e^{−aμ}, with μ the Airy saddle point.
- Compare to: the worldsheet and membrane instanton actions in the Fermi-gas literature (exponents not verified this session).
- New: a calibrated resurgence pipeline, to be reused on BFSS bootstrap data and on SFT g_s series.
