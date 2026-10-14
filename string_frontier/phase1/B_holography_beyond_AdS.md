# Front B — Holography beyond AdS (Phase 1 frontier map)

Conventions: **[V]** means the arXiv abstract page was fetched in this session and the title and authors match. Equation numbers are given only where I read them in fetched text. Checks I ran myself are in `phase1/B_checks/` (scripts plus `.out` files).
ID corrections found while checking: Maldacena's wavefunction paper is **astro-ph/0210603** (not hep-th/0210603). The Costello–Paquette–Sharma PRL is **2208.14233**. (2208.10567 is an unrelated neutrino paper.)

## 1. State of the art

**dS/CFT (future boundary)**
- ESTABLISHED (computation), CONJECTURE (duality). The asymptotic symmetries of dS3 give a Euclidean Virasoro×Virasoro algebra with c = 3ℓ/2G (hep-th/0106113 [V], eq. 4.2). Scalar weights are h± = 1 ± √(1−m²ℓ²) (eq. 5.3), which are complex for m²ℓ² > 1, so the dual is non-unitary.
- CONJECTURE. Asymptotically dS space has no precisely measurable observables, only "meta-observables" (hep-th/0106109 [V]). dS/CFT computes the late-time wavefunction, Ψ ≈ Z_CFT (astro-ph/0210603 [V]) (eqs. not checked).
- CONJECTURE (best-checked dS4 example). Vasiliev dS4 gravity is dual to the Sp(N) anticommuting-scalar CFT3. It is obtained from O(N)/AdS4 by N → −N, since N ∼ 1/(ΛG_N) (1108.5735 [V], introduction text). The boundary theory is non-unitary, and the bulk inner product is not fixed by the proposal.
- CONJECTURE. An SU(2)_k WZW model with k → −2 reproduces the classical S³ partition function (2110.03197 [V]).
- CONJECTURE with strong internal checks (2024–25). Canonical quantization of dS3 gives states that transform like CFT correlators with c ∈ 13+iℝ, where Im c = 3ℓ/2G. The inner product (1.2) is a string path integral, with g_s⁻² ∼ Z_grav(S³) (1.3). Integrated cosmological correlators equal amplitudes of the complex Liouville string (2409.17246 [V]), which is dual to a two-matrix integral. The claim S_dS = log|Z_S³| = log N_eff² (eq. 1.1) uses N_eff = ∫₂^{E₀} e^{S₀}ρ₀ with E₀ = 2cos(πb⁻²) (4.32–4.33). By the authors' own statement, the match "works regardless of the actual value of Z_S³" (text near 1.8 and 4.36). The choice of wavefunction is "the most speculative part" (2501.01486 [V]).

**Sphere partition function and dS entropy**
- ESTABLISHED. There is an exact one-loop formula for any effective field theory, written as a bulk character minus an edge character (2009.12464 [V], eqs. 1.4–1.11). Results:
  - dS3 gravity: S⁽¹⁾ = −3 log S⁽⁰⁾ + 5 log 2π.
  - dS4 gravity: S⁽¹⁾ = −5 log S⁽⁰⁾ − (571/45) log(ℓ/L) + … (eq. 1.12).
  - The pair (a, b) in S = S₀ − a log S₀ + b + … (eq. 1.3) is proposed as invariant data that any microscopic model must reproduce.
  - **Checked here:** from (1.4), (1.9), (1.11) and (1.13) I recomputed α₄ = −571/45 for the graviton (bulk −331/45, edge −16/3) and −31/45 for the photon (bulk −16/45, edge −1/3) (`check_adls_log.py`).
  - The dS3 result agrees with log Z_S³ = S_GH − 3 log S_GH + 5 log 2π ± 5πi/2 + … quoted in 2501.01486 (eq. 1.8).
- ESTABLISHED computation, CONJECTURE interpretation. The one-loop phase is Z_grav(S^D) = i^{D+2} × (positive) (2412.14014 [V], eq. 1).
  - Adding an observer gives i^{D+2}(−i)^{D−1} = i³. The latest version of 2412.14014 says an overall minus sign remains.
  - 2511.01400 [V] (Chen–Stanford–Tang–Yang) uses charged-black-hole observers. Its table (3.8) gives phase −i for lukewarm black holes, 1 for Nariai with Q < 1/4, and −i for Nariai with Q > 1/4. With a revised prescription, Z_Count is positive (4.7–4.8), but this "relies on … (−i)^{−4} = 1", which the authors say is not fully understood.
- ESTABLISHED. dS JT gravity gives no-boundary wavefunctions to all orders, plus topology sums (1904.01911 [V]).

**Observables and algebras**
- ESTABLISHED (semiclassical). The static-patch algebra, with operators dressed to an observer's worldline, is Type II₁. Its maximum-entropy state is empty dS. The entropy equals S_gen = A/4G + S_out up to a constant (2206.10780 [V], eq. 1.1), with A_hor/4G = A_dS/4G − β_dS E and p(E) = e^{−β_dS E} (eqs. 1.2–1.3).
- CONJECTURE. Proposed extensions:
  - a background-independent worldline algebra, with the Hartle–Hawking state as maximum entropy (2308.03663 [V]);
  - observer and clock algebras in cosmology (2406.02116 [V]);
  - the closed-universe Hilbert space is one-dimensional, and observers get an effective description of dimension ≈ e^{S_Ob} (2501.02359 [V]);
  - the need for observers comes from spontaneous breaking of time reversal (2601.01666 [V]).
- ESTABLISHED (large-volume Wheeler–DeWitt). Solutions with Λ > 0 are a universal phase times functionals that obey CFT Ward identities. The Euclidean vacuum is only one choice (2303.16315 [V]).

**Static-patch holography and DSSYK**
- ESTABLISHED. The chord Hamiltonian is H|n⟩ = |n+1⟩ + [n]_q|n−1⟩ (2402.00635 [V], eq. 2; 1811.02584 [V]).
  - **Checked here:** the spectrum edge is 2/√(1−q). The first three moments match chord counting (1, 2+q, 5+6q+3q²+q³) and the q-Hermite measure (`check_dssyk_transfer.py`).
  - This means E = 2cosθ/√(1−q). The rendered eq. (3) of 2402.00635 reads cosθ/√(1−q), a normalization or typo issue. 2310.16994 eq. (23) is consistent with my result.
- ESTABLISHED (correlator match), CONJECTURE (duality). In 2310.16994 [V], two infinite-temperature DSSYK models are coupled with H_L = H_R.
  - The exact two-point function is given in (42)–(43).
  - As λ → 0 it becomes G = 2 sinh(μτ)/(π sinh τ) with μ = 1−2Δ (eq. 45). This is the antipodal Green function in dS3 with m² = 4Δ(1−Δ).
  - The dictionary is T = J/2π (eq. 8) and R_dS/G_N = 8π/λ = 4πN/p² (eq. 11).
  - Remark: S_GH = πR/2G = 4π²/λ = 2π²N/p², which is not the full DSSYK entropy (≈ (N/2) log 2). Only entropy *differences* are matched (eqs. 9–10).
- CONJECTURE. Quantizing the phase space of Schwarzschild–dS3 and setting H = L_A/√(1−q) with L_A = 2 sin πα (eqs. 4–5, 25) reproduces the chord rule through skein relations (2402.00635). This requires gauging CPT, which breaks SL(2,ℂ) to SU(1,1).
- CONJECTURE. Competing dictionaries:
  - Sine-dilaton gravity: chord number = geodesic length; the "fake temperature" puzzle; discrete spacetime (2404.03535 [V]). The one-loop log correction matches "up to potential ordering ambiguities" (2411.15957 [V]).
  - Upper-edge scaling: scaling the ETH matrix model at the upper edge E = +E₀ gives the dS-JT disk (eq. 65) as an *unstable* saddle; the sine-dilaton solution at θ = π is −AdS₂ = dS₂ (eqs. 88–89) (2505.08116 [V]).
  - Earlier proposals and the temperature puzzle: 2109.14104, 2209.09997, 2401.08555; observer entropy matched to the DSSYK density: 2502.03869 (all [V]).
  - 2026: the dS3 static-patch density comes from an SL(2,ℤ) sum of Kerr-lens spacetimes; dS3 is dual to two copies of the complex Liouville string; in the simplest case this is dual to two copies of the DSSYK GΣ theory living on the observer's worldline (2606.26241 [V]).

**Cosmological bootstrap**
- ESTABLISHED. dS conformal Ward identities make the exchange 4-point function satisfy (Δ_u + μ² + ¼)F = Ĉ(u,v), with Δ_u = u²(1−u²)∂²_u − 2u³∂_u (1811.00024 [V], eqs. 2.30 and 3.11). Explicit solution: eq. 3.22. The particle-production term ∝ e^{−πμ} is fixed by eq. 3.21. Review: 2203.08121 [V].

**Flat space: celestial and Carrollian holography**
- ESTABLISHED. The infrared triangle links soft theorems, asymptotic symmetries and memory (1703.05448 [V]). Conformal primary wavefunctions on the principal series Δ = 1+iλ are complete, and massless amplitudes become SL(2,ℂ) correlators by Mellin transform (1705.01027 [V], 1701.00049 [V]).
- ESTABLISHED (tree level). Positive-helicity soft gravitons form the wedge of w₁₊∞: [w^p_m, w^q_n] = [m(q−1) − n(p−1)] w^{p+q−2}_{m+n} with 1−p ≤ m ≤ p−1 (2105.14346 [V], eqs. 2.5–2.9), building on 2103.03961 [V].
  - **Checked here:** 1864 structure constants of (2.5) rescaled by (2.6) reproduce (2.7), and 9261 Jacobi triples close (`check_w1inf.py`).
- ESTABLISHED. Celestial chiral algebras generally fail to be associative once loop effects or higher-derivative couplings are included. Associativity is guaranteed when the 4D theory lifts to a local holomorphic, anomaly-free theory on twistor space (2208.14233 [V], introduction; 2204.05301 [V]).
- CONJECTURE with checks. Top-down examples:
  - SO(8) WZW₄ plus scalar-flat Kähler gravity on Burns space is dual to the large-N chiral algebra on D1 branes, checked with 2- and 3-point functions (2208.14233 [V], 2201.02595 [V]);
  - all loop amplitudes on Burns space are chiral-algebra correlators (2306.00940 [V]);
  - self-dual QCD on self-dual backgrounds, matched to one loop, giving two-loop all-plus amplitudes (2412.02680 [V]).
  - **All of these are self-dual sectors.**
- CONJECTURE. A 2D CFT (Liouville with central charge c, plus a level-1 rank-N Kac–Moody algebra, plus a weight −3/2 fermion) generates tree-level MHV leaf amplitudes at large N and large c (2403.18896 [V]).
- ESTABLISHED (kinematics/perturbative), CONJECTURE (that the dual exists). Carrollian holography:
  - a sourced conformal Carrollian theory on I reproduces the celestial and BMS flux-balance Ward identities (2202.04702 [V], 2212.12553 [V]);
  - Carrollian amplitudes (2312.10138, 2202.08438 [V]); the flat limit of Lorentzian AdS correlators gives Carrollian amplitudes with a bulk-point singularity (2406.19343 [V]); lectures 2511.10162 [V]; Mellin flat-space limit 1011.1485 [V].
- ESTABLISHED. Celestial string amplitudes:
  - Type-I and heterotic 4-point Mellin integrals are finite because of exponential UV softness. α′ enters only as an overall (α′)^β, and the expansion parameter is 1/r (small angle) (1806.05688 [V], eqs. 3.10, 3.15–3.16, 3.26, 3.30–3.31).
  - Exponential softness implies no poles at Re β_APRS > 0. Poles at β = −2n have residues equal to low-energy coefficients (2012.04208 [V], toy eq. 24).
- OUR ATTEMPT (preliminary; I ran the numerics). I integrated eq. 3.15 of 1806.05688 directly at Re β < 0, contour-independent to ≥8 digits, and compared it with the mass-level representation (3.30)/(3.31) as transcribed.
  - The residue sum reproduces **exactly** the discontinuity across the s-channel poles: (above − below) = −iπS to 10 digits.
  - It does **not** reproduce the full integral. Example at r = 2, β = −1.3: direct = −0.035898 + 0.198276i, while (3.31) gives +0.144056 ∓ 0.198276i.
  - Also, the residue sum is nonzero at β = 1 (r = 2: S = −0.96242), although F_I has no O(w) term. So (3.31) cannot be the analytic continuation as written.
  - Possible causes: a transcription or branch issue, or a missing contribution from the contour at infinity or the u-channel poles (whose residues grow exponentially). Needs resolution (`check_celestial_string_direct*.py`).

## 2. The three sharpest open questions

1. **Is dS entropy a state count in D = 4?** Find a microscopic system (ideally from string theory) whose state count, or Type II₁ trace normalization with an observer, is real and positive and reproduces log Z(S⁴) through one loop. That means matching −5 log S⁽⁰⁾, the non-local constant terms of 2009.12464 (eq. 1.12), and the phase prescription of 2412.14014 / 2511.01400. *Answered if* a model reproduces these numbers without fitting them, or a no-go shows that no positive count can exist. (In dS3, the 2501.01486 count reproduces S_dS only structurally, independent of Z_S³.)
2. **Is DSSYK a hologram of dS, and of which dS?** Does the exact doubled-DSSYK two-point function (2310.16994 eqs. 42–43) agree at O(λ) with the one-loop-corrected dS3 propagator, using G/R = λ/8π? And can one dictionary hold all three pictures together: dS3 from the doubled model, dS₂-JT at the upper spectral edge as an unstable saddle (2505.08116), and sine-dilaton gravity (2404.03535)? This includes which temperature is the Gibbons–Hawking temperature. *Answered if* an explicit O(λ) bulk computation matches or fails to match, and if the entropy-difference identification (eqs. 9–11) survives one loop.
3. **Is there flat-space holography beyond self-dual sectors?** Is there a top-down dual, an independently defined and associative chiral or Carrollian operator algebra, for a 4D asymptotically flat theory containing non-self-dual Einstein gravity? It should reproduce at least tree-level MHV and NMHV 4-graviton amplitudes and the one-loop collinear splitting. *Answered if* such a construction passes these checks, or if an associativity or anomaly obstruction is proven (cf. 2208.14233, 2204.05301).

## 3. Key papers (all IDs verified via abstract fetch in this session)

| arXiv ID | Authors | Year | Title | Verified? | Why it matters |
|---|---|---|---|---|---|
| hep-th/0106113 | Strominger | 2001 | The dS/CFT Correspondence | yes | c = 3ℓ/2G; complex weights / non-unitarity |
| hep-th/0106109 | Witten | 2001 | Quantum Gravity in de Sitter Space | yes | Observables problem; finite Hilbert space |
| astro-ph/0210603 | Maldacena | 2002 | Non-Gaussian features of primordial fluctuations… | yes | dS/CFT as the wavefunction of the universe |
| 1108.5735 | Anninos, Hartman, Strominger | 2011 | Higher Spin Realization of the dS/CFT Correspondence | yes | Only explicit dS4 dual (Sp(N)) |
| 2009.12464 | Anninos, Denef, Law, Sun | 2020 | Quantum de Sitter horizon entropy from quasicanonical bulk, edge, sphere and topological string partition functions | yes | Exact one-loop S(S^D); invariant (a, b) test data |
| 1904.01911 | Maldacena, Turiaci, Yang | 2019 | Two dimensional Nearly de Sitter gravity | yes | dS JT gravity, no-boundary wavefunction |
| 2206.10780 | Chandrasekaran, Longo, Penington, Witten | 2022 | An Algebra of Observables for de Sitter Space | yes | Type II₁ static patch; S_gen |
| 2412.14014 | Maldacena | 2024 | Real observers solving imaginary problems | yes | Observer cancels most of the i^{D+2} phase |
| 2511.01400 | Chen, Stanford, Tang, Yang | 2025 | On the phase of the de Sitter density of states | yes | Positive density of states with black-hole observer |
| 2501.01486 | Collier, Eberhardt, Mühlmann | 2025 | A microscopic realization of dS3 | yes | dS3 ↔ matrix integral; GH entropy as an eigenvalue count |
| 1811.02584 | Berkooz, Isachenkov, Narovlansky, Torrents | 2018 | Towards a full solution of the large N double-scaled SYK model | yes | Chord-diagram solution of DSSYK |
| 2310.16994 | Narovlansky, Verlinde | 2023 | Double-scaled SYK and de Sitter Holography | yes | Exact dS3 Green function from doubled DSSYK |
| 2402.00635 | Verlinde | 2024 | Double-scaled SYK, Chords and de Sitter Gravity | yes | Chord rules from SdS3 Wilson loops |
| 2404.03535 | Blommaert, Mertens, Papalini | 2024 | The dilaton gravity hologram of double-scaled SYK | yes | Sine-dilaton bulk; discrete spacetime |
| 2505.08116 | Okuyama | 2025 | de Sitter JT gravity from double-scaled SYK | yes | dS-JT as an unstable upper-edge saddle |
| 2606.26241 | Blommaert, Tietto, Verlinde | 2026 | An observer's quantization of 3d de Sitter | yes | Static-patch density; CLS⊗CLS; worldline hologram |
| 1811.00024 | Arkani-Hamed, Baumann, Lee, Pimentel | 2018 | The Cosmological Bootstrap… | yes | Boundary-only dS correlators |
| 2105.14346 | Strominger | 2021 | w(1+infinity) and the Celestial Sphere | yes | Soft-graviton symmetry = wedge of w₁₊∞ |
| 2208.14233 | Costello, Paquette, Sharma | 2022 | Top-down holography in an asymptotically flat spacetime | yes | First top-down flat-space dual (Burns) |
| 2412.02680 | Bittleston, Costello, Zeng | 2024 | Self-Dual Gauge Theory from the Top Down | yes | Family of top-down self-dual duals, loop checks |
| 2403.18896 | Melton, Sharma, Strominger, Wang | 2024 | A Celestial Dual for MHV Amplitudes | yes | Explicit 2D CFT for MHV |
| 2406.19343 | Alday, Nocchi, Ruzziconi, Yelleshpur Srikant | 2024 | Carrollian Amplitudes from Holographic Correlators | yes | Flat limit of AdS = Carrollian limit |
| 2212.12553 | Donnay, Fiorucci, Herfray, Ruzziconi | 2022 | Bridging Carrollian and Celestial Holography | yes | Sourced Carrollian CFT ↔ BMS fluxes |
| 1806.05688 | Stieberger, Taylor | 2018 | Strings on Celestial Sphere | yes | Celestial superstring amplitudes |
| 2012.04208 | Arkani-Hamed, Pate, Raclariu, Strominger | 2020 | Celestial Amplitudes from UV to IR | yes | β-plane analyticity from UV softness |

## 4. Computable handles (Python, a few hours each)

**H1. Doubled-DSSYK → dS3 beyond leading order** (Q2)
- *Done:* the chord transfer matrix spectrum and moments (`check_dssyk_transfer.py`).
- *Compute:* implement 2310.16994 eqs. (40)–(43), taking the ϑ₁ conventions from its Appendix A. Evaluate G_Δ(E) and G_Δ(τ) at λ = 0.05–0.5.
  - Confirm the limit (44)–(45).
  - Confirm that the poles of (44) sit at the dS3 ℓ = 0 quasinormal frequencies ω = ±i(2Δ+2k) and ±i(2−2Δ+2k).
  - Fit the O(λ) correction δG(τ).
- *Compare:* the bulk one-loop matter-propagator correction of sine-dilaton gravity (2411.15957, AdS-side method), continued to θ ≈ π/2. Also compare the upper-edge density with Okuyama's dS-JT disk (2505.08116 eq. 65).
- *New:* the first O(G/R) test of the Narovlansky–Verlinde dictionary; a mismatch would falsify R/G = 8π/λ.

**H2. One-loop sphere data as a filter for microscopic models** (Q1)
- *Done:* α₄ = −571/45 (graviton) and −31/45 (photon) from the 2009.12464 character integral.
- *Compute, step 1:* run (1.4) for dS3 and dS4 with massive spin-s fields.
- *Compute, step 2:* verify numerically 2501.01486 (4.32)–(4.36), i.e. that (T̂₁,₁)² ∝ C_S² with T̂₁,₁ = 8b² sin(πb²) sin(πb⁻²)/(1−b⁴) and C_S² = 32π⁴[sin(πb²) sin(πb⁻²)/(b²−b⁻²)]². Check whether −3 log S_GH + 5 log 2π is an output of the count or only an input.
- *Compute, step 3:* compute log[ρ_DSSYK(θ)/ρ(π/2)] exactly in q and expand at small λ and fixed v = 1−2θ/π. Compare with −S_GH·α (2310.16994 eqs. 9–11), and extract the log λ and O(1) terms.
- *New:* a yes/no on whether either dS3 microscopic proposal encodes the invariant one-loop data (a, b) = (3, 5 log 2π), rather than having them absorbed into definitions.

**H3. Analytic structure of celestial string amplitudes** (flat space; Q3 groundwork)
- *Done (OUR ATTEMPT):* the direct iε evaluation of 1806.05688 eq. (3.15) disagrees with (3.30)/(3.31) except for the discontinuity (§1).
- *Next:* settle the representation, either by continuing from r < 0 where the formula was derived, or by tracking the u-channel poles and the contributions at infinity. Then repeat for Virasoro–Shapiro (type II, 4 gravitons).
- *Map:* all poles in the β-plane, and check the 2012.04208 claims: residues at β_APRS = −2n should equal the α′-expansion coefficients (ζ(3), ζ(5), …), and there should be no poles at Re β_APRS > 0.
- *New:* correct closed forms for the celestial type-I and Virasoro–Shapiro amplitudes, and how real-axis resonances interact with the "no UV poles" claim.

**H4. Do α′ corrections deform w₁₊∞?** (Q3)
- *Done:* (2.5) → (2.7) and Jacobi closure (`check_w1inf.py`).
- *Compute:* take collinear (OPE) limits of the tree-level 4- and 5-graviton amplitudes with α′³R⁴ (from the expansion of Virasoro–Shapiro). Extract the corrected conformally-soft current commutators, then test closure and Jacobi with the same code.
- *Compare:* (2.7). Check the literature on deformed w₁₊∞ before claiming novelty; no IDs were verified here.
- *New:* whether string theory preserves, deforms (e.g. to W₁₊∞), or breaks the celestial symmetry.

**H5. Phase of the dS density of states with observers** (Q1)
- *Compute:* reproduce the 2511.01400 table (3.8) by counting negative modes of the quadratic fluctuation operator of the 2D dilaton reduction (conformal gauge, finite differences). Then extend to a dS3 conical-defect observer and a rotating observer (cf. the Kerr-lens saddles in 2606.26241).
- *Compare:* i^{D+2} (2412.14014 eq. 1); and lukewarm −i, Nariai Q < 1/4: 1, Q > 1/4: −i.
- *New:* whether Z_Count > 0 is robust across observer types, i.e. whether the static patch has a positive density of states at all.
