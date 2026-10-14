# R2-T4: Testing the CEM dS₃/matrix-integral proposal against exact S³ data

**Target.** Collier, Eberhardt and Mühlmann (CEM), arXiv:2501.01486 [abstract fetched; published as SciPost Phys. 18, 131].

**Equation numbers.** Each CEM equation is cited twice: arXiv-HTML numbering / SciPost numbering. Our notes' "4.32–4.36" is the SciPost numbering. The two systems map as follows:

| arXiv-HTML | SciPost |
|---|---|
| (73)–(78) | (4.32)–(4.37) |
| (44)–(46) | (4.3)–(4.5) |
| (83)–(87) | (5.5)–(5.9) |

**Files.** Scripts and outputs are in this folder (`r2t4_A…D_*.py` / `.out`). Fetched texts are in `fetched/`.

## 0. Verdict (short)

1. **The microstate count passes only by construction.**
   - It reproduces S_GH, −3 log S_GH and 5 log 2π only because these are put in through the string coupling.
   - It does not reproduce the phase; the phase is removed by hand with |Z|.
   - Our script feeds the count deliberately wrong sphere partition functions, and it "reproduces" every one of them exactly.
   - CEM say this themselves ("works regardless of the actual value of Z_grav^{S³}").
   - The only content that depends on b is an identity internal to the complex Liouville string (CLS) / matrix model.
2. **CEM's separate exact proposal (85)/(5.7)** is a genuine but conditional prediction.
   - **Log coefficient −3:** an output, but only because of an ad hoc h-basis Jacobian. In the Liouville-momentum (p) basis it becomes −2.
   - **5 log 2π:** an input (the prefactor was fitted).
   - **Phase ±i:** an output, again tied to the h-basis.
   - **No test available yet:** the first coefficient that could discriminate is c₁, at order S₀⁻² (three-loop). There, (85) disagrees with both Chern–Simons all-loop answers:
     - c₁(CEM) = 4π²(π²/3−2) ≈ 50.92;
     - c₁(CCFZ, δ̂=0) = 4π⁴/3 ≈ 129.88;
     - c₁(ADLS, l=0) = 4π²(π²/3−6) ≈ −106.99.

   No independent three-loop computation exists, so this cannot yet be decided.

## 1. What 2501.01486 computes: inputs vs outputs

**ESTABLISHED (2501.01486, read in fetched text).**
- **Duality claim.** Integrated cosmological correlators (sums of norms of the Liouville wavefunction) = CLS amplitudes = resolvents of a double-scaled two-matrix integral, via (43)–(47)/(4.2)–(4.8). The genus counting parameter is e^{S₀}.
- **String coupling from Z(S³).** Equations (25)–(27)/(2.17) fix the string coupling by Z_grav^{S³} ∼ g_s⁻², "up to an overall b-independent constant". The worldsheet relation is g_s⁻² ∼ e^{2S₀} C_{S²} (44)/(4.3), with C_{S²} = 32π⁴[sin(πb²) sin(πb⁻²)/(b²−b⁻²)]² (45)/(4.4).
- **Z(S³) is taken as input.** Text above (28): "We take in the following Z_grav^{S³} = e^{S_dS} as an input."
- **The count.**
  - N_eff = ∫₂^{E₀} e^{S₀} ρ₀ dE, with E₀ = 2cos(πb⁻²) (73)/(4.32) and ρ₀ = (2/π) sinh(−iπb²) sin(−ib² arccosh(E/2)) (135).
  - S_micro = 2 log N_eff = 2 log(T_{1,1}/2πi) (74)/(4.33), where T_{1,1} = e^{S₀}·8b² sin(πb²) sin(πb⁻²)/(1−b⁴) (75)/(4.34).
  - The match condition reduces to T̂²_{1,1} ∼ C_{S²} (77)/(4.36).
  - The O(1) constant is then fixed by g_s⁻¹ = T_{1,1} = 2πi|Z|^{1/2} (78)/(4.37).
- **Their own caveat** (text after (77)): the matching "doesn't depend on the specific value of Z_grav^{S³}". The phase is handled by "the pragmatic solution of including an absolute value" (text after (8)/(1.8)).
- **The one-loop data are quoted, not derived.** CEM eq. (8)/(1.8) quotes log Z = S_GH − 3 log S_GH + 5 log 2π ± 5πi/2 with S_GH = πℓ/2G. Eq. (28)/(2.20) gives the same in b variables. The trace-mode phase comes from ζ_{S³}(0) = −5, (66)–(70).
- **The separate conjecture**, eqs. (83)–(87)/(5.5)–(5.9):
  - Z = S₁₁² = −16b² sin²(πb²) sin²(πb⁻²)/(1−b⁴)², using the S-kernel in the h basis. CEM write: "we do not have a better justification for it other than…"
  - The prefactor is "fixed … such that it will agree with the one-loop result".
  - This gives log Z = S_GH − 3 log S_GH + 5 log 2π − πi/2, with e^{S₀} = π/|b|.

**OUR ATTEMPT (`r2t4_A_counting_inputs.py`).**
- **(77) is an exact identity.** T̂²_{1,1}/C_{S²} = 2/π⁴, independent of b (sympy).
- **(74) checked by quadrature.** Integrating (135) at b² = iβ, β ∈ [0.05, 5], gives T̂/(2πi) to 40 digits, and E₀ is the first sign change of ρ₀.
- **Input test.** We fix e^{S₀} by (78) and count eigenvalues numerically for five inputs: ADLS one-loop, CEM (85), ADLS (8.18), and two deliberately wrong inputs ("a = −7, phase 0.3" and "S₀/2 + 11 log S₀").
  - For all five, log N_eff² = log|Z_in| holds to working precision (0 at 40 digits).
  - So (i) S_GH, (ii) −3 log S_GH and (iii) 5 log 2π are **inputs**.
  - (iv) the phase is **discarded**, not reproduced.
- **Where the large entropy comes from.** With (85), e^{S₀} = π/|b| ∝ S₀^{1/2}, so the genus-counting parameter grows only like √S_GH. The e^{S_GH/2} in N_eff comes from the width of the spectral window E₀ = 2cosh(π/β).

| quantity | role in the counting argument (73)–(78) |
|---|---|
| S_GH (leading) | input, via e^{S₀} fixed by Z (44)+(78) |
| −3 log S_GH | input (same route) |
| 5 log 2π | input; the O(1) constant is fixed by hand in (78) |
| phase ±5πi/2 | discarded (absolute value in (1)/(1.1)) |
| b-independence of T̂²/C_{S²} | the only output; internal to CLS / matrix model |

## 2. Independent data on Z(S³)

**ESTABLISHED (2009.12464, Anninos–Denef–Law–Sun (ADLS), read in fetched text).**
- **d = 2 has no particle contribution.** For massless spin s ≥ 2 in d = 2, χ_bulk = χ_edge = 0 (H.1).
- **One-loop formula.** Z⁽¹⁾_PI = i^{−P} γ^{dim G}/vol(G)_c · Z_char, with Z_char,fin = (2π)^{dim G} (H.5), (C.37). Here P₂ = 5 (H.4) and Z_fin = i⁻⁵(2πγ)⁶/(2π)⁴ (H.7).
- **Result.** log Z⁽¹⁾ = −3 log(2πℓ/4G) + 5 log 2π (1.23), (8.17).
- **Invariant data.** The pair (a, b) in S = S₀ − a log S₀ + b + Σ c_n S₀^{−2n} (1.3) is invariant once S₀ is chosen so that no odd powers of 1/S₀ appear; odd powers are local counterterms. ADLS state (a, b) = (3, 5 log 2π) is "macroscopic data" any microscopic model must reproduce (text near (1.3)).
- **Phase.** i^{−5} in the metric formulation (5.19), (H.7). In the Chern–Simons formulation the phase depends on the framing: canonical framing gives a real Z, and odd framing r gives (−i)^r at weak coupling (H.28)–(H.30).
- **All-loop CS answer.** (H.24) gives |√(2/K) sin(π/K) e^{πκ}|², with K = 2 + l + iκ and κ = ℓ/4G. Pure gravity with l = 0 and SO(4) normalization gives (8.18): S = S₀ + log|√(4/(2+iS₀/2π)) sin(π/(2+iS₀/2π))|². Footnote 17: the choice of l "does not affect the 1-loop based conclusions … but does affect the c_n."

**ESTABLISHED (2302.12281, Castro–Coman–Fliss–Zukowski (CCFZ), [abstract verified]).**
- The renormalized levels are k + 2 (2.9) and r = δ̂ ± is with s = ℓ/4G (2.83).
- The exact result is (2.82): Z = e^{2πs} (i e^{…}) (2/√(r_L r_R)) sin(π/r_L) sin(π/r_R).
- For δ̂ = 0, eq. (2.84) gives log Z = (π/2)(ℓ/G) + 3 log(4G/ℓ) + log(2π²i) + (16π²/3)(G/ℓ)² + …. The two-loop term vanishes, and the i comes from a Polchinski-type contour rotation.

**ESTABLISHED (other phase results).**
- Z(S^D) ∝ i^{D+2} at one loop; with an observer this is mostly cancelled, but an overall minus sign remains (2412.14014 v3 [abstract verified]).
- A positive density of states is found with a revised prescription (2511.01400 [abstract verified]).

**OUR ATTEMPT (`r2t4_B_adls_oneloop.py`).** We re-derived every ADLS ingredient independently.
- **Character finite part.** The finite part of ∫^× dt/2t (1+q)/(1−q) under prescription (5.12) is −0.918952 at ε = L⁻¹ = 10⁻⁶. This converges to −½ log 2π = ζ′(0), so log Z_char = 6 log 2π.
- **Group volume.** vol SO(4) = vol(S¹)vol(S²)vol(S³) = (2π)⁴.
- **γ factor.** γ² = 2π/S₀.
- **Total.** 5 log 2π = (+6 − 4 + 3) log 2π.
- **Phase.** P₂ = 5. Independently, by our own analytic continuation (checked against the direct sum at s = 2.5), ζ_{S³}(0) = −5.000 for the Wick-rotated trace modes. So the one-loop phase is e^{∓5πi/2} = ∓i.
- **Result.** log Z⁽¹⁾ + 3 log S₀ = 9.1893853320 − 1.5707963i at S₀ = 10, 10³ and 10⁶, i.e. 5 log 2π − iπ/2. This confirms the B-map entry and CEM (8).

## 3. Comparison of exact candidates (invariant ADLS form)

**OUR ATTEMPT (`r2t4_C_exact_candidates.py`).**

**Conventions.**
- CEM: b² = iβ, c = 13 + 6i(β − 1/β), and S_GH = (π/3) Im c = 2π(1/β − β).
- The unique S₀ with no odd powers is S₀ = 2π/β; this is proved by the exact series, in which d₁ = d₃ = d₅ = 0.
- ADLS: S₀ = 2πκ. CCFZ: S₀ = 2πs, multiplied by 2 for SO(4) vs SU(2)×SU(2) (ADLS text below (8.17)).

**Method.**
- Exact sympy series in 1/S₀.
- Independent numerical fits of log|Z| − S₀ = −a log S₀ + b + Σ_k d_k S₀^{−k}, using 60-digit evaluation including e^{−S₀} terms. The error bar is the spread over windows [40, 4000], [80, 8000], [160, 16000] and 6, 8 or 10 fit terms.

| candidate | a | b | odd d₁, d₃ | c₁ = d₂ (exact) | fit d₂ | phase arg Z |
|---|---|---|---|---|---|---|
| ADLS one-loop (1.23) | 3 | 5 log 2π | n/a | n/a | n/a | ∓π/2 (metric) |
| CEM (85), h-basis | 3.0 ± 10⁻¹² | 9.18938533205 ± 10⁻¹¹ (fitted by hand) | 0 | 4π²(π²/3−2) = 50.9220 | 50.92195 ± 10⁻⁶ | −π/2 (\|b\|<1), +π/2 (\|b\|>1), exactly |
| CEM construction, p-basis | 2.0 | log 128π⁴ | 0 | 4π⁴/3 = 129.879 | 129.8788 ± 2·10⁻⁴ | 0 |
| CCFZ (2.84), δ̂=0, ×2 | 3 | 5 log 2π | 0 | 4π⁴/3 = 129.879 | 129.8788 ± 2·10⁻⁴ | +π/2 (constant) |
| ADLS (8.18), l=0 | 3 | 5 log 2π | 0 | 4π²(π²/3−6) = −106.992 | −107.06 ± 0.55 | 0 (canonical framing) |

**Further exact results (OUR ATTEMPT, scripts C and D).**
- **Exact ratio to CCFZ.** |Z_CEM(85)| / |Z_CCFZ,δ̂=0| = (1 − e^{−S₀})²/(1 + 4π²/S₀²)² exactly; checked to 20 digits at S₀ = 3–300.
  - Both proposals have Re c = 13. CEM (126) gives b² = 1/k with k ∈ iℝ, which corresponds to δ̂ = 0; the published version states Re c = 13 is one-loop exact.
  - So their perturbative difference, c₁ difference = −8π², is entirely the h-basis Jacobian (1−b⁴)⁻².
- **(85) is not in the CS family.** Let L be the real part of the renormalized level, so that ADLS l = 0 is L = 2 and CCFZ δ̂ = 0 is L = 0. Then:
  - c₁(L) = 4π⁴/3 − 6π²L² and c₂(L) = (4π⁴/45)(135L⁴ − 180π²L² − 2π⁴).
  - Matching CEM's c₁ requires L² = 4/3. That L then gives c₂ = −20118, against CEM's c₂ = 16π⁴ − 8π⁸/45 = −128.30.
  - So (85) is not the continuation of SU(2)×SU(2) Chern–Simons at any level, already at three loops.
- **The h-basis choice drives both a and the phase.**
  - CEM's own S-kernel (83)/(5.5), continued to h = 0, 1, reproduces (84)/(5.6) up to the hand-fixed factor 2√2. The two branches of (84) correspond to S_{1,0} ∓ S_{1,1}, i.e. the sign choice of p ∝ (b − 1/b).
  - The same construction in the p basis gives S₁₁ = −8√2 sin(πb²) sin(πb⁻²). That yields a = 2 and no phase.
  - So the one-loop log coefficient and the ±i phase of (85) are outputs only of the h-basis Jacobian, which CEM did not justify.
- **CEM's own dictionary has a two-loop local term.** In the variable S_GH = (π/3) Im c, log|Z| − [S_GH − 3 log S_GH + 5 log 2π] → 4π²/S_GH (fitted ratio 0.99983 at S_GH = 10⁴). This odd (local) term is consistent with the published remark that Im c is renormalized at higher loops. It is harmless.

## 4. Verdict, precisely

**(a) The counting argument (73)–(78)/(4.32)–(4.37).**
- It passes on S_GH, −3 log S_GH and 5 log 2π **by construction**: all three are inputs through e^{S₀}.
- It **does not address** the phase (absolute value by fiat), and g_s is imaginary with "the sign … ambiguous" (text after (78)).
- Its non-trivial content, T̂² ∝ C_{S²}, is a consistency condition internal to CLS / the matrix model. CEM note it fails for the Virasoro minimal string. It is not a test against gravitational data.
- Status: ESTABLISHED by CEM's own text, and confirmed numerically by us (script A).

**(b) The exact conjecture (85)/(5.7).**
- **Leading term S_GH:** a genuine output, given Im c = 3ℓ/2G.
- **−3 log S_GH:** a genuine output **conditional** on the h-basis Jacobian; otherwise it is −2.
- **5 log 2π:** an input (prefactor fixed by matching).
- **Phase ∓i:** an output conditional on the same h-basis choice. It matches the metric one-loop phase i^{±5} up to the usual direction ambiguity, but phase data are themselves framing/contour dependent (ADLS (H.28)–(H.30); CCFZ (2.82)).
- **First real test: c₁ (order S₀⁻²).** Three published exact candidates give three different values:
  - CEM 50.92;
  - CCFZ 129.88;
  - ADLS −106.99.

  They agree on all one-loop invariants. Status: OUR ATTEMPT (scripts C, D).
- **Bottom line.** The invariant one-loop data (a, b, phase) cannot discriminate CEM from the Chern–Simons proposals, and the counting argument cannot fail them. Our verdict is **"passes only by construction"** for the count, and **"untested beyond one loop"** for (85).
- **What would decide it.** An independent three-loop computation of Euclidean S³ gravity, giving c₁. Equivalently, settling (i) the CS level shift (δ̂ = 0 vs 2) and (ii) the h- vs p-basis measure in complex Virasoro TQFT.

## 5. Literature, 2025–2026 (searched with Exa; IDs verified by abstract fetch)

- **SciPost referee report 3 on 2501.01486v1** (`fetched/scipost_…report3.txt`). It asks exactly our question: "It appears that the matching relies on Z taking the form (4.5)… would not hold if (4.5) is not correct?" It also notes that (5.7) "should be verifiable through an independent higher-loop semiclassical calculation".
  - The authors' reply could not be retrieved (bot check).
  - The published text keeps "regardless of the latter's specific value" and adds that Im c "can receive scheme dependent renormalizations at higher loops" while Re c = 13 is one-loop exact (SciPost §2, near (2.8)).
- **Blommaert–Tietto–Verlinde, 2606.26241 [verified].** They model dS₃ as CLS⊗CLS for an observer.
  - Their (5.32) quotes ρ₍₀,₁₎(α±→1) = 16 sinh²(ℏ/4) sinh²(4π²/ℏ) = S₁₁² as "completely consistent" with CEM (5.7).
  - OUR ATTEMPT (script D): with their dictionary (ℏ = 4πb², c = 13 ± 6i(b² − b⁻²)), the ratio is BTV/|CEM| = (1+β²)²/β, so a = 2 rather than 3. Their ∫∫ρ_disk ρ_disk gives a = 4.
  - The agreement therefore holds only at order e^{S_GH}. Since BTV's object is an observer density at fixed (E, J), this is not a contradiction. It does reinforce that a = 3 lives in normalization prefactors.
- **Chakravarty–Maloney–Namjou–Ross, 2505.06420 [verified].** Static-patch spectrum from a Lorentzian path integral. The abstract reports Bekenstein–Hawking agreement at the de Sitter temperature. It does not test CEM's count.
- **Collier–Eberhardt–Mühlmann, 2501.10265 [verified].** The CLS as sine-dilaton gravity, including sphere and disk partition functions. This is 2D; not checked in detail.
- **Phase follow-ups:** Maldacena 2412.14014 v3 (Nov 2025) and Chen–Stanford–Tang–Yang 2511.01400 [both verified]. Anninos et al. 2505.11330 [verified] treats 4D only.
- **Not found:** any paper that computes S³ gravity beyond one loop independently, or that tests (5.7)'s c₁. Three searches returned none. This is a search result, not a proof of absence.

## 6. Convention notes

- **Dictionary.** Im c = 3ℓ/2G and S_GH = πℓ/2G = (π/3) Im c.
- **Branches.** The |b|<1 branch (β → 0) and the |b|>1 branch (β → ∞) give the same |Z| and opposite phases.
- **Gauge group normalization.** SO(4) = (SU(2)×SU(2))/ℤ₂ multiplies Z by 2, i.e. shifts b by log 2. With SO(4), all candidates have b = 5 log 2π; the SU(2)×SU(2) form of CCFZ gives log 16π⁵.
- **Matrix measure.** e^{S₀} is a free genus-counting parameter in the double-scaled integral (2409.17246 eq. (3)). CEM fix it as a function of b only through (44)+(78). This is the precise sense in which the gravitational data are inputs.

## Files

- `r2t4_A_counting_inputs.py/.out`: identity (77), quadrature of (74), input test, and (87).
- `r2t4_B_adls_oneloop.py/.out`: independent recomputation of ADLS (H.5)–(H.7) and of ζ_{S³}(0) = −5.
- `r2t4_C_exact_candidates.py/.out`: series and fits for CEM-h, CEM-p, CCFZ and ADLS; phases; exact ratio; check of (84) from (83).
- `r2t4_D_cs_family.py/.out`: c₁(L), c₂(L) for the CS family and the BTV (5.32) comparison.
- `fetched/`: texts of 2501.01486 (arXiv HTML and SciPost published), 2009.12464, 2302.12281, 2409.17246, 2606.26241, and the referee report excerpt.
