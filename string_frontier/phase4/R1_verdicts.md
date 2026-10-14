# Phase 4 — Reviewer R1 verdicts

Reviewer: R1 (adversarial). Date: 2026-10-14. All code is in `phase4/R1_code/`. Each script was written from the primary-source formulas before I read the authors' scripts, and each has a `.out` file next to it.
In this report, "fetched" means I fetched the arXiv abstract or full text in this session with `mcp__Exa__web_fetch_exa`.

## Verdict table

| Claim | Verdict | One-line reason | R1 script |
|---|---|---|---|
| A1 one-loop D0 coefficients + BBPT citation | **FAILED (minor, one sign)** | 15/16, 0 and 9009/4096 are re-derived exactly. The v¹⁰/r¹⁹ coefficient is **+**182325/16384 in the Minkowski-Lagrangian convention the report uses for 15/16; the −182325/16384 is the Euclidean-signature coefficient, mislabelled. BBPT eqs. 13, 14 and 15 and Table 1 are verified. | `a1_a2_d0_phase_shift.py` |
| A2 three-loop target (16875/1024)N₁³N₂v⁸/(R⁷M²⁷r²¹); "no published check" | **VERIFIED-KNOWN** | The coefficient is exact, but it is just the next Taylor term of BBPT eq. (12), and BBPT explicitly name the v⁸/r²¹ three-loop term as the next test. "No published check" is plausible: I found none either. The proposed SU(2) comparison is ill-posed (see below). | `a1_a2_d0_phase_shift.py` |
| A3 BMN fuzzy-sphere spectrum; the "R/2 trap" | **VERIFIED** | DSV Table 1 is reproduced for N=2,3,4, with Σm_B = Σm_F = (2/9)N(8N²+1). The R/2 note is a *correct reading*, not an error in DSV: only the literal (1.3) coefficient iμ/8 reproduces the spectrum, and it agrees with DSV (4.4). | `a3_bmn_spectrum.py` |
| A4 Firat 1.62187e8 and λ→∞ limit 6561/4096 | **VERIFIED** (and the agent's branch choice is now justified) | Re-derived from the Fuchsian monodromy problem without Firat's eq. (5): 1.621870335e8 at λ*, tending to 6561/4096. Eq. (5) as printed in arXiv HTML gives the wrong v mod π and cannot reproduce Firat's own eq. (8) for any branch. The agent's inverted γ-ratio is the correct one. | `a4_hyperbolic_vertex_monodromy.py`, `a4_hyperbolic_vertex.py` |
| A5 ABJM Z(1)=1/4, Z(2)=1/(16π), Z/Z_Airy−1 | **VERIFIED** | Independent Nyström + elementary-symmetric (cancellation-free) evaluation. Z(3)=(π−3)/(64π) and Z(4)=(10−π²)/(1024π²) are also reproduced. The ratios agree with the agent's to ≲4·10⁻¹³ absolute. | `a5_abjm_fermigas.py` |
| A6 IKKT status claims / 2604.19836 | **VERIFIED (citations)** with a **label flag** | 1904.05914, 2307.01681, 2604.19836 and 1108.1540 were fetched and say what is claimed. Eq. (3.12) (ξ anisotropy), m_f, N=128 and "carefully examined" are all confirmed. However, "ESTABLISHED (claim of the KEK group)" is self-contradictory and should be CONJECTURE / numerical claim, and the "PTEP version" of 2604.19836 is unsupported. | — |
| B1 −571/45 and −31/45 (ADLS 2009.12464) | **VERIFIED** (reproduces the paper's own numbers) | The small-t Laurent coefficient of (1.4) with the (1.9)/(1.11) characters gives graviton −331/45 − 16/3 = −571/45 and photon −16/45 − 1/3 = −31/45. Control: a conformal scalar gives −1/90. | `b1_adls_log.py` |
| B2 w₁₊∞ bracket (2105.14346 eqs. 2.5–2.7), Jacobi | **VERIFIED** (known: Bakas 1989 / the paper itself) | 687 structure constants follow exactly from (2.5)+(2.6). Out-of-wedge targets have zero coefficient. Jacobi holds on 455 triples. | `b2_w1inf.py` |
| B3 DSSYK spectrum 2cosθ/√(1−q); Verlinde eq. (3) "omits factor 2" | **VERIFIED** (typo-level) | Spectral edge 2/√(1−q) to 10⁻⁵. Moments match q-weighted chord counting through H⁸. 2402.00635 eqs. (4)+(5) with 2πα = π−2θ give 2cosθ/√(1−q), contradicting its own eq. (3). 2505.08116 eq. (6) uses E₀=2𝒥/√(1−q). | `b3_dssyk_transfer.py` |
| **B4 Stieberger–Taylor 1806.05688: direct (3.15) ≠ (3.30)/(3.31)** | **VERIFIED — genuine error in the paper**, now resolved | The agent's numbers are reproduced to 12 digits. Cause: the u-channel poles at w = −m are missing. I derived the correct closed form and confirmed it to ≤10⁻¹⁵ (mostly ~10⁻²³) at 14 (r, β) points, including physical β=0.7i. A domain-of-validity explanation is excluded. Further consequences: the residue sign in (3.32) is wrong, and the paper's small-angle 1/r expansion misses the dominant Regge term ~½Γ(−β)(ln r)^β. | `b4_celestial_type1.py`, `b4_large_r.py` |
| B5 2025–26 citations | **VERIFIED** | 2501.01486, 2412.14014 (v3 "overall minus sign remains"), 2511.01400, 2505.08116 and 2606.26241 were all fetched; titles, authors and abstracts match. Minor slip: in 2505.08116, −AdS₂ = dS₂ is eqs. (92)–(93), not (88)–(89). Table (3.8) of 2511.01400 was not checked (fetch timed out). | — |
| S1 four derivations of D=26 | **VERIFIED-KNOWN** (textbook) | Operator implementations are correct: normal ordering, ghost vacuum, and truncation, which is exact by level counting. I hand-derived c_gh = −26. The Q² coefficients match Σ[(c−26)(m³−m)/12 − 2m(A+1)]. Parts 3 and 5 take a=1 as an input, which the script states. | (hand check below) |

---

## B4 in detail (priority claim)

**Setup, from the fetched full text of 1806.05688.**
- F_I(s,u) = Γ(1−s)Γ(1−u)/Γ(1−s−u) (3.10).
- I(r,β) = ½∫₀^∞ w^{−β−1}F_I(rw,−w)dw, with β = −(i/2)Σλ (3.15).
- The physical prescription is s → s+i0, i.e. w → w+i0 at fixed r>1.
- (3.30): I = 2πδ + (π/(2 sin πβ)) (−r)^β S(β), where S(β) = Σ_n Γ(n−n/r)/(n! n^β Γ(−n/r)).
- The paper prints (iπ/2)/sinh(½Σλ). This must be read as sinh(π/2·Σλ) to agree with its own (3.22) prefactor ½Γ(β)Γ(1−β). This is a harmless misprint.

**1. Direct evaluation (independent code).**
- I rotated the contour onto rays w = t e^{iθ} (θ = 0.15 and 0.6, giving identical results). I used the subtracted form, so the result is the meromorphic continuation for Re β < 2, and evaluated f−1 near 0 by its ζ-series.
- At r=2, β=−1.3: I₊ = −0.0358978933 + 0.1982756919i. The paper's formula (3.30), on the branch with the right discontinuity, gives +0.1440557225 + 0.1982756919i. Both numbers are identical to the agent's, so this is not a transcription or branch error by the agent.
- The discontinuity I₊ − I₋ = −iπ r^β S(β) is reproduced exactly. The formula is right only in its s-channel discontinuity.

**2. Why it fails.**
- Rewriting the Mellin integral as residues (the Hankel-contour step behind 3.31) needs the integrand to be closable in the w-plane.
- Stirling gives f(w) ≈ exp{w[iπ − K]} in the upper half plane, with K = r ln r − (r−1)ln(r−1). So f grows exponentially in a wedge |arg w − π| < arctan(π/K) around the **negative** real axis.
- That wedge contains the u-channel poles of Γ(1+w) at w = −m. Their residues R_m = (−1)^{m−1}Γ(1+rm)/((m−1)!Γ(1+(r−1)m)) grow like e^{mK}.
- APRS 2012.04208 (fetched, section 2.5) make the same general observation: "the residues of the poles at t → m … grow exponentially … which is why we can't express M(s,z) as a sum over its poles".
- In the paper's own derivation, the binomial expansion (3.19) of [r ln x − ln(1−x)]^{β−1} diverges near x→1. That is exactly the u-channel region.

**3. Corrected closed form (R1).** Solving the four linear relations between I_±, the u-side integrals J_±, and the s- and u-discontinuities gives

  I_±(r,β) = (π / (2 sin πβ)) · [ (r e^{∓iπ})^β S(β) − U_r(β) ],

  U_r(β) = Σ_{m≥1} (−1)^{m−1} m^{−β} Γ(1+rm)/(m! Γ(1+(r−1)m)).

- U is divergent but alternating. I defined it by the Lindelöf integral (1/2i)∫_{½−i∞}^{½+i∞} φ(z)/sin(πz) dz, which is analytic continuation in a fugacity x from small x to x=1.
- Checks at r=2: U(0) = 1−1/√5 and U(1) = −2 ln(2/(1+√5)). These are the closed forms from the generating function of the binomial coefficients C(2m, m), and they match to 15 digits.
- Agreement with direct integration:

| r | β | direct I₊ | (3.30) | R1 corrected | |direct−R1| |
|---|---|---|---|---|---|
| 2 | −1.3 | −0.0358979+0.1982757i | +0.1440557+0.1982757i | −0.0358979+0.1982757i | 8·10⁻²⁸ |
| 2 | 0.7i (physical) | −0.8840189+1.4339630i | −0.7798706+1.2415765i | identical | 3·10⁻²⁶ |
| 2 | 0.5 | −1.1941675+1.0404825i | 0+1.0404825i | identical | 4·10⁻²⁶ |
| 3.7 | −0.5+0.2i | 0.3061218+0.6794940i | −0.1223265+0.5574128i | identical | 9·10⁻²³ |
| 3.7 | 1.5 | 2.5055128+3.5532557i | 0+3.5532557i | identical | 4·10⁻¹⁵ |

  There are 14 points in total in `b4_celestial_type1.out`, covering both prescriptions.

- **Pole structure, a decisive test.** I_± must be analytic for Re β < 0, where the integral converges absolutely. It may have poles only at β = 0, 2, 3, … from the Taylor coefficients of f, and f has no O(w) term.
  - (3.30) has spurious poles. At β = −2 (inside the convergent region), the residue factor is (−r)^β S = −0.0358 for r=2. At β = 1 it is 0.9624.
  - In the corrected form, (−r)^βS − U = 0 to 10⁻²⁵ at β = 1 and β = −2.
  - At β=0, S(0) − U(0) = −1 exactly for r = 2, 3, 3.7 and 10. That is the required residue −½.

**4. Domain-of-validity defence excluded.** The paper derives (3.18) for r<0 and says it continues analytically to r>1.
- Along r = 1+2e^{iψ}, ψ∈[0,π] (from r=3 to r=−1 through the upper half plane), Re K(r) ≥ 1.386.
- So the real-w integral converges uniformly along the whole path. The paper's own continuation therefore equals the +i0 integral I computed.
- There is no range of r>1 or β in which (3.30) holds. Two meromorphic functions that differ on an open set differ everywhere.

**5. Consequence beyond the agent's claim (new, R1).** The paper's headline physics claim is: "QFT recovered at large r; the string corrections are a small-angle expansion in r⁻¹", i.e. I = 2πδ + O(r^{β−1}). This is false pointwise.
- Numerically, the large-r behaviour is Regge-dominated: I ≈ ½∫w^{−β−1}Γ(1+w)(rwe^{−iπ})^{−w}dw ≈ ½Γ(−β)(ln r − ln ln r − iπ)^β.
- At β=0.7i, r=10⁶: direct −0.555249−0.284371i, and the Regge integral agrees to 6 digits.
- By contrast, the paper's leading k=1 term ½Γ(β)Γ(1−β)(−r)^{β−1}ζ(β) has modulus 1.4·10⁻⁴ already at r=10⁴. Its modulus decays like 1/r.
- On the principal series, |(ln r)^{iy}| = 1, so the correction is O(1) and only oscillates as e^{iy ln ln r}. The field-theory limit can hold at most distributionally.
- A related consequence: the delta-function coefficient of the true amplitude on Re β=0 is π, not 2π, because F→0 in the UV. Equivalently, the "string correction" carries its own pole −1/(2β).
- The incorrect claim has propagated: Castiblanco–Giribet–Marin–Rojas 2405.01643 (fetched) restate "I(r,β)=2πδ+O(r^{β−1})" as their eq. (3.7). I did not check equation numbers beyond the search excerpt; (eq. not checked in full text).

**6. Novelty.**
- APRS 2012.04208 discuss the s/t split of the Euler integral for large-β asymptotics only. They do not correct the ST closed form.
- Casali–Pozzi 2602.03573 (fetched abstract) give a general "celestial dispersion relation" over poles and discontinuities, which is close in spirit to the corrected formula. I did not check whether they discuss ST.
- I found no paper stating that (3.30)/(3.31) is incomplete. This is a candidate for a genuine new correction.
- Label recommendation: B's "OUR ATTEMPT (preliminary)" can be upgraded to "OUR ATTEMPT, independently confirmed by R1, with resolution".

**7. Smaller errors in 1806.05688.**
- (3.32) has the wrong sign. Numerically, Res_{s=1}B(−s,1+s/2) = +0.5, while (3.32) gives −0.5. With the correct sign, (3.31) is not equal to (3.30).
- The sinh(½Σλ) misprint noted above.

---

## A2 in detail

- Expanding BBPT eq. (12), L′ = −p₋ẋ⁻ with ẋ⁻ = (√(1−hv²)−1)/h, gives (5/128)·h³v⁸·N₂/R = (16875/1024) N₁³N₂ v⁸/(R⁷M²⁷r²¹). This is exact (sympy). The next term is 354375/4096.
- Not new as a target. BBPT (fetched full text) write "+O(v⁸/r²¹)" in eq. (13) and state: "the skeptic is not willing to bet that the v⁸/r²¹ term, a three-loop effect, will show a discrepancy". The coefficient is a one-line Taylor expansion.
- "No published check exists": I also found no three-loop two-body computation. The searches surfaced Dine–Echols–Gray hep-th/9810021 (three-loop four-graviton v⁶ terms) and hep-th/9805007 (evidence *against* non-renormalization at v⁸). Verdict on this sub-claim: PLAUSIBLE.
- **Flaw in the proposed test (H1 "New").** At three loops the matrix-model N-dependence has two independent symmetric structures, (N₁³N₂+N₁N₂³) and N₁²N₂², plus non-planar O(N²) pieces.
  - At two loops BBPT eq. (15) used symmetry to fix a *single* structure. At three loops this is impossible.
  - An SU(2) three-loop number therefore cannot be compared with the probe coefficient. One needs general (N₁,N₂) or N₁≫N₂.
  - Also, without a v⁸ non-renormalization theorem (DEG 9805007), finite-N agreement is not expected; only leading large N is a test.

## A1 sign error (precise)

- From the background-gauge spectrum: ln Z = ∫ds/s e^{−sb²}[6+2cosh2sv−8cosh sv]/(2 sinh sv). The series is x³/2 + x⁷/480 − x⁹/6048 + x¹¹/57600 − …
- Wick rotation v_M = i v_E maps a Euclidean coefficient a_m to the Minkowski coefficient (−1)^k a_m Γ(m+½)/√π, with k=(m+1)/2.
- So the Minkowski coefficients of L = v²/2 + Σ c_k v^{2k}/r^{4k−1} are 15/16, 0, 9009/4096, **+182325/16384**, 61108047/524288, … All are positive.
- The agent's script identifies the Euclidean coefficient with the potential without the (−1)^k factor. This is harmless for k=2 and k=4 but flips k=5.
- BBPT also say the whole one-loop series is known (their ref. [13]), so these coefficients are not new.

## A3 detail

- I built real symmetric SO(9) γ's (16×16) and expanded around X_i=(μ/3R)J_i. The bosonic Hessian uses central differences. Fermions use A = −Rγ^i⊗ad(X_i) + 2c·iμγ¹²³⊗1.
- With the literal (1.3) coefficient (R/2)·(μ/4R) = μ/8, i.e. the term (iμ/8)Ψᵀγ¹²³Ψ (c = 1/8 in the script):
  - the fermion masses are (J+1/4)/3 and (J+3/4)/3 exactly, as in DSV sec. 5;
  - the bosons are (2j+1)/6, J/3, (J+1)/3, and N²−1 zero modes;
  - Σm_B = Σm_F = 14.667, 48.667 and 114.667 for N = 2, 3, 4.
- With "iμ/4" (no R/2), Σm_F = 18.67, 54.67 and 122.67, which is not supersymmetric. (The other orientation sign gives 16.0, 44.0 and 104.0.)
- So the literal reading is the only consistent one. This is the same as DSV (4.4), −(μ/4)ψ†ψ for complex ψ. Calling it a "trap" is fine; calling it an error in DSV would be wrong.

## A4 detail

- I solved the Fuchsian equation ψ''+½Tψ=0 by Frobenius series at 0 and at 1, and matched them at z=½.
- Real SL(2,ℝ) monodromy around z=1 fixes v ≡ −arg(αβ) (mod π) = 2.50610 at λ* = arcsinh(1)/π.
- I took the geodesic boundary of hole 1 as the waist λ log|ρ₁| = π(l̃+½) of the annulus containing the seam z=½. This gives l̃=−1, consistent with Firat's statement for λ≤10.
- Mapping radius r = 0.0356569. Then (1/3)r⁻⁶ = 1.621870335×10⁸, matching Firat eq. (8).
- For λ = 1, 2, 5, 10, 30 and 100 the coefficient is 9.6077, 2.4604, 1.7127, 1.6288, 1.60478 and 1.60207, tending to 6561/4096 = 1.601807.
- Firat eq. (5) as rendered gives v ≡ 1.24922 (mod π). With that value, eq. (7) gives 5.19×10²⁵ (principal branch), and no integer shift of v or l̃ gives 1.62×10⁸. The agent's variant, with the γ-ratio inverted, equals the monodromy value.
- So the printed (5) appears to carry an inverted γ-ratio. It could also be an HTML-rendering artefact; I could not see the PDF typography.

## S1 detail (operator checks)

- **Ghosts.** By hand, L₋₂|0⟩ = −(b₋₃c₁ + 2b₋₂c₀)|0⟩. L₂ returns −5−8 = −13 times |0⟩, and L₂|0⟩=0, so c_gh = 2(−13 − 4·0) = −26. L₀c₁|0⟩ = −c₁|0⟩.
- **Kac determinant.** The Gram matrix [[4h+c/2, 6h],[6h, 4h(2h+1)]] gives det = 2h(16h²+2h(c−5)+c). With L₁ψ = (3+g(4h+2)), L₂ψ = c/2+4h+6gh and h=−1: g=3/2, c=26.
- **BRST.** The script's output (−2(A+1), (c−26)/2−4(A+1), 2(c−26)−6(A+1)) equals (c−26)(m³−m)/12 − 2m(A+1) for m=1, 2, 3.
- **Truncation.** Q preserves total level, so terms with |n| > K vanish identically on the test states, and K=6 and K=8 agree.
- Conclusions follow. The script's references and labels are honest (textbook ESTABLISHED).

## Label audit (A and B reports)

- A: "ESTABLISHED (claim of the KEK group)" for the Lorentzian ≡ Euclidean equivalence should be CONJECTURE / reported numerical result.
- A: "Higher one-loop coefficients" are labelled OUR ATTEMPT but are known (BBPT ref. [13]). The −182325/16384 sign must be fixed.
- A: The three-loop "target" should be presented as BBPT's eq. (12) expanded, not as a derived new target. The SU(2) comparison plan must be revised (two N-structures).
- B: The B4 item is honestly labelled preliminary, and it is correct. B's H3 "settle the representation" is now done (section B4 above).
- B: The B3 remark ("normalization or typo") is appropriately modest.
