# RB4 — independent second referee: Stieberger–Taylor arXiv:1806.05688, celestial four-gluon type-I amplitude

Referee RB4. I wrote all code in `phase4/RB4_code/` from the fetched paper text **before** opening
`phase4/R1_code/`, `phase1/B_checks/` or `R1_verdicts.md`. The comparison with those files is in §6.

Sources fetched in this session with `mcp__Exa__web_fetch_exa`. Title and authors were checked on the abstract page.
- **1806.05688.** Stieberger, Taylor, "Strings on Celestial Sphere", v2 of 1 Sep 2018. Full text saved as `RB4_code/paper_1806.05688.txt`.
- **2405.01643.** Castiblanco, Giribet, Marin, Rojas. Full text fetched.
- **2012.04208.** Arkani-Hamed, Pate, Raclariu, Strominger. Full text fetched.
- **2106.11948.** Chang, Huang, Huang, Li. Full text fetched.
- **2408.02609.** Bockisch. Full text fetched.
- **2602.03573.** Casali, Pozzi. Full text fetched.

## 0. Verdict table

| sub-claim | verdict | one-line reason | script |
|---|---|---|---|
| (a) The pole sum (3.30)/(3.31) ≠ the direct integral (3.15), because the u-channel poles are dropped | **VERIFIED** (genuine error in the paper, not a branch, convention or domain slip) | (3.30) fails at every tested point for every reading: both iε, both branches of (−r)^β, sinh read as printed or with π, and (3.31) with either residue sign. That includes the physical β = 0.7i. The difference equals the u-channel term to 10⁻³⁰. | `rb4_direct.py`, `rb4_formulas.py` |
| (b) Corrected formula I_± = π/(2 sin πβ)[(r e^{∓iπ})^β S(β) − U_r(β)] | **VERIFIED** | I derived it myself (keyhole contour). I continued U_r by a *different* method (Abel/Mellin generating function, not Lindelöf). It matches direct integration to 10⁻¹⁷ … 10⁻³¹ at 8 points, including β = 0.7i. The residue at β=0 is −½: numerically exact, and **proved analytically** below. Every spurious pole cancels (β = 1, −1, −2, −3 checked). | `rb4_formulas.py`, `rb4_poles_larger.py` |
| (c) The paper's "small-angle expansion in r⁻¹" (2405.01643 eq. 3.7: I = 2πδ + O(r^{β−1})) is false; I ~ ½Γ(−β)(ln r)^β | **VERIFIED, with a domain qualifier**: true for Re β < 1, which includes the physical principal series. For Re β > 1 the r^{β−1}ζ(β) term does dominate. | At β = 0.7i, \|I\| stays ≈ 0.5 from r = 10 to 10³², while the paper's leading term is 10⁻³². The ratio I/[½Γ(−β)(ln r)^β] → 1 (1.04 at r = 10³²). Beyond the claim: even distributionally, I(r,·) → 0, not 2πδ (see §4.3). | `rb4_poles_larger.py`, `rb4_distributional.py` |
| (d) Residue sign in (3.32); sinh prefactor misprint | **VERIFIED** (both minor) | Numerically Res_{s=n}B(−s,1+s/r) = −Γ(n−n/r)/(n!Γ(−n/r)): +0.4 vs the printed −0.4 at r=2.5, n=1. (3.31) carries a compensating overall sign, so (3.31)+(3.32) = (3.30) as printed. sinh(½Σλ) must be sinh(π/2 Σλ) to match (3.22). | `rb4_formulas.py` |

Labels. The claim is correctly labelled **OUR ATTEMPT**. After R1 and RB4 it can be called "independently verified twice".
- The *ingredients* are ESTABLISHED in the literature (§5): cross-channel residues grow exponentially, and the Im part is the s-channel residue sum.
- I found no paper that states the ST formula is wrong or gives the corrected closed form. So VERIFIED, not VERIFIED-KNOWN.

---

## 1. Definitions read from the fetched 1806.05688 text (equation numbers as cited in the paper's own prose)

**Kinematics and form factor.**
- (3.2) r = z₁₂z₃₄/(z₂₃z₄₁).
- (3.3) 1/r = −u/s = sin²(θ/2), so the physical region is r > 1, with s > 0 and u = −s/r < 0.
- (3.10) F_I(s,u) = Γ(1−s)Γ(1−u)/Γ(1−s−u).
- The paper notes the s-channel poles at s = n "can be handled by using the iε prescription".

**The integral.**
- (3.15) I(r,β) := ½∫₀^∞ w^{−β−1}F_I(rw,−w)dw, with β := −(i/2)Σλ_n.
- On the principal series Δ_n = 1+iλ_n, β is purely imaginary.
- Integration variable: w = s/r = −u ∝ ω₄². The contour is the positive real w axis; the s-channel poles lie on it at w = n/r.

**Convergence.**
- (3.15) converges at w→0 iff Re β < 0, since F→1.
- At w→∞ it converges for every β, by the exponential decay (3.12).
- On Re β = 0 it is a distribution. The "2πδ(Σλ)" is the w→0 end.
- Since F − 1 = ζ(2) r w² + O(w³), the continuation in β has poles only at β = 0, 2, 3, 4, …

**Derivation route in the paper.**
- (3.17)/(3.18): the Euler integral is valid for r < 0 only. The paper says it will "perform analytic continuation to r>1 at the end".
- (3.19): binomial expansion of [r ln x − ln(1−x)]^{β−1} in powers of ln(1−x)/(r ln x).
- (3.22) then (3.26), giving the MZV series in (−r)^{−k}.
- (3.27)–(3.29) give (3.30): I = 2πδ + (iπ/2)(−r)^β sinh(½Σλ)^{−1} S_r(β), where S_r(β) = Σ_n Γ(n−n/r)/(n! n^β Γ(−n/r)).
- (3.31): Hankel form with Σ_n Res_{s=n}{(s/r)^{−β}B(−s,1+s/r)}.
- (3.32): the residue.
- No restriction on r>1 or on β is stated for (3.30)/(3.31). (3.24) is derived for Re β > 1. The result is applied on the principal series.

**My conventions.** I_± means s → s ± i0, i.e. r → r ± i0.
- Continuing the r<0 Euclidean integral through Im r > 0 gives exactly I_+. Checked: method B below is the paper's own (3.18) continued this way.

## 2. Direct evaluation (two independent methods), `rb4_direct.py`

**Method A: rotated w contour.** w = t e^{iθ}, subtracted so that it is valid for Re β < 2. F−1 near 0 is computed from the ζ-series of ln F plus expm1.
- At r=2, β=−1.3 the result is identical to 20 digits for θ = 0.05, 0.3, π/4, π/2, 2.0, 2.4 rad.
- θ = 2.4 lies in the second quadrant. This is allowed because the growth wedge is only 23.8° wide (§3).

**Method B: the paper's own (3.18), continued to r>1 with A − i0.** Here A = r ln x − ln(1−x). The path is indented below the zero x* of A, and x = e^{−t} is used for x < x*.
- Valid for Re β < 0.

| r | β | Method A (I₊) | Method B | \|A−B\| |
|---|---|---|---|---|
| 2 | −1.3 | −0.035897893322497804692 + 0.19827569193216443668i | same | <10⁻²⁰ |
| 2 | −0.5 | 0.555449581288710692 + 0.458477198302025087i | same | 4·10⁻¹⁸ |
| 1.3 | −2.2 | −0.15408634147384454 − 0.014606511315265778i | same | 5·10⁻³¹ |
| 3.7 | −0.4+1.1i | −0.728028265082892603 + 0.374432845900555414i | same | 8·10⁻¹⁴ |
| 10 | −0.25 | 1.404727287203646 + 0.367764141398259i | same | 3·10⁻⁹ (slow algebraic tail in B) |

- I₋ = conj(I₊) for real β and real r.
- The number quoted in claim (a), −0.035898 + 0.198276i, is confirmed.

## 3. Sub-claim (a): is the pole sum wrong, or misread?

### 3.1 The printed formula fails under every reading (`rb4_formulas.py`, section 3)

I evaluated (3.30) in every reading:
- sinh(½Σλ) = sinh(iβ) = i sin β, as printed, and the intended sinh(iπβ) = i sin πβ;
- (−r)^β = r^β e^{±iπβ};
- (3.31) evaluated literally with (3.32) and with the true residue.

None matches I₊ or I₋ at any tested point:

| r | β | direct I₊ | (3.30) intended, e^{−iπβ} branch (closest) |
|---|---|---|---|
| 2 | −1.3 | −0.0358979 + 0.1982757i | +0.1440557 + 0.1982757i |
| 2 | −0.5 | 0.5554496 + 0.4584772i | 0 + 0.4584772i |
| 3.7 | −0.4+1.1i | −0.7280283 + 0.3744328i | −0.7312981 + 0.3136471i |
| 5 | 0.7i (physical) | −1.0499931 + 0.9546326i | −0.8817859 + 0.7625807i |
| 2.5 | 0.3 | −2.0335117 + 0.8893919i | −0.6461810 + 0.8893919i |
| 10 | −0.25 | 1.4047273 + 0.3677641i | 0.3677641 + 0.3677641i |
| 1.5 | 1.5 | 1.5326686 + 1.7074340i | 0 + 1.7074340i |
| 4 | −2.7 | −0.0500628 + 0.0112011i | −0.0081381 + 0.0112011i |

- For real β the paper's formula reproduces **exactly the imaginary part**. That is the s-channel discontinuity.
- I checked I₊ − I₋ = −πi r^β S_r(β) to 10⁻³⁰ at four points, including β = 0.7i and complex β.
- The real part is wrong.

### 3.2 Why it fails (mechanism, re-derived)

1. **Growth wedge.** F(rw,−w) grows exponentially in a wedge around the **negative** w axis, where the u-channel poles of Γ(1+w) sit at w = −m.
   - Stirling gives |F| ~ exp[R(K cos φ − π sin φ)] on w = R e^{i(π−φ)}, with K = r ln r − (r−1)ln(r−1).
   - Checked numerically in `rb4_wedge.py`. The wedge half-angle is arctan(K/π): 23.8° at r=2 and 46.0° at r=10.
   - So the Hankel contour behind (3.31) cannot be closed. The big arc does not vanish, and the u-channel residues grow like e^{mK}.
2. **Keyhole bookkeeping.** With the w^{−β−1} cut on the positive axis:
   - 2I₊ − 2e^{−2πiβ}I₋ = 2πi Σ_u Res − (arc);
   - I₊ − I₋ = −πi Σ_s Res;
   - so I₊ = π/(2 sin πβ)[(r e^{−iπ})^β S_r(β) + S_{1/r}(β)] + (arc term), where S_{1/r}(β) = Σ_n (−1)^n C(rn,n) n^{−β} is the formal u-channel residue sum.
   - (3.30)/(3.31) is the first term only.
3. **The same error in the paper's own derivation.** The binomial expansion (3.19) converges only for x < x*, where r ln x = ln(1−x).
   - Near x → 1, the u-channel region, the term-wise integrals S(β−k,k) diverge for k > Re β − 1. They are defined only by continuation, which loses the x > x* contribution.
   - APRS (2012.04208, §2.5) split the same Euler integral at exactly this y* for that reason.
4. **No domain-of-validity defence.** I₊ is the continuation from r<0 (Method B agrees). S_r(β) diverges for r < 1, so (3.30) is defined only on r > 1, where it disagrees with the integral.
   - Two functions analytic on r > 1 that differ at points of that set cannot be rescued by a domain restriction.
   - The pole structure in β (§4.2) also differs, which no convention can fix.

**Verdict (a): VERIFIED.**

## 4. Sub-claims (b) and (c)

### 4.1 Corrected formula, my independent continuation (`rb4_formulas.py`)

I define S_{1/r}(β), which equals −U_r(β) in the claim's notation, by Abel/Mellin continuation of the generating function:
- G_r(z) = Σ_{n≥0} C(rn,n) zⁿ = B/(r − (r−1)B), with B = 1 + zB^r (GKP 5.61). Checked against the series at z = −0.05 to 25 digits.
- S_{1/r}(β) = Γ(β)^{−1}∫₀^∞ t^{β−1}[G_r(−e^{−t}) − 1]dt.
- On [0,1] I integrate the Taylor series of G at t=0 term by term. It has radius ≥ π, because the only singularity sits at z = 1/C > 0. The Taylor series comes from a power-series Newton solve, which makes it analytic in all β.
- This is a different continuation from R1's Lindelöf integral. Both give the unique analytic continuation.
- Equivalent closed form for Re β > 0:

  I_±(r,β) = π (r e^{∓iπ})^β S_r(β)/(2 sin πβ) + ½Γ(1−β) ∫₀^∞ t^{β−1}[G_r(−e^{−t}) − 1] dt.

**Agreement with direct integration**, |corrected − direct|:

| r | β | \|corrected − direct\| |
|---|---|---|
| 2 | −1.3 | 1·10⁻³¹ |
| 2 | −0.5 | 1·10⁻³¹ |
| 3.7 | −0.4+1.1i | 1·10⁻³¹ (I₊), 3·10⁻³² (I₋) |
| 5 | 0.7i | 1·10⁻³⁰ (I₊), 3·10⁻³² (I₋) |
| 2.5 | 0.3 | 5·10⁻³¹ |
| 10 | −0.25 | <10⁻³⁰ |
| 1.5 | 1.5 | 9·10⁻¹⁸ (limited by the direct integral near Re β → 2) |
| 4 | −2.7 | 1·10⁻³¹ |

**Pole structure** (r = 2.5 and r = 7; `rb4_poles_larger_out.txt`).
- The true I has poles only at β = 0 (residue −½) and at β = 2, 3, …; the residue at β=2 is −ζ(2)r/2.
- Paper (3.30): nonzero residues at β = −3, −2, −1 and +1. For r = 2.5 these are −0.0104, −0.0192, +0.0807 and +0.5307.
- Paper (3.30): wrong residue at β = 0 (−0.2153) and at β = 2 (−1.2975, against the true −2.0562 = −ζ(2)·2.5/2).
- Corrected formula: residue 0 (to 10⁻³¹) at β = −3, −2, −1, +1; −0.5 at β = 0; −2.0561676 at β = 2.
- So the paper's spurious poles occur at **every** integer β ≤ 1 other than 0, not only at 1 and −2. Its residues at the true poles are also wrong.

**Analytic proof of the residue −½ at β=0** (new, RB4). Let x ∈ (0,1) solve x + x^r = 1. Then:
- B_r(−1) = x, and S_{1/r}(0) + 1 = x/(r − (r−1)x).
- For t = 1 − 1/r, B_t(1) = 1/(1−x), because y = 1−x = x^r gives 1/y = 1 + y^{−t}. GKP 5.61 with exponent −1 then gives S_r(0) + 1 = r(1−x)/(r − (r−1)x).
- Summing, S_r(0) + S_{1/r}(0) + 2 = [r − (r−1)x]/[r − (r−1)x] = 1.
- Hence S_r(0) + S_{1/r}(0) = −1 for every r > 1, and the corrected residue ½[S_r(0) + S_{1/r}(0)] = −½ exactly.
- The paper's residue ½S_r(0) = −½·r(1−x)/(r−(r−1)x) is never −½.

**Verdict (b): VERIFIED**, by an independent derivation and an independent continuation.

### 4.2 Large r (small angle), `rb4_poles_larger.py`

Direct I₊ is computed with working precision 30 + 2·log₁₀ r digits. Note that F involves loggamma(1−rw) − loggamma(1−(r−1)w), so the precision must exceed log₁₀ r.

| r | β = 0.7i: I₊ | I₊/[½Γ(−β)(ln r)^β] | \|paper lead\| ∝ r^{β−1} | β = −0.5: I₊/[½Γ(−β)(ln r)^β] |
|---|---|---|---|---|
| 10 | −1.05367 + 0.67123i | 2.469 + 0.315i | 0.142 | 0.824 + 0.449i |
| 10² | −0.87860 + 0.18119i | 1.787 + 0.010i | 1.4·10⁻² | 0.979 + 0.331i |
| 10⁴ | −0.67700 − 0.13387i | 1.370 − 0.111i | 1.4·10⁻⁴ | 1.051 + 0.191i |
| 10⁸ | −0.45960 − 0.37187i | 1.174 − 0.096i | 1.4·10⁻⁸ | 1.055 + 0.097i |
| 10¹⁶ | −0.20287 − 0.50507i | 1.083 − 0.063i | 1.4·10⁻¹⁶ | 1.040 + 0.047i |
| 10³² | +0.06494 − 0.51829i | 1.040 − 0.038i | 1.4·10⁻³² | 1.025 + 0.023i |

At β = 0.5 the ratio goes 0.68 → 0.96 (r = 10 → 10³²), and |I₊| grows like (ln r)^{1/2} (14.6 at 10³²).

- **Principal series β = 0.7i.**
  - |I₊| ≈ 0.5–1.2 from r = 10 up to 10³² (0.522 at 10³²). It does not decay.
  - The paper's / 2405.01643's leading term ∝ r^{β−1} has modulus 0.142/r.
  - I₊/[½Γ(−β)(ln r)^β] → 1, slowly, like O(ln ln r/ln r).
  - The Regge form ½∫w^{−β−1}Γ(1+w)(rwe^{−iπ})^{−w}dw reproduces I₊ for β = −0.5 to 10⁻⁸ at r = 10⁸ and to ≥8 printed digits (imaginary part of the ratio 3·10⁻²⁸) at r = 10³².
- **Mechanism.** Analytically, the u-channel term alone gives h(t) ≈ −W(re^{−t})/(1+W). That is ≈ −1 for t < ln r, so the u-term is ≈ ½Γ(−β)(ln r)^β. The s-term is O((ln r)^{β−1}). At r = 10³ and β = 0.7i, the s-term is −0.3699−0.0339i and the u-term −0.3866+0.0213i.
- **Qualifier.** For Re β > 1 the paper's leading term π(−r)^{β−1}ζ(β)/(2 sin πβ) does dominate. At β = 1.5, Im I₊ = 4.1035·10⁸ (10¹⁶) and 4.1035·10¹⁶ (10³²), equal to that term, while Re I₊ (233.6, 692.8) tracks ½Γ(−β)(ln r)^β (ratio 0.88, 0.93). The O(r^{β−1}) statement fails exactly for Re β < 1, which is where the physical β ∈ iℝ lies.
- **The statement under review.** 2405.01643 eq. (3.7), fetched full text: "I(r,β) = 2πδ(λ₁+λ₂+λ₃+λ₄) + O(r^{β−1})". 1806.05688 §3: "small scattering angle expansion in r⁻¹ … quantum field theory is recovered in the kinematic limit of large r".

### 4.3 Beyond the claim: the forward limit fails even distributionally (`rb4_distributional.py`)

I pair the full I(r, iy) with φ(y) = e^{−y²/2}:

  ⟨I_r,φ⟩ = ½∫dw/w F(rw+i0,−w) φ̂(ln w).

- The contour is rotated to w = it. The pairing was cross-checked at r=10 against the y-integral on Re β = −0.3: 0.20432372+0.51885857i from both routes.
- For Yang–Mills (F ≡ 1) the pairing is exactly π, i.e. 2πδ(Σλ).

| r | ⟨I_r,φ⟩ | \|⟨I_r,φ⟩\|/π (YM/paper prediction: 1) |
|---|---|---|
| 10^{0.5} | 0.17832 + 0.98286i | 0.318 |
| 10 | 0.20432 + 0.51886i | 0.178 |
| 10² | 0.16858 + 0.21184i | 0.086 |
| 10⁴ | 0.07945 + 0.05604i | 0.031 |
| 10⁸ | 0.02238 + 0.00922i | 7.7·10⁻³ |
| 10¹⁶ | 0.00409 + 0.00099i | 1.3·10⁻³ |
| 10³² | 5.09·10⁻⁴ + 7.2·10⁻⁵i | 1.6·10⁻⁴ |
| 10⁶⁴ | 4.38·10⁻⁵ + 3.6·10⁻⁶i | 1.4·10⁻⁵ |

- The pairing goes to **0**, not π.
- Reason: at fixed w, F(rw,−w) ≈ Γ(1+w)(rw)^{−w} → 0 (Regge suppression). F ≈ 1 survives only for w ≲ 1/ln r, which shrinks away on the log-scale relevant to the Mellin transform.
- So the stringy factor does **not** reduce to the YM δ-function at r→∞ in any sense. The forward limit of the celestial string amplitude is Regge-dominated, not massless-exchange dominated.
- This is OUR ATTEMPT (RB4), verified numerically. It is stronger than R1's remark that the limit "can hold at most distributionally".
- Caveat: only one test function was used. The analytic argument (dominated convergence on the rotated contour) is general.

**Verdict (c): VERIFIED for Re β < 1, which includes the physical principal series.** The phrase "at large r the result behaves as ½Γ(−β)(ln r)^β" needs the qualifier Re β < 1.

## 5. Sub-claim (d), and the literature

**(d) Residue sign.** At r=2.5, n=1 the numerical contour residue of B(−s,1+s/r) at s=n is +0.4, against the printed −0.4. At n=2 it is +0.08 against −0.08; at r=3.3, n=3 it is +0.015026 against −0.015026.
- The true residue is −Γ(n−n/r)/(n!Γ(−n/r)).
- Inserting the true residue into (3.31) gives −(3.30). So (3.31) also has an overall sign slip, which compensates the slip in (3.32).

**(d) sinh misprint.** With ½Σλ = iβ, matching (3.22)'s prefactor ½Γ(β)Γ(1−β) = π/(2 sin πβ) requires sinh(π/2·Σλ). The same form appears in (3.26) and (5.12).
- The printed sinh(½Σλ) = i sin β gives yet another (wrong) function (§3.1, "printed" rows).
- **VERIFIED, minor.**

**Literature and novelty** (all fetched this session):
- **APRS, 2012.04208, §2.5.** Open-string residues in the cross channel "grow exponentially at large m, which is why we can't express M(s,z) as a sum over its poles". They split the Euler integral at y*. They cite ST as "a related analytic continuation" but do not say ST's closed form is wrong.
- **Chang–Huang–Huang–Li, 2106.11948.** Eqs. (4.11)/(4.14): Im Ψ is the s-channel residue sum, checked numerically for open strings. The u-channel celestial string amplitude "cannot be expressed as a sum over residues". The convergence fan has an r-dependent angle. ST is cited only generically, as one of [1–21]; there is no statement that (3.30) is incomplete.
- **Bockisch, 2408.02609, §6.1.** "Confirms" ST at the poles β = 2, 3, 4 by *truncating* the MZV sum to its absolutely convergent terms (his eq. 80). That truncation silently discards exactly the pieces the u-channel term cancels. He did not notice the error. The full (3.30) gives the residue −1.2975 at β=2 (r=2.5), against the true −2.0562.
- **Castiblanco et al., 2405.01643.** Eq. (3.7) repeats the incorrect O(r^{β−1}) statement.
- **Donnay et al., 2307.03551.** Repeats ST's tree-level field-theory-limit narrative; does not re-derive it.
- **Casali–Pozzi, 2602.03573.** Celestial dispersion relation and Regge limit. Does not discuss 1806.05688; the only Stieberger–Taylor work it cites is FFSTZ 2108.10337.
- I found no erratum. The arXiv v2 (1 Sep 2018) contains (3.30)–(3.32) as analysed.
- **Conclusion: not previously reported.** The ingredients are ESTABLISHED (APRS, CHHL); the specific correction is not in the literature I could find.

## 6. Comparison with R1 / B files (read only after the above)

**Agreement.** All of R1's and B's numbers that I spot-checked agree with mine to every printed digit:
- r=2, β=−1.3;
- r=2 and r=3.7 at β=0.7i;
- the u-term at r = 10, 100, 1000 (β = 0.7i);
- direct I₊ at r = 10⁴ (β = 0.7i): −0.67699569 − 0.13386917i;
- R1's U(0) = 1 − 1/√5 at r = 2, consistent with my sum rule.

Two different continuations (Lindelöf in R1, Abel generating function in RB4) agree, as they must.

**Discrepancies (all minor; none affects a verdict):**
1. R1 writes the growth wedge as |arg w − π| < arctan(π/K). It is arctan(**K/π**): 23.8° at r=2, not 66.2° (`rb4_wedge_out.txt`).
   - R1's exponent exp{w[iπ − K]} is right; the angle is inverted.
   - Evidence: my contour at arg w = 2.4 rad, i.e. 42.5° from the negative axis, converges and gives the same value to 20 digits.
2. R1 lists spurious poles at β = 1 and −2. In fact they occur at every integer β ≤ 1 except 0, and the residues at the true poles β ≥ 2 are also wrong.
3. R1 says the δ coefficient "is π, not 2π". This depends on convention:
   - 2πδ(Σλ) plus the analytic remainder taken as a boundary value from Re β > 0 equals the true distribution;
   - πδ(Σλ) plus the boundary value from Re β < 0 also equals it.
   - It is not an independent error of the paper. The real defect is the wrong residue at β=0 in (3.30).
4. R1 says the field-theory limit "can hold at most distributionally". RB4 finds it fails distributionally too (§4.3).

## 7. Files

All under `phase4/RB4_code/`:
- `rb4_direct.py` (+`_out.txt`): direct integral, two methods.
- `rb4_formulas.py` (+`_out.txt`): (3.30)/(3.31)/(3.32) in every reading, the discontinuity check, the corrected formula.
- `rb4_poles_larger.py` (+`_out.txt`): pole structure; large r up to 10³²; s/u decomposition.
- `rb4_distributional.py` (+`_out.txt`): weak limit at r→∞.
- `rb4_wedge.py` (+`_out.txt`): growth wedge.
- `paper_*.txt`: fetched source texts.
