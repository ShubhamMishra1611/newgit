# RT3 — Adversarial review of T3 (BP discretuum + MSW anthropics; type II minimal slopes vs DESI)

Reviewer: RT3. Target: `phase3/targets/T3/T3_report.md`, its scripts, and its outputs.
Protocol: `phase4/REVIEW_PROTOCOL.md`.

**Method.**
- I wrote all of my code from scratch, from the primary formulas, before reading T3's scripts. My enumeration library is `RT3_code/rt3_bp_core.py`.
- I used T3's scripts and outputs only afterwards, to find where results differed.
- Code and outputs are in `phase4/RT3_code/`. Fetched texts are in `phase4/RT3_code/sources/`.

## Verdict table

| # | Claim (T3) | Verdict | One-line reason | RT3 script |
|---|---|---|---|---|
| 1 | BP flatness: N(Λ∈[0,ε]) ∝ ε^p, p = 0.994–1.008, no gaps, χ²/dof ≈ 0.9 | **VERIFIED** | Completeness: an independent algorithm reproduces T3's J = 4 window counts exactly. Flatness: new seeds give p = 0.993–1.012 and pooled χ²/dof = 0.88–1.00. Power: the test has power (χ²/dof = 171 → 3.6×10⁵ for near-commensurate charges). | `rt3_bp1_identity_completeness.py`, `rt3_bp2_flatness_power.py`, `rt3_bp2c_seeds.py` |
| 2a | Exact identity #orthant = 2^−J Σ_U #(ℤ^U∩E); face-corrected density matches to ≤ 0.2 % | **VERIFIED** | The identity holds exactly by brute force (J = 2–5). Measured density / face formula = 0.998–1.005 on 2×10⁵ levels, i.e. within 2σ statistics. | `rt3_bp1_…`, `rt3_bp2_…` |
| 2b | BP (2.26) with D = 2^J underestimates the distinct density: 2.2× at n_max ≈ 30, 10^2.2–10^4.2 at J = 100 | **VERIFIED** (real effect, not a convention mismatch) | BP define D as the *typical* degeneracy and set it to 2^J. Points on coordinate faces have only 2^{#nonzero} sign images. Exact generating-function lattice counts at J = 100 agree with T3's face formula to 10^0.00. Caveat: BP's own J = 100 headline (q = 10^−1.6) used "D small", not 2^J; against that headline the true density is about 10^27 *lower*. | `rt3_bp3_largeJ_faces.py` |
| 2c | "…and by 10^60.7 in BP's large-J examples" (J = 300, q = 10^−0.87; §9 "their large-J examples") | **FAILED** | (i) BP contains no J = 300 example. (ii) The case has n_max = 10.5 < √J, which is outside the validity BP state for (2.24)–(2.26); BP say such fluxes "should simply be ignored". (iii) T3's face formula is itself wrong there: the exact lattice count gives a gain of 10^62.35, not 10^60.7. That is an error of 10^1.69. | `rt3_bp3_largeJ_faces.py` |
| 3 | Poisson spacings: var s ≈ 1, P(s<0.1) = 0.095, Σ²(L) ≈ L | **VERIFIED** | New seeds give var s = 0.979–1.017, P(s<0.1) = 0.094–0.096, KS p mostly > 0.1, and Σ²(L) = L within about 12 % up to L = 10³. This is expected Berry–Tabor behaviour, and T3 does not claim it as new. | `rt3_bp2_…`, `rt3_bp2c_…` |
| 4 | min\|Λ\| ~ Exp, mean 1/(2ρ(0)); pooled ⟨z⟩ = 0.943 ± 0.030, KS p = 0.28 | **VERIFIED** | 360 new exactly-checked shells give ⟨z⟩ = 0.932 ± 0.047, var 0.81, KS p = 0.51. Both runs sit about 6 % low, which is not significant. The R^−(J−2) scaling holds only for n_max ≫ J, as claimed: the exact lattice count matches the face formula at n_max ≥ 14. | `rt3_bp2_…`, `rt3_bp3_…` |
| 5 | Near-commensurate charges give gaps of up to 2×10⁶ spacings; repeated charges give J_eff | **VERIFIED** | Max gap 2.01×10⁶ at δS = 1.2; min\|Λ\| = 0.136 at δS = 0.12. Repeated charges: 315 and 1059 distinct values, reproduced exactly by an independent simplex count. One nuance is added in the details section. | `rt3_bp2_…`, `rt3_bp4_repeated.py` |
| 6 | AD/DD factor of 2 is a convention; #(fluxes, τ∈F) = 2Σσ + elliptic terms exactly | **VERIFIED** | Brute force is exact for L = 10–80. DD's own recipe gives 120 orbits at τ = 2i, so DD's "240" counts N and −N separately. The ρ-term is vacuous (ρ ∉ ℚ(i)). T3's disc counts at R = 0.20 and 0.25 are ill-posed: there are lattice points exactly on the circle. | `rt3_fac2.py` |
| 7 | MSW with Planck 2018: P = 0.086 at R_G = 1 Mpc, median 13; 0.043–0.29 over 0.5–3 Mpc; MSW reproduced first | **VERIFIED** | I re-derived MSW from the collapse model rather than their closed forms. This reproduces eq. 25 (2.5788/3.0103/3.2235), the quantiles (0.05897/1.4852/11.113) and all 18 σ entries of Table 1. My own CAMB run gives σ₀ = 2.1055 and P = 0.0857. EH98 and BBKS transfer functions move P by ≤ 5–7 %. Earlier MSW updates (Weinberg 2005, WMAP-1; Pogosian–Vilenkin 2007) are not cited. | `rt3_msw.py` |
| 8 | Min slopes: √(18/7) (IIB, R6+D9), √(50/19) (IIA, R6+F0), √2 with Q; pointwise | **VERIFIED** (math) — with a physics caveat | Kinetic terms come out exactly (3/4)(∂lnρ)² + (∂lnτ)² (sympy, 10d → 4d). All exponent vectors equal HKTT (2.5, 3.3, 4.1, 4.2, 4.4) and STT 2306.07327 §IV.1. The hull minima are reproduced exactly, and an ODE gives late-time ε = 9/7, 25/19 and 1. Caveat: the IIB minimiser needs a net positive D9-like tension, which IIB proper forbids (C₁₀ tadpole). Without it the sharp IIB value is √(8/3) = 1.633. 1.604 remains a valid lower bound. The pointwise bound is the standard convex-hull lemma and is not new. | `rt3_slopes.py`, `rt3_slopes_ode.py` |
| 9 | Mixed signs: no stable 2-term attractor; collinear (+,−,+) reach √(2/3) only toward strong coupling or zero volume; formal dS points with a flat direction | **VERIFIED** (2-term part: VERIFIED-KNOWN, STT 2303.03418 p. 3) | An independent line enumeration finds R6(+)+O4(−)+F4(+) and Q+O8+F2 at √(2/3). ODE check: ε = 1/3, with d ln g_s/dN = +1. Restricted to weak coupling and non-shrinking volume, the minimum is 1.622 with geometric sources and √2 with Q. The flat-direction lines are R6+O3+(F5/Cas) and F1+O9+Q, the no-scale (GKP-type) structure. | `rt3_slopes.py`, `rt3_slopes_ode.py` |
| 10 | DESI: λ = 0.698 is "4.1σ below √2, 5.2σ below 1.604" | **PLAUSIBLE** | The arithmetic is right ((√2−0.698)/0.173 = 4.14, and so on). But each σ is a Gaussian extrapolation of a 68 % interval into a 4–5σ tail: a λ²-Gaussian reading gives 6.3σ for √2. The fit uses DESY5 supernovae and compressed CMB priors. The headline labels these σ values "ESTABLISHED + OUR ATTEMPT"; they are OUR ATTEMPT only. | arithmetic in the details section |
| L1 | Literature IDs and quoted numbers (2411.07301, 1607.00180, 2303.03418, 2306.07327, 2209.11821, 2504.04226, 2101.11617, 1909.11063, 2309.03938, …) | **VERIFIED** | All 26 IDs I fetched exist and say what is claimed: P = 0.5 % and median 539 Λ_obs; 6.7/9.0/9.7 %; STT (II.3); CRV γ < √2 F-theory candidates; λ = 0.698 (+0.173/−0.202). The one exception is the "λ ≲ √3" attributed to 2405.09323, whose content I did not check. | Exa fetches |
| L2 | "D = 2^J … (text after (2.29))" in BP | **CITATION PROBLEM** (minor) | The D = 2^J sentence follows BP eq. (2.37), not (2.29). | `sources` (T3's BP text, line 272) |

## Details

### Claim 1 — BP flatness, completeness, genericity, power

**Completeness.**
- I rebuilt T3's J = 4 charges (seed 11, R = 1500) and listed every orthant point with |Λ| < 0.0226 by a *non-MITM* algorithm. It loops over n₁, vectorises over (n₂, n₃), and solves n₄ analytically, with a ±1 guard against floor errors.
- All 18 counts N[0,ε] and N[−ε,0] (2, 10, 24, …, 22550 | 5, 15, …, 22579) are identical to T3's.
- One trap: de-duplicating the levels by float64 value loses about 3 of 45 129 genuine levels, because distinct n-vectors can round to the same float64 at Q ≈ 2.25×10⁶. T3 counts points, not unique floats, so T3 is unaffected.
- The MITM is equal to brute force for J = 3–5.

**Genericity.**
- The charges are q = 1 + U[0,1). As float64 they are dyadic rationals with common denominator 2⁻⁵². That commensurability only matters at Λ ~ 10⁻¹⁶, which is far below every mean spacing tested (≥ 3×10⁻¹¹).
- Random reals satisfy every Diophantine condition almost surely, so the charges are generic.

**Flatness with new seeds** (J = 4/6/8/10/12; R = 800/200/70/40/26; 2×10⁵ levels each).
- Fitted p = 1.009, 1.001, 0.993, 1.012, 1.000.
- On the first seeds χ²/dof = 1.53, 1.11, 1.08, 1.54, 0.89, with two at p ≈ 0.01. These were fluctuations: six further seeds each at J = 4, 8, 10 give pooled χ²/dof = 1.000, 0.877, 0.968 over 264 dof.
- T3's quoted 0.89–0.95 are therefore representative.

**Power.**
- The same statistic, applied to near-commensurate charges c_i = 1 + δu_i (J = 8, S = 400.37), gives:

  | δS | χ²/dof | fitted p |
  |---|---|---|
  | 120 | 0.91 | 1.003 |
  | 12 | 171 | 0.997 |
  | 4 | 3.7×10³ | 0.991 |
  | 1.2 | 3.6×10⁵ | 1.33 |

- So the test would certainly have detected the known near-commensurate structure.
- Nuance: at δS = 12 and 4 the failure comes from a %-level offset and modulation of the density relative to the volume formula. The density near Λ = 0 is 3.5 % low at δS = 12 and 5 % high at δS = 4. Local flatness (p ≈ 1) and Poisson spacings survive at those δS. Genuine local non-flatness and gaps appear only at δS ≲ 1, as T3 says.

### Claim 2 — identity, face formula, and "BP underestimates"

**Identity.**
- A point with support T has 2^|T| sign images, and 2^(J−|T|) subsets U contain T. So it is counted exactly 2^J times on the right-hand side.
- Brute force (J = 2–5, 12 random cases) is exact.
- The smooth version is the saddle-point evaluation of Π_i ½(θ(t_i) + 1).

**Convention question.**
- BP (2.24) says, verbatim, "let the typical degeneracy be D … the volume per D grid points must then be less than the volume of the shell". They then set D = 2^J (paragraph after eq. 2.37).
- Points with n_i = 0 have fewer than 2^J images, so D = 2^J over-counts the degeneracy whenever faces matter. The underestimate of distinct values is real; it is not a convention mismatch.

**Exact large-J test.** Coefficients of f(x)^J with f = Σ_{n≥0} x^{n²}, averaged over ±8 integers around the shell, in BP units (|2λ_bare| = 2):

| J | log₁₀ q | n_max | BP (D = 2^J) | face formula | exact lattice | exact/face | exact/BP |
|---|---|---|---|---|---|---|---|
| 100 | −1.60 | 56 | 107.02 | 109.93 | 109.93 | 10^0.00 | 10^2.90 |
| 100 | −1.43 | 38 | 90.02 | 94.22 | 94.22 | 10^0.00 | 10^4.20 |
| 100 | −1.73 | 76 | 120.02 | 122.20 | 122.20 | 10^0.00 | 10^2.18 |
| 300 | −0.87 | 10.5 | 29.84 | 90.50 | **92.18** | **10^1.69** | **10^62.35** |
| 300 | −1.171 | 21 | 120.14 | 155.20 | 155.47 | 10^0.27 | 10^35.33 |

(The last four numeric columns before the ratios are log₁₀ of the density.)

- **2b** holds exactly at J = 100. But the context is missing: BP's own J = 100 number is the "D small" case (R2 verified this earlier). With D small, BP *over*state the distinct density by 2^100/10^2.9 ≈ 10^27.2.
- **2c fails as stated.** J = 300 appears nowhere in hep-th/0004134; I grepped the full text. At n_max < √J, BP explicitly say (2.24)–(2.26) "break down … the flux … should simply be ignored". T3's own continuum face formula also misses the true count by a factor of 49 there.
- T3's §2(c) does call this "an extrapolation" and "the breakdown regime BP flagged". §9 then reports it as "off … by 10^2–10^60 in their large-J examples", which misattributes the example and overstates the precision.

### Claims 3–5
- See the table.
- Repeated charges: J = 12 with three charges, each used 4 times (T3's seed 400).
  - Distinct values in |Λ| < 0.5: 315 at R = 8 and 1059 at R = 11. These are reproduced both by my MITM (keyed on (m_a) = sums of the four squares) and by a direct simplex count #{m ∈ ℤ³≥0 : Σc_a m_a ∈ [S−1, S+1)}.
  - The mean degeneracy is 862 and 3280.
  - Agreement with T3 is exact.

### Claim 6 — factor of 2
- **AD hep-th/0307049 §3.2** (text checked): "We can furthermore take a₁ > 0, taking into account a₁ < 0 by multiplying this result by 2 … N_vac(L) = 2σ(L)". The continuum volume π²L/3 matches *that* count.
- **DD hep-th/0404116 §5.1** prescribes "a₂ = 0, 0 ≤ b₁ < a₁ and a₁b₂ ≤ L". I implemented that recipe, mapped τ into F, and get 18 604 = Σσ representatives, of which **120** land at τ = 2i. DD state 240, which is #{A : 0 < 2|A|² ≤ 150} with A and −A counted separately. So DD's figures also count N and −N.
- Brute force over all fluxes with τ in a half-open F gives 192, 712, 2748, 6120, 10752 at L = 10, 20, 40, 60, 80. These equal 2Σσ + 2·#orbits(i) exactly.
- No flux ever stabilises at ρ = e^{2πi/3}, because τ = −B̄/Ā ∈ ℚ(i). T3's "4·#(orbits at ρ)" term is therefore identically 0.
- DD fig. 7 counts with N and −N separate:
  - R = 0.15: 632. R = 0.30: 2592. These match T3.
  - R = 0.20: strict 1156, closed 1252. R = 0.25: strict 1712, closed 1784.
  - T3's 1202 and 1734 lie between, decided by float rounding of points exactly on the circle. These two comparisons are ill-posed at about ±4 %. This is harmless to the conclusion.
- The statement that −1 ∈ SL(2,ℤ) acts as (−1)^{F_L}Ω is standard; T3's CONJECTURE label is over-cautious, not wrong.

### Claim 7 — MSW + Planck 2018
- **MSW model.** I re-derived it from the physics of MSW §2: threshold δ_c = (729ρ_V/500ρ̄)^{1/3}, infall fraction sδ/(δ_c + sδ), Gaussian N(δ). I did not use their closed forms.
- **Eq. (25):** 2.5788 / 3.0103 / 3.2235 (MSW 2.579 / 3.010 / 3.224).
- **Quantiles at s = 1:** 0.05897 / 1.4852 / 11.113 (MSW 0.0590 / 1.49 / 11.1).
- **Table 1** (BBKS, eq. 43 Γ, COBE δ_H with a = −0.95 and b = −0.169, f(λ₀) from eq. 48, z_rec = 1000): all 18 σ entries agree to the printed digits. At λ₀ = 0.7 I get P = 0.0530/0.1212, against MSW's 0.054/0.12.
- **Planck update.**
  - Inputs equal Planck 2018 Table 2, column TT,TE,EE+lowE+lensing (ω_b, ω_c, H₀, n_s, ln10¹⁰A_s, τ, Ω_m, Ω_Λ, σ₈; checked against T3's copy of 1807.06209).
  - My own CAMB 2.0.4 run gives σ₈ = 0.8109.
  - σ₀(R_G) for CDM+baryons with a Gaussian filter: 2.7920 / 2.1055 / 1.5100 / 1.2079. These are T3's values.
  - g₀ = 0.78799, ρ_obs/ρ_m0 = 2.1716, and P(s=1) = 0.0427 / 0.0857 / 0.1842 / 0.2929, with medians 30.4 / 13.0 / 4.8 / 2.5.
  - EH98 no-wiggle and BBKS, both normalised to σ₈, give P(1 Mpc) = 0.0887 and 0.0901.
  - Prior tilt ρ^α with α = ±0.05 gives 0.102 / 0.072, as T3 states.
  - MSW's L* filter formula with Planck values gives R_G = 0.9985 Mpc.
- **Minor points.**
  - T3's "M_f" column is the *total-matter* mass (2π)^{3/2}ρ_m0R_G³, whereas MSW define M_f with ρ_B0 (a factor of 6.4 smaller).
  - Earlier MSW-type updates exist and are uncited. Pogosian–Vilenkin (astro-ph/0611573, seen via ar5iv in search, not fetched from arxiv.org) quote Weinberg 2005 (WMAP-1, R_G = 2 Mpc) as P = 15.6 %, and report a WMAP-3 typicality of 22 %. T3's 18.4 % at R_G = 2 Mpc is consistent with those.
  - T3's own framing ("like-for-like update", "expected results confirmed") is honest.

### Claim 8 — minimal slopes
- **Reduction.** A sympy 10d → 4d minisuperspace reduction of ∫√−G e^{−2φ}(R + 4(∂φ)²), with τ = e^{−φ}ρ^{3/2} and an Einstein-frame rescaling, gives
  L ∝ −6a_E ȧ_E²/N_E + a_E³[(3/2)(lnρ)˙² + 2(lnτ)˙²]/N_E.
  That is K = (3/4)(∂lnρ)² + (∂lnτ)² relative to ½R, which is HKTT eq. (2.5).
- **Exponent vectors.** These come from explicit substitution and agree with HKTT (3.3), (3.4), (4.1), (4.2), (4.4) and with STT 2306.07327 §IV.1, where γ_δ̃ = d/√(d−2) − ½χ√(d−2).
- **Hull minima** (exact):

  | source set | min² |
  |---|---|
  | IIB, T3's set | 18/7 (R6–D9, t = 1/7) |
  | IIA | 50/19 (R6–F0) |
  | +Q | 2 |
  | IIB without D9 | **8/3**, i.e. R6 alone |
  | IIB without D9, plus a non-SUSY 10d one-loop Λ | 50/19 |

- **Checks.**
  - The pointwise random test over 2×10⁴ points shows no violation.
  - The ODE late-time ε is 9/7, 25/19, 1 and 4/3.
- **Physics caveat (beyond T3's "formal" remark).**
  - In IIB proper the C₁₀ tadpole forces zero net D9 charge. Positive net D9 tension then needs anti-D9s. With D9s that means a tachyon. With O9⁺ it means a *non-supersymmetric* model, such as the tachyon-free Sugimoto USp(32) string with its positive e^{−φ} tadpole.
  - So "IIB √(18/7)" is not attained within type IIB. Its physical realisation is a non-SUSY open-string model.
  - For IIB proper the sharp all-positive value is √(8/3) = 1.633.
  - 1.604 is still a correct *lower bound*, so T3's conclusion survives.
- **Novelty.** The pointwise statement is the convex-hull lemma (STT 2306.07327 "convex-hull criterion"; CRV 2209.11821), and T3 itself says "explicit instances of STT/CRV". The novelty bullet in §9 ("closed forms … with the pointwise bound") should drop the pointwise bound.

### Claim 9 — mixed signs
- **Two-term case.** The transverse potential f = Ae^{−ax} − Be^{−bx} at a critical point satisfies f'' = −ab·f. So V > 0 forces a maximum when ab > 0, and the critical point has V < 0 when ab < 0.
  - STT 2303.03418 (fetched) states this: "an accelerating universe involving a rolling dilaton minimally requires at least three terms in the potential, not all positive".
  - This part is VERIFIED-KNOWN.
- **Collinear (+,−,+) triples.** I enumerated every line that carries ≥ 3 source vectors.
  - The minimum |f| is √(2/3), on R6(+)/O4(−)/F4(+) and on Q/O8/F2. Both feet lie outside the segment, so the transverse minimum is only local (metastable).
  - Numeric example: e^{−x} − 1.9e^{−1.5x} + e^{−2x} has a local minimum at x = 0.446 with V = 0.077 > 0.
  - ODE from that minimum: ε = 1/3 and dφ₁₀/dN = +1, so the flow runs to strong coupling.
  - Restricted to weak coupling and non-shrinking string-frame volume (f_ρ ≥ 0 and f_τ ≥ √3 f_ρ), the minimum is 50/19 (R6/O8/F0) for geometric sources and 2 with Q. Both agree with T3.
- **Lines through the origin** carrying a (+,−,+) pattern: R6/O3/(F5 or Casimir) and F1/O9/Q. On these, V depends only on ρτ², which gives the GKP-type no-scale flat direction (the known O3/F5 exception, e.g. Wrase–Zagermann 1003.0029, seen in search highlights only).

### Claim 10 — DESI phrasing
- 2504.04226 (fetched abstract and HTML) gives λ = 0.698 +0.173/−0.202, or 0.722 +0.182/−0.208 with Ω_K free. The data are DESI DR2 + DESY5 + Planck *distance priors*.
- T3's numbers are correct arithmetic: 4.14, 5.24, 5.34, 3.80 and 0.68σ.
- But a 68 % half-width is extrapolated to the 4–5σ tail without the posterior. Two other readings:
  - In λ² (≈ (1+w₀)-like, σ_λ² ≈ 2λσ_λ = 0.242) the √2 point sits at 6.3σ.
  - With Pantheon+ instead of DESY5, DR1 gave λ = 0.48 +0.28 (2405.17396), i.e. about 3.3σ.
- The direction of the claim is robust; the specific σ values are not.
- The report text says "one-sided, using the DESI upper error"; only T3's .out file adds that this is a Gaussian reading of an asymmetric posterior.
- The headline labels the σ values "ESTABLISHED (…) + OUR ATTEMPT". They should be OUR ATTEMPT, with the Gaussian-extrapolation caveat stated in the report.

## Label check
- **CONJECTURE presented as ESTABLISHED:** none found.
- **Mislabels:**
  - (a) The DESI σ distances are tagged ESTABLISHED + OUR ATTEMPT (claim 10).
  - (b) The −1 = (−1)^{F_L}Ω identification is tagged CONJECTURE although it is standard (over-cautious, harmless).
- **Overclaims in §9 "Novelty":**
  - "10^2–10^60 in their [BP's] large-J examples" (2c, FAILED).
  - "with the pointwise bound" (known lemma).
- **Literature gaps, context only, journal references not fetched this session:**
  - Davenport–Lewis gap theorem for irrational positive forms (Bentkus–Götze 1999; Götze 2004, d ≥ 5). This already makes "no gaps at fixed shell width" a theorem for J ≥ 5.
  - Eskin–Margulis–Mozes (2005) for Poisson pair correlation of rectangular 2-tori.
  - Weinberg 2005 and Pogosian–Vilenkin 2007 for earlier MSW updates.
  - Sugimoto/brane-SUSY-breaking for the D9-like positive tadpole.

## Files
- `RT3_code/rt3_bp_core.py`: my BP enumeration library (MITM, exact Fraction checks, face formula).
- `rt3_bp1_identity_completeness.py`: identity brute force; independent J = 4 completeness check against T3's table.
- `rt3_bp2_flatness_power.py`, `rt3_bp2b_pulls.py`, `rt3_bp2c_seeds.py`: flatness, power, near-commensurate charges, min-z statistics, extra seeds.
- `rt3_bp3_largeJ_faces.py`: exact large-J lattice counts.
- `rt3_bp4_repeated.py`: J_eff check.
- `rt3_fac2.py`: AD/DD brute force.
- `rt3_msw.py`: MSW re-derivation, Table 1, Planck/CAMB/EH/BBKS.
- `rt3_slopes.py`, `rt3_slopes_ode.py`: reduction, exponent vectors, hulls, mixed-sign lines, ODE checks.
- Each script has a matching `.out` file.
- `sources/`: STT 2306.07327v2 and Akrami et al. 2504.04226v2 full texts, fetched via Exa.
