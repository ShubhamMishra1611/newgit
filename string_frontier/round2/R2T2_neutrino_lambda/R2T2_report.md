# R2-T2: Does consistency tie Λ to Standard Model data? (Standard Model on a circle + non-SUSY AdS conjecture)

Done by the lead agent. Scripts: `r2t2_casimir.py` (core + reproduction, output `r2t2_casimir_partA_output.txt`) and `r2t2_update.py` (update, output `r2t2_update_output.txt`).

## 0. Bottom line

- **Reproduced.** The Ibáñez–Martín-Lozano–Valenzuela (IMV) 3D results are reproduced to ≤ 0.2% for the AdS₃-onset masses, and within 2% for the dS₃-onset masses.
- **Updated.** With NuFIT 6.0 and Planck 2018, the non-SUSY AdS conjecture applied to the Standard Model on a circle requires:
  - Dirac neutrinos;
  - for normal ordering: m₁ < 7.64 meV (7.54–7.75 meV over the 3σ oscillation ranges), i.e. Σm_ν < 0.0699 eV;
  - for inverted ordering: m₃ < 2.53 meV, i.e. Σm_ν ∈ [0.0989, 0.1016] eV.
- **Majorana neutrinos** would need ρ_Λ ≥ 9.9 ρ_obs (normal ordering) or ≥ 389 ρ_obs (inverted ordering).
- **This is not a selection principle for Λ.**
  - It gives a lower bound ρ_Λ ≥ ρ_min(m_lightest), which vanishes as m_lightest → 0.
  - ρ_obs sits at that bound only if m₁ ≈ 7.6 meV.
  - It also rests on a conjecture plus an unproven stability assumption.

## 1. Setup (ESTABLISHED)

- **Radion potential.** Compactify the Standard Model coupled to gravity on a circle of radius R (hep-th/0703067; IMV 1706.05392 eqs. 3.1–3.12). In the 3D Einstein frame:
  - V(R) = (2πr³/R²)[ρ_Λ + Σ_i (−1)^{s_i} n_i ρ(R; m_i)].
  - Casimir density per degree of freedom (IMV eq. 3.10): ρ(R; m) = Σ_k 2m⁴/(2π)² · K₂(2πRmk)/(2πRmk)², which tends to π²/(90(2πR)⁴) as m → 0.
  - Bosons enter with −, fermions with periodic boundary conditions with +.
- **Field content.**
  - Massless bosons: photon + graviton, 4 degrees of freedom.
  - Neutrinos: 2 (Majorana) or 4 (Dirac) degrees of freedom each.
  - Heavier states are exponentially negligible at R ~ 1/m_ν.
- **Conjecture (CONJECTURE).** Non-SUSY AdS conjecture (Ooguri–Vafa 1610.01533): no stable non-SUSY AdS vacuum.
  - It is applied here as in IMV, to the 3D vacua.
  - IMV §2 discuss the stability caveat: the 3D AdS vacua could inherit the 4D dS instability.

## 2. Checks and reproductions (OUR ATTEMPT, `r2t2_casimir.py`)

- **Casimir density.** The massless limit reproduces π²/90 exactly. IMV's small-mR series (eq. 3.11) agrees with the Bessel sum to 6×10⁻⁵ at mR = 0.01, degrading to 5% at mR = 0.1, as expected for a truncated series.
- **Maximum with bosons and Λ only.** R_max = 75.47 eV⁻¹ numerically, against IMV eq. 3.6/3.7: (120π²Λ)^{−1/4} = 75.49 eV⁻¹ (7.55×10¹⁰ GeV⁻¹).
- **IMV inputs.** Λ = 2.6×10⁻¹¹ eV⁴ and PDG-2016-type splittings. The exact values IMV used were lost in text extraction; we assumed Δm²₂₁ = 7.53×10⁻⁵ and |Δm²₃₂| = 2.44 (NO) / 2.51 (IO) ×10⁻³ eV².

| Case | ours: dS onset | ours: AdS onset | IMV (sec. 3.4, Table 1) | ADNV (hep-th/0703067) |
|---|---|---|---|---|
| Majorana NO / IO | AdS at m_light = 0 | AdS at m_light = 0 | AdS always | AdS |
| Dirac NO | 6.585 meV | 7.724 meV | 6.7 / 7.73 meV | 7.1 / 8.3 meV |
| Dirac IO | 2.082 meV | 2.562 meV | 2.1 / 2.56 meV | 2.5 / 3.1 meV |

- The ADNV numbers differ because of their 2007 inputs (Λ and splittings); see the Λ scaling in §3.
- **Wording slip in IMV §3.5.** IMV's text calls m_ν₁ > 7.7 meV a "lower bound". Their own Table 1 and abstract make clear the conjecture gives an *upper* bound. Our computation confirms the upper-bound reading.

## 3. Update (OUR ATTEMPT, `r2t2_update.py`)

**Inputs**
- NuFIT 6.0 (2410.05380, Table 1, IC24 + SK-atm):
  - Δm²₂₁ = 7.49 (3σ: 6.92–8.05) × 10⁻⁵ eV²;
  - Δm²₃₁(NO) = 2.513 (2.451–2.578) × 10⁻³ eV²;
  - Δm²₃₂(IO) = −2.484 (−2.547 to −2.421) × 10⁻³ eV².
- Planck 2018 (1807.06209, Table 2): h = 0.6736 and Ω_Λ = 0.6847, giving ρ_Λ = 2.515×10⁻¹¹ eV⁴ = (2.239 meV)⁴ = 7.15×10⁻¹²¹ M_P⁴ (reduced Planck mass).

**Results**

| | NO (Dirac) | IO (Dirac) |
|---|---|---|
| dS₃ minimum appears at m_light > | 6.514 meV | 2.060 meV |
| AdS₃ appears at m_light > | **7.644 meV** (3σ corners: 7.538–7.754) | **2.534 meV** (2.507–2.562) |
| Allowed Σm_ν | 0.0588–0.0699 eV (3σ upper end 0.0692–0.0706) | 0.0989–0.1016 eV (upper end 0.1003–0.1029) |

- **Majorana, m_light = 0.** AdS₃ disappears only for ρ_Λ > 9.92 ρ_obs (NO) or > 389 ρ_obs (IO). IMV quoted "≈7×" for NO with older inputs.
- **Λ scaling (Dirac NO).** The onset follows m₁* ≈ (3.2–3.8) ρ_Λ^{1/4} over ρ_Λ ∈ [0.5, 5] ρ_obs.
- **Lower bound ρ_min(m₁) on ρ_Λ (Dirac NO)** needed to avoid AdS₃:

| m₁ | ρ_min |
|---|---|
| 1 meV | 7.3×10⁻³ ρ_obs |
| 4 meV | 0.18 ρ_obs |
| 7.7 meV | 1.02 ρ_obs |
| 23 meV (DESI m_l bound) | 37 ρ_obs |

  The ratio ρ_min^{1/4}/m₁ runs from 0.65 at 1 meV to 0.24 at 23 meV.

**Comparison with DESI DR2 (2503.14744; ESTABLISHED)**
- ΛCDM: Σm_ν < 0.0642 eV (95%, eq. 15) and m_l < 0.023 eV. w₀wₐ: Σm_ν < 0.163 eV.
- **NO.** The conjecture's bound (Σ < 0.0699 eV) is weaker than DESI-ΛCDM (0.0642) but much stronger than DESI-w₀wₐ (0.163). In m_l, 7.6 meV is about 3× stronger than DESI's 23 meV.
- **IO.** The conjecture confines Σm_ν to within 3% of its oscillation minimum, a sharp prediction. DESI-ΛCDM already disfavours IO; DESI-w₀wₐ does not.

## 4. Falsifiable consequences of "Standard Model + this conjecture" (OUR ATTEMPT; the logic is in IMV/ADNV, the numbers are updated here)

1. **Neutrinoless double-beta decay.** A detection (Majorana) falsifies the conjecture plus the minimal Standard Model, unless light BSM fermions exist; with Majorana neutrinos the vacuum energy would have to be ≥ 9.9 ρ_obs.
2. **Normal ordering.** A cosmological measurement Σm_ν > 0.070 eV (3σ envelope 0.0706 eV) falsifies it.
3. **Inverted ordering.** IO combined with Σm_ν > 0.103 eV falsifies it.
- **Caveat.** These hold only if the 3D AdS vacua are stable and the conjecture applies to them; IMV §2 and later work discuss this. Wilson lines and the QCD-scale potential at small R (IMV note added) are also not included.
- **The 2D torus** gives stronger bounds: m₁ < 4.12 meV (NO) and m₃ < 1.0 meV (IO), per IMV §4 with their inputs, ESTABLISHED. We did **not** recompute the 2D case. With NuFIT 6.0, m₁ = 4.1 meV corresponds to Σm_ν ≈ 0.064 eV, numerically the DESI-ΛCDM bound. This is a coincidence of two independent numbers and nothing more.

## 5. Dark energy that evolves (literature, ESTABLISHED)
- If dark energy is runaway quintessence (DESI's w₀wₐ preference), there is no 4D vacuum, and the 3D analysis changes.
  - Gonzalo–Ibáñez–Valenzuela (2109.10961) treat a simple quintessence potential.
  - They also show that the AdS distance conjecture and the Lüst–Palti–Vafa dS-type criterion, which need only local IR data, give similar constraints.
- We did not redo that analysis.
- In our potential, ρ_Λ → 0 with any m₁ > 0 produces an AdS₃ minimum: F(R) → −4ρ₀ < 0 at large R and is positive at small R. So the status of the asymptotic Λ → 0 endpoint under the conjecture matters, and it is open.

## 6. Novelty (honest)
- **Already published.** The mechanism, the Dirac requirement and the bound structure (ADNV 2007; IMV 2017; Hamada–Shiu 1707.06326; Gonzalo–Ibáñez–Valenzuela 2021).
- **New here, as far as our search shows.** Only the update to NuFIT 6.0, Planck 2018 and DESI DR2 with 3σ propagation (Σ bounds 0.0699 eV NO and 0.0989–0.1016 eV IO; Majorana factor 9.9×), and the explicit ρ_min(m₁) table.
- **Mission (2).** This is a consistency *relation* between Λ and the neutrino sector, not a principle that selects Λ ≈ 10⁻¹²².
