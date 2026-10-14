## 2. Frontier map (Phase 1)

Five parallel agents mapped the fronts. The full maps, with verified arXiv IDs and equation numbers where they were read, are in `phase1/`. Each agent verified 25–45 IDs by fetching the abstract page.

### A. Non-perturbative formulations (BFSS, IKKT, closed SFT, M-theory): `phase1/A_nonperturbative.md`

**State of the art**
- **ESTABLISHED.**
  - BFSS at finite N is DLCQ M-theory (Susskind; Seiberg; Sen).
  - Eikonal two-graviton checks hold: 15/16·v⁴/r⁷ at one loop and 225/64·v⁶/r¹⁴ at two loops (hep-th/9706072 eq. 13).
  - Lattice simulations reproduce E/N² = 7.41 T^{14/5}, but only in the IIA/'t Hooft regime.
  - Closed SFT is perturbatively complete and unitary. Hyperbolic vertices plus a 2024 recursion fix all interactions from the cubic one.
  - M-theory has a non-perturbative definition only in AdS or DLCQ settings, e.g. the ABJM Airy-function partition function.
- **CONJECTURE.**
  - The flat-space large-N BFSS limit at non-zero longitudinal momentum transfer.
  - Lorentzian IKKT emergent (3+1)-d expansion: the 2011 result was a Pauli-matrix artifact, and the 2026 version needs mass and anisotropy deformations.
  - Steinacker's emergent Einstein–Hilbert action.
  - SFT background independence beyond local.

**Sharpest open questions**
1. Does the large-N BFSS 2→2 amplitude with longitudinal momentum exchanged equal the 11d supergravity amplitude?
2. Does Lorentzian IKKT have a regulator-free phase with three non-Pauli-like growing directions?
3. Does hyperbolic closed SFT have a finite-distance, zero-action vacuum with empty cohomology?

### B. Holography beyond AdS (dS/CFT, DSSYK, celestial/Carrollian): `phase1/B_holography_beyond_AdS.md`

**State of the art**
- **ESTABLISHED.**
  - The exact one-loop sphere partition functions (2009.12464).
  - The Type II₁ observer algebra in dS (2206.10780).
  - Tree-level w₁₊∞ (2105.14346).
  - Self-dual top-down flat-space duals (2208.14233).
- **CONJECTURE.**
  - dS/CFT has only non-unitary or partial examples.
  - In dS₃, a matrix-integral match to pure gravity (2501.01486).
  - DSSYK ↔ dS has three competing dictionaries.
  - Whether dS entropy is a positive state count is contested (2412.14014 vs 2511.01400).

**Sharpest open questions**
1. Can a microscopic model reproduce the one-loop dS₄ entropy as a real, positive count?
2. Does doubled DSSYK match the dS₃ propagator at O(λ)?
3. Is there a top-down flat-space dual containing non-self-dual Einstein gravity?

### C. Vacuum selection (landscape counting, swampland): `phase1/C_vacuum_selection.md`

**Theorems, each with its scope**
- Harlow–Ooguri: no global symmetries, completeness and compactness. Proven within AdS/CFT only (1810.05337/05338).
- Finiteness of F-theory flux vacua at bounded tadpole: complex structure on one fixed CY₄ (2112.06995).
- Finding a Bousso–Polchinski vacuum with a target Λ is NP-complete (hep-th/0602072).
- 10d string universality (1006.1352).

**Conjectures**
- The full swampland program: WGC, Distance and Emergent String, dS/refined-dS/TCC, tadpole, cobordism, festina lente, and the dark dimension.
- The DESI DR2 exponential-quintessence fit gives λ = 0.698 (+0.17/−0.20) (2504.04226).

**Sharpest open questions**
1. Does controlled metastable dS exist, outliving the TCC bound?
2. Does any controlled asymptotic region have slope < √2 in d = 4?
3. Is the discretuum's Λ distribution flat near 0, the flat prior Weinberg needs?

### D. dS vacua (KKLT, LVS, criticisms): `phase1/D_dS_vacua.md`

**State of the art**
- **ESTABLISHED (leading order).** Explicit KKLT dS vacua (2406.13751: 5 compactifications, 30 vacua) and LVS constructions (2010.15903, 2206.12889).
- **Contested on four fronts.**
  - Anti-D3/KPV control at g_sM ≈ 1.
  - Conifold destabilization: g_sM² ≳ 46 (1809.06861), disputed by Lüst–Randall (2206.04708).
  - Odd flux quanta in CY orientifolds.
  - Warping corrections. 2512.17995 calls KKLT "effectively uncontrollable with current methods" and finds LVS corrections volume-suppressed.

**Sharpest open questions**
1. What are the true anti-D3 and conifold thresholds at the candidates' parameters?
2. Are odd fluxes allowed?
3. Is there an LVS geometry where every control inequality holds within a known tadpole?

### E. S-matrix bootstrap: `phase1/E_bootstrap.md`

**State of the art**
- **ESTABLISHED.** Every "string = unique solution" theorem has one UV-shaped input not derived from IR consistency:
  - level truncation (2406.02665, 2408.03362);
  - ultrasoftness plus minimal zeros (2508.09246);
  - monodromy or splitting (2605.11084);
  - maximal SUSY plus multi-point factorization (2601.11705, 2610.02124);
  - for gravity, finite spin at the first massive level (2607.14230).
- **ESTABLISHED (numerical).** String theory sits inside the allowed regions, not at a kink: 0 ≤ g₀M⁶/8πG ≤ 3 vs 2ζ(3) ≈ 2.40 (2102.08951 eqs. 3.41–3.42), and α_min ≈ 0.13 ± 0.02 vs 0.1389 (2102.02847 eq. 10).

**Sharpest open questions**
1. Can the infinite-spin-tower corner in gravity be excluded without assuming finite spin?
2. Is α_min = 0.138949 exactly?
3. Do the surviving 4-point deformations extend to consistent 5- and 6-point amplitudes?
