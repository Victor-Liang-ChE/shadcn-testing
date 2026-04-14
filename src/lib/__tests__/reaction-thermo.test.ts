/**
 * Reaction Thermochemistry Tests
 * Run: npx tsx src/lib/__tests__/reaction-thermo.test.ts
 *
 * Tests the core thermochemistry calculations used in the reaction-thermochemistry page.
 * Uses DIPPR Aly-Lee (eqno=107) coefficients sourced directly from the Supabase DB
 * (queried once and hardcoded here for reproducibility / offline use).
 *
 * KEY NOTE — "Frozen Composition" vs. Equilibrium AFT:
 * ─────────────────────────────────────────────────────
 * This page computes the FROZEN-COMPOSITION adiabatic flame temperature (AFT):
 * no dissociation of products is assumed, even at extreme temperatures.
 * For combustion reactions, FROZEN T_ad is significantly higher than the
 * EQUILIBRIUM T_ad (which accounts for CO₂ ⇌ CO+½O₂, H₂O ⇌ H₂+½O₂, etc.):
 *
 *   CH4 + 2O2 → CO2 + 2H2O
 *     Frozen T_ad (code):       ~5 100 K    (this code, extrapolated Aly-Lee)
 *     Equilibrium T_ad (NIST):  ~2 800–3 050 K
 *
 * If you want the equilibrium value you need minimization of Gibbs energy or a
 * NASA-polynomial equilibrium solver — this tool does not do that.
 */
import assert from 'assert';

// ─── Assertion helpers ─────────────────────────────────────────────────────────
let passed = 0; let failed = 0;
function check(desc: string, ok: boolean, detail = '') {
  if (ok) { console.log(`  ✓ ${desc}`); passed++; }
  else     { console.error(`  ✗ ${desc}${detail ? ': ' + detail : ''}`); failed++; }
}
function approx(desc: string, actual: number, expected: number, tolPct = 1) {
  const pct = Math.abs(actual - expected) / Math.abs(expected) * 100;
  check(desc + ` (got ${actual.toFixed(3)}, exp ${expected.toFixed(3)}, err ${pct.toFixed(2)}%)`, pct <= tolPct);
}
function within(desc: string, actual: number, lo: number, hi: number) {
  check(desc + ` (got ${actual.toFixed(1)}, range [${lo},${hi}])`, actual >= lo && actual <= hi);
}

// ─── DIPPR Aly-Lee (eqno 107) implementation ──────────────────────────────────
// Matches dipprCpDeltaH in page.tsx exactly.
// Cp [J/kmol/K] = A + B*(C/T/sinh(C/T))^2 + D*(E/T/cosh(E/T))^2
// H(T) [J/kmol]  = A*T + B*C*coth(C/T) - D*E*tanh(E/T)
function alyLeeH(A: number, B: number, C: number, D: number, E: number, T: number): number {
  const bTerm = (B !== 0 && C !== 0) ? B * C * (1 / Math.tanh(C / T)) : 0;
  const dTerm = (D !== 0 && E !== 0) ? D * E * Math.tanh(E / T)       : 0;
  return A * T + bTerm - dTerm; // J/kmol
}
function alyLeeDeltaH(A: number, B: number, C: number, D: number, E: number, T1: number, T2: number): number {
  return (alyLeeH(A, B, C, D, E, T2) - alyLeeH(A, B, C, D, E, T1)) / 1e6; // J/kmol → kJ/mol
}

// ─── NIST DB coefficients (verified against DB 2026-04-10) ────────────────────
// All gas-phase Aly-Lee coefficients (eqno=107, range 50–1000 K)
//
// Source: nist_shomate table  SELECT nist_id,a,b,c,d,e FROM nist_shomate WHERE phase='gas'
//         nist_compounds table: hf_gas_298_kjmol
//
// Species    formula   nist_id     hf_gas_298_kjmol (kJ/mol)
// CH4        C1        C74828      -74.6107
// CO2        CO2       C124389     -393.51
// H2O        H2O       C7732185    -241.818
// CO         CO        C630080     -110.546
// O2         O2        C7782447    null (element, 0 by convention)
// H2         H2        C1333740    null (element, 0 by convention)
// C3H8       C3H8      C74986      -104.455

interface AlyLeeCoeffs { A: number; B: number; C: number; D: number; E: number; }

const SH: Record<string, AlyLeeCoeffs> = {
  CH4: { A: 33239.6, B: 76765.2, C: 1992.37, D: 38824.5, E: 969.49   },
  CO2: { A: 29180.2, B: 32004.5, C: 1148.34, D: 19867.8, E: 508.32   },
  H2O: { A: 33274.8, B: 23863.8, C: 2359.34, D:  7351.22, E: 1075.62  },
  CO:  { A: 29103.3, B:  7607.94, C: 2874.37,  D:  8105.58, E: 1522.63  },
  O2:  { A: 29111.9, B:  8863.2, C: 2187.19, D:  8196.58, E: 1107.10  },
  H2:  { A: 41227.0, B:-11437.5, C:  396.298, D:-31098.4, E:  147.084  },
  C3H8:{ A: 34102.5, B:166702.0, C:  850.063, D: 48387.4, E:  276.694  },
};

const HF: Record<string, number> = {
  CH4:  -74.6107,
  CO2:  -393.51,
  H2O:  -241.818,
  CO:   -110.546,
  O2:    0,
  H2:    0,
  C3H8: -104.455,
};

function dH_species(sp: string, T1: number, T2: number): number {
  const c = SH[sp]; if (!c) return 0;
  return alyLeeDeltaH(c.A, c.B, c.C, c.D, c.E, T1, T2);
}

// ─── 1. Sensible heat integrals (NIST Shomate reference values) ───────────────
console.log('\n=== Sensible heat ΔH integrals vs NIST WebBook ===');
// CO2: NIST Shomate table H°(T)−H°(298.15 K) [kJ/mol]
//   500 K: 8.311   800 K: 22.811   1000 K: 33.397
approx('CO2  298→500 K',  dH_species('CO2', 298.15,  500), 8.311,  1.0);
approx('CO2  298→800 K',  dH_species('CO2', 298.15,  800), 22.811, 1.0);
approx('CO2  298→1000 K', dH_species('CO2', 298.15, 1000), 33.397, 1.5);

// H2O(g): NIST Shomate table H°(T)−H°(298.15 K) [kJ/mol]
//   500 K: 6.922   800 K: 17.991   1000 K: 26.000
approx('H2O  298→500 K',  dH_species('H2O', 298.15,  500), 6.922,  1.0);
approx('H2O  298→800 K',  dH_species('H2O', 298.15,  800), 17.991, 1.5);
approx('H2O  298→1000 K', dH_species('H2O', 298.15, 1000), 26.000, 1.5);

// CO: NIST table H°(T)−H°(298.15 K) [kJ/mol]
//   500 K: 5.930   800 K: 15.175
approx('CO   298→500 K',  dH_species('CO', 298.15,  500), 5.930,  1.0);
approx('CO   298→800 K',  dH_species('CO', 298.15,  800), 15.175, 1.5);

// H2: NIST table H°(T)−H°(298.15 K) [kJ/mol]
//   500 K: 5.882   800 K: 14.741
approx('H2   298→500 K',  dH_species('H2', 298.15,  500), 5.882,  1.0);
approx('H2   298→800 K',  dH_species('H2', 298.15,  800), 14.741, 1.5);

// ─── 2. ΔH_rxn at 298 K (standard Hess's law) ─────────────────────────────────
console.log('\n=== ΔH°_rxn at 298 K (Hess\'s law, gas-phase products) ===');
function dHrxn(products: [string, number][], reactants: [string, number][]): number {
  const s = (pairs: [string, number][]) => pairs.reduce((acc, [sp, n]) => acc + n * HF[sp], 0);
  return s(products) - s(reactants);
}

// CH4 + 2O2 → CO2 + 2H2O(g)  Literature (NIST LHV path, gas H2O): −802.3 kJ/mol
const dH_CH4 = dHrxn([['CO2',1],['H2O',2]], [['CH4',1],['O2',2]]);
approx('CH4 + 2O2 → CO2 + 2H2O(g)', dH_CH4, -802.3, 0.1);

// CO + ½O2 → CO2   Literature: −283.0 kJ/mol (NIST: −282.984)
const dH_CO = dHrxn([['CO2',1]], [['CO',1],['O2',0.5]]);
approx('CO + 0.5 O2 → CO2', dH_CO, -282.984, 0.1);

// H2 + ½O2 → H2O(g)   Literature: −241.8 kJ/mol
const dH_H2 = dHrxn([['H2O',1]], [['H2',1],['O2',0.5]]);
approx('H2 + 0.5 O2 → H2O(g)', dH_H2, -241.818, 0.1);

// C3H8 + 5O2 → 3CO2 + 4H2O(g)  Literature (LHV): −2043.3 kJ/mol
const dH_C3H8 = dHrxn([['CO2',3],['H2O',4]], [['C3H8',1],['O2',5]]);
approx('C3H8 + 5O2 → 3CO2 + 4H2O(g)', dH_C3H8, -2043.3, 0.1);

// CO + H2O(g) → CO2 + H2  (water-gas shift, gas phase)  Literature: −41.15 kJ/mol
const dH_WGS = dHrxn([['CO2',1],['H2',1]], [['CO',1],['H2O',1]]);
approx('CO + H2O → CO2 + H2 (WGS)', dH_WGS, -41.15, 0.1);

// ─── 3. Adiabatic flame temperature — iterative solver ─────────────────────────
console.log('\n=== Adiabatic flame temperature (frozen composition, no dissociation) ===');

/**
 * Solve for T_ad iteratively (bisection).
 * T_ad is the temperature where:
 *   Σ n_i * ΔH_sens,product(298→T_ad) = -ΔH_rxn(298K)
 * (all energy of reaction heats up products from 298→T_ad)
 */
function solveAdiabaticT(
  productSensHeat: (T: number) => number, // total product sensible heat sum [kJ/mol]
  dH_rxn: number,                         // kJ/mol (negative = exothermic)
  T_lo = 300, T_hi = 10000, tol = 0.5    // K
): number {
  const target = -dH_rxn; // must equal product sensible heat (exothermic → positive target)
  if (target <= 0) throw new Error('Reaction is endothermic; no adiabatic flame temp');
  let lo = T_lo, hi = T_hi;
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    if (productSensHeat(mid) < target) lo = mid; else hi = mid;
    if (hi - lo < tol) return 0.5 * (lo + hi);
  }
  return 0.5 * (lo + hi);
}

// ── Test A: CO + 0.5 O2 → CO2 ─────────────────────────────────────────────────
// Literature (frozen composition, Aly-Lee extrapolated to 6000 K):
//   T_ad_frozen ≈ 5180 K (our analytical estimate; matches DIPPR-based calculation)
//   T_ad_equilibrium ≈ 2973 K (NIST, accounting for CO₂ dissociation at high T)
{
  const T_ad_CO = solveAdiabaticT(
    T => 1 * dH_species('CO2', 298.15, T),
    dH_CO, // -282.964 kJ/mol
    300, 10000
  );
  console.log(`  CO + 0.5 O2 → CO2`);
  console.log(`    Frozen T_ad (code logic): ${T_ad_CO.toFixed(0)} K`);
  console.log(`    Equilibrium T_ad (NIST, with dissociation): ~2973 K`);
  console.log(`    Note: Difference is due to CO₂ dissociation → CO + ½O₂ at >2500 K,`);
  console.log(`          which absorbs energy and limits the real temperature.`);
  // The frozen T_ad should be well above 5000 K
  within('CO + ½O₂ → CO₂: frozen T_ad in expected range', T_ad_CO, 4800, 5500);
}

// ── Test B: H2 + 0.5 O2 → H2O(g) ─────────────────────────────────────────────
// Literature: frozen T_ad ≈ 5250 K; equilibrium T_ad ≈ 3080 K
{
  const T_ad_H2 = solveAdiabaticT(
    T => 1 * dH_species('H2O', 298.15, T),
    dH_H2,
    300, 10000
  );
  console.log(`\n  H2 + 0.5 O2 → H2O(g)`);
  console.log(`    Frozen T_ad (code logic): ${T_ad_H2.toFixed(0)} K`);
  console.log(`    Equilibrium T_ad (NIST): ~3080 K`);
  within('H₂ + ½O₂ → H₂O: frozen T_ad in expected range', T_ad_H2, 4900, 5600);
}

// ── Test C: CH4 + 2O2 → CO2 + 2H2O(g) ────────────────────────────────────────
// Literature: frozen T_ad ≈ 5100−5400 K; equilibrium T_ad ≈ 2800−3050 K
{
  const T_ad_CH4 = solveAdiabaticT(
    T => 1 * dH_species('CO2', 298.15, T) + 2 * dH_species('H2O', 298.15, T),
    dH_CH4,
    300, 10000
  );
  console.log(`\n  CH4 + 2O2 → CO2 + 2H2O(g)`);
  console.log(`    Frozen T_ad (code logic): ${T_ad_CH4.toFixed(0)} K`);
  console.log(`    Equilibrium T_ad (NIST): ~2800−3050 K`);
  console.log(`    Shomate data range: 50−1000 K. Extrapolation used above 1000 K.`);
  within('CH₄ + 2O₂: frozen T_ad in expected range', T_ad_CH4, 4500, 6000);
}

// ── Test D: C3H8 + 5O2 → 3CO2 + 4H2O(g) ─────────────────────────────────────
{
  const T_ad_C3H8 = solveAdiabaticT(
    T => 3 * dH_species('CO2', 298.15, T) + 4 * dH_species('H2O', 298.15, T),
    dH_C3H8,
    300, 10000
  );
  console.log(`\n  C3H8 + 5O2 → 3CO2 + 4H2O(g)`);
  console.log(`    Frozen T_ad (code logic): ${T_ad_C3H8.toFixed(0)} K`);
  console.log(`    Equilibrium T_ad (NIST): ~2267 K`);
  within('C₃H₈ + 5O₂: frozen T_ad in expected range', T_ad_C3H8, 4500, 6000);
}

// ── Test E: Verify Cp values at 298 K are physically reasonable ───────────────
console.log('\n=== Cp at 298.15 K [J/mol·K] vs NIST ===');
function alyLeeCp(A: number, B: number, C: number, D: number, E: number, T: number): number {
  const bTerm = B !== 0 && C !== 0 ? B * Math.pow(C / T / Math.sinh(C / T), 2) : 0;
  const dTerm = D !== 0 && E !== 0 ? D * Math.pow(E / T / Math.cosh(E / T), 2) : 0;
  return (A + bTerm + dTerm) / 1000; // J/kmol/K → J/mol/K
}
function Cp298(sp: string): number {
  const c = SH[sp]; if (!c) return 0;
  return alyLeeCp(c.A, c.B, c.C, c.D, c.E, 298.15);
}
// NIST values at 298.15 K [J/mol·K]:
approx('Cp(CH4, 298K)',  Cp298('CH4'), 35.7,  1.5);  // NIST: 35.69
approx('Cp(CO2, 298K)',  Cp298('CO2'), 37.1,  1.5);  // NIST: 37.11
approx('Cp(H2O, 298K)',  Cp298('H2O'), 33.6,  2.0);  // NIST: 33.58
approx('Cp(CO, 298K)',   Cp298('CO'),  29.1,  1.5);  // NIST: 29.14
approx('Cp(H2, 298K)',   Cp298('H2'),  28.8,  2.0);  // NIST: 28.84
approx('Cp(O2, 298K)',   Cp298('O2'),  29.4,  1.5);  // NIST: 29.38

// ─── Summary ──────────────────────────────────────────────────────────────────
console.log(`\n${'─'.repeat(60)}`);
console.log(`Results: ${passed} pass, ${failed} fail`);
if (failed > 0) process.exit(1);
