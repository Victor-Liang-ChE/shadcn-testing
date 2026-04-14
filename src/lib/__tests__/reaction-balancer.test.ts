/**
 * Tests for reaction-balancer.ts
 * Run: npx tsx src/lib/__tests__/reaction-balancer.test.ts
 */
import { parseFormula, checkBalance, autoBalance } from '../reaction-balancer';

// ─── Tiny assertion helper ────────────────────────────────────────────────────
let passed = 0; let failed = 0;
function expect(desc: string, actual: unknown, expected: unknown) {
  const ok = JSON.stringify(actual) === JSON.stringify(expected);
  if (ok) { console.log(`  ✓ ${desc}`); passed++; }
  else {
    console.error(`  ✗ ${desc}`);
    console.error(`      expected: ${JSON.stringify(expected)}`);
    console.error(`      actual:   ${JSON.stringify(actual)}`);
    failed++;
  }
}
function expectDeep(desc: string, actual: unknown, check: (v: unknown) => boolean) {
  const ok = check(actual);
  if (ok) { console.log(`  ✓ ${desc}`); passed++; }
  else { console.error(`  ✗ ${desc}: got ${JSON.stringify(actual)}`); failed++; }
}

// ─── parseFormula tests ───────────────────────────────────────────────────────
console.log('\n=== parseFormula ===');

expect('H2O',        parseFormula('H2O'),        { H: 2, O: 1 });
expect('CO2',        parseFormula('CO2'),         { C: 1, O: 2 });
expect('O2',         parseFormula('O2'),          { O: 2 });
expect('N2',         parseFormula('N2'),          { N: 2 });
expect('CH4',        parseFormula('CH4'),         { C: 1, H: 4 });
expect('C2H6',       parseFormula('C2H6'),        { C: 2, H: 6 });
expect('C6H12O6',    parseFormula('C6H12O6'),     { C: 6, H: 12, O: 6 });
expect('Fe2O3',      parseFormula('Fe2O3'),       { Fe: 2, O: 3 });
expect('NaCl',       parseFormula('NaCl'),        { Na: 1, Cl: 1 });
expect('Ca(OH)2',    parseFormula('Ca(OH)2'),     { Ca: 1, O: 2, H: 2 });
expect('CaCO3',      parseFormula('CaCO3'),       { Ca: 1, C: 1, O: 3 });
expect('Al2(SO4)3',  parseFormula('Al2(SO4)3'),   { Al: 2, S: 3, O: 12 });
expect('Mg3(PO4)2',  parseFormula('Mg3(PO4)2'),   { Mg: 3, P: 2, O: 8 });
expect('NH4NO3',     parseFormula('NH4NO3'),      { N: 2, H: 4, O: 3 });
expect('C3H8',       parseFormula('C3H8'),        { C: 3, H: 8 });
expect('C8H18',      parseFormula('C8H18'),       { C: 8, H: 18 });
expect('HCl',        parseFormula('HCl'),         { H: 1, Cl: 1 });
expect('H2SO4',      parseFormula('H2SO4'),       { H: 2, S: 1, O: 4 });
expect('KMnO4',      parseFormula('KMnO4'),       { K: 1, Mn: 1, O: 4 });
expect('K2Cr2O7',    parseFormula('K2Cr2O7'),     { K: 2, Cr: 2, O: 7 });
// Hydrate
expect('CuSO4·5H2O', parseFormula('CuSO4·5H2O'), { Cu: 1, S: 1, O: 9, H: 10 });
expect('CuSO4.5H2O', parseFormula('CuSO4.5H2O'), { Cu: 1, S: 1, O: 9, H: 10 });
// NIST suffix stripping
expect('CH4-N2 strip', parseFormula('CH4-N2'),    { C: 1, H: 4 }); // artifact stripped
// Charges ignored
expect('H+ charge',  parseFormula('H+'),  { H: 1 });
expect('OH- charge', parseFormula('OH-'), { O: 1, H: 1 });
// Nested parens
expect('Ca3(PO4)2',  parseFormula('Ca3(PO4)2'),  { Ca: 3, P: 2, O: 8 });
// Empty
expect('empty',      parseFormula(''),            {});
// Single element
expect('O',          parseFormula('O'),           { O: 1 });
// Complex organic
expect('C6H5OH',     parseFormula('C6H5OH'),      { C: 6, H: 6, O: 1 });
// Bracket notation
expect('[Fe(CN)6]',  parseFormula('[Fe(CN)6]'),   { Fe: 1, C: 6, N: 6 });

// ─── checkBalance tests ───────────────────────────────────────────────────────
console.log('\n=== checkBalance ===');

// CH4 + 2O2 → CO2 + 2H2O  (correctly balanced)
expect('CH4 combustion balanced', checkBalance(
  [{ formula: 'CH4', coefficient: 1 }, { formula: 'O2', coefficient: 2 }],
  [{ formula: 'CO2', coefficient: 1 }, { formula: 'H2O', coefficient: 2 }],
).status, 'balanced');

// CH4 + O2 → CO2 + H2O  (unbalanced but fixable)
expect('CH4 combustion fixable', checkBalance(
  [{ formula: 'CH4', coefficient: 1 }, { formula: 'O2', coefficient: 1 }],
  [{ formula: 'CO2', coefficient: 1 }, { formula: 'H2O', coefficient: 1 }],
).status, 'fixable');

// C3H8 + 5O2 → 3CO2 + 4H2O  (balanced)
expect('propane combustion balanced', checkBalance(
  [{ formula: 'C3H8', coefficient: 1 }, { formula: 'O2', coefficient: 5 }],
  [{ formula: 'CO2', coefficient: 3 }, { formula: 'H2O', coefficient: 4 }],
).status, 'balanced');

// Fe + O2 → Fe2O3 (unbalanced, fixable: 4Fe + 3O2 → 2Fe2O3)
expect('Fe oxidation fixable', checkBalance(
  [{ formula: 'Fe', coefficient: 1 }, { formula: 'O2', coefficient: 1 }],
  [{ formula: 'Fe2O3', coefficient: 1 }],
).status, 'fixable');

// Impossible: H2 + O2 → NaCl (different elements)
expect('impossible reaction', checkBalance(
  [{ formula: 'H2', coefficient: 1 }, { formula: 'O2', coefficient: 1 }],
  [{ formula: 'NaCl', coefficient: 1 }],
).status, 'impossible');

// N2 + H2 → NH3 (fixable: 1:3:2)
expect('Haber process fixable', checkBalance(
  [{ formula: 'N2', coefficient: 1 }, { formula: 'H2', coefficient: 1 }],
  [{ formula: 'NH3', coefficient: 1 }],
).status, 'fixable');

// N2 + 3H2 → 2NH3 (balanced)
expect('Haber process balanced', checkBalance(
  [{ formula: 'N2', coefficient: 1 }, { formula: 'H2', coefficient: 3 }],
  [{ formula: 'NH3', coefficient: 2 }],
).status, 'balanced');

// Photosynthesis: 6CO2 + 6H2O → C6H12O6 + 6O2  (balanced)
expect('photosynthesis balanced', checkBalance(
  [{ formula: 'CO2', coefficient: 6 }, { formula: 'H2O', coefficient: 6 }],
  [{ formula: 'C6H12O6', coefficient: 1 }, { formula: 'O2', coefficient: 6 }],
).status, 'balanced');

// Fractional: CO2 has no C on product side
expect('impossible missing C', checkBalance(
  [{ formula: 'CO2', coefficient: 1 }],
  [{ formula: 'H2O', coefficient: 1 }],
).status, 'impossible');

// Al + O2 → Al2O3 (fixable: 4Al + 3O2 → 2Al2O3)
expect('Al oxidation fixable', checkBalance(
  [{ formula: 'Al', coefficient: 1 }, { formula: 'O2', coefficient: 1 }],
  [{ formula: 'Al2O3', coefficient: 1 }],
).status, 'fixable');

// ─── autoBalance tests ────────────────────────────────────────────────────────
console.log('\n=== autoBalance ===');

function verifyBalanced(
  reactants: { formula: string; coefficient: number }[],
  products:  { formula: string; coefficient: number }[],
  result: { reactantCoeffs: number[]; productCoeffs: number[] } | null,
): boolean {
  if (!result) return false;
  const rAtoms = reactants.map(r => parseFormula(r.formula));
  const pAtoms = products.map(p => parseFormula(p.formula));
  const allEls = new Set<string>();
  for (const a of [...rAtoms, ...pAtoms]) for (const el of Object.keys(a)) allEls.add(el);
  for (const el of allEls) {
    let lhs = 0, rhs = 0;
    for (let i = 0; i < reactants.length; i++) lhs += (rAtoms[i][el] ?? 0) * result.reactantCoeffs[i];
    for (let i = 0; i < products.length;  i++) rhs += (pAtoms[i][el]  ?? 0) * result.productCoeffs[i];
    if (Math.abs(lhs - rhs) > 1e-6) return false;
  }
  return true;
}

// CH4 + O2 → CO2 + H2O  →  1:2:1:2
{
  const rr = [{ formula: 'CH4', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const pp = [{ formula: 'CO2', coefficient: 1 }, { formula: 'H2O', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('CH4 combustion balanced', res, v => verifyBalanced(rr, pp, v as any));
  expectDeep('CH4 coefficients correct', res, v =>
    JSON.stringify((v as any)?.reactantCoeffs) === '[1,2]' &&
    JSON.stringify((v as any)?.productCoeffs)  === '[1,2]',
  );
}

// N2 + H2 → NH3  →  1:3:2
{
  const rr = [{ formula: 'N2', coefficient: 1 }, { formula: 'H2', coefficient: 1 }];
  const pp = [{ formula: 'NH3', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('Haber balanced', res, v => verifyBalanced(rr, pp, v as any));
  expectDeep('Haber coefficients', res, v =>
    JSON.stringify((v as any)?.reactantCoeffs) === '[1,3]' &&
    JSON.stringify((v as any)?.productCoeffs)  === '[2]',
  );
}

// Fe + O2 → Fe2O3  →  4:3:2
{
  const rr = [{ formula: 'Fe', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const pp = [{ formula: 'Fe2O3', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('Fe oxidation balanced', res, v => verifyBalanced(rr, pp, v as any));
  expectDeep('Fe oxidation coefficients', res, v =>
    JSON.stringify((v as any)?.reactantCoeffs) === '[4,3]' &&
    JSON.stringify((v as any)?.productCoeffs)  === '[2]',
  );
}

// Al + O2 → Al2O3  →  4:3:2
{
  const rr = [{ formula: 'Al', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const pp = [{ formula: 'Al2O3', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('Al oxidation balanced', res, v => verifyBalanced(rr, pp, v as any));
}

// C3H8 + O2 → CO2 + H2O  →  1:5:3:4
{
  const rr = [{ formula: 'C3H8', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const pp = [{ formula: 'CO2', coefficient: 1 }, { formula: 'H2O', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('propane balanced', res, v => verifyBalanced(rr, pp, v as any));
  expectDeep('propane coefficients', res, v =>
    JSON.stringify((v as any)?.reactantCoeffs) === '[1,5]' &&
    JSON.stringify((v as any)?.productCoeffs)  === '[3,4]',
  );
}

// CaCO3 → CaO + CO2  →  1:1:1
{
  const rr = [{ formula: 'CaCO3', coefficient: 1 }];
  const pp = [{ formula: 'CaO', coefficient: 1 }, { formula: 'CO2', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('CaCO3 decomposition balanced', res, v => verifyBalanced(rr, pp, v as any));
}

// KMnO4 + HCl → KCl + MnCl2 + H2O + Cl2  (6-species, fixable: 2:16:2:2:8:5)
{
  const rr = [{ formula: 'KMnO4', coefficient: 1 }, { formula: 'HCl', coefficient: 1 }];
  const pp = [
    { formula: 'KCl', coefficient: 1 },
    { formula: 'MnCl2', coefficient: 1 },
    { formula: 'H2O', coefficient: 1 },
    { formula: 'Cl2', coefficient: 1 },
  ];
  const res = autoBalance(rr, pp);
  expectDeep('KMnO4+HCl balanced', res, v => verifyBalanced(rr, pp, v as any));
}

// Impossible: reactants H2+O2, products NaCl
{
  const rr = [{ formula: 'H2', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const pp = [{ formula: 'NaCl', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expect('impossible returns null', res, null);
}

// Already balanced stays the same (verifies it returns a solution anyway)
{
  const rr = [{ formula: 'CH4', coefficient: 1 }, { formula: 'O2', coefficient: 2 }];
  const pp = [{ formula: 'CO2', coefficient: 1 }, { formula: 'H2O', coefficient: 2 }];
  const res = autoBalance(rr, pp);
  expectDeep('already-balanced still returns solution', res, v => verifyBalanced(rr, pp, v as any));
}

// Photosynthesis: 6CO2 + 6H2O → C6H12O6 + 6O2
{
  const rr = [{ formula: 'CO2', coefficient: 1 }, { formula: 'H2O', coefficient: 1 }];
  const pp = [{ formula: 'C6H12O6', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('photosynthesis balanced', res, v => verifyBalanced(rr, pp, v as any));
  expectDeep('photosynthesis coefficients', res, v =>
    JSON.stringify((v as any)?.reactantCoeffs) === '[6,6]' &&
    JSON.stringify((v as any)?.productCoeffs)  === '[1,6]',
  );
}

// Single-reactant decomposition: H2O2 → H2O + O2  →  2:2:1
{
  const rr = [{ formula: 'H2O2', coefficient: 1 }];
  const pp = [{ formula: 'H2O', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('H2O2 decomposition balanced', res, v => verifyBalanced(rr, pp, v as any));
  expectDeep('H2O2 decomposition coefficients', res, v =>
    JSON.stringify((v as any)?.reactantCoeffs) === '[2]' &&
    JSON.stringify((v as any)?.productCoeffs)  === '[2,1]',
  );
}

// Edge: single species each side (trivial, same formula)
{
  const rr = [{ formula: 'H2O', coefficient: 1 }];
  const pp = [{ formula: 'H2O', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('H2O → H2O trivial', res, v => verifyBalanced(rr, pp, v as any));
}

// Edge: large coefficients needed — P4 + O2 → P4O10 (1:5:1) — formerly needs 5
{
  const rr = [{ formula: 'P4', coefficient: 1 }, { formula: 'O2', coefficient: 1 }];
  const pp = [{ formula: 'P4O10', coefficient: 1 }];
  const res = autoBalance(rr, pp);
  expectDeep('P4+O2 balanced', res, v => verifyBalanced(rr, pp, v as any));
}

// ─── Summary ─────────────────────────────────────────────────────────────────
console.log(`\n=== Results: ${passed} passed, ${failed} failed ===\n`);
if (failed > 0) process.exit(1);
