/**
 * reaction-balancer.ts
 *
 * Chemical formula parser and reaction balancer.
 *
 * parseFormula(str) → Record<string, number>  — atom counts
 * checkBalance(reactants, products) → BalanceResult
 * autoBalance(reactants, products) → number[] | null  — new coefficients or null if impossible
 */

// ─── Element detection ────────────────────────────────────────────────────────
// All 2-letter element symbols
const ELEMENTS_2 = new Set([
  'He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti','Cr','Mn',
  'Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Zr','Nb',
  'Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','Xe','Cs','Ba','La',
  'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
  'Ta','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
  'Ac','Th','Pa','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',
  'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og',
]);
const ELEMENTS_1 = new Set([
  'H','B','C','N','O','F','P','S','K','V','W','Y','I','U',
]);
function isElement(sym: string): boolean {
  return sym.length === 1 ? ELEMENTS_1.has(sym) : ELEMENTS_2.has(sym);
}

// Uppercase-letter pairs that represent 2-letter element symbols in NIST Hill notation.
// Both individual letters happen to be valid 1-letter element symbols too, so the
// ambiguity is resolved by Hill ordering: the second-letter element always comes
// BEFORE the first-letter element alphabetically, so they'd never appear adjacent
// as two separate atoms (e.g. I precedes N → N+I in Hill = "...I...N...", never "NI").
// CO and NO are intentionally excluded: they almost always mean C+O / N+O in org-chemistry.
const NIST_FORCE_2LETTER = new Set([
  'NI','CU','SN','SI','PB','PT','PD','AU','HG','TL','IN','BI',
]);

/** Convert a NIST all-uppercase formula (e.g. "C2H3CL", "NAOH", "FE") to
 * proper chemical notation ("C2H3Cl", "NaOH", "Fe") before parsing. */
function normalizeNistUppercase(f: string): string {
  let result = '';
  let i = 0;
  while (i < f.length) {
    const ch = f[i];
    if (/[A-Z]/.test(ch) && i + 1 < f.length && /[A-Z]/.test(f[i + 1])) {
      const second = f[i + 1];
      const candidate = ch + second.toLowerCase();
      if (ELEMENTS_2.has(candidate)) {
        // Safe to use the 2-letter symbol when:
        //   a) The second letter alone is NOT a valid element (no ambiguity), OR
        //   b) Hill ordering guarantees this pair = the 2-letter element.
        if (!ELEMENTS_1.has(second) || NIST_FORCE_2LETTER.has(ch + second)) {
          result += candidate;
          i += 2;
          continue;
        }
      }
    }
    result += ch;
    i++;
  }
  return result;
}

// ─── Formula parser ───────────────────────────────────────────────────────────
/**
 * Parse a chemical formula string into an atom-count map.
 *
 * Handles:
 *   - Standard formulas:     H2O, C6H12O6, Fe2O3
 *   - Parentheses/brackets:  Ca(OH)2, Al2(SO4)3, [Fe(CN)6]4-
 *   - Dots (hydrates):       CuSO4·5H2O, Na2SO4.10H2O
 *   - Charges (ignored):     SO4²⁻, H+, OH-, NH4+, PO4-3
 *   - Radical marker:        CH3• (bullet ignored)
 *   - NIST suffixes:         CH4-N2 (the "-Xn" artifact is stripped)
 *   - Isotopes not supported (treated as unknown element)
 *
 * Returns {} on empty/invalid input without throwing.
 */
export function parseFormula(raw: string): Record<string, number> {
  if (!raw || typeof raw !== 'string') return {};

  // Strip NIST artifact suffixes like "-N3" at end
  let f = raw.replace(/-[A-Z][a-z]?\d*$/, '');

  // Strip all-digit NIST suffixes like "-2" (e.g. isobutyl acetate "C6H12O2-2")
  // Must be done BEFORE the charge regex below, which would otherwise consume
  // the digit before the dash too (e.g. "O2-2" → would strip "2-2" leaving "O").
  f = f.replace(/-\d+$/, '');

  // Strip charge indicators at the end: +, -, +2, -3, 2+, 3-, ²⁻, etc.
  f = f.replace(/[²³⁴⁺⁻]+$/, '');         // unicode superscripts
  f = f.replace(/[\d]*[+-][\d]*$/, '');     // ASCII charges like 2+, -3, 4+2-

  // Replace dots and middle-dots (hydrate connector) with a left paren-style
  // by turning "·" into a multiply — we simply replace ·N with a multiplier
  // by inserting a group. Easiest: split on · and multiply.
  // But dots can also appear mid-formula. Simplest approach: replace · with +
  // and parse each part, summing results.
  // Normalize NIST all-uppercase formulas → proper chemical notation.
  // NIST stores formulas fully uppercase (C2H3CL, C2H3BR, FE, NAOH, NI…).
  // Must run after suffix/charge stripping so there are no stray hyphens.
  if (!/[a-z]/.test(f)) f = normalizeNistUppercase(f);

  if (f.includes('·') || f.includes('.')) {
    const parts = f.split(/[·.]/);
    const result: Record<string, number> = {};
    for (const part of parts) {
      const trimmed = part.trim();
      if (!trimmed) continue;
      // Check if part starts with a multiplier number e.g. "5H2O"
      const match = trimmed.match(/^(\d+)([A-Z].*)$/);
      const [, mult, rest] = match ? match : ['', '1', trimmed];
      const n = parseInt(mult || '1', 10);
      const sub = parseFormulaInner(rest);
      for (const [el, cnt] of Object.entries(sub)) {
        result[el] = (result[el] ?? 0) + cnt * n;
      }
    }
    return result;
  }

  return parseFormulaInner(f);
}

/** Core recursive parser — handles letters, digits, parens/brackets. */
function parseFormulaInner(f: string): Record<string, number> {
  const stack: Record<string, number>[] = [{}];
  let i = 0;

  while (i < f.length) {
    const ch = f[i];

    // Open group
    if (ch === '(' || ch === '[' || ch === '{') {
      stack.push({});
      i++;
      continue;
    }

    // Close group
    if (ch === ')' || ch === ']' || ch === '}') {
      i++;
      // Consume following number
      let numStr = '';
      while (i < f.length && /\d/.test(f[i])) { numStr += f[i++]; }
      const mult = numStr ? parseInt(numStr, 10) : 1;
      const closed = stack.pop() ?? {};
      const top = stack[stack.length - 1];
      for (const [el, cnt] of Object.entries(closed)) {
        top[el] = (top[el] ?? 0) + cnt * mult;
      }
      continue;
    }

    // Element symbol: uppercase letter, optionally followed by lowercase letters
    if (/[A-Z]/.test(ch)) {
      let sym = ch;
      // Try 2-char element first
      if (i + 1 < f.length && /[a-z]/.test(f[i + 1])) {
        const candidate = ch + f[i + 1];
        if (isElement(candidate)) {
          sym = candidate;
          i += 2;
        } else {
          i++;
        }
      } else {
        i++;
      }

      // Consume following digits
      let numStr = '';
      while (i < f.length && /\d/.test(f[i])) { numStr += f[i++]; }
      const count = numStr ? parseInt(numStr, 10) : 1;

      const top = stack[stack.length - 1];
      if (isElement(sym)) {
        top[sym] = (top[sym] ?? 0) + count;
      }
      // Non-element uppercase sequences (e.g. from malformed input) are silently skipped
      continue;
    }

    // Skip lowercase (shouldn't appear outside elements), digits consumed above,
    // charges, radical markers *, •, spaces, etc.
    i++;
  }

  return stack[0] ?? {};
}

// ─── Balance check ────────────────────────────────────────────────────────────
export type BalanceStatus =
  | 'balanced'      // atom counts match with current coefficients
  | 'fixable'       // same element set — integer coefficients exist that would balance it
  | 'impossible';   // element sets differ — can't be balanced by changing coefficients alone

export interface BalanceResult {
  status: BalanceStatus;
  /** Set of element symbols present on left but not right (or vice versa) */
  missingElements: string[];
  /** Suggested balanced coefficients [reactants..., products...] or null */
  suggestedCoeffs: number[] | null;
}

/**
 * Check whether the reaction is balanced.
 *
 * @param reactants  Array of { formula, coefficient } for the left side
 * @param products   Array of { formula, coefficient } for the right side
 */
export function checkBalance(
  reactants: { formula: string; coefficient: number }[],
  products:  { formula: string; coefficient: number }[],
): BalanceResult {
  if (reactants.length === 0 || products.length === 0) {
    return { status: 'balanced', missingElements: [], suggestedCoeffs: null };
  }

  // Parse all formulas
  const rAtoms = reactants.map(r => parseFormula(r.formula));
  const pAtoms = products.map(p => parseFormula(p.formula));

  // Collect all elements
  const allEls = new Set<string>();
  for (const a of [...rAtoms, ...pAtoms]) for (const el of Object.keys(a)) allEls.add(el);

  // Check which elements appear on only one side
  const rEls = new Set<string>();
  const pEls = new Set<string>();
  for (const a of rAtoms) for (const el of Object.keys(a)) rEls.add(el);
  for (const a of pAtoms) for (const el of Object.keys(a)) pEls.add(el);
  const missing: string[] = [];
  for (const el of allEls) {
    if (!rEls.has(el) || !pEls.has(el)) missing.push(el);
  }

  if (missing.length > 0) {
    return { status: 'impossible', missingElements: missing, suggestedCoeffs: null };
  }

  // Check if currently balanced
  const isBalancedNow = isCurrentlyBalanced(reactants, products, rAtoms, pAtoms, allEls);
  if (isBalancedNow) {
    return { status: 'balanced', missingElements: [], suggestedCoeffs: null };
  }

  // Try to find integer coefficients via null-space method
  const solution = solveBalancing(rAtoms, pAtoms, [...allEls]);
  if (solution) {
    return { status: 'fixable', missingElements: [], suggestedCoeffs: solution };
  }

  // Null-space method failed — try brute force for small systems
  const bruteForce = bruteForceBalance(rAtoms, pAtoms, [...allEls]);
  if (bruteForce) {
    return { status: 'fixable', missingElements: [], suggestedCoeffs: bruteForce };
  }

  return { status: 'impossible', missingElements: [], suggestedCoeffs: null };
}

function isCurrentlyBalanced(
  reactants: { formula: string; coefficient: number }[],
  products:  { formula: string; coefficient: number }[],
  rAtoms: Record<string, number>[],
  pAtoms: Record<string, number>[],
  allEls: Set<string>,
): boolean {
  const EPSILON = 1e-9;
  for (const el of allEls) {
    let lhs = 0, rhs = 0;
    for (let i = 0; i < reactants.length; i++) lhs += (rAtoms[i][el] ?? 0) * reactants[i].coefficient;
    for (let i = 0; i < products.length;  i++) rhs += (pAtoms[i][el]  ?? 0) * products[i].coefficient;
    if (Math.abs(lhs - rhs) > EPSILON * Math.max(1, Math.abs(lhs))) return false;
  }
  return true;
}

// ─── Balancing via null-space of stoichiometry matrix ────────────────────────
/**
 * Build the stoichiometry matrix:
 *   rows = elements, columns = species (reactants negative, products positive)
 * Find an integer null-space vector (all positive).
 *
 * Uses exact rational arithmetic (fractions as [num, den] pairs) to avoid fp errors.
 */
function solveBalancing(
  rAtoms: Record<string, number>[],
  pAtoms: Record<string, number>[],
  elements: string[],
): number[] | null {
  const nR = rAtoms.length;
  const nP = pAtoms.length;
  const nSpecies = nR + nP;
  const nEl = elements.length;

  // Build augmented matrix [A | 0] as rational fractions [num, den]
  // A[el][species]: reactants are negative (moving to RHS of balance), products positive
  type Frac = [number, number]; // [numerator, denominator]
  const gcd = (a: number, b: number): number => b === 0 ? Math.abs(a) : gcd(b, a % b);
  const frac = (n: number, d = 1): Frac => {
    if (d < 0) { n = -n; d = -d; }
    const g = gcd(Math.abs(n), d);
    return [n / g, d / g];
  };
  const fadd = (a: Frac, b: Frac): Frac => frac(a[0] * b[1] + b[0] * a[1], a[1] * b[1]);
  const fmul = (a: Frac, b: Frac): Frac => frac(a[0] * b[0], a[1] * b[1]);
  const fsub = (a: Frac, b: Frac): Frac => fadd(a, frac(-b[0], b[1]));
  const fdiv = (a: Frac, b: Frac): Frac => fmul(a, frac(b[1], b[0]));
  const fzero: Frac = [0, 1];
  const isZ = (f: Frac) => f[0] === 0;

  // Matrix rows = elements, columns = species + 1 identity augment
  // We solve A*x=0 via row reduction and read off free variables
  const rows = nEl;
  const cols = nSpecies;

  // Build matrix
  const M: Frac[][] = Array.from({ length: rows }, (_, r) => {
    const el = elements[r];
    return Array.from({ length: cols }, (_, c) => {
      if (c < nR) return frac(-(rAtoms[c][el] ?? 0)); // reactant: negative
      return frac(pAtoms[c - nR][el] ?? 0);           // product: positive
    });
  });

  // Gaussian elimination to row echelon form
  const pivotCols: number[] = [];
  let pivotRow = 0;
  for (let col = 0; col < cols && pivotRow < rows; col++) {
    // Find pivot in this column
    let found = -1;
    for (let r = pivotRow; r < rows; r++) {
      if (!isZ(M[r][col])) { found = r; break; }
    }
    if (found === -1) continue;
    // Swap
    [M[pivotRow], M[found]] = [M[found], M[pivotRow]];
    // Eliminate
    const pivot = M[pivotRow][col];
    for (let r = 0; r < rows; r++) {
      if (r === pivotRow) continue;
      if (isZ(M[r][col])) continue;
      const factor = fdiv(M[r][col], pivot);
      for (let c = 0; c < cols; c++) {
        M[r][c] = fsub(M[r][c], fmul(factor, M[pivotRow][c]));
      }
    }
    pivotCols.push(col);
    pivotRow++;
  }

  const rank = pivotCols.length;
  const freeCols = Array.from({ length: cols }, (_, i) => i).filter(c => !pivotCols.includes(c));

  // System is underdetermined; we need exactly 1 free variable for a unique solution (up to scale)
  if (freeCols.length !== 1) {
    // Multiple free variables: try setting all free vars = 1 (works for many common cases)
    if (freeCols.length > 1) {
      return tryFreeVarsAll1(M, pivotCols, freeCols, nSpecies);
    }
    return null;
  }

  const freeCol = freeCols[0];
  // Back-substitute: set x[freeCol] = 1 and solve for pivot variables
  const x: Frac[] = new Array(nSpecies).fill(fzero);
  x[freeCol] = frac(1);

  for (let i = pivotCols.length - 1; i >= 0; i--) {
    const pc = pivotCols[i];
    const row = M[i]; // after elimination the pivot row i has pivot at pc
    // Find which row has pivot at pc
    const pivRow = M.findIndex(r => !isZ(r[pc]));
    if (pivRow === -1) continue;
    let sum: Frac = fzero;
    for (let c = 0; c < nSpecies; c++) {
      if (c === pc) continue;
      sum = fadd(sum, fmul(M[pivRow][c], x[c]));
    }
    x[pc] = fdiv(frac(-sum[0], sum[1]), M[pivRow][pc]);
  }

  return normalizeSolution(x, nSpecies);
}

function tryFreeVarsAll1(
  M: [number, number][][],
  pivotCols: number[],
  freeCols: number[],
  nSpecies: number,
): number[] | null {
  type Frac = [number, number];
  const gcd = (a: number, b: number): number => b === 0 ? Math.abs(a) : gcd(b, a % b);
  const frac = (n: number, d = 1): Frac => {
    if (d < 0) { n = -n; d = -d; }
    const g = gcd(Math.abs(n), d);
    return [n / g, d / g];
  };
  const fmul = (a: Frac, b: Frac): Frac => frac(a[0] * b[0], a[1] * b[1]);
  const fadd = (a: Frac, b: Frac): Frac => frac(a[0] * b[1] + b[0] * a[1], a[1] * b[1]);
  const fdiv = (a: Frac, b: Frac): Frac => fmul(a, frac(b[1], b[0]));
  const fzero: Frac = [0, 1];
  const isZ = (f: Frac) => f[0] === 0;

  const x: Frac[] = new Array(nSpecies).fill(fzero);
  for (const fc of freeCols) x[fc] = frac(1);

  for (let i = pivotCols.length - 1; i >= 0; i--) {
    const pc = pivotCols[i];
    const pivRow = M.findIndex(r => !isZ(r[pc]));
    if (pivRow === -1) continue;
    let sum: Frac = fzero;
    for (let c = 0; c < nSpecies; c++) {
      if (c === pc) continue;
      sum = fadd(sum, fmul(M[pivRow][c], x[c]));
    }
    x[pc] = fdiv(frac(-sum[0], sum[1]), M[pivRow][pc]);
  }
  return normalizeSolution(x, nSpecies);
}

function normalizeSolution(x: [number, number][], nSpecies: number): number[] | null {
  const gcd = (a: number, b: number): number => b === 0 ? Math.abs(a) : gcd(b, a % b);

  // Check all defined
  if (x.some(f => f[1] === 0)) return null;

  // Multiply through by LCM of denominators to get integers
  let lcm = 1;
  for (const [, d] of x) {
    if (d === 0) return null;
    const g = gcd(Math.abs(lcm), Math.abs(d));
    lcm = Math.abs(lcm * d / g);
    if (lcm > 1e9) return null; // overflow guard
  }

  const ints = x.map(([n, d]) => Math.round((n * lcm) / d));

  // All must be non-zero
  if (ints.some(v => v === 0)) return null;

  // If any negative, flip all signs
  const anyNeg = ints.some(v => v < 0);
  const signs = anyNeg ? ints.map(v => -v) : ints;

  // Check all positive after flip
  if (signs.some(v => v <= 0)) return null;

  // Reduce to lowest terms
  let g = signs[0];
  for (const v of signs) { const gg = gcd(Math.abs(g), Math.abs(v)); if (gg > 0) g = gg; }
  return signs.map(v => v / g);
}

// ─── Brute-force for small systems ───────────────────────────────────────────
/**
 * For reactions with ≤ 6 species and small coefficients,
 * enumerate integer coefficients 1..12 to find a balanced solution.
 * Much slower but catches edge cases the algebraic method misses.
 */
function bruteForceBalance(
  rAtoms: Record<string, number>[],
  pAtoms: Record<string, number>[],
  elements: string[],
): number[] | null {
  const nR = rAtoms.length;
  const nP = pAtoms.length;
  const nSpecies = nR + nP;
  if (nSpecies > 6) return null; // too slow

  const MAX_COEFF = 12;

  function check(coeffs: number[]): boolean {
    for (const el of elements) {
      let lhs = 0, rhs = 0;
      for (let i = 0; i < nR; i++) lhs += (rAtoms[i][el] ?? 0) * coeffs[i];
      for (let i = 0; i < nP; i++) rhs += (pAtoms[i][el]  ?? 0) * coeffs[nR + i];
      if (lhs !== rhs) return false;
    }
    return true;
  }

  // Fix first coefficient = 1 (by convention) and enumerate the rest
  function recurse(pos: number, current: number[]): number[] | null {
    if (pos === nSpecies) {
      return check(current) ? [...current] : null;
    }
    const start = pos === 0 ? 1 : 1;
    for (let v = start; v <= MAX_COEFF; v++) {
      current[pos] = v;
      const result = recurse(pos + 1, current);
      if (result) return result;
    }
    return null;
  }

  return recurse(0, new Array(nSpecies).fill(1));
}

// ─── Public convenience wrapper ───────────────────────────────────────────────
/**
 * Auto-balance a reaction given arrays of formula strings and current coefficients.
 * Returns { reactantCoeffs, productCoeffs } or null if impossible.
 */
export function autoBalance(
  reactants: { formula: string; coefficient: number }[],
  products:  { formula: string; coefficient: number }[],
): { reactantCoeffs: number[]; productCoeffs: number[] } | null {
  const rAtoms = reactants.map(r => parseFormula(r.formula));
  const pAtoms = products.map(p => parseFormula(p.formula));
  const allEls = new Set<string>();
  for (const a of [...rAtoms, ...pAtoms]) for (const el of Object.keys(a)) allEls.add(el);

  const solution = solveBalancing(rAtoms, pAtoms, [...allEls])
    ?? bruteForceBalance(rAtoms, pAtoms, [...allEls]);

  if (!solution) return null;

  return {
    reactantCoeffs: solution.slice(0, reactants.length),
    productCoeffs:  solution.slice(reactants.length),
  };
}
