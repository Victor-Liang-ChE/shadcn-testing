import type { CanopyCompound, FluidPackage } from './types';
import type { InteractionParams, NRTLParams, UNIQUACParams, WilsonParams } from './thermo';
import { liquidMolarVolume_m3pmol } from './thermo';

export interface AspenBinaryParamRow {
  Compound1Name: string;
  Compound2Name: string;
  ElementLabel: string;
  Value: string | number;
}

function createMatrix(size: number, fill = 0): number[][] {
  return Array.from({ length: size }, () => new Array(size).fill(fill));
}

function buildIndex(names: string[]): Map<string, number> {
  return new Map(names.map((name, index) => [name.toUpperCase(), index]));
}

function parseValue(raw: string | number): number | null {
  const value = typeof raw === 'number' ? raw : parseFloat(raw);
  return Number.isFinite(value) ? value : null;
}

function getPair(
  row: AspenBinaryParamRow,
  indexByName: Map<string, number>,
): { i: number; j: number } | null {
  const i = indexByName.get(row.Compound1Name.toUpperCase());
  const j = indexByName.get(row.Compound2Name.toUpperCase());
  if (i === undefined || j === undefined || i === j) return null;
  return { i, j };
}

export function buildNrtlParamsFromAspenRows(
  compounds: CanopyCompound[],
  rows: AspenBinaryParamRow[],
): NRTLParams {
  const n = compounds.length;
  const names = compounds.map(c => c.name);
  const indexByName = buildIndex(names);
  const A = createMatrix(n);
  const B = createMatrix(n);
  const alpha = createMatrix(n, 0.3);
  const a_ext = createMatrix(n);
  const b_ext = createMatrix(n);
  const e_ext = createMatrix(n);
  const f_ext = createMatrix(n);
  const d_alpha = createMatrix(n);
  let hasExtended = false;

  for (const row of rows) {
    const pair = getPair(row, indexByName);
    if (!pair) continue;
    const value = parseValue(row.Value);
    if (value === null) continue;
    const { i, j } = pair;
    const label = row.ElementLabel.toLowerCase();
    if (label === 'aij') { a_ext[i][j] = value; hasExtended = true; }
    else if (label === 'aji') { a_ext[j][i] = value; hasExtended = true; }
    else if (label === 'bij') { b_ext[i][j] = value; hasExtended = true; }
    else if (label === 'bji') { b_ext[j][i] = value; hasExtended = true; }
    else if (label === 'cij') { alpha[i][j] = value; alpha[j][i] = value; }
    else if (label === 'dij') { d_alpha[i][j] = value; d_alpha[j][i] = value; hasExtended = true; }
    else if (label === 'eij') { e_ext[i][j] = value; hasExtended = true; }
    else if (label === 'eji') { e_ext[j][i] = value; hasExtended = true; }
    else if (label === 'fij') { f_ext[i][j] = value; hasExtended = true; }
    else if (label === 'fji') { f_ext[j][i] = value; hasExtended = true; }
  }

  for (let i = 0; i < n; i++) alpha[i][i] = 0;
  return {
    A,
    B,
    alpha,
    ...(hasExtended ? { a_ext, b_ext, e_ext, f_ext, d_alpha } : {}),
  };
}

export function buildUniquacParamsFromAspenRows(
  compounds: CanopyCompound[],
  rows: AspenBinaryParamRow[],
): UNIQUACParams {
  const n = compounds.length;
  const names = compounds.map(c => c.name);
  const indexByName = buildIndex(names);
  const a = createMatrix(n);
  const b = createMatrix(n);
  const a_ext = createMatrix(n);
  const b_ext = createMatrix(n);
  const c_ext = createMatrix(n);
  const d_ext = createMatrix(n);
  const e_ext = createMatrix(n);
  let hasExtended = false;
  const r = compounds.map(c => c.uniquac_r ?? 1);
  const q = compounds.map(c => c.uniquac_q ?? 1);

  for (const row of rows) {
    const pair = getPair(row, indexByName);
    if (!pair) continue;
    const value = parseValue(row.Value);
    if (value === null) continue;
    const { i, j } = pair;
    const label = row.ElementLabel.toLowerCase();
    if (label === 'aij') { a_ext[i][j] = value; hasExtended = true; }
    else if (label === 'aji') { a_ext[j][i] = value; hasExtended = true; }
    else if (label === 'bij') { b_ext[i][j] = value; hasExtended = true; }
    else if (label === 'bji') { b_ext[j][i] = value; hasExtended = true; }
    else if (label === 'cij') { c_ext[i][j] = value; hasExtended = true; }
    else if (label === 'cji') { c_ext[j][i] = value; hasExtended = true; }
    else if (label === 'dij') { d_ext[i][j] = value; hasExtended = true; }
    else if (label === 'dji') { d_ext[j][i] = value; hasExtended = true; }
    else if (label === 'eij') { e_ext[i][j] = value; hasExtended = true; }
    else if (label === 'eji') { e_ext[j][i] = value; hasExtended = true; }
  }

  return {
    a,
    b,
    r,
    q,
    ...(hasExtended ? { a_ext, b_ext, c_ext, d_ext, e_ext } : {}),
  };
}

export function buildWilsonParamsFromAspenRows(
  compounds: CanopyCompound[],
  rows: AspenBinaryParamRow[],
): WilsonParams {
  const n = compounds.length;
  const names = compounds.map(c => c.name);
  const indexByName = buildIndex(names);
  const A = createMatrix(n);
  const a_ext = createMatrix(n);
  const b_ext = createMatrix(n);
  const c_ext = createMatrix(n);
  const d_ext = createMatrix(n);
  const e_ext = createMatrix(n);
  let hasExtended = false;
  const molarVolumes = compounds.map(c => {
    const vM3PerMol = liquidMolarVolume_m3pmol(c, 298.15);
    return Number.isFinite(vM3PerMol) && vM3PerMol > 0 ? vM3PerMol * 1e6 : 100;
  });

  for (const row of rows) {
    const pair = getPair(row, indexByName);
    if (!pair) continue;
    const value = parseValue(row.Value);
    if (value === null) continue;
    const { i, j } = pair;
    const label = row.ElementLabel.toLowerCase();
    if (label === 'aij') { a_ext[i][j] = value; hasExtended = true; }
    else if (label === 'aji') { a_ext[j][i] = value; hasExtended = true; }
    else if (label === 'bij') { b_ext[i][j] = value; hasExtended = true; }
    else if (label === 'bji') { b_ext[j][i] = value; hasExtended = true; }
    else if (label === 'cij') { c_ext[i][j] = value; hasExtended = true; }
    else if (label === 'cji') { c_ext[j][i] = value; hasExtended = true; }
    else if (label === 'dij') { d_ext[i][j] = value; hasExtended = true; }
    else if (label === 'dji') { d_ext[j][i] = value; hasExtended = true; }
    else if (label === 'eij') { e_ext[i][j] = value; hasExtended = true; }
    else if (label === 'eji') { e_ext[j][i] = value; hasExtended = true; }
  }

  return {
    A,
    molarVolumes,
    ...(hasExtended ? { a_ext, b_ext, c_ext, d_ext, e_ext } : {}),
  };
}

export function buildEosInteractionParamsFromAspenRows(
  compounds: CanopyCompound[],
  rows: AspenBinaryParamRow[],
): Pick<InteractionParams, 'kij' | 'kij_a' | 'kij_b' | 'kij_c'> {
  const n = compounds.length;
  const names = compounds.map(c => c.name);
  const indexByName = buildIndex(names);
  const kij = createMatrix(n);
  const kij_a = createMatrix(n);
  const kij_b = createMatrix(n);
  const kij_c = createMatrix(n);
  let hasTempDependent = false;

  for (const row of rows) {
    const pair = getPair(row, indexByName);
    if (!pair) continue;
    const value = parseValue(row.Value);
    if (value === null) continue;
    const { i, j } = pair;
    const label = row.ElementLabel.toLowerCase();
    if (label === 'kij') {
      kij[i][j] = value;
      kij[j][i] = value;
    } else if (label === 'aij' || label === 'aji') {
      kij_a[i][j] = value;
      kij_a[j][i] = value;
      hasTempDependent = true;
    } else if (label === 'bij' || label === 'bji') {
      kij_b[i][j] = value;
      kij_b[j][i] = value;
      hasTempDependent = true;
    } else if (label === 'cij' || label === 'cji') {
      kij_c[i][j] = value;
      kij_c[j][i] = value;
      hasTempDependent = true;
    }
  }

  return {
    kij,
    ...(hasTempDependent ? { kij_a, kij_b, kij_c } : {}),
  };
}

export function buildInteractionParamsForPackage(
  pkg: FluidPackage,
  compounds: CanopyCompound[],
  rows: AspenBinaryParamRow[],
): InteractionParams | null {
  if (pkg === 'NRTL') {
    return { nrtl: buildNrtlParamsFromAspenRows(compounds, rows) };
  }
  if (pkg === 'UNIQUAC') {
    return { uniquac: buildUniquacParamsFromAspenRows(compounds, rows) };
  }
  if (pkg === 'Wilson') {
    return { wilson: buildWilsonParamsFromAspenRows(compounds, rows) };
  }
  if (pkg === 'Peng-Robinson' || pkg === 'SRK') {
    return buildEosInteractionParamsFromAspenRows(compounds, rows);
  }
  return null;
}
