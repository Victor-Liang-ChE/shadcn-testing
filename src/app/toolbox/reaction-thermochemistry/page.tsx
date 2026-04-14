'use client';

import React, { useState, useRef, useCallback, useMemo, useEffect } from 'react';
import { supabase } from '@/lib/supabaseClient';
import { useTheme } from 'next-themes';

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import { BarChart, LineChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import type { EChartsOption } from 'echarts';

echarts.use([TitleComponent, TooltipComponent, GridComponent, LegendComponent, BarChart, LineChart, CanvasRenderer]);

import { Card, CardContent, CardHeader, CardTitle } from '@/components/ui/card';
import { Input } from '@/components/ui/input';
import { Label } from '@/components/ui/label';
import { Button } from '@/components/ui/button';
import { Tabs, TabsList, TabsTrigger } from '@/components/ui/tabs';
import {
  Select,
  SelectContent,
  SelectItem,
  SelectTrigger,
  SelectValue,
} from '@/components/ui/select';
import { TooltipProvider, Tooltip, TooltipTrigger, TooltipContent } from '@/components/ui/tooltip';
import { Skeleton } from '@/components/ui/skeleton';

import { PlusCircle, Trash2, Flame, Snowflake, AlertCircle, Download } from 'lucide-react';
import { checkBalance, autoBalance, parseFormula, type BalanceResult } from '@/lib/reaction-balancer';

// ─── Chart color palette (shared between useMemo and legend JSX) ────────────
const CHR_FORM_R  = '#a78bfa'; // ΔHf reactant
const CHR_FORM_P  = '#2dd4bf'; // ΔHf product
const CHR_SENS_R  = '#c4b5fd'; // sensible reactant
const CHR_SENS_P  = '#99f6e4'; // sensible product
const CHR_LAT_R   = '#7c3aed'; // latent reactant
const CHR_LAT_P   = '#0d9488'; // latent product

/** Returns true when all atoms in a formula are the same element (standard-state element). */
function isMonoElemental(formula: string): boolean {
  const atoms = parseFormula(formula);
  return Object.keys(atoms).length === 1;
}


interface NistCompound {

  nist_id: string;
  name: string;
  formula: string;
  cas: string | null;
  mol_weight: number | null;
  hf_gas_298_kjmol: number | null;
  cp_gas_298_jmolk: number | null;
  tboil_k: number | null;
  tfus_k: number | null;
  tc_k: number | null;
  has_shomate: boolean;
  has_phase_change: boolean;
}

interface ShomateRow {
  phase: string;
  eqno: number;
  t_min_k: number;
  t_max_k: number;
  a: number | null; b: number | null; c: number | null; d: number | null;
  e: number | null; f: number | null; g: number | null;
}

interface PhaseTransitionRow {
  tboil_k: number | null;
  tfus_k: number | null;
  hvap_kjmol: number | null;
  hfus_kjmol: number | null;
}

/** Precomputed Joback group-contribution estimates stored in nist_joback_estimates. */
interface JobackEstimate {
  nist_id: string;
  hf_gas_298_kjmol: number | null;
  dgf_gas_298_kjmol: number | null;
  hfus_kjmol: number | null;
  hvap_kjmol: number | null;
  hvap_reference: string | null;   // 'watson_298K' | 'joback_Tb'
  cp_a: number | null;
  cp_b: number | null;
  cp_c: number | null;
  cp_d: number | null;
  cp_t_min_k: number | null;       // valid lower T limit (typically 300 K)
  cp_t_max_k: number | null;       // valid upper T limit (typically 1000 K)
  joback_group_count: number;
  solid_correction_incomplete: boolean;
  joback_success: boolean;
}

interface SpeciesRow {
  id: string;
  coefficient: string;
  query: string;
  selected: NistCompound | null;
  suggestions: NistCompound[];
  showSuggestions: boolean;
  quantity: string;
}

interface HeatSegment {
  label: string;
  T1_K: number;
  T2_K: number;
  dH_kjmol: number;
  type: 'sensible' | 'latent';
  phase?: string;
  warning?: string;
}

interface CalcResult {
  dH_rxn_kjmol: number;
  sensibleSegments: {
    nist_id: string;
    compound: NistCompound;
    speciesName: string;
    formula: string;
    coefficient: number;
    side: 'reactant' | 'product';
    segments: HeatSegment[];
    total_kjmol: number;
    missingData: boolean;
    approximated: boolean; // true when Joback estimates were used for any value
  }[];
  dH_total_kjmol: number;
  excessSegments: {
    nist_id: string;
    compound: NistCompound;
    speciesName: string;
    formula: string;
    excessCoefficient: number;
    stoichCoefficient: number; // raw stoichiometric coeff; used for live slider re-scaling
    segments: HeatSegment[];
    total_kjmol: number;
    missingData: boolean;
    approximated: boolean;
  }[];
  contributions: {
    formula: string;
    name: string;
    coefficient: number;
    hof_kjmol: number;
    contribution_kjmol: number;
    side: 'reactant' | 'product';
    phase298?: 'gas' | 'liquid' | 'solid';
  }[];
  equation: string;
  isExothermic: boolean;
  T1_K: number | null;
  T2_K: number | null;
  hasSensibleHeat: boolean;
  missingHessData: string[];
  missingSensibleData: boolean;
  limitingFormula: string | null;
  adiabaticFlameT_K: number | null;
  dGf_rxn_kjmol: number | null;
}

// ─── Constants ────────────────────────────────────────────────────────────────

const UNIT_CONVERSIONS: Record<string, { label: string; factor: number }> = {
  'kJ/mol':    { label: 'kJ/mol',     factor: 1 },
  'kcal/mol':  { label: 'kcal/mol',   factor: 0.239006 },
  'J/mol':     { label: 'J/mol',      factor: 1000 },
  'BTU/lbmol': { label: 'BTU/lbmol', factor: 429.923 },
};

const UNIT_CONVERSIONS_KG: Record<string, { label: string; factor: number }> = {
  'kJ/kg':   { label: 'kJ/kg',   factor: 1 },
  'J/kg':    { label: 'J/kg',    factor: 1000 },
  'kcal/kg': { label: 'kcal/kg', factor: 0.239006 },
  'BTU/lb':  { label: 'BTU/lb',  factor: 0.429923 },
};

// Absolute (total quantity) unit conversions — no /mol or /kg suffix
const UNIT_CONVERSIONS_ABS: Record<string, { label: string; factor: number }> = {
  'kJ':   { label: 'kJ',   factor: 1 },
  'kcal': { label: 'kcal', factor: 0.239006 },
  'J':    { label: 'J',    factor: 1000 },
  'BTU':  { label: 'BTU',  factor: 0.947817 },
};

// ─── Pure Helpers ─────────────────────────────────────────────────────────────

let _idCounter = 0;
function makeId() { return `sr-${++_idCounter}`; }
function makeSpecies(): SpeciesRow {
  return { id: makeId(), coefficient: '1', query: '', selected: null, suggestions: [], showSuggestions: false, quantity: '' };
}

function fmtNum(val: number): string {
  if (val === 0) return '0';
  const abs = Math.abs(val);
  if (abs < 1e-4) return val.toExponential(3);
  return parseFloat(val.toPrecision(4)).toString();
}

// ─── Formula Rendering ────────────────────────────────────────────────────────
const ELEMENTS_2 = new Set(['He','Li','Be','Ne','Na','Mg','Al','Si','Cl','Ar','Ca','Sc','Ti','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr']);
// Single-char element symbols — used to detect ambiguous pairs like "CO" (C then O, not Co cobalt)
const ALL_ELEMENTS_1 = new Set(['H','B','C','N','O','F','P','S','K','V','W','Y','I','U']);
/**
 * Normalise a chemical formula string:
 *  - Strip NIST artifact suffixes like "-N3" at the end.
 *  - Fix capitalisation (e.g. "HBR" → "HBr").
 *  - Never merge two consecutive uppercase letters into a 2-char element when the
 *    second character is a valid single-letter element (e.g. "CO2" stays "CO2",
 *    not "Co2"), because C and O are both valid 1-letter elements.
 */
function normalizeFormula(raw: string): string {
  // Strip NIST artifact suffixes: '-' followed by optional uppercase+lowercase+digits (e.g. "-N3", "-D2")
  // OR '-' followed by only digits (e.g. "-2" in isobutyl acetate C6H12O2-2)
  const trimmed = raw.replace(/-[A-Z][a-z]?\d*$/, '').replace(/-\d+$/, '');
  let out = ''; let i = 0;
  while (i < trimmed.length) {
    const ch = trimmed[i];
    if (/[A-Za-z]/.test(ch)) {
      const up = ch.toUpperCase();
      if (i + 1 < trimmed.length && /[A-Za-z]/.test(trimmed[i+1])) {
        const two = up + trimmed[i+1].toLowerCase();
        // Use 2-char element only when it is a known symbol AND either:
        //   (a) the second char is already lowercase (properly written, e.g. "Fe"), or
        //   (b) the second char alone is NOT a 1-letter element (e.g. "BR" → Br, since R is not an element)
        if (ELEMENTS_2.has(two) && (!ALL_ELEMENTS_1.has(trimmed[i+1].toUpperCase()) || /[a-z]/.test(trimmed[i+1]))) {
          out += two; i += 2; continue;
        }
      }
      out += up; i += 1;
    } else { out += ch; i += 1; }
  }
  return out;
}
/** Render a formula string as React nodes with <sub> for numbers */
function renderFormula(raw: string): React.ReactNode {
  const f = normalizeFormula(raw);
  const parts: React.ReactNode[] = []; let i = 0, k = 0;
  while (i < f.length) {
    if (/\d/.test(f[i])) {
      let n = ''; while (i < f.length && /\d/.test(f[i])) n += f[i++];
      parts.push(<sub key={k++} style={{ fontSize: '0.78em', lineHeight: 1 }}>{n}</sub>);
    } else {
      let s = ''; while (i < f.length && !/\d/.test(f[i])) s += f[i++];
      parts.push(<span key={k++}>{s}</span>);
    }
  }
  return <>{parts}</>;
}
/** ECharts axis label with rich-text subscript marker so digits render in the same font */
function echartsFormulaRich(raw: string): string {
  return normalizeFormula(raw).replace(/(\d+)/g, n => `{sub|${n}}`);
}

function toKelvin(val: number, unit: 'C' | 'K'): number {
  return unit === 'C' ? val + 273.15 : val;
}

/**
 * Compute enthalpy difference ΔH [kJ/mol] for a single Shomate/DIPPR segment.
 * eq107  (Aly-Lee gas):    Cp = A + B*(C/T/sinh(C/T))^2 + D*(E/T/cosh(E/T))^2  [J/kmol/K]
 *                           H(T) = A*T + B*C*coth(C/T) - D*E*tanh(E/T)          [J/kmol]
 * eq100  (polynomial solid/liq): Cp = A + B*T + C*T^2 + D*T^3 + E*T^4          [J/kmol/K]
 *                           H(T) = A*T + B*T^2/2 + … + E*T^5/5                  [J/kmol]
 * eq16   (liquid near Tc): Cp = A + B*T + C*T^2 + D*T^3 + E/(Tc-T)^g           [J/kmol/K]
 *                           H(T) antiderivative uses 1/((g-1)*(Tc-T)^(g-1))
 */
function dipprCpDeltaH(sh: ShomateRow, T1_K: number, T2_K: number): number {
  const { eqno } = sh;
  const a = sh.a ?? 0, b = sh.b ?? 0, c = sh.c ?? 0, d = sh.d ?? 0, e = sh.e ?? 0;

  if (eqno === 107) {
    // Aly-Lee: integral H(T) = A*T + B*C*coth(C/T) - D*E*tanh(E/T)  [J/kmol]
    const H = (T: number) => {
      const bTerm = (b !== 0 && c !== 0) ? b * c * (1 / Math.tanh(c / T)) : 0;
      const dTerm = (d !== 0 && e !== 0) ? d * e * Math.tanh(e / T) : 0;
      return a * T + bTerm - dTerm;
    };
    return (H(T2_K) - H(T1_K)) / 1e6; // J/kmol → kJ/mol

  } else if (eqno === 100) {
    // Standard polynomial
    const H = (T: number) => a*T + b*T**2/2 + c*T**3/3 + d*T**4/4 + e*T**5/5;
    return (H(T2_K) - H(T1_K)) / 1e6;

  } else if (eqno === 16 && sh.f != null) {
    // Liquid Cp with singular term near Tc: Cp = A+BT+CT^2+DT^3 + E/(Tc-T)^g
    const Tc = sh.f;
    const gv = sh.g ?? 3; // exponent (integer ≥ 1)
    const polyH = (T: number) => a*T + b*T**2/2 + c*T**3/3 + d*T**4/4;
    const singH = (T: number) => {
      const dT = Tc - T;
      if (gv === 1) return -e * Math.log(Math.abs(Math.max(dT, 0.2))); // ∫e/(Tc-T)dT = -e·ln|Tc-T|
      if (Math.abs(dT) < 0.2) return e / ((gv - 1) * (0.2 ** (gv - 1)));
      return e / ((gv - 1) * (dT ** (gv - 1)));
    };
    return ((polyH(T2_K) + singH(T2_K)) - (polyH(T1_K) + singH(T1_K))) / 1e6;

  } else {
    // Fallback: treat as polynomial
    const H = (T: number) => a*T + b*T**2/2 + c*T**3/3 + d*T**4/4 + e*T**5/5;
    return (H(T2_K) - H(T1_K)) / 1e6;
  }
}

/**
 * Integrate Joback group-contribution Cp polynomial from T1 to T2.
 * Cp [J/mol·K] = (Σa−37.93) + (Σb+0.210)T + (Σc−3.91×10⁻⁴)T² + (Σd+2.06×10⁻⁷)T³
 * Returns kJ/mol.
 */
function jobackCpDeltaH(j: JobackEstimate, T1_K: number, T2_K: number): number {
  // Clamp to Joback's validated Cp range (typically 300–1000 K) to avoid polynomial blow-up.
  const tMin = j.cp_t_min_k ?? 300;
  const tMax = j.cp_t_max_k ?? 1000;
  const t1 = Math.max(T1_K, tMin);
  const t2 = Math.min(T2_K, tMax);
  if (t2 <= t1) return 0;
  const A = (j.cp_a ?? 0) - 37.93;
  const B = (j.cp_b ?? 0) + 0.210;
  const C = (j.cp_c ?? 0) - 3.91e-4;
  const D = (j.cp_d ?? 0) + 2.06e-7;
  const dH = A * (t2 - t1)
    + B / 2 * (t2 * t2 - t1 * t1)
    + C / 3 * (t2 ** 3 - t1 ** 3)
    + D / 4 * (t2 ** 4 - t1 ** 4);
  return dH / 1000; // J/mol → kJ/mol
}

function findShomate(params: ShomateRow[], phase: string, T_K: number): ShomateRow | null {
  const dist = (v: ShomateRow) => v.t_max_k < T_K ? T_K - v.t_max_k : v.t_min_k > T_K ? v.t_min_k - T_K : 0;
  const closest = (rows: ShomateRow[]) => {
    const exact = rows.find(s => T_K >= s.t_min_k - 5 && T_K <= s.t_max_k + 5);
    return exact ?? rows.reduce((best, sh) => dist(sh) < dist(best) ? sh : best);
  };
  const byPhase = params.filter(s => s.phase === phase);
  if (byPhase.length > 0) return closest(byPhase);
  // Phase transition data missing (e.g. CO2 sublimes — no tboil_k) — fall back to gas › liquid › solid
  for (const fp of ['gas', 'liquid', 'solid'] as const) {
    if (fp === phase) continue;
    const fb = params.filter(s => s.phase === fp);
    if (fb.length > 0) return closest(fb);
  }
  return null;
}

function phaseAt(T_K: number, tfus: number | null, tboil: number | null): 'solid' | 'liquid' | 'gas' {
  if (tfus !== null && T_K < tfus) return 'solid';
  if (tboil !== null && T_K > tboil) return 'gas';
  return 'liquid';
}

/**
 * Compute the standard enthalpy of formation for the STABLE PHASE of a compound at T_ref=298.15 K.
 * The database always stores the ideal-gas value (hf_gas_298_kjmol). For species that are liquid
 * or solid at 298 K (e.g. H₂O, benzene), this differs from the stable-phase ΔHf by ΔHvap/ΔHfus.
 * Using a gas value WITH a sensible-heat path that crosses phase transitions double-counts the
 * phase-change enthalpy, so we correct here using Kirchhoff's law.
 *
 * hf_liq(T_ref) = hf_gas + ∫[T_ref→Tboil]Cp_gas − ΔHvap(Tboil) − ∫[T_ref→Tboil]Cp_liq
 */
function effectiveHf(
  compound: NistCompound,
  shParams: ShomateRow[],
  pt: PhaseTransitionRow | null,
  joback?: JobackEstimate | null,
): number {
  const hf_gas = compound.hf_gas_298_kjmol ?? 0;
  if (!shParams.length) return hf_gas;
  const T_ref = 298.15;
  const tboil = pt?.tboil_k ?? compound.tboil_k ?? null;
  const tfus  = pt?.tfus_k  ?? compound.tfus_k  ?? null;
  const hvap  = pt?.hvap_kjmol ?? (joback?.hvap_kjmol ?? null);
  const hfus  = pt?.hfus_kjmol ?? (joback?.hfus_kjmol ?? null);
  // Solid correction requires both hvap AND hfus; if Joback solid_correction_incomplete is set,
  // the group-contribution estimate may be partial — we still use it but mark hf_gas fallback.
  const solidIncomplete = joback?.solid_correction_incomplete ?? false;
  const ph = phaseAt(T_ref, tfus, tboil);
  if (ph === 'gas') return hf_gas; // already gas at T_ref — ideal-gas ΔHf is correct

  if (ph === 'liquid') {
    if (hvap === null || tboil === null) return hf_gas;
    const mid = (T_ref + tboil) / 2;
    const sh_gas = findShomate(shParams, 'gas',    mid);
    const sh_liq = findShomate(shParams, 'liquid', mid);
    if (!sh_gas || !sh_liq) return hf_gas;
    // Kirchhoff path: gas(T_ref)→gas(Tboil)→condense→liquid(Tboil)→liquid(T_ref)
    return hf_gas
      + dipprCpDeltaH(sh_gas, T_ref, tboil)  // heat gas T_ref → Tboil
      - hvap                                  // condense at Tboil
      - dipprCpDeltaH(sh_liq, T_ref, tboil); // cool liquid Tboil → T_ref (sign: subtract heating)
  }

  if (ph === 'solid') {
    if (solidIncomplete || hvap === null || hfus === null || tboil === null || tfus === null) return hf_gas;
    const midBoil = (T_ref + tboil) / 2;
    const midFus  = (Math.max(T_ref, tfus) + tboil) / 2;
    const sh_gas = findShomate(shParams, 'gas',    midBoil);
    const sh_liq = findShomate(shParams, 'liquid', midFus);
    const sh_sol = findShomate(shParams, 'solid',  (T_ref + Math.min(tfus, tboil)) / 2);
    if (!sh_gas) return hf_gas;
    const hf_liq_Tref = hf_gas
      + dipprCpDeltaH(sh_gas, T_ref, tboil)                           // gas: T_ref → Tboil
      - hvap                                                           // condense at Tboil
      + (sh_liq ? dipprCpDeltaH(sh_liq, tboil, tfus) : 0);           // cool liquid: Tboil → Tfus
    const hfus_at_Tref = hfus
      + (sh_liq ? dipprCpDeltaH(sh_liq, tfus, T_ref) : 0)            // cool liquid: Tfus → T_ref
      - (sh_sol ? dipprCpDeltaH(sh_sol, tfus, T_ref) : 0);           // cool solid:  Tfus → T_ref
    return hf_liq_Tref - hfus_at_Tref;
  }
  return hf_gas;
}

function computeSpeciesSensibleHeat(
  shParams: ShomateRow[],
  pt: PhaseTransitionRow | null,
  compound: NistCompound,
  T1_K: number,
  T2_K: number,
  joback?: JobackEstimate | null,
): { segments: HeatSegment[]; total_kjmol: number; missingData: boolean; approximated: boolean } {
  const tboil = pt?.tboil_k ?? compound.tboil_k ?? null;
  const tfus  = pt?.tfus_k  ?? compound.tfus_k  ?? null;
  // Joback fallback for latent heats when not in nist_phase_transitions
  const hvap  = pt?.hvap_kjmol ?? (joback?.hvap_kjmol ?? null);
  const hfus  = pt?.hfus_kjmol ?? (joback?.hfus_kjmol ?? null);
  const hvapIsJoback = pt?.hvap_kjmol == null && hvap != null;
  const hfusIsJoback = pt?.hfus_kjmol == null && hfus != null;
  // Describe the hvap source for the warning tooltip
  const hvapRef = joback?.hvap_reference;  // 'watson_298K' | 'joback_Tb' | null
  const hvapWarn = !hvapIsJoback ? undefined
    : hvapRef === 'watson_298K'
      ? '\u0394Hvap from Joback (Watson-corrected to 298 K)'
      : '\u0394Hvap from Joback estimate (at T\u1d47 — ~5–10% uncertainty)';
  // Check Joback Cp range against the segment temperatures
  const cpTmax = joback?.cp_t_max_k ?? 1000;

  const tMin = Math.min(T1_K, T2_K);
  const tMax = Math.max(T1_K, T2_K);
  const direction = T2_K >= T1_K ? 1 : -1;

  const checkpoints: number[] = [tMin];
  if (tfus  !== null && tfus  > tMin + 0.5 && tfus  < tMax - 0.5) checkpoints.push(tfus);
  if (tboil !== null && tboil > tMin + 0.5 && tboil < tMax - 0.5) checkpoints.push(tboil);
  checkpoints.push(tMax);
  checkpoints.sort((a, b) => a - b);

  const segments: HeatSegment[] = [];
  let missingData = false;
  let approximated = false;
  const hasJobackCp = !!(joback?.joback_success && joback.cp_a != null && joback.cp_b != null && joback.cp_c != null && joback.cp_d != null);
  // Warn when any part of the requested T range is outside the Joback Cp polynomial range
  const jobackCpOutOfRange = (tLo: number, tHi: number) => hasJobackCp && (tHi > cpTmax + 1 || tLo < (joback?.cp_t_min_k ?? 300) - 1);

  for (let i = 0; i < checkpoints.length - 1; i++) {
    const segT1 = checkpoints[i];
    const segT2 = checkpoints[i + 1];
    const midT  = (segT1 + segT2) / 2;
    const phase = phaseAt(midT, tfus, tboil);
    const sh1 = findShomate(shParams, phase, segT1);
    const sh2 = findShomate(shParams, phase, segT2);

    let dH = 0;
    let warning: string | undefined;

    if (sh1 && sh2) {
      if (sh1.t_min_k === sh2.t_min_k) {
        dH = dipprCpDeltaH(sh1, segT1, segT2);
      } else {
        const splitT = sh1.t_max_k;
        dH = dipprCpDeltaH(sh1, segT1, splitT) + dipprCpDeltaH(sh2, splitT, segT2);
      }
    } else if (sh1) {
      const clampT2 = Math.min(segT2, sh1.t_max_k);
      dH = dipprCpDeltaH(sh1, segT1, clampT2);
      if (clampT2 < segT2 - 1) {
        if (hasJobackCp) {
          dH += jobackCpDeltaH(joback!, clampT2, segT2);
          warning = `Cp extended ${clampT2.toFixed(0)}–${segT2.toFixed(0)} K via Joback estimate`;
          approximated = true;
        } else {
          warning = `Cp(T) data ends at ${sh1.t_max_k.toFixed(0)} K — truncated`;
          missingData = true;
        }
      }
    } else {
      if (hasJobackCp) {
        dH = jobackCpDeltaH(joback!, segT1, segT2);
        warning = jobackCpOutOfRange(segT1, segT2)
          ? `Cp from Joback estimate (${phase}, capped at ${cpTmax.toFixed(0)} K)`
          : `Cp from Joback estimate (${phase})`;
        approximated = true;
      } else {
        warning = `No Cp(T) data for ${phase} phase — assumed 0`;
        missingData = true;
      }
    }

    segments.push({
      label: `${phase.charAt(0).toUpperCase() + phase.slice(1)} (${segT1.toFixed(0)}\u2013${segT2.toFixed(0)} K)`,
      T1_K: direction < 0 ? segT2 : segT1, T2_K: direction < 0 ? segT1 : segT2, dH_kjmol: dH * direction,
      type: 'sensible', phase, warning,
    });

    if (i < checkpoints.length - 2) {
      const transT = checkpoints[i + 1];
      if (tfus !== null && Math.abs(transT - tfus) < 1) {
        segments.push({
          label: `Fusion at ${tfus.toFixed(1)} K`,
          T1_K: transT, T2_K: transT,
          dH_kjmol: hfus !== null ? hfus * direction : 0,
          type: 'latent',
          warning: hfus === null
            ? '\u0394Hfus not in database \u2014 assumed 0'
            : hfusIsJoback ? '\u0394Hfus from Joback estimate' : undefined,
        });
        if (hfus === null) missingData = true;
        if (hfusIsJoback) approximated = true;
      }
      if (tboil !== null && Math.abs(transT - tboil) < 1) {
        segments.push({
          label: `Vaporization at ${tboil.toFixed(1)} K`,
          T1_K: transT, T2_K: transT,
          dH_kjmol: hvap !== null ? hvap * direction : 0,
          type: 'latent',
          warning: hvap === null
            ? '\u0394Hvap not in database \u2014 assumed 0'
            : hvapIsJoback ? hvapWarn : undefined,
        });
        if (hvap === null) missingData = true;
        if (hvapIsJoback) approximated = true;
      }
    }
  }

  return { segments, total_kjmol: segments.reduce((s, x) => s + x.dH_kjmol, 0), missingData, approximated };
}
// ─── Name Formatting ──────────────────────────────────────────────────────────
const FILLER_WORDS = new Set(['and','of','the','in','at','on','a','an','with','for','to','by','from','into','per','bis','tris']);
// IUPAC locant prefixes that should stay lower-case (e.g. n-Butane, iso-Propanol, trans-2-Butene)
const IUPAC_LOWERCASE_PREFIXES = /^(n|iso|sec|tert|neo|cis|trans|alpha|beta|gamma|delta|ortho|meta|para|di|tri|tetra|penta|hexa|hepta|octa|bis|tris|tetrakis)-/;
function formatCompoundName(raw: string): string {
  return raw.trim().replace(/\s+/g, ' ').split(' ').map((word, i) => {
    if (/^\([a-z]+\)$/.test(word)) return word; // (s), (g), (aq) unchanged
    const lower = word.toLowerCase();
    // IUPAC prefix like "n-", "iso-", "trans-" stays lowercase; capitalise first letter after it
    const prefixMatch = lower.match(IUPAC_LOWERCASE_PREFIXES);
    if (prefixMatch) {
      const prefix = prefixMatch[0];
      const rest   = lower.slice(prefix.length);
      // Capitalise the first alphabetic character after the prefix (skipping digits/commas)
      return prefix + rest.replace(/^([^a-zA-Z]*)([a-zA-Z])/, (_, pre, ch) => pre + ch.toUpperCase());
    }
    if (i === 0 || !FILLER_WORDS.has(lower)) return lower.charAt(0).toUpperCase() + lower.slice(1);
    return lower;
  }).join(' ');
}

export default function ReactionThermochemistryPage() {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';

  const [reactants,   setReactants]   = useState<SpeciesRow[]>([makeSpecies()]);
  const [products,    setProducts]    = useState<SpeciesRow[]>([makeSpecies()]);
  const [unit,        setUnit]        = useState<string>('kJ/mol');
  const [heatUnit,    setHeatUnit]    = useState<'kJ/mol' | 'kJ/kg'>('kJ/mol');
  const [tempUnit,    setTempUnit]    = useState<'C' | 'K'>('C');
  const [enableTemp,  setEnableTemp]  = useState(false);
  const [t1Str,       setT1Str]       = useState('30');
  const [t2Str,       setT2Str]       = useState('100');
  const [result,      setResult]      = useState<CalcResult | null>(null);
  const [calcError,   setCalcError]   = useState<string | null>(null);
  const [loading,     setLoading]     = useState(false);
  const [conversion,  setConversion]  = useState(100);
  const [mounted,     setMounted]     = useState(false);
  const [kirchhoffOpen, setKirchhoffOpen] = useState(true);
  const [phaseChangeNotice, setPhaseChangeNotice] = useState<{formula: string; atK: number; transition: 'vap' | 'fus'}[] | null>(null);
  const [sensMode, setSensMode] = useState(false);
  const [sensVar, setSensVar] = useState<'T_in' | 'T_out' | 'conversion' | null>(null);
  const [sensOpts, setSensOpts] = useState<EChartsOption | null>(null);
  const [sensRunning, setSensRunning] = useState(false);
  const [pendingAdiabaticCalc, setPendingAdiabaticCalc] = useState(false);
  const [eqConvOpen, setEqConvOpen] = useState(false);
  const [perMoleMode, setPerMoleMode] = useState(true); // true = per mol/kg basis; false = absolute (total heat for given quantity)
  // Quantities committed on Enter/Calculate — so live typing doesn't update the chart
  const [committedQtys, setCommittedQtys] = useState<Map<string, string>>(new Map());
  // Confirmed balance snapshot — updated only when user explicitly confirms (Enter or
  // Calculate button). Null = no confirmation yet. Only non-null state shows warnings.
  const [confirmedBalanceCheck, setConfirmedBalanceCheck] = useState<BalanceResult | null>(null);
  const confirmBalanceRef = useRef(false); // set before calculate(); cleared inside it
  const [pendingAutoCalc, setPendingAutoCalc] = useState(false);
  useEffect(() => { setMounted(true); }, []);

  // Clear stale balance-check banner whenever the set of selected species changes.
  // This prevents showing an old "fixable" result after the user swaps a compound.
  const speciesKey = useMemo(() => {
    const r = reactants.filter(x => x.selected).map(x => x.selected!.formula).sort().join(',');
    const p = products .filter(x => x.selected).map(x => x.selected!.formula).sort().join(',');
    return `${r}\u2192${p}`;
  }, [reactants, products]);
  useEffect(() => { setConfirmedBalanceCheck(null); }, [speciesKey]);

  // Auto-dismiss phase change notice after 8 seconds
  useEffect(() => {
    if (!phaseChangeNotice || phaseChangeNotice.length === 0) return;
    const timer = setTimeout(() => setPhaseChangeNotice(null), 8000);
    return () => clearTimeout(timer);
  }, [phaseChangeNotice]);

  // Re-apply theme-appropriate colors to sensitivity chart when switching modes
  useEffect(() => {
    const chart = sensEchartsRef.current?.getEchartsInstance?.();
    if (!chart || !sensOpts) return;
    const textColor = isDark ? '#e2e8f0' : '#1a202c';
    const tooltipBg = isDark ? '#1e293b' : '#ffffff';
    const tooltipBorder = isDark ? '#475569' : '#e2e8f0';
    const axisLineColor = isDark ? '#e2e8f0' : '#1a202c';
    const splitLineColor = isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.08)';
    chart.setOption({
      tooltip: {
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        extraCssText: `box-shadow:0 4px 12px rgba(0,0,0,${isDark ? '0.4' : '0.1'});`,
        textStyle: { color: textColor },
      },
      xAxis: { nameTextStyle: { color: textColor }, axisLabel: { color: textColor }, axisLine: { lineStyle: { color: axisLineColor } } },
      yAxis: { nameTextStyle: { color: textColor }, axisLabel: { color: textColor }, axisLine: { lineStyle: { color: axisLineColor } }, splitLine: { lineStyle: { color: splitLineColor } } },
      title: { textStyle: { color: textColor, rich: { nm: { color: textColor }, sb: { color: textColor } } } },
    }, false);
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [isDark]);

  // After each calculation the chart re-renders; measure the actual ECharts grid rect
  // so the custom JSX legend is precisely centred over the bar plot area.
  useEffect(() => {
    if (!result) return;
    try {
      const inst = echartsRef.current?.getEchartsInstance?.();
      if (!inst) return;
      const rect = inst.getModel?.()?.getComponent?.('grid', 0)?.coordinateSystem?.getRect?.();
      const w = inst.getWidth?.();
      if (rect && w && w > 0) {
        setLegendOffset({ left: Math.round(rect.x), right: Math.round(w - rect.x - rect.width) });
      }
    } catch { /* ECharts not yet rendered */ }
  }, [result]);

  const wrapperRefs  = useRef<Map<string, HTMLDivElement>>(new Map());
  const searchTimers = useRef<Map<string, ReturnType<typeof setTimeout>>>(new Map());
  const autoCalcInit = useRef(false);
  const autoCalcDone = useRef(false);
  // Persistent Shomate+PT cache — keyed by nist_id, never evicted (data is immutable in DB)
  const shomatePersistCache = useRef<Map<string, { shParams: ShomateRow[]; pt: PhaseTransitionRow | null }>>(new Map());
  // Persistent Joback estimates cache — null value means "no Joback entry for this nist_id"
  const jobackPersistCache = useRef<Map<string, JobackEstimate | null>>(new Map());
  // ECharts instance ref — used to read exact grid rect for legend centering
  const echartsRef = useRef<any>(null);
  const sensEchartsRef = useRef<any>(null);
  const sliderCalcRef = useRef(false);
  const [legendOffset, setLegendOffset] = useState({ left: 64, right: 20 });

  /** Max temperature (K) with reliable Cp data for a species (from Shomate or Joback). */
  const getMaxCpT = useCallback((nistId: string): number => {
    const entry = shomatePersistCache.current.get(nistId);
    const shoParams = entry?.shParams ?? [];
    const maxSho = shoParams.length > 0 ? Math.max(...shoParams.map(s => s.t_max_k)) : 0;
    const jb = jobackPersistCache.current.get(nistId);
    const jbMax = jb?.cp_t_max_k ?? 0;
    return Math.max(maxSho, jbMax);
  }, []);

  const perKgFactor = useMemo(() => {
    if (heatUnit !== 'kJ/kg') return 1;
    const r = reactants.find(row => row.selected?.mol_weight);
    return r?.selected?.mol_weight ? 1000 / r.selected.mol_weight : 1;
  }, [heatUnit, reactants]);



  useEffect(() => {
    function onMouseDown(e: MouseEvent) {
      wrapperRefs.current.forEach((el, id) => {
        if (el && !el.contains(e.target as Node)) {
          setReactants(prev => prev.map(r => r.id === id ? { ...r, showSuggestions: false } : r));
          setProducts(prev  => prev.map(r => r.id === id ? { ...r, showSuggestions: false } : r));
        }
      });
    }
    document.addEventListener('mousedown', onMouseDown);
    return () => document.removeEventListener('mousedown', onMouseDown);
  }, []);

  // Preload CH4 combustion example with temperature range enabled (shows sensible + latent)
  useEffect(() => {
    const ids = ['C74828', 'C7782447', 'C124389', 'C7732185'];
    supabase
      .from('nist_compounds')
      .select('nist_id,name,formula,cas,mol_weight,hf_gas_298_kjmol,cp_gas_298_jmolk,tboil_k,tfus_k,tc_k,has_shomate,has_phase_change')
      .in('nist_id', ids)
      .then(({ data }) => {
        if (!data || data.length === 0) return;
        const byId = new Map<string, NistCompound>((data as NistCompound[]).map(c => [c.nist_id, c]));
        const makeRow = (id: string, coeff: string): SpeciesRow => {
          const c = byId.get(id) ?? null;
          return { id: makeId(), coefficient: coeff, query: c ? formatCompoundName(c.name) : '', selected: c, suggestions: [], showSuggestions: false, quantity: '' };
        };
        const r1 = makeRow('C74828', '1'); r1.quantity = '1';
        const r2 = makeRow('C7782447', '2'); r2.quantity = '2';
        setReactants([r1, r2]);
        setProducts ([makeRow('C124389', '1'), makeRow('C7732185', '2')]);
        setEnableTemp(true);
        setT1Str('30');
        setT2Str('150');
      });
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const patchRow = useCallback((side: 'reactants' | 'products', id: string, patch: Partial<SpeciesRow>) => {
    (side === 'reactants' ? setReactants : setProducts)(prev => prev.map(r => r.id === id ? { ...r, ...patch } : r));
  }, []);

  const fetchSuggestions = useCallback(async (side: 'reactants' | 'products', id: string, q: string) => {
    if (q.trim().length < 1) { patchRow(side, id, { suggestions: [], showSuggestions: false }); return; }
    const normKey = (s: string) => s.normalize('NFD').replace(/[\u0300-\u036f]/g, '').replace(/[.'\u2011]/g, '').toLowerCase().trim();
    const cols = 'nist_id,name,formula,cas,mol_weight,hf_gas_298_kjmol,cp_gas_298_jmolk,tboil_k,tfus_k,tc_k,has_shomate,has_phase_change';
    // Include compounds with null hf_gas_298_kjmol — elements in standard state (O2, N2, H2…)
    // store null because their ΔHf = 0 by convention; processHess handles null as 0.
    const base = () => supabase.from('nist_compounds').select(cols)
      .eq('has_shomate', true).eq('has_phase_change', true);
    const [prefixRes, containsRes] = await Promise.all([
      base().or(`name.ilike.${q}%,formula.ilike.${q}%`).limit(30),
      base().or(`name.ilike.%${q}%,formula.ilike.%${q}%`).limit(50),
    ]);
    const seen = new Set<string>();
    const all: NistCompound[] = [];
    for (const d of [...(prefixRes.data ?? []), ...(containsRes.data ?? [])]) {
      const c = d as NistCompound;
      if (!seen.has(c.nist_id)) { seen.add(c.nist_id); all.push(c); }
    }
    const normQ = normKey(q);
    // NIST names use inverted-CAS style: "BUTANE,-1,1-DIBROMO-" (parent, modifiers).
    // Strip common IUPAC/trivial prefixes (n-, iso-, sec-, cyclo-, etc.) from the base
    // name (before the first comma) so "n-butane" scores at the same tier as "butane"
    // but ranks above "butane,-1,1-dibromo-" by shorter total name length.
    const IUPAC_PREFIX = /^(n|iso|sec|tert|di|tri|tetra|neo|cyclo|bi|trans|cis|alpha|beta|gamma|delta|ortho|meta|para)-/;
    const score = (c: NistCompound): number => {
      const nn = normKey(c.name);
      const nf = normKey(c.formula);
      // Base name = text before first NIST comma-modifier (e.g. "butane" from "butane,-1,1-dibromo-")
      const baseName = nn.split(',')[0].trim();
      // Strip leading IUPAC prefix: "n-butane" → "butane", "iso-butane" → "butane"
      const coreName = baseName.replace(IUPAC_PREFIX, '');
      const coreWords = coreName.split(/[\s\-]+/).filter(Boolean);
      const allWords = nn.split(/[\s,()[\]\-]+/).filter(Boolean);

      if (nf === normQ) return 600;                                           // exact formula
      if (nn === normQ) return 500;                                           // exact full name
      if (nf.startsWith(normQ)) return 400 - nf.length;                      // formula prefix
      if (coreName === normQ) return 490 - nn.length;                        // core name exact
      if (coreName.startsWith(normQ)) return 300 - nn.length;               // core name prefix
      if (coreWords.some(w => w === normQ)) return 250 - nn.length;         // core word exact
      if (coreWords.some(w => w.startsWith(normQ))) return 200 - nn.length; // core word prefix
      if (allWords.some(w => w === normQ)) return 150 - nn.length;          // any word exact
      if (allWords.some(w => w.startsWith(normQ))) return 100 - nn.length;  // any word prefix
      return 10;
    };
    const sorted = all.sort((a, b) => score(b) - score(a) || a.formula.length - b.formula.length || a.name.length - b.name.length).slice(0, 12);
    patchRow(side, id, { suggestions: sorted, showSuggestions: true });
  }, [patchRow]);

  const handleQueryChange = useCallback((side: 'reactants' | 'products', id: string, value: string) => {
    patchRow(side, id, { query: value, selected: null, showSuggestions: false });
    const existing = searchTimers.current.get(id);
    if (existing) clearTimeout(existing);
    searchTimers.current.set(id, setTimeout(() => fetchSuggestions(side, id, value), 220));
  }, [patchRow, fetchSuggestions]);

  const handleSelect = useCallback((side: 'reactants' | 'products', id: string, row: NistCompound) => {
    patchRow(side, id, { selected: row, query: formatCompoundName(row.name), suggestions: [], showSuggestions: false });
  }, [patchRow]);

  const addRow = useCallback((side: 'reactants' | 'products') => {
    (side === 'reactants' ? setReactants : setProducts)(prev => [...prev, makeSpecies()]);
  }, []);

  const removeRow = useCallback((side: 'reactants' | 'products', id: string) => {
    (side === 'reactants' ? setReactants : setProducts)(prev => prev.filter(r => r.id !== id));
  }, []);

  const handleTempUnitToggle = useCallback((newUnit: string) => {
    const nu = newUnit as 'C' | 'K';
    const adj = nu === 'K' ? 273.15 : -273.15;
    const t1n = parseFloat(t1Str); const t2n = parseFloat(t2Str);
    setT1Str(!isNaN(t1n) ? String(+(t1n + adj).toFixed(2)) : t1Str);
    setT2Str(!isNaN(t2n) ? String(+(t2n + adj).toFixed(2)) : t2Str);
    setTempUnit(nu);
  }, [t1Str, t2Str]);

  // ─── Main Calculation ──────────────────────────────────────────────────────

  const calculate = useCallback(async () => {
    setCalcError(null);
    // Auto-select first suggestion for rows that have suggestions but no selection
    const applyAutoSelect = (rows: SpeciesRow[]): SpeciesRow[] =>
      rows.map(r => (!r.selected && r.suggestions.length > 0)
        ? { ...r, selected: r.suggestions[0], query: formatCompoundName(r.suggestions[0].name), suggestions: [], showSuggestions: false }
        : r);
    const resolvedReactants = applyAutoSelect(reactants);
    const resolvedProducts  = applyAutoSelect(products);
    if (resolvedReactants !== reactants) setReactants(resolvedReactants);
    if (resolvedProducts  !== products)  setProducts(resolvedProducts);
    const validReactants = resolvedReactants.filter(r => r.selected);
    const validProducts  = resolvedProducts.filter(p => p.selected);

    // Snapshot balance state if this run was triggered by Enter or Calculate button.
    // Auto-runs (initial load, slider effects) must NOT update this.
    if (confirmBalanceRef.current) {
      confirmBalanceRef.current = false;
      const rSpec = validReactants.map(r => ({ formula: r.selected!.formula, coefficient: Math.max(0.001, parseFloat(r.coefficient) || 1) }));
      const pSpec = validProducts .map(r => ({ formula: r.selected!.formula, coefficient: Math.max(0.001, parseFloat(r.coefficient) || 1) }));
      if (rSpec.length > 0 && pSpec.length > 0) setConfirmedBalanceCheck(checkBalance(rSpec, pSpec));
      else setConfirmedBalanceCheck(null);
    }

    if (validReactants.length === 0 || validProducts.length === 0) {
      setCalcError('Please select at least one reactant and one product from the autocomplete dropdown.');
      return;
    }
    setLoading(true);
    // Commit current qty inputs so chart + scale use the same snapshot calculate() sees
    setCommittedQtys(new Map(reactants.map(r => [r.id, r.quantity])));
    try {
      const contributions: CalcResult['contributions'] = [];
      const missingHessData: string[] = [];

      const processHess = (rows: SpeciesRow[], side: 'reactant' | 'product') => {
        for (const row of rows) {
          if (!row.selected) continue;
          const coeff = parseFloat(row.coefficient);
          if (isNaN(coeff) || coeff <= 0) continue;
          // Elements in standard state have ΔH°f = 0 by convention (NIST stores null for them)
          const hof = row.selected.hf_gas_298_kjmol ?? 0;
          contributions.push({
            formula: row.selected.formula, name: row.selected.name, coefficient: coeff,
            hof_kjmol: hof, contribution_kjmol: (side === 'product' ? 1 : -1) * coeff * hof, side,
            phase298: phaseAt(298.15, row.selected.tfus_k, row.selected.tboil_k),
          });
        }
      };
      processHess(validReactants, 'reactant');
      processHess(validProducts,  'product');

      if (contributions.length === 0) {
        setCalcError('None of the selected species have \u0394H\u00b0f(298 K) in the NIST database.');
        setLoading(false); return;
      }

      let dH_rxn_kjmol = contributions.reduce((s, c) => s + c.contribution_kjmol, 0);

      // ── Always fetch Shomate+PT (persistent cache) + apply effectiveHf correction ────
      // effectiveHf corrects NIST gas-phase ΔHf to stable-phase (liquid/solid) at 298 K.
      // This must run unconditionally — not only in sensible-heat mode — to avoid a
      // ~84 kJ/mol discontinuity when T_in≈T_out (equal temps skip sensible path but
      // still need the phase correction, e.g. 2×ΔHvap(H₂O) ≈ 88 kJ/mol jump).
      const allSpeciesForHf = [
        ...validReactants.map(r => ({ ...r, side: 'reactant' as const })),
        ...validProducts.map(r =>  ({ ...r, side: 'product'  as const })),
      ];
      // Fetch any species not yet in persistent cache
      await Promise.all(allSpeciesForHf.map(async (sp) => {
        if (!sp.selected || !sp.selected.has_shomate) return;
        if (shomatePersistCache.current.has(sp.selected.nist_id)) return; // already cached
        const nid = sp.selected.nist_id;
        const [shRes, ptRes, jbRes] = await Promise.all([
          supabase.from('nist_shomate').select('phase,eqno,t_min_k,t_max_k,a,b,c,d,e,f,g').eq('nist_id', nid),
          supabase.from('nist_phase_transitions').select('tboil_k,tfus_k,hvap_kjmol,hfus_kjmol').eq('nist_id', nid).maybeSingle(),
          jobackPersistCache.current.has(nid)
            ? Promise.resolve({ data: null, error: null })
            : supabase.from('nist_joback_estimates').select('nist_id,hf_gas_298_kjmol,dgf_gas_298_kjmol,hfus_kjmol,hvap_kjmol,hvap_reference,cp_a,cp_b,cp_c,cp_d,cp_t_min_k,cp_t_max_k,joback_group_count,solid_correction_incomplete,joback_success').eq('nist_id', nid).maybeSingle(),
        ]);
        shomatePersistCache.current.set(nid, {
          shParams: (shRes.data ?? []) as ShomateRow[],
          pt: (ptRes.data ?? null) as PhaseTransitionRow | null,
        });
        if (!jobackPersistCache.current.has(nid)) {
          jobackPersistCache.current.set(nid, (jbRes.data ?? null) as JobackEstimate | null);
        }
      }));
      // Apply effectiveHf correction to contributions + compute phase298
      for (const contrib of contributions) {
        const sp = allSpeciesForHf.find(r => r.selected?.formula === contrib.formula);
        if (!sp?.selected) continue;
        const cached = shomatePersistCache.current.get(sp.selected.nist_id);
        if (!cached) continue;
        const corrHof = effectiveHf(sp.selected, cached.shParams, cached.pt, jobackPersistCache.current.get(sp.selected.nist_id));
        contrib.hof_kjmol = corrHof;
        contrib.contribution_kjmol = (contrib.side === 'product' ? 1 : -1) * contrib.coefficient * corrHof;
        const tboil298 = cached.pt?.tboil_k ?? sp.selected.tboil_k ?? null;
        const tfus298  = cached.pt?.tfus_k  ?? sp.selected.tfus_k  ?? null;
        contrib.phase298 = phaseAt(298.15, tfus298, tboil298);
      }
      dH_rxn_kjmol = contributions.reduce((s, c) => s + c.contribution_kjmol, 0);

      // ΔG°rxn — only if ALL non-element species have Joback dgf_gas_298 coverage
      const _allHaveDgf = [...validReactants, ...validProducts].filter(sp => sp.selected).every(sp => {
        if (isMonoElemental(sp.selected!.formula)) return true;
        const jb = jobackPersistCache.current.get(sp.selected!.nist_id);
        return jb?.dgf_gas_298_kjmol != null;
      });
      let dGf_rxn_kjmol: number | null = null;
      if (_allHaveDgf) {
        let _dgf = 0;
        for (const sp of validReactants.filter(r => r.selected)) {
          const coeff = parseFloat(sp.coefficient); if (isNaN(coeff) || coeff <= 0) continue;
          const jb = jobackPersistCache.current.get(sp.selected!.nist_id);
          _dgf -= coeff * (isMonoElemental(sp.selected!.formula) ? 0 : (jb?.dgf_gas_298_kjmol ?? 0));
        }
        for (const sp of validProducts.filter(p => p.selected)) {
          const coeff = parseFloat(sp.coefficient); if (isNaN(coeff) || coeff <= 0) continue;
          const jb = jobackPersistCache.current.get(sp.selected!.nist_id);
          _dgf += coeff * (isMonoElemental(sp.selected!.formula) ? 0 : (jb?.dgf_gas_298_kjmol ?? 0));
        }
        dGf_rxn_kjmol = _dgf;
      }

      let adiabaticFlameT_K: number | null = null;
      let T1_K: number | null = null;
      let T2_K: number | null = null;
      const sensibleSegments: CalcResult['sensibleSegments'] = [];
      const excessSegments:   CalcResult['excessSegments']   = [];
      let missingSensibleData = false;

      const t1n = parseFloat(t1Str); const t2n = parseFloat(t2Str);
      // No |T1−T2| threshold — sensible corrections must run even when T_in = T_out
      // because the reactant-cooling and product-heating legs are both non-zero whenever
      // T_op ≠ T_ref (298.15 K), regardless of whether the inlet and outlet are equal.
      const wantSensible = enableTemp && !isNaN(t1n) && !isNaN(t2n);

      if (wantSensible) {
        T1_K = toKelvin(t1n, tempUnit);
        T2_K = toKelvin(t2n, tempUnit);
        if (T1_K <= 0 || T2_K <= 0) { setCalcError('Temperatures must be above 0 K.'); setLoading(false); return; }

        const T_ref = 298.15;
        const allSpecies = [
          ...validReactants.map(r => ({ ...r, side: 'reactant' as const })),
          ...validProducts.map(r => ({  ...r, side: 'product'  as const })),
        ];
        // Cache Shomate + PT data — use persistent cache, only fetch if new species.
        const shomateCache = new Map<string, { shParams: ShomateRow[]; pt: PhaseTransitionRow | null }>();

        await Promise.all(allSpecies.map(async (sp) => {
          if (!sp.selected || !sp.selected.has_shomate) return;
          const coeff = parseFloat(sp.coefficient);
          if (isNaN(coeff) || coeff <= 0) return;

          let entry = shomatePersistCache.current.get(sp.selected.nist_id);
          if (!entry) {
            const nid = sp.selected.nist_id;
            const [shRes, ptRes, jbRes] = await Promise.all([
              supabase.from('nist_shomate').select('phase,eqno,t_min_k,t_max_k,a,b,c,d,e,f,g').eq('nist_id', nid),
              supabase.from('nist_phase_transitions').select('tboil_k,tfus_k,hvap_kjmol,hfus_kjmol').eq('nist_id', nid).maybeSingle(),
              jobackPersistCache.current.has(nid)
                ? Promise.resolve({ data: null, error: null })
                : supabase.from('nist_joback_estimates').select('nist_id,hf_gas_298_kjmol,dgf_gas_298_kjmol,hfus_kjmol,hvap_kjmol,hvap_reference,cp_a,cp_b,cp_c,cp_d,cp_t_min_k,cp_t_max_k,joback_group_count,solid_correction_incomplete,joback_success').eq('nist_id', nid).maybeSingle(),
            ]);
            entry = {
              shParams: (shRes.data ?? []) as ShomateRow[],
              pt: (ptRes.data ?? null) as PhaseTransitionRow | null,
            };
            shomatePersistCache.current.set(nid, entry);
            if (!jobackPersistCache.current.has(nid)) {
              jobackPersistCache.current.set(nid, (jbRes.data ?? null) as JobackEstimate | null);
            }
          }
          shomateCache.set(sp.selected.nist_id, entry);

          const [sT1, sT2] = sp.side === 'reactant'
            ? [T1_K!, T_ref]   // actual direction: inlet → reference
            : [T_ref, T2_K!];  // actual direction: reference → outlet

          const joback = jobackPersistCache.current.get(sp.selected.nist_id) ?? null;
          const { segments, total_kjmol, missingData, approximated } = computeSpeciesSensibleHeat(entry.shParams, entry.pt, sp.selected, sT1, sT2, joback);
          if (missingData) missingSensibleData = true;

          sensibleSegments.push({
            nist_id: sp.selected.nist_id, compound: sp.selected,
            speciesName: sp.selected.name, formula: sp.selected.formula,
            coefficient: coeff, side: sp.side, segments, total_kjmol, missingData, approximated,
          });
        }));

        // ── Correct Hess ΔHf is already done above (always applied) ─────────────
        // The effectiveHf correction was applied unconditionally before this block.
        // Re-read dH_rxn_kjmol from current contributions (already corrected).

        // ── Physical + conversion-based excess ─────────────────────────────
        // Find limiting basis (min qty/stoich across reactants with quantities)
        const cf = conversion / 100;
        let limitingBasis = 1;
        const hasEnteredQtys = validReactants.some(r => r.quantity.trim() !== '' && r.selected);
        if (hasEnteredQtys) {
          let minBasis = Infinity;
          for (const r of validReactants) {
            if (!r.selected || r.quantity.trim() === '') continue;
            const qty = parseFloat(r.quantity);
            const coeff2 = parseFloat(r.coefficient);
            if (isNaN(qty) || qty <= 0 || isNaN(coeff2) || coeff2 <= 0) continue;
            let molsPerStoich: number;
            if (heatUnit === 'kJ/kg') {
              const mw = r.selected.mol_weight;
              if (!mw) continue;
              molsPerStoich = (qty * 1000 / mw) / coeff2;
            } else {
              molsPerStoich = qty / coeff2;
            }
            if (molsPerStoich < minBasis) minBasis = molsPerStoich;
          }
          if (isFinite(minBasis)) limitingBasis = minBasis;
        }

        // Helper: convert entered qty to moles for a reactant row
        const toActualMols = (sp: typeof validReactants[0]): number | null => {
          if (!sp.selected || sp.quantity.trim() === '') return null;
          const qty = parseFloat(sp.quantity);
          if (isNaN(qty) || qty <= 0) return null;
          if (heatUnit === 'kJ/kg') {
            const mw = sp.selected.mol_weight;
            return mw ? (qty * 1000 / mw) : null;
          }
          return qty;
        };

        // Fix reactant sensible coefficients: use all-entered moles normalised to per-stoich-unit
        // so that coeff * total_kjmol is in kJ per stoich-unit (same convention as products)
        if (hasEnteredQtys) {
          for (const seg of sensibleSegments) {
            if (seg.side !== 'reactant') continue;
            const r = validReactants.find(row => row.selected?.formula === seg.formula);
            if (!r) continue;
            const actualMols = toActualMols(r);
            if (actualMols === null) continue;
            // Store per-stoich-unit so that * df = * conv * limitingBasis gives actual kJ
            seg.coefficient = actualMols / limitingBasis;
          }
        }

        // Compute excess for each reactant
        const hasPhysicalExcess = hasEnteredQtys && validReactants.some(r => {
          const mols = toActualMols(r);
          if (mols === null) return false;
          const coeff2 = parseFloat(r.coefficient);
          if (isNaN(coeff2) || coeff2 <= 0) return false;
          return mols > coeff2 * limitingBasis + 1e-9;
        });

        if (cf < 1 || hasPhysicalExcess) {
          await Promise.all(validReactants.map(async (sp) => {
            if (!sp.selected || !sp.selected.has_shomate) return;
            const coeff = parseFloat(sp.coefficient);
            if (isNaN(coeff) || coeff <= 0) return;

            let excessCoeff: number;
            if (hasEnteredQtys) {
              const actualMols = toActualMols(sp);
              if (actualMols === null) {
                // No qty for this reactant — use conversion-based with limiting basis
                excessCoeff = coeff * limitingBasis * (1 - cf) / limitingBasis;
              } else {
                const reactedMols = coeff * limitingBasis * cf;
                const totalUnreacted = Math.max(0, actualMols - reactedMols);
                excessCoeff = totalUnreacted / limitingBasis;
              }
            } else {
              excessCoeff = coeff * (1 - cf);
            }

            if (excessCoeff <= 0) return;

            // Use persistent cache — avoid redundant DB fetch for reactant species
            let entry2 = shomatePersistCache.current.get(sp.selected.nist_id);
            if (!entry2) {
              const nid = sp.selected.nist_id;
              const [shRes, ptRes, jbRes] = await Promise.all([
                supabase.from('nist_shomate').select('phase,eqno,t_min_k,t_max_k,a,b,c,d,e,f,g').eq('nist_id', nid),
                supabase.from('nist_phase_transitions').select('tboil_k,tfus_k,hvap_kjmol,hfus_kjmol').eq('nist_id', nid).maybeSingle(),
                jobackPersistCache.current.has(nid)
                  ? Promise.resolve({ data: null, error: null })
                  : supabase.from('nist_joback_estimates').select('nist_id,hf_gas_298_kjmol,dgf_gas_298_kjmol,hfus_kjmol,hvap_kjmol,hvap_reference,cp_a,cp_b,cp_c,cp_d,cp_t_min_k,cp_t_max_k,joback_group_count,solid_correction_incomplete,joback_success').eq('nist_id', nid).maybeSingle(),
              ]);
              entry2 = {
                shParams: (shRes.data ?? []) as ShomateRow[],
                pt: (ptRes.data ?? null) as PhaseTransitionRow | null,
              };
              shomatePersistCache.current.set(nid, entry2);
              if (!jobackPersistCache.current.has(nid)) {
                jobackPersistCache.current.set(nid, (jbRes.data ?? null) as JobackEstimate | null);
              }
            }
            const joback2 = jobackPersistCache.current.get(sp.selected.nist_id) ?? null;
            const { segments: exSegs, total_kjmol: exTotal, missingData: exMissing, approximated: exApprox } =
              computeSpeciesSensibleHeat(entry2.shParams, entry2.pt, sp.selected, T_ref, T2_K!, joback2);
            if (exMissing) missingSensibleData = true;
            excessSegments.push({
              nist_id: sp.selected.nist_id, compound: sp.selected,
              speciesName: sp.selected.name, formula: sp.selected.formula,
              excessCoefficient: excessCoeff, stoichCoefficient: coeff,
              segments: exSegs, total_kjmol: exTotal, missingData: exMissing, approximated: exApprox,
            });
          }));
        }
        // ── Adiabatic flame temperature (bisection) ──────────────────────────────────────────────────────
        // Find T_out_ad: step1(T_in→298.15) + ΔHrxn + step3(298.15→T_out_ad) = 0
        {
          const _s1 = sensibleSegments.filter(s => s.side === 'reactant')
            .reduce((acc, sp) => acc + sp.coefficient * sp.total_kjmol, 0);
          const _target = -(_s1 + dH_rxn_kjmol);
          const _computeStep3 = (T_out: number): number => {
            let sum = 0;
            for (const sp of validProducts) {
              if (!sp.selected) continue;
              const coeff2 = parseFloat(sp.coefficient);
              if (isNaN(coeff2) || coeff2 <= 0) continue;
              const entry = shomatePersistCache.current.get(sp.selected.nist_id);
              if (!entry) continue;
              const jb2 = jobackPersistCache.current.get(sp.selected.nist_id) ?? null;
              const { total_kjmol: tot } = computeSpeciesSensibleHeat(entry.shParams, entry.pt, sp.selected, T_ref, T_out, jb2);
              sum += coeff2 * tot;
            }
            return sum;
          };
          const T_AD_MAX = 6000;
          const step3AtMax = _computeStep3(T_AD_MAX);
          if (_target >= 0 && _target <= step3AtMax) {
            let _lo = T_ref, _hi = T_AD_MAX;
            for (let _i = 0; _i < 60; _i++) {
              const _mid = (_lo + _hi) / 2;
              if (_computeStep3(_mid) < _target) _lo = _mid; else _hi = _mid;
            }
            adiabaticFlameT_K = (_lo + _hi) / 2;
          }
        }
      }

      let dH_sensible_correction = 0;
      for (const s of sensibleSegments) {
        dH_sensible_correction += s.coefficient * s.total_kjmol;
      }
      // Excess (unreacted) reactants go T_ref → T2 — heat added on product side
      for (const s of excessSegments) {
        dH_sensible_correction += s.excessCoefficient * s.total_kjmol;
      }
      const dH_total_kjmol = dH_rxn_kjmol + dH_sensible_correction;

      const fmtEq = (rows: SpeciesRow[]) =>
        rows.filter(r => r.selected).map(r => {
          const c = parseFloat(r.coefficient);
          return `${c !== 1 ? `${c} ` : ''}${r.selected!.formula}`;
        }).join(' + ');

      const equation = `${fmtEq(validReactants)} \u2192 ${fmtEq(validProducts)}`;
      // Correct total: reactant cooling (step 1) scales with ALL processed (no cf),
      // reaction + product heating scale with cf fraction, excess already has (1-cf) in coefficient.
      const _cf = wantSensible ? (conversion / 100) : 1;
      const _s1 = wantSensible ? sensibleSegments.filter(s => s.side === 'reactant').reduce((acc, sp) => acc + sp.coefficient * sp.total_kjmol, 0) : 0;
      const _s3 = wantSensible ? sensibleSegments.filter(s => s.side === 'product').reduce((acc, sp) => acc + sp.coefficient * sp.total_kjmol, 0) : 0;
      const _ex = wantSensible ? excessSegments.reduce((acc, sp) => acc + sp.excessCoefficient * sp.total_kjmol, 0) : 0;
      const displayDH = wantSensible ? (_cf * dH_rxn_kjmol + _s1 + _cf * _s3 + _ex) : dH_rxn_kjmol;

      // Compute limiting reactant (only when quantities provided)
      const calcLimiting = (): string | null => {
        // Limiting reactant only applies in absolute-quantity mode.
        // In per-mol/kg mode there is no physical quantity to compare.
        if (perMoleMode) return null;
        const withQty = validReactants.filter(r => r.selected && r.quantity.trim() !== '');
        if (withQty.length < 2) return null;
        const bases = withQty.flatMap(r => {
          const q = parseFloat(r.quantity);
          const coeff = parseFloat(r.coefficient);
          if (isNaN(q) || q <= 0 || isNaN(coeff) || coeff <= 0) return [];
          return [{ formula: r.selected!.formula, basis: q / coeff }];
        });
        if (bases.length < 2) return null;
        const min = Math.min(...bases.map(b => b.basis));
        const atMin = bases.filter(b => Math.abs(b.basis - min) < 1e-12);
        if (atMin.length === bases.length) return null; // all stoichiometric
        return atMin[0]?.formula ?? null;
      };

      // Detect phase changes for ephemeral notice (suppress when triggered by conversion slider)
      if (wantSensible && !sliderCalcRef.current) {
        const phaseChanges: {formula: string; atK: number; transition: 'vap' | 'fus'}[] = [];
        for (const sp of [...sensibleSegments, ...excessSegments]) {
          for (const seg of sp.segments) {
            if (seg.type === 'latent') {
              phaseChanges.push({ formula: sp.formula, atK: seg.T1_K, transition: seg.label.includes('Vapor') ? 'vap' : 'fus' });
            }
          }
        }
        setPhaseChangeNotice(phaseChanges.length > 0 ? phaseChanges : null);
      } else if (!wantSensible) {
        setPhaseChangeNotice(null);
      }
      sliderCalcRef.current = false;

      setResult({
        dH_rxn_kjmol, sensibleSegments, excessSegments, dH_total_kjmol, contributions, equation,
        isExothermic: displayDH < 0, T1_K, T2_K,
        hasSensibleHeat: wantSensible && sensibleSegments.length > 0,
        missingHessData, missingSensibleData,
        limitingFormula: calcLimiting(),
        adiabaticFlameT_K, dGf_rxn_kjmol,
      });
    } catch {
      setCalcError('An unexpected error occurred. Please try again.');
    }
    setLoading(false);
  }, [reactants, products, unit, heatUnit, enableTemp, t1Str, t2Str, tempUnit, conversion]);

  // Helper: trigger calculate if enough species are selected
  const tryCalc = useCallback(() => {
    if (reactants.some(r => r.selected) && products.some(p => p.selected)) calculate();
  }, [reactants, products, calculate]);

  const applyAutoBalance = useCallback(() => {
    const rSpecies = reactants.filter(r => r.selected).map(r => ({
      formula: r.selected!.formula,
      coefficient: Math.max(0.001, parseFloat(r.coefficient) || 1),
    }));
    const pSpecies = products.filter(p => p.selected).map(p => ({
      formula: p.selected!.formula,
      coefficient: Math.max(0.001, parseFloat(p.coefficient) || 1),
    }));
    const solution = autoBalance(rSpecies, pSpecies);
    if (!solution) return;
    const rSelected = reactants.filter(r => r.selected);
    const pSelected = products.filter(p => p.selected);
    setReactants(prev => prev.map(r => {
      const idx = rSelected.findIndex(s => s.id === r.id);
      if (idx === -1 || solution.reactantCoeffs[idx] === undefined) return r;
      return { ...r, coefficient: String(solution.reactantCoeffs[idx]) };
    }));
    setProducts(prev => prev.map(p => {
      const idx = pSelected.findIndex(s => s.id === p.id);
      if (idx === -1 || solution.productCoeffs[idx] === undefined) return p;
      return { ...p, coefficient: String(solution.productCoeffs[idx]) };
    }));
    // After auto-balance, the equation is guaranteed to be balanced — update the
    // confirmed snapshot immediately so the warning banner clears at once.
    setConfirmedBalanceCheck({ status: 'balanced', missingElements: [], suggestedCoeffs: null });
    // Signal that a recalculation is needed.
    setPendingAutoCalc(true);
  }, [reactants, products]);

  // Trigger recalculation after auto-balance has updated state.
  // Using useEffect ensures tryCalc() is called with the freshly re-rendered
  // reactants/products, avoiding the stale-closure pitfall of setTimeout.
  useEffect(() => {
    if (!pendingAutoCalc) return;
    setPendingAutoCalc(false);
    tryCalc();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pendingAutoCalc, tryCalc]);

  // Trigger recalculation after the adiabatic solve fills in T_out.
  useEffect(() => {
    if (!pendingAdiabaticCalc) return;
    setPendingAdiabaticCalc(false);
    tryCalc();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [pendingAdiabaticCalc, tryCalc]);

  const handleAdiabaticSolve = useCallback(() => {
    if (!result?.adiabaticFlameT_K) return;
    const T_ad = result.adiabaticFlameT_K;
    setT2Str(tempUnit === 'C' ? (T_ad - 273.15).toFixed(1) : T_ad.toFixed(0));
    setPendingAdiabaticCalc(true);
  }, [result, tempUnit]);

  // Trigger recalculation after the adiabatic solve fills in T_out.

  // Enter key in any input on the left side triggers calculate.
  // If the user typed something but didn't click a suggestion, auto-select the top suggestion.
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      // Auto-select first suggestion for any rows that have suggestions but no selection
      const autoSelect = (rows: SpeciesRow[], set: React.Dispatch<React.SetStateAction<SpeciesRow[]>>) => {
        set(prev => prev.map(r => (!r.selected && r.suggestions.length > 0)
          ? { ...r, selected: r.suggestions[0], query: formatCompoundName(r.suggestions[0].name), suggestions: [], showSuggestions: false }
          : r));
      };
      autoSelect(reactants, setReactants);
      autoSelect(products, setProducts);
      // Mark that this is a user-initiated run so calculate() will snapshot balance.
      confirmBalanceRef.current = true;
      // tryCalc will re-read state after the above setStates; defer one tick
      setTimeout(() => tryCalc(), 0);
    }
  }, [tryCalc, reactants, products]);

  useEffect(() => {
    if (autoCalcDone.current) return;
    if (reactants.some(r => r.selected) && products.some(p => p.selected)) {
      autoCalcDone.current = true;
      calculate();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [reactants, products]);

  // When the conversion slider moves below 100% and we don't yet have excess-segment
  // Shomate data (result was computed at CF=100%), trigger a recalculation once so
  // that section 1b appears immediately without the user pressing Calculate.
  useEffect(() => {
    if (!result || conversion >= 100) return;
    if (result.excessSegments.length === 0) {
      sliderCalcRef.current = true;
      tryCalc();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [conversion, result]);

  // When enableTemp changes, re-run calc so the right side updates immediately.
  const enableTempMounted = useRef(false);
  useEffect(() => {
    if (!enableTempMounted.current) { enableTempMounted.current = true; return; }
    if (autoCalcDone.current) tryCalc();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [enableTemp]);

  // ─── ECharts — combined chart below ────────────────────────────────────────

  // Absolute-mode quantity scale factor: use the LIMITING reactant (minimum qty/coeff ratio)
  const absQtyScale = useMemo(() => {
    if (perMoleMode) return 1;
    let minScale = Infinity;
    for (const r of reactants) {
      const qty_str = committedQtys.get(r.id) ?? '';
      if (!r.selected || qty_str.trim() === '') continue;
      const qty = parseFloat(qty_str);
      const coeff = parseFloat(r.coefficient);
      if (isNaN(qty) || qty <= 0 || isNaN(coeff) || coeff <= 0) continue;
      let scale: number;
      if (heatUnit === 'kJ/kg') {
        const mw = r.selected.mol_weight;
        if (!mw) continue;
        scale = (qty * 1000 / mw) / coeff;
      } else {
        scale = qty / coeff;
      }
      if (scale < minScale) minScale = scale;
    }
    return minScale === Infinity ? 1 : minScale;
  }, [perMoleMode, reactants, heatUnit, committedQtys]);

  // ─── Reactive excess coefficients (update live when conversion slider moves) ─

  const reactiveExcessCoeffs = useMemo((): Map<string, number> => {
    if (!result || !result.hasSensibleHeat || result.excessSegments.length === 0) return new Map();
    const cf = conversion / 100;
    const map = new Map<string, number>();
    for (const ex of result.excessSegments) {
      let newCoeff: number;
      const rRow = !perMoleMode ? reactants.find(r => r.selected?.formula === ex.formula) : undefined;
      const committedQty = rRow ? (committedQtys.get(rRow.id) ?? '') : '';
      const hasQty = !!rRow && committedQty.trim() !== '';
      if (hasQty) {
        const qty    = parseFloat(committedQty);
        const stoich = parseFloat(rRow!.coefficient);
        if (!isNaN(qty) && qty > 0 && !isNaN(stoich) && stoich > 0) {
          const actualMols = heatUnit === 'kJ/kg'
            ? (rRow!.selected?.mol_weight ? qty * 1000 / rRow!.selected.mol_weight : qty)
            : qty;
          const reactedMols = stoich * absQtyScale * cf;
          newCoeff = Math.max(0, actualMols - reactedMols) / absQtyScale;
        } else {
          newCoeff = Math.max(0, ex.stoichCoefficient * (1 - cf));
        }
      } else {
        // Use stored stoichCoefficient so we're not dependent on sensibleSegments lookup.
        // (sensibleSegments may have a different coefficient when hasEnteredQtys modified it.)
        newCoeff = Math.max(0, ex.stoichCoefficient * (1 - cf));
      }
      map.set(ex.formula, Math.max(0, newCoeff)); // always set — including 0, so fallback ?? sp.excessCoefficient never uses stale value
    }
    return map;
  }, [result, conversion, reactants, perMoleMode, heatUnit, absQtyScale, committedQtys]);

  // Fixed y-axis bounds: computed at CF=1 so the axis doesn't move as conversion changes.
  const yAxisBounds = useMemo(() => {
    if (!result || result.contributions.length === 0) return null;
    const map = perMoleMode ? (heatUnit === 'kJ/kg' ? UNIT_CONVERSIONS_KG : UNIT_CONVERSIONS) : UNIT_CONVERSIONS_ABS;
    const conv = map[unit] ?? { factor: 1, label: 'kJ' };
    const pkf = (perMoleMode && heatUnit === 'kJ/kg')
      ? (() => { const r = reactants.find(row => row.selected?.mol_weight); return r?.selected?.mol_weight ? 1000 / r.selected.mol_weight : 1; })()
      : 1;
    const db = conv.factor * pkf * absQtyScale; // dfBase at CF=1
    const vals: number[] = [0, result.dH_rxn_kjmol * db]; // use reaction enthalpy at CF=1 for stable axis bounds
    for (const c of result.contributions) vals.push(c.contribution_kjmol * db);
    if (enableTemp && result.hasSensibleHeat) {      for (const s of result.sensibleSegments) {
        vals.push(s.coefficient * s.total_kjmol * db);
        // at CF=0 all of this reactant is excess — use sensible segment as upper bound
        if (s.side === 'reactant') {
          const ex = result.excessSegments.find(e => e.formula === s.formula);
          if (ex) vals.push(s.coefficient * ex.total_kjmol * db);
        }
      }
      for (const ex of result.excessSegments) vals.push(ex.excessCoefficient * ex.total_kjmol * db);
    }
    const hi = Math.max(...vals, 0); const lo = Math.min(...vals, 0);
    const pad = (hi - lo || 1) * 0.18;
    return { min: lo - pad, max: hi + pad };
  }, [result, perMoleMode, heatUnit, unit, absQtyScale, enableTemp, reactants]); // tempMode needed to recompute bounds when mode changes

  // ─── Combined Enthalpy Chart ──────────────────────────────────────────────
  // One chart: grouped bars — "ΔH°f contribution" (Hess) and "Sensible/latent" per species.
  // Reactants group → ΔH°rxn bar → Products group → Excess group (if any)
  const combinedChartOption = useMemo((): EChartsOption => {
    if (!result || result.contributions.length === 0) return {};
    const activeConvMapC = perMoleMode
      ? (heatUnit === 'kJ/kg' ? UNIT_CONVERSIONS_KG : UNIT_CONVERSIONS)
      : UNIT_CONVERSIONS_ABS;
    const convC = activeConvMapC[unit] ?? { factor: 1, label: perMoleMode ? (heatUnit === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol') : 'kJ' };
    const pkf = (perMoleMode && heatUnit === 'kJ/kg')
      ? (() => { const r = reactants.find(row => row.selected?.mol_weight); return r?.selected?.mol_weight ? 1000 / r.selected.mol_weight : 1; })()
      : 1;
    const cf = conversion / 100;
    const dfBase = convC.factor * pkf * absQtyScale;
    const dfv   = dfBase * cf;

    const showSens = result.hasSensibleHeat && enableTemp;

    const lc = isDark ? '#cbd5e1' : '#475569';
    const sc = isDark ? 'rgba(255,255,255,0.07)' : 'rgba(0,0,0,0.07)';
    const ac = isDark ? '#ffffff' : '#000000';
    const ff = "'Merriweather Sans', sans-serif";
    const yFmt = (v: number) => { const r2 = Math.round(v); const a = Math.abs(r2); return a >= 1e6 ? `${(r2/1e6).toFixed(1)}M` : a >= 1e3 ? `${(r2/1e3).toFixed(1)}k` : `${r2}`; };

    type Entry = { formula: string; hessVal: number|null; sensRVal: number|null; latRVal: number|null; exSensRVal: number|null; sensPVal: number|null; latPVal: number|null; side: 'reactant'|'product' };
    const entries: Entry[] = [];
    const rContribs = result.contributions.filter(c => c.side === 'reactant');
    const pContribs = result.contributions.filter(c => c.side === 'product');

    rContribs.forEach(c => {
      const sensRow = result.sensibleSegments.find(s => s.formula === c.formula && s.side === 'reactant');
      const exRow   = result.excessSegments.find(s => s.formula === c.formula);
      let sensR: number|null = null, latR: number|null = null, exSensR: number|null = null;

      if (showSens && sensRow) {
        const sensOnly = sensRow.segments.filter(sg => sg.type !== 'latent').reduce((s,sg) => s+sg.dH_kjmol, 0);
        const latOnly  = sensRow.segments.filter(sg => sg.type === 'latent').reduce((s,sg) => s+sg.dH_kjmol, 0);

        if (!perMoleMode) {
          const rRow = reactants.find(r => r.selected?.formula === c.formula);
          const qty_i_str = rRow ? (committedQtys.get(rRow.id) ?? '') : '';
          const qty_i = qty_i_str.trim() !== '' ? parseFloat(qty_i_str) : NaN;
          if (!isNaN(qty_i) && qty_i > 0) {
            if (sensOnly !== 0) sensR = qty_i * sensOnly * convC.factor * pkf;
            if (latOnly  !== 0) latR  = qty_i * latOnly  * convC.factor * pkf;
            const reactedMoles = sensRow.coefficient * absQtyScale * cf;
            const excessMoles  = Math.max(0, qty_i - reactedMoles);
            if (excessMoles > 0 && exRow) exSensR = excessMoles * exRow.total_kjmol * convC.factor * pkf;
          } else {
            if (sensOnly !== 0) sensR = sensRow.coefficient * sensOnly * dfv;
            if (latOnly  !== 0) latR  = sensRow.coefficient * latOnly  * dfv;
            if (exRow) {
              const exCoeff = reactiveExcessCoeffs.get(c.formula) ?? exRow.excessCoefficient;
              if (exCoeff > 0) exSensR = exCoeff * exRow.total_kjmol * dfBase;
            }
          }
        } else {
          if (sensOnly !== 0) sensR = sensRow.coefficient * sensOnly * dfv;
          if (latOnly  !== 0) latR  = sensRow.coefficient * latOnly  * dfv;
          if (exRow) {
            const exCoeff = reactiveExcessCoeffs.get(c.formula) ?? exRow.excessCoefficient;
            if (exCoeff > 0) exSensR = exCoeff * exRow.total_kjmol * dfBase;
          }
        }
      }
      entries.push({ formula: c.formula, hessVal: c.contribution_kjmol*dfv,
        sensRVal: sensR, latRVal: latR, exSensRVal: exSensR,
        sensPVal: null, latPVal: null, side: 'reactant' });
    });

    pContribs.forEach(c => {
      const sensRow = result.sensibleSegments.find(s => s.formula === c.formula && s.side === 'product');
      let sensP: number|null = null, latP: number|null = null;
      if (showSens && sensRow) {
        const sensOnly = sensRow.segments.filter(sg => sg.type !== 'latent').reduce((s,sg) => s+sg.dH_kjmol, 0);
        const latOnly  = sensRow.segments.filter(sg => sg.type === 'latent').reduce((s,sg) => s+sg.dH_kjmol, 0);
        if (sensOnly !== 0) sensP = sensRow.coefficient * sensOnly * dfv;
        if (latOnly  !== 0) latP  = sensRow.coefficient * latOnly  * dfv;
      }
      entries.push({ formula: c.formula, hessVal: c.contribution_kjmol*dfv,
        sensRVal: null, latRVal: null, exSensRVal: null,
        sensPVal: sensP, latPVal: latP, side: 'product' });
    });

    const hessRSeries  = entries.map(e => e.side === 'reactant' ? e.hessVal  : null);
    const hessPSeries  = entries.map(e => e.side === 'product'  ? e.hessVal  : null);
    const sensRRaw     = entries.map(e => e.sensRVal);
    const latRSeries   = entries.map(e => e.latRVal);
    // Merge excess reactant sensible heat into the regular reactant sensible bar
    const sensRSeries  = entries.map((e, i) => {
      const s = sensRRaw[i];
      const ex = e.exSensRVal;
      if (s === null && ex === null) return null;
      return (s ?? 0) + (ex ?? 0);
    });
    const sensPSeries  = entries.map(e => e.sensPVal);
    const latPSeries   = entries.map(e => e.latPVal);

    const hasSensR   = sensRSeries.some(v => v !== null);
    const hasLatR    = latRSeries.some(v => v !== null);
    const hasSensP   = sensPSeries.some(v => v !== null);
    const hasLatP    = latPSeries.some(v => v !== null);
    const hasThermal = hasSensR || hasLatR || hasSensP || hasLatP;

    // Colors defined at module level as CHR_* constants.

    const mutedGroupColor = isDark ? '#94a3b8' : '#64748b';
    const xLabels = entries.map(e => echartsFormulaRich(e.formula));

    const markAreas: any[] = [];
    if (rContribs.length > 0) markAreas.push([
      { xAxis: 0,
        itemStyle: { color: 'rgba(167,139,250,0.06)' },
        label: { show: true, position: 'insideTop', formatter: '\u2014 Reactants \u2014',
          color: mutedGroupColor, fontSize: 8, fontFamily: ff, align: 'center' } },
      { xAxis: rContribs.length - 1 },
    ]);
    if (pContribs.length > 0) markAreas.push([
      { xAxis: rContribs.length,
        itemStyle: { color: 'rgba(45,212,191,0.06)' },
        label: { show: true, position: 'insideTop', formatter: '\u2014 Products \u2014',
          color: mutedGroupColor, fontSize: 8, fontFamily: ff, align: 'center' } },
      { xAxis: rContribs.length + pContribs.length - 1 },
    ]);

    const allVals = entries.flatMap(e => [e.hessVal, e.sensRVal, e.latRVal, e.exSensRVal, e.sensPVal, e.latPVal])
      .filter((v): v is number => v !== null);
    const maxAbs = Math.max(...allVals.map(Math.abs), 1);

    // Legend items — rich-text names for subscripted display
    // Series names use plain keys (for internal ECharts matching); legend.formatter
    // converts them to rich-text with subscripted f / sens / lat.
    const legendItems: { name: string; itemStyle: { color: string } }[] = [
      { name: 'dHf-R',    itemStyle: { color: CHR_FORM_R } },
      { name: 'dHf-P',    itemStyle: { color: CHR_FORM_P } },
      ...(hasSensR ? [{ name: 'dHsens-R', itemStyle: { color: CHR_SENS_R } }] : []),
      ...(hasLatR  ? [{ name: 'dHlat-R',  itemStyle: { color: CHR_LAT_R  } }] : []),
      ...(hasSensP ? [{ name: 'dHsens-P', itemStyle: { color: CHR_SENS_P } }] : []),
      ...(hasLatP  ? [{ name: 'dHlat-P',  itemStyle: { color: CHR_LAT_P  } }] : []),
    ];
    void legendItems; // rendered by custom JSX legend outside ECharts
    const legendFormatter = (name: string): string => {
      const side = name.endsWith('-R') ? ' (Reactants)' : ' (Products)';
      if (name.startsWith('dHf'))    return `ΔH°{sub|f}${side}`;
      if (name.startsWith('dHsens')) return `ΔH{sub|sens}${side}`;
      if (name.startsWith('dHlat'))  return `ΔH{sub|lat}${side}`;
      return name;
    };

    // Convert a chemical formula to HTML with subscript digits (for tooltip)
    const fmtFormulaHtml = (raw: string) =>
      normalizeFormula(raw).replace(/(\d+)/g, n => `<sub>${n}</sub>`);
    const tooltipSeriesNameHtml = (name: string): string => {
      if (name.startsWith('dHf'))    return 'ΔH°<sub>f</sub>';
      if (name.startsWith('dHsens')) return 'ΔH<sub>sens</sub>';
      if (name.startsWith('dHlat'))  return 'ΔH<sub>lat</sub>';
      return name;
    };

    const tooltipFn = (p: any) => {
      if (!entries[p.dataIndex]) return '';
      const val = p.value;
      if (val === null || val === undefined) return '';
      return `<b>${fmtFormulaHtml(entries[p.dataIndex].formula)}</b><br/>${tooltipSeriesNameHtml(p.seriesName)}: <b>${val >= 0 ? '+' : ''}${fmtNum(val)} ${convC.label}</b>`;
    };

    // All bars use tight grouped stacks: ΔHf, ΔHsens, ΔHlat each share one stack
    // across R and P. Since at any x-position only one side (R or P) is non-null,
    // stacking R+P of the same type gives exactly ONE bar per species per type.
    // This eliminates the empty-slot gap that appeared when they were separate series.
    // All bars use tight grouped stacks: ΔHf, ΔHsens, ΔHlat each share one stack
    // across R and P. Since at any x-position only one side (R or P) is non-null,
    // stacking R+P gives exactly ONE bar per species per type — no empty-slot gaps.
    const mkBar = (data: (number|null)[], color: string, stack: string) =>
      ({ type: 'bar' as const, barMaxWidth: 28, barGap: '0%', color, stack,
         data: data.map(v => v === null ? null : { value: v, itemStyle: { color, borderRadius: v >= 0 ? [4,4,0,0] : [0,0,4,4] } }) as any });

    void hasThermal; // no longer needed — legend is JSX

    return {
      animation: false,
      backgroundColor: 'transparent',
      textStyle: { fontFamily: ff },
      legend: undefined, // legend rendered as custom JSX outside the chart
      grid: { left: 64, right: 20, top: 28, bottom: 24, containLabel: true },
      xAxis: {
        type: 'category',
        data: xLabels,
        axisLabel: {
          color: lc, interval: 0, fontSize: 10, fontFamily: ff, lineHeight: 14,
          rich: { sub: { fontSize: 7.5, fontFamily: ff, verticalAlign: 'bottom', lineHeight: 14 } },
          formatter: (v: string) => v,
        },
        axisLine: { lineStyle: { color: ac } },
        axisTick: { lineStyle: { color: ac } },
      },
      yAxis: {
        type: 'value',
        name: `\u0394H (${convC.label})`,
        nameLocation: 'middle', nameRotate: 90,
        nameGap: Math.max(48, yFmt(-(yAxisBounds?.min ?? maxAbs)).length * 7 + 24),
        nameTextStyle: { color: lc, fontSize: 10, fontFamily: ff },
        axisLabel: { color: (value: number) => value > 0 ? '#38bdf8' : value < 0 ? '#fb923c' : lc, fontSize: 10, fontFamily: ff, formatter: (value: number) => value > 0 ? `+${yFmt(value)}` : yFmt(value), showMinLabel: false, showMaxLabel: false },
        axisLine: { show: true, lineStyle: { color: ac } },
        splitLine: { lineStyle: { color: sc } },
        min: yAxisBounds?.min,
        max: yAxisBounds?.max,
      },
      series: [
        {
          name: 'dHf-R',
          ...mkBar(hessRSeries, CHR_FORM_R, 'form'),
          markArea: { silent: true, data: markAreas as any },
          markLine: {
            symbol: ['none','none'], silent: true,
            // Correct total for markLine: step 1 scales with dfBase; reaction + step 3 scale with dfv; excess with dfBase
            ...((() => {
              const mS1 = result.hasSensibleHeat ? result.sensibleSegments.filter(s => s.side === 'reactant').reduce((a, sp) => a + sp.coefficient * sp.total_kjmol, 0) : 0;
              const mS3 = result.hasSensibleHeat ? result.sensibleSegments.filter(s => s.side === 'product').reduce((a, sp) => a + sp.coefficient * sp.total_kjmol, 0) : 0;
              const mEx = result.hasSensibleHeat ? result.excessSegments.reduce((a, sp) => a + (reactiveExcessCoeffs.get(sp.formula) ?? sp.excessCoefficient) * sp.total_kjmol, 0) : 0;
              const mTot = result.hasSensibleHeat ? mS1 * dfBase + result.dH_rxn_kjmol * dfv + mS3 * dfv + mEx * dfBase : result.dH_rxn_kjmol * dfv;
              const mColor = mTot < 0 ? '#fb923c' : mTot > 0 ? '#38bdf8' : lc;
              const mLineColor = mTot < 0 ? 'rgba(251,146,60,0.55)' : mTot > 0 ? 'rgba(56,189,248,0.55)' : (isDark ? 'rgba(255,255,255,0.35)' : 'rgba(0,0,0,0.35)');
              return {
                lineStyle: { type: 'dashed' as any, color: mLineColor, width: 1.5 },
                label: { show: true, formatter: `\u0394H\u00b0tot = ${fmtNum(mTot)} ${convC.label}`, position: 'insideStartTop' as any, color: mColor, fontSize: 9, fontFamily: ff },
                data: [{ yAxis: mTot }],
              };
            })()),
          },
        } as any,
        { name: 'dHf-P',    ...mkBar(hessPSeries, CHR_FORM_P, 'form') },
        ...(hasSensR ? [{ name: 'dHsens-R', ...mkBar(sensRSeries, CHR_SENS_R, 'sens') }] : []),
        ...(hasLatR  ? [{ name: 'dHlat-R',  ...mkBar(latRSeries,  CHR_LAT_R,  'lat')  }] : []),
        ...(hasSensP ? [{ name: 'dHsens-P', ...mkBar(sensPSeries, CHR_SENS_P, 'sens') }] : []),
        ...(hasLatP  ? [{ name: 'dHlat-P',  ...mkBar(latPSeries,  CHR_LAT_P,  'lat')  }] : []),
      ],
      tooltip: {
        trigger: 'item',
        backgroundColor: isDark ? '#1e293b' : '#fff',
        borderColor: isDark ? '#334155' : '#e2e8f0',
        textStyle: { color: isDark ? '#e2e8f0' : '#1e293b', fontSize: 12, fontFamily: ff },
        formatter: tooltipFn,
      },
    } as EChartsOption;
  }, [result, isDark, unit, heatUnit, reactants, conversion, enableTemp, perMoleMode, absQtyScale, reactiveExcessCoeffs, committedQtys, yAxisBounds]);

  // Chart legend meta — which thermal series are present (for custom JSX legend)
  const chartLegendMeta = useMemo(() => {
    if (!result || !result.hasSensibleHeat || !enableTemp) return { hasSens: false, hasLat: false };
    const segs = result.sensibleSegments;
    return {
      hasSens: segs.some(s => s.segments.some(sg => sg.type !== 'latent')),
      hasLat:  segs.some(s => s.segments.some(sg => sg.type === 'latent')),
    };
  }, [result, enableTemp]);
  const renderSide = (side: 'reactants' | 'products') => {
    const rows  = side === 'reactants' ? reactants : products;
    const label = side === 'reactants' ? 'Reactants' : 'Products';
    const showQty = side === 'reactants' && !perMoleMode;
    return (
      <div className="space-y-2">
        <Label className="text-xs font-semibold text-muted-foreground uppercase tracking-wide">{label}</Label>
        <div className="flex gap-2 px-0.5">
          <span className="w-16 text-center text-xs text-muted-foreground shrink-0">Coefficient</span>
          <span className="flex-1 text-xs text-muted-foreground">Species</span>
          {showQty && (
            <span className="w-20 text-center text-xs text-muted-foreground shrink-0">{heatUnit === 'kJ/kg' ? 'kg' : 'mol'}</span>
          )}
          <span className="w-9 shrink-0" />
        </div>
        {rows.map(row => {
          // Warn when a selected compound has null ΔHf but is NOT mono-elemental
          // (null on a multi-element compound means missing data, not "= 0 by convention").
          const warnMissingHf = !!row.selected && row.selected.hf_gas_298_kjmol === null &&
            !isMonoElemental(row.selected.formula);
          // Warn when the compound has no Shomate Cp data at all.
          const warnNoCp = !!row.selected && !row.selected.has_shomate;
          const showWarn = warnMissingHf || warnNoCp;
          return (
          <div key={row.id} className="relative flex gap-2 items-center" ref={el => { if (el) wrapperRefs.current.set(row.id, el); }}>
            {showWarn && (
              <Tooltip>
                <TooltipTrigger asChild>
                  <div className="absolute -left-5 top-1/2 -translate-y-1/2 text-yellow-500 cursor-help select-none">
                    <svg width="12" height="12" viewBox="0 0 12 12" aria-hidden="true" style={{ display: 'block', pointerEvents: 'none', userSelect: 'none' }}>
                      <path d="M6 1L11.2 10H0.8L6 1Z" fill="currentColor" />
                      <circle cx="6" cy="8" r="0.9" fill="black" />
                      <rect x="5.35" y="4.2" width="1.3" height="2.6" rx="0.65" fill="black" />
                    </svg>
                  </div>
                </TooltipTrigger>
                <TooltipContent side="left" className="max-w-52 text-xs">
                  {warnMissingHf && <p>&Delta;H&deg;<sub>f</sub> is missing for this compound and is treated as 0, which may be incorrect.</p>}
                  {warnNoCp && <p>No Cp / Shomate data &mdash; sensible heat cannot be computed for this species.</p>}
                </TooltipContent>
              </Tooltip>
            )}
            <Input className="w-16 text-center h-9 shrink-0" placeholder="1" value={row.coefficient}
              onChange={e => patchRow(side, row.id, { coefficient: e.target.value })}
              onKeyDown={handleKeyDown} />
            <div className="relative flex-1">
              <Input className="h-9" placeholder="Name or formula" value={row.query}
                onChange={e => handleQueryChange(side, row.id, e.target.value)}
                onKeyDown={handleKeyDown}
                onFocus={() => { if (row.suggestions.length > 0) patchRow(side, row.id, { showSuggestions: true }); }} />
              {row.showSuggestions && row.suggestions.length > 0 && (
                <div className="absolute z-50 top-full left-0 right-0 mt-1 rounded-md border bg-popover shadow-md max-h-64 overflow-y-auto">
                  {row.suggestions.map(s => (
                    <div key={s.nist_id} className="px-3 py-2 hover:bg-accent cursor-pointer"
                      onMouseDown={e => { e.preventDefault(); handleSelect(side, row.id, s); }}>
                      <div className="flex items-center justify-between gap-3">
                        <p className="text-sm font-medium leading-tight truncate">{formatCompoundName(s.name)}</p>
                        <span className="text-xs font-mono text-muted-foreground shrink-0">{renderFormula(s.formula)}</span>
                      </div>
                    </div>
                  ))}
                </div>
              )}
            </div>
            {showQty && (
              <Input className="w-20 text-center h-9 shrink-0"
                type="number" min="0" step="any"
                placeholder={heatUnit === 'kJ/kg' ? 'kg' : 'mol'}
                value={row.quantity}
                onChange={e => patchRow(side, row.id, { quantity: e.target.value })}
                onKeyDown={handleKeyDown} />
            )}
            {rows[0].id !== row.id ? (
              <Button size="icon" variant="ghost"
                className="h-9 w-9 shrink-0 text-red-400 hover:text-red-500 hover:bg-red-400/10"
                onClick={() => removeRow(side, row.id)}>
                <Trash2 className="h-4 w-4" />
              </Button>
            ) : <div className="h-9 w-9 shrink-0" />}
          </div>
          );
        })}
        <Button variant="outline" size="sm" className="w-full h-8 mt-1" onClick={() => addRow(side)}>
          <PlusCircle className="h-3.5 w-3.5 mr-1.5" />Add {label.slice(0, -1)}
        </Button>
      </div>
    );
  };

  const activeConvMap = perMoleMode
    ? (heatUnit === 'kJ/kg' ? UNIT_CONVERSIONS_KG : UNIT_CONVERSIONS)
    : UNIT_CONVERSIONS_ABS;
  const conv = activeConvMap[unit] ?? { factor: 1, label: perMoleMode ? (heatUnit === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol') : 'kJ' };

  const conversionFactor = conversion / 100;
  const perKgFactor2 = (perMoleMode && heatUnit === 'kJ/kg') ? perKgFactor : 1;
  const df  = conv.factor * perKgFactor2 * absQtyScale;
  const dfc = conv.factor * perKgFactor2 * conversionFactor * absQtyScale;
  // Correct display total: step 1 (reactant cooling) uses df — ALL reactants pass through regardless of conversion;
  // step 2 (reaction) and step 3 (product heating) use dfc — only cf fraction reacts/forms;
  // excess coefficient already encodes (1-cf) so it uses df.
  const _dStep1 = result?.hasSensibleHeat ? result.sensibleSegments.filter(s => s.side === 'reactant').reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0) : 0;
  const _dStep3 = result?.hasSensibleHeat ? result.sensibleSegments.filter(s => s.side === 'product').reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0) : 0;
  const _dExcess = result?.hasSensibleHeat ? result.excessSegments.reduce((s, sp) => s + (reactiveExcessCoeffs.get(sp.formula) ?? sp.excessCoefficient) * sp.total_kjmol, 0) : 0;
  const displayDH = result
    ? (result.hasSensibleHeat
        ? _dStep1 * df + result.dH_rxn_kjmol * dfc + _dStep3 * dfc + _dExcess * df
        : result.dH_rxn_kjmol * dfc)
    : 0;
  // Derived live from displayDH so the card updates as the conversion slider moves
  const liveIsExothermic = result ? displayDH < 0 : false;
  const fmtT = (T_K: number) =>
    tempUnit === 'C' ? `${(T_K - 273.15).toFixed(1)} \u00b0C` : `${T_K.toFixed(0)} K`;

  const runSensitivity = (varName: 'T_in' | 'T_out' | 'conversion') => {
    if (!result) return;
    setSensRunning(true);
    setSensVar(varName);
    const T_ref = 298.15;
    const pts: number[][] = [];
    const computeStep3 = (T_out: number): number => {
      let sum = 0;
      for (const sp of result.sensibleSegments.filter(s => s.side === 'product')) {
        const entry = shomatePersistCache.current.get(sp.nist_id);
        if (!entry) continue;
        const jb = jobackPersistCache.current.get(sp.nist_id) ?? null;
        const { total_kjmol } = computeSpeciesSensibleHeat(entry.shParams, entry.pt, sp.compound, T_ref, T_out, jb);
        sum += sp.coefficient * total_kjmol;
      }
      return sum;
    };
    const computeStep1 = (T_in: number): number => {
      let sum = 0;
      for (const sp of result.sensibleSegments.filter(s => s.side === 'reactant')) {
        const entry = shomatePersistCache.current.get(sp.nist_id);
        if (!entry) continue;
        const jb = jobackPersistCache.current.get(sp.nist_id) ?? null;
        const { total_kjmol } = computeSpeciesSensibleHeat(entry.shParams, entry.pt, sp.compound, T_in, T_ref, jb);
        sum += sp.coefficient * total_kjmol;
      }
      return sum;
    };
    const step1Fixed = result.sensibleSegments.filter(s => s.side === 'reactant')
      .reduce((acc, sp) => acc + sp.coefficient * sp.total_kjmol, 0);
    const step3Fixed = result.T2_K ? computeStep3(result.T2_K) : 0;
    const dH_rxn = result.dH_rxn_kjmol;

    // Compute max reliable Cp temperature from Shomate / Joback caches
    const productIds = products.filter(p => p.selected).map(p => p.selected!.nist_id);
    const reactantIds = reactants.filter(r => r.selected).map(r => r.selected!.nist_id);
    const maxCpProduct  = productIds.length  > 0 ? Math.min(...productIds.map(id  => getMaxCpT(id)).filter(v => v > 0)) : 3500;
    const maxCpReactant = reactantIds.length > 0 ? Math.min(...reactantIds.map(id => getMaxCpT(id)).filter(v => v > 0)) : 2000;

    // Data-point step: 5 K (or 1% for conversion) — gives smooth tooltips at round values
    const dataStep = 5; // K

    if (varName === 'conversion') {
      const _exStoich = result.excessSegments.reduce((acc, sp) => acc + sp.stoichCoefficient * sp.total_kjmol, 0);
      for (let pct = 0; pct <= 100; pct++) {
        const cf = pct / 100;
        pts.push([pct, (step1Fixed + cf * dH_rxn + cf * step3Fixed + (1 - cf) * _exStoich) * conv.factor]);
      }
    } else if (varName === 'T_out' && result.T1_K != null) {
      const maxT = Math.min(maxCpProduct, 6000);
      // Sweep from well below current T_out down to ~200 K so cold settings are always visible.
      // Use min(200, currentT-50) so we start 50 K below the current value, floored at 200 K max / 50 K min.
      const sweepLow_K = Math.max(50, Math.min(200, result.T2_K != null ? result.T2_K - 50 : 200));
      const startDisplay = tempUnit === 'C'
        ? Math.floor((sweepLow_K - 273.15) / dataStep) * dataStep
        : Math.floor(sweepLow_K / dataStep) * dataStep;
      const endDisplay = tempUnit === 'C' ? maxT - 273.15 : maxT;
      for (let d = startDisplay; d <= endDisplay; d += dataStep) {
        const T_K = tempUnit === 'C' ? d + 273.15 : d;
        pts.push([d, (step1Fixed + dH_rxn + computeStep3(T_K)) * conv.factor]);
      }
    } else if (varName === 'T_in' && result.T2_K != null) {
      const maxT = Math.min(maxCpReactant, 6000);
      const sweepLow_K = Math.max(50, Math.min(200, result.T1_K != null ? result.T1_K - 50 : 200));
      const startDisplay = tempUnit === 'C'
        ? Math.floor((sweepLow_K - 273.15) / dataStep) * dataStep
        : Math.floor(sweepLow_K / dataStep) * dataStep;
      const endDisplay = tempUnit === 'C' ? maxT - 273.15 : maxT;
      for (let d = startDisplay; d <= endDisplay; d += dataStep) {
        const T_K = tempUnit === 'C' ? d + 273.15 : d;
        pts.push([d, (computeStep1(T_K) + dH_rxn + step3Fixed) * conv.factor]);
      }
    }

    const font = "'Merriweather Sans', sans-serif";
    const textColor = isDark ? '#e2e8f0' : '#1a202c';
    const tooltipBg = isDark ? '#1e293b' : '#ffffff';
    const tooltipBorder = isDark ? '#475569' : '#e2e8f0';
    const axisLineColor = isDark ? '#e2e8f0' : '#1a202c';
    const splitLineColor = isDark ? 'rgba(255,255,255,0.08)' : 'rgba(0,0,0,0.08)';

    // ECharts rich-text axis names with proper subscripts
    const richSub = { name: { fontSize: 13, fontFamily: font }, sub: { fontSize: 9, fontFamily: font, verticalAlign: 'bottom' as const } };
    const xAxisName = varName === 'conversion' ? 'Conversion (%)'
      : varName === 'T_out' ? `{name|T}{sub|out}{name| (${tempUnit === 'C' ? '\u00b0C' : 'K'})}`
      : `{name|T}{sub|in}{name| (${tempUnit === 'C' ? '\u00b0C' : 'K'})}`;
    const yAxisName = `{name|\u0394H\u00b0}{sub|tot}{name| (${conv.label})}`;

    // Compute nice x-axis tick interval based on data range
    const xVals = pts.map(p => p[0]);
    const xRange = xVals.length > 1 ? xVals[xVals.length - 1] - xVals[0] : 100;
    const xMinInterval = varName === 'conversion'
      ? (xRange <= 20 ? 2 : xRange <= 50 ? 5 : 10)
      : (xRange <= 50 ? 5 : xRange <= 100 ? 10 : xRange <= 250 ? 25 : xRange <= 500 ? 50 : xRange <= 1000 ? 100 : 250);

    // Dynamic y-axis nameGap: must exceed tick-label pixel width so the name doesn't overlap.
    // With containLabel:true and grid.left = yNameGap+15, axis-line sits at
    // max(yNameGap+15, tickLabelWidth) from canvas left, so the rotated axis name
    // always lands ≥ 15 px inside the canvas edge.
    const yMax = Math.max(...pts.map(p => Math.abs(p[1])), 1);
    const yMin = Math.min(...pts.map(p => p[1]));
    const yDigits = Math.floor(Math.log10(yMax)) + 1;
    const hasNeg = yMin < 0;
    // Estimate tick-label width: each digit/char ≈ 7.5 px in 12 px font.
    const tickLabelWidth = Math.max(35, (yDigits + (hasNeg ? 2 : 1)) * 7.5 + 12);
    const yNameGap = Math.min(90, Math.max(55, tickLabelWidth + 12));
    const gridLeft = yNameGap + 15;

    // Annotate zero crossing (adiabatic operating point)
    let zeroCross: number | null = null;
    for (let i = 1; i < pts.length; i++) {
      if ((pts[i-1][1] >= 0 && pts[i][1] < 0) || (pts[i-1][1] <= 0 && pts[i][1] > 0)) {
        zeroCross = pts[i-1][0] + (0 - pts[i-1][1]) / (pts[i][1] - pts[i-1][1]) * (pts[i][0] - pts[i-1][0]);
        break;
      }
    }

    // Phase transition annotations for T_out and T_in sweeps
    // Collect tboil / tfus for all product (T_out) or reactant (T_in) species that fall in range.
    const phaseMarkLines: any[] = [];
    if (varName === 'T_out' || varName === 'T_in') {
      const sweepSpecies = varName === 'T_out'
        ? products.filter(p => p.selected).map(p => p.selected!)
        : reactants.filter(r => r.selected).map(r => r.selected!);
      const xLo = pts.length > 0 ? pts[0][0] : -Infinity;
      const xHi = pts.length > 0 ? pts[pts.length - 1][0] : Infinity;
      const seenT = new Set<number>();
      for (const sp of sweepSpecies) {
        const ptRow = shomatePersistCache.current.get(sp.nist_id)?.pt;
        const tboilK = ptRow?.tboil_k ?? sp.tboil_k ?? null;
        const tfusK  = ptRow?.tfus_k  ?? sp.tfus_k  ?? null;
        const richFormula = echartsFormulaRich(sp.formula);
        for (const [T_K, prefix, color] of [
          [tboilK, 'BP', isDark ? '#fb923c' : '#ea580c'] as [number | null, string, string],
          [tfusK,  'MP', isDark ? '#c084fc' : '#9333ea'] as [number | null, string, string],
        ]) {
          if (T_K == null) continue;
          const xVal = tempUnit === 'C' ? T_K - 273.15 : T_K;
          if (xVal < xLo || xVal > xHi) continue;
          if (seenT.has(xVal)) continue;
          seenT.add(xVal);
          phaseMarkLines.push({
            xAxis: xVal,
            lineStyle: { color, type: 'dashed', width: 1 },
            label: {
              formatter: `${prefix}(${richFormula})`,
              color,
              fontFamily: font,
              fontSize: 9,
              position: 'insideEndTop',
              align: 'left',
              offset: [4, 4],
              rich: { sub: { fontSize: 7, fontFamily: font, verticalAlign: 'bottom' } },
            },
          });
        }
      }
    }
    setSensOpts({
      backgroundColor: 'transparent',
      animation: false,
      tooltip: {
        trigger: 'axis',
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        extraCssText: `box-shadow:0 4px 12px rgba(0,0,0,${isDark ? '0.4' : '0.1'});`,
        textStyle: { color: textColor, fontFamily: font, fontSize: 11 },
        formatter: (p: any) => {
          // Read theme at render time so tooltip text color tracks theme switches
          const tc = typeof document !== 'undefined' && document.documentElement.classList.contains('dark')
            ? '#e2e8f0' : '#1a202c';
          const xVal = +p[0].data[0];
          const x = xVal.toFixed(varName === 'conversion' ? 0 : 1);
          const y = fmtNum(p[0].data[1]);
          const xLabel = varName === 'conversion' ? 'Conversion'
            : varName === 'T_out' ? 'T<sub style="font-size:0.8em">out</sub>'
            : 'T<sub style="font-size:0.8em">in</sub>';
          const xUnit = varName === 'conversion' ? '%' : (tempUnit === 'C' ? '\u00b0C' : 'K');
          return `<span style="font-family:'Merriweather Sans',sans-serif;font-size:11px;color:${tc}">${xLabel}: <b>${x} ${xUnit}</b><br/>\u0394H\u00b0<sub style="font-size:0.75em">tot</sub>: <b>${y}\u2009${conv.label}</b></span>`;
        },
      },
      xAxis: {
        type: 'value',
        name: xAxisName,
        nameLocation: 'middle',
        nameGap: 35,
        nameTextStyle: { color: textColor, fontFamily: font, rich: richSub },
        axisLabel: { color: textColor, fontFamily: font },
        axisLine: { lineStyle: { color: axisLineColor } },
        splitLine: { show: false },
        minInterval: xMinInterval,
      },
      yAxis: {
        type: 'value',
        name: yAxisName,
        nameLocation: 'middle',
        nameGap: yNameGap,
        nameTextStyle: { color: textColor, fontFamily: font, rich: richSub, padding: [0, 0, 0, 0] },
        axisLabel: { color: textColor, fontFamily: font, fontSize: 11 },
        axisLine: { lineStyle: { color: axisLineColor } },
        splitLine: { lineStyle: { color: splitLineColor } },
      },
      series: [{
        type: 'line',
        data: pts,
        smooth: false,
        lineStyle: { color: '#60a5fa', width: 2 },
        itemStyle: { color: '#60a5fa' },
        showSymbol: false,
        animation: false,
        markLine: (zeroCross != null || phaseMarkLines.length > 0) ? {
          silent: true, symbol: 'none',
          animation: false,
          data: [
            ...(zeroCross != null ? [{
              xAxis: zeroCross,
              lineStyle: { color: isDark ? '#4b5563' : '#9ca3af', type: 'dashed', width: 1 },
              label: {
                formatter: varName === 'conversion' ? 'Adiabatic\nConversion' : 'Adiabatic',
                color: isDark ? '#9ca3af' : '#6b7280',
                fontFamily: font,
                fontSize: 10,
                position: 'middle' as any,
                align: 'center',
                offset: [0, 0],
              },
            }] : []),
            ...phaseMarkLines,
          ],
        } : undefined,
      }],
      grid: { left: gridLeft, right: 30, top: 52, bottom: 60, containLabel: true },
      title: {
        text: varName === 'conversion'
          ? '\u0394H\u00b0 Sensitivity \u2014 Conversion'
          : `\u0394H\u00b0 Sensitivity \u2014 {nm|T}{sb|${varName === 'T_out' ? 'out' : 'in'}}{nm| (${tempUnit === 'C' ? '\u00b0C' : 'K'})}`,
        textStyle: {
          color: textColor, fontFamily: font, fontSize: 12, fontWeight: 'bold' as const,
          rich: {
            nm: { fontSize: 12, fontFamily: font, fontWeight: 'bold' as const, color: textColor },
            sb: { fontSize: 9, fontFamily: font, verticalAlign: 'bottom' as const, fontWeight: 'bold' as const, color: textColor },
          },
        },
        left: 'center',
        top: 6,
      },
    } as EChartsOption);
    setSensRunning(false);
  };

  const downloadPathwayCSV = () => {
    if (!result) return;
    const escape = (s: string) => `"${s.replace(/"/g, '""')}"`;
    const segStage = (seg: HeatSegment) =>
      seg.type === 'latent'
        ? `${seg.label.includes('Vapor') ? '\u0394Hvap' : '\u0394Hfus'} at ${seg.T1_K.toFixed(1)} K`
        : `${seg.T1_K.toFixed(1)}\u2192${seg.T2_K.toFixed(1)} K`;
    const csvRows: string[] = [
      ['Section', 'Formula', 'Stage', 'Type', 'T1 (K)', 'T2 (K)', 'dH (kJ/mol)'].join(','),
    ];
    for (const sp of result.sensibleSegments.filter(s => s.side === 'reactant')) {
      for (const seg of sp.segments) {
        csvRows.push([
          escape('Reactant Sensible Heat'), escape(normalizeFormula(sp.formula)),
          escape(segStage(seg)), seg.type,
          seg.T1_K.toFixed(2), seg.T2_K.toFixed(2),
          (seg.dH_kjmol * sp.coefficient).toFixed(4),
        ].join(','));
      }
    }
    csvRows.push([
      escape('Reaction at 298.15 K'), escape(result.equation),
      escape("Hess's Law"), 'reaction',
      '298.15', '298.15', result.dH_rxn_kjmol.toFixed(4),
    ].join(','));
    for (const c of result.contributions) {
      csvRows.push([
        escape('Reaction at 298.15 K'), escape(normalizeFormula(c.formula)),
        escape(`\u0394H\u00b0f (${c.side})`), 'formation',
        '298.15', '298.15',
        c.contribution_kjmol.toFixed(4),
      ].join(','));
    }
    for (const sp of result.sensibleSegments.filter(s => s.side === 'product')) {
      for (const seg of sp.segments) {
        csvRows.push([
          escape('Product Sensible Heat'), escape(normalizeFormula(sp.formula)),
          escape(segStage(seg)), seg.type,
          seg.T1_K.toFixed(2), seg.T2_K.toFixed(2),
          (seg.dH_kjmol * sp.coefficient).toFixed(4),
        ].join(','));
      }
    }
    for (const sp of result.excessSegments.filter(s => reactiveExcessCoeffs.has(s.formula))) {
      const exCoeff = reactiveExcessCoeffs.get(sp.formula) ?? sp.excessCoefficient;
      for (const seg of sp.segments) {
        csvRows.push([
          escape('Excess Reactant Heat'), escape(normalizeFormula(sp.formula)),
          escape(segStage(seg)), seg.type,
          seg.T1_K.toFixed(2), seg.T2_K.toFixed(2),
          (seg.dH_kjmol * exCoeff).toFixed(4),
        ].join(','));
      }
    }
    const blob = new Blob([csvRows.join('\n')], { type: 'text/csv;charset=utf-8;' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'enthalpy-pathway.csv'; a.click();
    URL.revokeObjectURL(url);
  };

  // ─── JSX ──────────────────────────────────────────────────────────────────

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">

          {/* Left: Controls */}
          <div className="lg:col-span-1 space-y-4">

            <Card>
              <CardHeader className="pb-2">
                <div className="flex items-center justify-between">
                  <CardTitle className="text-base">Build Reaction</CardTitle>
                  <div className="flex items-center gap-2">
                    <span className="text-xs text-muted-foreground">Per mol/kg</span>
                    <button
                      className={`relative inline-flex h-5 w-9 items-center rounded-full transition-colors shrink-0 ${perMoleMode ? 'bg-primary' : 'bg-muted'}`}
                      onClick={() => {
                        const next = !perMoleMode;
                        setPerMoleMode(next);
                        setUnit(next ? (heatUnit === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol') : 'kJ');
                      }}
                      role="switch" aria-checked={perMoleMode}>
                      <span className={`inline-block h-3.5 w-3.5 transform rounded-full bg-white transition-transform ${perMoleMode ? 'translate-x-4' : 'translate-x-0.5'}`} />
                    </button>
                  </div>
                </div>
              </CardHeader>
              <CardContent className="space-y-3">
                <div className="grid grid-cols-2 gap-3">
                  <div className="space-y-1.5">
                    <Label className="text-xs text-muted-foreground">&Delta;H units</Label>
                    <Select value={unit} onValueChange={setUnit}>
                      <SelectTrigger className="h-9 w-full"><SelectValue /></SelectTrigger>
                      <SelectContent>
                        {!perMoleMode ? (
                          <>
                            <SelectItem value="kJ">kJ</SelectItem>
                            <SelectItem value="kcal">kcal</SelectItem>
                            <SelectItem value="J">J</SelectItem>
                            <SelectItem value="BTU">BTU</SelectItem>
                          </>
                        ) : heatUnit === 'kJ/kg' ? (
                          <>
                            <SelectItem value="kJ/kg">kJ/kg</SelectItem>
                            <SelectItem value="J/kg">J/kg</SelectItem>
                            <SelectItem value="kcal/kg">kcal/kg</SelectItem>
                            <SelectItem value="BTU/lb">BTU/lb</SelectItem>
                          </>
                        ) : (
                          <>
                            <SelectItem value="kJ/mol">kJ/mol</SelectItem>
                            <SelectItem value="kcal/mol">kcal/mol</SelectItem>
                            <SelectItem value="J/mol">J/mol</SelectItem>
                            <SelectItem value="BTU/lbmol">BTU/lbmol</SelectItem>
                          </>
                        )}
                      </SelectContent>
                    </Select>
                  </div>
                  {perMoleMode && (
                  <div className="space-y-1.5">
                    <Label className="text-xs text-muted-foreground">Per</Label>
                    <Tabs value={heatUnit} onValueChange={v => {
                      const nu = v as 'kJ/mol' | 'kJ/kg';
                      setHeatUnit(nu);
                      setUnit(nu === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol');
                    }}>
                      <TabsList className="grid w-full grid-cols-2">
                        <TabsTrigger value="kJ/mol">mol</TabsTrigger>
                        <TabsTrigger value="kJ/kg">kg</TabsTrigger>
                      </TabsList>
                    </Tabs>
                  </div>
                  )}
                  {!perMoleMode && (
                  <div className="space-y-1.5">
                    <Label className="text-xs text-muted-foreground">Qty unit</Label>
                    <Tabs value={heatUnit} onValueChange={v => {
                      const nu = v as 'kJ/mol' | 'kJ/kg';
                      setHeatUnit(nu);
                    }}>
                      <TabsList className="grid w-full grid-cols-2">
                        <TabsTrigger value="kJ/mol">mol</TabsTrigger>
                        <TabsTrigger value="kJ/kg">kg</TabsTrigger>
                      </TabsList>
                    </Tabs>
                  </div>
                  )}
                </div>
                <div className="space-y-1.5">
                  <div className="flex justify-between items-center">
                    <Label className="text-xs text-muted-foreground">Conversion</Label>
                    <span className="text-xs font-mono tabular-nums text-muted-foreground">{conversion}%</span>
                  </div>
                  <input
                    type="range" min={0} max={100} step={1}
                    value={conversion}
                    onChange={e => setConversion(Number(e.target.value))}
                    className="w-full h-1.5 accent-primary cursor-pointer rounded-full"
                  />
                </div>
                <div className="border-t" />
                {renderSide('reactants')}
                <div className="border-t" />
                {renderSide('products')}
              </CardContent>
            </Card>

            <Card>
              <CardHeader className="pb-2">
                <div className="flex items-center justify-between gap-2">
                  <div className="flex items-center gap-2 min-w-0">
                    <CardTitle className="text-base shrink-0">Temperature Change</CardTitle>
                    <Tabs value={tempUnit} onValueChange={handleTempUnitToggle}>
                      <TabsList className="h-6 p-0.5">
                        <TabsTrigger value="C" className="h-5 px-2 text-xs">&deg;C</TabsTrigger>
                        <TabsTrigger value="K" className="h-5 px-2 text-xs">K</TabsTrigger>
                      </TabsList>
                    </Tabs>
                  </div>
                  <button
                    className={`relative inline-flex h-5 w-9 items-center rounded-full transition-colors shrink-0 ${enableTemp ? 'bg-primary' : 'bg-muted'}`}
                    onClick={() => setEnableTemp(v => !v)}
                    role="switch" aria-checked={enableTemp}>
                    <span className={`inline-block h-3.5 w-3.5 transform rounded-full bg-white transition-transform ${enableTemp ? 'translate-x-4' : 'translate-x-0.5'}`} />
                  </button>
                </div>
              </CardHeader>
              <CardContent className={`space-y-3 transition-opacity ${(!mounted || enableTemp) ? 'opacity-100' : 'opacity-40 pointer-events-none'}`}>
                <div className="grid grid-cols-2 gap-2">
                  <div className="flex items-center gap-2">
                    <Label className="text-xs text-muted-foreground shrink-0"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">in</sub>:</span></Label>
                    <Input
                      className={`h-9 flex-1 ${enableTemp && toKelvin(parseFloat(t1Str), tempUnit) <= 0 && t1Str !== '' ? 'border-destructive' : ''}`}
                      value={t1Str} onChange={e => setT1Str(e.target.value)} onKeyDown={handleKeyDown}
                      placeholder={tempUnit === 'C' ? '30' : '303'} />
                  </div>
                  <div className="flex items-center gap-2">
                    <Label className="text-xs text-muted-foreground shrink-0"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">out</sub>:</span></Label>
                    <Input
                      className={`h-9 flex-1 ${enableTemp && toKelvin(parseFloat(t2Str), tempUnit) <= 0 && t2Str !== '' ? 'border-destructive' : ''}`}
                      value={t2Str} onChange={e => setT2Str(e.target.value)} onKeyDown={handleKeyDown}
                      placeholder={tempUnit === 'C' ? '100' : '373'} />
                  </div>
                </div>
                {enableTemp && (() => {
                  const t1K = toKelvin(parseFloat(t1Str), tempUnit);
                  const t2K = toKelvin(parseFloat(t2Str), tempUnit);
                  const t1Err = !isNaN(t1K) && t1K <= 0;
                  const t2Err = !isNaN(t2K) && t2K <= 0;
                  return <>
                    {t1Err && <p className="text-xs text-destructive">T_in must be above 0 K ({tempUnit==='C'?'min: \u2212273.15 \u00b0C':'min: 0 K'})</p>}
                    {t2Err && <p className="text-xs text-destructive">T_out must be above 0 K ({tempUnit==='C'?'min: \u2212273.15 \u00b0C':'min: 0 K'})</p>}
                    {!t1Err && !t2Err && (
                      <>
                        <Button variant="outline" className="w-full h-9 text-xs"
                          onClick={handleAdiabaticSolve} disabled={loading || !result?.adiabaticFlameT_K}>
                          <span className="inline-flex items-baseline">Solve Adiabatic T<sub className="text-[0.65em] leading-none">out</sub></span>
                        </Button>
                      </>
                    )}
                  </>;
                })()}
              </CardContent>
            </Card>

            <Button className="w-full h-10" onClick={() => { confirmBalanceRef.current = true; calculate(); }} disabled={loading}>
              {loading ? 'Calculating\u2026' : <><span className="inline-flex items-baseline">Calculate &Delta;H&deg;<sub className="text-[0.65em] leading-none">tot</sub></span></>}
            </Button>

            {/* Sensitivity Analysis */}
            <Button
              variant={sensMode ? 'default' : 'outline'}
              className="w-full h-9 text-xs"
              onClick={() => { if (sensMode) { setSensMode(false); setSensVar(null); setSensOpts(null); } else { setSensMode(true); } }}
              disabled={!result}
            >
              {sensMode ? 'Exit Sensitivity Analysis' : 'Sensitivity Analysis'}
            </Button>
            {sensMode && result && (
              <div className="grid gap-1.5" style={{ gridTemplateColumns: `repeat(${[result.hasSensibleHeat && result.T1_K != null, result.hasSensibleHeat && result.T2_K != null, true].filter(Boolean).length}, 1fr)` }}>
                {result.hasSensibleHeat && result.T1_K != null && (
                  <Button size="sm" variant={sensVar === 'T_in' ? 'default' : 'outline'} className="text-xs h-8 px-2"
                    onClick={() => runSensitivity('T_in')} disabled={sensRunning}>
                    <span className="inline-flex items-baseline gap-[1px]">Sweep T<sub className="text-[0.6em] leading-none">in</sub></span>
                  </Button>
                )}
                {result.hasSensibleHeat && result.T2_K != null && (
                  <Button size="sm" variant={sensVar === 'T_out' ? 'default' : 'outline'} className="text-xs h-8 px-2"
                    onClick={() => runSensitivity('T_out')} disabled={sensRunning}>
                    <span className="inline-flex items-baseline gap-[1px]">Sweep T<sub className="text-[0.6em] leading-none">out</sub></span>
                  </Button>
                )}
                <Button size="sm" variant={sensVar === 'conversion' ? 'default' : 'outline'} className="text-xs h-8 px-2"
                  onClick={() => runSensitivity('conversion')} disabled={sensRunning}>
                  Sweep X
                </Button>
              </div>
            )}

            {/* Balance check banner — only shown after Enter or Calculate press */}
            {confirmedBalanceCheck && confirmedBalanceCheck.status !== 'balanced' && (
              <div className={`flex flex-col gap-2 rounded-md border p-3 text-xs ${
                confirmedBalanceCheck.status === 'fixable'
                  ? 'border-yellow-500/50 bg-yellow-500/10 text-yellow-600 dark:text-yellow-400'
                  : 'border-red-500/50 bg-red-500/10 text-red-600 dark:text-red-400'
              }`}>
                <div className="flex items-start gap-2">
                  <AlertCircle className="h-3.5 w-3.5 mt-0.5 shrink-0" />
                  {confirmedBalanceCheck.status === 'fixable' ? (
                    <span>
                      This reaction doesn&apos;t appear to be balanced.<br />
                      Adjusting the coefficients can fix it.
                    </span>
                  ) : confirmedBalanceCheck.missingElements.length > 0 ? (
                    <span>
                      This reaction cannot be balanced as&nbsp;
                      {confirmedBalanceCheck.missingElements.join(', ')}
                      {confirmedBalanceCheck.missingElements.length === 1 ? ' appears' : ' appear'} on only one side.<br />
                      Please check that all species are correct.
                    </span>
                  ) : (
                    <span>
                      This reaction cannot be balanced &mdash; no valid positive integer coefficients exist for this combination of species.<br />
                      Try adding a missing co-product (e.g. HCl, H₂O) or check that all species are correct.
                    </span>
                  )}
                </div>
                {confirmedBalanceCheck.status === 'fixable' && confirmedBalanceCheck.suggestedCoeffs && (
                  <button
                    onClick={applyAutoBalance}
                    className="self-start rounded border border-yellow-500/60 bg-yellow-500/10 px-3 py-1 text-xs font-medium text-yellow-700 dark:text-yellow-300 hover:bg-yellow-500/20 transition-colors"
                  >
                    Auto-balance coefficients
                  </button>
                )}
              </div>
            )}

            {phaseChangeNotice && phaseChangeNotice.length > 0 && (
              <div className="flex items-start justify-between gap-2 rounded-md border border-blue-500/40 bg-blue-500/10 p-3 text-xs text-blue-400">
                <span>
                  <span className="font-semibold">Phase change detected</span> during heat path:{' '}
                  {phaseChangeNotice.map((pc, i) => (
                    <span key={i}>{i > 0 ? ', ' : ''}
                      <span className="font-mono">{renderFormula(pc.formula)}</span>
                      {' '}({pc.transition === 'vap' ? 'vaporization' : 'fusion'} at {pc.atK.toFixed(1)} K)
                    </span>
                  ))}.
                </span>
                <button onClick={() => setPhaseChangeNotice(null)} className="shrink-0 text-blue-400/60 hover:text-blue-400 text-sm leading-none">&times;</button>
              </div>
            )}

            {calcError && (
              <div className="flex items-start gap-2 rounded-md border border-destructive/50 bg-destructive/10 p-3 text-sm text-destructive">
                <AlertCircle className="h-4 w-4 mt-0.5 shrink-0" /><span>{calcError}</span>
              </div>
            )}
          </div>

          {/* Right: Results */}
          <div className="lg:col-span-2 space-y-4">
            {sensMode && sensOpts ? (
              <Card>
                <CardContent className="pt-4">
                  <div className="relative">
                    <ReactECharts ref={sensEchartsRef} option={sensOpts} style={{ height: 560 }} notMerge />
                    <button
                      onClick={() => {
                        const chart = sensEchartsRef.current?.getEchartsInstance?.();
                        if (!chart) return;
                        const url = chart.getDataURL({ type: 'png', pixelRatio: 2, backgroundColor: isDark ? '#0f172a' : '#ffffff' });
                        const a = document.createElement('a');
                        a.href = url;
                        a.download = `sensitivity-${sensVar ?? 'plot'}.png`;
                        a.click();
                      }}
                      className="absolute bottom-2 right-2 z-10 p-1.5 rounded bg-background/80 hover:bg-background border border-border text-muted-foreground hover:text-foreground transition-colors"
                      title="Download plot as PNG"
                    >
                      <Download className="h-4 w-4" />
                    </button>
                  </div>
                </CardContent>
              </Card>
            ) : result ? (
              <div className={`space-y-4 ${loading ? 'opacity-60 pointer-events-none transition-opacity duration-200' : 'transition-opacity duration-200'}`}>
                <Card className={`border-2 transition-colors ${liveIsExothermic ? 'border-orange-400/50' : 'border-sky-400/50'}`}>
                  <CardContent className="pt-5 pb-5">
                    <div className="flex items-center justify-between gap-4">
                      <div className="flex-1 min-w-0 text-center">
                        <p className="text-2xl font-mono font-semibold leading-snug break-words mb-3">
                          {result.contributions.filter(c => c.side === 'reactant').map((c, i) => (
                            <React.Fragment key={i}>
                              {i > 0 && <span className="mx-1.5 opacity-50">+</span>}
                              {c.coefficient !== 1 && <span>{c.coefficient} </span>}
                              {result.limitingFormula === c.formula ? (
                                <Tooltip>
                                  <TooltipTrigger asChild>
                                    <span className="text-yellow-400 cursor-default">{renderFormula(c.formula)}</span>
                                  </TooltipTrigger>
                                  <TooltipContent side="top">Limiting Reactant</TooltipContent>
                                </Tooltip>
                              ) : <span>{renderFormula(c.formula)}</span>}
                            </React.Fragment>
                          ))}
                          <span className="mx-2.5 opacity-50">&rarr;</span>
                          {result.contributions.filter(c => c.side === 'product').map((c, i) => (
                            <React.Fragment key={i}>
                              {i > 0 && <span className="mx-1.5 opacity-50">+</span>}
                              {c.coefficient !== 1 && <span>{c.coefficient} </span>}
                              {renderFormula(c.formula)}
                            </React.Fragment>
                          ))}
                        </p>
                        <div className="flex items-baseline gap-2 justify-center">
                          <span className="text-xl font-mono">&Delta;H&deg;<sub className="text-[0.65em] leading-none">tot</sub> =</span>
                          <span className="text-2xl font-bold tabular-nums tracking-tight">
                            {displayDH >= 0 ? '+' : ''}{fmtNum(displayDH)}
                          </span>
                          <span className="text-xl">{conv.label}</span>
                        </div>
                      </div>
                      <div className={`flex flex-col items-center justify-center rounded-xl px-6 py-4 shrink-0 ${liveIsExothermic ? 'bg-orange-400/10 text-orange-400' : 'bg-sky-400/10 text-sky-400'}`}>
                        {liveIsExothermic ? <Flame className="h-9 w-9 mb-1" /> : <Snowflake className="h-9 w-9 mb-1" />}
                        <span className="font-semibold text-sm">{liveIsExothermic ? 'Exothermic' : 'Endothermic'}</span>
                        <span className="text-xs opacity-70 mt-0.5">{liveIsExothermic ? 'Releases heat' : 'Absorbs heat'}</span>
                      </div>
                    </div>
                    {result.missingHessData.length > 0 && (
                      <div className="mt-3 flex items-start gap-2 rounded-md border border-yellow-500/50 bg-yellow-500/10 p-2 text-xs text-yellow-600 dark:text-yellow-400">
                        <AlertCircle className="h-3.5 w-3.5 mt-0.5 shrink-0" />
                        <span>&Delta;H&deg;f missing for: {result.missingHessData.join(', ')}. These species are excluded from &Delta;H&deg;rxn.</span>
                      </div>
                    )}
                    {result.missingSensibleData && (
                      <div className="mt-2 flex items-start gap-2 rounded-md border border-yellow-500/50 bg-yellow-500/10 p-2 text-xs text-yellow-600 dark:text-yellow-400">
                        <AlertCircle className="h-3.5 w-3.5 mt-0.5 shrink-0" />
                        <span>Cp(T) or latent heat data missing for one or more species; correction is partial.</span>
                      </div>
                    )}
                  </CardContent>
                </Card>

                {result.hasSensibleHeat ? (
                  /* ── 3-step Kirchhoff story ── */
                  <Card>
                    <CardHeader className="py-0 cursor-pointer select-none" onClick={() => setKirchhoffOpen(v => !v)}>
                      <div className="flex items-center gap-2">
                        <CardTitle className="text-sm">Thermochemical Enthalpy Pathway</CardTitle>
                        {kirchhoffOpen && (
                          <button
                            onClick={e => { e.stopPropagation(); downloadPathwayCSV(); }}
                            className="flex items-center gap-1 text-[10px] text-muted-foreground hover:text-foreground border border-muted-foreground/30 hover:border-foreground/50 rounded px-1.5 py-0.5 transition-colors"
                          >
                            <Download className="w-3 h-3" />CSV
                          </button>
                        )}
                        <span className="ml-auto text-muted-foreground text-xs">{kirchhoffOpen ? '▲' : '▼'}</span>
                      </div>
                    </CardHeader>
                    <CardContent className={`space-y-5 pb-4 ${kirchhoffOpen ? '' : 'hidden'}`}>

                      {/* Shared table renderer */}
                      {(() => {
                        const renderSegTable = (
                          rows: { coeff: number; formula: string; colorClass: string; segments: HeatSegment[]; missingData: boolean; approximated?: boolean }[],
                          emptyMsg: React.ReactNode,
                          scale = dfc
                        ) => {
                          const hasRows = rows.some(r => r.segments.some(seg => seg.type === 'latent' || Math.abs(seg.T2_K - seg.T1_K) > 0.5));
                          if (!hasRows) return <p className="text-xs text-muted-foreground italic pl-6">{emptyMsg}</p>;
                          return (
                            <div className="overflow-x-auto">
                              <table className="w-full text-xs table-fixed">
                                <colgroup>
                                  <col style={{ width: '20%' }} />
                                  <col style={{ width: '35%' }} />
                                  <col style={{ width: '13%' }} />
                                  <col style={{ width: '32%' }} />
                                </colgroup>
                                <thead>
                                  <tr className="border-b text-muted-foreground bg-muted/10">
                                    <th className="px-2 py-1 text-left font-medium">Species</th>
                                    <th className="px-2 py-1 text-left font-medium">Stage</th>
                                    <th className="px-2 py-1 text-center font-medium">Type</th>
                                    <th className="px-2 py-1 text-right font-medium">&Delta;H ({conv.label})</th>
                                  </tr>
                                </thead>
                                <tbody>
                                  {rows.flatMap((r, ri) => {
                                    const visSegs = r.segments.filter(seg => seg.type === 'latent' || Math.abs(seg.T2_K - seg.T1_K) > 0.5);
                                    // Phase label is per-segment: (g)/(l)/(s) for sensible, (l→g)/(s→l) for latent
                                    return visSegs.map((seg, si) => {
                                      const segPhaseLabel = seg.type === 'latent'
                                        ? (seg.label.includes('Vapor') ? '(l\u2192g)' : '(s\u2192l)')
                                        : seg.phase === 'gas' ? '(g)' : seg.phase === 'liquid' ? '(l)' : seg.phase === 'solid' ? '(s)' : null;
                                      return (
                                      <tr key={`${ri}-${si}`} className="border-b last:border-0 hover:bg-muted/20">
                                        <td className={`px-2 py-1.5 font-mono ${r.colorClass} whitespace-nowrap overflow-hidden text-ellipsis`}>
                                            {renderFormula(r.formula)}
                                            {segPhaseLabel && <span className="ml-1 text-muted-foreground font-normal text-[10px]">{segPhaseLabel}</span>}
                                            {(r.missingData || r.approximated) && (() => {
                                              const notes = r.segments.flatMap(s => s.warning ? [s.warning] : []).filter((v, i, a) => a.indexOf(v) === i);
                                              const isApproxOnly = r.approximated && !r.missingData;
                                              const tipText = notes.length ? notes.join(' \u00b7 ') : (isApproxOnly ? 'Joback group-contribution estimate (approximate)' : 'Missing data \u2014 assumed 0');
                                              return (
                                                <Tooltip>
                                                  <TooltipTrigger asChild>
                                                    <span className={`ml-1 cursor-help ${isApproxOnly ? 'text-cyan-400' : 'text-yellow-500'}`}>&#x26A0;</span>
                                                  </TooltipTrigger>
                                                  <TooltipContent side="top" className="max-w-[260px] text-xs whitespace-pre-wrap">{tipText}</TooltipContent>
                                                </Tooltip>
                                              );
                                            })()}
                                          </td>
                                          <td className="px-2 py-1.5 text-muted-foreground whitespace-nowrap overflow-hidden text-ellipsis">
                                            {seg.type === 'latent'
                                              ? <>{seg.label.includes('Vapor') ? <>&Delta;H<sub style={{fontSize:'0.78em',lineHeight:1}}>vap</sub></> : <>&Delta;H<sub style={{fontSize:'0.78em',lineHeight:1}}>fus</sub></>}{' at '}{fmtT(seg.T1_K)}</>
                                              : `${fmtT(seg.T1_K)} \u2192 ${fmtT(seg.T2_K)}`}
                                            {seg.warning && (
                                              <Tooltip>
                                                <TooltipTrigger asChild>
                                                  <span className="ml-1 text-yellow-500 cursor-help">&#x26A0;</span>
                                                </TooltipTrigger>
                                                <TooltipContent side="top" className="max-w-[260px] text-xs whitespace-pre-wrap">{seg.warning}</TooltipContent>
                                              </Tooltip>
                                            )}
                                          </td>
                                          <td className="px-2 py-1.5 text-center">
                                            <span className="text-[10px] text-muted-foreground">{seg.type === 'latent' ? 'Latent' : 'Sensible'}</span>
                                          </td>
                                          <td className={`px-2 py-1.5 text-right tabular-nums font-medium whitespace-nowrap ${seg.dH_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                            {seg.dH_kjmol >= 0 ? '+' : ''}{fmtNum(seg.dH_kjmol * r.coeff * scale)}
                                          </td>
                                      </tr>
                                      );
                                    });
                                  })}
                                </tbody>
                              </table>
                            </div>
                          );
                        };

                        const rSegs = result.sensibleSegments.filter(sp => sp.side === 'reactant');
                        const pSegs = result.sensibleSegments.filter(sp => sp.side === 'product');
                        // Only show excess section when there is actually non-zero excess
                        // (at 100% CF with stoichiometric reactants all coefficients are 0)
                        const showExcess = result.hasSensibleHeat && reactiveExcessCoeffs.size > 0
                          && [...reactiveExcessCoeffs.values()].some(v => v > 1e-9);
                        const hasReactantHeat = rSegs.some(sp => sp.segments.some(seg => seg.type === 'latent' || Math.abs(seg.T2_K - seg.T1_K) > 0.5));

                        return (
                          <>
                            {/* Step 1 — split into 1a/1b when excess exists */}
                            {showExcess ? (
                              <div>
                                <p className="text-xs font-semibold text-violet-400 mb-2 flex items-center gap-1.5">
                                  <span className="inline-flex items-center justify-center w-5 h-4 rounded-full bg-violet-400/20 text-[10px] font-bold shrink-0">1a</span>
                                  Reactants sensible + latent heat ({fmtT(result.T1_K!)} &rarr; {fmtT(298.15)})
                                </p>
                                {renderSegTable(
                                  rSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-violet-400', segments: sp.segments, missingData: sp.missingData, approximated: sp.approximated })),
                                  <span className="text-xs text-muted-foreground italic"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">in</sub></span> = {fmtT(298.15)} &mdash; no sensible heat needed for reactants.</span>
                                )}
                                {hasReactantHeat && (() => {
                                  const total1 = rSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                                  return (
                                    <div className="mt-1 flex justify-end text-xs pr-2">
                                      <span className={`tabular-nums font-bold ${total1 >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                        {total1 >= 0 ? '+' : ''}{fmtNum(total1 * df)} {conv.label}
                                      </span>
                                    </div>
                                  );
                                })()}
                                {/* 1b — excess */}
                                <div className="mt-4">
                                  <p className="text-xs font-semibold text-lime-400 mb-2 flex items-center gap-1.5">
                                    <span className="inline-flex items-center justify-center w-5 h-4 rounded-full bg-lime-400/20 text-[10px] font-bold shrink-0">1b</span>
                                    Excess / unreacted sensible + latent heat &mdash; {fmtT(298.15)} &rarr; {fmtT(result.T2_K!)}
                                  </p>
                                  {result.excessSegments.some(sp => reactiveExcessCoeffs.has(sp.formula)) ? (
                                    <>
                                      {renderSegTable(
                                        result.excessSegments
                                          .filter(sp => reactiveExcessCoeffs.has(sp.formula))
                                          .map(sp => ({ coeff: +(reactiveExcessCoeffs.get(sp.formula)!).toFixed(4), formula: sp.formula, colorClass: 'text-violet-400', segments: sp.segments, missingData: sp.missingData, approximated: sp.approximated })),
                                        <span className="text-xs text-muted-foreground italic">No excess heat needed.</span>,
                                        df
                                      )}
                                      {(() => {
                                        const totalEx = result.excessSegments
                                          .filter(sp => reactiveExcessCoeffs.has(sp.formula))
                                          .reduce((s, sp) => s + (reactiveExcessCoeffs.get(sp.formula) ?? 0) * sp.total_kjmol, 0);
                                        return (
                                          <div className="mt-1 flex justify-end text-xs pr-2">
                                            <span className={`tabular-nums font-bold ${totalEx >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                              {totalEx >= 0 ? '+' : ''}{fmtNum(totalEx * df)} {conv.label}
                                            </span>
                                          </div>
                                        );
                                      })()}
                                    </>
                                  ) : (
                                    <p className="text-xs text-muted-foreground italic px-1">Press &ldquo;Calculate &Delta;H&deg;<sub className="text-[0.65em] leading-none">tot</sub>&rdquo; to compute excess heat contributions.</p>
                                  )}
                                </div>
                              </div>
                            ) : (
                              <div>
                                <p className="text-xs font-semibold text-violet-400 mb-2 flex items-center gap-1.5">
                                  <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-violet-400/20 text-[10px] font-bold shrink-0">1</span>
                                  Reactants sensible + latent heat ({fmtT(result.T1_K!)} &rarr; {fmtT(298.15)})
                                </p>
                                {renderSegTable(
                                  rSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-violet-400', segments: sp.segments, missingData: sp.missingData, approximated: sp.approximated })),
                                  <span className="text-xs text-muted-foreground italic"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">in</sub></span> = {fmtT(298.15)} &mdash; no sensible heat needed for reactants.</span>
                                )}
                                {hasReactantHeat && (() => {
                                  const total1 = rSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                                  return (
                                    <div className="mt-1 flex justify-end text-xs pr-2">
                                      <span className={`tabular-nums font-bold ${total1 >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                        {total1 >= 0 ? '+' : ''}{fmtNum(total1 * df)} {conv.label}
                                      </span>
                                    </div>
                                  );
                                })()}
                              </div>
                            )}

                            {/* Step 2 */}
                            <div>
                              <p className="text-xs font-semibold text-amber-400 mb-2 flex items-center gap-1.5">
                                <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-amber-400/20 text-[10px] font-bold shrink-0">2</span>
                                <span className="inline-flex items-baseline">&Delta;H&deg;<sub className="text-[0.65em] leading-none">rxn</sub></span> at {fmtT(298.15)} &mdash; Hess&apos;s Law
                              </p>
                              <div className="overflow-x-auto">
                                <table className="w-full text-xs table-fixed">
                                  <colgroup>
                                    <col style={{ width: '20%' }} />
                                    <col style={{ width: '35%' }} />
                                    <col style={{ width: '13%' }} />
                                    <col style={{ width: '32%' }} />
                                  </colgroup>
                                  <thead>
                                    <tr className="border-b text-muted-foreground bg-muted/10">
                                      <th className="px-2 py-1 text-left font-medium">Species</th>
                                      <th className="px-2 py-1 text-left font-medium">&Delta;H&deg;<sub>f</sub> ({conv.label})</th>
                                      <th className="px-2 py-1 text-center font-medium">Coefficient</th>
                                      <th className="px-2 py-1 text-right font-medium">&Delta;H ({conv.label})</th>
                                    </tr>
                                  </thead>
                                  <tbody>
                                    {result.contributions.map((c, i) => (
                                      <tr key={i} className="border-b last:border-0 hover:bg-muted/20">
                                        <td className={`px-2 py-1.5 font-mono whitespace-nowrap overflow-hidden text-ellipsis ${c.side === 'reactant' ? 'text-violet-400' : 'text-teal-400'}`}>
                                          {renderFormula(c.formula)}
                                          {c.phase298 && <span className="ml-0.5 text-muted-foreground font-normal text-[9px] font-sans">({c.phase298 === 'gas' ? 'g' : c.phase298 === 'liquid' ? 'l' : 's'})</span>}
                                        </td>
                                        <td className="px-2 py-1.5 tabular-nums text-muted-foreground whitespace-nowrap overflow-hidden text-ellipsis">{fmtNum(c.hof_kjmol * df)}</td>
                                        <td className="px-2 py-1.5 text-center tabular-nums text-muted-foreground">{c.side === 'reactant' ? '\u2212' : '+'}{c.coefficient}</td>
                                        <td className={`px-2 py-1.5 text-right tabular-nums font-medium whitespace-nowrap ${c.contribution_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                          {c.contribution_kjmol >= 0 ? '+' : ''}{fmtNum(c.contribution_kjmol * dfc)}
                                        </td>
                                      </tr>
                                    ))}
                                    <tr className="border-t-2 bg-muted/20">
                                      <td className="px-2 py-1.5 font-bold" colSpan={3}><span className="inline-flex items-baseline">&Delta;H&deg;<sub className="text-[0.65em] leading-none">rxn</sub></span> ({fmtT(298.15)})</td>
                                      <td className={`px-2 py-1.5 text-right tabular-nums font-bold whitespace-nowrap ${result.dH_rxn_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                        {result.dH_rxn_kjmol >= 0 ? '+' : ''}{fmtNum(result.dH_rxn_kjmol * dfc)} {conv.label}
                                      </td>
                                    </tr>
                                  </tbody>
                                </table>
                              </div>
                            </div>

                            {/* Step 3 */}
                            <div>
                              <p className="text-xs font-semibold text-teal-400 mb-2 flex items-center gap-1.5">
                                <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-teal-400/20 text-[10px] font-bold shrink-0">3</span>
                                Products sensible + latent heat ({fmtT(298.15)} &rarr; {fmtT(result.T2_K!)})
                              </p>
                              {renderSegTable(
                                pSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-teal-400', segments: sp.segments, missingData: sp.missingData, approximated: sp.approximated })),
                                <span className="text-xs text-muted-foreground italic"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">out</sub></span> = {fmtT(298.15)} &mdash; no sensible heat needed for products.</span>
                              )}
                              {pSegs.length > 0 && (() => {
                                const total3 = pSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                                return (
                                  <div className="mt-1 flex justify-end text-xs pr-2">
                                    <span className={`tabular-nums font-bold ${total3 >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                      {total3 >= 0 ? '+' : ''}{fmtNum(total3 * dfc)} {conv.label}
                                    </span>
                                  </div>
                                );
                              })()}
                            </div>

                            {/* Total — use correctly weighted sum */}
                            {(() => {
                              const kT1 = rSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                              const kT3 = pSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                              const kEx = result.excessSegments.reduce((s, sp) => s + (reactiveExcessCoeffs.get(sp.formula) ?? 0) * sp.total_kjmol, 0);
                              const kTot = kT1 * df + result.dH_rxn_kjmol * dfc + kT3 * dfc + kEx * df;
                              return (
                                <div className="border-t-2 mt-2 pt-3 flex items-center justify-between pr-2">
                                  <span className="text-xs font-bold"><span className="inline-flex items-baseline">&Delta;H&deg;<sub className="text-[0.65em] leading-none">tot</sub></span> ({fmtT(result.T1_K!)} &rarr; {fmtT(result.T2_K!)})</span>
                                  <span className={`tabular-nums font-bold text-sm whitespace-nowrap ${kTot >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                    {kTot >= 0 ? '+' : ''}{fmtNum(kTot)} {conv.label}
                                  </span>
                                </div>
                              );
                            })()}
                          </>
                        );
                      })()}
                    </CardContent>
                  </Card>
                ) : (
                  /* ── Standalone Hess's Law table ── */
                  <Card>
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm">&Delta;H&deg;<sub>f</sub> Contributions &mdash; Hess&apos;s Law (298.15 K)</CardTitle>
                    </CardHeader>
                    <CardContent className="p-0 pb-4">
                      <div className="overflow-x-auto">
                        <table className="w-full text-sm">
                          <thead>
                            <tr className="border-b text-xs text-foreground bg-muted/20">
                              <th className="px-4 py-2 text-center font-semibold">Species</th>
                              <th className="px-3 py-2 text-center font-semibold">Formula</th>
                              <th className="px-3 py-2 text-center font-semibold">Coefficient</th>
                              <th className="px-3 py-2 text-center font-semibold whitespace-nowrap">&Delta;H&deg;<sub>f</sub> ({conv.label})</th>
                              <th className="px-4 py-2 text-center font-semibold whitespace-nowrap">Contribution ({conv.label})</th>
                            </tr>
                          </thead>
                          <tbody>
                            {result.contributions.map((c, i) => (
                              <tr key={i} className="border-b last:border-0 hover:bg-muted/30 transition-colors">
                                <td className="px-4 py-2 text-xs max-w-[140px] truncate text-center" title={c.name}>{formatCompoundName(c.name)}</td>
                                <td className="px-3 py-2 font-mono text-xs text-center">
                                  {renderFormula(c.formula)}
                                  {c.phase298 && <span className="ml-0.5 text-muted-foreground font-normal text-[9px] font-sans">({c.phase298 === 'gas' ? 'g' : c.phase298 === 'liquid' ? 'l' : 's'})</span>}
                                </td>
                                <td className="px-3 py-2 text-center text-xs tabular-nums">{c.side === 'reactant' ? '\u2212' : '+'}{c.coefficient}</td>
                                <td className="px-3 py-2 text-center tabular-nums text-xs">{fmtNum(c.hof_kjmol * df)}</td>
                                <td className={`px-4 py-2 text-center tabular-nums text-xs font-medium ${c.contribution_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                  {c.contribution_kjmol >= 0 ? '+' : ''}{fmtNum(c.contribution_kjmol * dfc)}
                                </td>
                              </tr>
                            ))}
                            <tr className="border-t-2 bg-muted/20">
                              <td className="px-4 py-3 font-bold text-sm text-center">
                                <span className="inline-flex items-baseline">&Delta;H&deg;<sub className="text-[0.65em] leading-none">rxn</sub></span>
                                {tempUnit === 'C' ? ' (25 \u00b0C)' : ' (298 K)'}
                              </td>
                              <td colSpan={3} />
                              <td className={`px-4 py-3 text-center tabular-nums font-bold text-sm ${result.dH_rxn_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                {result.dH_rxn_kjmol >= 0 ? '+' : ''}{fmtNum(result.dH_rxn_kjmol * dfc)} {conv.label}
                              </td>
                            </tr>
                          </tbody>
                        </table>
                      </div>
                    </CardContent>
                  </Card>
                )}

                {result.dGf_rxn_kjmol != null && (() => {
                  const R_kJ = 8.314e-3; // kJ/(mol·K)
                  const T0 = 298.15;
                  const dG = result.dGf_rxn_kjmol;
                  const Kp298 = Math.exp(-dG / (R_kJ * T0));
                  const T2 = result.T2_K;
                  const Kp_T2 = T2 && Math.abs(T2 - T0) > 1
                    ? Kp298 * Math.exp(-(result.dH_rxn_kjmol / R_kJ) * (1/T2 - 1/T0))
                    : null;
                  const fmtK = (K: number) => K > 1e6 || K < 1e-6 ? K.toExponential(2) : K.toPrecision(3);
                  return (
                    <Card>
                      <CardHeader className="pb-1 cursor-pointer select-none" onClick={() => setEqConvOpen(v => !v)}>
                        <div className="flex items-center justify-between">
                          <CardTitle className="text-sm">&Delta;G&deg;<sub>rxn</sub> &amp; Equilibrium Constants</CardTitle>
                          <span className="text-muted-foreground text-xs">{eqConvOpen ? '\u25b2' : '\u25bc'}</span>
                        </div>
                      </CardHeader>
                      {eqConvOpen && (
                        <CardContent className="pb-4">
                          <div className="grid grid-cols-2 gap-2.5 text-xs">
                            <div className="rounded-md bg-muted/30 p-2.5">
                              <p className="text-muted-foreground mb-0.5">&Delta;G&deg;<sub>rxn</sub>(298 K)</p>
                              <p className={`font-mono font-semibold ${dG < 0 ? 'text-green-400' : 'text-red-400'}`}>{dG >= 0 ? '+' : ''}{fmtNum(dG)} kJ/mol</p>
                              <p className="text-muted-foreground text-[10px] mt-0.5">{dG < 0 ? 'Spontaneous at 298 K' : 'Non-spontaneous at 298 K'}</p>
                            </div>
                            <div className="rounded-md bg-muted/30 p-2.5">
                              <p className="text-muted-foreground mb-0.5">K<sub>p</sub>(298 K)</p>
                              <p className="font-mono font-semibold">{fmtK(Kp298)}</p>
                              <p className="text-muted-foreground text-[10px] mt-0.5">{Kp298 > 100 ? 'Strongly favors products' : Kp298 < 0.01 ? 'Strongly favors reactants' : 'Mixed equilibrium'}</p>
                            </div>
                            {Kp_T2 != null && (
                              <>
                                <div className="rounded-md bg-muted/30 p-2.5">
                                  <p className="text-muted-foreground mb-0.5">K<sub>p</sub>({fmtT(T2!)})</p>
                                  <p className="font-mono font-semibold">{fmtK(Kp_T2)}</p>
                                  <p className="text-muted-foreground text-[10px] mt-0.5">Van&apos;t Hoff (const &Delta;H)</p>
                                </div>
                                <div className="rounded-md bg-muted/30 p-2.5">
                                  <p className="text-muted-foreground mb-0.5">&Delta;G&deg;<sub>rxn</sub>({fmtT(T2!)})</p>
                                  <p className="font-mono font-semibold">{fmtNum(-R_kJ * 1000 * T2! * Math.log(Kp_T2))} J/mol</p>
                                  <p className="text-muted-foreground text-[10px] mt-0.5">Approx. (const &Delta;H)</p>
                                </div>
                              </>
                            )}
                          </div>
                          <p className="text-[10px] text-muted-foreground italic mt-2">Joback ideal-gas &Delta;G&deg;<sub>f</sub> values. Gas-phase only; liquid/solid-phase requires solvation corrections.</p>
                        </CardContent>
                      )}
                    </Card>
                  );
                })()}

                {result.dGf_rxn_kjmol != null && (() => {
                  const R_kJ = 8.314e-3; // kJ/(mol·K)
                  const T0 = 298.15;
                  const dG = result.dGf_rxn_kjmol;
                  const Kp298 = Math.exp(-dG / (R_kJ * T0));
                  const T2 = result.T2_K;
                  const Kp_T2 = T2 && Math.abs(T2 - T0) > 1
                    ? Kp298 * Math.exp(-(result.dH_rxn_kjmol / R_kJ) * (1/T2 - 1/T0))
                    : null;
                  const fmtK = (K: number) => K > 1e6 || K < 1e-6 ? K.toExponential(2) : K.toPrecision(3);
                  return (
                    <Card>
                      <CardHeader className="pb-1 cursor-pointer select-none" onClick={() => setEqConvOpen(v => !v)}>
                        <div className="flex items-center justify-between">
                          <CardTitle className="text-sm">&Delta;G&deg;<sub>rxn</sub> &amp; Equilibrium Constants</CardTitle>
                          <span className="text-muted-foreground text-xs">{eqConvOpen ? '▲' : '▼'}</span>
                        </div>
                      </CardHeader>
                      {eqConvOpen && (
                        <CardContent className="pb-4">
                          <div className="grid grid-cols-2 gap-2.5 text-xs">
                            <div className="rounded-md bg-muted/30 p-2.5">
                              <p className="text-muted-foreground mb-0.5">&Delta;G&deg;<sub>rxn</sub>(298 K)</p>
                              <p className={`font-mono font-semibold ${dG < 0 ? 'text-green-400' : 'text-red-400'}`}>{dG >= 0 ? '+' : ''}{fmtNum(dG)} kJ/mol</p>
                              <p className="text-muted-foreground text-[10px] mt-0.5">{dG < 0 ? 'Spontaneous at 298 K' : 'Non-spontaneous at 298 K'}</p>
                            </div>
                            <div className="rounded-md bg-muted/30 p-2.5">
                              <p className="text-muted-foreground mb-0.5">K<sub>p</sub>(298 K)</p>
                              <p className="font-mono font-semibold">{fmtK(Kp298)}</p>
                              <p className="text-muted-foreground text-[10px] mt-0.5">{Kp298 > 100 ? 'Strongly favors products' : Kp298 < 0.01 ? 'Strongly favors reactants' : 'Mixed equilibrium'}</p>
                            </div>
                            {Kp_T2 != null && (
                              <>
                                <div className="rounded-md bg-muted/30 p-2.5">
                                  <p className="text-muted-foreground mb-0.5">K<sub>p</sub>({fmtT(T2!)})</p>
                                  <p className="font-mono font-semibold">{fmtK(Kp_T2)}</p>
                                  <p className="text-muted-foreground text-[10px] mt-0.5">Van&apos;t Hoff (const &Delta;H)</p>
                                </div>
                                <div className="rounded-md bg-muted/30 p-2.5">
                                  <p className="text-muted-foreground mb-0.5">&Delta;G&deg;<sub>rxn</sub>({fmtT(T2!)})</p>
                                  <p className="font-mono font-semibold">{fmtNum(-R_kJ * 1000 * T2! * Math.log(Kp_T2))} J/mol</p>
                                  <p className="text-muted-foreground text-[10px] mt-0.5">Approx. (const &Delta;H)</p>
                                </div>
                              </>
                            )}
                          </div>
                          <p className="text-[10px] text-muted-foreground italic mt-2">Joback ideal-gas &Delta;G&deg;<sub>f</sub> values. Gas-phase only; liquid/solid requires solvation corrections.</p>
                        </CardContent>
                      )}
                    </Card>
                  );
                })()}

                <Card>
                  <CardHeader className="pb-1"><CardTitle className="text-sm">Enthalpy Component Breakdown</CardTitle></CardHeader>
                  <CardContent className="pt-4 pb-2">
                    <ReactECharts ref={echartsRef} option={combinedChartOption} style={{ height: 300 }} />
                    {/* Custom split-color legend */}
                    {(chartLegendMeta.hasSens || chartLegendMeta.hasLat) && (
                      <div className="flex flex-wrap justify-center gap-x-4 gap-y-1.5 pt-1" style={{ fontFamily: "'Merriweather Sans', sans-serif", fontSize: 10, paddingLeft: legendOffset.left, paddingRight: legendOffset.right }}>
                        {/* Formation */}
                        <span className="flex items-center gap-1.5">
                          <svg width="20" height="10" style={{ borderRadius: 2, display: 'inline-block', flexShrink: 0 }}>
                            <rect x="0" y="0" width="10" height="10" fill={CHR_FORM_R} />
                            <rect x="10" y="0" width="10" height="10" fill={CHR_FORM_P} />
                          </svg>
                          <span className="text-muted-foreground">&Delta;H&deg;<sub style={{ fontSize: '0.65em' }}>f</sub></span>
                        </span>
                        {/* Sensible */}
                        {chartLegendMeta.hasSens && (
                          <span className="flex items-center gap-1.5">
                            <svg width="20" height="10" style={{ borderRadius: 2, display: 'inline-block', flexShrink: 0 }}>
                              <rect x="0" y="0" width="10" height="10" fill={CHR_SENS_R} />
                              <rect x="10" y="0" width="10" height="10" fill={CHR_SENS_P} />
                            </svg>
                            <span className="text-muted-foreground">&Delta;H<sub style={{ fontSize: '0.65em' }}>sens</sub></span>
                          </span>
                        )}
                        {/* Latent */}
                        {chartLegendMeta.hasLat && (
                          <span className="flex items-center gap-1.5">
                            <svg width="20" height="10" style={{ borderRadius: 2, display: 'inline-block', flexShrink: 0 }}>
                              <rect x="0" y="0" width="10" height="10" fill={CHR_LAT_R} />
                              <rect x="10" y="0" width="10" height="10" fill={CHR_LAT_P} />
                            </svg>
                            <span className="text-muted-foreground">&Delta;H<sub style={{ fontSize: '0.65em' }}>lat</sub></span>
                          </span>
                        )}

                      </div>
                    )}
                  </CardContent>
                </Card>
              </div>
            ) : loading ? (
              <div className="space-y-4">
                <Card className="border-2"><CardContent className="pt-5 pb-5">
                  <Skeleton className="h-6 w-2/3 mx-auto mb-4" />
                  <Skeleton className="h-10 w-1/2 mx-auto mb-2" />
                  <Skeleton className="h-4 w-1/3 mx-auto" />
                </CardContent></Card>
                <Card><CardContent className="pt-4 pb-4 space-y-2">
                  <Skeleton className="h-4 w-1/2" />
                  <Skeleton className="h-32 w-full" />
                </CardContent></Card>
                <Card><CardContent className="pt-4 pb-4">
                  <Skeleton className="h-48 w-full" />
                </CardContent></Card>
              </div>
            ) : (
              <div className="flex h-full min-h-[400px] items-center justify-center">
                <p className="text-sm text-muted-foreground">Fill in reactants and products, then calculate.</p>
              </div>
            )}
          </div>

        </div>
      </div>
    </TooltipProvider>
  );
}

