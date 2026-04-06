'use client';

import React, { useState, useRef, useCallback, useMemo, useEffect } from 'react';
import { supabase } from '@/lib/supabaseClient';
import { useTheme } from 'next-themes';

import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import { BarChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';
import type { EChartsOption } from 'echarts';

echarts.use([TitleComponent, TooltipComponent, GridComponent, LegendComponent, BarChart, CanvasRenderer]);

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

import { PlusCircle, Trash2, Flame, Snowflake, AlertCircle } from 'lucide-react';

// ─── Types ────────────────────────────────────────────────────────────────────
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
    speciesName: string;
    formula: string;
    coefficient: number;
    side: 'reactant' | 'product';
    segments: HeatSegment[];
    total_kjmol: number;
    missingData: boolean;
  }[];
  dH_total_kjmol: number;
  excessSegments: {
    speciesName: string;
    formula: string;
    excessCoefficient: number;
    segments: HeatSegment[];
    total_kjmol: number;
    missingData: boolean;
  }[];
  contributions: {
    formula: string;
    name: string;
    coefficient: number;
    hof_kjmol: number;
    contribution_kjmol: number;
    side: 'reactant' | 'product';
  }[];
  equation: string;
  isExothermic: boolean;
  T1_K: number | null;
  T2_K: number | null;
  hasSensibleHeat: boolean;
  missingHessData: string[];
  missingSensibleData: boolean;
  limitingFormula: string | null;
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
  const trimmed = raw.replace(/-[A-Z][a-z]?\d*$/, '');
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

function computeSpeciesSensibleHeat(
  shParams: ShomateRow[],
  pt: PhaseTransitionRow | null,
  compound: NistCompound,
  T1_K: number,
  T2_K: number
): { segments: HeatSegment[]; total_kjmol: number; missingData: boolean } {
  const tboil = pt?.tboil_k ?? compound.tboil_k ?? null;
  const tfus  = pt?.tfus_k  ?? compound.tfus_k  ?? null;
  const hvap  = pt?.hvap_kjmol ?? null;
  const hfus  = pt?.hfus_kjmol ?? null;

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
      if (clampT2 < segT2 - 1) { warning = `Cp(T) data ends at ${sh1.t_max_k.toFixed(0)} K`; missingData = true; }
    } else {
      warning = `No Cp(T) data for ${phase} phase`; missingData = true;
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
          warning: hfus === null ? '\u0394Hfus not in database' : undefined,
        });
        if (hfus === null) missingData = true;
      }
      if (tboil !== null && Math.abs(transT - tboil) < 1) {
        segments.push({
          label: `Vaporization at ${tboil.toFixed(1)} K`,
          T1_K: transT, T2_K: transT,
          dH_kjmol: hvap !== null ? hvap * direction : 0,
          type: 'latent',
          warning: hvap === null ? '\u0394Hvap not in database' : undefined,
        });
        if (hvap === null) missingData = true;
      }
    }
  }

  return { segments, total_kjmol: segments.reduce((s, x) => s + x.dH_kjmol, 0), missingData };
}
// ─── Name Formatting ──────────────────────────────────────────────────────────
const FILLER_WORDS = new Set(['and','of','the','in','at','on','a','an','with','for','to','by','from','into','per','bis','tris']);
function formatCompoundName(raw: string): string {
  return raw.trim().replace(/\s+/g, ' ').split(' ').map((word, i) => {
    if (/^\([a-z]+\)$/.test(word)) return word; // (s), (g), (aq) unchanged
    const lower = word.toLowerCase();
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
  const [perMoleMode, setPerMoleMode] = useState(true); // true = per mol/kg basis; false = absolute (total heat for given quantity)
  // Quantities committed on Enter/Calculate — so live typing doesn't update the chart
  const [committedQtys, setCommittedQtys] = useState<Map<string, string>>(new Map());
  useEffect(() => { setMounted(true); }, []);

  const wrapperRefs  = useRef<Map<string, HTMLDivElement>>(new Map());
  const searchTimers = useRef<Map<string, ReturnType<typeof setTimeout>>>(new Map());
  const autoCalcInit = useRef(false);
  const autoCalcDone = useRef(false);

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
    const base = () => supabase.from('nist_compounds').select(cols)
      .eq('has_shomate', true).eq('has_phase_change', true).not('hf_gas_298_kjmol', 'is', null);
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
    const score = (c: NistCompound): number => {
      const nn = normKey(c.name);
      const nf = normKey(c.formula);
      if (nf === normQ) return 400;             // exact formula match — highest priority
      if (nn === normQ) return 300;             // exact name match
      if (nf.startsWith(normQ)) return 250 - nf.length; // formula prefix: shorter wins
      if (nn.startsWith(normQ)) return 200;
      const words = nn.split(/[\s,()[\]\-]+/);
      if (words.some(w => w.startsWith(normQ))) return 100;
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
          });
        }
      };
      processHess(validReactants, 'reactant');
      processHess(validProducts,  'product');

      if (contributions.length === 0) {
        setCalcError('None of the selected species have \u0394H\u00b0f(298 K) in the NIST database.');
        setLoading(false); return;
      }

      const dH_rxn_kjmol = contributions.reduce((s, c) => s + c.contribution_kjmol, 0);

      let T1_K: number | null = null;
      let T2_K: number | null = null;
      const sensibleSegments: CalcResult['sensibleSegments'] = [];
      const excessSegments:   CalcResult['excessSegments']   = [];
      let missingSensibleData = false;

      const t1n = parseFloat(t1Str); const t2n = parseFloat(t2Str);
      const wantSensible = enableTemp && !isNaN(t1n) && !isNaN(t2n) && Math.abs(t1n - t2n) > 0.01;

      if (wantSensible) {
        T1_K = toKelvin(t1n, tempUnit);
        T2_K = toKelvin(t2n, tempUnit);
        if (T1_K <= 0 || T2_K <= 0) { setCalcError('Temperatures must be above 0 K.'); setLoading(false); return; }

        const T_ref = 298.15;
        const allSpecies = [
          ...validReactants.map(r => ({ ...r, side: 'reactant' as const })),
          ...validProducts.map(r => ({  ...r, side: 'product'  as const })),
        ];

        await Promise.all(allSpecies.map(async (sp) => {
          if (!sp.selected || !sp.selected.has_shomate) return;
          const coeff = parseFloat(sp.coefficient);
          if (isNaN(coeff) || coeff <= 0) return;

          const [shRes, ptRes] = await Promise.all([
            supabase.from('nist_shomate').select('phase,eqno,t_min_k,t_max_k,a,b,c,d,e,f,g').eq('nist_id', sp.selected.nist_id),
            supabase.from('nist_phase_transitions').select('tboil_k,tfus_k,hvap_kjmol,hfus_kjmol').eq('nist_id', sp.selected.nist_id).maybeSingle(),
          ]);

          const shParams = (shRes.data ?? []) as ShomateRow[];
          const pt = (ptRes.data ?? null) as PhaseTransitionRow | null;

          const [sT1, sT2] = sp.side === 'reactant'
            ? [T1_K!, T_ref]   // actual direction: inlet → reference
            : [T_ref, T2_K!];  // actual direction: reference → outlet

          const { segments, total_kjmol, missingData } = computeSpeciesSensibleHeat(shParams, pt, sp.selected, sT1, sT2);
          if (missingData) missingSensibleData = true;

          sensibleSegments.push({
            speciesName: sp.selected.name, formula: sp.selected.formula,
            coefficient: coeff, side: sp.side, segments, total_kjmol, missingData,
          });
        }));

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

            const [shRes, ptRes] = await Promise.all([
              supabase.from('nist_shomate').select('phase,eqno,t_min_k,t_max_k,a,b,c,d,e,f,g').eq('nist_id', sp.selected.nist_id),
              supabase.from('nist_phase_transitions').select('tboil_k,tfus_k,hvap_kjmol,hfus_kjmol').eq('nist_id', sp.selected.nist_id).maybeSingle(),
            ]);
            const shParams2 = (shRes.data ?? []) as ShomateRow[];
            const pt2 = (ptRes.data ?? null) as PhaseTransitionRow | null;
            const { segments: exSegs, total_kjmol: exTotal, missingData: exMissing } =
              computeSpeciesSensibleHeat(shParams2, pt2, sp.selected, T_ref, T2_K!);
            if (exMissing) missingSensibleData = true;
            excessSegments.push({
              speciesName: sp.selected.name, formula: sp.selected.formula,
              excessCoefficient: excessCoeff, segments: exSegs, total_kjmol: exTotal, missingData: exMissing,
            });
          }));
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

      setResult({
        dH_rxn_kjmol, sensibleSegments, excessSegments, dH_total_kjmol, contributions, equation,
        isExothermic: displayDH < 0, T1_K, T2_K,
        hasSensibleHeat: wantSensible && sensibleSegments.length > 0,
        missingHessData, missingSensibleData,
        limitingFormula: calcLimiting(),
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
      tryCalc();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [conversion, result]);

  // When enableTemp changes (user toggles temp switch), re-run calc so the right side
  // updates immediately without needing to press Calculate.
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
          const sensRow = result.sensibleSegments.find(s => s.formula === ex.formula && s.side === 'reactant');
          newCoeff = sensRow ? Math.max(0, sensRow.coefficient * (1 - cf)) : 0;
        }
      } else {
        const sensRow = result.sensibleSegments.find(s => s.formula === ex.formula && s.side === 'reactant');
        newCoeff = sensRow ? Math.max(0, sensRow.coefficient * (1 - cf)) : 0;
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
    if (enableTemp && result.hasSensibleHeat) {
      for (const s of result.sensibleSegments) {
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
  }, [result, perMoleMode, heatUnit, unit, absQtyScale, enableTemp, reactants]);

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
    const sensRSeries  = entries.map(e => e.sensRVal);
    const latRSeries   = entries.map(e => e.latRVal);
    const exSensRSeries = entries.map(e => e.exSensRVal);
    const sensPSeries  = entries.map(e => e.sensPVal);
    const latPSeries   = entries.map(e => e.latPVal);

    const hasSensR   = sensRSeries.some(v => v !== null);
    const hasLatR    = latRSeries.some(v => v !== null);
    const hasExSensR = exSensRSeries.some(v => v !== null);
    const hasSensP   = sensPSeries.some(v => v !== null);
    const hasLatP    = latPSeries.some(v => v !== null);
    const hasThermal = hasSensR || hasLatR || hasExSensR || hasSensP || hasLatP;

    const rPurple   = '#a78bfa';
    const pTeal     = '#2dd4bf';
    const rPurpleL  = '#c4b5fd'; // sensible (reacted reactant)
    const rAmber    = '#fbbf24'; // sensible (excess reactant) — stacks on rPurpleL
    const pTealL    = '#99f6e4'; // sensible (product)
    const latentColor = '#f472b6';

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

    // Legend items with explicit colors
    const legendItems: { name: string; itemStyle: { color: string } }[] = [
      { name: '\u0394H\u00b0f (R)', itemStyle: { color: rPurple } },
      { name: '\u0394H\u00b0f (P)', itemStyle: { color: pTeal } },
      ...(hasSensR   ? [{ name: 'Sensible (R)',      itemStyle: { color: rPurpleL } }] : []),
      ...(hasExSensR ? [{ name: 'Sensible-Excess (R)', itemStyle: { color: rAmber } }] : []),
      ...(hasLatR    ? [{ name: 'Latent (R)',         itemStyle: { color: latentColor } }] : []),
      ...(hasSensP   ? [{ name: 'Sensible (P)',       itemStyle: { color: pTealL } }] : []),
      ...(hasLatP    ? [{ name: 'Latent (P)',         itemStyle: { color: latentColor } }] : []),
    ];

    const tooltipFn = (p: any) => {
      if (!entries[p.dataIndex]) return '';
      const val = p.value;
      if (val === null || val === undefined) return '';
      return `<b>${entries[p.dataIndex].formula}</b><br/>${p.seriesName}: <b>${val >= 0 ? '+' : ''}${fmtNum(val)} ${convC.label}</b>`;
    };

    const mkBar = (data: (number|null)[], color: string, stack?: string) =>
      ({ type: 'bar' as const, barMaxWidth: 20, color, stack,
         data: data.map(v => v === null ? null : { value: v, itemStyle: { color, borderRadius: v >= 0 ? [4,4,0,0] : [0,0,4,4] } }) as any });

    const bottomPad = hasThermal ? 56 : 32;

    return {
      animation: false,
      backgroundColor: 'transparent',
      textStyle: { fontFamily: ff },
      legend: hasThermal ? {
        data: legendItems,
        bottom: 4, left: 'center',
        textStyle: { color: lc, fontSize: 10, fontFamily: ff },
        itemWidth: 12, itemHeight: 8,
        type: 'scroll' as any,
      } : undefined,
      grid: { left: 64, right: 20, top: 28, bottom: bottomPad, containLabel: true },
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
        axisLabel: { color: (value: number) => value > 0 ? '#fb923c' : value < 0 ? '#38bdf8' : lc, fontSize: 10, fontFamily: ff, formatter: yFmt },
        axisLine: { show: true, lineStyle: { color: ac } },
        splitLine: { lineStyle: { color: sc } },
        min: yAxisBounds?.min,
        max: yAxisBounds?.max,
      },
      series: [
        {
          name: '\u0394H\u00b0f (R)',
          ...mkBar(hessRSeries, rPurple),
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
        { name: '\u0394H\u00b0f (P)', ...mkBar(hessPSeries, pTeal) },
        ...(hasSensR   ? [{ name: 'Sensible (R)',       ...mkBar(sensRSeries,   rPurpleL, 'thermalR') }] : []),
        ...(hasExSensR ? [{ name: 'Sensible-Excess (R)', ...mkBar(exSensRSeries, rAmber,   'thermalR') }] : []),
        ...(hasLatR    ? [{ name: 'Latent (R)',          ...mkBar(latRSeries,    latentColor) }] : []),
        ...(hasSensP   ? [{ name: 'Sensible (P)',        ...mkBar(sensPSeries,   pTealL,   'thermalP') }] : []),
        ...(hasLatP    ? [{ name: 'Latent (P)',          ...mkBar(latPSeries,    latentColor) }] : []),
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

  // ─── Species Row Renderer ─────────────────────────────────────────────────

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
        {rows.map(row => (
          <div key={row.id} className="flex gap-2 items-center" ref={el => { if (el) wrapperRefs.current.set(row.id, el); }}>
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
                      <p className="text-sm font-medium leading-tight">{formatCompoundName(s.name)}</p>
                      <div className="flex justify-between text-xs text-muted-foreground mt-0.5 gap-2">
                        <span className="font-mono">{renderFormula(s.formula)}</span>
                        <span className="flex gap-2 items-center">
                            {s.hf_gas_298_kjmol !== null && (() => {
                            // ΔHf is always molar — convert with /mol units; fallback to kJ/mol if user is in /kg mode
                            const molConv = UNIT_CONVERSIONS[unit] ?? { factor: 1, label: 'kJ/mol' };
                            return <span>&#x394;H&#xb0;<sub>f</sub> = {fmtNum(s.hf_gas_298_kjmol * molConv.factor)} {molConv.label}</span>;
                          })()}

                        </span>
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
        ))}
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
                  if (!isNaN(t1K) && t1K <= 0) return <p className="text-xs text-destructive">T_in must be above 0 K ({tempUnit==='C'?'min: −273.15 °C':'min: 0 K'})</p>;
                  if (!isNaN(t2K) && t2K <= 0) return <p className="text-xs text-destructive">T_out must be above 0 K ({tempUnit==='C'?'min: −273.15 °C':'min: 0 K'})</p>;
                  return null;
                })()}
              </CardContent>
            </Card>

            <Button className="w-full h-10" onClick={calculate} disabled={loading}>
              {loading ? 'Calculating\u2026' : <><span className="inline-flex items-baseline">Calculate &Delta;H&deg;<sub className="text-[0.65em] leading-none">tot</sub></span></>}
            </Button>

            {calcError && (
              <div className="flex items-start gap-2 rounded-md border border-destructive/50 bg-destructive/10 p-3 text-sm text-destructive">
                <AlertCircle className="h-4 w-4 mt-0.5 shrink-0" /><span>{calcError}</span>
              </div>
            )}
          </div>

          {/* Right: Results */}
          <div className="lg:col-span-2 space-y-4">
            {result ? (
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
                                    <span className="text-yellow-400 underline underline-offset-2 cursor-default">{renderFormula(c.formula)}</span>
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
                          <span className="text-xl text-muted-foreground font-mono">&Delta;H&deg;<sub className="text-[0.65em] leading-none">tot</sub> =</span>
                          <span className="text-2xl font-bold tabular-nums tracking-tight">
                            {displayDH >= 0 ? '+' : ''}{fmtNum(displayDH)}
                          </span>
                          <span className="text-xl text-muted-foreground">{conv.label}</span>
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
                      <div className="flex items-center justify-between">
                        <CardTitle className="text-sm">Thermochemical Enthalpy Pathway</CardTitle>
                        <span className="text-muted-foreground text-xs ml-2">{kirchhoffOpen ? '▲' : '▼'}</span>
                      </div>
                    </CardHeader>
                    <CardContent className={`space-y-5 pb-4 ${kirchhoffOpen ? '' : 'hidden'}`}>

                      {/* Shared table renderer */}
                      {(() => {
                        const renderSegTable = (
                          rows: { coeff: number; formula: string; colorClass: string; segments: HeatSegment[]; missingData: boolean }[],
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
                                    return visSegs.map((seg, si) => (
                                      <tr key={`${ri}-${si}`} className="border-b last:border-0 hover:bg-muted/20">
                                        <td className={`px-2 py-1.5 font-mono ${r.colorClass} whitespace-nowrap overflow-hidden text-ellipsis`}>
                                            {renderFormula(r.formula)}
                                            {r.missingData && <span className="ml-1 text-yellow-500" title="Missing data">&#x26A0;</span>}
                                          </td>
                                          <td className="px-2 py-1.5 text-muted-foreground whitespace-nowrap overflow-hidden text-ellipsis">
                                            {seg.type === 'latent'
                                              ? `${seg.label.includes('Vapor') ? 'Vaporization' : 'Fusion'} at ${fmtT(seg.T1_K)}`
                                              : `${fmtT(seg.T1_K)} \u2192 ${fmtT(seg.T2_K)}`}
                                            {seg.warning && <span className="ml-1 text-yellow-500" title={seg.warning}>&#x26A0;</span>}
                                          </td>
                                          <td className="px-2 py-1.5 text-center">
                                            <span className={`text-[10px] ${seg.type === 'latent' ? 'text-purple-400' : 'text-green-400'}`}>{seg.type === 'latent' ? 'Latent' : 'Sensible'}</span>
                                          </td>
                                          <td className={`px-2 py-1.5 text-right tabular-nums font-medium whitespace-nowrap ${seg.dH_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                            {seg.dH_kjmol >= 0 ? '+' : ''}{fmtNum(seg.dH_kjmol * r.coeff * scale)}
                                          </td>
                                      </tr>
                                    ));
                                  })}
                                </tbody>
                              </table>
                            </div>
                          );
                        };

                        const rSegs = result.sensibleSegments.filter(sp => sp.side === 'reactant');
                        const pSegs = result.sensibleSegments.filter(sp => sp.side === 'product');
                        const showExcess = reactiveExcessCoeffs.size > 0 && result.hasSensibleHeat;
                        const hasReactantHeat = rSegs.some(sp => sp.segments.some(seg => seg.type === 'latent' || Math.abs(seg.T2_K - seg.T1_K) > 0.5));

                        return (
                          <>
                            {/* Step 1 — split into 1a/1b when excess exists */}
                            {showExcess ? (
                              <div>
                                <p className="text-xs font-semibold text-violet-400 mb-2 flex items-center gap-1.5">
                                  <span className="inline-flex items-center justify-center w-5 h-4 rounded-full bg-violet-400/20 text-[10px] font-bold shrink-0">1a</span>
                                  Reactants sensible heat ({fmtT(result.T1_K!)} &rarr; {fmtT(298.15)})
                                </p>
                                {renderSegTable(
                                  rSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-violet-400', segments: sp.segments, missingData: sp.missingData })),
                                  <span className="text-xs text-muted-foreground italic"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">in</sub></span> = {fmtT(298.15)} &mdash; no sensible heat needed for reactants.</span>
                                )}
                                {hasReactantHeat && (() => {
                                  const total1 = rSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                                  return (
                                    <div className="mt-1 flex justify-end text-xs">
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
                                    Excess / unreacted &mdash; {fmtT(298.15)} &rarr; {fmtT(result.T2_K!)}
                                  </p>
                                  {result.excessSegments.some(sp => reactiveExcessCoeffs.has(sp.formula)) ? (
                                    <>
                                      {renderSegTable(
                                        result.excessSegments
                                          .filter(sp => reactiveExcessCoeffs.has(sp.formula))
                                          .map(sp => ({ coeff: +(reactiveExcessCoeffs.get(sp.formula)!).toFixed(4), formula: sp.formula, colorClass: 'text-lime-400', segments: sp.segments, missingData: sp.missingData })),
                                        <span className="text-xs text-muted-foreground italic">No excess heat needed.</span>,
                                        df
                                      )}
                                      {(() => {
                                        const totalEx = result.excessSegments
                                          .filter(sp => reactiveExcessCoeffs.has(sp.formula))
                                          .reduce((s, sp) => s + (reactiveExcessCoeffs.get(sp.formula) ?? 0) * sp.total_kjmol, 0);
                                        return (
                                          <div className="mt-1 flex justify-end text-xs">
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
                                  Reactants sensible heat ({fmtT(result.T1_K!)} &rarr; {fmtT(298.15)})
                                </p>
                                {renderSegTable(
                                  rSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-violet-400', segments: sp.segments, missingData: sp.missingData })),
                                  <span className="text-xs text-muted-foreground italic"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">in</sub></span> = {fmtT(298.15)} &mdash; no sensible heat needed for reactants.</span>
                                )}
                                {hasReactantHeat && (() => {
                                  const total1 = rSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                                  return (
                                    <div className="mt-1 flex justify-end text-xs">
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
                                        <td className={`px-2 py-1.5 font-mono whitespace-nowrap overflow-hidden text-ellipsis ${c.side === 'reactant' ? 'text-violet-400' : 'text-teal-400'}`}>{renderFormula(c.formula)}</td>
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
                                pSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-teal-400', segments: sp.segments, missingData: sp.missingData })),
                                <span className="text-xs text-muted-foreground italic"><span className="inline-flex items-baseline">T<sub className="text-[0.65em] leading-none">out</sub></span> = {fmtT(298.15)} &mdash; no sensible heat needed for products.</span>
                              )}
                              {pSegs.length > 0 && (() => {
                                const total3 = pSegs.reduce((s, sp) => s + sp.coefficient * sp.total_kjmol, 0);
                                return (
                                  <div className="mt-1 flex justify-end text-xs">
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
                                <div className="border-t-2 mt-2 pt-3 flex items-center justify-between">
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
                                <td className="px-3 py-2 font-mono text-xs text-center">{renderFormula(c.formula)}</td>
                                <td className="px-3 py-2 text-center text-xs tabular-nums">{c.side === 'reactant' ? '\u2212' : '+'}{c.coefficient}</td>
                                <td className="px-3 py-2 text-center tabular-nums text-xs">{fmtNum(c.hof_kjmol * df)}</td>
                                <td className={`px-4 py-2 text-center tabular-nums text-xs font-medium ${c.contribution_kjmol >= 0 ? 'text-sky-400' : 'text-orange-400'}`}>
                                  {c.contribution_kjmol >= 0 ? '+' : ''}{fmtNum(c.contribution_kjmol * dfc)}
                                </td>
                              </tr>
                            ))}
                            <tr className="border-t-2 bg-muted/20">
                              <td className="px-4 py-3 font-bold text-sm text-center" colSpan={3}><span className="inline-flex items-baseline">&Delta;H&deg;<sub className="text-[0.65em] leading-none">rxn</sub></span> (298 K)</td>
                              <td className="px-3 py-2 text-center text-xs font-semibold">Total:</td>
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

                <Card>
                  <CardHeader className="pb-1"><CardTitle className="text-sm">Enthalpy Component Breakdown</CardTitle></CardHeader>
                  <CardContent><ReactECharts option={combinedChartOption} style={{ height: 300 }} /></CardContent>
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

