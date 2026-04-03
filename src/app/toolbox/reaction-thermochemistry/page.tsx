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
import { PlusCircle, Trash2, Flame, Snowflake, AlertCircle } from 'lucide-react';

// ─── Types ────────────────────────────────────────────────────────────────────

interface NistCompound {
  nist_id: string;
  name: string;
  formula: string;
  cas: string | null;
  mol_weight: number | null;
  hf_gas_298_kjmol: number | null;
  s_gas_298_jmolk: number | null;
  cp_gas_298_jmolk: number | null;
  tboil_k: number | null;
  tfus_k: number | null;
  has_shomate: boolean;
  has_phase_change: boolean;
}

interface ShomateRow {
  phase: string;
  t_min_k: number;
  t_max_k: number;
  a: number; b: number; c: number; d: number;
  e: number; f: number; g: number; h: number;
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

// ─── Pure Helpers ─────────────────────────────────────────────────────────────

let _idCounter = 0;
function makeId() { return `sr-${++_idCounter}`; }
function makeSpecies(): SpeciesRow {
  return { id: makeId(), coefficient: '1', query: '', selected: null, suggestions: [], showSuggestions: false, quantity: '' };
}

function fmtNum(val: number): string {
  const abs = Math.abs(val);
  if (abs === 0) return '0.000';
  if (abs < 0.001 || abs >= 1e6) return val.toExponential(3);
  if (abs < 1) return val.toFixed(4);
  return val.toFixed(3);
}

function toKelvin(val: number, unit: 'C' | 'K'): number {
  return unit === 'C' ? val + 273.15 : val;
}

function shomateH(sh: ShomateRow, T_K: number): number {
  const t = T_K / 1000;
  return sh.a * t + (sh.b * t * t) / 2 + (sh.c * t * t * t) / 3 +
    (sh.d * t * t * t * t) / 4 - sh.e / t + sh.f - sh.h;
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
        dH = shomateH(sh1, segT2) - shomateH(sh1, segT1);
      } else {
        const splitT = sh1.t_max_k;
        dH = (shomateH(sh1, splitT) - shomateH(sh1, segT1)) +
             (shomateH(sh2, segT2)  - shomateH(sh2, splitT));
      }
    } else if (sh1) {
      const clampT2 = Math.min(segT2, sh1.t_max_k);
      dH = shomateH(sh1, clampT2) - shomateH(sh1, segT1);
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
export default function ReactionThermochemistryPage() {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';

  const [reactants,   setReactants]   = useState<SpeciesRow[]>([makeSpecies()]);
  const [products,    setProducts]    = useState<SpeciesRow[]>([makeSpecies()]);
  const [unit,        setUnit]        = useState<string>('kJ/mol');
  const [heatUnit,    setHeatUnit]    = useState<'kJ/mol' | 'kJ/kg'>('kJ/mol');
  const [tempUnit,    setTempUnit]    = useState<'C' | 'K'>('C');
  const [enableTemp,  setEnableTemp]  = useState(false);
  const [t1Str,       setT1Str]       = useState('25');
  const [t2Str,       setT2Str]       = useState('100');
  const [result,      setResult]      = useState<CalcResult | null>(null);
  const [calcError,   setCalcError]   = useState<string | null>(null);
  const [loading,     setLoading]     = useState(false);
  const [conversion,  setConversion]  = useState(100);
  const [mounted,     setMounted]     = useState(false);
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

  const limitingReactantFormula = useMemo<string | null>(() => {
    const withQty = reactants.filter(r => r.selected && r.quantity.trim() !== '');
    if (withQty.length < 2) return null;
    const bases = withQty.flatMap(r => {
      const q = parseFloat(r.quantity);
      const coeff = parseFloat(r.coefficient);
      if (isNaN(q) || q <= 0 || isNaN(coeff) || coeff <= 0) return [];
      return [{ formula: r.selected!.formula, basis: q / coeff }];
    });
    if (bases.length < 2) return null;
    const min = Math.min(...bases.map(b => b.basis));
    return bases.find(b => b.basis === min)?.formula ?? null;
  }, [reactants]);

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
      .select('nist_id,name,formula,cas,mol_weight,hf_gas_298_kjmol,s_gas_298_jmolk,cp_gas_298_jmolk,tboil_k,tfus_k,has_shomate,has_phase_change')
      .in('nist_id', ids)
      .then(({ data }) => {
        if (!data || data.length === 0) return;
        const byId = new Map<string, NistCompound>((data as NistCompound[]).map(c => [c.nist_id, c]));
        const makeRow = (id: string, coeff: string): SpeciesRow => {
          const c = byId.get(id) ?? null;
          return { id: makeId(), coefficient: coeff, query: c?.name ?? '', selected: c, suggestions: [], showSuggestions: false, quantity: '' };
        };
        const r1 = makeRow('C74828', '1'); r1.quantity = '1';
        const r2 = makeRow('C7782447', '2'); r2.quantity = '2';
        setReactants([r1, r2]);
        setProducts ([makeRow('C124389', '1'), makeRow('C7732185', '2')]);
        setEnableTemp(true);
        setT1Str('25');
        setT2Str('150');
      });
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const patchRow = useCallback((side: 'reactants' | 'products', id: string, patch: Partial<SpeciesRow>) => {
    (side === 'reactants' ? setReactants : setProducts)(prev => prev.map(r => r.id === id ? { ...r, ...patch } : r));
  }, []);

  const fetchSuggestions = useCallback(async (side: 'reactants' | 'products', id: string, q: string) => {
    if (q.trim().length < 1) { patchRow(side, id, { suggestions: [], showSuggestions: false }); return; }
    const { data } = await supabase
      .from('nist_compounds')
      .select('nist_id,name,formula,cas,mol_weight,hf_gas_298_kjmol,s_gas_298_jmolk,cp_gas_298_jmolk,tboil_k,tfus_k,has_shomate,has_phase_change')
      .or(`name.ilike.%${q}%,formula.ilike.%${q}%`)
      // Only show compounds with ΔH°f data (or element convention), Shomate Cp(T), AND phase-change data
      .eq('has_shomate', true)
      .eq('has_phase_change', true)
      .limit(12);
    if (data) {
      const ql = q.toLowerCase();
      const sorted = [...(data as NistCompound[])].sort((a, b) => {
        const aP = a.name.toLowerCase().startsWith(ql) || a.formula.toLowerCase().startsWith(ql) ? 0 : 1;
        const bP = b.name.toLowerCase().startsWith(ql) || b.formula.toLowerCase().startsWith(ql) ? 0 : 1;
        return aP - bP;
      });
      patchRow(side, id, { suggestions: sorted, showSuggestions: true });
    }
  }, [patchRow]);

  const handleQueryChange = useCallback((side: 'reactants' | 'products', id: string, value: string) => {
    patchRow(side, id, { query: value, selected: null, showSuggestions: false });
    const existing = searchTimers.current.get(id);
    if (existing) clearTimeout(existing);
    searchTimers.current.set(id, setTimeout(() => fetchSuggestions(side, id, value), 220));
  }, [patchRow, fetchSuggestions]);

  const handleSelect = useCallback((side: 'reactants' | 'products', id: string, row: NistCompound) => {
    patchRow(side, id, { selected: row, query: row.name, suggestions: [], showSuggestions: false });
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
    const validReactants = reactants.filter(r => r.selected);
    const validProducts  = products.filter(p => p.selected);

    if (validReactants.length === 0 || validProducts.length === 0) {
      setCalcError('Please select at least one reactant and one product from the autocomplete dropdown.');
      return;
    }
    setLoading(true);
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
            supabase.from('nist_shomate').select('phase,t_min_k,t_max_k,a,b,c,d,e,f,g,h').eq('nist_id', sp.selected.nist_id),
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

        // Excess (unreacted) reactants travel T_ref → T2 at the same conditions as products
        const cf = conversion / 100;
        if (cf < 1) {
          await Promise.all(validReactants.map(async (sp) => {
            if (!sp.selected || !sp.selected.has_shomate) return;
            const coeff = parseFloat(sp.coefficient);
            if (isNaN(coeff) || coeff <= 0) return;

            const [shRes, ptRes] = await Promise.all([
              supabase.from('nist_shomate').select('phase,t_min_k,t_max_k,a,b,c,d,e,f,g,h').eq('nist_id', sp.selected.nist_id),
              supabase.from('nist_phase_transitions').select('tboil_k,tfus_k,hvap_kjmol,hfus_kjmol').eq('nist_id', sp.selected.nist_id).maybeSingle(),
            ]);
            const shParams2 = (shRes.data ?? []) as ShomateRow[];
            const pt2 = (ptRes.data ?? null) as PhaseTransitionRow | null;
            const excessCoeff = coeff * (1 - cf);
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
      const displayDH = wantSensible ? dH_total_kjmol : dH_rxn_kjmol;

      setResult({
        dH_rxn_kjmol, sensibleSegments, excessSegments, dH_total_kjmol, contributions, equation,
        isExothermic: displayDH < 0, T1_K, T2_K,
        hasSensibleHeat: wantSensible && sensibleSegments.length > 0,
        missingHessData, missingSensibleData,
      });
    } catch {
      setCalcError('An unexpected error occurred. Please try again.');
    }
    setLoading(false);
  }, [reactants, products, unit, enableTemp, t1Str, t2Str, tempUnit, conversion]);

  // Helper: trigger calculate if enough species are selected
  const tryCalc = useCallback(() => {
    if (reactants.some(r => r.selected) && products.some(p => p.selected)) calculate();
  }, [reactants, products, calculate]);

  // Enter key in any input on the left side triggers calculate
  const handleKeyDown = useCallback((e: React.KeyboardEvent) => {
    if (e.key === 'Enter') { e.preventDefault(); tryCalc(); }
  }, [tryCalc]);

  useEffect(() => {
    if (!autoCalcInit.current) { autoCalcInit.current = true; return; }
    if (reactants.some(r => r.selected) && products.some(p => p.selected)) calculate();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [unit]);

  useEffect(() => {
    if (!autoCalcInit.current) return;
    tryCalc();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [enableTemp]);

  useEffect(() => {
    if (autoCalcDone.current) return;
    if (reactants.some(r => r.selected) && products.some(p => p.selected)) {
      autoCalcDone.current = true;
      calculate();
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [reactants, products]);

  // ─── ECharts ──────────────────────────────────────────────────────────────

  const chartOption = useMemo((): EChartsOption => {
    if (!result || result.contributions.length === 0) return {};
    const activeConvMapC = heatUnit === 'kJ/kg' ? UNIT_CONVERSIONS_KG : UNIT_CONVERSIONS;
    const conv = activeConvMapC[unit] ?? { factor: 1, label: heatUnit === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol' };
    const pkf = heatUnit === 'kJ/kg'
      ? (() => { const r = reactants.find(row => row.selected?.mol_weight); return r?.selected?.mol_weight ? 1000 / r.selected.mol_weight : 1; })()
      : 1;
    const cf = conversion / 100;
    const numR  = result.contributions.filter(c => c.side === 'reactant').length;
    const hasGap = numR > 0 && numR < result.contributions.length;

    const chartLabels: string[] = [];
    const chartData: { value: number | null; itemStyle: { color: string; borderRadius?: number[] } }[] = [];
    result.contributions.forEach((c, i) => {
      if (i === numR && hasGap) { chartLabels.push(''); chartData.push({ value: null, itemStyle: { color: 'transparent' } }); }
      chartLabels.push(c.formula);
      const v = c.contribution_kjmol * conv.factor * pkf * cf;
      chartData.push({ value: v, itemStyle: { color: c.side === 'product' ? '#38bdf8' : '#EE6666', borderRadius: v >= 0 ? [4, 4, 0, 0] : [0, 0, 4, 4] } });
    });

    const values = result.contributions.map(c => c.contribution_kjmol * conv.factor * pkf * cf);
    const lc = isDark ? '#cbd5e1' : '#475569';
    const sc = isDark ? 'rgba(255,255,255,0.07)' : 'rgba(0,0,0,0.07)';
    const ac = isDark ? '#ffffff' : '#000000';
    const ff = "'Merriweather Sans', sans-serif";
    const yFmt = (v: number) => { const r = Math.round(v); const a = Math.abs(r); return a >= 1e6 ? `${(r/1e6).toFixed(1)}M` : a >= 1e3 ? `${(r/1e3).toFixed(1)}k` : `${r}`; };
    const maxAbs = Math.max(...values.map(Math.abs), 1);

    return {
      animation: false,
      backgroundColor: 'transparent', textStyle: { fontFamily: ff },
      grid: { left: 10, right: 30, top: 28, bottom: 10, containLabel: true },
      xAxis: { type: 'category', data: chartLabels,
        axisLabel: { color: lc, interval: 0, fontSize: 10, fontFamily: ff },
        axisLine: { lineStyle: { color: ac } },
        axisTick: { lineStyle: { color: ac }, interval: (_: number, v: string) => v !== '' } },
      yAxis: { type: 'value', name: `\u0394H\u00b0f contribution (${conv.label})`,
        nameLocation: 'middle', nameRotate: 90, nameGap: Math.max(50, yFmt(-maxAbs).length * 9 + 15),
        nameTextStyle: { color: lc, fontSize: 11, fontFamily: ff },
        axisLabel: { color: lc, fontSize: 11, fontFamily: ff, formatter: yFmt },
        axisLine: { show: false }, splitLine: { lineStyle: { color: sc } } },
      series: [{ type: 'bar', barMaxWidth: 48, data: chartData as any,
        markLine: hasGap ? { symbol: ['none','none'], silent: true,
          lineStyle: { type: 'dashed' as any, color: isDark ? 'rgba(255,255,255,0.3)' : 'rgba(0,0,0,0.3)', width: 1 },
          label: { show: false }, data: [{ xAxis: numR }] } : undefined,
        markArea: { silent: true, data: [
          [{ xAxis: 0, label: { show: true, formatter: 'Reactants', position: 'insideBottom' as any, color: '#EE6666', fontSize: 9, fontFamily: ff }, itemStyle: { color: 'rgba(238,102,102,0.06)' } }, { xAxis: numR }],
          [{ xAxis: numR, label: { show: true, formatter: 'Products', position: 'insideBottom' as any, color: '#38bdf8', fontSize: 9, fontFamily: ff }, itemStyle: { color: 'rgba(56,189,248,0.06)' } }, { xAxis: chartLabels.length - 1 }],
        ] as any } }],
      tooltip: { trigger: 'item',
        backgroundColor: isDark ? '#1e293b' : '#fff', borderColor: isDark ? '#334155' : '#e2e8f0',
        textStyle: { color: isDark ? '#e2e8f0' : '#1e293b', fontSize: 12, fontFamily: ff },
        formatter: (p: any) => {
          const ri: number = p.dataIndex;
          if (hasGap && ri === numR) return '';
          const ai = hasGap && ri > numR ? ri - 1 : ri;
          if (ai < 0 || ai >= result.contributions.length) return '';
          const c = result.contributions[ai];
          const sgn = c.side === 'product' ? '+' : '\u2212';
          const hof = fmtNum(c.hof_kjmol * conv.factor * pkf);
          return `<b>${c.name}</b><br/>${c.formula}<br/>\u0394H\u00b0<sub>f</sub> = ${hof} ${conv.label}<br/>${sgn}${c.coefficient} \u00d7 ${hof} = <b>${fmtNum(c.contribution_kjmol * conv.factor * pkf * cf)}</b>`;
        } },
    };
  }, [result, isDark, unit, heatUnit, reactants, conversion]);

  const sensibleChartOption = useMemo((): EChartsOption | null => {
    if (!result?.hasSensibleHeat) return null;
    const activeConvMapS = heatUnit === 'kJ/kg' ? UNIT_CONVERSIONS_KG : UNIT_CONVERSIONS;
    const conv = activeConvMapS[unit] ?? { factor: 1, label: heatUnit === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol' };
    const pkf = heatUnit === 'kJ/kg'
      ? (() => { const r = reactants.find(row => row.selected?.mol_weight); return r?.selected?.mol_weight ? 1000 / r.selected.mol_weight : 1; })()
      : 1;
    const cf = conversion / 100;
    const dfc = conv.factor * pkf * cf;
    const bars: { label: string; value: number; color: string }[] = [];
    for (const s of result.sensibleSegments.filter(x => x.side === 'reactant')) {
      const v = s.coefficient * s.total_kjmol * dfc;
      if (Math.abs(v) > 0.0001) bars.push({ label: s.formula, value: v, color: '#EE6666' });
    }
    bars.push({ label: '\u0394H\u00b0rxn', value: result.dH_rxn_kjmol * dfc, color: '#fbbf24' });
    for (const s of result.sensibleSegments.filter(x => x.side === 'product')) {
      const v = s.coefficient * s.total_kjmol * dfc;
      if (Math.abs(v) > 0.0001) bars.push({ label: s.formula, value: v, color: '#38bdf8' });
    }
    for (const s of result.excessSegments) {
      const v = s.excessCoefficient * s.total_kjmol * dfc;
      if (Math.abs(v) > 0.0001) bars.push({ label: `${s.formula}*`, value: v, color: '#f97316' });
    }
    const labels = bars.map(b => b.label);
    const chartData = bars.map(b => ({ value: b.value, itemStyle: { color: b.color, borderRadius: b.value >= 0 ? [4, 4, 0, 0] : [0, 0, 4, 4] as any } }));
    const maxAbs = Math.max(...bars.map(b => Math.abs(b.value)), 1);
    const lc = isDark ? '#cbd5e1' : '#475569';
    const sc = isDark ? 'rgba(255,255,255,0.07)' : 'rgba(0,0,0,0.07)';
    const ff = "'Merriweather Sans', sans-serif";
    const yFmt = (v: number) => { const r = Math.round(v); const a = Math.abs(r); return a >= 1e6 ? `${(r/1e6).toFixed(1)}M` : a >= 1e3 ? `${(r/1e3).toFixed(1)}k` : `${r}`; };
    return {
      animation: false,
      backgroundColor: 'transparent', textStyle: { fontFamily: ff },
      grid: { left: 10, right: 30, top: 36, bottom: 10, containLabel: true },
      xAxis: { type: 'category', data: labels,
        axisLabel: { color: lc, interval: 0, fontSize: 10, fontFamily: ff },
        axisLine: { lineStyle: { color: isDark ? '#ffffff' : '#000000' } },
        axisTick: { lineStyle: { color: isDark ? '#ffffff' : '#000000' } } },
      yAxis: { type: 'value', name: `\u0394H contribution (${conv.label})`,
        nameLocation: 'middle', nameRotate: 90, nameGap: Math.max(60, yFmt(-maxAbs).length * 9 + 15),
        nameTextStyle: { color: lc, fontSize: 11, fontFamily: ff },
        axisLabel: { color: lc, fontSize: 11, fontFamily: ff, formatter: yFmt },
        axisLine: { show: false }, splitLine: { lineStyle: { color: sc } } },
      series: [{ type: 'bar', barMaxWidth: 48, data: chartData as any,
        markLine: { symbol: ['none', 'none'], silent: true,
          lineStyle: { type: 'dashed' as any, color: isDark ? 'rgba(255,255,255,0.35)' : 'rgba(0,0,0,0.35)', width: 1.5 },
          label: { show: true, formatter: `\u0394H\u00b0tot = ${fmtNum(result.dH_total_kjmol * dfc)} ${conv.label}`, position: 'insideStartTop' as any, color: lc, fontSize: 9, fontFamily: ff },
          data: [{ yAxis: result.dH_total_kjmol * dfc }] } }],
      tooltip: { trigger: 'item',
        backgroundColor: isDark ? '#1e293b' : '#fff', borderColor: isDark ? '#334155' : '#e2e8f0',
        textStyle: { color: isDark ? '#e2e8f0' : '#1e293b', fontSize: 12, fontFamily: ff },
        formatter: (p: any) => {
          const b = bars[p.dataIndex];
          if (!b) return '';
          return `<b>${b.label}</b><br/>${fmtNum(b.value)} ${conv.label}`;
        } },
    };
  }, [result, isDark, unit, heatUnit, reactants, conversion]);

  // ─── Species Row Renderer ─────────────────────────────────────────────────

  const renderSide = (side: 'reactants' | 'products') => {
    const rows  = side === 'reactants' ? reactants : products;
    const label = side === 'reactants' ? 'Reactants' : 'Products';
    const showQty = side === 'reactants';
    return (
      <div className="space-y-2">
        <Label className="text-xs font-semibold text-muted-foreground uppercase tracking-wide">{label}</Label>
        <div className="flex gap-2 px-0.5">
          <span className="w-16 text-center text-xs text-muted-foreground shrink-0">Coeff.</span>
          <span className="flex-1 text-xs text-muted-foreground">Species</span>
          {showQty && <span className="w-20 text-center text-xs text-muted-foreground shrink-0">{heatUnit === 'kJ/kg' ? 'kg' : 'mol'}</span>}
          <span className="w-9 shrink-0" />
        </div>
        {rows.map(row => (
          <div key={row.id} className="flex gap-2 items-center" ref={el => { if (el) wrapperRefs.current.set(row.id, el); }}>
            <Input className="w-16 text-center h-9 shrink-0" placeholder="1" value={row.coefficient}
              onChange={e => patchRow(side, row.id, { coefficient: e.target.value })}
              onKeyDown={handleKeyDown} />
            <div className="relative flex-1">
              <Input className="h-9" placeholder="Name or formula\u2026" value={row.query}
                onChange={e => handleQueryChange(side, row.id, e.target.value)}
                onKeyDown={handleKeyDown}
                onFocus={() => { if (row.suggestions.length > 0) patchRow(side, row.id, { showSuggestions: true }); }} />
              {row.showSuggestions && row.suggestions.length > 0 && (
                <div className="absolute z-50 top-full left-0 right-0 mt-1 rounded-md border bg-popover shadow-md max-h-64 overflow-y-auto">
                  {row.suggestions.map(s => (
                    <div key={s.nist_id} className="px-3 py-2 hover:bg-accent cursor-pointer"
                      onMouseDown={e => { e.preventDefault(); handleSelect(side, row.id, s); }}>
                      <p className="text-sm font-medium leading-tight">{s.name}</p>
                      <div className="flex justify-between text-xs text-muted-foreground mt-0.5 gap-2">
                        <span className="font-mono">{s.formula}</span>
                        <span className="flex gap-2 items-center">
                          {s.hf_gas_298_kjmol !== null
                            ? <span>&#x394;H&#xb0;<sub>f</sub> {s.hf_gas_298_kjmol.toFixed(1)} kJ/mol</span>
                            : <span className="text-yellow-500 dark:text-yellow-400">element &#x394;H&#xb0;<sub>f</sub>=0</span>}

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

  const activeConvMap = heatUnit === 'kJ/kg' ? UNIT_CONVERSIONS_KG : UNIT_CONVERSIONS;
  const conv = activeConvMap[unit] ?? { factor: 1, label: heatUnit === 'kJ/kg' ? 'kJ/kg' : 'kJ/mol' };
  const conversionFactor = conversion / 100;
  const df  = conv.factor * perKgFactor;          // display factor for intrinsic values (ΔHf per species)
  const dfc = conv.factor * perKgFactor * conversionFactor; // display factor for extensive values
  const displayDH = result
    ? (result.hasSensibleHeat ? result.dH_total_kjmol : result.dH_rxn_kjmol) * dfc
    : 0;
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
                <CardTitle className="text-base">Build Reaction</CardTitle>
              </CardHeader>
              <CardContent className="space-y-3">
                <div className="grid grid-cols-2 gap-3">
                  <div className="space-y-1.5">
                    <Label className="text-xs text-muted-foreground">&Delta;H units</Label>
                    <Select value={unit} onValueChange={setUnit}>
                      <SelectTrigger className="h-9 w-full"><SelectValue /></SelectTrigger>
                      <SelectContent>
                        {heatUnit === 'kJ/kg' ? (
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
                  <CardTitle className="text-base">Temperature Change</CardTitle>
                  <div className="flex items-center gap-2 shrink-0">
                    <Tabs value={tempUnit} onValueChange={handleTempUnitToggle}>
                      <TabsList className="h-6 p-0.5">
                        <TabsTrigger value="C" className="h-5 px-2 text-xs">&deg;C</TabsTrigger>
                        <TabsTrigger value="K" className="h-5 px-2 text-xs">K</TabsTrigger>
                      </TabsList>
                    </Tabs>
                    <button
                      className={`relative inline-flex h-5 w-9 items-center rounded-full transition-colors ${enableTemp ? 'bg-primary' : 'bg-muted'}`}
                      onClick={() => setEnableTemp(v => !v)}
                      role="switch" aria-checked={enableTemp}>
                      <span className={`inline-block h-3.5 w-3.5 transform rounded-full bg-white transition-transform ${enableTemp ? 'translate-x-4' : 'translate-x-0.5'}`} />
                    </button>
                  </div>
                </div>
              </CardHeader>
              <CardContent className={`space-y-3 transition-opacity ${(!mounted || enableTemp) ? 'opacity-100' : 'opacity-40 pointer-events-none'}`}>
                <div className="grid grid-cols-2 gap-2">
                  <div className="space-y-1">
                    <Label className="text-xs text-muted-foreground">T<sub style={{verticalAlign:'-0.1em',fontSize:'0.8em'}}>in</sub></Label>
                    <Input className="h-9" value={t1Str} onChange={e => setT1Str(e.target.value)} onKeyDown={handleKeyDown} placeholder={tempUnit === 'C' ? '25' : '298'} />
                  </div>
                  <div className="space-y-1">
                    <Label className="text-xs text-muted-foreground">T<sub style={{verticalAlign:'-0.1em',fontSize:'0.8em'}}>out</sub></Label>
                    <Input className="h-9" value={t2Str} onChange={e => setT2Str(e.target.value)} onKeyDown={handleKeyDown} placeholder={tempUnit === 'C' ? '100' : '373'} />
                  </div>
                </div>
              </CardContent>
            </Card>

            <Button className="w-full h-10" onClick={calculate} disabled={loading}>
              {loading ? 'Calculating\u2026' : <>Calculate &Delta;H&deg;<sub style={{fontSize:'0.75em'}}>tot</sub></>}
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
                <Card className={`border-2 transition-colors ${result.isExothermic ? 'border-orange-400/50' : 'border-sky-400/50'}`}>
                  <CardContent className="pt-5 pb-5">
                    <div className="flex items-start justify-between gap-4">
                      <div className="flex-1 min-w-0 text-center">
                        <p className="text-xl font-mono font-semibold leading-snug break-words mb-3">
                          {result.contributions.filter(c => c.side === 'reactant').map((c, i) => (
                            <React.Fragment key={i}>
                              {i > 0 && <span className="mx-1.5 opacity-50">+</span>}
                              {c.coefficient !== 1 && <span className="opacity-60">{c.coefficient} </span>}
                              {limitingReactantFormula === c.formula ? (
                                <Tooltip>
                                  <TooltipTrigger asChild>
                                    <span className="text-yellow-400 underline underline-offset-2 cursor-default">{c.formula}</span>
                                  </TooltipTrigger>
                                  <TooltipContent side="top">Limiting Reactant</TooltipContent>
                                </Tooltip>
                              ) : <span>{c.formula}</span>}
                            </React.Fragment>
                          ))}
                          <span className="mx-2.5 opacity-50">&rarr;</span>
                          {result.contributions.filter(c => c.side === 'product').map((c, i) => (
                            <React.Fragment key={i}>
                              {i > 0 && <span className="mx-1.5 opacity-50">+</span>}
                              {c.coefficient !== 1 && <span className="opacity-60">{c.coefficient} </span>}
                              {c.formula}
                            </React.Fragment>
                          ))}
                        </p>
                        <div className="flex items-baseline gap-2 justify-center">
                          <span className="text-4xl font-bold tabular-nums tracking-tight">
                            {displayDH >= 0 ? '+' : ''}{fmtNum(displayDH)}
                          </span>
                          <span className="text-xl text-muted-foreground">{conv.label}</span>
                        </div>
                        <p className="text-xs text-muted-foreground mt-1.5">
                          {result.hasSensibleHeat
                            ? <>&Delta;H&deg;<sub>tot</sub> &mdash; {fmtT(result.T1_K!)} &rarr; {fmtT(result.T2_K!)}</>
                            : <>&Delta;H&deg;<sub>rxn</sub> at 298.15 K</>}
                        </p>
                      </div>
                      <div className={`flex flex-col items-center justify-center rounded-xl px-6 py-4 shrink-0 ${result.isExothermic ? 'bg-orange-400/10 text-orange-400' : 'bg-sky-400/10 text-sky-400'}`}>
                        {result.isExothermic ? <Flame className="h-9 w-9 mb-1" /> : <Snowflake className="h-9 w-9 mb-1" />}
                        <span className="font-semibold text-sm">{result.isExothermic ? 'Exothermic' : 'Endothermic'}</span>
                        <span className="text-xs opacity-70 mt-0.5">{result.isExothermic ? 'Releases heat' : 'Absorbs heat'}</span>
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
                    <CardHeader className="pb-2">
                      <CardTitle className="text-sm">Enthalpy Path &mdash; Kirchhoff + Hess&apos;s Law</CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-5 pb-4">

                      {/* Shared table renderer */}
                      {(() => {
                        const renderSegTable = (
                          rows: { coeff: number; formula: string; colorClass: string; segments: HeatSegment[]; missingData: boolean }[],
                          emptyMsg: string
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
                                  {rows.flatMap((r, ri) =>
                                    r.segments
                                      .filter(seg => seg.type === 'latent' || Math.abs(seg.T2_K - seg.T1_K) > 0.5)
                                      .map((seg, si) => (
                                        <tr key={`${ri}-${si}`} className="border-b last:border-0 hover:bg-muted/20">
                                          <td className={`px-2 py-1.5 font-mono ${r.colorClass} whitespace-nowrap overflow-hidden text-ellipsis`}>
                                            {r.formula}
                                            {r.missingData && <span className="ml-1 text-yellow-500" title="Missing data">&#x26A0;</span>}
                                          </td>
                                          <td className="px-2 py-1.5 text-muted-foreground whitespace-nowrap overflow-hidden text-ellipsis">
                                            {seg.type === 'latent' ? seg.label : `${fmtT(seg.T1_K)} \u2192 ${fmtT(seg.T2_K)}`}
                                            {seg.warning && <span className="ml-1 text-yellow-500" title={seg.warning}>&#x26A0;</span>}
                                          </td>
                                          <td className="px-2 py-1.5 text-center">
                                            <span className={`text-[10px] ${seg.type === 'latent' ? 'text-purple-400' : 'text-blue-400'}`}>{seg.type}</span>
                                          </td>
                                          <td className={`px-2 py-1.5 text-right tabular-nums font-medium whitespace-nowrap ${seg.dH_kjmol >= 0 ? 'text-orange-400' : 'text-sky-400'}`}>
                                            {seg.dH_kjmol >= 0 ? '+' : ''}{fmtNum(seg.dH_kjmol * dfc)}
                                          </td>
                                        </tr>
                                      ))
                                  )}
                                </tbody>
                              </table>
                            </div>
                          );
                        };

                        const rSegs = result.sensibleSegments.filter(sp => sp.side === 'reactant');
                        const pSegs = result.sensibleSegments.filter(sp => sp.side === 'product');
                        const hasExcess = result.excessSegments.length > 0 && result.excessSegments.some(sp => sp.segments.some(seg => seg.type === 'latent' || Math.abs(seg.T2_K - seg.T1_K) > 0.5));

                        return (
                          <>
                            {/* Step 1 */}
                            <div>
                              <p className="text-xs font-semibold text-[#EE6666] mb-2 flex items-center gap-1.5">
                                <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-[#EE6666]/20 text-[10px] font-bold shrink-0">1</span>
                                Reactants sensible heat ({fmtT(result.T1_K!)} &rarr; {fmtT(298.15)})
                              </p>
                              {renderSegTable(
                                rSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-[#EE6666]', segments: sp.segments, missingData: sp.missingData })),
                                `T\u1d35\u208f = ${fmtT(298.15)} \u2014 no sensible heat needed for reactants.`
                              )}
                            </div>

                            {/* Step 2 */}
                            <div>
                              <p className="text-xs font-semibold text-amber-400 mb-2 flex items-center gap-1.5">
                                <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-amber-400/20 text-[10px] font-bold shrink-0">2</span>
                                &Delta;H&deg;<sub>rxn</sub> at {fmtT(298.15)} &mdash; Hess&apos;s Law
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
                                      <th className="px-2 py-1 text-center font-medium">&nu;</th>
                                      <th className="px-2 py-1 text-right font-medium">&Delta;H ({conv.label})</th>
                                    </tr>
                                  </thead>
                                  <tbody>
                                    {result.contributions.map((c, i) => (
                                      <tr key={i} className="border-b last:border-0 hover:bg-muted/20">
                                        <td className={`px-2 py-1.5 font-mono whitespace-nowrap overflow-hidden text-ellipsis ${c.side === 'reactant' ? 'text-[#EE6666]' : 'text-sky-400'}`}>{c.formula}</td>
                                        <td className="px-2 py-1.5 tabular-nums text-muted-foreground whitespace-nowrap overflow-hidden text-ellipsis">{fmtNum(c.hof_kjmol * df)}</td>
                                        <td className="px-2 py-1.5 text-center tabular-nums text-muted-foreground">{c.side === 'reactant' ? '\u2212' : '+'}{c.coefficient}</td>
                                        <td className={`px-2 py-1.5 text-right tabular-nums font-medium whitespace-nowrap ${c.contribution_kjmol >= 0 ? 'text-sky-400' : 'text-red-400'}`}>
                                          {c.contribution_kjmol >= 0 ? '+' : ''}{fmtNum(c.contribution_kjmol * dfc)}
                                        </td>
                                      </tr>
                                    ))}
                                    <tr className="border-t-2 bg-muted/20">
                                      <td className="px-2 py-1.5 font-bold" colSpan={3}>&Delta;H&deg;<sub>rxn</sub> ({fmtT(298.15)})</td>
                                      <td className={`px-2 py-1.5 text-right tabular-nums font-bold whitespace-nowrap ${result.dH_rxn_kjmol >= 0 ? 'text-sky-400' : 'text-red-400'}`}>
                                        {result.dH_rxn_kjmol >= 0 ? '+' : ''}{fmtNum(result.dH_rxn_kjmol * dfc)}
                                      </td>
                                    </tr>
                                  </tbody>
                                </table>
                              </div>
                            </div>

                            {/* Step 3 */}
                            <div>
                              <p className="text-xs font-semibold text-sky-400 mb-2 flex items-center gap-1.5">
                                <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-sky-400/20 text-[10px] font-bold shrink-0">3</span>
                                Products sensible + latent heat ({fmtT(298.15)} &rarr; {fmtT(result.T2_K!)})
                              </p>
                              {renderSegTable(
                                pSegs.map(sp => ({ coeff: sp.coefficient, formula: sp.formula, colorClass: 'text-sky-400', segments: sp.segments, missingData: sp.missingData })),
                                `T\u1d37\u1d57 = ${fmtT(298.15)} \u2014 no sensible heat needed for products.`
                              )}
                            </div>

                            {/* Excess reactants (if conversion < 100%) */}
                            {hasExcess && (
                              <div>
                                <p className="text-xs font-semibold text-amber-500 mb-2 flex items-center gap-1.5">
                                  <span className="inline-flex items-center justify-center w-4 h-4 rounded-full bg-amber-500/20 text-[10px] font-bold shrink-0">+</span>
                                  Unreacted excess ({(100 - conversion).toFixed(0)}%) &mdash; {fmtT(298.15)} &rarr; {fmtT(result.T2_K!)}
                                </p>
                                {renderSegTable(
                                  result.excessSegments.map(sp => ({ coeff: +sp.excessCoefficient.toFixed(4), formula: sp.formula, colorClass: 'text-amber-500', segments: sp.segments, missingData: sp.missingData })),
                                  'No excess heat needed.'
                                )}
                              </div>
                            )}

                            {/* Total */}
                            <div className="border-t-2 pt-3 flex items-center justify-between">
                              <span className="text-xs font-bold">&Delta;H&deg;<sub>tot</sub> ({fmtT(result.T1_K!)} &rarr; {fmtT(result.T2_K!)})</span>
                              <span className={`tabular-nums font-bold text-sm whitespace-nowrap ${result.dH_total_kjmol >= 0 ? 'text-orange-400' : 'text-sky-400'}`}>
                                {result.dH_total_kjmol >= 0 ? '+' : ''}{fmtNum(result.dH_total_kjmol * dfc)} {conv.label}
                              </span>
                            </div>
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
                              <th className="px-3 py-2 text-center font-semibold">&nu;</th>
                              <th className="px-3 py-2 text-center font-semibold whitespace-nowrap">&Delta;H&deg;<sub>f</sub> ({conv.label})</th>
                              <th className="px-4 py-2 text-center font-semibold whitespace-nowrap">Contribution ({conv.label})</th>
                            </tr>
                          </thead>
                          <tbody>
                            {result.contributions.map((c, i) => (
                              <tr key={i} className="border-b last:border-0 hover:bg-muted/30 transition-colors">
                                <td className="px-4 py-2 text-xs max-w-[140px] truncate text-center" title={c.name}>{c.name}</td>
                                <td className="px-3 py-2 font-mono text-xs text-center">{c.formula}</td>
                                <td className="px-3 py-2 text-center text-xs tabular-nums">{c.side === 'reactant' ? '\u2212' : '+'}{c.coefficient}</td>
                                <td className="px-3 py-2 text-center tabular-nums text-xs">{fmtNum(c.hof_kjmol * df)}</td>
                                <td className={`px-4 py-2 text-center tabular-nums text-xs font-medium ${c.contribution_kjmol >= 0 ? 'text-sky-400' : 'text-red-400'}`}>
                                  {c.contribution_kjmol >= 0 ? '+' : ''}{fmtNum(c.contribution_kjmol * dfc)}
                                </td>
                              </tr>
                            ))}
                            <tr className="border-t-2 bg-muted/20">
                              <td className="px-4 py-3 font-bold text-sm text-center" colSpan={3}>&Delta;H&deg;<sub>rxn</sub> (298 K)</td>
                              <td className="px-3 py-2 text-center text-xs font-semibold">Total:</td>
                              <td className={`px-4 py-3 text-center tabular-nums font-bold text-sm ${result.dH_rxn_kjmol >= 0 ? 'text-sky-400' : 'text-red-400'}`}>
                                {result.dH_rxn_kjmol >= 0 ? '+' : ''}{fmtNum(result.dH_rxn_kjmol * dfc)}
                              </td>
                            </tr>
                          </tbody>
                        </table>
                      </div>
                    </CardContent>
                  </Card>
                )}

                <Card>
                  <CardHeader className="pb-1"><CardTitle className="text-sm">Enthalpy Breakdown (Hess&apos;s Law)</CardTitle></CardHeader>
                  <CardContent><ReactECharts option={chartOption} style={{ height: 280 }} /></CardContent>
                </Card>
                {result.hasSensibleHeat && sensibleChartOption && (
                  <Card>
                    <CardHeader className="pb-1"><CardTitle className="text-sm">Sensible &amp; Latent Heat Pathway</CardTitle></CardHeader>
                    <CardContent><ReactECharts option={sensibleChartOption} style={{ height: 240 }} /></CardContent>
                  </Card>
                )}
              </div>
            ) : (
              <div className="flex h-full min-h-[400px] items-center justify-center">
                <p className="text-sm text-muted-foreground animate-pulse">Calculating&#x2026;</p>
              </div>
            )}
          </div>

        </div>
      </div>
    </TooltipProvider>
  );
}

