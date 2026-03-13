'use client';

import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { supabase } from '@/lib/supabaseClient';
import { type SupabaseClient as _SupabaseClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";

// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
// Import specific types from the main echarts package
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart, ScatterChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

// Register necessary ECharts components
echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent,
  ToolboxComponent,
  LineChart,
  ScatterChart,
  CanvasRenderer
]);

import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip as ShadTooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { ArrowLeftRight, Download } from 'lucide-react'; // Import swap icon

// Import VLE calculation logic and types
import {
  calculatePsat_Pa as libCalculatePsat_Pa,
  // UNIFAC
  calculateBubbleTemperatureUnifac,
  calculateBubblePressureUnifac,
  fetchUnifacInteractionParams,
  type UnifacParameters,
  // NRTL
  fetchNrtlParameters,
  calculateBubbleTemperatureNrtl,
  calculateBubblePressureNrtl,
  type NrtlInteractionParams,
  // Peng–Robinson
  fetchPrInteractionParams,
  calculateBubbleTemperaturePr,
  calculateBubblePressurePr,
  type PrInteractionParams,
  // SRK
  fetchSrkInteractionParams,
  calculateBubbleTemperatureSrk,
  calculateBubblePressureSrk,
  type SrkInteractionParams,
  // UNIQUAC
  fetchUniquacInteractionParams,
  calculateBubbleTemperatureUniquac,
  calculateBubblePressureUniquac,
  type UniquacInteractionParams as LibUniquacInteractionParams,
  // Wilson
  fetchWilsonInteractionParams,
  calculateBubbleTemperatureWilson,
  calculateBubblePressureWilson,
  type WilsonInteractionParams as LibWilsonInteractionParams,
  antoineBoilingPointSolverLocal
} from '@/lib/vle-calculations';

// Import Shared VLE Types
import type {
    CompoundData,
    BubbleDewResult,
} from '@/lib/vle-types';
import { fetchCompoundDataFromHysys, fetchCompoundSuggestions, resolveSimName, formatCompoundName } from '@/lib/antoine-utils';
import type { CompoundAlias } from '@/lib/antoine-utils';


type EChartsPoint = [number, number] | [number | null, number | null];
type FluidPackageTypeMcCabe = 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac' | 'wilson';
type SensitivityVar = 'tc_or_p' | 'xd' | 'xf' | 'xb' | 'q' | 'r' | 'murphree';

// Round a slider minimum UP to the nearest step increment so the min tick aligns
// with the rest of the slider's grid.  Falls back to 3-sig-fig rounding when
// no step is provided (e.g. when used outside slider context).
function niceMin(v: number, step?: number): number {
    if (!isFinite(v)) return v;
    if (step && step > 0) {
        return Math.ceil(v / step) * step;
    }
    if (v === 0) return 0;
    if (v < 0) {
        const absV = Math.abs(v);
        const mag = Math.pow(10, Math.floor(Math.log10(absV)) - 2);
        return -Math.floor(absV / mag) * mag;
    }
    const mag = Math.pow(10, Math.floor(Math.log10(v)) - 2);
    return Math.ceil(v / mag) * mag;
}

// Pure stage-counting function used by sensitivity sweep (no state side-effects).
function countStagesForSweep(
    xValues: number[], yValues: number[],
    params: { xd: number; xb: number; xf: number; q: number; r: number; murphreeEfficiency: number; buffer: number }
): { stages: number; feedStage: number } {
    const { xd, xb, xf, q, r, murphreeEfficiency, buffer } = params;
    const rectifyingSlope = r / (r + 1);
    const rectifyingIntercept = xd / (r + 1);
    const feedSlope = Math.abs(q - 1) < 1e-6 ? Infinity : q / (q - 1);
    const feedIntercept = feedSlope === Infinity ? NaN : xf - feedSlope * xf;
    let xIntersect: number, yIntersect: number;
    if (feedSlope === Infinity) {
        xIntersect = xf;
        yIntersect = rectifyingSlope * xf + rectifyingIntercept;
    } else {
        xIntersect = (feedIntercept - rectifyingIntercept) / (rectifyingSlope - feedSlope);
        yIntersect = rectifyingSlope * xIntersect + rectifyingIntercept;
    }
    xIntersect = Math.max(0, Math.min(1, xIntersect));
    yIntersect = Math.max(0, Math.min(1, yIntersect));
    const strippingSlope = (xIntersect === xb) ? Infinity : (yIntersect - xb) / (xIntersect - xb);
    const strippingIntercept = strippingSlope === Infinity ? NaN : xb - strippingSlope * xb;
    let yDataForStepping = yValues;
    if (murphreeEfficiency < 1.0) {
        yDataForStepping = xValues.map((x, i) => {
            const y_eq = yValues[i];
            const y_op = x >= xIntersect
                ? rectifyingSlope * x + rectifyingIntercept
                : (strippingSlope === Infinity ? xb : strippingSlope * x + strippingIntercept);
            return y_op + murphreeEfficiency * (y_eq - y_op);
        });
    }
    let stageCount = 0, feedStageCount = 0;
    let currentX = xd, currentY = xd;
    let previousSectionIsRectifying = true;
    const MAX_STAGES = 300;
    for (let i = 0; i < MAX_STAGES; i++) {
        if (currentX <= xb + buffer) break;
        let intersectX = NaN;
        for (let j = 0; j < xValues.length - 1; j++) {
            if ((yDataForStepping[j] <= currentY && yDataForStepping[j + 1] >= currentY) ||
                (yDataForStepping[j] >= currentY && yDataForStepping[j + 1] <= currentY)) {
                intersectX = Math.abs(yDataForStepping[j + 1] - yDataForStepping[j]) < 1e-9
                    ? (yDataForStepping[j] === currentY ? xValues[j] : xValues[j + 1])
                    : xValues[j] + ((currentY - yDataForStepping[j]) / (yDataForStepping[j + 1] - yDataForStepping[j])) * (xValues[j + 1] - xValues[j]);
                break;
            }
        }
        if (isNaN(intersectX)) break;
        // Detect near-pinch / azeotrope — stepping is making no progress toward xb
        if (i > 0 && Math.abs(intersectX - currentX) < 1e-5) break;
        stageCount++;
        if (intersectX <= xb + buffer) break;
        const currentSectionIsRectifying = intersectX >= xIntersect;
        let nextY = currentSectionIsRectifying
            ? rectifyingSlope * intersectX + rectifyingIntercept
            : (strippingSlope === Infinity ? xb : strippingSlope * intersectX + strippingIntercept);
        nextY = Math.max(0, Math.min(1, nextY));
        if (feedStageCount === 0 && currentSectionIsRectifying !== previousSectionIsRectifying) feedStageCount = stageCount;
        previousSectionIsRectifying = currentSectionIsRectifying;
        currentX = intersectX;
        currentY = nextY;
    }
    // If we exhausted MAX_STAGES without converging, treat as infeasible (null)
    const converged = currentX <= xb + buffer * 10 || stageCount < MAX_STAGES;
    if (!converged) return { stages: 0, feedStage: 0 };
    return { stages: stageCount, feedStage: (feedStageCount === 0 && stageCount > 0) ? stageCount : feedStageCount };
}

// Pure synchronous VLE computation used by T/P sensitivity sweep.
function computeEquilibriumPointsSync(
    data1: CompoundData, data2: CompoundData,
    fluidPackage: FluidPackageTypeMcCabe,
    activityParams: any,
    useTemperature: boolean,
    value: number, // tempC when useTemperature, pressureBar when !useTemperature
    N: number = 80
): { x: number[]; y: number[] } | null {
    try {
        const components: CompoundData[] = [data1, data2];
        const stepSize = 1 / (N - 1);
        const xFeedValues = Array.from({ length: N }, (_, i) => parseFloat((i * stepSize).toFixed(4)));
        const xOut: number[] = [0.0], yOut: number[] = [0.0];
        const pressurePa = !useTemperature ? value * 1e5 : null;
        const Tbp1 = (!useTemperature && data1.antoine) ? antoineBoilingPointSolverLocal(data1.antoine, pressurePa!) : null;
        const Tbp2 = (!useTemperature && data2.antoine) ? antoineBoilingPointSolverLocal(data2.antoine, pressurePa!) : null;
        let prevBubbleT: number | null = null;
        for (const x1 of xFeedValues.slice(1, -1)) {
            let result: any = null;
            if (useTemperature) {
                const Tk = value + 273.15;
                const pGuess = (data1.antoine && data2.antoine)
                    ? (libCalculatePsat_Pa(data1.antoine, Tk) + libCalculatePsat_Pa(data2.antoine, Tk)) / 2
                    : 101325;
                if (fluidPackage === 'unifac') result = calculateBubblePressureUnifac(components, x1, Tk, activityParams as UnifacParameters);
                else if (fluidPackage === 'nrtl') result = calculateBubblePressureNrtl(components, x1, Tk, activityParams as NrtlInteractionParams);
                else if (fluidPackage === 'pr') result = calculateBubblePressurePr(components, x1, Tk, activityParams as PrInteractionParams, pGuess);
                else if (fluidPackage === 'srk') result = calculateBubblePressureSrk(components, x1, Tk, activityParams as SrkInteractionParams, pGuess);
                else if (fluidPackage === 'uniquac') result = calculateBubblePressureUniquac(components, x1, Tk, activityParams as LibUniquacInteractionParams);
                else if (fluidPackage === 'wilson') result = calculateBubblePressureWilson(components, x1, Tk, activityParams as LibWilsonInteractionParams);
            } else {
                let tGuessDefault = (Tbp1 !== null && Tbp2 !== null) ? x1 * Tbp1 + (1 - x1) * Tbp2
                    : (Tbp1 ?? Tbp2 ?? 373.15);
                const tGuess = Math.max(150, Math.min(prevBubbleT ?? tGuessDefault, 1000));
                if (fluidPackage === 'unifac') result = calculateBubbleTemperatureUnifac(components, x1, pressurePa!, activityParams as UnifacParameters, tGuess);
                else if (fluidPackage === 'nrtl') result = calculateBubbleTemperatureNrtl(components, x1, pressurePa!, activityParams as NrtlInteractionParams, tGuess);
                else if (fluidPackage === 'pr') result = calculateBubbleTemperaturePr(components, x1, pressurePa!, activityParams as PrInteractionParams, tGuess);
                else if (fluidPackage === 'srk') result = calculateBubbleTemperatureSrk(components, x1, pressurePa!, activityParams as SrkInteractionParams, tGuess);
                else if (fluidPackage === 'uniquac') result = calculateBubbleTemperatureUniquac(components, x1, pressurePa!, activityParams as LibUniquacInteractionParams, tGuess);
                else if (fluidPackage === 'wilson') result = calculateBubbleTemperatureWilson(components, x1, pressurePa!, activityParams as LibWilsonInteractionParams, tGuess);
                if (result?.T_K) prevBubbleT = result.T_K;
            }
            if (result && result.error === undefined && typeof result.comp1_equilibrium === 'number' && !isNaN(result.comp1_equilibrium)) {
                xOut.push(x1);
                yOut.push(result.comp1_equilibrium);
            }
        }
        xOut.push(1.0); yOut.push(1.0);
        return { x: xOut, y: yOut };
    } catch { return null; }
}

export default function McCabeThielePage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  
  // Input States - Updated defaults for initial load
  const [comp1Name, setComp1Name] = useState('Methanol');
  const [comp2Name, setComp2Name] = useState('H2O');
  const [temperatureC, setTemperatureC] = useState<number | null>(60);
  const [pressureBar, setPressureBar] = useState<number | null>(1); // Default, but will be overridden by useTemperature=true
  
  // New string states for input fields
  const [_temperatureInput, _setTemperatureInput] = useState<string>(temperatureC !== null ? String(temperatureC) : '');
  const [_pressureInput, _setPressureInput] = useState<string>(pressureBar !== null ? String(pressureBar) : '');
  const [tempMax, setTempMax] = useState<string>('100');
  const [pressureMax, setPressureMax] = useState<string>('10');
  const [qMinSlider, setQMinSlider] = useState<string>('-1');
  const [qMaxSlider, setQMaxSlider] = useState<string>('2');
  const [rMinSlider, setRMinSlider] = useState<string>('0.1');
  const [rMaxSlider, setRMaxSlider] = useState<string>('10');

  const [useTemperature, setUseTemperature] = useState(true);
  const [fluidPackage, setFluidPackage] = useState<FluidPackageTypeMcCabe>('uniquac');

  // Parameter States
  const [xd, setXd] = useState(0.9);
  const [xb, setXb] = useState(0.1);
  const [xf, setXf] = useState(0.5);
  const [q, setQ] = useState(1.0);
  const [r, setR] = useState(1.5);

  const buffer = 1e-4;

  // Data & Control States
  const [equilibriumData, setEquilibriumData] = useState<{ x: number[], y: number[] } | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // --- NEW: Caching states for component data and parameters ---
  const [comp1Data, setComp1Data] = useState<CompoundData | null>(null);
  const [comp2Data, setComp2Data] = useState<CompoundData | null>(null);
  const [interactionParams, setInteractionParams] = useState<any>(null);
  // Track the exact input names+package that produced the current cached data.
  // Prevents API refetches on every temperature/pressure slider move.
  const lastFetchKeyRef = useRef<string | null>(null);

  // Result States
  const [stages, setStages] = useState<number | null>(null);
  const [feedStage, setFeedStage] = useState<number | null>(null);
  const [rMin, setRMin] = useState<number | null>(null);
  const [murphreeEfficiency, setMurphreeEfficiency] = useState(1.0);

  // State for ECharts options - Use the imported EChartsOption type
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null); // Ref for ECharts instance
  // State for displayed parameters in the title
  const [displayedComp1, setDisplayedComp1] = useState('');
  const [displayedComp2, setDisplayedComp2] = useState('');
  const [displayedTemp, setDisplayedTemp] = useState<number | null>(null);
  const [displayedPressure, setDisplayedPressure] = useState<number | null>(null);
  const [displayedUseTemp, setDisplayedUseTemp] = useState(true);
  const [displayedFluidPackage, setDisplayedFluidPackage] = useState<FluidPackageTypeMcCabe>('unifac');

  // Suggestion states
  const [comp1Suggestions, setComp1Suggestions] = useState<CompoundAlias[]>([]);
  const [comp2Suggestions, setComp2Suggestions] = useState<CompoundAlias[]>([]);
  const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
  const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
  const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);

  // Display names (human-readable FullName from HYSYS ALIASES)
  const [comp1DisplayName, setComp1DisplayName] = useState('Methanol');
  const [comp2DisplayName, setComp2DisplayName] = useState('Water');

  const input1Ref = useRef<HTMLInputElement>(null);
  const input2Ref = useRef<HTMLInputElement>(null);
  const activeComponentRef = useRef<HTMLDivElement>(null);

  // Flag to trigger automatic graph regeneration after UI actions like swap, toggle, or fluid-package change
  const [autoGeneratePending, setAutoGeneratePending] = useState(false);
  const [autoGenerateOnCompChange, setAutoGenerateOnCompChange] = useState(true);

  // Sensitivity analysis state
  const [sensitivityMode, setSensitivityMode] = useState(false);
  const [activeSensVar, setActiveSensVar] = useState<SensitivityVar | null>(null);
  const [sensitivityRunning, setSensitivityRunning] = useState(false);
  const [sensitivityEchartsOptions, setSensitivityEchartsOptions] = useState<EChartsOption | null>(null);
  const sensitivityDebounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  // Derive valid pressure/temperature bounds from Antoine Tmin/Tmax AND critical properties.
  // Same logic as binary-phase-diagrams so sliders can't reach unphysical conditions.
  const antoineLimits = useMemo(() => {
      if (!comp1Data?.antoine || !comp2Data?.antoine) return null;
      const a1 = comp1Data.antoine, a2 = comp2Data.antoine;
      const Tc1 = comp1Data.Tc_K ?? comp1Data.prParams?.Tc_K ?? null;
      const Tc2 = comp2Data.Tc_K ?? comp2Data.prParams?.Tc_K ?? null;
      const Pc1 = comp1Data.Pc_Pa ?? comp1Data.prParams?.Pc_Pa ?? null;
      const Pc2 = comp2Data.Pc_Pa ?? comp2Data.prParams?.Pc_Pa ?? null;
      const antoineTmaxK = Math.min(a1.Tmax_K, a2.Tmax_K);
      const critTmaxK    = (Tc1 != null && Tc2 != null) ? Math.min(Tc1, Tc2) : (Tc1 ?? Tc2 ?? Infinity);
      const maxTK = Math.min(antoineTmaxK, critTmaxK);
      const pMaxPa_antoine = Math.min(
          libCalculatePsat_Pa(a1, a1.Tmax_K), libCalculatePsat_Pa(a2, a2.Tmax_K)
      );
      const critPmaxPa = (Pc1 != null && Pc2 != null) ? Math.min(Pc1, Pc2) : (Pc1 ?? Pc2 ?? Infinity);
      const maxPPa = Math.min(
          isFinite(pMaxPa_antoine) && pMaxPa_antoine > 0 ? pMaxPa_antoine : Infinity, critPmaxPa
      );
      const pMinPa = Math.max(libCalculatePsat_Pa(a1, a1.Tmin_K), libCalculatePsat_Pa(a2, a2.Tmin_K));
      return {
          maxPBar: isFinite(maxPPa)  && maxPPa  > 0 ? maxPPa  / 1e5 : null,
          minPBar: isFinite(pMinPa)  && pMinPa  > 0 ? Math.max(0.1, pMinPa / 1e5) : 0.1,
          maxTC:   isFinite(maxTK) ? maxTK - 273.15 : null,
          minTC:   isFinite(Math.max(a1.Tmin_K, a2.Tmin_K)) ? Math.max(a1.Tmin_K, a2.Tmin_K) - 273.15 : null,
      };
  }, [comp1Data, comp2Data]);

  // When compound limits change, clamp the Max-input boxes and current slider values.
  useEffect(() => {
      if (!antoineLimits) return;
      // maxPBar clamping removed — users may explore supercritical/high-pressure conditions.
      if (antoineLimits.maxTC !== null) {
          const cur = parseFloat(tempMax);
          if (isFinite(cur) && cur > antoineLimits.maxTC)
              setTempMax(String(Math.floor(antoineLimits.maxTC)));
      }
      if (antoineLimits.minPBar !== null)
          setPressureBar(prev => (prev !== null && prev < antoineLimits.minPBar!) ? antoineLimits.minPBar! : prev);
      if (antoineLimits.minTC !== null)
          setTemperatureC(prev => (prev !== null && prev < antoineLimits.minTC!) ? antoineLimits.minTC! : prev);
  }, [antoineLimits]); // eslint-disable-line react-hooks/exhaustive-deps

  // Fugacity-deviation diagnostic using the Pitzer-Tsonopoulos second-virial correlation
  // (Tsonopoulos 1974, AIChE J. 20:263). Computes |φ−1| per compound at its estimated
  // boiling-point temperature under the system pressure, then reports the worst case.
  // Activity-coefficient models ignore this correction entirely (φᵢ = 1 assumed).
  const modelPressureWarning = useMemo(() => {
      if (useTemperature) return null;
      const activityModels: FluidPackageTypeMcCabe[] = ['nrtl', 'wilson', 'uniquac', 'unifac'];
      if (!activityModels.includes(fluidPackage)) return null;
      if (!pressureBar || !comp1Data || !comp2Data) return null;
      const P_Pa = pressureBar * 1e5;
      const pitzDev = (comp: typeof comp1Data): number | null => {
          const Tc = comp.Tc_K ?? comp.prParams?.Tc_K ?? null;
          const Pc = comp.Pc_Pa ?? comp.prParams?.Pc_Pa ?? null;
          const omega = comp.prParams?.omega ?? comp.srkParams?.omega ?? null;
          if (!Tc || !Pc || omega == null || Tc <= 0 || Pc <= 0) return null;
          const Tbp = (comp.antoine ? antoineBoilingPointSolverLocal(comp.antoine, P_Pa) : null) ?? (Tc * 0.65);
          const Tr = Tbp / Tc;
          if (Tr <= 0 || !isFinite(Tr)) return null;
          const B0 = 0.083 - 0.422 / Math.pow(Tr, 1.6);
          const B1 = 0.139 - 0.172 / Math.pow(Tr, 4.2);
          const B_m3mol = (B0 + omega * B1) * (8.314 * Tc / Pc);
          const lnPhi = B_m3mol * P_Pa / (8.314 * Tbp);
          return Math.abs(Math.exp(lnPhi) - 1);
      };
      const maxDev = Math.max(pitzDev(comp1Data) ?? 0, pitzDev(comp2Data) ?? 0);
      if (maxDev === 0) return null;
      return null; // Force-disable all UI warnings for fugacity.
      // const pct = maxDev * 100;
      // if (pct > 15) return { level: 'error' as const, pct };
      // if (pct > 5)  return { level: 'warn'  as const, pct };
      // return null;
  }, [useTemperature, fluidPackage, pressureBar, comp1Data, comp2Data]);

  // --- useEffect for initial graph load ---
  useEffect(() => {
    // Automatically calculate and display the graph on initial component mount
    // with the default states (Methanol/Water, 60C, UNIQUAC).
    if (supabase && autoGenerateOnCompChange) { // Ensure supabase is initialized before first calculation
        calculateEquilibriumCurve();
        setAutoGenerateOnCompChange(false);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase]); // Add supabase as a dependency if it can be initially undefined

  useEffect(() => {
    if (autoGeneratePending) {
        calculateEquilibriumCurve();
        setAutoGeneratePending(false);
    }
  }, [autoGeneratePending]); // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    if (autoGenerateOnCompChange) return;
    calculateEquilibriumCurve();
  }, [temperatureC, pressureBar, fluidPackage]); // eslint-disable-line react-hooks/exhaustive-deps

  // --- Data Fetching (Updated to use HYSYS tables) ---
  async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
    if (!supabase) { throw new Error("Supabase client not initialized."); }
    try {
        const result = await fetchCompoundDataFromHysys(supabase, compoundName);
        return result;
    } catch (err: any) {
        console.error(`Error fetching data for ${compoundName} in McCabe:`, err.message);
        setError(`Data fetch failed for ${compoundName}: ${err.message}`);
        return null;
    }
  }

  // --- Fetch Suggestions (from HYSYS ALIASES table) ---
  const fetchSuggestions = useCallback(async (inputValue: string, inputTarget: 'comp1' | 'comp2') => {
    if (!inputValue || inputValue.length < 2 || !supabase) {
      if (inputTarget === 'comp1') {
        setComp1Suggestions([]);
        setShowComp1Suggestions(false);
      } else {
        setComp2Suggestions([]);
        setShowComp2Suggestions(false);
      }
      return;
    }

    try {
      const suggestions = await fetchCompoundSuggestions(supabase, inputValue);
      if (inputTarget === 'comp1') {
        setComp1Suggestions(suggestions);
        setShowComp1Suggestions(suggestions.length > 0);
      } else {
        setComp2Suggestions(suggestions);
        setShowComp2Suggestions(suggestions.length > 0);
      }
    } catch (err) {
      console.error("Error fetching suggestions:", err);
      if (inputTarget === 'comp1') setComp1Suggestions([]); else setComp2Suggestions([]);
    }
  }, [supabase]);

  const handleComp1NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp1Name(newValue);
    setComp1DisplayName(newValue);
    setActiveSuggestionInput('comp1');
    if (newValue.trim() === "") {
      setShowComp1Suggestions(false);
      setComp1Suggestions([]);
    } else {
      fetchSuggestions(newValue, 'comp1');
    }
  };

  const handleComp2NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp2Name(newValue);
    setComp2DisplayName(newValue);
    setActiveSuggestionInput('comp2');
    if (newValue.trim() === "") {
      setShowComp2Suggestions(false);
      setComp2Suggestions([]);
    } else {
      fetchSuggestions(newValue, 'comp2');
    }
  };

  const handleSuggestionClick = (alias: CompoundAlias, inputTarget: 'comp1' | 'comp2') => {
    if (inputTarget === 'comp1') {
      setComp1Name(alias.simName);
      setComp1DisplayName(formatCompoundName(alias.fullName));
      setShowComp1Suggestions(false);
      setComp1Suggestions([]);
    } else {
      setComp2Name(alias.simName);
      setComp2DisplayName(formatCompoundName(alias.fullName));
      setShowComp2Suggestions(false);
      setComp2Suggestions([]);
    }
    // setActiveSuggestionInput(null); // Don't reset active input here, let blur handle it or next focus
  };
  
  // Effect to handle clicks outside of the suggestion boxes
  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      const target = event.target as Node;
      
      // If there's an active component with suggestions showing, check if click is outside
      if (activeComponentRef.current && !activeComponentRef.current.contains(target)) {
        setShowComp1Suggestions(false);
        setShowComp2Suggestions(false);
      }
    }
    
    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [activeSuggestionInput]);

  const calculateEquilibriumCurve = useCallback(async (showLoading: boolean = true, pointsCount: number = 101) => {
    if (showLoading) {
      setLoading(true);
      // setEquilibriumData(null); // This line is removed to keep the old graph visible.
    }
    setError(null);

    if (!comp1Name || !comp2Name) {
        setError("Please enter valid compound names.");
        setLoading(false);
        return;
    }
    if (!supabase) {
        setError("Supabase client not available.");
        setLoading(false);
        return;
    }

    // Validate temperature or pressure based on mode
    if (useTemperature && temperatureC === null) {
        setError("Temperature must be provided to calculate.");
        setLoading(false);
        return;
    }
    if (!useTemperature && pressureBar === null) {
        setError("Pressure must be provided to calculate.");
        setLoading(false);
        return;
    }

    try {
        let data1 = comp1Data;
        let data2 = comp2Data;
        let activityParameters = interactionParams;

        // Check if we need to re-fetch data
        const currentFetchKey = `${comp1Name}|${comp2Name}|${fluidPackage}`;
        const needsFetching = !data1 || !data2 || !activityParameters || lastFetchKeyRef.current !== currentFetchKey;

        if (needsFetching) {
            // Resolve display names to SimNames (handles case where user typed a FullName without selecting a suggestion)
            const [sim1, sim2] = await Promise.all([resolveSimName(supabase, comp1Name), resolveSimName(supabase, comp2Name)]);
            const resolvedComp1 = sim1 || comp1Name;
            const resolvedComp2 = sim2 || comp2Name;
            const [fetchedData1, fetchedData2] = await Promise.all([fetchCompoundDataLocal(resolvedComp1), fetchCompoundDataLocal(resolvedComp2)]);
            if (!fetchedData1 || !fetchedData2) {
                setLoading(false);
                return;
            }
            
            const fetchedComponents: CompoundData[] = [fetchedData1, fetchedData2];
            let fetchedActivityParameters: any;

            if (fluidPackage === 'unifac') {
                if (!fetchedData1.unifacGroups || !fetchedData2.unifacGroups) throw new Error("UNIFAC groups missing.");
                const allSubgroupIds = new Set<number>();
                fetchedComponents.forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                fetchedActivityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
            } else if (fluidPackage === 'nrtl') {
                fetchedActivityParameters = await fetchNrtlParameters(supabase, fetchedData1.name, fetchedData2.name);
                if ((fetchedActivityParameters as any)?._usedUnifacFallback) {
                    setError(`Warning: NRTL parameters not found in database for ${fetchedData1.name}–${fetchedData2.name}. Using UNIFAC-based estimate.`);
                }
            } else if (fluidPackage === 'pr') {
                if (!fetchedData1.prParams || !fetchedData2.prParams) throw new Error("PR pure component params missing.");
                fetchedActivityParameters = await fetchPrInteractionParams(supabase, fetchedData1.name, fetchedData2.name);
            } else if (fluidPackage === 'srk') {
                if (!fetchedData1.srkParams || !fetchedData2.srkParams) throw new Error("SRK pure component params missing.");
                fetchedActivityParameters = await fetchSrkInteractionParams(supabase, fetchedData1.name, fetchedData2.name);
            } else if (fluidPackage === 'uniquac') {
                if (!fetchedData1.uniquacParams || !fetchedData2.uniquacParams) throw new Error("UNIQUAC pure component params missing.");
                fetchedActivityParameters = await fetchUniquacInteractionParams(supabase, fetchedData1.name, fetchedData2.name);
                if ((fetchedActivityParameters as any)?._usedUnifacFallback) {
                    setError(`Warning: UNIQUAC parameters not found in database for ${fetchedData1.name}–${fetchedData2.name}. Using UNIFAC-based estimate.`);
                }
            } else if (fluidPackage === 'wilson') {
                if (!fetchedData1.wilsonParams || !fetchedData2.wilsonParams) throw new Error("Wilson pure component params missing.");
                fetchedActivityParameters = await fetchWilsonInteractionParams(supabase, fetchedData1.name, fetchedData2.name);
            } else {
                throw new Error(`Unsupported fluid package: ${fluidPackage}`);
            }
            
            // Update cache and use the new data for this run
            setComp1Data(fetchedData1);
            setComp2Data(fetchedData2);
            setInteractionParams(fetchedActivityParameters);
            lastFetchKeyRef.current = currentFetchKey;
            data1 = fetchedData1;
            data2 = fetchedData2;
            activityParameters = fetchedActivityParameters;
        }

        if (!data1 || !data2 || !activityParameters) {
             setError("Component data or parameters not available for calculation.");
             if (showLoading) setLoading(false);
             return;
        }

        const components: CompoundData[] = [data1, data2];

        const stepSize = 1 / (pointsCount - 1);
        const x_feed_values = Array.from({ length: pointsCount }, (_, i) => parseFloat((i * stepSize).toFixed(4)));
        const calculated_x_values: number[] = [0.0];
        const calculated_y_values: number[] = [0.0];

        // Pre-compute boiling points for pressure mode (computed once, not per-composition)
        const mcTbp1: number | null = !useTemperature ? antoineBoilingPointSolverLocal(data1.antoine, pressureBar! * 1e5) : null;
        const mcTbp2: number | null = !useTemperature ? antoineBoilingPointSolverLocal(data2.antoine, pressureBar! * 1e5) : null;
        // Sequential continuation: seed each composition's bubble-T solver from the previous result
        let mcPrevBubbleT: number | null = null;

        // Loop over interior points only to avoid issues at pure component conditions
        for (const x1_val of x_feed_values.slice(1, -1)) {
            let resultPoint: BubbleDewResult | null = null;
            if (useTemperature) {
                const currentFixedTempK = temperatureC! + 273.15;
                const initialPressureGuess = (libCalculatePsat_Pa(data1.antoine!, currentFixedTempK) + libCalculatePsat_Pa(data2.antoine!, currentFixedTempK)) / 2 || 101325;
                
                if (fluidPackage === 'unifac') resultPoint = calculateBubblePressureUnifac(components, x1_val, currentFixedTempK, activityParameters as UnifacParameters);
                else if (fluidPackage === 'nrtl') resultPoint = calculateBubblePressureNrtl(components, x1_val, currentFixedTempK, activityParameters as NrtlInteractionParams);
                else if (fluidPackage === 'pr') resultPoint = calculateBubblePressurePr(components, x1_val, currentFixedTempK, activityParameters as PrInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'srk') resultPoint = calculateBubblePressureSrk(components, x1_val, currentFixedTempK, activityParameters as SrkInteractionParams, initialPressureGuess);
                else if (fluidPackage === 'uniquac') resultPoint = calculateBubblePressureUniquac(components, x1_val, currentFixedTempK, activityParameters as LibUniquacInteractionParams);
                else if (fluidPackage === 'wilson') resultPoint = calculateBubblePressureWilson(components, x1_val, currentFixedTempK, activityParameters as LibWilsonInteractionParams);
            } else {
                const currentFixedPressurePa = pressureBar! * 1e5;

                let defaultTGuess: number;
                if (mcTbp1 !== null && mcTbp2 !== null) {
                    defaultTGuess = x1_val * mcTbp1 + (1 - x1_val) * mcTbp2;
                } else if (mcTbp1 !== null) {
                    defaultTGuess = mcTbp1;
                } else if (mcTbp2 !== null) {
                    defaultTGuess = mcTbp2;
                } else {
                    defaultTGuess = 373.15;
                }
                const initialTempGuess = Math.max(150, Math.min(mcPrevBubbleT ?? defaultTGuess, 1000));

                if (fluidPackage === 'unifac') resultPoint = calculateBubbleTemperatureUnifac(components, x1_val, currentFixedPressurePa, activityParameters as UnifacParameters, initialTempGuess);
                else if (fluidPackage === 'nrtl') resultPoint = calculateBubbleTemperatureNrtl(components, x1_val, currentFixedPressurePa, activityParameters as NrtlInteractionParams, initialTempGuess);
                else if (fluidPackage === 'pr') resultPoint = calculateBubbleTemperaturePr(components, x1_val, currentFixedPressurePa, activityParameters as PrInteractionParams, initialTempGuess);
                else if (fluidPackage === 'srk') resultPoint = calculateBubbleTemperatureSrk(components, x1_val, currentFixedPressurePa, activityParameters as SrkInteractionParams, initialTempGuess);
                else if (fluidPackage === 'uniquac') resultPoint = calculateBubbleTemperatureUniquac(components, x1_val, currentFixedPressurePa, activityParameters as LibUniquacInteractionParams, initialTempGuess);
                else if (fluidPackage === 'wilson') resultPoint = calculateBubbleTemperatureWilson(components, x1_val, currentFixedPressurePa, activityParameters as LibWilsonInteractionParams, initialTempGuess);

                if (resultPoint?.T_K) mcPrevBubbleT = resultPoint.T_K;
            }

            if (resultPoint && resultPoint.error === undefined && typeof resultPoint.comp1_equilibrium === 'number' && !isNaN(resultPoint.comp1_equilibrium)) {
                calculated_x_values.push(x1_val);
                calculated_y_values.push(resultPoint.comp1_equilibrium);
            } else {
                console.warn(`McCabe: Calculation failed for x1=${x1_val} with ${fluidPackage}. Error: ${resultPoint?.error}`);
            }
        }
        
        calculated_x_values.push(1.0);
        calculated_y_values.push(1.0);

        const sortedPairs = calculated_x_values.map((x_val, i) => ({ x: x_val, y: calculated_y_values[i] }))
            .sort((a, b) => a.x - b.x);

        setEquilibriumData({ x: sortedPairs.map(p => p.x), y: sortedPairs.map(p => p.y) });

        setDisplayedComp1(comp1DisplayName || comp1Name);
        setDisplayedComp2(comp2DisplayName || comp2Name);
        setDisplayedUseTemp(useTemperature);
        setDisplayedFluidPackage(fluidPackage);
        if (useTemperature) { setDisplayedTemp(temperatureC); setDisplayedPressure(null); }
        else { setDisplayedTemp(null); setDisplayedPressure(pressureBar); }

    } catch (err: any) {
        console.error("McCabe: Error calculating equilibrium curve:", err);
        setError(`Calculation failed: ${err.message}`);
        // if (showLoading) setEquilibriumData(null); // This line is removed to keep the graph visible on error.
    } finally {
        if (showLoading) setLoading(false);
    }
  }, [comp1Name, comp2Name, useTemperature, temperatureC, pressureBar, fluidPackage, comp1Data, comp2Data, interactionParams, displayedFluidPackage]);

  // REMOVED: useEffect that automatically called calculateEquilibriumCurve
  // useEffect(() => {
  //   calculateEquilibriumCurve();
  // }, [calculateEquilibriumCurve]); 

  // Initial load with default values can be triggered manually by the user via the button.
  // If an initial automatic load with defaults is strictly required *once* on page load,
  // a separate useEffect with an empty dependency array could call calculateEquilibriumCurve(),
  // but the current request is to stop fetching until the button is pressed.

  const interpolateY = (xTarget: number, xData: number[], yData: number[]): number | null => {
    if (!xData || !yData || xData.length === 0 || xData.length !== yData.length) return null;
    
    // Handle cases where xTarget is outside the data range (clamp to ends)
    if (xTarget <= xData[0]) return yData[0];
    if (xTarget >= xData[xData.length - 1]) return yData[yData.length - 1];
  
    for (let i = 0; i < xData.length - 1; i++) {
      if (xTarget >= xData[i] && xTarget <= xData[i + 1]) {
        const x1 = xData[i], y1 = yData[i];
        const x2 = xData[i + 1], y2 = yData[i + 1];
        if (x2 === x1) return y1; // Avoid division by zero if points are identical
        return y1 + (y2 - y1) * (xTarget - x1) / (x2 - x1);
      }
    }
    return null; // Should ideally be covered by boundary checks
  };

  const formatNumberToPrecision = (num: any, precision: number = 3): string => {
    if (typeof num === 'number') {
      if (num === 0) return '0';
      const fixed = num.toPrecision(precision);
      // Remove trailing zeros after decimal point, but keep integer part if it ends in zero
      if (fixed.includes('.')) {
        return parseFloat(fixed).toString(); 
      }
      return fixed;
    }
    return String(num); // Return as string if not a number
  };


    const generateEChartsOptions = useCallback((xValues: number[], yValues: number[]) => {
    if (!xValues || !yValues || xValues.length === 0) {
        setEchartsOptions({});
        setRMin(null);
        return;
    }
    const series: SeriesOption[] = [];
    
    // --- Define Operating Lines ---
    const rectifyingSlope = r / (r + 1);
    const rectifyingIntercept = xd / (r + 1);
    let feedSlope: number;
    if (Math.abs(q - 1) < 1e-6) { feedSlope = Infinity; } else { feedSlope = q / (q - 1); }
    const feedIntercept = (feedSlope === Infinity) ? NaN : xf - feedSlope * xf;
    let xIntersect: number, yIntersect: number;
    if (feedSlope === Infinity) { xIntersect = xf; yIntersect = rectifyingSlope * xf + rectifyingIntercept; }
    else { xIntersect = (feedIntercept - rectifyingIntercept) / (rectifyingSlope - feedSlope); yIntersect = rectifyingSlope * xIntersect + rectifyingIntercept; }
    xIntersect = Math.max(0, Math.min(1, xIntersect));
    yIntersect = Math.max(0, Math.min(1, yIntersect));
    const strippingSlope = (xIntersect === xb) ? Infinity : (yIntersect - xb) / (xIntersect - xb);
    const strippingIntercept = (strippingSlope === Infinity) ? NaN : xb - strippingSlope * xb;

    // --- Define Data for Calculations ---
    let yDataForStepping = yValues;
    if (murphreeEfficiency < 1.0) {
        yDataForStepping = xValues.map((x, i) => {
            const y_eq = yValues[i];
            let y_op: number;
            if (x >= xIntersect) {
                y_op = rectifyingSlope * x + rectifyingIntercept;
            } else {
                y_op = (strippingSlope === Infinity) ? xb : strippingSlope * x + strippingIntercept;
            }
            return y_op + murphreeEfficiency * (y_eq - y_op);
        });
    }

    // --- PASS 1: Calculate Stage Data and Find Final X Position ---
    let stageCount = 0;
    let feedStageCount = 0;
    let currentX = xd;
    let currentY = xd;
    let finalStageX = xb;
    const stageLineData: EChartsPoint[] = [];
    let previousSectionIsRectifying = true;

    for(let i = 0; i < 300; i++) {
        if (currentX <= xb + buffer) break;

        let intersectX = NaN;
        for (let j = 0; j < xValues.length - 1; j++) {
            if (yDataForStepping[j] === null || yDataForStepping[j+1] === null) continue;
            if ((yDataForStepping[j] <= currentY && yDataForStepping[j+1] >= currentY) || (yDataForStepping[j] >= currentY && yDataForStepping[j+1] <= currentY)) {
                 if (Math.abs(yDataForStepping[j+1] - yDataForStepping[j]) < 1e-9) {
                     intersectX = (yDataForStepping[j] === currentY) ? xValues[j] : xValues[j+1];
                 } else {
                     const fraction = (currentY - yDataForStepping[j]) / (yDataForStepping[j+1] - yDataForStepping[j]);
                     intersectX = xValues[j] + fraction * (xValues[j+1] - xValues[j]);
                 }
                 break;
            }
        }
        if (isNaN(intersectX)) break;
        // Near-pinch / azeotrope detection: stepping has stalled
        if (i > 0 && Math.abs(intersectX - currentX) < 1e-5) break;

        finalStageX = intersectX;
        stageCount++;

        stageLineData.push([currentX, currentY]);
        // Remove the extra horizontal push here since the blocks handle it

        const isFinalStage = intersectX <= xb + buffer;
        if (isFinalStage) {
            stageLineData.push([intersectX, currentY]); // Final horizontal
            // Drop to the operating line (not y=x) to avoid visually crossing the drawn operating line
            let finalOpY: number;
            if (intersectX >= xIntersect) {
                finalOpY = rectifyingSlope * intersectX + rectifyingIntercept;
            } else if (intersectX >= xb) {
                finalOpY = (strippingSlope === Infinity) ? xb : strippingSlope * intersectX + strippingIntercept;
            } else {
                finalOpY = intersectX; // past xb: draw to diagonal
            }
            stageLineData.push([intersectX, finalOpY]); // Final vertical
            stageLineData.push([null, null]);
            break;
        } else {
            const currentSectionIsRectifying = intersectX >= xIntersect;
            let nextY;
            if (currentSectionIsRectifying) {
                nextY = rectifyingSlope * intersectX + rectifyingIntercept;
            } else {
                nextY = (strippingSlope === Infinity) ? xb : strippingSlope * intersectX + strippingIntercept;
            }
            nextY = Math.max(0, Math.min(1, nextY));

            stageLineData.push([intersectX, currentY]); // Horizontal to equilibrium
            stageLineData.push([intersectX, nextY]);   // Vertical to operating line
            stageLineData.push([null, null]);
            
            if (feedStageCount === 0 && currentSectionIsRectifying !== previousSectionIsRectifying) {
                feedStageCount = stageCount;
            }
            previousSectionIsRectifying = currentSectionIsRectifying;
            currentX = intersectX;
            currentY = nextY;
        }
    }
    
    setStages(stageCount);
    if (feedStageCount === 0 && stageCount > 0) setFeedStage(stageCount);
    else setFeedStage(feedStageCount);
    
    // --- PASS 2: Assemble All Chart Series With Finalized Data ---
    series.push({ name: 'Equilibrium Line', type: 'line', data: xValues.map((x, i) => [x, yValues[i]]), color: 'cyan', symbol: 'none', lineStyle: { width: 3.5 }, z: 5, animation: false });
    series.push({ name: 'y = x Line', type: 'line', data: [[0, 0], [1, 1]], color: 'blue', symbol: 'none', lineStyle: { width: 3.5, type: 'dotted' }, animation: false });

    if (murphreeEfficiency < 1.0) {
        const plotData = xValues
            .map((x, i) => (x >= finalStageX && x <= xd) ? [x, yDataForStepping[i]] as EChartsPoint : null)
            .filter(p => p !== null);

        const yAtFinalStageX = interpolateY(finalStageX, xValues, yDataForStepping);
        if (yAtFinalStageX !== null) {
            plotData.unshift([finalStageX, yAtFinalStageX]);
        }
        
        series.push({ name: 'Effective Equilibrium', type: 'line', data: plotData, connectNulls: false, color: 'cyan', symbol: 'none', lineStyle: { width: 3.5, type: 'dashed' }, z: 4, animation: false });
    }

    const feedLineData: EChartsPoint[] = [[xf, xf]];
    if (feedSlope === Infinity) { feedLineData.push([xf, yIntersect]); } else { feedLineData.push([xIntersect, yIntersect]); }
    series.push({ name: 'Rectifying Section', type: 'line', data: [[xd, xd], [xIntersect, yIntersect]], color: 'orange', symbol: 'none', lineStyle: { width: 3.5 }, animation: false });
    series.push({ name: 'Feed Section', type: 'line', data: feedLineData, color: 'red', symbol: 'none', lineStyle: { width: 3.5 }, animation: false });
    series.push({ name: 'Stripping Section', type: 'line', data: [[xIntersect, yIntersect], [xb, xb]], color: 'green', symbol: 'none', lineStyle: { width: 3.5 }, animation: false });
    series.push({ name: 'Key Points', type: 'scatter', data: [{ value: [xd, xd], itemStyle: { color: 'orange' }, name: 'Distillate' },{ value: [xb, xb], itemStyle: { color: 'green' }, name: 'Bottoms' },{ value: [xf, xf], itemStyle: { color: 'red' }, name: 'Feed' }], symbolSize: 8, label: { show: false }, tooltip: { formatter: (params: any) => `${params.name}: (${params.value[0].toFixed(3)}, ${params.value[1].toFixed(3)})` }, animation: false });
    series.push({ name: 'Stages', type: 'line', data: stageLineData, color: 'black', symbol: 'none', lineStyle: { width: 3 }, connectNulls: false, legendHoverLink: false, animation: false });

    // --- FINAL Rmin Calculation ---
    let rMinValue: number | null = null;
    let maxOperatingSlope = -Infinity;

    // Step 1: Calculate q-line/equilibrium intersection
    let intersectionPinch: { x: number; y: number } | null = null;
    if (Math.abs(q - 1) < 1e-9) {
        const y_val = interpolateY(xf, xValues, yValues);
        if (y_val !== null) {
            intersectionPinch = { x: xf, y: y_val };
        }
    } else {
        for (let i = 0; i < xValues.length - 1; i++) {
            const x1 = xValues[i], y1 = yValues[i];
            const x2 = xValues[i + 1], y2 = yValues[i + 1];
            if (x1 === x2) continue;
            const m_eq = (y2 - y1) / (x2 - x1);
            const c_eq = y1 - m_eq * x1;
            const denominator = (q - 1) * m_eq - q;
            if (Math.abs(denominator) < 1e-9) continue;
            const x_int = (-(q - 1) * c_eq - xf) / denominator;
            if (x_int >= Math.min(x1, x2) - 1e-6 && x_int <= Math.max(x1, x2) + 1e-6) {
                const y_int = m_eq * x_int + c_eq;
                intersectionPinch = { x: x_int, y: y_int };
                break;
            }
        }
    }

    if (intersectionPinch && intersectionPinch.x <= xd) {
        if (Math.abs(xd - intersectionPinch.x) > 1e-9) {
            maxOperatingSlope = (xd - intersectionPinch.y) / (xd - intersectionPinch.x);
        }
    }

    // Step 2: Check for a more limiting tangency
    if (q >= 1 || q <= 0) {
        // ✅ FIXED: Start loop at i = 1 to skip the invalid (0, 0) point.
        for (let i = 1; i < xValues.length; i++) {
            const x_eq = xValues[i];
            const y_eq = yValues[i];

            if (x_eq > xd + 1e-6) continue;

            let isCandidate = false;
            if (q >= 1) {
                if (x_eq <= xf + 1e-6) isCandidate = true;
            } else { // q <= 0
                if (x_eq >= xf - 1e-6) isCandidate = true;
            }
            
            if (isCandidate) {
                if (Math.abs(xd - x_eq) < 1e-9) continue;
                const slope = (xd - y_eq) / (xd - x_eq);
                if (slope > maxOperatingSlope) {
                    maxOperatingSlope = slope;
                }
            }
        }
    }

    // Step 3: Final Calculation
    if (maxOperatingSlope > -Infinity && maxOperatingSlope < 1) {
        rMinValue = maxOperatingSlope / (1 - maxOperatingSlope);
    }
    setRMin(rMinValue !== null && rMinValue >= 0 ? rMinValue : null);

    const capitalizeFirst = (str: string) => str ? str.charAt(0).toUpperCase() + str.slice(1) : '';
    const dispComp1Cap = capitalizeFirst(displayedComp1);
    const dispComp2Cap = capitalizeFirst(displayedComp2);
    const titleCondition = displayedUseTemp ? 
        (displayedTemp !== null ? `${Number(displayedTemp.toFixed(1))} °C` : 'N/A Temp') : 
        (displayedPressure !== null ? `${Number(displayedPressure.toFixed(2))} bar` : 'N/A Pressure');
    const titleText = displayedComp1 && displayedComp2 ?
        `McCabe-Thiele: ${dispComp1Cap}/${dispComp2Cap} at ${titleCondition}` :
        'McCabe-Thiele Diagram';
    
    const axisCompLabel = dispComp1Cap;
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';

    setEchartsOptions({
        backgroundColor: 'transparent',
        title: { text: titleText, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
        grid: { left: '5%', right: '5%', bottom: '5%', top: '5%', containLabel: true },
        xAxis: {
            type: 'value', min: 0, max: 1, interval: 0.1,
            name: `Liquid Mole Fraction ${axisCompLabel} (x)`,
            nameLocation: 'middle', nameGap: 30, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
        },
        yAxis: {
            type: 'value', min: 0, max: 1, interval: 0.1,
            name: `Vapor Mole Fraction ${axisCompLabel} (y)`,
            nameLocation: 'middle', nameGap: 40, nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } }, axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans', formatter: '{value}' }, splitLine: { show: false },
        },
        legend: {
            orient: 'vertical', right: '2%', top: 'center',
            data: (() => {
                const names = series.map(s => s.name).filter((n): n is string => typeof n === 'string');
                const filtered = names.filter(name => name !== 'Key Points' && name !== 'Stages');
                let ordered: string[] = [];
                if (filtered.includes('Equilibrium Line')) ordered.push('Equilibrium Line');
                if (filtered.includes('Effective Equilibrium') && murphreeEfficiency < 1.0) ordered.push('Effective Equilibrium');
                filtered.forEach(name => {
                  if (name !== 'Equilibrium Line' && name !== 'Effective Equilibrium') ordered.push(name);
                });
                return ordered as string[];
            })(),
            textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
            itemWidth: 25, itemHeight: 2,
        },        tooltip: {
            show: true,
            trigger: 'axis',
            triggerOn: 'mousemove',
            backgroundColor: isDark ? '#1a1a2e' : '#ffffff',
            borderColor: isDark ? '#55aaff' : '#d1d5db',
            borderWidth: 1,
            padding: [8, 12],
            textStyle: {
                color: textColor,
                fontSize: 12,
                fontFamily: 'Merriweather Sans'
            },
            axisPointer: {
                type: 'line',
                lineStyle: {
                    color: isDark ? '#55aaff' : '#999',
                    type: 'dashed',
                    width: 1,
                }
            },
            formatter: function (params: any) {
                if (!Array.isArray(params) || params.length === 0) return '';
                const xVal = params[0].axisValue;
                if (xVal === undefined) return '';

                const lines: string[] = [];
                lines.push(`<span style="font-weight:600;">x = ${formatNumberToPrecision(xVal, 4)}</span>`);

                // Series to show: VLE curve, operating lines only (not Stages, not y=x)
                const SHOW_SERIES = new Set(['VLE Curve', 'Rectifying Line', 'Stripping Line', 'Feed Line (q-line)']);

                // Key points (scatter) — only if exactly on them
                params.forEach((param: any) => {
                    if (param.seriesName === 'Key Points' && param.data?.name) {
                        lines.push(`<span style="color:${param.color};">● ${param.data.name}: (${formatNumberToPrecision(param.value[0], 3)}, ${formatNumberToPrecision(param.value[1], 3)})</span>`);
                    }
                });

                // Curve/line y-values
                params.forEach((param: any) => {
                    if (!SHOW_SERIES.has(param.seriesName)) return;
                    if (!param.value || !Array.isArray(param.value) || param.value.length < 2) return;
                    const y = param.value[1];
                    if (y === null || y === undefined || isNaN(y)) return;
                    lines.push(`<span style="color:${param.color};">● ${param.seriesName}: ${formatNumberToPrecision(y, 4)}</span>`);
                });

                return lines.join('<br/>');
            }
        }, 
        animation: false,
        toolbox: {
            show: false
        },
        series: series,
    });
  }, [xd, xb, xf, q, r, murphreeEfficiency, displayedComp1, displayedComp2, displayedTemp, displayedPressure, displayedUseTemp, displayedFluidPackage, buffer, resolvedTheme]);

  useEffect(() => {
    if (equilibriumData?.x && equilibriumData?.y) {
      generateEChartsOptions(equilibriumData.x, equilibriumData.y);
    } else {
        setEchartsOptions({}); setStages(null); setFeedStage(null);
    }
  }, [equilibriumData, generateEChartsOptions]);

  const handleUpdateGraphClick = () => {
    calculateEquilibriumCurve();
  };


  
  const getFeedQualityState = () => {
    if (q === 1) return "Saturated Liquid";
    if (q === 0) return "Saturated Vapor";
    if (q > 1) return "Subcooled Liquid";
    if (q > 0 && q < 1) return "Partially Vaporized";
    return "Super Heated Vapor";
  };

  const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>, inputTarget?: 'comp1' | 'comp2') => {
    if (event.key === 'Enter') {
      event.preventDefault();
      if (inputTarget === 'comp1' && showComp1Suggestions && comp1Suggestions.length > 0) {
        // Potentially select first suggestion or handle differently
        // For now, just close suggestions and update graph
        setShowComp1Suggestions(false);
      } else if (inputTarget === 'comp2' && showComp2Suggestions && comp2Suggestions.length > 0) {
        setShowComp2Suggestions(false);
      }
      // Trigger graph update on Enter key, which now internally syncs numeric states
      handleUpdateGraphClick();
    } else if (event.key === 'Escape') {
        if (inputTarget === 'comp1') setShowComp1Suggestions(false);
        if (inputTarget === 'comp2') setShowComp2Suggestions(false);
    }
  };
  
  const handleSwapComponents = () => {
    const tempName = comp1Name;
    const tempDisplay = comp1DisplayName;
    setComp1Name(comp2Name);
    setComp1DisplayName(comp2DisplayName);
    setComp2Name(tempName);
    setComp2DisplayName(tempDisplay);
    setAutoGeneratePending(true); // trigger auto regenerate after swap
  };

  const updateCompositions = (type: 'xd' | 'xf' | 'xb', value: number) => {
    // Get the current state values to check against
    const currentXd = xd;
    const currentXf = xf;
    const currentXb = xb;

    let targetValue = value; // Start with the value from the slider event

    // Check constraints and adjust the targetValue if needed (snap to boundary)
    if (type === 'xd') {
        // Prevent xd from going below or equal to xf + buffer
        if (targetValue <= currentXf + buffer) {
            targetValue = currentXf + buffer; // Snap to the boundary
        }
        // Ensure targetValue doesn't exceed absolute max
        targetValue = Math.min(targetValue, 0.99);
        setXd(targetValue);
    } else if (type === 'xf') {
        // Prevent xf from going above or equal to xd - buffer
        if (targetValue >= currentXd - buffer) {
            targetValue = currentXd - buffer; // Snap to the boundary
        }
        // Prevent xf from going below or equal to xb + buffer
        if (targetValue <= currentXb + buffer) {
            targetValue = currentXb + buffer; // Snap to the boundary
        }
         // Ensure targetValue is within absolute bounds
        targetValue = Math.max(0.01, Math.min(targetValue, 0.99));
        setXf(targetValue);
    } else if (type === 'xb') {
        // Prevent xb from going above or equal to xf - buffer
        if (targetValue >= currentXf - buffer) {
            targetValue = currentXf - buffer; // Snap to the boundary
        }
        // Ensure targetValue doesn't go below absolute min
        targetValue = Math.max(targetValue, 0.01);
        setXb(targetValue);
    }

    // Note: We directly set the state for the changed slider.
    // The constraints prevent invalid values based on the *current* state of the other sliders.
    // This avoids complex cascading updates and relies on the user not being able to move sliders past the calculated boundaries.
  };

  // Fluid packages available for the currently-loaded compound pair.
  // Packages whose pure-component params are missing are hidden from the dropdown.
  const availablePackages = useMemo(() => {
    if (!comp1Data || !comp2Data) return ['wilson', 'uniquac', 'nrtl', 'pr', 'srk'];
    const pkgs: FluidPackageTypeMcCabe[] = [];
    if (comp1Data.wilsonParams && comp2Data.wilsonParams) pkgs.push('wilson');
    if (comp1Data.uniquacParams && comp2Data.uniquacParams) pkgs.push('uniquac');
    pkgs.push('nrtl'); // always available (UNIFAC fallback)
    if (comp1Data.prParams && comp2Data.prParams) pkgs.push('pr');
    if (comp1Data.srkParams && comp2Data.srkParams) pkgs.push('srk');
    return pkgs;
  }, [comp1Data, comp2Data]);

  // Auto-switch away from an unavailable package (e.g. compound changed).
  useEffect(() => {
    if (!availablePackages.includes(fluidPackage as FluidPackageTypeMcCabe)) {
      setFluidPackage('nrtl');
    }
  }, [availablePackages]); // eslint-disable-line react-hooks/exhaustive-deps

  // Helper to compute a 'nice' slider step given the maximum value
  const computeStep = (maxVal: number): number => {
    const desiredSteps = 100;
    const raw = maxVal / desiredSteps;
    const exponent = Math.floor(Math.log10(raw));
    const pow10 = Math.pow(10, exponent);
    const mant = raw / pow10;
    let niceMant: number;
    if (mant <= 1) niceMant = 1;
    else if (mant <= 2) niceMant = 2;
    else if (mant <= 5) niceMant = 5;
    else niceMant = 10;
    return niceMant * pow10;
  };

  // Silent recalculation when fixed condition slider changes
  const hasMounted = useRef(false);
  useEffect(() => {
      if (!hasMounted.current) {
          hasMounted.current = true;
          return;
      }
      // Use higher resolution (101 points => 0.01 spacing) for real-time slider updates
      calculateEquilibriumCurve(false, 101);
  }, [temperatureC, pressureBar]);

  // Enforce a minimum pressure of 0.1 bar when in pressure mode
  useEffect(() => {
    if (!useTemperature && pressureBar !== null && pressureBar < 0.1) {
      setPressureBar(0.1);
    }
  }, [useTemperature, pressureBar]);

  // When autoGeneratePending is set, wait for the relevant state to update (next render) then regenerate the graph
  useEffect(() => {
    if (autoGeneratePending) {
      handleUpdateGraphClick();
      setAutoGeneratePending(false);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [comp1Name, comp2Name, useTemperature, fluidPackage, autoGeneratePending]);

  const handleDownloadChart = () => {
    const chart = echartsRef.current?.getEchartsInstance();
    if (!chart) return;
    const url = chart.getDataURL({ type: 'png', pixelRatio: 2, backgroundColor: resolvedTheme === 'dark' ? '#0f172a' : '#ffffff' });
    const a = document.createElement('a');
    a.href = url;
    a.download = 'chart.png';
    a.click();
  };

  // Build and display a sensitivity chart for the given variable.
  const runSensitivity = useCallback(async (varName: SensitivityVar) => {
    setActiveSensVar(varName);
    setSensitivityRunning(true);
    // Yield so React paints the "running" state before the heavy synchronous computation.
    await new Promise(resolve => setTimeout(resolve, 0));

    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';
    const stageParams = { xd, xb, xf, q, r, murphreeEfficiency, buffer };

    let minVal: number, maxVal: number, xLabel: string;
    if (varName === 'tc_or_p') {
      if (useTemperature) { minVal = niceMin(antoineLimits?.minTC ?? 0); maxVal = parseFloat(tempMax) || 200; xLabel = 'Temperature (°C)'; }
      else { minVal = niceMin(antoineLimits?.minPBar ?? 0.1); maxVal = parseFloat(pressureMax) || 10; xLabel = 'Pressure (bar)'; }
    } else if (varName === 'xd') { minVal = 0.01; maxVal = 0.99; xLabel = 'Distillate (xD)'; }
    else if (varName === 'xf') { minVal = 0.01; maxVal = 0.99; xLabel = 'Feed (xF)'; }
    else if (varName === 'xb') { minVal = 0.01; maxVal = 0.99; xLabel = 'Bottoms (xB)'; }
    else if (varName === 'q') { minVal = parseFloat(qMinSlider) || -1; maxVal = parseFloat(qMaxSlider) || 2; xLabel = 'Feed Quality (q)'; }
    else if (varName === 'r') { minVal = parseFloat(rMinSlider) || 0.1; maxVal = parseFloat(rMaxSlider) || 10; xLabel = 'Reflux Ratio (R)'; }
    else { minVal = 1; maxVal = 100; xLabel = 'Murphree Efficiency (%)'; }

    // Compute a nice round step size targeting ~200 points across the range.
    const TARGET_POINTS = 200;
    const rawStep = (maxVal - minVal) / (TARGET_POINTS - 1);
    const stepExp = Math.floor(Math.log10(rawStep));
    const stepPow = Math.pow(10, stepExp);
    const stepMant = rawStep / stepPow;
    let niceStep: number;
    if (stepMant <= 1) niceStep = 1 * stepPow;
    else if (stepMant <= 2) niceStep = 2 * stepPow;
    else if (stepMant <= 5) niceStep = 5 * stepPow;
    else niceStep = 10 * stepPow;
    // Build x list: always include minVal, then every niceStep from the snapped start.
    const snappedStart = Math.ceil((minVal + niceStep * 0.0001) / niceStep) * niceStep;
    const xAxisValues: number[] = [parseFloat(minVal.toPrecision(8))];
    for (let v = snappedStart; v <= maxVal + niceStep * 0.0001; v = parseFloat((v + niceStep).toPrecision(14))) {
      xAxisValues.push(parseFloat(v.toPrecision(8)));
    }
    if (xAxisValues[xAxisValues.length - 1] < maxVal - niceStep * 0.001) {
      xAxisValues.push(parseFloat(maxVal.toPrecision(8)));
    }

    const stagesArr: (number | null)[] = [];
    const feedStagesArr: (number | null)[] = [];

    for (const v of xAxisValues) {
      let eqX: number[], eqY: number[];
      if (varName === 'tc_or_p') {
        if (!comp1Data || !comp2Data || !interactionParams) { stagesArr.push(null); feedStagesArr.push(null); continue; }
        const pts = computeEquilibriumPointsSync(comp1Data, comp2Data, fluidPackage, interactionParams, useTemperature, v, 80);
        if (!pts) { stagesArr.push(null); feedStagesArr.push(null); continue; }
        eqX = pts.x; eqY = pts.y;
      } else {
        if (!equilibriumData) { stagesArr.push(null); feedStagesArr.push(null); continue; }
        eqX = equilibriumData.x; eqY = equilibriumData.y;
      }
      const overrideParams = { ...stageParams };
      if (varName === 'xd') overrideParams.xd = v;
      else if (varName === 'xf') overrideParams.xf = v;
      else if (varName === 'xb') overrideParams.xb = v;
      else if (varName === 'q') overrideParams.q = v;
      else if (varName === 'r') overrideParams.r = v;
      else if (varName === 'murphree') overrideParams.murphreeEfficiency = v / 100;
      const result = countStagesForSweep(eqX, eqY, overrideParams);
      stagesArr.push(result.stages > 0 ? result.stages : null);
      feedStagesArr.push(result.feedStage > 0 ? result.feedStage : null);
    }

    const sensOptions: EChartsOption = {
      backgroundColor: 'transparent',
      animation: false,
      title: { text: `Sensitivity: Stages vs ${xLabel}`, left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
      legend: { bottom: 5, data: ['Stages', 'Feed Stage'], textStyle: { color: textColor, fontFamily: 'Merriweather Sans' } },
      grid: { left: '5%', right: '5%', bottom: '10%', top: '12%', containLabel: true },
      xAxis: { type: 'value' as const, name: xLabel, nameLocation: 'middle' as const, nameGap: 30, nameTextStyle: { color: textColor, fontFamily: 'Merriweather Sans' }, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans' }, axisLine: { lineStyle: { color: textColor } }, splitLine: { show: false } },
      yAxis: { type: 'value' as const, name: 'Stages / Feed Stage', nameLocation: 'middle' as const, nameGap: 40, minInterval: 1, nameTextStyle: { color: textColor, fontFamily: 'Merriweather Sans' }, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans' }, axisLine: { lineStyle: { color: textColor } }, splitLine: { show: false } },
      series: [
        { name: 'Stages', type: 'line' as const, data: xAxisValues.map((v, i) => [v, stagesArr[i]]), color: 'cyan', symbol: 'circle', symbolSize: 4, lineStyle: { width: 3 }, connectNulls: false },
        { name: 'Feed Stage', type: 'line' as const, data: xAxisValues.map((v, i) => [v, feedStagesArr[i]]), color: 'orange', symbol: 'circle', symbolSize: 4, lineStyle: { width: 3 }, connectNulls: false },
      ],
      tooltip: { trigger: 'axis' as const, backgroundColor: isDark ? '#08306b' : '#ffffff', borderColor: isDark ? '#55aaff' : '#333333', textStyle: { color: textColor, fontFamily: 'Merriweather Sans' } },
    };
    setSensitivityEchartsOptions(sensOptions);
    setSensitivityRunning(false);
  }, [equilibriumData, comp1Data, comp2Data, interactionParams, fluidPackage, useTemperature,
      xd, xb, xf, q, r, murphreeEfficiency, buffer, antoineLimits, tempMax, pressureMax,
      qMinSlider, qMaxSlider, rMinSlider, rMaxSlider, resolvedTheme]);

  // Auto-rerun sensitivity when sliders change while a sweep is active.
  useEffect(() => {
    if (activeSensVar === null) return;
    runSensitivity(activeSensVar);
  }, [activeSensVar, xd, xb, xf, q, r, murphreeEfficiency, equilibriumData]); // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-4 md:px-16">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardContent className="space-y-6 p-4">
                <div className="space-y-4">
                  {/* Temperature/Pressure Toggle and Input */}
                  <div>
                    <Tabs value={useTemperature ? "temperature" : "pressure"} onValueChange={(value) => { 
                        setUseTemperature(value === "temperature");
                        setAutoGeneratePending(true); // trigger auto regenerate after mode change
                    }}>
                      <TabsList className="flex w-full">
                        <TabsTrigger value="temperature" className="flex-1 text-center">Temperature</TabsTrigger>
                        <TabsTrigger value="pressure" className="flex-1 text-center">Pressure</TabsTrigger>
                      </TabsList>
                    </Tabs>
                  </div>

                  {/* Component Inputs - Side by side with placeholders and swap button */}
                  <div className="flex items-center gap-2">
                    <div className="relative flex-1">
                      <Input
                        ref={input1Ref}
                        id="comp1Name"
                        value={comp1DisplayName}
                        onChange={handleComp1NameChange}
                        onKeyDown={(e) => handleKeyDown(e, 'comp1')}
                        onFocus={() => {
                            setActiveSuggestionInput('comp1');
                            if (comp1DisplayName.trim() !== "" && comp1Suggestions.length > 0) setShowComp1Suggestions(true);
                            else if (comp1DisplayName.trim() !== "") fetchSuggestions(comp1DisplayName, 'comp1');
                        }}
                        placeholder="Methanol"
                        required
                        className="w-full"
                        autoComplete="off"
                      />
                      {showComp1Suggestions && comp1Suggestions.length > 0 && (
                        <div ref={activeSuggestionInput === 'comp1' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                          {comp1Suggestions.map((alias, index) => (
                            <div
                              key={index}
                              onClick={() => handleSuggestionClick(alias, 'comp1')}
                              className="px-3 py-2 hover:bg-accent cursor-pointer text-sm"
                            >
                              {formatCompoundName(alias.fullName)}
                            </div>
                          ))}
                        </div>
                      )}
                    </div>
                    <Button variant="ghost" size="icon" onClick={handleSwapComponents} title="Swap Components">
                        <ArrowLeftRight className="h-4 w-4" />
                    </Button>
                    <div className="relative flex-1">
                      <Input
                        ref={input2Ref}
                        id="comp2Name"
                        value={comp2DisplayName}
                        onChange={handleComp2NameChange}
                        onKeyDown={(e) => handleKeyDown(e, 'comp2')}
                        onFocus={() => {
                            setActiveSuggestionInput('comp2');
                            if (comp2DisplayName.trim() !== "" && comp2Suggestions.length > 0) setShowComp2Suggestions(true);
                            else if (comp2DisplayName.trim() !== "") fetchSuggestions(comp2DisplayName, 'comp2');
                        }}
                        placeholder="Water"
                        required
                        className="w-full"
                        autoComplete="off"
                      />
                      {showComp2Suggestions && comp2Suggestions.length > 0 && (
                        <div ref={activeSuggestionInput === 'comp2' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                          {comp2Suggestions.map((alias, index) => (
                            <div
                              key={index}
                              onClick={() => handleSuggestionClick(alias, 'comp2')}
                              className="px-3 py-2 hover:bg-accent cursor-pointer text-sm"
                            >
                              {formatCompoundName(alias.fullName)}
                            </div>
                          ))}
                        </div>
                      )}
                    </div>
                  </div>
                  {/* Fluid Package Selector - Label side-by-side, reordered */}
                  <div className="flex items-center gap-2">
                      <Label htmlFor="fluidPackageMcCabe" className="text-sm font-medium whitespace-nowrap">Fluid Package:</Label>
                      <Select value={fluidPackage} onValueChange={(value) => {
                          setFluidPackage(value as FluidPackageTypeMcCabe);
                          setAutoGeneratePending(true); // trigger auto regenerate after fluid package change
                      }}>
                          <SelectTrigger id="fluidPackageMcCabe" className="flex-1"><SelectValue placeholder="Select fluid package" /></SelectTrigger>
                          <SelectContent>
                              {availablePackages.includes('wilson') && <SelectItem value="wilson">Wilson</SelectItem>}
                              {availablePackages.includes('uniquac') && <SelectItem value="uniquac">UNIQUAC</SelectItem>}
                              <SelectItem value="nrtl">NRTL</SelectItem>
                              {availablePackages.includes('pr') && <SelectItem value="pr">Peng-Robinson</SelectItem>}
                              {availablePackages.includes('srk') && <SelectItem value="srk">SRK</SelectItem>}
                          </SelectContent>
                      </Select>
                  </div>
                </div>
                <Button onClick={handleUpdateGraphClick} disabled={loading} className="w-full">
                  {loading ? 'Calculating...' : 'Update Graph'}
                </Button>
                {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
              </CardContent>
            </Card>
            {/* Slider Card */}
            <Card>
               <CardContent className="space-y-8 pt-2 pb-2">
                {/* Sensitivity Analysis Toggle */}
                <Button
                  variant={sensitivityMode ? 'default' : 'outline'}
                  className="w-full"
                  onClick={() => {
                    if (sensitivityMode) {
                      setSensitivityMode(false);
                      setActiveSensVar(null);
                      setSensitivityEchartsOptions(null);
                    } else {
                      setSensitivityMode(true);
                      setActiveSensVar(null);
                    }
                  }}
                >
                  {sensitivityMode ? 'Exit Sensitivity Analysis' : 'Run Sensitivity Analysis'}
                </Button>
                {/* Fixed condition slider */}
                {activeSensVar === 'tc_or_p' ? null : (
                <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                  <div className={`space-y-2${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                  <div className="flex items-center justify-between">
                    <Label htmlFor="fixedCondition">
                      {useTemperature ? `Temperature (°C): ${formatNumberToPrecision(temperatureC ?? 0, 4)}` : `Pressure (bar): ${formatNumberToPrecision(pressureBar ?? 0, 4)}`}
                    </Label>
                    <div className="flex items-center gap-1">
                      <Label className="text-xs text-muted-foreground">Max:</Label>
                      <Input
                        type="text"
                        value={useTemperature ? tempMax : pressureMax}
                        onChange={(e: React.ChangeEvent<HTMLInputElement>) => {
                          const raw = e.target.value;
                          const num = parseFloat(raw);
                          if (useTemperature) {
                            const CLAMP_MAX = antoineLimits?.maxTC ?? 2000;
                            if (raw === "" || isNaN(num)) setTempMax(raw);
                            else setTempMax(String(Math.min(num, CLAMP_MAX)));
                          } else {
                            const CLAMP_MAX = antoineLimits?.maxPBar ?? 2000;
                            if (raw === "" || isNaN(num)) setPressureMax(raw);
                            else setPressureMax(String(Math.min(num, CLAMP_MAX)));
                          }
                        }}
                        className="w-16 h-8 text-xs"
                      />
                    </div>
                  </div>
                  <Slider
                    id="fixedCondition"
                    min={useTemperature ? -100 : 0.1}
                    max={useTemperature ? (parseFloat(tempMax) || 200) : (parseFloat(pressureMax) || 10)}
                    step={computeStep(useTemperature ? (parseFloat(tempMax) || 200) : (parseFloat(pressureMax) || 10))}
                    value={[useTemperature ? (temperatureC || 0) : (pressureBar || 0)]}
                    onValueChange={([v]: number[]) => {
                      if (useTemperature) setTemperatureC(v);
                      else setPressureBar(v);
                    }}
                    className="w-full"
                  />
                  </div>
                  {sensitivityMode && activeSensVar === null && (
                    <div className="absolute inset-0 flex">
                      <Button variant="outline" className="w-full" onClick={() => runSensitivity('tc_or_p')} disabled={sensitivityRunning}>
                        {`Sweep ${useTemperature ? 'Temperature' : 'Pressure'}`}
                      </Button>
                    </div>
                  )}
                </div>
                )}
                 {/* xd Slider */}
                 {activeSensVar === 'xd' ? null : (
                 <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                   <div className={`space-y-3${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                   <Label htmlFor="xd">
                     <span dangerouslySetInnerHTML={{ __html: `Distillate Composition (x<sub>D</sub>): ${xd.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xd" min={0.01} max={0.99} step={0.01} value={[xd]} onValueChange={(value) => updateCompositions('xd', value[0])} style={{ '--primary': 'hsl(38 92% 50%)' } as React.CSSProperties}/>
                   </div>
                   {sensitivityMode && activeSensVar === null && (
                     <div className="absolute inset-0 flex">
                       <Button variant="outline" className="w-full" onClick={() => runSensitivity('xd')} disabled={sensitivityRunning}>
                         Sweep Distillate (xD)
                       </Button>
                     </div>
                   )}
                 </div>
                 )}
                 {/* xf Slider */}
                 {activeSensVar === 'xf' ? null : (
                 <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                   <div className={`space-y-3${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                   <Label htmlFor="xf">
                     <span dangerouslySetInnerHTML={{ __html: `Feed Composition (x<sub>F</sub>): ${xf.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xf" min={0.01} max={0.99} step={0.01} value={[xf]} onValueChange={(value) => updateCompositions('xf', value[0])} style={{ '--primary': 'hsl(0 84% 60%)' } as React.CSSProperties}/>
                   </div>
                   {sensitivityMode && activeSensVar === null && (
                     <div className="absolute inset-0 flex">
                       <Button variant="outline" className="w-full" onClick={() => runSensitivity('xf')} disabled={sensitivityRunning}>
                         Sweep Feed (xF)
                       </Button>
                     </div>
                   )}
                 </div>
                 )}
                 {/* xb Slider */}
                 {activeSensVar === 'xb' ? null : (
                 <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                   <div className={`space-y-3${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                   <Label htmlFor="xb">
                     <span dangerouslySetInnerHTML={{ __html: `Bottoms Composition (x<sub>B</sub>): ${xb.toFixed(2)}` }} />
                   </Label>
                   <Slider id="xb" min={0.01} max={0.99} step={0.01} value={[xb]} onValueChange={(value) => updateCompositions('xb', value[0])} style={{ '--primary': 'hsl(142 71% 45%)' } as React.CSSProperties}/>
                   </div>
                   {sensitivityMode && activeSensVar === null && (
                     <div className="absolute inset-0 flex">
                       <Button variant="outline" className="w-full" onClick={() => runSensitivity('xb')} disabled={sensitivityRunning}>
                         Sweep Bottoms (xB)
                       </Button>
                     </div>
                   )}
                 </div>
                 )}
                 {/* q Slider */}
                 {activeSensVar === 'q' ? null : (
                 <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                   <div className={`space-y-3${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                   <div className="flex items-center justify-between">
                     <Label htmlFor="q" className="flex items-center">Feed Quality (q): {q.toFixed(2)}
                         <ShadTooltip><TooltipTrigger asChild><Button variant="ghost" size="icon" className="h-5 w-5 rounded-full"><span className="text-xs">ⓘ</span></Button></TooltipTrigger><TooltipContent><p>{getFeedQualityState()}</p></TooltipContent></ShadTooltip>
                     </Label>
                     <div className="flex items-center gap-1">
                       <Label className="text-xs text-muted-foreground">Min:</Label>
                       <Input
                         type="text"
                         value={qMinSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setQMinSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                       <Label className="text-xs text-muted-foreground">Max:</Label>
                       <Input
                         type="text"
                         value={qMaxSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setQMaxSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                     </div>
                   </div>
                   <Slider id="q" min={parseFloat(qMinSlider)} max={parseFloat(qMaxSlider)} step={0.05} value={[q]} onValueChange={(value) => setQ(value[0])} style={{ '--primary': 'hsl(60 100% 50%)' } as React.CSSProperties}/>
                   </div>
                   {sensitivityMode && activeSensVar === null && (
                     <div className="absolute inset-0 flex">
                       <Button variant="outline" className="w-full" onClick={() => runSensitivity('q')} disabled={sensitivityRunning}>
                         Sweep Feed Quality (q)
                       </Button>
                     </div>
                   )}
                 </div>
                 )}
                 {/* r Slider */}
                 {activeSensVar === 'r' ? null : (
                 <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                   <div className={`space-y-3${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                   <div className="flex justify-between items-center">
                     <Label htmlFor="r" className="flex justify-between items-center">
                       <span>Reflux Ratio (R): {r.toFixed(2)}</span>
                       {rMin !== null && (
                          <span className="text-xs text-muted-foreground ml-2">
                              R<sub>min</sub>: {rMin.toFixed(2)}
                          </span>
                       )}
                     </Label>
                     <div className="flex items-center gap-1">
                       <Label className="text-xs text-muted-foreground">Min:</Label>
                       <Input
                         type="text"
                         value={rMinSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setRMinSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                       <Label className="text-xs text-muted-foreground">Max:</Label>
                       <Input
                         type="text"
                         value={rMaxSlider}
                         onChange={(e: React.ChangeEvent<HTMLInputElement>) => setRMaxSlider(e.target.value)}
                         className="w-12 h-8 text-xs"
                       />
                     </div>
                   </div>
                   <Slider id="r" min={parseFloat(rMinSlider)} max={parseFloat(rMaxSlider)} step={0.05} value={[r]} onValueChange={(value) => setR(value[0])} style={{ '--primary': 'hsl(262 84% 58%)' } as React.CSSProperties}/>
                   </div>
                   {sensitivityMode && activeSensVar === null && (
                     <div className="absolute inset-0 flex">
                       <Button variant="outline" className="w-full" onClick={() => runSensitivity('r')} disabled={sensitivityRunning}>
                         Sweep Reflux Ratio (R)
                       </Button>
                     </div>
                   )}
                 </div>
                 )}
                 {/* Murphree Efficiency Slider */}
                 {activeSensVar === 'murphree' ? null : (
                 <div className={`relative${sensitivityMode && activeSensVar === null ? ' h-10 overflow-hidden' : ''}`}>
                   <div className={`space-y-3${sensitivityMode && activeSensVar === null ? ' invisible' : ''}`}>
                   <Label htmlFor="murphree" className="flex justify-between items-center">
                     <span dangerouslySetInnerHTML={{ __html: `Murphree Efficiency (E<sub>M</sub>): ${(murphreeEfficiency * 100).toFixed(0)}%` }} />
                   </Label>
                   <Slider
                     id="murphree"
                     min={0.01}
                     max={1}
                     step={0.01}
                     value={[murphreeEfficiency]}
                     onValueChange={(value) => setMurphreeEfficiency(value[0])}
                     style={{ '--primary': 'hsl(180 100% 50%)' } as React.CSSProperties}
                   />
                   </div>
                   {sensitivityMode && activeSensVar === null && (
                     <div className="absolute inset-0 flex">
                       <Button variant="outline" className="w-full" onClick={() => runSensitivity('murphree')} disabled={sensitivityRunning}>
                         Sweep Murphree Efficiency
                       </Button>
                     </div>
                   )}
                 </div>
                 )}
               </CardContent>
            </Card>
            {stages !== null && feedStage !== null && !sensitivityMode && (
              <Card>
                <CardContent>
                  <div className="grid grid-cols-2 gap-4 text-center">
                    <div className="p-4 bg-muted rounded-md">
                      <p className="text-sm font-medium">Number of Stages</p>
                      <p className="text-2xl font-bold">{stages}</p>
                    </div>
                    <div className="p-4 bg-muted rounded-md">
                      <p className="text-sm font-medium">Feed Stage</p>
                      <p className="text-2xl font-bold">{feedStage}</p>
                    </div>
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
          {/* Column 2: Plot and Results Cards */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2">
                <div className="relative aspect-square rounded-md">
                   {/* Sensitivity Analysis Chart */}
                   {activeSensVar !== null && sensitivityEchartsOptions && (
                     <ReactECharts
                       ref={echartsRef}
                       echarts={echarts}
                       option={sensitivityEchartsOptions}
                       style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }}
                       notMerge={true}
                       lazyUpdate={false}
                     />
                   )}
                   {/* Normal McCabe-Thiele Chart */}
                   {activeSensVar === null && (
                     <>
                   {!equilibriumData && !error && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Please provide inputs and update graph.</div> )}
                   {error && !loading && ( <div className="absolute inset-0 flex items-center justify-center bg-background/80 text-red-400 p-4 text-center rounded-md">{error}</div> )}
                  {equilibriumData && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts 
                  key={`echarts-${murphreeEfficiency === 1.0 ? 'no-eff' : 'with-eff'}`}
                  ref={echartsRef} 
                  echarts={echarts} 
                  option={echartsOptions} 
                  style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }} 
                  notMerge={true} 
                  lazyUpdate={false} 
                />
                  )}
                     </>
                   )}
                  {(activeSensVar !== null ? (sensitivityEchartsOptions !== null) : (Object.keys(echartsOptions).length > 0)) && (
                    <button onClick={handleDownloadChart} className="absolute bottom-2 right-2 z-10 p-1.5 rounded bg-background/80 hover:bg-background border border-border text-muted-foreground hover:text-foreground transition-colors" title="Download chart as PNG"><Download className="h-4 w-4" /></button>
                  )}
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}