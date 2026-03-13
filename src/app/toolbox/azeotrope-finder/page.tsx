'use client';

import React, { useState, useEffect, useCallback, useRef, useMemo } from 'react';
import { supabase } from '@/lib/supabaseClient';
import { type SupabaseClient as _SupabaseClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";

// ECharts imports (ensure all necessary components are registered as in other pages)
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption } from 'echarts';
import { LineChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  DataZoomComponent,
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([
  TitleComponent, TooltipComponent, GridComponent, LegendComponent, ToolboxComponent, DataZoomComponent, LineChart, CanvasRenderer
]);

// Shadcn UI imports
import { Card, CardContent } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { TooltipProvider } from "@/components/ui/tooltip"; // Assuming you might use Shadcn tooltips
import { ArrowLeftRight, Terminal, Download } from 'lucide-react';
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs"; // Added Tabs

// VLE Calculation Library – consolidated
import {
    calculatePsat_Pa as libCalculatePsat_Pa,
    // UNIFAC
    calculateBubbleTemperatureUnifac,
    calculateBubblePressureUnifac,
    fetchUnifacInteractionParams,
    // NRTL
    fetchNrtlParameters,
    calculateBubbleTemperatureNrtl,
    calculateBubblePressureNrtl,
    // Peng–Robinson
    fetchPrInteractionParams,
    calculateBubbleTemperaturePr,
    calculateBubblePressurePr,
    // SRK
    fetchSrkInteractionParams,
    calculateBubbleTemperatureSrk,
    calculateBubblePressureSrk,
    // UNIQUAC
    fetchUniquacInteractionParams,
    calculateBubbleTemperatureUniquac,
    calculateBubblePressureUniquac,
    // Wilson
    fetchWilsonInteractionParams,
    calculateBubbleTemperatureWilson,
    calculateBubblePressureWilson,
    antoineBoilingPointSolverLocal
} from '@/lib/vle-calculations';

// Shared VLE Types
import type {
    CompoundData, BubbleDewResult
} from '@/lib/vle-types';
import { fetchCompoundDataFromHysys, fetchCompoundSuggestions, resolveSimName, formatCompoundName } from '@/lib/antoine-utils';
import type { CompoundAlias } from '@/lib/antoine-utils';

type FluidPackageTypeAzeotrope = 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac' | 'wilson';
type AzeotropeScanType = 'vs_P_find_T' | 'vs_T_find_P'; // Scan P, find T_az OR Scan T, find P_az


// Debounce function removed for instantaneous suggestions

// formatNumberToPrecision function (copy from mccabe-thiele)
const formatNumberToPrecision = (num: any, precision: number = 3): string => {
    if (typeof num === 'number') {
      if (num === 0) return '0';
      const fixed = num.toPrecision(precision);
      if (fixed.includes('.')) { return parseFloat(fixed).toString(); }
      return fixed;
    }
    return String(num);
};


export default function AzeotropeFinderPage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  // NOTE: Caching is disabled to avoid stale data carrying over between
  // fluid-package switches.
  
  // Input States
  const [comp1Name, setComp1Name] = useState('Acetone'); // Default to acetone
  const [comp2Name, setComp2Name] = useState('WATER');   // Default to water
  const [fluidPackage, setFluidPackage] = useState<FluidPackageTypeAzeotrope>('uniquac'); // Default to uniquac
  const [azeotropeScanType, setAzeotropeScanType] = useState<AzeotropeScanType>('vs_P_find_T'); // Default to Scan Pressure
  // const scanSteps = 50; // Always use 50 scan points, no UI to change this // REMOVED

  // Data & Control States
  const [azeotropeScanData, setAzeotropeScanData] = useState<{ scanVal: number, x_az: number, dependentVal: number, branchIndex: number }[]>([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  // const [logMessages, setLogMessages] = useState<string[]>([]); // Removed logMessages

  // ECharts State
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null);

  // Suggestion states (copy from mccabe-thiele)
  const [comp1Suggestions, setComp1Suggestions] = useState<CompoundAlias[]>([]);
  const [comp2Suggestions, setComp2Suggestions] = useState<CompoundAlias[]>([]);
  const [showComp1Suggestions, setShowComp1Suggestions] = useState(false);
  const [showComp2Suggestions, setShowComp2Suggestions] = useState(false);
  const [activeSuggestionInput, setActiveSuggestionInput] = useState<'comp1' | 'comp2' | null>(null);
  const input1Ref = useRef<HTMLInputElement>(null);
  const input2Ref = useRef<HTMLInputElement>(null);
  const activeComponentRef = useRef<HTMLDivElement>(null);

  // Display names (human-readable FullName from HYSYS ALIASES)
  const [comp1DisplayName, setComp1DisplayName] = useState('Acetone');
  const [comp2DisplayName, setComp2DisplayName] = useState('Water');
  
  // Displayed parameters for title
  const [displayedComp1, setDisplayedComp1] = useState('');
  const [displayedComp2, setDisplayedComp2] = useState('');
  const [displayedFluidPackage, setDisplayedFluidPackage] = useState<FluidPackageTypeAzeotrope | ''>('');
  const [displayedScanType, setDisplayedScanType] = useState<AzeotropeScanType | ''>('');
  // Concurrency guard so stale scans don't overwrite newer results
  const scanGenerationRef = useRef(0);
  // Cached compound data + interaction params — avoids re-fetching on scan-type or mode changes
  const [cachedComp1Data, setCachedComp1Data] = useState<CompoundData | null>(null);
  const [cachedComp2Data, setCachedComp2Data] = useState<CompoundData | null>(null);
  const [cachedInteractionParams, setCachedInteractionParams] = useState<any>(null);
  const lastFetchKeyRef = useRef<string | null>(null);
  // Pre-fetched compound data used to determine which fluid packages have the
  // required parameters for the current compound pair.
  const [comp1PreData, setComp1PreData] = useState<CompoundData | null>(null);
  const [comp2PreData, setComp2PreData] = useState<CompoundData | null>(null);

  // This new useEffect will handle re-scanning when the fluid package changes.
  const didMountFluidPkg = useRef(false);
  useEffect(() => {
    if (didMountFluidPkg.current) {
      handleAzeotropeScan();
    } else {
      didMountFluidPkg.current = true;
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [fluidPackage]);

  // ADD THIS NEW EFFECT
  // This new useEffect will handle re-scanning when the scan type changes.
  const didMountScanType = useRef(false);
  useEffect(() => {
    if (didMountScanType.current) {
        handleAzeotropeScan();
    } else {
        didMountScanType.current = true;
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [azeotropeScanType]);

  // Derive which fluid packages have the necessary parameters for the loaded pair.
  // Only packages whose required data is present in both compounds are shown.
  const availablePackages = useMemo((): Set<FluidPackageTypeAzeotrope> => {
    const d1 = comp1PreData, d2 = comp2PreData;
    // Before either compound is resolved, allow everything so the UI isn't blank.
    if (!d1 || !d2) return new Set(['unifac', 'nrtl', 'wilson', 'uniquac', 'pr', 'srk']);
    const pkgs = new Set<FluidPackageTypeAzeotrope>();
    pkgs.add('nrtl'); // NRTL always falls back to UNIFAC estimate
    if (d1.unifacGroups && d2.unifacGroups) pkgs.add('unifac');
    if (d1.uniquacParams && d2.uniquacParams) pkgs.add('uniquac');
    if (d1.wilsonParams && d2.wilsonParams) pkgs.add('wilson');
    if (d1.prParams && d2.prParams) pkgs.add('pr');
    if (d1.srkParams && d2.srkParams) pkgs.add('srk');
    return pkgs;
  }, [comp1PreData, comp2PreData]);

  // Auto-switch to the best available fluid package when the current one becomes unsupported.
  useEffect(() => {
    if (comp1PreData && comp2PreData && !availablePackages.has(fluidPackage)) {
      const preference: FluidPackageTypeAzeotrope[] = ['uniquac', 'nrtl', 'wilson', 'unifac', 'pr', 'srk'];
      const fallback = preference.find(p => availablePackages.has(p));
      if (fallback) setFluidPackage(fallback);
    }
  // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [availablePackages]);

  // --- Logging Hook (copy from test/page.tsx or mccabe-thiele) ---
  // Removed logging useEffect and related state

  // --- useEffect for initial graph load ---
  useEffect(() => {
    if (supabase && comp1Name && comp2Name) {
        handleAzeotropeScan();
        // Pre-fetch compound data so the fluid-package dropdown can be filtered immediately.
        fetchCompoundDataFromHysys(supabase, comp1Name).then(setComp1PreData).catch(() => {});
        fetchCompoundDataFromHysys(supabase, comp2Name).then(setComp2PreData).catch(() => {});
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [supabase]); // Run once when supabase client is ready

  // --- Data Fetching (fetchCompoundDataLocal - uses HYSYS tables) ---
  async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
    console.log(`AzeotropeFinder: Fetching data for ${compoundName} from DB...`);
    if (!supabase) { throw new Error("Supabase client not initialized."); }
    try {
        const result = await fetchCompoundDataFromHysys(supabase, compoundName);
        return result;
    } catch (err: any) {
        console.error(`Error fetching data for ${compoundName} in AzeotropeFinder:`, err.message);
        setError(`Data fetch failed for ${compoundName}: ${err.message}`);
        return null;
    }
  }

  // --- Suggestion Handlers (using HYSYS ALIASES table) ---
  const fetchSuggestions = useCallback(async (inputValue: string, inputTarget: 'comp1' | 'comp2') => {
    if (!inputValue || inputValue.length < 2 || !supabase) {
      if (inputTarget === 'comp1') { setComp1Suggestions([]); setShowComp1Suggestions(false); }
      else { setComp2Suggestions([]); setShowComp2Suggestions(false); }
      return;
    }
    try {
      const suggestions = await fetchCompoundSuggestions(supabase, inputValue);
      if (inputTarget === 'comp1') { setComp1Suggestions(suggestions); setShowComp1Suggestions(suggestions.length > 0); }
      else { setComp2Suggestions(suggestions); setShowComp2Suggestions(suggestions.length > 0); }
    } catch (err) { console.error("Azeo: Error fetching suggestions:", err); if (inputTarget === 'comp1') setComp1Suggestions([]); else setComp2Suggestions([]); }
  }, [supabase]);

  // Removed debouncedFetchSuggestions for instantaneous suggestions

  const handleComp1NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp1Name(newValue);
    setComp1DisplayName(newValue);
    setActiveSuggestionInput('comp1');
    if (newValue.trim() === "") { setShowComp1Suggestions(false); setComp1Suggestions([]); }
    else { fetchSuggestions(newValue, 'comp1'); }
  };

  const handleComp2NameChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newValue = e.target.value;
    setComp2Name(newValue);
    setComp2DisplayName(newValue);
    setActiveSuggestionInput('comp2');
    if (newValue.trim() === "") { setShowComp2Suggestions(false); setComp2Suggestions([]); }
    else { fetchSuggestions(newValue, 'comp2'); }
  };

  const handleSuggestionClick = (alias: CompoundAlias, inputTarget: 'comp1' | 'comp2') => {
    if (inputTarget === 'comp1') {
      setComp1Name(alias.simName); setComp1DisplayName(formatCompoundName(alias.fullName)); setShowComp1Suggestions(false); setComp1Suggestions([]);
      if (supabase) fetchCompoundDataFromHysys(supabase, alias.simName).then(setComp1PreData).catch(() => {});
    } else {
      setComp2Name(alias.simName); setComp2DisplayName(formatCompoundName(alias.fullName)); setShowComp2Suggestions(false); setComp2Suggestions([]);
      if (supabase) fetchCompoundDataFromHysys(supabase, alias.simName).then(setComp2PreData).catch(() => {});
    }
  };

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
    return () => { document.removeEventListener("mousedown", handleClickOutside); };
  }, [activeSuggestionInput]);

  const handleSwapComponents = () => { 
    const tempName = comp1Name;
    const tempDisplay = comp1DisplayName;
    setComp1Name(comp2Name);
    setComp1DisplayName(comp2DisplayName);
    setComp2Name(tempName);
    setComp2DisplayName(tempDisplay);
    // Auto-trigger azeotrope scan after swap
    setTimeout(() => handleAzeotropeScan(), 100);
  };

  // --- Enter / Escape key handler for compound inputs ---
  const handleKeyDown = (event: React.KeyboardEvent<HTMLInputElement>, inputTarget: 'comp1' | 'comp2') => {
    if (event.key === 'Enter') {
      event.preventDefault();
      setShowComp1Suggestions(false);
      setShowComp2Suggestions(false);
      handleAzeotropeScan();
    } else if (event.key === 'Escape') {
      if (inputTarget === 'comp1') setShowComp1Suggestions(false);
      if (inputTarget === 'comp2') setShowComp2Suggestions(false);
    }
  };

  // --- CSV Export Function ---
  const handleExportToCSV = () => {
    if (azeotropeScanData.length === 0) {
      alert("No data to export.");
      return;
    }

    let csvContent = "";
    let scanValHeader = "";
    let xAzHeader = `Azeotropic ${displayedComp1 || 'Comp1'} Composition (x1)`;
    let dependentValHeader = "";

    if (displayedScanType === 'vs_P_find_T') {
      scanValHeader = "Pressure (bar)";
      dependentValHeader = "Azeotropic Temperature (°C)";
    } else { // vs_T_find_P
      scanValHeader = "Temperature (°C)";
      dependentValHeader = "Azeotropic Pressure (bar)";
    }

    const headers = [scanValHeader, xAzHeader, dependentValHeader, "Branch"];
    csvContent += headers.join(",") + "\r\n";

    azeotropeScanData.forEach(row => {
      const scanValFormatted = formatNumberToPrecision(row.scanVal, 4);
      const xAzFormatted = formatNumberToPrecision(row.x_az, 4);
      const dependentValFormatted = (displayedScanType === 'vs_P_find_T')
        ? formatNumberToPrecision(row.dependentVal, 2)
        : formatNumberToPrecision(row.dependentVal, 4);

      const csvRow = [scanValFormatted, xAzFormatted, dependentValFormatted, String(row.branchIndex + 1)].join(",");
      csvContent += csvRow + "\r\n";
    });

    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement("a");
    if (link.download !== undefined) {
      const url = URL.createObjectURL(blob);
      const filename = `azeotrope_locus_${displayedComp1 || 'Comp1'}_${displayedComp2 || 'Comp2'}_${displayedFluidPackage || 'model'}.csv`;
      link.setAttribute("href", url);
      link.setAttribute("download", filename);
      link.style.visibility = 'hidden';
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
    }
  };

  // --- Bisection Solver (copy from test/page.tsx) ---
  function bisectionSolve(
    func: (x: number) => number | null, a: number, b: number, tolerance: number = 1e-5, maxIter: number = 50
  ): number | null {
    // ... (Full implementation from test/page.tsx)
    let fa = func(a); let fb = func(b);
    if (fa === null || fb === null || fa * fb >=0) {
        if (fa !== null && Math.abs(fa) < tolerance) return a;
        if (fb !== null && Math.abs(fb) < tolerance) return b;
        return null;
    }
    let c = a;
    for (let i = 0; i < maxIter; i++) {
        c = (a+b)/2; const fc = func(c);
        if (fc === null || Math.abs(fc) < tolerance || (b-a)/2 < tolerance) return (fc === null && Math.abs(func(c-tolerance) || 1) > tolerance && Math.abs(func(c+tolerance) || 1) > tolerance) ? null : c; // check if fc is null due to bad point
        if (fa * fc < 0) { b = c; } else { a = c; fa = fc; }
    }
    return null;
  }




  // --- Main Azeotrope Scan Logic ---
  const handleAzeotropeScan = useCallback(async () => {
    const myScanId = ++scanGenerationRef.current;
    setLoading(true); setError(null);

    if (!comp1Name || !comp2Name) { setError("Please enter compound names."); setLoading(false); return; }
    if (!supabase) { setError("Supabase not initialized."); setLoading(false); return; }

    try {
        // Resolve display names to SimNames
        const [sim1, sim2] = await Promise.all([resolveSimName(supabase, comp1Name), resolveSimName(supabase, comp2Name)]);
        const resolvedComp1 = sim1 || comp1Name;
        const resolvedComp2 = sim2 || comp2Name;

        // Use cache when compound names + fluid package haven't changed
        const fetchKey = `${resolvedComp1}|${resolvedComp2}|${fluidPackage}`;
        let data1 = cachedComp1Data;
        let data2 = cachedComp2Data;
        let activityParameters: any = cachedInteractionParams;
        const needsFetch = !data1 || !data2 || !activityParameters || lastFetchKeyRef.current !== fetchKey;

        if (needsFetch) {
            const [fetched1, fetched2] = await Promise.all([fetchCompoundDataLocal(resolvedComp1), fetchCompoundDataLocal(resolvedComp2)]);
            if (!fetched1 || !fetched2) { setLoading(false); return; }
            data1 = fetched1; data2 = fetched2;

            // Check for pure component system
            if (fetched1.name.toLowerCase() === fetched2.name.toLowerCase()) {
                setError("Cannot find an azeotrope for a pure component. Please select two different compounds.");
                setLoading(false);
                return;
            }

            // Validate required parameters based on fluidPackage
            if (fluidPackage === 'unifac' && (!fetched1.unifacGroups || !fetched2.unifacGroups)) throw new Error("UNIFAC groups missing.");
            if (fluidPackage === 'pr' && (!fetched1.prParams || !fetched2.prParams)) throw new Error("PR parameters missing.");
            if (fluidPackage === 'srk' && (!fetched1.srkParams || !fetched2.srkParams)) throw new Error("SRK parameters missing.");
            if (fluidPackage === 'uniquac' && (!fetched1.uniquacParams || !fetched2.uniquacParams)) throw new Error("UNIQUAC parameters missing.");
            if (fluidPackage === 'wilson' && (!fetched1.wilsonParams || !fetched2.wilsonParams)) throw new Error("Wilson parameters missing.");

            // Fetch interaction parameters
            if (fluidPackage === 'unifac') {
                const allSubgroupIds = new Set<number>();
                [fetched1, fetched2].forEach(comp => { if (comp.unifacGroups) Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id))); });
                activityParameters = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds));
            } else if (fluidPackage === 'nrtl') {
                activityParameters = await fetchNrtlParameters(supabase, fetched1.name, fetched2.name);
                if ((activityParameters as any)?._usedUnifacFallback) {
                    setError(`NRTL parameters not found in database for ${fetched1.name}–${fetched2.name}. Using UNIFAC-based estimate.`);
                }
            } else if (fluidPackage === 'pr') {
                activityParameters = await fetchPrInteractionParams(supabase, fetched1.name, fetched2.name);
            } else if (fluidPackage === 'srk') {
                activityParameters = await fetchSrkInteractionParams(supabase, fetched1.name, fetched2.name);
            } else if (fluidPackage === 'uniquac') {
                activityParameters = await fetchUniquacInteractionParams(supabase, fetched1.name, fetched2.name);
                if ((activityParameters as any)?._usedUnifacFallback) {
                    setError(`UNIQUAC parameters not found in database for ${fetched1.name}–${fetched2.name}. Using UNIFAC-based estimate.`);
                }
            } else if (fluidPackage === 'wilson') {
                activityParameters = await fetchWilsonInteractionParams(supabase, fetched1.name, fetched2.name);
            } else {
                throw new Error(`Unsupported fluid package: ${fluidPackage}`);
            }

            // Update cache
            setCachedComp1Data(fetched1);
            setCachedComp2Data(fetched2);
            setCachedInteractionParams(activityParameters);
            lastFetchKeyRef.current = fetchKey;
            // Also keep PreData in sync for fluid-package dropdown
            setComp1PreData(fetched1);
            setComp2PreData(fetched2);
        } else {
            // Still need to validate for pure component check when using cache
            if (data1!.name.toLowerCase() === data2!.name.toLowerCase()) {
                setError("Cannot find an azeotrope for a pure component.");
                setLoading(false);
                return;
            }
        }

        // data1 and data2 are guaranteed non-null after the cache/fetch block above
        const d1 = data1!;
        const d2 = data2!;
        const components: CompoundData[] = [d1, d2];

        const scanResultsArray: { scanVal: number, x_az: number, dependentVal: number, branchIndex: number }[] = [];
        
        // Derive scan bounds from Antoine parameters and critical properties
        let scanPMinBar = 0.01, scanPMaxBar = 20;
        let scanTMinK = 230, scanTMaxK = 550;
        if (d1.antoine && d2.antoine) {
            const Pc1 = d1.Pc_Pa ?? d1.prParams?.Pc_Pa ?? Infinity;
            const Pc2 = d2.Pc_Pa ?? d2.prParams?.Pc_Pa ?? Infinity;
            const Tc1 = d1.Tc_K ?? d1.prParams?.Tc_K ?? Infinity;
            const Tc2 = d2.Tc_K ?? d2.prParams?.Tc_K ?? Infinity;
            const Pmin1 = libCalculatePsat_Pa(d1.antoine, d1.antoine.Tmin_K);
            const Pmin2 = libCalculatePsat_Pa(d2.antoine, d2.antoine.Tmin_K);
            const Pmax1 = Math.min(libCalculatePsat_Pa(d1.antoine, d1.antoine.Tmax_K), isFinite(Pc1) ? Pc1 : Infinity);
            const Pmax2 = Math.min(libCalculatePsat_Pa(d2.antoine, d2.antoine.Tmax_K), isFinite(Pc2) ? Pc2 : Infinity);
            const computedPMin = Math.max(Pmin1, Pmin2) / 1e5;
            const computedPMax = Math.min(isFinite(Pmax1) ? Pmax1 : 20e5, isFinite(Pmax2) ? Pmax2 : 20e5) / 1e5;
            if (computedPMax > computedPMin) { scanPMinBar = Math.max(computedPMin, 0.01); scanPMaxBar = computedPMax; }

            const computedTMin = Math.max(d1.antoine.Tmin_K, d2.antoine.Tmin_K);
            const computedTMax = Math.min(
                Math.min(d1.antoine.Tmax_K, isFinite(Tc1) ? Tc1 : Infinity),
                Math.min(d2.antoine.Tmax_K, isFinite(Tc2) ? Tc2 : Infinity)
            );
            if (computedTMax > computedTMin) { scanTMinK = computedTMin; scanTMaxK = computedTMax; }
        }

        let currentScanStartValue: number, currentScanEndValue: number, stepValue: number;
        if (azeotropeScanType === 'vs_P_find_T') {
            currentScanStartValue = scanPMinBar;
            currentScanEndValue = scanPMaxBar;
            stepValue = 0.1; // bar increment
        } else { // vs_T_find_P
            currentScanStartValue = scanTMinK;
            currentScanEndValue = scanTMaxK;
            stepValue = 1; // K increment
        }

        // Use index-based loop to avoid floating-point accumulation in the loop variable
        const totalSteps = Math.ceil((currentScanEndValue - currentScanStartValue) / stepValue);
        const scanValues = Array.from({ length: totalSteps + 1 }, (_, i) =>
            parseFloat((currentScanStartValue + i * stepValue).toPrecision(10))
        ).filter(v => v <= currentScanEndValue + stepValue * 0.5);

        for (const currentScanVal of scanValues) {
            // (iterates clean rounded values with no FP drift)

            const objectiveFunction = (x1_guess: number): number | null => {
                if (x1_guess <= 1e-5 || x1_guess >= 1.0 - 1e-5) return 1.0; // Avoid edges for stability

                let bubbleResult: BubbleDewResult | null = null;
                if (azeotropeScanType === 'vs_P_find_T') {
                    const P_system_Pa = currentScanVal * 100000; // Convert bar to Pa
                    const T_bp1_est = antoineBoilingPointSolverLocal(d1.antoine, P_system_Pa) || 350;
                    const T_bp2_est = antoineBoilingPointSolverLocal(d2.antoine, P_system_Pa) || 350;
                    const initialTempGuess = x1_guess * T_bp1_est + (1 - x1_guess) * T_bp2_est;

                    if (fluidPackage === 'unifac') bubbleResult = calculateBubbleTemperatureUnifac(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'nrtl') bubbleResult = calculateBubbleTemperatureNrtl(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'pr') bubbleResult = calculateBubbleTemperaturePr(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'srk') bubbleResult = calculateBubbleTemperatureSrk(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'uniquac') bubbleResult = calculateBubbleTemperatureUniquac(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'wilson') bubbleResult = calculateBubbleTemperatureWilson(components, x1_guess, P_system_Pa, activityParameters, initialTempGuess);

                } else { // vs_T_find_P
                    const T_system_K = currentScanVal;
                    const P_sat1_est = libCalculatePsat_Pa(d1.antoine!, T_system_K);
                    const P_sat2_est = libCalculatePsat_Pa(d2.antoine!, T_system_K);
                    const initialPressureGuess = (P_sat1_est && P_sat2_est) ? (x1_guess * P_sat1_est + (1 - x1_guess) * P_sat2_est) : 101325;
                    
                    if (fluidPackage === 'unifac') bubbleResult = calculateBubblePressureUnifac(components, x1_guess, T_system_K, activityParameters);
                    else if (fluidPackage === 'nrtl') bubbleResult = calculateBubblePressureNrtl(components, x1_guess, T_system_K, activityParameters);
                    else if (fluidPackage === 'pr') bubbleResult = calculateBubblePressurePr(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'srk') bubbleResult = calculateBubblePressureSrk(components, x1_guess, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'uniquac') bubbleResult = calculateBubblePressureUniquac(components, x1_guess, T_system_K, activityParameters);
                    else if (fluidPackage === 'wilson') bubbleResult = calculateBubblePressureWilson(components, x1_guess, T_system_K, activityParameters);
                }

                if (bubbleResult && bubbleResult.error === undefined && typeof bubbleResult.comp1_equilibrium === 'number' && !isNaN(bubbleResult.comp1_equilibrium)) {
                    return x1_guess - bubbleResult.comp1_equilibrium;
                }
                return null; // Indicate failure to solver
            };

            // Sample across composition space to find ALL sign changes (multiple azeotropes)
            const NUM_GRID = 100;
            const xGrid = Array.from({ length: NUM_GRID }, (_, k) => 0.0001 + k * (0.9998 / (NUM_GRID - 1)));
            const fGrid = xGrid.map(x => objectiveFunction(x));

            const foundRoots: number[] = [];
            for (let k = 0; k < NUM_GRID - 1; k++) {
                const fa = fGrid[k];
                const fb = fGrid[k + 1];
                if (fa === null || fb === null) continue;
                if (fa * fb < 0) {
                    const root = bisectionSolve(objectiveFunction, xGrid[k], xGrid[k + 1], 1e-5, 50);
                    if (root !== null && root > 0.001 && root < 0.999) {
                        // Deduplicate roots that are numerically identical
                        if (!foundRoots.some(r => Math.abs(r - root) < 0.005)) {
                            foundRoots.push(root);
                        }
                    }
                }
            }
            foundRoots.sort((a, b) => a - b);

            for (let branchIdx = 0; branchIdx < foundRoots.length; branchIdx++) {
                const x_az_found = foundRoots[branchIdx];
                let dependent_val_found: number | null = null;
                let finalBubbleResult: BubbleDewResult | null = null;

                if (azeotropeScanType === 'vs_P_find_T') {
                    const P_system_Pa = currentScanVal * 100000;
                    const T_bp1_est = antoineBoilingPointSolverLocal(d1.antoine, P_system_Pa) || 350;
                    const T_bp2_est = antoineBoilingPointSolverLocal(d2.antoine, P_system_Pa) || 350;
                    const initialTempGuess = x_az_found * T_bp1_est + (1 - x_az_found) * T_bp2_est;

                    if (fluidPackage === 'unifac') finalBubbleResult = calculateBubbleTemperatureUnifac(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'nrtl') finalBubbleResult = calculateBubbleTemperatureNrtl(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'pr') finalBubbleResult = calculateBubbleTemperaturePr(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'srk') finalBubbleResult = calculateBubbleTemperatureSrk(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'uniquac') finalBubbleResult = calculateBubbleTemperatureUniquac(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);
                    else if (fluidPackage === 'wilson') finalBubbleResult = calculateBubbleTemperatureWilson(components, x_az_found, P_system_Pa, activityParameters, initialTempGuess);

                    if (finalBubbleResult && finalBubbleResult.T_K) {
                        dependent_val_found = finalBubbleResult.T_K - 273.15;
                    }
                } else {
                    const T_system_K = currentScanVal;
                    const P_sat1_est = libCalculatePsat_Pa(d1.antoine!, T_system_K);
                    const P_sat2_est = libCalculatePsat_Pa(d2.antoine!, T_system_K);
                    const initialPressureGuess = (P_sat1_est && P_sat2_est) ? (x_az_found * P_sat1_est + (1 - x_az_found) * P_sat2_est) : 101325;

                    if (fluidPackage === 'unifac') finalBubbleResult = calculateBubblePressureUnifac(components, x_az_found, T_system_K, activityParameters);
                    else if (fluidPackage === 'nrtl') finalBubbleResult = calculateBubblePressureNrtl(components, x_az_found, T_system_K, activityParameters);
                    else if (fluidPackage === 'pr') finalBubbleResult = calculateBubblePressurePr(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'srk') finalBubbleResult = calculateBubblePressureSrk(components, x_az_found, T_system_K, activityParameters, initialPressureGuess);
                    else if (fluidPackage === 'uniquac') finalBubbleResult = calculateBubblePressureUniquac(components, x_az_found, T_system_K, activityParameters);
                    else if (fluidPackage === 'wilson') finalBubbleResult = calculateBubblePressureWilson(components, x_az_found, T_system_K, activityParameters);

                    if (finalBubbleResult && finalBubbleResult.P_Pa) dependent_val_found = finalBubbleResult.P_Pa / 100000;
                }

                if (dependent_val_found !== null) {
                    scanResultsArray.push({ scanVal: azeotropeScanType === 'vs_T_find_P' ? currentScanVal - 273.15 : currentScanVal, x_az: x_az_found, dependentVal: dependent_val_found, branchIndex: branchIdx });
                }
            }
        }
        if (myScanId === scanGenerationRef.current) {
            setAzeotropeScanData(scanResultsArray);
            setDisplayedComp1(comp1DisplayName || comp1Name);
            setDisplayedComp2(comp2DisplayName || comp2Name);
            setDisplayedFluidPackage(fluidPackage);
            setDisplayedScanType(azeotropeScanType);
        }

    } catch (err: any) {
        console.error("Azeotrope Scan failed:", err);
        setError(`Scan failed: ${err.message}`);
    } finally {
        if (myScanId === scanGenerationRef.current) {
            setLoading(false);
            // setLogMessages(prev => ["--- Azeotrope Scan Complete ---", ...prev]); // Logging removed
        }
    }
  }, [comp1Name, comp2Name, fluidPackage, azeotropeScanType, supabase]);

  // --- ECharts Plotting ---
  const generateAzeotropeEChartsOptions = useCallback(() => {
    if (azeotropeScanData.length === 0) {
      setEchartsOptions({});
      return;
    }

    const comp1DisplayCap = displayedComp1 ? displayedComp1.charAt(0).toUpperCase() + displayedComp1.slice(1) : "Comp1";
    const scanParamName = displayedScanType === 'vs_P_find_T' ? "Azeotropic Pressure (bar)" : "Azeotropic Temperature (°C)";
    const dependentParamName = displayedScanType === 'vs_P_find_T' ? `Azeotropic Temperature (°C)` : `Azeotropic Pressure (bar)`;
    const compositionLegendName = `Azeotropic ${comp1DisplayCap} Comp. (x₁)`;
    const compositionAxisName = `Azeotropic ${comp1DisplayCap} Composition (x₁)`;
    const titleText = displayedComp1 && displayedComp2 ? `Azeotrope Locus: ${displayedComp1}/${displayedComp2}` : 'Azeotrope Locus';


    // Group data by branch index for separate series per azeotrope locus
    const branchIndices = [...new Set(azeotropeScanData.map(d => d.branchIndex))].sort((a, b) => a - b);
    const BRANCH_COLORS: [string, string][] = [
        ['#22c55e', '#f59e0b'], // branch 0: green / amber
        ['#3b82f6', '#ef4444'], // branch 1: blue / red
        ['#a855f7', '#f97316'], // branch 2: purple / orange
        ['#06b6d4', '#ec4899'], // branch 3: cyan / pink
    ];
    const chartSeries: any[] = [];
    for (const branchIdx of branchIndices) {
        const branchData = azeotropeScanData
            .filter(d => d.branchIndex === branchIdx)
            .sort((a, b) => a.scanVal - b.scanVal);
        const [compColor, depColor] = BRANCH_COLORS[branchIdx % BRANCH_COLORS.length];
        const branchLabel = branchIndices.length > 1 ? ` (Branch ${branchIdx + 1})` : '';
        chartSeries.push({
            name: compositionLegendName + branchLabel,
            type: 'line', yAxisIndex: 0,
            data: branchData.map(d => [d.scanVal, d.x_az]),
            smooth: true, lineStyle: { color: compColor, width: 2.5 }, itemStyle: { color: compColor }, symbol: 'none'
        });
        chartSeries.push({
            name: dependentParamName + branchLabel,
            type: 'line', yAxisIndex: 1,
            data: branchData.map(d => [d.scanVal, d.dependentVal]),
            smooth: true, lineStyle: { color: depColor, width: 2.5 }, itemStyle: { color: depColor }, symbol: 'none'
        });
    }

    // Theme-dependent colors
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';

    const options: EChartsOption = {
      title: { 
        text: titleText, 
        left: 'center', 
        top: '10',
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' },
      },
      backgroundColor: 'transparent',
      tooltip: {
        trigger: 'axis',
        axisPointer: {
          type: 'cross',
          label: {
            backgroundColor: resolvedTheme === 'dark' ? '#1e293b' : '#ffffff',
            color: textColor,
            fontFamily: 'Merriweather Sans',
            borderColor: resolvedTheme === 'dark' ? '#3b82f6' : '#333333',
            borderWidth: 1,
          }
        },
        backgroundColor: resolvedTheme === 'dark' ? '#1e293b' : '#ffffff',
        borderColor: resolvedTheme === 'dark' ? '#3b82f6' : '#333333',
        textStyle: { color: textColor, fontFamily: 'Merriweather Sans' },
        formatter: (params: any) => {
            if (!Array.isArray(params) || params.length === 0) return '';
            let tooltipHtml = `${scanParamName}: <b>${formatNumberToPrecision(params[0].axisValue, 4)}</b><br/>`;
            params.forEach((param: any) => {
                const seriesName = param.seriesName;
                const value = param.value[1];
                let precision = 3;
                if (seriesName === compositionLegendName) precision = 3;
                else if (displayedScanType === 'vs_T_find_P') precision = 4; 
                else if (displayedScanType === 'vs_P_find_T') precision = 3; 
                
                const formattedValue = (displayedScanType === 'vs_P_find_T' && seriesName === dependentParamName)
                    ? parseFloat(value.toFixed(1)).toString()
                    : formatNumberToPrecision(value, precision);

                tooltipHtml += `<span style="color: ${param.color};"><b>${seriesName}: ${formattedValue}</b></span><br/>`;
            });
            return tooltipHtml;
        }
      },
      legend: {
        show: true,
        orient: 'horizontal', bottom: '0', left: 'center',
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        itemWidth: 25, itemHeight: 4,
      },
      grid: { left: '8%', right: '10%', bottom: '8%', top: '12%', containLabel: true },
      xAxis: {
        show: true, 
        type: 'value', name: scanParamName, nameLocation: 'middle', nameGap: 30,
        axisLabel: { color: textColor, fontFamily: 'Merriweather Sans', fontSize: 16 },
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor, width: 2.5 } },
        splitLine: { show: false },
        min: azeotropeScanData.length > 0 ? parseFloat(Math.min(...azeotropeScanData.map(d => d.scanVal)).toPrecision(4)) : undefined,
        max: azeotropeScanData.length > 0 ? parseFloat(Math.max(...azeotropeScanData.map(d => d.scanVal)).toPrecision(4)) : undefined,
      },
      yAxis: [
        {
          show: true, 
          type: 'value', name: compositionAxisName, nameLocation: 'middle', nameGap: 50, min: 0, max: 1,
          position: 'left', 
          axisLabel: { color: textColor, formatter: (v: number) => formatNumberToPrecision(v,3), fontFamily: 'Merriweather Sans', fontSize: 16 }, // Changed color, Added fontSize
          nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' }, // Added fontSize
          axisLine: { lineStyle: { color: textColor, width: 2.5 } },
          splitLine: { show: false }, // Removed grid lines
        },
        {
          show: true, 
          type: 'value', name: dependentParamName, nameLocation: 'middle', nameGap: 60, position: 'right',
          axisLabel: { 
            color: textColor, // Changed color
            formatter: (v: number) => (displayedScanType === 'vs_P_find_T') ? parseFloat(v.toFixed(1)).toString() : formatNumberToPrecision(v, 3), 
            fontFamily: 'Merriweather Sans',
            fontSize: 16 // Added fontSize
          },
          nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' }, // Added fontSize
          axisLine: { lineStyle: { color: textColor, width: 2.5 } },
          splitLine: { show: false }, // Removed grid lines
          scale: true,
        }
      ],
      series: chartSeries,
      animation: false,
      toolbox: { show: false },
    };
    setEchartsOptions(options);
  }, [azeotropeScanData, displayedComp1, displayedComp2, displayedFluidPackage, displayedScanType, resolvedTheme]);

  useEffect(() => {
    generateAzeotropeEChartsOptions();
  }, [azeotropeScanData, generateAzeotropeEChartsOptions]);

  const handleDownloadChart = () => {
    const chart = echartsRef.current?.getEchartsInstance();
    if (!chart) return;
    const url = chart.getDataURL({ type: 'png', pixelRatio: 2, backgroundColor: resolvedTheme === 'dark' ? '#0f172a' : '#ffffff' });
    const a = document.createElement('a');
    a.href = url;
    a.download = 'chart.png';
    a.click();
  };

  // --- JSX Structure ---
  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Left Column: Inputs */}
          <div className="lg:col-span-1 space-y-6">
            <Card>
              {/* <CardHeader>{/* <CardTitle>Azeotrope Finder Inputs</CardTitle> * /}</CardHeader> */} {/* Removed CardHeader */}
              <CardContent className="p-4 space-y-4"> {/* Added p-4 for consistent padding */}
                <Tabs value={azeotropeScanType} onValueChange={(value) => {
                    setAzeotropeScanType(value as AzeotropeScanType);
                    // The direct call to handleAzeotropeScan has been removed.
                    // The new useEffect will trigger the scan.
                }} className="mb-4">
                    <TabsList className="grid w-full grid-cols-2">
                        <TabsTrigger value="vs_P_find_T">Scan Pressure</TabsTrigger>
                        <TabsTrigger value="vs_T_find_P">Scan Temperature</TabsTrigger>
                    </TabsList>
                </Tabs>
                {/* Component Inputs with Suggestions (similar to McCabe-Thiele) */}
                <div className="flex items-center gap-2">
                    <div className="relative flex-1">
                        {/* <Label htmlFor="comp1NameAzeo">Component 1</Label> */}
                        <Input ref={input1Ref} id="comp1NameAzeo" value={comp1DisplayName} onChange={handleComp1NameChange}
                               onKeyDown={(e) => handleKeyDown(e, 'comp1')}
                               onFocus={() => { setActiveSuggestionInput('comp1'); if (comp1DisplayName.trim()) fetchSuggestions(comp1DisplayName, 'comp1');}}
                               placeholder="e.g., Ethanol" autoComplete="off" />
                        {showComp1Suggestions && comp1Suggestions.length > 0 && (
                            <div ref={activeSuggestionInput === 'comp1' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                {comp1Suggestions.map((alias, i) => <div key={i} onClick={() => handleSuggestionClick(alias, 'comp1')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{formatCompoundName(alias.fullName)}</div>)}
                            </div>
                        )}
                    </div>
                    <Button variant="ghost" size="icon" onClick={handleSwapComponents} title="Swap Components" className="self-end mb-1"><ArrowLeftRight className="h-4 w-4" /></Button>
                    <div className="relative flex-1">
                        {/* <Label htmlFor="comp2NameAzeo">Component 2</Label> */}
                        <Input ref={input2Ref} id="comp2NameAzeo" value={comp2DisplayName} onChange={handleComp2NameChange}
                               onKeyDown={(e) => handleKeyDown(e, 'comp2')}
                               onFocus={() => { setActiveSuggestionInput('comp2'); if (comp2DisplayName.trim()) fetchSuggestions(comp2DisplayName, 'comp2');}}
                               placeholder="e.g., Water" autoComplete="off" />
                        {showComp2Suggestions && comp2Suggestions.length > 0 && (
                            <div ref={activeSuggestionInput === 'comp2' ? activeComponentRef : null} className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                {comp2Suggestions.map((alias, i) => <div key={i} onClick={() => handleSuggestionClick(alias, 'comp2')} className="px-3 py-2 hover:bg-accent cursor-pointer text-sm">{formatCompoundName(alias.fullName)}</div>)}
                            </div>
                        )}
                    </div>
                </div>

                <div className="flex items-center gap-2">
                  <Label htmlFor="fluidPackageAzeo" className="whitespace-nowrap">Fluid Package:</Label>
                  <Select value={fluidPackage} onValueChange={(value) => {
                    setFluidPackage(value as FluidPackageTypeAzeotrope);
                  }}>
                    <SelectTrigger id="fluidPackageAzeo" className="flex-1"><SelectValue /></SelectTrigger>
                                              <SelectContent>
                              {availablePackages.has('wilson') && <SelectItem value="wilson">Wilson</SelectItem>}
                              {availablePackages.has('uniquac') && <SelectItem value="uniquac">UNIQUAC</SelectItem>}
                              {availablePackages.has('nrtl') && <SelectItem value="nrtl">NRTL</SelectItem>}
                              {availablePackages.has('pr') && <SelectItem value="pr">Peng-Robinson</SelectItem>}
                              {availablePackages.has('srk') && <SelectItem value="srk">SRK</SelectItem>}
                          </SelectContent>
                  </Select>
                </div>
                {/* Scan Type Select removed */}
                {/* 
                <div>
                  <Label htmlFor="azeotropeScanTypeAzeo">Scan Type:</Label>
                  <Select value={azeotropeScanType} onValueChange={(value) => setAzeotropeScanType(value as AzeotropeScanType)}>
                    <SelectTrigger id="azeotropeScanTypeAzeo"><SelectValue /></SelectTrigger>
                    <SelectContent>
                      <SelectItem value="vs_P_find_T">Scan Pressure (find x_az, T_az)</SelectItem>
                      <SelectItem value="vs_T_find_P">Scan Temperature (find x_az, P_az)</SelectItem>
                    </SelectContent>
                  </Select>
                </div>
                */}
                {/* Scan Steps input removed - defaults to 50 */}
                {/* 
                <div className="flex items-center gap-2">
                  <Label htmlFor="scanStepsAzeo" className="whitespace-nowrap">Scan Points:</Label>
                  <Input id="scanStepsAzeo" type="number" value={scanSteps} onChange={(e) => setScanSteps(Math.max(5, parseInt(e.target.value)))} min="5" max="200" className="flex-1"/>
                </div>
                */}
                <Button onClick={handleAzeotropeScan} disabled={loading} className="w-full">
                  {loading ? 'Scanning...' : 'Find Azeotropes'}
                </Button>
                <Button 
                  onClick={handleExportToCSV} 
                  disabled={loading || azeotropeScanData.length === 0} 
                  className="w-full mt-2"
                  variant="outline"
                >
                  Export to CSV
                </Button>
                {error && <Alert variant="destructive" className="mt-2"><Terminal className="h-4 w-4" /><AlertTitle>Error</AlertTitle><AlertDescription>{error}</AlertDescription></Alert>}
              </CardContent>
            </Card>
            {/* Log Card Removed */}
            {/* 
            <Card>
                <CardHeader><CardTitle>Calculation Log</CardTitle></CardHeader>
                <CardContent><div className="h-48 overflow-y-auto bg-muted p-2 rounded text-xs font-mono">{logMessages.length === 0 ? <p>No messages yet...</p> : logMessages.map((msg, index) => <div key={index}>{msg}</div>)}</div></CardContent>
            </Card>
            */}
          </div>

          {/* Right Column: Plot */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardContent className="py-2"> {/* Minimal padding for graph card content */}
                <div className="relative aspect-square rounded-md"> {/* Plot area */}
                  {loading && (<div className="absolute inset-0 flex items-center justify-center text-muted-foreground"><div className="text-center"><div className="mb-2">Scanning for Azeotropes...</div><div className="text-sm text-muted-foreground/70">Using {fluidPackage.toUpperCase()} model.</div></div></div>)}
                  {!loading && !displayedComp1 && !error && (<div className="absolute inset-0 flex items-center justify-center text-muted-foreground">Please provide inputs and find azeotropes.</div>)}
                  {error && !loading && (<div className="absolute inset-0 flex items-center justify-center text-red-400">Error: {error}</div>)}
                  {!loading && displayedComp1 && !error && azeotropeScanData.length === 0 && (
                    <div className="absolute inset-0 flex items-center justify-center text-muted-foreground text-center px-8">
                      No azeotrope detected in the scanned range for these conditions.
                    </div>
                  )}
                  {!loading && azeotropeScanData.length > 0 && Object.keys(echartsOptions).length > 0 && (
                    <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={true} />
                  )}
                  {Object.keys(echartsOptions).length > 0 && <button onClick={handleDownloadChart} className="absolute bottom-2 right-2 z-10 p-1.5 rounded bg-background/80 hover:bg-background border border-border text-muted-foreground hover:text-foreground transition-colors" title="Download chart as PNG"><Download className="h-4 w-4" /></button>}
                </div>
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}
