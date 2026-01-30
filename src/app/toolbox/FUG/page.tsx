'use client';

import React, { useState, useEffect, useRef } from 'react';
import { createClient } from '@supabase/supabase-js';
import { useTheme } from "next-themes";
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import { LineChart, BarChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  MarkLineComponent,
  MarkPointComponent
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

echarts.use([
  TitleComponent, TooltipComponent, GridComponent, LegendComponent, 
  MarkLineComponent, MarkPointComponent, LineChart, BarChart, CanvasRenderer
]);

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { TooltipProvider } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Trash2, PlusCircle, Calculator, Search } from 'lucide-react';

// --- Types ---
type Compound = {
  id: number;
  name: string; // "CompoundID.value"
  formula: string; 
  antoine: { A: number; B: number; C: number }; // Parsed from "AntoineVaporPressure"
  z: number; // Mole fraction (user input)
  Tb: number; // Normal Boiling Point (K) for sorting/filtering
  alpha?: number; // Relative volatility (calculated)
  K?: number; // K-value (calculated)
  xD?: number; // Distillate composition (calculated)
  xB?: number; // Bottoms composition (calculated)
};

// --- Supabase Client ---
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || '';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || '';
const supabase = createClient(supabaseUrl, supabaseAnonKey);

export default function FUGSimulation() {
  const { theme, resolvedTheme } = useTheme(); // Use resolvedTheme for accurate dark/light detection
  const textColor = resolvedTheme === 'dark' ? '#ffffff' : '#000000'; // Dynamic text color
  const echartsRef = useRef<ReactECharts>(null);
  const [loading, setLoading] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');
  
  // Simulation State
  const [compounds, setCompounds] = useState<Compound[]>([]);
  const [pressure, setPressure] = useState<number>(1.01325); // Default to ~1 atm in Bar
  const [qValue, setQValue] = useState<number>(1.0); // Feed Quality (1=Liquid, 0=Vapor)
  const [lightKeyIdx, setLightKeyIdx] = useState<string>("");
  const [heavyKeyIdx, setHeavyKeyIdx] = useState<string>("");
  const [recLK, setRecLK] = useState<number>(0.99); // Recovery of LK in Distillate
  const [recHK, setRecHK] = useState<number>(0.99); // Recovery of HK in Bottoms
  const [R_factor, setRFactor] = useState<number>(1.2); // R / Rmin multiplier

  // Autocomplete suggestions for component search
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const suggestionsRef = useRef<HTMLDivElement>(null);

  // Ref to ensure initial auto-calc runs only once on mount when defaults loaded
  const initialCalcRef = useRef(false);

  // Results State
  const [results, setResults] = useState<{
    N_min: number;
    R_min: number;
    R_actual: number;
    N_actual: number;
    FeedStage: number;
    FeedTemp: number;
    GillilandData: number[][];
    OperatingPoint: number[];
  } | null>(null);

  const fetchSuggestions = async (query: string) => {
    if (!query || query.length < 2) {
      setSuggestions([]);
      setShowSuggestions(false);
      return;
    }
    const { data } = await supabase
      .from('compound_properties')
      .select('name')
      .ilike('name', `${query}%`)
      .limit(5);
    
    if (data) {
      setSuggestions(data.map((d: any) => d.name));
      setShowSuggestions(true);
    }
  };

  const handleSearchChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const val = e.target.value;
    setSearchQuery(val);
    fetchSuggestions(val);
  };

  const handleSuggestionClick = (name: string) => {
    setSearchQuery(name);
    setShowSuggestions(false);
    handleAddCompound(name); // FIX: Immediately add the specific compound
  };

  // Close dropdown when clicking outside
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (suggestionsRef.current && !suggestionsRef.current.contains(event.target as Node)) {
        setShowSuggestions(false);
      }
    };
    document.addEventListener("mousedown", handleClickOutside);
    return () => document.removeEventListener("mousedown", handleClickOutside);
  }, []);

  // Calculate normal boiling point from Antoine at 1 atm (Pa)
  const calculateTb = (A: number, B: number, C: number) => {
    // ln(P) = A - B/(T+C)  =>  T = B/(A - ln(P)) - C
    const P = 101325;
    const val = A - Math.log(P);
    if (val === 0) return 373.15; // fallback
    return (B / val) - C;
  };

  // --- Helper: Fetch Compound ---
  const handleAddCompound = async (specificName?: string) => {
    // FIX: Limit to 7 compounds
    if (compounds.length >= 7) {
        alert("Maximum of 7 compounds allowed.");
        return null;
    }

    const query = specificName || searchQuery;
    if (!query) return null;
    setLoading(true);

    // Step 1: Fetch potential matches (fetch a few to sort through them)
    // If a specific name was clicked, we search for THAT exact name.
    // If "Enter" was pressed (fuzzy), we search generally.
    let dbQuery: any = supabase
      .from('compound_properties') 
      .select('id, name, properties');

    if (specificName) {
        dbQuery = dbQuery.eq('name', specificName); // Exact match for dropdown clicks
    } else {
        dbQuery = dbQuery.ilike('name', `%${query}%`).limit(10); // Fuzzy for typing
    }

    const { data, error } = await dbQuery;

    if (error) {
        console.error("Supabase error:", error);
        setLoading(false);
        return null;
    }

    if (data && data.length > 0) {
        // FIX: Intelligent Selection Logic
        // If we did a fuzzy search, find the best match:
        // 1. Exact match (e.g. "Benzene" vs "1,2,3...Benzene")
        // 2. Starts with query
        // 3. Shortest string (usually the simplest compound)
        let row = data[0];
        
        if (!specificName) {
             const lowerQ = query.toLowerCase();
             row = data.find((d: any) => d.name.toLowerCase() === lowerQ) || 
                   data.find((d: any) => d.name.toLowerCase().startsWith(lowerQ)) ||
                   data.sort((a: any, b: any) => a.name.length - b.name.length)[0]; // Fallback to shortest
        }

        const props = row.properties; 

        if (props && props.AntoineVaporPressure) {
             const A = parseFloat(props.AntoineVaporPressure.A?.value ?? props.AntoineVaporPressure.A);
             const B = parseFloat(props.AntoineVaporPressure.B?.value ?? props.AntoineVaporPressure.B);
             const C = parseFloat(props.AntoineVaporPressure.C?.value ?? props.AntoineVaporPressure.C);
             
             const newCompound: Compound = {
                id: row.id, 
                name: row.name,
                formula: props.StructureFormula?.value || "N/A",
                antoine: { A, B, C },
                z: 0.1,
                Tb: calculateTb(A, B, C) // Calculate boiling point immediately
            };
            setCompounds(prev => [...prev, newCompound]);
            setSearchQuery('');
            setSuggestions([]); // Clear suggestions
            setLoading(false);
            return newCompound;
        } else {
            alert(`Compound "${row.name}" lacks Antoine Coefficients.`);
            setLoading(false);
            return null;
        }
    } else {
        // Silent fail for defaults, alert for user
        if (!specificName) alert("Compound not found.");
        setLoading(false);
        return null;
    }
  };

  // Auto-load defaults on mount
  useEffect(() => {
    const loadDefaults = async () => {
        // Prevent double loading if already populated
        if (compounds.length > 0) return;
        // Load sequentially to maintain order and capture IDs
        const m = await handleAddCompound("Methanol");
        const e = await handleAddCompound("Ethanol");
        const w = await handleAddCompound("Water");

        // Set default keys: Methanol (light) and Water (heavy) by their IDs if present
        if (m) setLightKeyIdx(m.id.toString());
        if (w) setHeavyKeyIdx(w.id.toString());

        // Trigger initial calculation after state settles
        setTimeout(() => {
          calculateFUG();
        }, 60);
    };
    loadDefaults();
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Run initial calculation once defaults are loaded and keys are set
  useEffect(() => {
    if (initialCalcRef.current) return;
    if (compounds.length >= 3 && lightKeyIdx && heavyKeyIdx) {
      initialCalcRef.current = true;
      // Slight delay so React state updates finish
      setTimeout(() => calculateFUG(), 40);
    }
  }, [compounds, lightKeyIdx, heavyKeyIdx]);

  // When there are fewer than 3 components, clear results so the graph hides
  useEffect(() => {
    if (compounds.length < 3) {
      setResults(null);
    }
  }, [compounds.length]);

  const removeCompound = (index: number) => {
    const newArr = [...compounds];
    const removed = newArr.splice(index, 1)[0];
    setCompounds(newArr);

    // If removed compound was selected as a key, clear that selection
    if (removed) {
      const removedId = removed.id;
      if (lightKeyIdx && parseInt(lightKeyIdx) === removedId) setLightKeyIdx("");
      if (heavyKeyIdx && parseInt(heavyKeyIdx) === removedId) setHeavyKeyIdx("");
    }
  }; 

  const updateZ = (index: number, val: number) => {
    const newArr = [...compounds];
    newArr[index].z = val;
    setCompounds(newArr);
  };

  // --- Core Calc: Antoine Pressure ---
  // ln(P_sat [Pa]) = A - B / (T [K] + C)  (Standard form for SI units usually found in these DBs)
  // NOTE: Check your specific DB units. Assuming Pa and K based on "7732-18-5" JSON example.
  // The JSON provided showed A=23.4, B=3987. This fits ln(P) = A - B/(T+C) for Pa.
  const getPsat = (c: Compound, T: number) => {
    return Math.exp(c.antoine.A - c.antoine.B / (T + c.antoine.C));
  };

  // --- Core Calc: Bubble Point Temperature ---
  // Find T where Sum(z_i * K_i) = 1  => Sum(z_i * P_sat_i / P) = 1
  const solveBubblePoint = (comps: Compound[], P: number) => {
    let T = 350; // Initial Guess (K)
    let error = 1.0;
    let iter = 0;
    
    while (Math.abs(error) > 1e-4 && iter < 50) {
        let K_sum = 0;
        let dSum_dT = 0; // Derivative for Newton Raphson

        comps.forEach(c => {
            const Psat = getPsat(c, T);
            const K = Psat / P;
            K_sum += c.z * K;
            
            // Derivative of Psat w.r.t T: dP/dT = P * [ B / (T+C)^2 ]
            const dP_dT = Psat * (c.antoine.B / Math.pow(T + c.antoine.C, 2));
            const dK_dT = dP_dT / P;
            dSum_dT += c.z * dK_dT;
        });

        error = K_sum - 1;
        // Newton Step: T_new = T_old - f(T)/f'(T)
        if (dSum_dT !== 0) {
            T = T - error / dSum_dT;
        } else {
            T += 1; // Fallback
        }
        iter++;
    }
    return T;
  };

  // --- FUG Solver ---
  const calculateFUG = () => {
    // Require at least 3 compounds to match UI and ensure meaningful calculation
    if (compounds.length < 3 || !lightKeyIdx || !heavyKeyIdx) return;

    // CONVERT BAR TO PASCALS FOR MATH
    const P_Pa = pressure * 100000;

    // 1. Normalize Composition
    const totalZ = compounds.reduce((sum, c) => sum + c.z, 0);
    const normalizedComps = compounds.map(c => ({ ...c, z: c.z / totalZ }));

    // 2. Determine Bubble Point of Feed to get average Alpha
    const T_feed = solveBubblePoint(normalizedComps, P_Pa);

    // 3. Calculate Alphas relative to Heavy Key
    const hkIndex = compounds.findIndex(c => c.id === parseInt(heavyKeyIdx));
    const lkIndex = compounds.findIndex(c => c.id === parseInt(lightKeyIdx));

    if (hkIndex < 0 || lkIndex < 0) {
      alert('Invalid key selection. Please select valid Light and Heavy Keys.');
      return;
    }

    const P_hk = getPsat(normalizedComps[hkIndex], T_feed);
    const K_hk = P_hk / P_Pa;

    const compsWithAlpha = normalizedComps.map(c => {
      const Psat = getPsat(c, T_feed);
      const K = Psat / P_Pa;
      return { ...c, K, alpha: K / K_hk };
    });

    const alphaLK = compsWithAlpha[lkIndex].alpha!;

    // --- CRITICAL VALIDATION ---
    // If the user picked a Light Key that is actually heavier than the Heavy Key, 
    // alphaLK will be < 1, and the math breaks.
    if (alphaLK <= 1.0001) {
      alert("Invalid Selection: The Light Key must be more volatile (higher Alpha) than the Heavy Key. Please swap your selection.");
      return;
    }

    // 4. Mass Balance (Distribute components)
    const FeedFlow = 100; // Basis: 100 kmol/hr
    let D_total = 0;

    const compsWithDist = compsWithAlpha.map((c, idx) => {
      let d_i = 0; // Flow in Distillate
      const f_i = c.z * FeedFlow;

      if (idx === lkIndex) {
        d_i = f_i * recLK;
      } else if (idx === hkIndex) {
        d_i = f_i * (1 - recHK);
      } else if (c.alpha! > alphaLK) {
        d_i = f_i * 0.999; // Lighters go up
      } else {
        d_i = f_i * 0.001; // Heavies go down
      }
      D_total += d_i;
      return { ...c, d_i, b_i: f_i - d_i };
    });

    const compsFinal = compsWithDist.map(c => ({
      ...c,
      xD: c.d_i / D_total,
      xB: c.b_i / (FeedFlow - D_total)
    }));

    // 5. Fenske Equation (Minimum Stages)
    const dLK = compsFinal[lkIndex].d_i;
    const bLK = compsFinal[lkIndex].b_i;
    const dHK = compsFinal[hkIndex].d_i;
    const bHK = compsFinal[hkIndex].b_i;

    const separationFactor = (dLK / bLK) * (bHK / dHK);
    const N_min = Math.log10(separationFactor) / Math.log10(alphaLK);

    // 6. Underwood Equations (Minimum Reflux)
    // Step A: Robust Theta Search (Bisection)
    // The root theta lies between alpha_HK (1.0) and alpha_LK.
    let low = 1.000001;
    let high = alphaLK - 0.000001;
    let theta = (low + high) / 2;
    const target = 1 - qValue;

    for (let i = 0; i < 100; i++) {
      let f = 0;
      compsFinal.forEach(c => {
        // Guard against division by zero if theta lands exactly on an alpha
        const denom = (c.alpha! - theta);
        if (Math.abs(denom) > 1e-9) {
          f += (c.alpha! * c.z) / denom;
        }
      });

      if (Math.abs(f - target) < 1e-6) break;

      // Since f(theta) is increasing between poles:
      // If f > target, we are too high -> move High down
      // If f < target, we are too low -> move Low up
      if (f > target) {
        high = theta;
      } else {
        low = theta;
      }
      theta = (low + high) / 2;
    }

    // Step B: Calculate Rmin
    let sumUnderwood = 0;
    compsFinal.forEach(c => {
      sumUnderwood += (c.alpha! * c.xD!) / (c.alpha! - theta);
    });
    const R_min = sumUnderwood - 1;

    // 7. Gilliland Correlation
    const R_actual = R_min * R_factor;
    const X = (R_actual - R_min) / (R_actual + 1);
    const Y = 0.75 * (1 - Math.pow(X, 0.5668));
    const N_actual = (N_min + Y) / (1 - Y);

    // 8. Kirkbride (Feed Stage)
    const B_total = FeedFlow - D_total;
    const zHK = compsFinal[hkIndex].z;
    const zLK = compsFinal[lkIndex].z;
    const xBLK = compsFinal[lkIndex].xB!;
    const xDHK = compsFinal[hkIndex].xD!;

    const kirkArg = (B_total / D_total) * (zHK / zLK) * Math.pow(xBLK / xDHK, 2);
    const logRatio = 0.206 * Math.log10(kirkArg);
    const ratio = Math.pow(10, logRatio);

    const N_strip = N_actual / (1 + ratio);
    const N_rect = N_actual - N_strip;

    const curveData = [];
    for (let x = 0; x <= 1; x += 0.01) {
      const y = 0.75 * (1 - Math.pow(x, 0.5668));
      curveData.push([x, y]);
    }

    setResults({
      N_min,
      R_min,
      R_actual,
      N_actual,
      FeedStage: N_rect,
      FeedTemp: T_feed,
      GillilandData: curveData,
      OperatingPoint: [X, Y]
    });
  };

  // Live update when the user changes the reflux multiplier slider: update only
  // the derived operating point and stage count so the plot updates in real-time.
  useEffect(() => {
    // Only depend on R_factor and update results safely using a functional update
    // to avoid an infinite loop caused by writing to `results` and then having
    // `results` trigger the effect again.
    if (!results) return;
    setResults(prev => {
      if (!prev || !prev.N_actual) return prev;
      const R_actual = prev.R_min * R_factor;
      const X = (R_actual - prev.R_min) / (R_actual + 1);
      const Y = 0.75 * (1 - Math.pow(X, 0.5668));
      const N_actual = (prev.N_min + Y) / (1 - Y);

      const eps = 1e-9;
      const sameR = Math.abs((prev.R_actual ?? 0) - R_actual) < eps;
      const sameN = Math.abs((prev.N_actual ?? 0) - N_actual) < eps;
      const sameX = Math.abs((prev.OperatingPoint?.[0] ?? 0) - X) < eps;
      const sameY = Math.abs((prev.OperatingPoint?.[1] ?? 0) - Y) < eps;

      // Scale FeedStage proportionally to change in N_actual (keep ratio rect/strip constant)
      const scale = prev.N_actual ? (N_actual / prev.N_actual) : 1;
      const newFeedStage = (prev.FeedStage ?? 0) * scale;
      const sameFeed = Math.abs((prev.FeedStage ?? 0) - newFeedStage) < eps;

      if (sameR && sameN && sameX && sameY && sameFeed) return prev;
      return { ...prev, R_actual, N_actual, OperatingPoint: [X, Y], FeedStage: newFeedStage };
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [R_factor]);

  // Auto-run full FUG calculation when core inputs change (run immediately to keep red dot responsive)
  useEffect(() => {
    // Require keys and at least 3 components
    if (compounds.length < 3 || !lightKeyIdx || !heavyKeyIdx) return;

    // Run calculation immediately on any input change for instant visual feedback
    calculateFUG();

    // No cleanup necessary
  }, [pressure, qValue, recLK, recHK, compounds, lightKeyIdx, heavyKeyIdx]);

  // --- ECharts Options ---
  const gillilandOption = results ? {
    backgroundColor: 'transparent',
    // Global text style to enforce Merriweather Sans across all chart text
    textStyle: { fontFamily: 'Merriweather Sans, sans-serif', color: textColor },
    title: { 
        text: 'Gilliland Correlation', 
        left: 'center', 
        textStyle: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans, sans-serif' } 
    },
    tooltip: {
        trigger: 'axis',
        renderMode: 'html',
        backgroundColor: resolvedTheme === 'dark' ? '#2563eb' : 'rgba(255,255,255,0.95)',
        borderColor: resolvedTheme === 'dark' ? '#2563eb' : '#ddd',
        padding: 8,
        textStyle: { fontFamily: 'Merriweather Sans, sans-serif', color: resolvedTheme === 'dark' ? '#ffffff' : '#000000' },
        formatter: (params: any) => {
            if (!params || !params.length) return '';
            const xVal = params[0].axisValue;
            const color = resolvedTheme === 'dark' ? '#ffffff' : '#000000';
            let html = `<div style="color:${color}; font-family:'Merriweather Sans', 'Merriweather Sans', sans-serif; font-size:14px; font-weight:400; line-height:1.2; margin:0; padding:6px 4px;">Reflux Param (X): <b>${Number(xVal).toPrecision(3)}</b><br/>`;
            params.forEach((p: any) => {
                if (p.value && p.value.length === 2) {
                    const yVal = p.value[1];
                    html += `${p.seriesName}: <b>${Number(yVal).toPrecision(3)}</b><br/>`;
                }
            });
            html += '</div>';
            return html;
        }
    },
    xAxis: { 
        name: '(R-Rmin)/(R+1)', 
        min: 0, max: 1, 
        nameLocation: 'middle', 
        nameGap: 30,
        axisLabel: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' },
        nameTextStyle: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' },
        axisLine: { lineStyle: { color: textColor } },
        splitLine: { show: false }
    },
    yAxis: { 
        name: '(N-Nmin)/(N+1)', 
        min: 0, max: 1,
        nameLocation: 'middle',
        nameRotate: 90,
        nameGap: 40,
        axisLabel: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' },
        nameTextStyle: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' },
        axisLine: { lineStyle: { color: textColor } },
        splitLine: { show: false }
    },
    grid: { right: '10%', top: '12%', bottom: '10%', left: '10%' },
    series: [
      {
        name: 'Correlation',
        type: 'line',
        z: 1,
        zlevel: 1,
        showSymbol: false,
        data: results.GillilandData,
        smooth: true,
        lineStyle: { width: 3, color: '#3b82f6' }
      },
      {
        name: 'Design Point',
        type: 'scatter',
        z: 10,
        zlevel: 10,
        data: [results.OperatingPoint],
        itemStyle: { color: '#ef4444' },
        symbolSize: 9,
        tooltip: {
             renderMode: 'html',
             backgroundColor: resolvedTheme === 'dark' ? '#2563eb' : 'rgba(255,255,255,0.95)',
             borderColor: resolvedTheme === 'dark' ? '#2563eb' : '#ddd',
             padding: 8,
             textStyle: { fontFamily: 'Merriweather Sans, sans-serif', color: resolvedTheme === 'dark' ? '#ffffff' : '#000000' },
             formatter: (p: any) => {
                 const color = resolvedTheme === 'dark' ? '#ffffff' : '#000000';
                 return `<div style="color:${color}; font-family:'Merriweather Sans', 'Merriweather Sans', sans-serif; font-size:14px; font-weight:500; line-height:1.2; margin:0; padding:6px 4px;">Design Point<br/>X: ${p.value[0].toPrecision(3)}<br/>Y: ${p.value[1].toPrecision(3)}</div>`;
             }
        }
      }
    ]
  } : {};  

  return (
    <TooltipProvider>
      <div className="p-8 space-y-8 max-w-7xl mx-auto">
        <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
            
          {/* LEFT PANEL: INPUTS */}
          <div className="lg:col-span-4 space-y-6">
            
            {/* 1. Add Compounds */}
            <Card>
              <CardHeader className="pb-3">
                <CardTitle className="text-lg font-medium flex items-center">
                   <Search className="w-4 h-4 mr-2" /> Add Components
                </CardTitle>
              </CardHeader>
              <CardContent>
                {/* FIX: Removed the Button square, kept only Input */}
                <div className="relative mb-4">
                    <Input 
                        placeholder="Search chemical..." 
                        value={searchQuery}
                        onChange={handleSearchChange}
                        onFocus={() => { if(suggestions.length > 0) setShowSuggestions(true); }}
                        // FIX: Enter key now triggers auto-select logic in handleAddCompound
                        onKeyDown={(e) => e.key === 'Enter' && handleAddCompound()}
                        autoComplete="off"
                    />
                    {showSuggestions && suggestions.length > 0 && (
                        <div 
                            ref={suggestionsRef} 
                            className="absolute z-20 w-full bg-background border border-input rounded-md shadow-lg mt-1 top-full max-h-48 overflow-y-auto"
                        >
                            {suggestions.map((s, i) => (
                                <div 
                                    key={i} 
                                    onClick={() => handleSuggestionClick(s)} 
                                    className="px-3 py-2 hover:bg-accent cursor-pointer text-sm"
                                >
                                    {s}
                                </div>
                            ))}
                        </div>
                    )}
                </div>
                
                <div className="space-y-3 max-h-[300px] overflow-y-auto pr-2">
                    {compounds.map((c, i) => (
                        <div key={i} className="flex items-center space-x-4 text-sm bg-secondary/50 p-3 rounded-md">
                            <div className="w-24 font-semibold truncate" title={c.name}>
                                {c.name}
                            </div>
                            
                            <div className="flex-1 flex flex-col justify-center">
                                <div className="flex justify-between text-[10px] text-muted-foreground mb-1">
                                    <span>Mole Frac</span>
                                    <span>{c.z.toFixed(2)}</span>
                                </div>
                                <input 
                                    type="range" 
                                    min="0.01" max="0.99" step="0.01"
                                    value={c.z} 
                                    onChange={(e) => updateZ(i, parseFloat(e.target.value))}
                                    className="w-full h-2 bg-secondary rounded-lg appearance-none cursor-pointer accent-primary"
                                />
                            </div>

                            <Button variant="ghost" size="sm" onClick={() => removeCompound(i)}>
                                <Trash2 className="w-4 h-4 text-red-400" />
                            </Button> 
                        </div>
                    ))}
                    {compounds.length === 0 && <div className="text-sm text-center text-muted-foreground py-4">No components added.</div>}
                </div>
                {/* FIX: Removed "Total Composition" text div here */}
              </CardContent>
            </Card>

            {/* 2. Process Specs */}
            <Card>
                <CardHeader className="pb-3">
                    <CardTitle className="text-lg font-medium">Design Specifications</CardTitle>
                </CardHeader>
                <CardContent className="space-y-6">
                    {/* Key Selection (Dropdowns) */}
                    <div className="grid grid-cols-2 gap-4">
                        <div className="space-y-1">
                            <Label>Light Key (LK)</Label>
                            <Select value={lightKeyIdx} onValueChange={setLightKeyIdx}>
                                <SelectTrigger><SelectValue placeholder="Select LK" /></SelectTrigger>
                                <SelectContent>
                                            {compounds.map((c) => {
                                        let isDisabled = false;
                                        const hkId = heavyKeyIdx ? parseInt(heavyKeyIdx) : NaN;
                                        if (!Number.isNaN(hkId)) {
                                            const hkC = compounds.find(cc => cc.id === hkId);
                                            if (hkC) {
                                                if (c.id === hkId) isDisabled = true; // Can't be same
                                                if (c.Tb >= hkC.Tb) isDisabled = true; // Can't be heavier than HK
                                            }
                                        }
                                        if (isDisabled) return null;

                                        return <SelectItem key={c.id} value={c.id.toString()}>{c.name}</SelectItem>;
                                    })}
                                </SelectContent>
                            </Select>
                        </div>
                        <div className="space-y-1">
                            <Label>Heavy Key (HK)</Label>
                            <Select value={heavyKeyIdx} onValueChange={setHeavyKeyIdx}>
                                <SelectTrigger><SelectValue placeholder="Select HK" /></SelectTrigger>
                                <SelectContent>
                                    {compounds.map((c) => {
                                        let isDisabled = false;
                                        const lkId = lightKeyIdx ? parseInt(lightKeyIdx) : NaN;
                                        if (!Number.isNaN(lkId)) {
                                            const lkC = compounds.find(cc => cc.id === lkId);
                                            if (lkC) {
                                                if (c.id === lkId) isDisabled = true;
                                                if (c.Tb <= lkC.Tb) isDisabled = true; // Can't be lighter than LK
                                            }
                                        }
                                        if (isDisabled) return null;

                                        return <SelectItem key={c.id} value={c.id.toString()}>{c.name}</SelectItem>;
                                    })}
                                </SelectContent>
                            </Select>
                        </div>
                    </div>

                    {/* Recovery Sliders */}
                    <div className="space-y-4 border-t pt-4">
                        <div className="space-y-2">
                            <div className="flex justify-between">
                                <Label>LK Recovery (Distillate)</Label>
                                <span className="text-xs font-mono text-muted-foreground">{(recLK * 100).toFixed(1)}%</span>
                            </div>
                            <input 
                                type="range" min="0.50" max="1.00" step="0.001"
                                value={recLK} 
                                onChange={e=>setRecLK(parseFloat(e.target.value))}
                                className="w-full h-2 bg-secondary rounded-lg appearance-none cursor-pointer accent-blue-500"
                            />
                        </div>
                        <div className="space-y-2">
                            <div className="flex justify-between">
                                <Label>HK Recovery (Bottoms)</Label>
                                <span className="text-xs font-mono text-muted-foreground">{(recHK * 100).toFixed(1)}%</span>
                            </div>
                            <input 
                                type="range" min="0.50" max="1.00" step="0.001"
                                value={recHK} 
                                onChange={e=>setRecHK(parseFloat(e.target.value))}
                                className="w-full h-2 bg-secondary rounded-lg appearance-none cursor-pointer accent-blue-500"
                            />
                        </div>
                    </div>

                    {/* Pressure & Feed Q */}
                    <div className="grid grid-cols-2 gap-4 border-t pt-4">
                        <div className="space-y-1">
                            <Label>Pressure (Bar)</Label>
                            <Input 
                                type="number" 
                                value={pressure} 
                                onChange={e=>setPressure(parseFloat(e.target.value))} 
                                step={0.1}
                            />
                        </div>
                        <div className="space-y-1">
                            <Label>Feed Quality (q)</Label>
                            <Input 
                                type="number" 
                                value={qValue} 
                                onChange={e=>setQValue(parseFloat(e.target.value))} 
                                step={0.1} 
                                placeholder="1.0 = Liquid"
                            />
                        </div>
                    </div>

                    {/* Reflux Multiplier */}
                    <div className="space-y-3 pt-2 border-t">
                        <div className="flex justify-between">
                            <Label>Reflux Multiplier (R / Rmin)</Label>
                            <span className="text-sm font-mono text-blue-400">{R_factor}x</span>
                        </div>
                        <input 
                            type="range" 
                            min="1.05" max="3.0" step="0.05" 
                            value={R_factor} 
                            onChange={(e) => setRFactor(parseFloat(e.target.value))}
                            className="w-full h-2 bg-secondary rounded-lg appearance-none cursor-pointer accent-green-500"
                        />
                    </div>

                    <Button className="w-full" onClick={calculateFUG} disabled={compounds.length < 3}>
                        <Calculator className="w-4 h-4 mr-2" /> Calculate Column
                    </Button>
                </CardContent>
            </Card>
          </div>

          {/* RIGHT PANEL: RESULTS */}
          <div className="lg:col-span-8 space-y-6">
            {results ? (
                <>
                 <Card className="h-[560px]">
                    <CardContent className="h-[520px]">
                        <ReactECharts 
                            ref={echartsRef} 
                            option={gillilandOption} 
                            style={{ height: '100%', width: '100%' }} 
                        />
                    </CardContent>
                 </Card>

                 <div className="grid grid-cols-1 md:grid-cols-5 gap-4">
                    <Card className="bg-primary/5 border-primary/20">
                        <CardContent className="p-4 text-center">
                            <div className="text-sm text-muted-foreground">Minimum Stages</div>
                            <div className="text-2xl font-bold">{results.N_min.toFixed(1)}</div>
                            <div className="text-xs text-muted-foreground">Fenske Eq</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-primary/5 border-primary/20">
                        <CardContent className="p-4 text-center">
                            <div className="text-sm text-muted-foreground">Min Reflux Ratio</div>
                            <div className="text-2xl font-bold">{results.R_min.toFixed(2)}</div>
                            <div className="text-xs text-muted-foreground">Underwood Eq</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-primary/5 border-primary/20">
                        <CardContent className="p-4 text-center">
                            <div className="text-sm text-muted-foreground">Feed Bubble Pt</div>
                            {/* Converted to Celsius for readability */}
                            <div className="text-2xl font-bold">{(results.FeedTemp - 273.15).toFixed(1)} °C</div>
                            <div className="text-xs text-muted-foreground">{(results.FeedTemp).toFixed(1)} K</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-blue-500/10 border-blue-500/20">
                        <CardContent className="p-4 text-center">
                            <div className="text-sm text-muted-foreground">Actual Stages</div>
                            <div className="text-3xl font-bold text-blue-500">{results.N_actual.toFixed(1)}</div>
                            <div className="text-xs text-muted-foreground">@ R = {results.R_actual.toFixed(2)}</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-blue-500/10 border-blue-500/20">
                        <CardContent className="p-4 text-center">
                            <div className="text-sm text-muted-foreground">Feed Location</div>
                            <div className="text-3xl font-bold text-blue-500">{results.FeedStage.toFixed(1)}</div>
                            <div className="text-xs text-muted-foreground">Stages from Top</div>
                        </CardContent>
                    </Card>
                 </div> 
                </>
            ) : (
                <div className="h-full flex flex-col items-center justify-center text-muted-foreground border-2 border-dashed rounded-lg p-12">
                    <Calculator className="w-12 h-12 mb-4 opacity-50" />
                    <p>Add components and set keys to run the FUG design simulation.</p>
                </div>
            )}
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}