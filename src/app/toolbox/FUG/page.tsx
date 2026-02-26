'use client';

import React, { useState, useEffect, useRef, useLayoutEffect } from 'react';
import { supabase } from '@/lib/supabaseClient';
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
import { Tabs, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip as ShadTooltip, TooltipContent, TooltipTrigger } from "@/components/ui/tooltip";
import { TooltipProvider } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Trash2, Calculator, Search, Lock, Unlock } from 'lucide-react';

// --- Types ---
type Compound = {
  id: number;
  name: string; // "CompoundID.value"
  formula: string; 
  antoine: { A: number; B: number; C: number }; // Parsed from "AntoineVaporPressure"
  z: number; // Mole fraction (user input)
  Tb: number; // Normal Boiling Point (K) for sorting/filtering
  locked?: boolean; // UI lock to hold mole fraction fixed
  alpha?: number; // Relative volatility (calculated)
  K?: number; // K-value (calculated)
  xD?: number; // Distillate composition (calculated)
  xB?: number; // Bottoms composition (calculated)
};

// Add this type for detailed flash results (vapor fraction + phase compositions)
type FlashResult = {
  psi: number;          // Vapor Fraction
  stateMsg: string;     // "Subcooled", "Two-Phase", etc.
  phases: { 
    name: string; 
    z: number; // Feed
    x: number; // Liquid
    y: number; // Vapor
  }[];
};

// Helper: smart number formatting (up to 3 significant figures, remove trailing zeros)
const formatNum = (num: number) => {
  if (!Number.isFinite(num) || Math.abs(num) < 1e-6) return "0";
  return parseFloat(num.toPrecision(3)).toString();
};

// Slider-specific formatting: show up to 2 decimal places; if smaller than 0.01 display '0'
const formatSliderValue = (num: number) => {
  if (!Number.isFinite(num) || Math.abs(num) < 0.01) return '0';
  // Show up to 2 decimals but strip trailing zeros (e.g., 0.10 -> '0.1')
  const s = num.toFixed(2).replace(/\.?(0+)$/,'');
  return s;
};

// Phase fraction formatting for flash table: show up to thousandths (3 decimal places).
// Any value with magnitude smaller than 0.001 is shown as '0' to avoid noisy sub-thousandth precision.
const formatPhaseFraction = (num: number) => {
  if (!Number.isFinite(num) || Math.abs(num) < 1e-6) return '0';
  if (Math.abs(num) < 0.001) return '0';
  // Round to 3 decimal places and remove unnecessary trailing zeros
  const rounded = Math.round(num * 1000) / 1000;
  return parseFloat(rounded.toFixed(3)).toString();
};

// Fixed 3-decimal formatter for phase bar chart (always show 3 decimal places unless <0.001 -> '0')
const formatPhaseFractionFixed3 = (num: number) => {
  if (!Number.isFinite(num) || Math.abs(num) < 1e-6) return '0';
  if (Math.abs(num) < 0.001) return '0';
  return num.toFixed(3);
};

// Small helper component: shrink text to fit the container without ellipses
function ShrinkToFit({ text }: { text: string }) {
  const ref = React.useRef<HTMLDivElement | null>(null);
  const [fontSize, setFontSize] = React.useState<number>(14);
  useLayoutEffect(() => {
    const el = ref.current;
    if (!el) return;
    // Start from a reasonable max and shrink until it fits or min font reached
    let fs = 14;
    const minFs = 10;
    el.style.fontSize = fs + 'px';
    // Measure and shrink
    while (el.scrollWidth > el.clientWidth && fs > minFs) {
      fs -= 1;
      el.style.fontSize = fs + 'px';
    }
    setFontSize(fs);
  }, [text]);

  return (
    <div ref={ref} style={{ fontSize: `${fontSize}px`, lineHeight: 1 }} className="whitespace-nowrap overflow-hidden">
      {text}
    </div>
  );
}

// --- FUG Simulation Component ---
export default function FUGSimulation() {
  const { resolvedTheme } = useTheme(); // Use resolvedTheme for accurate dark/light detection
  const textColor = resolvedTheme === 'dark' ? '#ffffff' : '#000000'; // Dynamic text color
  const [showFlashBarChart, setShowFlashBarChart] = useState<boolean>(false);


  const echartsRef = useRef<ReactECharts>(null);
  const [_loading, setLoading] = useState(false);
  const [searchQuery, setSearchQuery] = useState('');
  
  // Simulation State
  const [compounds, setCompounds] = useState<Compound[]>([]);
  const [pressure, setPressure] = useState<number>(1); // Default to 1 bar
  const [qValue, setQValue] = useState<number>(1.0); // Feed Quality (1=Liquid, 0=Vapor)
  const [inputTemp, setInputTemp] = useState<number>(350); // Default Kelvin
  // Feed Spec mode: 'temp' = isothermal feed spec, 'q' = quality-specified flash
  const [feedSpecMode, setFeedSpecMode] = useState<'temp'|'q'>('temp');
  // Controlled string inputs so the user can clear the field while we keep the last valid numeric value
  const [pressureInput, setPressureInput] = useState<string>(pressure.toString());
  const [qInput, setQInput] = useState<string>(qValue.toString());
  const [inputTempInput, setInputTempInput] = useState<string>(inputTemp.toString());

  // Handler: keeps the visible string in sync while only updating numeric state when value is a valid number
  // Added optional `minAllowed` to prevent values <= minAllowed (useful to disallow 0 or negative inputs)
  const handleNumericInputString = (
    str: string,
    setStr: (s: string) => void,
    setNum: (n: number) => void,
    _prevNum: number,
    minAllowed?: number
  ) => {
    setStr(str);
    // Allow user to clear the box or type a lone '-' or '.' while editing without changing numeric state
    if (str === '' || str === '-' || str === '.' || str === '-.') return;
    const num = parseFloat(str);
    if (!isNaN(num)) {
      if (typeof minAllowed === 'number' && num <= minAllowed) {
        // Clamp to the minimum allowed and reflect immediately in both string and numeric state
        setNum(minAllowed);
        setStr(String(minAllowed));
      } else {
        setNum(num);
      }
    }
  };

  // Keep string views in sync when numeric values change programmatically
  useEffect(() => setPressureInput(pressure.toString()), [pressure]);
  useEffect(() => setQInput(qValue.toString()), [qValue]);
  // Keep inputTempInput user-friendly and avoid long floating-point noise that can trigger
  // the 'text-transparent' clipping. Round to 3 decimal places and trim trailing zeros.
  useEffect(() => setInputTempInput(() => {
    const v = Number(inputTemp);
    if (!Number.isFinite(v)) return '';
    return parseFloat(v.toFixed(3)).toString();
  }), [inputTemp]);

  const [flashState, setFlashState] = useState<{ msg: string; psi?: number } | null>(null);
  // Flash result details (x, y for each component)
  const [flashDetails, setFlashDetails] = useState<FlashResult | null>(null);
  const [lightKeyIdx, setLightKeyIdx] = useState<string>("");
  const [heavyKeyIdx, setHeavyKeyIdx] = useState<string>("");
  const [recLK, setRecLK] = useState<number>(0.99); // Recovery of LK in Distillate
  const [recHK, setRecHK] = useState<number>(0.99); // Recovery of HK in Bottoms
  const [R_factor, setRFactor] = useState<number>(1.2); // R / Rmin multiplier

  // Informational message when automatic key swaps occur
  const [keySwapMsg, setKeySwapMsg] = useState<string | null>(null);

  // Autocomplete suggestions for component search
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [showSuggestions, setShowSuggestions] = useState(false);
  const suggestionsRef = useRef<HTMLDivElement>(null);

  // Ref to ensure initial auto-calc runs only once on mount when defaults loaded
  const initialCalcRef = useRef(false);
  // Guard to prevent double-loading defaults (React Strict Mode / multiple mounts)
  const defaultsLoadedRef = useRef(false);
  // Track previous compound count to detect transitions (e.g., 2 -> 3)
  const prevCompCountRef = useRef<number>(compounds.length);
  // Suppress compound-change effect when we trigger immediate calculations (avoid duplicate runs)
  const suppressCompoundEffectRef = useRef(false);

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
      .limit(10); // Fetch a few more to allow for filtering
    
    if (data) {
      // Filter out compounds that are already added (case-insensitive)
      const existingNames = new Set(compounds.map(c => c.name.toLowerCase()));
      const filtered = data
        .map((d: any) => d.name)
        .filter((name: string) => !existingNames.has(name.toLowerCase()))
        .slice(0, 5); // Take top 5 unique

      setSuggestions(filtered);
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

  // --- Core Calc: Dew Point Temperature ---
  // Find T where Sum(z_i / K_i) = 1
  const solveDewPoint = (comps: Compound[], P: number) => {
    let T = 350; // Initial Guess (K)
    let error = 1.0;
    let iter = 0;

    while (Math.abs(error) > 1e-4 && iter < 50) {
      let invK_sum = 0;
      let dSum_dT = 0;

      comps.forEach(c => {
        const Psat = getPsat(c, T);
        const K = Psat / P;
        invK_sum += c.z / K;

        // Derivative d(1/K)/dT = -1/K^2 * dK/dT
        const dP_dT = Psat * (c.antoine.B / Math.pow(T + c.antoine.C, 2));
        const dK_dT = dP_dT / P;
        const dInvK_dT = -(1 / (K * K)) * dK_dT;

        dSum_dT += c.z * dInvK_dT;
      });

      error = invK_sum - 1;
      if (dSum_dT !== 0) {
        T = T - error / dSum_dT;
      } else {
        T += 1;
      }
      iter++;
    }
    return T;
  };

  // --- Core Calc: Rachford-Rice (Isothermal Flash) ---
  // Solve for psi (Vapor Fraction) where Sum( z_i * (K_i - 1) / (1 + psi*(K_i - 1)) ) = 0
  const solveRachfordRice = (comps: Compound[], P: number, T: number) => {
    // 1. Calculate K-values at current T
    const compsWithK = comps.map(c => {
      const Psat = getPsat(c, T);
      return { ...c, K: Psat / P };
    });

    // Newton-Raphson for Psi
    let psi = 0.5; // Guess 50% vapor
    let iter = 0;
    let error = 1.0;

    while (Math.abs(error) > 1e-6 && iter < 50) {
      let f = 0;
      let df = 0;

      compsWithK.forEach(c => {
        const K = c.K!;
        const num = c.z * (K - 1);
        const den = 1 + psi * (K - 1);

        f += num / den;
        // Derivative: - (z * (K-1)^2) / (1 + psi*(K-1))^2
        df -= (num * (K - 1)) / (den * den);
      });

      error = f;
      if (Math.abs(df) > 1e-9) {
        let newPsi = psi - f / df;
        // Dampen or clamp if it shoots out of bounds [0,1]
        if (newPsi < 0) newPsi = 0.0001;
        if (newPsi > 1) newPsi = 0.9999;
        psi = newPsi;
      } else {
        break;
      }
      iter++;
    }

    return { psi, compsWithK };
  };

  // --- Helper: Fetch Compound ---
  const handleAddCompound = async (specificName?: string): Promise<Compound | null> => {
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

        // Prevent adding duplicates (case-insensitive)
        if (compounds.some(c => c.name.toLowerCase() === row.name.toLowerCase())) {
            alert(`Compound "${row.name}" is already added.`);
            setLoading(false);
            return null;
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
            // Atomically insert and renormalize so mole fractions sum to 1; respect locked components
            // Use updater form to avoid race conditions when handleAddCompound is called concurrently
            let addedCompound: Compound | null = null;
            setCompounds(prev => {
              // Prevent duplicates atomically by checking against the current prev
              if (prev.some(c => c.id === row.id || c.name.toLowerCase() === row.name.toLowerCase())) {
                return prev;
              }

              const arr = prev.map(p => ({ ...p }));
              const newZ = newCompound.z;
              const lockedSum = arr.reduce((s, c) => s + (c.locked ? c.z : 0), 0);
              const unlocked = arr.filter(c => !c.locked);
              const unlockedSum = unlocked.reduce((s, c) => s + c.z, 0);

              // Ensure the new compound doesn't exceed available room
              const availableForNew = Math.max(0, 1 - lockedSum);
              const finalNewZ = Math.min(newZ, availableForNew);

              const remaining = Math.max(0, 1 - lockedSum - finalNewZ);

              if (unlockedSum > 0) {
                const scale = remaining / unlockedSum;
                for (let i = 0; i < arr.length; i++) {
                  if (!arr[i].locked) arr[i].z = arr[i].z * scale;
                }
              } else if (arr.length > 0) {
                const unlockedCount = arr.filter(c => !c.locked).length;
                if (unlockedCount > 0) {
                  const share = remaining / unlockedCount;
                  for (let i = 0; i < arr.length; i++) if (!arr[i].locked) arr[i].z = share;
                }
              }

              const compiled = { ...newCompound, z: finalNewZ };
              addedCompound = compiled;
              return [...arr, compiled];
            });

            // If another concurrent call already added the same compound, alert and abort
            if (!addedCompound) {
              alert(`Compound "${row.name}" is already added.`);
              setLoading(false);
              return null;
            }

            setSearchQuery('');
            setSuggestions([]); // Clear suggestions
            setLoading(false);
            return addedCompound;
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
        // Prevent double loading in StrictMode or if already loaded
        if (defaultsLoadedRef.current) return;
        defaultsLoadedRef.current = true;
        if (compounds.length > 0) return;
        // Load sequentially to maintain order and capture IDs
        const m = await handleAddCompound("Methanol");
        await handleAddCompound("Ethanol");
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

  // If user just added a third component (2 -> 3), auto-select keys if missing and run calculation
  useEffect(() => {
    const prev = prevCompCountRef.current;
    const curr = compounds.length;
    // Only act on a transition from <3 to >=3
    if (prev < 3 && curr >= 3) {
      // If keys are not set, pick first and last as sensible defaults
      if (!lightKeyIdx && compounds.length >= 1) {
        setLightKeyIdx(String(compounds[0].id));
      }
      if (!heavyKeyIdx && compounds.length >= 1) {
        setHeavyKeyIdx(String(compounds[compounds.length - 1].id));
      }
      // Give React a tick to apply key changes and then calculate
      setTimeout(() => {
        if (compounds.length >= 3) calculateFUG();
      }, 40);
    }
    prevCompCountRef.current = curr;
  }, [compounds.length, lightKeyIdx, heavyKeyIdx]);

  const removeCompound = (index: number) => {
    const newArr = [...compounds];
    const removed = newArr.splice(index, 1)[0];

    // Re-normalize unlocked components so total mole fraction = 1 (locked ones hold their values)
    const lockedSum = newArr.reduce((s, c) => s + (c.locked ? c.z : 0), 0);
    const unlocked = newArr.filter(c => !c.locked);
    const unlockedSum = unlocked.reduce((s, c) => s + c.z, 0);
    const remaining = Math.max(0, 1 - lockedSum);
    if (unlockedSum > 0) {
      const scale = remaining / unlockedSum;
      for (let i = 0; i < newArr.length; i++) {
        if (!newArr[i].locked) newArr[i].z = newArr[i].z * scale;
      }
    } else {
      // If no unlocked components, distribute the remaining equally (rare)
      const count = newArr.length - (newArr.filter(c => c.locked).length);
      if (count > 0) {
        const share = remaining / count;
        for (let i = 0; i < newArr.length; i++) if (!newArr[i].locked) newArr[i].z = share;
      }
    }

    setCompounds(newArr);

    // If removed compound was selected as a key, clear that selection
    if (removed) {
      const removedId = removed.id;
      if (lightKeyIdx && parseInt(lightKeyIdx) === removedId) setLightKeyIdx("");
      if (heavyKeyIdx && parseInt(heavyKeyIdx) === removedId) setHeavyKeyIdx("");
    }
  }; 

  const updateZ = (index: number, val: number) => {
    // Only apply updates to unlocked components; sliders are disabled for locked ones
    const clamped = Math.max(0, val);
    const newArr = compounds.map((c) => ({ ...c }));
    // If the target is locked, ignore
    if (newArr[index].locked) return;

    // Calculate sums
    const lockedSum = newArr.reduce((s, c) => s + (c.locked ? c.z : 0), 0);
    const targetNew = Math.min(clamped, Math.max(0, 1 - lockedSum));

    // Sum of unlocked others (excluding index)
    const otherUnlocked = newArr.map((c, i) => ({ c, i })).filter(item => !item.c.locked && item.i !== index);
    const sumOtherUnlocked = otherUnlocked.reduce((s, item) => s + item.c.z, 0);

    const remaining = Math.max(0, 1 - lockedSum - targetNew);

    if (sumOtherUnlocked > 0) {
      const scale = remaining / sumOtherUnlocked;
      for (const item of otherUnlocked) {
        newArr[item.i].z = newArr[item.i].z * scale;
      }
    } else {
      // Distribute remaining equally among other unlocked (if any)
      const count = otherUnlocked.length;
      if (count > 0) {
        const share = remaining / count;
        for (const item of otherUnlocked) newArr[item.i].z = share;
      }
    }

    // Set target value
    newArr[index].z = targetNew;

    // Immediately reflect composition change and trigger calculations for instant feedback
    suppressCompoundEffectRef.current = true;
    setCompounds(newArr);

    // Run calculations using the new array so the graphs update without waiting for React state
    calculateFUG(newArr);
    if (feedSpecMode === 'temp') runFlashCalc(newArr);
    else computeFlashFromQ(parseFloat(qInput), newArr);
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

  const runFlashCalc = (compsOverride?: Compound[], tempOverride?: number) => {
    const comps = compsOverride ?? compounds;
    const T_use = (typeof tempOverride === 'number') ? tempOverride : inputTemp;
    if (comps.length === 0) return;

    const P_Pa = pressure * 100000;
    
    // 1. Check Boundaries
    const T_bub = solveBubblePoint(comps, P_Pa);
    const T_dew = solveDewPoint(comps, P_Pa);

    let calculatedQ = 0;
    let psi = 0;
    let stateMsg = "";
    let currentComps = comps.map(c => ({...c, K: getPsat(c, T_use)/P_Pa}));

    if (T_use < T_bub) {
      stateMsg = "Subcooled Liquid";
      psi = 0;
      calculatedQ = 1.05; // Just an indicator for "Liquid"
    } else if (T_use > T_dew) {
      stateMsg = "Superheated Vapor";
      psi = 1;
      calculatedQ = -0.05; // Indicator for "Vapor"
    } else {
      stateMsg = "Two-Phase Mixture";
      const rr = solveRachfordRice(comps, P_Pa, T_use);
      psi = rr.psi;
      currentComps = rr.compsWithK;
      calculatedQ = 1 - psi;
    }

    // 2. Calculate Phase Compositions (x and y)
    // x_i = z_i / (1 + psi(K_i - 1))
    // y_i = K_i * x_i
    const phases = currentComps.map((c: any) => {
      const denom = 1 + psi * (c.K! - 1);
      const x = c.z / denom;
      const y = c.K! * x;
      return { name: c.name, z: c.z, x, y };
    });

    setQValue(parseFloat(calculatedQ.toFixed(3)));
    setFlashDetails({ psi, stateMsg, phases });
    setFlashState({ msg: stateMsg, psi });
  };

  // Helper: find feed temperature that yields desired quality Q via bisection between bubble and dew
  const findTemperatureForQuality = (desiredQ: number, P_Pa: number, compsOverride?: Compound[]) => {
    // desiredQ: 1 = liquid, 0 = vapor
    const comps = compsOverride ?? compounds;
    const desiredPsi = 1 - desiredQ;
    const T_bub = solveBubblePoint(comps, P_Pa);
    const T_dew = solveDewPoint(comps, P_Pa);
    if (T_bub === null || T_dew === null) return { success: false, msg: 'Unable to compute bubble/dew bounds at this pressure.' };

    // Handle boundary qualities quickly
    if (desiredQ >= 1) return { success: true, T: T_bub, psi: 0, stateMsg: 'Subcooled Liquid' };
    if (desiredQ <= 0) return { success: true, T: T_dew, psi: 1, stateMsg: 'Superheated Vapor' };

    let low = T_bub, high = T_dew;

    for (let i = 0; i < 60; i++) {
      const mid = 0.5 * (low + high);
      try {
        const rr = solveRachfordRice(comps, P_Pa, mid);
        const psiMid = rr.psi;
        if (!Number.isFinite(psiMid)) break;
        if (Math.abs(psiMid - desiredPsi) < 1e-6) return { success: true, T: mid, psi: psiMid, compsWithK: rr.compsWithK };
        if (psiMid < desiredPsi) low = mid; else high = mid;
      } catch (err) {
        break;
      }
    }
    return { success: false, msg: `Unable to converge to the requested quality at this pressure.` };
  };

  const computeFlashFromQ = (desiredQ: number, compsOverride?: Compound[]) => {
    setFlashState(null);
    const comps = compsOverride ?? compounds;
    if (comps.length === 0) return;
    if (!Number.isFinite(desiredQ)) {
      setFlashState({ msg: 'Quality must be a number.' });
      setTimeout(() => setFlashState(null), 4000);
      return;
    }
    const P_Pa = pressure * 100000;
    const res = findTemperatureForQuality(desiredQ, P_Pa, comps);
    if (!res.success) {
      setFlashState({ msg: res.msg ?? 'Quality computation failed.' });
      setTimeout(() => setFlashState(null), 5000);
      return;
    }

    const T_found = res.T as number;
    const psi = typeof res.psi === 'number' ? res.psi : (1 - desiredQ);
    // Compute K-values and phases similar to runFlashCalc
    let currentComps = comps.map(c => ({ ...c, K: getPsat(c, T_found) / P_Pa }));
    if (psi === 0) {
      const phases = currentComps.map(c => ({ name: c.name, z: c.z, x: c.z, y: c.K! * c.z }));
      setInputTemp(T_found);
      setQValue(parseFloat(desiredQ.toFixed(3)));
      setFlashDetails({ psi: 0, stateMsg: 'Subcooled Liquid', phases });
      setFlashState({ msg: 'Subcooled Liquid', psi: 0 });
      return;
    }
    if (psi === 1) {
      const phases = currentComps.map(c => ({ name: c.name, z: c.z, x: c.z / (c.K! || 1), y: c.z }));
      setInputTemp(T_found);
      setQValue(parseFloat(desiredQ.toFixed(3)));
      setFlashDetails({ psi: 1, stateMsg: 'Superheated Vapor', phases });
      setFlashState({ msg: 'Superheated Vapor', psi: 1 });
      return;
    }

    // Two-phase: use compsWithK from the bisection result when available
    const compsWithK = (res as any).compsWithK ?? currentComps;
    const phases = compsWithK.map((c: any) => {
      const denom = 1 + psi * (c.K! - 1);
      const x = c.z / denom;
      const y = c.K! * x;
      return { name: c.name, z: c.z, x, y };
    });
    setInputTemp(T_found);
    setQValue(parseFloat(desiredQ.toFixed(3)));
    setFlashDetails({ psi, stateMsg: 'Two-Phase Mixture', phases });
    setFlashState({ msg: 'Two-Phase Mixture', psi });
  };

  // --- FUG Solver ---
  const calculateFUG = (compsOverride?: Compound[]) => {
    const comps = compsOverride ?? compounds;
    // Require at least 3 compounds to match UI and ensure meaningful calculation
    if (comps.length < 3 || !lightKeyIdx || !heavyKeyIdx) return;

    // CONVERT BAR TO PASCALS FOR MATH
    const P_Pa = pressure * 100000;

    // 1. Normalize Composition
    const totalZ = comps.reduce((sum, c) => sum + c.z, 0);
    const normalizedComps = comps.map(c => ({ ...c, z: c.z / totalZ }));

    // 2. Determine Bubble Point of Feed to get average Alpha
    const T_feed = solveBubblePoint(normalizedComps, P_Pa);

    // 3. Calculate Alphas relative to Heavy Key
    let hkIndex = compounds.findIndex(c => c.id === parseInt(heavyKeyIdx));
    let lkIndex = compounds.findIndex(c => c.id === parseInt(lightKeyIdx));

    if (hkIndex < 0 || lkIndex < 0) {
      alert('Invalid key selection. Please select valid Light and Heavy Keys.');
      return;
    }

    const P_hk = getPsat(normalizedComps[hkIndex], T_feed);
    const K_hk = P_hk / P_Pa;

    let compsWithAlpha = normalizedComps.map(c => {
      const Psat = getPsat(c, T_feed);
      const K = Psat / P_Pa;
      return { ...c, K, alpha: K / K_hk };
    });

    let alphaLK = compsWithAlpha[lkIndex].alpha!;

    // --- CRITICAL VALIDATION ---
    // If the user picked a Light Key that is actually heavier than the Heavy Key,
    // alphaLK will be <= 1. Automatically swap the keys and proceed if that fixes the issue.
    if (alphaLK <= 1.0001) {
      // Auto-swap light/heavy
      const oldLK = lightKeyIdx;
      const oldHK = heavyKeyIdx;
      setLightKeyIdx(oldHK);
      setHeavyKeyIdx(oldLK);
      setKeySwapMsg('Light/Heavy Keys were swapped to maintain volatility ordering.');
      setTimeout(() => setKeySwapMsg(null), 5000);

      // Recompute indices using swapped roles
      const newHkIndex = compounds.findIndex(c => c.id === parseInt(oldLK)); // previous LK is now HK
      const newLkIndex = compounds.findIndex(c => c.id === parseInt(oldHK)); // previous HK is now LK

      if (newHkIndex < 0 || newLkIndex < 0) {
        alert('Key swap failed due to invalid selection. Please re-select keys.');
        return;
      }

      const P_hk_new = getPsat(normalizedComps[newHkIndex], T_feed);
      const K_hk_new = P_hk_new / P_Pa;
      const compsWithAlphaNew = normalizedComps.map(c => {
        const Psat = getPsat(c, T_feed);
        const K = Psat / P_Pa;
        return { ...c, K, alpha: K / K_hk_new };
      });
      const alphaLKNew = compsWithAlphaNew[newLkIndex].alpha!;
      if (alphaLKNew <= 1.0001) {
        // Still invalid after swap
        alert('Unable to find a valid Light/Heavy key pairing at this pressure. Please reselect keys.');
        return;
      }
      // Use the recomputed arrays/values moving forward
      // Overwrite locals so rest of calculation uses updated arrays
      // (reassign variables used below)
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      hkIndex = newHkIndex;
      // eslint-disable-next-line @typescript-eslint/no-unused-vars
      // @ts-ignore
      lkIndex = newLkIndex;
      // Replace compsWithAlpha and alphaLK for the rest of the calculation
      // @ts-ignore
      compsWithAlpha = compsWithAlphaNew;
      // @ts-ignore
      alphaLK = alphaLKNew;
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
    // If a previous immediate calc suppressed this effect, clear the flag and skip
    if (suppressCompoundEffectRef.current) { suppressCompoundEffectRef.current = false; return; }

    // Require keys and at least 3 components
    if (compounds.length < 3 || !lightKeyIdx || !heavyKeyIdx) return;

    // Run calculation immediately on any input change for instant visual feedback
    calculateFUG();

    // No cleanup necessary
  }, [pressure, qValue, recLK, recHK, compounds, lightKeyIdx, heavyKeyIdx]);

  // Re-run flash calc when inputs affecting it change (By Temp or By Q)
  useEffect(() => {
    // If a previous immediate calc suppressed this effect, clear it and skip
    if (suppressCompoundEffectRef.current) { suppressCompoundEffectRef.current = false; return; }

    if (feedSpecMode === 'temp') {
      runFlashCalc();
      return;
    }
    // feedSpecMode === 'q' -> debounce compute from Q input
    const num = parseFloat(qInput);
    const t = setTimeout(() => {
      if (!isNaN(num)) computeFlashFromQ(num);
    }, 300);
    return () => clearTimeout(t);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [inputTemp, compounds, pressure, feedSpecMode, qInput]);

  // Clean up flash chart toggle when flash data is cleared
  useEffect(() => {
    if (!flashDetails) setShowFlashBarChart(false);
  }, [flashDetails]);

  // Force chart resize/re-render when switching chart mode to avoid stale state
  useEffect(() => {
    // If we have a ref to the chart, resize to ensure correct rendering
    const inst = echartsRef.current?.getEchartsInstance && echartsRef.current.getEchartsInstance();
    if (inst) {
      // Use a small timeout to allow DOM changes to settle
      setTimeout(() => inst.resize(), 50);
    }
  }, [showFlashBarChart, flashDetails]);

  // --- ECharts Options ---
  const gillilandOption = results ? {
    backgroundColor: 'transparent',
    // Global text style to enforce Merriweather Sans across all chart text
    textStyle: { fontFamily: 'Merriweather Sans, sans-serif', color: textColor },
    title: { 
        text: 'Gilliland Correlation', 
        left: 'center', 
        top: '4%',
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
            let html = `<div style="color:${color}; font-family:'Merriweather Sans', 'Merriweather Sans', sans-serif; font-size:14px; font-weight:400; line-height:1.2; margin:0; padding:6px 4px;">Reflux Parameter: <b>${formatNum(Number(xVal))}</b><br/>`;
            params.forEach((p: any) => {
                if (p.value && p.value.length === 2) {
                    const yVal = p.value[1];
                    html += `${p.seriesName}: <b>${formatNum(Number(yVal))}</b><br/>`;
                }
            });
            html += '</div>';
            return html;
        }
    },
    xAxis: { 
        name: '{main|(R-R}{sub|min}{main|)/(R+1)}',
        min: 0, max: 1, 
        nameLocation: 'middle', 
        nameGap: 30,
        axisLabel: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' },
        nameTextStyle: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif', rich: { main: { fontSize: 12, fontFamily: 'Merriweather Sans, sans-serif' }, sub: { fontSize: 10, padding: [6,0,0,0], fontFamily: 'Merriweather Sans, sans-serif' } } },
        axisLine: { lineStyle: { color: textColor } },
        splitLine: { show: false }
    },
    yAxis: { 
        name: '{main|(N-N}{sub|min}{main|)/(N+1)}', 
        min: 0, max: 1,
        nameLocation: 'middle',
        nameRotate: 90,
        nameGap: 40,
        axisLabel: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' },
        nameTextStyle: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif', rich: { main: { fontSize: 12, fontFamily: 'Merriweather Sans, sans-serif' }, sub: { fontSize: 10, padding: [6,0,0,0], fontFamily: 'Merriweather Sans, sans-serif' } } },
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
                 return `<div style="color:${color}; font-family:'Merriweather Sans', 'Merriweather Sans', sans-serif; font-size:14px; font-weight:500; line-height:1.2; margin:0; padding:6px 4px;">Design Point<br/>X: ${formatNum(p.value[0])}<br/>Y: ${formatNum(p.value[1])}</div>`;
             }
        }
      }
    ]
  } : {};  

  // Flash mode bar chart (phase composition)
  const flashBarOption = (flashDetails && results) ? {
    backgroundColor: 'transparent',
    textStyle: { fontFamily: 'Merriweather Sans, sans-serif', color: textColor },
    title: { text: 'Phase Composition', left: 'center', top: '4%', textStyle: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans, sans-serif' } },
    tooltip: {
      trigger: 'axis',
      axisPointer: { type: 'shadow' },
      renderMode: 'html',
      backgroundColor: resolvedTheme === 'dark' ? '#2563eb' : 'rgba(255,255,255,0.95)',
      borderColor: resolvedTheme === 'dark' ? '#2563eb' : '#ddd',
      padding: 8,
      textStyle: { fontFamily: 'Merriweather Sans, sans-serif', color: resolvedTheme === 'dark' ? '#ffffff' : '#000000' },
      formatter: (params: any) => {
        if (!Array.isArray(params) || params.length === 0) return '';
        const comp = params[0].axisValue;
        const color = resolvedTheme === 'dark' ? '#ffffff' : '#000000';
        let html = `<div style="color:${color}; font-family:'Merriweather Sans', 'Merriweather Sans', sans-serif; font-size:14px; font-weight:400; line-height:1.2; margin:0; padding:6px 4px;"><b>${comp}</b><br/>`;
        params.forEach((p: any) => {
          html += `${p.seriesName}: <b>${formatPhaseFractionFixed3(Number(p.value))}</b><br/>`;
        });
        html += '</div>';
        return html;
      }
    },
    legend: { data: ['Liquid', 'Vapor'], bottom: 4, left: 'center', orient: 'horizontal', textStyle: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' } },
    xAxis: { type: 'category', data: flashDetails.phases.map(p => p.name), axisLabel: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif' } },
    yAxis: { name: 'Mole Fraction', nameLocation: 'middle', nameRotate: 90, nameGap: 44, type: 'value', min: 0, max: 1, axisLabel: { color: textColor, fontFamily: 'Merriweather Sans, sans-serif', formatter: (val: any) => { const n = Number(val); if (!Number.isFinite(n)) return '0.0'; return n.toFixed(1); } }, axisLine: { show: true, lineStyle: { color: textColor, width: 1 } }, axisTick: { show: true, lineStyle: { color: textColor } }, splitLine: { show: false } },
    grid: { left: '8%', right: '8%', top: '14%', bottom: '12%' },
    series: [
      { name: 'Liquid', type: 'bar', data: flashDetails.phases.map(p => Number(formatPhaseFractionFixed3(p.x))), itemStyle: { color: '#60a5fa' } },
      { name: 'Vapor', type: 'bar', data: flashDetails.phases.map(p => Number(formatPhaseFractionFixed3(p.y))), itemStyle: { color: '#fb923c' } }
    ]
  } : {};

  return (
    <TooltipProvider>
      <div className="p-8 space-y-8 max-w-7xl mx-auto fug-range">
        <style>{`
          /* Page-scoped slider styling for native range inputs (improves contrast in light mode) */
          .fug-range input[type="range"] { -webkit-appearance: none; appearance: none; background: transparent; }
          .fug-range input[type="range"]::-webkit-slider-runnable-track {
            height: 8px; background: #cbd5e1; border-radius: 9999px; border: 1px solid #cbd5e1;
          }
          .fug-range input[type="range"]::-webkit-slider-thumb {
            -webkit-appearance: none; width: 18px; height: 18px; margin-top: -5px; border-radius: 9999px; background: #0f172a; border: none;
            box-shadow: 0 0 0 4px rgba(15,23,42,0.05);
          }
          .fug-range input[type="range"]::-moz-range-track { height: 8px; background: #cbd5e1; border-radius: 9999px; border: 1px solid #cbd5e1; }
          .fug-range input[type="range"]::-moz-range-thumb { width: 18px; height: 18px; background: #0f172a; border-radius: 9999px; border: none; }
          /* Dark mode overrides (keep existing dark appearance) */
          .dark .fug-range input[type="range"]::-webkit-slider-runnable-track { background: #2b3a49; border-color: transparent; }
          .dark .fug-range input[type="range"]::-moz-range-track { background: #2b3a49; border-color: transparent; }
        `}</style>
        <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
            
          {/* LEFT PANEL: INPUTS */}
          <div className="lg:col-span-4 space-y-6">
            
            {/* 1. Add Compounds Card (simplified layout, no scrolling titles) */}
            <Card>
              <CardContent className="pt-4 pb-4">
                <div className="relative mb-4">
                    <Search className="absolute left-3 top-2.5 h-4 w-4 text-muted-foreground" />
                    <Input 
                        placeholder="Search chemical..." 
                        value={searchQuery}
                        onChange={handleSearchChange}
                        onFocus={() => { if(suggestions.length > 0) setShowSuggestions(true); }}
                        onKeyDown={(e) => e.key === 'Enter' && handleAddCompound()}
                        autoComplete="off"
                        className="pl-9"
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
                
                {/* Removed scroll class to allow expansion */}
                <div className="space-y-3 pr-2">
                    {compounds.map((c, i) => (
                        <div key={i} className="flex flex-col space-y-2 text-sm bg-secondary/50 p-3 rounded-md">

                          <div className="flex items-center justify-between">
                            <div className="flex items-baseline min-w-0">
                              <div className="font-semibold min-w-0" title={c.name}>
                                <ShrinkToFit text={c.name} />
                              </div>
                              <div className="text-[12px] text-muted-foreground ml-2">Mole Fraction</div>
                            </div>

                            {/* Right: Icons */}
                            <div className="flex items-center space-x-2">
                              <ShadTooltip>
                                <TooltipTrigger asChild>
                                  <Button
                                    variant="ghost"
                                    size="sm"
                                    className="p-1"
                                    onClick={() => {
                                      const unlockedCount = compounds.filter(x => !x.locked).length;
                                      if (!c.locked && unlockedCount <= 2) return;
                                      setCompounds(prev => prev.map((p, idx) => idx === i ? { ...p, locked: !p.locked } : p));
                                    }}
                                    disabled={!c.locked && compounds.filter(x => !x.locked).length <= 2}
                                    aria-label={c.locked ? 'Unlock mole fraction' : 'Lock mole fraction'}
                                  >
                                    {c.locked ? <Lock className="w-3 h-3 text-green-500" /> : <Unlock className="w-3 h-3 text-muted-foreground" />}
                                  </Button>
                                </TooltipTrigger>
                                <TooltipContent>
                                  <div className="text-sm">{c.locked ? `Unlock: allow ${c.name} mole fraction to change` : `Hold ${c.name} mole fraction constant`}</div>
                                </TooltipContent>
                              </ShadTooltip>

                              <ShadTooltip>
                                <TooltipTrigger asChild>
                                  <Button
                                    variant="ghost"
                                    size="sm"
                                    className="p-1"
                                    onClick={() => removeCompound(i)}
                                    aria-label={`Remove ${c.name}`}
                                  >
                                    <Trash2 className="w-3 h-3 text-red-400" />
                                  </Button>
                                </TooltipTrigger>
                                <TooltipContent>
                                  <div className="text-sm">Remove {c.name}</div>
                                </TooltipContent>
                              </ShadTooltip>
                            </div>
                          </div>

                          {/* Slider row with reserved right space for numeric value */}
                          <div className="flex items-center gap-3">
                              <div className="flex-1">
                                <input 
                                    type="range" 
                                    min="0" max="1" step="0.01"
                                    value={c.z} 
                                    onChange={(e) => updateZ(i, parseFloat(e.target.value))}
                                    className={`w-full h-2 bg-slate-300 border border-slate-300 dark:bg-secondary dark:border-transparent rounded-lg appearance-none ${c.locked ? 'opacity-60 cursor-not-allowed' : 'cursor-pointer'} accent-primary`}
                                    disabled={c.locked}
                                />
                              </div>

                              <div className="w-8 text-right font-mono text-sm text-muted-foreground">
                                {formatSliderValue(c.z)}
                              </div>
                          </div>

                        </div>
                    ))}
                    {compounds.length === 0 && <div className="text-sm text-center text-muted-foreground py-4">No components added.</div>}
                </div>
              </CardContent>
            </Card>

            {/* 2. Process Specs */}
            <Card>
                <CardContent className="space-y-6 pt-4 pb-4">
  {/* Key Selection (Dropdowns) */}
  <div className="grid grid-cols-2 gap-4 min-h-[56px] items-center">
    <div className="flex flex-col justify-center space-y-1">
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
                if (c.id === hkId) isDisabled = true;
                if (c.Tb >= hkC.Tb) isDisabled = true;
              }
            }
            if (isDisabled) return null;
            return <SelectItem key={c.id} value={c.id.toString()}>{c.name}</SelectItem>;
          })}
        </SelectContent>
      </Select>
    </div>
    <div className="flex flex-col justify-center space-y-1">
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
                if (c.Tb <= lkC.Tb) isDisabled = true;
              }
            }
            if (isDisabled) return null;
            return <SelectItem key={c.id} value={c.id.toString()}>{c.name}</SelectItem>;
          })}
        </SelectContent>
      </Select>
    </div>
  </div>

  {keySwapMsg && (
    <div className="text-xs text-amber-400 mt-1">{keySwapMsg}</div>
  )}

  {/* Recovery Sliders */}
  <div className="space-y-4 pt-4">
    <div className="space-y-2">
      <div className="flex justify-between">
        <Label>LK Recovery (Distillate)</Label>
        <span className="text-xs font-mono text-muted-foreground">{formatNum(recLK * 100)}%</span>
      </div>
      <input 
        type="range" min="0.50" max="1.00" step="0.001"
        value={recLK} 
        onChange={e => setRecLK(parseFloat(e.target.value))}
        className="w-full h-2 bg-slate-300 border border-slate-300 dark:bg-secondary dark:border-transparent rounded-lg appearance-none cursor-pointer accent-blue-500"
      />
    </div>
    <div className="space-y-2">
      <div className="flex justify-between">
        <Label>HK Recovery (Bottoms)</Label>
        <span className="text-xs font-mono text-muted-foreground">{formatNum(recHK * 100)}%</span>
      </div>
      <input 
        type="range" min="0.50" max="1.00" step="0.001"
        value={recHK} 
        onChange={e => setRecHK(parseFloat(e.target.value))}
        className="w-full h-2 bg-slate-300 border border-slate-300 dark:bg-secondary dark:border-transparent rounded-lg appearance-none cursor-pointer accent-blue-500"
      />
    </div>
  </div>

  {/* REORDERED: Reflux Multiplier appears here first */}
  <div className="space-y-3 pt-4">
    <div className="flex justify-between">
      <Label>Reflux Multiplier <span className="text-xs">(R / R<sub>min</sub>)</span></Label>
      <span className="text-sm font-mono text-blue-400">{formatNum(R_factor)}</span>
    </div>
    <input 
      type="range" 
      min="1.05" max="3.0" step="0.05" 
      value={R_factor} 
      onChange={(e) => setRFactor(parseFloat(e.target.value))}
      className="w-full h-2 bg-slate-300 border border-slate-300 dark:bg-secondary dark:border-transparent rounded-lg appearance-none cursor-pointer accent-green-500"
    />
  </div>

  {/* Pressure & Flash Calculation */}
  <div className="grid grid-cols-3 gap-4 pt-4">
    <div className="space-y-1 col-span-1">
      <Label>Pressure</Label>
      <div className="relative flex items-center">
        <Input 
          type="number" 
          value={pressureInput} 
          onChange={e => handleNumericInputString(e.target.value, setPressureInput, setPressure, pressure, 0.01)} 
          step={0.01}
          title={pressureInput}
          className={`w-24 px-1 pr-12 overflow-hidden whitespace-nowrap ${pressureInput && pressureInput.length > 7 ? 'text-transparent caret-transparent' : 'text-center font-mono'}`}
        />
        <span className="absolute right-3 inset-y-0 flex items-center text-xs text-muted-foreground pointer-events-none">bar</span>
      </div>
    </div>

    {/* RACHFORD-RICE TOGGLE AREA */}
    <div className="space-y-1 col-span-2 relative pt-3.5">
      <div className="mb-1">
        <div className="absolute -top-8 left-0 z-20">
          <Tabs value={feedSpecMode} onValueChange={(v) => setFeedSpecMode(v === 'q' ? 'q' : 'temp')}>
            <TabsList className="inline-flex rounded-full bg-secondary p-1 shadow-sm">
              <TabsTrigger value="temp" className="text-xs px-2">Feed Temp</TabsTrigger>
              <TabsTrigger value="q" className="text-xs px-2">Feed Quality</TabsTrigger>
            </TabsList>
          </Tabs>
        </div>
      </div>

      {feedSpecMode === 'temp' ? (
        <div className="flex items-center space-x-2">
            <div className="relative flex-1">
                <Input 
                  type="number" 
                  value={inputTempInput} 
                  onChange={(e) => handleNumericInputString(e.target.value, setInputTempInput, setInputTemp, inputTemp, 1)} 
                  title={inputTempInput}
                  className={`pr-12 min-w-[88px] overflow-hidden whitespace-nowrap ${inputTempInput && inputTempInput.length > 7 ? 'text-transparent caret-transparent' : 'text-center font-mono'}`}
                />
                <span className="absolute right-3 inset-y-0 flex items-center text-xs text-muted-foreground">K</span>
            </div>
            {/* Arrow UI */}
            <div className="text-muted-foreground">→</div>
            <div className="flex items-center min-w-[80px] justify-center">
                <span className="font-mono text-sm">Q = {formatNum(qValue)}</span>
            </div>
        </div>
        ) : (
        <div className="flex items-center space-x-2">
            <div className="relative flex-1">
                <Input 
                  type="number" 
                  value={qInput} 
                  onChange={e => handleNumericInputString(e.target.value, setQInput, setQValue, qValue)} 
                  step={0.01}
                  title={qInput}
                  className={`pr-12 min-w-[88px] overflow-hidden whitespace-nowrap ${qInput && qInput.length > 7 ? 'text-transparent caret-transparent' : 'text-center font-mono'}`}
                />

            </div>
            {/* Arrow UI */}
            <div className="text-muted-foreground">→</div>
            <div className="flex items-center min-w-[80px] justify-center">
                <span className="font-mono text-sm">T = {inputTemp ? formatNum(inputTemp) : '—'} K</span>
            </div>
        </div>
        )}

    </div>
  </div>


</CardContent>
            </Card>
          </div>

          {/* RIGHT PANEL: RESULTS */}
          <div className="lg:col-span-8 space-y-6">
            {results ? (
                <>
                 <Card className="h-[560px]">
                    <CardContent className="h-[520px] relative">
                          <div className="absolute right-4 top-4 z-20">
                            <Tabs value={showFlashBarChart ? 'bar' : 'gill'} onValueChange={(v) => setShowFlashBarChart(v === 'bar')}>
                              <TabsList className="inline-flex rounded-full bg-secondary p-1 shadow-sm">
                                <TabsTrigger value="gill" className="text-xs px-3">Gilliland</TabsTrigger>
                                <TabsTrigger value="bar" className="text-xs px-3">Phase Comp</TabsTrigger>
                              </TabsList>
                            </Tabs>
                          </div>
                        <ReactECharts 
                            key={showFlashBarChart && flashDetails ? 'flashBar' : 'gilliland'}
                            ref={echartsRef} 
                            option={showFlashBarChart && flashDetails ? flashBarOption : gillilandOption} 
                            style={{ height: '100%', width: '100%' }} 
                        />
                    </CardContent>
                 </Card>

                 <div className="grid grid-cols-1 md:grid-cols-5 gap-4">
                    <Card className="bg-primary/5 border-primary/20">
                        <CardContent className="p-0 text-center">
                            <div className="text-sm text-muted-foreground">Minimum Stages</div>
                            <div className="text-xl font-bold">{formatNum(results.N_min)}</div>
                            <div className="text-xs text-muted-foreground">Fenske Eq.</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-primary/5 border-primary/20">
                        <CardContent className="p-0 text-center">
                            <div className="text-sm text-muted-foreground">Min Reflux Ratio.</div>
                            <div className="text-xl font-bold">{formatNum(results.R_min)}</div>
                            <div className="text-xs text-muted-foreground">Underwood Eq.</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-primary/5 border-primary/20">
                        <CardContent className="p-0 text-center">
                            <div className="text-sm text-muted-foreground">Feed Bubble Pt.</div>
                            {/* Converted to Celsius for readability */}
                            <div className="text-xl font-bold">{formatNum(results.FeedTemp - 273.15)} °C</div>
                            <div className="text-xs text-muted-foreground">{formatNum(results.FeedTemp)} K</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-blue-500/10 border-blue-500/20">
                        <CardContent className="p-0 text-center">
                            <div className="text-sm text-muted-foreground">Actual Stages</div>
                            <div className="text-2xl font-bold text-blue-500">{formatNum(results.N_actual)}</div>
                            <div className="text-xs text-muted-foreground">at R = {formatNum(results.R_actual)}</div>
                        </CardContent>
                    </Card>
                    <Card className="bg-blue-500/10 border-blue-500/20">
                        <CardContent className="p-0 text-center">
                            <div className="text-sm text-muted-foreground">Feed Location</div>
                            <div className="text-2xl font-bold text-blue-500">{formatNum(results.FeedStage)}</div>
                            <div className="text-xs text-muted-foreground">Stages from Top</div>
                        </CardContent>
                    </Card>
                 </div>

        {/* 3. NEW LOCATION: Phase Composition Table (Full Width) */}
        <Card>
                <CardHeader className="pb-2">
                    <CardTitle className="text-sm font-medium flex justify-between items-center">
                        <span>Phase Equilibrium (Flash Calculation)</span>
                        {flashDetails ? (
                        <span className={`text-xs px-2 py-1 rounded ${
                            flashDetails.psi > 0 && flashDetails.psi < 1 
                            ? "bg-blue-100 text-blue-700 dark:bg-blue-900 dark:text-blue-100" 
                            : "bg-secondary text-secondary-foreground"
                        }`}>
                            {flashDetails.stateMsg} — {formatNum(flashDetails.psi * 100)}% Vapor
                        </span>
                        ) : (
                          <span className="text-xs text-muted-foreground">Feed spec will generate phase data here</span>
                        )}
                        {/* Show flash error messages (non-state messages) here without affecting feed input layout */}
                        {flashState && flashState.msg && !["Subcooled Liquid","Superheated Vapor","Two-Phase Mixture"].includes(flashState.msg) && (
                          <div className="text-xs text-amber-400 ml-3">{flashState.msg}</div>
                        )}
                    </CardTitle>
                </CardHeader>
                <CardContent>
                    <div className="overflow-x-auto">
                    {flashDetails ? (
                        <table className="w-full text-sm text-left">
                            <thead className="text-xs text-muted-foreground uppercase bg-secondary/50">
                                <tr>
                                    <th className="px-4 py-2 rounded-l">Component</th>
                                    <th className="px-4 py-2 text-right">Feed (z)</th>
                                    <th className="px-4 py-2 text-right text-blue-500">Liquid (x)</th>
                                    <th className="px-4 py-2 text-right text-orange-500 rounded-r">Vapor (y)</th>
                                </tr>
                            </thead>
                            <tbody>
                                {flashDetails.phases.map((p, i) => (
                                    <tr key={i} className="border-b last:border-0 hover:bg-muted/50">
                                        <td className="px-4 py-2 font-medium">{p.name}</td>
                                        <td className="px-4 py-2 text-right font-mono">{formatNum(p.z)}</td>
                                        <td className="px-4 py-2 text-right font-mono text-blue-600 dark:text-blue-400">{formatPhaseFraction(p.x)}</td>
                                        <td className="px-4 py-2 text-right font-mono text-orange-600 dark:text-orange-400">{formatPhaseFraction(p.y)}</td>
                                    </tr>
                                ))}
                            </tbody>
                        </table>
                    ) : (
                        <div className="text-sm text-muted-foreground py-6 text-center">No phase data yet. Select a feed spec or enter a feed quality.</div>
                    )}
                    </div>
                </CardContent>
            </Card>
                </>
            ) : (
                <div className="h-full flex flex-col items-center justify-center text-muted-foreground border-2 border-dashed rounded-lg p-12">
                    <Calculator className="w-12 h-12 mb-4 opacity-50" />
                    <p>Add at least 3 components to run the simulation.</p>
                </div>
            )}
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}