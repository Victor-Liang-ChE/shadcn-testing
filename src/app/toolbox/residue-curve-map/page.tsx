'use client';

import React, { useState, useEffect, useRef, useCallback, useMemo } from 'react';
import dynamic from 'next/dynamic';
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { createClient } from '@supabase/supabase-js';
import type { SupabaseClient } from '@supabase/supabase-js';
import { 
    CompoundData, 
    WilsonPureComponentParams as WilsonPureParams, 
    UnifacGroupComposition,
    PrPureComponentParams, 
    SrkPureComponentParams, 
    UniquacPureComponentParams,
    AntoineParams, 
} from '@/lib/vle-types';
import { 
    simulateResidueCurveODE,
    type FluidPackageResidue
} from '@/lib/residue-curves-ode';
import { 
    TernaryWilsonParams,
    TernaryNrtlParams,
    TernaryPrParams,
    TernarySrkParams,
    TernaryUniquacParams,
    findAzeotropeNRTL,
    findAzeotropePr,
    findAzeotropeSrk,
    findAzeotropeUniquac,
    findAzeotropeWilson,
    findAzeotropeUnifac,
    systematicBinaryAzeotropeSearch,
    estimateTernaryGuessFromBinary,
    systematicTernaryAzeotropeSearch,
    type AzeotropeResult
} from '@/lib/residue-curves-ode';
import {
  R_gas_const_J_molK,
  calculatePsat_Pa,
  calculateSaturationTemperaturePurePr,
  calculateSaturationTemperaturePureSrk,
  fetchWilsonInteractionParams,
  fetchUnifacInteractionParams,
  fetchNrtlParameters,
  fetchPrInteractionParams,
  fetchSrkInteractionParams,
  fetchUniquacInteractionParams,
  type UnifacParameters,
  type NrtlInteractionParams,
  type PrInteractionParams,
  type SrkInteractionParams,
  type UniquacInteractionParams,
} from '@/lib/vle-calculations';

import { fetchAndConvertThermData, FetchedCompoundThermData } from '@/lib/antoine-utils'; 
import type { Data, Layout } from 'plotly.js'; 
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"; 
import { Table, TableBody, TableCell, TableHead, TableHeader, TableRow } from "@/components/ui/table";
import {
  evaluateResidueODE
} from '@/lib/residue-curves-ode';

// Dynamically import Plotly to avoid SSR issues
const Plot = dynamic(
  () => import('react-plotly.js').then((mod) => mod.default),
  {
    ssr: false,
    loading: () => <div className="flex justify-center items-center h-[600px]"><p>Loading plot...</p></div>
  }
);

// Initialize Supabase client (replace with your actual URL and anon key)
// Ensure these are environment variables in a real application
const SUPABASE_URL = process.env.NEXT_PUBLIC_SUPABASE_URL || "your_supabase_url";
const SUPABASE_ANON_KEY = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || "your_supabase_anon_key";

let supabase: SupabaseClient | null = null;
if (SUPABASE_URL !== "your_supabase_url" && SUPABASE_ANON_KEY !== "your_supabase_anon_key") {
    supabase = createClient(SUPABASE_URL, SUPABASE_ANON_KEY);
} else {
    console.warn("Supabase client not initialized. Please provide NEXT_PUBLIC_SUPABASE_URL and NEXT_PUBLIC_SUPABASE_ANON_KEY environment variables.");
}

// Debounce function removed for instantaneous suggestions

interface ComponentInputState {
    name: string; 
}

// To store component data along with its original index and NBP for sorting
interface ProcessedComponentData {
    name: string;
    casNumber: string;
    thermData: FetchedCompoundThermData; // This now contains criticalProperties
    originalIndex: number; // 0, 1, or 2 based on input order
    bp_at_Psys_K?: number | null; // Boiling point at system pressure
    isSCF?: boolean;
}

const initialComponentState = (): ComponentInputState => ({
    name: '',
});

// Default average boiling point for initial guess if NBP is missing
const DEFAULT_INITIAL_T_K = 350; // Kelvin, e.g., around 77°C

// Helper function to convert ternary coordinates to paper coordinates for annotations
const convertTernaryToPaperCoordinates = (
    p_x_array: number[], // Original mole fractions [x_comp0, x_comp1, x_comp2]
    plotSortedCompDefs: ProcessedComponentData[] // [Lightest, Intermediate, Heaviest]
): { x: number; y: number } | null => {
    if (plotSortedCompDefs.length !== 3) return null;

    const x_light = p_x_array[plotSortedCompDefs[0].originalIndex];
    const x_intermediate = p_x_array[plotSortedCompDefs[1].originalIndex];
    const x_heavy = p_x_array[plotSortedCompDefs[2].originalIndex];

    const x_paper = x_light * 1 + x_intermediate * 0.5; // More direct: x_paper = x_c_axis_data + x_a_axis_data * 0.5
    const y_paper = x_intermediate * (Math.sqrt(3) / 2);   // y_paper = x_a_axis_data * Math.sqrt(3)/2

    return { x: x_paper, y: y_paper };
};

export type FluidPackageTypeResidue = 'wilson' | 'unifac' | 'nrtl' | 'pr' | 'srk' | 'uniquac';

// define a tiny generic curve type that all models satisfy
interface ResidueCurvePoint { 
  x: number[]; 
  T_K: number; 
  step?: number; 
}
type ResidueCurve = ResidueCurvePoint[];


// Define a common interface for displaying azeotropes on the page
interface AzeotropeDisplayInfo {
    x: number[];
    T_K: number;
    errorNorm?: number;
    // Other fields like 'converged', 'iterations', 'message' can be part of specific results
    // but are not strictly needed for basic plotting.
}

/**
 * Analyzes residue curves to find approximate azeotrope locations.
 * Azeotropes are the start/end points of curves that are not pure components.
 * @param curves The array of residue curves from the ODE simulation.
 * @returns An array of unique, approximate azeotrope compositions.
 */
function findApproximateAzeotropesFromCurves(
  curves: ResidueCurve[],
): number[][] {
  const potentialAzeotropes: number[][] = [];
  const PURE_COMPONENT_THRESHOLD = 0.999; 

  for (const curve of curves) {
    if (curve.length < 2) continue;

    const startPoint = curve[0].x;
    const endPoint = curve[curve.length - 1].x;

    if (!startPoint.some(frac => frac > PURE_COMPONENT_THRESHOLD)) {
      potentialAzeotropes.push(startPoint);
    }
    if (!endPoint.some(frac => frac > PURE_COMPONENT_THRESHOLD)) {
      potentialAzeotropes.push(endPoint);
    }
  }

  // --- REVISED DE-DUPLICATION LOGIC ---
  const uniqueAzeotropes: number[][] = [];
  const binaryAzeotropes: number[][] = [];
  const ternaryAzeotropes: number[][] = [];
  
  const BINARY_THRESHOLD = 0.015; // Mole fraction below this is considered 'on the edge'
  const MIN_SQ_DIST_TERNARY = 0.02 * 0.02; // Wider tolerance for ternary azeotropes
  const MIN_SQ_DIST_BINARY = 0.02 * 0.02;  // Separate tolerance for binary azeotropes

  for (const p of potentialAzeotropes) {
    if (p.some(frac => frac < BINARY_THRESHOLD)) {
      binaryAzeotropes.push(p);
    } else {
      ternaryAzeotropes.push(p);
    }
  }

  // De-duplicate ternary points
  const uniqueTernary: number[][] = [];
  for (const p of ternaryAzeotropes) {
    if (!uniqueTernary.some(up => ((p[0]-up[0])**2 + (p[1]-up[1])**2 + (p[2]-up[2])**2) < MIN_SQ_DIST_TERNARY)) {
      uniqueTernary.push(p);
    }
  }

  // De-duplicate binary points
  const uniqueBinary: number[][] = [];
  for (const p of binaryAzeotropes) {
    if (!uniqueBinary.some(up => ((p[0]-up[0])**2 + (p[1]-up[1])**2 + (p[2]-up[2])**2) < MIN_SQ_DIST_BINARY)) {
      uniqueBinary.push(p);
    }
  }

  return [...uniqueTernary, ...uniqueBinary];
}

// Map per-model azeotrope result names to unified interface
type NrtlAzeotropeResult = AzeotropeResult;
type PrAzeotropeResult = AzeotropeResult;
type SrkAzeotropeResult = AzeotropeResult & { y: number[] };
type UniquacAzeotropeResult = AzeotropeResult;

export default function TernaryResidueMapPage() {
    // Add cache reference for efficient data handling
    const componentDataCache = useRef(new Map<string, { casNumber: string, thermData: FetchedCompoundThermData } | null>());
    
    const [componentsInput, setComponentsInput] = useState<ComponentInputState[]>([
        { name: 'Acetone' },
        { name: 'Chloroform' },
        { name: 'Methanol' },
    ]);
    const [componentSuggestions, setComponentSuggestions] = useState<string[][]>([[], [], []]); // Suggestions for each component
    const [showSuggestions, setShowSuggestions] = useState<boolean[]>([false, false, false]); // Control dropdown visibility
    const [activeSuggestionIndex, setActiveSuggestionIndex] = useState<number | null>(null);
    const inputRefs = useRef<(HTMLInputElement | null)[]>([]);
    const suggestionsContainerRefs = useRef<(HTMLUListElement | null)[]>([]);
    const [systemPressure, setSystemPressure] = useState<string>('1'); // System pressure in bar, default to 1
    const [fluidPackage, setFluidPackage] = useState<FluidPackageTypeResidue>('wilson'); // New state for fluid package
    const [displayedFluidPackage, setDisplayedFluidPackage] = useState<FluidPackageTypeResidue>(fluidPackage); // State for title
    const [residueCurves, setResidueCurves] = useState<ResidueCurve[]>([]);
    const [cleanedResidueCurves, setCleanedResidueCurves] = useState<ResidueCurve[]>([]); // For filtered curves
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [plotlyData, setPlotlyData] = useState<Data[]>([]);
    const [plotlyLayout, setPlotlyLayout] = useState<Partial<Layout>>({});
    const [plotSortedComponents, setPlotSortedComponents] = useState<ProcessedComponentData[]>([]);
    const [plotAxisTitles, setPlotAxisTitles] = useState({ a: 'Comp 2 (M)', b: 'Comp 3 (H)', c: 'Comp 1 (L)' });
    const initialMapGenerated = useRef(false); // Ref to track initial generation
    const [currentTheme, setCurrentTheme] = useState<'dark' | 'light'>('dark'); // Default to dark
    const [plotContainerRef, setPlotContainerRef] = React.useState<HTMLDivElement | null>(null); // State to hold the div element
    const [directAzeotropes, setDirectAzeotropes] = useState<AzeotropeDisplayInfo[]>([]); // Use common display info
    const [scfMessage, setScfMessage] = useState<string | null>(null);
    const [displayedPressure, setDisplayedPressure] = useState<string>(systemPressure);
    // Names frozen after last successful generation to prevent flicker while typing
    const [displayedComponentNames, setDisplayedComponentNames] = useState<string[]>(componentsInput.map(c=>c.name));
    // Track which azeotrope (if any) is being hovered so we can highlight it on the plot
    const [highlightedAzeoIdx, setHighlightedAzeoIdx] = useState<number | null>(null);
    const [backendComps, setBackendComps] = useState<CompoundData[]>([]);
    const [backendPkgParams, setBackendPkgParams] = useState<any>(null);
    const [backendPressurePa, setBackendPressurePa] = useState<number>(0);

    // Removed debounced version for direct calls
    const setHighlightedAzeoIdxDirect = setHighlightedAzeoIdx;

    // Trigger re-generation automatically when the fluid-package selection changes
    const didMountFluidPkg = useRef(false);
    useEffect(() => {
        if (didMountFluidPkg.current) {
            handleGenerateMap();
        } else {
            didMountFluidPkg.current = true;
        }
        // eslint-disable-next-line react-hooks/exhaustive-deps
    }, [fluidPackage]);

    useEffect(() => {
        // Function to check and set theme
        const updateTheme = () => {
            if (typeof window !== 'undefined') {
                const isDarkMode = document.documentElement.classList.contains('dark');
                setCurrentTheme(isDarkMode ? 'dark' : 'light');
            }
        };

        updateTheme(); // Initial check

        // Observe class changes on <html> for dynamic theme switching
        if (typeof window !== 'undefined' && window.MutationObserver) {
            const observer = new MutationObserver(mutationsList => {
                for (const mutation of mutationsList) {
                    if (mutation.type === 'attributes' && mutation.attributeName === 'class') {
                        updateTheme();
                        break;
                    }
                }
            });
            observer.observe(document.documentElement, { attributes: true, attributeFilter: ['class'] });
            return () => observer.disconnect(); // Cleanup observer on unmount
        }
    }, []);

    const fetchComponentSuggestions = async (query: string, index: number) => {
        if (!supabase || !query.trim()) {
            setComponentSuggestions(prev => {
                const updated = [...prev];
                updated[index] = [];
                return updated;
            });
            setShowSuggestions(prev => {
                const newShow = [...prev];
                newShow[index] = false;
                return newShow;
            });
            return;
        }
        try {
            const { data, error } = await supabase
                .from('compounds')
                .select('name')
                .ilike('name', `${query.trim()}%`) // Changed to "starts with" pattern
                .limit(5);

            if (error) {
                console.error(`Error fetching suggestions for "${query}":`, error.message);
                setComponentSuggestions(prev => {
                    const updated = [...prev];
                    updated[index] = [];
                    return updated;
                });
                setShowSuggestions(prev => {
                    const newShow = [...prev];
                    newShow[index] = false;
                    return newShow;
                });
                return;
            }

            if (data) {
                const suggestions = data.map((item: { name: string }) => item.name);
                setComponentSuggestions(prev => {
                    const updated = [...prev];
                    updated[index] = suggestions;
                    return updated;
                });
                setShowSuggestions(prev => {
                    const newShow = [...prev];
                    newShow[index] = suggestions.length > 0;
                    return newShow;
                });
            }
        } catch (err) {
            console.error(`Error fetching suggestions for "${query}":`, err);
            setComponentSuggestions(prev => {
                const updated = [...prev];
                updated[index] = [];
                return updated;
            });
            setShowSuggestions(prev => {
                const newShow = [...prev];
                newShow[index] = false;
                return newShow;
            });
        }
    };
    
    // Removed debounced version for instantaneous suggestions

    const handleComponentInputChange = (index: number, field: keyof ComponentInputState, value: string) => {
        const newInputs = [...componentsInput];
        newInputs[index] = { ...newInputs[index], [field]: value };
        setComponentsInput(newInputs);

        if (field === 'name') {
            setActiveSuggestionIndex(index);
            if (value.trim() === "") {
                setShowSuggestions(prev => {
                    const updated = [...prev];
                    updated[index] = false;
                    return updated;
                });
                setComponentSuggestions(prev => {
                    const updated = [...prev];
                    updated[index] = [];
                    return updated;
                });
            } else {
                // Direct call instead of debounced call
                fetchComponentSuggestions(value, index);
            }
        }
    };

    const handleSuggestionClick = (index: number, suggestion: string) => {
        const newInputs = [...componentsInput];
        newInputs[index].name = suggestion;
        setComponentsInput(newInputs);

        setShowSuggestions(prev => {
            const updated = [...prev];
            updated[index] = false;
            return updated;
        });
        setComponentSuggestions(prev => {
            const updated = [...prev];
            updated[index] = [];
            return updated;
        });
        inputRefs.current[index]?.focus();
    };

    useEffect(() => {
        function handleClickOutside(event: MouseEvent) {
            if (activeSuggestionIndex !== null) {
                const inputEl = inputRefs.current[activeSuggestionIndex];
                const suggestionsEl = suggestionsContainerRefs.current[activeSuggestionIndex];
                if (
                    inputEl && !inputEl.contains(event.target as Node) &&
                    suggestionsEl && !suggestionsEl.contains(event.target as Node)
                ) {
                    setShowSuggestions(prev => {
                        const updated = [...prev];
                        if (activeSuggestionIndex !== null) {
                           updated[activeSuggestionIndex] = false;
                        }
                        return updated;
                    });
                }
            }
        }
        document.addEventListener("mousedown", handleClickOutside);
        return () => {
            document.removeEventListener("mousedown", handleClickOutside);
        }
    }, [activeSuggestionIndex]);

    const generateStartingPoints = (
        numPerEdge: number = 8,
        numInternalDivisions: number = 10, // New parameter for systematic internal points
        fluidPackageType?: FluidPackageTypeResidue // To handle specific seeds like PR/SRK corners
    ): number[][] => {
        const points: number[][] = [];
        const edgeEpsilon = 0.005; // Start slightly off the actual edge for numerical stability
    
        // Helper to add point if valid and normalized
        const addNormalizedPoint = (p: number[]) => {
            const sum = p[0] + p[1] + p[2];
            if (sum <= 1e-6) return; // Avoid division by zero or near-zero sum
            const normP = [p[0] / sum, p[1] / sum, p[2] / sum];
            
            // Filters points too close to pure (e.g., >0.999 or <0.001 for a component)
            // This is generally okay as curves should naturally reach these extremes.
            const minFrac = 1e-3; 
            if (normP.every(val => val >= minFrac && val <= 1.0 - minFrac)) {
                // Avoid adding effectively duplicate points
                if (!points.some(existingP => 
                    Math.abs(existingP[0] - normP[0]) < 1e-4 &&
                    Math.abs(existingP[1] - normP[1]) < 1e-4 &&
                    Math.abs(existingP[2] - normP[2]) < 1e-4
                )) {
                    points.push(normP);
                }
            } else {
                // If point was filtered by minFrac, log it for debugging if necessary,
                // but generally this filter is intended.
                // console.log("Filtered point by minFrac: ", normP, " from original: ", p);
            }
        };
        
        // 1. Points along edges
        for (let i = 1; i <= numPerEdge; i++) {
            const x_ratio = i / (numPerEdge + 1.0);
            addNormalizedPoint([x_ratio * (1.0 - edgeEpsilon), (1.0 - x_ratio) * (1.0 - edgeEpsilon), edgeEpsilon]); // Edge 0-1
            addNormalizedPoint([x_ratio * (1.0 - edgeEpsilon), edgeEpsilon, (1.0 - x_ratio) * (1.0 - edgeEpsilon)]); // Edge 0-2
            addNormalizedPoint([edgeEpsilon, x_ratio * (1.0 - edgeEpsilon), (1.0 - x_ratio) * (1.0 - edgeEpsilon)]); // Edge 1-2
        }
    
        // 2. Systematic Internal Points
        const N = numInternalDivisions;
        for (let i = 1; i < N; i++) { // i corresponds to x0 * N
            const x0 = i / N;
            for (let j = 1; j < N - i; j++) { // j corresponds to x1 * N
                const x1 = j / N;
                const x2 = 1.0 - x0 - x1; // k = N - i - j, so x2 = k/N
                // By loop construction, x0, x1, x2 are all >= 1/N.
                // addNormalizedPoint will handle normalization and its own minFrac check.
                addNormalizedPoint([x0, x1, x2]);
            }
        }
    
        // 3. User-defined specific internal points (from the original list)
        const internalPointsSeed = [
            [0.333, 0.333, 0.334],
            [0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5],
            [0.1, 0.1, 0.8], [0.1, 0.8, 0.1], [0.8, 0.1, 0.1],
            [0.2, 0.2, 0.6], [0.2, 0.6, 0.2], [0.6, 0.2, 0.2],
            [0.4, 0.4, 0.2], [0.4, 0.2, 0.4], [0.2, 0.4, 0.4],
            [0.1, 0.2, 0.7], [0.1, 0.7, 0.2], [0.2, 0.1, 0.7], [0.2, 0.7, 0.1], [0.7, 0.1, 0.2], [0.7, 0.2, 0.1],
            [0.15, 0.35, 0.5], [0.15, 0.5, 0.35], [0.35, 0.15, 0.5], [0.35, 0.5, 0.15], [0.5, 0.15, 0.35], [0.5, 0.35, 0.15],
            [0.2, 0.3, 0.5], [0.2, 0.5, 0.3], [0.3, 0.2, 0.5], [0.3, 0.5, 0.2], [0.5, 0.2, 0.3], [0.5, 0.3, 0.2],
            [0.3, 0.3, 0.4], [0.3, 0.4, 0.3], [0.4, 0.3, 0.3],
        ];
        internalPointsSeed.forEach(p => addNormalizedPoint(p));
        
        // 4. Special near-corner seeds for PR/SRK, added via addNormalizedPoint
        if (fluidPackageType === 'pr' || fluidPackageType === 'srk') {
            const cornerEps = 0.005; // This is the same as edgeEpsilon.
                                     // addNormalizedPoint will ensure they are valid and not overly redundant.
            // These points like [0.99, 0.005, 0.005] should pass the current addNormalizedPoint minFrac filter.
            addNormalizedPoint([1 - 2 * cornerEps, cornerEps, cornerEps]); // e.g. [0.99, 0.005, 0.005] - slightly more internal
            addNormalizedPoint([cornerEps, 1 - 2 * cornerEps, cornerEps]);
            addNormalizedPoint([cornerEps, cornerEps, 1 - 2 * cornerEps]);
        }
        
        return points;
    };

    // Cached component data fetching function
    async function fetchComponentData(compoundName: string): Promise<{ casNumber: string, thermData: FetchedCompoundThermData } | null> {
        const cacheKey = compoundName.toLowerCase();
        
        // 1. Check the cache first
        if (componentDataCache.current.has(cacheKey)) {
            console.log(`ResidueCurveMap: Cache HIT for ${compoundName}.`);
            const cachedData = componentDataCache.current.get(cacheKey);
            if (cachedData === null) {
                setError(`Data fetch previously failed for ${compoundName}. Please check the name.`);
            }
            return cachedData || null;
        }

        console.log(`ResidueCurveMap: Cache MISS for ${compoundName}. Fetching from DB.`);
        if (!supabase) throw new Error("Supabase client is not available.");

        try {
            // 2. Fetch CAS number
            const casNumber = await fetchCasNumberByName(supabase, compoundName);
            
            // 3. Fetch thermodynamic data
            const thermData = await fetchAndConvertThermData(supabase, casNumber);
            
            const finalComponentData = { casNumber, thermData };

            // 4. Store the successful result in the cache
            componentDataCache.current.set(cacheKey, finalComponentData);
            return finalComponentData;

        } catch (err: any) {
            console.error(`Error fetching data for ${compoundName} in ResidueCurveMap:`, err.message);
            setError(err.message);
            
            // 5. Store failure in the cache to prevent re-fetching invalid names
            componentDataCache.current.set(cacheKey, null);
            return null;
        }
    }

    const fetchCasNumberByName = async (supabaseClient: SupabaseClient, name: string): Promise<string> => {
        if (!name.trim()) {
            console.error("fetchCasNumberByName: Name is empty.");
            throw new Error("Component name cannot be empty.");
        }
        const trimmedName = name.trim();
        console.log(`fetchCasNumberByName: Fetching CAS for name: "${trimmedName}"`); // Added log
        const { data, error } = await supabaseClient
            .from('compounds')
            .select('cas_number, name') // Also select name for logging
            .ilike('name', trimmedName) // Changed from %${trimmedName}% to trimmedName for exact match (case-insensitive)
            .limit(1)
            .single();

        if (error) {
            console.error(`fetchCasNumberByName: Error fetching CAS for "${trimmedName}":`, error);
            throw new Error(`Could not fetch CAS number for ${trimmedName}: ${error.message}`);
        }
        if (!data || !data.cas_number) {
            console.error(`fetchCasNumberByName: No CAS number found for "${trimmedName}". Data received:`, data);
            throw new Error(`No CAS number found for ${trimmedName}.`);
        }
        console.log(`fetchCasNumberByName: Found CAS: ${data.cas_number} for DB name: "${data.name}" (queried as "${trimmedName}")`); // Added log
        return data.cas_number;
    };

    // Helper function to calculate boiling point of a pure component at P_sys_Pa using Secant method
    function calculateBoilingPointAtPressureSecant(
        antoine: AntoineParams,
        P_sys_Pa: number,
        initialT_K: number,
        psatFunc: (antoine: AntoineParams, T_K: number) => number,
        maxIter: number = 30,
        relTolerance: number = 1e-5, // Relative tolerance for f(T)/P_sys
        absTolerancePa: number = 0.1, // Absolute tolerance for f(T) in Pa
        tempStepTol: number = 0.01 // Kelvin tolerance for T_new - T_curr
    ): number | null {
        let T_curr = initialT_K > 0 ? initialT_K : (antoine.Tmin_K || 273.15);
        let f_curr_val = psatFunc(antoine, T_curr) - P_sys_Pa;

        if (isNaN(f_curr_val)) return null;

        // Initial perturbation for T_prev
        let T_prev = T_curr * (f_curr_val > 0 ? 1.05 : 0.95); // Perturb based on f_curr sign
        if (T_prev <= 0 || Math.abs(T_prev - T_curr) < tempStepTol / 2) { // Ensure T_prev is valid and distinct
            T_prev = T_curr + (T_curr > antoine.Tmin_K + tempStepTol ? -tempStepTol : tempStepTol); 
        }
        if (T_prev <=0) T_prev = T_curr * 0.9; // Final fallback for T_prev

        let f_prev_val = psatFunc(antoine, T_prev) - P_sys_Pa;

        if (isNaN(f_prev_val)) { // If first perturbation failed, try another direction
            T_prev = T_curr + (T_prev > T_curr ? -2 * (T_prev - T_curr) : 2 * (T_curr - T_prev)); // Reverse and double perturbation
            if (T_prev <=0) T_prev = T_curr * 1.1; // Fallback
            f_prev_val = psatFunc(antoine, T_prev) - P_sys_Pa;
            if (isNaN(f_prev_val)) return null; // Still NaN, give up
        }

        for (let i = 0; i < maxIter; i++) {
            if (Math.abs(f_curr_val / P_sys_Pa) < relTolerance || Math.abs(f_curr_val) < absTolerancePa) {
                return T_curr;
            }

            const denominator = f_curr_val - f_prev_val;
            if (Math.abs(denominator) < 1e-6) { // Denominator too small
                return (Math.abs(f_curr_val / P_sys_Pa) < relTolerance * 10) ? T_curr : null;
            }

            const T_next = T_curr - f_curr_val * (T_curr - T_prev) / denominator;

            if (isNaN(T_next) || !isFinite(T_next) || T_next <= 0 || T_next < (antoine.Tmin_K || 0) * 0.8 || T_next > (antoine.Tmax_K || Infinity) * 1.2 ) {
                T_prev = T_curr;
                f_prev_val = f_curr_val;
                T_curr = T_curr + (f_curr_val > 0 ? -0.5 : 0.5) * (T_curr * 0.005 + 0.05);
                if (T_curr <=0) T_curr = initialT_K * (0.5 + Math.random()*0.5);
                f_curr_val = psatFunc(antoine, T_curr) - P_sys_Pa;
                if(isNaN(f_curr_val)) return null;
                continue;
            }
            
            if (Math.abs(T_next - T_curr) < tempStepTol) {
                 return T_next;
            }

            T_prev = T_curr;
            f_prev_val = f_curr_val;
            T_curr = T_next;
            f_curr_val = psatFunc(antoine, T_curr) - P_sys_Pa;
            if (isNaN(f_curr_val)) return null;
        }
        return (Math.abs(f_curr_val / P_sys_Pa) < relTolerance * 10) ? T_curr : null;
    }

    const handleGenerateMap = useCallback(async () => {
        if (!supabase) {
            setError("Supabase client is not initialized. Cannot fetch parameters.");
            return;
        }
        const generatingFluidPackage = fluidPackage; // Capture fluid package for this generation run
        setIsLoading(true);
        setError(null);

        try {
            const P_system_Pa = parseFloat(systemPressure) * 100000; // Convert bar to Pa
            if (isNaN(P_system_Pa) || P_system_Pa <= 0) {
                throw new Error("Invalid system pressure input.");
            }

            // Use cached data fetching
            const componentDataPromises = componentsInput.map(input => 
                fetchComponentData(input.name)
            );
            const fetchedComponentDataArray = await Promise.all(componentDataPromises);
            
            // Extract cas numbers and therm data from cached results
            const fetchedCasNumbers = fetchedComponentDataArray.map(data => data?.casNumber || '');
            const fetchedThermDataArray = fetchedComponentDataArray.map(data => data?.thermData || null);

            const processedComponents: ProcessedComponentData[] = componentsInput.map((input, index) => {
                const thermData = fetchedThermDataArray[index];
                if (!thermData) {
                    throw new Error(`Failed to fetch thermodynamic data for ${input.name}`);
                }
                if (thermData.criticalProperties) {
                    const { Tc_K, Pc_Pa, omega } = thermData.criticalProperties;
                    // console.log(`EOS BP attempt → ${input.name}: Tc=${Tc_K.toFixed(2)} K, Pc=${(Pc_Pa/1e5).toFixed(2)} bar, ω=${omega.toFixed(3)}`);
                } else {
                    // console.log(`EOS BP attempt → ${input.name}: critical properties NOT found`);
                }
                let bp_at_Psys_K_val: number | null = null;

                // Attempt EOS-based boiling point first if critical properties are available
                if (thermData.criticalProperties) {
                    const compObj: CompoundData = {
                        name: input.name,
                        antoine: thermData.antoine,
                        prParams: thermData.criticalProperties as PrPureComponentParams,
                        srkParams: thermData.criticalProperties as SrkPureComponentParams,
                        unifacGroups: null,
                        uniquacParams: null,
                        wilsonParams: null,
                    } as unknown as CompoundData;
                    bp_at_Psys_K_val = calculateSaturationTemperaturePurePr(compObj, P_system_Pa);
                    if (bp_at_Psys_K_val !== null && !isNaN(bp_at_Psys_K_val)) {
                        // console.log(`EOS-PR BP → ${input.name}: T_sat@${(P_system_Pa/1e5).toFixed(2)} bar = ${bp_at_Psys_K_val.toFixed(2)} K`);
                    } else {
                        // console.warn(`EOS-PR BP failed for ${input.name}`);
                    }

                    const { Tc_K, Pc_Pa } = thermData.criticalProperties;
                    // If operating above BOTH Pc and Tc – treat as supercritical fluid
                    const isSCF = (P_system_Pa > Pc_Pa) && (bp_at_Psys_K_val !== null ? bp_at_Psys_K_val > Tc_K : true);
                    if (isSCF) {
                        // console.log(`Supercritical conditions detected for ${input.name} (Pc=${(Pc_Pa/1e5).toFixed(2)} bar, Tc=${Tc_K.toFixed(1)} K). Marking as SCF.`);
                        bp_at_Psys_K_val = Number.NaN; // Sentinel for SCF
                    }
                } else if (thermData.antoine) {
                    // Fallback to Antoine
                    bp_at_Psys_K_val = calculateBoilingPointAtPressureSecant(
                        thermData.antoine,
                        P_system_Pa,
                        thermData.nbp_K || DEFAULT_INITIAL_T_K,
                        calculatePsat_Pa
                    );
                    // console.log(`Antoine BP fallback → ${input.name}: T_sat@${(P_system_Pa/1e5).toFixed(2)} bar = ${bp_at_Psys_K_val?.toFixed(2) ?? 'NaN'} K`);
                }

                return {
                    name: input.name,
                    casNumber: fetchedCasNumbers[index],
                    thermData: thermData,
                    originalIndex: index,
                    bp_at_Psys_K: bp_at_Psys_K_val,
                    isSCF: Number.isNaN(bp_at_Psys_K_val),
                };
            });

            // Validate pure component parameters based on selected fluid package
            if (fluidPackage === 'wilson' && processedComponents.some(pc => pc.thermData.V_L_m3mol == null || pc.thermData.V_L_m3mol <= 0)) {
                throw new Error("Wilson: Liquid molar volume (V_L_m3mol) missing or invalid.");
            }
            if (fluidPackage === 'unifac' && processedComponents.some(pc => !pc.thermData.unifacGroups || Object.keys(pc.thermData.unifacGroups).length === 0)) {
                throw new Error("UNIFAC: UNIFAC groups missing for one or more components.");
            }
            if ((fluidPackage === 'pr' || fluidPackage === 'srk') && processedComponents.some(pc => !pc.thermData.criticalProperties)) {
                throw new Error(`${fluidPackage.toUpperCase()}: Critical parameters (Tc, Pc, omega) missing.`);
            }
            if (fluidPackage === 'uniquac' && processedComponents.some(pc => !pc.thermData.uniquacParams)) {
                throw new Error("UNIQUAC: r and q parameters missing.");
            }
            if (processedComponents.some(pc => !pc.thermData.antoine)) {
                throw new Error("Antoine parameters missing for one or more components.");
            }

            // Log calculated BPs for components before sorting
            processedComponents.forEach(pc=>{
                // console.log('[BP] Component', pc.name, 'bp_at_Psys_K=', pc.bp_at_Psys_K?.toFixed(2), 'nbp_K=', pc.thermData.nbp_K?.toFixed(2));
            });

            const sortedComponents = [...processedComponents].sort((a, b) => {
                const bpA = a.bp_at_Psys_K ?? a.thermData.nbp_K ?? Infinity;
                const bpB = b.bp_at_Psys_K ?? b.thermData.nbp_K ?? Infinity;
                return bpA - bpB;
            });

            // console.log('[SORT] L/M/H order names:', sortedComponents.map(c=>c.name));

            setPlotSortedComponents(sortedComponents);

            // --- NEW: Directly compute and set axis titles here ---
            const formatBp = (bp_K: number | null | undefined) => 
                (bp_K && !Number.isNaN(bp_K) ? `${(bp_K - 273.15).toFixed(1)}°C` : 'N/A');

            const labelByOrig: Record<number,string> = {};
            if (sortedComponents.length === 3){
                labelByOrig[sortedComponents[0].originalIndex] = 'L';
                labelByOrig[sortedComponents[1].originalIndex] = 'M';
                labelByOrig[sortedComponents[2].originalIndex] = 'H';
            }
            const getBp = (orig:number)=> sortedComponents.find(c=>c.originalIndex===orig)?.bp_at_Psys_K ?? sortedComponents.find(c=>c.originalIndex===orig)?.thermData.nbp_K;

            setPlotAxisTitles({
                c: `${componentsInput[0]?.name || 'Comp 1'} (${labelByOrig[0] || ''}${sortedComponents.length===3?`, ${formatBp(getBp(0))}`:''})`,
                a: `${componentsInput[1]?.name || 'Comp 2'} (${labelByOrig[1] || ''}${sortedComponents.length===3?`, ${formatBp(getBp(1))}`:''})`,
                b: `${componentsInput[2]?.name || 'Comp 3'} (${labelByOrig[2] || ''}${sortedComponents.length===3?`, ${formatBp(getBp(2))}`:''})`
            });

            // Supercritical component check
            const scfComponents = sortedComponents.filter(pc => pc.isSCF);

            if (scfComponents.length > 0) {
                const numScf = scfComponents.length;
                let scfNamesWithValues = '';
                if (numScf === 1) {
                    const comp = scfComponents[0];
                    const pcBar = comp.thermData.criticalProperties ? (comp.thermData.criticalProperties.Pc_Pa / 1e5).toFixed(2) : 'N/A';
                    scfNamesWithValues = `${comp.name} (Pc=${pcBar} bar)`;
                    setScfMessage(`${scfNamesWithValues} is above its Pc at ${systemPressure} bar.\nA ternary residue-curve diagram cannot be generated under super-critical conditions.`);
                } else {
                    scfNamesWithValues = scfComponents.map((comp, index) => {
                        const pcBar = comp.thermData.criticalProperties ? (comp.thermData.criticalProperties.Pc_Pa / 1e5).toFixed(2) : 'N/A';
                        let nameStr = `${comp.name} (Pc=${pcBar} bar)`;
                        if (numScf > 1 && index === numScf - 1) {
                            nameStr = `and ${nameStr}`;
                        }
                        return nameStr;
                    }).join(numScf > 2 ? ', ' : ' ');
                    
                    setScfMessage(`${scfNamesWithValues} are above their Pc at ${systemPressure} bar.\nA ternary residue-curve diagram cannot be generated under super-critical conditions.`);
                }
                setIsLoading(false);
                setError("Supercritical component detected.");
                setResidueCurves([]);
                setDirectAzeotropes([]);
                return;
            } else {
                setScfMessage(null);
            }

            // Important: create a shallow copy before sorting by originalIndex so we don't mutate `sortedComponents`.
            const compoundsForBackend: CompoundData[] = [...sortedComponents]
                .sort((a,b) => a.originalIndex - b.originalIndex)
                .map(pc => ({
                    name: pc.name,
                    cas_number: pc.casNumber,
                    antoine: pc.thermData.antoine!,
                    wilsonParams: fluidPackage === 'wilson' ? { V_L_m3mol: pc.thermData.V_L_m3mol! } as WilsonPureParams : undefined,
                    unifacGroups: fluidPackage === 'unifac' ? pc.thermData.unifacGroups! : undefined,
                    prParams: fluidPackage === 'pr' ? pc.thermData.criticalProperties! as PrPureComponentParams : undefined,
                    srkParams: fluidPackage === 'srk' ? pc.thermData.criticalProperties! as SrkPureComponentParams : undefined,
                    uniquacParams: fluidPackage === 'uniquac' ? pc.thermData.uniquacParams! : undefined,
            }));

            const cas = compoundsForBackend.map(c => c.cas_number!);
            let activityModelParams: any;

            if (fluidPackage === 'wilson') {
                const fetchWilsonBothWays = async (c1:string,c2:string) => {
                    let p = await fetchWilsonInteractionParams(supabase, c1, c2);
                    if ((p.a12_J_mol === 0 && p.a21_J_mol === 0)) {
                        // Try the opposite orientation – many tables are asymmetric.
                        const pRev = await fetchWilsonInteractionParams(supabase, c2, c1);
                        if (pRev.a12_J_mol !== 0 || pRev.a21_J_mol !== 0) {
                            p = { a12_J_mol: pRev.a21_J_mol, a21_J_mol: pRev.a12_J_mol };
                        }
                    }
                    return p;
                };

                const params01 = await fetchWilsonBothWays(cas[0], cas[1]);
                const params02 = await fetchWilsonBothWays(cas[0], cas[2]);
                const params12 = await fetchWilsonBothWays(cas[1], cas[2]);

                const symmetric = (a_forward:number|undefined, a_reverse:number|undefined):[number,number] => {
                    const f = (a_forward !== undefined && !Number.isNaN(a_forward)) ? a_forward : undefined;
                    const r = (a_reverse !== undefined && !Number.isNaN(a_reverse)) ? a_reverse : undefined;
                    if (f !== undefined && r !== undefined) return [f,r];
                    if (f !== undefined) return [f,f];
                    if (r !== undefined) return [r,r];
                    return [0,0];
                };

                const [a01,a10] = symmetric(params01.a12_J_mol, params01.a21_J_mol);
                const [a02,a20] = symmetric(params02.a12_J_mol, params02.a21_J_mol);
                const [a12,a21] = symmetric(params12.a12_J_mol, params12.a21_J_mol);

                activityModelParams = {
                    a01_J_mol: a01,
                    a10_J_mol: a10,
                    a02_J_mol: a02,
                    a20_J_mol: a20,
                    a12_J_mol: a12,
                    a21_J_mol: a21,
                } as TernaryWilsonParams;
            } else if (fluidPackage === 'unifac') {
                const allSubgroupIds = new Set<number>();
                compoundsForBackend.forEach(comp => {
                    if (comp.unifacGroups) {
                        Object.keys(comp.unifacGroups).forEach(id => allSubgroupIds.add(parseInt(id)));
                    }
                });
                if (allSubgroupIds.size === 0 && compoundsForBackend.some(c => !c.unifacGroups || Object.keys(c.unifacGroups).length === 0)) {
                     throw new Error("UNIFAC selected, but no UNIFAC subgroup IDs found for any component after fetching data.");
                }
                activityModelParams = await fetchUnifacInteractionParams(supabase, Array.from(allSubgroupIds)) as UnifacParameters;
                if (!activityModelParams) {
                    throw new Error("Failed to fetch UNIFAC interaction parameters. The parameters object is null.");
                }
            } else if (fluidPackage === 'nrtl') {
                const params01 = await fetchNrtlParameters(supabase, cas[0], cas[1]);
                const params02 = await fetchNrtlParameters(supabase, cas[0], cas[2]);
                const params12 = await fetchNrtlParameters(supabase, cas[1], cas[2]);
                if (!params01 || !params02 || !params12) {
                    throw new Error("NRTL parameters for one or more binary pairs could not be fetched.");
                }
                activityModelParams = {
                    g01_J_mol: (params01 as any).g_ij_J_mol ?? 0, g10_J_mol: (params01 as any).g_ji_J_mol ?? 0, alpha01: (params01 as any).alpha_ij ?? 0,
                    g02_J_mol: (params02 as any).g_ij_J_mol ?? 0, g20_J_mol: (params02 as any).g_ji_J_mol ?? 0, alpha02: (params02 as any).alpha_ij ?? 0,
                    g12_J_mol: (params12 as any).g_ij_J_mol ?? 0, g21_J_mol: (params12 as any).g_ji_J_mol ?? 0, alpha12: (params12 as any).alpha_ij ?? 0,
                } as TernaryNrtlParams;
            } else if (fluidPackage === 'pr') {
                const params01 = await fetchPrInteractionParams(supabase, cas[0], cas[1]); 
                const params02 = await fetchPrInteractionParams(supabase, cas[0], cas[2]);
                const params12 = await fetchPrInteractionParams(supabase, cas[1], cas[2]);
                activityModelParams = {
                    k01: params01.k_ij ?? 0, k10: params01.k_ji ?? 0,
                    k02: params02.k_ij ?? 0, k20: params02.k_ji ?? 0,
                    k12: params12.k_ij ?? 0, k21: params12.k_ji ?? 0,
                } as TernaryPrParams;
            } else if (fluidPackage === 'srk') {
                const params01 = await fetchSrkInteractionParams(supabase, cas[0], cas[1]); 
                const params02 = await fetchSrkInteractionParams(supabase, cas[0], cas[2]);
                const params12 = await fetchSrkInteractionParams(supabase, cas[1], cas[2]);
                activityModelParams = {
                    k01: params01.k_ij ?? 0, k10: params01.k_ji ?? 0,
                    k02: params02.k_ij ?? 0, k20: params02.k_ji ?? 0,
                    k12: params12.k_ij ?? 0, k21: params12.k_ji ?? 0,
                } as TernarySrkParams;
            } else if (fluidPackage === 'uniquac') {
                const params01 = await fetchUniquacInteractionParams(supabase, cas[0], cas[1]);
                const params02 = await fetchUniquacInteractionParams(supabase, cas[0], cas[2]);
                const params12 = await fetchUniquacInteractionParams(supabase, cas[1], cas[2]);
                if (!params01 || !params02 || !params12) {
                    throw new Error("UNIQUAC parameters for one or more binary pairs could not be fetched.");
                }
                activityModelParams = {
                    a01_J_mol: ((params01 as any).a_ij_K ?? 0) * R_gas_const_J_molK, a10_J_mol: ((params01 as any).a_ji_K ?? 0) * R_gas_const_J_molK,
                    a02_J_mol: ((params02 as any).a_ij_K ?? 0) * R_gas_const_J_molK, a20_J_mol: ((params02 as any).a_ji_K ?? 0) * R_gas_const_J_molK,
                    a12_J_mol: ((params12 as any).a_ij_K ?? 0) * R_gas_const_J_molK, a21_J_mol: ((params12 as any).a_ji_K ?? 0) * R_gas_const_J_molK,
                } as TernaryUniquacParams;
            } else {
                throw new Error(`Unsupported fluid package: ${fluidPackage}`);
            }
            
            const avgNBP_at_Psys = sortedComponents
                .map(p => p.bp_at_Psys_K) // Use BP at system pressure
                .filter(bp => bp !== null && bp !== undefined) as number[];
            
            const initialTempGuess_K = avgNBP_at_Psys.length > 0 
                ? avgNBP_at_Psys.reduce((sum, val) => sum + val, 0) / avgNBP_at_Psys.length 
                : (sortedComponents
                    .map(p => p.thermData.nbp_K)
                    .filter(nbp => nbp !== null) as number[]
                  ).length > 0 
                    ? (sortedComponents
                        .map(p => p.thermData.nbp_K)
                        .filter(nbp => nbp !== null) as number[]
                      ).reduce((sum, val) => sum + val, 0) / (sortedComponents
                        .map(p => p.thermData.nbp_K)
                        .filter(nbp => nbp !== null) as number[]
                      ).length
                    : DEFAULT_INITIAL_T_K;

            // -------------------------------------------------------------
            // Helper: convert a composition expressed in Light/Medium/Heavy
            // (i.e. sortedComponents order) back to the user-input order that
            // `compoundsForBackend` and all thermodynamic routines expect.
            // -------------------------------------------------------------
            const toInputOrder = (xSorted: number[]): number[] => {
                const xInput = [0, 0, 0];
                sortedComponents.forEach((comp, sIdx) => {
                    xInput[comp.originalIndex] = xSorted[sIdx] ?? 0;
                });
                // console.debug('[toInputOrder] L/M/H -> input order:', xSorted.map(v=>v.toFixed(3)), '=>', xInput.map(v=>v.toFixed(3)));
                return xInput;
            };

            // Placeholder – will be populated later after identifying unstable nodes
            let startingCompositions: number[][] = [];

            // --- Integrate residue curves from each seed ---
            const curvePromises: Promise<ResidueCurve | null>[] = [];

            const d_xi_step_for_model =
                (fluidPackage === 'pr' || fluidPackage === 'srk') ? 0.002 : 0.03;
            const max_steps = (fluidPackage === 'unifac') ? 600 :
                              (fluidPackage === 'pr' || fluidPackage === 'srk') ? 1200 : 800;

            for (const start_x of startingCompositions) {
                switch (fluidPackage) {
                    case 'wilson':
                        curvePromises.push(simulateResidueCurveODE('wilson', start_x, P_system_Pa, compoundsForBackend, activityModelParams as TernaryWilsonParams, d_xi_step_for_model, max_steps, initialTempGuess_K));
                        break;
                    case 'unifac':
                        curvePromises.push(simulateResidueCurveODE('unifac', start_x, P_system_Pa, compoundsForBackend, activityModelParams as UnifacParameters, d_xi_step_for_model, max_steps, initialTempGuess_K));
                        break;
                    case 'nrtl':
                        curvePromises.push(simulateResidueCurveODE('nrtl', start_x, P_system_Pa, compoundsForBackend, activityModelParams as TernaryNrtlParams, d_xi_step_for_model, max_steps, initialTempGuess_K));
                        break;
                    case 'pr':
                        curvePromises.push(simulateResidueCurveODE('pr', start_x, P_system_Pa, compoundsForBackend, activityModelParams as TernaryPrParams, d_xi_step_for_model, max_steps, initialTempGuess_K));
                        break;
                    case 'srk':
                        curvePromises.push(simulateResidueCurveODE('srk', start_x, P_system_Pa, compoundsForBackend, activityModelParams as TernarySrkParams, d_xi_step_for_model, max_steps, initialTempGuess_K));
                        break;
                    case 'uniquac':
                        curvePromises.push(simulateResidueCurveODE('uniquac', start_x, P_system_Pa, compoundsForBackend, activityModelParams as TernaryUniquacParams, d_xi_step_for_model, max_steps, initialTempGuess_K));
                        break;
                    default:
                        break;
                }
            }

            const curveResults = await Promise.all(curvePromises);
            const validCurves = curveResults.filter(c => c && c.length > 1) as ResidueCurve[];
            setResidueCurves(validCurves);

            // --- Call Direct Azeotrope Solver using Hybrid RCM Approach ---
            let foundAzeotropes: AzeotropeDisplayInfo[] = [];
            let binaryHits: AzeotropeResult[] = [];
            const BINARY_ERR_TOL = 0.01; // global threshold for binary acceptance

            if (compoundsForBackend.length === 3) {
                // --- STEP 1: Generate a diverse set of initial guesses ---

                // a) Start with the original hardcoded guesses for broad coverage
                const initialGuessPool: { x: number[], T_K: number }[] = [
                    { x: [0.333, 0.333, 0.334], T_K: initialTempGuess_K }, // Center guess
                ];
                
                // Add other binary/vertex guesses from the original code
                if (sortedComponents.length === 3) {
                    initialGuessPool.push({ x: toInputOrder([0.9, 0.05, 0.05]), T_K: sortedComponents[0].bp_at_Psys_K || sortedComponents[0].thermData.nbp_K || initialTempGuess_K });
                    initialGuessPool.push({ x: toInputOrder([0.05, 0.9, 0.05]), T_K: sortedComponents[1].bp_at_Psys_K || sortedComponents[1].thermData.nbp_K || initialTempGuess_K });
                    initialGuessPool.push({ x: toInputOrder([0.05, 0.05, 0.9]), T_K: sortedComponents[2].bp_at_Psys_K || sortedComponents[2].thermData.nbp_K || initialTempGuess_K });
                }
                if (compoundsForBackend.length === 3) {
                    const T0 = sortedComponents.find(c => c.originalIndex === 0)?.bp_at_Psys_K || sortedComponents.find(c => c.originalIndex === 0)?.thermData.nbp_K || initialTempGuess_K;
                    const T1 = sortedComponents.find(c => c.originalIndex === 1)?.bp_at_Psys_K || sortedComponents.find(c => c.originalIndex === 1)?.thermData.nbp_K || initialTempGuess_K;
                    const T2 = sortedComponents.find(c => c.originalIndex === 2)?.bp_at_Psys_K || sortedComponents.find(c => c.originalIndex === 2)?.thermData.nbp_K || initialTempGuess_K;
                    initialGuessPool.push({ x: toInputOrder([0.5, 0.5, 0.001]), T_K: (T0+T1)/2 });
                    initialGuessPool.push({ x: toInputOrder([0.5, 0.001, 0.5]), T_K: (T0+T2)/2 });
                    initialGuessPool.push({ x: toInputOrder([0.001, 0.5, 0.5]), T_K: (T1+T2)/2 });
                }

                // d) Systematic Binary Search – derive additional initial guesses
                try {
                  // Run binary search in Light/Medium/Heavy order to make results
                  // independent of user input order, then map compositions back to
                  // the original order with `toInputOrder`.
                  const compsForBinary = [...sortedComponents]
                        .map((sc, idx) => {
                            const found = compoundsForBackend.find(c => c.cas_number === sc.casNumber)!;
                            // console.log(`[MapLBH->Binary] idx ${idx} LBH=${sc.name} maps to backend pos=${compoundsForBackend.indexOf(found)} name=${found.name}`);
                            return found;
                        }) as CompoundData[];

                  // Re-map Wilson interaction parameters to match L/M/H order so they
                  // stay consistent with `compsForBinary`.
                  let paramsForBinary = activityModelParams;
                  if (fluidPackage === 'wilson') {
                      const perm: number[] = sortedComponents.map(sc => sc.originalIndex); // newIndex -> origIndex
                      const a = (i:number,j:number) => {
                        const pair = `${i}${j}`;
                        switch(pair){
                          case '01': return activityModelParams.a01_J_mol;
                          case '10': return activityModelParams.a10_J_mol;
                          case '02': return activityModelParams.a02_J_mol;
                          case '20': return activityModelParams.a20_J_mol;
                          case '12': return activityModelParams.a12_J_mol;
                          case '21': return activityModelParams.a21_J_mol;
                          default: return 0;
                        }
                      };
                      // If a specific orientation value is missing (0), fall back to the reverse
                      const aSafe = (i:number,j:number):number => {
                        const val = a(i,j);
                        if (val !== 0 && !Number.isNaN(val)) return val;
                        const rev = a(j,i);
                        return (rev !== 0 && !Number.isNaN(rev)) ? rev : 0;
                      };
                      const aLM = (i:number,j:number)=> aSafe(perm[i],perm[j]);
                      paramsForBinary = {
                        a01_J_mol: aLM(0,1), a10_J_mol: aLM(1,0),
                        a02_J_mol: aLM(0,2), a20_J_mol: aLM(2,0),
                        a12_J_mol: aLM(1,2), a21_J_mol: aLM(2,1)
                      } as TernaryWilsonParams;
                      // console.log('[ParamReorder] Wilson params remapped for L/M/H order');
                  }

                  const binaryHitsSorted = systematicBinaryAzeotropeSearch(
                        generatingFluidPackage as FluidPackageResidue, P_system_Pa, compsForBinary, paramsForBinary, initialTempGuess_K);

                  binaryHits = binaryHitsSorted.map(h => ({ ...h, x: toInputOrder(h.x) }));

                  // console.log(`[BinarySearch] ${binaryHits.length} binary azeotrope candidates found (order-independent).`);
                  // console.log('[BinarySearch] compsForBinary order (names):', compsForBinary.map(c=>c.name));
                  for (const hit of binaryHits) {
                    if (hit.converged && hit.errorNorm !== undefined && hit.errorNorm < BINARY_ERR_TOL) {
                      // add as starting guess
                      initialGuessPool.push({ x: hit.x, T_K: hit.T_K });
                      // also record directly for display so that even if the
                      // subsequent ternary solver fails we still list the
                      // binary azeotrope detected on the edge.
                      foundAzeotropes.push({ x: [...hit.x], T_K: hit.T_K, errorNorm: hit.errorNorm });
                    }
                    // console.log(`[BinaryHit ${binaryHits.indexOf(hit)}] x(input order)=`, hit.x.map(v=>v.toFixed(3)), 'errorNorm=', hit.errorNorm?.toExponential(2));
                  }
                  const ternarySeed = estimateTernaryGuessFromBinary(binaryHits);
                  if (ternarySeed) {
                    initialGuessPool.push({ x: ternarySeed.x, T_K: ternarySeed.T_K }); // ternarySeed.x already in input order
                    // console.log('[BinarySearch] Added averaged ternary seed.');
                  }
                } catch (err) {
                  console.warn('[BinarySearch] Error during systematic scan:', err);
                }

                // --- NEW systematic interior ternary scan ---
                try {
                  const ternaryHits = systematicTernaryAzeotropeSearch(generatingFluidPackage as FluidPackageResidue, P_system_Pa, compoundsForBackend, activityModelParams, initialTempGuess_K);
                  // console.log(`[TernarySearch] ${ternaryHits.length} interior ternary azeotrope candidates found.`);
                  for (const th of ternaryHits) {
                    if (th.converged) {
                      initialGuessPool.push({ x: th.x, T_K: th.T_K });
                    }
                  }
                } catch (err) {
                  console.warn('[TernarySearch] Error during interior scan:', err);
                }

                // b) **NEW**: Use RCM topology to find high-quality, targeted guesses
                const approxAzeotropesFromRCM = findApproximateAzeotropesFromCurves(validCurves);
                // console.log(`[Azeotrope Solver] Found ${approxAzeotropesFromRCM.length} potential azeotropes from RCM topology.`);

                for (const approx_x of approxAzeotropesFromRCM) {
                    // Use the temperature from the curve start/end point as the initial T guess
                    let T_guess_for_x = initialTempGuess_K;
                    let found_T_guess = false;
                    for (const curve of validCurves) {
                        const firstPoint = curve[0];
                        const lastPoint = curve[curve.length - 1];
                        
                        // Check if the guess matches the start point (potential min-boiling azeotrope)
                        if (firstPoint && !found_T_guess && firstPoint.x.every((val, i) => Math.abs(val - approx_x[i]) < 1e-4)) {
                            T_guess_for_x = firstPoint.T_K;
                            found_T_guess = true;
                        }
                        
                        // Check if the guess matches the end point (potential max-boiling azeotrope)
                        if (lastPoint && !found_T_guess && lastPoint.x.every((val, i) => Math.abs(val - approx_x[i]) < 1e-4)) {
                            T_guess_for_x = lastPoint.T_K;
                            found_T_guess = true;
                        }

                        if (found_T_guess) break; // Found a good T, move to next azeotrope guess
                    }
                    initialGuessPool.push({ x: approx_x, T_K: T_guess_for_x });
                }

                // c) De-duplicate the final list of guesses to avoid redundant solver calls
                const uniqueGuesses: { x: number[], T_K: number }[] = [];
                const MIN_SQ_DIST_GUESSES = 0.005 * 0.005;

                for (const guess of initialGuessPool) {
                    if (!uniqueGuesses.some(uniqueGuess => {
                        const distSq = (guess.x[0] - uniqueGuess.x[0])**2 + (guess.x[1] - uniqueGuess.x[1])**2 + (guess.x[2] - uniqueGuess.x[2])**2;
                        return distSq < MIN_SQ_DIST_GUESSES;
                    })) {
                        uniqueGuesses.push(guess);
                    }
                }
                // console.log(`[Azeotrope Solver] Refining ${uniqueGuesses.length} unique initial guesses.`);
                // uniqueGuesses.forEach((g,idx)=>console.log(`[UniqueGuess ${idx}] x=`, g.x.map(v=>v.toFixed(3)), 'T=', (g.T_K-273.15).toFixed(1)));

                // --- STEP 2: Run the local solver with the improved guess list ---
                let results: (AzeotropeResult | null)[] = [];

                // --- Debugging: validate parameter completeness ---
                // --- End Debugging validation block ---

                const solverPromises = uniqueGuesses.map(guess => { // Use the enhanced 'uniqueGuesses' list
                    switch (generatingFluidPackage) {
                        case 'wilson':
                            return findAzeotropeWilson(P_system_Pa, compoundsForBackend, activityModelParams as TernaryWilsonParams, guess.x, guess.T_K);
                        case 'nrtl':
                            return findAzeotropeNRTL(P_system_Pa, compoundsForBackend, activityModelParams as TernaryNrtlParams, guess.x, guess.T_K);
                        case 'unifac':
                             return findAzeotropeUnifac(P_system_Pa, compoundsForBackend, activityModelParams as UnifacParameters, guess.x, guess.T_K);
                        case 'pr':
                            return findAzeotropePr(P_system_Pa, compoundsForBackend, activityModelParams as TernaryPrParams, guess.x, guess.T_K);
                        case 'srk':
                            return findAzeotropeSrk(P_system_Pa, compoundsForBackend, activityModelParams as TernarySrkParams, guess.x, guess.T_K);
                        case 'uniquac':
                            return findAzeotropeUniquac(P_system_Pa, compoundsForBackend, activityModelParams as TernaryUniquacParams, guess.x, guess.T_K);
                        default:
                            return Promise.resolve(null);
                    }
                });
                
                results = await Promise.all(solverPromises);
                
                // console.log(`[Azeotrope Solver] Raw results for ${generatingFluidPackage}:`, JSON.parse(JSON.stringify(results)));

                const convergedResults = results.filter(r => {
                    if (!r) return false;
                    if (generatingFluidPackage === 'srk') {
                        const srkRes = r as SrkAzeotropeResult;
                        // Ensure y exists and is an array of 3 numbers for SRK
                        if (!srkRes.y || srkRes.y.length !== 3 || !srkRes.x || srkRes.x.length !== 3) {
                            // console.warn("SRK Azeotrope result is missing x or y, or they are not 3-component arrays. Skipping.", srkRes);
                            return false;
                        }
                        const errorNormYX = Math.sqrt(
                            (srkRes.y[0] - srkRes.x[0])**2 +
                            (srkRes.y[1] - srkRes.x[1])**2 +
                            (srkRes.y[2] - srkRes.x[2])**2
                        );
                        (r as any).calculatedErrorNormYX = errorNormYX;
                        return errorNormYX < 0.01; // Threshold for considering SRK result an azeotrope
                    }
                    // For all other models, they should have a 'converged' property
                    const convFlag = (r as any).converged === true;
                    const errVal: number | undefined = (r as any).errorNorm;
                    return convFlag || (errVal !== undefined && errVal < BINARY_ERR_TOL);
                }) as Array<AzeotropeResult & { calculatedErrorNormYX?: number }>;

                // console.log(`[Azeotrope Solver] Converged results:`, JSON.parse(JSON.stringify(convergedResults)));

                const uniqueDirectAzeotropes: AzeotropeDisplayInfo[] = [];
                // Increase tolerance to merge near-identical solutions (binary vs ternary duplicates)
                const minSqDistAzeo = 0.09 * 0.09;

                for (const res of convergedResults) {
                    // Ensure x is a 3-component array. SRK result should already be.
                    if (!res.x || res.x.length !== 3) {
                        // console.warn("Azeotrope result 'x' is not a 3-component array. Skipping.", res);
                        continue;
                    }

                    // --- NEW: Filter out azeotropes too close to a pure component corner ---
                    const PURE_COMPONENT_CORNER_THRESHOLD = 0.98;
                    if (res.x.some(xi => xi > PURE_COMPONENT_CORNER_THRESHOLD)) {
                        // console.log(`[Azeotrope Filter] Rejecting azeotrope near corner: x=[${res.x.map(v=>v.toFixed(3)).join(',')}]`);
                        continue; // Skip this result as it's essentially a pure component
                    }

                    let displayErrorNorm: number | undefined;
                    let res_T_K = res.T_K;

                    if (generatingFluidPackage === 'srk') {
                        displayErrorNorm = (res as any).calculatedErrorNormYX;
                        // If T_K from SRK result is 0 or invalid, use initialTempGuess_K as a fallback.
                        // This is a workaround; findAzeotropeSrk should ideally return a correct T_K.
                        if (!res_T_K || res_T_K <= 0) {
                            // console.warn(`SRK Azeotrope result has T_K = ${res_T_K}. Using initialTempGuess_K (${initialTempGuess_K.toFixed(2)} K) as fallback.`);
                            res_T_K = initialTempGuess_K; 
                        }
                    } else {
                        // For all other models
                        displayErrorNorm = (res as AzeotropeResult).errorNorm;
                    }
                    
                    let isTooClose = false;
                    for (const uniqueAzeo of uniqueDirectAzeotropes) {
                        const distSq = (res.x[0] - uniqueAzeo.x[0]) ** 2 +
                                       (res.x[1] - uniqueAzeo.x[1]) ** 2 +
                                       (res.x[2] - uniqueAzeo.x[2]) ** 2;
                        if (distSq < minSqDistAzeo) {
                            isTooClose = true;
                            // If this result is better (smaller error norm), replace the existing one
                            if (displayErrorNorm !== undefined && (uniqueAzeo.errorNorm === undefined || displayErrorNorm < uniqueAzeo.errorNorm)) {
                                 uniqueAzeo.x = [...res.x];
                                 uniqueAzeo.T_K = res_T_K; // Use T_K from the current result
                                 uniqueAzeo.errorNorm = displayErrorNorm;
                            }
                            break;
                        }
                    }
                    if (!isTooClose) {
                        uniqueDirectAzeotropes.push({
                            x: [...res.x],
                            T_K: res_T_K, // Use T_K from the current result
                            errorNorm: displayErrorNorm
                        });
                    }
                }
                // console.log(`[Azeotrope Solver] Unique azeotropes for state:`, JSON.parse(JSON.stringify(uniqueDirectAzeotropes)));

                // --- NEW: Merge original binary hits that may have been lost if the
                //     local ternary solver did not re-converge on the edge.
                const mergedAzeos: AzeotropeDisplayInfo[] = [...uniqueDirectAzeotropes];
                // Use same tolerance for merging edge hits to avoid duplicates
                const EDGE_MIN_DIST_SQ = minSqDistAzeo;

                for(const bh of binaryHits){
                    if(!bh.converged) continue;
                    if(bh.errorNorm === undefined || bh.errorNorm >= BINARY_ERR_TOL) continue; // skip if not precise enough
                    // skip if third component exceeds 0.02 (should already be small)
                    if(bh.x[0] < 0.02 && bh.x[1] < 0.02 && bh.x[2] < 0.02) continue; // safety
                    const exists = mergedAzeos.some(az=>{
                        const dSq = (az.x[0]-bh.x[0])**2 + (az.x[1]-bh.x[1])**2 + (az.x[2]-bh.x[2])**2;
                        return dSq < EDGE_MIN_DIST_SQ;
                    });
                    if(!exists){
                        mergedAzeos.push({ x:[...bh.x], T_K: bh.T_K, errorNorm: bh.errorNorm });
                    }
                }

                foundAzeotropes = mergedAzeos;
            }
            
            setDirectAzeotropes(foundAzeotropes);
            if (foundAzeotropes.length > 0) {
                // console.log(`Directly found ${foundAzeotropes.length} potential azeotrope(s) using ${generatingFluidPackage}:`, foundAzeotropes.map(az => ({x: az.x, T_K: az.T_K, err: az.errorNorm})));
            } else {
                // console.log(`No azeotropes found directly using ${generatingFluidPackage} with the given initial guesses.`);
            }
            // --- End Direct Azeotrope Solver Call ---

            // store backend info for eigenvalue classification
            setBackendComps(compoundsForBackend);
            setBackendPkgParams(activityModelParams);
            setBackendPressurePa(P_system_Pa);

            // --- NEW OCCUPANCY GRID FILL LOGIC ---

            // 1. Generate a very dense, shuffled grid of candidate seeds
            const candidateSeeds = generateStartingPoints(10, 12).map(toInputOrder);
            shuffleArray(candidateSeeds);

            const allCurves: ResidueCurve[] = [];
            const GRID_SIZE = 120; // Finer grid for more precise checking
            const occupancyGrid = Array(GRID_SIZE + 1).fill(0).map(() => Array(GRID_SIZE + 1).fill(false));
            const Y_SCALE = 1 / (Math.sqrt(3) / 2); // Pre-calculate for normalization

            const d_xi_step = (fluidPackage === 'pr' || fluidPackage === 'srk') ? 0.002 : 0.03;
            const m_steps = (fluidPackage === 'unifac') ? 600 : (fluidPackage === 'pr' || fluidPackage === 'srk') ? 1200 : 800;

            const getGridCoords = (ternaryPoint: number[]): { j: number, i: number } | null => {
                if (sortedComponents.length !== 3) return null;
                const x_light = ternaryPoint[sortedComponents[0].originalIndex];
                const x_intermediate = ternaryPoint[sortedComponents[1].originalIndex];
                
                const x_paper = x_light + x_intermediate * 0.5;
                const y_paper = x_intermediate * (Math.sqrt(3) / 2);

                const j = Math.floor(x_paper * GRID_SIZE);
                const i = Math.floor(y_paper * Y_SCALE * GRID_SIZE);

                if (i >= 0 && i <= GRID_SIZE && j >= 0 && j <= GRID_SIZE) {
                    return { j, i };
                }
                return null;
            }

            // 2. Iteratively process shuffled seeds, adding curves to empty grid cells
            for (const seed of candidateSeeds) {
                const coords = getGridCoords(seed);
                if (!coords || occupancyGrid[coords.i][coords.j]) {
                    continue; // Skip if outside grid or cell is occupied
                }

                const newCurve = await simulateResidueCurveODE(
                    generatingFluidPackage,
                    seed,
                    P_system_Pa,
                    compoundsForBackend,
                    activityModelParams,
                    d_xi_step,
                    m_steps,
                    initialTempGuess_K
                );

                if (newCurve && newCurve.length > 1) {
                    allCurves.push(newCurve);
                    // Mark all cells this new curve passes through as occupied
                    for (const p of newCurve) {
                        const curveCoords = getGridCoords(p.x);
                        if (curveCoords) {
                            occupancyGrid[curveCoords.i][curveCoords.j] = true;
                        }
                    }
                }
            }
            setResidueCurves(allCurves);

            // store backend info for eigenvalue classification
            setBackendComps(compoundsForBackend);
            setBackendPkgParams(activityModelParams);
            setBackendPressurePa(P_system_Pa);

        } catch (err: any) {
            console.error("Error generating map:", err);
            setError(err.message || "An unknown error occurred.");
        } finally {
            setIsLoading(false);
            setDisplayedFluidPackage(generatingFluidPackage); // Update displayed fluid package when generation finishes
            setDisplayedPressure(systemPressure);
            setDisplayedComponentNames(componentsInput.map(ci=>ci.name));
        }
    }, [componentsInput, systemPressure, supabase, fetchCasNumberByName, fetchAndConvertThermData, fluidPackage]); // Removed initialTempGuess_K

    useEffect(() => {
        if (supabase && !initialMapGenerated.current) {
            handleGenerateMap();
            initialMapGenerated.current = true;
        }
    }, [supabase, handleGenerateMap]);

    useEffect(() => {
        if (!residueCurves || residueCurves.length === 0) {
            setCleanedResidueCurves([]);
            // console.log(`Plotting useEffect: No residueCurves or empty. Setting cleanedResidueCurves to [].`);
            return;
        }

        // Use a much tighter threshold for PR and SRK (small steps) to keep more points for line drawing.
        // Other models use a more standard threshold.
        const moleFractionThreshold = (fluidPackage === 'pr' || fluidPackage === 'srk') ? 1e-9 : 1e-5; 

        const filterCurves = (curves: ResidueCurve[]): ResidueCurve[] => {
            return curves.map(curve => {
                if (curve.length < 2) return curve; // Keep short curves (0 or 1 point) as is

                const cleanedCurve: ResidueCurve = [curve[0]]; // Always include the first point
                for (let i = 1; i < curve.length; i++) {
                    const currentPoint = curve[i];
                    const lastCleanedPoint = cleanedCurve[cleanedCurve.length - 1];
                    
                    let diffSq = 0;
                    for (let j = 0; j < currentPoint.x.length; j++) {
                        const d = currentPoint.x[j] - lastCleanedPoint.x[j];
                        diffSq += d * d;
                    }

                    if (Math.sqrt(diffSq) > moleFractionThreshold) {
                        cleanedCurve.push(currentPoint);
                    }
                }
                return cleanedCurve;
            }).filter(curve => curve.length > 0); // Ensure no empty curves if all points were identical
        };

        const filteredForPlotting = filterCurves(residueCurves);
        // Verbose plotting logs removed for production.

        if (residueCurves.length > 0 && filteredForPlotting.length === 0) {
            console.warn(`Plotting useEffect: All curves were filtered out by the ${moleFractionThreshold.toExponential()} threshold. Original curves might be too short or stagnant.`);
        }
        // Removed verbose curve preview log.

        setCleanedResidueCurves(filteredForPlotting);
    }, [residueCurves, fluidPackage]); // Ensure fluidPackage is a dependency

    const memoizedPlotData = useMemo(() => {
        const modelName = displayedFluidPackage.charAt(0).toUpperCase() + displayedFluidPackage.slice(1); // Use displayedFluidPackage

        let traces: Data[] = []; 

        // plotSortedComponents is [Lightest, Intermediate, Heaviest]
        // User's desired visual: Right=Lightest, Top=Intermediate, Left=Heaviest
        // User's interpretation of Plotly axes: a=Top, b=Left, c=Right

        // Titles now come directly from state
        const { a: titleA, b: titleB, c: titleC } = plotAxisTitles;
        
        const plotTitleColor = currentTheme === 'dark' ? '#e5e7eb' : '#1f2937'; 
        // Define font family string and color for consistent use
        const merriweatherFamilyString = 'Merriweather Sans, Arial, sans-serif';
        const plotFontColor = currentTheme === 'dark' ? '#e5e7eb' : '#1f2937';

        // Base font object for global settings and easy reuse
        const basePlotFontObject = { family: merriweatherFamilyString, color: plotFontColor };

        const axisTitleFont = { ...basePlotFontObject, size: 15 };
        const tickFont = { ...basePlotFontObject, size: 14 };

        const baseLayout: Partial<Layout> = {
            title: {
                text: `Ternary Residue Curve Map (ODE Simulation)`,
                font: { family: merriweatherFamilyString, size: 18, color: plotTitleColor }, 
            },
            ternary: {
                sum: 1,
                aaxis: { title: { text: titleA, font: axisTitleFont, standoff: 35 }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: 'transparent' }, // Top (Intermediate)
                baxis: { title: { text: titleB, font: axisTitleFont, standoff: 35 }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: 'transparent' }, // Left (Heaviest)
                caxis: { title: { text: titleC, font: axisTitleFont, standoff: 35 }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: 'transparent' }, // Right (Lightest)
                bgcolor: 'transparent',
            },
            margin: { l: 90, r: 90, b: 90, t: 110, pad: 10 }, // Increased margins for better text display
            autosize: true,
            paper_bgcolor: 'transparent', 
            font: basePlotFontObject, // Global font setting
            plot_bgcolor: 'transparent', 
            annotations: [], 
            shapes: [], 
        };

        // Get container dimensions for coordinate transformation
        let containerWidthPx = 800; // Default/fallback
        let containerHeightPx = 600; // Default/fallback
        if (plotContainerRef) { // Check if plotContainerRef (the div element) exists
            containerWidthPx = plotContainerRef.offsetWidth;
            containerHeightPx = plotContainerRef.offsetHeight;
        }

        const marginLeftPx = baseLayout.margin?.l ?? 70;
        const marginRightPx = baseLayout.margin?.r ?? 70;
        const marginBottomPx = baseLayout.margin?.b ?? 70;
        const marginTopPx = baseLayout.margin?.t ?? 90;

        // Plot area in pixel dimensions (inside margins)
        const plotAreaClientWidthPx = containerWidthPx - marginLeftPx - marginRightPx;
        const plotAreaClientHeightPx = containerHeightPx - marginTopPx - marginBottomPx;

        // Plot area origin and dimensions in paper units
        const paperMarginL = marginLeftPx / containerWidthPx;
        const paperMarginR = marginRightPx / containerWidthPx;
        const paperMarginB = marginBottomPx / containerHeightPx;
        const paperMarginT = marginTopPx / containerHeightPx;

        const plotAreaPaperOriginX = paperMarginL;
        const plotAreaPaperOriginY = paperMarginB;
        const plotAreaPaperWidth = 1.0 - paperMarginL - paperMarginR;
        const plotAreaPaperHeight = 1.0 - paperMarginT - paperMarginB;
        
        let finalTriOriginX_paper = plotAreaPaperOriginX;
        let finalTriOriginY_paper = plotAreaPaperOriginY;
        let actualRenderedTriangleBase_paper = plotAreaPaperWidth; 
        let actualRenderedTriangleHeight_paper = plotAreaPaperHeight; 

        // Declare triangleActualBasePx and triangleActualHeightPx here for wider scope
        let triangleActualBasePx: number = plotAreaClientWidthPx; 
        let triangleActualHeightPx: number = plotAreaClientHeightPx; 

        if (plotAreaClientWidthPx > 0 && plotAreaClientHeightPx > 0) {
            // Determine if width or height is the limiting dimension for an equilateral triangle
            if (plotAreaClientWidthPx * (Math.sqrt(3)/2) <= plotAreaClientHeightPx) {
                triangleActualBasePx = plotAreaClientWidthPx;
                triangleActualHeightPx = triangleActualBasePx * (Math.sqrt(3)/2);
            } else {
                triangleActualHeightPx = plotAreaClientHeightPx;
                triangleActualBasePx = triangleActualHeightPx / (Math.sqrt(3)/2);
            }

            actualRenderedTriangleBase_paper = triangleActualBasePx / containerWidthPx;
            actualRenderedTriangleHeight_paper = triangleActualHeightPx / containerHeightPx;

            const triangleOffsetX_paper = (plotAreaPaperWidth - actualRenderedTriangleBase_paper) / 2;
            const triangleOffsetY_paper = (plotAreaPaperHeight - actualRenderedTriangleHeight_paper) / 2;

            finalTriOriginX_paper = plotAreaPaperOriginX + triangleOffsetX_paper;
            finalTriOriginY_paper = plotAreaPaperOriginY + triangleOffsetY_paper;
        }

        if (cleanedResidueCurves.length === 0 || plotSortedComponents.length !== 3) { // Use cleanedResidueCurves
            const layoutUpdate: Partial<Layout> = {
                ...baseLayout, 
                title: {
                    ...(typeof baseLayout.title === 'object' ? baseLayout.title : {}), // Keep base title properties
                    text: `Ternary Residue Curve Map (${modelName} - No data or NBP info pending)` // Specific title for this case
                },
                annotations: [{ 
                    text: isLoading ? "Generating initial map..." : "Enter component data and click 'Generate Residue Map'.",
                    showarrow: false,
                    xref: 'paper',
                    yref: 'paper',
                    x: 0.5,
                    y: 0.5,
                    font: basePlotFontObject, // Use defined base font object
                }],
                shapes: [] 
            };
            return { baseTraces: [], baseLayout: layoutUpdate, minAz: [], maxAz: [], sadAz: [] };
        }

        traces = cleanedResidueCurves.map((curve, index) => {
            const compNameA = componentsInput[1]?.name || 'Comp 2'; // top axis
            const compNameB = componentsInput[2]?.name || 'Comp 3'; // left axis
            const compNameC = componentsInput[0]?.name || 'Comp 1'; // right axis

            const labelByOrig: Record<number,string> = {};
            if(plotSortedComponents.length===3){
                labelByOrig[plotSortedComponents[0].originalIndex] = 'L';
                labelByOrig[plotSortedComponents[1].originalIndex] = 'M';
                labelByOrig[plotSortedComponents[2].originalIndex] = 'H';
            }
            const labelA = labelByOrig[1] ?? '';
            const labelB = labelByOrig[2] ?? '';
            const labelC = labelByOrig[0] ?? '';

            return {
                type: 'scatterternary',
                mode: curve.length > 1 ? 'lines' : 'markers',
                a: curve.map(p => p.x[1]), // component 2
                b: curve.map(p => p.x[2]), // component 3
                c: curve.map(p => p.x[0]), // component 1
                name: `Curve ${index + 1}`,
                showlegend: false,
                line: { width: 1.5, color: '#60a5fa' }, // Example color
                hoverinfo: 'skip'
            };
        });

        const arrowMarkerTraceData = {
            a: [] as number[],
            b: [] as number[],
            c: [] as number[],
            angles: [] as number[],
            text: [] as string[], 
        };
        // --- Globals for arrow placement ---
        const allPlacedArrowIdealCoords: { x: number; y: number }[] = [];
        // ADJUST THESE:
        const minSqDistBetweenPlacedArrows_ideal = (0.025) * (0.025); // Smaller value for denser arrows
        const maxTotalArrowsOnPlot = 500;                             // Allow more arrows overall
        let placedArrowCountOverall = 0;


        let firstArrowLogged = false; 
        const targetArrowPixelSize = 10; 
        const arrowMarkerColor = currentTheme === 'dark' ? '#FFFFFF' : '#333333'; // Adjusted for theme

        const plotAreaPxWidth = containerWidthPx - (baseLayout.margin?.l ?? 70) - (baseLayout.margin?.r ?? 70);
        const plotAreaPxHeight = containerHeightPx - (baseLayout.margin?.t ?? 90) - (baseLayout.margin?.b ?? 70);

        if (plotSortedComponents.length === 3 && plotContainerRef && actualRenderedTriangleBase_paper > 0 && actualRenderedTriangleHeight_paper > 0 && plotAreaPxWidth > 0 && plotAreaPxHeight > 0) { // Check plotContainerRef
            cleanedResidueCurves.forEach((curve, curveIndex) => { 
                if (curve.length < 2) return;
                if (placedArrowCountOverall >= maxTotalArrowsOnPlot) return;

                let lastCheckedPointIdealCoords: { x: number; y: number } | null = null;
                const minIdealDistSqForCheckingNextPointOnCurve = (0.02) * (0.02); 

                for (let headIndex = 0; headIndex < curve.length; headIndex++) {
                    if (placedArrowCountOverall >= maxTotalArrowsOnPlot) break;

                    const currentPointTernary = curve[headIndex].x;
                    const currentPointIdealCoords = convertTernaryToPaperCoordinates(currentPointTernary, plotSortedComponents);

                    if (!currentPointIdealCoords) continue;

                    let shouldCheckThisPointForArrow = false;
                    if (!lastCheckedPointIdealCoords) {
                        if (headIndex >= Math.floor(curve.length * 0.05) && curve.length > 1) { 
                            shouldCheckThisPointForArrow = true;
                        }
                    } else {
                        const distSqFromLastCheck =
                            (currentPointIdealCoords.x - lastCheckedPointIdealCoords.x)**2 +
                            (currentPointIdealCoords.y - lastCheckedPointIdealCoords.y)**2;
                        if (distSqFromLastCheck >= minIdealDistSqForCheckingNextPointOnCurve) {
                            shouldCheckThisPointForArrow = true;
                        }
                    }

                    if (shouldCheckThisPointForArrow) {
                        lastCheckedPointIdealCoords = currentPointIdealCoords; 

                        let tooCloseToExistingArrow = false;
                        for (const placedCoord of allPlacedArrowIdealCoords) {
                            const distSqToPlaced = (currentPointIdealCoords.x - placedCoord.x)**2 +
                                                   (currentPointIdealCoords.y - placedCoord.y)**2;
                            if (distSqToPlaced < minSqDistBetweenPlacedArrows_ideal) {
                                tooCloseToExistingArrow = true;
                                break;
                            }
                        }

                        if (tooCloseToExistingArrow) continue;

                        const centroidPointData = curve[headIndex]; // currentPoint is the centroid for this arrow
                        let point1_for_direction_ternary: number[], point2_for_direction_ternary: number[];
                        
                        if (headIndex + 1 < curve.length) {
                            point1_for_direction_ternary = curve[headIndex].x;
                            point2_for_direction_ternary = curve[headIndex + 1].x;
                        } else if (headIndex > 0) {
                            point1_for_direction_ternary = curve[headIndex - 1].x;
                                                       point2_for_direction_ternary = curve[headIndex].x;
                        } else {
                            continue; 
                        }
                        
                        const arePointsEffectivelyIdentical =
                            Math.abs(point1_for_direction_ternary[0] - point2_for_direction_ternary[0]) < 1e-7 &&
                            Math.abs(point1_for_direction_ternary[1] - point2_for_direction_ternary[1]) < 1e-7 &&

                            Math.abs(point1_for_direction_ternary[2] - point2_for_direction_ternary[2]) < 1e-7;

                        if (arePointsEffectivelyIdentical) continue;
                        
                        const localPrevDirCoords = convertTernaryToPaperCoordinates(point1_for_direction_ternary, plotSortedComponents);
                        const localNextDirCoords = convertTernaryToPaperCoordinates(point2_for_direction_ternary, plotSortedComponents);

                        if (!localPrevDirCoords || !localNextDirCoords) continue;

                        const dx_ideal_CALCULATED = localNextDirCoords.x - localPrevDirCoords.x;
                        const dy_ideal_CALCULATED = localNextDirCoords.y - localPrevDirCoords.y;

                        if (Math.abs(dx_ideal_CALCULATED) < 1e-9 && Math.abs(dy_ideal_CALCULATED) < 1e-9) {
                            continue;
                        }
                        
                        const triangleActualPixelBase = actualRenderedTriangleBase_paper * containerWidthPx;
                        const pixel_dx = dx_ideal_CALCULATED * triangleActualPixelBase;
                        const pixel_dy = dy_ideal_CALCULATED * triangleActualPixelBase;
                        
                        if (Math.abs(pixel_dx) < 1e-6 && Math.abs(pixel_dy) < 1e-6) {
                            continue; 
                        }
                            
                        // Calculate arrow rotation so that the triangle symbol points along the flow direction.
                        // The base 'triangle-up' in Plotly points along +Y (upward). To align it with the
                        // vector (dx,dy) in screen space, rotate the symbol by (raw_angle − 90°).
                        const raw_screen_angle_deg = Math.atan2(-pixel_dy, pixel_dx) * (180 / Math.PI);
                        const final_angle_deg = raw_screen_angle_deg + 90;

                        arrowMarkerTraceData.a.push(centroidPointData.x[1]);
                        arrowMarkerTraceData.b.push(centroidPointData.x[2]);
                        arrowMarkerTraceData.c.push(centroidPointData.x[0]);
                        arrowMarkerTraceData.angles.push(final_angle_deg);
                        // No tooltip text for arrow markers
                        arrowMarkerTraceData.text.push('');
                        
                        allPlacedArrowIdealCoords.push(currentPointIdealCoords); 
                        placedArrowCountOverall++;

                        if (!firstArrowLogged && curveIndex === 0) { 
                            // Arrow marker debug disabled.
                        }
                    } // Closes if (shouldCheckThisPointForArrow)
                } // Closes for (let headIndex ...) loop
            }); // Closes cleanedResidueCurves.forEach
        } // Closes if (plotSortedComponents.length === 3 && plotContainerRef && ...)
        
        const allTraces = [...traces];
        if (arrowMarkerTraceData.a.length > 0) {
            allTraces.push({
                type: 'scatterternary',
                mode: 'markers',
                a: arrowMarkerTraceData.a,
                b: arrowMarkerTraceData.b,
                c: arrowMarkerTraceData.c,
                text: arrowMarkerTraceData.text,
                hoverinfo: 'skip',
                hoverlabel: {
                    bgcolor: currentTheme === 'dark' ? 'rgba(50,50,50,0.85)' : 'rgba(250,250,250,0.85)', // Adjusted for theme
                    font: {
                        family: merriweatherFamilyString, // Explicitly set family
                        size: 14,
                        color: plotFontColor // Use defined color
                    },
                    bordercolor: currentTheme === 'dark' ? '#777' : '#ccc'
                },
                marker: {
                    symbol: 'triangle-up',
                    size: targetArrowPixelSize,
                    color: arrowMarkerColor, // Theme-aware
                    angle: arrowMarkerTraceData.angles,
                    standoff: 5,
                } as any,
                cliponaxis: false,
                showlegend: false,
            } as Data);
        }

        const minAz: AzeotropeDisplayInfo[] = [];
        const maxAz: AzeotropeDisplayInfo[] = [];
        const sadAz: AzeotropeDisplayInfo[] = [];

        // --- Directly Found Azeotropes to Plot (Two-Trace Method for Animation) ---
        if (directAzeotropes.length > 0 && plotSortedComponents.length === 3) {
            // --- Robust Classification Logic ---
            const classifyEigen = (az:AzeotropeDisplayInfo):'min'|'max'|'saddle'|'unknown' => {
               const h = 1e-4;
               if(backendComps.length!==3 || !backendPkgParams) return 'unknown';
               const base = evaluateResidueODE(displayedFluidPackage, az.x, az.T_K, backendPressurePa, backendComps as unknown as CompoundData[], backendPkgParams);
               if(!base) return 'unknown';
               const J:number[][]=[ [0,0,0],[0,0,0],[0,0,0] ];
               for(let j=0;j<3;j++){
                  const xPert=[...az.x]; xPert[j]+=h; 
                  const r = evaluateResidueODE(displayedFluidPackage, xPert, az.T_K, backendPressurePa, backendComps as unknown as CompoundData[], backendPkgParams);
                  if(!r) return 'unknown';
                  for(let i=0;i<3;i++) J[i][j]=(r.d[i]-base.d[i])/h;
               }
               const A=[
                 [J[0][0]-J[2][0], J[0][1]-J[2][1]],
                 [J[1][0]-J[2][0], J[1][1]-J[2][1]]
               ];
               const tr=A[0][0]+A[1][1];
               const det=A[0][0]*A[1][1]-A[0][1]*A[1][0];
               const disc=tr*tr-4*det;
               if(disc<0) return 'saddle';
               const sqrtDisc=Math.sqrt(disc);
               const l1=0.5*(tr+sqrtDisc);
               const l2=0.5*(tr-sqrtDisc);
               if(l1<0 && l2<0) return 'max';
               if(l1>0 && l2>0) return 'min';
               return 'saddle';
            };

            for(const az of directAzeotropes){
               const cls=classifyEigen(az);
               switch(cls){
                   case 'min': minAz.push(az); break;
                   case 'max': maxAz.push(az); break;
                   case 'saddle': default: sadAz.push(az); break;
               }
            }

            const fmt3 = (v:number)=>{
                if(Math.abs(v) < 0.0005) return '0'; // treat very small values as zero
                const s = Number(v).toPrecision(3);
                // Convert to number then back to string to avoid exponential notation if possible
                const num = parseFloat(s);
                // For values that are still in exponential form after parseFloat (e.g., 1e-5), force fixed
                if(num.toString().includes('e')){
                    return num.toFixed(3).replace(/\.0+$/,'');
                }
                return num.toString();
            };
            const mkText = (az:AzeotropeDisplayInfo)=>
                `${displayedComponentNames[1] || 'Comp 2'}: ${fmt3(az.x[1])}<br>`+
                `${displayedComponentNames[2] || 'Comp 3'}: ${fmt3(az.x[2])}<br>`+
                `${displayedComponentNames[0] || 'Comp 1'}: ${fmt3(az.x[0])}<br>`+
                `T: ${fmt3(az.T_K-273.15)}°C`;

            const pushAzeoTrace = (dataArr:typeof minAz, color:string, name:string) => {
                if (dataArr.length === 0) return;
                allTraces.push({
                    type: 'scatterternary',
                    mode: 'markers',
                    a: dataArr.map(az=>az.x[1]),
                    b: dataArr.map(az=>az.x[2]),
                    c: dataArr.map(az=>az.x[0]),
                    cliponaxis:false,
                    text: dataArr.map(az=>mkText(az)),
                    hoverinfo:'text',
                    marker:{symbol:'star',size:14,color:color,opacity:1,line:{color: color, width:1.5}},
                    name:name,
                    legendgroup:'azeotropes',
                    hoverlabel: {
                        font: {
                            family: merriweatherFamilyString,
                            size: 14,
                            color: plotFontColor
                        }
                    }
                } as Data);
            };

            pushAzeoTrace(minAz, '#00C000', 'Min-boiling'); // brighter green
            pushAzeoTrace(maxAz, '#9900FF', 'Max-boiling'); // saturated purple
            pushAzeoTrace(sadAz, 'red', 'Saddle');
        }
        // --- End Add Directly Found Azeotropes to Plot ---


        const finalLayout: Partial<Layout> = {
            ...baseLayout,
            shapes: [], // Clear old shapes if any, or manage them if other shapes are needed
            annotations: [...(baseLayout.annotations || [])], 
            title: { 
                ...(typeof baseLayout.title === 'object' ? baseLayout.title : { text: '', font: { family: merriweatherFamilyString, size: 18, color: plotTitleColor} } ), 
                text: `Ternary Residue Curve Map at ${displayedPressure} bar` 
            },
            legend: {
                x: 0.8,
                xanchor: 'left' as const,
                y: 1,
                font: basePlotFontObject,
                bgcolor: currentTheme === 'dark' ? 'rgba(8, 48, 107, 0.8)' : 'rgba(220, 230, 240, 0.8)',
                bordercolor: currentTheme === 'dark' ? '#4b5563' : '#9ca3af',
                borderwidth: 1,
                tracegroupgap: 5,
            },
            // Removed transition to eliminate marker animation
        };

        return { baseTraces: allTraces, baseLayout: finalLayout, minAz, maxAz, sadAz };

    }, [cleanedResidueCurves, componentsInput, displayedPressure, plotSortedComponents, isLoading, currentTheme, displayedFluidPackage, directAzeotropes, plotAxisTitles, plotContainerRef, backendComps, backendPkgParams, backendPressurePa]);

    useEffect(() => {
        const { baseTraces, baseLayout, minAz, maxAz, sadAz } = memoizedPlotData;

        const finalTraces = [...baseTraces];

        // Highlight trace reinstated (only when row hovered)
        if (highlightedAzeoIdx !== null && highlightedAzeoIdx < directAzeotropes.length) {
            const highlightedAzeo = directAzeotropes[highlightedAzeoIdx];

            // --- Determine color based on earlier classification ---
            let highlightColor = 'red'; // default to saddle color
            if (minAz.includes(highlightedAzeo)) {
                highlightColor = '#00C000';
            } else if (maxAz.includes(highlightedAzeo)) {
                highlightColor = '#9900FF';
            }

            finalTraces.push({
                type: 'scatterternary',
                mode: 'markers',
                a: [highlightedAzeo.x[1]],
                b: [highlightedAzeo.x[2]],
                c: [highlightedAzeo.x[0]],
                cliponaxis:false,
                hoverinfo:'skip',
                marker:{symbol:'star',size:28,color:highlightColor,opacity:1,line:{color: currentTheme==='dark'? '#fff':'#000', width:2}},
                showlegend:false,
                legendgroup:'azeotropes'
            } as Data);
        }
        
        setPlotlyData(finalTraces);
        setPlotlyLayout(baseLayout);

    }, [memoizedPlotData, highlightedAzeoIdx, directAzeotropes, currentTheme]);

    const handleGenerateClick = () => {
        handleGenerateMap();
    };

    return (
        <div className="container mx-auto p-4">
            <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
                {/* Left panel: Inputs */}
                <div className="lg:col-span-1 space-y-6">
                    <Card>
                        <CardHeader>
                            <CardTitle>Residue Map Inputs</CardTitle>
                        </CardHeader>
                        <CardContent className="space-y-4">
                            {/* Component names */}
                            <div className="space-y-3"> {/* Grouping div for component inputs */}
                                {componentsInput.map((component, index) => (
                                    <div className="flex items-center space-x-2" key={index}> {/* Flex container for label and input */}
                                        <Label htmlFor={`name-${index}`} className="whitespace-nowrap">{`Component ${index + 1}:`}</Label>
                                        <div className="relative w-full"> {/* Relative container for input and suggestions */}
                                            <Input
                                                id={`name-${index}`}
                                                ref={el => { inputRefs.current[index] = el; }}
                                                value={component.name}
                                                onChange={e => handleComponentInputChange(index, 'name', e.target.value)}
                                                placeholder={`e.g. Acetone`}
                                                autoComplete="off"
                                                onKeyDown={e => {
                                                    if (e.key === 'Enter') {
                                                        e.preventDefault();
                                                        handleGenerateClick();
                                                    }
                                                }}
                                                onFocus={() => {
                                                    setActiveSuggestionIndex(index);
                                                    if (componentsInput[index].name.trim()) {
                                                        fetchComponentSuggestions(componentsInput[index].name, index);
                                                    }
                                                }}
                                            />
                                            {showSuggestions[index] && componentSuggestions[index].length > 0 && (
                                                <ul
                                                    ref={el => { suggestionsContainerRefs.current[index] = el; }}
                                                    className={`absolute z-10 mt-1 w-full max-h-40 overflow-y-auto rounded-md shadow-md
                                                      ${currentTheme === 'dark'
                                                        ? 'bg-gray-800 border border-gray-700 text-gray-200'
                                                        : 'bg-white border border-gray-300 text-gray-900'}`}
                                                >
                                                    {componentSuggestions[index].map((suggestion, si) => (
                                                        <li
                                                            key={si}
                                                            className={`px-4 py-2 cursor-pointer
                                                              ${currentTheme === 'dark' ? 'hover:bg-gray-700' : 'hover:bg-gray-100'}`}
                                                            onMouseDown={() => handleSuggestionClick(index, suggestion)}
                                                        >
                                                            {suggestion}
                                                        </li>
                                                    ))}
                                                </ul>
                                            )}
                                        </div> {/* End of relative container for input and suggestions */}
                                    </div>
                                ))}

                            </div> {/* End of grouping div for component inputs */}
                            {/* Pressure on same line */}
                            <div className="flex items-center space-x-2">
                                <Label htmlFor="systemPressure" className="whitespace-nowrap">
                                    Pressure:
                                </Label>
                                <div className="flex items-center w-full"> {/* Wrapper for input and unit */}
                                    <Input
                                        id="systemPressure"
                                        type="number"
                                        value={systemPressure}
                                        onChange={e => setSystemPressure(e.target.value)}
                                        onKeyDown={e => {
                                            if (e.key === 'Enter') {
                                                e.preventDefault();
                                                handleGenerateClick();
                                            }
                                        }}
                                        placeholder="1.0"
                                        className="flex-grow"
                                    />
                                    <span className="ml-2 text-muted-foreground">bar</span> {/* Unit display */}
                                </div>
                            </div>
                            {/* Fluid package on same line */}
                            <div className="flex items-center space-x-2">
                                <Label htmlFor="fluidPackage" className="whitespace-nowrap">
                                    Fluid Package:
                                </Label>
                                <Select
                                    value={fluidPackage}
                                    onValueChange={v => {
                                        setFluidPackage(v as FluidPackageTypeResidue);
                                        // Auto-trigger map generation after fluid package change
                                        setTimeout(() => handleGenerateMap(), 100);
                                    }}
                                >
                                    <SelectTrigger id="fluidPackage">
                                        <SelectValue placeholder="Select model" />
                                    </SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="wilson">Wilson</SelectItem>
                                        <SelectItem value="unifac">UNIFAC</SelectItem>
                                        <SelectItem value="nrtl">NRTL</SelectItem>
                                        <SelectItem value="pr">Peng–Robinson</SelectItem>
                                        <SelectItem value="srk">SRK</SelectItem>
                                        <SelectItem value="uniquac">UNIQUAC</SelectItem>
                                    </SelectContent>
                                </Select>
                            </div>
                            <Button onClick={handleGenerateClick} disabled={isLoading} className="w-full">
                                {isLoading ? "Generating..." : "Generate Map"}
                            </Button>
                        </CardContent>
                    </Card>

                    {directAzeotropes.length > 0 && (
                        <Card>
                            <CardHeader>
                                <CardTitle>Azeotropic Composition and Temperature</CardTitle>
                            </CardHeader>
                            <CardContent>
                                <Table>
                                    <TableHeader>
                                        <TableRow className="hover:bg-transparent">
                                            <TableHead className="px-2 text-center">{displayedComponentNames[0] || 'Comp 1'}</TableHead>
                                            <TableHead className="px-2 text-center">{displayedComponentNames[1] || 'Comp 2'}</TableHead>
                                            <TableHead className="px-2 text-center">{displayedComponentNames[2] || 'Comp 3'}</TableHead>
                                            <TableHead className="px-2 text-center">T(°C)</TableHead>
                                        </TableRow>
                                    </TableHeader>
                                    <TableBody>
                                        {directAzeotropes.map((az, index) => {
                                            return (
                                                <TableRow
                                                    key={index}
                                                    onMouseEnter={() => setHighlightedAzeoIdx(index)}
                                                    onMouseLeave={() => setHighlightedAzeoIdx(null)}
                                                >
                                                 <TableCell className="px-2 text-center">{az.x[0].toFixed(3)}</TableCell>
                                                 <TableCell className="px-2 text-center">{az.x[1].toFixed(3)}</TableCell>
                                                 <TableCell className="px-2 text-center">{az.x[2].toFixed(3)}</TableCell>
                                                 <TableCell className="px-2 text-center">{(az.T_K - 273.15).toFixed(1)}</TableCell>
                                                 </TableRow>
                                            );
                                        })}
                                    </TableBody>
                                </Table>
                            </CardContent>
                        </Card>
                    )}
                </div>

                {/* Right panel: Plot */}
                <div className="lg:col-span-2">
                    <Card>
                        <CardContent className="p-2">
                            <div 
                                className="relative w-full h-[600px] md:h-[700px]" 
                                ref={setPlotContainerRef}
                            >
                                {scfMessage && !isLoading ? (
                                    <div className="flex items-center justify-center w-full h-full text-center px-4">
                                        <p className="text-lg md:text-2xl font-semibold whitespace-pre-line">{scfMessage}</p>
                                    </div>
                                ) : error && !isLoading ? (
                                    <div className="absolute inset-0 flex items-center justify-center bg-background/80 z-10">
                                        <Card className="border-destructive">
                                            <CardHeader>
                                                <CardTitle>Error</CardTitle>
                                            </CardHeader>
                                            <CardContent>
                                                <p>{error}</p>
                                            </CardContent>
                                        </Card>
                                    </div>
                                ) : (
                                    <Plot
                                        data={plotlyData}
                                        layout={plotlyLayout}
                                        style={{ width: '100%', height: '100%' }}
                                        config={{ displayModeBar: false }}
                                        useResizeHandler
                                    />
                                )}
                            </div>
                        </CardContent>
                    </Card>
                </div>
            </div>
        </div>
    );
}

/**
 * Shuffles an array in-place using the Fisher-Yates algorithm.
 * @param array The array to shuffle.
 */
function shuffleArray<T>(array: T[]): void {
    for (let i = array.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [array[i], array[j]] = [array[j], array[i]];
    }
}

/**
 * Filters a set of residue curves to remove ones that are too close to each other
 * in the central region of the ternary diagram, preventing visual clutter.
 * @param curves The raw array of residue curves.
 * @param distanceThreshold The minimum mole fraction distance between curves in the central area.
 * @returns A filtered array of residue curves.
 */
function filterOverlappingCurves(curves: ResidueCurve[], distanceThreshold: number): ResidueCurve[] {
    if (curves.length === 0) return [];

    const centralRegionMin = 0.05; // Not on an edge
    const centralRegionMax = 0.95; // Not in a corner

    const isCentral = (x: number[]) => {
        return x.every(xi => xi > centralRegionMin && xi < centralRegionMax);
    };

    // Pre-filter points in the central region for each curve to optimize.
    const centralPointsByCurve = curves.map(curve => curve.filter(p => isCentral(p.x)));

    const acceptedCurvesIndices: number[] = [];
    if (curves.length > 0) {
        acceptedCurvesIndices.push(0); // Always accept the first curve.
    }

    for (let i = 1; i < curves.length; i++) {
        const candidateCentralPoints = centralPointsByCurve[i];
        if (candidateCentralPoints.length === 0) {
            // If curve has no central points, it doesn't contribute to central clutter, so accept.
            acceptedCurvesIndices.push(i);
            continue;
        }

        let isTooClose = false;
        for (const acceptedIndex of acceptedCurvesIndices) {
            const acceptedCentralPoints = centralPointsByCurve[acceptedIndex];
            if (acceptedCentralPoints.length === 0) continue;

            // Check for proximity
            for (const p_candidate of candidateCentralPoints) {
                for (const p_accepted of acceptedCentralPoints) {
                    const distSq =
                        (p_candidate.x[0] - p_accepted.x[0]) ** 2 +
                        (p_candidate.x[1] - p_accepted.x[1]) ** 2 +
                        (p_candidate.x[2] - p_accepted.x[2]) ** 2;

                    if (distSq < distanceThreshold ** 2) {
                        isTooClose = true;
                        break;
                    }
                }
                if (isTooClose) break;
            }
            if (isTooClose) break;
        }

        if (!isTooClose) {
            acceptedCurvesIndices.push(i);
        }
    }
    return acceptedCurvesIndices.map(i => curves[i]);
}

// Classify an azeotrope using the Jacobian eigenvalues of the residue-curve ODE
function classifyAzeotropeEigen(
    az: AzeotropeDisplayInfo,
    pkg: FluidPackageTypeResidue,
    P_Pa: number,
    comps: CompoundData[],
    pkgParams: any,
): 'min' | 'max' | 'saddle' | 'unknown' {
    const h = 1e-4;
    const base = evaluateResidueODE(pkg, az.x, az.T_K, P_Pa, comps, pkgParams);
    if (!base) return 'unknown';
    const J: number[][] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (let j = 0; j < 3; j++) {
        const xPert = [...az.x];
        xPert[j] += h;
        const r = evaluateResidueODE(pkg, xPert, az.T_K, P_Pa, comps, pkgParams);
        if (!r) return 'unknown';
        for (let i = 0; i < 3; i++) {
            J[i][j] = (r.d[i] - base.d[i]) / h;
        }
    }
    // Project 3×3 Jacobian onto 2D tangent plane to assess stability
    const A = [
        [J[0][0] - J[2][0], J[0][1] - J[2][1]],
        [J[1][0] - J[2][0], J[1][1] - J[2][1]],
    ];
    const tr = A[0][0] + A[1][1];
    const det = A[0][0] * A[1][1] - A[0][1] * A[1][0];
    const disc = tr * tr - 4 * det;
    if (disc < 0) return 'saddle';
    const sqrtDisc = Math.sqrt(disc);
    const l1 = 0.5 * (tr + sqrtDisc);
    const l2 = 0.5 * (tr - sqrtDisc);
    if (l1 < 0 && l2 < 0) return 'max';
    if (l1 > 0 && l2 > 0) return 'min';
    return 'saddle';
}

