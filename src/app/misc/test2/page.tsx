'use client';

import React, { useState, useEffect, useRef, useCallback } from 'react';
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
    findAzeotropeNRTL,
    findAzeotropePr,
    findAzeotropeSrk,
    findAzeotropeUniquac,
    type TernaryWilsonParams,
    type TernaryNrtlParams,
    type TernaryPrParams,
    type TernarySrkParams,
    type TernaryUniquacParams,
    type AzeotropeResult,
} from '@/lib/residue-curves-ode';
import { fetchAndConvertThermData, FetchedCompoundThermData } from '@/lib/antoine-utils'; 
import type { Data, Layout } from 'plotly.js'; 
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select"; 
import { fetchWilsonInteractionParams, R_gas_const_J_molK } from '@/lib/vle-calculations-wilson';
import { fetchUnifacInteractionParams, UnifacParameters, calculatePsat_Pa } from '@/lib/vle-calculations-unifac';
import { fetchNrtlParameters, type NrtlInteractionParams } from '@/lib/vle-calculations-nrtl';
import { fetchPrInteractionParams, type PrInteractionParams } from '@/lib/vle-calculations-pr';
import { fetchSrkInteractionParams, type SrkInteractionParams } from '@/lib/vle-calculations-srk';
import { fetchUniquacInteractionParams, type UniquacInteractionParams } from '@/lib/vle-calculations-uniquac';

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

// Debounce function
function debounce<F extends (...args: any[]) => any>(func: F, waitFor: number) {
  let timeout: ReturnType<typeof setTimeout> | null = null;
  const debounced = (...args: Parameters<F>) => {
    if (timeout !== null) { clearTimeout(timeout); timeout = null; }
    timeout = setTimeout(() => func(...args), waitFor);
  };
  return debounced as (...args: Parameters<F>) => void;
}

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

// Alias specific azeotrope result type names used later in the file to the unified AzeotropeResult
type NrtlAzeotropeResult = AzeotropeResult;
type PrAzeotropeResult = AzeotropeResult;
type SrkAzeotropeResult = AzeotropeResult & { y?: number[] };
type UniquacAzeotropeResult = AzeotropeResult;

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

// Placeholder type for UniquacAzeotropeResult if not properly imported.
// Ensure this matches the actual structure from residue-curves-ode-uniquac.ts
// type UniquacAzeotropeResult = NrtlAzeotropeResult; // Example if similar to NRTL

export default function TernaryResidueMapPage() {
    const [componentsInput, setComponentsInput] = useState<ComponentInputState[]>([
        { name: 'Acetone' },
        { name: 'Methanol' },
        { name: 'Water' },
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
    const initialMapGenerated = useRef(false); // Ref to track initial generation
    const [currentTheme, setCurrentTheme] = useState<'dark' | 'light'>('dark'); // Default to dark
    const [plotContainerRef, setPlotContainerRef] = React.useState<HTMLDivElement | null>(null); // State to hold the div element
    const [directAzeotropes, setDirectAzeotropes] = useState<AzeotropeDisplayInfo[]>([]); // Use common display info

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
    
    const debouncedFetchComponentSuggestions = useCallback(debounce(fetchComponentSuggestions, 300), [supabase]);

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
                debouncedFetchComponentSuggestions(value, index);
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
        setDirectAzeotropes([]); // Clear previous direct azeotropes

        try {
            const P_system_Pa = parseFloat(systemPressure) * 100000; // Convert bar to Pa
            if (isNaN(P_system_Pa) || P_system_Pa <= 0) {
                throw new Error("Invalid system pressure input.");
            }

            const casNumbersPromises = componentsInput.map(input => 
                fetchCasNumberByName(supabase!, input.name)
            );
            const fetchedCasNumbers = await Promise.all(casNumbersPromises);

            const thermoDataPromises = fetchedCasNumbers.map(casn => 
                fetchAndConvertThermData(supabase!, casn)
            );
            const fetchedThermDataArray = await Promise.all(thermoDataPromises);

            const processedComponents: ProcessedComponentData[] = componentsInput.map((input, index) => {
                const thermData = fetchedThermDataArray[index];
                let bp_at_Psys_K_val: number | null = null;
                if (thermData.antoine) {
                    bp_at_Psys_K_val = calculateBoilingPointAtPressureSecant(
                        thermData.antoine,
                        P_system_Pa,
                        thermData.nbp_K || DEFAULT_INITIAL_T_K, // Use NBP as initial guess, or default
                        calculatePsat_Pa // Pass the imported Psat function
                    );
                }
                return {
                    name: input.name,
                    casNumber: fetchedCasNumbers[index],
                    thermData: thermData,
                    originalIndex: index,
                    bp_at_Psys_K: bp_at_Psys_K_val,
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

            const sortedForPlot = [...processedComponents].sort((a, b) => {
                // Sorting for plot axes (L, M, H) should still be based on NBP @ 1atm for consistency
                if (a.thermData.nbp_K === null && b.thermData.nbp_K === null) return 0;
                if (a.thermData.nbp_K === null) return 1;
                if (b.thermData.nbp_K === null) return -1;
                return a.thermData.nbp_K - b.thermData.nbp_K;
            });
            setPlotSortedComponents(sortedForPlot);

            const compoundsForBackend: CompoundData[] = processedComponents
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
                const params01 = await fetchWilsonInteractionParams(supabase, cas[0], cas[1]);
                const params02 = await fetchWilsonInteractionParams(supabase, cas[0], cas[2]);
                const params12 = await fetchWilsonInteractionParams(supabase, cas[1], cas[2]);
                activityModelParams = {
                    a01_J_mol: params01.a12_J_mol, a10_J_mol: params01.a21_J_mol,
                    a02_J_mol: params02.a12_J_mol, a20_J_mol: params02.a21_J_mol,
                    a12_J_mol: params12.a12_J_mol, a21_J_mol: params12.a21_J_mol,
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
            
            const avgNBP_at_Psys = processedComponents
                .map(p => p.bp_at_Psys_K) // Use BP at system pressure
                .filter(bp => bp !== null && bp !== undefined) as number[];
            
            const initialTempGuess_K = avgNBP_at_Psys.length > 0 
                ? avgNBP_at_Psys.reduce((sum, val) => sum + val, 0) / avgNBP_at_Psys.length 
                : (processedComponents
                    .map(p => p.thermData.nbp_K)
                    .filter(nbp => nbp !== null) as number[]
                  ).length > 0 
                    ? (processedComponents
                        .map(p => p.thermData.nbp_K)
                        .filter(nbp => nbp !== null) as number[]
                      ).reduce((sum, val) => sum + val, 0) / (processedComponents
                        .map(p => p.thermData.nbp_K)
                        .filter(nbp => nbp !== null) as number[]
                      ).length
                    : DEFAULT_INITIAL_T_K;

            // Determine parameters for starting point generation
            // Use different (potentially fewer) points for UNIFAC if it's much slower
            const numPerEdgeUnifac = 6; // Reduced for UNIFAC
            const numInternalDivisionsUnifac = 7; // Reduced for UNIFAC

            const numPerEdge = (fluidPackage === 'unifac') ? numPerEdgeUnifac : 
                               (fluidPackage === 'pr' || fluidPackage === 'srk') ? 10 : 8; // Adjusted non-UNIFAC from 6 to 8
            const numInternalDivisions = (fluidPackage === 'unifac') ? numInternalDivisionsUnifac :
                                         (fluidPackage === 'pr' || fluidPackage === 'srk') ? 10 : 8; // Adjusted non-UNIFAC from 7 to 8
            
            // Call the modified generateStartingPoints function
            let startingCompositions = generateStartingPoints(numPerEdge, numInternalDivisions, fluidPackage);
            
            const allCurvesPromises: Promise<ResidueCurve | null>[] = [];

            // Adjust d_xi_step for each model
            const base_d_xi_step_activity = 0.03; // Wilson, NRTL, UNIQUAC, UNIFAC - larger step
            const base_d_xi_step_eos = 0.035;    // Larger step for PR/SRK

            const d_xi_step_for_model =
                (fluidPackage === 'pr' || fluidPackage === 'srk') ? base_d_xi_step_eos
                : base_d_xi_step_activity; // UNIFAC now uses this

            const max_steps_unifac = 600; // Reduced max steps for UNIFAC
            const max_steps = (fluidPackage === 'unifac') ? max_steps_unifac :
                              (fluidPackage === 'pr' || fluidPackage === 'srk') ? 8000 : 800;

            for (const start_x_3comp of startingCompositions) {
                let promise: Promise<ResidueCurve | null>;
                switch (fluidPackage) {
                    case 'wilson':
                        promise = simulateResidueCurveODE('wilson', start_x_3comp, P_system_Pa, compoundsForBackend, activityModelParams as TernaryWilsonParams, d_xi_step_for_model, max_steps, initialTempGuess_K);
                        break;
                    case 'unifac':
                        promise = simulateResidueCurveODE('unifac', start_x_3comp, P_system_Pa, compoundsForBackend, activityModelParams as UnifacParameters, d_xi_step_for_model, max_steps, initialTempGuess_K);
                        break;
                    case 'nrtl':
                        promise = simulateResidueCurveODE('nrtl', start_x_3comp, P_system_Pa, compoundsForBackend, activityModelParams as TernaryNrtlParams, d_xi_step_for_model, max_steps, initialTempGuess_K);
                        break;
                    case 'pr':
                        promise = simulateResidueCurveODE(
                            'pr',
                            start_x_3comp,
                            P_system_Pa,
                            compoundsForBackend,
                            activityModelParams as TernaryPrParams,
                            d_xi_step_for_model,
                            max_steps,
                            initialTempGuess_K
                        );
                        break;
                    case 'srk':
                        promise = simulateResidueCurveODE(
                            'srk',
                            start_x_3comp,
                            P_system_Pa,
                            compoundsForBackend,
                            activityModelParams as TernarySrkParams,
                            d_xi_step_for_model,
                            max_steps,
                            initialTempGuess_K
                        );
                        break;
                    case 'uniquac':
                        promise = simulateResidueCurveODE('uniquac', start_x_3comp, P_system_Pa, compoundsForBackend, activityModelParams as TernaryUniquacParams, d_xi_step_for_model, max_steps, initialTempGuess_K);
                        break;
                    default:
                        throw new Error(`Simulation function not implemented for ${fluidPackage}`);
                }
                allCurvesPromises.push(promise);
            }
            
            const resolvedCurves = await Promise.all(allCurvesPromises);

            const validCurves = resolvedCurves.filter(curve => curve !== null && curve.length > 1) as ResidueCurve[];

            console.log(`Residue Curve Generation (${fluidPackage}): Total starting points: ${startingCompositions.length}. Resolved curves (pre-filter): ${resolvedCurves.length}. Null/Short curves: ${resolvedCurves.filter(c => c === null || c.length <=1).length}. Valid curves (length > 1): ${validCurves.length}`);
            if (validCurves.length > 0) {
                console.log(`First valid curve (first 5 points): `, validCurves[0].slice(0, 5).map(p => ({x: p.x.map(val => val.toFixed(4)), T: p.T_K.toFixed(2)})));
                 if (validCurves[0].length > 5) {
                    console.log(`First valid curve (last 5 points): `, validCurves[0].slice(-5).map(p => ({x: p.x.map(val => val.toFixed(4)), T: p.T_K.toFixed(2)})));
                }
            }

            setResidueCurves(validCurves);

            // --- Call Direct Azeotrope Solver ---
            let foundAzeotropes: AzeotropeDisplayInfo[] = [];

            if (compoundsForBackend.length === 3) {
                const initialGuessesRaw: { x: number[], T_K: number }[] = [ // Added type for initialGuessesRaw elements
                    { x: [0.333, 0.333, 0.334], T_K: initialTempGuess_K }, // Center
                ];
                if (plotSortedComponents.length === 3) {
                    initialGuessesRaw.push({ x: [0.9, 0.05, 0.05], T_K: plotSortedComponents[0].bp_at_Psys_K || plotSortedComponents[0].thermData.nbp_K || initialTempGuess_K });
                    initialGuessesRaw.push({ x: [0.05, 0.9, 0.05], T_K: plotSortedComponents[1].bp_at_Psys_K || plotSortedComponents[1].thermData.nbp_K || initialTempGuess_K });
                    initialGuessesRaw.push({ x: [0.05, 0.05, 0.9], T_K: plotSortedComponents[2].bp_at_Psys_K || plotSortedComponents[2].thermData.nbp_K || initialTempGuess_K });
                }
                if (compoundsForBackend.length === 3) {
                    const T0 = plotSortedComponents.find(c => c.originalIndex === 0)?.bp_at_Psys_K || plotSortedComponents.find(c => c.originalIndex === 0)?.thermData.nbp_K || initialTempGuess_K;
                    const T1 = plotSortedComponents.find(c => c.originalIndex === 1)?.bp_at_Psys_K || plotSortedComponents.find(c => c.originalIndex === 1)?.thermData.nbp_K || initialTempGuess_K;
                    const T2 = plotSortedComponents.find(c => c.originalIndex === 2)?.bp_at_Psys_K || plotSortedComponents.find(c => c.originalIndex === 2)?.thermData.nbp_K || initialTempGuess_K;
                    initialGuessesRaw.push({ x: [0.5, 0.5, 0.001], T_K: (T0+T1)/2 });
                    initialGuessesRaw.push({ x: [0.5, 0.001, 0.5], T_K: (T0+T2)/2 });
                    initialGuessesRaw.push({ x: [0.001, 0.5, 0.5], T_K: (T1+T2)/2 });
                }

                let results: (NrtlAzeotropeResult | PrAzeotropeResult | SrkAzeotropeResult | UniquacAzeotropeResult | null)[] = [];

                if (generatingFluidPackage === 'nrtl') {
                    const foundAzeotropesPromises = initialGuessesRaw.map(guess => // 'guess' type is inferred from initialGuessesRaw
                        findAzeotropeNRTL(
                            P_system_Pa,
                            compoundsForBackend,
                            activityModelParams as TernaryNrtlParams,
                            guess.x,
                            guess.T_K,
                            150, 1e-7
                        )
                    );
                    results = await Promise.all(foundAzeotropesPromises);
                } else if (generatingFluidPackage === 'pr') {
                    results = initialGuessesRaw.map(guess => // 'guess' type is inferred
                        findAzeotropePr(
                            P_system_Pa,
                            compoundsForBackend,
                            activityModelParams as TernaryPrParams,
                            guess.x,
                            guess.T_K,
                            150, 1e-7
                        )
                    );
                } else if (generatingFluidPackage === 'srk') {
                    results = initialGuessesRaw.map(guess => // 'guess' type is inferred
                        findAzeotropeSrk(
                            guess.x, 
                            P_system_Pa,
                            compoundsForBackend,
                            activityModelParams as TernarySrkParams,
                            [1, 1, 1] 
                        )
                    );
                } else if (generatingFluidPackage === 'uniquac') {
                    // Assuming findAzeotropeUniquac and UniquacAzeotropeResult are properly defined and imported
                    // from '@/lib/residue-curves-ode-uniquac'
                    if (typeof findAzeotropeUniquac === 'function') {
                        results = initialGuessesRaw.map(guess =>
                            findAzeotropeUniquac(
                                P_system_Pa,
                                compoundsForBackend,
                                activityModelParams as TernaryUniquacParams,
                                guess.x, 
                                guess.T_K,
                                150, 1e-7 // Assuming similar parameters to NRTL/PR
                            )
                        );
                    } else {
                        console.warn("UNIQUAC direct azeotrope solver (findAzeotropeUniquac) is not available. Skipping.");
                        results = [];
                    }
                }

                const convergedResults = results.filter(r => {
                    if (!r) return false;
                    if (generatingFluidPackage === 'srk') {
                        const srkRes = r as SrkAzeotropeResult;
                        // Ensure y exists and is an array of 3 numbers for SRK
                        if (!srkRes.y || srkRes.y.length !== 3 || !srkRes.x || srkRes.x.length !== 3) {
                            console.warn("SRK Azeotrope result is missing x or y, or they are not 3-component arrays. Skipping.", srkRes);
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
                    // For PR, NRTL, and now UNIQUAC, they should have a 'converged' property
                    return (r as NrtlAzeotropeResult | PrAzeotropeResult | UniquacAzeotropeResult).converged;
                }) as Array<NrtlAzeotropeResult | PrAzeotropeResult | UniquacAzeotropeResult | (SrkAzeotropeResult & { calculatedErrorNormYX?: number })>;


                const uniqueDirectAzeotropes: AzeotropeDisplayInfo[] = [];
                const minSqDistAzeo = 0.02 * 0.02;

                for (const res of convergedResults) {
                    // Ensure x is a 3-component array. SRK result should already be.
                    if (!res.x || res.x.length !== 3) {
                        console.warn("Azeotrope result 'x' is not a 3-component array. Skipping.", res);
                        continue;
                    }

                    let displayErrorNorm: number | undefined;
                    let res_T_K = res.T_K;

                    if (generatingFluidPackage === 'srk') {
                        displayErrorNorm = (res as any).calculatedErrorNormYX;
                        // If T_K from SRK result is 0 or invalid, use initialTempGuess_K as a fallback.
                        // This is a workaround; findAzeotropeSrk should ideally return a correct T_K.
                        if (!res_T_K || res_T_K <= 0) {
                            console.warn(`SRK Azeotrope result has T_K = ${res_T_K}. Using initialTempGuess_K (${initialTempGuess_K.toFixed(2)} K) as fallback.`);
                            res_T_K = initialTempGuess_K; 
                        }
                    } else {
                        // For PR/NRTL/UNIQUAC
                        displayErrorNorm = (res as NrtlAzeotropeResult | PrAzeotropeResult | UniquacAzeotropeResult).errorNorm;
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
                foundAzeotropes = uniqueDirectAzeotropes;
            }
            
            setDirectAzeotropes(foundAzeotropes);
            if (foundAzeotropes.length > 0) {
                console.log(`Directly found ${foundAzeotropes.length} potential azeotrope(s) using ${generatingFluidPackage}:`, foundAzeotropes.map(az => ({x: az.x, T_K: az.T_K, err: az.errorNorm})));
            } else {
                console.log(`No azeotropes found directly using ${generatingFluidPackage} with the given initial guesses.`);
            }
            // --- End Direct Azeotrope Solver Call ---

        } catch (err: any) {
            console.error("Error generating map:", err);
            setError(err.message || "An unknown error occurred.");
        } finally {
            setIsLoading(false);
            setDisplayedFluidPackage(generatingFluidPackage); // Update displayed fluid package when generation finishes
        }
    }, [componentsInput, systemPressure, supabase, fetchCasNumberByName, fetchAndConvertThermData, fluidPackage, plotSortedComponents]); // Removed initialTempGuess_K

    useEffect(() => {
        if (supabase && !initialMapGenerated.current) {
            handleGenerateMap();
            initialMapGenerated.current = true;
        }
    }, [supabase, handleGenerateMap]);

    useEffect(() => {
        if (!residueCurves || residueCurves.length === 0) {
            setCleanedResidueCurves([]);
            console.log(`Plotting useEffect: No residueCurves or empty. Setting cleanedResidueCurves to [].`);
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
        console.log(`Plotting useEffect: residueCurves length: ${residueCurves.length}. Filtered for plotting (${moleFractionThreshold.toExponential()} threshold for ${fluidPackage}): ${filteredForPlotting.length}`);
        if (residueCurves.length > 0 && filteredForPlotting.length === 0) {
            console.warn(`Plotting useEffect: All curves were filtered out by the ${moleFractionThreshold.toExponential()} threshold. Original curves might be too short or stagnant.`);
        }
        if (filteredForPlotting.length > 0 && filteredForPlotting[0].length > 0) { // Added check for filteredForPlotting[0].length
            console.log(`Plotting useEffect: First filtered curve for plotting (first 5 points): `, filteredForPlotting[0].slice(0,5).map(p => ({x: p.x.map(val => typeof val === 'number' ? val.toFixed(4) : 'N/A'), T_K: typeof p.T_K === 'number' ? p.T_K.toFixed(2) : 'N/A'})));
        }

        setCleanedResidueCurves(filteredForPlotting);
    }, [residueCurves, fluidPackage]); // Ensure fluidPackage is a dependency

    useEffect(() => {
        const pressureInUnits = parseFloat(systemPressure); 
        const pressureDisplay = isNaN(pressureInUnits) ? systemPressure : pressureInUnits.toFixed(3);
        const modelName = displayedFluidPackage.charAt(0).toUpperCase() + displayedFluidPackage.slice(1); // Use displayedFluidPackage

        let traces: Data[] = []; 

        // plotSortedComponents is [Lightest, Intermediate, Heaviest]
        // User's desired visual: Right=Lightest, Top=Intermediate, Left=Heaviest
        // User's interpretation of Plotly axes: a=Top, b=Left, c=Right

        // Titles for Plotly's a, b, c axes based on user's interpretation to achieve desired visual
        let titleA = componentsInput[1]?.name || 'Comp 2 (M)'; // a-axis (Top) = Intermediate
        let titleB = componentsInput[2]?.name || 'Comp 3 (H)'; // b-axis (Left) = Heaviest
        let titleC = componentsInput[0]?.name || 'Comp 1 (L)'; // c-axis (Right) = Lightest

        if (plotSortedComponents.length === 3) {
            const formatBp = (bp_K: number | null | undefined) => 
                bp_K ? `${(bp_K - 273.15).toFixed(1)}°C` : 'N/A';

            // a-axis (Top) = Intermediate component (plotSortedComponents[1])
            // Display BP at system pressure if available, else NBP
            const bpA = plotSortedComponents[1].bp_at_Psys_K ?? plotSortedComponents[1].thermData.nbp_K;
            titleA = `${plotSortedComponents[1].name} (M, ${formatBp(bpA)})`;
            
            // b-axis (Left) = Heaviest component (plotSortedComponents[2])
            const bpB = plotSortedComponents[2].bp_at_Psys_K ?? plotSortedComponents[2].thermData.nbp_K;
            titleB = `${plotSortedComponents[2].name} (H, ${formatBp(bpB)})`;

            // c-axis (Right) = Lightest component (plotSortedComponents[0])
            const bpC = plotSortedComponents[0].bp_at_Psys_K ?? plotSortedComponents[0].thermData.nbp_K;
            titleC = `${plotSortedComponents[0].name} (L, ${formatBp(bpC)})`;
        }
        
        const plotTitleColor = currentTheme === 'dark' ? '#e5e7eb' : '#1f2937'; 
        const plotFont = { family: 'Merriweather Sans, Arial, sans-serif', color: currentTheme === 'dark' ? '#e5e7eb' : '#1f2937' }; 
        const axisTitleFont = { ...plotFont, size: 15 };
        const tickFont = { ...plotFont, size: 14 };

        const baseLayout: Partial<Layout> = {
            title: {
                text: `Ternary Residue Curve Map (${modelName} Model, ODE Simulation)`,
                font: { family: 'Merriweather Sans, Arial, sans-serif', size: 18, color: plotTitleColor }, 
            },
            ternary: {
                sum: 1,
                aaxis: { title: { text: titleA, font: axisTitleFont }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: currentTheme === 'dark' ? '#2d3748' : '#cbd5e1' }, // Top (Intermediate)
                baxis: { title: { text: titleB, font: axisTitleFont }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: currentTheme === 'dark' ? '#2d3748' : '#cbd5e1' }, // Left (Heaviest)
                caxis: { title: { text: titleC, font: axisTitleFont }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: currentTheme === 'dark' ? '#2d3748' : '#cbd5e1' }, // Right (Lightest)
                bgcolor: currentTheme === 'dark' ? '#08306b' : '#f0f4f8',
            },
            margin: { l: 70, r: 70, b: 70, t: 90, pad: 5 }, // Standard margins
            autosize: true,
            paper_bgcolor: currentTheme === 'dark' ? '#08306b' : '#f0f4f8', 
            font: plotFont, 
            plot_bgcolor: currentTheme === 'dark' ? '#08306b' : '#f0f4f8', 
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
            setPlotlyData([]); 
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
                    font: plotFont,
                }],
                shapes: [] 
            };
            setPlotlyLayout(layoutUpdate);
            return;
        }

        traces = cleanedResidueCurves.map((curve, index) => {
            const compNameA = plotSortedComponents[1]?.name || 'Comp A';
            const compNameB = plotSortedComponents[2]?.name || 'Comp B';
            const compNameC = plotSortedComponents[0]?.name || 'Comp C';
            // Convert T_K to Celsius directly for customdata for hovertemplate
            const tempsCelsius = curve.map(p => (p.T_K - 273.15).toFixed(1));

            return {
                type: 'scatterternary',
                mode: curve.length > 1 ? 'lines' : 'markers',
                a: curve.map(p => p.x[plotSortedComponents[1].originalIndex]), 
                b: curve.map(p => p.x[plotSortedComponents[2].originalIndex]), 
                c: curve.map(p => p.x[plotSortedComponents[0].originalIndex]), 
                name: `Curve ${index + 1}`,
                showlegend: false,
                line: { width: 1.5, color: '#60a5fa' }, // Example color
                customdata: tempsCelsius, // Store Celsius strings for hover
                hovertemplate:
                    `<b>${compNameA} (M)</b>: %{a:.3f}<br>` + 
                    `<b>${compNameB} (H)</b>: %{b:.3f}<br>` + 
                    `<b>${compNameC} (L)</b>: %{c:.3f}<br>` + 
                    `<b>T</b>: %{customdata} °C<extra></extra>`, // Use customdata directly
                hoverlabel: {
                    font: {
                        family: 'Merriweather Sans, Arial, sans-serif',
                        size: 14,
                        color: plotFont.color
                    }
                }
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
        const plotAreaPxHeight = containerHeightPx - (baseLayout.margin?.t ?? 90) - (baseLayout.margin?.b ?? 90);

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
                            
                        let raw_screen_angle_deg = Math.atan2(-pixel_dy, pixel_dx) * (180 / Math.PI);
                        let final_angle_deg = raw_screen_angle_deg + 90; 

                        // Convert T_K to Celsius for display in hover text
                        const temp_C = centroidPointData.T_K - 273.15;
                        arrowMarkerTraceData.a.push(centroidPointData.x[plotSortedComponents[1].originalIndex]);
                        arrowMarkerTraceData.b.push(centroidPointData.x[plotSortedComponents[2].originalIndex]);
                        arrowMarkerTraceData.c.push(centroidPointData.x[plotSortedComponents[0].originalIndex]);
                        arrowMarkerTraceData.angles.push(final_angle_deg);
                        arrowMarkerTraceData.text.push(`Arrow on Curve ${curveIndex + 1}<br>T=${temp_C.toFixed(1)}°C`);
                        
                        allPlacedArrowIdealCoords.push(currentPointIdealCoords); 
                        placedArrowCountOverall++;

                        if (!firstArrowLogged && curveIndex === 0) { 
                            console.log(`First arrow MARKER (curve ${curveIndex}, arrow ${placedArrowCountOverall}): 
                                CentroidIndex: ${headIndex}, P1_tern: ${point1_for_direction_ternary.map(v => v.toFixed(3)).join(',')}, P2_tern: ${point2_for_direction_ternary.map(v => v.toFixed(3)).join(',')}
                                Ideal dx=${dx_ideal_CALCULATED.toFixed(3)}, dy=${dy_ideal_CALCULATED.toFixed(3)}
                                Pixel dx=${pixel_dx.toFixed(3)}, dy=${pixel_dy.toFixed(3)}, Angle=${final_angle_deg.toFixed(1)}`);
                            firstArrowLogged = true; 
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
                hoverinfo: 'text',
                hoverlabel: {
                    bgcolor: currentTheme === 'dark' ? 'rgba(50,50,50,0.85)' : 'rgba(250,250,250,0.85)', // Adjusted for theme
                    font: {
                        family: 'Merriweather Sans, Arial, sans-serif',
                        size: 14,
                        color: plotFont.color // Already theme-aware
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
                showlegend: false,
            } as Data);
        }

        // --- Directly Found Azeotropes to Plot ---
        if (directAzeotropes.length > 0 && plotSortedComponents.length === 3) {
            const directAzeotropePointsData = {
                a: [] as number[],
                b: [] as number[],
                c: [] as number[],
                text: [] as string[],
            };

            directAzeotropes.forEach(az => {
                if (!az.x || az.x.length !== 3) return;
                directAzeotropePointsData.a.push(az.x[plotSortedComponents[1].originalIndex]);
                directAzeotropePointsData.b.push(az.x[plotSortedComponents[2].originalIndex]);
                directAzeotropePointsData.c.push(az.x[plotSortedComponents[0].originalIndex]);
                directAzeotropePointsData.text.push(
                    `<b>Azeotrope</b><br>` +
                    `${plotSortedComponents[1].name}: ${az.x[plotSortedComponents[1].originalIndex].toFixed(4)}<br>` +
                    `${plotSortedComponents[2].name}: ${az.x[plotSortedComponents[2].originalIndex].toFixed(4)}<br>` +
                    `${plotSortedComponents[0].name}: ${az.x[plotSortedComponents[0].originalIndex].toFixed(4)}<br>` +
                    `T: ${(az.T_K - 273.15).toFixed(2)}°C<br>` +
                    `||y-x||: ${az.errorNorm ? az.errorNorm.toExponential(2) : 'N/A'}`
                );
            });

            if (directAzeotropePointsData.a.length > 0) {
                allTraces.push({
                    type: 'scatterternary',
                    mode: 'markers',
                    a: directAzeotropePointsData.a,
                    b: directAzeotropePointsData.b,
                    c: directAzeotropePointsData.c,
                    text: directAzeotropePointsData.text,
                    hoverinfo: 'text',
                    marker: {
                        symbol: 'star', // Different symbol
                        size: 14,
                        color: '#FFD700', // Gold color
                        line: { color: currentTheme === 'dark' ? '#FFFFFF' : '#000000', width: 1.5 }
                    },
                    name: 'Azeotropes', 
                    showlegend: true,
                    legendgroup: 'azeotropes', // Use a common group name
                    hoverlabel: { 
                        font: { family: 'Merriweather Sans, Arial, sans-serif', size: 12, color: plotFont.color},
                        bgcolor: currentTheme === 'dark' ? 'rgba(40,40,40,0.9)' : 'rgba(240,240,240,0.9)',
                        bordercolor: currentTheme === 'dark' ? '#e5e7eb' : '#1f2937'
                    }
                } as Data);
            }
        }
        // --- End Add Directly Found Azeotropes to Plot ---


        const finalLayout = {
            ...baseLayout,
            shapes: [], // Clear old shapes if any, or manage them if other shapes are needed
            annotations: [...(baseLayout.annotations || [])], 
            title: { 
                ...(typeof baseLayout.title === 'object' ? baseLayout.title : { text: '', font: { family: 'Merriweather Sans, Arial, sans-serif', size: 18, color: plotTitleColor} } ), 
                text: `Ternary Residue Curve Map @ ${pressureDisplay} bar (${modelName} Model)` 
            },
            legend: { // Ensure legend is configured
                font: plotFont,
                bgcolor: currentTheme === 'dark' ? 'rgba(8, 48, 107, 0.8)' : 'rgba(220, 230, 240, 0.8)', // Semi-transparent background
                bordercolor: currentTheme === 'dark' ? '#4b5563' : '#9ca3af',
                borderwidth: 1,
                tracegroupgap: 5, // Gap between legend groups
            }
        };

        setPlotlyData(allTraces);
        setPlotlyLayout(finalLayout);

    }, [cleanedResidueCurves, componentsInput, systemPressure, plotSortedComponents, isLoading, currentTheme, displayedFluidPackage, directAzeotropes]); // Depend on directAzeotropes

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
                            <div className="space-y-2">
                                {componentsInput.map((component, index) => (
                                    <Card className="flex-1" key={index}>
                                        <CardContent className="space-y-3">
                                            <div className="relative">
                                                <Label htmlFor={`name-${index}`}>Name</Label>
                                                <Input
                                                    id={`name-${index}`}
                                                    ref={el => { inputRefs.current[index] = el; }}
                                                    value={component.name}
                                                    onChange={e => handleComponentInputChange(index, 'name', e.target.value)}
                                                    placeholder={`e.g. Acetone`}
                                                    autoComplete="off"
                                                    onFocus={() => {
                                                        setActiveSuggestionIndex(index);
                                                        if (componentsInput[index].name.trim()) {
                                                            debouncedFetchComponentSuggestions(componentsInput[index].name, index);
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
                                            </div>
                                        </CardContent>
                                    </Card>
                                ))}
                            </div>
                            {/* Pressure on same line */}
                            <div className="flex items-center space-x-2">
                                <Label htmlFor="systemPressure" className="whitespace-nowrap">
                                    Pressure (bar)
                                </Label>
                                <Input
                                    id="systemPressure"
                                    type="number"
                                    value={systemPressure}
                                    onChange={e => setSystemPressure(e.target.value)}
                                    placeholder="1.0"
                                />
                            </div>
                            {/* Activity model on same line */}
                            <div className="flex items-center space-x-2">
                                <Label htmlFor="fluidPackage" className="whitespace-nowrap">
                                    Activity Model
                                </Label>
                                <Select
                                    value={fluidPackage}
                                    onValueChange={v => setFluidPackage(v as FluidPackageTypeResidue)}
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
                </div>

                {/* Right panel: Plot */}
                <div className="lg:col-span-2">
                    <Card>
                        <CardHeader>
                            {/* <CardTitle>Ternary Residue Curve Map</CardTitle> */} {/* Removed this line */}
                        </CardHeader>
                        <CardContent className="bg-[#08306b] p-1 rounded-md">
                            <div 
                                className="relative h-[600px]" 
                                ref={setPlotContainerRef} // Assign the state setter to the ref prop
                            >
                                <Plot 
                                    data={plotlyData} 
                                    layout={plotlyLayout} 
                                    useResizeHandler 
                                    style={{ width: '100%', height: '100%' }} 
                                />
                            </div>
                        </CardContent>
                    </Card>
                </div>
            </div>
        </div>
    );
}

