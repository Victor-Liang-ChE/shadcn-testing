'use client';

import React, { useState, useEffect, useRef, useCallback } from 'react';
import dynamic from 'next/dynamic'; // Import dynamic
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { createClient } from '@supabase/supabase-js';
import type { SupabaseClient } from '@supabase/supabase-js';
import { CompoundData, WilsonPureComponentParams as WilsonPureParams } from '@/lib/vle-types';
import { 
    simulateSingleResidueCurveODE, 
    ResidueCurveODE, 
    TernaryWilsonParams
} from '@/lib/residue-curves-ode-wilson';
import { fetchWilsonInteractionParams } from '@/lib/vle-calculations-wilson'; 

import { fetchAndConvertThermData, FetchedCompoundThermData } from '@/lib/antoine-utils'; 
import type { Data, Layout } from 'plotly.js'; 

// Dynamically import Plotly to avoid SSR issues
const Plot = dynamic(() => import('react-plotly.js'), {
  ssr: false,
  loading: () => <div className="flex justify-center items-center h-[600px]"><p>Loading plot...</p></div>
});

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
    thermData: FetchedCompoundThermData; // This should include Antoine, NBP, and V_L_m3mol
    originalIndex: number; // 0, 1, or 2 based on input order
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
    const [residueCurves, setResidueCurves] = useState<ResidueCurveODE[]>([]); // NEW: Using ResidueCurveODE
    const [cleanedResidueCurves, setCleanedResidueCurves] = useState<ResidueCurveODE[]>([]); // For filtered curves
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [plotlyData, setPlotlyData] = useState<Data[]>([]);
    const [plotlyLayout, setPlotlyLayout] = useState<Partial<Layout>>({});
    const [plotSortedComponents, setPlotSortedComponents] = useState<ProcessedComponentData[]>([]);
    const initialMapGenerated = useRef(false); // Ref to track initial generation
    const [currentTheme, setCurrentTheme] = useState<'dark' | 'light'>('dark'); // Default to dark
    const plotContainerRef = useRef<HTMLDivElement>(null); // Ref for the plot container div

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
                .ilike('name', `%${query.trim()}%`)
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
        };
    }, [activeSuggestionIndex]);

    const generateStartingPoints = (numPerEdge: number = 8): number[][] => {
        const points: number[][] = [];
        const edgeEpsilon = 0.005; // Start slightly off the actual edge for numerical stability

        // Helper to add point if valid and normalized
        const addNormalizedPoint = (p: number[]) => {
            const sum = p[0] + p[1] + p[2];
            if (sum <= 1e-6) return; // Avoid division by zero or near-zero sum
            const normP = [p[0] / sum, p[1] / sum, p[2] / sum];
            // Ensure points are not extremely close to pure components after normalization
            const minFrac = 1e-3;
            if (normP.every(val => val >= minFrac && val <= 1.0 - minFrac)) {
                 // Avoid adding exact duplicates based on a reasonable precision
                if (!points.some(existingP => 
                    Math.abs(existingP[0] - normP[0]) < 1e-4 &&
                    Math.abs(existingP[1] - normP[1]) < 1e-4 &&
                    Math.abs(existingP[2] - normP[2]) < 1e-4
                )) {
                    points.push(normP);
                }
            }
        };
        
        // Points along edges (slightly offset)
        // Edge 0-1 (x2 is small)
        for (let i = 1; i <= numPerEdge; i++) {
            const x0_ratio = i / (numPerEdge + 1.0);
            addNormalizedPoint([x0_ratio * (1.0 - edgeEpsilon), (1.0 - x0_ratio) * (1.0 - edgeEpsilon), edgeEpsilon]);
        }
        // Edge 0-2 (x1 is small)
        for (let i = 1; i <= numPerEdge; i++) {
            const x0_ratio = i / (numPerEdge + 1.0);
            addNormalizedPoint([x0_ratio * (1.0 - edgeEpsilon), edgeEpsilon, (1.0 - x0_ratio) * (1.0 - edgeEpsilon)]);
        }
        // Edge 1-2 (x0 is small)
        for (let i = 1; i <= numPerEdge; i++) {
            const x1_ratio = i / (numPerEdge + 1.0);
            addNormalizedPoint([edgeEpsilon, x1_ratio * (1.0 - edgeEpsilon), (1.0 - x1_ratio) * (1.0 - edgeEpsilon)]);
        }

        // Add a more comprehensive set of internal points
        const internalPointsSeed = [
            [0.333, 0.333, 0.334],
            [0.5, 0.25, 0.25], [0.25, 0.5, 0.25], [0.25, 0.25, 0.5],
            [0.1, 0.1, 0.8], [0.1, 0.8, 0.1], [0.8, 0.1, 0.1],
            [0.2, 0.2, 0.6], [0.2, 0.6, 0.2], [0.6, 0.2, 0.2],
            [0.4, 0.4, 0.2], [0.4, 0.2, 0.4], [0.2, 0.4, 0.4],
            [0.1, 0.2, 0.7], [0.1, 0.7, 0.2], [0.2, 0.1, 0.7], [0.2, 0.7, 0.1], [0.7, 0.1, 0.2], [0.7, 0.2, 0.1],
            [0.15, 0.35, 0.5], [0.15, 0.5, 0.35], [0.35, 0.15, 0.5], [0.35, 0.5, 0.15], [0.5, 0.15, 0.35], [0.5, 0.35, 0.15],
            [0.2, 0.3, 0.5], [0.2, 0.5, 0.3], [0.3, 0.2, 0.5], [0.3, 0.5, 0.2], [0.5, 0.2, 0.3], [0.5, 0.3, 0.2],
            // Add points closer to the center if needed
            [0.3, 0.3, 0.4], [0.3, 0.4, 0.3], [0.4, 0.3, 0.3],
        ];
        internalPointsSeed.forEach(p => addNormalizedPoint(p));
        
        return points;
    };

    const fetchCasNumberByName = async (supabaseClient: SupabaseClient, name: string): Promise<string> => {
        if (!name.trim()) {
            throw new Error("Component name cannot be empty.");
        }
        const { data, error } = await supabaseClient
            .from('compounds')
            .select('cas_number')
            .ilike('name', `%${name.trim()}%`)
            .limit(1)
            .single();

        if (error || !data || !data.cas_number) {
            console.error(`Error fetching CAS for name '${name}':`, error);
            throw new Error(`CAS number not found for component "${name}". Ensure the name is correct and exists in the database.`);
        }
        return data.cas_number;
    };
    
    const handleGenerateMap = useCallback(async () => {
        if (!supabase) {
            setError("Supabase client is not initialized. Cannot fetch parameters.");
            return;
        }
        setIsLoading(true);
        setError(null);

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

            const processedComponents: ProcessedComponentData[] = componentsInput.map((input, index) => ({
                name: input.name,
                casNumber: fetchedCasNumbers[index],
                thermData: fetchedThermDataArray[index],
                originalIndex: index,
            }));
            
            if (processedComponents.some(pc => pc.thermData.V_L_m3mol == null || pc.thermData.V_L_m3mol <= 0)) {
                throw new Error("Liquid molar volume (V_L_m3mol) missing or invalid for one or more components.");
            }
             if (processedComponents.some(pc => !pc.thermData.antoine)) {
                throw new Error("Antoine parameters missing for one or more components.");
            }

            const sortedForPlot = [...processedComponents].sort((a, b) => {
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
                    wilsonParams: { V_L_m3mol: pc.thermData.V_L_m3mol! } as WilsonPureParams,
                    unifacGroups: [], 
            }));

            const cas = compoundsForBackend.map(c => c.cas_number!);
            const params01 = await fetchWilsonInteractionParams(supabase, cas[0], cas[1]);
            const params02 = await fetchWilsonInteractionParams(supabase, cas[0], cas[2]);
            const params12 = await fetchWilsonInteractionParams(supabase, cas[1], cas[2]);

            const ternaryWilsonParams: TernaryWilsonParams = {
                a01_J_mol: params01.a12_J_mol, a10_J_mol: params01.a21_J_mol,
                a02_J_mol: params02.a12_J_mol, a20_J_mol: params02.a21_J_mol,
                a12_J_mol: params12.a12_J_mol, a21_J_mol: params12.a21_J_mol,
            };
            
            const avgNBP = processedComponents
                .map(p => p.thermData.nbp_K)
                .filter(nbp => nbp !== null) as number[];
            const initialTempGuess_K = avgNBP.length > 0 
                ? avgNBP.reduce((sum, val) => sum + val, 0) / avgNBP.length 
                : DEFAULT_INITIAL_T_K;

            const startingCompositions = generateStartingPoints(8); // Increased numPerEdge, internal points are fixed in the function
            const allCurvesPromises: Promise<ResidueCurveODE | null>[] = [];

            for (const start_x_3comp of startingCompositions) {
                allCurvesPromises.push(
                    simulateSingleResidueCurveODE(
                        start_x_3comp,
                        P_system_Pa,
                        compoundsForBackend,
                        ternaryWilsonParams,
                        0.02, // dt_frac
                        1000,  // max_steps increased from 300
                        initialTempGuess_K 
                    )
                );
            }
            
            const resolvedCurves = await Promise.all(allCurvesPromises);
            const validCurves = resolvedCurves.filter(curve => curve !== null && curve.length > 1) as ResidueCurveODE[];

            // Apply new cleaning logic to validCurves before setting state
            const cleanedCurves = validCurves.map(curve => {
                if (curve.length < 2) return curve;
                const uniquePointsCurve: ResidueCurveODE = [curve[0]]; // Start with the first point
                for (let k = 1; k < curve.length; k++) {
                    const prevP = uniquePointsCurve[uniquePointsCurve.length - 1].x;
                    const currP = curve[k].x;
                    const isDuplicate =
                        Math.abs(prevP[0] - currP[0]) < 1e-7 &&
                        Math.abs(prevP[1] - currP[1]) < 1e-7 &&
                        Math.abs(prevP[2] - currP[2]) < 1e-7;
                    if (!isDuplicate) {
                        uniquePointsCurve.push(curve[k]);
                    }
                }
                return uniquePointsCurve;
            }).filter(curve => curve.length > 0); // Ensure no empty curves if all points were identical


            if (cleanedCurves.length === 0 && startingCompositions.length > 0) {
                const errorMsg = "No valid residue curves generated (or all points were duplicates). This might be due to: \n1. Calculation failures for all starting points (e.g., bubble point errors). \n2. Missing or incorrect interaction parameters (e.g., Aij=0 for some pairs). \n3. Issues with pure component parameters. \nPlease check console logs and database for parameter completeness.";
                console.warn(errorMsg);
                setError(errorMsg); // Display a more informative error to the user
            } else if (cleanedCurves.length < startingCompositions.length * 0.5) { // If many curves failed
                 console.warn(`Many residue curve calculations failed or resulted in very short/duplicate curves (${startingCompositions.length - cleanedCurves.length} out of ${startingCompositions.length}). Check parameters and console for errors.`);
                 // Optionally, append a warning to the existing error or set a non-critical warning message
            }
            setResidueCurves(cleanedCurves); // Set the 1e-7 cleaned curves

        } catch (err: any) {
            console.error("Error generating map:", err);
            setError(err.message || "An unknown error occurred.");
        } finally {
            setIsLoading(false);
        }
    }, [componentsInput, systemPressure, supabase, fetchCasNumberByName, fetchAndConvertThermData]);

    useEffect(() => {
        if (supabase && !initialMapGenerated.current) {
            handleGenerateMap();
            initialMapGenerated.current = true;
        }
    }, [supabase, handleGenerateMap]);

    useEffect(() => {
        if (!residueCurves || residueCurves.length === 0) {
            setCleanedResidueCurves([]);
            return;
        }

        const moleFractionThreshold = 1e-5; // Threshold for considering points different

        const filterCurves = (curves: ResidueCurveODE[]): ResidueCurveODE[] => {
            return curves.map(curve => {
                if (curve.length < 2) return curve; // Keep short curves (0 or 1 point) as is

                const cleanedCurve: ResidueCurveODE = [curve[0]]; // Always include the first point
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

        setCleanedResidueCurves(filterCurves(residueCurves));
    }, [residueCurves]);

    useEffect(() => {
        const pressureInUnits = parseFloat(systemPressure); 
        const pressureDisplay = isNaN(pressureInUnits) ? systemPressure : pressureInUnits.toFixed(3);

        let traces: Data[] = []; 

        // plotSortedComponents is [Lightest, Intermediate, Heaviest]
        // User's desired visual: Right=Lightest, Top=Intermediate, Left=Heaviest
        // User's interpretation of Plotly axes: a=Top, b=Left, c=Right

        // Titles for Plotly's a, b, c axes based on user's interpretation to achieve desired visual
        let titleA = componentsInput[1]?.name || 'Comp 2 (M)'; // a-axis (Top) = Intermediate
        let titleB = componentsInput[2]?.name || 'Comp 3 (H)'; // b-axis (Left) = Heaviest
        let titleC = componentsInput[0]?.name || 'Comp 1 (L)'; // c-axis (Right) = Lightest

        if (plotSortedComponents.length === 3) {
            const formatNbp = (nbp_K: number | null | undefined) => 
                nbp_K ? `${(nbp_K - 273.15).toFixed(1)}°C` : 'N/A';

            // a-axis (Top) = Intermediate component (plotSortedComponents[1])
            titleA = `${plotSortedComponents[1].name} (M, ${formatNbp(plotSortedComponents[1].thermData.nbp_K)})`;
            // b-axis (Left) = Heaviest component (plotSortedComponents[2])
            titleB = `${plotSortedComponents[2].name} (H, ${formatNbp(plotSortedComponents[2].thermData.nbp_K)})`;
            // c-axis (Right) = Lightest component (plotSortedComponents[0])
            titleC = `${plotSortedComponents[0].name} (L, ${formatNbp(plotSortedComponents[0].thermData.nbp_K)})`;
        }
        
        const plotTitleColor = currentTheme === 'dark' ? '#e5e7eb' : '#1f2937'; 
        const plotFont = { family: 'Merriweather Sans, Arial, sans-serif', color: '#e5e7eb' }; 
        const axisTitleFont = { ...plotFont, size: 15 };
        const tickFont = { ...plotFont, size: 14 };

        const baseLayout: Partial<Layout> = {
            title: {
                text: 'Ternary Residue Curve Map (Wilson Model, ODE Simulation)',
                font: { family: 'Merriweather Sans, Arial, sans-serif', size: 18, color: plotTitleColor }, 
            },
            ternary: {
                sum: 1,
                aaxis: { title: { text: titleA, font: axisTitleFont }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: '#2d3748' }, // Top (Intermediate)
                baxis: { title: { text: titleB, font: axisTitleFont }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: '#2d3748' }, // Left (Heaviest)
                caxis: { title: { text: titleC, font: axisTitleFont }, min: 0, max: 1, tickformat: '.1f', tickfont: tickFont, linecolor: '#4b5563', gridcolor: '#2d3748' }, // Right (Lightest)
                bgcolor: '#08306b',
            },
            margin: { l: 70, r: 70, b: 70, t: 90, pad: 5 }, // Standard margins
            autosize: true,
            paper_bgcolor: '#08306b', 
            font: plotFont, 
            plot_bgcolor: '#08306b', 
            annotations: [], 
            shapes: [], 
        };

        // Get container dimensions for coordinate transformation
        let containerWidthPx = 800; // Default/fallback
        let containerHeightPx = 600; // Default/fallback
        if (plotContainerRef.current) {
            containerWidthPx = plotContainerRef.current.offsetWidth;
            containerHeightPx = plotContainerRef.current.offsetHeight;
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
                    text: 'Ternary Residue Curve Map (No data or NBP info pending)' // Specific title for this case
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

        traces = cleanedResidueCurves.map((curve, index) => { // Use cleanedResidueCurves
            const getFractionForAxis = (point_x_array: number[], targetOriginalIndex: number): number => {
                return point_x_array[targetOriginalIndex];
            };

            return {
                type: 'scatterternary',
                mode: 'lines',
                a: curve.map(p => getFractionForAxis(p.x, plotSortedComponents[1].originalIndex)), 
                b: curve.map(p => getFractionForAxis(p.x, plotSortedComponents[2].originalIndex)), 
                c: curve.map(p => getFractionForAxis(p.x, plotSortedComponents[0].originalIndex)), 
                name: `Curve ${index + 1}`,
                showlegend: false,          // <--- ADDED THIS LINE
                line: { width: 1.5, color: '#60a5fa' }, // Restored color
                customdata: curve.map(p => p.T_K),  // Restored customdata
                hovertemplate:  // Restored hovertemplate
                    `<b>${plotSortedComponents[1].name} (M)</b>: %{a:.3f}<br>` + 
                    `<b>${plotSortedComponents[2].name} (H)</b>: %{b:.3f}<br>` + 
                    `<b>${plotSortedComponents[0].name} (L)</b>: %{c:.3f}<br>` + 
                    `<b>T</b>: %{customdata:.1f} K<extra></extra>`,
            };
        });

        const arrowMarkerTraceData = { // This is already being reset
            a: [] as number[],
            b: [] as number[],
            c: [] as number[],
            angles: [] as number[],
            text: [] as string[], 
        };
        // --- Globals for arrow placement ---
        const allPlacedArrowIdealCoords: { x: number; y: number }[] = [];
        // ADJUST THESE:
        const minSqDistBetweenPlacedArrows_ideal = (0.05) * (0.05); // Smaller value for denser arrows
        const maxTotalArrowsOnPlot = 250;                             // Allow more arrows overall
        let placedArrowCountOverall = 0;


        let firstArrowLogged = false; 
        const targetArrowPixelSize = 10; 
        const arrowMarkerColor = '#FFFFFF';

        const plotAreaPxWidth = containerWidthPx - (baseLayout.margin?.l ?? 70) - (baseLayout.margin?.r ?? 70);
        const plotAreaPxHeight = containerHeightPx - (baseLayout.margin?.t ?? 90) - (baseLayout.margin?.b ?? 90);

        if (plotSortedComponents.length === 3 && plotContainerRef.current && actualRenderedTriangleBase_paper > 0 && actualRenderedTriangleHeight_paper > 0 && plotAreaPxWidth > 0 && plotAreaPxHeight > 0) {
            cleanedResidueCurves.forEach((curve, curveIndex) => { 
                if (curve.length < 2) return;
                if (placedArrowCountOverall >= maxTotalArrowsOnPlot) return;

                let lastCheckedPointIdealCoords: { x: number; y: number } | null = null;
                const minIdealDistSqForCheckingNextPointOnCurve = (0.03) * (0.03); 

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

                        arrowMarkerTraceData.a.push(centroidPointData.x[plotSortedComponents[1].originalIndex]);
                        arrowMarkerTraceData.b.push(centroidPointData.x[plotSortedComponents[2].originalIndex]);
                        arrowMarkerTraceData.c.push(centroidPointData.x[plotSortedComponents[0].originalIndex]);
                        arrowMarkerTraceData.angles.push(final_angle_deg);
                        arrowMarkerTraceData.text.push(`Arrow on Curve ${curveIndex + 1}<br>T=${centroidPointData.T_K.toFixed(1)}K`);
                        
                        allPlacedArrowIdealCoords.push(currentPointIdealCoords); 
                        placedArrowCountOverall++;

                        if (!firstArrowLogged && curveIndex === 0) { 
                            console.log(`First arrow MARKER (curve ${curveIndex}, arrow ${placedArrowCountOverall}): 
                                CentroidIndex: ${headIndex}, P1_tern: ${point1_for_direction_ternary.map(v => v.toFixed(3)).join(',')}, P2_tern: ${point2_for_direction_ternary.map(v => v.toFixed(3)).join(',')}
                                Ideal dx=${dx_ideal_CALCULATED.toFixed(3)}, dy=${dy_ideal_CALCULATED.toFixed(3)}
                                Pixel dx=${pixel_dx.toFixed(3)}, dy=${pixel_dy.toFixed(3)}, Angle=${final_angle_deg.toFixed(1)}`);
                            firstArrowLogged = true; 
                        }
                    }
                }
            });
        }
        
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
                marker: {
                    symbol: 'triangle-up',
                    color: arrowMarkerColor,
                    size: targetArrowPixelSize,
                    angle: arrowMarkerTraceData.angles,
                    line: { width: 0 }
                },
                name: 'Direction Arrows',
                showlegend: false,
            } as Data);
        }

        const finalLayout = {
            ...baseLayout,
            shapes: [], // Clear old shapes if any, or manage them if other shapes are needed
            annotations: [...(baseLayout.annotations || [])], 
            title: { 
                ...(typeof baseLayout.title === 'object' ? baseLayout.title : { text: '', font: { family: 'Merriweather Sans, Arial, sans-serif', size: 18, color: plotTitleColor} } ), 
                text: `Ternary Residue Curve Map @ ${pressureDisplay} bar` 
            },
            legend: { font: plotFont }
        };

        setPlotlyData(allTraces);
        setPlotlyLayout(finalLayout);

    }, [cleanedResidueCurves, componentsInput, systemPressure, plotSortedComponents, isLoading, currentTheme]); // Depend on cleanedResidueCurves

    return (
        <div className="container mx-auto p-4">
            <Card className="mb-6">
                <CardHeader><CardTitle>Ternary Residue Curve Map Inputs (Wilson ODE)</CardTitle></CardHeader>
                <CardContent className="space-y-6">
                    <div className="flex flex-col md:flex-row items-stretch gap-2 md:gap-4">
                        {componentsInput.map((component, index) => (
                            <Card className="flex-1" key={index}>
                                <CardHeader className="pb-2 pt-4">
                                    <CardTitle className="text-base">Component {index + 1}</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-3">
                                    <div className="relative">
                                        <Label htmlFor={`name-${index}`}>Name</Label>
                                        <Input
                                            id={`name-${index}`}
                                            ref={el => { inputRefs.current[index] = el; }}
                                            value={component.name}
                                            onChange={(e) => handleComponentInputChange(index, 'name', e.target.value)}
                                            placeholder={`e.g., Acetone`}
                                            autoComplete="off"
                                            onFocus={() => {
                                                setActiveSuggestionIndex(index);
                                                // Fetch suggestions if input is not empty and suggestions are not already shown
                                                if (component.name.trim() && (!showSuggestions[index] || componentSuggestions[index].length === 0) ) {
                                                    debouncedFetchComponentSuggestions(component.name, index);
                                                } else if (component.name.trim() && componentSuggestions[index].length > 0) {
                                                     setShowSuggestions(prev => {
                                                        const updated = [...prev];
                                                        updated[index] = true; 
                                                        return updated;
                                                    });
                                                }
                                            }}
                                        />
                                        {showSuggestions[index] && componentSuggestions[index].length > 0 && (
                                            <ul 
                                                ref={el => { suggestionsContainerRefs.current[index] = el; }}
                                                className={`absolute z-10 rounded-md shadow-md mt-1 max-h-40 overflow-y-auto w-full
                                                ${currentTheme === 'dark' 
                                                    ? 'bg-gray-800 border border-gray-700 text-gray-200' 
                                                    : 'bg-white border border-gray-300 text-gray-900'}`}
                                            >
                                                {componentSuggestions[index].map((suggestion, suggestionIndex) => (
                                                    <li
                                                        key={suggestionIndex}
                                                        className={`px-4 py-2 cursor-pointer 
                                                            ${currentTheme === 'dark' 
                                                                ? 'hover:bg-gray-700' 
                                                                : 'hover:bg-gray-100'}`}
                                                        onMouseDown={() => handleSuggestionClick(index, suggestion)} // Use onMouseDown to fire before onBlur
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
                    <div>
                        <Label htmlFor="systemPressure">System Pressure (bar)</Label>
                        <Input 
                            id="systemPressure" 
                            value={systemPressure} 
                            onChange={(e) => setSystemPressure(e.target.value)} 
                            placeholder="e.g., 1.01325" 
                            type="number"
                        />
                    </div>
                    <Button onClick={handleGenerateMap} disabled={isLoading || !supabase}>
                        {isLoading ? 'Generating (Wilson ODE)...' : 'Generate Residue Map (Wilson ODE)'}
                    </Button>
                    {error && <p className="text-red-500">{error}</p>}
                    {!supabase && <p className="text-orange-500">Supabase client not configured. Parameters cannot be fetched.</p>}
                </CardContent>
            </Card>

            <Card>
                <CardHeader>
                    <CardTitle style={{ fontFamily: 'Merriweather Sans', color: '#08306b' }}>
                        Residue Curve Map (Wilson Model, ODE Simulation)
                    </CardTitle>
                </CardHeader>
                <CardContent className="bg-[#08306b] p-1 rounded-md">
                  <div ref={plotContainerRef} style={{ width: '100%', height: '600px' }}> {/* Added ref here */}
                    {(isLoading && cleanedResidueCurves.length === 0) ? ( // Use cleanedResidueCurves
                        <div className="text-center p-10 text-slate-300 min-h-[600px] flex items-center justify-center" style={{fontFamily: 'Merriweather Sans'}}>
                            <div>Generating initial map...</div>
                        </div>
                    ) : (cleanedResidueCurves.length > 0 || isLoading) ? ( // Use cleanedResidueCurves
                        <Plot
                            data={plotlyData}
                            layout={plotlyLayout}
                            style={{ width: '100%', height: '100%' }} // Changed to 100% to fill parent div
                            useResizeHandler={true}
                            config={{ responsive: true, displaylogo: false }}
                        />
                    ) : (
                        <div className="text-center p-10 text-slate-300 min-h-[600px] flex items-center justify-center" style={{fontFamily: 'Merriweather Sans'}}>
                            <div>
                                <p>Enter component names and system pressure, then click "Generate Residue Map" to see the plot.</p>
                                <p>CAS, Antoine, V<sub>L</sub> (molar vol), NBP will be fetched. Wilson binary params also fetched.</p>
                                <p className="mt-2">
                                    Plotly Ternary Axes: a (Intermediate NBP), b (Highest NBP), c (Lowest NBP).
                                </p>
                                <p className="text-xs mt-1">Ensure `fetchAndConvertThermData` provides liquid molar volumes (V_L_m3mol).</p>
                            </div>
                        </div>
                    )}
                   </div> {/* Closing plotContainerRef div */}
                </CardContent>
            </Card>
        </div>
    );
}

const MIN_MOLE_FRACTION = 1e-9;

