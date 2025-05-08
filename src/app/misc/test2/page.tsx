'use client';

import React, { useState, useEffect, useRef } from 'react';
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { createClient } from '@supabase/supabase-js';
import type { SupabaseClient } from '@supabase/supabase-js';
import { CompoundData } from '@/lib/vle-types'; // Adjust path as needed
import {
    fetchTernaryNrtlParameters,
    traceResidueCurve,
    ResidueCurve,
    ResidueCurvePoint,
    TernaryNrtlParameters
} from '@/lib/residue-curves-nrtl'; // Adjust path as needed
import { fetchAndConvertThermData, FetchedCompoundThermData } from '@/lib/antoine-utils'; // Updated import
import Plot from 'react-plotly.js'; // Import Plotly
import type { Data, Layout } from 'plotly.js'; // Import Plotly types

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

interface ComponentInputState {
    name: string; 
}

// To store component data along with its original index and NBP for sorting
interface ProcessedComponentData {
    name: string;
    casNumber: string;
    thermData: FetchedCompoundThermData;
    originalIndex: number; // 0, 1, or 2 based on input order
}

const initialComponentState = (): ComponentInputState => ({
    name: '',
});

export default function TernaryResidueMapPage() {
    const [componentsInput, setComponentsInput] = useState<ComponentInputState[]>([
        { name: 'Acetone' },
        { name: 'Methanol' },
        { name: 'Water' },
    ]);
    const [systemPressure, setSystemPressure] = useState<string>('101.325'); // System pressure in kPa
    const [residueCurves, setResidueCurves] = useState<ResidueCurve[]>([]);
    const [isLoading, setIsLoading] = useState<boolean>(false);
    const [error, setError] = useState<string | null>(null);
    const [plotlyData, setPlotlyData] = useState<Data[]>([]);
    const [plotlyLayout, setPlotlyLayout] = useState<Partial<Layout>>({});
    // State to store NBP-sorted component data for plot labeling and data mapping
    const [plotSortedComponents, setPlotSortedComponents] = useState<ProcessedComponentData[]>([]);

    const handleComponentInputChange = (index: number, field: keyof ComponentInputState, value: string) => {
        const newInputs = [...componentsInput];
        newInputs[index] = { ...newInputs[index], [field]: value };
        setComponentsInput(newInputs);
    };

    const generateStartingPoints = (numPerEdge: number = 5, numInternal: number = 3): [number, number][] => {
        const points: [number, number][] = [];
        const step = 1.0 / (numPerEdge + 1);

        // Edge 1-3 (x2 = 0)
        for (let i = 1; i <= numPerEdge; i++) points.push([i * step, 0]);
        // Edge 2-3 (x1 = 0)
        for (let i = 1; i <= numPerEdge; i++) points.push([0, i * step]);
        // Edge 1-2 (x3 = 0 => x1 + x2 = 1)
        for (let i = 1; i <= numPerEdge; i++) points.push([i * step, 1 - (i * step)]);
        
        // Some internal points (crude generation)
        // Example: a few points on lines from vertices to opposite midpoints or centroid
        if (numInternal > 0) {
            points.push([0.333, 0.333]); // Centroid-ish
            points.push([0.5, 0.25]);
            points.push([0.25, 0.5]);
            points.push([0.1, 0.45]);
            points.push([0.45, 0.1]);
        }
        // Filter out points too close to vertices or duplicate
        return points.filter(p => p[0] > 1e-3 && p[1] > 1e-3 && (p[0] + p[1]) < (1 - 1e-3));
    };

    const fetchCasNumberByName = async (supabaseClient: SupabaseClient, name: string): Promise<string> => {
        if (!name.trim()) {
            throw new Error("Component name cannot be empty.");
        }
        const { data, error } = await supabaseClient
            .from('compounds')
            .select('cas_number')
            .ilike('name', `%${name.trim()}%`) // Use ilike for case-insensitive partial match, or eq for exact match
            .limit(1)
            .single();

        if (error || !data || !data.cas_number) {
            console.error(`Error fetching CAS for name '${name}':`, error);
            throw new Error(`CAS number not found for component "${name}". Ensure the name is correct and exists in the database.`);
        }
        return data.cas_number;
    };
    
    const handleGenerateMap = async () => {
        if (!supabase) {
            setError("Supabase client is not initialized. Cannot fetch parameters.");
            return;
        }
        setIsLoading(true);
        setError(null);
        setResidueCurves([]);

        try {
            const P_system_Pa = parseFloat(systemPressure) * 1000; // Convert kPa to Pa
            if (isNaN(P_system_Pa) || P_system_Pa <= 0) {
                throw new Error("Invalid system pressure input. Must be a positive number.");
            }

            // Fetch CAS numbers from names
            const casNumbersPromises = componentsInput.map(input => 
                fetchCasNumberByName(supabase!, input.name)
            );
            const fetchedCasNumbers = await Promise.all(casNumbersPromises);

            // Fetch Thermo data (Antoine + NBP) for each component using fetched CAS numbers
            const thermoDataPromises = fetchedCasNumbers.map(casn => 
                fetchAndConvertThermData(supabase!, casn) // Use new function
            );
            const fetchedThermDataArray = await Promise.all(thermoDataPromises);

            // Prepare component data with NBP for sorting and general use
            const processedComponents: ProcessedComponentData[] = componentsInput.map((input, index) => ({
                name: input.name,
                casNumber: fetchedCasNumbers[index],
                thermData: fetchedThermDataArray[index],
                originalIndex: index,
            }));

            // Sort components by NBP for plot axes (lowest to highest)
            // Handle cases where NBP might be null
            const sortedForPlot = [...processedComponents].sort((a, b) => {
                if (a.thermData.nbp_K === null && b.thermData.nbp_K === null) return 0;
                if (a.thermData.nbp_K === null) return 1; // Put nulls at the end (higher NBP)
                if (b.thermData.nbp_K === null) return -1;
                return a.thermData.nbp_K - b.thermData.nbp_K;
            });
            
            if (sortedForPlot.some(c => c.thermData.nbp_K === null)) {
                setError("One or more components missing Normal Boiling Point data. Plot axes may not represent boiling point order.");
                // Proceeding, but the user is warned. Or, could throw new Error here.
            }
            setPlotSortedComponents(sortedForPlot);

            // Create CompoundData array in the original input order for backend functions
            const compoundsForBackend: CompoundData[] = processedComponents
                .sort((a,b) => a.originalIndex - b.originalIndex) // Ensure original order
                .map(pc => ({
                    name: pc.name,
                    cas_number: pc.casNumber,
                    antoine: pc.thermData.antoine,
                    unifacGroups: [], // Or other properties as needed by backend
            }));

            const nrtlParams: TernaryNrtlParameters = await fetchTernaryNrtlParameters(
                supabase,
                compoundsForBackend[0].cas_number!, // cas_number is now guaranteed by fetchCasNumberByName
                compoundsForBackend[1].cas_number!,
                compoundsForBackend[2].cas_number!
            );

            const startingPoints = generateStartingPoints(7, 5); 
            const allCurves: ResidueCurve[] = [];

            for (const startPair of startingPoints) {
                // traceResidueCurve expects compounds in a specific order (matching x1, x2, x3)
                // and P_system_Pa (assuming backend is updated for constant P)
                const curve = traceResidueCurve(startPair, compoundsForBackend, P_system_Pa, nrtlParams, 0.005, 2000);
                if (curve.length > 1) allCurves.push(curve);
            }
            setResidueCurves(allCurves);

        } catch (err: any) {
            console.error("Error generating map:", err);
            setError(err.message || "An unknown error occurred.");
        } finally {
            setIsLoading(false);
        }
    };

    useEffect(() => {
        const pressureInUnits = parseFloat(systemPressure); 
        const pressureDisplay = isNaN(pressureInUnits) ? systemPressure : pressureInUnits.toFixed(1);

        // Default titles if sorted data is not yet available
        let titleA = componentsInput[0].name || 'Comp 1';
        let titleB = componentsInput[1].name || 'Comp 2';
        let titleC = componentsInput[2].name || 'Comp 3';

        if (plotSortedComponents.length === 3) {
            titleA = `${plotSortedComponents[0].name} (Lowest NBP)`;
            titleB = `${plotSortedComponents[1].name} (Mid NBP)`;
            titleC = `${plotSortedComponents[2].name} (Highest NBP)`;
        }

        const baseLayout: Partial<Layout> = {
            title: 'Ternary Residue Curve Map',
            ternary: {
                sum: 1,
                aaxis: { title: titleA, min: 0, max: 1, tickformat: '.1f' },
                baxis: { title: titleB, min: 0, max: 1, tickformat: '.1f' },
                caxis: { title: titleC, min: 0, max: 1, tickformat: '.1f' },
            },
            margin: { l: 90, r: 90, b: 90, t: 90, pad: 5 }, // Increased margins for longer titles
            autosize: true,
        };

        if (residueCurves.length === 0 || plotSortedComponents.length !== 3) {
            setPlotlyData([]);
            setPlotlyLayout({
                ...baseLayout,
                title: 'Ternary Residue Curve Map (No data or NBP info pending)',
                annotations: [{
                    text: "Enter component data and click 'Generate Residue Map'.",
                    showarrow: false,
                    xref: 'paper',
                    yref: 'paper',
                    x: 0.5,
                    y: 0.5,
                }]
            });
            return;
        }

        // Map residue curve data (x1,x2,x3 for original components 0,1,2)
        // to Plotly axes (a,b,c for lowest, mid, highest NBP components)
        const traces: Data[] = residueCurves.map((curve, index) => {
            const getFraction = (p: ResidueCurvePoint, originalIndex: number): number => {
                if (originalIndex === 0) return p.x1;
                if (originalIndex === 1) return p.x2;
                if (originalIndex === 2) return p.x3;
                return 0; // Should not happen
            };

            return {
                type: 'scatterternary',
                mode: 'lines',
                a: curve.map(p => getFraction(p, plotSortedComponents[0].originalIndex)), // Lowest NBP
                b: curve.map(p => getFraction(p, plotSortedComponents[1].originalIndex)), // Mid NBP
                c: curve.map(p => getFraction(p, plotSortedComponents[2].originalIndex)), // Highest NBP
                name: `Curve ${index + 1}`,
                line: { width: 1.5 },
                hovertemplate: 
                    `${plotSortedComponents[0].name}: %{a:.3f}<br>` +
                    `${plotSortedComponents[1].name}: %{b:.3f}<br>` +
                    `${plotSortedComponents[2].name}: %{c:.3f}<extra></extra>`,
            };
        });

        setPlotlyData(traces);
        setPlotlyLayout({
            ...baseLayout,
            title: `Ternary Residue Curve Map @ ${pressureDisplay} kPa`,
        });

    }, [residueCurves, componentsInput, systemPressure, plotSortedComponents]);

    return (
        <div className="container mx-auto p-4">
            <Card className="mb-6">
                <CardHeader><CardTitle>Ternary Residue Curve Map Inputs</CardTitle></CardHeader>
                <CardContent className="space-y-6">
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-4"> {/* Horizontal layout for components */}
                        {componentsInput.map((comp, index) => (
                            <Card key={index} className="flex-1"> {/* Use Card for each component */}
                                <CardHeader className="pb-2 pt-4">
                                    <CardTitle className="text-base">Component {index + 1}</CardTitle>
                                </CardHeader>
                                <CardContent className="space-y-3">
                                    <div>
                                        <Label htmlFor={`name-${index}`}>Name</Label>
                                        <Input 
                                            id={`name-${index}`} 
                                            value={comp.name} 
                                            onChange={(e) => handleComponentInputChange(index, 'name', e.target.value)} 
                                            placeholder={`e.g., ${index === 0 ? 'Acetone' : index === 1 ? 'Methanol' : 'Water'}`} 
                                        />
                                    </div>
                                </CardContent>
                            </Card>
                        ))}
                    </div>
                    <div>
                        <Label htmlFor="systemPressure">System Pressure (kPa)</Label>
                        <Input 
                            id="systemPressure" 
                            value={systemPressure} 
                            onChange={(e) => setSystemPressure(e.target.value)} 
                            placeholder="e.g., 101.325" 
                            type="number"
                        />
                    </div>
                    <Button onClick={handleGenerateMap} disabled={isLoading || !supabase}>
                        {isLoading ? 'Generating...' : 'Generate Residue Map'}
                    </Button>
                    {error && <p className="text-red-500">{error}</p>}
                    {!supabase && <p className="text-orange-500">Supabase client not configured. Parameters cannot be fetched.</p>}
                </CardContent>
            </Card>

            <Card>
                <CardHeader><CardTitle>Residue Curve Map</CardTitle></CardHeader>
                <CardContent>
                    {(residueCurves.length > 0 || isLoading) ? (
                        <Plot
                            data={plotlyData}
                            layout={plotlyLayout}
                            style={{ width: '100%', height: '600px' }}
                            useResizeHandler={true}
                            config={{ responsive: true }}
                        />
                    ) : (
                        <div className="text-center p-10 text-muted-foreground min-h-[600px] flex items-center justify-center">
                            <div>
                                <p>Enter component names and system pressure, then click "Generate Residue Map" to see the plot.</p>
                                <p>CAS numbers, Antoine parameters, and Normal Boiling Points will be fetched automatically.</p>
                                <p className="mt-2">
                                    Plotly Ternary Axes: a (Lowest NBP), b (Intermediate NBP), c (Highest NBP).
                                </p>
                            </div>
                        </div>
                    )}
                </CardContent>
            </Card>
        </div>
    );
}

