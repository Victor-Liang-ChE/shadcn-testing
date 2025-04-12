"use client";
import React, { useState, useEffect, useCallback, KeyboardEvent, useRef } from 'react';

// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption, SeriesOption } from 'echarts';
import { LineChart } from 'echarts/charts';
import {
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent // For save image feature
} from 'echarts/components';
import { CanvasRenderer } from 'echarts/renderers';

// Register necessary ECharts components
echarts.use([
  TitleComponent,
  TooltipComponent,
  GridComponent,
  LegendComponent,
  ToolboxComponent,
  LineChart,
  CanvasRenderer
]);

// Import Shadcn UI components
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";

// Define a color palette for species
const colorPalette = [
  '#5470C6', '#91CC75', '#FAC858', '#EE6666', '#73C0DE',
  '#3BA272', '#FC8452', '#9A60B4', '#EA7CCC'
];

// --- Helper Functions ... ---
// Function to detect unique species in the order they appear
function detectUniqueSpeciesOrdered(reactions: string[]): string[] {
  const orderedSpecies: string[] = [];
  const uniqueSpecies = new Set<string>();

  for (const reaction of reactions) {
    if (!reaction) continue;

    // Split the reaction into reactants and products
    const parts = reaction.split('->');
    if (parts.length !== 2) continue; // Skip invalid format
    const [reactants, products] = parts;

    // Extract species from reactants and products
    const reactantSpecies = reactants.split('+').map(species =>
      species.trim().replace(/^\d+\*?/, '').trim()).filter(s => s); // Filter empty strings

    const productSpecies = products.split('+').map(species =>
      species.trim().replace(/^\d+\*?/, '').trim()).filter(s => s); // Filter empty strings

    // Add species to the ordered list if they are not already in the set
    for (const species of [...reactantSpecies, ...productSpecies]) {
      if (species && !uniqueSpecies.has(species)) {
        uniqueSpecies.add(species);
        orderedSpecies.push(species);
      }
    }
  }

  return orderedSpecies;
}

// Function to format species names with subscripted numbers (returns string for ECharts)
function formatSpeciesEcharts(species: string): string {
  return species.replace(/(\d+)/g, '<sub>$1</sub>');
}

// Adaptive RK45 solver (solveODERK45 function remains the same as in the original file)
function solveODERK45(
  dydt: (t: number, y: number[]) => number[],
  y0: number[],
  tSpan: [number, number],
  tol: number = 1e-6,
  h_initial: number = 0.1,
  maxSteps: number = 10000
): { t: number[]; y: number[][] } {
  const [t0, tf] = tSpan;
  let t = t0;
  let y = y0.slice(); // copy initial state
  const t_arr = [t];
  const y_arr = [y.slice()]; // y_arr stores states at each time step

  // Adaptive step parameters
  let h = h_initial;
  const safety = 0.84;
  const minFactor = 0.2;
  const maxFactor = 5.0;

  // Dormand–Prince coefficients (RK45)
  // c values
  const c2 = 1/5, c3 = 3/10, c4 = 4/5, c5 = 8/9, c6 = 1, c7 = 1;
  // a coefficients (lower-triangular matrix)
  const a21 = 1/5;
  const a31 = 3/40,  a32 = 9/40;
  const a41 = 44/45, a42 = -56/15, a43 = 32/9;
  const a51 = 19372/6561, a52 = -25360/2187, a53 = 64448/6561, a54 = -212/729;
  const a61 = 9017/3168,  a62 = -355/33,    a63 = 46732/5247, a64 = 49/176,  a65 = -5103/18656;
  const a71 = 35/384,   a72 = 0,          a73 = 500/1113,  a74 = 125/192, a75 = -2187/6784, a76 = 11/84;
  // b coefficients for the 5th‑order solution
  const b1 = 35/384, b2 = 0, b3 = 500/1113, b4 = 125/192, b5 = -2187/6784, b6 = 11/84, b7 = 0;
  // b* coefficients for the 4th‑order solution
  const b1s = 5179/57600, b2s = 0, b3s = 7571/16695, b4s = 393/640, b5s = -92097/339200, b6s = 187/2100, b7s = 1/40;

  let stepCount = 0;
  while (t < tf && stepCount < maxSteps) {
    // Adjust final step if necessary
    if (t + h > tf) {
      h = tf - t;
    }
    if (h <= 1e-15) break; // Prevent infinite loop if h becomes too small

    // Compute the RK stages:
    const k1 = dydt(t, y);
    const y2 = y.map((yi, i) => yi + h * a21 * k1[i]);
    const k2 = dydt(t + c2 * h, y2);
    const y3 = y.map((yi, i) => yi + h * (a31 * k1[i] + a32 * k2[i]));
    const k3 = dydt(t + c3 * h, y3);
    const y4 = y.map((yi, i) => yi + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]));
    const k4 = dydt(t + c4 * h, y4);
    const y5 = y.map((yi, i) => yi + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]));
    const k5 = dydt(t + c5 * h, y5);
    const y6 = y.map((yi, i) => yi + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]));
    const k6 = dydt(t + c6 * h, y6);
    const y7 = y.map((yi, i) => yi + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]));
    const k7 = dydt(t + c7 * h, y7);

    // Compute 5th-order solution
    const y5th = y.map((yi, i) =>
      yi + h * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i])
    );
    // Compute 4th-order solution
    const y4th = y.map((yi, i) =>
      yi + h * (b1s * k1[i] + b2s * k2[i] + b3s * k3[i] + b4s * k4[i] + b5s * k5[i] + b6s * k6[i] + b7s * k7[i])
    );

    // Estimate the error using the maximum absolute difference relative to tolerance
    let error_ratio = 0;
    for (let i = 0; i < y.length; i++) {
        const scale = tol + Math.abs(y[i]) * tol; // Relative tolerance scaling
        error_ratio = Math.max(error_ratio, Math.abs(y5th[i] - y4th[i]) / scale);
    }
    error_ratio = Math.max(error_ratio, 1e-10); // Avoid division by zero

    // If error is within tolerance (error_ratio <= 1), accept the step
    if (error_ratio <= 1) {
      t += h;
      y = y5th.slice();
      t_arr.push(t);
      y_arr.push(y.slice());
    }

    // Adjust the step size for the next iteration.
    const factor = safety * Math.pow(error_ratio, -1 / 5);
    const factorClamped = Math.min(maxFactor, Math.max(minFactor, factor));
    h = h * factorClamped;
    stepCount++;
  }

  if (stepCount >= maxSteps) {
    console.warn(`RK45 solver reached maximum steps (${maxSteps})`);
  }

  // Transpose y_arr to match the previous format (rows = species, cols = time points)
  const numSpecies = y0.length;
  const numTimePoints = t_arr.length;
  const y_transposed: number[][] = Array.from({ length: numSpecies }, () => Array(numTimePoints).fill(0));
  for (let i = 0; i < numTimePoints; i++) {
    for (let j = 0; j < numSpecies; j++) {
      y_transposed[j][i] = y_arr[i][j];
    }
  }

  return { t: t_arr, y: y_transposed };
}


// --- KineticsPage Component ---
export default function KineticsPage() {
  // Default reaction setup
  const defaultReaction = '2H2 + O2 -> 2H2O';
  const defaultReactants = '2H2 + O2';
  const defaultProducts = '2H2O';
  const defaultRateConstant = 1;
  const defaultSpecies = detectUniqueSpeciesOrdered([defaultReaction]);
  const defaultConcentrations: Record<string, number | ''> = {};
  defaultSpecies.forEach(sp => { defaultConcentrations[sp] = 1; });
  const initialSpeciesColors: Record<string, string> = {};
   defaultSpecies.forEach((sp, index) => {
     initialSpeciesColors[sp] = colorPalette[index % colorPalette.length];
   });


  const [reactionInputs, setReactionInputs] = useState<{
    reaction: string;
    reactants: string;
    products: string;
    rateConstant: number | ''
  }[]>([{
    reaction: defaultReaction,
    reactants: defaultReactants,
    products: defaultProducts,
    rateConstant: defaultRateConstant
  }]);
  const [confirmedReactions, setConfirmedReactions] = useState<boolean>(true); // Start as confirmed
  const [species, setSpecies] = useState<string[]>(defaultSpecies);
  const [concentrations, setConcentrations] = useState<Record<string, number | ''>>(defaultConcentrations);
  const [speciesColors, setSpeciesColors] = useState<Record<string, string>>(initialSpeciesColors); // State for species colors
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const [showGraph, setShowGraph] = useState<boolean>(true);
  const [loading, setLoading] = useState<boolean>(false);
  const echartsRef = useRef<ReactECharts | null>(null); // Ref for ECharts instance

  // Function to generate ECharts options (modified reactionGraphing)
  const generateEchartsOptions = useCallback((
    reactions: string[],
    ks: number[],
    C0: Record<string, number>,
    currentSpeciesColors: Record<string, string> // Pass colors
  ): EChartsOption => {
    setLoading(true);
    try {
      if (ks.length !== reactions.length) {
        throw new Error("The number of rate constants does not match the number of reactions.");
      }
      const orderedSpecies = detectUniqueSpeciesOrdered(reactions);
      const missingSpecies = orderedSpecies.filter(species => !C0.hasOwnProperty(species));
      if (missingSpecies.length > 0) {
        throw new Error(`The following species are missing in C0: ${missingSpecies.join(', ')}`);
      }

      function odes(t: number, y: number[]): number[] {
        const dydt = new Array(orderedSpecies.length).fill(0);
        const currentConcentrations: Record<string, number> = {};
        orderedSpecies.forEach((species, i) => {
          currentConcentrations[species] = Math.max(0, y[i]); // Ensure non-negative concentrations
        });

        reactions.forEach((reaction, i) => {
          if (!reaction) return;
          const parts = reaction.split('->');
          if (parts.length !== 2) return;
          const [reactantsStr, productsStr] = parts;

          const reactantSpecies: [string, number][] = [];
          reactantsStr.split('+').forEach(speciesPart => {
            const match = speciesPart.trim().match(/^(?:(\d+)\*?)?(.+)$/);
            if (!match) return;
            const [, coeffStr, sp] = match;
            const coefficient = coeffStr ? parseInt(coeffStr) : 1;
            if (sp && orderedSpecies.includes(sp.trim())) {
              reactantSpecies.push([sp.trim(), coefficient]);
            }
          });

          const productSpecies: [string, number][] = [];
           productsStr.split('+').forEach(speciesPart => {
            const match = speciesPart.trim().match(/^(?:(\d+)\*?)?(.+)$/);
            if (!match) return;
            const [, coeffStr, sp] = match;
            const coefficient = coeffStr ? parseInt(coeffStr) : 1;
            if (sp && orderedSpecies.includes(sp.trim())) {
              productSpecies.push([sp.trim(), coefficient]);
            }
          });

          // Calculate the reaction rate
          let rate = ks[i];
          for (const [sp, coeff] of reactantSpecies) {
             rate *= Math.pow(currentConcentrations[sp], coeff);
          }

          // Update differential forms
          reactantSpecies.forEach(([sp, coeff]) => {
            const index = orderedSpecies.indexOf(sp);
            if (index !== -1) dydt[index] -= rate * coeff;
          });
          productSpecies.forEach(([sp, coeff]) => {
             const index = orderedSpecies.indexOf(sp);
             if (index !== -1) dydt[index] += rate * coeff;
          });
        });
        return dydt;
      }

      const y0 = orderedSpecies.map(species => C0[species]);
      const tSpan: [number, number] = [0, 10];
      const initial_h = 0.01;
      const tolerance = 1e-5;

      const solution = solveODERK45(odes, y0, tSpan, tolerance, initial_h);

      // Restore steady state time calculation for axis limit
      const maxConcentration = Math.max(1e-9, ...solution.y.flat());
      const relativeTolerance = maxConcentration * 1e-4;
      let steadyStateTime = tSpan[1];
      for (let i = 1; i < solution.t.length; i++) {
        const concentrationDiff = orderedSpecies.map((_, j) =>
          Math.abs(solution.y[j][i] - solution.y[j][i - 1]));
        if (concentrationDiff.every(diff => diff < relativeTolerance)) {
          steadyStateTime = solution.t[i];
          break;
        }
      }
      steadyStateTime = steadyStateTime < tSpan[1] ? steadyStateTime * 1.1 : steadyStateTime;
      steadyStateTime = Math.min(steadyStateTime, tSpan[1]);

      // Restore dynamic yRange calculation for axis limit
      const allYValues = solution.y.flat().filter(y => !isNaN(y) && isFinite(y));
      const minY = Math.min(0, ...allYValues);
      const maxY = Math.max(1e-9, ...allYValues);
      const yRange = [minY, maxY * 1.1]; // Dynamic range with buffer

      const seriesData: SeriesOption[] = orderedSpecies.map((species, i) => ({
        name: formatSpeciesEcharts(species), // Should return 'H<sub>2</sub>O' etc.
        type: 'line',
        smooth: true,
        symbol: 'none',
        data: solution.t.map((time, timeIndex) => [time, solution.y[i][timeIndex]]),
        itemStyle: {
          color: currentSpeciesColors[species] || colorPalette[i % colorPalette.length]
        },
        emphasis: {
            focus: 'series'
        },
        // animation: false // Keep animation disabled per series if needed
      }));

      return {
        backgroundColor: '#08306b',
        animation: false, // Keep animations disabled globally
        title: {
            text: 'Concentration Profiles',
            left: 'center',
            textStyle: {
                color: '#fff',
                fontSize: 18,
                fontFamily: 'sans-serif'
            }
        },
        tooltip: {
          trigger: 'axis',
          formatter: (params: any) => {
             let tooltipText = `Time: ${params[0].axisValueLabel}<br/>`;
             params.forEach((param: any) => {
                 tooltipText += `${param.marker}${param.seriesName}: ${param.value[1].toPrecision(4)}<br/>`;
             });
             return tooltipText;
          }
        },
        legend: {
          // Ensure data uses objects with 'name' containing rich text
          data: orderedSpecies.map(sp => ({ name: formatSpeciesEcharts(sp) })),
          textStyle: {
             color: '#fff',
             fontSize: 12,
             fontFamily: 'sans-serif',
             // Ensure rich text configuration is correct
             rich: {
                sub: {
                    fontSize: 9, // Smaller font size for subscript
                    verticalAlign: 'bottom' // Align subscript to bottom
                }
             }
          },
          top: 'bottom',
          type: 'scroll',
        },
        grid: {
          left: '5%',
          right: '5%',
          bottom: '12%', // Match top padding
          top: '12%',    // Match bottom padding
          containLabel: true
        },
        toolbox: {
          feature: {
            saveAsImage: {
                name: 'kinetics-plot',
                backgroundColor: '#08306b' // Match background
            }
          },
          iconStyle: {
             borderColor: '#fff' // White icons
          },
          orient: 'vertical', // Keep vertical orientation
          right: 10, // Position from right
          bottom: 40 // Adjust position if needed
        },
        xAxis: {
          type: 'value',
          name: 'Time',
          nameLocation: 'middle',
          nameGap: 30,
          min: 0,
          max: steadyStateTime, // Use dynamic max
          axisLabel: {
              color: '#fff',
              fontSize: 14,
              fontFamily: 'sans-serif',
              // Apply consistent formatting to x-axis labels
              formatter: (value: number) => {
                  if (Math.abs(value) < 1e-6) return '0';
                  // Use toFixed for time for better readability if appropriate, or toPrecision
                  return value.toFixed(2); // Example: 2 decimal places for time
              }
          },
          nameTextStyle: {
              color: '#fff',
              fontSize: 15,
              fontFamily: 'sans-serif'
          },
          axisLine: { lineStyle: { color: '#fff' } },
          axisTick: { lineStyle: { color: '#fff' } },
          splitLine: { show: false }
        },
        yAxis: {
          type: 'value',
          name: 'Concentration',
          nameLocation: 'middle',
          nameGap: 50,
          min: yRange[0], // Use dynamic min
          max: yRange[1], // Use dynamic max
          axisLabel: {
              color: '#fff',
              fontSize: 14,
              fontFamily: 'sans-serif',
              formatter: (value: number) => { // Consistent formatter for y-axis
                  if (Math.abs(value) < 1e-6) return '0';
                  // Apply precision consistently, including max value
                  return value.toPrecision(3); // Example: 3 significant figures
              }
          },
          nameTextStyle: {
              color: '#fff',
              fontSize: 15,
              fontFamily: 'sans-serif'
          },
          axisLine: { lineStyle: { color: '#fff' } },
          axisTick: { lineStyle: { color: '#fff' } },
          splitLine: { show: false }
        },
        series: seriesData,
      };
    } catch (error) {
      console.error("Error generating ECharts options:", error);
      // Return an empty options object or handle error state appropriately
      return {};
    } finally {
      setLoading(false);
    }
  }, []); // Dependencies remain empty

  // Function to trigger graph update - Modified to accept optional arguments
  const submitForm = useCallback((
    reactionsOverride?: string[],
    ksOverride?: number[],
    concentrationsOverride?: Record<string, number>
  ) => {
    const validReactions = reactionsOverride ?? reactionInputs
      .filter(input => input.reactants.trim() && input.products.trim())
      .map(input => `${input.reactants.trim()} -> ${input.products.trim()}`);

    const validRateConstants = ksOverride ?? reactionInputs
      .filter(input => input.reactants.trim() && input.products.trim())
      .map(input => (typeof input.rateConstant === 'number' ? input.rateConstant : 1)); // Default k=1 if empty

    const validConcentrations: Record<string, number> = {};
    let allConcValid = true;

    if (concentrationsOverride) {
        // If overrides are provided, use them directly (assuming they are valid numbers)
        Object.assign(validConcentrations, concentrationsOverride);
    } else {
        // Otherwise, derive from state and validate
        for (const sp of species) {
            const concValue = concentrations[sp];
            if (concValue === '' || isNaN(Number(concValue))) {
                allConcValid = false;
                break; // Exit loop early if any concentration is invalid
            }
            validConcentrations[sp] = Number(concValue);
        }
    }

    // Check if all species needed for the *current* valid reactions have concentrations
    const currentSpeciesForReactions = detectUniqueSpeciesOrdered(validReactions);
    const missingSpeciesInConcentrations = currentSpeciesForReactions.filter(sp => !(sp in validConcentrations));

    if (!allConcValid || missingSpeciesInConcentrations.length > 0) {
        console.warn("Some initial concentrations are invalid or missing for the current reactions.", { missingSpeciesInConcentrations });
        // Maybe show a subtle warning instead of alert?
        return;
    }

    if (validReactions.length === 0) {
        console.warn("No valid reactions entered.");
        setEchartsOptions({}); // Clear graph
        setShowGraph(false);
        return;
    }

    // Pass the current speciesColors to the generation function
    const options = generateEchartsOptions(validReactions, validRateConstants, validConcentrations, speciesColors);
    setEchartsOptions(options);
    setShowGraph(true);

  }, [reactionInputs, concentrations, species, generateEchartsOptions, speciesColors]); // Keep dependencies

  // Generate initial plot on mount
  useEffect(() => {
    submitForm(); // Call submitForm directly to generate initial plot
  }, []); // Run only once on mount


  // --- Event Handlers ---
  const addReactionInput = () => {
    setReactionInputs([...reactionInputs, {
      reaction: '',
      reactants: '',
      products: '',
      rateConstant: 1
    }]);
    setConfirmedReactions(false); // Require confirmation after adding
  };

  const removeReactionInput = () => {
    if (reactionInputs.length > 1) {
      setReactionInputs(reactionInputs.slice(0, -1));
      setConfirmedReactions(false); // Require confirmation after removing
    }
  };

  const handleReactionChange = (index: number, field: 'reactants' | 'products' | 'rateConstant', value: string | number) => {
    const newInputs = [...reactionInputs];
    const currentInput = { ...newInputs[index] };

    if (field === 'reactants') {
      currentInput.reactants = value as string;
    } else if (field === 'products') {
      currentInput.products = value as string;
    } else {
      currentInput.rateConstant = value === '' ? '' : Number(value);
    }

    // Reconstruct reaction string
    currentInput.reaction = `${currentInput.reactants} -> ${currentInput.products}`;
    newInputs[index] = currentInput;

    setReactionInputs(newInputs);
    setConfirmedReactions(false); // Require confirmation after change
  };

  const confirmReactions = useCallback(() => {
    const validReactions = reactionInputs
      .filter(input => input.reactants.trim() && input.products.trim())
      .map(input => `${input.reactants.trim()} -> ${input.products.trim()}`);

    const validRateConstants = reactionInputs
      .filter(input => input.reactants.trim() && input.products.trim())
      .map(input => (typeof input.rateConstant === 'number' ? input.rateConstant : 1));

    if (validReactions.length > 0) {
      const detectedSpecies = detectUniqueSpeciesOrdered(validReactions);
      setSpecies(detectedSpecies);

      // Update species colors
      const newSpeciesColors: Record<string, string> = {};
      detectedSpecies.forEach((sp, index) => {
        newSpeciesColors[sp] = speciesColors[sp] || colorPalette[index % colorPalette.length]; // Preserve existing or assign new
      });
      setSpeciesColors(newSpeciesColors);

      // Prepare the next concentration state for both UI and immediate plotting
      const currentConcentrationsState = { ...concentrations };
      const nextConcentrationsState: Record<string, number | ''> = {};
      const concentrationsForPlot: Record<string, number> = {}; // For immediate plot

      detectedSpecies.forEach(sp => {
        const currentValue = currentConcentrationsState.hasOwnProperty(sp)
                                     ? currentConcentrationsState[sp]
                                     : 1; // Default to 1 if new
        nextConcentrationsState[sp] = currentValue; // Update state for UI
        // Ensure a valid number for plotting, defaulting to 1 if empty/invalid
        concentrationsForPlot[sp] = (typeof currentValue === 'number' && !isNaN(currentValue)) ? currentValue : 1;
      });

      setConcentrations(nextConcentrationsState); // Update UI state
      setConfirmedReactions(true);

      // Trigger graph update immediately with the calculated data
      // Pass overrides to submitForm
      submitForm(validReactions, validRateConstants, concentrationsForPlot);
      // Removed setTimeout(submitForm, 0);

    } else {
      setSpecies([]);
      setConcentrations({});
      setSpeciesColors({});
      setEchartsOptions({});
      setShowGraph(false);
      setConfirmedReactions(false);
    }
  }, [reactionInputs, concentrations, submitForm, speciesColors]); // submitForm is now a dependency

  const handleConcentrationChange = (species: string, value: number | '') => {
    setConcentrations(prev => ({ ...prev, [species]: value === '' ? '' : Number(value) }));
    // Trigger graph update immediately after slider change if reactions are confirmed
    if (confirmedReactions) {
        setTimeout(submitForm, 0);
    }
  };

  // Trigger confirmation on Enter key in reaction inputs
  const handleReactionKeyDown = (e: KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      confirmReactions();
    }
  };

  // Trigger submit on Enter key in concentration inputs
  const handleConcentrationKeyDown = (e: KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter') {
      e.preventDefault();
      if (confirmedReactions) { // Only submit if reactions are confirmed
          submitForm();
      }
    }
  };

  // Format slider value display
  const formatSliderValue = (value: number | ''): string => {
    return typeof value === 'number' ? value.toFixed(2) : '0.00';
  };

  // Handle slider changes and trigger graph update
  const handleSliderChange = (type: 'rateConstant' | 'concentration', index: number | string, value: number) => {
    let shouldSubmit = false;
    if (type === 'rateConstant') {
      // Update reaction input state
      const newInputs = [...reactionInputs];
      newInputs[index as number] = { ...newInputs[index as number], rateConstant: value };
      // Reconstruct reaction string (optional, but good practice if reaction string is used elsewhere)
      const currentInput = newInputs[index as number];
      currentInput.reaction = `${currentInput.reactants} -> ${currentInput.products}`;
      setReactionInputs(newInputs);
      // Mark for submission only if reactions were already confirmed
      shouldSubmit = confirmedReactions;
    } else { // concentration
      // Update concentration state directly using the callback form for safety
      setConcentrations(prev => ({ ...prev, [index as string]: value }));
      // Mark for submission only if reactions were already confirmed
      shouldSubmit = confirmedReactions;
    }

    // Trigger submit immediately after state update if needed
    // Use setTimeout to allow state update to likely flush before submitForm reads it
    if (shouldSubmit) {
        setTimeout(submitForm, 0);
    }
  };


  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">

          {/* Column 1: Controls */}
          <div className="lg:col-span-1 space-y-6">
            {/* Reactions Card */}
            <Card>
              <CardContent className="space-y-6 pt-4">
                {reactionInputs.map((input, index) => (
                  <div key={index} className="space-y-6 border-b pb-6 last:border-b-0 last:pb-0">
                    <div className="flex items-center gap-2">
                      <Input
                        value={input.reactants}
                        onChange={(e) => handleReactionChange(index, 'reactants', e.target.value)}
                        onKeyDown={handleReactionKeyDown}
                        placeholder="e.g., 2H2 + O2"
                        className="flex-1"
                      />
                      <span className="font-bold text-lg">→</span>
                      <Input
                        value={input.products}
                        onChange={(e) => handleReactionChange(index, 'products', e.target.value)}
                        onKeyDown={handleReactionKeyDown}
                        placeholder="e.g., 2H2O"
                        className="flex-1"
                      />
                    </div>
                    <div className="space-y-2">
                      <Label htmlFor={`rate-${index}`}>
                        {/* Use dangerouslySetInnerHTML for precise subscript rendering */}
                        <span dangerouslySetInnerHTML={{ __html: `Rate Constant (k<sub>${index + 1}</sub>): ${formatSliderValue(input.rateConstant)}` }} />
                      </Label>
                      <Slider
                        id={`rate-${index}`}
                        min={0}
                        max={10}
                        step={0.01}
                        value={[Number(input.rateConstant || 0)]}
                        // Pass index and value directly
                        onValueChange={(val) => handleSliderChange('rateConstant', index, val[0])}
                        style={{ '--primary': 'hsl(var(--primary))' } as React.CSSProperties}
                      />
                    </div>
                  </div>
                ))}
                <div className="flex gap-2">
                  <Button variant="outline" onClick={addReactionInput} className="flex-1">
                    Add Reaction
                  </Button>
                  <Button variant="outline" onClick={removeReactionInput} disabled={reactionInputs.length <= 1} className="flex-1">
                    Remove Last Reaction
                  </Button>
                </div>
                 <Button onClick={confirmReactions} className="w-full" disabled={confirmedReactions}>
                    Confirm Reactions & Update Species
                 </Button>
              </CardContent>
            </Card>

            {/* Concentrations Card */}
            {species.length > 0 && (
              <Card>
                <CardHeader>
                  <CardTitle>Initial Concentrations</CardTitle>
                </CardHeader>
                <CardContent className="space-y-6">
                  {species.map((sp) => (
                    <div key={sp} className="space-y-2">
                      <Label htmlFor={`conc-${sp}`}>
                        <span dangerouslySetInnerHTML={{ __html: `[${formatSpeciesEcharts(sp)}]₀: ${formatSliderValue(concentrations[sp])}` }} />
                      </Label>
                      <Slider
                        id={`conc-${sp}`}
                        min={0}
                        max={10}
                        step={0.01}
                        value={[Number(concentrations[sp] || 0)]}
                        onValueChange={(val) => handleSliderChange('concentration', sp, val[0])}
                        onKeyDown={handleConcentrationKeyDown}
                        style={{ '--primary': speciesColors[sp] || 'hsl(var(--primary))' } as React.CSSProperties}
                      />
                    </div>
                  ))}
                </CardContent>
              </Card>
            )}
          </div>

          {/* Column 2: Plot */}
          <div className="lg:col-span-2">
            <Card>
              <CardContent className="pt-6">
                <div className="relative h-[500px] md:h-[600px] rounded-md overflow-hidden" style={{ backgroundColor: '#08306b' }}>
                  {loading && (
                     <div className="absolute inset-0 flex items-center justify-center text-white z-10 bg-opacity-50 backdrop-blur-sm">
                       Generating plot...
                     </div>
                  )}
                  {showGraph && Object.keys(echartsOptions).length > 0 ? (
                    <ReactECharts
                      ref={echartsRef}
                      echarts={echarts}
                      option={echartsOptions}
                      style={{ height: '100%', width: '100%' }}
                      notMerge={true} // Important for applying options like animation: false
                      lazyUpdate={false}
                    />
                  ) : (
                    <div className="absolute inset-0 flex items-center justify-center text-gray-400"> {/* Adjusted text color for dark bg */}
                      {species.length > 0 ? "Adjust parameters to see the concentration profiles." : "Enter reactions and click 'Confirm Reactions'."}
                    </div>
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
