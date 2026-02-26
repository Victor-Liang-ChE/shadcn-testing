"use client";
import { useState, useEffect, useCallback, KeyboardEvent, useRef } from 'react';
import { useTheme } from "next-themes";

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
import { Tooltip as _Tooltip, TooltipProvider } from "@/components/ui/tooltip";

// Define a color palette for species
const colorPalette = [
  '#5470C6', '#91CC75', '#FAC858', '#EE6666', '#73C0DE',
  '#3BA272', '#FC8452', '#9A60B4', '#EA7CCC'
];

// --- Helper Functions --- (detectUniqueSpeciesOrdered, solveODERK45 remain the same)
// Function to detect unique species in the order they appear
function detectUniqueSpeciesOrdered(reactions: string[]): string[] {
  const orderedSpecies: string[] = [];
  const uniqueSpecies = new Set<string>();

  for (const reaction of reactions) {
    if (!reaction) continue;
    const parts = reaction.split('->');
    if (parts.length !== 2) continue;
    const [reactants, products] = parts;
    const reactantSpecies = reactants.split('+').map(species => species.trim().replace(/^\d+\*?/, '').trim()).filter(s => s);
    const productSpecies = products.split('+').map(species => species.trim().replace(/^\d+\*?/, '').trim()).filter(s => s);
    for (const species of [...reactantSpecies, ...productSpecies]) {
      if (species && !uniqueSpecies.has(species)) {
        uniqueSpecies.add(species);
        orderedSpecies.push(species);
      }
    }
  }
  return orderedSpecies;
}

// Adaptive RK45 solver
function solveODERK45(
  dydt: (t: number, y: number[]) => number[],
  y0: number[],
  tSpan: [number, number],
  tol: number = 1e-6,
  h_initial: number = 0.1,
  maxSteps: number = 10000
): { t: number[]; y: number[][] } {
  const [t0, tf] = tSpan;
  let t = t0; let y = y0.slice();
  const t_arr = [t]; const y_arr = [y.slice()];
  let h = h_initial; const safety = 0.84; const minFactor = 0.2; const maxFactor = 5.0;
  const c2=1/5, c3=3/10, c4=4/5, c5=8/9, c6=1, c7=1;
  const a21=1/5; const a31=3/40, a32=9/40; const a41=44/45, a42=-56/15, a43=32/9;
  const a51=19372/6561, a52=-25360/2187, a53=64448/6561, a54=-212/729;
  const a61=9017/3168, a62=-355/33, a63=46732/5247, a64=49/176, a65=-5103/18656;
  const a71=35/384, a72=0, a73=500/1113, a74=125/192, a75=-2187/6784, a76=11/84;
  const b1=35/384, b2=0, b3=500/1113, b4=125/192, b5=-2187/6784, b6=11/84, b7=0;
  const b1s=5179/57600, b2s=0, b3s=7571/16695, b4s=393/640, b5s=-92097/339200, b6s=187/2100, b7s=1/40;
  let stepCount = 0;
  while (t < tf && stepCount < maxSteps) {
    if (t + h > tf) h = tf - t;
    if (h <= 1e-15) break;
    const k1 = dydt(t, y); const y2 = y.map((yi, i) => yi + h * a21 * k1[i]);
    const k2 = dydt(t + c2 * h, y2); const y3 = y.map((yi, i) => yi + h * (a31 * k1[i] + a32 * k2[i]));
    const k3 = dydt(t + c3 * h, y3); const y4 = y.map((yi, i) => yi + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]));
    const k4 = dydt(t + c4 * h, y4); const y5 = y.map((yi, i) => yi + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]));
    const k5 = dydt(t + c5 * h, y5); const y6 = y.map((yi, i) => yi + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]));
    const k6 = dydt(t + c6 * h, y6); const y7 = y.map((yi, i) => yi + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]));
    const k7 = dydt(t + c7 * h, y7);
    const y5th = y.map((yi, i) => yi + h * (b1 * k1[i] + b2 * k2[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i]));
    const y4th = y.map((yi, i) => yi + h * (b1s * k1[i] + b2s * k2[i] + b3s * k3[i] + b4s * k4[i] + b5s * k5[i] + b6s * k6[i] + b7s * k7[i]));
    let error_ratio = 0;
    for (let i = 0; i < y.length; i++) { const scale = tol + Math.abs(y[i]) * tol; error_ratio = Math.max(error_ratio, Math.abs(y5th[i] - y4th[i]) / scale); }
    error_ratio = Math.max(error_ratio, 1e-10);
    if (error_ratio <= 1) { t += h; y = y5th.slice(); t_arr.push(t); y_arr.push(y.slice()); }
    const factor = safety * Math.pow(error_ratio, -1 / 5); const factorClamped = Math.min(maxFactor, Math.max(minFactor, factor));
    h = h * factorClamped; stepCount++;
  }
  if (stepCount >= maxSteps) console.warn(`RK45 solver reached maximum steps (${maxSteps})`);
  const numSpecies = y0.length; const numTimePoints = t_arr.length;
  const y_transposed: number[][] = Array.from({ length: numSpecies }, () => Array(numTimePoints).fill(0));
  for (let i = 0; i < numTimePoints; i++) { for (let j = 0; j < numSpecies; j++) { y_transposed[j][i] = y_arr[i][j]; } }
  return { t: t_arr, y: y_transposed };
}


// --- KineticsPage Component ---
export default function KineticsPage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  
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
  
  // Add state for max slider values
  const [maxRateConstantSlider, setMaxRateConstantSlider] = useState<string>('10');
  const [maxConcentrationSlider, setMaxConcentrationSlider] = useState<string>('10');
  const [confirmedReactions, setConfirmedReactions] = useState<boolean>(true);
  const [species, setSpecies] = useState<string[]>(defaultSpecies);
  const [concentrations, setConcentrations] = useState<Record<string, number | ''>>(defaultConcentrations);
  const [speciesColors, setSpeciesColors] = useState<Record<string, string>>(initialSpeciesColors);
  const [echartsOptions, setEchartsOptions] = useState<EChartsOption>({});
  const [showGraph, setShowGraph] = useState<boolean>(true);
  const [loading, setLoading] = useState<boolean>(false);
  const echartsRef = useRef<ReactECharts | null>(null);

  // Function to generate ECharts options (modified reactionGraphing)
  const generateEchartsOptions = useCallback((
    reactions: string[],
    ks: number[],
    C0: Record<string, number>,
    currentSpeciesColors: Record<string, string> // Pass colors
  ): EChartsOption => {
    setLoading(true);
    try {
      // Theme-aware colors
      const textColor = resolvedTheme === 'dark' ? '#ffffff' : '#000000';
      if (ks.length !== reactions.length) {
        throw new Error("The number of rate constants does not match the number of reactions.");
      }
      const orderedSpecies = detectUniqueSpeciesOrdered(reactions);
      const missingSpecies = orderedSpecies.filter(species => !C0.hasOwnProperty(species));
      if (missingSpecies.length > 0) {
        throw new Error(`The following species are missing in C0: ${missingSpecies.join(', ')}`);
      }

      function odes(_t: number, y: number[]): number[] {
        const dydt = new Array(orderedSpecies.length).fill(0);
        const currentConcentrations: Record<string, number> = {};
        orderedSpecies.forEach((species, i) => { currentConcentrations[species] = Math.max(0, y[i]); });
        reactions.forEach((reaction, i) => {
          if (!reaction) return;
          const parts = reaction.split('->'); if (parts.length !== 2) return;
          const [reactantsStr, productsStr] = parts;
          const reactantSpecies: [string, number][] = [];
          reactantsStr.split('+').forEach(speciesPart => { const match = speciesPart.trim().match(/^(?:(\d+)\*?)?(.+)$/); if (!match) return; const [, coeffStr, sp] = match; const coefficient = coeffStr ? parseInt(coeffStr) : 1; if (sp && orderedSpecies.includes(sp.trim())) reactantSpecies.push([sp.trim(), coefficient]); });
          const productSpecies: [string, number][] = [];
          productsStr.split('+').forEach(speciesPart => { const match = speciesPart.trim().match(/^(?:(\d+)\*?)?(.+)$/); if (!match) return; const [, coeffStr, sp] = match; const coefficient = coeffStr ? parseInt(coeffStr) : 1; if (sp && orderedSpecies.includes(sp.trim())) productSpecies.push([sp.trim(), coefficient]); });
          let rate = ks[i];
          for (const [sp, coeff] of reactantSpecies) rate *= Math.pow(currentConcentrations[sp], coeff);
          reactantSpecies.forEach(([sp, coeff]) => { const index = orderedSpecies.indexOf(sp); if (index !== -1) dydt[index] -= rate * coeff; });
          productSpecies.forEach(([sp, coeff]) => { const index = orderedSpecies.indexOf(sp); if (index !== -1) dydt[index] += rate * coeff; });
        });
        return dydt;
      }

      const y0 = orderedSpecies.map(species => C0[species]);
      const tSpan: [number, number] = [0, 10];
      const initial_h = 0.01;
      const tolerance = 1e-5;

      const solution = solveODERK45(odes, y0, tSpan, tolerance, initial_h);

      // --- Calculate dynamic limits ---
      const maxConcentration = Math.max(1e-9, ...solution.y.flat().filter(y => isFinite(y)));
      const relativeTolerance = maxConcentration * 1e-4;
      let steadyStateTime = tSpan[1];
      for (let i = 1; i < solution.t.length; i++) {
        const concentrationDiff = orderedSpecies.map((_, j) => Math.abs(solution.y[j][i] - solution.y[j][i - 1]));
        if (concentrationDiff.every(diff => diff < relativeTolerance)) {
          steadyStateTime = solution.t[i];
          break;
        }
      }
      steadyStateTime = steadyStateTime < tSpan[1] ? steadyStateTime * 1.1 : steadyStateTime;
      steadyStateTime = Math.min(steadyStateTime, tSpan[1]);

      const allYValues = solution.y.flat().filter(y => !isNaN(y) && isFinite(y));
      const maxY = Math.max(1e-9, ...allYValues);
      const yAxisMax = maxY * 1.1;
      const xAxisMax = steadyStateTime;
      // -------------------------------

      // Interpolate data for smoother tooltips
      const interpolateData = (t: number[], y: number[], targetPoints: number = 1000) => {
        if (t.length < 2) return { t, y };
        
        const tMin = t[0];
        const tMax = t[t.length - 1];
        const newT: number[] = [];
        const newY: number[] = [];
        
        for (let i = 0; i < targetPoints; i++) {
          const currentT = tMin + (tMax - tMin) * i / (targetPoints - 1);
          newT.push(currentT);
          
          // Linear interpolation
          let j = 0;
          while (j < t.length - 1 && t[j + 1] < currentT) j++;
          
          if (j >= t.length - 1) {
            newY.push(y[y.length - 1]);
          } else if (j === 0 && currentT <= t[0]) {
            newY.push(y[0]);
          } else {
            const t1 = t[j], t2 = t[j + 1];
            const y1 = y[j], y2 = y[j + 1];
            const interpolatedY = y1 + (y2 - y1) * (currentT - t1) / (t2 - t1);
            newY.push(interpolatedY);
          }
        }
        
        return { t: newT, y: newY };
      };

      const seriesData: SeriesOption[] = orderedSpecies.map((species, i) => {
        const interpolated = interpolateData(solution.t, solution.y[i], 1000);
        return {
          name: species, // Use plain species name
          type: 'line', smooth: true, symbol: 'none',
          data: interpolated.t.map((time, timeIndex) => [time, interpolated.y[timeIndex]]),
          itemStyle: { color: currentSpeciesColors[species] || colorPalette[i % colorPalette.length] },
          emphasis: { focus: 'series' },
          lineStyle: { width: 3.5 }
        };
      });

      return {
        backgroundColor: 'transparent',
        animation: false,
        title: {
            text: 'Concentration Profiles', left: 'center', top: '0%',
            textStyle: { 
              fontSize: 18, 
              fontFamily: 'Merriweather Sans',
              color: textColor
            }
        },
        tooltip: {
          trigger: 'axis',
          backgroundColor: resolvedTheme === 'dark' ? '#08306b' : '#ffffff',
          borderColor: resolvedTheme === 'dark' ? '#55aaff' : '#333333',
          borderWidth: 1,
          textStyle: {
            color: textColor,
            fontSize: 12,
            fontFamily: 'Merriweather Sans'
          },
          axisPointer: {
            type: 'cross',
            label: {
              show: true,
              backgroundColor: resolvedTheme === 'dark' ? '#08306b' : '#ffffff',
              color: textColor,
              borderColor: resolvedTheme === 'dark' ? '#55aaff' : '#333333',
              borderWidth: 1,
              fontFamily: 'Merriweather Sans',
              formatter: function (params: any) {
                if (params.axisDimension === 'x') {
                  return `Time: ${params.value.toFixed(2)}`;
                } else {
                  return `Conc: ${params.value.toFixed(3)}`;
                }
              }
            }
          },
          formatter: (params: any) => {
             let tooltipText = `Time: ${params[0].axisValueLabel}<br/>`;
             params.forEach((param: any) => {
                 // Manually format the plain seriesName to HTML for the tooltip
                 const formattedName = param.seriesName.replace(/(\d+)/g, '<sub>$1</sub>');
                 tooltipText += `<span style="color: ${param.color};">${formattedName}: ${param.value[1].toPrecision(3)}</span><br/>`;
             });
             return tooltipText;
          }
        },
        legend: {
          // Use plain species names for data matching.
          data: orderedSpecies.map(sp => ({ name: sp })), // Use plain species name here
          // Add a formatter to convert plain name to rich text for display
          formatter: (name: string) => {
              return name.replace(/(\d+)/g, '{sub|$1}');
          },
          textStyle: {
             color: textColor,
             fontSize: 12,
             fontFamily: 'Merriweather Sans',
             rich: { // Keep rich text config for legend rendering
                sub: {
                    fontSize: 9,
                    verticalAlign: 'bottom',
                    fontFamily: 'Merriweather Sans'
                }
             }
          },
          top: '96%',
          type: 'scroll',
          itemWidth: 25,
          itemHeight: 2
        },
        grid: { // Reduce top spacing to bring graph closer to title
          left: '5%', right: '5%', bottom: '10%', top: '3%', containLabel: true
        },
        toolbox: {
          show: false
        },
        xAxis: {
          type: 'value', name: 'Time', nameLocation: 'middle', nameGap: 30,
          min: 0,
          max: xAxisMax,
          axisLabel: {
            showMaxLabel: false, // <<< HIDE MAX LABEL for X-axis
            color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans',
            formatter: (value: number) => { // Keep formatter for non-max labels
                if (Math.abs(value) < 1e-6) return '0';
                return value.toFixed(2);
            }
          },
          nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
          axisLine: { lineStyle: { color: textColor, width: 2 } },
          axisTick: { lineStyle: { color: textColor, width: 2 } },
          splitLine: { show: false }
        },
        yAxis: {
          type: 'value', name: 'Concentration', nameLocation: 'middle', nameGap: 40,
          min: 0,
          max: yAxisMax,
          axisLabel: {
            showMaxLabel: false, // <<< HIDE MAX LABEL for Y-axis
            color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans',
            formatter: (value: number) => { // Keep formatter for non-max labels
                if (Math.abs(value) < 1e-6) return '0';
                return value.toPrecision(3);
            }
          },
          nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
          axisLine: { lineStyle: { color: textColor, width: 2 } },
          axisTick: { lineStyle: { color: textColor, width: 2 } },
          splitLine: { show: false }
        },
        series: seriesData, // series.name is now plain
      };
    } catch (error) {
      console.error("Error generating ECharts options:", error);
      return {};
    } finally {
      setLoading(false);
    }
  }, [resolvedTheme]); // Dependencies include resolvedTheme

  // --- submitForm, useEffect, and event handlers remain the same ---
  const submitForm = useCallback((
    reactionsOverride?: string[], ksOverride?: number[], concentrationsOverride?: Record<string, number>
  ) => {
    const validReactions = reactionsOverride ?? reactionInputs.filter(input => input.reactants.trim() && input.products.trim()).map(input => `${input.reactants.trim()} -> ${input.products.trim()}`);
    const validRateConstants = ksOverride ?? reactionInputs.filter(input => input.reactants.trim() && input.products.trim()).map(input => (typeof input.rateConstant === 'number' ? input.rateConstant : 1));
    const validConcentrations: Record<string, number> = {}; let allConcValid = true;
    if (concentrationsOverride) { Object.assign(validConcentrations, concentrationsOverride); }
    else { for (const sp of species) { const concValue = concentrations[sp]; if (concValue === '' || isNaN(Number(concValue))) { allConcValid = false; break; } validConcentrations[sp] = Number(concValue); } }
    const currentSpeciesForReactions = detectUniqueSpeciesOrdered(validReactions);
    const missingSpeciesInConcentrations = currentSpeciesForReactions.filter(sp => !(sp in validConcentrations));
    if (!allConcValid || missingSpeciesInConcentrations.length > 0) { console.warn("Concentrations invalid/missing.", { missingSpeciesInConcentrations }); return; }
    if (validReactions.length === 0) { console.warn("No valid reactions."); setEchartsOptions({}); setShowGraph(false); return; }
    const options = generateEchartsOptions(validReactions, validRateConstants, validConcentrations, speciesColors);
    setEchartsOptions(options); setShowGraph(true);
  }, [reactionInputs, concentrations, species, generateEchartsOptions, speciesColors]);

  useEffect(() => { submitForm(); // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  // Re-generate chart when theme changes (if reactions are confirmed)
  useEffect(() => {
    if (confirmedReactions && species.length > 0) {
      submitForm();
    }
  }, [resolvedTheme, confirmedReactions, species.length, submitForm]);

  const addReactionInput = () => { setReactionInputs([...reactionInputs, { reaction: '', reactants: '', products: '', rateConstant: 1 }]); setConfirmedReactions(false); };
  const removeReactionInput = () => { if (reactionInputs.length > 1) { setReactionInputs(reactionInputs.slice(0, -1)); setConfirmedReactions(false); } };
  const handleReactionChange = (index: number, field: 'reactants' | 'products' | 'rateConstant', value: string | number) => { const newInputs = [...reactionInputs]; const currentInput = { ...newInputs[index] }; if (field === 'reactants') currentInput.reactants = value as string; else if (field === 'products') currentInput.products = value as string; else currentInput.rateConstant = value === '' ? '' : Number(value); currentInput.reaction = `${currentInput.reactants} -> ${currentInput.products}`; newInputs[index] = currentInput; setReactionInputs(newInputs); setConfirmedReactions(false); };
  const confirmReactions = useCallback(() => { const validReactions = reactionInputs.filter(input => input.reactants.trim() && input.products.trim()).map(input => `${input.reactants.trim()} -> ${input.products.trim()}`); const validRateConstants = reactionInputs.filter(input => input.reactants.trim() && input.products.trim()).map(input => (typeof input.rateConstant === 'number' ? input.rateConstant : 1)); if (validReactions.length > 0) { const detectedSpecies = detectUniqueSpeciesOrdered(validReactions); setSpecies(detectedSpecies); const newSpeciesColors: Record<string, string> = {}; detectedSpecies.forEach((sp, index) => { newSpeciesColors[sp] = speciesColors[sp] || colorPalette[index % colorPalette.length]; }); setSpeciesColors(newSpeciesColors); const currentConcentrationsState = { ...concentrations }; const nextConcentrationsState: Record<string, number | ''> = {}; const concentrationsForPlot: Record<string, number> = {}; detectedSpecies.forEach(sp => { const currentValue = currentConcentrationsState.hasOwnProperty(sp) ? currentConcentrationsState[sp] : 1; nextConcentrationsState[sp] = currentValue; concentrationsForPlot[sp] = (typeof currentValue === 'number' && !isNaN(currentValue)) ? currentValue : 1; }); setConcentrations(nextConcentrationsState); setConfirmedReactions(true); submitForm(validReactions, validRateConstants, concentrationsForPlot); } else { setSpecies([]); setConcentrations({}); setSpeciesColors({}); setEchartsOptions({}); setShowGraph(false); setConfirmedReactions(false); } }, [reactionInputs, concentrations, submitForm, speciesColors]);
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
  const handleReactionKeyDown = (e: KeyboardEvent<HTMLInputElement>) => { if (e.key === 'Enter') { e.preventDefault(); confirmReactions(); } };
  const handleConcentrationKeyDown = (e: KeyboardEvent<HTMLInputElement>) => { if (e.key === 'Enter') { e.preventDefault(); if (confirmedReactions) { submitForm(); } } };
  const formatSliderValue = (value: number | ''): string => {
    return typeof value === 'number' ? value.toFixed(2) : '0.00';
  };

  // Helper for formatting label text in JSX (replaces direct use of formatSpeciesEcharts)
  const formatLabelHtml = (species: string): string => {
      return species.replace(/(\d+)/g, '<sub>$1</sub>');
  };


  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-8 md:px-16">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">

          {/* Column 1: Controls */}
          <div className="lg:col-span-1 space-y-6">
            {/* Reactions Card */}
            <Card>
              <CardContent className="space-y-8 pt-6 pb-3">
                {reactionInputs.map((input, index) => (
                  <div key={index} className="space-y-6 border-b pb-6 last:border-b-0 last:pb-0">
                    <div className="flex items-center gap-2">
                      <Input value={input.reactants} onChange={(e) => handleReactionChange(index, 'reactants', e.target.value)} onKeyDown={handleReactionKeyDown} placeholder="e.g., 2H2 + O2" className="flex-1"/>
                      <span className="font-bold text-lg">→</span>
                      <Input value={input.products} onChange={(e) => handleReactionChange(index, 'products', e.target.value)} onKeyDown={handleReactionKeyDown} placeholder="e.g., 2H2O" className="flex-1"/>
                    </div>
                    <div className="space-y-3">
                      <div className="flex items-center gap-3">
                        <Label htmlFor={`rate-${index}`} className="flex-1">
                          <span dangerouslySetInnerHTML={{ __html: `Rate Constant (k<sub>${index + 1}</sub>): ${formatSliderValue(input.rateConstant)}` }} />
                        </Label>
                      </div>
                      <div className="flex items-center gap-3">
                        <Slider 
                          id={`rate-${index}`} 
                          min={0} 
                          max={parseFloat(maxRateConstantSlider)} 
                          step={0.01} 
                          value={[Number(input.rateConstant || 0)]} 
                          onValueChange={(val) => handleSliderChange('rateConstant', index, val[0])} 
                          style={{ '--primary': 'hsl(var(--primary))' } as React.CSSProperties}
                          className="flex-1"
                        />
                        <div className="flex items-center gap-1">
                          <Label htmlFor={`maxRateInput-${index}`} className="text-xs text-muted-foreground">Max:</Label>
                          <Input
                            id={`maxRateInput-${index}`}
                            type="number"
                            value={maxRateConstantSlider}
                            onChange={(e) => setMaxRateConstantSlider(e.target.value)}
                            className="w-20 h-8 text-xs"
                            min="0.01"
                          />
                        </div>
                      </div>
                    </div>
                  </div>
                ))}
                <div className="flex gap-2">
                  <Button variant="outline" onClick={addReactionInput} className="flex-1">Add Reaction</Button>
                  <Button variant="outline" onClick={removeReactionInput} disabled={reactionInputs.length <= 1} className="flex-1">Remove Last</Button>
                </div>
                 <Button onClick={confirmReactions} className="w-full" disabled={confirmedReactions}>Confirm Reactions & Update Species</Button>
              </CardContent>
            </Card>

            {/* Concentrations Card */}
            {species.length > 0 && (
              <Card>
                <CardHeader><CardTitle>Initial Concentrations</CardTitle></CardHeader>
                <CardContent className="space-y-8 pt-2 pb-3">
                  {species.map((sp) => (
                    <div key={sp} className="space-y-3">
                      <div className="flex items-center gap-3">
                        <Label htmlFor={`conc-${sp}`} className="flex-1">
                          <span dangerouslySetInnerHTML={{ __html: `[${formatLabelHtml(sp)}]₀: ${formatSliderValue(concentrations[sp])}` }} />
                        </Label>
                      </div>
                      <div className="flex items-center gap-3">
                        <Slider
                          id={`conc-${sp}`}
                          min={0}
                          max={parseFloat(maxConcentrationSlider)}
                          step={0.01}
                          value={[Number(concentrations[sp] || 0)]}
                          onValueChange={(val) => handleSliderChange('concentration', sp, val[0])}
                          onKeyDown={handleConcentrationKeyDown}
                          style={{ '--primary': speciesColors[sp] || 'hsl(var(--primary))' } as React.CSSProperties}
                          className="flex-1"
                        />
                        <div className="flex items-center gap-1">
                          <Label htmlFor={`maxConcInput-${sp}`} className="text-xs text-muted-foreground">Max:</Label>
                          <Input
                            id={`maxConcInput-${sp}`}
                            type="number"
                            value={maxConcentrationSlider}
                            onChange={(e) => setMaxConcentrationSlider(e.target.value)}
                            className="w-20 h-8 text-xs"
                            min="0.01"
                          />
                        </div>
                      </div>
                    </div>
                  ))}
                </CardContent>
              </Card>
            )}
          </div>

          {/* Column 2: Plot */}
          <div className="lg:col-span-2">
            <Card>
              <CardContent className="p-1">
                <div className="relative aspect-square rounded-md overflow-hidden">
                  {loading && ( <div className="absolute inset-0 flex items-center justify-center text-muted-foreground z-10 bg-opacity-50 backdrop-blur-sm">Generating plot...</div> )}
                  {showGraph && Object.keys(echartsOptions).length > 0 ? (
                    <ReactECharts ref={echartsRef} echarts={echarts} option={echartsOptions} style={{ height: '100%', width: '100%' }} notMerge={true} lazyUpdate={false}/>
                  ) : (
                    <div className="absolute inset-0 flex items-center justify-center text-muted-foreground"> {species.length > 0 ? "Adjust parameters..." : "Enter reactions..."} </div>
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