'use client';

import React, { useState, useEffect, useCallback, useMemo, useRef } from 'react';
import { useTheme } from "next-themes";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { TooltipProvider } from "@/components/ui/tooltip";
import { Slider } from "@/components/ui/slider";
import { AlertCircle, PlusCircle, Trash2, BarChart3 } from 'lucide-react';


// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
import type { EChartsOption } from 'echarts';
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

// Import from lib files
import {
  Component,
  Kinetics,
  ReactionData,
  parseCompoundsFromReaction,
  parseParallelReactions
} from '@/lib/reaction-parser';
import {
  CalculationResults,
  solveCSTRParallel,
  solvePFR_ODE_System, // Use the new ODE system solver
  PFR_ODE_ProfilePoint, // Use the new profile point type
  findLimitingReactant
} from '@/lib/reactor-solver';
import {
  solveRecycleSystem,
  RecycleCalculationResult
} from '@/lib/recycle-solver';

// Helper types
type ReactorType = 'PFR' | 'CSTR';
type ReactionPhase = 'Liquid' | 'Gas';

const R_GAS_CONSTANT = 8.314; // J/(mol*K) or m³*kPa/(mol*K) depending on units

// Interface for a single point on the graph
interface DataPoint {
  volume: number;
  conversion: number;
  selectivity?: number;
  flowRates: { [key: string]: number };
  compositions: { [key: string]: number };
}

interface ExtendedCalculationResults extends CalculationResults {
  dataPoints: DataPoint[];
}

export default function ReactorDesignPage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  
  const [reactorType, setReactorType] = useState<ReactorType>('CSTR');
  const [reactionPhase, setReactionPhase] = useState<ReactionPhase>('Liquid');
  const [operationMode, setOperationMode] = useState<'SinglePass' | 'Recycle'>('SinglePass');
  const [reactionString, setReactionString] = useState<string>('A + B -> C');
  const [reactions, setReactions] = useState<ReactionData[]>([
    { id: '1', reactants: 'A + B', products: 'C', AValue: '1e6', EaValue: '50000', AValueBackward: '1e4', EaValueBackward: '60000', isEquilibrium: false },
    { id: '2', reactants: 'C', products: 'D', AValue: '5e5', EaValue: '45000', AValueBackward: '1e3', EaValueBackward: '55000', isEquilibrium: false }
  ]); // Parallel reactions: A+B->C and C->D
  
  // State for components (auto-generated but with editable reaction orders)
  const [components, setComponents] = useState<Component[]>([]);
  
  const [kinetics, setKinetics] = useState<Kinetics>({
    rateConstantInputMode: 'directK',
    kValue: '0.1', // Units depend on reaction order and concentration units
    AValue: '1e6',
    EaValue: '50000', // e.g., J/mol
    reactionTempK: '350', // Kelvin
  });

  const [reactorVolume, setReactorVolume] = useState<string>('100'); // e.g., cubic meters
  const [totalPressure, setTotalPressure] = useState<string>('1'); // e.g., bar (for gas phase)
  const [volumetricFlowRate, setVolumetricFlowRate] = useState<string>('1'); // v0 (for liquid phase, m³/s)
  const [maxVolumeSlider, setMaxVolumeSlider] = useState<string>('1000'); // User-defined max for volume slider
  const [maxTemperatureSlider, setMaxTemperatureSlider] = useState<string>('500'); // User-defined max for temperature slider
  const [maxFlowRateSlider, setMaxFlowRateSlider] = useState<string>('10'); // User-defined max for flow rate slider
  const [maxPressureSlider, setMaxPressureSlider] = useState<string>('50'); // User-defined max for pressure slider
  const [recycleFraction, setRecycleFraction] = useState<string>('50'); // Recycle fraction percentage (0-100%)

  const [calculationResults, setCalculationResults] = useState<CalculationResults | null>(null);
  const [recycleResults, setRecycleResults] = useState<RecycleCalculationResult | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [calculationError, setCalculationError] = useState<string | null>(null);

  // Debounce timer reference
  const debounceTimer = useRef<NodeJS.Timeout | null>(null);

  // Graph-related state
  const [showGraph, setShowGraph] = useState(false);
  const [graphType, setGraphType] = useState<'conversion' | 'selectivity' | 'flowrates' | 'composition'>('conversion');
  const [graphOptions, setGraphOptions] = useState<EChartsOption>({});
  const echartsRef = useRef<ReactECharts | null>(null);

  // Handle Enter key press to trigger calculation
  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleCalculate();
    }
  };

  // Component Management
  const addReaction = () => {
    // Cap at 6 reactions maximum
    if (reactions.length >= 6) {
      return;
    }
    const newId = (reactions.length + 1).toString();
    setReactions([...reactions, { id: newId, reactants: '', products: '', AValue: '1e6', EaValue: '50000', AValueBackward: '1e4', EaValueBackward: '60000', isEquilibrium: false }]);
  };

  const removeReaction = (id: string) => {
    const updatedReactions = reactions.filter(reaction => reaction.id !== id);
    setReactions(updatedReactions);
    
    // Update the main reaction string to reflect the first remaining reaction
    if (updatedReactions.length > 0) {
      const firstReaction = updatedReactions[0];
      const arrow = firstReaction.isEquilibrium ? '<=>' : '->';
      setReactionString(`${firstReaction.reactants} ${arrow} ${firstReaction.products}`);
    } else {
      // If no reactions left, set a default
      setReactionString('A + B -> C');
    }
  };

  const updateReaction = (id: string, field: 'reactants' | 'products' | 'AValue' | 'EaValue' | 'AValueBackward' | 'EaValueBackward' | 'isEquilibrium', value: string | boolean) => {
    setReactions(reactions.map(reaction => 
      reaction.id === id ? { ...reaction, [field]: value } : reaction
    ));
    
    // Update the main reaction string based on the updated reaction
    if (field === 'reactants' || field === 'products' || field === 'isEquilibrium') {
      const updatedReaction = reactions.find(r => r.id === id);
      if (updatedReaction) {
        const newReactants = field === 'reactants' ? value as string : updatedReaction.reactants;
        const newProducts = field === 'products' ? value as string : updatedReaction.products;
        const reactionIsEquilibrium = field === 'isEquilibrium' ? value as boolean : updatedReaction.isEquilibrium;
        const arrow = reactionIsEquilibrium ? '<=>' : '->';
        setReactionString(`${newReactants} ${arrow} ${newProducts}`);
      }
    }
  };

  const handleComponentChange = (id: string, field: keyof Component, value: string | boolean | { [reactionId: string]: string }) => {
    if (field === 'reactionOrders') {
      // Update the reaction orders for a specific component
      setComponents(prevComponents => 
        prevComponents.map(comp => 
          comp.id === id ? { ...comp, reactionOrders: value as { [reactionId: string]: string } } : comp
        )
      );
    } else if (field === 'initialFlowRate') {
      // Update the initial flow rate for a specific component
      setComponents(prevComponents => 
        prevComponents.map(comp => 
          comp.id === id ? { ...comp, initialFlowRate: value as string } : comp
        )
      );
    }
  };

  const handleKineticsChange = (field: keyof Kinetics, value: string) => {
    setKinetics(prev => ({ ...prev, [field]: value }));
  };

  // Auto-detect components from reactions and update state
  useEffect(() => {
    const allDetectedNames = new Set<string>();
    const reactantNames = new Set<string>();
    const productNames = new Set<string>();
    
    reactions.forEach(reaction => {
      const reactionStr = `${reaction.reactants} -> ${reaction.products}`;
      const names = parseCompoundsFromReaction(reactionStr);
      names.forEach(name => allDetectedNames.add(name));
      
      // Track which components are reactants vs products
      if (reaction.reactants) {
        const reactants = parseCompoundsFromReaction(reaction.reactants + ' -> dummy');
        reactants.forEach(name => reactantNames.add(name));
      }
      if (reaction.products) {
        const products = parseCompoundsFromReaction('dummy -> ' + reaction.products);
        products.forEach(name => productNames.add(name));
      }
    });
    
    const detectedNames = Array.from(allDetectedNames);
    const newComponents = detectedNames.map(name => {
      // Find existing component to preserve user input
      const existingComponent = components.find(comp => comp.name === name);
      
      // Set initial flow rate: 1 if reactant in any reaction, 0 if only product
      // Special case: C should start with 0 flow rate (intermediate product)
      const initialFlowRate = name === 'C' ? '0' : (reactantNames.has(name) ? '1' : '0');
      
      // Initialize reaction orders for each reaction
      const reactionOrders: { [reactionId: string]: string } = {};
      reactions.forEach(reaction => {
        const isReactantInReaction = reaction.reactants && reaction.reactants.toLowerCase().includes(name.toLowerCase());
        const isProductInReaction = reaction.products && reaction.products.toLowerCase().includes(name.toLowerCase());
        
        if (reaction.isEquilibrium) {
          // For equilibrium reactions, add both forward and backward keys
          const forwardKey = `${reaction.id}-fwd`;
          const backwardKey = `${reaction.id}-rev`;
          
          if (isReactantInReaction) {
            reactionOrders[forwardKey] = existingComponent?.reactionOrders?.[forwardKey] || '1';
          }
          if (isProductInReaction) {
            reactionOrders[backwardKey] = existingComponent?.reactionOrders?.[backwardKey] || '1';
          }
        } else {
          // For forward-only reactions
          if (isReactantInReaction) {
            reactionOrders[reaction.id] = existingComponent?.reactionOrders?.[reaction.id] || '1';
          }
        }
      });
      
      return {
        id: `comp-${name}`,
        name,
        initialFlowRate: existingComponent?.initialFlowRate || initialFlowRate,
        reactionOrders
      };
    });
    
    // Only update if components actually changed to prevent loops
    if (JSON.stringify(newComponents) !== JSON.stringify(components)) {
      setComponents(newComponents);
    }
  }, [reactions]); // Only depend on reactions, not components to avoid loops

  // Memoize parsed parallel reactions to prevent unnecessary recalculations
  const parsedParallelReactions = useMemo(() => {
    const parallelResult = parseParallelReactions(reactionString, components, kinetics, reactions, false);
    return parallelResult;
  }, [reactionString, components, kinetics, reactions]);

  // Update reaction string when reactions change (but don't trigger infinite loops)
  useEffect(() => {
    if (reactions.length > 0) {
      const firstReaction = reactions[0];
      const arrow = firstReaction.isEquilibrium ? '<=>' : '->';
      const newReactionString = `${firstReaction.reactants} ${arrow} ${firstReaction.products}`;
      
      // Only update if different to prevent loops
      if (newReactionString !== reactionString) {
        setReactionString(newReactionString);
      }
    }
  }, [reactions, reactionString]);

  // Main Calculation Logic - Updated for single pass and recycle modes
  const handleCalculate = useCallback(async () => {
    setIsLoading(true);
    setCalculationResults(null);
    setRecycleResults(null);
    setCalculationError(null);

    if (parsedParallelReactions.error || parsedParallelReactions.reactions.length === 0) {
      setCalculationError(parsedParallelReactions.error || "No valid reactions to process.");
      setIsLoading(false);
      return;
    }
    
    // Common parameters
    const V = parseFloat(reactorVolume);
    const T_K = parseFloat(kinetics.reactionTempK);
    const v0_input_Ls = parseFloat(volumetricFlowRate);
    const initialFlowRates = components.reduce((acc, comp) => {
      acc[comp.name] = parseFloat(comp.initialFlowRate) || 0;
      return acc;
    }, {} as { [key: string]: number });

    let v0_calc: number;
    if (reactionPhase === 'Liquid') {
      if (isNaN(v0_input_Ls) || v0_input_Ls <= 0) {
        setCalculationError("Liquid phase requires a positive volumetric flow rate.");
        setIsLoading(false);
        return;
      }
      v0_calc = v0_input_Ls;
    } else { // Gas Phase
      const P_total_kPa = parseFloat(totalPressure) * 100;
      if (isNaN(P_total_kPa) || P_total_kPa <= 0) {
        setCalculationError("Gas phase requires a positive total pressure.");
        setIsLoading(false);
        return;
      }
      const F_total0 = Object.values(initialFlowRates).reduce((sum, f) => sum + f, 0);
      if (T_K <= 0) {
        setCalculationError("Gas phase requires a positive absolute temperature.");
        setIsLoading(false);
        return;
      }
      const R_GAS_CONSTANT = 8.314;
      v0_calc = (F_total0 * R_GAS_CONSTANT * T_K) / P_total_kPa;
    }

    try {
      if (operationMode === 'Recycle') {
        // Solve recycle system
        const recycleResult = await solveRecycleSystem(
          parsedParallelReactions,
          V,
          initialFlowRates,
          v0_input_Ls,
          components,
          reactorType,
          reactionPhase,
          kinetics,
          parseFloat(recycleFraction),
          totalPressure
        );
        
        if (recycleResult.calculationError || recycleResult.error) {
          setCalculationError(recycleResult.calculationError || recycleResult.error || "Recycle calculation failed.");
        } else {
          setRecycleResults(recycleResult);
        }
      } else {
        // Single pass operation
        let finalFlowRates: { [key: string]: number };

        if (reactorType === 'PFR') {
          const pfrProfile = solvePFR_ODE_System(
            parsedParallelReactions, V, initialFlowRates, v0_calc, components.map(c => c.name)
          );
          finalFlowRates = pfrProfile.length > 0 ? pfrProfile[pfrProfile.length - 1].flowRates : initialFlowRates;
        } else { // CSTR
          finalFlowRates = solveCSTRParallel(
              parsedParallelReactions, V, initialFlowRates, v0_calc
          );
        }
        
        const limitingReactant = findLimitingReactant(parsedParallelReactions, components);
        let conversion: { reactantName: string; value: number } | undefined;
        if (limitingReactant) {
          const F_A0 = limitingReactant.initialFlow;
          const F_A = finalFlowRates[limitingReactant.name];
          const X = F_A0 > 0 ? (F_A0 - F_A) / F_A0 : 0;
          conversion = { reactantName: limitingReactant.name, value: Math.max(0, Math.min(X, 1)) };
        }

        setCalculationResults({ conversion, outletFlowRates: finalFlowRates });
      }

    } catch (e: any) {
      setCalculationError(`Calculation Error: ${e.message}`);
    } finally {
      setIsLoading(false);
    }
  }, [parsedParallelReactions, reactorType, reactorVolume, kinetics.reactionTempK, components, reactionPhase, totalPressure, volumetricFlowRate, operationMode, recycleFraction]);

  // Add helper function for 3 significant figures formatting
  const formatToSigFigs = (num: number, sigFigs: number = 3): string => {
    if (num === 0) return '0.00';
    const magnitude = Math.floor(Math.log10(Math.abs(num)));
    const factor = Math.pow(10, sigFigs - 1 - magnitude);
    return (Math.round(num * factor) / factor).toString();
  };

  // Debounced calculation trigger
  const triggerDebouncedCalculation = useCallback(() => {
    if (debounceTimer.current) {
      clearTimeout(debounceTimer.current);
    }
    debounceTimer.current = setTimeout(() => {
      handleCalculate();
    }, 1); // 1ms debounce delay
  }, [handleCalculate]);

  // Effect to run calculation when primary inputs change (debounced)
  useEffect(() => {
    triggerDebouncedCalculation();
  }, [
    reactorVolume, 
    kinetics.reactionTempK, 
    volumetricFlowRate, 
    totalPressure,
    recycleFraction,
    triggerDebouncedCalculation
  ]);


  // --- Graphing Logic ---
  // Generate graph data when results are available
  const generateGraphData = useCallback(() => {
    if (!parsedParallelReactions || parsedParallelReactions.error || parsedParallelReactions.reactions.length === 0) {
      setGraphOptions({});
      return;
    }

    const maxVol = parseFloat(maxVolumeSlider);
    const T_K = parseFloat(kinetics.reactionTempK);
    const v0_input_Ls = parseFloat(volumetricFlowRate);
    const initialFlowRates = components.reduce((acc, comp) => {
        acc[comp.name] = parseFloat(comp.initialFlowRate) || 0;
        return acc;
    }, {} as { [key: string]: number });
    
    let v0_calc: number;
    if (reactionPhase === 'Liquid') {
        if (isNaN(v0_input_Ls) || v0_input_Ls <= 0) return;
        v0_calc = v0_input_Ls;
    } else {
        const P_total_kPa = parseFloat(totalPressure) * 100;
        if (isNaN(P_total_kPa) || P_total_kPa <= 0) return;
        const F_total0 = Object.values(initialFlowRates).reduce((sum, f) => sum + f, 0);
        if (T_K <= 0) return;
        const R_GAS_CONSTANT = 8.314;
        v0_calc = (F_total0 * R_GAS_CONSTANT * T_K) / P_total_kPa;
    }

    let profile: PFR_ODE_ProfilePoint[];
    if (reactorType === 'PFR') {
        profile = solvePFR_ODE_System(parsedParallelReactions, maxVol, initialFlowRates, v0_calc, components.map(c => c.name));
    } else { // CSTR
        profile = [];
        const steps = 50; // Reduced from 100 for better performance
        for (let i = 1; i <= steps; i++) {
            const vol = (i / steps) * maxVol;
            const flowRates = solveCSTRParallel(parsedParallelReactions, vol, initialFlowRates, v0_calc);
            profile.push({ volume: vol, flowRates });
        }
    }

    const limitingReactant = findLimitingReactant(parsedParallelReactions, components);
    if (!limitingReactant) return;
    const F_A0 = limitingReactant.initialFlow;
    
    const dataPoints: DataPoint[] = profile.map(point => {
        const { volume, flowRates } = point;
        const F_A = flowRates[limitingReactant!.name];
        const conversion = F_A0 > 0 ? (F_A0 - F_A) / F_A0 : 0;
        
        // --- This is the new block for Overall Selectivity ---
        // Note: This calculation is for the Overall Selectivity of component 'C'.
        
        // 1. Calculate the net moles of the desired product ('C') formed up to this point.
        const moles_product_formed = (flowRates['C'] || 0) - (initialFlowRates['C'] || 0);

        // 2. Calculate the moles of the limiting reactant consumed up to this point.
        const moles_reactant_consumed = F_A0 - flowRates[limitingReactant!.name];

        let selectivity: number;
        // 3. Calculate selectivity, avoiding division by zero at the reactor inlet.
        if (moles_reactant_consumed < 1e-9) {
          selectivity = 0; // At the very beginning, selectivity can be considered zero.
        } else {
          selectivity = moles_product_formed / moles_reactant_consumed;
        }
        // --- End of new block ---

        // Calculate compositions (mol%)
        const totalFlow = Object.values(flowRates).reduce((sum, flow) => sum + flow, 0);
        const compositions: { [key: string]: number } = {};
        if (totalFlow > 1e-9) {
            Object.keys(flowRates).forEach(componentName => {
                compositions[componentName] = (flowRates[componentName] / totalFlow) * 100;
            });
        } else {
            Object.keys(flowRates).forEach(componentName => {
                compositions[componentName] = 0;
            });
        }

        return { volume, conversion, selectivity, flowRates, compositions };
    });

    // Charting logic
    // Theme-dependent colors
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';
    
    let chartOptions: EChartsOption = {
        backgroundColor: 'transparent',
        animation: false,
        animationDuration: 0,
        title: { text: 'Conversion vs Volume', left: 'center', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } },
        grid: { left: '10%', right: '10%', bottom: '15%', top: '15%', containLabel: true },
        xAxis: { 
          type: 'value', 
          name: 'Reactor Volume (m³)', 
          nameLocation: 'middle', 
          nameGap: 30, 
          nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
          axisLine: { lineStyle: { color: textColor } }, 
          axisTick: { lineStyle: { color: textColor } },
          axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          splitLine: { show: false }
        },
        yAxis: { 
          type: 'value', 
          name: 'Percentage (%)', 
          nameLocation: 'middle', 
          nameGap: 40, 
          nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
          axisLine: { lineStyle: { color: textColor } }, 
          axisTick: { lineStyle: { color: textColor } },
          axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          splitLine: { show: false }
        },
        legend: { 
          orient: 'horizontal',
          bottom: '5%', 
          left: 'center',
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          data: ['Conversion', 'Selectivity'],
          itemWidth: 25,
          itemHeight: 2
        },
        tooltip: { 
          trigger: 'axis', 
          backgroundColor: resolvedTheme === 'dark' ? '#08306b' : '#ffffff', 
          borderColor: resolvedTheme === 'dark' ? '#55aaff' : '#333333', 
          borderWidth: 1,
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
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
                  return `Volume: ${params.value.toFixed(2)} m³`;
                } else {
                  return `${params.value.toFixed(1)}%`;
                }
              }
            }
          },
          formatter: function(params: any) {
            if (!Array.isArray(params)) return '';
            const volume = params[0]?.axisValue || 0;
            let tooltipContent = `<div style="color: ${textColor};">Volume: ${formatToSigFigs(volume)} m³<br/>`;
            
            params.forEach((param: any) => {
              if (param.seriesName === 'Conversion' || param.seriesName === 'Selectivity') {
                const value = formatToSigFigs(param.value[1]);
                const color = param.color;
                tooltipContent += `<span style="color: ${color};">${param.seriesName}: ${value}%</span><br/>`;
              }
            });
            
            tooltipContent += '</div>';
            return tooltipContent;
          }
        },
        series: []
    };

    if (graphType === 'conversion') {
        // Find the data point closest to the current reactor volume
        const currentVolume = parseFloat(reactorVolume);
        const closestPoint = dataPoints.reduce((closest, point) => 
            Math.abs(point.volume - currentVolume) < Math.abs(closest.volume - currentVolume) ? point : closest
        );

        chartOptions.title = { text: 'Conversion vs Volume', left: 'center', top: '3%', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } };
        chartOptions.legend = { 
          orient: 'horizontal',
          bottom: '5%', 
          left: 'center',
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          data: ['Conversion']
        };

        chartOptions.series = [
            { 
              name: 'Conversion', 
              type: 'line', 
              showSymbol: false, 
              data: dataPoints.map(p => [p.volume, p.conversion * 100]), 
              color: '#FF6347', // Add this for legend and tooltip colors
              lineStyle: { color: '#FF6347', width: 2 },
              emphasis: { lineStyle: { width: 3 } },
              markPoint: {
                data: [{
                  name: 'Current Volume',
                  coord: [closestPoint.volume, closestPoint.conversion * 100],
                  symbol: 'circle',
                  symbolSize: 8,
                  itemStyle: { color: '#FF0000', borderColor: '#FFFFFF', borderWidth: 2 },
                  label: { show: false }
                }]
              }
            }
        ];
    } else if (graphType === 'selectivity') {
        // Find the data point closest to the current reactor volume
        const currentVolume = parseFloat(reactorVolume);
        const closestPoint = dataPoints.reduce((closest, point) => 
            Math.abs(point.volume - currentVolume) < Math.abs(closest.volume - currentVolume) ? point : closest
        );

        chartOptions.title = { text: 'Selectivity vs Volume', left: 'center', top: '3%', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } };
        chartOptions.legend = { 
          orient: 'horizontal',
          bottom: '5%', 
          left: 'center',
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          data: ['Selectivity']
        };

        chartOptions.series = [
            { 
              name: 'Selectivity', 
              type: 'line', 
              showSymbol: false, 
              data: dataPoints.filter(p => p.volume > 0).map(p => [p.volume, (p.selectivity ?? 0) * 100]), 
              color: '#00BFFF', // Add this for legend and tooltip colors
              lineStyle: { color: '#00BFFF', width: 2 },
              emphasis: { lineStyle: { width: 3 } },
              markPoint: {
                data: [{
                  name: 'Current Volume',
                  coord: [closestPoint.volume, (closestPoint.selectivity ?? 0) * 100],
                  symbol: 'circle',
                  symbolSize: 8,
                  itemStyle: { color: '#FF0000', borderColor: '#FFFFFF', borderWidth: 2 },
                  label: { show: false }
                }]
              }
            }
        ];
    } else if (graphType === 'flowrates') {
        chartOptions.title = { text: 'Flow Rates vs Volume', left: 'center', top: '3%', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } };
        chartOptions.yAxis = { 
          type: 'value', 
          name: 'Flow Rate (mol/s)', 
          nameLocation: 'middle', 
          nameGap: 50, 
          nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
          axisLine: { lineStyle: { color: textColor } }, 
          axisTick: { lineStyle: { color: textColor } },
          axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          splitLine: { show: false }
        };
        chartOptions.legend = {
          orient: 'horizontal',
          bottom: '5%', 
          left: 'center',
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          data: components.map(c => c.name)
        };
        chartOptions.tooltip = { 
          trigger: 'axis', 
          backgroundColor: isDark ? '#08306b' : '#ffffff', 
          borderColor: isDark ? '#55aaff' : '#333333', 
          borderWidth: 1,
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          formatter: function(params: any) {
            if (!Array.isArray(params)) return '';
            const volume = params[0]?.axisValue || 0;
            let tooltipContent = `<div style="color: ${textColor};">Volume: ${formatToSigFigs(volume)} m³<br/>`;
            
            params.forEach((param: any) => {
              const value = formatToSigFigs(param.value[1]);
              const color = param.color;
              tooltipContent += `<span style="color: ${color};">● ${param.seriesName}: ${value} mol/s</span><br/>`;
            });
            
            tooltipContent += '</div>';
            return tooltipContent;
          }
        };
        const colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#feca57', '#ff9ff3'];
        
        // Find the data point closest to the current reactor volume for flow rates chart
        const currentVolume = parseFloat(reactorVolume);
        const closestPoint = dataPoints.reduce((closest, point) => 
            Math.abs(point.volume - currentVolume) < Math.abs(closest.volume - currentVolume) ? point : closest
        );
        
        chartOptions.series = components.map((comp, i) => ({
            name: comp.name,
            type: 'line',
            showSymbol: false,
            data: dataPoints.map(p => [p.volume, p.flowRates[comp.name]]),
            color: colors[i % colors.length], // Add this for legend and tooltip colors
            lineStyle: { color: colors[i % colors.length], width: 2 },
            emphasis: { lineStyle: { width: 3 } },
            markPoint: {
              data: [{
                name: 'Current Volume',
                coord: [closestPoint.volume, closestPoint.flowRates[comp.name]],
                symbol: 'circle',
                symbolSize: 8,
                itemStyle: { color: '#FF0000', borderColor: '#FFFFFF', borderWidth: 2 },
                label: { show: false }
              }]
            }
        }));
    } else { // 'composition'
        chartOptions.title = { text: 'Composition vs Volume', left: 'center', top: '3%', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } };
        chartOptions.yAxis = { 
          type: 'value', 
          name: 'Composition (mol%)', 
          nameLocation: 'middle', 
          nameGap: 50, 
          nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
          axisLine: { lineStyle: { color: textColor } }, 
          axisTick: { lineStyle: { color: textColor } },
          axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          splitLine: { show: false }
        };
        chartOptions.legend = {
          orient: 'horizontal',
          bottom: '5%', 
          left: 'center',
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          data: components.map(c => c.name)
        };
        chartOptions.tooltip = { 
          trigger: 'axis', 
          backgroundColor: isDark ? '#08306b' : '#ffffff', 
          borderColor: isDark ? '#55aaff' : '#333333', 
          borderWidth: 1,
          textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
          formatter: function(params: any) {
            if (!Array.isArray(params)) return '';
            const volume = params[0]?.axisValue || 0;
            let tooltipContent = `<div style="color: ${textColor};">Volume: ${formatToSigFigs(volume)} m³<br/>`;
            
            params.forEach((param: any) => {
              const value = formatToSigFigs(param.value[1]);
              const color = param.color;
              tooltipContent += `<span style="color: ${color};">● ${param.seriesName}: ${value}%</span><br/>`;
            });
            
            tooltipContent += '</div>';
            return tooltipContent;
          }
        };
        const colors = ['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#feca57', '#ff9ff3'];
        
        // Find the data point closest to the current reactor volume for composition chart
        const currentVolume = parseFloat(reactorVolume);
        const closestPoint = dataPoints.reduce((closest, point) => 
            Math.abs(point.volume - currentVolume) < Math.abs(closest.volume - currentVolume) ? point : closest
        );
        
        chartOptions.series = components.map((comp, i) => ({
            name: comp.name,
            type: 'line',
            showSymbol: false,
            data: dataPoints.map(p => [p.volume, p.compositions[comp.name]]),
            color: colors[i % colors.length], // Add this for legend and tooltip colors
            lineStyle: { color: colors[i % colors.length], width: 2 },
            emphasis: { lineStyle: { width: 3 } },
            markPoint: {
              data: [{
                name: 'Current Volume',
                coord: [closestPoint.volume, closestPoint.compositions[comp.name]],
                symbol: 'circle',
                symbolSize: 8,
                itemStyle: { color: '#FF0000', borderColor: '#FFFFFF', borderWidth: 2 },
                label: { show: false }
              }]
            }
        }));
    }

    setGraphOptions(chartOptions);

  }, [parsedParallelReactions, reactorType, reactorVolume, maxVolumeSlider, kinetics.reactionTempK, volumetricFlowRate, components, reactionPhase, totalPressure, graphType, resolvedTheme]);


  // Update graph when calculation results change - immediate updates
  useEffect(() => {
    if (showGraph && calculationResults && operationMode === 'SinglePass') {
      generateGraphData();
    }
  }, [showGraph, calculationResults, graphType, reactorType, maxVolumeSlider, operationMode, generateGraphData]);

  // Reactor SVG Visualization
  const renderReactorSVG = () => {
    // ... SVG code is unchanged ...
    if (reactorType === 'PFR') {
      return (
        <div className="w-full h-64 bg-card rounded flex items-center justify-center">
          <svg viewBox="0 0 400 200" className="w-full h-full">
            {/* PFR Tube - slightly smaller */}
            <rect x="60" y="70" width="280" height="60" fill="lightblue" stroke="currentColor" strokeWidth="2" rx="30"/>
            
            {/* Feed Arrow - longer tail */}
            <defs>
              <marker id="arrowhead" markerWidth="7" markerHeight="5" refX="6" refY="2.5" orient="auto">
                <polygon points="0 0, 7 2.5, 0 5" fill="currentColor"/>
              </marker>
            </defs>
            <line x1="10" y1="100" x2="55" y2="100" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="32" y="85" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Feed</text>
            
            {/* Product Arrow - longer tail */}
            <line x1="345" y1="100" x2="390" y2="100" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="367" y="85" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Product</text>
            
            {/* PFR Label */}
            <text x="200" y="110" fontSize="20" fill="currentColor" className="text-foreground" textAnchor="middle" fontWeight="bold">PFR</text>
          </svg>
        </div>
      );
    } else {
      return (
        <div className="w-full h-64 bg-card rounded flex items-center justify-center">
          <svg viewBox="0 0 400 240" className="w-full h-full">
            {/* CSTR Tank - even larger rectangular shape matching CSTRVisualization.tsx */}
            <rect x="80" y="20" width="240" height="180" fill="lightblue" stroke="currentColor" strokeWidth="4" rx="25"/>
            
            {/* Feed Arrow - longer tail */}
            <defs>
              <marker id="arrowhead" markerWidth="7" markerHeight="5" refX="6" refY="2.5" orient="auto">
                <polygon points="0 0, 7 2.5, 0 5" fill="currentColor"/>
              </marker>
            </defs>
            <line x1="20" y1="105" x2="75" y2="105" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="47" y="90" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Feed</text>
            
            {/* Product Arrow - longer tail */}
            <line x1="325" y1="105" x2="380" y2="105" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="352" y="90" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Product</text>
            
            {/* Stirrer shaft */}
            <line x1="200" y1="0" x2="200" y2="155" stroke="currentColor" strokeWidth="8"/>
            
            {/* Stirrer blades with spinning animation - matching CSTRVisualization.tsx */}
            <g className="animate-spin" style={{transformOrigin: '200px 155px'}}>
              {/* A horizontal blade centered at (200, 155) */}
              <line x1="160" y1="155" x2="240" y2="155" stroke="currentColor" strokeWidth="8"/>
              {/* A vertical blade of the same length, also centered at (200, 155) */}
              <line x1="200" y1="115" x2="200" y2="195" stroke="currentColor" strokeWidth="8"/>
            </g>
            
            {/* CSTR Label */}
            <text x="200" y="230" fontSize="28" fill="currentColor" className="text-foreground" textAnchor="middle" fontWeight="bold">CSTR</text>
          </svg>
        </div>
      );
    }
  };


  return (
    // ... The entire JSX return block is unchanged ...
    <TooltipProvider>
      <div className="container mx-auto p-4 space-y-6">
        <div className="grid grid-cols-1 md:grid-cols-5 gap-6"> {/* Changed md:grid-cols-2 to md:grid-cols-5 */}
          {/* Left Column: Inputs */}
          <div className="space-y-6 md:col-span-2"> {/* Added md:col-span-2 */}
            {/* Reactor Type and Phase Selection */}
            <Card>
              <CardContent className="space-y-6 pt-3 relative">
                {/* Phase Selection - Top Left */}
                <div className="absolute top-3 left-3">
                  <div className="flex items-center gap-2 bg-muted rounded-lg p-1">
                    <Button
                      variant={reactionPhase === 'Liquid' ? 'default' : 'ghost'}
                      size="sm"
                      onClick={() => setReactionPhase('Liquid')}
                      className="h-8 px-3"
                    >
                      Liquid
                    </Button>
                    <Button
                      variant={reactionPhase === 'Gas' ? 'default' : 'ghost'}
                      size="sm"
                      onClick={() => setReactionPhase('Gas')}
                      className="h-8 px-3"
                    >
                      Gas
                    </Button>
                  </div>
                </div>
                
                {/* Operation Mode Selection - Top Right */}
                <div className="absolute top-3 right-3">
                  <div className="flex items-center gap-2 bg-muted rounded-lg p-1">
                    <Button
                      variant={operationMode === 'SinglePass' ? 'default' : 'ghost'}
                      size="sm"
                      onClick={() => setOperationMode('SinglePass')}
                      className="h-8 px-3"
                    >
                      Single Pass
                    </Button>
                    <Button
                      variant={operationMode === 'Recycle' ? 'default' : 'ghost'}
                      size="sm"
                      onClick={() => {
                        setOperationMode('Recycle');
                        setShowGraph(false); // Turn off graph mode when switching to recycle
                      }}
                      className="h-8 px-3"
                    >
                      Recycle
                    </Button>
                  </div>
                </div>
                
                {/* Add top margin to account for absolute positioned buttons */}
                <div className="mt-12"></div>
                
                {/* Reactor Volume Slider */}
                <div className="space-y-3">
                  <div className="flex items-center gap-3">
                    <div className="flex items-center gap-1">
                      <Label htmlFor="reactorVolumeSlider" className="text-sm font-medium whitespace-nowrap">Reactor Volume:</Label>
                      <span className="text-sm font-medium w-12 text-right">{reactorVolume}</span>
                      <span className="text-xs text-muted-foreground">m³</span>
                    </div>
                    <Slider
                      id="reactorVolumeSlider"
                      min={1}
                      max={parseFloat(maxVolumeSlider)}
                      step={1}
                      value={[parseFloat(reactorVolume)]}
                      onValueChange={(value) => setReactorVolume(value[0].toString())}
                      className="flex-1"
                    />
                    <div className="flex items-center gap-1">
                      <Label htmlFor="maxVolumeInput" className="text-xs text-muted-foreground">Max:</Label>
                      <Input
                        id="maxVolumeInput"
                        type="number"
                        value={maxVolumeSlider}
                        onChange={(e) => setMaxVolumeSlider(e.target.value)}
                        className="w-20 h-8 text-xs"
                        min="1"
                      />
                    </div>
                  </div>
                </div>

                {/* Temperature Slider */}
                <div className="space-y-3">
                  <div className="flex items-center gap-3">
                    <div className="flex items-center gap-1">
                      <Label htmlFor="temperatureSlider" className="text-sm font-medium whitespace-nowrap">Temperature:</Label>
                      <span className="text-sm font-medium w-12 text-right">{kinetics.reactionTempK}</span>
                      <span className="text-xs text-muted-foreground">K</span>
                    </div>
                    <Slider
                      id="temperatureSlider"
                      min={250}
                      max={parseFloat(maxTemperatureSlider)}
                      step={1}
                      value={[parseFloat(kinetics.reactionTempK)]}
                      onValueChange={(value) => handleKineticsChange('reactionTempK', value[0].toString())}
                      className="flex-1"
                    />
                    <div className="flex items-center gap-1">
                      <Label htmlFor="maxTempInput" className="text-xs text-muted-foreground">Max:</Label>
                      <Input
                        id="maxTempInput"
                        type="number"
                        value={maxTemperatureSlider}
                        onChange={(e) => setMaxTemperatureSlider(e.target.value)}
                        className="w-20 h-8 text-xs"
                        min="250"
                      />
                    </div>
                  </div>
                </div>

                {/* Volumetric Flow Rate / Pressure */}
                <div className="space-y-3">
                  {reactionPhase === 'Gas' && (
                    <div className="flex items-center gap-3">
                      <div className="flex items-center gap-1">
                        <Label htmlFor="pressureSlider" className="text-sm font-medium whitespace-nowrap">Total Pressure:</Label>
                        <span className="text-sm font-medium w-12 text-right">{totalPressure}</span>
                        <span className="text-xs text-muted-foreground">bar</span>
                      </div>
                      <Slider
                        id="pressureSlider"
                        min={0.1}
                        max={parseFloat(maxPressureSlider)}
                        step={0.1}
                        value={[parseFloat(totalPressure)]}
                        onValueChange={(value) => setTotalPressure(value[0].toString())}
                        className="flex-1"
                      />
                      <div className="flex items-center gap-1">
                        <Label htmlFor="maxPressureInput" className="text-xs text-muted-foreground">Max:</Label>
                        <Input
                          id="maxPressureInput"
                          type="number"
                          value={maxPressureSlider}
                          onChange={(e) => setMaxPressureSlider(e.target.value)}
                          className="w-20 h-8 text-xs"
                          min="0.1"
                        />
                      </div>
                    </div>
                  )}
                  {reactionPhase === 'Liquid' && (
                    <div className="flex items-center gap-3">
                      <div className="flex items-center gap-1">
                        <Label htmlFor="volumetricFlowRateSlider" className="text-sm font-medium whitespace-nowrap">Vol Flow Rate:</Label>
                        <span className="text-sm font-medium w-12 text-right">{volumetricFlowRate}</span>
                        <span className="text-xs text-muted-foreground">m³/s</span>
                      </div>
                      <Slider
                        id="volumetricFlowRateSlider"
                        min={0.1}
                        max={parseFloat(maxFlowRateSlider)}
                        step={0.1}
                        value={[parseFloat(volumetricFlowRate)]}
                        onValueChange={(value) => setVolumetricFlowRate(value[0].toString())}
                        className="flex-1"
                      />
                      <div className="flex items-center gap-1">
                        <Label htmlFor="maxFlowRateInput" className="text-xs text-muted-foreground">Max:</Label>
                        <Input
                          id="maxFlowRateInput"
                          type="number"
                          value={maxFlowRateSlider}
                          onChange={(e) => setMaxFlowRateSlider(e.target.value)}
                          className="w-20 h-8 text-xs"
                          min="0.1"
                          step="0.1"
                        />
                      </div>
                    </div>
                  )}
                </div>

                {/* Recycle Fraction Slider - Only for Recycle Mode */}
                {operationMode === 'Recycle' && (
                  <div className="space-y-3">
                    <div className="flex items-center gap-3">
                      <div className="flex items-center gap-1">
                        <Label htmlFor="recycleFractionSlider" className="text-sm font-medium whitespace-nowrap">Recycle Fraction:</Label>
                        <span className="text-sm font-medium w-12 text-right">{recycleFraction}</span>
                        <span className="text-xs text-muted-foreground">%</span>
                      </div>
                      <Slider
                        id="recycleFractionSlider"
                        min={0}
                        max={100}
                        step={1}
                        value={[parseFloat(recycleFraction)]}
                        onValueChange={(value) => setRecycleFraction(value[0].toString())}
                        className="flex-1"
                      />
                      <div className="flex items-center gap-1">
                        <span className="text-xs text-muted-foreground w-20">0-100%</span>
                      </div>
                    </div>
                  </div>
                )}
              </CardContent>
            </Card>

            {/* Reaction Stoichiometry */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center justify-between">
                  Reaction Stoichiometry
                  <Button 
                    variant="outline" 
                    size="sm" 
                    onClick={addReaction}
                    disabled={reactions.length >= 6}
                  >
                    <PlusCircle className="mr-2 h-4 w-4" />Add Reaction
                  </Button>
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                {reactions.map((reaction, index) => (
                  <div key={reaction.id}>
                    <div className="grid grid-cols-12 gap-3 items-center">
                      <div className="col-span-2 flex items-center justify-start">
                        <Label className="text-sm font-medium whitespace-nowrap">Reaction {index + 1}:</Label>
                      </div>
                      <div className="col-span-3">
                        <Input
                          placeholder="A + B"
                          value={reaction.reactants}
                          onChange={(e) => updateReaction(reaction.id, 'reactants', e.target.value)}
                        />
                      </div>
                      <div className="col-span-2 text-center flex items-center justify-center">
                        <div className="flex flex-col items-center gap-2">
                          <div className="flex items-center gap-2 bg-muted rounded-lg p-1">
                            <Button
                              variant={!reaction.isEquilibrium ? 'default' : 'ghost'}
                              size="sm"
                              onClick={() => updateReaction(reaction.id, 'isEquilibrium', false)}
                              className="h-8 px-3"
                            >
                              →
                            </Button>
                            <Button
                              variant={reaction.isEquilibrium ? 'default' : 'ghost'}
                              size="sm"
                              onClick={() => updateReaction(reaction.id, 'isEquilibrium', true)}
                              className="h-8 px-3"
                            >
                              ⇌
                            </Button>
                          </div>
                        </div>
                      </div>
                      <div className="col-span-3">
                        <Input
                          placeholder="C"
                          value={reaction.products}
                          onChange={(e) => updateReaction(reaction.id, 'products', e.target.value)}
                        />
                      </div>
                      <div className="col-span-1"></div>
                      <div className="col-span-1 flex justify-start">
                        {index > 0 && (
                          <button
                            onClick={() => removeReaction(reaction.id)}
                            className="text-red-500 hover:text-red-700"
                          >
                            <Trash2 className="h-4 w-4" />
                          </button>
                        )}
                      </div>
                    </div>
                    
                    {/* A and Ea values for this reaction */}
                    <div className="mt-2 pl-4 space-y-3">
                      {/* Forward reaction parameters */}
                      <div className="grid grid-cols-2 gap-4">
                        <div className="flex items-center gap-2">
                          <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                            {reaction.isEquilibrium ? (
                              <span>A<sub style={{fontSize: '0.75em', lineHeight: '1'}}>forward</sub>:</span>
                            ) : (
                              'A:'
                            )}
                          </Label>
                          <Input
                            type="number"
                            value={reaction.AValue}
                            onChange={(e) => updateReaction(reaction.id, 'AValue', e.target.value)}
                            placeholder="1e6"
                            className="flex-1"
                          />
                        </div>
                        <div className="flex items-center gap-2">
                          <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                            {reaction.isEquilibrium ? (
                              <span>Ea<sub style={{fontSize: '0.75em', lineHeight: '1'}}>forward</sub>:</span>
                            ) : (
                              'Ea:'
                            )}
                          </Label>
                          <div className="flex items-center gap-1 flex-1">
                            <Input
                              type="number"
                              value={reaction.EaValue}
                              onChange={(e) => updateReaction(reaction.id, 'EaValue', e.target.value)}
                              placeholder="50000"
                            />
                            <span className="text-xs text-muted-foreground">J/mol</span>
                          </div>
                        </div>
                      </div>
                      
                      {/* Backward reaction parameters (only shown for equilibrium) */}
                      {reaction.isEquilibrium && (
                        <div className="grid grid-cols-2 gap-4">
                          <div className="flex items-center gap-2">
                            <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                              <span>A<sub style={{fontSize: '0.75em', lineHeight: '1'}}>backward</sub>:</span>
                            </Label>
                            <Input
                              type="number"
                              value={reaction.AValueBackward || ''}
                              onChange={(e) => updateReaction(reaction.id, 'AValueBackward', e.target.value)}
                              placeholder="1e4"
                              className="flex-1"
                            />
                          </div>
                          <div className="flex items-center gap-2">
                            <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                              <span>Ea<sub style={{fontSize: '0.75em', lineHeight: '1'}}>backward</sub>:</span>
                            </Label>
                            <div className="flex items-center gap-1 flex-1">
                              <Input
                                type="number"
                                value={reaction.EaValueBackward || ''}
                                onChange={(e) => updateReaction(reaction.id, 'EaValueBackward', e.target.value)}
                                placeholder="60000"
                              />
                              <span className="text-xs text-muted-foreground">J/mol</span>
                            </div>
                          </div>
                        </div>
                      )}
                    </div>
                    
                    {index < reactions.length - 1 && <hr className="my-4" />}
                  </div>
                ))}
              </CardContent>
            </Card>

            {/* Component Feed Composition */}
            <Card>
              <CardHeader>
                <CardTitle className="font-bold">Components</CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="mb-2">
                  <div className={`grid gap-2 items-center text-xs font-medium text-muted-foreground`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                    <div></div>
                    <div></div>
                    <div className="text-center flex justify-center" style={{gridColumn: `span ${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}`}}>
                      <Label className="text-sm font-medium">Reaction Order</Label>
                    </div>
                  </div>
                </div>
                <div className={`grid gap-2 items-center text-xs font-medium text-muted-foreground border-b pb-1`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                  <div className="text-center flex justify-center">
                    <Label className="text-sm font-medium">Name</Label>
                  </div>
                  <div className="text-center flex justify-center">
                    <Label className="text-sm font-medium">Initial Flow Rate (mol/s)</Label>
                  </div>
                  {reactions.map((reaction, index) => (
                    reaction.isEquilibrium ? (
                      <React.Fragment key={reaction.id}>
                        <div className="text-center flex justify-center">
                          <Label className="text-xs font-medium">Rxn {index + 1} Fwd</Label>
                        </div>
                        <div className="text-center flex justify-center">
                          <Label className="text-xs font-medium">Rxn {index + 1} Rev</Label>
                        </div>
                      </React.Fragment>
                    ) : (
                      <div key={reaction.id} className="text-center flex justify-center">
                        <Label className="text-xs font-medium">Rxn {index + 1}</Label>
                      </div>
                    )
                  ))}
                </div>
                {components.map((comp, index) => {
                  return (
                    <div key={comp.id} className={`grid gap-2 items-center`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                      <div className="flex justify-center">
                        <div className="text-sm font-medium p-2">
                          {comp.name}
                        </div>
                      </div>
                      <div className="flex justify-center">
                        <Input
                          type="number"
                          value={comp.initialFlowRate}
                          onChange={(e) => handleComponentChange(comp.id, 'initialFlowRate', e.target.value)}
                          onKeyDown={handleKeyDown}
                          placeholder="1.0"
                          className="w-full text-center"
                        />
                      </div>
                      {reactions.map((reaction, reactionIndex) => {
                        const isReactantInReaction = reaction.reactants && reaction.reactants.toLowerCase().includes(comp.name.toLowerCase());
                        const isProductInReaction = reaction.products && reaction.products.toLowerCase().includes(comp.name.toLowerCase());
                        
                        if (reaction.isEquilibrium) {
                          const forwardKey = `${reaction.id}-fwd`;
                          const backwardKey = `${reaction.id}-rev`;
                          
                          return (
                            <React.Fragment key={reaction.id}>
                              {/* Forward reaction order */}
                              <div className="flex items-center justify-center">
                                {isReactantInReaction ? (
                                  <Input
                                    type="number"
                                    value={comp.reactionOrders?.[forwardKey] || ''}
                                    onChange={(e) => {
                                      const newReactionOrders = { ...comp.reactionOrders, [forwardKey]: e.target.value };
                                      handleComponentChange(comp.id, 'reactionOrders', newReactionOrders);
                                    }}
                                    placeholder="0"
                                    className="w-16 h-8 text-center text-xs"
                                    step="0.1"
                                    min="0"
                                  />
                                ) : (
                                  <Input
                                    type="text"
                                    value="-"
                                    disabled
                                    className="w-16 h-8 text-center text-xs bg-muted text-muted-foreground cursor-not-allowed"
                                  />
                                )}
                              </div>
                              {/* Reverse reaction order */}
                              <div className="flex items-center justify-center">
                                {isProductInReaction ? (
                                  <Input
                                    type="number"
                                    value={comp.reactionOrders?.[backwardKey] || ''}
                                    onChange={(e) => {
                                      const newReactionOrders = { ...comp.reactionOrders, [backwardKey]: e.target.value };
                                      handleComponentChange(comp.id, 'reactionOrders', newReactionOrders);
                                    }}
                                    placeholder="0"
                                    className="w-16 h-8 text-center text-xs"
                                    step="0.1"
                                    min="0"
                                  />
                                ) : (
                                  <Input
                                    type="text"
                                    value="-"
                                    disabled
                                    className="w-16 h-8 text-center text-xs bg-muted text-muted-foreground cursor-not-allowed"
                                  />
                                )}
                              </div>
                            </React.Fragment>
                          );
                        } else {
                          return (
                            <div key={reaction.id} className="flex items-center justify-center">
                              {isReactantInReaction ? (
                                <Input
                                  type="number"
                                  value={comp.reactionOrders?.[reaction.id] || ''}
                                  onChange={(e) => {
                                    const newReactionOrders = { ...comp.reactionOrders, [reaction.id]: e.target.value };
                                    handleComponentChange(comp.id, 'reactionOrders', newReactionOrders);
                                  }}
                                  placeholder="0"
                                  className="w-16 h-8 text-center text-xs"
                                  step="0.1"
                                  min="0"
                                />
                              ) : (
                                <Input
                                  type="text"
                                  value="-"
                                  disabled
                                  className="w-16 h-8 text-center text-xs bg-muted text-muted-foreground cursor-not-allowed"
                                />
                              )}
                            </div>
                          );
                        }
                      })}
                    </div>
                  );
                })}
              </CardContent>
            </Card>
          </div>

          {/* Right Column: Visualization and Results */}
          <div className="space-y-6 md:col-span-3"> {/* Added md:col-span-3 */}
            {!showGraph ? (
              <>
                {/* Reactor Visualization and Type Selection */}
                <Card>
              <CardContent className="p-0">
                <div className="flex justify-center pt-6">
                  <div className="flex items-center gap-2 bg-muted rounded-lg p-1">
                    <Button
                      variant={reactorType === 'CSTR' ? 'default' : 'ghost'}
                      size="sm"
                      onClick={() => setReactorType('CSTR')}
                      className="h-8 px-3"
                    >
                      CSTR
                    </Button>
                    <Button
                      variant={reactorType === 'PFR' ? 'default' : 'ghost'}
                      size="sm"
                      onClick={() => setReactorType('PFR')}
                      className="h-8 px-3"
                    >
                      PFR
                    </Button>
                  </div>
                </div>
                <div className="p-6 pb-2">
                  {renderReactorSVG()}
                </div>
              </CardContent>
            </Card>

            {/* Results Display */}
            {(calculationResults || recycleResults) && (
              <Card>
                <CardHeader>
                  <div className="flex items-center justify-between">
                    <CardTitle>Results {operationMode === 'Recycle' ? '- Recycle System' : '- Single Pass'}</CardTitle>
                    <div className="flex items-center gap-2 text-xs text-muted-foreground">
                      <div className="flex items-center gap-1">
                        <div className="px-2 py-1 bg-yellow-100 dark:bg-yellow-900/50 rounded-sm text-[10px] font-medium text-gray-700 dark:text-gray-300">
                          Limiting Reactant
                        </div>
                      </div>
                    </div>
                  </div>
                </CardHeader>
                <CardContent>
                  <div className="overflow-x-auto">
                    {operationMode === 'SinglePass' && calculationResults ? (
                      // Single Pass Results Table
                      (<table className="w-full border-collapse border border-border">
                        <thead>
                          <tr className="bg-muted/50">
                            <th className="border border-border px-3 py-2 text-center font-medium">Comp.</th>
                            <th className="border border-border px-3 py-2 text-center font-medium">Conversion (%)</th>
                            <th className="border border-border px-3 py-2 text-center font-medium">Selectivity (%)</th>
                            <th className="border border-border px-3 py-2 text-center font-medium">Outlet Flow Rate (mol/s)</th>
                            <th className="border border-border px-3 py-2 text-center font-medium">Composition (mol%)</th>
                          </tr>
                        </thead>
                        <tbody>
                          {components.map((comp) => {
                            const isReactantInAnyReaction = reactions.some(r =>
                              (r.reactants || '').split(/[+\s]+/).map(s => s.trim().replace(/^\d+/, '')).includes(comp.name)
                            );
                            const isProductInAnyReaction = reactions.some(r =>
                              (r.products || '').split(/[+\s]+/).map(s => s.trim().replace(/^\d+/, '')).includes(comp.name)
                            );

                            let conversionValue: string | number = '-';
                            if (isReactantInAnyReaction && calculationResults?.outletFlowRates) {
                              const F_R_in = parseFloat(comp.initialFlowRate || '0');
                              const F_R_out = calculationResults.outletFlowRates?.[comp.name] ?? 0;
                              if (F_R_in > 1e-9) {
                                const X_R = (F_R_in - F_R_out) / F_R_in;
                                conversionValue = formatToSigFigs(Math.max(0, Math.min(1, X_R)) * 100);
                              } else {
                                conversionValue = F_R_out > 1e-9 ? '0.00' : '-';
                              }
                            }

                            let selectivityValue: string | number = '-';
                            if (isProductInAnyReaction && calculationResults?.conversion && calculationResults?.outletFlowRates) {
                              const F_P_out = calculationResults.outletFlowRates?.[comp.name] ?? 0;
                              const F_P_in = parseFloat(comp.initialFlowRate || '0');
                              const moles_P_formed = Math.max(0, F_P_out - F_P_in);

                              const keyReactantName = calculationResults.conversion.reactantName;
                              const F_key_in_comp = components.find(c => c.name === keyReactantName);
                              const F_key_in = parseFloat(F_key_in_comp?.initialFlowRate || '0');
                              const X_key = calculationResults.conversion.value;
                              const moles_key_reacted = F_key_in * X_key;

                              if (moles_key_reacted > 1e-9) {
                                const rawSelectivity = moles_P_formed / moles_key_reacted;
                                selectivityValue = formatToSigFigs(Math.max(0, rawSelectivity) * 100);
                              } else {
                                selectivityValue = (moles_P_formed > 1e-9) ? '-' : '0.00';
                              }
                            }
                            
                            const isLimitingReactant = comp.name === calculationResults?.conversion?.reactantName;

                            // Calculate composition (mol%)
                            let compositionValue: string = '-';
                            if (calculationResults?.outletFlowRates) {
                              const totalOutletFlow = Object.values(calculationResults.outletFlowRates).reduce((sum, flow) => sum + flow, 0);
                              const compOutletFlow = calculationResults.outletFlowRates[comp.name] ?? 0;
                              if (totalOutletFlow > 1e-9) {
                                const molPercent = (compOutletFlow / totalOutletFlow) * 100;
                                compositionValue = formatToSigFigs(molPercent);
                              }
                            }

                            return (
                              <tr 
                                key={comp.id} 
                                className={`${isLimitingReactant ? 'bg-yellow-100 dark:bg-yellow-900/30' : ''}`}
                              >
                                <td className={`border border-border px-3 py-2 text-center ${isLimitingReactant ? 'font-bold' : ''}`}>
                                  {comp.name}
                                </td>
                                <td className="border border-border px-3 py-2 text-center">
                                  {conversionValue}
                                </td>
                                <td className="border border-border px-3 py-2 text-center">
                                  {selectivityValue}
                                </td>
                                <td className="border border-border px-3 py-2 text-center">
                                  {calculationResults?.outletFlowRates?.[comp.name] !== undefined ? formatToSigFigs(calculationResults.outletFlowRates[comp.name]) : '-'}
                                </td>
                                <td className="border border-border px-3 py-2 text-center">
                                  {compositionValue}
                                </td>
                              </tr>
                            );
                          })}
                        </tbody>
                      </table>)
                    ) : operationMode === 'Recycle' && recycleResults ? (
                      // Recycle Results Table
                      (<div className="space-y-4">
                        {/* Summary Information */}
                        <div className="grid grid-cols-2 gap-4 p-4 bg-muted/50 rounded-lg">
                          <div className="text-center">
                            <div className="text-sm font-medium text-muted-foreground">Overall Conv.</div>
                            <div className="text-lg font-bold">
                              {recycleResults.overallConversion ? `${formatToSigFigs(recycleResults.overallConversion.value * 100)}%` : '-'}
                            </div>
                            <div className="text-xs text-muted-foreground">
                              {recycleResults.overallConversion?.reactantName || ''}
                            </div>
                          </div>
                          <div className="text-center">
                            <div className="text-sm font-medium text-muted-foreground">Single Pass Conversion</div>
                            <div className="text-lg font-bold">
                              {recycleResults.singlePassConversion ? `${formatToSigFigs(recycleResults.singlePassConversion.value * 100)}%` : '-'}
                            </div>
                            <div className="text-xs text-muted-foreground">
                              {recycleResults.singlePassConversion?.reactantName || ''} (Reactor Only)
                            </div>
                          </div>
                        </div>
                        {/* Convergence Issues Warning */}
                        {!recycleResults.converged && (
                          <div className="p-4 bg-red-50 dark:bg-red-900/20 border border-red-200 dark:border-red-800 rounded-lg">
                            <h4 className="font-semibold text-red-800 dark:text-red-200 mb-2">⚠️ System Not Converged</h4>
                            <div className="space-y-2 text-sm text-red-700 dark:text-red-300">
                              {recycleResults.singlePassConversion && recycleResults.singlePassConversion.value < 0.05 && (
                                <div className="p-2 bg-blue-100 dark:bg-blue-900/50 rounded">
                                  <div className="font-medium text-blue-800 dark:text-blue-200">Low Single Pass Conversion Detected!</div>
                                  <div className="text-blue-700 dark:text-blue-300 text-sm">
                                    Single pass conversion is only {formatToSigFigs(recycleResults.singlePassConversion.value * 100)}%. 
                                    Try:
                                    <br />• Increasing reactor volume from {reactorVolume}m³ to {Math.ceil(parseFloat(reactorVolume) * 5)}m³+
                                    <br />• Increasing temperature from {kinetics.reactionTempK}K to {Math.ceil(parseFloat(kinetics.reactionTempK) + 50)}K+
                                  </div>
                                </div>
                              )}
                              <div className="mt-2">
                                <div className="font-medium">General Solutions:</div>
                                <div className="ml-4">- Increase reactor volume to allow more reaction</div>
                                <div className="ml-4">- Increase temperature to increase reaction rate</div>
                                <div className="ml-4">- Check that reaction kinetics parameters are reasonable</div>
                              </div>
                            </div>
                          </div>
                        )}
                        {/* Detailed Results Table */}
                        <table className="w-full border-collapse border border-border">
                          <thead>
                            <tr className="bg-muted/50">
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Component</th>
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Overall Conv (%)</th>
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Sel (%)</th>
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Outlet (mol/s)</th>
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Outlet Comp (%)</th>
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Rec (mol/s)</th>
                              <th className="border border-border px-2 py-2 text-center font-medium text-xs">Rec Comp (%)</th>
                            </tr>
                          </thead>
                          <tbody>
                            {components.map((comp) => {
                              const isReactantInAnyReaction = reactions.some(r =>
                                (r.reactants || '').split(/[+\s]+/).map(s => s.trim().replace(/^\d+/, '')).includes(comp.name)
                              );
                              const isProductInAnyReaction = reactions.some(r =>
                                (r.products || '').split(/[+\s]+/).map(s => s.trim().replace(/^\d+/, '')).includes(comp.name)
                             
                              );

                              // Overall conversion (fresh feed to product stream)
                              let conversionValue: string = '-';
                              if (isReactantInAnyReaction) {
                                const F_fresh = parseFloat(comp.initialFlowRate || '0');
                                const F_product = recycleResults.productStream[comp.name] || 0;
                                if (F_fresh > 1e-9) {
                                  const X_overall = (F_fresh - F_product) / F_fresh;
                                  conversionValue = formatToSigFigs(Math.max(0, Math.min(1, X_overall)) * 100);
                                } else {
                                  conversionValue = F_product > 1e-9 ? '0.00' : '-';
                                }
                              }

                              // Overall selectivity (product formed / reactant consumed from fresh feed)
                              let selectivityValue: string = '-';
                              if (isProductInAnyReaction && recycleResults.overallConversion) {
                                const F_product_out = recycleResults.productStream[comp.name] || 0;
                                const F_product_in = parseFloat(comp.initialFlowRate || '0'); // Fresh feed of this product
                                const moles_product_formed = Math.max(0, F_product_out - F_product_in);

                                const keyReactantName = recycleResults.overallConversion.reactantName;
                                const F_key_fresh_comp = components.find(c => c.name === keyReactantName);
                                const F_key_fresh = parseFloat(F_key_fresh_comp?.initialFlowRate || '0');
                                const X_key_overall = recycleResults.overallConversion.value;
                                const moles_key_consumed = F_key_fresh * X_key_overall;

                                if (moles_key_consumed > 1e-9) {
                                  const overallSelectivity = moles_product_formed / moles_key_consumed;
                                  selectivityValue = formatToSigFigs(Math.max(0, overallSelectivity) * 100);
                                } else {
                                  selectivityValue = (moles_product_formed > 1e-9) ? '-' : '0.00';
                                }
                              }

                              const isLimitingReactant = comp.name === recycleResults.overallConversion?.reactantName;

                              // Calculate outlet composition (mol%)
                              let outletCompositionValue: string = '-';
                              const totalOutletFlow = Object.values(recycleResults.reactorOutletStreamFinal).reduce((sum, flow) => sum + flow, 0);
                              const compOutletFlow = recycleResults.reactorOutletStreamFinal[comp.name] || 0;
                              if (totalOutletFlow > 1e-9) {
                                const molPercent = (compOutletFlow / totalOutletFlow) * 100;
                                outletCompositionValue = formatToSigFigs(molPercent);
                              }

                              // Calculate recycle composition (mol%)
                              let recycleCompositionValue: string = '-';
                              const totalRecycleFlow = Object.values(recycleResults.recycleStreamFinal).reduce((sum, flow) => sum + flow, 0);
                              const compRecycleFlow = recycleResults.recycleStreamFinal[comp.name] || 0;
                              if (totalRecycleFlow > 1e-9) {
                                const molPercent = (compRecycleFlow / totalRecycleFlow) * 100;
                                recycleCompositionValue = formatToSigFigs(molPercent);
                              }

                              return (
                                <tr 
                                  key={comp.id} 
                                  className={`${isLimitingReactant ? 'bg-yellow-100 dark:bg-yellow-900/30' : ''}`}
                                >
                                  <td className={`border border-border px-2 py-2 text-center text-xs ${isLimitingReactant ? 'font-bold' : ''}`}>
                                    {comp.name}
                                  </td>
                                  <td className="border border-border px-2 py-2 text-center text-xs">
                                    {conversionValue}
                                  </td>
                                  <td className="border border-border px-2 py-2 text-center text-xs">
                                    {selectivityValue}
                                  </td>
                                  <td className="border border-border px-2 py-2 text-center text-xs">
                                    {formatToSigFigs(recycleResults.reactorOutletStreamFinal[comp.name] || 0)}
                                  </td>
                                  <td className="border border-border px-2 py-2 text-center text-xs">
                                    {outletCompositionValue}
                                  </td>
                                  <td className="border border-border px-2 py-2 text-center text-xs">
                                    {formatToSigFigs(recycleResults.recycleStreamFinal[comp.name] || 0)}
                                  </td>
                                  <td className="border border-border px-2 py-2 text-center text-xs">
                                    {recycleCompositionValue}
                                  </td>
                                </tr>
                              );
                            })}
                          </tbody>
                        </table>
                        {/* System Performance Analysis - Moved below the table */}
                        <div className="p-4 bg-gray-50 dark:bg-gray-900/20 border border-gray-200 dark:border-gray-800 rounded-lg">
                          <h4 className="font-semibold text-gray-800 dark:text-gray-200 mb-2">System Performance Analysis</h4>
                          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
                            <div>
                              <div className="text-muted-foreground">Reactor Inlet Total Flow</div>
                              <div className="font-medium">
                                {formatToSigFigs(Object.values(recycleResults.reactorInletStreamFinal).reduce((sum, flow) => sum + flow, 0))} mol/s
                              </div>
                            </div>
                            <div>
                              <div className="text-muted-foreground">Reactor Outlet Total Flow</div>
                              <div className="font-medium">
                                {formatToSigFigs(Object.values(recycleResults.reactorOutletStreamFinal).reduce((sum, flow) => sum + flow, 0))} mol/s
                              </div>
                            </div>
                            <div>
                              <div className="text-muted-foreground">Recycle Ratio</div>
                              <div className="font-medium">
                                {(() => {
                                  const totalRecycle = Object.values(recycleResults.recycleStreamFinal).reduce((sum, flow) => sum + flow, 0);
                                  const totalFresh = components.reduce((sum, comp) => sum + (parseFloat(comp.initialFlowRate) || 0), 0);
                                  return totalFresh > 1e-9 ? formatToSigFigs(totalRecycle / totalFresh) : '-';
                                })()}
                              </div>
                            </div>
                            <div>
                              <div className="text-muted-foreground">Material Balance</div>
                              <div className="font-medium">
                                {(() => {
                                  const totalIn = components.reduce((sum, comp) => sum + (parseFloat(comp.initialFlowRate) || 0), 0);
                                  const totalOut = Object.values(recycleResults.productStream).reduce((sum, flow) => sum + flow, 0);
                                  const balance = totalIn > 1e-9 ? (totalOut / totalIn) * 100 : 0;
                                  return `${formatToSigFigs(balance)}%`;
                                })()}
                              </div>
                            </div>
                          </div>
                        </div>
                      </div>)
                    ) : null}
                  </div>
                </CardContent>
              </Card>
            )}
              </>
            ) : (
              /* Graph Display - Replaces visualization and results */
              (<Card>
                <CardContent>
                  <div className="relative w-full aspect-square rounded-md">
                    {/* PFR/CSTR Switch on top-left corner of graph */}
                    <div className="absolute top-4 left-4 z-10">
                      <div className="flex items-center gap-2 bg-muted rounded-lg p-1">
                        <Button
                          variant={reactorType === 'PFR' ? 'default' : 'ghost'}
                          size="sm"
                          onClick={() => setReactorType('PFR')}
                          className="h-8 px-3"
                        >
                          PFR
                        </Button>
                        <Button
                          variant={reactorType === 'CSTR' ? 'default' : 'ghost'}
                          size="sm"
                          onClick={() => setReactorType('CSTR')}
                          className="h-8 px-3"
                        >
                          CSTR
                        </Button>
                      </div>
                    </div>
                    
                    {/* Graph Type Selector on top-right corner of graph */}
                    <div className="absolute top-4 right-4 z-10">
                      <Select value={graphType} onValueChange={(value: 'conversion' | 'selectivity' | 'flowrates' | 'composition') => setGraphType(value)}>
                        <SelectTrigger className="w-56">
                          <SelectValue placeholder="Select graph type" />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="conversion">Conversion vs Volume</SelectItem>
                          <SelectItem value="selectivity">Selectivity vs Volume</SelectItem>
                          <SelectItem value="flowrates">Flow Rates vs Volume</SelectItem>
                          <SelectItem value="composition">Composition vs Volume</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                    {Object.keys(graphOptions).length > 0 && graphOptions.series && (graphOptions.series as any[]).length > 0 ? (
                      <ReactECharts 
                        ref={echartsRef} 
                        echarts={echarts} 
                        option={graphOptions} 
                        style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }} 
                        notMerge={true} 
                        lazyUpdate={false}
                      />
                    ) : (
                      <div className="absolute inset-0 flex items-center justify-center text-white">
                        <div className="text-center">
                          <div className="mb-2">Generating analysis data...</div>
                          <div className="text-sm text-gray-300">This may take a moment</div>
                        </div>
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>)
            )}

            {/* Graph Visualization Button */}
            {(calculationResults || recycleResults) && operationMode === 'SinglePass' && (
              <div className="flex justify-center">
                <Button 
                  onClick={() => setShowGraph(!showGraph)}
                  variant="outline"
                  className="flex items-center gap-2"
                >
                  <BarChart3 className="h-4 w-4" />
                  {showGraph ? 'Hide Analysis Graph' : 'Show Analysis Graph'}
                </Button>
              </div>
            )}

            {/* Error Display */}
            {calculationError && (
              <Alert variant="destructive">
                <AlertCircle className="h-4 w-4" />
                <AlertTitle>Error</AlertTitle>
                <AlertDescription>{calculationError}</AlertDescription>
              </Alert>
            )}
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}