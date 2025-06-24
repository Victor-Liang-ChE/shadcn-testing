'use client';

import React, { useState, useMemo, useCallback, memo, useRef, useEffect } from 'react';
import { useTheme } from "next-themes";

// Import ECharts components
import ReactECharts from 'echarts-for-react';
import * as echarts from 'echarts/core';
// Import specific types from the main echarts package
import type { EChartsOption, SeriesOption } from 'echarts';
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

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Info } from 'lucide-react';

// --- Type Definitions ---
type TuningMethod = 'ziegler-nichols' | 'itae' | 'amigo' | 'imc';
type ControllerType = 'P' | 'PI' | 'PID';
type ItaeInputType = 'disturbance' | 'setpoint';
type ImcModelCase = 'A' | 'B' | 'C' | 'D' | 'E' | 'F' | 'G' | 'H' | 'I' | 'J' | 'K' | 'L' | 'M' | 'N' | 'O'; // All IMC model cases

interface TuningResult {
  kc: number | null;
  tauI: number | null;
  tauD: number | null;
}

// --- Helper Data for Tuning Rules ---

// For ITAE Rule: Y = A * (theta/tau)^B
const itaeParams = {
  disturbance: {
    PI: { P: { A: 0.859, B: -0.977 }, I: { A: 0.674, B: -0.680 } },
    PID: { P: { A: 1.357, B: -0.947 }, I: { A: 0.842, B: -0.738 }, D: { A: 0.381, B: 0.995 } },
  },
  setpoint: {
    // Note: Integral mode has a different formula for setpoint changes
    PI: { P: { A: 0.586, B: -0.916 }, I: { A: 1.03, B: -0.165 } },
    PID: { P: { A: 0.965, B: -0.85 }, I: { A: 0.796, B: -0.1465 }, D: { A: 0.308, B: 0.929 } },
  },
};

// --- SliderGroup Component Definition ---
interface SliderGroupProps {
  label: React.ReactNode;
  tooltip: string;
  value: string;
  onChange: (v: string) => void;
  min: number;
  max: string;
  step: number;
  unit?: string;
  maxValue: string;
  onMaxChange: (v: string) => void;
}

const SliderGroup = ({ label, tooltip, value, onChange, min, max, step, unit, maxValue, onMaxChange }: SliderGroupProps) => {
  return (
    <div className="space-y-2">
      {/* Label and value above the slider */}
      <div className="flex items-center gap-2">
        <Tooltip>
          <TooltipTrigger asChild>
            <Label className="text-sm font-medium cursor-help">
              {label}
            </Label>
          </TooltipTrigger>
          <TooltipContent><p>{tooltip}</p></TooltipContent>
        </Tooltip>
        <div className="flex items-center gap-1">
          <span className="text-sm font-medium">{value}</span>
          {unit && <span className="text-xs text-muted-foreground">{unit}</span>}
        </div>
      </div>
      {/* Slider and max input on the same row */}
      <div className="flex items-center gap-3">
        <Slider
          min={min}
          max={parseFloat(max)}
          step={step}
          value={[parseFloat(value)]}
          onValueChange={([v]) => {
            onChange(v.toString());
          }}
          className="flex-1"
        />
        <div className="flex items-center gap-1">
          <Label className="text-xs text-muted-foreground">Max:</Label>
          <Input
            type="number"
            value={maxValue}
            onChange={e => onMaxChange(e.target.value)}
            className="w-20 h-8 text-xs"
            min={min.toString()}
            step={step.toString()}
          />
        </div>
      </div>
    </div>
  );
};

// --- Main Component ---
export default function PidTuningPage() {
  const { resolvedTheme } = useTheme(); // Get the resolved theme ('light' or 'dark')
  const purpleColor = '#9333ea'; // A nice purple for the controller output
  
  // Input States
  const [tuningMethod, setTuningMethod] = useState<TuningMethod>('imc');
  const [controllerType, setControllerType] = useState<ControllerType>('PI');
  
  // Model Parameter String States
  const [k, setK] = useState('1.0');
  const [tau, setTau] = useState('1.0');
  const [tau2, setTau2] = useState('5.0'); // For Case B
  const [tau3, setTau3] = useState('2.0'); // For Case I
  const [beta, setBeta] = useState('1.0'); // For Case D
  const [theta, setTheta] = useState('0.1');
  const [kcu, setKcu] = useState('2.0');
  const [pu, setPu] = useState('5.0');
  const [tauC, setTauC] = useState('2.0'); // For IMC method
  const [zeta, setZeta] = useState('0.5'); // For SOPTD model

  // Specific Method States
  const [itaeInputType, setItaeInputType] = useState<ItaeInputType>('disturbance');
  const [imcModelCase, setImcModelCase] = useState<ImcModelCase>('A');

  // Maximum value states for sliders
  const [maxK, setMaxK] = useState<string>('10');
  const [maxTau, setMaxTau] = useState<string>('10');
  const [maxTau2, setMaxTau2] = useState<string>('25');
  const [maxTau3, setMaxTau3] = useState<string>('10');
  const [maxBeta, setMaxBeta] = useState<string>('5');
  const [maxTheta, setMaxTheta] = useState<string>('20');
  const [maxKcu, setMaxKcu] = useState<string>('10');
  const [maxPu, setMaxPu] = useState<string>('25');
  const [maxTauC, setMaxTauC] = useState<string>('10');
  const [maxZeta, setMaxZeta] = useState<string>('2');

  // Control & Result States
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<TuningResult | null>(null);
  const [displayedMethod, setDisplayedMethod] = useState<string>('');
  
  // Simulation Control States
  const [currentSetpoint, setCurrentSetpoint] = useState<number>(1);
  const [simulationTime, setSimulationTime] = useState<number>(0);
  const [isSimulationRunning, setIsSimulationRunning] = useState<boolean>(false);
  const [currentPV, setCurrentPV] = useState<number>(0);
  const [controlError, setControlError] = useState<number>(0);
  const [divergenceWarning, setDivergenceWarning] = useState<boolean>(false);
  const [resetCountdown, setResetCountdown] = useState<number>(0);
  const [showErrorGraph, setShowErrorGraph] = useState<boolean>(false);
  const [graphMode, setGraphMode] = useState<'process' | 'error' | 'bode' | 'nyquist' | 'power'>('process');
  const [graphInitialized, setGraphInitialized] = useState<boolean>(false);
  
  // Add a ref for immediate setpoint updates that the simulation can access
  const currentSetpointRef = useRef<number>(1);
  // Add a ref for the error graph state that the simulation can access
  const showErrorGraphRef = useRef<boolean>(false);
  const graphModeRef = useRef<'process' | 'error' | 'bode' | 'nyquist' | 'power'>('process');
  // Add a ref for the PID result that the simulation can access
  const resultRef = useRef<TuningResult | null>(null);
  // Add refs for warning states to prevent graph flickering
  const divergenceWarningRef = useRef<boolean>(false);
  const resetCountdownRef = useRef<number>(0);

  // Add these states and refs for the new multi-view functionality
  const [isMultiView, setIsMultiView] = useState<boolean>(false);
  const [multiViewData, setMultiViewData] = useState<{ time: number; pv: number; sp: number; output: number }[]>([]);
  const isMultiViewRef = useRef<boolean>(false);

  // Add state for the new overlapped graph view
  const [isOverlappedView, setIsOverlappedView] = useState<boolean>(true);
  const isOverlappedViewRef = useRef<boolean>(false);
  useEffect(() => {
    isOverlappedViewRef.current = isOverlappedView;
  }, [isOverlappedView]);

  const multiViewProcessRef = useRef<ReactECharts | null>(null);
  const multiViewErrorRef = useRef<ReactECharts | null>(null);
  const multiViewOutputRef = useRef<ReactECharts | null>(null);
  const overlappedChartRef = useRef<ReactECharts | null>(null); // Ref for the new overlapped chart

  const echartsRef = useRef<ReactECharts | null>(null); // Ref for ECharts instance
  const simulationIntervalRef = useRef<NodeJS.Timeout | null>(null);
  const countdownIntervalRef = useRef<NodeJS.Timeout | null>(null); // Ref for countdown timer
  
  // Add simulation state refs for real-time simulation
  const simulationStateRef = useRef({
    pv: 0,
    prev_pv: 0,
    integral_error: 0,
    filtered_derivative: 0,
    y1: 0, // For second-order models
    y2: 0, // For second-order models
    delayBuffer: [] as number[],
  });

  const simulationDataRef = useRef<{ time: number; pv: number; sp: number; output: number }[]>([]);

  const multiViewProcessOptions = useMemo(() => {
    return getMultiViewProcessOptions(resolvedTheme === 'dark', []);
  }, [resolvedTheme]);

  const multiViewErrorOptions = useMemo(() => {
    return getMultiViewErrorOptions(resolvedTheme === 'dark', []);
  }, [resolvedTheme]);

  const multiViewOutputOptions = useMemo(() => {
    return getMultiViewOutputOptions(resolvedTheme === 'dark', []);
  }, [resolvedTheme]);

  function getOverlappedGraphOptions(isDark: boolean): EChartsOption {
    const textColor = isDark ? 'white' : '#000000';
    const tooltipBg = isDark ? '#08306b' : '#f8f9fa';
    const tooltipBorder = isDark ? '#55aaff' : '#dee2e6';

    return {
      backgroundColor: 'transparent',
      title: { text: 'Combined Simulation View', textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' }, left: 'center', top: '2%' },
      grid: { left: '8%', right: '8%', bottom: '15%', top: '12%', containLabel: true },
      tooltip: {
        trigger: 'axis',
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        borderWidth: 1,
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        axisPointer: { type: 'cross', label: { backgroundColor: tooltipBg, color: textColor, borderColor: tooltipBorder, fontFamily: 'Merriweather Sans' } },
        formatter: function(params: any) {
          if (!params || !Array.isArray(params) || params.length === 0) return '';
          const time = params[0]?.value[0];
          if (time === undefined) return '';

          let tooltip = `Time: ${parseFloat(time).toPrecision(3)} s<br/>`;
          const seriesMap = new Map<string, { value: number, color: string }>();

          params.forEach((param: any) => {
            seriesMap.set(param.seriesName, {
              value: parseFloat(param.value[1]),
              color: param.color
            });
          });

          const formatLine = (name: string, color: string) => {
            if (seriesMap.has(name)) {
              const series = seriesMap.get(name)!;
              tooltip += `<span style="color: ${color};">${name}: ${series.value.toPrecision(3)}</span><br/>`;
            }
          };

          formatLine('Process Variable', '#00ff51ff');
          formatLine('Setpoint', '#ef4444');
          formatLine('Control Error', '#3b82f6');
          formatLine('Controller Output', purpleColor);

          return tooltip;
        }
      },
      legend: {
        data: [
          { name: 'Process Variable', icon: 'rect', itemStyle: { color: '#00ff51ff' } },
          { name: 'Setpoint', icon: 'rect', itemStyle: { color: '#ef4444' } },
          { name: 'Control Error', icon: 'rect', itemStyle: { color: '#3b82f6' } },
          { name: 'Controller Output', icon: 'rect', itemStyle: { color: purpleColor } }
        ],
        bottom: '2%',
        left: 'center',
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        itemWidth: 25,
        itemHeight: 2
      },
      xAxis: {
        type: 'value',
        name: 'Time (s)',
        min: 0,
        max: 30,
        nameLocation: 'middle',
        nameGap: 25,
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { 
          color: textColor, 
          fontSize: 16, 
          fontFamily: 'Merriweather Sans',
          showMinLabel: false,
          showMaxLabel: false
        },
        splitLine: { show: false } 
      },
      yAxis: {
        type: 'value',
        name: 'Value',
        min: -0.5,
        max: 2.5,
        nameLocation: 'middle',
        nameGap: 50,
        axisLine: { lineStyle: { color: textColor } },
        axisLabel: {
          color: textColor,
          fontSize: 16,
          fontFamily: 'Merriweather Sans'
        },
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false }
      },
      series: [
        { name: 'Process Variable', type: 'line', data: [], showSymbol: false, lineStyle: { color: '#00ff51ff', width: 3 } },
        { name: 'Setpoint', type: 'line', data: [], showSymbol: false, lineStyle: { color: '#ef4444', type: 'dashed', width: 3 } },
        { name: 'Control Error', type: 'line', data: [], showSymbol: false, lineStyle: { color: '#3b82f6', width: 3 } },
        { name: 'Controller Output', type: 'line', data: [], showSymbol: false, lineStyle: { color: purpleColor, width: 3 } }
      ]
    };
  }

  const overlappedGraphOptions = useMemo(() => {
    return getOverlappedGraphOptions(resolvedTheme === 'dark');
  }, [resolvedTheme]);

  // Single-step simulation function for real-time updates
  const runSingleSimulationStep = (
    time: number,
    dt: number,
    processParams: { K: number; tau: number; theta: number; tau2?: number; zeta?: number },
    controller: { Kc: number; tauI: number | null; tauD: number | null },
    modelCase: ImcModelCase,
    setpointValue: number
  ) => {
    const { K, tau, theta, tau2, zeta } = processParams;
    const { Kc, tauI, tauD } = controller;

    // Get current state from the ref
    const state = simulationStateRef.current;
    let { pv, prev_pv, integral_error, filtered_derivative, y1, y2, delayBuffer } = state;

    // --- Derivative Filter ---
    const alpha = 0.1;
    const T_filter = (tauD && tauD > 0) ? Math.max(alpha * tauD, dt) : dt;

    const error = setpointValue - pv;

    // --- PID Controller with Derivative Filter ---
    const p_term = Kc * error;
    if (tauI && tauI > 0) {
      integral_error += (Kc / tauI) * error * dt;
    }

    const raw_derivative = (prev_pv - pv) / dt;
    filtered_derivative += (dt / T_filter) * (raw_derivative - filtered_derivative);
    const derivative_term = (tauD && tauD > 0) ? Kc * tauD * filtered_derivative : 0;
    
    const controller_output = p_term + integral_error + derivative_term;
    
    // --- Process Simulation (Handles different models) ---
    // If theta is very small (< 0.1), skip delay to minimize startup delay
    let delayed_input;
    if (theta < 0.1) {
      delayed_input = controller_output;
    } else {
      delayed_input = delayBuffer.shift() || 0;
      delayBuffer.push(controller_output);
    }

    // Select the correct process model for simulation
    if (modelCase === 'B' || modelCase === 'K' || modelCase === 'I') {
      const tau1_sim = tau;
      const tau2_sim = tau2 || tau;
      const y1_deriv = (K * delayed_input - y1) / tau1_sim;
      y1 += y1_deriv * dt;
      const y2_deriv = (y1 - y2) / tau2_sim;
      y2 += y2_deriv * dt;
      pv = y2;
    } else {
      const pv_deriv = tau > 0 ? (K * delayed_input - pv) / tau : 0;
      pv += pv_deriv * dt;
    }

    prev_pv = pv;
    
    // Update the state in the ref for the next tick
    simulationStateRef.current = { pv, prev_pv, integral_error, filtered_derivative, y1, y2, delayBuffer };

    // Return the data point for this step including controller output
    return { time, pv, sp: setpointValue, output: controller_output };
  };

  // Simulation function for closed-loop response (kept for compatibility but not used in real-time mode)
  const simulateClosedLoopResponse = (
    processParams: { K: number; tau: number; theta: number; tau2?: number; zeta?: number },
    controller: { Kc: number; tauI: number | null; tauD: number | null },
    modelCase: ImcModelCase,
    setpointValue: number = 1,
    disturbanceTime?: number
  ): { time: number; pv: number; sp: number }[] => {
    
    const { K, tau, theta, tau2, zeta } = processParams;
    const { Kc, tauI, tauD } = controller;

    // Simulation parameters
    const dt = 0.1; // Time step
    const simulationTime = Math.max(tau * 15, theta * 5, 80);
    const n_steps = Math.floor(simulationTime / dt);
    
    // --- Derivative Filter ---
    const alpha = 0.1; // Filter constant (0.1 is a common value)
    const T_filter = (tauD && tauD > 0) ? Math.max(alpha * tauD, dt) : dt;

    // Initial conditions
    let pv = 0;
    let prev_pv = 0;
    let integral_error = 0;
    let filtered_derivative = 0;

    // Second-order model states (if needed)
    let y1 = 0; // Output of first lag
    let y2 = 0; // Output of second lag

    // Time delay buffer
    const delaySteps = Math.floor(theta / dt);
    const delayBuffer = new Array(delaySteps).fill(0);
    
    const results = [{ time: 0, pv: 0, sp: 0 }];

    for (let i = 1; i <= n_steps; i++) {
      const time = i * dt;
      
      // Setpoint logic - step change at t > dt
      let setpoint = time > dt ? setpointValue : 0;
      
      const error = setpoint - pv;

      // --- PID Controller with Derivative Filter ---
      const p_term = Kc * error;
      if (tauI && tauI > 0) {
        integral_error += (Kc / tauI) * error * dt;
      }

      // 1. Calculate the raw derivative on PV
      const raw_derivative = (prev_pv - pv) / dt;
      // 2. Pass it through the low-pass filter
      filtered_derivative += (dt / T_filter) * (raw_derivative - filtered_derivative);
      // 3. The final D-term uses the filtered value
      const derivative_term = (tauD && tauD > 0) ? Kc * tauD * filtered_derivative : 0;
      
      const controller_output = p_term + integral_error + derivative_term;
      
      // --- Process Simulation (Handles different models) ---
      const delayed_input = delayBuffer.shift() || 0;
      delayBuffer.push(controller_output); // No continuous disturbance

      // Select the correct process model for simulation
      if (modelCase === 'B' || modelCase === 'K' || modelCase === 'I') {
        // Second-Order Overdamped (Two Lags in series)
        const tau1_sim = tau;
        const tau2_sim = tau2 || tau; // Fallback if tau2 not provided
        const y1_deriv = (K * delayed_input - y1) / tau1_sim;
        y1 += y1_deriv * dt;
        const y2_deriv = (y1 - y2) / tau2_sim;
        y2 += y2_deriv * dt;
        pv = y2;
      } else {
        // Default to FOPTD for all other relevant cases (G, H, ITAE, AMIGO etc.)
        const pv_deriv = tau > 0 ? (K * delayed_input - pv) / tau : 0;
        pv += pv_deriv * dt;
      }

      prev_pv = pv;
      results.push({ time, pv, sp: setpoint });
    }
    return results;
  };

  // Helper function to get available controller types
  const getAvailableControllerTypes = (method: TuningMethod): ControllerType[] => {
    switch (method) {
      case 'ziegler-nichols':
        return ['P', 'PI', 'PID'];
      case 'itae':
        return ['PI', 'PID'];
      case 'amigo':
        return ['PI', 'PID'];
      case 'imc':
        return ['PI', 'PID']; // Allow PI controllers as some IMC rules result in PI
      default:
        return ['P', 'PI', 'PID'];
    }
  };

  // Helper function to get available IMC cases for controller type
  const getAvailableImcCases = (controllerType: ControllerType): ImcModelCase[] => {
    // These lists are now correct based on our previous discussion
    const piOnlyCases: ImcModelCase[] = ['A', 'E', 'G', 'M'];
    const pidCapableCases: ImcModelCase[] = ['B', 'C', 'D', 'F', 'H', 'I', 'J', 'K', 'L', 'N', 'O'];
    
    if (controllerType === 'PID') {
      return pidCapableCases;
    } else if (controllerType === 'PI') {
      // CHANGE THIS LINE: Instead of returning all models, return only the ones
      // that result in a PI controller according to the table.
      return piOnlyCases;
    }
    
    return []; // No P-only IMC rules
  };

  // Helper function to determine controller type from IMC model case
  const getControllerTypeFromImcCase = (modelCase: ImcModelCase): ControllerType => {
    const piOnlyCases: ImcModelCase[] = ['A', 'E', 'G', 'M'];
    return piOnlyCases.includes(modelCase) ? 'PI' : 'PID';
  };

  // --- Calculation Logic ---
  const handleCalculate = React.useCallback(() => {
    setLoading(true);
    setError(null);
    setResult(null);

    try {
        let calcResult: TuningResult | null = null;
        
        // Parse common parameters
        const numK = parseFloat(k);
        const numTau = parseFloat(tau);
        const numTheta = parseFloat(theta);
        
        switch (tuningMethod) {
            case 'ziegler-nichols':
                const numKcu = parseFloat(kcu);
                const numPu = parseFloat(pu);
                if (isNaN(numKcu) || isNaN(numPu) || numKcu <= 0 || numPu <= 0) {
                    throw new Error("K_cu and P_u must be positive numbers.");
                }
                calcResult = calculateZieglerNichols(controllerType, numKcu, numPu);
                setDisplayedMethod('Ziegler-Nichols Continuous Cycling');
                break;
            
            case 'itae':
                if (isNaN(numK) || isNaN(numTau) || isNaN(numTheta) || numTau <= 0 || numTheta <= 0) {
                    throw new Error("K, τ, and θ must be valid positive numbers.");
                }
                calcResult = calculateItae(controllerType, itaeInputType, numK, numTau, numTheta);
                setDisplayedMethod(`ITAE for ${itaeInputType === 'setpoint' ? 'Set-Point Changes' : 'Disturbance Rejection'}`);
                break;
            
            case 'amigo':
                if (isNaN(numK) || isNaN(numTheta) || numTheta <= 0) {
                    throw new Error("K and θ must be valid positive numbers.");
                }
                const numTauAmigo = isNaN(numTau) ? 0 : numTau; // Tau can be 0 for integrator model
                calcResult = calculateAmigo(controllerType, numK, numTauAmigo, numTheta);
                setDisplayedMethod('AMIGO Tuning Rules');
                break;

            case 'imc':
                const numTauC = parseFloat(tauC);
                 if (isNaN(numK) || isNaN(numTau) || isNaN(numTheta) || isNaN(numTauC) || numTauC <= 0) {
                    throw new Error("K, τ, θ, and τ_c must be valid positive numbers.");
                }
                calcResult = calculateImc(imcModelCase, controllerType, numK, numTau, numTheta, numTauC);
                setDisplayedMethod('IMC');
                break;
        }
        
        if (!calcResult) {
            throw new Error("Calculation could not be completed for the selected options.");
        }

        setResult(calcResult);

    } catch (e: any) {
        setError(e.message);
        setResult(null);
    } finally {
        setLoading(false);
    }
  }, [k, tau, tau2, tau3, beta, theta, kcu, pu, tauC, zeta, tuningMethod, controllerType, itaeInputType, imcModelCase]);

  // Initial calculation on mount
  React.useEffect(() => {
    handleCalculate();
  }, []);

  // Effect to run calculation when inputs change
  React.useEffect(() => {
    handleCalculate();
  }, [
    k, tau, tau2, tau3, beta, theta, kcu, pu, tauC, zeta,
    tuningMethod, controllerType, itaeInputType, imcModelCase,
    handleCalculate
  ]);

  // Effect to reset simulation only when model-related parameters change (not slider values)
  React.useEffect(() => {
    if (result && isSimulationRunning) {
      resetSimulation();
    }
  }, [
    tuningMethod, controllerType, itaeInputType, imcModelCase
  ]);

  // Auto-start simulation when result is first calculated
  React.useEffect(() => {
    if (result && !isSimulationRunning && !loading) {
      startSimulation();
    }
  }, [result, loading]);

  // Update controller type when tuning method changes
  React.useEffect(() => {
    if (tuningMethod === 'imc') {
      // For IMC, auto-determine controller type from selected model case
      const requiredControllerType = getControllerTypeFromImcCase(imcModelCase);
      setControllerType(requiredControllerType);
    } else {
      // For other methods, check if current controller type is available
      const availableTypes = getAvailableControllerTypes(tuningMethod);
      if (!availableTypes.includes(controllerType)) {
        setControllerType(availableTypes[availableTypes.length - 1]);
      }
    }
  }, [tuningMethod]);

  // Auto-determine controller type from IMC model case for IMC method
  React.useEffect(() => {
    if (tuningMethod === 'imc') {
      const requiredControllerType = getControllerTypeFromImcCase(imcModelCase);
      if (controllerType !== requiredControllerType) {
        setControllerType(requiredControllerType);
      }
    }
  }, [imcModelCase, tuningMethod]);

  // Keep the refs in sync with the state
  React.useEffect(() => {
    showErrorGraphRef.current = showErrorGraph;
    graphModeRef.current = graphMode;
  }, [showErrorGraph, graphMode]);

  // Keep the result ref in sync with the result state
  React.useEffect(() => {
    resultRef.current = result;
  }, [result]);

  // Keep the warning refs in sync with the warning states
  React.useEffect(() => {
    divergenceWarningRef.current = divergenceWarning;
    resetCountdownRef.current = resetCountdown;
  }, [divergenceWarning, resetCountdown]);

  // Keep the setpoint ref in sync with the setpoint state
  React.useEffect(() => {
    currentSetpointRef.current = currentSetpoint;
  }, [currentSetpoint]);

  // ADD THIS NEW useEffect to sync the multi-view ref and data
  React.useEffect(() => {
    isMultiViewRef.current = isMultiView;
    // When switching to multi-view, copy the current simulation data to populate the charts
    if (isMultiView) {
      const currentData = simulationDataRef.current.length > 0 
        ? simulationDataRef.current.map(p => ({ ...p, output: 0 })) // Add output field with default value
        : [{ time: 0, pv: 0, sp: currentSetpointRef.current, output: 0 }];
      setMultiViewData(currentData);
    } else {
      // Clear data when leaving multi-view to save memory
      setMultiViewData([]);
    }
  }, [isMultiView]);

  // Sync showErrorGraph with graphMode for backward compatibility
  React.useEffect(() => {
    setShowErrorGraph(graphMode === 'error');
  }, [graphMode]);

  // --- Individual Tuning Rule Functions ---
  
  const calculateZieglerNichols = (ct: ControllerType, kcu_val: number, pu_val: number): TuningResult => {
      let kc: number | null = null, tauI: number | null = null, tauD: number | null = null;
      switch(ct) {
          case 'P':
              kc = 0.5 * kcu_val;
              break;
          case 'PI':
              kc = 0.45 * kcu_val;
              tauI = pu_val / 1.2;
              break;
          case 'PID':
              kc = 0.6 * kcu_val;
              tauI = pu_val / 2.0;
              tauD = pu_val / 8.0;
              break;
      }
      return { kc, tauI, tauD };
  };

  const calculateItae = (ct: ControllerType, it: ItaeInputType, k_val: number, tau_val: number, theta_val: number): TuningResult => {
      if (ct === 'P') {
          throw new Error("ITAE rules are not provided for P-only controllers in this set.");
      }
      
      const ratio = theta_val / tau_val;
      let kc: number | null = null, tauI: number | null = null, tauD: number | null = null;
      
      const rules = itaeParams[it][ct];
      
      if (!rules) throw new Error(`ITAE rules not available for ${ct} controller with ${it} input.`);

      const p_params = rules.P as { A: number; B: number };
      const y_p = p_params.A * Math.pow(ratio, p_params.B);
      kc = y_p / k_val;

      if (ct === 'PI' || ct === 'PID') {
          const i_params = rules.I as { A: number; B: number };
          if (it === 'setpoint') {
              const val = i_params.A + i_params.B * ratio;
              tauI = tau_val / val;
          } else {
              const y_i = i_params.A * Math.pow(ratio, i_params.B);
              tauI = y_i * tau_val;
          }
      }

      // FIX 1: Add a type guard `&& 'D' in rules` to ensure the 'D' property exists before access.
      if (ct === 'PID' && 'D' in rules) {
          const d_params = rules.D; // No `as` type assertion needed now
          const y_d = d_params.A * Math.pow(ratio, d_params.B);
          tauD = y_d * tau_val;
      }

      return { kc, tauI, tauD };
  };
    
  const calculateAmigo = (ct: ControllerType, k_val: number, tau_val: number, theta_val: number): TuningResult => {
      let kc: number | null = null, tauI: number | null = null, tauD: number | null = null;
      const isFOPTD = tau_val > 0;

      if (ct === 'PI') {
          if (isFOPTD) {
              const tau_norm = tau_val / theta_val;
              kc = (0.15 / k_val) + (0.35 - (tau_norm * theta_val) / Math.pow(theta_val + tau_val, 2)) * (tau_val / (k_val * theta_val));
              tauI = 0.35 * theta_val + (13 * Math.pow(tau_val, 2)) / (Math.pow(tau_val, 2) + 12 * tau_val * theta_val + 7 * Math.pow(theta_val, 2));
          } else { 
              kc = 0.35 / (k_val * theta_val);
              tauI = 13.4 * theta_val;
          }
      } else if (ct === 'PID') {
           if (isFOPTD) {
              kc = (1/k_val) * (0.2 + 0.45 * (tau_val / theta_val));
              tauI = (0.4 * theta_val + 0.8 * tau_val) / (theta_val + 0.1 * tau_val);
              tauD = (0.5 * tau_val * theta_val) / (0.3 * tau_val + theta_val);
           } else { 
              kc = 0.45 / (k_val * theta_val);
              tauI = 8.0 * theta_val;
              tauD = 0.5 * theta_val;
           }
      } else {
        throw new Error("AMIGO rules are only provided for PI and PID controllers.");
      }
      return { kc, tauI, tauD };
  };

  const calculateImc = (mcase: ImcModelCase, ct: ControllerType, k_val: number, tau_val: number, theta_val: number, tauc_val: number): TuningResult => {
      let kc: number | null = null, tauI: number | null = null, tauD: number | null = null;

      switch (mcase) {
          case 'A': { // Table Case A: First-Order
              if (ct === 'PID') throw new Error("IMC Case A is for PI controllers only.");
              kc = (1 / k_val) * (tau_val / tauc_val);
              tauI = tau_val;
              tauD = null; // No derivative action
              break;
          }

          case 'B': { // Table Case B: Second-Order (2 time constants)
              const numTau2 = parseFloat(tau2);
              if (isNaN(numTau2) || numTau2 <= 0) throw new Error("τ₂ must be a positive number.");
              kc = (1 / k_val) * ((tau_val + numTau2) / tauc_val);
              tauI = tau_val + numTau2;
              tauD = (tau_val * numTau2) / (tau_val + numTau2);
              break;
          }

          case 'C': { // Table Case C: Second-Order (underdamped)
              const numZeta = parseFloat(zeta);
              if (isNaN(numZeta) || numZeta <= 0) throw new Error("Zeta (ζ) must be a positive number.");
              kc = (1 / k_val) * ((2 * numZeta * tau_val) / tauc_val);
              tauI = 2 * numZeta * tau_val;
              tauD = tau_val / (2 * numZeta);
              break;
          }

          case 'D': { // Table Case D: SOPTD + RHP Zero
              const numZeta = parseFloat(zeta);
              if (isNaN(numZeta) || numZeta <= 0) throw new Error("Zeta (ζ) must be a positive number.");
              const numBeta = parseFloat(beta);
              if (isNaN(numBeta) || numBeta <= 0) throw new Error("Beta (β) must be a positive number.");
              
              kc = (1 / k_val) * ( (2 * numZeta * tau_val) / (tauc_val + numBeta) );
              tauI = 2 * numZeta * tau_val;
              tauD = tau_val / (2 * numZeta);
              break;
          }

          case 'E': { // Table Case E: Integrator
              if (ct === 'PID') throw new Error("IMC Case E is for PI controllers only.");
              kc = (1 / k_val) * (2 / tauc_val);
              tauI = 2 * tauc_val;
              tauD = null;
              break;
          }

          case 'F': { // Table Case F: Integrator + Lag
              kc = (1 / k_val) * ((2 * tauc_val + tau_val) / Math.pow(tauc_val, 2));
              tauI = 2 * tauc_val + tau_val;
              tauD = (2 * tauc_val * tau_val) / (2 * tauc_val + tau_val);
              break;
          }

          case 'G': { // Table Case G: FOPTD
              if (ct === 'PID') throw new Error("IMC Case G from the table is for PI controllers only.");
              if (tauc_val + theta_val <= 0) throw new Error("τc + θ must be positive.");
              
              kc = (1 / k_val) * (tau_val / (tauc_val + theta_val));
              tauI = tau_val;
              tauD = null;
              break;
          }

          case 'H': { // Table Case H: FOPTD (PID Controller)
              if (tauc_val + theta_val / 2 <= 0) throw new Error("τc + θ/2 must be positive.");
              if (2 * tau_val + theta_val === 0) throw new Error("Denominator (2τ + θ) cannot be zero.");

              kc = (1 / k_val) * ((tau_val + theta_val / 2) / (tauc_val + theta_val / 2));
              tauI = tau_val + theta_val / 2;
              tauD = (tau_val * theta_val) / (2 * tau_val + theta_val);
              break;
          }

          case 'I': { // Table Case I: 2nd-Order + Zero + Delay
              const numTau2 = parseFloat(tau2);
              if (isNaN(numTau2) || numTau2 <= 0) throw new Error("τ₂ must be positive.");
              const numTau3 = parseFloat(tau3);
              if (isNaN(numTau3) || numTau3 < 0) throw new Error("τ₃ must be non-negative.");
              if (tauc_val + theta_val <= 0) throw new Error("τc + θ must be positive.");

              const tau_sum = tau_val + numTau2 - numTau3;
              if (tau_sum <= 0) throw new Error("For Case I, (τ₁ + τ₂ - τ₃) must be positive.");

              kc = (1 / k_val) * (tau_sum / (tauc_val + theta_val));
              tauI = tau_sum;
              tauD = (tau_val * numTau2 - tau_sum * numTau3) / tau_sum;
              break;
          }

          case 'J': { // Table Case J: SOPTD + Zero + Delay
              const numZeta = parseFloat(zeta);
              const numTau3 = parseFloat(tau3);
              if (isNaN(numZeta) || isNaN(numTau3) || numZeta <= 0 || numTau3 < 0) throw new Error("τ, ζ, and τ₃ must be valid positive numbers.");
              if (tauc_val + theta_val <= 0) throw new Error("τc + θ must be positive.");

              const term = 2 * numZeta * tau_val - numTau3;
              if (term <= 0) throw new Error("For Case J, (2ζτ - τ₃) must be positive.");

              kc = (1 / k_val) * (term / (tauc_val + theta_val));
              tauI = term;
              tauD = (Math.pow(tau_val, 2) - term * numTau3) / term;
              break;
          }

          case 'K': { // Table Case K: 2nd-Order + RHP Zero + Delay
              const numTau2 = parseFloat(tau2);
              const numTau3 = parseFloat(tau3);
              if (isNaN(numTau2) || isNaN(numTau3) || numTau2 <= 0 || numTau3 <= 0) throw new Error("τ₁, τ₂, and τ₃ must be positive.");
              if (tauc_val + numTau3 + theta_val === 0) throw new Error("Denominator (τc + τ₃ + θ) cannot be zero.");

              const common_term_1 = (numTau3 * theta_val) / (tauc_val + numTau3 + theta_val);
              const tau_i_val = tau_val + numTau2 + common_term_1;
              if (tau_i_val === 0) throw new Error("Denominator for second term of τD cannot be zero.");
              
              const kc_k_num = tau_i_val; // Per table, numerator of Kc*K is the same as τI 
              const kc_k_den = tauc_val + numTau3 + theta_val;

              kc = (1 / k_val) * (kc_k_num / kc_k_den);
              tauI = tau_i_val; // Formula for τI from table 
              tauD = common_term_1 + (tau_val * numTau2) / tau_i_val; // Formula for τD from table 
              break;
          }

          case 'L': { // Table Case L: SOPTD + RHP Zero + Delay
              const numZeta = parseFloat(zeta);
              const numTau3 = parseFloat(tau3);
              if (isNaN(numZeta) || isNaN(numTau3) || numZeta <= 0 || numTau3 <= 0) throw new Error("τ, ζ, and τ₃ must be positive.");
              if (tauc_val + numTau3 + theta_val === 0) throw new Error("Denominator (τc + τ₃ + θ) cannot be zero.");

              const common_term_1 = (numTau3 * theta_val) / (tauc_val + numTau3 + theta_val);
              const tau_i_val = (2 * numZeta * tau_val) + common_term_1;
              if (tau_i_val === 0) throw new Error("Denominator for second term of τD cannot be zero.");

              const kc_k_num = tau_i_val; // Per table, numerator of Kc*K is the same as τI 
              const kc_k_den = tauc_val + numTau3 + theta_val;

              kc = (1 / k_val) * (kc_k_num / kc_k_den);
              tauI = tau_i_val; // Formula for τI from table 
              tauD = common_term_1 + (Math.pow(tau_val, 2)) / tau_i_val; // Formula for τD from table 
              break;
          }

          case 'M': { // Table Case M: Integrator + Delay
              if (ct === 'PID') throw new Error("IMC Case M from the table is for PI controllers only.");
              if (Math.pow(tauc_val + theta_val, 2) === 0) throw new Error("Denominator cannot be zero.");

              const kc_k_term = (2 * tauc_val + theta_val) / Math.pow(tauc_val + theta_val, 2);
              kc = (1 / k_val) * kc_k_term;
              tauI = 2 * tauc_val + theta_val;
              tauD = null; // Set to null as per the table.
              break;
          }

          case 'N': { // Model: Ke⁻θs/s
              const term1 = tauc_val + theta_val / 2;
              if (term1 <= 0) throw new Error("τc + θ/2 must be positive.");
              const term2 = 2 * tauc_val + theta_val;
              
              kc = (1 / k_val) * (term2 / Math.pow(term1, 2));
              tauI = term2;
              tauD = (tauc_val * theta_val + Math.pow(theta_val, 2) / 4) / term2;
              break;
          }

          case 'O': { // Model: Ke⁻θs/(s(τs+1))
              const term1 = tauc_val + theta_val;
              if (term1 <= 0) throw new Error("τc + θ must be positive.");
              const term2 = 2 * tauc_val + tau_val + theta_val;
              
              kc = (1 / k_val) * (term2 / Math.pow(term1, 2));
              tauI = term2;
              tauD = ((2 * tauc_val + theta_val) * tau_val) / term2;
              break;
          }
          
          default:
              throw new Error(`IMC model case ${mcase} not implemented.`);
      }
      return { kc, tauI, tauD };
  };

  const calculateFrequencyResponse = (
    processParams: { K: number; tau: number; theta: number; },
    controller: { Kc: number; tauI: number | null; tauD: number | null }
  ) => {
    const { K, tau, theta } = processParams;
    const { Kc, tauI, tauD } = controller;
    const bodeData = { magnitude: [] as [number, number][], phase: [] as [number, number][] };
    const nyquistData = [] as [number, number][];
    const nyquistDataNegative = [] as [number, number][]; // For negative frequencies
    const points = 200;

    for (let i = 0; i <= points; i++) {
      const w = Math.pow(10, -2 + (i * 4) / points); // Freq range 10^-2 to 10^2

      // Gp(s) = K*e^(-s*theta)/(tau*s+1)
      const gp_real = K * (Math.cos(w*theta) - w*tau*Math.sin(w*theta)) / (1 + w*w*tau*tau);
      const gp_imag = K * (-Math.sin(w*theta) - w*tau*Math.cos(w*theta)) / (1 + w*w*tau*tau);
      
      // Gc(s) = Kc*(1 + 1/(tauI*s) + tauD*s)
      const gc_real = Kc;
      let gc_imag = (Kc * (tauD || 0) * w) - (tauI ? Kc / (tauI * w) : 0);

      // Gol(s) = Gc(s) * Gp(s)
      const gol_real = gc_real * gp_real - gc_imag * gp_imag;
      const gol_imag = gc_real * gp_imag + gc_imag * gp_real;

      const magnitude = Math.sqrt(gol_real**2 + gol_imag**2);
      bodeData.magnitude.push([w, 20 * Math.log10(magnitude)]);
      bodeData.phase.push([w, Math.atan2(gol_imag, gol_real) * (180 / Math.PI)]);

      // --- Add this for Nyquist Plot ---
      nyquistData.push([gol_real, gol_imag]);
      nyquistDataNegative.push([gol_real, -gol_imag]); // Complex conjugate for w < 0
    }
    
    // The full Nyquist contour is the negative frequency data reversed, then the positive data
    const fullNyquistData = [...nyquistDataNegative.reverse(), ...nyquistData];

    return { bodeData, nyquistData: fullNyquistData };
  };

  // Simulation Control Functions - Live changes without reset
  const handleSetpointChange = (direction: 'up' | 'down') => {
    setCurrentSetpoint(prevSetpoint => {
      let newSetpoint = prevSetpoint;
      
      if (direction === 'up' && prevSetpoint < 2) {
        newSetpoint = prevSetpoint + 1;
      } else if (direction === 'down' && prevSetpoint > 0) {
        newSetpoint = prevSetpoint - 1;
      }
      
      currentSetpointRef.current = newSetpoint;
      
      simulationDataRef.current = simulationDataRef.current.map(point => ({
        ...point,
        sp: newSetpoint
      }));
      
      const echartsInstance = echartsRef.current?.getEchartsInstance();
      if (echartsInstance) {
        if (graphModeRef.current === 'bode' || graphModeRef.current === 'nyquist' || graphModeRef.current === 'power') {
          // No update needed for these modes when setpoint changes,
          // as their data is not directly tied to the setpoint value itself
          // in a way that requires an immediate series data update.
        } else if (graphModeRef.current === 'error') {
          echartsInstance.setOption({
            series: [
              { name: 'Control Error', data: simulationDataRef.current.map(p => [p.time, p.sp - p.pv]) },
              { name: 'Zero Line', data: simulationDataRef.current.map(p => [p.time, 0]) }
            ]
          });
        } else { // 'process'
          echartsInstance.setOption({
            series: [
              { name: 'Process Variable', data: simulationDataRef.current.map(p => [p.time, p.pv]) },
              { name: 'Setpoint', data: simulationDataRef.current.map(p => [p.time, p.sp]) },
              { name: 'Control Error Range', data: [simulationDataRef.current[simulationDataRef.current.length - 1]?.time || 0] }
            ]
          });
        }
      }
      
      return newSetpoint;
    });
  };

  const triggerDisturbance = (direction: 'up' | 'down') => {
    const currentPV = simulationStateRef.current.pv;
    
    if (direction === 'up' && currentPV < 2.0) {
      simulationStateRef.current.pv = currentPV + 0.5;
    } else if (direction === 'down' && currentPV > 0.0) {
      simulationStateRef.current.pv = currentPV - 0.5;
    }
    
    simulationStateRef.current.prev_pv = simulationStateRef.current.pv;
  };

  // Start simulation automatically and provide reset functionality
  const startSimulation = (forceStart = false) => {
    if (isSimulationRunning && !forceStart) return; // Already running, unless forced
    
    const dt = 0.1;
    const numTheta = parseFloat(theta);
    // If theta is very small, use minimal delay buffer to reduce startup delay
    const delaySteps = numTheta < 0.1 ? 0 : Math.floor(numTheta / dt);

    simulationStateRef.current = {
      pv: 0,
      prev_pv: 0,
      integral_error: 0,
      filtered_derivative: 0,
      y1: 0,
      y2: 0,
      delayBuffer: new Array(delaySteps).fill(0),
    };
    
    // Reset data array using the current ref value
    simulationDataRef.current = [{ time: 0, pv: 0, sp: currentSetpointRef.current, output: 0 }];
    setSimulationTime(0);
    setIsSimulationRunning(true);
    
    const startTime = Date.now();
    
    simulationIntervalRef.current = setInterval(() => {
        if (!resultRef.current) {
            return;
        }

        const elapsedTime = (Date.now() - startTime) / 1000;
        const time = parseFloat(elapsedTime.toFixed(1));
        const dt = 0.05;

        // --- Get parameters for this tick ---
        const processParams = {
            K: parseFloat(k),
            tau: parseFloat(tau),
            theta: parseFloat(theta),
            tau2: parseFloat(tau2),
            zeta: parseFloat(zeta)
        };
        const controllerParams = {
            Kc: resultRef.current.kc || 0,
            tauI: resultRef.current.tauI,
            tauD: resultRef.current.tauD
        };
        const simModelCase = tuningMethod === 'imc' ? imcModelCase : 'G';
        
        // --- Run one step of the simulation ---
        const newDataPoint = runSingleSimulationStep(
            time,
            dt,
            processParams,
            controllerParams,
            simModelCase,
            currentSetpointRef.current
        );
        
        simulationDataRef.current.push(newDataPoint);

        // Update current PV and error for the UI
        setCurrentPV(newDataPoint.pv);
        setControlError(newDataPoint.sp - newDataPoint.pv);

        // Check for divergence - if control error exceeds 3, trigger warning
        const currentError = Math.abs(newDataPoint.sp - newDataPoint.pv);
        if (currentError > 3 && !divergenceWarning && !countdownIntervalRef.current) {
          setDivergenceWarning(true);
          setResetCountdown(3);
          
          // Start countdown timer (only if one isn't already running)
          countdownIntervalRef.current = setInterval(() => {
            setResetCountdown(prev => {
              if (prev <= 1) {
                clearInterval(countdownIntervalRef.current!);
                countdownIntervalRef.current = null;
                // Reset the simulation
                resetSimulation();
                setDivergenceWarning(false);
                setResetCountdown(0);
                return 0;
              }
              return prev - 1;
            });
          }, 1000);
        }

        // Keep only the last 2000 data points to prevent memory issues
        if (simulationDataRef.current.length > 2000) {
          simulationDataRef.current = simulationDataRef.current.slice(-2000);
        }

        // Calculate dynamic x-axis range to keep current time on left/middle of screen
        const currentTime = time;
        const timeWindow = 30; // Show 30 seconds worth of data
        // Continuous sliding window – current time appears ~75% across
        const xMin = Math.max(0, currentTime - timeWindow * 0.75);
        const xMax = Math.max(timeWindow, currentTime + timeWindow * 0.25);

        // --- ROBUST LIVE UPDATE LOGIC ---
        if (isMultiViewRef.current) {
          if (isOverlappedViewRef.current) {
            // --- Update logic for the new OVERLAPPED graph ---
            const overlappedChart = overlappedChartRef.current?.getEchartsInstance();
            if (overlappedChart) {
              const allData = simulationDataRef.current;
              overlappedChart.setOption({
                animationDuration: 0,
                animationDurationUpdate: 0,
                xAxis: { min: xMin, max: xMax, animation: false },
                series: [
                  { name: 'Process Variable', data: allData.map(p => [p.time, p.pv]) },
                  { name: 'Setpoint', data: allData.map(p => [p.time, p.sp]) },
                  { name: 'Control Error', data: allData.map(p => [p.time, p.sp - p.pv]) },
                  { name: 'Controller Output', data: allData.map(p => [p.time, p.output]) }
                ]
              });
            }
          } else {
            // --- This block will now run correctly for MULTI-VIEW (separated) ---
            const allData = simulationDataRef.current;
            const processChart = multiViewProcessRef.current?.getEchartsInstance();
            const errorChart = multiViewErrorRef.current?.getEchartsInstance();
            const outputChart = multiViewOutputRef.current?.getEchartsInstance();

            if (processChart) {
              processChart.setOption({
                animationDuration: 0,
                animationDurationUpdate: 0,
                xAxis: { min: xMin, max: xMax, animation: false },
                series: [{ data: allData.map(p => [p.time, p.pv]) }, { data: allData.map(p => [p.time, p.sp]) }]
              });
            }
            if (errorChart) {
              errorChart.setOption({
                animationDuration: 0,
                animationDurationUpdate: 0,
                xAxis: { min: xMin, max: xMax, animation: false },
                series: [{ data: allData.map(p => [p.time, p.sp - p.pv]) }, { data: allData.map(p => [p.time, 0]) }]
              });
            }
            if (outputChart) {
              outputChart.setOption({
                animationDuration: 0,
                animationDurationUpdate: 0,
                xAxis: { min: xMin, max: xMax, animation: false },
                yAxis: {
                  animationDuration: 0,
                  animationDurationUpdate: 0
                },
                series: [{ data: allData.map(p => [p.time, p.output]) }]
              });
            }
          }
        } else {
          // --- This is the original, working logic for SINGLE-VIEW ---
          const singleChart = echartsRef.current?.getEchartsInstance();
          if (singleChart) {
            if (graphModeRef.current === 'bode' || graphModeRef.current === 'nyquist') {
              // No update needed
            } else if (graphModeRef.current === 'error') {
              singleChart.setOption({
                xAxis: { min: xMin, max: xMax },
                series: [
                  { name: 'Control Error', data: simulationDataRef.current.map(p => [p.time, p.sp - p.pv]) },
                  { name: 'Zero Line', data: simulationDataRef.current.map(p => [p.time, 0]) }
                ]
              });
            } else if (graphModeRef.current === 'power') {
              singleChart.setOption({
                xAxis: { min: xMin, max: xMax },
                yAxis: {
                  min: -0.5,
                  max: 2.5,
                  animation: false, // ensures y-axis does not glide
                },
                series: [
                  { name: 'Controller Output', data: simulationDataRef.current.map(p => [p.time, p.output]) }
                ]
              });
            } else { // 'process'
              singleChart.setOption({
                xAxis: { min: xMin, max: xMax },
                series: [
                  { name: 'Process Variable', data: simulationDataRef.current.map(p => [p.time, p.pv]) },
                  { name: 'Setpoint', data: simulationDataRef.current.map(p => [p.time, p.sp]) },
                  { name: 'Control Error Range', data: [simulationDataRef.current[simulationDataRef.current.length - 1]?.time || 0] }
                ]
              });
            }
          }
        }

        setSimulationTime(time);
        if (time >= 1000) resetSimulation();
    }, 100);
  };

  // Reset simulation - stops current simulation and restarts it
  const resetSimulation = () => {
    if (simulationIntervalRef.current) {
      clearInterval(simulationIntervalRef.current);
      simulationIntervalRef.current = null;
    }
    if (countdownIntervalRef.current) {
      clearInterval(countdownIntervalRef.current);
      countdownIntervalRef.current = null;
    }
    setIsSimulationRunning(false);
    setDivergenceWarning(false);
    setResetCountdown(0);
    
    // Force start simulation immediately
    startSimulation(true);
  };

  // Cleanup simulation interval on unmount
  React.useEffect(() => {
    return () => {
      if (simulationIntervalRef.current) {
        clearInterval(simulationIntervalRef.current);
      }
      if (countdownIntervalRef.current) {
        clearInterval(countdownIntervalRef.current);
      }
    };
  }, []);

  // Add these helper functions for multi-view chart options
  function getMultiViewProcessOptions(isDark: boolean, data: any[]): EChartsOption {
    const textColor = isDark ? 'white' : '#000000';
    const tooltipBg = isDark ? '#08306b' : '#f8f9fa';
    const tooltipBorder = isDark ? '#55aaff' : '#dee2e6';
    
    return {
      backgroundColor: 'transparent',
      title: { 
        text: 'Process Simulation', 
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' }, 
        left: 'center',
        top: '2%'
      },
      grid: { left: '8%', right: '8%', bottom: '15%', top: '12%', containLabel: true },
      tooltip: {
        show: true,
        trigger: 'axis',
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        borderWidth: 1,
        textStyle: { 
          color: textColor,
          fontSize: 12,
          fontFamily: 'Merriweather Sans'
        },
        formatter: function(params: any) {
          if (!params || !Array.isArray(params) || params.length < 1) return '';
          const time = params[0]?.value[0];
          if (time === undefined) return '';

          let tooltip = `Time: ${parseFloat(time).toPrecision(3)} s<br/>`;

          params.forEach((param: any) => {
            const value = parseFloat(param.value[1]).toPrecision(3);
            const color = param.seriesName === 'Process Variable' ? '#00ff51ff' : '#ef4444';
            tooltip += `<span style="color: ${color};">${param.seriesName}: ${value}</span><br/>`;
          });

          return tooltip;
        },
        axisPointer: {
          type: 'cross',
          label: {
            show: true,
            backgroundColor: tooltipBg,
            color: textColor,
            borderColor: tooltipBorder,
            borderWidth: 1,
            shadowBlur: 0,
            shadowColor: 'transparent',
            fontFamily: 'Merriweather Sans'
          },
          crossStyle: { color: '#999' }
        }
      },
      xAxis: {
        type: 'value',
        name: 'Time (s)',
        min: 0,
        max: 30,
        nameLocation: 'middle',
        nameGap: 25,
        nameTextStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: {
          color: textColor,
          fontSize: 12,
          fontFamily: 'Merriweather Sans',
          showMinLabel: false,
          showMaxLabel: false
        },
        splitLine: { show: false }
      },
      yAxis: {
        type: 'value',
        name: 'Process Variable',
        nameLocation: 'middle',
        nameGap: 45,
        nameTextStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        min: -0.5,
        max: 2.5,
        scale: true,
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false }
      },
      legend: {
        data: [
          { name: 'Process Variable', icon: 'rect', itemStyle: { color: '#00ff51ff' } },
          { name: 'Setpoint', icon: 'rect', itemStyle: { color: '#ef4444' } }
        ],
        bottom: '2%',
        left: 'center',
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        itemWidth: 25,
        itemHeight: 2
      },
      series: [
        { 
          name: 'Process Variable', 
          type: 'line', 
          data: data.map(p => [p.time, p.pv]), 
          showSymbol: false, 
          lineStyle: { color: '#00ff51ff', width: 3 },
          smooth: true
        },
        { 
          name: 'Setpoint', 
          type: 'line', 
          data: data.map(p => [p.time, p.sp]), 
          showSymbol: false, 
          lineStyle: { color: '#ef4444', type: 'dashed', width: 2 } 
        }
      ]
    };
  };

  function getMultiViewErrorOptions(isDark: boolean, data: any[]): EChartsOption {
    const textColor = isDark ? 'white' : '#000000';
    const tooltipBg = isDark ? '#08306b' : '#f8f9fa';
    const tooltipBorder = isDark ? '#55aaff' : '#dee2e6';
    
    return {
      backgroundColor: 'transparent',
      title: { 
        text: 'Control Error', 
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' }, 
        left: 'center',
        top: '2%'
      },
      grid: { left: '8%', right: '8%', bottom: '15%', top: '12%', containLabel: true },
      tooltip: {
        show: true,
        trigger: 'axis',
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        borderWidth: 1,
        textStyle: { 
          color: textColor,
          fontSize: 12,
          fontFamily: 'Merriweather Sans'
        },
        formatter: function(params: any) {
          if (!params || !Array.isArray(params)) return '';
          const time = params[0]?.value[0];
          if (time === undefined) return '';
          let tooltip = `Time: ${parseFloat(time).toPrecision(3)} s<br/>`;
          params.forEach((param: any) => {
            if (param.seriesName === 'Control Error') {
              const value = parseFloat(param.value[1]).toPrecision(3);
              tooltip += `<span style=\"color: #3b82f6;\">Control Error: ${value}</span><br/>`;
            }
          });
          return tooltip;
        },
        axisPointer: {
          type: 'cross',
          label: {
            show: true,
            backgroundColor: tooltipBg,
            color: textColor,
            borderColor: tooltipBorder,
            borderWidth: 1,
            shadowBlur: 0,
            shadowColor: 'transparent',
            fontFamily: 'Merriweather Sans'
          },
          crossStyle: { color: '#999' }
        }
      },
      xAxis: {
        type: 'value',
        name: 'Time (s)',
        min: 0,
        max: 30,
        nameLocation: 'middle',
        nameGap: 25,
        nameTextStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: {
          color: textColor,
          fontSize: 12,
          fontFamily: 'Merriweather Sans',
          showMinLabel: false,
          showMaxLabel: false
        },
        splitLine: { show: false }
      },
      yAxis: {
        type: 'value',
        name: 'Control Error',
        nameLocation: 'middle',
        nameGap: 45,
        nameTextStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        min: -3,
        max: 3,
        scale: true,
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false }
      },
      legend: {
        data: [
          { 
            name: 'Control Error', 
            icon: 'rect',
            itemStyle: { color: '#3b82f6' }
          }
        ],
        bottom: '2%',
        left: 'center',
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        itemWidth: 25,
        itemHeight: 2
      },
      series: [
        { 
          name: 'Control Error', 
          type: 'line', 
          data: data.map(p => [p.time, p.sp - p.pv]), 
          showSymbol: false, 
          lineStyle: { color: '#3b82f6', width: 3 },
          smooth: true
        },
        { 
          name: 'Zero Line', 
          type: 'line', 
          data: data.map(p => [p.time, 0]), 
          showSymbol: false, 
          lineStyle: { color: '#888888', type: 'dashed', width: 1 },
          silent: true
        }
      ]
    };
  };

  function getMultiViewOutputOptions(isDark: boolean, data: any[]): EChartsOption {
    const textColor = isDark ? 'white' : '#000000';
    const tooltipBg = isDark ? '#08306b' : '#f8f9fa';
    const tooltipBorder = isDark ? '#55aaff' : '#dee2e6';
    
    return {
      backgroundColor: 'transparent',
      title: { 
        text: 'Controller Output', 
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' }, 
        left: 'center',
        top: '2%'
      },
      grid: { left: '8%', right: '8%', bottom: '15%', top: '12%', containLabel: true },
      tooltip: {
        show: true,
        trigger: 'axis',
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        borderWidth: 1,
        textStyle: { 
          color: textColor,
          fontSize: 12,
          fontFamily: 'Merriweather Sans'
        },
        formatter: function(params: any) {
          if (!params || !Array.isArray(params) || params.length === 0) return '';
          const time = params[0]?.value[0];
          if (time === undefined) return '';
          let tooltip = `Time: ${parseFloat(time).toPrecision(3)} s<br/>`;
          params.forEach((param: any) => {
            if (param.seriesName === 'Controller Output') {
              const value = parseFloat(param.value[1]).toPrecision(3);
              tooltip += `<span style="color: ${purpleColor};">${param.seriesName}: ${value}</span><br/>`;
            }
          });
          return tooltip;
        },
        axisPointer: {
          type: 'cross',
          label: {
            show: true,
            backgroundColor: tooltipBg,
            color: textColor,
            borderColor: tooltipBorder,
            borderWidth: 1,
            shadowBlur: 0,
            shadowColor: 'transparent',
            fontFamily: 'Merriweather Sans'
          },
          crossStyle: { color: '#999' }
        }
      },
      xAxis: {
        type: 'value',
        name: 'Time (s)',
        min: 0,
        max: 30,
        nameLocation: 'middle',
        nameGap: 25,
        nameTextStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: {
          color: textColor,
          fontSize: 12,
          fontFamily: 'Merriweather Sans',
          showMinLabel: false,
          showMaxLabel: false
        },
        splitLine: { show: false }
      },
      yAxis: {
        type: 'value',
        name: 'Controller Output',
        nameLocation: 'middle',
        nameGap: 45,
        nameTextStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        min: -0.5,
        max: 2.5,
        scale: true,
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false }
      },
      legend: {
        data: [{ name: 'Controller Output', icon: 'rect', itemStyle: { color: purpleColor } }],
        bottom: '2%',
        left: 'center',
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        itemWidth: 25,
        itemHeight: 2
      },
      series: [
        { 
          name: 'Controller Output', 
          type: 'line', 
          data: data.map(p => [p.time, p.output]), 
          showSymbol: false, 
          lineStyle: { color: purpleColor, width: 3 },
          smooth: true
        }
      ]
    };
  };

  const formatResult = (value: number | null) => {
    if (value === null || isNaN(value)) return 'N/A';
    if (Math.abs(value) < 0.001 && value !== 0) return value.toExponential(2);
    if (value === 0) return '0.00';
    
    // Always format to exactly 3 significant figures
    return value.toPrecision(3);
  }

  const frequencyResponseData = useMemo(() => {
    if (!result) return null;
    const processParams = { K: parseFloat(k), tau: parseFloat(tau), theta: parseFloat(theta) };
    if (isNaN(processParams.K) || isNaN(processParams.tau) || isNaN(processParams.theta)) return null;
    const controller = { Kc: result.kc || 0, tauI: result.tauI, tauD: result.tauD };
    return calculateFrequencyResponse(processParams, controller);
  }, [result, k, tau, theta]);

  // Graph options for ECharts - Real-time simulation
  const graphOptions = useMemo((): EChartsOption => {
    // Theme-dependent colors
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';
    const tooltipBg = isDark ? '#08306b' : '#f8f9fa';
    const tooltipBorder = isDark ? '#55aaff' : '#dee2e6';
    const axisLineStyle = { lineStyle: { color: textColor } };

    // --- BODE PLOT CONFIGURATION ---
    if (graphMode === 'bode') {
      if (!frequencyResponseData) return {};
      return {
        backgroundColor: 'transparent',
        animation: false,
        title: { 
          text: `${controllerType} Controller Bode Plot`, 
          left: 'center',
          top: '0%',
          textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } 
        },
        tooltip: {
          show: true,
          trigger: 'axis',
          backgroundColor: tooltipBg,
          borderColor: tooltipBorder,
          borderWidth: 1,
          textStyle: { 
            color: textColor,
            fontSize: 12,
            fontFamily: 'Merriweather Sans'
          },
          formatter: function(params: any) {
            if (!params || !Array.isArray(params) || params.length === 0) return '';
            
            const frequency = params[0]?.value[0];
            if (frequency === undefined) return '';
            
            let tooltip = `<span style="color: ${textColor};">Frequency: ${parseFloat(frequency).toPrecision(3)} rad/s</span><br/>`;
            
            params.forEach((param: any) => {
              if (param.seriesName === 'Magnitude') {
                const magnitude = param.value[1];
                tooltip += `<span style="color: #00ff51ff;">Magnitude: ${parseFloat(magnitude).toPrecision(3)} dB</span><br/>`;
              } else if (param.seriesName === 'Phase') {
                const phase = param.value[1];
                tooltip += `<span style="color: #3b82f6;">Phase: ${parseFloat(phase).toPrecision(3)} deg</span><br/>`;
              }
            });
            
            return tooltip;
          },
          axisPointer: {
            type: 'cross',
            label: {
              show: true,
              backgroundColor: tooltipBg,
              color: textColor,
              borderColor: tooltipBorder,
              borderWidth: 1,
              shadowBlur: 0,
              shadowColor: 'transparent',
              fontFamily: 'Merriweather Sans',
              formatter: function(params: any) {
                // Format axis pointer labels to 3 significant figures
                if (params.axisDimension === 'x') {
                  // For frequency (x-axis), format to 3 sig figs
                  return parseFloat(params.value).toPrecision(3);
                } else {
                  // For magnitude/phase (y-axis), format to 3 sig figs
                  return parseFloat(params.value).toPrecision(3);
                }
              }
            },
            crossStyle: {
              color: '#999'
            }
          }
        },
        axisPointer: { link: [{ xAxisIndex: 'all' }] },
        grid: [
          { left: '8%', right: '5%', top: '15%', height: '35%', containLabel: true }, // Magnitude plot
          { left: '8%', right: '5%', bottom: '12%', height: '35%', containLabel: true }  // Phase plot
        ],
        xAxis: [
          { 
            type: 'log', 
            gridIndex: 0, 
            axisLine: { lineStyle: { color: textColor } },
            axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { show: false },
            splitLine: { show: false }
          },
          { 
            type: 'log', 
            gridIndex: 1, 
            name: 'Frequency (rad/s)', 
            nameLocation: 'middle', 
            nameGap: 30,
            nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } },
            axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
            splitLine: { show: false }
          }
        ],
        yAxis: [
          { 
            type: 'value', 
            gridIndex: 0, 
            name: 'Magnitude (dB)', 
            nameLocation: 'middle', 
            nameGap: 50,
            nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } },
            axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
            splitLine: { show: false }
          },
          { 
            type: 'value', 
            gridIndex: 1, 
            name: 'Phase (deg)', 
            nameLocation: 'middle', 
            nameGap: 50,
            nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
            axisLine: { lineStyle: { color: textColor } },
            axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
            axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
            splitLine: { show: false }
          }
        ],
        series: [
          { 
            name: 'Magnitude', 
            type: 'line', 
            xAxisIndex: 0, 
            yAxisIndex: 0, 
            showSymbol: false, 
            data: frequencyResponseData.bodeData.magnitude, 
            lineStyle: { color: '#00ff51ff', width: 2 } 
          },
          { 
            name: 'Phase', 
            type: 'line', 
            xAxisIndex: 1, 
            yAxisIndex: 1, 
            showSymbol: false, 
            data: frequencyResponseData.bodeData.phase, 
            lineStyle: { color: '#3b82f6', width: 2 } 
          }
        ]
      };
    }
    
    // --- NYQUIST PLOT CONFIGURATION ---
    if (graphMode === 'nyquist') {
      if (!frequencyResponseData) return {};
      return {
        backgroundColor: 'transparent',
        animation: false,
        title: { 
          text: `${controllerType} Controller Nyquist Plot`, 
          left: 'center',
          top: '0%',
          textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } 
        },
        grid: { left: '5%', right: '5%', bottom: '12%', top: '8%', containLabel: true },
        tooltip: {
          show: true,
          trigger: 'item',
          backgroundColor: tooltipBg,
          borderColor: tooltipBorder,
          borderWidth: 1,
          textStyle: { 
            color: textColor,
            fontSize: 12,
            fontFamily: 'Merriweather Sans'
          },
          formatter: function(params: any) {
            if (!params || !params.value || !Array.isArray(params.value)) return '';
            
            const real = params.value[0];
            const imag = params.value[1];
            
            return `<span style="color: #00ff51ff;">Real: ${parseFloat(real).toPrecision(3)}</span><br/>` +
                   `<span style="color: #00ff51ff;">Imaginary: ${parseFloat(imag).toPrecision(3)}</span>`;
          },
          axisPointer: {
            type: 'cross',
            label: {
              show: true,
              backgroundColor: tooltipBg,
              color: textColor,
              borderColor: tooltipBorder,
              borderWidth: 1,
              shadowBlur: 0,
              shadowColor: 'transparent',
              fontFamily: 'Merriweather Sans'
            },
            crossStyle: {
              color: '#999'
            }
          }
        },
        xAxis: { 
          type: 'value', 
          name: 'Real Axis', 
          nameLocation: 'middle',
          nameGap: 30,
          nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
          scale: true, 
          axisLine: { lineStyle: { color: textColor } },
          axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
          axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
          splitLine: { show: false }
        },
        yAxis: { 
          type: 'value', 
          name: 'Imaginary Axis', 
          nameLocation: 'middle',
          nameGap: 50,
          nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
          scale: true, 
          axisLine: { lineStyle: { color: textColor } },
          axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
          axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
          splitLine: { show: false }
        },
        series: [
          {
            name: 'Open Loop Response',
            type: 'line',
            showSymbol: false,
            data: frequencyResponseData.nyquistData,
            lineStyle: { color: '#00ff51ff', width: 2 },
            markPoint: {
              // These are the new, combined mark points
              data: [
                // 1. The existing critical point marker
                {
                  name: 'Critical Point',
                  coord: [-1, 0],
                  symbol: 'cross',
                  symbolSize: 12,
                  label: {
                    show: true,
                    formatter: '(-1, 0)',
                    color: textColor,
                    fontFamily: 'Merriweather Sans',
                    fontSize: 14
                  }
                },
                // 2. The new direction arrows, generated dynamically
                ...(() => {
                  const arrows = [];
                  const nyquistPoints = frequencyResponseData.nyquistData;
                  if (!nyquistPoints || nyquistPoints.length < 2) return [];

                  // Place ~6 arrows, avoiding the start/end where the plot can be chaotic
                  const numArrows = 6;
                  const startIndex = Math.floor(nyquistPoints.length * 0.1);
                  const endIndex = Math.floor(nyquistPoints.length * 0.9);
                  const step = Math.floor((endIndex - startIndex) / numArrows);

                  if (step <= 0) return []; // Not enough data points

                  for (let i = startIndex; i < endIndex; i += step) {
                    const currentPoint = nyquistPoints[i];
                    // Look ahead a few points to get a stable direction vector
                    const nextPoint = nyquistPoints[Math.min(i + 5, nyquistPoints.length - 1)];

                    if (currentPoint && nextPoint) {
                      const dx = nextPoint[0] - currentPoint[0];
                      const dy = nextPoint[1] - currentPoint[1];

                      // Only add an arrow if there's noticeable movement
                      if (Math.abs(dx) > 1e-6 || Math.abs(dy) > 1e-6) {
                        const angle = Math.atan2(dy, dx) * (180 / Math.PI);
                        arrows.push({
                          name: `Direction Arrow ${i}`,
                          coord: currentPoint,
                          // FIX: Use a custom SVG path for the arrow symbol.
                          // This places the tip of the arrow directly on the curve.
                          symbol: 'path://M 0 0 L -10 5 L -10 -5 Z',
                          symbolSize: 12,
                          symbolRotate: angle,
                          itemStyle: {
                            color: '#00ff51ff'
                          },
                          label: { show: false } // Hide labels for arrows
                        });
                      }
                    }
                  }
                  return arrows;
                })()
              ]
            }
          } as any
        ]
      };
    }
    
    // --- EXISTING PROCESS & ERROR PLOT CONFIGURATION ---
    // Only initialize chart if we have calculated results
    if (!result) {
      console.log('Graph not showing - no results:', { result });
      return {};
    }

    // Mark graph as initialized once we have results
    if (!graphInitialized) {
      setGraphInitialized(true);
    }

    return {
      backgroundColor: 'transparent',
      animation: false, // Disable animation for live updates
      title: { 
        text: graphMode === 'error' ? `${controllerType} Controller Error` : 
              graphMode === 'power' ? `${controllerType} Controller Output` : 
              `${controllerType} Controller Simulation`, 
        left: 'center',
        top: '2%',
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } 
      },
      grid: { left: '5%', right: '5%', bottom: '15%', top: '8%', containLabel: true },
      xAxis: {
        type: 'value',
        name: 'Time (s)',
        nameLocation: 'middle',
        nameGap: 35,
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { 
          color: textColor, 
          fontSize: 16, 
          fontFamily: 'Merriweather Sans',
          showMinLabel: false,
          showMaxLabel: false
        },
        splitLine: { show: false },
        min: 0,
        max: 30 // Dynamic range will be set during simulation
      },
      yAxis: {
        type: 'value',
        name: graphMode === 'error' ? 'Control Error' : 
              graphMode === 'power' ? 'Controller Output' : 
              'Process Variable',
        nameLocation: 'middle',
        nameGap: 40,
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        min: graphMode === 'error' ? -3 : graphMode === 'power' ? (currentSetpointRef.current || 1) - 1 : -0.5,
        max: graphMode === 'error' ? 3 : graphMode === 'power' ? (currentSetpointRef.current || 1) + 1 : 2.5,
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false }
      },
      tooltip: {
        show: true,
        trigger: 'axis',
        backgroundColor: tooltipBg,
        borderColor: tooltipBorder,
        borderWidth: 1,
        textStyle: { 
          color: isDark ? 'white' : '#000000',
          fontSize: 12,
          fontFamily: 'Merriweather Sans'
        },
        formatter: function(params: any) {
          if (!params || !Array.isArray(params)) return '';
          const time = params[0]?.value[0]; // Get time from first series data point
          let tooltip = `Time: ${parseFloat(time).toPrecision(3)} s<br/>`;
          
          if (graphMode === 'error') {
            // Error graph tooltip
            params.forEach((param: any) => {
              if (param.seriesName === 'Control Error') {
                const value = parseFloat(param.value[1]).toPrecision(3);
                tooltip += `<span style="color: #3b82f6;">Control Error: ${value}</span><br/>`;
              }
            });
          } else if (graphMode === 'power') {
            // Controller Output graph tooltip
            params.forEach((param:any)=>{
              const value = parseFloat(param.value[1]).toPrecision(3);
              tooltip += `<span style="color: ${purpleColor};">Controller Output: ${value}</span><br/>`;
            });
          } else {
            // Main graph tooltip
            let pvValue = 0;
            let spValue = 0;
            
            params.forEach((param: any) => {
              const value = parseFloat(param.value[1]).toPrecision(3);
              // Map series names to their correct colors
              const color = param.seriesName === 'Process Variable' ? '#00ff51ff' : 
                           param.seriesName === 'Setpoint' ? '#ef4444' : param.color;
              
              // Skip the Control Error Range series in tooltip
              if (param.seriesName !== 'Control Error Range') {
                tooltip += `<span style="color: ${color};">${param.seriesName}: ${value}</span><br/>`;
                
                // Store values for error calculation
                if (param.seriesName === 'Process Variable') {
                  pvValue = parseFloat(param.value[1]);
                } else if (param.seriesName === 'Setpoint') {
                  spValue = parseFloat(param.value[1]);
                }
              }
            });
            
            // Calculate and display error
            const error = spValue - pvValue;
            tooltip += `<span style="color: #3b82f6;">Control Error: ${error.toFixed(2)}</span><br/>`;
          }
          
          return tooltip;
        },
        axisPointer: {
          type: 'cross',
          label: {
            show: true,
            backgroundColor: tooltipBg,
            color: isDark ? 'white' : '#000000',
            borderColor: tooltipBorder,
            borderWidth: 1,
            shadowBlur: 0,
            shadowColor: 'transparent',
            fontFamily: 'Merriweather Sans'
          },
          crossStyle: {
            color: '#999'
          }
        }
      },
      legend: {
        data: graphMode === 'error' 
          ? [{ name: 'Control Error', itemStyle: { color: '#3b82f6' } }]
          : graphMode === 'power'
          ? [{ name: 'Controller Output', itemStyle: { color: purpleColor } }]
          : [
              { name: 'Process Variable', itemStyle: { color: '#00ff51ff' } },
              { name: 'Setpoint', itemStyle: { color: '#ef4444' } }
            ],
        bottom: '2%',
        left: 'center',
        textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        itemWidth: 40,
        itemHeight: 3
      },
      series: [
        // Main series based on graph mode
        ...(graphMode === 'error' ? [
          {
            name: 'Control Error',
            type: 'line',
            data: simulationDataRef.current.map(p => [p.time, p.sp - p.pv]), // Populate with current error data
            showSymbol: false,
            lineStyle: { width: 3, color: '#ffffff' },
            smooth: true
          },
          {
            name: 'Zero Line',
            type: 'line',
            data: simulationDataRef.current.map(p => [p.time, 0]), // Zero reference line
            showSymbol: false,
            lineStyle: { type: 'dashed', color: '#888888', width: 1 },
            silent: true // Don't show in tooltip
          }
        ] : graphMode === 'power' ? [
          {
            name: 'Controller Output',
            type: 'line',
            data: simulationDataRef.current.map(p => [p.time, p.output]), // Controller output data
            showSymbol: false,
            lineStyle: { width: 3, color: purpleColor },
            smooth: true
          }
        ] : [
          {
            name: 'Process Variable',
            type: 'line',
            data: [], // Start with empty data - will be populated by simulation
            showSymbol: false,
            lineStyle: { width: 3, color: '#00ff51ff' },
            smooth: true
          },
          {
            name: 'Setpoint',
            type: 'line',
            data: [], // Start with empty data - will be populated by simulation  
            showSymbol: false,
            lineStyle: { type: 'dashed', color: '#ef4444', width: 2 }
          },
          {
            name: 'Control Error Range',
            type: 'custom',
            renderItem: (params: any, api: any) => {
              if (!simulationDataRef.current.length) return null;
              
              const currentData = simulationDataRef.current[simulationDataRef.current.length - 1];
              if (!currentData) return null;
              
              const currentTime = currentData.time;
              const currentPV = currentData.pv;
              const currentSP = currentData.sp;
              const error = currentSP - currentPV;
              
              // Only show if we have valid data and the current time is visible
              const timeRange = 30; // Current visible time window
              const currentVisibleTime = Math.max(0, currentTime - timeRange);
              
              if (currentTime < currentVisibleTime) return null;
              
              // Calculate positions
              const x = api.coord([currentTime, 0])[0];
              const y1 = api.coord([0, currentSP])[1];
              const y2 = api.coord([0, currentPV])[1];
              const yMid = (y1 + y2) / 2;
              
              // Theme-aware colors
              const isDark = resolvedTheme === 'dark';
              const lineColor = Math.abs(error) > 3 ? '#ef4444' : (isDark ? '#ffffff' : '#000000');
              const textColor = Math.abs(error) > 3 ? '#ef4444' : (isDark ? '#ffffff' : '#000000');
              
              return {
                type: 'group',
                children: [
                  // Vertical line showing the error range
                  {
                    type: 'line',
                    shape: {
                      x1: x,
                      y1: y1,
                      x2: x,
                      y2: y2
                    },
                    style: {
                      stroke: lineColor,
                      lineWidth: 2,
                      opacity: 0.7
                    }
                  },
                  // Top bracket
                  {
                    type: 'line',
                    shape: {
                      x1: x - 5,
                      y1: y1,
                      x2: x + 5,
                      y2: y1
                    },
                    style: {
                      stroke: lineColor,
                      lineWidth: 2
                    }
                  },
                  // Bottom bracket
                  {
                    type: 'line',
                    shape: {
                      x1: x - 5,
                      y1: y2,
                      x2: x + 5,
                      y2: y2
                    },
                    style: {
                      stroke: lineColor,
                      lineWidth: 2
                    }
                  },
                  // Error value text
                  {
                    type: 'text',
                    style: {
                      text: error.toFixed(2),
                      x: x + 10,
                      y: yMid,
                      textVerticalAlign: 'middle',
                      textAlign: 'left',
                      fill: textColor,
                      fontSize: 12,
                      fontFamily: 'monospace',
                      fontWeight: 'bold'
                    }
                  }
                ]
              };
            },
            data: [0] // Dummy data to trigger renderItem
          } as any
        ]),
        // Warning overlay (always present, shows for both graph modes)
        {
          name: 'Warning Overlay',
          type: 'custom',
          renderItem: (params: any, api: any) => {
            if (!divergenceWarningRef.current) return null;
            
            // Get the current visible x-axis range to properly center the warning
            const currentTime = simulationDataRef.current.length > 0 
              ? simulationDataRef.current[simulationDataRef.current.length - 1].time 
              : 50;
            const timeWindow = 30; // Match the timeWindow from the simulation
            const xMin = Math.max(0, currentTime - timeWindow * 0.7);
            const xMax = Math.max(timeWindow, currentTime + timeWindow * 0.3);
            const xMidValue = (xMin + xMax) / 2; // Center of current visible range
            
            // Use coordinate conversion to get center of the actual plot area
            const yMidValue = graphMode === 'error' ? 0 : 1; // Middle of y-axis range
            const centerCoords = api.coord([xMidValue, yMidValue]); // Use center of current visible x-range
            const centerX = centerCoords[0];
            const centerY = centerCoords[1];
            
            const isDark = resolvedTheme === 'dark';
            
            return {
              type: 'group',
              children: [
                // Warning background
                {
                  type: 'rect',
                  shape: {
                    x: centerX - 150,
                    y: centerY - 40,
                    width: 300,
                    height: 80
                  },
                  style: {
                    fill: isDark ? 'rgba(127, 29, 29, 0.95)' : 'rgba(254, 226, 226, 0.95)',
                    stroke: '#ef4444',
                    lineWidth: 2
                  }
                },
                // Warning title
                {
                  type: 'text',
                  style: {
                    text: '⚠ Warning: Poor Tuning Parameters!',
                    x: centerX,
                    y: centerY - 15,
                    textAlign: 'center',
                    textVerticalAlign: 'middle',
                    fill: isDark ? '#fecaca' : '#dc2626',
                    fontSize: 14,
                    fontWeight: 'bold',
                    fontFamily: 'Merriweather Sans'
                  }
                },
                // Countdown text
                {
                  type: 'text',
                  style: {
                    text: `Process diverging! Resetting in ${resetCountdownRef.current}s...`,
                    x: centerX,
                    y: centerY + 10,
                    textAlign: 'center',
                    textVerticalAlign: 'middle',
                    fill: isDark ? '#f87171' : '#b91c1c',
                    fontSize: 12,
                    fontFamily: 'Merriweather Sans'
                  }
                }
              ]
            };
          },
          data: [0], // Dummy data to trigger renderItem
          silent: true // Don't show in tooltip
        } as any
      ].filter(Boolean)
    };

  }, [!!result, controllerType, resolvedTheme, graphMode]); // Reduced dependencies to prevent flickering

  // Effect to handle graph switching between process and error views
  React.useEffect(() => {
    const echartsInstance = echartsRef.current?.getEchartsInstance();
    if (!echartsInstance || !simulationDataRef.current.length) return;

    // Force a complete chart re-render when switching graph types
    // This ensures the correct series structure is used
    echartsInstance.setOption(graphOptions, true); // true = not merge, replace completely
    
    // Then update with current data
    if (graphMode === 'bode') {
      // Bode plot doesn't need data updates - it's static based on calculated parameters
      // The graph options already contain all the necessary data
    } else if (graphMode === 'nyquist') {
      // Nyquist plot doesn't need data updates - it's static based on calculated parameters
      // The graph options already contain all the necessary data
    } else if (graphMode === 'error') {
      // Switch to error graph - provide complete series configurations
      const errorData = simulationDataRef.current.map(p => [p.time, p.sp - p.pv]);
      const zeroData = simulationDataRef.current.map(p => [p.time, 0]);
      echartsInstance.setOption({
        series: [
          {
            name: 'Control Error',
            type: 'line',
            showSymbol: false,
            lineStyle: { width: 3, color: '#3b82f6' },
            smooth: true,
            data: errorData,
          },
          {
            name: 'Zero Line',
            type: 'line',
            showSymbol: false,
            lineStyle: { type: 'dashed', color: '#888888', width: 1 },
            silent: true,
            data: zeroData,
          }
        ]
      });
    } else if (graphMode === 'power') {
      // Switch to power graph - provide complete series configurations
      echartsInstance.setOption({
        series: [
          {
            name: 'Controller Output',
            type: 'line',
            showSymbol: false,
            lineStyle: { width: 3, color: purpleColor },
            smooth: true,
            data: simulationDataRef.current.map(p => [p.time, p.output]),
          }
        ]
      });
    } else {
      // Switch to process graph - provide complete series configurations
      echartsInstance.setOption({
        series: [
          {
            name: 'Process Variable',
            type: 'line',
            showSymbol: false,
            lineStyle: { width: 3, color: '#00ff51ff' },
            smooth: true,
            data: simulationDataRef.current.map(p => [p.time, p.pv]),
          },
          {
            name: 'Setpoint',
            type: 'line',
            showSymbol: false,
            lineStyle: { type: 'dashed', color: '#ef4444', width: 2 },
            data: simulationDataRef.current.map(p => [p.time, p.sp]),
          },
          {
            name: 'Control Error Range',
            type: 'custom',
            renderItem: (params: any, api: any) => {
              if (!simulationDataRef.current.length) return null;
              const currentData = simulationDataRef.current[simulationDataRef.current.length - 1];
              if (!currentData) return null;
              const currentTime = currentData.time;
              const currentPV = currentData.pv;
              const currentSP = currentData.sp;
              const error = currentSP - currentPV;
              const timeRange = 30;
              const currentVisibleTime = Math.max(0, currentTime - timeRange);
              if (currentTime < currentVisibleTime) return null;
              const x = api.coord([currentTime, 0])[0];
              const y1 = api.coord([0, currentSP])[1];
              const y2 = api.coord([0, currentPV])[1];
              const yMid = (y1 + y2) / 2;
              const isDark = resolvedTheme === 'dark';
              const lineColor = Math.abs(error) > 3 ? '#ef4444' : (isDark ? '#ffffff' : '#000000');
              const textColor = Math.abs(error) > 3 ? '#ef4444' : (isDark ? '#ffffff' : '#000000');
              return {
                type: 'group',
                children: [
                  { type: 'line', shape: { x1: x, y1: y1, x2: x, y2: y2 }, style: { stroke: lineColor, lineWidth: 2, opacity: 0.7 } },
                  { type: 'line', shape: { x1: x - 5, y1: y1, x2: x + 5, y2: y1 }, style: { stroke: lineColor, lineWidth: 2 } },
                  { type: 'line', shape: { x1: x - 5, y1: y2, x2: x + 5, y2: y2 }, style: { stroke: lineColor, lineWidth: 2 } },
                  { type: 'text', style: { text: error.toFixed(2), x: x + 10, y: yMid, textVerticalAlign: 'middle', textAlign: 'left', fill: textColor, fontSize: 12, fontFamily: 'monospace', fontWeight: 'bold' } }
                ]
              };
            },
            data: [simulationDataRef.current[simulationDataRef.current.length - 1]?.time || 0] // Trigger custom series re-render
          }
        ]
      });
    }
  }, [graphMode, graphOptions]);

  // Effect to update Bode and Nyquist plots when frequency response data changes
  React.useEffect(() => {
    if (frequencyResponseData) {
      const echartsInstance = echartsRef.current?.getEchartsInstance();
      if (echartsInstance) {
        if (graphMode === 'bode') {
          // Complete series configuration for Bode plot
          echartsInstance.setOption({
            series: [
              {
                name: 'Magnitude',
                type: 'line',
                xAxisIndex: 0,
                yAxisIndex: 0,
                showSymbol: false,
                data: frequencyResponseData.bodeData.magnitude,
                lineStyle: { color: '#00ff51ff', width: 2 }
              },
              {
                name: 'Phase',
                type: 'line',
                xAxisIndex: 1,
                yAxisIndex: 1,
                showSymbol: false,
                data: frequencyResponseData.bodeData.phase,
                lineStyle: { color: '#3b82f6', width: 2 }
              }
            ]
          });
        } else if (graphMode === 'nyquist') {
          // Complete series configuration for Nyquist plot
          const isDark = resolvedTheme === 'dark';
          const textColor = isDark ? 'white' : '#000000';

          echartsInstance.setOption({
            series: [
              {
                name: 'Open Loop Response',
                type: 'line',
                showSymbol: false,
                data: frequencyResponseData.nyquistData,
                lineStyle: { color: '#00ff51ff', width: 2 },
                markPoint: {
                  data: [
                    // 1. Critical Point
                    {
                      name: 'Critical Point',
                      coord: [-1, 0],
                      symbol: 'cross',
                      symbolSize: 12,
                      label: {
                        show: true,
                        formatter: '(-1, 0)',
                        color: textColor,
                        fontFamily: 'Merriweather Sans',
                        fontSize: 14
                      }
                    },
                    // 2. Direction Arrows, regenerated with the new data
                    ...(() => {
                      const arrows = [];
                      const nyquistPoints = frequencyResponseData.nyquistData;
                      if (!nyquistPoints || nyquistPoints.length < 2) return [];

                      const numArrows = 6;
                      const startIndex = Math.floor(nyquistPoints.length * 0.1);
                      const endIndex = Math.floor(nyquistPoints.length * 0.9);
                      const step = Math.floor((endIndex - startIndex) / numArrows);

                      if (step <= 0) return [];

                      for (let i = startIndex; i < endIndex; i += step) {
                        const currentPoint = nyquistPoints[i];
                        const nextPoint = nyquistPoints[Math.min(i + 5, nyquistPoints.length - 1)];

                        if (currentPoint && nextPoint) {
                          const dx = nextPoint[0] - currentPoint[0];
                          const dy = nextPoint[1] - currentPoint[1];

                          if (Math.abs(dx) > 1e-6 || Math.abs(dy) > 1e-6) {
                            const angle = Math.atan2(dy, dx) * (180 / Math.PI);
                            arrows.push({
                              coord: currentPoint,
                              // FIX: Use the same custom SVG path here as well.
                              symbol: 'path://M 0 0 L -10 5 L -10 -5 Z',
                              symbolSize: 12,
                              symbolRotate: angle,
                              itemStyle: { color: '#00ff51ff' },
                              label: { show: false }
                            });
                          }
                        }
                      }
                      return arrows;
                    })()
                  ]
                }
              }
            ]
          });
        }
      }
    }
  }, [frequencyResponseData, graphMode, resolvedTheme]);

  // Cleanup simulation on component unmount
  useEffect(() => {
    return () => {
      if (simulationIntervalRef.current) {
        clearInterval(simulationIntervalRef.current);
      }
    };
  }, []);

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 sm:p-6 lg:p-8">
        {isMultiView && (multiViewData.length > 0 || simulationDataRef.current.length > 0) ? (
          // --- NEW: MULTI-GRAPH VIEW ---
          <div className="space-y-4">
            <div className="flex justify-between items-center">
              <Button variant="outline" onClick={() => setIsMultiView(false)}>
                ← Back to Single Graph View
              </Button>
              <Button variant="outline" onClick={() => setIsOverlappedView(!isOverlappedView)}>
                {isOverlappedView ? 'Separate Graphs' : 'Overlap Graphs'}
              </Button>
            </div>

            {isOverlappedView ? (
              // --- Overlapped View (Corrected) ---
              <div className="grid grid-cols-1 lg:grid-cols-5 gap-6">
                <div className="lg:col-span-2 space-y-6">
                  {/* PID Tuning Parameters Card */}
                  <Card>
                    <CardHeader className="pb-1">
                      <CardTitle>{controllerType} Tuning Parameters</CardTitle>
                    </CardHeader>
                    <CardContent className="space-y-6 p-4 pt-1">
                      <div className="space-y-4">
                        {/* Tuning Method Selector */}
                        <div className="flex items-center gap-2">
                            <Label htmlFor="tuningMethod" className="text-sm font-medium whitespace-nowrap">Tuning Method:</Label>
                            <Select value={tuningMethod} onValueChange={(value) => setTuningMethod(value as TuningMethod)}>
                                <SelectTrigger id="tuningMethod" className="flex-1"><SelectValue /></SelectTrigger>
                                <SelectContent className="max-h-[400px]">
                                    <SelectItem value="imc">IMC</SelectItem>
                                    <SelectItem value="itae">ITAE</SelectItem>
                                    <SelectItem value="amigo">AMIGO</SelectItem>
                                    <SelectItem value="ziegler-nichols">Ziegler-Nichols</SelectItem>
                                </SelectContent>
                            </Select>
                        </div>

                        {/* Controller Type Selector (hidden for IMC method) */}
                        {tuningMethod !== 'imc' && (
                          <div className="flex items-center gap-2">
                              <Label htmlFor="controllerType" className="text-sm font-medium whitespace-nowrap">Controller Type:</Label>
                              <Select value={controllerType} onValueChange={(value) => setControllerType(value as ControllerType)}>
                                  <SelectTrigger id="controllerType" className="flex-1"><SelectValue /></SelectTrigger>
                                  <SelectContent className="max-h-[400px]">
                                      {getAvailableControllerTypes(tuningMethod).map(type => (
                                        <SelectItem key={type} value={type}>{type}</SelectItem>
                                      ))}
                                  </SelectContent>
                              </Select>
                          </div>
                        )}

                        {/* Model Case Selector (only shown for IMC) */}
                        {tuningMethod === 'imc' && (
                          <div className="flex items-center gap-2">
                              <Label htmlFor="imcModelCase" className="text-sm font-medium whitespace-nowrap">Process Model:</Label>
                              <Select value={imcModelCase} onValueChange={(v) => setImcModelCase(v as ImcModelCase)}>
                                  <SelectTrigger id="imcModelCase" className="flex-1"><SelectValue /></SelectTrigger>
                                  <SelectContent className="max-h-[400px]">
                                      {/* Show all IMC cases A-O */}
                                      {(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'] as ImcModelCase[]).map(caseVal => {
                                        const descriptions = {
                                          'A': 'A: First-Order',
                                          'B': 'B: Second-Order Overdamped',
                                          'C': 'C: Second-Order Underdamped',
                                          'D': 'D: Second-Order Underdamped + RHP Zero',
                                          'E': 'E: Integrator',
                                          'F': 'F: First-Order + Integrator',
                                          'G': 'G: FOPTD',
                                          'H': 'H: FOPTD',
                                          'I': 'I: SOPTD Overdamped + LHP Zero',
                                          'J': 'J: SOPTD Underdamped + LHP Zero',
                                          'K': 'K: SOPTD Overdamped + RHP Zero',
                                          'L': 'L: SOPTD Underdamped + RHP Zero',
                                          'M': 'M: Integrator + Delay',
                                          'N': 'N: Integrator + Delay',
                                          'O': 'O: FOPTD + Integrator',
                                        };
                                        return (
                                          <SelectItem key={caseVal} value={caseVal}>
                                            {descriptions[caseVal]}
                                          </SelectItem>
                                        );
                                      })}
                                  </SelectContent>
                              </Select>
                          </div>
                        )}

                        {/* ITAE Input Type Selector (only shown for ITAE) */}
                        {tuningMethod === 'itae' && (
                          <div className="flex items-center gap-2">
                              <Label htmlFor="itaeInputType" className="text-sm font-medium whitespace-nowrap">Input Type:</Label>
                              <Select value={itaeInputType} onValueChange={(v) => setItaeInputType(v as ItaeInputType)}>
                                  <SelectTrigger id="itaeInputType" className="flex-1"><SelectValue /></SelectTrigger>
                                  <SelectContent className="max-h-[400px]">
                                      <SelectItem value="disturbance">Disturbance Rejection</SelectItem>
                                      <SelectItem value="setpoint">Set-Point Tracking</SelectItem>
                                  </SelectContent>
                              </Select>
                          </div>
                        )}
                        <hr />
                        {/* Inline Input Logic with Sliders */}
                        {tuningMethod === 'ziegler-nichols' && (
                          <>
                            <SliderGroup 
                                label={<span>Ultimate Gain (K<sub>cu</sub>):</span>} 
                                tooltip="The proportional gain at which the system has sustained oscillations."
                                value={kcu} onChange={setKcu} min={0.1} max={maxKcu} step={0.01}
                                maxValue={maxKcu} onMaxChange={setMaxKcu}
                            />
                            <SliderGroup 
                                label={<span>Ultimate Period (P<sub>u</sub>):</span>} 
                                tooltip="The period of the sustained oscillations at the ultimate gain (in time units)."
                                value={pu} onChange={setPu} min={0.1} max={maxPu} step={0.01}
                                maxValue={maxPu} onMaxChange={setMaxPu}
                            />
                          </>
                        )}

                        {tuningMethod === 'itae' && (
                           <>
                              <SliderGroup 
                                  label="Process Gain (K):" tooltip="Ratio of change in output to the change in input at steady state."
                                  value={k} onChange={setK} min={0.1} max={maxK} step={0.01}
                                  maxValue={maxK} onMaxChange={setMaxK}
                              />
                              <SliderGroup 
                                  label="Time Constant (τ):" tooltip="Time it takes for the process to reach 63.2% of its final value."
                                  value={tau} onChange={setTau} min={0.1} max={maxTau} step={0.01}
                                  maxValue={maxTau} onMaxChange={setMaxTau}
                              />
                              <SliderGroup 
                                  label="Dead Time (θ):" tooltip="Delay before the process output starts to respond to an input change."
                                  value={theta} onChange={setTheta} min={0.1} max={maxTheta} step={0.01}
                                  maxValue={maxTheta} onMaxChange={setMaxTheta}
                              />
                          </>
                        )}
                        
                        {tuningMethod === 'amigo' && (
                          <>
                            <SliderGroup 
                                label="Process Gain (K):" tooltip="For FOPTD, standard gain. For Integrator model, K is the integrator gain."
                                value={k} onChange={setK} min={0.1} max={maxK} step={0.01}
                                maxValue={maxK} onMaxChange={setMaxK}
                            />
                            <SliderGroup 
                                label="Time Constant (τ):" tooltip="Time constant for FOPTD model. Set to 0 for integrator model."
                                value={tau} onChange={setTau} min={0} max={maxTau} step={0.01}
                                maxValue={maxTau} onMaxChange={setMaxTau}
                            />
                            <SliderGroup 
                                label="Dead Time (θ):" tooltip="Delay before the process output starts to respond."
                                value={theta} onChange={setTheta} min={0.1} max={maxTheta} step={0.01}
                                maxValue={maxTheta} onMaxChange={setMaxTheta}
                            />
                          </>
                        )}

                        {tuningMethod === 'imc' && (() => {
                            const showTau = ['A', 'C', 'D', 'F', 'G', 'H', 'J', 'L', 'O'].includes(imcModelCase);
                            const showTau1 = ['B', 'I', 'K'].includes(imcModelCase);
                            const showTau2 = ['B', 'I', 'K'].includes(imcModelCase);
                            const showTau3 = ['I', 'J', 'K', 'L'].includes(imcModelCase);
                            const showZeta = ['C', 'D', 'J', 'L'].includes(imcModelCase);
                            const showTheta = ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'].includes(imcModelCase);
                            const showBeta = imcModelCase === 'D';
                            let tau3Label = <span>Zero Time Constant (τ₃):</span>;
                            if (['K', 'L'].includes(imcModelCase)) { tau3Label = <span>RHP Zero Time Constant (τ₃):</span>; }

                            return (
                                <>
                                    <SliderGroup label="Process Gain (K):" tooltip="Ratio of change in output to the change in input at steady state." value={k} onChange={setK} min={0.1} max={maxK} step={0.01} maxValue={maxK} onMaxChange={setMaxK} />
                                    {showTau && <SliderGroup label={<span>Time Constant (τ):</span>} tooltip="Primary time constant of the process model." value={tau} onChange={setTau} min={0.1} max={maxTau} step={0.01} maxValue={maxTau} onMaxChange={setMaxTau} />}
                                    {showTau1 && <SliderGroup label={<span>Time Constant 1 (τ₁):</span>} tooltip="First time constant of the second-order model." value={tau} onChange={setTau} min={0.1} max={maxTau} step={0.01} maxValue={maxTau} onMaxChange={setMaxTau} />}
                                    {showTau2 && <SliderGroup label={<span>Time Constant 2 (τ₂):</span>} tooltip="Second time constant of the second-order model." value={tau2} onChange={setTau2} min={0.1} max={maxTau2} step={0.01} maxValue={maxTau2} onMaxChange={setMaxTau2} />}
                                    {showTau3 && <SliderGroup label={tau3Label} tooltip="The time constant of the numerator zero term." value={tau3} onChange={setTau3} min={0.1} max={maxTau3} step={0.01} maxValue={maxTau3} onMaxChange={setMaxTau3} />}
                                    {showZeta && <SliderGroup label={<span>Damping Factor (ζ):</span>} tooltip="The damping factor of the second-order model." value={zeta} onChange={setZeta} min={0.05} max={maxZeta} step={0.01} maxValue={maxZeta} onMaxChange={setMaxZeta} />}
                                    {showBeta && <SliderGroup label={<span>RHP Zero (β):</span>} tooltip="The right-hand plane zero term. Must be positive." value={beta} onChange={setBeta} min={0.1} max={maxBeta} step={0.01} maxValue={maxBeta} onMaxChange={setMaxBeta} />}
                                    {showTheta && <SliderGroup label="Dead Time (θ):" tooltip="Delay before the process output starts to respond." value={theta} onChange={setTheta} min={0.1} max={maxTheta} step={0.01} maxValue={maxTheta} onMaxChange={setMaxTheta} />}
                                    <SliderGroup label={<span>IMC Tuning Parameter (τ<sub>c</sub>):</span>} tooltip="The desired closed-loop time constant. A larger value gives a slower, more robust response." value={tauC} onChange={setTauC} min={0.1} max={maxTauC} step={0.01} maxValue={maxTauC} onMaxChange={setMaxTauC} />
                                </>
                            );
                        })()}
                      </div>
                      {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
                    </CardContent>
                  </Card>
                  
                  {/* Simulation Control Card */}
                  <Card>
                    <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
                      <CardTitle>Simulation Control</CardTitle>
                      <Button variant="outline" size="sm" onClick={resetSimulation} disabled={!result || loading || divergenceWarning}>
                        Reset
                      </Button>
                    </CardHeader>
                    <CardContent className="space-y-4">
                      <div className="space-y-3">
                        <div className="flex items-center gap-4">
                          <div className="flex items-center gap-2">
                            <Label className="text-sm font-medium">Setpoint Value: {currentSetpoint}</Label>
                            <Button variant="outline" size="sm" onClick={() => handleSetpointChange('down')} disabled={currentSetpoint <= 0 || divergenceWarning} className="w-8 h-6 p-0 text-xs">-</Button>
                            <Button variant="outline" size="sm" onClick={() => handleSetpointChange('up')} disabled={currentSetpoint >= 2 || divergenceWarning} className="w-8 h-6 p-0 text-xs">+</Button>
                          </div>
                          <div className="flex items-center gap-2">
                            <Label className="text-sm font-medium">Process Disturbance:</Label>
                            <Button variant="outline" size="sm" onClick={() => triggerDisturbance('down')} disabled={!isSimulationRunning} className="w-12 h-6 p-0 text-xs">-0.5</Button>
                            <Button variant="outline" size="sm" onClick={() => triggerDisturbance('up')} disabled={!isSimulationRunning} className="w-12 h-6 p-0 text-xs">+0.5</Button>
                          </div>
                        </div>
                      </div>
                    </CardContent>
                  </Card>

                  {/* Controller Settings Card */}
                  <Card>
                    <CardHeader>
                      <CardTitle>Calculated Controller Settings</CardTitle>
                    </CardHeader>
                    <CardContent>
                      {loading && ( <div className="flex items-center justify-center text-muted-foreground p-8">Calculating...</div> )}
                      {!loading && !result && ( <div className="flex items-center justify-center text-muted-foreground p-8">Please provide inputs and start tuning.</div> )}
                      {result && !loading && (
                          <div className={`grid grid-cols-1 gap-4 ${controllerType === 'P' ? 'sm:grid-cols-1' : controllerType === 'PI' ? 'sm:grid-cols-2' : 'sm:grid-cols-3'}`}>
                              <Tooltip><TooltipTrigger asChild><div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2"><p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'K<sub>c</sub>:' }} /><p className="text-2xl font-bold font-mono">{formatResult(result.kc)}</p></div></TooltipTrigger><TooltipContent><p>Proportional Gain</p></TooltipContent></Tooltip>
                              {controllerType !== 'P' && (<Tooltip><TooltipTrigger asChild><div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2"><p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'τ<sub>I</sub>:' }} /><p className="text-2xl font-bold font-mono">{formatResult(result.tauI)}</p></div></TooltipTrigger><TooltipContent><p>Integral Time</p></TooltipContent></Tooltip>)}
                              {controllerType !== 'PI' && controllerType !== 'P' && (<Tooltip><TooltipTrigger asChild><div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2"><p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'τ<sub>D</sub>:' }} /><p className="text-2xl font-bold font-mono">{formatResult(result.tauD)}</p></div></TooltipTrigger><TooltipContent><p>Derivative Time</p></TooltipContent></Tooltip>)}
                          </div>
                      )}
                    </CardContent>
                  </Card>
                </div>

                <div className="lg:col-span-3">
                  <Card>
                    <CardContent className="py-2">
                      <div className="relative aspect-square rounded-md pt-0 px-2 pb-0">
                        <ReactECharts
                          ref={overlappedChartRef}
                          echarts={echarts}
                          option={overlappedGraphOptions}
                          style={{ height: '100%', width: '100%' }}
                        />
                      </div>
                    </CardContent>
                  </Card>
                </div>
              </div>
            ) : (
              // --- Separated View ---
              <>
                <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">
                  {/* Process Chart */}
                  <Card>
                    <CardContent className="p-2">
                      <div className="h-[400px]">
                        <ReactECharts
                          ref={multiViewProcessRef}
                          echarts={echarts}
                          option={multiViewProcessOptions}
                          style={{ height: '100%', width: '100%' }}
                        />
                      </div>
                    </CardContent>
                  </Card>

                  {/* Error Chart */}
                  <Card>
                    <CardContent className="p-2">
                      <div className="h-[400px]">
                        <ReactECharts
                          ref={multiViewErrorRef}
                          echarts={echarts}
                          option={multiViewErrorOptions}
                          style={{ height: '100%', width: '100%' }}
                        />
                      </div>
                    </CardContent>
                  </Card>
                  
                  {/* Controller Output Chart */}
                  <Card>
                    <CardContent className="p-2">
                      <div className="h-[400px]">
                        <ReactECharts
                          ref={multiViewOutputRef}
                          echarts={echarts}
                          option={multiViewOutputOptions}
                          style={{ height: '100%', width: '100%' }}
                        />
                      </div>
                    </CardContent>
                  </Card>
                </div>
                
                {/* Shared Simulation Controls */}
                <Card>
                  <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
                    <CardTitle>Simulation Control</CardTitle>
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={resetSimulation}
                      disabled={!result || loading || divergenceWarning}
                    >
                      Reset
                    </Button>
                  </CardHeader>
                  <CardContent className="space-y-4">
                    <div className="space-y-3">
                      {/* Setpoint Control */}
                      <div className="flex items-center gap-2">
                        <Label className="text-sm font-medium">Setpoint Value: {currentSetpoint}</Label>
                        <Button variant="outline" size="sm" onClick={() => handleSetpointChange('down')} disabled={currentSetpoint <= 0 || divergenceWarning} className="w-8 h-6 p-0 text-xs">-</Button>
                        <Button variant="outline" size="sm" onClick={() => handleSetpointChange('up')} disabled={currentSetpoint >= 2 || divergenceWarning} className="w-8 h-6 p-0 text-xs">+</Button>
                      </div>

                      {/* Process Disturbance Control */}
                      <div className="flex items-center gap-2">
                        <Label className="text-sm font-medium">Process Disturbance:</Label>
                        <Button variant="outline" size="sm" onClick={() => triggerDisturbance('down')} disabled={!isSimulationRunning} className="w-12 h-6 p-0 text-xs">-0.5</Button>
                        <Button variant="outline" size="sm" onClick={() => triggerDisturbance('up')} disabled={!isSimulationRunning} className="w-12 h-6 p-0 text-xs">+0.5</Button>
                      </div>
                    </div>
                  </CardContent>
                </Card>
              </>
            )}
          </div>
        ) : (
          // --- NEW SINGLE GRAPH VIEW STRUCTURE ---
          <div className="space-y-4">
            {/* 1. Consistent header bar */}
            <div className="flex justify-between items-center">
              <div></div>
              <Button variant="secondary" onClick={() => setIsMultiView(true)} disabled={!result}>
                View All Dynamic Graphs
              </Button>
            </div>

            {/* 2. Main grid layout */}
            <div className="grid grid-cols-1 lg:grid-cols-5 gap-6">
              {/* Column 1: Controls */}
              <div className="lg:col-span-2 space-y-6">
                {/* PID Tuning Parameters card */}
                <Card>
                  <CardHeader className="pb-1">
                    <CardTitle>{controllerType} Tuning Parameters</CardTitle>
                  </CardHeader>
                  <CardContent className="space-y-6 p-4 pt-1">
                    <div className="space-y-4">
                      {/* Tuning Method Selector */}
                      <div className="flex items-center gap-2">
                        <Label htmlFor="tuningMethod" className="text-sm font-medium whitespace-nowrap">Tuning Method:</Label>
                        <Select value={tuningMethod} onValueChange={(value) => setTuningMethod(value as TuningMethod)}>
                          <SelectTrigger id="tuningMethod" className="flex-1"><SelectValue /></SelectTrigger>
                          <SelectContent className="max-h-[400px]">
                            <SelectItem value="imc">IMC</SelectItem>
                            <SelectItem value="itae">ITAE</SelectItem>
                            <SelectItem value="amigo">AMIGO</SelectItem>
                            <SelectItem value="ziegler-nichols">Ziegler-Nichols</SelectItem>
                          </SelectContent>
                        </Select>
                      </div>
                      {/* Controller Type Selector (hidden for IMC method) */}
                      {tuningMethod !== 'imc' && (
                        <div className="flex items-center gap-2">
                          <Label htmlFor="controllerType" className="text-sm font-medium whitespace-nowrap">Controller Type:</Label>
                          <Select value={controllerType} onValueChange={(value) => setControllerType(value as ControllerType)}>
                            <SelectTrigger id="controllerType" className="flex-1"><SelectValue /></SelectTrigger>
                            <SelectContent className="max-h-[400px]">
                              {getAvailableControllerTypes(tuningMethod).map(type => (
                                <SelectItem key={type} value={type}>{type}</SelectItem>
                              ))}
                            </SelectContent>
                          </Select>
                        </div>
                      )}
                      {/* Model Case Selector (only shown for IMC) */}
                      {tuningMethod === 'imc' && (
                        <div className="flex items-center gap-2">
                          <Label htmlFor="imcModelCase" className="text-sm font-medium whitespace-nowrap">Process Model:</Label>
                          <Select value={imcModelCase} onValueChange={(v) => setImcModelCase(v as ImcModelCase)}>
                            <SelectTrigger id="imcModelCase" className="flex-1"><SelectValue /></SelectTrigger>
                            <SelectContent className="max-h-[400px]">
                              {(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'] as ImcModelCase[]).map(caseVal => {
                                const descriptions = {
                                  'A': 'A: First-Order',
                                  'B': 'B: Second-Order Overdamped',
                                  'C': 'C: Second-Order Underdamped',
                                  'D': 'D: Second-Order Underdamped + RHP Zero',
                                  'E': 'E: Integrator',
                                  'F': 'F: First-Order + Integrator',
                                  'G': 'G: FOPTD',
                                  'H': 'H: FOPTD',
                                  'I': 'I: SOPTD Overdamped + LHP Zero',
                                  'J': 'J: SOPTD Underdamped + LHP Zero',
                                  'K': 'K: SOPTD Overdamped + RHP Zero',
                                  'L': 'L: SOPTD Underdamped + RHP Zero',
                                  'M': 'M: Integrator + Delay',
                                  'N': 'N: Integrator + Delay',
                                  'O': 'O: FOPTD + Integrator',
                                };
                                return (
                                  <SelectItem key={caseVal} value={caseVal}>
                                    {descriptions[caseVal]}
                                  </SelectItem>
                                );
                              })}
                            </SelectContent>
                          </Select>
                        </div>
                      )}
                      {/* ITAE Input Type Selector (only shown for ITAE) */}
                      {tuningMethod === 'itae' && (
                        <div className="flex items-center gap-2">
                          <Label htmlFor="itaeInputType" className="text-sm font-medium whitespace-nowrap">Input Type:</Label>
                          <Select value={itaeInputType} onValueChange={(v) => setItaeInputType(v as ItaeInputType)}>
                            <SelectTrigger id="itaeInputType" className="flex-1"><SelectValue /></SelectTrigger>
                            <SelectContent className="max-h-[400px]">
                              <SelectItem value="disturbance">Disturbance Rejection</SelectItem>
                              <SelectItem value="setpoint">Set-Point Tracking</SelectItem>
                            </SelectContent>
                          </Select>
                        </div>
                      )}
                      <hr />
                      {/* Inline Input Logic */}
                      {tuningMethod === 'ziegler-nichols' && (
                        <>
                          <SliderGroup label={<span>Ultimate Gain (K<sub>cu</sub>):</span>} tooltip="The proportional gain at which the system has sustained oscillations." value={kcu} onChange={setKcu} min={0.1} max={maxKcu} step={0.01} maxValue={maxKcu} onMaxChange={setMaxKcu} />
                          <SliderGroup label={<span>Ultimate Period (P<sub>u</sub>):</span>} tooltip="The period of the sustained oscillations at the ultimate gain (in time units)." value={pu} onChange={setPu} min={0.1} max={maxPu} step={0.01} maxValue={maxPu} onMaxChange={setMaxPu} />
                        </>
                      )}
                      {tuningMethod === 'itae' && (
                        <>
                          <SliderGroup label="Process Gain (K):" tooltip="Ratio of change in output to the change in input at steady state." value={k} onChange={setK} min={0.1} max={maxK} step={0.01} maxValue={maxK} onMaxChange={setMaxK} />
                          <SliderGroup label="Time Constant (τ):" tooltip="Time it takes for the process to reach 63.2% of its final value." value={tau} onChange={setTau} min={0.1} max={maxTau} step={0.01} maxValue={maxTau} onMaxChange={setMaxTau} />
                          <SliderGroup label="Dead Time (θ):" tooltip="Delay before the process output starts to respond to an input change." value={theta} onChange={setTheta} min={0.1} max={maxTheta} step={0.01} maxValue={maxTheta} onMaxChange={setMaxTheta} />
                        </>
                      )}
                      {tuningMethod === 'amigo' && (
                        <>
                          <SliderGroup label="Process Gain (K):" tooltip="For FOPTD, standard gain. For Integrator model, K is the integrator gain." value={k} onChange={setK} min={0.1} max={maxK} step={0.01} maxValue={maxK} onMaxChange={setMaxK} />
                          <SliderGroup label="Time Constant (τ):" tooltip="Time constant for FOPTD model. Set to 0 for integrator model." value={tau} onChange={setTau} min={0} max={maxTau} step={0.01} maxValue={maxTau} onMaxChange={setMaxTau} />
                          <SliderGroup label="Dead Time (θ):" tooltip="Delay before the process output starts to respond." value={theta} onChange={setTheta} min={0.1} max={maxTheta} step={0.01} maxValue={maxTheta} onMaxChange={setMaxTheta} />
                        </>
                      )}
                      {tuningMethod === 'imc' && (() => {
                        const showTau = ['A', 'C', 'D', 'F', 'G', 'H', 'J', 'L', 'O'].includes(imcModelCase);
                        const showTau1 = ['B', 'I', 'K'].includes(imcModelCase);
                        const showTau2 = ['B', 'I', 'K'].includes(imcModelCase);
                        const showTau3 = ['I', 'J', 'K', 'L'].includes(imcModelCase);
                        const showZeta = ['C', 'D', 'J', 'L'].includes(imcModelCase);
                        const showTheta = ['G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O'].includes(imcModelCase);
                        const showBeta = imcModelCase === 'D';
                        let tau3Label = <span>Zero Time Constant (τ₃):</span>;
                        if (['K', 'L'].includes(imcModelCase)) { tau3Label = <span>RHP Zero Time Constant (τ₃):</span>; }
                        return (
                          <>
                            <SliderGroup label="Process Gain (K):" tooltip="Ratio of change in output to the change in input at steady state." value={k} onChange={setK} min={0.1} max={maxK} step={0.01} maxValue={maxK} onMaxChange={setMaxK} />
                            {showTau && <SliderGroup label={<span>Time Constant (τ):</span>} tooltip="Primary time constant of the process model." value={tau} onChange={setTau} min={0.1} max={maxTau} step={0.01} maxValue={maxTau} onMaxChange={setMaxTau} />}
                            {showTau1 && <SliderGroup label={<span>Time Constant 1 (τ₁):</span>} tooltip="First time constant of the second-order model." value={tau} onChange={setTau} min={0.1} max={maxTau} step={0.01} maxValue={maxTau} onMaxChange={setMaxTau} />}
                            {showTau2 && <SliderGroup label={<span>Time Constant 2 (τ₂):</span>} tooltip="Second time constant of the second-order model." value={tau2} onChange={setTau2} min={0.1} max={maxTau2} step={0.01} maxValue={maxTau2} onMaxChange={setMaxTau2} />}
                            {showTau3 && <SliderGroup label={tau3Label} tooltip="The time constant of the numerator zero term." value={tau3} onChange={setTau3} min={0.1} max={maxTau3} step={0.01} maxValue={maxTau3} onMaxChange={setMaxTau3} />}
                            {showZeta && <SliderGroup label={<span>Damping Factor (ζ):</span>} tooltip="The damping factor of the second-order model." value={zeta} onChange={setZeta} min={0.05} max={maxZeta} step={0.01} maxValue={maxZeta} onMaxChange={setMaxZeta} />}
                            {showBeta && <SliderGroup label={<span>RHP Zero (β):</span>} tooltip="The right-hand plane zero term. Must be positive." value={beta} onChange={setBeta} min={0.1} max={maxBeta} step={0.01} maxValue={maxBeta} onMaxChange={setMaxBeta} />}
                            {showTheta && <SliderGroup label="Dead Time (θ):" tooltip="Delay before the process output starts to respond." value={theta} onChange={setTheta} min={0.1} max={maxTheta} step={0.01} maxValue={maxTheta} onMaxChange={setMaxTheta} />}
                            <SliderGroup label={<span>IMC Tuning Parameter (τ<sub>c</sub>):</span>} tooltip="The desired closed-loop time constant. A larger value gives a slower, more robust response." value={tauC} onChange={setTauC} min={0.1} max={maxTauC} step={0.01} maxValue={maxTauC} onMaxChange={setMaxTauC} />
                          </>
                        );
                      })()}
                    </div>
                    {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
                  </CardContent>
                </Card>
                {/* Simulation Control card */}
                <Card>
                  <CardHeader className="flex flex-row items-center justify-between space-y-0 pb-2">
                    <CardTitle>Simulation Control</CardTitle>
                    <Button variant="outline" size="sm" onClick={resetSimulation} disabled={!result || loading || divergenceWarning}>
                      Reset
                    </Button>
                  </CardHeader>
                  <CardContent className="space-y-4">
                    <div className="space-y-3">
                      <div className="flex items-center gap-4">
                        <div className="flex items-center gap-2">
                          <Label className="text-sm font-medium">Setpoint Value: {currentSetpoint}</Label>
                          <Button variant="outline" size="sm" onClick={() => handleSetpointChange('down')} disabled={currentSetpoint <= 0 || divergenceWarning} className="w-8 h-6 p-0 text-xs">-</Button>
                          <Button variant="outline" size="sm" onClick={() => handleSetpointChange('up')} disabled={currentSetpoint >= 2 || divergenceWarning} className="w-8 h-6 p-0 text-xs">+</Button>
                        </div>
                        <div className="flex items-center gap-2">
                          <Label className="text-sm font-medium">Process Disturbance:</Label>
                          <Button variant="outline" size="sm" onClick={() => triggerDisturbance('down')} disabled={!isSimulationRunning} className="w-12 h-6 p-0 text-xs">-0.5</Button>
                          <Button variant="outline" size="sm" onClick={() => triggerDisturbance('up')} disabled={!isSimulationRunning} className="w-12 h-6 p-0 text-xs">+0.5</Button>
                        </div>
                      </div>
                    </div>
                  </CardContent>
                </Card>
                {/* Controller Settings card */}
                <Card>
                  <CardHeader>
                    <CardTitle>Calculated Controller Settings</CardTitle>
                  </CardHeader>
                  <CardContent>
                    {loading && ( <div className="flex items-center justify-center text-muted-foreground p-8">Calculating...</div> )}
                    {!loading && !result && ( <div className="flex items-center justify-center text-muted-foreground p-8">Please provide inputs and start tuning.</div> )}
                    {result && !loading && (
                      <div className={`grid grid-cols-1 gap-4 ${controllerType === 'P' ? 'sm:grid-cols-1' : controllerType === 'PI' ? 'sm:grid-cols-2' : 'sm:grid-cols-3'}`}>
                        <Tooltip><TooltipTrigger asChild><div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2"><p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'K<sub>c</sub>:' }} /><p className="text-2xl font-bold font-mono">{formatResult(result.kc)}</p></div></TooltipTrigger><TooltipContent><p>Proportional Gain</p></TooltipContent></Tooltip>
                        {controllerType !== 'P' && (<Tooltip><TooltipTrigger asChild><div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2"><p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'τ<sub>I</sub>:' }} /><p className="text-2xl font-bold font-mono">{formatResult(result.tauI)}</p></div></TooltipTrigger><TooltipContent><p>Integral Time</p></TooltipContent></Tooltip>)}
                        {controllerType !== 'PI' && controllerType !== 'P' && (<Tooltip><TooltipTrigger asChild><div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2"><p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'τ<sub>D</sub>:' }} /><p className="text-2xl font-bold font-mono">{formatResult(result.tauD)}</p></div></TooltipTrigger><TooltipContent><p>Derivative Time</p></TooltipContent></Tooltip>)}
                      </div>
                    )}
                  </CardContent>
                </Card>
              </div>

              {/* Column 2: Graph */}
              <div className="lg:col-span-3">
                <Card>
                  <CardContent className="py-2">
                    <div className="relative aspect-square rounded-md pt-0 px-2 pb-0">
                      <ReactECharts
                        ref={echartsRef}
                        echarts={echarts}
                        option={graphOptions}
                        style={{ height: '100%', width: '100%' }}
                      />
                    </div>
                  </CardContent>
                </Card>
              </div>
            </div>
          </div>
        )}
      </div>
    </TooltipProvider>
  );
}