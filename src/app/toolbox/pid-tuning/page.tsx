'use client';

import React, { useState, useMemo, useCallback, memo, useRef } from 'react';
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
  
  // Input States
  const [tuningMethod, setTuningMethod] = useState<TuningMethod>('imc');
  const [controllerType, setControllerType] = useState<ControllerType>('PI');
  
  // Model Parameter String States
  const [k, setK] = useState('1.0');
  const [tau, setTau] = useState('1.0');
  const [tau2, setTau2] = useState('5.0'); // For Case B
  const [tau3, setTau3] = useState('2.0'); // For Case I
  const [beta, setBeta] = useState('1.0'); // For Case D
  const [theta, setTheta] = useState('2.0');
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
  
  // Debounce timer for calculations
  const debounceTimer = React.useRef<NodeJS.Timeout | null>(null);
  const echartsRef = useRef<ReactECharts | null>(null); // Ref for ECharts instance

  // Simulation function for closed-loop response
  const simulateClosedLoopResponse = (
    processParams: { K: number; tau: number; theta: number; tau2?: number; zeta?: number },
    controller: { Kc: number; tauI: number | null; tauD: number | null },
    modelCase: ImcModelCase
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
      const setpoint = time > dt ? 1 : 0;
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
      const delayed_output = delayBuffer.shift() || 0;
      delayBuffer.push(controller_output);

      // Select the correct process model for simulation
      if (modelCase === 'B' || modelCase === 'K' || modelCase === 'I') {
        // Second-Order Overdamped (Two Lags in series)
        const tau1_sim = tau;
        const tau2_sim = tau2 || tau; // Fallback if tau2 not provided
        const y1_deriv = (K * delayed_output - y1) / tau1_sim;
        y1 += y1_deriv * dt;
        const y2_deriv = (y1 - y2) / tau2_sim;
        y2 += y2_deriv * dt;
        pv = y2;
      } else {
        // Default to FOPTD for all other relevant cases (G, H, ITAE, AMIGO etc.)
        const pv_deriv = tau > 0 ? (K * delayed_output - pv) / tau : 0;
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

  // Effect to run calculation when inputs change (debounced)
  React.useEffect(() => {
    if (debounceTimer.current) {
      clearTimeout(debounceTimer.current);
    }
    debounceTimer.current = setTimeout(() => {
      handleCalculate();
    }, 1); // 1ms debounce delay - same as reactor design
  }, [
    k, tau, tau2, tau3, beta, theta, kcu, pu, tauC, zeta,
    tuningMethod, controllerType, itaeInputType, imcModelCase,
    handleCalculate
  ]);

  // Update controller type when tuning method changes
  React.useEffect(() => {
    const availableTypes = getAvailableControllerTypes(tuningMethod);
    if (!availableTypes.includes(controllerType)) {
      setControllerType(availableTypes[availableTypes.length - 1]);
    }
  }, [tuningMethod]);

  // Update IMC model case when controller type changes for IMC method
  React.useEffect(() => {
    if (tuningMethod === 'imc') {
      const availableCases = getAvailableImcCases(controllerType);
      if (!availableCases.includes(imcModelCase)) {
        // Set to a default compatible case
        setImcModelCase(availableCases[0]);
      }
    }
  }, [controllerType, tuningMethod]);

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

  const formatResult = (value: number | null) => {
    if (value === null || isNaN(value)) return '---';
    if (value === 0) return '0.00';
    
    // Always format to exactly 3 significant figures
    return value.toPrecision(3);
  }

  // Graph options for ECharts
  const graphOptions = useMemo((): EChartsOption => {
    // Only run simulation if we have calculated results
    if (!result) {
      console.log('Graph not showing - no results:', { result });
      return {};
    }

    // --- FIX: Assemble the correct process parameters based on the selected model ---
    const processParams = {
      K: parseFloat(k),
      tau: parseFloat(tau),
      theta: parseFloat(theta),
      tau2: parseFloat(tau2), // for Case B, I, K
      zeta: parseFloat(zeta)  // for Case C, D, H, etc.
    };

    const controllerParams = {
      Kc: result.kc || 0,
      tauI: result.tauI,
      tauD: result.tauD
    };

    console.log('Process params:', processParams);
    console.log('Controller params:', controllerParams);

    // Check if parameters are valid for simulation
    if (isNaN(processParams.K) || isNaN(processParams.tau) || isNaN(processParams.theta) || !controllerParams.Kc) {
      console.log('Invalid parameters for simulation');
      return {};
    }

    // Determine which IMC model to use for the simulation.
    // Default to 'G' (FOPTD) if not IMC, as it's the most common.
    const simModelCase = tuningMethod === 'imc' ? imcModelCase : 'G';

    const simulationData = simulateClosedLoopResponse(processParams, controllerParams, simModelCase);
    console.log('Simulation data points:', simulationData.length);

    // Theme-dependent colors
    const isDark = resolvedTheme === 'dark';
    const textColor = isDark ? 'white' : '#000000';
    const tooltipBg = isDark ? '#08306b' : '#f8f9fa';
    const tooltipBorder = isDark ? '#55aaff' : '#dee2e6';

    return {
      backgroundColor: 'transparent',
      title: { 
        text: `${controllerType} Controller Step Response`, 
        left: 'center',
        top: '0%',
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' } 
      },
      grid: { left: '5%', right: '5%', bottom: '12%', top: '5%', containLabel: true },
      xAxis: {
        type: 'value',
        name: 'Time (s)',
        nameLocation: 'middle',
        nameGap: 30,
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        axisLine: { lineStyle: { color: textColor } },
        axisTick: { lineStyle: { color: textColor }, length: 5, inside: false },
        axisLabel: { color: textColor, fontSize: 16, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false }
      },
      yAxis: {
        type: 'value',
        name: 'Process Variable',
        nameLocation: 'middle',
        nameGap: 40,
        nameTextStyle: { color: textColor, fontSize: 15, fontFamily: 'Merriweather Sans' },
        min: 0,
        max: 2,
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
          const time = params[0]?.axisValue;
          let tooltip = `Time: ${parseFloat(time).toPrecision(3)} s<br/>`;
          params.forEach((param: any) => {
            const value = parseFloat(param.value[1]).toPrecision(3);
            // Map series names to their correct colors
            const color = param.seriesName === 'Process Variable' ? '#3b82f6' : 
                         param.seriesName === 'Setpoint' ? '#ef4444' : param.color;
            tooltip += `<span style="color: ${color};">${param.seriesName}: ${value}</span><br/>`;
          });
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
        data: [
          { name: 'Process Variable', itemStyle: { color: '#3b82f6' } },
          { name: 'Setpoint', itemStyle: { color: '#ef4444' } }
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
          data: simulationData.map(p => [p.time, p.pv]),
          showSymbol: false,
          lineStyle: { width: 3, color: '#3b82f6' },
          smooth: true
        },
        {
          name: 'Setpoint',
          type: 'line',
          data: simulationData.map(p => [p.time, p.sp]),
          showSymbol: false,
          lineStyle: { type: 'dashed', color: '#ef4444', width: 2 }
        }
      ]
    };

  }, [result, k, tau, theta, tau2, zeta, imcModelCase, tuningMethod, controllerType, resolvedTheme]);

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-6 px-4 md:px-16">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Column 1: Controls (50% width on large screens) */}
          <div className="space-y-6">
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

                  {/* Controller Type Selector */}
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

                  {/* Model Case Selector (only shown for IMC) */}
                  {tuningMethod === 'imc' && (
                    <div className="flex items-center gap-2">
                        <Label htmlFor="imcModelCase" className="text-sm font-medium whitespace-nowrap">Model Case:</Label>
                        <Select value={imcModelCase} onValueChange={(v) => setImcModelCase(v as ImcModelCase)}>
                            <SelectTrigger id="imcModelCase" className="flex-1"><SelectValue /></SelectTrigger>
                            <SelectContent className="max-h-[400px]">
                                {getAvailableImcCases(controllerType).map(caseVal => {
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
                      <SliderGroup 
                          label={<span>Ultimate Gain (K<sub>cu</sub>):</span>} 
                          tooltip="The proportional gain at which the system has sustained oscillations."
                          value={kcu}
                          onChange={setKcu}
                          min={0.1}
                          max={maxKcu}
                          step={0.01}
                          maxValue={maxKcu}
                          onMaxChange={setMaxKcu}
                      />
                      <SliderGroup 
                          label={<span>Ultimate Period (P<sub>u</sub>):</span>} 
                          tooltip="The period of the sustained oscillations at the ultimate gain (in time units)."
                          value={pu}
                          onChange={setPu}
                          min={0.1}
                          max={maxPu}
                          step={0.01}
                          maxValue={maxPu}
                          onMaxChange={setMaxPu}
                      />
                    </>
                  )}

                  {tuningMethod === 'itae' && (
                     <>
                        <SliderGroup 
                            label="Process Gain (K):" 
                            tooltip="Ratio of change in output to the change in input at steady state."
                            value={k}
                            onChange={setK}
                            min={0.1}
                            max={maxK}
                            step={0.01}
                            maxValue={maxK}
                            onMaxChange={setMaxK}
                        />
                        <SliderGroup 
                            label="Time Constant (τ):" 
                            tooltip="Time it takes for the process to reach 63.2% of its final value."
                            value={tau}
                            onChange={setTau}
                            min={0.1}
                            max={maxTau}
                            step={0.01}
                            maxValue={maxTau}
                            onMaxChange={setMaxTau}
                        />
                        <SliderGroup 
                            label="Dead Time (θ):" 
                            tooltip="Delay before the process output starts to respond to an input change."
                            value={theta}
                            onChange={setTheta}
                            min={0.1}
                            max={maxTheta}
                            step={0.01}
                            maxValue={maxTheta}
                            onMaxChange={setMaxTheta}
                        />
                    </>
                  )}
                  
                  {tuningMethod === 'amigo' && (
                    <>
                      <SliderGroup 
                          label="Process Gain (K):" 
                          tooltip="For FOPTD, standard gain. For Integrator model, K is the integrator gain."
                          value={k}
                          onChange={setK}
                          min={0.1}
                          max={maxK}
                          step={0.01}
                          maxValue={maxK}
                          onMaxChange={setMaxK}
                      />
                      <SliderGroup 
                          label="Time Constant (τ):" 
                          tooltip="Time constant for FOPTD model. Set to 0 for integrator model."
                          value={tau}
                          onChange={setTau}
                          min={0}
                          max={maxTau}
                          step={0.01}
                          maxValue={maxTau}
                          onMaxChange={setMaxTau}
                      />
                      <SliderGroup 
                          label="Dead Time (θ):" 
                          tooltip="Delay before the process output starts to respond."
                          value={theta}
                          onChange={setTheta}
                          min={0.1}
                          max={maxTheta}
                          step={0.01}
                          maxValue={maxTheta}
                          onMaxChange={setMaxTheta}
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
                      if (['K', 'L'].includes(imcModelCase)) {
                          tau3Label = <span>RHP Zero Time Constant (τ₃):</span>;
                      }

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
          </div>
           {/* Column 2: Results (50% width on large screens) */}
          <div className="space-y-6">
            {/* Graph Card */}
            <Card>
              <CardContent className="p-0 pb-0">
                <div className="aspect-square rounded-md pt-0 px-2 pb-0">
                   {Object.keys(graphOptions).length > 0 ? (
                     <ReactECharts 
                       ref={echartsRef} 
                       echarts={echarts} 
                       option={graphOptions} 
                       style={{ height: '100%', width: '100%', borderRadius: '0.375rem', overflow: 'hidden' }} 
                       notMerge={false} 
                       lazyUpdate={false} 
                     />
                   ) : (
                     <div className="flex items-center justify-center h-full text-muted-foreground">
                       Adjust parameters above to see controller response simulation
                     </div>
                   )}
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
                        <Tooltip>
                            <TooltipTrigger asChild>
                                <div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2">
                                    <p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'K<sub>c</sub>:' }} />
                                    <p className="text-2xl font-bold font-mono">{formatResult(result.kc)}</p>
                                </div>
                            </TooltipTrigger>
                            <TooltipContent><p>Proportional Gain</p></TooltipContent>
                        </Tooltip>
                        {controllerType !== 'P' && (
                            <Tooltip>
                                <TooltipTrigger asChild>
                                    <div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2">
                                        <p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'τ<sub>I</sub>:' }} />
                                        <p className="text-2xl font-bold font-mono">{formatResult(result.tauI)}</p>
                                    </div>
                                </TooltipTrigger>
                                <TooltipContent><p>Integral Time</p></TooltipContent>
                            </Tooltip>
                        )}
                        {controllerType !== 'PI' && controllerType !== 'P' && (
                            <Tooltip>
                                <TooltipTrigger asChild>
                                    <div className="p-4 bg-muted rounded-md cursor-help flex items-center justify-center gap-2">
                                        <p className="text-2xl font-medium" dangerouslySetInnerHTML={{ __html: 'τ<sub>D</sub>:' }} />
                                        <p className="text-2xl font-bold font-mono">{formatResult(result.tauD)}</p>
                                    </div>
                                </TooltipTrigger>
                                <TooltipContent><p>Derivative Time</p></TooltipContent>
                            </Tooltip>
                        )}
                    </div>
                )}
              </CardContent>
            </Card>
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}