'use client';

import React, { useState, useMemo } from 'react';

import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Tooltip, TooltipContent, TooltipProvider, TooltipTrigger } from "@/components/ui/tooltip";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Info } from 'lucide-react';

// --- Type Definitions ---
type TuningMethod = 'ziegler-nichols' | 'itae' | 'amigo' | 'imc';
type ControllerType = 'P' | 'PI' | 'PID';
type ItaeInputType = 'disturbance' | 'setpoint';
type ImcModelCase = 'G' | 'H'; // Starting with common FOPTD models from the table

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

// --- Main Component ---
export default function PidTuningPage() {
  // Input States
  const [tuningMethod, setTuningMethod] = useState<TuningMethod>('imc');
  const [controllerType, setControllerType] = useState<ControllerType>('PID');
  
  // Model Parameter String States
  const [k, setK] = useState('1.0');
  const [tau, setTau] = useState('10.0');
  const [theta, setTheta] = useState('2.0');
  const [kcu, setKcu] = useState('2.0');
  const [pu, setPu] = useState('5.0');
  const [tauC, setTauC] = useState('2.0'); // For IMC method

  // Specific Method States
  const [itaeInputType, setItaeInputType] = useState<ItaeInputType>('disturbance');
  const [imcModelCase, setImcModelCase] = useState<ImcModelCase>('G');

  // Control & Result States
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [result, setResult] = useState<TuningResult | null>(null);
  const [displayedMethod, setDisplayedMethod] = useState<string>('');

  // Initial calculation on mount
  React.useEffect(() => {
    handleCalculate();
  }, []);

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
        return ['PID'];
      default:
        return ['P', 'PI', 'PID'];
    }
  };

  // Update controller type when tuning method changes
  React.useEffect(() => {
    const availableTypes = getAvailableControllerTypes(tuningMethod);
    if (!availableTypes.includes(controllerType)) {
      setControllerType(availableTypes[availableTypes.length - 1]);
    }
  }, [tuningMethod]);

  // --- Calculation Logic ---
  const handleCalculate = () => {
    setLoading(true);
    setError(null);
    setResult(null);

    // Simulate calculation delay
    setTimeout(() => {
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
                    setDisplayedMethod(`IMC for Model Case ${imcModelCase}`);
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
    }, 500); // 500ms delay to show loading state
  };
  
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
      
      if (ct !== 'PID') throw new Error("IMC rules from this table are for PID controllers.");

      switch (mcase) {
          case 'G':
              kc = (1/k_val) * ( (tau_val + 0.5 * theta_val) / (tauc_val + 0.5 * theta_val) );
              tauI = tau_val + 0.5 * theta_val;
              tauD = (tau_val * theta_val) / (2 * tau_val + theta_val);
              break;
          
          case 'H':
              kc = (1/k_val) * ( (theta_val/2) / (tauc_val + theta_val/2) );
              tauI = theta_val/2;
              tauD = 0;
              break;
          default:
              throw new Error(`IMC model case ${mcase} not implemented.`);
      }
      return { kc, tauI, tauD };
  };

  // --- Render Helper for Inputs ---
  const renderInputs = () => {
    switch(tuningMethod) {
        case 'ziegler-nichols':
            return (
                <>
                    <InputGroup label="Ultimate Gain (K_cu)" tooltip="The proportional gain at which the system has sustained oscillations.">
                        <Input id="kcu" value={kcu} onChange={(e) => setKcu(e.target.value)} placeholder="e.g., 2.0" />
                    </InputGroup>
                    <InputGroup label="Ultimate Period (P_u)" tooltip="The period of the sustained oscillations at the ultimate gain (in time units).">
                        <Input id="pu" value={pu} onChange={(e) => setPu(e.target.value)} placeholder="e.g., 5.0" />
                    </InputGroup>
                </>
            );
        case 'itae':
            return (
                <>
                    <InputGroup label="Process Gain (K)" tooltip="Ratio of change in output to the change in input at steady state.">
                        <Input id="k" value={k} onChange={(e) => setK(e.target.value)} placeholder="e.g., 1.0" />
                    </InputGroup>
                    <InputGroup label="Time Constant (τ)" tooltip="Time it takes for the process to reach 63.2% of its final value.">
                        <Input id="tau" value={tau} onChange={(e) => setTau(e.target.value)} placeholder="e.g., 10.0" />
                    </InputGroup>
                    <InputGroup label="Dead Time (θ)" tooltip="Delay before the process output starts to respond to an input change.">
                        <Input id="theta" value={theta} onChange={(e) => setTheta(e.target.value)} placeholder="e.g., 2.0" />
                    </InputGroup>
                    <div className="flex items-center gap-2">
                        <Label htmlFor="itaeInputType" className="text-sm font-medium whitespace-nowrap">Input Type:</Label>
                        <Select value={itaeInputType} onValueChange={(v) => setItaeInputType(v as ItaeInputType)}>
                            <SelectTrigger id="itaeInputType" className="flex-1"><SelectValue /></SelectTrigger>
                            <SelectContent>
                                <SelectItem value="disturbance">Disturbance Rejection</SelectItem>
                                <SelectItem value="setpoint">Set-Point Tracking</SelectItem>
                            </SelectContent>
                        </Select>
                    </div>
                </>
            );
        case 'amigo':
             return (
                <>
                    <InputGroup label="Process Gain (K)" tooltip="For FOPTD, standard gain. For Integrator model, K is the integrator gain.">
                        <Input id="k" value={k} onChange={(e) => setK(e.target.value)} placeholder="e.g., 1.0" />
                    </InputGroup>
                    <InputGroup label="Time Constant (τ)" tooltip="Time constant for FOPTD model. Set to 0 for integrator model.">
                        <Input id="tau" value={tau} onChange={(e) => setTau(e.target.value)} placeholder="e.g., 10.0 or 0" />
                    </InputGroup>
                    <InputGroup label="Dead Time (θ)" tooltip="Delay before the process output starts to respond.">
                        <Input id="theta" value={theta} onChange={(e) => setTheta(e.target.value)} placeholder="e.g., 2.0" />
                    </InputGroup>
                </>
             );
        case 'imc':
            return (
                <>
                    <InputGroup label="Process Gain (K)" tooltip="Ratio of change in output to the change in input at steady state.">
                        <Input id="k" value={k} onChange={(e) => setK(e.target.value)} placeholder="e.g., 1.0" />
                    </InputGroup>
                    <InputGroup label="Time Constant (τ)" tooltip="Primary time constant of the process model.">
                        <Input id="tau" value={tau} onChange={(e) => setTau(e.target.value)} placeholder="e.g., 10.0" />
                    </InputGroup>
                    <InputGroup label="Dead Time (θ)" tooltip="Delay before the process output starts to respond.">
                        <Input id="theta" value={theta} onChange={(e) => setTheta(e.target.value)} placeholder="e.g., 2.0" />
                    </InputGroup>
                     <InputGroup label={<span>IMC Tuning Parameter (τ<sub>c</sub>)</span>} tooltip="The desired closed-loop time constant. A larger value gives a slower, more robust response.">
                        <Input id="tauC" value={tauC} onChange={(e) => setTauC(e.target.value)} placeholder="e.g., 2.0" />
                    </InputGroup>
                </>
            );
    }
  };
  
  // FIX 2: Define a specific props interface for the child component.
  interface InputGroupChildProps {
    id: string; // The child must have an id prop.
  }

  // Use the specific interface to type the `children` prop.
  const InputGroup = ({ label, tooltip, children }: { label: string | React.ReactNode, tooltip: string, children: React.ReactElement<InputGroupChildProps> }) => (
    <div className="space-y-2">
      <Label htmlFor={children.props.id} className="flex items-center">
        {label}
        <Tooltip>
          <TooltipTrigger asChild>
            <Button variant="ghost" size="icon" className="h-5 w-5 rounded-full ml-1">
              <Info className="h-3 w-3" />
            </Button>
          </TooltipTrigger>
          <TooltipContent><p>{tooltip}</p></TooltipContent>
        </Tooltip>
      </Label>
      {children}
    </div>
  );

  const formatResult = (value: number | null) => {
    if (value === null || isNaN(value)) return '---';
    if (Math.abs(value) < 1e-4 && value !== 0) return value.toExponential(3);
    return value.toFixed(4);
  }

  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8 px-8 md:px-32">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Column 1: Controls */}
          <div className="lg:col-span-1 space-y-6">
            <Card>
              <CardHeader>
                <CardTitle>PID Tuning Parameters</CardTitle>
              </CardHeader>
              <CardContent className="space-y-6 p-4">
                <div className="space-y-4">
                  {/* Controller Type Selector */}
                  <div className="flex items-center gap-2">
                      <Label htmlFor="controllerType" className="text-sm font-medium whitespace-nowrap">Controller Type:</Label>
                      <Select value={controllerType} onValueChange={(value) => setControllerType(value as ControllerType)}>
                          <SelectTrigger id="controllerType" className="flex-1"><SelectValue /></SelectTrigger>
                          <SelectContent>
                              {getAvailableControllerTypes(tuningMethod).map(type => (
                                <SelectItem key={type} value={type}>{type}</SelectItem>
                              ))}
                          </SelectContent>
                      </Select>
                  </div>

                  {/* Tuning Method Selector */}
                  <div className="flex items-center gap-2">
                      <Label htmlFor="tuningMethod" className="text-sm font-medium whitespace-nowrap">Tuning Method:</Label>
                      <Select value={tuningMethod} onValueChange={(value) => setTuningMethod(value as TuningMethod)}>
                          <SelectTrigger id="tuningMethod" className="flex-1"><SelectValue /></SelectTrigger>
                          <SelectContent>
                              <SelectItem value="itae">ITAE</SelectItem>
                              <SelectItem value="ziegler-nichols">Ziegler-Nichols</SelectItem>
                              <SelectItem value="amigo">AMIGO</SelectItem>
                              <SelectItem value="imc">IMC</SelectItem>
                          </SelectContent>
                      </Select>
                  </div>

                  {/* Model Case Selector (only shown for IMC) */}
                  {tuningMethod === 'imc' && (
                    <div className="flex items-center gap-2">
                        <Label htmlFor="imcModelCase" className="text-sm font-medium whitespace-nowrap">Model Case:</Label>
                        <Select value={imcModelCase} onValueChange={(v) => setImcModelCase(v as ImcModelCase)}>
                            <SelectTrigger id="imcModelCase" className="flex-1"><SelectValue /></SelectTrigger>
                            <SelectContent>
                                <SelectItem value="G">G: FOPTD</SelectItem>
                                <SelectItem value="H">H: Pure Time Delay</SelectItem>
                            </SelectContent>
                        </Select>
                    </div>
                  )}
                  <hr />
                  {/* Dynamically Rendered Inputs */}
                  {renderInputs()}
                </div>
                <Button onClick={handleCalculate} disabled={loading} className="w-full">
                  {loading ? 'Calculating...' : 'Calculate Tuning'}
                </Button>
                {error && <p className="text-sm text-red-500 mt-2">{error}</p>}
              </CardContent>
            </Card>
          </div>
          
          {/* Column 2: Results */}
          <div className="lg:col-span-2 space-y-6">
            <Card>
              <CardHeader>
                <CardTitle>Calculated Controller Settings</CardTitle>
                <p className="text-sm text-muted-foreground pt-1">{displayedMethod || 'Results will appear here'}</p>
              </CardHeader>
              <CardContent>
                {loading && ( <div className="flex items-center justify-center text-muted-foreground p-8">Calculating...</div> )}
                {!loading && !result && ( <div className="flex items-center justify-center text-muted-foreground p-8">Please provide inputs and click calculate.</div> )}
                {result && !loading && (
                    <div className="grid grid-cols-1 sm:grid-cols-3 gap-4 text-center">
                        <div className="p-4 bg-muted rounded-md">
                            <p className="text-sm font-medium" dangerouslySetInnerHTML={{ __html: 'Proportional Gain (K<sub>c</sub>)' }} />
                            <p className="text-2xl font-bold font-mono">{formatResult(result.kc)}</p>
                        </div>
                        <div className="p-4 bg-muted rounded-md">
                            <p className="text-sm font-medium" dangerouslySetInnerHTML={{ __html: 'Integral Time (τ<sub>I</sub>)' }} />
                            <p className="text-2xl font-bold font-mono">{formatResult(result.tauI)}</p>
                        </div>
                        <div className="p-4 bg-muted rounded-md">
                            <p className="text-sm font-medium" dangerouslySetInnerHTML={{ __html: 'Derivative Time (τ<sub>D</sub>)' }} />
                            <p className="text-2xl font-bold font-mono">{formatResult(result.tauD)}</p>
                        </div>
                    </div>
                )}
              </CardContent>
            </Card>

            <Card>
                <CardHeader>
                    <CardTitle>Controller Equation (Standard Form)</CardTitle>
                </CardHeader>
                <CardContent className="text-center font-mono text-lg bg-muted p-6 rounded-md">
                    G<sub>c</sub>(s) = K<sub>c</sub> (1 + <sup>1</sup>⁄<sub>τ<sub>I</sub>s</sub> + τ<sub>D</sub>s)
                </CardContent>
            </Card>

          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}