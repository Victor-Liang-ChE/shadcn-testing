'use client';

import React, { useState, useEffect, useCallback, useMemo } from 'react';
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { TooltipProvider, Tooltip, TooltipTrigger, TooltipContent } from "@/components/ui/tooltip";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { AlertCircle, PlusCircle, RefreshCw, Trash2 } from 'lucide-react';

// Helper types
type ReactorType = 'PFR' | 'CSTR';
type ReactionPhase = 'Liquid' | 'Gas';
type RateConstantInputMode = 'directK' | 'arrhenius';

interface Component {
  id: string;
  name: string;
  initialFlowRate: string; // F_i0 (e.g., mol/s)
  // For stoichiometry parsing results:
  isReactant?: boolean;
  isProduct?: boolean;
  stoichiometricCoefficient?: number;
  reactionOrder?: string; // For reactants
}

interface Kinetics {
  rateConstantInputMode: RateConstantInputMode;
  kValue: string; // Rate constant k
  AValue: string; // Pre-exponential factor A
  EaValue: string; // Activation energy Ea (e.g., J/mol)
  reactionTempK: string; // Reaction temperature in Kelvin
}

interface ParsedComponent extends Component {
  stoichiometricCoefficient: number; // Now mandatory after parsing
  role: 'reactant' | 'product' | 'inert'; // Inert if present in feed but not in reaction
  reactionOrderNum?: number; // Parsed reaction order
}

interface ParsedReaction {
  reactants: ParsedComponent[];
  products: ParsedComponent[];
  allInvolved: ParsedComponent[]; // All components from stoichiometry
  rateConstantAtT?: number; // Calculated k at reaction temperature
}

interface CalculationResults {
  conversion?: { reactantName: string; value: number }; // X for a key reactant
  outletFlowRates?: { [componentName: string]: number };
  error?: string;
  // Add selectivity, concentrations etc. later
}

const R_GAS_CONSTANT = 8.314; // J/(mol*K) or L*kPa/(mol*K) depending on units

export default function ReactorDesignPage() {
  const [reactorType, setReactorType] = useState<ReactorType>('PFR');
  const [reactionPhase, setReactionPhase] = useState<ReactionPhase>('Liquid');
  const [reactionString, setReactionString] = useState<string>('A + B -> C');
  const [reactions, setReactions] = useState<Array<{id: string, reactants: string, products: string}>>([
    { id: '1', reactants: 'A + B', products: 'C' }
  ]); // Example: 2 H2 + O2 -> 2 H2O
  
  const [components, setComponents] = useState<Component[]>([
    { id: 'comp1', name: 'A', initialFlowRate: '10', reactionOrder: '1' },
    { id: 'comp2', name: 'B', initialFlowRate: '10', reactionOrder: '1' },
    { id: 'comp3', name: 'C', initialFlowRate: '0' },
  ]);

  const [kinetics, setKinetics] = useState<Kinetics>({
    rateConstantInputMode: 'directK',
    kValue: '0.1', // Units depend on reaction order and concentration units
    AValue: '1e6',
    EaValue: '50000', // e.g., J/mol
    reactionTempK: '300', // Kelvin
  });

  const [reactorVolume, setReactorVolume] = useState<string>('100'); // e.g., Liters
  const [totalPressure, setTotalPressure] = useState<string>('101.325'); // e.g., kPa (for gas phase)
  const [volumetricFlowRate, setVolumetricFlowRate] = useState<string>('1'); // v0 (for liquid phase, L/s)


  const [parsedReaction, setParsedReaction] = useState<ParsedReaction | null>(null);
  const [calculationResults, setCalculationResults] = useState<CalculationResults | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [calculationError, setCalculationError] = useState<string | null>(null);

  // Handle Enter key press to trigger calculation
  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleCalculate();
    }
  };

  // Component Management
  const addReaction = () => {
    const newId = (reactions.length + 1).toString();
    setReactions([...reactions, { id: newId, reactants: '', products: '' }]);
  };

  const removeReaction = (id: string) => {
    if (reactions.length > 1) {
      setReactions(reactions.filter(reaction => reaction.id !== id));
    }
  };

  const updateReaction = (id: string, field: 'reactants' | 'products', value: string) => {
    setReactions(reactions.map(reaction => 
      reaction.id === id ? { ...reaction, [field]: value } : reaction
    ));
    // Update the main reaction string based on the first reaction
    if (id === reactions[0]?.id || reactions.length === 1) {
      const updatedReaction = reactions.find(r => r.id === id);
      if (updatedReaction) {
        const newReactants = field === 'reactants' ? value : updatedReaction.reactants;
        const newProducts = field === 'products' ? value : updatedReaction.products;
        setReactionString(`${newReactants} -> ${newProducts}`);
      }
    }
  };

  const addComponent = () => {
    setComponents([...components, { id: `comp${Date.now()}`, name: '', initialFlowRate: '', reactionOrder: '' }]);
  };

  const removeComponent = (id: string) => {
    setComponents(components.filter(c => c.id !== id));
  };

  const handleComponentChange = (id: string, field: keyof Component, value: string) => {
    setComponents(components.map(c => c.id === id ? { ...c, [field]: value } : c));
  };

  const handleKineticsChange = (field: keyof Kinetics, value: string) => {
    setKinetics(prev => ({ ...prev, [field]: value }));
  };
  
  // Stoichiometry Parsing Logic
  const parseStoichiometry = useCallback(() => {
    setCalculationError(null);
    setParsedReaction(null);
    if (!reactionString.trim()) {
      setCalculationError("Reaction stoichiometry cannot be empty.");
      return;
    }

    const arrow = reactionString.includes('->') ? '->' : reactionString.includes('<=>') ? '<=>' : null;
    if (!arrow) {
      setCalculationError("Invalid reaction string: missing '->' or '<=>'.");
      return;
    }

    const [reactantsStr, productsStr] = reactionString.split(arrow).map(s => s.trim());
    if (!reactantsStr || !productsStr) {
      setCalculationError("Invalid reaction string: reactants or products side is empty.");
      return;
    }
    
    const parseSide = (sideStr: string, role: 'reactant' | 'product'): ParsedComponent[] => {
      return sideStr.split('+').map(termStr => {
        const term = termStr.trim();
        const match = term.match(/^(\d*)\s*([A-Za-z0-9]+)$/); // Matches "2 H2O" or "N2"
        if (!match) {
          throw new Error(`Invalid term in reaction: ${term}`);
        }
        const coeff = match[1] ? parseInt(match[1]) : 1;
        const name = match[2];

        const feedComponent = components.find(c => c.name.trim().toLowerCase() === name.trim().toLowerCase());
        if (!feedComponent) {
          throw new Error(`Component '${name}' in reaction not found in feed components list. Please add it.`);
        }
        
        const reactionOrderNum = role === 'reactant' && feedComponent.reactionOrder ? parseFloat(feedComponent.reactionOrder) : undefined;
        if (role === 'reactant' && feedComponent.reactionOrder && reactionOrderNum !== undefined && (isNaN(reactionOrderNum) || reactionOrderNum < 0)) {
            throw new Error(`Invalid reaction order for ${name}. Must be a non-negative number.`);
        }

        return {
          ...feedComponent,
          stoichiometricCoefficient: coeff,
          role,
          reactionOrderNum: reactionOrderNum
        };
      });
    };

    try {
      const parsedReactants = parseSide(reactantsStr, 'reactant');
      const parsedProducts = parseSide(productsStr, 'product');
      
      const allInvolvedNames = new Set([...parsedReactants.map(r => r.name), ...parsedProducts.map(p => p.name)]);
      const allInvolvedFromFeed = components
        .filter(fc => allInvolvedNames.has(fc.name))
        .map(fc => {
            const reactantInfo = parsedReactants.find(r => r.name === fc.name);
            const productInfo = parsedProducts.find(p => p.name === fc.name);
            if (reactantInfo) return reactantInfo;
            if (productInfo) return productInfo;
            return { // Should not happen if logic is correct
                ...fc, 
                stoichiometricCoefficient: 0, 
                role: 'inert' as 'inert', // Should be caught by earlier check
            };
        });


      // Calculate k at T using Arrhenius equation
      let kAtT: number | undefined;
      const T_K = parseFloat(kinetics.reactionTempK);
      if (isNaN(T_K) || T_K <= 0) throw new Error("Reaction temperature must be a positive number.");

      const A = parseFloat(kinetics.AValue);
      const Ea = parseFloat(kinetics.EaValue);
      if (isNaN(A) || isNaN(Ea)) throw new Error("Invalid A or Ea value for Arrhenius equation.");
      kAtT = A * Math.exp(-Ea / (R_GAS_CONSTANT * T_K));
      if (kAtT <=0) throw new Error("Rate constant (k) must be positive.");


      setParsedReaction({ reactants: parsedReactants, products: parsedProducts, allInvolved: allInvolvedFromFeed, rateConstantAtT: kAtT });
    } catch (e: any) {
      setCalculationError(`Stoichiometry/Kinetics Error: ${e.message}`);
      setParsedReaction(null);
    }
  }, [reactionString, components, kinetics]);

  // Effect to auto-parse stoichiometry when inputs change
  useEffect(() => {
    parseStoichiometry();
  }, [parseStoichiometry]);


  // Main Calculation Logic
  const handleCalculate = useCallback(async () => {
    setIsLoading(true);
    setCalculationResults(null);
    setCalculationError(null); // Clear previous specific calculation errors

    if (!parsedReaction || !parsedReaction.rateConstantAtT) {
      if (!calculationError) setCalculationError("Reaction stoichiometry or kinetics not parsed correctly. Please check inputs.");
      setIsLoading(false);
      return;
    }
    
    const { reactants, products, rateConstantAtT } = parsedReaction;
    const V = parseFloat(reactorVolume);
    const T_K = parseFloat(kinetics.reactionTempK); // Already validated in parseStoichiometry

    if (isNaN(V) || V <= 0) {
      setCalculationError("Reactor volume must be a positive number.");
      setIsLoading(false);
      return;
    }
    if (reactants.length === 0) {
        setCalculationError("No reactants identified in the reaction.");
        setIsLoading(false);
        return;
    }

    // For simplicity, assume the first reactant is the limiting one or key reactant
    // A more robust solution would identify the limiting reactant based on stoichiometry and feed.
    const keyReactant = reactants[0];
    const F_key0 = parseFloat(keyReactant.initialFlowRate);
    if (isNaN(F_key0) || F_key0 <= 0) {
        setCalculationError(`Initial flow rate for key reactant ${keyReactant.name} must be positive.`);
        setIsLoading(false);
        return;
    }

    // Initial concentrations and total volumetric flow rate
    let C_key0: number; // Concentration of key reactant at inlet
    let v0_calc: number; // Total volumetric flow rate at inlet (L/s or m3/s)

    const P_total_kPa = parseFloat(totalPressure); // For gas phase
    const v0_input_Ls = parseFloat(volumetricFlowRate); // For liquid phase if provided

    // Calculate initial concentrations based on phase
    const initialConcentrations: { [name: string]: number } = {};
    const F_total0 = components.reduce((sum, c) => sum + (parseFloat(c.initialFlowRate) || 0), 0);

    if (reactionPhase === 'Liquid') {
        if (!isNaN(v0_input_Ls) && v0_input_Ls > 0) {
            v0_calc = v0_input_Ls;
        } else {
            // Try to estimate v0 if not given, e.g. if one component's C0 is known, or sum of partial molar volumes.
            // This is a simplification. Often v0 is a direct input for liquid systems.
            setCalculationError("For liquid phase, please provide total volumetric flow rate (v0).");
            setIsLoading(false); return;
        }
        components.forEach(c => {
            const Fi0 = parseFloat(c.initialFlowRate) || 0;
            initialConcentrations[c.name] = Fi0 / v0_calc;
        });
        C_key0 = initialConcentrations[keyReactant.name];

    } else { // Gas Phase
        if (isNaN(P_total_kPa) || P_total_kPa <= 0) {
            setCalculationError("For gas phase, total pressure must be a positive number.");
            setIsLoading(false); return;
        }
        // Ideal Gas Law: C_i = P_i / (R*T) = (y_i * P_total) / (R*T)
        // y_i = F_i0 / F_total0
        // R in L*kPa/(mol*K) is 8.314
        components.forEach(c => {
            const Fi0 = parseFloat(c.initialFlowRate) || 0;
            const yi0 = F_total0 > 0 ? Fi0 / F_total0 : 0;
            initialConcentrations[c.name] = (yi0 * P_total_kPa) / (R_GAS_CONSTANT * T_K);
        });
        C_key0 = initialConcentrations[keyReactant.name];
        // v0 for gas phase: F_total0 * R * T_K / P_total_kPa
        v0_calc = (F_total0 * R_GAS_CONSTANT * T_K) / P_total_kPa; 
    }
    
    if (C_key0 <= 0 && F_key0 > 0) { // Check if key reactant concentration is valid
        setCalculationError(`Initial concentration of key reactant ${keyReactant.name} is zero or negative. Check feed conditions.`);
        setIsLoading(false); return;
    }


    // Placeholder for rate law function
    // Example: -rA = k * C_A^orderA * C_B^orderB
    // Concentrations here are current concentrations, not initial.
    const rateLaw = (currentConcentrations: { [name: string]: number }): number => {
        let rate = rateConstantAtT!;
        for (const r of reactants) {
            const conc = currentConcentrations[r.name] || 0;
            const order = r.reactionOrderNum || 0; // Default to 0 if not specified, though validated earlier
            if (conc < 0) return -1; // Invalid concentration
            rate *= Math.pow(conc, order);
        }
        return rate; // This is the rate of disappearance of the key reactant (e.g. -r_A)
    };


    try {
        let conversion = 0;
        const outletFlowRates: { [componentName: string]: number } = {};

        if (reactorType === 'CSTR') {
            // CSTR: V = F_A0 * X / (-r_A_exit)
            // Need to solve for X. This often requires a numerical solver for non-linear -r_A.
            // For simplicity, let's assume a first-order reaction A -> Products, -rA = k * CA
            // CA = CA0 * (1-X)
            // V = F_A0 * X / (k * CA0 * (1-X))
            // Let tau = V / v0_calc. Da = k * tau (Damkohler number for 1st order)
            // X = Da / (1 + Da)
            // This is a MAJOR simplification. A general CSTR solver is complex.

            // Iterative solver for X (e.g., bisection or Newton-Raphson)
            let X_guess = 0.5; // Initial guess for conversion
            const maxIter = 100;
            let iter = 0;
            const tol = 1e-6;

            while (iter < maxIter) {
                const currentConcentrations_exit: { [name: string]: number } = {};
                // Calculate exit concentrations based on X_guess for all reactants
                reactants.forEach(r => {
                    const Fi0_r = parseFloat(r.initialFlowRate);
                    // Stoichiometric ratio relative to keyReactant
                    const nu_r_over_nu_key = (r.stoichiometricCoefficient / keyReactant.stoichiometricCoefficient);
                    const Ci0_r = initialConcentrations[r.name];
                    
                    if (reactionPhase === 'Liquid') {
                         // Theta_i = Ci0_r / C_key0
                        const theta_i = Ci0_r / C_key0;
                        currentConcentrations_exit[r.name] = C_key0 * (theta_i - nu_r_over_nu_key * X_guess);
                    } else { // Gas phase - more complex due to epsilon
                        // For simplicity, approximate with liquid phase equations or use molar flows
                        // F_i = F_i0 - (nu_i / nu_key) * F_key0 * X
                        // C_i = F_i / v_exit. v_exit = v0 * (1 + epsilon * X) * (P0/P) * (T/T0)
                        // This gets very involved. Let's use a simplified liquid-phase like concentration calc for now.
                        const Fi_exit = Fi0_r - nu_r_over_nu_key * F_key0 * X_guess;
                        // This needs a proper gas phase concentration model at exit
                        currentConcentrations_exit[r.name] = Math.max(0, initialConcentrations[r.name] * (1 - (nu_r_over_nu_key / (Fi0_r/F_key0)) * X_guess) ); // Very rough approx.
                        if (r.name === keyReactant.name) currentConcentrations_exit[r.name] = C_key0 * (1 - X_guess);
                    }
                     if (currentConcentrations_exit[r.name] < 0) currentConcentrations_exit[r.name] = 0; // Cannot be negative
                });

                const r_A_exit = rateLaw(currentConcentrations_exit);
                if (r_A_exit <= 0) { // Rate cannot be zero or negative if reaction is proceeding
                     // If rate is zero, it means X might be too high (reactant consumed) or too low (no reaction)
                    if (X_guess > 1e-3 && Object.values(currentConcentrations_exit).some(c => c < 1e-9)) { // Reactant nearly consumed
                        // X_guess might be the solution or slightly overestimate
                        break; 
                    }
                    // If X_guess is small and rate is zero, something is wrong or k is zero
                    // setCalculationError("CSTR: Calculated reaction rate at exit is zero or negative. Check kinetics or orders.");
                    // For now, break and use current X_guess or set to 0 if rate is always 0.
                    if (X_guess < tol && r_A_exit <=0) X_guess = 0; // No reaction
                    break; 
                }

                const X_calc = (r_A_exit * V) / (F_key0); // From V * r_A = F_A0 * X
                
                if (Math.abs(X_calc - X_guess) < tol) {
                    X_guess = X_calc;
                    break;
                }
                X_guess = (X_calc + X_guess) / 2; // Simple averaging, not robust.
                // A better root finding method (e.g. f(X) = F_A0*X - (-r_A(X))*V = 0) is needed.
                iter++;
            }
            if (iter === maxIter) {
                 console.warn("CSTR solver did not converge. Result might be inaccurate.");
            }
            conversion = Math.max(0, Math.min(X_guess, 1.0)); // Ensure conversion is between 0 and 1

        } else { // PFR
            // PFR: dX/dV = -r_A / F_A0
            // Requires numerical integration (e.g., Euler, Runge-Kutta)
            // This is a placeholder for a proper ODE solver.
            // Simple Euler method for demonstration (highly inaccurate for many cases):
            const dV = V / 100; // Number of steps
            let X_current = 0;
            let currentVolume = 0;

            for (let i = 0; i < 100; i++) {
                const currentConcentrations_step: { [name: string]: number } = {};
                 reactants.forEach(r => {
                    const Ci0_r = initialConcentrations[r.name];
                    const nu_r_over_nu_key = (r.stoichiometricCoefficient / keyReactant.stoichiometricCoefficient);
                    if (reactionPhase === 'Liquid') {
                        const theta_i = Ci0_r / C_key0;
                        currentConcentrations_step[r.name] = C_key0 * (theta_i - nu_r_over_nu_key * X_current);
                    } else { // Gas phase - simplified
                         currentConcentrations_step[r.name] = Math.max(0, initialConcentrations[r.name] * (1 - (nu_r_over_nu_key / ( (parseFloat(r.initialFlowRate)/F_key0) ) ) * X_current) );
                         if (r.name === keyReactant.name) currentConcentrations_step[r.name] = C_key0 * (1-X_current);
                    }
                    if (currentConcentrations_step[r.name] < 0) currentConcentrations_step[r.name] = 0;
                });

                const r_A_step = rateLaw(currentConcentrations_step);
                if (r_A_step <= 0) { // Reaction stopped or invalid rate
                    if (X_current > 1e-3 && Object.values(currentConcentrations_step).some(c => c < 1e-9)) {
                         break; // Reactant consumed
                    }
                    // If rate is zero early, something is wrong.
                    break;
                }
                
                const dX = (r_A_step / F_key0) * dV;
                X_current += dX;
                currentVolume += dV;
                if (X_current >= 1.0) { X_current = 1.0; break; }
                if (X_current < 0) { X_current = 0; break; } // Should not happen with positive rate
            }
            conversion = Math.max(0, Math.min(X_current, 1.0));
        }

        // Calculate outlet flow rates for all components (reactants, products, inerts)
        components.forEach(c => {
            const F_i0 = parseFloat(c.initialFlowRate) || 0;
            const parsedCompInfo = parsedReaction.allInvolved.find(pc => pc.name === c.name);

            if (parsedCompInfo) { // Component is involved in the reaction
                const nu_i = parsedCompInfo.role === 'reactant' ? -parsedCompInfo.stoichiometricCoefficient : parsedCompInfo.stoichiometricCoefficient;
                // F_i = F_i0 + (nu_i / nu_key_abs) * F_key0 * X
                // nu_key_abs is absolute stoich coeff of key reactant
                const nu_key_abs = keyReactant.stoichiometricCoefficient; 
                outletFlowRates[c.name] = F_i0 + (nu_i / nu_key_abs) * F_key0 * conversion;
            } else { // Inert component
                outletFlowRates[c.name] = F_i0;
            }
            if (outletFlowRates[c.name] < 0) outletFlowRates[c.name] = 0; // Flow rate cannot be negative
        });
        
        setCalculationResults({
            conversion: { reactantName: keyReactant.name, value: conversion },
            outletFlowRates,
        });

    } catch (e: any) {
        setCalculationError(`Calculation Error: ${e.message}`);
    } finally {
        setIsLoading(false);
    }
  }, [parsedReaction, reactorType, reactorVolume, kinetics, components, reactionPhase, totalPressure, volumetricFlowRate, calculationError]);

  // Remove automatic calculation - calculations will only run when button is pressed
  // useEffect(() => {
  //   calculateReactorPerformance();
  // }, [reactorType, reactionPhase, reactorVolume, volumetricFlowRate, totalPressure, kinetics, components, reactionString]);


  // Reactor SVG Visualization
  const ReactorVisualization = useMemo(() => {
    const commonArrowProps = { stroke: "currentColor", strokeWidth: "2", fill: "currentColor" };
    const reactorFill = "fill-gray-300 dark:fill-gray-700";
    const reactorStroke = "stroke-gray-500 dark:stroke-gray-400";

    if (reactorType === 'PFR') {
      return (
        <svg viewBox="0 0 220 100" className="w-full h-auto max-h-[200px] text-gray-800 dark:text-gray-200">
          <title>Plug Flow Reactor (PFR)</title>
          {/* Main Body */}
          <rect x="20" y="30" width="180" height="40" rx="15" ry="15" className={`${reactorFill} ${reactorStroke}`} strokeWidth="1"/>
          {/* Inlet */}
          {/* Feed inlet */}
          <path d="M25 50 L45 50" {...commonArrowProps} />
          <polygon points="40,45 45,50 40,55" {...commonArrowProps} />
          <text x="10" y="40" fontSize="10" textAnchor="middle" fill="currentColor">Feed</text>
          {/* Outlet */}
          {/* Shaft of the arrow */}
          <path d="M200 50 L215 50" {...commonArrowProps} />
          {/* Arrowhead, tip at x=220, base at x=215. No transform needed. */}
          <polygon points="215,45 220,50 215,55" {...commonArrowProps} />
          {/* Text label for Product (position can remain as is or be adjusted as preferred) */}
          <text x="207" y="40" fontSize="10" textAnchor="middle" fill="currentColor">Product</text>
        </svg>
      );
    } else { // CSTR
      return (
        <svg viewBox="0 0 100 170" className="w-full h-auto max-h-[250px] text-gray-800 dark:text-gray-200">
          <title>Continuous Stirred-Tank Reactor (CSTR)</title>
          {/* Tank Body */}
          <rect x="15" y="20" width="70" height="120" rx="5" ry="5" className={`${reactorFill} ${reactorStroke}`} strokeWidth="1"/>
          {/* Stirrer Shaft */}
          <line x1="50" y1="5" x2="50" y2="85" stroke="gray" strokeWidth="3" />
          {/* Stirrer Blades (Animated) */}
          <g style={{ transformOrigin: '50px 85px' }} className="animate-spin-slow"> {/* CSS animation defined in global styles or here */}
            <line x1="30" y1="85" x2="70" y2="85" stroke="gray" strokeWidth="5" />
            <line x1="50" y1="70" x2="50" y2="100" stroke="gray" strokeWidth="5" />
          </g>
          {/* Inlet */}
          <path d="M0 40 L15 40" {...commonArrowProps} />
          <polygon points="10,35 15,40 10,45" {...commonArrowProps} />
          <text x="2" y="30" fontSize="10" fill="currentColor">Feed</text>
          {/* Outlet */}
          <path d="M50 140 L50 160" {...commonArrowProps} />
          <polygon points="45,155 50,160 55,155" {...commonArrowProps} />
          <text x="53" y="155" fontSize="10" fill="currentColor">Product</text>
        </svg>
      );
    }
  }, [reactorType]);


  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 space-y-6">
        <div className="grid grid-cols-2 gap-6">
          {/* Left Panel - Input Controls (1/2 of width) */}
          <div className="space-y-6">
            {/* System Configuration & Kinetics */}
            <Card>
              <CardContent className="space-y-4 relative">
                {/* Continuous vertical divider */}
                <div className="absolute left-1/2 top-0 h-full w-px bg-border transform -translate-x-1/2"></div>
                
                <div className="grid grid-cols-2 gap-4">
                  <div className="flex items-center gap-2 pr-4">
                    <Label htmlFor="reactorType" className="text-sm whitespace-nowrap">Reactor Type:</Label>
                    <div className="flex-1 min-w-0">
                      <Select value={reactorType} onValueChange={(value) => setReactorType(value as ReactorType)}>
                        <SelectTrigger id="reactorType" className="w-full">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="PFR">PFR</SelectItem>
                          <SelectItem value="CSTR">CSTR</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                  </div>
                  <div className="flex items-center gap-2 pl-4">
                    <Label htmlFor="reactionPhase" className="text-sm whitespace-nowrap">Phase:</Label>
                    <div className="flex-1 min-w-0">
                      <Select value={reactionPhase} onValueChange={(value) => setReactionPhase(value as ReactionPhase)}>
                        <SelectTrigger id="reactionPhase" className="w-full">
                          <SelectValue />
                        </SelectTrigger>
                        <SelectContent>
                          <SelectItem value="Liquid">Liquid</SelectItem>
                          <SelectItem value="Gas">Gas</SelectItem>
                        </SelectContent>
                      </Select>
                    </div>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="flex items-center gap-2 pr-4">
                    <Label htmlFor="reactorVolume" className="text-sm whitespace-nowrap">Reactor Volume:</Label>
                    <div className="flex items-center gap-1 flex-1 min-w-0">
                      <Input
                        id="reactorVolume"
                        type="number"
                        value={reactorVolume}
                        onChange={(e) => setReactorVolume(e.target.value)}
                        onKeyDown={handleKeyDown}
                        placeholder="100"
                      />
                      <span className="text-xs text-muted-foreground">L</span>
                    </div>
                  </div>
                  <div className="flex items-center gap-2 pl-4">
                    {reactionPhase === 'Gas' && (
                      <>
                        <Label htmlFor="totalPressure" className="text-sm whitespace-nowrap">Total Pressure:</Label>
                        <div className="flex items-center gap-1 flex-1 min-w-0">
                          <Input
                            id="totalPressure"
                            type="number"
                            value={totalPressure}
                            onChange={(e) => setTotalPressure(e.target.value)}
                            onKeyDown={handleKeyDown}
                            placeholder="101.325"
                          />
                          <span className="text-xs text-muted-foreground">kPa</span>
                        </div>
                      </>
                    )}
                    {reactionPhase === 'Liquid' && (
                      <>
                        <Label htmlFor="volumetricFlowRate" className="text-sm whitespace-nowrap">Vol Flow Rate:</Label>
                        <div className="flex items-center gap-1 flex-1 min-w-0">
                          <Input
                            id="volumetricFlowRate"
                            type="number"
                            value={volumetricFlowRate}
                            onChange={(e) => setVolumetricFlowRate(e.target.value)}
                            onKeyDown={handleKeyDown}
                            placeholder="1.0"
                          />
                          <span className="text-xs text-muted-foreground">L/s</span>
                        </div>
                      </>
                    )}
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="flex items-center gap-2 pr-4">
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <Label htmlFor="preExponentialFactor" className="text-sm whitespace-nowrap cursor-help">A:</Label>
                      </TooltipTrigger>
                      <TooltipContent>
                        <p>Pre-exponential factor</p>
                      </TooltipContent>
                    </Tooltip>
                    <Input
                      id="preExponentialFactor"
                      type="number"
                      value={kinetics.AValue}
                      onChange={(e) => handleKineticsChange('AValue', e.target.value)}
                      onKeyDown={handleKeyDown}
                      placeholder="1e6"
                      className="flex-1 min-w-0"
                    />
                  </div>
                  <div className="flex items-center gap-2 pl-4">
                    <Tooltip>
                      <TooltipTrigger asChild>
                        <Label htmlFor="activationEnergy" className="text-sm whitespace-nowrap cursor-help">Ea:</Label>
                      </TooltipTrigger>
                      <TooltipContent>
                        <p>Activation energy</p>
                      </TooltipContent>
                    </Tooltip>
                    <div className="flex items-center gap-1 flex-1 min-w-0">
                      <Input
                        id="activationEnergy"
                        type="number"
                        value={kinetics.EaValue}
                        onChange={(e) => handleKineticsChange('EaValue', e.target.value)}
                        onKeyDown={handleKeyDown}
                        placeholder="50000"
                      />
                      <span className="text-xs text-muted-foreground">J/mol</span>
                    </div>
                  </div>
                </div>

                <div className="grid grid-cols-2 gap-4">
                  <div className="flex items-center gap-2 pr-4">
                    <Label htmlFor="reactionTemp" className="text-sm whitespace-nowrap">Temperature:</Label>
                    <div className="flex items-center gap-1 flex-1 min-w-0">
                      <Input
                        id="reactionTemp"
                        type="number"
                        value={kinetics.reactionTempK}
                        onChange={(e) => handleKineticsChange('reactionTempK', e.target.value)}
                        onKeyDown={handleKeyDown}
                        placeholder="300"
                      />
                      <span className="text-xs text-muted-foreground">K</span>
                    </div>
                  </div>
                  <div className="pl-4"></div>
                </div>
              </CardContent>
            </Card>

            {/* Reaction Stoichiometry */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center justify-between">
                  Reaction Stoichiometry
                  <Button variant="outline" size="sm" onClick={addReaction}>
                    <PlusCircle className="mr-2 h-4 w-4" />Add Reaction
                  </Button>
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                {reactions.map((reaction, index) => (
                  <div key={reaction.id}>
                    <div className="grid grid-cols-12 gap-2 items-center">
                      <div className="col-span-4">
                        <Input
                          placeholder="A + B"
                          value={reaction.reactants}
                          onChange={(e) => updateReaction(reaction.id, 'reactants', e.target.value)}
                        />
                      </div>
                      <div className="col-span-2 text-center">
                        <span className="text-lg">→</span>
                      </div>
                      <div className="col-span-4">
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
                    {index < reactions.length - 1 && <hr className="my-4" />}
                  </div>
                ))}
              </CardContent>
            </Card>

            {/* Component Feed Composition */}
            <Card>
              <CardContent className="space-y-4">
                <div className="flex justify-between items-center">
                  <span className="text-sm font-medium">Components</span>
                  <Button variant="outline" size="sm" onClick={addComponent}><PlusCircle className="mr-2 h-4 w-4" />Add Component</Button>
                </div>
                <div className="grid grid-cols-12 gap-2 items-center text-xs font-medium text-muted-foreground border-b pb-1">
                  <div className="col-span-3">Name</div>
                  <div className="col-span-4">Initial Flow Rate</div>
                  <div className="col-span-3">Reaction Order</div>
                  <div className="col-span-2"></div>
                </div>
                {components.map((comp, index) => (
                  <div key={comp.id} className="grid grid-cols-12 gap-2 items-center border-b pb-2">
                    <div className="col-span-3">
                      <Input
                        placeholder="A"
                        value={comp.name}
                        onChange={(e) => handleComponentChange(comp.id, 'name', e.target.value)}
                      />
                    </div>
                    <div className="col-span-4">
                      <div className="flex items-center gap-1">
                        <Input
                          type="number"
                          placeholder="10.0"
                          value={comp.initialFlowRate}
                          onChange={(e) => handleComponentChange(comp.id, 'initialFlowRate', e.target.value)}
                          onKeyDown={handleKeyDown}
                        />
                        <span className="text-xs text-muted-foreground whitespace-nowrap">mol/s</span>
                      </div>
                    </div>
                    <div className="col-span-3">
                      <div className="flex items-center gap-1">
                        <Input
                          type="number"
                          placeholder="1.0"
                          value={comp.reactionOrder}
                          onChange={(e) => handleComponentChange(comp.id, 'reactionOrder', e.target.value)}
                          onKeyDown={handleKeyDown}
                        />
                      </div>
                    </div>
                    <div className="col-span-2 flex justify-end">
                      <Button
                        variant="ghost"
                        size="icon"
                        onClick={() => removeComponent(comp.id)}
                        disabled={components.length <= 1}
                        className="text-red-500 hover:text-red-700 hover:bg-red-50"
                      >
                        <Trash2 className="h-4 w-4" />
                      </Button>
                    </div>
                  </div>
                ))}
              </CardContent>
            </Card>
          </div>

          {/* Right Panel - Visualization and Results (1/2 of width) */}
          <div className="space-y-6">
            {/* Reactor Visualization */}
            <Card>
              <CardContent>
                <div className="flex justify-center items-center h-64 bg-gray-50 dark:bg-gray-900 rounded-lg">
                  {ReactorVisualization}
                </div>
              </CardContent>
            </Card>

            {/* Calculate Button - Always below visualization */}
            <div className="w-full">
              <Button
                onClick={handleCalculate}
                disabled={isLoading || !parsedReaction}
                className="w-full"
                size="lg"
              >
                {isLoading ? 'Calculating...' : 'Calculate Reactor Performance'}
              </Button>
            </div>

            {/* Error Display */}
            {calculationError && (
              <Alert variant="destructive">
                <AlertCircle className="h-4 w-4" />
                <AlertTitle>Error</AlertTitle>
                <AlertDescription>{calculationError}</AlertDescription>
              </Alert>
            )}

            {/* Results Display */}
            {calculationResults && (
              <Card>
                <CardContent>
                  <div className="space-y-4">
                    {calculationResults.conversion && (
                      <div className="flex items-center gap-2">
                      <span className="font-semibold">Conversion (X<sub>{calculationResults.conversion.reactantName}</sub>):</span>
                      <span>{(calculationResults.conversion.value * 100).toFixed(2)}%</span>
                    </div>
                    )}

                    {calculationResults.outletFlowRates && (
                      <div>
                        <h4 className="font-semibold">Outlet Flow Rates:</h4>
                        <div className="space-y-1">
                          {Object.entries(calculationResults.outletFlowRates).map(([comp, flowRate]) => (
                            <div key={comp} className="flex items-center gap-1">
                              <span>{comp}:</span>
                              <span>{flowRate.toFixed(4)} mol/s</span>
                            </div>
                          ))}
                        </div>
                      </div>
                    )}
                  </div>
                </CardContent>
              </Card>
            )}
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}
