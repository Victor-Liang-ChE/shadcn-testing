"use client"

import React, { useState, useMemo, useEffect, useCallback } from 'react'
import ReactECharts from 'echarts-for-react'
import { EChartsOption } from 'echarts'
import { PlusCircle, Trash2, ArrowLeft } from 'lucide-react'
import { Label } from '@/components/ui/label'
import { Input } from '@/components/ui/input'
import { Slider } from '@/components/ui/slider'
import { useTheme } from 'next-themes';
import { Button } from '@/components/ui/button';
import { CardTitle } from '@/components/ui/card';
import { Card, CardHeader, CardContent } from '@/components/ui/card';
import { solvePFR_ODE_System } from '@/lib/reactor-solver'
import { parseParallelReactions } from '@/lib/reaction-parser';
import { Select } from '@/components/ui/select';
import { SelectItem } from '@/components/ui/select';
import { SelectTrigger } from '@/components/ui/select';
import { SelectValue } from '@/components/ui/select';
import { SelectContent } from '@/components/ui/select';

// --- Type Definitions ---
interface ReactionSetup {
  id: string
  AValue: string
  EaValue: string
  isEquilibrium: boolean
  AValueBackward?: string
  EaValueBackward?: string
}

interface ComponentSetup {
  id: string
  name: string
  reactionData: {
    [reactionId: string]: {
      stoichiometry: string
      order: string
      orderReverse?: string // Added for reversible reactions
    }
  }
  initialFlowRate?: string
}

type ReactorType = 'PFR' | 'CSTR'
type GraphType = 'selectivity' | 'volume' | 'conversion' | 'flowrates' | 'composition' | 'selectivityVsConversion';

// Default Data (now used as initial state for the editable form)
const initialReactions: ReactionSetup[] = [
  { id: '1', AValue: '3.66e14', EaValue: '101600', isEquilibrium: false },
  { id: '2', AValue: '4.77e16', EaValue: '110850', isEquilibrium: false },
]

const initialComponents: ComponentSetup[] = [
  { id: 'comp-A', name: '1-Butene', reactionData: { '1': { stoichiometry: '-1', order: '1' }, '2': { stoichiometry: '-1', order: '1' } } },
  { id: 'comp-B', name: 'Isobutane', reactionData: { '1': { stoichiometry: '-1', order: '1' }, '2': { stoichiometry: '0', order: '0' } } },
  { id: 'comp-C', name: 'Isooctane', reactionData: { '1': { stoichiometry: '1', order: '0' }, '2': { stoichiometry: '-1', order: '1' } } },
  { id: 'comp-D', name: 'Dodecane', reactionData: { '1': { stoichiometry: '0', order: '0' }, '2': { stoichiometry: '1', order: '0' } } },
]


// --- Calculation Engine (Now more generic) ---
const R = 8.314

// Helper function to perform linear interpolation
const interpolate = (
  targetX: number,
  p1: { x: number; y: number },
  p2: { x: number; y: number }
): number => {
  if (p1.x === p2.x) {
    return p1.y; // Avoid division by zero
  }
  // Standard linear interpolation formula
  return p1.y + ((targetX - p1.x) * (p2.y - p1.y)) / (p2.x - p1.x);
};

const calculateKinetics = (reactions: ReactionSetup[], tempK: number) => {
  // This assumes the first two reactions are the key ones for the model.
  const k1 = parseFloat(reactions[0]?.AValue || '0') * Math.exp(-parseFloat(reactions[0]?.EaValue || '0') / (R * tempK))
  const k2 = parseFloat(reactions[1]?.AValue || '0') * Math.exp(-parseFloat(reactions[1]?.EaValue || '0') / (R * tempK))
  return { k1, k2 }
}

// This is the new, debug-friendly solver function
const solveForInitialFlows = (
  targetProductionRate: number,
  targetComponentName: string,
  reactorVolume: number,
  reactorType: ReactorType,
  molarRatios: { numeratorId: string; value: number }[],
  reactions: ReactionSetup[],
  components: ComponentSetup[],
  tempK: number,
  simBasis: { limitingReactantId: string; desiredProductId: string }
): { [key: string]: number } | null => {
  const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
  if (!limitingReactant) {
    console.error("DEBUG: Limiting reactant not found in component list.");
    return null;
  }

  // --- Start of Debug Logging ---
  console.log("--- [Debug] Starting solveForInitialFlows ---");
  console.log(`Target: ${targetProductionRate.toFixed(4)} mol/s of ${targetComponentName}`);
  console.log(`At Reactor Volume: ${reactorVolume} m³ | Temp: ${tempK} K`);

  const errorFunction = (limitingReactantFlowGuess: number): number => {
    if (limitingReactantFlowGuess <= 0) return 1e6;

    const initialFlowRates = components.reduce((acc, comp) => {
        if (comp.id === limitingReactant.id) { acc[comp.name] = limitingReactantFlowGuess; } 
        else {
            const ratioInfo = molarRatios.find(r => r.numeratorId === comp.id);
            acc[comp.name] = ratioInfo ? limitingReactantFlowGuess * ratioInfo.value : 0;
        }
        return acc;
    }, {} as {[key: string]: number});
    
    // Self-contained, single-point calculation
    let outletFlows: { [key: string]: number };
    const v0 = 1.0; // Basis for liquid phase
    // Map ComponentSetup[] to Component[] for parseParallelReactions
    const componentsForParsing = components.map(comp => ({
      ...comp,
      initialFlowRate: (comp as any).initialFlowRate !== undefined ? (comp as any).initialFlowRate : '0'
    }));
    const kineticsForParsing = {
      rateConstantInputMode: 'arrhenius' as const,
      kValue: '0',
      AValue: '0',
      EaValue: '0',
      reactionTempK: tempK.toString()
    };
    // Convert ReactionSetup[] to ReactionData[]
    const reactionsData = reactions.map(r => {
      const reactants = components
        .filter(c => parseFloat(c.reactionData[r.id]?.stoichiometry || '0') < 0)
        .map(c => {
          const s = Math.abs(parseFloat(c.reactionData[r.id]?.stoichiometry || '0'));
          return s === 1 ? c.name : `${s} ${c.name}`;
        })
        .join(' + ');
      const products = components
        .filter(c => parseFloat(c.reactionData[r.id]?.stoichiometry || '0') > 0)
        .map(c => {
          const s = parseFloat(c.reactionData[r.id]?.stoichiometry || '0');
          return s === 1 ? c.name : `${s} ${c.name}`;
        })
        .join(' + ');
      return {
        id: r.id,
        reactants,
        products,
        AValue: r.AValue,
        EaValue: r.EaValue,
        AValueBackward: r.AValueBackward,
        EaValueBackward: r.EaValueBackward,
        isEquilibrium: r.isEquilibrium,
      };
    });
    const parsed = parseParallelReactions('', componentsForParsing, kineticsForParsing, reactionsData, false);

    if (reactorType === 'CSTR') {
        outletFlows = solveCSTRParallel(parsed, reactorVolume, initialFlowRates, v0, 'Liquid', 1, tempK);
    } else {
        const componentNames = parsed.allInvolvedComponents.map(c => c.name);
        const profile = solvePFR_ODE_System(parsed, reactorVolume, initialFlowRates, v0, componentNames);
        outletFlows = profile.length > 0 ? profile[profile.length - 1].flowRates : initialFlowRates;
    }
    
    const calculatedProduction = outletFlows[targetComponentName] || 0;
    return calculatedProduction - targetProductionRate;
  };

  // --- Robust Secant Method with Logging ---
  let x0 = 0.1;
  let x1 = 100.0; 
  let f0 = errorFunction(x0);
  let f1 = errorFunction(x1);
  const maxIter = 50;
  const tol = 1e-6;

  console.log(`Initial Guesses: f(${x0.toFixed(2)}) = ${f0.toFixed(4)}, f(${x1.toFixed(2)}) = ${f1.toFixed(4)}`);

  if (f0 * f1 > 0) {
      console.warn("Solver Warning: The target production rate may be unachievable. The error is the same sign at both low and high feed rates.");
  }

  for (let i = 0; i < maxIter; i++) {
    if (Math.abs(f1) < tol) {
        console.log(`%c[Debug] Solver CONVERGED in ${i + 1} iterations.`, 'color: #4CAF50; font-weight: bold;');
        const finalFlowRate = x1;
        return components.reduce((acc, comp) => {
            if (comp.id === limitingReactant.id) { acc[comp.name] = finalFlowRate; } 
            else {
                const ratioInfo = molarRatios.find(r => r.numeratorId === comp.id);
                acc[comp.name] = ratioInfo ? finalFlowRate * ratioInfo.value : 0;
            }
            return acc;
        }, {} as {[key: string]: number});
    }

    const denominator = f1 - f0;
    if (Math.abs(denominator) < 1e-12) {
        console.warn(`[Debug] Solver STOPPED: Denominator is too small. The function is too flat.`);
        return null;
    }
    
    const x2 = x1 - f1 * (x1 - x0) / denominator;
    if (x2 <= 0) {
        console.warn(`[Debug] Solver STOPPED: Next guess for feed rate was zero or negative (${x2.toFixed(4)}).`);
        return null;
    }
    
    x0 = x1; f0 = f1; x1 = x2; f1 = errorFunction(x1);
  }

  console.error(`[Debug] Solver FAILED to converge after ${maxIter} iterations.`);
  console.log(`Final State: Guess = ${x1.toFixed(4)}, Error = ${f1.toExponential(3)}`);
  return null;
};

// --- Update calculateCstrData to return outletFlowsAtV ---
const calculateCstrData = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number },
    v0: number,
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    V_max?: number
) => {
    // 1. Parse reactions into a generic structure the solver can use
    const allInvolvedComponents = components.filter(c => Object.values(c.reactionData).some(d => d.stoichiometry && parseFloat(d.stoichiometry) !== 0));
    const allComponentNames = allInvolvedComponents.map(c => c.name);
    const parsedReactions = {
        reactions: reactions.map(reactionInfo => {
            const reactants: any[] = [];
            const products: any[] = [];
            allInvolvedComponents.forEach(comp => {
                const reactionData = comp.reactionData?.[reactionInfo.id];
                if (reactionData?.stoichiometry && parseFloat(reactionData.stoichiometry) !== 0) {
                    const parsedComp = {
                        name: comp.name,
                        stoichiometricCoefficient: Math.abs(parseFloat(reactionData.stoichiometry)),
                        reactionOrderNum: parseFloat(reactionData.order || '0')
                    };
                    if (parseFloat(reactionData.stoichiometry) < 0) reactants.push(parsedComp);
                    else products.push(parsedComp);
                }
            });
            const A = parseFloat(reactionInfo.AValue);
            const Ea = parseFloat(reactionInfo.EaValue);
            return { reactants, products, rateConstantAtT: A * Math.exp(-Ea / (R * tempK)) };
        }),
        allInvolvedComponents: allInvolvedComponents.map(c => ({ name: c.name, id: c.id })),
    };

    let solverGuess: number[] | undefined = undefined;
    const nPoints = 100;
    const maxVolume = 2000;
    const minVolume = 0.01;
    const logMinV = Math.log10(minVolume);
    const logMaxV = Math.log10(maxVolume);
    const logStep = (logMaxV - logMinV) / (nPoints - 1);
    const selectivityData: { x: number; y: number }[] = [];
    const volumeData: { x: number; y: number }[] = [];
    let outletFlowsAtV: { [key: string]: number } = {};
    let closestDiff = Infinity;
    
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);

    if (!limitingReactant || !desiredProduct) return { selectivityData, volumeData, outletFlowsAtV: {} };

    const F_A0 = initialFlowRates[limitingReactant.name] || 0;
    if (F_A0 <= 1e-9) return { selectivityData, volumeData, outletFlowsAtV: {} };

    // 2. Loop over a logarithmic range of reactor volumes
    for (let i = 0; i < nPoints; i++) {
        const logV = logMinV + i * logStep;
        const V = Math.pow(10, logV);
        
        // 3. Call the solver, using the last result as the initial guess
        const outletFlowRates = solveCSTRParallel(
            parsedReactions, V, initialFlowRates, v0, 'Liquid', 1, tempK, solverGuess
        );
        
        solverGuess = parsedReactions.allInvolvedComponents.map((c:any) => outletFlowRates[c.name] || 0);

        // 4. Calculate metrics for the graph
        const F_A_final = outletFlowRates[limitingReactant.name] || 0;
        const conversion = (F_A0 - F_A_final) / F_A0;

        const molesReactantConsumed = F_A0 - F_A_final;
        const molesProductFormed = (outletFlowRates[desiredProduct.name] || 0) - (initialFlowRates[desiredProduct.name] || 0);
        
        // This simplified selectivity is valid for many common reaction schemes
        const selectivity = molesReactantConsumed > 1e-9 ? molesProductFormed / molesReactantConsumed : 0;

        if (conversion >= 0.001 && conversion <= 0.999) {
            selectivityData.push({ x: conversion, y: Math.max(0, selectivity) });
            volumeData.push({ x: conversion, y: V });
        }
    }
    
    selectivityData.sort((a, b) => a.x - b.x);
    volumeData.sort((a, b) => a.x - b.x);

    return { selectivityData, volumeData, outletFlowsAtV };
};

// --- Update calculatePfrData to return outletFlowsAtV ---
const calculatePfrData = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    V_max: number, // This will now be the dynamically calculated max volume
    initialFlowRates: { [key: string]: number },
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string }
) => {
    // 1. Parse reactions into a generic structure the solver can use
    const allInvolvedComponents = components.filter(c => Object.values(c.reactionData).some(d => d.stoichiometry && parseFloat(d.stoichiometry) !== 0));
    const allComponentNames = allInvolvedComponents.map(c => c.name);
    const parsedReactions = {
        reactions: reactions.map(reactionInfo => {
            const reactants: any[] = [];
            const products: any[] = [];
            allInvolvedComponents.forEach(comp => {
                const reactionData = comp.reactionData?.[reactionInfo.id];
                if (reactionData?.stoichiometry && parseFloat(reactionData.stoichiometry) !== 0) {
                    const parsedComp = {
                        name: comp.name,
                        stoichiometricCoefficient: Math.abs(parseFloat(reactionData.stoichiometry)),
                        reactionOrderNum: parseFloat(reactionData.order || '0')
                    };
                    if (parseFloat(reactionData.stoichiometry) < 0) reactants.push(parsedComp);
                    else products.push(parsedComp);
                }
            });
            const A = parseFloat(reactionInfo.AValue);
            const Ea = parseFloat(reactionInfo.EaValue);
            return { reactants, products, rateConstantAtT: A * Math.exp(-Ea / (R * tempK)) };
        }),
    };

    // 2. Define the PFR differential equations (dF/dV = r)
    const dFdV = (V: number, F: number[]): number[] => {
        const concentrations: { [key: string]: number } = {};
        const F_dict = allComponentNames.reduce((acc, name, i) => { acc[name] = F[i]; return acc; }, {} as {[key:string]:number});
        
        const v0 = 1.0; // Assuming constant liquid volumetric flow rate for rate calc
        allComponentNames.forEach(name => concentrations[name] = (F_dict[name] > 0 ? F_dict[name] : 0) / v0);

        const reactionRates = parsedReactions.reactions.map((reaction:any) => {
            let rate = reaction.rateConstantAtT || 0;
            reaction.reactants.forEach((reactant:any) => {
                rate *= Math.pow(concentrations[reactant.name], reactant.reactionOrderNum || 0);
            });
            return rate;
        });

        return allComponentNames.map(name => {
            let netRateOfFormation = 0; // This is r_j
            parsedReactions.reactions.forEach((reaction:any, j:number) => {
                const reactantInfo = reaction.reactants.find((r:any) => r.name === name);
                const productInfo = reaction.products.find((p:any) => p.name === name);
                if (reactantInfo) netRateOfFormation -= (reactantInfo.stoichiometricCoefficient || 0) * reactionRates[j];
                if (productInfo) netRateOfFormation += (productInfo.stoichiometricCoefficient || 0) * reactionRates[j];
            });
            return netRateOfFormation; // dFj/dV = r_j
        });
    };

    // MODIFICATION: Use a logarithmic scale for the volume span
    // The BDF solver will choose its own steps
    const V_span: [number, number] = [0, V_max];
    
    const F0 = allComponentNames.map(name => initialFlowRates[name] || 0);
    const results = solveODE_BDF(dFdV, F0, V_span);

    const selectivityData: { x: number; y: number }[] = [];
    const volumeData: { x: number; y: number }[] = [];
    let outletFlowsAtV: { [key: string]: number } = {};
    let closestDiff = Infinity;
    
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);
    if (!limitingReactant || !desiredProduct || !results.y) return { selectivityData, volumeData, outletFlowsAtV: {} };

    const F_A0 = initialFlowRates[limitingReactant.name] || 0;
    if (F_A0 <= 1e-9) return { selectivityData, volumeData, outletFlowsAtV: {} };

    results.t.forEach((v, i) => {
        const currentFlowRates = allComponentNames.reduce((acc, name, j) => {
            acc[name] = results.y[j][i];
            return acc;
        }, {} as { [key: string]: number });
        
        const F_A_current = currentFlowRates[limitingReactant.name] || 0;
        const conversion = (F_A0 - F_A_current) / F_A0;

        const molesReactantConsumed = F_A0 - F_A_current;
        const molesProductFormed = (currentFlowRates[desiredProduct.name] || 0) - (initialFlowRates[desiredProduct.name] || 0);
        
        const selectivity = molesReactantConsumed > 1e-9 ? molesProductFormed / molesReactantConsumed : 0;
        
        if (conversion >= 0.001 && conversion <= 0.999) {
            selectivityData.push({ x: conversion, y: Math.max(0, selectivity) });
            volumeData.push({ x: conversion, y: v });
        }
        // Track the closest volume to targetVolume
        // if (targetVolume !== undefined) {
        //     const diff = Math.abs(v - targetVolume);
        //     if (diff < closestDiff) {
        //         closestDiff = diff;
        //         outletFlowsAtV = currentFlowRates;
        //     }
        // }
    });

    return { selectivityData, volumeData, outletFlowsAtV };
};

// --- Linear System Solver ---
/**
 * Solves a linear system of equations Ax = b using Gaussian elimination.
 */
function solveLinearSystem(A: number[][], b: number[]): number[] | null {
    const n = A.length;
    // Create copies to avoid modifying the original arrays
    const A_copy = A.map(row => [...row]);
    const b_copy = [...b];

    for (let i = 0; i < n; i++) {
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(A_copy[k][i]) > Math.abs(A_copy[maxRow][i])) maxRow = k;
        }
        [A_copy[i], A_copy[maxRow]] = [A_copy[maxRow], A_copy[i]];
        [b_copy[i], b_copy[maxRow]] = [b_copy[maxRow], b_copy[i]];
        
        if (Math.abs(A_copy[i][i]) < 1e-12) return null; // System is singular

        for (let k = i + 1; k < n; k++) {
            const factor = A_copy[k][i] / A_copy[i][i];
            b_copy[k] -= factor * b_copy[i];
            for (let j = i; j < n; j++) {
                A_copy[k][j] -= factor * A_copy[i][j];
            }
        }
    }

    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        let sum = 0;
        for (let j = i + 1; j < n; j++) sum += A_copy[i][j] * x[j];
        x[i] = (b_copy[i] - sum) / A_copy[i][i];
    }
    return x;
}

/**
 * Solves a system of ODEs using a custom BDF method, ideal for stiff systems.
 */
function solveODE_BDF(
  derivatives: (t: number, y: number[]) => number[],
  y0: number[],
  tSpan: [number, number],
): { t: number[], y: number[][] } {
  const [t0, tf] = tSpan;
  let t = t0;
  let y = [...y0];
  const t_out = [t0];
  const y_out = [y0];
  const n = y0.length;

  // MODIFICATION: The initial step 'h' is now capped at 1.0 to prevent a giant first jump.
  let h = Math.min((tf - t0) / 1000, 1.0) || 1e-5; 
  
  const min_h = 1e-12;
  const newton_tol = 1e-8;
  const newton_max_iter = 20;
  
  while (t < tf) {
    if (t + h > tf) h = tf - t;
    if (h < min_h) break;

    let y_next = [...y];
    let converged = false;

    // Newton's method to solve the implicit BDF equation
    for (let iter = 0; iter < newton_max_iter; iter++) {
        const f_next = derivatives(t + h, y_next);
        const G = y_next.map((val, i) => val - y[i] - h * f_next[i]);
        
        if (Math.sqrt(G.reduce((sum, val) => sum + val*val, 0)) < newton_tol) {
            converged = true;
            break;
        }
        
        const J_f: number[][] = Array.from({length: n}, () => Array(n).fill(0));
        const h_eps = 1e-8;
        for(let j=0; j<n; j++){
            const y_h = [...y_next];
            y_h[j] += h_eps;
            const fx_h = derivatives(t+h, y_h);
            for(let i=0; i<n; i++){
                J_f[i][j] = (fx_h[i] - f_next[i]) / h_eps;
            }
        }
        
        const J_G = Array.from({length: n}, (_, i) => 
            Array.from({length: n}, (_, j) => (i === j ? 1 : 0) - h * J_f[i][j])
        );

        const delta_y = solveLinearSystem(J_G, G.map(val => -val));

        if (!delta_y) break;
        y_next = y_next.map((val, i) => val + delta_y[i]);
    }
    
    if (converged) {
        y = y_next.map(v => Math.max(0, v)); // Ensure no negative flows
        t += h;
        t_out.push(t);
        y_out.push([...y]);
        h = Math.min(h * 1.2, tf - t); // Increase step size
    } else {
        h *= 0.5; // Reduce step size if convergence fails
    }
  }

  // Transpose results for easier processing
  const y_transposed: number[][] = Array.from({ length: n }, () => []);
  for (let i = 0; i < t_out.length; i++) {
    for (let j = 0; j < n; j++) {
      y_transposed[j][i] = y_out[i][j];
    }
  }
  return { t: t_out, y: y_transposed };
}

// --- CSTR Parallel Solver ---
function solveCSTRParallel(
  parsedReactions: any,
  V: number,
  initialFlowRates: { [key: string]: number },
  v0: number,
  reactionPhase: 'Liquid' | 'Gas',
  totalPressure_bar: number,
  temp_K: number,
  initialGuess?: number[]
): { [key: string]: number } {
  // CORRECTED: Get component order directly from the parsed data
  const componentNames = parsedReactions.allInvolvedComponents.map((c: any) => c.name);
  
  let F = initialGuess ? [...initialGuess] : componentNames.map((name: string) => initialFlowRates[name] || 0);
  const n = F.length;
  
  const max_iter = 50;
  const tol = 1e-9;

  const G = (current_F: number[]): number[] => {
    let concentrations: { [key: string]: number } = {};
    const F_dict = componentNames.reduce((acc: { [key: string]: number }, name: string, i: number) => { acc[name] = current_F[i]; return acc; }, {} as {[key:string]:number});

    if (reactionPhase === 'Liquid') {
        componentNames.forEach((name: string) => concentrations[name] = (F_dict[name] > 0 ? F_dict[name] : 0) / v0);
    } else {
        const F_total = current_F.reduce((sum: number, f: number) => sum + f, 0);
        if (F_total < 1e-9) componentNames.forEach((name: string) => concentrations[name] = 0);
        else {
          const P_kPa = totalPressure_bar * 100;
          componentNames.forEach((name: string) => concentrations[name] = (F_dict[name] > 0 ? F_dict[name] : 0) / F_total * P_kPa / (R * temp_K));
        }
    }

    const rates = parsedReactions.reactions.map((reaction: any) => {
      let rate = reaction.rateConstantAtT || 0;
      reaction.reactants.forEach((reactant: any) => {
        rate *= Math.pow(concentrations[reactant.name], reactant.reactionOrderNum || 0);
      });
      return rate;
    });

    return componentNames.map((name: string, i: number) => {
      let netRateOfFormation = 0;
      parsedReactions.reactions.forEach((reaction: any, j: number) => {
        const reactantInfo = reaction.reactants.find((r: any) => r.name === name);
        const productInfo = reaction.products.find((p: any) => p.name === name);
        if (reactantInfo) netRateOfFormation -= (reactantInfo.stoichiometricCoefficient || 0) * rates[j];
        if (productInfo) netRateOfFormation += (productInfo.stoichiometricCoefficient || 0) * rates[j];
      });
      return (initialFlowRates[name] || 0) - current_F[i] + netRateOfFormation * V;
    });
  };

  for (let iter = 0; iter < max_iter; iter++) {
    const G_current = G(F);
    if (Math.sqrt(G_current.reduce((s: number, v: number) => s + v * v, 0)) < tol) break;
    const J_G: number[][] = Array(n).fill(0).map(() => Array(n).fill(0));
    const epsilon = 1e-8;
    for (let j = 0; j < n; j++) {
        const F_plus_h = [...F];
        F_plus_h[j] += epsilon;
        const G_plus_h = G(F_plus_h);
        for (let i = 0; i < n; i++) J_G[i][j] = (G_plus_h[i] - G_current[i]) / epsilon;
    }
    const delta_F = solveLinearSystem(J_G, G_current.map((val: number) => -val));
    if (!delta_F) break;
    const dampingFactor = 0.7;
    F = F.map((val: number, i: number) => val + dampingFactor * delta_F[i]);
  }

  return componentNames.reduce((acc: { [key: string]: number }, name: string, index: number) => {
    acc[name] = F[index] > 0 ? F[index] : 0;
    return acc;
  }, {} as { [key: string]: number });
}

// --- React Components ---

const KineticsInput = ({ onNext, reactionsSetup, setReactionsSetup, componentsSetup, setComponentsSetup, simBasis, setSimBasis }: { 
    onNext: () => void,
    reactionsSetup: ReactionSetup[],
    setReactionsSetup: React.Dispatch<React.SetStateAction<ReactionSetup[]>>,
    componentsSetup: ComponentSetup[],
    setComponentsSetup: React.Dispatch<React.SetStateAction<ComponentSetup[]>>,
    simBasis: any,
    setSimBasis: any
}) => {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';
  const cardBg = 'bg-card';
  const cardFg = 'text-card-foreground';
  const mainBg = 'bg-background';
  const mainFg = 'text-foreground';
  const mutedFg = 'text-muted-foreground';

  const addReactionSetup = () => {
    if (reactionsSetup.length >= 4) return
    const newId = (Math.max(0, ...reactionsSetup.map(r => parseInt(r.id))) + 1).toString()
    setReactionsSetup([...reactionsSetup, { id: newId, AValue: '1e6', EaValue: '50000', isEquilibrium: false }])
  }

  const removeReactionSetup = (idToRemove: string) => {
    if (reactionsSetup.length <= 1) return;
    setReactionsSetup(reactionsSetup.filter(r => r.id !== idToRemove))
    setComponentsSetup(componentsSetup.map(comp => {
      const newReactionData = { ...comp.reactionData }
      delete newReactionData[idToRemove]
      return { ...comp, reactionData: newReactionData }
    }))
  }

  // Update updateReactionSetup to handle boolean and string fields
  const updateReactionSetup = (id: string, field: keyof ReactionSetup, value: any) => {
    setReactionsSetup(reactionsSetup.map(reaction =>
      reaction.id === id ? { ...reaction, [field]: value } : reaction
    ))
  }

  const addComponentSetup = () => {
    const newId = `comp-${Date.now()}`
    setComponentsSetup([...componentsSetup, { id: newId, name: '', reactionData: {} }])
  }

  const removeComponentSetup = (idToRemove: string) => {
    if (componentsSetup.length <= 2) return;
    setComponentsSetup(componentsSetup.filter(c => c.id !== idToRemove))
  }

  const handleComponentSetupChange = (id: string, field: keyof ComponentSetup, value: any) => {
      setComponentsSetup(prev => prev.map(c => c.id === id ? { ...c, [field]: value } : c))
  }
  
  // --- MODIFICATION START ---
  // Validate that every reaction has at least one reactant and one product.
  const validationResult = useMemo(() => {
    for (const reaction of reactionsSetup) {
        let hasReactant = false;
        let hasProduct = false;

        // Check all components for their role in this specific reaction
        for (const comp of componentsSetup) {
            const stoichValue = parseFloat(comp.reactionData[reaction.id]?.stoichiometry || '0');
            if (stoichValue < 0) {
                hasReactant = true;
            } else if (stoichValue > 0) {
                hasProduct = true;
            }
        }

        // If a reaction is incomplete, invalidate the form.
        if (!hasReactant || !hasProduct) {
            return {
                isValid: false,
                message: `Reaction ${reaction.id} must have at least one reactant (negative stoichiometry) and one product (positive stoichiometry).`
            };
        }
    }

    // If all reactions are valid
    return { isValid: true, message: '' };
  }, [reactionsSetup, componentsSetup]);
  // --- MODIFICATION END ---
  
  const ReactionPreview = useMemo(() => {
    const generatePreview = (reaction: ReactionSetup) => {
        const reactants = componentsSetup
            .filter(c => c.name && parseFloat(c.reactionData[reaction.id]?.stoichiometry || '0') < 0)
            .map(c => {
                const stoich = Math.abs(parseFloat(c.reactionData[reaction.id]?.stoichiometry || '0'));
                return stoich === 1 ? c.name : `${stoich} ${c.name}`;
            })
            .join(' + ');

        const products = componentsSetup
            .filter(c => c.name && parseFloat(c.reactionData[reaction.id]?.stoichiometry || '0') > 0)
            .map(c => {
                const stoich = parseFloat(c.reactionData[reaction.id]?.stoichiometry || '0');
                return stoich === 1 ? c.name : `${stoich} ${c.name}`;
            })
            .join(' + ');

        const rateLawReactants = componentsSetup
            .filter(c => c.name && parseFloat(c.reactionData[reaction.id]?.order || '0') !== 0)
            .map(c => {
                const order = c.reactionData[reaction.id]?.order;
                const exponent = (order && order !== '1') ? <sup>{order}</sup> : null;
                return (
                    <React.Fragment key={c.id}>
                        {' '}C<sub>{c.name}</sub>
                        {exponent}
                    </React.Fragment>
                );
            });

        // MODIFICATION: Conditionally render the reverse exponent
        const rateLawProducts = reaction.isEquilibrium ? componentsSetup
            .filter(c => c.name && parseFloat(c.reactionData[reaction.id]?.orderReverse || '0') !== 0)
            .map(c => {
                const reverseOrder = c.reactionData[reaction.id]?.orderReverse;
                const exponent = (reverseOrder && reverseOrder !== '1') ? <sup>{reverseOrder}</sup> : null;
                return (
                    <React.Fragment key={c.id}>
                        {' '}C<sub>{c.name}</sub>
                        {exponent}
                    </React.Fragment>
                );
            }) : [];

        return {
            equation: `${reactants || '...'} ${reaction.isEquilibrium ? '⇌' : '→'} ${products || '...'}`,
            rateLaw: rateLawReactants,
            rateLawReverse: rateLawProducts
        };
    };
    
    return (
        <Card>
            <CardHeader>
                {/* MODIFICATION: Ensure title class matches others */}
                <CardTitle>Parsed Reactions & Rate Laws</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
                {reactionsSetup.map(r => {
                    const preview = generatePreview(r)
                    return (
                        <div key={r.id} className="p-3 bg-card rounded-md">
                          <p className="font-bold text-primary">Rxn {r.id}:</p>
                            <p className="font-medium text-center mb-1">{preview.equation}</p>
                          <p className="font-mono text-sm text-center text-muted-foreground w-full">
                            Rate = k
                            {/* Conditionally add ',f' to the subscript if reversible */}
                            <sub>{r.id}{r.isEquilibrium ? ',f' : ''}</sub>
                            {preview.rateLaw}
                            {/* Show the reverse rate term if reversible */}
                            {r.isEquilibrium && preview.rateLawReverse.length > 0 && (
                                <> - k<sub>{r.id},r</sub>{preview.rateLawReverse}</>
                            )}
                            </p>
                        </div>
                    )
                })}
            </CardContent>
        </Card>
    )
  }, [reactionsSetup, componentsSetup])

  return (
    <div className={`min-h-screen flex flex-col items-center justify-center p-4 sm:p-8 ${mainBg} ${mainFg}`}>
      <div className="container mx-auto space-y-6">
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            <div className={`p-6 rounded-lg shadow-lg ${cardBg} ${cardFg}`}> 
              <Card>
                <CardHeader>
                  <CardTitle className="flex items-center justify-between">
                    Reactions & Kinetics
                    {/* This button was missing */}
                    <Button variant="outline" size="sm" onClick={addReactionSetup} disabled={reactionsSetup.length >= 6}>
                        <PlusCircle className="mr-2 h-4 w-4" />
                        Add Reaction
                    </Button>
                  </CardTitle>
                </CardHeader>
                <CardContent className="space-y-6">
    {reactionsSetup.map((reaction, index) => (
      <div key={reaction.id} className={`space-y-3 ${index > 0 ? 'border-t pt-6' : ''}`}>
        <div className="flex items-center justify-between">
            <Label className="font-bold">Reaction {reaction.id}</Label>
            {/* This ensures the remove button shows for all but the first reaction */}
            {index > 0 && 
              <Button variant="destructive" size="sm" onClick={() => removeReactionSetup(reaction.id)} aria-label="Remove">
                <Trash2 className="h-4 w-4" />
              </Button>
            }
        </div>
        {/* Forward Reaction Inputs */}
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="flex items-center gap-2">
              <Label className="text-sm whitespace-nowrap" htmlFor={`a-value-${reaction.id}`}>
                {reaction.isEquilibrium ? <span>A<sub className="relative -left-px">f</sub>:</span> : 'A:'}
              </Label>
              <Input id={`a-value-${reaction.id}`} type="text" value={reaction.AValue} onChange={(e) => updateReactionSetup(reaction.id, 'AValue', e.target.value)} />
            </div>
            <div className="flex items-center gap-2">
              <Label className="text-sm whitespace-nowrap" htmlFor={`ea-value-${reaction.id}`}>
                {reaction.isEquilibrium ? <span>Ea<sub className="relative -left-px">f</sub>:</span> : 'Ea:'}
              </Label>
              <Input id={`ea-value-${reaction.id}`} type="text" value={reaction.EaValue} onChange={(e) => updateReactionSetup(reaction.id, 'EaValue', e.target.value)} />
              <span className="text-xs text-muted-foreground">J/mol</span>
            </div>
        </div>
        {/* Conditionally show reverse reaction inputs */}
        {reaction.isEquilibrium && (
            <div className="grid grid-cols-1 md:grid-cols-2 gap-4 pt-2 border-t">
                <div className="flex items-center gap-2">
                  <Label className="text-sm whitespace-now-wrap" htmlFor={`a-value-rev-${reaction.id}`}> 
                    <span>A<sub className="relative -left-px">r</sub>:</span>
                  </Label>
                  <Input id={`a-value-rev-${reaction.id}`} type="text" value={reaction.AValueBackward || ''} onChange={(e) => updateReactionSetup(reaction.id, 'AValueBackward', e.target.value)} />
                </div>
                <div className="flex items-center gap-2">
                  <Label className="text-sm whitespace-now-wrap" htmlFor={`ea-value-rev-${reaction.id}`}> 
                    <span>Ea<sub className="relative -left-px">r</sub>:</span>
                  </Label>
                  <Input id={`ea-value-rev-${reaction.id}`} type="text" value={reaction.EaValueBackward || ''} onChange={(e) => updateReactionSetup(reaction.id, 'EaValueBackward', e.target.value)} />
                  <span className="text-xs text-muted-foreground">J/mol</span>
                </div>
            </div>
        )}
      </div>
    ))}
</CardContent>
              </Card>
            </div>
            {ReactionPreview}
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          <div className={`lg:col-span-2 p-6 rounded-lg shadow-lg ${cardBg} ${cardFg}`}> 
              <div className="mb-4">
                <CardTitle className="flex items-center justify-between">
                  Components
                  <Button variant="secondary" size="sm" onClick={addComponentSetup}>
                      <PlusCircle className="mr-2 h-4 w-4" />Add Component
                  </Button>
                </CardTitle>
              </div>
              <div className="overflow-x-auto">
                  {/* First header row: Main column names and Reaction X */}
                  <div className="grid gap-2 items-center text-xs font-medium text-muted-foreground pb-2 border-b" style={{gridTemplateColumns: `40px 1fr 1fr repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`}}>
                      <div></div> {/* Spacer for trashcan icon */}
                      <div className="text-center">Name</div>
                      <div className="text-center">Initial Flow (mol/s)</div>
                      {reactionsSetup.map(r => (
                          <div key={r.id} className="text-center col-span-2">Reaction {r.id}</div>
                      ))}
                  </div>
                  {/* Second header row: Stoich. and Order */}
                  <div className="grid gap-2 items-center text-xs font-medium text-muted-foreground pb-2 border-b" style={{gridTemplateColumns: `40px 1fr 1fr repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`}}>
                      <div></div><div></div><div></div> {/* Spacers */}
                      {reactionsSetup.map(r => (
                          <React.Fragment key={r.id}>
                              <div className="text-center">Stoich.</div>
                              <div className="text-center">Order</div>
                          </React.Fragment>
                      ))}
                  </div>

                  {/* Data rows for each component */}
                  {componentsSetup.map((comp, index) => (
                      <div key={comp.id} className="grid gap-2 items-center py-2 relative" style={{gridTemplateColumns: `40px 1fr 1fr repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`}}>
                          {/* Trashcan button */}
                          {index > 0 ? (
                            <Button variant="ghost" size="icon" className="h-8 w-8" onClick={() => removeComponentSetup(comp.id)}>
                                <Trash2 className="h-4 w-4 text-red-500" />
                            </Button>
                          ) : (
                            <Button variant="ghost" size="icon" className="h-8 w-8" disabled>
                              <Trash2 className="h-4 w-4 text-muted-foreground" />
                            </Button>
                          )}

                          {/* Name Input */}
                          <div className="relative">
                              <Input 
                                  value={comp.name} 
                                  onChange={e => {
                                      handleComponentSetupChange(comp.id, 'name', e.target.value)
                                      // fetchSuggestions(e.target.value) // Removed undefined function
                                  }}
                                  // ... other props
                              />
                              {/* ... suggestion box logic ... */}
                          </div>

                          {/* Initial Flow Input */}
                          <Input type="number" value={comp.initialFlowRate} onChange={e => handleComponentSetupChange(comp.id, 'initialFlowRate', e.target.value)} placeholder="0.0" className="h-8 text-center"/>
                          
                          {/* Stoich. and Order inputs for each reaction */}
                          {reactionsSetup.map(r => (
                              <React.Fragment key={r.id}>
                                  <Input type="number" placeholder="0" value={comp.reactionData[r.id]?.stoichiometry || ''} onChange={e => handleComponentSetupChange(comp.id, 'reactionData', {...comp.reactionData, [r.id]: {...comp.reactionData[r.id], stoichiometry: e.target.value}})} className="h-8 text-center" step="0.1" />
                                  <Input type="number" placeholder="0" value={comp.reactionData[r.id]?.order || ''} onChange={e => handleComponentSetupChange(comp.id, 'reactionData', {...comp.reactionData, [r.id]: {...comp.reactionData[r.id], order: e.target.value}})} className="h-8 text-center" step="0.1" />
                              </React.Fragment>
                          ))}
                      </div>
                  ))}
              </div>
          </div>
          <div className={`p-6 rounded-lg shadow-lg ${cardBg} ${cardFg}`}>
            <Card>
                <CardHeader>
                    <CardTitle>Simulation Basis</CardTitle>
                </CardHeader>
                <CardContent className="space-y-4">
                    {/* Limiting Reactant Dropdown */}
                    <div>
                        <Label>Limiting Reactant</Label>
                        <Select 
                            value={simBasis.limitingReactantId} 
                            onValueChange={(value) => setSimBasis({...simBasis, limitingReactantId: value})}
                        >
                            <SelectTrigger className="w-full mt-1">
                                <SelectValue placeholder="Select a reactant" />
                            </SelectTrigger>
                            <SelectContent>
                                {componentsSetup
                                    .filter(c => {
                                        const isReactant = Object.values(c.reactionData).some(d => parseFloat(d.stoichiometry || '0') < 0);
                                        const isProduct = Object.values(c.reactionData).some(d => parseFloat(d.stoichiometry || '0') > 0);
                                        return isReactant && !isProduct;
                                    })
                                    .map(c => (
                                        <SelectItem key={c.id} value={c.id}>{c.name}</SelectItem>
                                    ))
                                }
                            </SelectContent>
                        </Select>
                    </div>
                    {/* Desired Product Dropdown */}
                    <div>
                        <Label>Desired Product</Label>
                        <Select 
                            value={simBasis.desiredProductId} 
                            onValueChange={(value) => setSimBasis({...simBasis, desiredProductId: value})}
                        >
                            <SelectTrigger className="w-full mt-1">
                                <SelectValue placeholder="Select a product" />
                            </SelectTrigger>
                            <SelectContent>
                                {componentsSetup
                                    .filter(c => Object.values(c.reactionData).some(d => parseFloat(d.stoichiometry || '0') > 0))
                                    .map(c => (
                                        <SelectItem key={c.id} value={c.id}>{c.name}</SelectItem>
                                    ))
                                }
                            </SelectContent>
                        </Select>
                    </div>
                </CardContent>
            </Card>
          </div>
        </div>

        <div className="text-right mt-2">
            <div title={!validationResult.isValid ? validationResult.message : undefined} className="inline-block">
                {/* Use the Button component for the primary action */}
                <Button onClick={onNext} disabled={!validationResult.isValid}>
                Next: Configure Simulation →
                </Button>
            </div>
        </div>
      </div>
    </div>
  )
}

// Interface for the new molarRatios prop
interface MolarRatio {
    numeratorId: string;
    value: number;
}

const ReactorSimulator = ({ onBack, reactions, components, simBasis, molarRatios, setMolarRatios }: { 
    onBack: () => void, 
    reactions: ReactionSetup[], 
    components: ComponentSetup[],
    simBasis: any,
    molarRatios: MolarRatio[],
    setMolarRatios: React.Dispatch<React.SetStateAction<MolarRatio[]>>
}) => {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';
  const textColor = isDark ? '#E5E7EB' : '#1F2937';
  const cardBg = isDark ? 'bg-card' : 'bg-card'; // Tailwind + theme
  const cardFg = isDark ? 'text-card-foreground' : 'text-card-foreground';
  const mainBg = isDark ? 'bg-background' : 'bg-background';
  const mainFg = isDark ? 'text-foreground' : 'text-foreground';

  const [reactorType, setReactorType] = useState<ReactorType>('PFR')
  const [graphType, setGraphType] = useState<GraphType>('selectivity')
  const [temperature, setTemperature] = useState(4)
  
  const [molarRatioMin, setMolarRatioMin] = useState('2');
  const [molarRatioMax, setMolarRatioMax] = useState('20');
  const [tempMin, setTempMin] = useState('4');
  const [tempMax, setTempMax] = useState('20');

  const [prodRate, setProdRate] = useState('100');
  const [prodMolarMass, setProdMolarMass] = useState('114.23');
  const [reactorVolume, setReactorVolume] = useState('100'); // New state for reactor volume

  const tempK = temperature + 273.15
  
  const generateGraphData = useCallback(() => {
    const v0 = 1.0; // Assume a constant volumetric flow rate for the liquid phase basis
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    if (!limitingReactant) return { series: [], xAxis: '', yAxis: '', legend: [] };

    // 1. Use a fixed basis for the initial simulation run (e.g., 1 mol/s of limiting reactant).
    const F_A0_BASIS = 1.0;
    const basisFlowRates = components.reduce((acc, comp) => {
        if (comp.id === limitingReactant.id) {
            acc[comp.name] = F_A0_BASIS;
        } else {
            const ratioInfo = molarRatios.find(r => r.numeratorId === comp.id);
            // As MR increases, this flow rate increases, correctly increasing total flow.
            acc[comp.name] = ratioInfo ? F_A0_BASIS * ratioInfo.value : 0;
        }
        return acc;
    }, {} as {[key: string]: number});

    // 2. Dynamically find the required V_max to ensure the graph reaches ~99% conversion on the basis run.
    let V_max;
    if (reactorType === 'PFR') {
        V_max = findVolumeForConversion_PFR(0.99, reactions, components, basisFlowRates, tempK, simBasis);
    } else { // CSTR
        V_max = findVolumeForConversion_CSTR(0.99, reactions, components, basisFlowRates, v0, tempK, simBasis);
    }
    
    // 3. Run the appropriate solver to get the performance curves for the basis feed rate.
    const basisResults = reactorType === 'CSTR' 
        ? calculateCstrData(reactions, components, basisFlowRates, v0, tempK, simBasis, V_max)
        : calculatePfrData(reactions, components, V_max, basisFlowRates, tempK, simBasis);
    
    // 4. Scale the volume data based on the desired production rate from the UI.
    const targetProduction_mol_s = (parseFloat(prodRate) * 1e6) / parseFloat(prodMolarMass) / (365 * 24 * 3600);
    
    const scaledVolumeData = basisResults.volumeData.map(point => {
        const conversion = point.x;
        const basisVolume = point.y;
        
        // Find corresponding selectivity by interpolating from the simulation results
        const sData = basisResults.selectivityData;
        let selectivity = 0;
        if (sData.length > 0) {
            const p1 = sData.findLast(p => p.x <= conversion) ?? sData[0];
            const p2 = sData.find(p => p.x >= conversion) ?? sData[sData.length - 1];
            selectivity = interpolate(conversion, p1, p2);
        }

        if (conversion < 1e-6 || selectivity < 1e-6) return null; // Avoid division by zero

        // Calculate the actual feed rate required for the target production
        const F_A0_actual = targetProduction_mol_s / (selectivity * conversion);
        // Determine how much to scale the basis volume by
        const scalingFactor = F_A0_actual / F_A0_BASIS;
        const actualVolume = basisVolume * scalingFactor;
        
        return { x: conversion, y: actualVolume };
    }).filter((p): p is {x: number, y: number} => p !== null && isFinite(p.y));

    // --- The rest of the function formats the data for the ECharts component ---
    const dataToShow = {
        selectivityData: basisResults.selectivityData,
        volumeData: scaledVolumeData,
    };
    
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);
    const xLabel = `Conversion of ${limitingReactant?.name || 'Limiting Reactant'}`;
    let yLabel = '';
    let legend: string[] = [];
    let series: any[] = [];
    
    if (graphType === 'volume') {
      yLabel = 'Reactor Volume (m³)';
      legend = ['Volume'];
      series = [{ name: legend[0], type: 'line', data: dataToShow.volumeData.map(d => [d.x, d.y]), smooth: true, showSymbol: false, lineStyle: { width: 2 } }];
    } else { 
      yLabel = `Selectivity to ${desiredProduct?.name || 'Product'}`;
      legend = [desiredProduct?.name || 'Product'];
      series = [{ name: legend[0], type: 'line', data: dataToShow.selectivityData.map(d => [d.x, d.y]), smooth: true, showSymbol: false, lineStyle: { width: 2 } }];
    }
    
    return { series, xAxis: xLabel, yAxis: yLabel, legend };
  }, [reactorType, molarRatios, temperature, graphType, tempK, reactions, components, simBasis, prodRate, prodMolarMass]);

  const graphData = generateGraphData();

  const chartOptions: EChartsOption = {
    animation: false,
    backgroundColor: 'transparent',
    title: {
        text: `${reactorType} Performance`,
        left: 'center',
        top: '2%',
        textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' }
    },
    grid: { left: '5%', right: '3%', bottom: '8%', top: '10%', containLabel: true },
    xAxis: {
        type: 'value',
        name: graphData.xAxis,
        nameLocation: 'middle', nameGap: 30,
        nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
        min: 0, max: 1,
        axisLine: { lineStyle: { color: textColor } },
        axisLabel: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false },
    },
    yAxis: {
        type: 'value',
        name: graphData.yAxis,
        nameLocation: 'middle',
        nameGap: 50,
        nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
        min: 0, 
        
        // MODIFICATION: Set max to 1 for selectivity/conversion, otherwise auto-scale
        max: ['selectivity', 'selectivityVsConversion', 'conversion'].includes(graphType) ? 1 : 'dataMax',
        
        axisLine: { lineStyle: { color: textColor } },
        axisLabel: { 
            color: textColor, 
            fontSize: 12, 
            fontFamily: 'Merriweather Sans',
            formatter: (value: number) => value.toPrecision(2), // Use toPrecision for better formatting
            showMaxLabel: true 
        },
        splitLine: { show: false },
    },
    legend: {
      orient: 'horizontal', bottom: -5, left: 'center',
      textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
      data: graphData.legend,
      itemWidth: 25, itemHeight: 2
    },
    tooltip: { 
      trigger: 'axis',
      backgroundColor: isDark ? 'rgba(31, 41, 55, 0.9)' : '#f9fafb',
      borderColor: isDark ? '#4B5563' : '#d1d5db',
      textStyle: { color: textColor, fontSize: 12, fontFamily: 'Merriweather Sans' },
      
      // 1. Format the tooltips that appear on the axes
      axisPointer: {
        type: 'cross',
        label: {
          backgroundColor: isDark ? '#374151' : '#e5e7eb', // A dark grey color
          formatter: function (params) {
            // Use toPrecision(3) for 3 significant figures
            return parseFloat(params.value as string).toPrecision(3);
          }
        }
      },

      // 2. Format the main data tooltip that appears on hover
      formatter: (params: any) => {
        if (!params || params.length === 0) {
            return '';
        }

        // Get the x-axis name (e.g., "Conversion") and format its value
        const xAxisLabel = graphData.xAxis.split(' ')[0]; 
        const xAxisValue = parseFloat(params[0].axisValue).toPrecision(3);
        
        let tooltipContent = `<strong>${xAxisLabel}</strong>: ${xAxisValue}<br/>`;

        // Add a line for each data series (e.g., Volume)
        params.forEach((p: any) => {
            const seriesName = p.seriesName;
            const seriesValue = parseFloat(p.value[1]).toPrecision(3);
            const marker = p.marker; // The colored dot icon
            tooltipContent += `${marker} ${seriesName}: <strong>${seriesValue}</strong>`;
        });

        return tooltipContent;
      }
    },
    series: graphData.series
  };
  
  const limitingReactantName = components.find(c => c.id === simBasis.limitingReactantId)?.name || 'Limiting';

  // Move handleMolarRatioChange here so it can access setMolarRatios from props
  const handleMolarRatioChange = (numeratorId: string, newValue: number) => {
    setMolarRatios((prevRatios: MolarRatio[]) => 
      prevRatios.map((ratio: MolarRatio) => 
        ratio.numeratorId === numeratorId ? { ...ratio, value: newValue } : ratio
      )
    );
  };

  return (
    <div className={`min-h-screen flex flex-col p-4 ${mainBg} ${mainFg}`}>
      {/* MODIFICATION: Wrap the grid in a snapping container */}
      <div className="container mx-auto flex-grow">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 h-full">
        
            {/* Control Panel Card */}
            <div className={`lg:col-span-1 p-6 rounded-lg shadow-lg flex flex-col space-y-6 ${cardBg} ${cardFg}`}>
            <div className="flex justify-start">
                 <button onClick={onBack} className="bg-slate-700 hover:bg-slate-600 text-white font-bold py-2 px-4 rounded-lg transition-colors flex items-center">
                    <ArrowLeft className="h-4 w-4 mr-2"/> Back
                </button>
            </div>
                
                {/* MODIFICATION: Add the input boxes back */}
            <div className="space-y-2">
                <label className="text-sm font-medium">Desired Production Rate (kta)</label>
                    <input type="number" value={prodRate} onChange={e => setProdRate(e.target.value)} className="w-full bg-background text-foreground p-2 rounded"/>
            </div>
            <div className="space-y-2">
                <label className="text-sm font-medium">Product Molar Mass (g/mol)</label>
                    <input type="number" value={prodMolarMass} onChange={e => setProdMolarMass(e.target.value)} className="w-full bg-background text-foreground p-2 rounded"/>
            </div>
                
                {/* Molar ratio sliders using ParameterSlider */}
                {molarRatios.map(ratio => {
                    const numeratorName = components.find(c => c.id === ratio.numeratorId)?.name || 'N/A';
                    return (
                        <div key={ratio.numeratorId} className="pt-4 border-t border-border">
                            <ParameterSlider
                                /* MODIFICATION: Parentheses removed from the label */
                                label={`${numeratorName} / ${limitingReactantName}`}
                                unit="" // Unit is now part of the label
                                value={ratio.value.toString()}
                                onValueChange={(value) => handleMolarRatioChange(ratio.numeratorId, value)}
                                min={parseFloat(molarRatioMin)}
                                max={parseFloat(molarRatioMax)}
                                step={1}
                                maxSliderValue={molarRatioMax}
                                onMaxSliderChange={setMolarRatioMax}
                                minSliderValue={molarRatioMin}
                                onMinSliderChange={setMolarRatioMin}
                            />
                </div>
                    );
                })}

                {/* Temperature slider using ParameterSlider */}
                <ParameterSlider
                    label="Temperature"
                    unit="°C"
                    value={temperature.toString()}
                    onValueChange={(value) => setTemperature(value)}
                    min={parseFloat(tempMin)}
                    max={parseFloat(tempMax)}
                    step={1}
                    maxSliderValue={tempMax}
                    onMaxSliderChange={setTempMax}
                    minSliderValue={tempMin}
                    onMinSliderChange={setTempMin}
                />
            </div>

            {/* Graph Card */}
            <div className={`lg:col-span-2 p-4 rounded-lg shadow-lg aspect-square ${cardBg} ${cardFg}`}> 
              {/* Chart section is unchanged */}
            <div className="relative mb-2">
                <div className="absolute top-0 left-0 z-10">
                        {/* Wrap the buttons in a muted background for a nice toggle effect */}
                        <div className="flex items-center gap-1 rounded-lg p-1 bg-muted">
                            <Button
                                onClick={() => setReactorType('PFR')}
                                variant={reactorType === 'PFR' ? 'default' : 'ghost'}
                                size="sm"
                                className="text-xs px-3 py-1 h-auto"
                            >
                                PFR
                            </Button>
                            <Button
                                onClick={() => setReactorType('CSTR')}
                                variant={reactorType === 'CSTR' ? 'default' : 'ghost'}
                                size="sm"
                                className="text-xs px-3 py-1 h-auto"
                            >
                                CSTR
                            </Button>
                    </div>
                </div>
                <div className="absolute top-0 right-0 z-10">
                        <select onChange={(e) => setGraphType(e.target.value as GraphType)} value={graphType} className={`p-1 rounded-md text-xs ${isDark ? 'bg-muted' : 'bg-muted'}`}> 
                        <option value="selectivity">Selectivity vs. Conversion</option>
                        <option value="volume">Volume vs. Conversion</option>
                    </select>
                </div>
            </div>
          <ReactECharts
            // MODIFICATION: Add reactorType to the key to force a full re-render on switch
            key={`${resolvedTheme}-${reactorType}`}
            // If you have echartsRef and echarts, include them as props
            // ref={echartsRef}
            // echarts={echarts}
            option={chartOptions}
            style={{ height: '100%', width: '100%', minHeight: '450px' }}
            notMerge={true}
          />
            </div>
            
        </div>
      </div>
    </div>
  )
}

// Place this new component before your main Page component
interface ParameterSliderProps {
  label: string;
  value: string;
  unit: string;
  onValueChange: (value: number) => void;
  min: number;
  max: number;
  step: number;
  maxSliderValue: string;
  onMaxSliderChange: (value: string) => void;
  minSliderValue: string;
  onMinSliderChange: (value: string) => void;
  isReadOnly?: boolean;
}

const ParameterSlider: React.FC<ParameterSliderProps> = ({
  label,
  value,
  unit,
  onValueChange,
  min,
  max,
  step,
  maxSliderValue,
  onMaxSliderChange,
  minSliderValue,
  onMinSliderChange,
  isReadOnly = false
}) => {
  return (
    <div className="space-y-2">
      <div className="flex items-center justify-between">
        <Label htmlFor={`${label}-slider`} className="text-sm font-medium">
          {label}: {value} {unit}
        </Label>
        <div className="flex items-center gap-1">
          <Label htmlFor={`${label}-min-input`} className="text-xs text-muted-foreground">Min:</Label>
          {/* MODIFICATION: Changed to your preferred width of w-14 */}
          <Input id={`${label}-min-input`} type="number" value={minSliderValue} onChange={(e) => onMinSliderChange(e.target.value)} className="w-14 h-8 text-xs" />
          <Label htmlFor={`${label}-max-input`} className="text-xs text-muted-foreground">Max:</Label>
          {/* MODIFICATION: Changed to your preferred width of w-14 */}
          <Input id={`${label}-max-input`} type="number" value={maxSliderValue} onChange={(e) => onMaxSliderChange(e.target.value)} className="w-14 h-8 text-xs" readOnly={isReadOnly} />
        </div>
      </div>
      <Slider
        id={`${label}-slider`}
        min={min}
        max={max}
        step={step}
        value={[parseFloat(value)]}
        onValueChange={(val) => onValueChange(val[0])}
        className="w-full"
      />
    </div>
  );
};

// This helper function finds the PFR volume required for a target conversion
const findVolumeForConversion_PFR = (
    targetConversion: number,
    // Pass in all the necessary data for a simulation run
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number },
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string }
): number => {
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    if (!limitingReactant) return 2000; // Default fallback

    // This "error function" will be used by our solver
    const errorFunction = (volumeGuess: number): number => {
        if (volumeGuess <= 0) return -targetConversion; // Return large error for invalid input

        // Run a PFR simulation up to the guessed volume
        const profile = calculatePfrData(reactions, components, volumeGuess, initialFlowRates, tempK, simBasis);
        
        // Find the maximum conversion achieved in the simulation
        let maxConversion = 0;
        if (profile.selectivityData.length > 0) {
            maxConversion = profile.selectivityData.reduce((max, p) => p.x > max ? p.x : max, 0);
        }

        return maxConversion - targetConversion;
    };

    // Simple solver (bisection method) to find the volume
    let low = 0.1;
    let high = 50000; // Start with a very large max volume
    const maxIter = 50;

    for (let i = 0; i < maxIter; i++) {
        const mid = (low + high) / 2;
        const error = errorFunction(mid);
        if (Math.abs(error) < 0.01) return mid; // Close enough
        if (error < 0) {
            low = mid; // Conversion is too low, need more volume
        } else {
            high = mid; // Conversion is too high, need less volume
        }
    }
    return high; // Return the best guess if not fully converged
};

// This helper function finds the CSTR volume required for a target conversion
const findVolumeForConversion_CSTR = (
    targetConversion: number,
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number },
    v0: number,
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string }
): number => {
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    if (!limitingReactant) return 2000; // Default fallback

    const errorFunction = (volumeGuess: number): number => {
        if (volumeGuess <= 0) return -targetConversion;

        // Use the same parsing logic as in calculateCstrData
        const allInvolvedComponents = components.filter(c => Object.values(c.reactionData).some(d => d.stoichiometry && parseFloat(d.stoichiometry) !== 0));
        const parsedReactions = {
            reactions: reactions.map(reactionInfo => {
                const reactants: any[] = [];
                const products: any[] = [];
                allInvolvedComponents.forEach(comp => {
                    const reactionData = comp.reactionData?.[reactionInfo.id];
                    if (reactionData?.stoichiometry && parseFloat(reactionData.stoichiometry) !== 0) {
                        const parsedComp = {
                            name: comp.name,
                            stoichiometricCoefficient: Math.abs(parseFloat(reactionData.stoichiometry)),
                            reactionOrderNum: parseFloat(reactionData.order || '0')
                        };
                        if (parseFloat(reactionData.stoichiometry) < 0) reactants.push(parsedComp);
                        else products.push(parsedComp);
                    }
                });
                const A = parseFloat(reactionInfo.AValue);
                const Ea = parseFloat(reactionInfo.EaValue);
                return { reactants, products, rateConstantAtT: A * Math.exp(-Ea / (8.314 * tempK)) };
            }),
            allInvolvedComponents: allInvolvedComponents.map(c => ({ name: c.name, id: c.id })),
        };
        const outletFlows = solveCSTRParallel(parsedReactions, volumeGuess, initialFlowRates, v0, 'Liquid', 1, tempK);
        
        const F_A0 = initialFlowRates[limitingReactant.name] || 0;
        const F_A_final = outletFlows[limitingReactant.name] || 0;
        const currentConversion = (F_A0 > 1e-9) ? (F_A0 - F_A_final) / F_A0 : 0;

        return currentConversion - targetConversion;
    };

    // Use a simple bisection solver to find the required volume
    let low = 0.1;
    let high = 50000; // Start with a very large max volume to ensure we find the target
    const maxIter = 50;

    for (let i = 0; i < maxIter; i++) {
        const mid = (low + high) / 2;
        if (mid === low || mid === high) break; // No progress
        const error = errorFunction(mid);
        if (Math.abs(error) < 0.01) return mid; // Target found
        if (error < 0) {
            low = mid; // Conversion is too low, need more volume
        } else {
            high = mid; // Conversion is too high, need less volume
        }
    }
    return high; // Return the best guess
};

export default function Home() {
  const [showSimulator, setShowSimulator] = useState(false)
  
  const [reactionsSetup, setReactionsSetup] = useState<ReactionSetup[]>(initialReactions);
  const [componentsSetup, setComponentsSetup] = useState<ComponentSetup[]>(initialComponents);
  
  const [simBasis, setSimBasis] = useState({
    limitingReactantId: 'comp-A',
    desiredProductId: 'comp-C',
  });

  const [molarRatios, setMolarRatios] = useState<Array<{numeratorId: string, value: number}>>([]);

  // --- MODIFICATION START ---
  // This function now correctly identifies ONLY primary reactants (chemicals that are
  // consumed but not also produced), which fixes the bug.
  const getReactants = useCallback((components: ComponentSetup[], reactions: ReactionSetup[]) => {
      const componentRoles = new Map<string, {isReactant: boolean, isProduct: boolean}>();

      // First, determine if each component is a reactant or a product in any reaction
      for (const comp of components) {
          const roles = { isReactant: false, isProduct: false };
          for (const reaction of reactions) {
              const stoich = parseFloat(comp.reactionData[reaction.id]?.stoichiometry || '0');
              if (stoich < 0) {
                  roles.isReactant = true;
              } 
              else if (stoich > 0) {
                  roles.isProduct = true;
              }
          }
          componentRoles.set(comp.id, roles);
      }

      // A true "reactant" for our purposes is one that is consumed but never produced (i.e., not an intermediate)
      const reactantIds = new Set<string>();
      for (const [id, roles] of componentRoles.entries()) {
          if (roles.isReactant && !roles.isProduct) {
              reactantIds.add(id);
          }
      }
      return components.filter(comp => reactantIds.has(comp.id));
  }, []);
  // --- MODIFICATION END ---

  // 4. Effect to automatically manage the molar ratios based on the limiting reactant
  useEffect(() => {
    const allReactants = getReactants(componentsSetup, reactionsSetup);
    const otherReactants = allReactants.filter(r => r.id !== simBasis.limitingReactantId);

    // Update the molarRatios state, preserving existing values where possible
    setMolarRatios(prevRatios => {
        const newRatios = otherReactants.map(reactant => {
            const existingRatio = prevRatios.find(r => r.numeratorId === reactant.id);
            return {
                numeratorId: reactant.id,
                value: existingRatio ? existingRatio.value : 12.0 // Default value for new ratios
            };
        });
        
        // Prevents unnecessary re-renders if the ratios haven't changed
        if (JSON.stringify(newRatios) !== JSON.stringify(prevRatios)) {
             return newRatios;
        }
        return prevRatios;
    });

  }, [componentsSetup, reactionsSetup, simBasis.limitingReactantId, getReactants]);


  // 5. Updated effect to handle component deletion
  useEffect(() => {
    const componentIds = new Set(componentsSetup.map(c => c.id));
    const allReactants = getReactants(componentsSetup, reactionsSetup);

    if (!componentIds.has(simBasis.limitingReactantId)) {
        setSimBasis(b => ({...b, limitingReactantId: allReactants[0]?.id || ''}))
    }
    if (!componentIds.has(simBasis.desiredProductId)) {
        // Fallback to the first component if the desired product is deleted
        setSimBasis(b => ({...b, desiredProductId: componentsSetup[0]?.id || ''}))
    }
  }, [componentsSetup, simBasis.limitingReactantId, simBasis.desiredProductId, getReactants, reactionsSetup]);

  if (showSimulator) {
    return <ReactorSimulator 
              onBack={() => setShowSimulator(false)} 
              reactions={reactionsSetup} 
              components={componentsSetup} 
              simBasis={simBasis}
              // Pass the new molar ratio state down
              molarRatios={molarRatios}
              setMolarRatios={setMolarRatios}
            />
  } else {
    return <KineticsInput 
              onNext={() => setShowSimulator(true)} 
              reactionsSetup={reactionsSetup}
              setReactionsSetup={setReactionsSetup}
              componentsSetup={componentsSetup}
              setComponentsSetup={setComponentsSetup}
              simBasis={simBasis}
              setSimBasis={setSimBasis}
            />
  }
}