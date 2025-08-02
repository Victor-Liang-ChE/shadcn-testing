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
import { Switch } from '@/components/ui/switch';

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
  molarMass?: string // New
  density?: string // New
}

type ReactorType = 'PFR' | 'CSTR'
type GraphType = 'selectivity' | 'volume' | 'conversion' | 'flowrates' | 'composition' | 'selectivityVsConversion';

// Default Data (now used as initial state for the editable form)
const initialReactions: ReactionSetup[] = [
  { id: '1', AValue: '3.66e14', EaValue: '101600', isEquilibrium: false },
  { id: '2', AValue: '4.77e16', EaValue: '110850', isEquilibrium: false },
]

const initialComponents: ComponentSetup[] = [
  { id: 'comp-A', name: '1-Butene', molarMass: '56.11', density: '595', reactionData: { '1': { stoichiometry: '-1', order: '1' }, '2': { stoichiometry: '-1', order: '1' } } },
  { id: 'comp-B', name: 'Isobutane', molarMass: '58.12', density: '593', reactionData: { '1': { stoichiometry: '-1', order: '1' }, '2': { stoichiometry: '0', order: '0' } } },
  { id: 'comp-C', name: 'Isooctane', molarMass: '114.23', density: '692', reactionData: { '1': { stoichiometry: '1', order: '0' }, '2': { stoichiometry: '-1', order: '1' } } },
  { id: 'comp-D', name: 'Dodecane', molarMass: '170.34', density: '750', reactionData: { '1': { stoichiometry: '0', order: '0' }, '2': { stoichiometry: '1', order: '0' } } },
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

// Helper function to perform linear interpolation from an array of data points
const interpolateFromData = (
  targetX: number,
  dataPoints: { x: number; y: number }[]
): number => {
  if (dataPoints.length < 2) {
    return NaN // Not enough data to interpolate
  }

  // Ensure data is sorted by x-value for reliable searching
  const sortedPoints = [...dataPoints].sort((a, b) => a.x - b.x)

  // Find the two points that bracket the targetX
  let p1: { x: number; y: number } | null = null
  let p2: { x: number; y: number } | null = null

  if (targetX < sortedPoints[0].x) {
    // If target is before the first point, we can't reliably get a value
    return NaN
  }
  
  if (targetX > sortedPoints[sortedPoints.length - 1].x) {
    // If target is after the last point, also return NaN
    return NaN
  }

  for (let i = 0; i < sortedPoints.length - 1; i++) {
    if (sortedPoints[i].x <= targetX && sortedPoints[i + 1].x >= targetX) {
      p1 = sortedPoints[i]
      p2 = sortedPoints[i + 1]
      break
    }
  }

  if (!p1 || !p2) {
    return NaN // Should not happen if logic is correct
  }

  // Handle vertical line case to prevent division by zero
  if (p1.x === p2.x) {
    return p1.y
  }

  // Standard linear interpolation formula
  return p1.y + ((targetX - p1.x) * (p2.y - p1.y)) / (p2.x - p1.x)
}

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
    return null;
  }

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
        outletFlows = solveCSTRParallel(parsed, components, reactorVolume, initialFlowRates, 'Liquid', 1, tempK);
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

  if (f0 * f1 > 0) {
      return null;
  }

  for (let i = 0; i < maxIter; i++) {
    if (Math.abs(f1) < tol) {
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
        return null;
    }
    
    const x2 = x1 - f1 * (x1 - x0) / denominator;
    if (x2 <= 0) {
        return null;
    }
    
    x0 = x1; f0 = f1; x1 = x2; f1 = errorFunction(x1);
  }

  return null;
};

// Replace your existing calculateCstrData function with this one.
const calculateCstrData = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number }, // This will be treated as a basis
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    reactionPhase: 'Liquid' | 'Gas',
    pressure_bar: number,
    targetProduction_kta: number
) => {
    // This function also follows the same pattern:
    // 1. Run the CSTR simulation loop with high resolution to generate raw data.
    // 2. Interpolate that raw data onto the common conversion grid.
    
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId)
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId)
    if (!limitingReactant || !desiredProduct) return { selectivityData: [], volumeData: [] }

    // --- 1. Generate high-resolution raw data from the simulation ---
    const rawSelectivityData: { x: number; y: number }[] = []
    const rawVolumeData: { x: number; y: number }[] = []
    
    const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
    const P_C_target_mols = (targetProduction_kta * 1e9) / productMolarMass / (365 * 24 * 3600);

    const parsedReactions = {
        reactions: reactions.map(reactionInfo => {
            const reactants: any[] = [];
            const products: any[] = [];
            components.forEach(comp => {
                const reactionData = comp.reactionData?.[reactionInfo.id];
                if (reactionData?.stoichiometry && parseFloat(reactionData.stoichiometry) !== 0) {
                    const parsedComp = { name: comp.name, stoichiometricCoefficient: Math.abs(parseFloat(reactionData.stoichiometry)), reactionOrderNum: parseFloat(reactionData.order || '0') };
                    if (parseFloat(reactionData.stoichiometry) < 0) reactants.push(parsedComp);
                    else products.push(parsedComp);
                }
            });
            const A = parseFloat(reactionInfo.AValue); const Ea = parseFloat(reactionInfo.EaValue);
            return { 
                id: reactionInfo.id,
                reactants, 
                products, 
                rateConstantAtT: A * Math.exp(-Ea / (8.314 * tempK)) 
            };
        }),
        allInvolvedComponents: components.map(c => ({ name: c.name, id: c.id })),
    };

    const getRateOfFormation = (current_F: { [key: string]: number }): number => {
        let concentrations: { [key: string]: number } = {};
        const F_vec = parsedReactions.allInvolvedComponents.map((c:any) => current_F[c.name] || 0)
        
        if (reactionPhase === 'Liquid') {
            let v = 0;
            components.forEach(comp => {
                const molarMass = parseFloat(comp.molarMass || '1');
                const density = parseFloat(comp.density || '1000');
                if (density > 0) v += (current_F[comp.name] * molarMass) / density;
            });
            if (v < 1e-9) v = 1e-9;
            parsedReactions.allInvolvedComponents.forEach((c:any) => concentrations[c.name] = Math.max(0, current_F[c.name]) / v);
        } else { // Gas Phase
            const F_total = F_vec.reduce((sum, f) => sum + f, 0);
            if (F_total < 1e-9) return 0;
            const P_Pa = pressure_bar * 1e5;
            parsedReactions.allInvolvedComponents.forEach((c:any) => concentrations[c.name] = Math.max(0, current_F[c.name]) / F_total * P_Pa / (R * tempK));
        }

        const reactionRates = parsedReactions.reactions.map((reaction: any) => {
            let rate = reaction.rateConstantAtT || 0;
            reaction.reactants.forEach((reactant: any) => {
                rate *= Math.pow(concentrations[reactant.name], reactant.reactionOrderNum || 0);
            });
            return rate;
        });

        let r_C = 0; // Net rate of formation for the desired product
        parsedReactions.reactions.forEach((reaction: any, j: number) => {
            const stoich = parseFloat(desiredProduct.reactionData[reaction.id]?.stoichiometry || '0');
            r_C += stoich * reactionRates[j];
        });
        return r_C;
    };

    // **FIX START**: Calculate initial selectivity and initialize data arrays correctly.
    const initialRates = calculateNetFormationRates(initialFlowRates, reactions, components, tempK, reactionPhase, pressure_bar);
    const r_P_initial = initialRates[desiredProduct.name] || 0;
    const r_A_initial = initialRates[limitingReactant.name] || 0;
    const initialSelectivity = (r_A_initial < -1e-9) ? (r_P_initial / -r_A_initial) : 1.0;

    // **FIX END**

    const F_A0_basis = initialFlowRates[limitingReactant.name] || 0;
    let solverGuess = parsedReactions.allInvolvedComponents.map((c:any) => initialFlowRates[c.name] || 0);

    const logV_min = -8;
    const logV_max = 8;
    const n_steps = 500; // Use more steps for higher resolution raw data
    const step_size = (logV_max - logV_min) / n_steps;

    for (let i = 0; i <= n_steps; i++) {
        const logV_basis = logV_min + i * step_size;
        const V_basis = Math.pow(10, logV_basis);

        const outletFlows_basis = solveCSTR_Robust(
            parsedReactions, components, V_basis, initialFlowRates, 
            reactionPhase, pressure_bar, tempK, solverGuess
        );
        
        solverGuess = parsedReactions.allInvolvedComponents.map((c:any) => outletFlows_basis[c.name] || 0);

        const F_A_final_basis = outletFlows_basis[limitingReactant.name] || 0;
        
        // **FIX START**: Filter out points with negligible conversion to prevent bad data.
        const molesReactantConsumed = F_A0_basis - F_A_final_basis;
        if (molesReactantConsumed < 1e-7) {
            continue;
        }
        // **FIX END**

        const conversion = molesReactantConsumed / F_A0_basis;
        const r_C = getRateOfFormation(outletFlows_basis);
        let V_true = (r_C > 1e-9) ? P_C_target_mols / r_C : Infinity;
        if (reactionPhase === 'Liquid') {
            V_true *= 0.001; // Convert L to m³
        }
        
        const molesProductFormed = (outletFlows_basis[desiredProduct.name] || 0) - (initialFlowRates[desiredProduct.name] || 0);
        
        // **FIX START**: Calculate selectivity directly now that we've filtered.
        const selectivity = molesProductFormed / molesReactantConsumed;
        // **FIX END**

        if (conversion >= 1e-6 && conversion <= 0.9999 && V_true < Infinity) {
             rawSelectivityData.push({ x: conversion, y: Math.max(0, selectivity) });
             rawVolumeData.push({ x: conversion, y: V_true });
        }

        if (conversion > 0.999) break;
    }

    // --- 2. Create the common grid and interpolate onto it ---
    const commonConversionGrid = Array.from({ length: 100 }, (_, i) => i * 0.01)
    
    const selectivityData = commonConversionGrid.map(conv => ({
        x: conv,
        y: interpolateFromData(conv, rawSelectivityData)
    })).filter(p => !isNaN(p.y))

    const volumeData = commonConversionGrid.map(conv => ({
        x: conv,
        y: interpolateFromData(conv, rawVolumeData)
    })).filter(p => !isNaN(p.y))

    return { selectivityData, volumeData }
}

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
  const newton_max_iter = 100;
  
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

// --- CSTR Solver (Updated with Debugging) ---
function solveCSTRParallel(
  parsedReactions: any,
  components: ComponentSetup[],
  V: number,
  initialFlowRates: { [key: string]: number },
  reactionPhase: 'Liquid' | 'Gas',
  totalPressure_bar: number,
  temp_K: number,
  initialGuess?: number[],
  debug = false // Debug flag
): { [key: string]: number } {
  const componentNames = parsedReactions.allInvolvedComponents.map((c: any) => c.name);
  let F = initialGuess ? [...initialGuess] : componentNames.map((name: string) => initialFlowRates[name] || 0);
  const n = F.length;
  const max_iter = 100;
  const tol = 1e-9;

  if (debug) {
    console.groupCollapsed(`[CSTR Solver Debug] V=${V.toExponential(2)}`);
    console.log("Initial Guess (F):", F.map((f: number) => f.toExponential(3)));
    console.log("Initial Flow Rates (F0):", initialFlowRates);
    console.log("Temperature (K):", temp_K);
  }

  const G = (current_F: number[]): number[] => {
    let concentrations: { [key: string]: number } = {};
    const F_dict = componentNames.reduce((acc: { [key: string]: number }, name: string, i: number) => { acc[name] = current_F[i]; return acc; }, {} as {[key:string]:number});

    if (reactionPhase === 'Liquid') {
        let v = 0;
        componentNames.forEach((name: string) => {
            const compInfo = components.find(c => c.name === name);
            const molarMass = parseFloat(compInfo?.molarMass || '1');
            const density = parseFloat(compInfo?.density || '1000');
            if (density > 0) {
                v += (F_dict[name] * molarMass) / density;
            }
        });
        if (v < 1e-9) v = 1e-9;
        componentNames.forEach((name: string) => concentrations[name] = (F_dict[name] > 0 ? F_dict[name] : 0) / v);
    } else {
        const F_total = current_F.reduce((sum: number, f: number) => sum + f, 0);
        if (F_total < 1e-9) {
            componentNames.forEach((name: string) => concentrations[name] = 0);
        } else {
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
    const errorNorm = Math.sqrt(G_current.reduce((s, v) => s + v * v, 0));

    if (debug) {
        console.log(`--- Iteration ${iter} ---`);
        console.log("Current F:", F.map((f: number) => f.toExponential(3)));
        console.log("Error Norm:", errorNorm.toExponential(3));
    }
    
    if (errorNorm < tol) {
        if (debug) {
            console.log(`%cCONVERGED in ${iter + 1} iterations.`, 'color: #4CAF50; font-weight: bold;');
            console.groupEnd();
        }
        break;
    }
    
    const J_G: number[][] = Array(n).fill(0).map(() => Array(n).fill(0));
    for (let j = 0; j < n; j++) {
        const F_plus_h = [...F];
        const h = (Math.abs(F[j]) || 1e-8) * 1e-7;
        F_plus_h[j] += h;
        const G_plus_h = G(F_plus_h);
        for (let i = 0; i < n; i++) J_G[i][j] = (G_plus_h[i] - G_current[i]) / h;
    }

    const delta_F = solveLinearSystem(J_G, G_current.map((val: number) => -val));
    
    if (!delta_F) {
        if (debug) {
            console.warn(`%cLinear solve failed at iter ${iter}. Jacobian may be singular.`, 'color: #F44336;');
            console.log("Jacobian (J_G):", J_G);
            console.groupEnd();
        }
        break;
    }
    
    const dampingFactor = 0.8;
    F = F.map((val: number, i: number) => val + dampingFactor * delta_F[i]);
    F = F.map((val: number) => Math.max(1e-12, val));
    
    if (iter === max_iter - 1) {
        if (debug) {
            console.error(`%cFAILED to converge after ${max_iter} iterations.`, 'color: #F44336; font-weight: bold;');
            console.groupEnd();
        }
    }
  }

  return componentNames.reduce((acc: { [key: string]: number }, name: string, index: number) => {
    acc[name] = F[index] > 0 ? F[index] : 0;
    return acc;
  }, {} as { [key: string]: number });
}

// Add this new, more robust solver function to your file.
const solveCSTR_Robust = (
  parsedReactions: any,
  components: ComponentSetup[],
  V: number,
  initialFlowRates: { [key: string]: number },
  reactionPhase: 'Liquid' | 'Gas',
  totalPressure_bar: number,
  temp_K: number,
  initialGuess?: number[]
): { [key: string]: number } => {
  const componentNames = parsedReactions.allInvolvedComponents.map((c: any) => c.name);
  let F = initialGuess ? [...initialGuess] : componentNames.map((name: string) => initialFlowRates[name] || 0);
  const n = F.length;
  const max_iter = 100;
  const tol = 1e-9;

  // --- Helper function G(F) = F0 - F + r*V ---
  const G = (current_F: number[]): number[] => {
    let concentrations: { [key: string]: number } = {};
    const F_dict = componentNames.reduce((acc: { [key: string]: number }, name: string, i: number) => { acc[name] = current_F[i]; return acc; }, {} as {[key:string]:number});

    if (reactionPhase === 'Liquid') {
        let v = 0;
        componentNames.forEach((name: string) => {
            const compInfo = components.find(c => c.name === name);
            const molarMass = parseFloat(compInfo?.molarMass || '1');
            const density = parseFloat(compInfo?.density || '1000');
            if (density > 0) v += (F_dict[name] * molarMass) / density;
        });
        if (v < 1e-9) v = 1e-9;
        componentNames.forEach((name: string) => concentrations[name] = Math.max(0, F_dict[name]) / v);
    } else { // Gas Phase
        const F_total = current_F.reduce((sum: number, f: number) => sum + Math.max(0, f), 0);
        if (F_total < 1e-9) {
            componentNames.forEach((name: string) => concentrations[name] = 0);
        } else {
          const P_Pa = totalPressure_bar * 1e5;
          componentNames.forEach((name: string) => concentrations[name] = (Math.max(0, F_dict[name])) / F_total * P_Pa / (R * temp_K));
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
  // --- End of G(F) helper ---

  for (let iter = 0; iter < max_iter; iter++) {
    const G_current = G(F);
    const errorNorm = Math.sqrt(G_current.reduce((s, v) => s + v * v, 0));
    
    if (errorNorm < tol) break;
    
    // --- Calculate Jacobian and Newton Step (same as before) ---
    const J_G: number[][] = Array(n).fill(0).map(() => Array(n).fill(0));
    for (let j = 0; j < n; j++) {
        const F_plus_h = [...F];
        const h = (Math.abs(F[j]) || 1e-8) * 1e-7;
        F_plus_h[j] += h;
        const G_plus_h = G(F_plus_h);
        for (let i = 0; i < n; i++) J_G[i][j] = (G_plus_h[i] - G_current[i]) / h;
    }

    const delta_F = solveLinearSystem(J_G, G_current.map((val: number) => -val));
    if (!delta_F) break; // Jacobian is singular, can't proceed
    
    // --- NEW: Backtracking Line Search to ensure stability ---
    let alpha = 1.0; // Initial step size (full Newton step)
    const beta = 0.5; // Step reduction factor
    const c = 1e-4; // Armijo condition constant
    let F_new = F;
    let step_accepted = false;

    for (let ls_iter = 0; ls_iter < 10; ls_iter++) { // Max 10 line search steps
      F_new = F.map((val: number, i: number) => val + alpha * delta_F[i]);

      // Ensure physical solution (no negative flows)
      if (F_new.some((val: number) => val < 0)) {
        alpha *= beta;
        continue;
      }
      
      const G_new = G(F_new);
      const newErrorNorm = Math.sqrt(G_new.reduce((s, v) => s + v * v, 0));

      // Check for sufficient decrease (Armijo condition)
      if (newErrorNorm < (1 - c * alpha) * errorNorm) {
        step_accepted = true;
        break;
      }
      
      alpha *= beta; // Reduce step size and try again
    }
    
    if (step_accepted) {
      F = F_new;
    } else {
      break; // Line search failed, exit main loop
    }
    // --- End of Line Search ---
  }

  return componentNames.reduce((acc: { [key: string]: number }, name: string, index: number) => {
    acc[name] = F[index] > 0 ? F[index] : 0;
    return acc;
  }, {} as { [key: string]: number });
};

const calculatePfrGasPhase = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    tempK: number,
    pressure_bar: number,
    molarRatios: Array<{ numeratorId: string; value: number }>,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    targetProduction_kta: number
) => {
    // This function follows the exact same pattern as the liquid phase function above:
    // 1. Run the simulation with high resolution to generate raw data.
    // 2. Interpolate that raw data onto the common conversion grid.
    
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId)
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId)
    if (!limitingReactant || !desiredProduct) return { selectivityData: [], volumeData: [] }

    // --- 1. Generate high-resolution raw data from the simulation ---
    const rawSelectivityData: { x: number; y: number }[] = []
    const rawVolumeData: { x: number; y: number }[] = []
    
    const P_Pa = pressure_bar * 1e5;
    const limitingReactantIndex = components.findIndex(c => c.id === simBasis.limitingReactantId);
    const desiredProductIndex = components.findIndex(c => c.id === simBasis.desiredProductId);

    // --- FIX START: Corrected dF/dV function for gas phase ---
    const dFdV = (V: number, F_vec: number[]) => {
        const F_dict = components.reduce((acc, comp, i) => { acc[comp.name] = F_vec[i]; return acc; }, {} as { [key: string]: number });
        const F_tot = F_vec.reduce((sum, f) => sum + f, 0);
        if (F_tot < 1e-12) return Array(components.length).fill(0);

        // First, calculate concentrations in standard SI units (mol/m^3).
        const concentrations_m3 = components.reduce((acc, comp) => {
            acc[comp.name] = (F_dict[comp.name] * P_Pa) / (F_tot * R * tempK);
            return acc;
        }, {} as { [key: string]: number });

        // Since rate constants are in L-based units, convert concentrations to mol/L for calculations.
        const concentrations_L: { [key: string]: number } = {};
        for (const key in concentrations_m3) {
            concentrations_L[key] = concentrations_m3[key] * 1e-3; // Convert mol/m³ to mol/L
        }

        const reactionRates = reactions.map(reactionInfo => {
            // k is in L-based units (e.g., L/mol/s)
            const k_Liters = parseFloat(reactionInfo.AValue) * Math.exp(-parseFloat(reactionInfo.EaValue) / (R * tempK));
            let rate = k_Liters;
            components.forEach(comp => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
                
                // Rate law must only use reactant concentrations and L-based units.
                if (order > 0 && stoich < 0) {
                     rate *= Math.pow(Math.max(0, concentrations_L[comp.name]), order);
                }
            });
            return rate; // The rate is now correctly in mol/(L*s).
        });
        
        return components.map(comp => {
            let netRateOfFormation = 0;
            reactions.forEach((reactionInfo, j) => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                netRateOfFormation += stoich * reactionRates[j];
            });
            // With rate in mol/(L*s), the integrated volume will be in Liters.
            return netRateOfFormation;
        });
    };
    // --- FIX END ---

    const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
    const P_C_mol_s = (targetProduction_kta * 1e6) / productMolarMass / (365 * 24 * 3600);

    const F0_basis_dict: { [key: string]: number } = { [limitingReactant.name]: 1.0 };
    molarRatios.forEach(ratio => {
        const compName = components.find(c => c.id === ratio.numeratorId)?.name;
        if (compName) F0_basis_dict[compName] = ratio.value;
    });
    const F0_basis_vec = components.map(c => F0_basis_dict[c.name] || 0);
    const F_tot0_basis = F0_basis_vec.reduce((sum, f) => sum + f, 0);

    // --- Calculate initial selectivity ---
    const initialRates = calculateNetFormationRates(F0_basis_dict, reactions, components, tempK, 'Gas', pressure_bar);
    const r_P_initial = initialRates[desiredProduct.name] || 0;
    const r_A_initial = initialRates[limitingReactant.name] || 0;
    const initialSelectivity = (r_A_initial < -1e-9) ? (r_P_initial / -r_A_initial) : 1.0;

    // --- Start of existing simulation loop ---
    let F_current = F0_basis_vec;
    let V_current = 0;

    const logV_min = -8;
    const logV_max = 8; // Previously 3
    const n_steps = 500; // Use more steps for higher resolution raw data
    const v_points: number[] = []
    for (let i = 0; i <= n_steps; i++) {
        const logV = logV_min + i * (logV_max - logV_min) / n_steps
        v_points.push(Math.pow(10, logV))
    }

    for (const V_next of v_points) {
        if (V_next <= V_current) continue;
        
        const odeResults = solveODE_BDF(dFdV, F_current, [V_current, V_next]);
        
        if (!odeResults.y || odeResults.y[0].length < 2) continue;
        
        F_current = odeResults.y.map(comp_y => comp_y[comp_y.length - 1]);
        V_current = V_next;
        
        const F_A_current_val = F_current[limitingReactantIndex];
        const F_P_current_val = F_current[desiredProductIndex];
        
        const conversion = (F0_basis_vec[limitingReactantIndex] - F_A_current_val) / F0_basis_vec[limitingReactantIndex];
        const molesConsumed = F0_basis_vec[limitingReactantIndex] - F_A_current_val;
        const selectivity = molesConsumed > 1e-9 ? F_P_current_val / molesConsumed : 0;
        
        if (conversion >= 1e-6 && conversion < 0.9999 && selectivity > 1e-6) {
            // --- correct Eq. F-35 scaling ---
            const F_tot0_basis = F0_basis_vec.reduce((s,f)=>s+f,0);          // ΣF₀,basis
            const F_A0_true   = P_C_mol_s / (selectivity*conversion);  // Eq.F-10
            const F_tot0_true = F_A0_true*(1 + molarRatios
                                                 .reduce((s,r)=>s+r.value,0)); // Eq.F-11 + others
            const V_true      = 1e-3 * V_current * (F_tot0_true / F_tot0_basis);    // L → m³, Eq.F-35
            // --- correct Eq. F-35 scaling ---
            
            rawSelectivityData.push({ x: conversion, y: selectivity });
            rawVolumeData.push({ x: conversion, y: V_true });
        }
        
        if (conversion >= 0.999) break;
    }
    // --- End of existing simulation loop ---

    // --- 2. Create the common grid and interpolate onto it ---
    const commonConversionGrid = Array.from({ length: 100 }, (_, i) => i * 0.01)
    
    const selectivityData = commonConversionGrid.map(conv => ({
        x: conv,
        y: interpolateFromData(conv, rawSelectivityData)
    })).filter(p => !isNaN(p.y))

    const volumeData = commonConversionGrid.map(conv => ({
        x: conv,
        y: interpolateFromData(conv, rawVolumeData)
    })).filter(p => !isNaN(p.y))

    return { selectivityData, volumeData }
}

const calculatePfrLiquidPhase = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    tempK: number,
    molarRatios: Array<{ numeratorId: string; value: number }>,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    targetProduction_kta: number
) => {
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);
    if (!limitingReactant || !desiredProduct) return { selectivityData: [], volumeData: [] };

    const rawSelectivityData: { x: number; y: number }[] = [];
    const rawVolumeData: { x: number; y: number }[] = [];
    
    const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
    const P_C_target_mols = (targetProduction_kta * 1e9) / productMolarMass / (365 * 24 * 3600);

    const limitingReactantIndex = components.findIndex(c => c.id === simBasis.limitingReactantId);
    const desiredProductIndex = components.findIndex(c => c.id === simBasis.desiredProductId);

    // --- Basis Calculation Setup ---
    const F0_basis_dict: { [key: string]: number } = { [limitingReactant.name]: 1.0 };
    molarRatios.forEach(ratio => {
        const compName = components.find(c => c.id === ratio.numeratorId)?.name;
        if (compName) F0_basis_dict[compName] = ratio.value;
    });
    const F0_basis_vec = components.map(c => F0_basis_dict[c.name] || 0);

    // --- Calculate initial selectivity ---
    const initialRates = calculateNetFormationRates(F0_basis_dict, reactions, components, tempK, 'Liquid', 1); // Pressure=1 is arbitrary for liquid
    const r_P_initial = initialRates[desiredProduct.name] || 0;
    const r_A_initial = initialRates[limitingReactant.name] || 0;
    const initialSelectivity = (r_A_initial < -1e-9) ? (r_P_initial / -r_A_initial) : 1.0;

    // --- Calculate volumetric flow 'q' for the basis run (L/s) ---
    let q_basis = 0;
    components.forEach(comp => {
        const flow = F0_basis_dict[comp.name] || 0;
        if (flow > 0) {
            const molarMass = parseFloat(comp.molarMass || '1');
            const density = parseFloat(comp.density || '1000'); // g/L
            q_basis += (flow * molarMass) / density;
        }
    });
    if (q_basis < 1e-9) q_basis = 1e-9;
    const q_basis_L  = q_basis;

    // --- PFR Differential Equations (dF/dV = r) ---
    const dFdV = (V: number, F_vec: number[]) => {
        const concentrations = F_vec.map(f => f / q_basis_L); // mol/L
        const C_dict = components.reduce((acc, comp, i) => { acc[comp.name] = concentrations[i]; return acc; }, {} as { [key: string]: number });

        const reactionRates = reactions.map(reactionInfo => {
            const k_f = parseFloat(reactionInfo.AValue) * Math.exp(-parseFloat(reactionInfo.EaValue) / (R * tempK));
            let rate = k_f;
            components.forEach(comp => {
                const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
                if (order > 0 && parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0') < 0) {
                     rate *= Math.pow(Math.max(0, C_dict[comp.name]), order);
                }
            });
            return rate; // mol/(L*s)
        });

        return components.map(comp => {
            let netRateOfFormation = 0;
            reactions.forEach((reactionInfo, j) => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                netRateOfFormation += stoich * reactionRates[j];
            });
            return netRateOfFormation; // r_j in mol/(L*s)
        });
    };
    
    // --- Simulation Loop ---
    let F_current = F0_basis_vec;
    let V_current = 0;
    const logV_min = -8;
    const logV_max = 8;
    const n_steps = 500;
    const v_points: number[] = [];
    for (let i = 0; i <= n_steps; i++) {
        const logV = logV_min + i * (logV_max - logV_min) / n_steps;
        v_points.push(Math.pow(10, logV));
    }

    for (const V_next of v_points) {
        if (V_next <= V_current) continue;
        const odeResults = solveODE_BDF(dFdV, F_current, [V_current, V_next]);
        if (!odeResults.y || odeResults.y[0].length < 2) continue;
        
        F_current = odeResults.y.map(comp_y => comp_y[comp_y.length - 1]);
        V_current = V_next;
        
        const F_A_out_basis = F_current[limitingReactantIndex];
        const P_C_out_basis = F_current[desiredProductIndex];
        const conversion = (F0_basis_vec[limitingReactantIndex] - F_A_out_basis) / F0_basis_vec[limitingReactantIndex];
        const P_C_basis_mols = P_C_out_basis - (F0_basis_vec[desiredProductIndex] || 0);
        
        const selectivity = (F0_basis_vec[limitingReactantIndex] - F_A_out_basis > 1e-9) ? P_C_basis_mols / (F0_basis_vec[limitingReactantIndex] - F_A_out_basis) : 0;
        
        // --- Correct Volume Scaling Logic ---
        // 1. Calculate the true molar flow rate of the limiting reactant into the reactor.
        const F_A0_true = (selectivity > 1e-9 && conversion > 1e-9) ? P_C_target_mols / (selectivity * conversion) : 0;

        // 2. Calculate the true total volumetric flow rate into the reactor.
        let q_true = 0;
        const F0_true_dict: {[key: string]: number} = { [limitingReactant.name]: F_A0_true };
         molarRatios.forEach(ratio => {
            const compName = components.find(c => c.id === ratio.numeratorId)?.name;
            if (compName) F0_true_dict[compName] = F_A0_true * ratio.value;
        });

        components.forEach(comp => {
            const flow = F0_true_dict[comp.name] || 0;
            if (flow > 0) {
                const molarMass = parseFloat(comp.molarMass || '1');
                const density = parseFloat(comp.density || '1000');
                q_true += (flow * molarMass) / density;
            }
        });

        // 3. Scale the volume by the ratio of volumetric flows and convert from Liters to m³.
        const V_true = (q_basis > 1e-9 && q_true > 0) ? (1e-3 * V_current * (q_true / q_basis)) : Infinity;
        
        if (conversion >= 1e-6 && conversion < 0.9999 && V_true < Infinity && selectivity > 0) {
            rawSelectivityData.push({ x: conversion, y: selectivity });
            rawVolumeData.push({ x: conversion, y: V_true });
        }
        if (conversion >= 0.999) break;
    }

    const commonConversionGrid = Array.from({ length: 101 }, (_, i) => i * 0.01);
    
    const selectivityData = commonConversionGrid.map(conv => ({
        x: conv,
        y: interpolateFromData(conv, rawSelectivityData)
    })).filter(p => !isNaN(p.y) && p.y >= 0);

    const volumeData = commonConversionGrid.map(conv => ({
        x: conv,
        y: interpolateFromData(conv, rawVolumeData)
    })).filter(p => !isNaN(p.y) && p.y >= 0);

    return { selectivityData, volumeData };
};

const calculateNetFormationRates = (
    current_F: { [key: string]: number },
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    tempK: number,
    reactionPhase: 'Liquid' | 'Gas',
    pressure_bar: number
): { [key: string]: number } => {
    let concentrations: { [key: string]: number } = {};
    const R = 8.314;

    if (reactionPhase === 'Liquid') {
        let v_Liters = 0;
        components.forEach(comp => {
            const molarMass = parseFloat(comp.molarMass || '1');
            const density = parseFloat(comp.density || '1000');
            if (density > 0 && current_F[comp.name]) {
                v_Liters += (current_F[comp.name] * molarMass) / density;
            }
        });
        if (v_Liters < 1e-9) v_Liters = 1e-9;
        const v_m3 = v_Liters / 1000;
        components.forEach(comp => {
            concentrations[comp.name] = Math.max(0, current_F[comp.name] || 0) / v_m3;
        });
    } else { // Gas Phase
        const F_total = Object.values(current_F).reduce((sum, f) => sum + (f || 0), 0);
        if (F_total < 1e-9) return components.reduce((acc, c) => ({...acc, [c.name]: 0}), {});
        const P_Pa = pressure_bar * 1e5;
        components.forEach(comp => {
            concentrations[comp.name] = Math.max(0, current_F[comp.name] || 0) / F_total * P_Pa / (R * tempK);
        });
    }

    const reactionNetRates = reactions.map(reactionInfo => {
        const k_f = parseFloat(reactionInfo.AValue) * Math.exp(-parseFloat(reactionInfo.EaValue) / (R * tempK));
        let rate_f = k_f;
        components.forEach(comp => {
            const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
            const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
            if (stoich < 0 && order > 0) {
                rate_f *= Math.pow(Math.max(0, concentrations[comp.name]), order);
            }
        });
        return rate_f;
    });

    const componentRates: { [key: string]: number } = {};
    components.forEach(comp => {
        let netRateOfFormation = 0;
        reactions.forEach((reactionInfo, j) => {
            const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
            netRateOfFormation += stoich * reactionNetRates[j];
        });
        componentRates[comp.name] = netRateOfFormation;
    });

    return componentRates;
};

// --- React Components ---

const KineticsInput = ({ 
    onNext, 
    reactionsSetup, setReactionsSetup, 
    componentsSetup, setComponentsSetup, 
    simBasis, setSimBasis,
    reactionPhase, setReactionPhase, // New props for reaction phase
    prodRate, setProdRate // New props for production rate
}: { 
    onNext: () => void,
    reactionsSetup: ReactionSetup[],
    setReactionsSetup: React.Dispatch<React.SetStateAction<ReactionSetup[]>>,
    componentsSetup: ComponentSetup[],
    setComponentsSetup: React.Dispatch<React.SetStateAction<ComponentSetup[]>>,
    simBasis: any,
    setSimBasis: any,
    reactionPhase: 'Liquid' | 'Gas',
    setReactionPhase: React.Dispatch<React.SetStateAction<'Liquid' | 'Gas'>>,
    prodRate: string,
    setProdRate: React.Dispatch<React.SetStateAction<string>>
}) => {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';
  const cardBg = 'bg-card';
  const cardFg = 'text-card-foreground';
  const mainBg = 'bg-background';
  const mainFg = 'text-foreground';
  
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

  const handleComponentSetupChange = (id: string, field: 'name' | 'reactionData' | 'molarMass' | 'density', value: any) => {
      setComponentsSetup(prev => prev.map(c => c.id === id ? { ...c, [field]: value } : c))
  }
  
  const validationResult = useMemo(() => {
    for (const reaction of reactionsSetup) {
        let hasReactant = false;
        let hasProduct = false;

        for (const comp of componentsSetup) {
            const stoichValue = parseFloat(comp.reactionData[reaction.id]?.stoichiometry || '0');
            if (stoichValue < 0) {
                hasReactant = true;
            } else if (stoichValue > 0) {
                hasProduct = true;
            }
        }

        if (!hasReactant || !hasProduct) {
            return {
                isValid: false,
                message: `Reaction ${reaction.id} must have at least one reactant (negative stoichiometry) and one product (positive stoichiometry).`
            };
        }
    }

    return { isValid: true, message: '' };
  }, [reactionsSetup, componentsSetup]);
  
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
                <CardTitle>Parsed Reactions & Rate Laws</CardTitle>
            </CardHeader>
            <CardContent className="space-y-4">
                {reactionsSetup.map(r => {
                    const preview = generatePreview(r)
                    return (
                        <div key={r.id} className="p-3 bg-card rounded-md">
                          <p className="font-bold text-primary text-center">Reaction {r.id}</p>
                            <p className="font-medium text-center mb-1">{preview.equation}</p>
                          <p className="font-mono text-sm text-center text-muted-foreground w-full">
                            Rate = k
                            <sub>{r.id}{r.isEquilibrium ? ',f' : ''}</sub>
                            {preview.rateLaw}
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

  const componentsGridCols = `40px 1fr 1fr 1fr repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`;

  return (
    <div className={`min-h-screen flex flex-col p-4 ${mainBg} ${mainFg}`}>
      <div className="container mx-auto space-y-6">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
            {/* Reactions & Kinetics Card */}
            <div className={`p-6 rounded-lg shadow-lg ${cardBg} ${cardFg}`}> 
                <div className="mb-6">
                  <CardTitle className="flex items-center justify-between">
                    Reactions & Kinetics
                    <Button variant="outline" size="sm" onClick={addReactionSetup} disabled={reactionsSetup.length >= 6}>
                        <PlusCircle className="mr-2 h-4 w-4" />
                        Add Reaction
                    </Button>
                  </CardTitle>
                </div>
                <div className="space-y-6">
                    {reactionsSetup.map((reaction, index) => (
                    <div key={reaction.id} className={`space-y-4 ${index > 0 ? 'border-t pt-6' : ''}`}> 
                        <div className="flex items-center justify-between">
                            <Label className="font-bold">Reaction {reaction.id}</Label>
                            <div className='flex items-center gap-4'>
                            <div className='flex items-center gap-2'>
                                <Label htmlFor={`reversible-switch-${reaction.id}`} className="text-sm">Reversible</Label>
                                <Switch
                                id={`reversible-switch-${reaction.id}`}
                                checked={reaction.isEquilibrium}
                                onCheckedChange={(checked) => updateReactionSetup(reaction.id, 'isEquilibrium', checked)}
                                />
                            </div>
                            {index > 0 ? (
                                <Button variant="destructive" size="sm" onClick={() => removeReactionSetup(reaction.id)} aria-label="Remove reaction">
                                <Trash2 className="h-4 w-4" />
                                </Button>
                            ) : (
                                <Button variant="outline" size="sm" disabled aria-label="Cannot remove first reaction">
                                <Trash2 className="h-4 w-4" />
                                </Button>
                            )}
                            </div>
                        </div>
                        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                            <div className="flex items-center gap-2">
                            <Label className="w-16 text-sm" htmlFor={`a-value-${reaction.id}`}> 
                                <span>A<sub className="font-medium relative -left-px">{reaction.id}{reaction.isEquilibrium ? ',f' : ''}</sub>:</span>
                            </Label>
                            <Input id={`a-value-${reaction.id}`} type="text" value={reaction.AValue} onChange={(e) => updateReactionSetup(reaction.id, 'AValue', e.target.value)} />
                            </div>
                            <div className="flex items-center gap-2">
                            <Label className="w-16 text-sm" htmlFor={`ea-value-${reaction.id}`}> 
                                <span>Ea<sub className="font-medium relative -left-px">{reaction.id}{reaction.isEquilibrium ? ',f' : ''}</sub>:</span>
                            </Label>
                            <Input id={`ea-value-${reaction.id}`} type="text" value={reaction.EaValue} onChange={(e) => updateReactionSetup(reaction.id, 'EaValue', e.target.value)} />
                            <span className="text-xs text-muted-foreground">J/mol</span>
                            </div>
                        </div>
                        {/* Conditionally show reverse reaction inputs */}
                        {reaction.isEquilibrium && (
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-4 pt-4 border-t border-dashed">
                                <div className="flex items-center gap-2">
                                <Label className="w-16 text-sm" htmlFor={`a-value-rev-${reaction.id}`}> 
                                    <span>A<sub className="font-medium relative -left-px">{reaction.id},r</sub>:</span>
                                </Label>
                                <Input id={`a-value-rev-${reaction.id}`} type="text" value={reaction.AValueBackward || ''} placeholder='e.g. 1e12' onChange={(e) => updateReactionSetup(reaction.id, 'AValueBackward', e.target.value)} />
                                </div>
                                <div className="flex items-center gap-2">
                                <Label className="w-16 text-sm" htmlFor={`ea-value-rev-${reaction.id}`}> 
                                    <span>Ea<sub className="font-medium relative -left-px">{reaction.id},r</sub>:</span>
                                </Label>
                                <Input id={`ea-value-rev-${reaction.id}`} type="text" value={reaction.EaValueBackward || ''} placeholder='e.g. 80000' onChange={(e) => updateReactionSetup(reaction.id, 'EaValueBackward', e.target.value)} />
                                <span className="text-xs text-muted-foreground">J/mol</span>
                                </div>
                            </div>
                        )}
                    </div>
                    ))}
                </div>
            </div>
            {/* Simulation Basis Card */}
            <div className={`p-6 rounded-lg shadow-lg ${cardBg} ${cardFg}`}>
                <div className="mb-6">
                    <CardTitle>Simulation Basis</CardTitle>
                </div>
                <div className="space-y-4">
                    {/* The Reaction Phase toggle is now at the top without a title */}
                    <div className="flex items-center gap-1 rounded-lg p-1 bg-muted">
                        <Button
                            onClick={() => setReactionPhase('Liquid')}
                            variant={reactionPhase === 'Liquid' ? 'default' : 'ghost'}
                            size="sm"
                            className="flex-1 text-xs px-3 py-1 h-auto"
                        >
                            Liquid
                        </Button>
                        <Button
                            onClick={() => setReactionPhase('Gas')}
                            variant={reactionPhase === 'Gas' ? 'default' : 'ghost'}
                            size="sm"
                            className="flex-1 text-xs px-3 py-1 h-auto"
                        >
                            Gas
                        </Button>
                    </div>

                    {/* The rest of the controls remain the same */}
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
                    <div>
                        <Label>Desired Production Rate (kta)</Label>
                        <Input
                            type="number"
                            value={prodRate}
                            onChange={e => setProdRate(e.target.value)}
                            className="w-full mt-1"
                        />
                    </div>
                </div>
            </div>
            
            {/* Parsed Reactions Card */}
            {ReactionPreview}
        </div>

        {/* Components Card */}
        <div className={`p-6 rounded-lg shadow-lg ${cardBg} ${cardFg}`}> 
              <div className="mb-4">
                <CardTitle className="flex items-center justify-between">
                  Components
                  <Button variant="secondary" size="sm" onClick={addComponentSetup}>
                      <PlusCircle className="mr-2 h-4 w-4" />Add Component
                  </Button>
                </CardTitle>
              </div>
              <div className="overflow-x-auto">
                  {(() => {
                    const componentsGridCols = `40px 1fr 1fr 1fr repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`;
                    return (
                        <>
                            {/* Header Row 1: Main Titles */}
                            <div className="grid items-center text-xs font-medium text-muted-foreground" style={{gridTemplateColumns: componentsGridCols}}>
                                <div/> {/* Spacer */}
                                <div className="text-center">Name</div>
                                <div className="text-center">Molar Mass</div>
                                <div className="text-center">Density</div>
                                {reactionsSetup.map(r => (
                                    <div key={r.id} className="text-center col-span-2">{`Reaction ${r.id}`}</div>
                                ))}
                            </div>
                            {/* Header Row 2: Units and Sub-headers */}
                            <div className="grid items-center text-xs font-medium text-muted-foreground pb-2 border-b" style={{gridTemplateColumns: componentsGridCols}}>
                                <div/> {/* Spacer */}
                                <div/> {/* Spacer for Name */}
                                <div className="text-center">(g/mol)</div>
                                <div className="text-center">(g/L)</div>
                                {reactionsSetup.map(r => (
                                    <React.Fragment key={r.id}>
                                        <div className="text-center">Stoich.</div>
                                        <div className="text-center">
                                            {r.isEquilibrium ? 'Order (fwd/rev)' : 'Order'}
                                        </div>
                                    </React.Fragment>
                                ))}
                            </div>

                            {componentsSetup.map((comp, index) => (
                                <div key={comp.id} className="grid gap-2 items-center py-2" style={{gridTemplateColumns: componentsGridCols}}>
                                    {index > 1 ? (
                                        <Button variant="ghost" size="icon" className="h-8 w-8" onClick={() => removeComponentSetup(comp.id)}>
                                            <Trash2 className="h-4 w-4 text-red-500" />
                                        </Button>
                                    ) : (
                                        <Button variant="ghost" size="icon" className="h-8 w-8" disabled>
                                        <Trash2 className="h-4 w-4 text-muted-foreground" />
                                        </Button>
                                    )}

                                    {/* Correctly ordered input fields */}
                                    <Input value={comp.name} onChange={e => handleComponentSetupChange(comp.id, 'name', e.target.value)} />
                                    <Input type="number" placeholder="g/mol" value={comp.molarMass || ''} onChange={e => handleComponentSetupChange(comp.id, 'molarMass', e.target.value)} className="h-8 text-center"/>
                                    <Input type="number" placeholder="g/L" value={comp.density || ''} onChange={e => handleComponentSetupChange(comp.id, 'density', e.target.value)} className="h-8 text-center"/>
                                    {reactionsSetup.map(r => (
                                        <React.Fragment key={r.id}>
                                            <Input type="number" placeholder="0" value={comp.reactionData[r.id]?.stoichiometry || ''} onChange={e => handleComponentSetupChange(comp.id, 'reactionData', {...comp.reactionData, [r.id]: {...comp.reactionData[r.id], stoichiometry: e.target.value}})} className="h-8 text-center" step="0.1" />
                                            <div className="flex items-center gap-1">
                                                <Input 
                                                type="number" 
                                                placeholder={r.isEquilibrium ? "fwd" : "0"}
                                                value={comp.reactionData[r.id]?.order || ''} 
                                                onChange={e => handleComponentSetupChange(comp.id, 'reactionData', {...comp.reactionData, [r.id]: {...comp.reactionData[r.id], order: e.target.value}})} 
                                                className="h-8 text-center" 
                                                step="0.1" 
                                                />
                                                {r.isEquilibrium && (
                                                <Input 
                                                    type="number" 
                                                    placeholder="rev"
                                                    value={comp.reactionData[r.id]?.orderReverse || ''} 
                                                    onChange={e => handleComponentSetupChange(comp.id, 'reactionData', {...comp.reactionData, [r.id]: {...comp.reactionData[r.id], orderReverse: e.target.value}})} 
                                                    className="h-8 text-center" 
                                                    step="0.1" 
                                                />
                                                )}
                                            </div>
                                        </React.Fragment>
                                    ))}
                                </div>
                            ))}
                        </>
                    )
                  })()}
              </div>
          </div>

        <div className="text-right mt-2">
            <div title={!validationResult.isValid ? validationResult.message : undefined} className="inline-block">
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

const ReactorSimulator = ({ 
    onBack, 
    reactions, components, 
    simBasis, 
    molarRatios, setMolarRatios,
    reactionPhase, // Receive reactionPhase as a prop
    prodRate // Receive prodRate as a prop
}: { 
    onBack: () => void, 
    reactions: ReactionSetup[], 
    components: ComponentSetup[],
    simBasis: any,
    molarRatios: MolarRatio[],
    setMolarRatios: React.Dispatch<React.SetStateAction<MolarRatio[]>>,
    reactionPhase: 'Liquid' | 'Gas',
    prodRate: string 
}) => {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';
  const textColor = isDark ? '#E5E7EB' : '#1F2937';
  const cardBg = 'bg-card';
  const cardFg = 'text-card-foreground';
  const mainBg = 'bg-background';
  const mainFg = 'text-foreground';

  const [reactorTypes, setReactorTypes] = useState({ pfr: true, cstr: false })

  // Add this handler function to manage toggling
  const handleReactorTypeChange = (type: 'pfr' | 'cstr') => {
      setReactorTypes(prev => {
          const newState = { ...prev, [type]: !prev[type] }
          
          // This logic prevents deselecting both buttons
          if (!newState.pfr && !newState.cstr) {
              return prev // Return the previous state if the user tries to deselect the last active button
          }
          
          return newState
      })
  }
  const [graphType, setGraphType] = useState<GraphType>('selectivity')
  const [temperature, setTemperature] = useState(4)
  
  const [pressure, setPressure] = useState(5); // Initial pressure in bar
  const [pressureMin, setPressureMin] = useState('1');
  const [pressureMax, setPressureMax] = useState('50');
  
  const [molarRatioMin, setMolarRatioMin] = useState('2');
  const [molarRatioMax, setMolarRatioMax] = useState('20');
  const [tempMin, setTempMin] = useState('4');
  const [tempMax, setTempMax] = useState('20');
  // Removed local prodRate state

  const generateGraphData = useCallback(() => {
    const tempK = temperature + 273.15
    
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId)
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId)

    if (!limitingReactant || !desiredProduct) {
        return { series: [], xAxis: '', yAxis: '', legend: [], yMax: 'dataMax' }
    }

    const xLabel = `Conversion of ${limitingReactant?.name || 'Limiting Reactant'}`
    let yLabel = ''
    const legend: string[] = []
    const series: any[] = []
    
    const colors = {
        pfr: '#22c55e', // green
        cstr: '#E53E3E'  // red
    }

    let overallReasonableMax = -1 // Used to find a good axis limit across all visible curves

    // Loop through the selected reactor types
    for (const type in reactorTypes) {
        if (reactorTypes[type as keyof typeof reactorTypes]) {
            let dataToShow: {
                selectivityData: { x: number; y: number }[];
                volumeData: { x: number; y: number }[];
            } = { selectivityData: [], volumeData: [] }

            const reactorType = type.toUpperCase() as 'PFR' | 'CSTR'

            // --- Calculation for each type (logic is unchanged) ---
            if (reactorType === 'PFR') {
                if (reactionPhase === 'Gas') {
                    dataToShow = calculatePfrGasPhase(
                        reactions, components, tempK, pressure, molarRatios, simBasis, parseFloat(prodRate)
                    )
                } else {
                    dataToShow = calculatePfrLiquidPhase(
                        reactions, components, tempK, molarRatios, simBasis, parseFloat(prodRate)
                    )
                }
            } else { // CSTR
                const F_A0_guess = 10.0
                const initialFlowRates = { [limitingReactant.name]: F_A0_guess }
                molarRatios.forEach(ratio => {
                    const compName = components.find(c => c.id === ratio.numeratorId)?.name
                    if (compName) initialFlowRates[compName] = F_A0_guess * ratio.value
                })
                components.forEach(c => {
                    if (initialFlowRates[c.name] === undefined) initialFlowRates[c.name] = 0
                })
                
                dataToShow = calculateCstrData(
                    reactions, components, initialFlowRates, tempK, simBasis, 
                    reactionPhase,
                    pressure,
                    parseFloat(prodRate)
                )
            }

            // --- **FIX START**: Determine a reasonable Y-axis max for volume graphs ---
            if (graphType === 'volume') {
                const volumeData = dataToShow.volumeData;
                const reasonablePoints = volumeData.filter(p => p.x < 0.95); // Find all points below 98% conversion
                if (reasonablePoints.length > 0) {
                    const lastReasonablePoint = reasonablePoints[reasonablePoints.length - 1];
                    // If this curve's reasonable max is higher, it becomes the new overall max
                    if (lastReasonablePoint.y > overallReasonableMax) {
                        overallReasonableMax = lastReasonablePoint.y;
                    }
                }
            }
            // --- **FIX END** ---
// This formatter function creates the custom tooltip content
      formatter: (params: any) => {
        if (!params || params.length === 0) {
            return ''
        }

        const xAxisLabel = graphData.xAxis.split(' ')[0]
        const xAxisValue = parseFloat(parseFloat(params[0].axisValue).toPrecision(6))
        let tooltipContent = `<strong>${xAxisLabel}</strong>: ${xAxisValue}`

        params.forEach((p: any) => {
            const seriesName = p.seriesName
            const rawValue = parseFloat(p.value[1])
            let formattedValue = ''

            if (graphType === 'selectivity') {
                formattedValue = rawValue.toPrecision(3)
            } else { // 'volume'
                if (rawValue === 0) {
                    formattedValue = '0';
                } else if (rawValue < 1) {
                    formattedValue = rawValue.toPrecision(3);
                } else {
                    const magnitude = Math.floor(Math.log10(rawValue));
                    const scale = Math.pow(10, magnitude - 2);
                    formattedValue = (Math.round(rawValue / scale) * scale).toLocaleString('en-US', {useGrouping: true, maximumFractionDigits: 20});
                }
                // --- MODIFICATION: Add units to the volume string ---
                formattedValue += ' m³';
            }
            
            const seriesColor = p.color

            tooltipContent += `<br/>`
            tooltipContent += `<span style="color: ${seriesColor}; font-weight: bold;">${seriesName}:</span> <strong>${formattedValue}</strong>`
        })

        return tooltipContent
      }
            // --- MODIFICATION: Check how many reactors are active ---
            const activeReactorCount = (reactorTypes.pfr ? 1 : 0) + (reactorTypes.cstr ? 1 : 0);

            // --- MODIFICATION: Set the series name conditionally ---
            const baseName = graphType === 'volume' ? 'Reactor Volume' : 'Selectivity';
            const seriesName = activeReactorCount > 1
                ? `${reactorType.toUpperCase()} - ${baseName}`
                : baseName;
            legend.push(seriesName)
            
            if (graphType === 'volume') {
                yLabel = 'Reactor Volume (m³)'

                const finalVolumeData = dataToShow.volumeData
                                            .map(d => [d.x, d.y]); // already in m³

                series.push({ 
                    name: seriesName, 
                    type: 'line', 
                    data: finalVolumeData, 
                    smooth: true, 
                    showSymbol: false, 
                    color: colors[type as keyof typeof colors],
                    lineStyle: { width: 4 } 
                })
            } else { 
                yLabel = `Selectivity to ${desiredProduct?.name || 'Product'}`
                series.push({ 
                    name: seriesName, 
                    type: 'line', 
                    data: dataToShow.selectivityData.map(d => [d.x, d.y]), 
                    smooth: true, 
                    showSymbol: false, 
                    color: colors[type as keyof typeof colors],
                    lineStyle: { width: 4 }
                })
            }
        }
    }
    
    // Set the final Y-axis max value
    let yAxisMax: number | 'dataMax' = 'dataMax';
    if (graphType === 'volume' && overallReasonableMax > 0) {
        // --- MODIFICATION: Set max directly. ECharts will handle the scaling. ---
        yAxisMax = overallReasonableMax;
    }

    return { series, xAxis: xLabel, yAxis: yLabel, legend, yMax: yAxisMax }

}, [reactorTypes, reactionPhase, pressure, molarRatios, temperature, graphType, reactions, components, simBasis, prodRate])

  const graphData = generateGraphData();
    // --- MODIFICATION: Set chart title based on selection ---
    let chartTitle = 'Reactor Performance'; // Default title
    if (reactorTypes.pfr && !reactorTypes.cstr) {
        chartTitle = 'PFR Performance';
    } else if (!reactorTypes.pfr && reactorTypes.cstr) {
        chartTitle = 'CSTR Performance';
    }

    const chartOptions: EChartsOption = {
        animation: false,
        backgroundColor: 'transparent',
        title: {
                text: chartTitle, // Use the new dynamic title
                left: 'center',
                top: '2%',
                textStyle: { color: textColor, fontSize: 18, fontFamily: 'Merriweather Sans' }
        },
    grid: { left: '5%', right: '3%', bottom: '8%', top: '10%', containLabel: true },
    xAxis: {
        type: 'value',
        name: graphData.xAxis,
    // --- MODIFICATION: Add the scale property ---
    scale: true,
        nameLocation: 'middle', nameGap: 30,
        nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
        min: 0, max: 1,
        axisLine: { lineStyle: { color: textColor } },
        axisLabel: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
        splitLine: { show: false },
    },
    yAxis: {
        type: 'value',
        name: graphData.yAxis,
        nameLocation: 'middle',
        nameGap: 50,
        nameTextStyle: { color: textColor, fontSize: 14, fontFamily: 'Merriweather Sans' },
        min: 0,
        
        // **FIX START**: Use the dynamic max value for volume graphs
        max: graphType === 'selectivity' ? 1 : graphData.yMax,
        // **FIX END**

        axisTick: {
            show: true
        },

        axisLine: { lineStyle: { color: textColor } },
        axisLabel: {
            color: textColor,
            fontSize: 14,
            fontFamily: 'Merriweather Sans',
            // This formatter hides the label for the maximum tick
            formatter: (value: number) => {
                if (
                    graphType === 'volume' &&
                    typeof graphData.yMax === 'number' &&
                    value >= graphData.yMax
                ) {
                    return '';
                }
                // Formatting for all other labels remains the same
                if (graphType === 'volume') {
                    if (value === 0) return '0';
                    if (value < 1) return value.toPrecision(3);
                    const magnitude = Math.floor(Math.log10(value));
                    const scale = Math.pow(10, magnitude - 2);
                    return (Math.round(value / scale) * scale).toString();
                }
                return value.toString();
            }
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
      
      axisPointer: {
        type: 'cross',
        label: {
          backgroundColor: isDark ? '#374151' : '#e5e7eb',
          color: textColor,
          fontFamily: 'Merriweather Sans',
          fontSize: 12,
          formatter: function (params) {
            return parseFloat(params.value as string).toPrecision(3)
          }
        }
      },

      // This formatter function creates the custom tooltip content 
       formatter: (params: any) => { 
         if (!params || params.length === 0) { 
             return '' 
         } 
 
         const xAxisLabel = graphData.xAxis.split(' ')[0] 
         const xAxisValue = parseFloat(parseFloat(params[0].axisValue).toPrecision(6)) 
         let tooltipContent = `<strong>${xAxisLabel}</strong>: ${xAxisValue}` 
 
         params.forEach((p: any) => { 
             const seriesName = p.seriesName // This name is now updated 
             const rawValue = parseFloat(p.value[1]) 
             let formattedValue = '' 
 
             // --- MODIFICATION: Conditional Value Formatting --- 
             if (graphType === 'selectivity') { 
                 // Format selectivity to 3 significant figures 
                 formattedValue = rawValue.toPrecision(3) 
             } else { // 'volume' 
                 if (rawValue === 0) { 
                     formattedValue = '0'; 
                 } else if (rawValue < 1) { 
                     // For values < 1, use scientific notation with 3 sig figs 
                     formattedValue = rawValue.toPrecision(3); 
                 } else { 
                     // For values >= 1, round to 3 sig figs without scientific notation 
                     const magnitude = Math.floor(Math.log10(rawValue)); 
                     const scale = Math.pow(10, magnitude - 2); // Adjust scale for 3 sig figs 
                     formattedValue = (Math.round(rawValue / scale) * scale).toLocaleString('en-US', {useGrouping: true, maximumFractionDigits: 20}); 
                 } 
                 // --- MODIFICATION: Add units to the volume string --- 
                 formattedValue += ' m³'; 
             } 
             // --- END MODIFICATION ---
             const seriesColor = p.color 
 
            tooltipContent += `<br/>` 
             tooltipContent += `<span style="color: ${seriesColor}; font-weight: bold;">${seriesName}:</span> <strong>${formattedValue}</strong>` 
         }) 
 
         return tooltipContent
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
      <div className="container mx-auto flex-grow">
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6 h-full">
            {/* Control Panel Card */}
            <div className={`lg:col-span-1 p-6 rounded-lg shadow-lg flex flex-col space-y-6 ${cardBg} ${cardFg}`}>
            <div className="flex justify-start">
                 <button onClick={onBack} className="bg-slate-700 hover:bg-slate-600 text-white font-bold py-2 px-4 rounded-lg transition-colors flex items-center">
                    <ArrowLeft className="h-4 w-4 mr-2"/> Back
                </button>
            </div>
                
            {molarRatios.map(ratio => {
                const numeratorName = components.find(c => c.id === ratio.numeratorId)?.name || 'N/A';
                return (
                    <div key={ratio.numeratorId} className="pt-4 border-t border-border">
                        <ParameterSlider
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



            {/* Pressure Slider (only for Gas phase) */}
            {reactionPhase === 'Gas' && (
                <ParameterSlider
                    label="Pressure"
                    unit="bar"
                    value={pressure.toString()}
                    onValueChange={(value) => setPressure(value)}
                    min={parseFloat(pressureMin)}
                    max={parseFloat(pressureMax)}
                    step={0.5}
                    maxSliderValue={pressureMax}
                    onMaxSliderChange={setPressureMax}
                    minSliderValue={pressureMin}
                    onMinSliderChange={setPressureMin}
                />
            )}
            </div>

            {/* Graph Card */}
            <div className={`lg:col-span-2 p-4 rounded-lg shadow-lg aspect-square ${cardBg} ${cardFg}`}> 
              {/* Chart section is unchanged */}
            <div className="relative mb-2">
                <div className="absolute top-0 left-0 z-10">
                        {/* Wrap the buttons in a muted background for a nice toggle effect */}
                        <div className="flex items-center gap-1 rounded-lg p-1 bg-muted">
                            <Button
                                onClick={() => handleReactorTypeChange('pfr')}
                                variant={reactorTypes.pfr ? 'default' : 'ghost'}
                                size="sm"
                                className="text-xs px-3 py-1 h-auto"
                            >
                                PFR
                            </Button>
                            <Button
                                onClick={() => handleReactorTypeChange('cstr')}
                                variant={reactorTypes.cstr ? 'default' : 'ghost'}
                                size="sm"
                                className="text-xs px-3 py-1 h-auto"
                            >
                                CSTR
                            </Button>
                    </div>
                </div>
                <div className="absolute top-0 right-0 z-10">
                    <Select value={graphType} onValueChange={(value) => setGraphType(value as GraphType)}>
                        <SelectTrigger className="w-[200px] h-8 text-xs">
                            <SelectValue placeholder="Select graph type" />
                        </SelectTrigger>
                        <SelectContent>
                            <SelectItem value="selectivity">Selectivity vs. Conversion</SelectItem>
                            <SelectItem value="volume">Volume vs. Conversion</SelectItem>
                        </SelectContent>
                    </Select>
                </div>
            </div>
          <ReactECharts
            key={`${resolvedTheme}-${reactorTypes.pfr ? 'pfr' : 'cstr'}`}
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
          <Input id={`${label}-min-input`} type="number" value={minSliderValue} onChange={(e) => onMinSliderChange(e.target.value)} className="w-14 h-8 text-xs" />
          <Label htmlFor={`${label}-max-input`} className="text-xs text-muted-foreground">Max:</Label>
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

// --- FIND VOLUME HELPER FUNCTIONS (Updated) ---

const findVolumeForConversion_PFR = (
    targetConversion: number,
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number },
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string }
): number => {
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    if (!limitingReactant) return 2000;

    const errorFunction = (volumeGuess: number): number => {
        if (volumeGuess <= 0) return -targetConversion;
        const profile = calculatePfrData(reactions, components, volumeGuess, initialFlowRates, tempK, simBasis);
        let maxConversion = 0;
        if (profile.volumeData.length > 0) {
            maxConversion = Math.max(...profile.volumeData.map(p => p.x));
        }
        return maxConversion - targetConversion;
    };

    let low = 0.1;
    let high = 100.0; // Start with a reasonable guess
    const maxBracketingIter = 10; // Failsafe to prevent infinite loops

    // Step 1: Bracket the root by expanding the search range if needed.
    for (let i = 0; i < maxBracketingIter; i++) {
        if (errorFunction(high) >= 0) break; // Root is bracketed
        low = high;
        high *= 10; // Exponentially increase search volume
        if (i === maxBracketingIter - 1) return high; // Return max if bracketing fails
    }

    // Step 2: Use bisection to refine the volume once bracketed.
    const maxIter = 50;
    for (let i = 0; i < maxIter; i++) {
        const mid = (low + high) / 2;
        if (mid === low || mid === high) break; // No progress
        const error = errorFunction(mid);
        if (Math.abs(error) < 0.01) return mid * 1.05; // Return result with a 5% buffer
        if (error < 0) {
            low = mid; // Conversion too low, need more volume
        } else {
            high = mid; // Conversion too high, need less volume
        }
    }
    return high * 1.1; // Return a safe upper bound if bisection doesn't converge
};



// Add this new helper function to your file
const calculatePfrDataByTau = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialConcentrations: { [key: string]: number },
    tempK: number,
    t_max: number
): Array<{ tau: number; concentrations: { [key: string]: number } }> => {
    const componentNames = components.map(c => c.name).filter(name => initialConcentrations[name] !== undefined);

    const dCdTau = (tau: number, C_vec: number[]): number[] => {
        const C_dict = componentNames.reduce((acc, name, i) => {
            acc[name] = C_vec[i];
            return acc;
        }, {} as { [key: string]: number });

        // --- CORRECTED RATE CALCULATION ---
        const reactionNetRates = reactions.map(reactionInfo => {
            // Forward reaction rate calculation
            const k_f = parseFloat(reactionInfo.AValue) * Math.exp(-parseFloat(reactionInfo.EaValue) / (R * tempK));
            let rate_f = k_f;
            components.forEach(comp => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
                if (stoich < 0 && order > 0) { // It's a reactant
                    rate_f *= Math.pow(Math.max(0, C_dict[comp.name]), order);
                }
            });

            if (!reactionInfo.isEquilibrium) {
                return rate_f; // Return only forward rate if not reversible
            }

            // Reverse reaction rate calculation
            const k_r = parseFloat(reactionInfo.AValueBackward || '0') * Math.exp(-parseFloat(reactionInfo.EaValueBackward || '0') / (R * tempK));
            let rate_r = k_r;
            components.forEach(comp => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                const reverseOrder = parseFloat(comp.reactionData[reactionInfo.id]?.orderReverse || '0');
                if (stoich > 0 && reverseOrder > 0) { // It's a product
                    rate_r *= Math.pow(Math.max(0, C_dict[comp.name]), reverseOrder);
                }
            });

            return rate_f - rate_r; // Return the net rate
        });
        // --- END CORRECTION ---

        return componentNames.map(name => {
            let netRateOfFormation = 0; // r_i
            reactions.forEach((reactionInfo, j) => {
                const stoich = parseFloat(components.find(c => c.name === name)?.reactionData[reactionInfo.id]?.stoichiometry || '0');
                netRateOfFormation += stoich * reactionNetRates[j];
            });
            return netRateOfFormation; // dC/dτ = r_i
        });
    };
    
    const C0 = componentNames.map(name => initialConcentrations[name] || 0);
    const odeResults = solveODE_BDF(dCdTau, C0, [0, t_max]);
    
    const results = [];
    if (odeResults && odeResults.t) {
        for (let i = 0; i < odeResults.t.length; i++) {
            const tau = odeResults.t[i];
            const concentrations = componentNames.reduce((acc, name, j) => {
                acc[name] = odeResults.y[j][i];
                return acc;
            }, {} as { [key: string]: number });
            results.push({ tau, concentrations });
        }
    }
    return results;
};

// CSTR solver by residence time
const calculateCstrDataByTau = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialConcentrations: { [key: string]: number },
    tempK: number,
    t_max: number
): Array<{ tau: number; concentrations: { [key: string]: number } }> => {
    const componentNames = components.map(c => c.name).filter(name => initialConcentrations[name] !== undefined);
    const C0_vec = componentNames.map(name => initialConcentrations[name] || 0);
    
    let solverGuess: number[] | undefined = C0_vec;
    const results = [];
    const n_steps = 100;

    // Use a logarithmic scale for tau and loop backwards (for continuation)
    const log_t_min = Math.log10(Math.max(0.01, t_max / 1e6));
    const log_t_max = Math.log10(t_max);
    
    for (let i = n_steps; i >= 0; i--) {
        const log_tau = log_t_min + (i/n_steps)*(log_t_max - log_t_min);
        const tau = Math.pow(10, log_tau);
        
        // This is a simplified CSTR root-finder. Your more complex `solveCSTRParallel` is more robust.
        // For simplicity here, we abstract the solver logic.
        const C_out_vec = solveCSTR_by_tau(tau, reactions, components, initialConcentrations, tempK, solverGuess);
        solverGuess = C_out_vec; // Use result as next guess

        const concentrations = componentNames.reduce((acc, name, j) => {
            acc[name] = C_out_vec[j];
            return acc;
        }, {} as { [key: string]: number });

        results.push({ tau, concentrations });
    }
    return results.reverse(); // Return in ascending order of tau
};

// You will also need a CSTR solver that works by tau, similar to the PFR one.
// Here is a simplified version.
const solveCSTR_by_tau = (
    tau: number,
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialConcentrations: { [key: string]: number },
    tempK: number,
    initialGuess?: number[]
): number[] => {
    const componentNames = components.map(c => c.name).filter(name => initialConcentrations[name] !== undefined);
    const C0_vec = componentNames.map(name => initialConcentrations[name] || 0);

    let C_vec = initialGuess ? [...initialGuess] : [...C0_vec];
    const n = C_vec.length;
    
    const G = (current_C_vec: number[]): number[] => {
        const C_dict = componentNames.reduce((acc, name, i) => {
            acc[name] = current_C_vec[i];
            return acc;
        }, {} as { [key: string]: number });

        // --- CORRECTED RATE CALCULATION ---
        const reactionNetRates = reactions.map(reactionInfo => {
            const k_f = parseFloat(reactionInfo.AValue) * Math.exp(-parseFloat(reactionInfo.EaValue) / (R * tempK));
            let rate_f = k_f;
            components.forEach(comp => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
                if (stoich < 0 && order > 0) {
                    rate_f *= Math.pow(Math.max(0, C_dict[comp.name]), order);
                }
            });

            if (!reactionInfo.isEquilibrium) {
                return rate_f;
            }

            const k_r = parseFloat(reactionInfo.AValueBackward || '0') * Math.exp(-parseFloat(reactionInfo.EaValueBackward || '0') / (R * tempK));
            let rate_r = k_r;
            components.forEach(comp => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                const reverseOrder = parseFloat(comp.reactionData[reactionInfo.id]?.orderReverse || '0');
                if (stoich > 0 && reverseOrder > 0) {
                    rate_r *= Math.pow(Math.max(0, C_dict[comp.name]), reverseOrder);
                }
            });

            return rate_f - rate_r;
        });
        // --- END CORRECTION ---

        return componentNames.map((name, i) => {
            let netRateOfFormation = 0; // r_i
            reactions.forEach((reactionInfo, j) => {
                const stoich = parseFloat(components.find(c => c.name === name)?.reactionData[reactionInfo.id]?.stoichiometry || '0');
                netRateOfFormation += stoich * reactionNetRates[j];
            });
            return C0_vec[i] - current_C_vec[i] + tau * netRateOfFormation;
        });
    };

    const max_iter = 50;
    const tol = 1e-9;
    for (let iter = 0; iter < max_iter; iter++) {
        const G_current = G(C_vec);
        if (Math.sqrt(G_current.reduce((s, v) => s + v * v, 0)) < tol) break;
        
        const J_G: number[][] = Array(n).fill(0).map(() => Array(n).fill(0));
        for (let j = 0; j < n; j++) {
            const C_plus_h = [...C_vec];
            const h = (Math.abs(C_vec[j]) || 1e-8) * 1e-7;
            C_plus_h[j] += h;
            const G_plus_h = G(C_plus_h);
            for (let i = 0; i < n; i++) J_G[i][j] = (G_plus_h[i] - G_current[i]) / h;
        }

        const delta_C = solveLinearSystem(J_G, G_current.map(val => -val));
        if (!delta_C) break;
        
        const dampingFactor = 0.8;
        C_vec = C_vec.map((val, i) => val + dampingFactor * delta_C[i]);
        C_vec = C_vec.map(val => Math.max(0, val));
    }

    return C_vec;
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
  
  // State for production rate now lives in the parent component
  const [prodRate, setProdRate] = useState('100');

  // Add reactionPhase to your state declarations
  const [reactionPhase, setReactionPhase] = useState<'Liquid' | 'Gas'>('Liquid'); 

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
              molarRatios={molarRatios}
              setMolarRatios={setMolarRatios}
              reactionPhase={reactionPhase} // Pass the phase state down
              prodRate={prodRate}
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
              reactionPhase={reactionPhase}      // Pass state down
              setReactionPhase={setReactionPhase} // Pass setter down
              prodRate={prodRate} // Pass state and setter down
              setProdRate={setProdRate}
            />
  }
}