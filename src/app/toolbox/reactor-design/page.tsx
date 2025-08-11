
"use client"

let concentrations: { [key: string]: number } = {};

import React, { useState, useMemo, useEffect, useCallback } from 'react'
import ReactECharts from 'echarts-for-react'
import { EChartsOption } from 'echarts'
import { PlusCircle, Trash2, ArrowLeft, ArrowRight } from 'lucide-react'
import { Label } from '@/components/ui/label'
import { Input } from '@/components/ui/input'
import { Slider } from '@/components/ui/slider'
import { useTheme } from 'next-themes';
import { Button } from '@/components/ui/button';
import { CardTitle } from '@/components/ui/card';
import { Card, CardHeader, CardContent } from '@/components/ui/card';
import { Select } from '@/components/ui/select';
import { SelectItem } from '@/components/ui/select';
import { SelectTrigger } from '@/components/ui/select';
import { SelectValue } from '@/components/ui/select';
import { SelectContent } from '@/components/ui/select';
import { Switch } from '@/components/ui/switch';

// ==================================================================
// ▼▼▼ PLACE THE UNIT CONVERSION UTILITY HERE ▼▼▼
// ==================================================================

const units = {
    YEAR_IN_SECONDS: 365 * 24 * 3600,
    PASCALS_PER_BAR: 1e5,
    CUBIC_METERS_PER_LITER: 0.001,
    GRAMS_PER_KILOTONNE: 1e9,
    KELVIN_OFFSET: 273.15,
}

const convertUnit = {
    kinetics: {
        kcToKp: (kc: number, tempK: number, totalOrder: number): number => {
            // Converts a rate constant from a concentration basis (kc) to a pressure basis (kp).
            // kc is based on mol/L, kp is based on Pascals.
            // The conversion is: kp = kc * (R*T)^(-n) where n is the total reactant order.
            // R must be in (L·Pa)/(mol·K) to match the units.
            // R = 8.314 m³·Pa/mol·K * 1000 L/m³ = 8314 L·Pa/mol·K.
            const R_LITERS_PASCALS = 8314;
            return kc * Math.pow(R_LITERS_PASCALS * tempK, -totalOrder);
        }
    },
    pressure: {
        barToPa: (bar: number): number => bar * units.PASCALS_PER_BAR,
    },
    volume: {
        litersToCubicMeters: (liters: number): number => liters * units.CUBIC_METERS_PER_LITER,
    },
    concentration: {
        molPerCubicMeterToMolPerLiter: (molPerM3: number): number => molPerM3 * units.CUBIC_METERS_PER_LITER,
    },
    temperature: {
        celsiusToKelvin: (celsius: number): number => celsius + units.KELVIN_OFFSET,
    },
    flowRate: {
        ktaToGramsPerSecond: (kta: number): number => (kta * units.GRAMS_PER_KILOTONNE) / units.YEAR_IN_SECONDS,
        ktaToMolPerSecond: (kta: number, molarMass_g_per_mol: number): number => {
            if (molarMass_g_per_mol <= 0) return 0
            const gramsPerSecond = convertUnit.flowRate.ktaToGramsPerSecond(kta)
            return gramsPerSecond / molarMass_g_per_mol
        }
    },
    energy: {
        convertEaToJules: (Ea: number, unit: string): number => {
            switch (unit) {
                case 'kcal/mol': return Ea * 4184;
                case 'cal/mol': return Ea * 4.184;
                case 'kJ/mol': return Ea * 1000;
                case 'J/mol':
                default:
                    return Ea;
            }
        }
    }
}

// ==================================================================
// ▲▲▲ END OF UTILITY SECTION ▲▲▲
// ==================================================================

// --- Type Definitions ---
interface ReactionSetup {
    id: string
    AValue: string
    AUnit: string // New
    EaValue: string
    EaUnit: string // New
    isEquilibrium: boolean
    AValueBackward?: string
    AUnitBackward?: string // New
    EaValueBackward?: string
    EaUnitBackward?: string // New
    rateLawBasis?: 'concentration' | 'partialPressure'
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
    // Butane Alkylation (Liquid phase -> simple units)
    { id: '1', AValue: '3.66e14', AUnit: 'mol,L,s', EaValue: '101600', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
    { id: '2', AValue: '4.77e16', AUnit: 'mol,L,s', EaValue: '110850', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
]

const initialComponents: ComponentSetup[] = [
  { id: 'comp-A', name: '1-Butene', molarMass: '56.11', density: '595', reactionData: { '1': { stoichiometry: '-1', order: '1' }, '2': { stoichiometry: '-1', order: '1' } } },
  { id: 'comp-B', name: 'Isobutane', molarMass: '58.12', density: '593', reactionData: { '1': { stoichiometry: '-1', order: '1' }, '2': { stoichiometry: '0', order: '0' } } },
  { id: 'comp-C', name: 'Isooctane', molarMass: '114.23', density: '692', reactionData: { '1': { stoichiometry: '1', order: '0' }, '2': { stoichiometry: '-1', order: '1' } } },
  { id: 'comp-D', name: 'Dodecane', molarMass: '170.34', density: '750', reactionData: { '1': { stoichiometry: '0', order: '0' }, '2': { stoichiometry: '1', order: '0' } } },
]


// --- Calculation Engine (Now more generic) ---
const R = 8.314

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

const calculateCstrData = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number }, // Molar flows based on F_A0=1.0
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    reactionPhase: 'Liquid' | 'Gas',
    pressure_bar: number,
    targetProduction_kta: number
) => {
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);

    if (!limitingReactant || !desiredProduct) return { selectivityData: [], volumeData: [] };

    // --- START: ADD THIS SECTION ---
    const primaryReactionId = reactions[0].id;
    const stoichCoeffStr = limitingReactant.reactionData[primaryReactionId]?.stoichiometry;
    const stoichFactor = Math.abs(parseFloat(stoichCoeffStr || '0'));
    // --- END: ADD THIS SECTION ---

    const rawSelectivityData: { x: number; y: number }[] = [];
    const rawVolumeData: { x: number; y: number }[] = [];

    // --- Prepare inputs for the robust solver ---
    const allInvolvedComponents = components.filter(c => 
        Object.values(c.reactionData).some(d => d.stoichiometry && parseFloat(d.stoichiometry) !== 0)
    );
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
            return { 
                reactants, 
                products, 
                rateConstantAtT: A * Math.exp(-Ea / (R * tempK)),
                rateLawBasis: reactionInfo.rateLawBasis // Pass basis to solver
            };
        }),
        allInvolvedComponents: allInvolvedComponents,
    };

    // --- Basis Calculation Setup ---
    const F_A0_basis = 1.0;
    const F_in_basis: { [key: string]: number } = { [limitingReactant.name]: F_A0_basis };
    const F_A_in_basis_for_ratio = initialFlowRates[limitingReactant.name] || 1.0;
    for (const comp of components) {
        if (comp.id !== limitingReactant.id) {
            const ratio = (initialFlowRates[comp.name] || 0) / F_A_in_basis_for_ratio;
            F_in_basis[comp.name] = F_A0_basis * ratio;
        }
    }
    components.forEach(c => { if (F_in_basis[c.name] === undefined) F_in_basis[c.name] = 0; });
    
    // --- Iterate across a logarithmic range of basis volumes ---
    const logV_min = -8; // from 1e-8 Liters
    const logV_max = 12;  // to 1e12 Liters
    // MODIFICATION: Increased steps from 100 to 500 for a smoother curve
    const n_steps = 500; 
    let last_F_out_basis: number[] | undefined = undefined;

    for (let i = 0; i <= n_steps; i++) {
        const logV = logV_min + i * (logV_max - logV_min) / n_steps;
        const V_basis_liters = Math.pow(10, logV);
        
        // 1. Solve the CSTR equations for the current basis volume
        const F_out_basis = solveCSTR_Robust(
            parsedReactions,
            components,
            V_basis_liters,
            F_in_basis,
            reactionPhase,
            pressure_bar,
            tempK,
            last_F_out_basis // Pass the last result as the initial guess
        );

        // After successfully getting F_out_basis, update the guess for the next iteration
        const F_out_vec = parsedReactions.allInvolvedComponents.map((c: any) => F_out_basis[c.name] || 0);
        if (F_out_vec.some(f => !isNaN(f))) {
            last_F_out_basis = F_out_vec;
        }

        // 2. Calculate conversion and selectivity from the solver's results
    const F_A_out = F_out_basis[limitingReactant.name] || 0;
    const F_P_out = F_out_basis[desiredProduct.name] || 0;

    const conversion = (F_A0_basis - F_A_out) / F_A0_basis;
    const molesConsumed_A = F_A0_basis - F_A_out;
    const molesFormed_P = F_P_out - (F_in_basis[desiredProduct.name] || 0);

    // --- START: MODIFY THIS LINE ---
    const selectivity = molesConsumed_A > 1e-9 ? (stoichFactor * molesFormed_P) / molesConsumed_A : 0;
    // --- END: MODIFY THIS LINE ---
        
        if (conversion < 0.001 || conversion > 0.999 || selectivity <= 0) continue;

        // 3. Scale the volume to meet the true production target
        const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
        const P_target_mols = convertUnit.flowRate.ktaToMolPerSecond(targetProduction_kta, productMolarMass);
        
        const F_A0_true = P_target_mols / (conversion * selectivity);
        const V_true_liters = V_basis_liters * (F_A0_true / F_A0_basis);
        const V_true_m3 = convertUnit.volume.litersToCubicMeters(V_true_liters);

        if (isFinite(V_true_m3)) {
            rawVolumeData.push({ x: conversion, y: V_true_m3 });
            rawSelectivityData.push({ x: conversion, y: selectivity });
        }
    }

    // 4. Interpolate results for a smooth curve (this part is unchanged)
    const commonConversionGrid = Array.from({ length: 100 }, (_, i) => i * 0.01);
    const selectivityData = commonConversionGrid.map(conv => ({x: conv, y: interpolateFromData(conv, rawSelectivityData)})).filter(p => !isNaN(p.y));
    const volumeData = commonConversionGrid.map(conv => ({x: conv, y: interpolateFromData(conv, rawVolumeData)})).filter(p => !isNaN(p.y));

    return { selectivityData, volumeData };
}// --- Update calculatePfrData to return outletFlowsAtV ---
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

// Add this new, more robust solver function to your file.
const solveCSTR_Robust = (
  parsedReactions: any,
  components: ComponentSetup[],
  V: number, // Liters
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

  const G = (current_F: number[]): number[] => {
    const F_dict = componentNames.reduce((acc: { [key: string]: number }, name: string, i: number) => { acc[name] = current_F[i]; return acc; }, {} as {[key:string]:number});
    let rates: number[];

    if (reactionPhase === 'Liquid') {
        let v = 0;
        componentNames.forEach((name: string) => {
            const compInfo = components.find(c => c.name === name);
            const molarMass = parseFloat(compInfo?.molarMass || '1');
            const density = parseFloat(compInfo?.density || '1000');
            if (density > 0) v += (F_dict[name] * molarMass) / density;
        });
        if (v < 1e-9) v = 1e-9; // v is in Liters
        
        let concentrations: { [key: string]: number } = {};
        componentNames.forEach((name: string) => concentrations[name] = Math.max(0, F_dict[name]) / v);
        
        rates = parsedReactions.reactions.map((reaction: any) => {
            let rate = reaction.rateConstantAtT || 0;
            reaction.reactants.forEach((reactant: any) => {
                rate *= Math.pow(concentrations[reactant.name], reactant.reactionOrderNum || 0);
            });
            return rate; // mol/(L*s)
        });

    } else { // Gas Phase
        const F_total = current_F.reduce((sum: number, f: number) => sum + Math.max(0, f), 0);
        let concentrations_L: { [key: string]: number } = {};
        let partialPressures_Pa: { [key: string]: number } = {};

        if (F_total < 1e-9) {
            componentNames.forEach((name: string) => {
                concentrations_L[name] = 0;
                partialPressures_Pa[name] = 0;
            });
        } else {
            const P_Pa = convertUnit.pressure.barToPa(totalPressure_bar);

            // First, calculate concentration in mol/m³
            const concentrations_m3 = componentNames.reduce((acc: { [key: string]: number }, name: string) => {
                acc[name] = (Math.max(0, F_dict[name]) / F_total) * P_Pa / (R * temp_K);
                return acc;
            }, {} as { [key: string]: number });

            // Now, correctly calculate and assign the values needed for the rate law
            for (const name of componentNames) {
                concentrations_L[name] = convertUnit.concentration.molPerCubicMeterToMolPerLiter(concentrations_m3[name]);
                partialPressures_Pa[name] = concentrations_m3[name] * R * temp_K;
            }
        }
        
        // --- START: Replace the 'rates' calculation with this corrected block ---
        rates = parsedReactions.reactions.map((reaction: any) => {
            const k_c = reaction.rateConstantAtT || 0;
            let rate: number;
            const usePartialPressure = reaction.rateLawBasis === 'partialPressure';

            if (usePartialPressure) {
                let totalOrder = 0;
                reaction.reactants.forEach((r: any) => { totalOrder += r.reactionOrderNum || 0; });
                const k_p = convertUnit.kinetics.kcToKp(k_c, temp_K, totalOrder);
                rate = k_p;
            } else {
                rate = k_c;
            }

            reaction.reactants.forEach((reactant: any) => {
                if (usePartialPressure) {
                    rate *= Math.pow(partialPressures_Pa[reactant.name], reactant.reactionOrderNum || 0);
                } else {
                    rate *= Math.pow(concentrations_L[reactant.name], reactant.reactionOrderNum || 0);
                }
            });
            return rate;
        });
        // --- END: Replacement block ---
    }

    return componentNames.map((name: string, i: number) => {
      let netRateOfFormation = 0; // r_j in mol/(L*s)
      parsedReactions.reactions.forEach((reaction: any, j: number) => {
        const reactantInfo = reaction.reactants.find((r: any) => r.name === name);
        const productInfo = reaction.products.find((p: any) => p.name === name);
        if (reactantInfo) netRateOfFormation -= (reactantInfo.stoichiometricCoefficient || 0) * rates[j];
        if (productInfo) netRateOfFormation += (productInfo.stoichiometricCoefficient || 0) * rates[j];
      });
      // F0 (mol/s) - F (mol/s) + r (mol/L/s) * V (L) = mol/s
      return (initialFlowRates[name] || 0) - current_F[i] + netRateOfFormation * V;
    });
  };

  // The rest of the function (Newton's method solver) is unchanged...
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

const calculateNetReactionRates = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    concentrations_L: { [key: string]: number },
    partialPressures_Pa: { [key: string]: number },
    tempK: number
) => {
    return reactions.map(reactionInfo => {
        // --- Time Unit Conversion for A-Value ---
        let A_f = parseFloat(reactionInfo.AValue);
        if (reactionInfo.AUnit.includes('h')) {
            A_f /= 3600.0; // Convert hours to seconds
        }
        
        // --- Forward Rate Calculation ---
        const Ea_f_J = convertUnit.energy.convertEaToJules(parseFloat(reactionInfo.EaValue), reactionInfo.EaUnit);
        const k_f_c = A_f * Math.exp(-Ea_f_J / (R * tempK));
        
        const usePartialPressure = reactionInfo.rateLawBasis === 'partialPressure';
        let k_f: number;

        // MODIFICATION: Skip kcToKp conversion if basis is already partial pressure
        if (usePartialPressure) {
            k_f = k_f_c;
        } else {
            k_f = k_f_c;
        }
        
        let rate_f = k_f;
        components.forEach(comp => {
            const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
            const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
            
            if (order > 0 && stoich < 0) {
                 if (usePartialPressure) {
                    rate_f *= Math.pow(Math.max(0, partialPressures_Pa[comp.name]), order);
                 } else {
                    rate_f *= Math.pow(Math.max(0, concentrations_L[comp.name]), order);
                 }
            }
        });

        // --- Reverse Rate Calculation (only if reaction is reversible) ---
        let rate_r = 0;
        if (reactionInfo.isEquilibrium) {
            // --- Time Unit Conversion for Reverse A-Value ---
            let A_r = parseFloat(reactionInfo.AValueBackward || '0');
            if (reactionInfo.AUnitBackward?.includes('h')) {
                A_r /= 3600.0; // Convert hours to seconds
            }
            const Ea_r_J = convertUnit.energy.convertEaToJules(parseFloat(reactionInfo.EaValueBackward || '0'), reactionInfo.EaUnitBackward || 'J/mol');
            const k_r_c = A_r * Math.exp(-Ea_r_J / (R * tempK));
            
            let k_r: number;
            if (usePartialPressure) {
                let totalOrderReverse = 0;
                components.forEach(comp => {
                    const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                    if (stoich > 0) { // It's a product
                        totalOrderReverse += parseFloat(comp.reactionData[reactionInfo.id]?.orderReverse || '0');
                    }
                });
                k_r = convertUnit.kinetics.kcToKp(k_r_c, tempK, totalOrderReverse);
            } else {
                k_r = k_r_c;
            }
            
            rate_r = k_r;
            components.forEach(comp => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                const orderReverse = parseFloat(comp.reactionData[reactionInfo.id]?.orderReverse || '0');
                
                if (orderReverse > 0 && stoich > 0) { // Product with non-zero reverse order
                     if (usePartialPressure) {
                        rate_r *= Math.pow(Math.max(0, partialPressures_Pa[comp.name]), orderReverse);
                     } else {
                        rate_r *= Math.pow(Math.max(0, concentrations_L[comp.name]), orderReverse);
                     }
                }
            });
        }

        return rate_f - rate_r; // Returns the net rate for the reaction
    });
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
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId)
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId)
    if (!limitingReactant || !desiredProduct) return { selectivityData: [], volumeData: [] }

    // --- START: ADD THIS SECTION ---
    const primaryReactionId = reactions[0].id;
    const stoichCoeffStr = limitingReactant.reactionData[primaryReactionId]?.stoichiometry;
    const stoichFactor = Math.abs(parseFloat(stoichCoeffStr || '0'));
    // --- END: ADD THIS SECTION ---

    const rawSelectivityData: { x: number; y: number }[] = []
    const rawVolumeData: { x: number; y: number }[] = []
    
    const P_Pa = pressure_bar * 1e5;
    const limitingReactantIndex = components.findIndex(c => c.id === simBasis.limitingReactantId);
    const desiredProductIndex = components.findIndex(c => c.id === simBasis.desiredProductId);

    const dFdV = (V: number, F_vec: number[]) => {
        const F_dict = components.reduce((acc, comp, i) => { acc[comp.name] = F_vec[i]; return acc; }, {} as { [key: string]: number });
        const F_tot = F_vec.reduce((sum, f) => sum + f, 0);
        if (F_tot < 1e-12) return Array(components.length).fill(0);

        const concentrations_m3 = components.reduce((acc, comp) => {
            acc[comp.name] = (F_dict[comp.name] * P_Pa) / (F_tot * R * tempK);
            return acc;
        }, {} as { [key: string]: number });

        const partialPressures_Pa = components.reduce((acc, comp) => {
            acc[comp.name] = concentrations_m3[comp.name] * R * tempK;
            return acc;
        }, {} as { [key: string]: number });
        
        const concentrations_L: { [key: string]: number } = {};
        for (const key in concentrations_m3) {
            concentrations_L[key] = concentrations_m3[key] * 1e-3;
        }

        // The function now calls the new helper to get rates
        const reactionNetRates = calculateNetReactionRates(reactions, components, concentrations_L, partialPressures_Pa, tempK);
        
        return components.map(comp => {
            let netRateOfFormation = 0; // This will be in mol/(L*s)
            reactions.forEach((reactionInfo, j) => {
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                netRateOfFormation += stoich * reactionNetRates[j];
            });
            
            // MODIFICATION: Convert rate from per-Liter to per-m³ for the solver
            return netRateOfFormation * 1000;
        });
    };

    // The rest of the function remains unchanged...
    const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
    const P_C_mol_s = convertUnit.flowRate.ktaToMolPerSecond(targetProduction_kta, productMolarMass);

    const F0_basis_dict: { [key: string]: number } = { [limitingReactant.name]: 1.0 };
    molarRatios.forEach(ratio => {
        const compName = components.find(c => c.id === ratio.numeratorId)?.name;
        if (compName) F0_basis_dict[compName] = ratio.value;
    });
    const F0_basis_vec = components.map(c => F0_basis_dict[c.name] || 0);

    let F_current = F0_basis_vec;
    let V_current = 0;

    const logV_min = -8;
    const logV_max = 8;
    const n_steps = 500;
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
    const molesFormed = F_P_current_val - (F0_basis_vec[desiredProductIndex] || 0)
    // --- START: MODIFY THIS LINE ---
    const selectivity = molesConsumed > 1e-9 ? (stoichFactor * molesFormed) / molesConsumed : 0;
    // --- END: MODIFY THIS LINE ---
        
        if (conversion >= 1e-6 && conversion < 0.9999 && selectivity > 1e-6) {
            const F_A0_true   = P_C_mol_s / (selectivity*conversion);
            const F_tot0_basis = F0_basis_vec.reduce((s,f)=>s+f,0);
            const F_tot0_true = F_A0_true / (F0_basis_vec[limitingReactantIndex] / F_tot0_basis);
            const V_true = V_current * (F_tot0_true / F_tot0_basis); // Corrected Line
            rawSelectivityData.push({ x: conversion, y: selectivity });
            rawVolumeData.push({ x: conversion, y: V_true });
        }
        if (conversion >= 0.999) break;
    }
    
    const commonConversionGrid = Array.from({ length: 100 }, (_, i) => i * 0.01)
    const selectivityData = commonConversionGrid.map(conv => ({x: conv, y: interpolateFromData(conv, rawSelectivityData)})).filter(p => !isNaN(p.y))
    const volumeData = commonConversionGrid.map(conv => ({x: conv, y: interpolateFromData(conv, rawVolumeData)})).filter(p => !isNaN(p.y))
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

    // --- START: ADD THIS SECTION ---
    // Find the stoichiometric coefficient of the limiting reactant in the primary reaction (reaction[0])
    const primaryReactionId = reactions[0].id;
    const stoichCoeffStr = limitingReactant.reactionData[primaryReactionId]?.stoichiometry;
    const stoichFactor = Math.abs(parseFloat(stoichCoeffStr || '0'));
    // --- END: ADD THIS SECTION ---

    const rawSelectivityData: { x: number; y: number }[] = [];
    const rawVolumeData: { x: number; y: number }[] = [];
    
    const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
    const P_C_target_mols = convertUnit.flowRate.ktaToMolPerSecond(targetProduction_kta, productMolarMass);

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
            const Ea_J_per_mol = convertUnit.energy.convertEaToJules(
                parseFloat(reactionInfo.EaValue), 
                reactionInfo.EaUnit
            );
            const k_f = parseFloat(reactionInfo.AValue) * Math.exp(-Ea_J_per_mol / (R * tempK));
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

    // --- START: MODIFY THIS LINE ---
    const molesConsumed = F0_basis_vec[limitingReactantIndex] - F_A_out_basis;
    const selectivity = (molesConsumed > 1e-9) ? (stoichFactor * P_C_basis_mols) / molesConsumed : 0;
    // --- END: MODIFY THIS LINE ---
        
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
    const scaled_volume_liters = V_current * (q_true / q_basis);
    const V_true = (q_basis > 1e-9 && q_true > 0) ? convertUnit.volume.litersToCubicMeters(scaled_volume_liters) : Infinity;
        
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
        const Ea_J_per_mol = convertUnit.energy.convertEaToJules(
            parseFloat(reactionInfo.EaValue), 
            reactionInfo.EaUnit
        );
        const k_f = parseFloat(reactionInfo.AValue) * Math.exp(-Ea_J_per_mol / (R * tempK));
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
    reactionPhase, setReactionPhase,
    prodRate, setProdRate,
    setTemperature, // New Prop
    setPressure,    // New Prop
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
    setProdRate: React.Dispatch<React.SetStateAction<string>>,
    setTemperature: React.Dispatch<React.SetStateAction<number>>, // New Prop Type
    setPressure: React.Dispatch<React.SetStateAction<number>>,    // New Prop Type
}) => {
  const { resolvedTheme } = useTheme();
  
  const addReactionSetup = () => {
    if (reactionsSetup.length >= 4) return
    const newId = (Math.max(0, ...reactionsSetup.map(r => parseInt(r.id))) + 1).toString()
    setReactionsSetup([...reactionsSetup, { id: newId, AValue: '1e6', AUnit: 'mol/(L*s)', EaValue: '50000', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' }])
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

  // --- HDA Preset Data ---
  const loadHdaPreset = () => {
    setReactionsSetup([
        // HDA Rxn 1 (Gas phase, Concentration basis -> simple units)
        { id: '1', AValue: '1.188e14', AUnit: 'mol,L,h', EaValue: '52.0', EaUnit: 'kcal/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
        // HDA Rxn 2 (Gas phase, Pressure basis -> pressure units)
        { id: '2', AValue: '1.1647e-5', AUnit: 'mol,L,Pa,h', EaValue: '30.19', EaUnit: 'kcal/mol', isEquilibrium: true, AValueBackward: '4.56e-7', AUnitBackward: 'mol,L,Pa,h', EaValueBackward: '23.35', EaUnitBackward: 'kcal/mol', rateLawBasis: 'partialPressure'}
    ]);
    setComponentsSetup([
        { id: 'comp-T', name: 'Toluene', molarMass: '92.14', density: '867', reactionData: { '1': { stoichiometry: '-1', order: '1' } } },
        { id: 'comp-H2', name: 'Hydrogen', molarMass: '2.016', density: '89.9', reactionData: { '1': { stoichiometry: '-1', order: '0.5' }, '2': { stoichiometry: '1', order: '0', orderReverse: '1' } } },
        { id: 'comp-B', name: 'Benzene', molarMass: '78.11', density: '876', reactionData: { '1': { stoichiometry: '1', order: '0' }, '2': { stoichiometry: '-2', order: '2' } } },
        { id: 'comp-CH4', name: 'Methane', molarMass: '16.04', density: '717', reactionData: { '1': { stoichiometry: '1', order: '0' } } },
        { id: 'comp-D', name: 'Diphenyl', molarMass: '154.21', density: '1080', reactionData: { '2': { stoichiometry: '1', order: '0', orderReverse: '1' } } },
    ]);
    setSimBasis({ limitingReactantId: 'comp-T', desiredProductId: 'comp-B' });
    setReactionPhase('Gas');
    setProdRate('100');
        setTemperature(650); // MODIFICATION: Set default temperature to 650 C
        setPressure(30);    // MODIFICATION: Set default pressure to 30 bar
  };

  const loadAlkylationPreset = () => {
    setReactionsSetup(initialReactions);
    setComponentsSetup(initialComponents);
    setSimBasis({ limitingReactantId: 'comp-A', desiredProductId: 'comp-C' });
    setReactionPhase('Liquid');
    setProdRate('100');
    setTemperature(4); // MODIFICATION: Corrected from 10 to 4
    setPressure(10);
  }
  
  const loadDmcPreset = () => {
      setReactionsSetup([
          // DMC Synthesis (Liquid phase -> simple units)
          { id: '1', AValue: '1.6e11', AUnit: 'mol,L,s', EaValue: '104600', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
          { id: '2', AValue: '6.7e12', AUnit: 'mol,L,s', EaValue: '98750', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
      ]);
      setComponentsSetup([
          { id: 'comp-MeOH', name: 'Methanol', molarMass: '32.04', density: '792', reactionData: { '1': { stoichiometry: '-2', order: '2' } } },
          { id: 'comp-CO', name: 'Carbon Monoxide', molarMass: '28.01', density: '789', reactionData: { '1': { stoichiometry: '-1', order: '0' }, '2': { stoichiometry: '-1', order: '0' } } },
          { id: 'comp-O2', name: 'Oxygen', molarMass: '32.00', density: '1141', reactionData: { '1': { stoichiometry: '-0.5', order: '0.5' }, '2': { stoichiometry: '-0.5', order: '0.5' } } },
          { id: 'comp-DMC', name: 'DMC', molarMass: '90.08', density: '1069', reactionData: { '1': { stoichiometry: '1', order: '0' } } },
          { id: 'comp-H2O', name: 'Water', molarMass: '18.02', density: '1000', reactionData: { '1': { stoichiometry: '1', order: '0' } } },
          { id: 'comp-CO2', name: 'Carbon Dioxide', molarMass: '44.01', density: '1101', reactionData: { '2': { stoichiometry: '1', order: '0' } } },
      ]);
      setSimBasis({ limitingReactantId: 'comp-O2', desiredProductId: 'comp-DMC' });
      setReactionPhase('Liquid');
      setProdRate('50');
      setTemperature(80); // Set default temperature for next screen
      setPressure(25);    // Set default pressure for next screen
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
        
        const rateSymbol = reactionPhase === 'Gas' && reaction.rateLawBasis === 'partialPressure' ? 'P' : 'C';

        const rateLawReactants = componentsSetup
            .filter(c => c.name && parseFloat(c.reactionData[reaction.id]?.order || '0') !== 0)
            .map(c => {
                const order = c.reactionData[reaction.id]?.order;
                const exponent = (order && order !== '1') ? <sup>{order}</sup> : null;
                return (
                    <React.Fragment key={c.id}>
                        {' '}{rateSymbol}<sub>{c.name}</sub>
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
                        {' '}{rateSymbol}<sub>{c.name}</sub>
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

  const componentsGridCols = `40px 1fr 1fr 1fr repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`;

  return (
    <div className={`min-h-screen flex flex-col p-4 bg-background text-foreground`}>
      <div className="container mx-auto space-y-6">
        
        {/* --- NEW HEADER ROW --- */}
        <div className="relative flex justify-center items-center py-2 mb-4">
            {/* Centered Preset Buttons */}
            <div className="flex justify-center items-center gap-4">
                <Label className="text-sm font-semibold text-muted-foreground">Load a Preset:</Label>
                <Button variant="outline" onClick={loadAlkylationPreset}>Butane Alkylation</Button>
                <Button variant="outline" onClick={loadHdaPreset}>HDA</Button>
                <Button variant="outline" onClick={loadDmcPreset}>DMC Synthesis</Button>
            </div>
            
            {/* Right-aligned Next Button */}
                    <div className="absolute right-0">
                        <div title={!validationResult.isValid ? validationResult.message : undefined} className="inline-block">
                            <Button 
                                onClick={onNext} 
                                disabled={!validationResult.isValid}
                                variant="secondary" // This variant is theme-aware
                                className="px-4 py-2 rounded-lg font-bold flex items-center"
                                aria-label="Next"
                            >
                                Next <ArrowRight className="h-4 w-4" />
                            </Button>
                        </div>
                    </div>
        </div>

        <div className="grid grid-cols-1 lg:grid-cols-12 gap-6">
            {/* Reactions & Kinetics Card */}
            <div className={`lg:col-span-6`}> 
                {/* MODIFICATION: Added h-full to the Card component */}
                <Card className="p-6 h-full">
                    <div className="mb-0">
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
                            {/* MODIFICATION: Displaying new reaction names */}
                            <Label className="font-bold">
                                {(() => {
                                    const numSideReactions = reactionsSetup.length - 1;
                                    if (index === 0) return 'Primary Reaction';
                                    return numSideReactions === 1 ? 'Side Reaction' : `Side Reaction ${index}`;
                                })()}
                            </Label>
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
                            {/* A Value Input with Unit Selector */}
                            <div className="flex items-center gap-2">
                                <Label className="w-16 text-sm">
                                    <span>A<sub className="font-medium relative -left-px">{reaction.id}{reaction.isEquilibrium ? ',f' : ''}</sub>:</span>
                                </Label>
                                <Input type="text" value={reaction.AValue} onChange={(e) => updateReactionSetup(reaction.id, 'AValue', e.target.value)} />
                                
                                {/* --- MODIFICATION: Conditional Rendering of the Select Menu --- */}
                                {(() => {
                                    const showPressureUnits = reactionPhase === 'Gas' && reaction.rateLawBasis === 'partialPressure';
                                    return showPressureUnits ? (
                                        <Select value={reaction.AUnit} onValueChange={(value) => updateReactionSetup(reaction.id, 'AUnit', value)}>
                                            <SelectTrigger className="w-[180px] text-xs"><SelectValue /></SelectTrigger>
                                            <SelectContent>
                                                {/* Pressure-based units */}
                                                <SelectItem value="mol,L,Pa,s">mol, L, Pa, s</SelectItem>
                                                <SelectItem value="mol,m^3,Pa,s">mol, m³, Pa, s</SelectItem>
                                                <SelectItem value="mol,L,Pa,h">mol, L, Pa, h</SelectItem>
                                                <SelectItem value="mol,m^3,Pa,h">mol, m³, Pa, h</SelectItem>
                                                <SelectItem value="mol,L,bar,s">mol, L, bar, s</SelectItem>
                                                <SelectItem value="mol,m^3,bar,s">mol, m³, bar, s</SelectItem>
                                                <SelectItem value="mol,L,bar,h">mol, L, bar, h</SelectItem>
                                                <SelectItem value="mol,m^3,bar,h">mol, m³, bar, h</SelectItem>
                                            </SelectContent>
                                        </Select>
                                    ) : (
                                        <Select value={reaction.AUnit} onValueChange={(value) => updateReactionSetup(reaction.id, 'AUnit', value)}>
                                            <SelectTrigger className="w-[180px] text-xs"><SelectValue /></SelectTrigger>
                                            <SelectContent>
                                                {/* Simple units */}
                                                <SelectItem value="mol,L,s">mol, L, s</SelectItem>
                                                <SelectItem value="mol,m^3,s">mol, m³, s</SelectItem>
                                                <SelectItem value="mol,L,h">mol, L, h</SelectItem>
                                                <SelectItem value="mol,m^3,h">mol, m³, h</SelectItem>
                                            </SelectContent>
                                        </Select>
                                    );
                                })()}
                            </div>
                            {/* Ea Value Input with Unit Selector */}
                            <div className="flex items-center gap-2">
                                <Label className="w-16 text-sm">
                                    <span>Ea<sub className="font-medium relative -left-px">{reaction.id}{reaction.isEquilibrium ? ',f' : ''}</sub>:</span>
                                </Label>
                                <Input type="text" value={reaction.EaValue} onChange={(e) => updateReactionSetup(reaction.id, 'EaValue', e.target.value)} />
                                <Select value={reaction.EaUnit} onValueChange={(value) => updateReactionSetup(reaction.id, 'EaUnit', value)}>
                                    <SelectTrigger className="w-[120px] text-xs"><SelectValue /></SelectTrigger>
                                    <SelectContent>
                                        <SelectItem value="J/mol">J/mol</SelectItem>
                                        <SelectItem value="kJ/mol">kJ/mol</SelectItem>
                                        <SelectItem value="cal/mol">cal/mol</SelectItem>
                                        <SelectItem value="kcal/mol">kcal/mol</SelectItem>
                                    </SelectContent>
                                </Select>
                            </div>
                        </div>
                        {/* Conditionally show reverse reaction inputs with unit selectors */}
                        {reaction.isEquilibrium && (
                            <div className="grid grid-cols-1 md:grid-cols-2 gap-4 pt-4 border-t border-dashed">
                                {/* A Value Reverse */}
                                <div className="flex items-center gap-2">
                                    <Label className="w-16 text-sm">
                                        <span>A<sub className="font-medium relative -left-px">{reaction.id},r</sub>:</span>
                                    </Label>
                                    <Input type="text" value={reaction.AValueBackward || ''} placeholder='e.g. 1e12' onChange={(e) => updateReactionSetup(reaction.id, 'AValueBackward', e.target.value)} />
                                    
                                    {/* --- MODIFICATION: Conditional Rendering of the Select Menu --- */}
                                    {(() => {
                                        const showPressureUnits = reactionPhase === 'Gas' && reaction.rateLawBasis === 'partialPressure';
                                        return showPressureUnits ? (
                                            <Select value={reaction.AUnitBackward} onValueChange={(value) => updateReactionSetup(reaction.id, 'AUnitBackward', value)}>
                                                <SelectTrigger className="w-[180px] text-xs"><SelectValue /></SelectTrigger>
                                                <SelectContent>
                                                    {/* Pressure-based units */}
                                                    <SelectItem value="mol,L,Pa,s">mol, L, Pa, s</SelectItem>
                                                    <SelectItem value="mol,m^3,Pa,s">mol, m³, Pa, s</SelectItem>
                                                    <SelectItem value="mol,L,Pa,h">mol, L, Pa, h</SelectItem>
                                                    <SelectItem value="mol,m^3,Pa,h">mol, m³, Pa, h</SelectItem>
                                                    <SelectItem value="mol,L,bar,s">mol, L, bar, s</SelectItem>
                                                    <SelectItem value="mol,m^3,bar,s">mol, m³, bar, s</SelectItem>
                                                    <SelectItem value="mol,L,bar,h">mol, L, bar, h</SelectItem>
                                                    <SelectItem value="mol,m^3,bar,h">mol, m³, bar, h</SelectItem>
                                                </SelectContent>
                                            </Select>
                                        ) : (
                                            <Select value={reaction.AUnitBackward} onValueChange={(value) => updateReactionSetup(reaction.id, 'AUnitBackward', value)}>
                                                <SelectTrigger className="w-[180px] text-xs"><SelectValue /></SelectTrigger>
                                                <SelectContent>
                                                    {/* Simple units */}
                                                    <SelectItem value="mol,L,s">mol, L, s</SelectItem>
                                                    <SelectItem value="mol,m^3,s">mol, m³, s</SelectItem>
                                                    <SelectItem value="mol,L,h">mol, L, h</SelectItem>
                                                    <SelectItem value="mol,m^3,h">mol, m³, h</SelectItem>
                                                </SelectContent>
                                            </Select>
                                        );
                                    })()}
                                </div>
                                {/* Ea Value Reverse */}
                                <div className="flex items-center gap-2">
                                    <Label className="w-16 text-sm">
                                        <span>Ea<sub className="font-medium relative -left-px">{reaction.id},r</sub>:</span>
                                    </Label>
                                    <Input type="text" value={reaction.EaValueBackward || ''} placeholder='e.g. 80000' onChange={(e) => updateReactionSetup(reaction.id, 'EaValueBackward', e.target.value)} />
                                    <Select value={reaction.EaUnitBackward} onValueChange={(value) => updateReactionSetup(reaction.id, 'EaUnitBackward', value)}>
                                        <SelectTrigger className="w-[120px] text-xs"><SelectValue /></SelectTrigger>
                                        <SelectContent>
                                            <SelectItem value="J/mol">J/mol</SelectItem>
                                            <SelectItem value="kJ/mol">kJ/mol</SelectItem>
                                            <SelectItem value="cal/mol">cal/mol</SelectItem>
                                            <SelectItem value="kcal/mol">kcal/mol</SelectItem>
                                        </SelectContent>
                                    </Select>
                                </div>
                            </div>
                        )}
                    </div>
                    ))}
                </div>
                </Card>
            </div>

            {/* Simulation Basis Card */}
            <div className={`lg:col-span-3`}>
                {/* MODIFICATION: Added h-full to the Card component */}
                <Card className="p-6 h-full">
                    <div className="mb-0">
                        <CardTitle>Simulation Basis</CardTitle>
                    </div>
                <div className="space-y-4">
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
                                {/* Only show reactants from the Primary Reaction */}
                                {componentsSetup
                                    .filter(c => parseFloat(c.reactionData[reactionsSetup[0]?.id]?.stoichiometry || '0') < 0)
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
                                {/* Only show products from the Primary Reaction */}
                                {componentsSetup
                                    .filter(c => parseFloat(c.reactionData[reactionsSetup[0]?.id]?.stoichiometry || '0') > 0)
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
                </Card>
            </div>

            {/* Parsed Reactions & Rate Laws Card */}
            <div className="lg:col-span-3 flex flex-col">
                {/* MODIFICATION: Added h-full to the Card component */}
                <Card className="h-full flex flex-col"> 
                    <CardHeader>
                        <CardTitle>Parsed Reactions & Rate Laws</CardTitle>
                    </CardHeader>
                    <CardContent className="flex-grow space-y-4">
                        {reactionsSetup.map((r, index) => {
                            const preview = generatePreview(r)
                            return (
                                                                <div key={r.id} className="p-3 bg-card rounded-md">
                                                                                                        <div className={`flex w-full ${reactionPhase === 'Gas' ? 'flex-row justify-between items-center' : 'flex-col items-center'}`}>
                                                                                                            <p className={`font-bold text-primary ${reactionPhase === 'Gas' ? 'text-left' : 'text-center'} mb-0`}>
                                                                                {(() => {
                                                                                        const numSideReactions = reactionsSetup.length - 1;
                                                                                        if (index === 0) return 'Primary Reaction';
                                                                                        return numSideReactions === 1 ? 'Side Reaction' : `Side Reaction ${index}`;
                                                                                })()}
                                                                        </p>
                                    
                                    {/* The switch container is now always rendered but conditionally hidden */}
                                    <div className={`${reactionPhase !== 'Gas' ? 'hidden' : 'flex items-center gap-2'}`}>
                                        <div className="flex items-center gap-1 rounded-lg p-1 bg-muted">
                                            <Button
                                                onClick={() => updateReactionSetup(r.id, 'rateLawBasis', 'concentration')}
                                                variant={(r.rateLawBasis === 'concentration' || !r.rateLawBasis) ? 'default' : 'ghost'}
                                                size="sm"
                                                className="text-xs px-2 py-1 h-auto"
                                            >
                                                Conc.
                                            </Button>
                                            <Button
                                                onClick={() => updateReactionSetup(r.id, 'rateLawBasis', 'partialPressure')}
                                                variant={r.rateLawBasis === 'partialPressure' ? 'default' : 'ghost'}
                                                size="sm"
                                                className="text-xs px-2 py-1 h-auto"
                                            >
                                                Pressure
                                            </Button>
                                        </div>
                                    </div>

                                  </div>
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
            </div>
        </div>

        {/* Components Card (unchanged) */}
        <div className={`p-6 rounded-lg shadow-lg bg-card text-card-foreground`}> 
              <div className="mb-4">
                <CardTitle className="flex items-center justify-between">
                  Components
                  <Button variant="outline" size="sm" onClick={addComponentSetup}>
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
                                {reactionsSetup.map((r, index) => {
                                    const numSideReactions = reactionsSetup.length - 1;
                                    let reactionTitle = 'Primary Reaction';
                                    if (index > 0) {
                                        reactionTitle = numSideReactions === 1 ? 'Side Reaction' : `Side Reaction ${index}`;
                                    }
                                    return (
                                        <div key={r.id} className="text-center col-span-2">
                                            {reactionTitle}
                                        </div>
                                    );
                                })}
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
      </div>
    </div>
  )
}

// Interface for the new molarRatios prop
interface MolarRatio {
    numeratorId: string;
    value: number;
}

const formatTooltipNumber = (num: number | string, sigFigs: number = 3): string => {
    const numFloat = typeof num === 'string' ? parseFloat(num) : num;
    if (isNaN(numFloat) || numFloat === 0) return '0';

    // Get the number to 3 significant figures
    const preparedNum = numFloat.toPrecision(sigFigs);

    // Convert back to a number to drop trailing zeros, then back to a string
    return Number(preparedNum).toString();
}

const ReactorSimulator = ({ 
    onBack, 
    reactions, components, 
    simBasis, 
    molarRatios, setMolarRatios,
    reactionPhase,
    prodRate,
    temperature, setTemperature, // Now receives state and setter
    pressure, setPressure,       // Now receives state and setter
}: { 
    onBack: () => void, 
    reactions: ReactionSetup[], 
    components: ComponentSetup[],
    simBasis: any,
    molarRatios: MolarRatio[],
    setMolarRatios: React.Dispatch<React.SetStateAction<MolarRatio[]>>,
    reactionPhase: 'Liquid' | 'Gas',
    prodRate: string,
    temperature: number,
    setTemperature: React.Dispatch<React.SetStateAction<number>>,
    pressure: number,
    setPressure: React.Dispatch<React.SetStateAction<number>>,
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
  
  const [pressureMin, setPressureMin] = useState('1');
  const [pressureMax, setPressureMax] = useState('50');
  
  const [molarRatioMin, setMolarRatioMin] = useState('1');
  const [molarRatioMax, setMolarRatioMax] = useState('25');
  const [tempMin, setTempMin] = useState('0');
  const [tempMax, setTempMax] = useState('100');
  // Removed local prodRate state

  const generateGraphData = useCallback(() => {
    const tempK = convertUnit.temperature.celsiusToKelvin(temperature)
    
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
                if (volumeData.length > 1) { // Need at least 2 points to interpolate
                    // Interpolate to get volume at both low (5%) and high (95%) conversion
                    const volAt5Percent = interpolateFromData(0.05, volumeData);
                    const volAt95Percent = interpolateFromData(0.95, volumeData);

                    // Determine the larger of the two volumes, ignoring any non-finite results
                    let curveReasonableMax = -1;
                    if (isFinite(volAt5Percent) && isFinite(volAt95Percent)) {
                        curveReasonableMax = Math.max(volAt5Percent, volAt95Percent);
                    } else if (isFinite(volAt5Percent)) {
                        curveReasonableMax = volAt5Percent;
                    } else if (isFinite(volAt95Percent)) {
                        curveReasonableMax = volAt95Percent;
                    }
                    
                    // If this curve's reasonable max is higher than the overall max we've seen so far,
                    // it becomes the new overall max for the chart's Y-axis.
                    if (curveReasonableMax > overallReasonableMax) {
                        overallReasonableMax = curveReasonableMax;
                    }
                }
            }
            // --- **FIX END** ---

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
    
    // Set the final Y-axis max value with a hard cap
    let yAxisMax: number | 'dataMax' = 'dataMax';
    if (graphType === 'volume') {
        const hardCap = 5000; // The maximum practical volume
        if (overallReasonableMax > 0) {
            // The axis maximum will be the smaller of the calculated max or the hard cap.
            // This keeps the graph nicely scaled for smaller volumes but caps it for large ones.
            yAxisMax = Math.min(overallReasonableMax, hardCap);
        } else {
            // If all calculated volumes are enormous, just default the axis to the cap.
            yAxisMax = hardCap;
        }
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

         // MODIFICATION: Using the new helper for all numbers
         const xAxisValue = formatTooltipNumber(params[0].axisValue);
         let tooltipContent = `<strong>Conversion</strong>: ${xAxisValue}`; 

         params.forEach((p: any) => { 
             const seriesName = p.seriesName;
             const rawValue = parseFloat(p.value[1]);
             let formattedValue = ''; 

             if (graphType === 'selectivity') { 
                 formattedValue = formatTooltipNumber(rawValue);
             } else { // 'volume' 
                 formattedValue = formatTooltipNumber(rawValue);
                 formattedValue += ' m³'; 
             }
             
             const seriesColor = p.color;

            tooltipContent += `<br/>`;
             tooltipContent += `<span style="color: ${seriesColor}; font-weight: bold;">${seriesName}:</span> <strong>${formattedValue}</strong>`;
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
                    // MODIFICATION: Removed pt-4 for consistent spacing
                    <div key={ratio.numeratorId}>
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

            {/* MODIFICATION: Removed wrapper divs with borders for the sliders below */}
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

export default function Home() {
  const [showSimulator, setShowSimulator] = useState(false)
  
  const [reactionsSetup, setReactionsSetup] = useState<ReactionSetup[]>(initialReactions);
  const [componentsSetup, setComponentsSetup] = useState<ComponentSetup[]>(initialComponents);
  
  const [simBasis, setSimBasis] = useState({
    limitingReactantId: 'comp-A',
    desiredProductId: 'comp-C',
  });

  const [molarRatios, setMolarRatios] = useState<Array<{numeratorId: string, value: number}>>([]);
  
  // State for parameters now lives in the parent component
  const [prodRate, setProdRate] = useState('100');
  const [reactionPhase, setReactionPhase] = useState<'Liquid' | 'Gas'>('Liquid'); 
  const [temperature, setTemperature] = useState(4); // State lifted here
  const [pressure, setPressure] = useState(25);       // State lifted here 

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
    if (!reactionsSetup || reactionsSetup.length === 0 || !componentsSetup || componentsSetup.length === 0) {
        return;
    }

    // Identify the Primary Reaction (always the first one)
    const primaryReactionId = reactionsSetup[0].id;

    // Find all reactants ONLY in the Primary Reaction
    const primaryReactants = componentsSetup.filter(comp => {
        const stoich = parseFloat(comp.reactionData[primaryReactionId]?.stoichiometry || '0');
        return stoich < 0;
    });

    // Filter out the selected limiting reactant to find the others
    const otherPrimaryReactants = primaryReactants.filter(r => r.id !== simBasis.limitingReactantId);

    // Create a slider for each remaining primary reactant
    const newRatios = otherPrimaryReactants.map(reactant => {
        const existingRatio = molarRatios.find(r => r.numeratorId === reactant.id);
        let defaultValue = 1.0;

        // --- MODIFICATION: Set default ratios based on preset ---
        const limitingReactantName = componentsSetup.find(c => c.id === simBasis.limitingReactantId)?.name;

        if (limitingReactantName === 'Toluene' && reactant.name === 'Hydrogen') {
            defaultValue = 5.0; // HDA Preset
        } else if (limitingReactantName === '1-Butene' && reactant.name === 'Isobutane') {
            defaultValue = 6.0; // Butane Alkylation Preset
        } else if (limitingReactantName === 'Oxygen') {
            if (reactant.name === 'Methanol') defaultValue = 12.0; // DMC Preset
            if (reactant.name === 'Carbon Monoxide') defaultValue = 24.0; // DMC Preset
        }

        return {
            numeratorId: reactant.id,
            value: existingRatio ? existingRatio.value : defaultValue 
        };
    });

    if (JSON.stringify(newRatios) !== JSON.stringify(molarRatios)) {
        setMolarRatios(newRatios);
    }

  }, [reactionsSetup, componentsSetup, simBasis.limitingReactantId, molarRatios]);


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

  // 6. Context-aware validation of AUnit selection based on phase and rate law basis
    // 6. Context-aware validation and auto-swap of AUnit selection based on rateLawBasis changes
    useEffect(() => {
        if (!reactionsSetup || reactionsSetup.length === 0) return;

        const swapUnit = (unit: string, toPressure: boolean) => {
            // Only swap if unit is in the expected set
            if (toPressure) {
                // concentration → pressure
                if (unit === 'mol,L,s') return 'mol,L,Pa,s';
                if (unit === 'mol,L,h') return 'mol,L,Pa,h';
                if (unit === 'mol,m^3,s') return 'mol,m^3,Pa,s';
                if (unit === 'mol,m^3,h') return 'mol,m^3,Pa,h';
            } else {
                // pressure → concentration
                if (unit === 'mol,L,Pa,s' || unit === 'mol,L,bar,s') return 'mol,L,s';
                if (unit === 'mol,L,Pa,h' || unit === 'mol,L,bar,h') return 'mol,L,h';
                if (unit === 'mol,m^3,Pa,s' || unit === 'mol,m^3,bar,s') return 'mol,m^3,s';
                if (unit === 'mol,m^3,Pa,h' || unit === 'mol,m^3,bar,h') return 'mol,m^3,h';
            }
            return unit;
        };

        const updatedReactions = reactionsSetup.map(reaction => {
            let updated = { ...reaction };
            // Only swap for gas phase
            if (reactionPhase === 'Gas') {
                // Forward
                if (reaction.rateLawBasis === 'partialPressure') {
                    updated.AUnit = swapUnit(reaction.AUnit, true);
                } else if (reaction.rateLawBasis === 'concentration') {
                    updated.AUnit = swapUnit(reaction.AUnit, false);
                }
                // Reverse (if equilibrium)
                        if (reaction.isEquilibrium) {
                            if (reaction.rateLawBasis === 'partialPressure') {
                                updated.AUnitBackward = swapUnit(reaction.AUnitBackward ?? '', true);
                            } else if (reaction.rateLawBasis === 'concentration') {
                                updated.AUnitBackward = swapUnit(reaction.AUnitBackward ?? '', false);
                            }
                        }
            }
            return updated;
        });
        // Only update if something changed
        if (JSON.stringify(updatedReactions) !== JSON.stringify(reactionsSetup)) {
            setReactionsSetup(updatedReactions);
        }
    }, [reactionPhase, reactionsSetup]);

  if (showSimulator) {
    return <ReactorSimulator 
              onBack={() => setShowSimulator(false)} 
              reactions={reactionsSetup} 
              components={componentsSetup} 
              simBasis={simBasis}
              molarRatios={molarRatios}
              setMolarRatios={setMolarRatios}
              reactionPhase={reactionPhase}
              prodRate={prodRate}
              temperature={temperature}     // Pass state down
              setTemperature={setTemperature} // Pass setter down
              pressure={pressure}           // Pass state down
              setPressure={setPressure}       // Pass setter down
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
              reactionPhase={reactionPhase}      
              setReactionPhase={setReactionPhase}
              prodRate={prodRate}
              setProdRate={setProdRate}
              setTemperature={setTemperature} // Pass setter down
              setPressure={setPressure}       // Pass setter down
            />
  }
}