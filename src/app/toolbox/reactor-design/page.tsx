
"use client"

let concentrations: { [key: string]: number } = {};

import React, { useState, useMemo, useEffect, useCallback, useRef } from 'react'
import ReactECharts from 'echarts-for-react'
import { EChartsOption } from 'echarts'
import { PlusCircle, Trash2, ArrowLeft, ArrowRight, Unlink2, Info } from 'lucide-react'
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

import {
    Tooltip,
    TooltipContent,
    TooltipProvider,
    TooltipTrigger,
} from '@/components/ui/tooltip'

import { createClient, SupabaseClient } from '@supabase/supabase-js';
import type { CompoundData } from '@/lib/vle-types';

// Import BDF solver from reactor-solver (aliased to avoid name conflict)
import { solveODE_BDF as solveODE_BDF_imported, solveLinearSystem as solveLinearSystem_imported } from '@/lib/reactor-solver';

import {
    // All VLE calculations and data fetching now come from the multicomponent file
    calculateNrtlGammaMulticomponent,
    calculateWilsonGammaMulticomponent,
    calculateUniquacGammaMulticomponent,
    calculatePrFugacityCoefficientsMulticomponent,
    calculateSrkFugacityCoefficientsMulticomponent,

    // Fetching functions (now moved into the multicomponent file)
    fetchNrtlParameters,
    fetchWilsonInteractionParams,
    fetchPrInteractionParams,
    fetchSrkInteractionParams,
    fetchUniquacInteractionParams,

    // Type definitions for parameter objects and matrices
    NrtlInteractionParams,
    WilsonInteractionParams,
    PrSrkInteractionParams, // Combined type for PR and SRK
    UniquacInteractionParams,
    NrtlParameterMatrix,
    WilsonParameterMatrix,
    UniquacParameterMatrix,
    PrSrkParameterMatrix,
} from '@/lib/vle-calculations-multicomponent';

// The generic InteractionParams type needs to be redefined here
type InteractionParams = NrtlInteractionParams | WilsonInteractionParams | PrSrkInteractionParams | UniquacInteractionParams;

// --- Supabase Client Initialization ---
const supabaseUrl = process.env.NEXT_PUBLIC_SUPABASE_URL || '';
const supabaseAnonKey = process.env.NEXT_PUBLIC_SUPABASE_ANON_KEY || '';
let supabase: SupabaseClient;
if (supabaseUrl && supabaseAnonKey) {
    try {
        supabase = createClient(supabaseUrl, supabaseAnonKey);
    } catch (error) {
        console.error('Error initializing Supabase client:', error);
    }
} else {
    console.error('Supabase URL or Anon Key is missing.');
}

// --- Property Calculation Utilities ---
// NOTE: Please replace this placeholder with your actual function from '@/lib/property-equations'.
// This example implements a simplified version for liquid density (DIPPR 105).
const calculatePropertyByEquation = (
    equationNumber: number,
    tempK: number,
    coeffs: { [key: string]: number | null },
    Tc?: number | null
): number | null => {
    if (equationNumber === 105) {
        // DIPPR Equation for Liquid Density (kg/m^3)
        const { A, B, C } = coeffs;
        if (
            typeof A !== 'number' ||
            typeof B !== 'number' ||
            typeof C !== 'number' ||
            typeof Tc !== 'number'
        )
            return null;
        const Tr = tempK / Tc;
        if (Tr >= 1) return null; // Not valid at or above critical temp
        return A / Math.pow(B, 1 + Math.pow(1 - Tr, C));
    }
    console.warn(
        `Equation number ${equationNumber} is not implemented in this placeholder.`
    );
    return null;
};

// Utility function to safely parse coefficients
const parseCoefficient = (value: any): number | null => {
    if (typeof value === 'number') return value;
    if (typeof value === 'string') {
        const num = parseFloat(value);
        return isNaN(num) ? null : num;
    }
    return null;
};

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

    // --- ADD THESE NEW FIELDS ---
    isDbLocked?: boolean // To disable inputs after fetching
    liquidDensityData?: any // To store raw Supabase data for density
    criticalTemp?: number | null // To store critical temperature for density calc
        // --- ADD THESE NEW FIELDS ---
        phase?: 'Gas' | 'Liquid'
    // --- UPDATE: henryData now stores constants keyed by temperature row ID ---
    henryData?: { [tempId: string]: string } // Stores constants keyed by temperature row ID
    // --- ADD THIS FIELD ---
    cas_number?: string;
}

// --- NEW: Define types for the fluid package and parameters ---
type VleFluidPackage = 'Ideal' | 'NRTL' | 'Wilson' | 'UNIQUAC' | 'Peng-Robinson' | 'SRK';

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

// Helper function to perform linear interpolation (with extrapolation) from an array of data points
const interpolateFromData = (
    targetX: number,
    dataPoints: { x: number; y: number }[]
): number => {
    if (dataPoints.length === 0) {
        return NaN; // No data to interpolate from
    }
    if (dataPoints.length === 1) {
        return dataPoints[0].y; // Only one point, return its value
    }

    // Ensure data is sorted by x-value for reliable searching
    const sortedPoints = [...dataPoints].sort((a, b) => a.x - b.x);
    const n = sortedPoints.length;

    // If target is before the first point, return the first point's value
    if (targetX <= sortedPoints[0].x) {
        return sortedPoints[0].y;
    }
  
    // --- START: MODIFICATION ---
    // If the target is after the last point, return NaN instead of extrapolating.
    // This causes the plotting function to filter out the point and end the line naturally.
    if (targetX > sortedPoints[n - 1].x) {
        return NaN;
    }
    // --- END: MODIFICATION ---

    // Find the two points that bracket the targetX for interpolation
    let p1: { x: number; y: number } | null = null;
    let p2: { x: number; y: number } | null = null;

    for (let i = 0; i < sortedPoints.length - 1; i++) {
        if (sortedPoints[i].x <= targetX && sortedPoints[i + 1].x >= targetX) {
            p1 = sortedPoints[i];
            p2 = sortedPoints[i + 1];
            break;
        }
    }

    if (!p1 || !p2) {
        return NaN; // Should not happen with the new logic, but keep as a fallback
    }

    // Handle vertical line case to prevent division by zero
    if (p1.x === p2.x) {
        return p1.y;
    }

    // Standard linear interpolation formula
    return p1.y + ((targetX - p1.x) * (p2.y - p1.y)) / (p2.x - p1.x);
};


/**
 * Calculates a Henry's Law constant at a target temperature using interpolation
 * or extrapolation based on a provided set of data points.
 * Uses a physically meaningful linear relationship between ln(H) and 1/T (Kelvin).
 *
 * @param targetTemp_C The temperature (°C) at which to find the constant.
 * @param henryData The component's map of { tempId: constantValue }.
 * @param definedTemps The array of { id: tempId, temp: tempValue_C } from the UI.
 * @returns The calculated Henry's constant, or null if calculation is not possible.
 */
const getHenryConstant = (
    targetTemp_C: number,
    henryData: { [tempId: string]: string } | undefined,
    definedTemps: { id: string; temp: string }[]
): number | null => {
    if (!henryData) return null;

    const targetTemp_K = targetTemp_C + 273.15;

    // 1. Create a clean, sorted list of valid [1/T_K, ln(H)] data points
    const dataPoints = definedTemps
        .map(t => {
            const tempVal_C = parseFloat(t.temp);
            const henryVal = parseFloat(henryData[t.id]);
            if (!isNaN(tempVal_C) && !isNaN(henryVal) && henryVal > 0) {
                return {
                    inv_T: 1 / (tempVal_C + 273.15),
                    ln_H: Math.log(henryVal),
                };
            }
            return null;
        })
        .filter((p): p is { inv_T: number; ln_H: number } => p !== null)
        .sort((a, b) => a.inv_T - b.inv_T); // Sorts from high T to low T

    if (dataPoints.length < 2) {
        // Cannot interpolate or extrapolate without at least two points
        return dataPoints.length === 1 ? Math.exp(dataPoints[0].ln_H) : null;
    }

    const target_inv_T = 1 / targetTemp_K;

    // 2. Determine points for interpolation or extrapolation
    let p1: { inv_T: number; ln_H: number };
    let p2: { inv_T: number; ln_H: number };

    if (target_inv_T <= dataPoints[0].inv_T) {
        // Extrapolate above the highest temperature
        p1 = dataPoints[1];
        p2 = dataPoints[0];
    } else if (target_inv_T >= dataPoints[dataPoints.length - 1].inv_T) {
        // Extrapolate below the lowest temperature
        p1 = dataPoints[dataPoints.length - 1];
        p2 = dataPoints[dataPoints.length - 2];
    } else {
        // Interpolate within the temperature range
        const upperIndex = dataPoints.findIndex(p => p.inv_T > target_inv_T);
        p1 = dataPoints[upperIndex - 1];
        p2 = dataPoints[upperIndex];
    }
  
    // 3. Perform linear interpolation/extrapolation
    if (p1.inv_T === p2.inv_T) {
        // Avoid division by zero if temperatures are identical
        return Math.exp(p1.ln_H);
    }

    const slope = (p2.ln_H - p1.ln_H) / (p2.inv_T - p1.inv_T);
    const ln_H_target = p1.ln_H + slope * (target_inv_T - p1.inv_T);

    return Math.exp(ln_H_target);
};

const calculateCstrData = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    initialFlowRates: { [key: string]: number }, // Molar flows based on F_A0=1.0
    tempK: number,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    reactionPhase: 'Liquid' | 'Gas' | 'Mixed',
    pressure_bar: number,
    targetProduction_kta: number,
    // --- ADD henryUnit PARAM ---
    henryUnit: 'Pa' | 'kPa' | 'bar',
    // VLE parameters
    fluidPackage?: VleFluidPackage,
    vleComponentData?: Map<string, CompoundData>,
    interactionParameters?: Map<string, InteractionParams | null>,
    henryTemperatures?: { id: string; temp: string }[] // 👈 ADD THIS
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
            // Preserve original reaction metadata (AUnit, EaUnit, etc.) while adding parsed reactants/products
            return {
                ...reactionInfo,
                reactants,
                products,
                rateConstantAtT: A * Math.exp(-Ea / (R * tempK)),
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
    console.log(`[CSTR Tracer] Starting loop step ${i + 1}/${n_steps + 1}, V_basis = ${V_basis_liters.toExponential()}`);
        
        // 1. Solve the CSTR equations for the current basis volume
        const F_out_basis = solveCSTR_Robust(
            parsedReactions,
            components,
            V_basis_liters,
            F_in_basis,
            reactionPhase,
            pressure_bar,
            tempK,
            last_F_out_basis, // Pass the last result as the initial guess
            henryUnit,
            fluidPackage,
            vleComponentData,
            interactionParameters,
            henryTemperatures // 👈 ADD THIS
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

    // --- START: MODIFICATION ---
    const productStoichFactor = parseFloat(desiredProduct.reactionData[primaryReactionId]?.stoichiometry || '0');
    const limitingReactantStoichFactor = parseFloat(limitingReactant.reactionData[primaryReactionId]?.stoichiometry || '0');

    // The ratio of stoichiometric coefficients is needed for an accurate selectivity calculation.
    const stoichRatio = Math.abs(productStoichFactor / limitingReactantStoichFactor);

    const selectivity = molesConsumed_A > 1e-9 ? (molesFormed_P / molesConsumed_A) / stoichRatio : 0;
    // --- END: MODIFICATION ---
        
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
    const componentNames = allComponentNames;
    const results = solveODE_BDF_imported(dFdV, F0, V_span, componentNames);

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

// Local linear solver and BDF solver removed — using shared implementations from '@/lib/reactor-solver'

const solveCSTR_Robust = (
    parsedReactions: any,
    components: ComponentSetup[],
    V: number, // Liters
    initialFlowRates: { [key: string]: number },
    reactionPhase: 'Liquid' | 'Gas' | 'Mixed',
    totalPressure_bar: number,
    temp_K: number,
    initialGuess?: number[],
    henryUnit?: 'Pa' | 'kPa' | 'bar',
    fluidPackage?: VleFluidPackage,
    vleComponentData?: Map<string, CompoundData>,
    interactionParameters?: Map<string, InteractionParams | null>,
    henryTemperatures?: { id: string; temp: string }[]
): { [key: string]: number } => {
    const componentNames = parsedReactions.allInvolvedComponents.map((c: any) => c.name);
    let F = initialGuess ? [...initialGuess] : componentNames.map((name: string) => initialFlowRates[name] || 0);
    const n = F.length;
    const max_iter = 100;
    const tol = 1e-9;

    const G = (current_F: number[]): number[] => {
        const F_dict = componentNames.reduce((acc: { [key: string]: number }, name: string, i: number) => { acc[name] = current_F[i]; return acc; }, {} as { [key: string]: number });
        let rates: number[];
        const F_total = current_F.reduce((sum, f) => sum + Math.max(0, f), 0);
        if (F_total < 1e-12) return Array(n).fill(0).map((_, i) => (initialFlowRates[componentNames[i]] || 0) - current_F[i]);

        const rateBasis: { [key: string]: number } = {};
        const allCompSetupData = parsedReactions.allInvolvedComponents.map((c: any) => components.find(comp => comp.name === c.name)).filter(Boolean) as ComponentSetup[];

        if (reactionPhase === 'Mixed') {
            const liquidComps = allCompSetupData.filter(c => c.phase === 'Liquid');
            const gasComps = allCompSetupData.filter(c => c.phase === 'Gas');
            const liquidCompNames = liquidComps.map(c => c.name);

            const F_liquid_total = liquidCompNames.reduce((sum, name) => sum + (F_dict[name] || 0), 0);
            const liquidMoleFractions = F_liquid_total > 0 ? liquidCompNames.map(name => (F_dict[name] || 0) / F_liquid_total) : [];

            let liquidGammas = new Map(liquidCompNames.map(name => [name, 1.0]));
            if (fluidPackage !== 'Ideal' && liquidComps.length > 0 && interactionParameters) {
                const liquidVleData = liquidComps.map(c => ({ ...(vleComponentData?.get(c.id)), phase: c.phase })) as any[];
                
                // --- START OF FIX: Remap interaction parameter keys to local indices ---
                const originalIndexToLocalIndex = new Map<number, number>();
                liquidComps.forEach((comp, localIndex) => {
                    const originalIndex = allCompSetupData.findIndex(c => c.id === comp.id);
                    if (originalIndex !== -1) {
                        originalIndexToLocalIndex.set(originalIndex, localIndex);
                    }
                });

                const localInteractionParams = new Map<string, InteractionParams>();
                for (const [key, value] of interactionParameters.entries()) {
                    const [i_orig_str, j_orig_str] = key.split('-');
                    const i_orig = parseInt(i_orig_str, 10);
                    const j_orig = parseInt(j_orig_str, 10);

                    if (originalIndexToLocalIndex.has(i_orig) && originalIndexToLocalIndex.has(j_orig) && value) {
                        const i_local = originalIndexToLocalIndex.get(i_orig)!;
                        const j_local = originalIndexToLocalIndex.get(j_orig)!;
                        localInteractionParams.set(`${i_local}-${j_local}`, value);
                    }
                }
                // --- END OF FIX ---

                let calculatedGammas: number[] | null = null;
                switch (fluidPackage) {
                    case 'NRTL': calculatedGammas = calculateNrtlGammaMulticomponent(liquidVleData, liquidMoleFractions, temp_K, localInteractionParams as any); break;
                    case 'Wilson': calculatedGammas = calculateWilsonGammaMulticomponent(liquidVleData, liquidMoleFractions, temp_K, localInteractionParams as any); break;
                    case 'UNIQUAC': calculatedGammas = calculateUniquacGammaMulticomponent(liquidVleData, liquidMoleFractions, temp_K, localInteractionParams as any); break;
                }
                if (calculatedGammas) {
                    liquidCompNames.forEach((name, i) => liquidGammas.set(name, calculatedGammas[i]));
                }
            }

            let q_liquid_L_per_s = 0;
            liquidComps.forEach(comp => {
                const flow = F_dict[comp.name] || 0;
                if (flow > 0) {
                    q_liquid_L_per_s += (flow * parseFloat(comp.molarMass || '1')) / parseFloat(comp.density || '1000');
                }
            });
            if (q_liquid_L_per_s < 1e-12) q_liquid_L_per_s = 1e-12;
            const C_total_liquid_mol_L = F_liquid_total > 0 ? F_liquid_total / q_liquid_L_per_s : 0;

            const concentrations_L: { [key: string]: number } = {};
            liquidComps.forEach(comp => {
                concentrations_L[comp.name] = (F_dict[comp.name] || 0) / q_liquid_L_per_s;
                rateBasis[comp.name] = (liquidGammas.get(comp.name) || 1.0) * Math.max(0, concentrations_L[comp.name]);
            });

            const definedTemps = (window as any).__HENRY_TEMPS || [];
            gasComps.forEach(comp => {
                const partial_pressure_bar = (F_dict[comp.name] / F_total) * totalPressure_bar;
                let H_val = getHenryConstant(temp_K - 273.15, comp.henryData, definedTemps);
                if (H_val !== null) {
                    if (henryUnit === 'Pa') H_val /= 1e5; else if (henryUnit === 'kPa') H_val /= 100;
                }
                if (H_val && H_val > 0) {
                    concentrations_L[comp.name] = (partial_pressure_bar / H_val) * C_total_liquid_mol_L;
                } else {
                    concentrations_L[comp.name] = 0;
                }
                rateBasis[comp.name] = 1.0 * Math.max(0, concentrations_L[comp.name]);
            });
            rates = calculateNetReactionRates(parsedReactions.reactions, components, rateBasis, {}, temp_K);

        } else if (reactionPhase === 'Liquid') {
            const moleFractions = componentNames.map((name: string) => (F_dict[name] || 0) / F_total);
            let gammas = Array(n).fill(1.0);
            const allVleData = allCompSetupData.map(c => vleComponentData?.get(c.id)) as any[];
            if (fluidPackage !== 'Ideal' && allVleData.every(Boolean) && interactionParameters) {
                switch (fluidPackage) {
                    case 'NRTL': gammas = calculateNrtlGammaMulticomponent(allVleData, moleFractions, temp_K, interactionParameters as any) || gammas; break;
                    case 'Wilson': gammas = calculateWilsonGammaMulticomponent(allVleData, moleFractions, temp_K, interactionParameters as any) || gammas; break;
                    case 'UNIQUAC': gammas = calculateUniquacGammaMulticomponent(allVleData, moleFractions, temp_K, interactionParameters as any) || gammas; break;
                }
            }
            let v_Liters = 0;
            componentNames.forEach((name: string) => {
                const compInfo = components.find(c => c.name === name);
                v_Liters += (F_dict[name] * parseFloat(compInfo?.molarMass || '1')) / parseFloat(compInfo?.density || '1000');
            });
            if (v_Liters < 1e-9) v_Liters = 1e-9;
            componentNames.forEach((name: string, i: number) => {
                rateBasis[name] = gammas[i] * Math.max(0, F_dict[name]) / v_Liters;
            });
            rates = calculateNetReactionRates(parsedReactions.reactions, components, rateBasis, {}, temp_K);

        } else { // Gas Phase
            const moleFractions = componentNames.map((name: string) => (F_dict[name] || 0) / F_total);
            let phis = Array(n).fill(1.0);
            const allVleData = allCompSetupData.map(c => vleComponentData?.get(c.id)) as any[];
            if (fluidPackage !== 'Ideal' && allVleData.every(Boolean) && interactionParameters) {
                 switch (fluidPackage) {
                    case 'Peng-Robinson': phis = calculatePrFugacityCoefficientsMulticomponent(allVleData, moleFractions, temp_K, totalPressure_bar * 1e5, interactionParameters as any, 'vapor') || phis; break;
                    case 'SRK': phis = calculateSrkFugacityCoefficientsMulticomponent(allVleData, moleFractions, temp_K, totalPressure_bar * 1e5, interactionParameters as any, 'vapor') || phis; break;
                }
            }
            const P_Pa = convertUnit.pressure.barToPa(totalPressure_bar);
            componentNames.forEach((name: string, i: number) => {
                rateBasis[name] = phis[i] * moleFractions[i] * P_Pa;
            });
            rates = calculateNetReactionRates(parsedReactions.reactions, components, {}, rateBasis, temp_K);
        }

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

    // --- Newton-Raphson solver part is unchanged ---
    for (let iter = 0; iter < max_iter; iter++) {
        const G_current = G(F);
        const errorNorm = Math.sqrt(G_current.reduce((s, v) => s + v * v, 0));
        if (errorNorm < tol) break;
        const J_G: number[][] = Array(n).fill(0).map(() => Array(n).fill(0));
        for (let j = 0; j < n; j++) {
            const F_plus_h = [...F];
            const h = (Math.abs(F[j]) || 1e-8) * 1e-7;
            F_plus_h[j] += h;
            const G_plus_h = G(F_plus_h);
            for (let i = 0; i < n; i++) J_G[i][j] = (G_plus_h[i] - G_current[i]) / h;
        }
        const delta_F = solveLinearSystem_imported(J_G, G_current.map((val: number) => -val));
        if (!delta_F) break;
        let alpha = 1.0;
        const beta = 0.5;
        const c = 1e-4;
        let F_new = F;
        let step_accepted = false;
        for (let ls_iter = 0; ls_iter < 10; ls_iter++) {
            F_new = F.map((val: number, i: number) => val + alpha * delta_F[i]);
            F_new = F_new.map((v: number) => Math.max(0, v));
            const G_new = G(F_new);
            const newErrorNorm = Math.sqrt(G_new.reduce((s, v) => s + v * v, 0));
            if (newErrorNorm < (1 - c * alpha) * errorNorm) {
                step_accepted = true;
                break;
            }
            alpha *= beta;
        }
        if (step_accepted) { F = F_new; }
        else { break; }
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
        // Convert from hours to seconds if needed
        if (reactionInfo.AUnit.includes('h')) {
            A_f /= 3600.0; 
        }
        
        // --- Forward Rate Calculation ---
        const Ea_f_J = convertUnit.energy.convertEaToJules(parseFloat(reactionInfo.EaValue), reactionInfo.EaUnit);
        const k_f_c = A_f * Math.exp(-Ea_f_J / (R * tempK));
        
        const usePartialPressure = reactionInfo.rateLawBasis === 'partialPressure';
        let k_f: number;

        const AUnit = reactionInfo.AUnit || '';
        const AUnitIsPressureBased = /Pa|bar/.test(AUnit);
        if (usePartialPressure && !AUnitIsPressureBased) {
            let totalOrder = 0;
            components.forEach(comp => {
                const order = parseFloat(comp.reactionData[reactionInfo.id]?.order || '0');
                const stoich = parseFloat(comp.reactionData[reactionInfo.id]?.stoichiometry || '0');
                if (stoich < 0 && order > 0) totalOrder += order;
            });
            k_f = convertUnit.kinetics.kcToKp(k_f_c, tempK, totalOrder);
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

        // --- START: MODIFICATION ---
        // --- Reverse Rate Calculation (only if reaction is reversible) ---
        let rate_r = 0;
        if (reactionInfo.isEquilibrium) {
            let A_r = parseFloat(reactionInfo.AValueBackward || '0');
            if (reactionInfo.AUnitBackward?.includes('h')) {
                A_r /= 3600.0;
            }
            const Ea_r_J = convertUnit.energy.convertEaToJules(parseFloat(reactionInfo.EaValueBackward || '0'), reactionInfo.EaUnitBackward || 'J/mol');
            const k_r_c = A_r * Math.exp(-Ea_r_J / (R * tempK));
            
            let k_r: number;
            const AUnitBackwardIsPressureBased = /Pa|bar/.test(reactionInfo.AUnitBackward || '');
            if (usePartialPressure && !AUnitBackwardIsPressureBased) {
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
        // --- END: MODIFICATION ---

        return rate_f - rate_r; // Returns the net rate for the reaction
    });
};

const calculatePfrData_Optimized = (
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    tempK: number,
    pressure_bar: number,
    molarRatios: Array<{ numeratorId: string; value: number }>,
    simBasis: { limitingReactantId: string; desiredProductId: string },
    targetProduction_kta: number,
    reactionPhase: 'Liquid' | 'Gas',
    // --- ADD THESE NEW PARAMETERS ---
    fluidPackage?: VleFluidPackage,
    vleComponentData?: Map<string, CompoundData>,
    interactionParameters?: Map<string, InteractionParams | null>
) => {
    // Add this log at the top
    console.log(`[DEBUG-PFR] 🧠 PFR solver using package: ${fluidPackage}. Received interaction params:`, interactionParameters);

    // 1. Initial Setup
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);
    if (!limitingReactant || !desiredProduct) {
        return { selectivityData: [], volumeData: [] };
    }
    const allComponentNames = components.map(c => c.name);
    const limitingReactantIndex = components.findIndex(c => c.id === simBasis.limitingReactantId);
    const desiredProductIndex = components.findIndex(c => c.id === simBasis.desiredProductId);

    const primaryReactionId = reactions[0].id;
    const limitingReactantStoichInReaction = Math.abs(parseFloat(limitingReactant.reactionData[primaryReactionId]?.stoichiometry || '1'));
    const desiredProductStoichInReaction = Math.abs(parseFloat(desiredProduct.reactionData[primaryReactionId]?.stoichiometry || '1'));
    // The stoichiometric ratio is crucial for accurate selectivity
    const stoichRatio = desiredProductStoichInReaction / limitingReactantStoichInReaction;

    // 2. Pre-calculate all temperature-dependent constants ONCE for performance.
    const parsedReactions = reactions.map(reactionInfo => {
        const Ea_f_J = convertUnit.energy.convertEaToJules(parseFloat(reactionInfo.EaValue), reactionInfo.EaUnit);
        let A_f = parseFloat(reactionInfo.AValue);
        if (reactionInfo.AUnit.includes('h')) A_f /= 3600.0;
        const k_f = A_f * Math.exp(-Ea_f_J / (R * tempK));

        let k_r = 0;
        if (reactionInfo.isEquilibrium) {
            const Ea_r_J = convertUnit.energy.convertEaToJules(parseFloat(reactionInfo.EaValueBackward || '0'), reactionInfo.EaUnitBackward || 'J/mol');
            let A_r = parseFloat(reactionInfo.AValueBackward || '0');
            if (reactionInfo.AUnitBackward?.includes('h')) A_r /= 3600.0;
            k_r = A_r * Math.exp(-Ea_r_J / (R * tempK));
        }

        const participants = allComponentNames.map(name => {
            const comp = components.find(c => c.name === name);
            const reactionData = comp?.reactionData[reactionInfo.id];
            return {
                name: name,
                stoichiometry: parseFloat(reactionData?.stoichiometry || '0'),
                order: parseFloat(reactionData?.order || '0'),
                orderReverse: parseFloat(reactionData?.orderReverse || '0'),
            };
        });

        return {
            rateLawBasis: reactionInfo.rateLawBasis,
            isEquilibrium: reactionInfo.isEquilibrium,
            k_f: k_f,
            k_r: k_r,
            participants: participants,
        };
    });

    // 3. Define a ROBUST Rate Law Function (dF/dV = r)
    const dFdV = (V: number, F_vec: number[]): number[] => {
        // This object will hold the basis for the rate law (activity or fugacity)
        const rateBasis: { [key: string]: number } = {};
        const F_total = F_vec.reduce((sum, f) => sum + Math.max(0, f), 0);
        if (F_total < 1e-12) return Array(allComponentNames.length).fill(0);

        const moleFractions = F_vec.map(f => Math.max(0, f) / F_total);

        if (reactionPhase === 'Liquid') {
            // --- Activity Coefficient Models ---
            let gammas: number[] = Array(allComponentNames.length).fill(1.0); // Default to ideal

            // Get component data for multicomponent calculations
            const allCompData = components.map(c => {
                const data = vleComponentData?.get(c.id);
                if (!data) return null;
                return { ...data, phase: c.phase }; // Merge phase property
            }).filter(Boolean) as (CompoundData & { phase?: 'Liquid' | 'Gas' })[];
            
            if (fluidPackage !== 'Ideal' && allCompData.length === allComponentNames.length && interactionParameters) {
                switch (fluidPackage) {
                    case 'NRTL':
                        gammas = calculateNrtlGammaMulticomponent(allCompData, moleFractions, tempK, interactionParameters as unknown as NrtlParameterMatrix) || gammas;
                        // Add this log to see the result
                        console.log('[DEBUG-PFR] NRTL gammas:', gammas);
                        break;
                    case 'Wilson':
                        gammas = calculateWilsonGammaMulticomponent(allCompData, moleFractions, tempK, interactionParameters as unknown as WilsonParameterMatrix) || gammas;
                        // Add this log to see the result
                        console.log('[DEBUG-PFR] Wilson gammas:', gammas);
                        break;
                    case 'UNIQUAC':
                        gammas = calculateUniquacGammaMulticomponent(allCompData, moleFractions, tempK, interactionParameters as unknown as UniquacParameterMatrix) || gammas;
                        // Add this log to see the result
                        console.log('[DEBUG-PFR] UNIQUAC gammas:', gammas);
                        break;
                }
            }

            // Calculate liquid volumetric flow rate and concentrations
            let q_Liters = 0;
            allComponentNames.forEach((name, i) => {
                const flow = F_vec[i];
                if (flow > 0) {
                    const comp = components[i];
                    q_Liters += (flow * parseFloat(comp.molarMass || '1')) / parseFloat(comp.density || '1000');
                }
            });
            if (q_Liters < 1e-9) q_Liters = 1e-9;
            
            // Calculate activity = gamma * concentration
            allComponentNames.forEach((name, i) => {
                const concentration = Math.max(0, F_vec[i]) / q_Liters;
                rateBasis[name] = gammas[i] * concentration;
            });

        } else { // Gas Phase
            // --- Fugacity Coefficient Models (Equations of State) ---
            let phis: number[] = Array(allComponentNames.length).fill(1.0); // Default to ideal

            // Get component data for multicomponent calculations
            const allCompData = components.map(c => vleComponentData?.get(c.id)).filter(Boolean) as CompoundData[];

            if (fluidPackage !== 'Ideal' && allCompData.length === allComponentNames.length && interactionParameters) {
                switch (fluidPackage) {
                    case 'Peng-Robinson':
                        phis = calculatePrFugacityCoefficientsMulticomponent(allCompData, moleFractions, tempK, pressure_bar * 1e5, interactionParameters as unknown as PrSrkParameterMatrix, 'vapor') || phis;
                        // Add this log to see the result
                        console.log('[DEBUG-PFR] Peng-Robinson phis:', phis);
                        break;
                    case 'SRK':
                        phis = calculateSrkFugacityCoefficientsMulticomponent(allCompData, moleFractions, tempK, pressure_bar * 1e5, interactionParameters as unknown as PrSrkParameterMatrix, 'vapor') || phis;
                        // Add this log to see the result
                        console.log('[DEBUG-PFR] SRK phis:', phis);
                        break;
                }
            }

            // Calculate fugacity = phi * partial pressure
            const P_Pa = convertUnit.pressure.barToPa(pressure_bar);
            allComponentNames.forEach((name, i) => {
                const partialPressure = moleFractions[i] * P_Pa;
                rateBasis[name] = phis[i] * partialPressure;
            });
        }
        
        // --- Calculate net rates using the new rateBasis (activity or fugacity) ---
        const reactionNetRates = parsedReactions.map(pReaction => {
            let rate_f = pReaction.k_f;
            let rate_r = pReaction.k_r;
            // Use pressure-based rate constants if the basis is fugacity
            const usePartialPressure = pReaction.rateLawBasis === 'partialPressure' || reactionPhase === 'Gas';

            pReaction.participants.forEach(p => {
                if (p.stoichiometry < 0 && p.order > 0) {
                    const basisValue = (reactionPhase === 'Liquid') ? rateBasis[p.name] : rateBasis[p.name];
                    rate_f *= Math.pow(Math.max(0, basisValue), p.order);
                }
                if (p.stoichiometry > 0 && p.orderReverse > 0) {
                    const basisValue = (reactionPhase === 'Liquid') ? rateBasis[p.name] : rateBasis[p.name];
                    rate_r *= Math.pow(Math.max(0, basisValue), p.orderReverse);
                }
            });
            return rate_f - (pReaction.isEquilibrium ? rate_r : 0);
        });

        const dF_vec = Array(allComponentNames.length).fill(0);
        parsedReactions.forEach((pReaction, i) => {
            pReaction.participants.forEach((p, j) => {
                dF_vec[j] += p.stoichiometry * reactionNetRates[i];
            });
        });

        // The ODE solver expects dF/dV in mol/L/s for liquid and mol/m³/s for gas.
        // Our rate laws are already in mol/L/s (for liquid) or based on Pa (for gas).
        // The gas phase rate calculation needs a conversion factor.
        if (reactionPhase === 'Gas') {
             for (let i = 0; i < dF_vec.length; i++) dF_vec[i] *= 1000;
        }

        return dF_vec;
    };

    // 4. Set up and run the ODE solver ONCE
    const F0_basis_dict: { [key: string]: number } = { [limitingReactant.name]: 1.0 };
    molarRatios.forEach(ratio => {
        const compName = components.find(c => c.id === ratio.numeratorId)?.name;
        if (compName) F0_basis_dict[compName] = ratio.value;
    });
    const F0_basis_vec = components.map(c => F0_basis_dict[c.name] || 0);
    const F_A_in_basis = F0_basis_vec[limitingReactantIndex];
    if (F_A_in_basis <= 0) return { selectivityData: [], volumeData: [] };

    const V_max_basis = 1e8; // A large basis volume to ensure high conversion
    const V_span: [number, number] = [0, V_max_basis];

    const rawSolution = solveODE_BDF_imported(dFdV, F0_basis_vec, V_span, allComponentNames);

    // 5. Process the high-quality solver output
    const rawSelectivityData: { x: number; y: number }[] = [];
    const rawVolumeData: { x: number; y: number }[] = [];

    const productMolarMass = parseFloat(desiredProduct.molarMass || '1.0');
    const P_target_mols = convertUnit.flowRate.ktaToMolPerSecond(targetProduction_kta, productMolarMass);
    const F_P_in_basis = F0_basis_vec[desiredProductIndex] || 0;

    for (let i = 0; i < rawSolution.t.length; i++) {
        const V_basis_current = rawSolution.t[i];
        const F_A_out_basis = rawSolution.y[limitingReactantIndex][i];
        const F_P_out_basis = rawSolution.y[desiredProductIndex][i];

        const conversion = (F_A_in_basis - F_A_out_basis) / F_A_in_basis;
        const molesReactantConsumed = F_A_in_basis - F_A_out_basis;
        const netMolesProductFormed = F_P_out_basis - F_P_in_basis;
        
        const selectivity = (molesReactantConsumed > 1e-9) ? (netMolesProductFormed / molesReactantConsumed) / stoichRatio : 0;
        
        let V_true_m3 = Infinity;

        if (selectivity > 1e-9 && conversion > 1e-9) {
            const F_A0_true = P_target_mols / (selectivity * conversion);
            const scaling_factor = F_A0_true / F_A_in_basis;
            const V_scaled = V_basis_current * scaling_factor;
            V_true_m3 = (reactionPhase === 'Liquid') ? convertUnit.volume.litersToCubicMeters(V_scaled) : V_scaled;
        }

        if (conversion >= 1e-6 && conversion < 1.0 && V_true_m3 < Infinity && selectivity >= 0) {
            rawSelectivityData.push({ x: conversion, y: selectivity });
            rawVolumeData.push({ x: conversion, y: V_true_m3 });
        }
    }

    // 6. Interpolate results for a perfectly smooth plot
    const commonConversionGrid = Array.from({ length: 101 }, (_, i) => i * 0.01);
    const selectivityData = commonConversionGrid.map(conv => ({
        x: conv, y: interpolateFromData(conv, rawSelectivityData)
    })).filter(p => !isNaN(p.y) && p.y >= 0);
    const volumeData = commonConversionGrid.map(conv => ({
        x: conv, y: interpolateFromData(conv, rawVolumeData)
    })).filter(p => !isNaN(p.y) && p.y >= 0);

    return { selectivityData, volumeData };
};const calculateNetFormationRates = (
    current_F: { [key: string]: number },
    reactions: ReactionSetup[],
    components: ComponentSetup[],
    tempK: number,
    reactionPhase: 'Liquid' | 'Gas' | 'Mixed',
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
  reactionsSetup,
  setReactionsSetup,
  componentsSetup,
  setComponentsSetup,
  simBasis,
  setSimBasis,
  reactionPhase,
  setReactionPhase,
  prodRate,
  setProdRate,
  setTemperature,
  setPressure,
    // ▼▼▼ ADD THIS PROP ▼▼▼
    setReactorTypes,
    // ▲▲▲ ADD THIS PROP ▲▲▲
  fetchComponentProperties,
  isFetching,
    fluidPackage,
    setFluidPackage,
    // Henry props passed from parent
    henryTemperatures,
    setHenryTemperatures,
    henryUnit,
    setHenryUnit,
}: {
  onNext: () => void
  reactionsSetup: ReactionSetup[]
  setReactionsSetup: React.Dispatch<React.SetStateAction<ReactionSetup[]>>
  componentsSetup: ComponentSetup[]
  setComponentsSetup: React.Dispatch<React.SetStateAction<ComponentSetup[]>>
  simBasis: any
  setSimBasis: any
    reactionPhase: 'Liquid' | 'Gas' | 'Mixed'
    setReactionPhase: React.Dispatch<React.SetStateAction<'Liquid' | 'Gas' | 'Mixed'>>
  prodRate: string
  setProdRate: React.Dispatch<React.SetStateAction<string>>
  setTemperature: React.Dispatch<React.SetStateAction<number>>
  setPressure: React.Dispatch<React.SetStateAction<number>>
    // ▼▼▼ ADD THIS PROP TYPE ▼▼▼
    setReactorTypes: React.Dispatch<React.SetStateAction<{ pfr: boolean; cstr: boolean }>>
    // ▲▲▲ ADD THIS PROP TYPE ▲▲▲
  fetchComponentProperties: (
    name: string,
    componentId: string
  ) => Promise<void>
  isFetching: boolean
    fluidPackage: VleFluidPackage;
    setFluidPackage: React.Dispatch<React.SetStateAction<VleFluidPackage>>;
    // Henry's Law table state (managed by parent)
    henryTemperatures: { id: string; temp: string }[];
    setHenryTemperatures: React.Dispatch<React.SetStateAction<{ id: string; temp: string }[]>>;
    henryUnit: 'Pa' | 'kPa' | 'bar';
    setHenryUnit: React.Dispatch<React.SetStateAction<'Pa' | 'kPa' | 'bar'>>;
    // interaction parameters are managed at the parent level now
}) => {
  const { resolvedTheme } = useTheme()

  const [suggestions, setSuggestions] = useState<string[]>([])
  const [showSuggestions, setShowSuggestions] = useState(false)
  const [activeInputId, setActiveInputId] = useState<string | null>(null)
  const containerRef = useRef<HTMLDivElement>(null)

  const handleFetchSuggestions = useCallback(async (value: string) => {
    if (!value || value.length < 2 || !supabase) {
      setSuggestions([])
      setShowSuggestions(false)
      return
    }
    try {
      const { data, error } = await supabase
        .from('compound_properties')
        .select('name')
        .ilike('name', `${value}%`)
        .limit(5)
      if (error) throw error
      const fetchedNames = data?.map(d => d.name) || []
      setSuggestions(fetchedNames)
      setShowSuggestions(fetchedNames.length > 0)
    } catch (err) {
      console.error('Suggestion fetch error:', err)
    }
  }, [])

  const handleNameChange = (
    e: React.ChangeEvent<HTMLInputElement>,
    id: string
  ) => {
    const value = e.target.value
    setComponentsSetup(prev =>
      prev.map(c => (c.id === id ? { ...c, name: value } : c))
    )
    setActiveInputId(id)
    handleFetchSuggestions(value)
  }

  const handleSuggestionClick = (name: string) => {
    if (activeInputId) {
      fetchComponentProperties(name, activeInputId)
    }
    setShowSuggestions(false)
    setActiveInputId(null)
  }

  const handleUnlockComponent = (id: string) => {
    setComponentsSetup(prev =>
      prev.map(c =>
        c.id === id
          ? {
              ...c,
              isDbLocked: false,
              liquidDensityData: undefined,
              criticalTemp: undefined,
            }
          : c
      )
    )
  }

  useEffect(() => {
    function handleClickOutside(event: MouseEvent) {
      if (
        containerRef.current &&
        !containerRef.current.contains(event.target as Node)
      ) {
        setShowSuggestions(false)
        setActiveInputId(null)
      }
    }
    document.addEventListener('mousedown', handleClickOutside)
    return () => document.removeEventListener('mousedown', handleClickOutside)
  }, [])

    // --- NEW: Generic component property change helper ---
    const handleComponentPropertyChange = (
        id: string,
        field: keyof ComponentSetup,
        value: any
    ) => {
        setComponentsSetup(prev =>
            prev.map(c => (c.id === id ? { ...c, [field]: value } : c))
        )
    }

        // --- NEW: Generate unique component pairs for the new card ---
        const componentPairs = useMemo(() => {
            const pairs: [ComponentSetup, ComponentSetup][] = [];
            for (let i = 0; i < componentsSetup.length; i++) {
                for (let j = i + 1; j < componentsSetup.length; j++) {
                    pairs.push([componentsSetup[i], componentsSetup[j]]);
                }
            }
            return pairs;
        }, [componentsSetup]);

    // Interaction parameter fetching is handled at the parent level now.

  // ... (rest of the component logic for adding/removing reactions, presets, etc.)

                // --- FIX #2: Corrected grid columns with fixed widths for consistent margins ---
                const reactionCols = `repeat(${reactionsSetup.length * 2}, minmax(0, 1fr))`;
                const componentsGridCols =
                    reactionPhase === 'Gas'
                        ? `40px minmax(0, 1.5fr) 120px ${reactionCols}`
                        : reactionPhase === 'Liquid'
                        ? `40px minmax(0, 1.5fr) 120px 120px ${reactionCols}`
                        : `40px minmax(0, 1.5fr) 120px 110px ${reactionCols}`; // Mixed

        // --- CHANGE #1: Logic to determine if Henry's Law card should be visible ---
        const hasGasComponent = useMemo(
            () =>
                reactionPhase === 'Mixed' &&
                componentsSetup.some(c => c.phase === 'Gas'),
            [reactionPhase, componentsSetup]
        );

  // Unit mapping function to convert between Gas and Liquid phase units
  const mapUnits = (currentUnit: string, toPhase: 'Gas' | 'Liquid'): string => {
    if (toPhase === 'Gas') {
      // Converting from Liquid to Gas - add pressure units (default to Pa)
      const liquidToGasMap: { [key: string]: string } = {
        'mol,L,s': 'mol,L,Pa,s',
        'mol,L,h': 'mol,L,Pa,h',
        'mol,m^3,s': 'mol,m^3,Pa,s',
        'mol,m^3,h': 'mol,m^3,Pa,h'
      };
      return liquidToGasMap[currentUnit] || currentUnit;
    } else {
      // Converting from Gas to Liquid - remove pressure units
      const gasToLiquidMap: { [key: string]: string } = {
        'mol,L,Pa,s': 'mol,L,s',
        'mol,L,Pa,h': 'mol,L,h', 
        'mol,m^3,Pa,s': 'mol,m^3,s',
        'mol,m^3,Pa,h': 'mol,m^3,h',
        'mol,L,bar,s': 'mol,L,s',
        'mol,L,bar,h': 'mol,L,h',
        'mol,m^3,bar,s': 'mol,m^3,s',
        'mol,m^3,bar,h': 'mol,m^3,h'
      };
      return gasToLiquidMap[currentUnit] || currentUnit;
    }
  };

    // Parameter UI moved to parent; this component is simplified.

    // Enhanced phase switching function with automatic unit mapping
    const handlePhaseChange = (newPhase: 'Liquid' | 'Gas' | 'Mixed') => {
        if (newPhase === reactionPhase) return // No change needed

        // For Mixed phase, don't remap units; rely on per-component phase selection
        if (newPhase === 'Mixed') {
            setReactionPhase(newPhase)
            return
        }

        // Map units for all reactions when toggling between Liquid and Gas
        const updatedReactions = reactionsSetup.map(reaction => ({
            ...reaction,
            AUnit: mapUnits(reaction.AUnit, newPhase),
            AUnitBackward: reaction.AUnitBackward
                ? mapUnits(reaction.AUnitBackward, newPhase)
                : reaction.AUnitBackward,
        }))

        setReactionsSetup(updatedReactions)
        setReactionPhase(newPhase)
    }
  
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
        // --- START: MODIFICATION ---
        // First set up the reaction data with the correct Gas phase units
        const hdaReactions: ReactionSetup[] = [
                // HDA Rxn 1 (Gas phase, Concentration basis -> simple units)
                // Python units were sqrt(L/gmol)/h. The closest available unit is mol,L,h. The value is correct for this time basis.
                { id: '1', AValue: '1.188e14', AUnit: 'mol,L,Pa,h', EaValue: '52.0', EaUnit: 'kcal/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
                // HDA Rxn 2 (Gas phase, Pressure basis -> pressure units)
                // Python units were gmol/(L*Pa^2*h). The closest available unit is mol,L,Pa,h. The rate law uses the correct squared order.
                { id: '2', AValue: '1.1647e-5', AUnit: 'mol,L,Pa,h', EaValue: '30.19', EaUnit: 'kcal/mol', isEquilibrium: true, AValueBackward: '4.56e-7', AUnitBackward: 'mol,L,Pa,h', EaValueBackward: '23.35', EaUnitBackward: 'kcal/mol', rateLawBasis: 'partialPressure'}
        ];
        
        setReactionsSetup(hdaReactions);
        // --- END: MODIFICATION ---
        setComponentsSetup([
                { id: 'comp-T', name: 'Toluene', molarMass: '92.14', density: '867', reactionData: { '1': { stoichiometry: '-1', order: '1' } } },
                { id: 'comp-H2', name: 'Hydrogen', molarMass: '2.016', density: '89.9', reactionData: { '1': { stoichiometry: '-1', order: '0.5' }, '2': { stoichiometry: '1', order: '0', orderReverse: '1' } } },
                { id: 'comp-B', name: 'Benzene', molarMass: '78.11', density: '876', reactionData: { '1': { stoichiometry: '1', order: '0' }, '2': { stoichiometry: '-2', order: '2' } } },
                { id: 'comp-CH4', name: 'Methane', molarMass: '16.04', density: '717', reactionData: { '1': { stoichiometry: '1', order: '0' } } },
                { id: 'comp-D', name: 'Diphenyl', molarMass: '154.21', density: '1080', reactionData: { '2': { stoichiometry: '1', order: '0', orderReverse: '1' } } },
        ]);
        setSimBasis({ limitingReactantId: 'comp-T', desiredProductId: 'comp-B' });
        setReactionPhase('Gas'); // Set phase directly since we already have the correct units
        setProdRate('100');
        setTemperature(650);
        setPressure(30);
        // ▼▼▼ ADD THIS LINE ▼▼▼
        setReactorTypes({ pfr: true, cstr: false }); // HDA defaults to PFR
    };

  const loadAlkylationPreset = () => {
    setReactionsSetup(initialReactions);
    setComponentsSetup(initialComponents);
    setSimBasis({ limitingReactantId: 'comp-A', desiredProductId: 'comp-C' });
    setReactionPhase('Liquid'); // Set phase directly since initialReactions already has correct Liquid units
    setProdRate('100');
    setTemperature(4); // MODIFICATION: Corrected from 10 to 4
    setPressure(10);
    // ▼▼▼ ADD THIS LINE ▼▼▼
    setReactorTypes({ pfr: false, cstr: true }); // Alkylation defaults to CSTR
  }
  
                // --- DMC Preset Loader with Default Henry's Law Data ---
    const loadDmcPreset = async () => {
        if (!supabase) {
            console.error('Supabase client not available.');
            return;
        }

        // --- Data from the provided image ---
        const presetTemperatures = ['80', '90', '100', '110', '120', '130'];
        const presetHenryConstants = {
            'Oxygen': ['1526', '1451', '1379', '1303', '1225', '1146'],
            'Carbon monoxide': ['4265', '3843', '3484', '3156', '2792', '2496'],
            'Carbon dioxide': ['206', '222', '237', '250', '262', '271'],
        };
        // ------------------------------------

        const componentNames = [
            'Methanol',
            'Carbon monoxide',
            'Oxygen',
            'Dimethyl carbonate',
            'Water',
            'Carbon dioxide',
        ];
        const gasPhaseNames = ['Carbon monoxide', 'Oxygen', 'Carbon dioxide'];

        // 1. Create the temperature rows for the UI
        const newHenryTemperatures = presetTemperatures.map(t => ({
            id: `temp-dmc-${t}`, // Give each row a unique, stable ID
            temp: t,
        }));

            // 2. Fetch all component data from Supabase
            const fetchPromises = componentNames.map(async name => {
                try {
                    const { data: compoundDbData, error: compoundError } = await supabase
                        .from('compound_properties')
                        .select('id, name, properties')
                        .ilike('name', name)
                        .limit(1)
                        .single();
                    if (compoundError) throw compoundError;

                    if (!compoundDbData || !compoundDbData.properties) {
                        throw new Error(`No properties found for compound '${name}'`);
                    }

                    // Add debugging
                    console.log('DMC preset - fetched compound:', name, compoundDbData.properties);
                    console.log('DMC preset - available keys:', Object.keys(compoundDbData.properties || {}));

                    // Extract molecular weight from properties JSON - handle if it's an object
                    let molecularWeight = compoundDbData.properties.MolecularWeight || compoundDbData.properties.molecular_weight;
                    
                    // Try other possible keys if the standard ones don't work
                    if (!molecularWeight) {
                      const possibleKeys = ['MW', 'mw', 'MolarMass', 'molar_mass', 'mol_weight'];
                      for (const key of possibleKeys) {
                        if (compoundDbData.properties[key]) {
                          molecularWeight = compoundDbData.properties[key];
                          break;
                        }
                      }
                    }
                    
                    if (typeof molecularWeight === 'object' && molecularWeight !== null) {
                      molecularWeight = molecularWeight.value || molecularWeight.Value || molecularWeight.val || Object.values(molecularWeight)[0];
                    }
                    
                    console.log('DMC preset - extracted molecular weight:', molecularWeight);

                    return {
                        id: `comp-${compoundDbData.name.replace(/\s/g, '')}`,
                        name: compoundDbData.name,
                        molarMass: String(molecularWeight),
                        density: '',
                        reactionData: {},
                        isDbLocked: true,
                        liquidDensityData: compoundDbData.properties['LiquidDensity'] || compoundDbData.properties['Liquid density'],
                        criticalTemp: parseCoefficient(compoundDbData.properties['CriticalTemperature'] || compoundDbData.properties['Critical temperature']),
                        phase: gasPhaseNames.includes(compoundDbData.name) ? 'Gas' : 'Liquid',
                    } as ComponentSetup;
                } catch (error) {
                    console.error(`Failed to fetch DMC component: '${name}'. Check database for this entry.`, error);
                    return null;
                }
            });

            const fetchedRaw = await Promise.all(fetchPromises);
            const fetchedComponents = fetchedRaw.filter(Boolean) as ComponentSetup[];

            if (fetchedComponents.length !== 6) {
                console.error('Could not fetch all DMC components. Check the log above for details.');
                return;
            }

            // 3. Add reaction and Henry's Law data to the fetched components
            const finalComponents = fetchedComponents.map(c => {
                switch (c.name) {
                    case 'Methanol':
                        c.reactionData = { '1': { stoichiometry: '-2', order: '2' } };
                        break;
                    case 'Carbon monoxide':
                        c.reactionData = {
                            '1': { stoichiometry: '-1', order: '0' },
                            '2': { stoichiometry: '-1', order: '0' },
                        };
                        break;
                    case 'Oxygen':
                        c.reactionData = {
                            '1': { stoichiometry: '-0.5', order: '0.5' },
                            '2': { stoichiometry: '-0.5', order: '0.5' },
                        };
                        break;
                    case 'Dimethyl carbonate':
                        c.id = 'comp-DMC';
                        c.reactionData = { '1': { stoichiometry: '1', order: '0' } };
                        break;
                    case 'Water':
                        c.reactionData = { '1': { stoichiometry: '1', order: '0' } };
                        break;
                    case 'Carbon dioxide':
                        c.reactionData = { '2': { stoichiometry: '1', order: '0' } };
                        break;
                }

                // Assign Henry's Law data if the component is a gas in the preset
                if (c.name in presetHenryConstants) {
                    const constants = presetHenryConstants[c.name as keyof typeof presetHenryConstants];
                    const newHenryData: { [key: string]: string } = {};
                    newHenryTemperatures.forEach((tempRow, index) => {
                        newHenryData[tempRow.id] = constants[index];
                    });
                    c.henryData = newHenryData;
                }
                return c;
            });


        // 4. Set all the states for the preset
        const dmcReactions: ReactionSetup[] = [
            { id: '1', AValue: '1.6e11', AUnit: 'mol,L,s', EaValue: '104600', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
            { id: '2', AValue: '6.7e12', AUnit: 'mol,L,s', EaValue: '98750', EaUnit: 'J/mol', isEquilibrium: false, rateLawBasis: 'concentration' },
        ];
    
        setReactionsSetup(dmcReactions);
        setHenryTemperatures(newHenryTemperatures); // Set the temperature rows
        setComponentsSetup(finalComponents);       // Set components with populated data
        setSimBasis({ limitingReactantId: 'comp-Oxygen', desiredProductId: 'comp-DMC' });
        setReactionPhase('Mixed');
        setHenryUnit('bar');
    // --- FIX: Ensure fluid package is set correctly ---
    setFluidPackage('NRTL');
        setProdRate('50');
        setTemperature(80);
        setPressure(25);
        // ▼▼▼ ADD THIS LINE ▼▼▼
        setReactorTypes({ pfr: false, cstr: true }); // DMC defaults to CSTR
    };

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
            let hasReactant = false
            let hasProduct = false

            for (const comp of componentsSetup) {
                const stoichValue = parseFloat(
                    comp.reactionData[reaction.id]?.stoichiometry || '0'
                )
                if (stoichValue < 0) {
                    hasReactant = true
                } else if (stoichValue > 0) {
                    hasProduct = true
                }
            }

            if (!hasReactant || !hasProduct) {
                return {
                    isValid: false,
                    message: `Reaction ${reaction.id} must have at least one reactant (negative stoichiometry) and one product (positive stoichiometry).`,
                }
            }
        }

        // --- NEW: Mixed phase requires all components from DB (locked) ---
        if (reactionPhase === 'Mixed') {
            const allComponentsLocked = componentsSetup.every(c => c.isDbLocked)
            if (!allComponentsLocked) {
                return {
                    isValid: false,
                    message:
                        'For Mixed Phase, all components must be selected from the database.',
                }
            }
        }

        return { isValid: true, message: '' }
    }, [reactionsSetup, componentsSetup, reactionPhase])
  
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

  // --- NEW: State to manage the temperature rows for Henry's Law ---
    // Henry's Law temps & unit are passed in from parent (Home) via props

  // --- Updated Handlers ---
  const handleAddTemperatureRow = () => {
    setHenryTemperatures([
      ...henryTemperatures,
      { id: `temp-${Date.now()}`, temp: '' },
    ]);
  };
  
  // --- Updated to prevent deleting the first two rows ---
  const handleRemoveTemperatureRow = (idToRemove: string) => {
    if (henryTemperatures.length <= 2) return; // Keep at least two rows
    setHenryTemperatures(henryTemperatures.filter(row => row.id !== idToRemove));
    // Also clean up the orphaned data from components
    setComponentsSetup(prev =>
      prev.map(comp => {
        if (comp.henryData && comp.henryData[idToRemove]) {
          const newHenryData = { ...comp.henryData };
          delete newHenryData[idToRemove];
          return { ...comp, henryData: newHenryData };
        }
        return comp;
      })
    );
  };

  const handleTemperatureChange = (id: string, value: string) => {
    setHenryTemperatures(
      henryTemperatures.map(row => (row.id === id ? { ...row, temp: value } : row))
    );
  };

  const handleHenryConstantChange = (
    componentId: string,
    tempId: string,
    value: string
  ) => {
    setComponentsSetup(prev =>
      prev.map(comp => {
        if (comp.id === componentId) {
          return {
            ...comp,
            henryData: {
              ...comp.henryData,
              [tempId]: value,
            },
          };
        }
        return comp;
      })
    );
  };

  // Memoize the list of gas components to avoid re-filtering on every render
  const gasComponents = useMemo(
    () => componentsSetup.filter(c => c.phase === 'Gas'),
    [componentsSetup]
  );

  // --- Updated grid definition for new trash can position ---
  const henrysLawGridCols = `40px 1fr ${'1fr '.repeat(gasComponents.length)}`;

  return (
    <TooltipProvider>
    <div className={`min-h-screen flex flex-col p-4 bg-background text-foreground`}>
      <div className="container mx-auto space-y-6">
        
        {/* --- NEW HEADER ROW --- */}
        <div className="relative flex justify-center items-center py-2 mb-4">
            {/* Centered Preset Buttons */}
            <div className="flex justify-center items-center gap-4">
                <Label className="text-sm font-semibold text-muted-foreground">Load a preset from my design courses:</Label>
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
                                                        onClick={() => handlePhaseChange('Liquid')}
                                                        variant={reactionPhase === 'Liquid' ? 'default' : 'ghost'}
                                                        size="sm"
                                                        className="flex-1 text-xs px-3 py-1 h-auto"
                                                >
                                                        Liquid
                                                </Button>
                                                <Button
                                                        onClick={() => handlePhaseChange('Gas')}
                                                        variant={reactionPhase === 'Gas' ? 'default' : 'ghost'}
                                                        size="sm"
                                                        className="flex-1 text-xs px-3 py-1 h-auto"
                                                >
                                                        Gas
                                                </Button>
                                                <Button
                                                        onClick={() => handlePhaseChange('Mixed')}
                                                        variant={reactionPhase === 'Mixed' ? 'default' : 'ghost'}
                                                        size="sm"
                                                        className="flex-1 text-xs px-3 py-1 h-auto"
                                                >
                                                        Mixed
                                                </Button>
                                        </div>

                                        {/* --- NEW: Fluid Package Dropdown --- */}
                                        {reactionPhase === 'Mixed' && (
                                            <div>
                                                <Label>Fluid Package</Label>
                                                <Select value={fluidPackage} onValueChange={v => setFluidPackage(v as VleFluidPackage)}>
                                                    <SelectTrigger className="w-full mt-1">
                                                        <SelectValue />
                                                    </SelectTrigger>
                                                    <SelectContent>
                                                        <SelectItem value="Ideal">Ideal</SelectItem>
                                                        <SelectItem value="NRTL">NRTL</SelectItem>
                                                        <SelectItem value="Wilson">Wilson</SelectItem>
                                                        <SelectItem value="UNIQUAC">UNIQUAC</SelectItem>
                                                        <SelectItem value="Peng-Robinson">Peng-Robinson</SelectItem>
                                                        <SelectItem value="SRK">SRK</SelectItem>
                                                    </SelectContent>
                                                </Select>
                                            </div>
                                        )}

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

                {/* Components Card */}
                <div
                    ref={containerRef}
                    className={`p-6 rounded-lg shadow-lg bg-card text-card-foreground`}
                >
                    <div className="mb-4">
                        <CardTitle className="flex items-center justify-between">
                            Components
                            <Button variant="outline" size="sm" onClick={addComponentSetup}>
                                <PlusCircle className="mr-2 h-4 w-4" />
                                Add Component
                            </Button>
                        </CardTitle>
                    </div>
                    {/* ▼▼▼ CHANGE IS HERE ▼▼▼ */}
                    {/* The <div> with "overflow-x-auto" has been removed to fix clipping */}
                    <>
                                                {/* --- UPDATED TWO-ROW HEADER --- */}
                                                <div
                                                    className="grid items-center text-xs font-medium text-muted-foreground"
                                                    style={{ gridTemplateColumns: componentsGridCols }}
                                                >
                                                    {/* Row 1: Main Labels */}
                                                    <div /> {/* Spacer for trash icon */}
                                                    <div className="text-center">Name</div>
                                                    <div className="text-center">Molar Weight</div>
                                                    {reactionPhase === 'Mixed' ? (
                                                        <div className="text-center">Phase</div>
                                                    ) : reactionPhase === 'Liquid' ? (
                                                        <div className="text-center">Density</div>
                                                    ) : null}

                                                    {/* Reaction Titles */}
                                                    {reactionsSetup.map((r, index) => {
                                                        const numSideReactions = reactionsSetup.length - 1
                                                        let reactionTitle = 'Primary Reaction'
                                                        if (index > 0) {
                                                            reactionTitle =
                                                                numSideReactions === 1
                                                                    ? 'Side Reaction'
                                                                    : `Side Reaction ${index}`
                                                        }
                                                        return (
                                                            <div key={r.id} className="text-center col-span-2">
                                                                {reactionTitle}
                                                            </div>
                                                        )
                                                    })}
                                                </div>

                                                <div
                                                    className="grid items-center text-xs font-medium text-muted-foreground pb-2 border-b"
                                                    style={{ gridTemplateColumns: componentsGridCols }}
                                                >
                                                    {/* Row 2: Units and Sub-headers */}
                                                    <div /> {/* Spacer for trash icon */}
                                                    <div /> {/* Spacer for Name */}
                                                    <div className="text-center">(g/mol)</div>
                                                    {reactionPhase === 'Mixed' ? (
                                                        <div /> 
                                                    ) : reactionPhase === 'Liquid' ? (
                                                        <div className="text-center">(g/L)</div>
                                                    ) : null}

                                                    {/* Reaction Sub-headers */}
                                                    {reactionsSetup.map(r => (
                                                        <React.Fragment key={r.id}>
                                                            <div className="text-center">Stoich.</div>
                                                            <div className="text-center">
                                                                {r.isEquilibrium ? 'Order (fwd/rev)' : 'Order'}
                                                            </div>
                                                        </React.Fragment>
                                                    ))}
                                                </div>

                        {/* The actual component input rows */}
                        {componentsSetup.map((comp, index) => (
                            <div
                                key={comp.id}
                                className="grid gap-2 items-center py-2"
                                style={{ gridTemplateColumns: componentsGridCols }}
                            >
                                {index > 1 ? (
                                    <Button
                                        variant="ghost"
                                        size="icon"
                                        className="h-8 w-8"
                                        onClick={() => removeComponentSetup(comp.id)}
                                    >
                                        <Trash2 className="h-4 w-4 text-red-500" />
                                    </Button>
                                ) : (
                                    <Button variant="ghost" size="icon" className="h-8 w-8" disabled>
                                        <Trash2 className="h-4 w-4 text-muted-foreground" />
                                    </Button>
                                )}

                                <div className="relative">
                                    <div className="flex items-center">
                                        <Input
                                            value={comp.name}
                                            onChange={e => handleNameChange(e, comp.id)}
                                            onFocus={() => {
                                                setActiveInputId(comp.id)
                                                if (comp.name && comp.name.length > 1) {
                                                    handleFetchSuggestions(comp.name)
                                                }
                                            }}
                                            className="h-8 text-center"
                                            autoComplete="off"
                                        />
                                        {comp.isDbLocked && (
                                            <button
                                                onClick={() => handleUnlockComponent(comp.id)}
                                                className="absolute right-2 text-muted-foreground hover:text-foreground transition-colors"
                                                title="Unlock to edit manually"
                                            >
                                                <Unlink2 className="h-4 w-4" />
                                            </button>
                                        )}
                                    </div>
                                    {/* Suggestions Dropdown */}
                                    {showSuggestions &&
                                        activeInputId === comp.id &&
                                        suggestions.length > 0 && (
                                            <div className="absolute z-50 w-full bg-background border border-input rounded-md shadow-lg mt-1">
                                                {suggestions.map((s, i) => (
                                                    <div
                                                        key={i}
                                                        onClick={() => handleSuggestionClick(s)}
                                                        className="px-3 py-2 hover:bg-accent cursor-pointer text-sm text-center"
                                                    >
                                                        {s}
                                                    </div>
                                                ))}
                                            </div>
                                        )}
                                </div>
                                <Input
                                    type="number"
                                    placeholder="g/mol"
                                    value={comp.molarMass || ''}
                                    onChange={e =>
                                        handleComponentSetupChange(comp.id, 'molarMass', e.target.value)
                                    }
                                    className={`h-8 text-center transition-colors ${
                                        comp.isDbLocked ? 'bg-muted opacity-75' : ''
                                    }`}
                                    readOnly={comp.isDbLocked}
                                />
                                {/* --- NEW: Conditional Phase Selector or Density Input --- */}
                                {/* --- Finalized Phase/Density Column --- */}
                                {reactionPhase === 'Mixed' ? (
                                    <Select
                                        value={comp.phase}
                                        onValueChange={value =>
                                            handleComponentPropertyChange(comp.id, 'phase', value as any)
                                        }
                                    >
                                        {/* --- CHANGE #3: Consistent Width --- */}
                                        <SelectTrigger className="h-8 w-[90px] text-xs mx-auto">
                                            {/* --- CHANGE #2: New Placeholder --- */}
                                            <SelectValue placeholder="L / G" />
                                        </SelectTrigger>
                                        <SelectContent>
                                            {/* --- CHANGE #2: New Item Labels --- */}
                                            <SelectItem value="Liquid">L</SelectItem>
                                            <SelectItem value="Gas">G</SelectItem>
                                        </SelectContent>
                                    </Select>
                                ) : (
                                    reactionPhase === 'Liquid' && (
                                        <Input
                                            type={comp.isDbLocked ? 'text' : 'number'}
                                            placeholder="g/L"
                                            value={comp.isDbLocked ? 'T. dependent' : comp.density || ''}
                                            onChange={e =>
                                                handleComponentPropertyChange(comp.id, 'density', e.target.value)
                                            }
                                            className={`h-8 text-center transition-colors ${
                                                comp.isDbLocked ? 'bg-muted opacity-75 italic' : ''
                                            }`}
                                            readOnly={comp.isDbLocked}
                                        />
                                    )
                                )}
                                {reactionsSetup.map(r => (
                                    <React.Fragment key={r.id}>
                                        <Input
                                            type="number"
                                            placeholder="0"
                                            value={comp.reactionData[r.id]?.stoichiometry || ''}
                                            onChange={e =>
                                                handleComponentSetupChange(comp.id, 'reactionData', {
                                                    ...comp.reactionData,
                                                    [r.id]: {
                                                        ...comp.reactionData[r.id],
                                                        stoichiometry: e.target.value,
                                                    },
                                                })
                                            }
                                            className="h-8 text-center"
                                            step="0.1"
                                        />
                                        <div className="flex items-center gap-1">
                                            <Input
                                                type="number"
                                                placeholder={r.isEquilibrium ? 'fwd' : '0'}
                                                value={comp.reactionData[r.id]?.order || ''}
                                                onChange={e =>
                                                    handleComponentSetupChange(comp.id, 'reactionData', {
                                                        ...comp.reactionData,
                                                        [r.id]: {
                                                            ...comp.reactionData[r.id],
                                                            order: e.target.value,
                                                        },
                                                    })
                                                }
                                                className="h-8 text-center"
                                                step="0.1"
                                            />
                                            {r.isEquilibrium && (
                                                <Input
                                                    type="number"
                                                    placeholder="rev"
                                                    value={comp.reactionData[r.id]?.orderReverse || ''}
                                                    onChange={e =>
                                                        handleComponentSetupChange(comp.id, 'reactionData', {
                                                            ...comp.reactionData,
                                                            [r.id]: {
                                                                ...comp.reactionData[r.id],
                                                                orderReverse: e.target.value,
                                                            },
                                                        })
                                                    }
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
                    {/* ▲▲▲ END OF CHANGES ▲▲▲ */}
                                </div>

                    {/* --- Corrected Henry's Law Card --- */}
                    {reactionPhase === 'Mixed' && hasGasComponent && (
                        <Card className="p-6">
                            <CardHeader className="p-0 pb-4 flex flex-row justify-between items-center">
                                <CardTitle className="flex items-center gap-2">
                                    Henry's Law Constants
                                    <Tooltip>
                                        <TooltipTrigger asChild>
                                            <button>
                                                <Info className="h-4 w-4 text-muted-foreground" />
                                            </button>
                                        </TooltipTrigger>
                                        <TooltipContent className="max-w-xs">
                                            <p className="text-sm font-medium">
                                                H = P*y/x
                                            </p>
                                            <p className="text-xs text-muted-foreground mt-1">
                                                The Henry volatility defined via aqueous-phase mixing ratio.
                                            </p>
                                        </TooltipContent>
                                    </Tooltip>
                                </CardTitle>
                                <div className="flex items-center gap-2">
                                    <Select value={henryUnit} onValueChange={(v) => setHenryUnit(v as any)}>
                                        <SelectTrigger className="w-[90px] h-9 text-xs">
                                            <SelectValue />
                                        </SelectTrigger>
                                        <SelectContent>
                                            <SelectItem value="Pa">Pa</SelectItem>
                                            <SelectItem value="kPa">kPa</SelectItem>
                                            <SelectItem value="bar">bar</SelectItem>
                                        </SelectContent>
                                    </Select>
                                    <Button
                                        variant="outline"
                                        size="sm"
                                        onClick={handleAddTemperatureRow}
                                    >
                                        <PlusCircle className="mr-2 h-4 w-4" />
                                        Add Temp
                                    </Button>
                                </div>
                            </CardHeader>
                            <CardContent className="space-y-2 p-0">
                                {/* Table Header */}
                                <div
                                    className="grid items-center text-xs font-medium text-muted-foreground pb-2 border-b"
                                    style={{ gridTemplateColumns: henrysLawGridCols }}
                                >
                                    <div /> {/* Spacer for delete button */}
                                    <div className="text-center">Temperature (°C)</div>
                                    {gasComponents.map(comp => (
                                        <div key={comp.id} className="text-center font-semibold">
                                            {comp.name}
                                        </div>
                                    ))}
                                </div>

                                {/* Table Rows */}
                                {henryTemperatures.map((row, index) => (
                                    <div
                                        key={row.id}
                                        className="grid items-center gap-2"
                                        style={{ gridTemplateColumns: henrysLawGridCols }}
                                    >
                                        <Button
                                            variant="ghost"
                                            size="icon"
                                            className="h-8 w-8"
                                            onClick={() => handleRemoveTemperatureRow(row.id)}
                                            // --- FIX: Trash can is disabled based on row index ---
                                            disabled={index < 2}
                                        >
                                            <Trash2
                                                className={`h-4 w-4 ${
                                                    // --- FIX: Icon is greyed out based on row index ---
                                                    index < 2
                                                        ? 'text-muted-foreground'
                                                        : 'text-red-500'
                                                }`}
                                            />
                                        </Button>
                                        <Input
                                            type="number"
                                            placeholder="e.g., 80"
                                            value={row.temp}
                                            onChange={e => handleTemperatureChange(row.id, e.target.value)}
                                            className="h-8 text-center"
                                            // --- FIX: The "disabled" prop has been removed ---
                                        />
                                        {gasComponents.map(comp => (
                                            <Input
                                                key={comp.id}
                                                type="number"
                                                placeholder={`H (${henryUnit}·m³/mol)`}
                                                value={comp.henryData?.[row.id] || ''}
                                                onChange={e =>
                                                    handleHenryConstantChange(
                                                        comp.id,
                                                        row.id,
                                                        e.target.value
                                                    )
                                                }
                                                className="h-8 text-center"
                                                // --- FIX: The "disabled" prop has been removed ---
                                            />
                                        ))}
                                    </div>
                                ))}
                            </CardContent>
                        </Card>
                    )}
                    {/* --- Interaction Parameters Card (shows when a non-ideal package is selected) --- */}
                    {/* Interaction Parameters card removed - handled at parent level */}
        </div>
      </div>
    </TooltipProvider>
  );
};

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
    // --- ADD THESE NEW PROPS ---
    reactorTypes, setReactorTypes,
    henryTemperatures,
    henryUnit,
    fluidPackage,
    vleComponentData,
    interactionParameters,
}: { 
    onBack: () => void, 
    reactions: ReactionSetup[], 
    components: ComponentSetup[],
    simBasis: any,
    molarRatios: MolarRatio[],
    setMolarRatios: React.Dispatch<React.SetStateAction<MolarRatio[]>>,
    reactionPhase: 'Liquid' | 'Gas' | 'Mixed',
    prodRate: string,
    temperature: number,
    setTemperature: React.Dispatch<React.SetStateAction<number>>,
    pressure: number,
    setPressure: React.Dispatch<React.SetStateAction<number>>,
    // --- ADD THESE NEW PROPS ---
    reactorTypes: { pfr: boolean; cstr: boolean },
    setReactorTypes: React.Dispatch<React.SetStateAction<{ pfr: boolean; cstr: boolean }>>,
    henryTemperatures: { id: string; temp: string }[],
    henryUnit: 'Pa' | 'kPa' | 'bar',
    fluidPackage: VleFluidPackage,
    vleComponentData: Map<string, CompoundData>,
    interactionParameters: Map<string, InteractionParams | null>,
}) => {
  const { resolvedTheme } = useTheme();
  const isDark = resolvedTheme === 'dark';
  const textColor = isDark ? '#E5E7EB' : '#1F2937';
  const cardBg = 'bg-card';
  const cardFg = 'text-card-foreground';
  const mainBg = 'bg-background';
  const mainFg = 'text-foreground';

  // reactorTypes/state are passed in from parent; use prop-based handler
  const handleReactorTypeChange = (type: 'pfr' | 'cstr') => {
      setReactorTypes(prev => {
          const newState = { ...prev, [type]: !prev[type] }
          if (!newState.pfr && !newState.cstr) {
              return prev // prevent deselecting both
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
    console.log("=============== NEW CALCULATION ===============");
    const tempK = convertUnit.temperature.celsiusToKelvin(temperature);
    
    const limitingReactant = components.find(c => c.id === simBasis.limitingReactantId);
    const desiredProduct = components.find(c => c.id === simBasis.desiredProductId);

    if (!limitingReactant || !desiredProduct) {
        return { series: [], xAxis: '', yAxis: '', legend: [], yMax: 'dataMax' }
    }

    const xLabel = `Conversion of ${limitingReactant?.name || 'Limiting Reactant'}`;
    let yLabel = '';
    
    // We will now build a list of ALL series, and toggle their visibility
    const allSeries: any[] = [];
    const allLegends: string[] = [];
    
    const hardCap = 5000; // The maximum practical volume
    const colors = {
        pfr: '#22c55e', // green
        cstr: '#E53E3E'  // red
    };

    // This will be calculated using ONLY data from visible series
    let overallReasonableMax = -1;
    const visibleVolumeDataSets: { x: number; y: number }[][] = [];


    for (const type in reactorTypes) {
        if (reactorTypes[type as keyof typeof reactorTypes]) {
            const reactorType = type.toUpperCase() as 'PFR' | 'CSTR';
            let dataToShow: {
                selectivityData: { x: number; y: number }[];
                volumeData: { x: number; y: number }[];
            } = { selectivityData: [], volumeData: [] };

            // --- Calculation for each type ---
            if (reactorType === 'PFR') {
                const pfrPhase: 'Liquid' | 'Gas' = (reactionPhase === 'Mixed') ? 'Gas' : reactionPhase;
                dataToShow = calculatePfrData_Optimized(
                    reactions, components, tempK, pressure, molarRatios, simBasis, parseFloat(prodRate), pfrPhase,
                    // --- ADD THESE ---
                    fluidPackage,
                    vleComponentData,
                    interactionParameters
                );
            } else { // CSTR
                const F_A0_guess = 10.0;
                const initialFlowRates = { [limitingReactant.name]: F_A0_guess };
                molarRatios.forEach(ratio => {
                    const compName = components.find(c => c.id === ratio.numeratorId)?.name;
                    if (compName) initialFlowRates[compName] = F_A0_guess * ratio.value;
                });
                components.forEach(c => {
                    if (initialFlowRates[c.name] === undefined) initialFlowRates[c.name] = 0;
                });
                
                dataToShow = calculateCstrData(
                    reactions, components, initialFlowRates, tempK, simBasis, 
                    reactionPhase, pressure, parseFloat(prodRate), henryUnit,
                    fluidPackage, vleComponentData, interactionParameters,
                    henryTemperatures // 👈 ADD THIS
                );
            }
            
            const activeReactorCount = (reactorTypes.pfr ? 1 : 0) + (reactorTypes.cstr ? 1 : 0);
            const seriesNamePrefix = activeReactorCount > 1 ? `${reactorType.toUpperCase()} - ` : '';
            
            let seriesObject: any = {
                type: 'line',
                smooth: true,
                showSymbol: false,
                color: colors[type as keyof typeof colors],
                lineStyle: { width: 4 }
            };

            if (graphType === 'volume') {
                yLabel = 'Reactor Volume (m³)';
                const volumePoints = dataToShow.volumeData;
                seriesObject.name = `${seriesNamePrefix}Reactor Volume`;
                seriesObject.data = volumePoints.map(d => [d.x, d.y]);
                
                // Check for visibility but don't filter the series
                const isVisible = volumePoints.some(p => p.y <= hardCap);
                
                // --- THIS IS THE KEY CHANGE ---
                // The series is ALWAYS included, but its visibility is toggled.
                seriesObject.show = isVisible;
                
                if (isVisible) {
                    visibleVolumeDataSets.push(volumePoints);
                }

            } else { // Selectivity graph is always visible
                yLabel = `Selectivity to ${desiredProduct?.name || 'Product'}`;
                seriesObject.name = `${seriesNamePrefix}Selectivity`;
                seriesObject.data = dataToShow.selectivityData.map(d => [d.x, d.y]);
                seriesObject.show = true;
            }
            
            allSeries.push(seriesObject);
            allLegends.push(seriesObject.name);
        }
    }
    
    // Calculate the Y-axis max based ONLY on the visible datasets
    if (graphType === 'volume') {
        for (const volumePoints of visibleVolumeDataSets) {
             if (volumePoints.length > 1) {
                const volAt5Percent = interpolateFromData(0.05, volumePoints);
                const volAt95Percent = interpolateFromData(0.95, volumePoints);
                let curveReasonableMax = -1;
                if (isFinite(volAt5Percent) && isFinite(volAt95Percent)) {
                    curveReasonableMax = Math.max(volAt5Percent, volAt95Percent);
                } else if (isFinite(volAt5Percent)) {
                    curveReasonableMax = volAt5Percent;
                } else if (isFinite(volAt95Percent)) {
                    curveReasonableMax = volAt95Percent;
                }

                if (curveReasonableMax > overallReasonableMax) {
                    overallReasonableMax = curveReasonableMax;
                }
            }
        }
    }

    let yAxisMax: number | 'dataMax' = 'dataMax';
    if (graphType === 'volume') {
        if (overallReasonableMax > 0) {
            yAxisMax = Math.min(overallReasonableMax, hardCap);
        } else {
            yAxisMax = 100; // Default max if no curves are visible
        }
    }
    
    // Return ALL series and legends. ECharts handles showing/hiding.
    return { series: allSeries, xAxis: xLabel, yAxis: yLabel, legend: allLegends, yMax: yAxisMax }

}, [reactorTypes, reactionPhase, pressure, molarRatios, temperature, graphType, reactions, components, simBasis, prodRate, henryTemperatures, henryUnit, fluidPackage, vleComponentData, interactionParameters])

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
        
        max: graphType === 'selectivity' ? 1 : graphData.yMax,

        axisTick: {
            show: true
        },

        axisLine: { lineStyle: { color: textColor } },
        axisLabel: {
            color: textColor,
            fontSize: 14,
            fontFamily: 'Merriweather Sans',
            formatter: (value: number) => {
                if (
                    graphType === 'volume' &&
                    typeof graphData.yMax === 'number' &&
                    value >= graphData.yMax
                ) {
                    return '';
                }
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
            // --- CHANGE #1: Use toFixed(2) for the axis pointer label ---
            return parseFloat(params.value as string).toFixed(2);
          }
        }
      },

      formatter: (params: any) => { 
         if (!params || params.length === 0) { 
             return '' 
         } 

         // --- CHANGE #2: Use toFixed(2) for the tooltip header ---
         const xAxisValue = parseFloat(params[0].axisValue).toFixed(2);
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
    const [reactionPhase, setReactionPhase] = useState<'Liquid' | 'Gas' | 'Mixed'>('Liquid'); 
    // --- ADD THIS NEW STATE ---
    const [fluidPackage, setFluidPackage] = useState<VleFluidPackage>('Ideal');
    const [interactionParameters, setInteractionParameters] = useState<Map<string, InteractionParams | null>>(new Map());
        // This map will store fetched pure component VLE data
        const [vleComponentData, setVleComponentData] = useState<Map<string, CompoundData>>(new Map());
  const [temperature, setTemperature] = useState(4); // State lifted here
  const [pressure, setPressure] = useState(25);       // State lifted here 
  const [isFetching, setIsFetching] = useState(false) // Add loading state
    // Loading and error states for VLE preparation
    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);

    // --- ADD THIS NEW STATE FOR HENRY'S LAW TABLE ---
    const [henryTemperatures, setHenryTemperatures] = useState<
        { id: string; temp: string }[]
    >([
        { id: `temp-${Date.now()}-1`, temp: '' },
        { id: `temp-${Date.now()}-2`, temp: '' },
    ]);
    const [henryUnit, setHenryUnit] = useState<'Pa' | 'kPa' | 'bar'>('Pa');
    // --- END ADD ---

        // --- NEW: reactorTypes state to control which reactors are active ---
        const [reactorTypes, setReactorTypes] = useState<{ pfr: boolean; cstr: boolean }>({ pfr: true, cstr: false });

  // --- ADD THIS FUNCTION TO FETCH COMPONENT DATA ---
  const fetchComponentProperties = useCallback(
    async (name: string, componentId: string) => {
      if (!supabase) {
        console.error('Supabase not initialized.')
        return
      }
      setIsFetching(true)
      try {
                const { data: compoundDbData, error: compoundError } = await supabase
                    .from('compound_properties')
                    .select('id, name, properties')
                    .ilike('name', name)
                    .limit(1)
                    .single()

        if (compoundError || !compoundDbData) {
          throw new Error(`Compound '${name}' not found.`)
        }

        if (!compoundDbData.properties) {
          throw new Error(`No properties found for ${compoundDbData.name}.`)
        }

        // Add debugging to see what we're getting
        console.log('Fetched compound data:', compoundDbData.name, compoundDbData.properties);
        console.log('Available property keys:', Object.keys(compoundDbData.properties || {}));

    const liquidDensityData = compoundDbData.properties['LiquidDensity'] || compoundDbData.properties['Liquid density']
    const cas_number = compoundDbData.properties['CasNumber'] || compoundDbData.properties['cas_number'];
        
        // Handle molecular weight - it might be an object or a simple value
        let molarMass = compoundDbData.properties['MolecularWeight'] || compoundDbData.properties['molecular_weight'];
        
        // Try other possible keys if the standard ones don't work
        if (!molarMass) {
          const possibleKeys = ['MW', 'mw', 'MolarMass', 'molar_mass', 'mol_weight'];
          for (const key of possibleKeys) {
            if (compoundDbData.properties[key]) {
              molarMass = compoundDbData.properties[key];
              break;
            }
          }
        }
        
        if (typeof molarMass === 'object' && molarMass !== null) {
          // If it's an object, try to extract the value
          molarMass = molarMass.value || molarMass.Value || molarMass.val || Object.values(molarMass)[0];
        }
        
        console.log('Extracted molar mass:', molarMass);
        
        const critTemp = parseCoefficient(
          compoundDbData.properties['CriticalTemperature'] || compoundDbData.properties['Critical temperature']
        )
        let densityValue = ''

        if (liquidDensityData && liquidDensityData.eqno && critTemp) {
          const tempK = convertUnit.temperature.celsiusToKelvin(temperature)
          const densityCoeffs = {
            A: parseCoefficient(liquidDensityData.A),
            B: parseCoefficient(liquidDensityData.B),
            C: parseCoefficient(liquidDensityData.C),
            D: parseCoefficient(liquidDensityData.D),
          }
          // Assuming DB density is kg/m^3 which is equivalent to g/L.
          const densityKgM3 = calculatePropertyByEquation(
            liquidDensityData.eqno,
            tempK,
            densityCoeffs,
            critTemp
          )
          if (densityKgM3 !== null) {
            densityValue = densityKgM3.toPrecision(4)
          }
        }

        setComponentsSetup(prev =>
          prev.map(c =>
                                c.id === componentId
                            ? {
                                    ...c,
                                    name: compoundDbData.name,
                                    molarMass: molarMass ? String(molarMass) : '',
                                    density: densityValue,
                                    liquidDensityData: liquidDensityData,
                                    criticalTemp: critTemp,
                                    cas_number: cas_number,
                                    isDbLocked: true,
                                }
              : c
          )
        )
      } catch (err: any) {
        console.error('Error fetching component properties:', err.message)
      } finally {
        setIsFetching(false)
      }
    },
    [temperature]
  ) // Depends on temperature for initial calculation

  // --- Data Fetching helper used across VLE tools (adapted from mccabe-thiele) ---
  async function fetchCompoundDataLocal(compoundName: string): Promise<CompoundData | null> {
    if (!supabase) { throw new Error('Supabase client not initialized.'); }
    try {
        const { data: compoundDbData, error: compoundError } = await supabase
            .from('compound_properties')
            .select('id, name, properties')
            .ilike('name', compoundName)
            .limit(1);

        if (compoundError) throw new Error(`Supabase compound query error: ${compoundError.message}`);
        if (!compoundDbData || compoundDbData.length === 0) throw new Error(`Compound '${compoundName}' not found.`);
        const compound = compoundDbData[0];
        const foundName = compound.name;
        const casNumber = compound.properties?.CasNumber || compound.properties?.cas_number;

        // Use the properties directly from the compound_properties table
        const properties = compound.properties;
        if (!properties) {
            throw new Error(`No properties found for ${foundName}.`);
        }

        // Extract Antoine params
        let antoine: any | null = null;
        const antoineChemsep = properties.Antoine || properties.AntoineVaporPressure;
        if (antoineChemsep?.A && antoineChemsep.B && antoineChemsep.C) {
            antoine = {
                A: parseFloat(antoineChemsep.A?.value ?? antoineChemsep.A),
                B: parseFloat(antoineChemsep.B?.value ?? antoineChemsep.B),
                C: parseFloat(antoineChemsep.C?.value ?? antoineChemsep.C),
                Tmin_K: parseFloat(antoineChemsep.Tmin?.value ?? antoineChemsep.Tmin ?? 0),
                Tmax_K: parseFloat(antoineChemsep.Tmax?.value ?? antoineChemsep.Tmax ?? 10000),
                Units: antoineChemsep.units || 'Pa',
                EquationNo: antoineChemsep.eqno
            };
        }
        if (!antoine || isNaN(antoine.A) || isNaN(antoine.B) || isNaN(antoine.C)) throw new Error(`Failed to extract valid Antoine params for ${foundName}.`);

        // UNIFAC groups (optional)
        let unifacGroups: any = null;
        if (properties.elements_composition?.UNIFAC) {
            unifacGroups = {};
            for (const key in properties.elements_composition.UNIFAC) {
                const subgroupId = parseInt(key); const count = parseInt(properties.elements_composition.UNIFAC[key]);
                if (!isNaN(subgroupId) && !isNaN(count) && count > 0) unifacGroups[subgroupId] = count;
            }
            if (Object.keys(unifacGroups).length === 0) unifacGroups = null;
        }

        // PR/SRK/Uniquac/Wilson params
        const tcPropObj = properties['Critical temperature'];
        const pcPropObj = properties['Critical pressure'];
        const omegaPropObj = properties['Acentric factor'];
        let prParams = null;
        let srkParams = null;
        if (tcPropObj && pcPropObj && omegaPropObj) {
            const Tc_K_val = parseFloat(tcPropObj.value);
            const pcValue = parseFloat(pcPropObj.value);
            const pcUnits = String(pcPropObj.units).toLowerCase();
            let Pc_Pa_val: number;
            if (pcUnits === 'pa') Pc_Pa_val = pcValue;
            else if (pcUnits === 'kpa') Pc_Pa_val = pcValue * 1e3;
            else if (pcUnits === 'mpa') Pc_Pa_val = pcValue * 1e6;
            else if (pcUnits === 'bar') Pc_Pa_val = pcValue * 1e5;
            else { Pc_Pa_val = pcValue; console.warn(`Unknown Pc units ('${pcUnits}') for ${foundName}. Assuming Pa.`); }
            prParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: parseFloat(omegaPropObj.value) };
            srkParams = { Tc_K: Tc_K_val, Pc_Pa: Pc_Pa_val, omega: parseFloat(omegaPropObj.value) };
        }

        // Wilson params (liquid molar volume)
        let wilsonParams = null;
        const vLObj = properties['Liquid molar volume'];
        if (vLObj && vLObj.value) {
            const vVal = parseFloat(vLObj.value);
            // Assume units are m^3/mol if small, otherwise convert if units provided (best-effort)
            wilsonParams = { V_L_m3mol: vVal };
        }

        const compoundData: CompoundData = {
            name: foundName,
            cas_number: casNumber,
            antoine,
            unifacGroups: unifacGroups,
            prParams: prParams,
            srkParams: srkParams,
            uniquacParams: null,
            wilsonParams: wilsonParams,
        };
        return compoundData;
    } catch (err: any) {
        console.error('fetchCompoundDataLocal error:', err.message || err);
        return null;
    }
  }

  // --- NEW: Prepare VLE data function (fetches pure component data and binary params) ---
    const prepareVleData = async () => {
        // Skip fetching VLE params when using Ideal model or when the reaction phase is not Mixed
        if (fluidPackage === 'Ideal' || reactionPhase !== 'Mixed' || !supabase) {
            setInteractionParameters(new Map()); // Clear params just in case
            return true; // Immediately exit with success
        }
    
    console.log(`[DEBUG] 🚀 Starting VLE data preparation for: ${fluidPackage}`);
    console.log(`🔵 [Debug] Starting to build ${fluidPackage} Parameter Map for ${componentsSetup.length} components...`);
    setLoading(true);
    setError(null);

    try {
      // 1. Fetch Pure Component VLE Data (this part is unchanged)
      const fetchedVleData = new Map(vleComponentData);
      const dataPromises = componentsSetup.map(async (comp) => {
        if (!comp.name || fetchedVleData.has(comp.id)) return;
        const data = await fetchCompoundDataLocal(comp.name);
        if (data) {
          fetchedVleData.set(comp.id, data);
        } else {
          throw new Error(`Failed to fetch VLE properties for ${comp.name}.`);
        }
      });
      await Promise.all(dataPromises);
      setVleComponentData(fetchedVleData);

      // 2. Fetch Binary Interaction Parameters and store by INDEX
      const newInteractionParams = new Map();
      const componentIdToIndexMap = new Map<string, number>();
      componentsSetup.forEach((comp, index) => {
          componentIdToIndexMap.set(comp.id, index);
      });

      const paramPromises = [];
      for (let i = 0; i < componentsSetup.length; i++) {
          for (let j = i + 1; j < componentsSetup.length; j++) {
              const comp1_setup = componentsSetup[i];
              const comp2_setup = componentsSetup[j];
              
              const comp1_vle = fetchedVleData.get(comp1_setup.id);
              const comp2_vle = fetchedVleData.get(comp2_setup.id);

              if (!comp1_vle?.cas_number || !comp2_vle?.cas_number) continue;

              const key = `${i}-${j}`;

              const promise = (async () => {
                  let params: InteractionParams | null = null;
                  let sourceType = 'Unknown';
                  try {
                      const isLiquidModel = fluidPackage === 'NRTL' || fluidPackage === 'Wilson' || fluidPackage === 'UNIQUAC';
                      const involvesGas = comp1_setup.phase === 'Gas' || comp2_setup.phase === 'Gas';

                      if (isLiquidModel && involvesGas) return; 

                      switch (fluidPackage) {
                          case 'NRTL':
                              params = await fetchNrtlParameters(supabase, comp1_vle.cas_number!, comp2_vle.cas_number!, (source) => { sourceType = source; });
                              break;
                          case 'Wilson':
                              params = await fetchWilsonInteractionParams(supabase, comp1_vle.cas_number!, comp2_vle.cas_number!);
                              sourceType = 'Database'; // Assuming no estimation for now
                              break;
                          case 'UNIQUAC':
                              params = await fetchUniquacInteractionParams(supabase, comp1_vle.cas_number!, comp2_vle.cas_number!);
                              sourceType = 'Database'; // Assuming no estimation for now
                              break;
                          case 'Peng-Robinson':
                              params = await fetchPrInteractionParams(supabase, comp1_vle.cas_number!, comp2_vle.cas_number!);
                              sourceType = 'Database'; // Assuming no estimation for now
                              break;
                          case 'SRK':
                              params = await fetchSrkInteractionParams(supabase, comp1_vle.cas_number!, comp2_vle.cas_number!);
                              sourceType = 'Database'; // Assuming no estimation for now
                              break;
                      }
                  } catch (err) {
                      console.error(`Error fetching params for pair ${i}-${j}:`, err);
                      params = null;
                  }

                  if (params) {
                      console.log(`- ✅ VLE params for pair ${i}-${j} (${comp1_setup.name} - ${comp2_setup.name}) [Source: ${sourceType}]:`, params);
                      newInteractionParams.set(key, params);
                  } else if (sourceType !== 'Unknown'){
                      console.error(`❌ FETCH FAILED for pair: ${comp1_setup.name} - ${comp2_setup.name}`);
                  }
              })();
              paramPromises.push(promise);
          }
      }
            
      await Promise.all(paramPromises);
      console.log('[DEBUG] 🏁 Final interaction parameter map:', newInteractionParams);
      setInteractionParameters(newInteractionParams as any);

      return true;

    } catch (err: any) {
        setError(err.message);
        return false;
    } finally {
        setLoading(false);
    }
  };

  // --- ADD THIS EFFECT TO UPDATE DENSITY ON TEMP CHANGE ---
  useEffect(() => {
    const tempK = convertUnit.temperature.celsiusToKelvin(temperature)

    const needsUpdate = componentsSetup.some(
      c => c.isDbLocked && c.liquidDensityData
    )

    if (needsUpdate) {
      setComponentsSetup(prev =>
        prev.map(c => {
          if (c.isDbLocked && c.liquidDensityData && c.criticalTemp) {
            const { eqno, ...coeffs } = c.liquidDensityData
            const densityKgM3 = calculatePropertyByEquation(
              eqno,
              tempK,
              coeffs,
              c.criticalTemp
            )
            if (densityKgM3 !== null) {
              return { ...c, density: densityKgM3.toPrecision(4) }
            }
          }
          return c
        })
      )
    }
  }, [temperature]) 

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
            if (reactant.name === 'Carbon monoxide') defaultValue = 24.0; // DMC Preset
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
                            onBack={() => {
                                // A trick to make Henry's Law data available to the CSTR solver
                                (window as any).__HENRY_TEMPS = null; 
                                setShowSimulator(false)
                            }} 
                            reactions={reactionsSetup} 
                            components={componentsSetup} 
                            simBasis={simBasis}
                            molarRatios={molarRatios}
                            setMolarRatios={setMolarRatios}
                            reactionPhase={reactionPhase}
                            prodRate={prodRate}
                            temperature={temperature}
                            setTemperature={setTemperature}
                            pressure={pressure}
                            setPressure={setPressure}
                            // --- PASS NEW PROPS ---
                            reactorTypes={reactorTypes}
                            setReactorTypes={setReactorTypes}
                            henryTemperatures={henryTemperatures}
                            henryUnit={henryUnit}
                            fluidPackage={fluidPackage}
                            vleComponentData={vleComponentData}
                            interactionParameters={interactionParameters}
                        />
  } else {
        return <KineticsInput 
                            onNext={async () => { // Make this function async
                                const success = await prepareVleData(); // Await the data fetching
                                if (success) { // Only proceed if data was fetched successfully
                                    (window as any).__HENRY_TEMPS = henryTemperatures;
                                    setShowSimulator(true)
                                } else {
                                    // Optionally, you can alert the user that something went wrong
                                    console.error("Failed to prepare VLE data. Cannot proceed to simulator.");
                                }
                            }} 
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
                            setTemperature={setTemperature}
                            setPressure={setPressure}
                            fetchComponentProperties={fetchComponentProperties}
                            isFetching={isFetching}
                            fluidPackage={fluidPackage}
                            setFluidPackage={setFluidPackage}
                            // --- PASS NEW PROPS ---
                            henryTemperatures={henryTemperatures}
                            setHenryTemperatures={setHenryTemperatures}
                            henryUnit={henryUnit}
                            setHenryUnit={setHenryUnit}
                            setReactorTypes={setReactorTypes}
                        />
  }
}