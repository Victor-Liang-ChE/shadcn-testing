import {
  Component,
  ParsedParallelReactions,
  // CalculationResults as SinglePassCalculationResults, // Not directly used, but good for context
  solveCSTRParallel,
  solvePFR_ODE_System,
  findLimitingReactant,
  PFR_ODE_ProfilePoint,
} from './reactor-solver';
import { Kinetics } from './reaction-parser';

const R_GAS_CONSTANT = 8.314; // J/(mol*K) or L*kPa/(mol*K)

// --- Type Definitions ---
export type Stream = { [componentName: string]: number };

// SeparatorRule interface is removed as per new logic

export interface RecycleCalculationResult {
  productStream: Stream;
  recycleStreamFinal: Stream;
  reactorInletStreamFinal: Stream;
  reactorOutletStreamFinal: Stream;
  overallConversion?: { reactantName: string; value: number };
  singlePassConversion?: { reactantName: string; value: number }; // Conversion in the reactor itself
  iterations: number;
  converged: boolean;
  error?: string; // For convergence errors
  calculationError?: string; // Error from the calculation process itself (e.g., bad inputs)
  debugInfo?: {
    iterationHistory: Array<{
      iteration: number;
      recycleStreamError: number;
      recycleStream: Stream;
      reactorOutlet: Stream;
    }>;
    finalError: number;
    tolerance: number;
    maxIterations: number;
    convergenceIssues: string[];
  };
}

// --- Helper Functions ---

/**
 * Mixes two streams (fresh feed and recycle stream).
 * Assumes component names in both streams are consistent.
 */
function mixStreams(freshFeed: Stream, recycleStream: Stream, allComponentNames: string[]): Stream {
  const mixedStream: Stream = {};
  for (const name of allComponentNames) {
    mixedStream[name] = (freshFeed[name] || 0) + (recycleStream[name] || 0);
  }
  return mixedStream;
}

/**
 * Separates a single stream into a product stream and a recycle stream.
 * Logic:
 * - All components are recycled at the specified recycle fraction.
 * - The remainder goes to the product stream.
 */
function separateStream(
  reactorOutletStream: Stream,
  allComponentNames: string[],
  recycleFractionPercent: number
): { productStream: Stream; recycleStream: Stream } {
  const productStream: Stream = {};
  const recycleStream: Stream = {};
  
  const recycleFraction = recycleFractionPercent / 100; // Convert percentage to fraction

  for (const name of allComponentNames) {
    const totalOut = reactorOutletStream[name] || 0;

    recycleStream[name] = totalOut * recycleFraction;
    productStream[name] = totalOut * (1 - recycleFraction);
  }
  return { productStream, recycleStream };
}

/**
 * Calculates the error between two streams to check for convergence.
 * Uses sum of absolute differences of molar flow rates for each component.
 */
function calculateStreamError(streamA: Stream, streamB: Stream, allComponentNames: string[]): number {
  let error = 0;
  for (const name of allComponentNames) {
    error += Math.abs((streamA[name] || 0) - (streamB[name] || 0));
  }
  return error;
}

/**
 * Creates an empty stream (all flow rates zero) for all known components.
 */
function createEmptyStream(allComponentNames: string[]): Stream {
  const stream: Stream = {};
  for (const name of allComponentNames) {
    stream[name] = 0;
  }
  return stream;
}

// --- Main Recycle System Solver ---
export async function solveRecycleSystem(
  parsedParallelReactions: ParsedParallelReactions,
  reactorVolume: number,
  freshFeedRates: Stream,
  volumetricFlowRateFresh: number, // v0 for the fresh feed stream (used for liquid phase)
  components: Component[], // Full component list from UI
  reactorType: 'PFR' | 'CSTR',
  reactionPhase: 'Liquid' | 'Gas',
  kinetics: Kinetics,
  recycleFractionPercent: number, // Recycle fraction percentage (0-100%)
  totalPressureInput?: string // Total pressure in bar (for gas phase)
): Promise<RecycleCalculationResult> {
  // Initial setup
  const allComponentNames = components.map(c => c.name);
  // const separatorRules = separatorRulesString.map(r => ({ // Removed
  //   componentName: r.componentName,
  //   recycleFraction: parseFloat(r.recycleFraction) || 0,
  // }));

  let recycleStream_guess: Stream = createEmptyStream(allComponentNames);
  let reactorInletStream_final: Stream = mixStreams(freshFeedRates, recycleStream_guess, allComponentNames);
  let reactorOutletStream_final: Stream = { ...reactorInletStream_final };
  let productStream_final: Stream = { ...freshFeedRates };

  let converged = false;
  const maxIterations = 1000;
  const tolerance = 1e-6; // Tolerance for sum of absolute flow rate differences
  let iterationCount = 0; // To store the number of iterations
  
  // Debug information tracking
  const iterationHistory: Array<{
    iteration: number;
    recycleStreamError: number;
    recycleStream: Stream;
    reactorOutlet: Stream;
  }> = [];
  const convergenceIssues: string[] = [];
  let finalError = 0;

  for (let i = 0; i < maxIterations; i++) {
    iterationCount = i + 1; // Update iteration count
    
    // 1. MIXING POINT BALANCE: Calculate the reactor's true inlet stream
    const reactorInletStream = mixStreams(freshFeedRates, recycleStream_guess, allComponentNames);

    // 2. CALCULATE VOLUMETRIC FLOW RATE (v0) FOR THE REACTOR INLET
    let v0_reactorInlet: number;
    const totalMolarFlowReactorInlet = Object.values(reactorInletStream).reduce((sum, f) => sum + f, 0);

    if (reactionPhase === 'Liquid') {
      const totalMolarFlowFresh = Object.values(freshFeedRates).reduce((sum, f) => sum + f, 0);
      // const totalMolarFlowRecycle = Object.values(recycleStream_guess).reduce((sum, f) => sum + f, 0); // Not directly needed for this v0 calc method
      if (totalMolarFlowFresh > 1e-9) {
        v0_reactorInlet = volumetricFlowRateFresh * (totalMolarFlowReactorInlet / totalMolarFlowFresh);
      } else if (totalMolarFlowReactorInlet > 1e-9) {
        v0_reactorInlet = volumetricFlowRateFresh; // Fallback if no fresh feed but recycle exists.
                                                  // This assumes recycle stream has similar density properties to what v0_fresh implies.
                                                  // A more rigorous approach might need density inputs.
      } else {
        v0_reactorInlet = 0;
      }
       if (v0_reactorInlet < 1e-9 && totalMolarFlowReactorInlet > 1e-9) {
        // This case might indicate an issue (e.g., v0_fresh was 0 but there's flow).
        // Proceeding, but concentrations could be very high.
      }
    } else { // Gas Phase
      const T_K = parseFloat(kinetics.reactionTempK);
      const P_total_kPa = (parseFloat(totalPressureInput || "1") || 1) * 100;
      if (T_K <= 0 || P_total_kPa <= 0) {
        return {
          productStream: createEmptyStream(allComponentNames),
          recycleStreamFinal: recycleStream_guess,
          reactorInletStreamFinal: reactorInletStream,
          reactorOutletStreamFinal: createEmptyStream(allComponentNames),
          iterations: iterationCount,
          converged: false,
          calculationError: "Invalid temperature or pressure for gas phase v0 calculation.",
        };
      }
      v0_reactorInlet = (totalMolarFlowReactorInlet * R_GAS_CONSTANT * T_K) / P_total_kPa;
    }
    
    if (i % 10 === 0 || i < 5) {
      console.log('💧 Flow Calculations:', {
        totalMolarFlowReactorInlet,
        v0_reactorInlet,
        phase: reactionPhase
      });
    }
    
    // Check for zero volume or zero flow conditions before solving reactor
    if (reactorVolume <= 0 || (v0_reactorInlet <= 1e-9 && totalMolarFlowReactorInlet > 1e-9)) {
       reactorOutletStream_final = { ...reactorInletStream };
       if (i % 10 === 0 || i < 5) {
         // Zero volume or flow condition - no reaction occurring
       }
    } else if (v0_reactorInlet < 1e-9 && totalMolarFlowReactorInlet < 1e-9) {
       reactorOutletStream_final = createEmptyStream(allComponentNames);
       if (i % 10 === 0 || i < 5) {
         // No flow detected - empty outlet stream
       }
    } else {
        // 3. SOLVE THE REACTOR (SINGLE PASS)
        let singlePassResultFlows: Stream;
        
        if (reactorType === 'CSTR') {
          singlePassResultFlows = await solveCSTRParallel(
            parsedParallelReactions,
            reactorVolume,
            reactorInletStream,
            v0_reactorInlet
          );
        } else { // PFR
          const pfrProfile: PFR_ODE_ProfilePoint[] = await solvePFR_ODE_System(
            parsedParallelReactions,
            reactorVolume,
            reactorInletStream,
            v0_reactorInlet,
            allComponentNames
          );
          singlePassResultFlows = pfrProfile.length > 0 ? pfrProfile[pfrProfile.length - 1].flowRates : { ...reactorInletStream };
        }
        reactorOutletStream_final = singlePassResultFlows;
        
    }

    // 4. SEPARATOR BALANCE (using new logic)
    const { productStream, recycleStream: recycleStream_calculated } = separateStream(
      reactorOutletStream_final,
      allComponentNames,
      recycleFractionPercent // Pass recycle fraction percentage here
    );

    // 5. CHECK FOR CONVERGENCE
    const error = calculateStreamError(recycleStream_calculated, recycleStream_guess, allComponentNames);
    finalError = error;

    // Store iteration history for debugging
    iterationHistory.push({
      iteration: iterationCount,
      recycleStreamError: error,
      recycleStream: { ...recycleStream_calculated },
      reactorOutlet: { ...reactorOutletStream_final }
    });

    reactorInletStream_final = reactorInletStream; // Store the inlet for this iteration
    productStream_final = productStream; // Store product stream for this iteration

    if (error < tolerance) {
      converged = true;
      recycleStream_guess = recycleStream_calculated; // Final converged recycle stream
      break;
    }

    // 6. UPDATE GUESS with convergence acceleration
    const dampingFactor = 0.5; // Use underrelaxation to improve stability
    
    // Apply damping: new_guess = damping * calculated + (1-damping) * old_guess
    const dampedRecycleStream: Stream = {};
    allComponentNames.forEach(name => {
      const calculated = recycleStream_calculated[name] || 0;
      const oldGuess = recycleStream_guess[name] || 0;
      dampedRecycleStream[name] = dampingFactor * calculated + (1 - dampingFactor) * oldGuess;
    });
    
    recycleStream_guess = dampedRecycleStream; // Use damped update
    
    // Add debugging information about convergence issues
    if (iterationCount > 10 && error > tolerance * 1000) {
      const issue = `Large error at iteration ${iterationCount}: ${error.toExponential(3)}`;
      convergenceIssues.push(issue);
      if (iterationCount % 20 === 0) {
        // Suggest reactor improvements
        const firstReactantName = parsedParallelReactions.reactions[0]?.reactants[0]?.name || 'A';
        const reactorInletFlow = reactorInletStream[firstReactantName] || 0;
        const reactorOutletFlow = reactorOutletStream_final[firstReactantName] || 0;
        const singlePassConv = reactorInletFlow > 1e-9 ? (reactorInletFlow - reactorOutletFlow) / reactorInletFlow : 0;
        
        if (singlePassConv < 0.05) {
          // Store suggestions in convergence issues instead of console
          convergenceIssues.push('Single pass conversion is very low (<5%). Consider increasing reactor volume or temperature.');
        }
      }
    }
    
    if (iterationCount > maxIterations * 0.3) { // Check after 30% of max iterations
      const totalReactorInlet = Object.values(reactorInletStream).reduce((sum, f) => sum + f, 0);
      const totalReactorOutlet = Object.values(reactorOutletStream_final).reduce((sum, f) => sum + f, 0);
      const totalRecycle = Object.values(recycleStream_calculated).reduce((sum, f) => sum + f, 0);
      const totalFresh = Object.values(freshFeedRates).reduce((sum, f) => sum + f, 0);
      
      // Check for very low conversion that makes recycle impractical
      if (totalReactorOutlet > totalReactorInlet * 0.98 && totalRecycle > totalFresh * 20) {
        const issue = `System appears unworkable: Very low reactor conversion with excessive recycle buildup`;
        convergenceIssues.push(issue);
        
        // Early termination for hopeless cases
        if (iterationCount > maxIterations * 0.5 && error > tolerance * 1000000) {
          break;
        }
      }
      
      if (totalReactorOutlet < totalReactorInlet * 0.01) {
        const issue = `Very low reactor conversion detected - reactor may be too small or conditions too mild`;
        convergenceIssues.push(issue);
      }
      
      if (totalRecycle > totalReactorInlet * 10) {
        const issue = `Very high recycle ratio detected - system may be unstable`;
        convergenceIssues.push(issue);
      }
    }
  }

  // Set finalError after loop completes
  if (iterationHistory.length > 0) {
    finalError = iterationHistory[iterationHistory.length - 1].recycleStreamError;
  }

  // --- Calculate Overall Conversion (based on fresh feed and final product stream) ---
  let overallConversionResult: { reactantName: string; value: number } | undefined = undefined;
  
  // Create a temporary Component array with fresh feed rates for findLimitingReactant
  const freshFeedComponents = components.map(c => ({
    ...c,
    initialFlowRate: (freshFeedRates[c.name] || 0).toString(),
  }));

  const overallLimitingReactant = findLimitingReactant(parsedParallelReactions, freshFeedComponents);

  if (overallLimitingReactant) {
    const F_A0_fresh = overallLimitingReactant.initialFlow; // This is from freshFeedRates
    const F_A_product = productStream_final[overallLimitingReactant.name] || 0;
    if (F_A0_fresh > 1e-9) { // Avoid division by zero
      overallConversionResult = {
        reactantName: overallLimitingReactant.name,
        value: Math.max(0, Math.min(1, (F_A0_fresh - F_A_product) / F_A0_fresh)),
      };
    } else if (F_A0_fresh <= 1e-9 && F_A_product <= 1e-9) { // No limiting reactant fed, none in product
        overallConversionResult = {
            reactantName: overallLimitingReactant.name,
            value: 0, // Or undefined, as conversion is ill-defined
        };
    }
  }
  
  // --- Calculate Single Pass Conversion (based on reactorInletStream_final and reactorOutletStream_final) ---
    let singlePassConversionResult: { reactantName: string; value: number } | undefined = undefined;
    const reactorInletComponents = components.map(c => ({
        ...c,
        initialFlowRate: (reactorInletStream_final[c.name] || 0).toString(),
    }));
    const singlePassLimitingReactant = findLimitingReactant(parsedParallelReactions, reactorInletComponents);

    if (singlePassLimitingReactant) {
        const F_A0_reactor = singlePassLimitingReactant.initialFlow;
        const F_A_reactor_outlet = reactorOutletStream_final[singlePassLimitingReactant.name] || 0;
        if (F_A0_reactor > 1e-9) {
            singlePassConversionResult = {
                reactantName: singlePassLimitingReactant.name,
                value: Math.max(0, Math.min(1, (F_A0_reactor - F_A_reactor_outlet) / F_A0_reactor)),
            };
        } else if (F_A0_reactor <= 1e-9 && F_A_reactor_outlet <= 1e-9) {
             singlePassConversionResult = {
                reactantName: singlePassLimitingReactant.name,
                value: 0,
            };
        }    }

  return {
    productStream: productStream_final,
    recycleStreamFinal: recycleStream_guess,
    reactorInletStreamFinal: reactorInletStream_final, 
    reactorOutletStreamFinal: reactorOutletStream_final,
    overallConversion: overallConversionResult,
    singlePassConversion: singlePassConversionResult,
    iterations: iterationCount,
    converged,
    error: converged ? undefined : `Recycle loop did not converge after ${iterationCount} iterations.`,
    debugInfo: {
      iterationHistory,
      finalError,
      tolerance,
      maxIterations,
      convergenceIssues
    }
  };
}

