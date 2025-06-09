// Import types from reaction-parser (TODO: fix import issue)
interface Component {
  id: string;
  name: string;
  initialFlowRate: string;
  isReactant?: boolean;
  isProduct?: boolean;
  stoichiometricCoefficient?: number;
  reactionOrders?: { [reactionId: string]: string };
}

interface ParsedComponent extends Component {
  stoichiometricCoefficient: number;
  role: 'reactant' | 'product' | 'inert';
  reactionOrderNum?: number;
  reactionOrderNumBackward?: number;
}

interface ParsedReaction {
  reactants: ParsedComponent[];
  products: ParsedComponent[];
  allInvolved: ParsedComponent[];
  rateConstantAtT?: number;
  rateConstantBackwardAtT?: number;
  isEquilibriumReaction?: boolean;
}

interface ParsedParallelReactions {
  reactions: ParsedReaction[];
  allInvolvedComponents: ParsedComponent[];
  error?: string;
}

export interface CalculationResults {
  conversion?: { reactantName: string; value: number }; // X for a key reactant
  outletFlowRates?: { [componentName: string]: number };
  error?: string;
  // Add selectivity, concentrations etc. later
}

// Calculate combined rate for a specific component across all parallel reactions
export const calculateCombinedRate = (
  componentName: string,
  concentrations: { [name: string]: number },
  parsedParallelReactions: ParsedParallelReactions
): number => {
  console.log(`[DEBUG calculateCombinedRate] Calculating for ${componentName} across ${parsedParallelReactions.reactions.length} reactions`);
  
  let totalRate = 0;

  parsedParallelReactions.reactions.forEach((reaction, index) => {
    // Find the component in this reaction
    const componentInReaction = reaction.allInvolved.find(comp => comp.name === componentName);
    
    if (!componentInReaction) {
      console.log(`[DEBUG] Component ${componentName} not involved in reaction ${index + 1}`);
      return; // Component not involved in this reaction
    }

    console.log(`[DEBUG] Component ${componentName} found in reaction ${index + 1} as ${componentInReaction.role}, stoich=${componentInReaction.stoichiometricCoefficient}`);

    // Calculate rate for this specific reaction
    const { reactants, products, rateConstantAtT, rateConstantBackwardAtT, isEquilibriumReaction } = reaction;
    
    console.log(`[DEBUG] Reaction ${index + 1} rate constants: forward=${rateConstantAtT}, backward=${rateConstantBackwardAtT}`);
    
    if (rateConstantAtT === undefined) {
      return; // Skip if no forward rate constant
    }

    // Calculate forward rate
    let forwardRate = rateConstantAtT;
    reactants.forEach(reactant => {
      const conc = concentrations[reactant.name] || 0;
      const order = reactant.reactionOrderNum !== undefined ? reactant.reactionOrderNum : 1;
      forwardRate *= Math.pow(Math.max(0, conc), order);
    });

    // Calculate backward rate (if equilibrium)
    let backwardRate = 0;
    if (isEquilibriumReaction && rateConstantBackwardAtT !== undefined && rateConstantBackwardAtT >= 0) {
      backwardRate = rateConstantBackwardAtT;
      products.forEach(product => {
        const conc = concentrations[product.name] || 0;
        const order = product.reactionOrderNumBackward !== undefined ? product.reactionOrderNumBackward : 1;
        backwardRate *= Math.pow(Math.max(0, conc), order);
      });
    }

    // Net rate for this reaction
    const netReactionRate = forwardRate - backwardRate;
    
    // Contribution to this component's rate based on stoichiometry
    let stoichCoeff = 0;
    if (componentInReaction.role === 'reactant') {
      stoichCoeff = -componentInReaction.stoichiometricCoefficient; // Negative for consumption
    } else if (componentInReaction.role === 'product') {
      stoichCoeff = componentInReaction.stoichiometricCoefficient; // Positive for production
    }
    
    totalRate += stoichCoeff * netReactionRate;
  });

  return totalRate;
};

// CSTR solver for parallel reactions: V = F_A0 * X / (-r_A_total)
export const solveCSTRParallel = (
  parsedParallelReactions: ParsedParallelReactions,
  F_A0: number, // Flow rate of key reactant
  C_A0: number, // Initial concentration of key reactant
  V: number,
  initialConcentrations: { [name: string]: number }, // Initial concentrations of ALL components
  keyReactantName: string,
  allInvolvedComponents: ParsedComponent[]
): number => {
  let X_guess = 0.5; // Initial guess for conversion
  const maxIter = 100;
  let iter = 0;
  const tol = 1e-6;

  const C_key0 = initialConcentrations[keyReactantName];
  
  // Find the key reactant's total stoichiometric coefficient across all reactions
  let totalKeyReactantStoich = 0;
  parsedParallelReactions.reactions.forEach(reaction => {
    const keyInReaction = reaction.allInvolved.find(comp => comp.name === keyReactantName);
    if (keyInReaction && keyInReaction.role === 'reactant') {
      totalKeyReactantStoich += keyInReaction.stoichiometricCoefficient;
    }
  });

  if (totalKeyReactantStoich === 0) {
    return 0; // Key reactant not involved in any reaction
  }

  while (iter < maxIter) {
    const currentConcentrations: { [name: string]: number } = {};
    
    // Calculate concentrations at exit based on current conversion guess
    allInvolvedComponents.forEach(comp => {
      const C_i0 = initialConcentrations[comp.name];
      let currentConc: number;

      if (comp.name === keyReactantName) {
        currentConc = C_key0 * (1 - X_guess);
      } else {
        // For parallel reactions, need to sum contributions from all reactions
        let totalChange = 0;
        
        parsedParallelReactions.reactions.forEach(reaction => {
          const compInReaction = reaction.allInvolved.find(c => c.name === comp.name);
          const keyInReaction = reaction.allInvolved.find(c => c.name === keyReactantName);
          
          if (compInReaction && keyInReaction && keyInReaction.role === 'reactant') {
            const nu_i_signed = compInReaction.role === 'reactant' 
              ? -compInReaction.stoichiometricCoefficient 
              : compInReaction.stoichiometricCoefficient;
            const nu_key_abs = keyInReaction.stoichiometricCoefficient;
            
            // Assume equal conversion across all reactions for simplification
            totalChange += (nu_i_signed / nu_key_abs) * C_key0 * X_guess;
          }
        });
        
        currentConc = C_i0 + totalChange;
      }
      currentConcentrations[comp.name] = Math.max(0, currentConc);
    });

    // Calculate combined rate for the key reactant
    const r_A_total = calculateCombinedRate(keyReactantName, currentConcentrations, parsedParallelReactions);
    
    if (r_A_total >= 0) break; // Rate should be negative for consumption

    const X_calc = (-r_A_total * V) / F_A0;
    
    if (Math.abs(X_calc - X_guess) < tol) {
      X_guess = X_calc;
      break;
    }
    X_guess = (X_calc + X_guess) / 2;
    iter++;
  }
  
  return Math.max(0, Math.min(X_guess, 1.0));
};

// PFR solver for parallel reactions: dX/dV = -r_A_total / F_A0
export const solvePFRParallel = (
  parsedParallelReactions: ParsedParallelReactions,
  F_A0: number, // Flow rate of key reactant
  C_A0: number, // Initial concentration of key reactant
  V: number,
  initialConcentrations: { [name: string]: number }, // Initial concentrations of ALL components
  keyReactantName: string,
  allInvolvedComponents: ParsedComponent[]
): number => {
  const steps = 100;
  const dV = V / steps;
  let X_current = 0;

  const C_key0 = initialConcentrations[keyReactantName];
  
  // Find the key reactant's total stoichiometric coefficient across all reactions
  let totalKeyReactantStoich = 0;
  parsedParallelReactions.reactions.forEach(reaction => {
    const keyInReaction = reaction.allInvolved.find(comp => comp.name === keyReactantName);
    if (keyInReaction && keyInReaction.role === 'reactant') {
      totalKeyReactantStoich += keyInReaction.stoichiometricCoefficient;
    }
  });

  if (totalKeyReactantStoich === 0) {
    return 0; // Key reactant not involved in any reaction
  }

  for (let i = 0; i < steps; i++) {
    const currentConcentrations: { [name: string]: number } = {};
    
    allInvolvedComponents.forEach(comp => {
      const C_i0 = initialConcentrations[comp.name];
      let currentConc: number;

      if (comp.name === keyReactantName) {
        currentConc = C_key0 * (1 - X_current);
      } else {
        // For parallel reactions, need to sum contributions from all reactions
        let totalChange = 0;
        
        parsedParallelReactions.reactions.forEach(reaction => {
          const compInReaction = reaction.allInvolved.find(c => c.name === comp.name);
          const keyInReaction = reaction.allInvolved.find(c => c.name === keyReactantName);
          
          if (compInReaction && keyInReaction && keyInReaction.role === 'reactant') {
            const nu_i_signed = compInReaction.role === 'reactant' 
              ? -compInReaction.stoichiometricCoefficient 
              : compInReaction.stoichiometricCoefficient;
            const nu_key_abs = keyInReaction.stoichiometricCoefficient;
            
            // Assume equal conversion across all reactions for simplification
            totalChange += (nu_i_signed / nu_key_abs) * C_key0 * X_current;
          }
        });
        
        currentConc = C_i0 + totalChange;
      }
      currentConcentrations[comp.name] = Math.max(0, currentConc);
    });

    // Calculate combined rate for the key reactant
    const r_A_total = calculateCombinedRate(keyReactantName, currentConcentrations, parsedParallelReactions);
    
    if (r_A_total >= 0) break; // Rate should be negative for consumption
    
    const dX = (-r_A_total / F_A0) * dV;
    X_current += dX;
    
    if (X_current >= 1.0) {
      X_current = 1.0;
      break;
    }
  }
  
  return Math.max(0, Math.min(X_current, 1.0));
};

// Calculate outlet flow rates for parallel reactions
export const calculateOutletFlowRatesParallel = (
  components: Component[],
  parsedParallelReactions: ParsedParallelReactions,
  keyReactantName: string,
  F_key0: number, // Initial molar flow rate of key reactant
  conversion: number, // Overall conversion of key reactant from solver
  reactorVolume: number,
  volumetricFlowRate: number, // Total volumetric flow rate (v0)
  reactorType: 'CSTR' | 'PFR'
): { [componentName: string]: number } => {
  console.log(`[DEBUG calculateOutletFlowRatesParallel] Start. KeyR: ${keyReactantName}, X: ${conversion}, V: ${reactorVolume}, v0: ${volumetricFlowRate}, Type: ${reactorType}`);
  const outletFlowRates: { [componentName: string]: number } = {};
  const numReactions = parsedParallelReactions.reactions.length;
  const reactionExtents: number[] = new Array(numReactions).fill(0);

  // Pass 1: Calculate all reaction extents sequentially
  parsedParallelReactions.reactions.forEach((reaction, reactionIndex) => {
    let currentReactionExtent = 0;
    const keyInThisReaction = reaction.allInvolved.find(comp => comp.name === keyReactantName);

    if (keyInThisReaction && keyInThisReaction.role === 'reactant') {
      // This reaction involves the key reactant.
      // Note: If keyReactant is in multiple reactions, 'conversion' is overall.
      // This calculation assumes this reaction takes a proportional share or is the main one.
      const nu_key_abs = keyInThisReaction.stoichiometricCoefficient;
      if (nu_key_abs > 0) {
        currentReactionExtent = (F_key0 * conversion) / nu_key_abs;
      } else {
        currentReactionExtent = 0; // Should not happen for a reactant
      }
      console.log(`[DEBUG Extent Pass] Primary reaction ${reactionIndex + 1} (key reactant ${keyReactantName} involved), nu_key: ${nu_key_abs}, extent: ${currentReactionExtent}`);
    } else {
      // Secondary reaction (does not involve the keyReactantName directly)
      // Calculate extent based on available reactants from prior reactions and its own kinetics.
      console.log(`[DEBUG Extent Pass] Secondary reaction ${reactionIndex + 1} (key reactant ${keyReactantName} NOT involved)`);
      let minAvailableExtentForThisReaction = Infinity;

      if (reaction.reactants.length === 0) {
        currentReactionExtent = 0; // No reactants, extent is 0
        console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1} has no reactants. Extent = 0.`);
      } else {
        reaction.reactants.forEach(reactant => {
          const reactantCompInfo = components.find(comp => comp.name === reactant.name);
          let availableFlowOfReactant = parseFloat(reactantCompInfo?.initialFlowRate || '0');

          // Add production of this reactant from all previous reactions
          for (let prevReactionIdx = 0; prevReactionIdx < reactionIndex; prevReactionIdx++) {
            const prevReaction = parsedParallelReactions.reactions[prevReactionIdx];
            const reactantAsProductInPrev = prevReaction.allInvolved.find(
              p => p.name === reactant.name && p.role === 'product'
            );
            if (reactantAsProductInPrev) {
              const productionFromPrevReaction = reactantAsProductInPrev.stoichiometricCoefficient * reactionExtents[prevReactionIdx];
              availableFlowOfReactant += productionFromPrevReaction;
              console.log(`[DEBUG Extent Pass] Reactant ${reactant.name} for reaction ${reactionIndex + 1}: +${productionFromPrevReaction} from reaction ${prevReactionIdx + 1} (extent ${reactionExtents[prevReactionIdx]})`);
            }
          }
          
          availableFlowOfReactant = Math.max(0, availableFlowOfReactant); // Ensure non-negative

          if (reactant.stoichiometricCoefficient > 0) {
            const maxExtentFromThisReactant = availableFlowOfReactant / reactant.stoichiometricCoefficient;
            minAvailableExtentForThisReaction = Math.min(minAvailableExtentForThisReaction, maxExtentFromThisReactant);
            console.log(`[DEBUG Extent Pass] Reactant ${reactant.name} for reaction ${reactionIndex + 1}: initial=${reactantCompInfo?.initialFlowRate || '0'}, total available=${availableFlowOfReactant}, stoichCoeff=${reactant.stoichiometricCoefficient}, max_extent_contrib=${maxExtentFromThisReactant}`);
          } else {
            // Stoichiometric coefficient is zero or negative, this reactant cannot limit extent.
             console.log(`[DEBUG Extent Pass] Reactant ${reactant.name} for reaction ${reactionIndex + 1} has non-positive stoich coeff ${reactant.stoichiometricCoefficient}.`);
          }
        });

        if (minAvailableExtentForThisReaction === Infinity || minAvailableExtentForThisReaction < 0) {
             minAvailableExtentForThisReaction = 0; // If no limiting reactant found or negative, extent is 0
        }
        
        console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1}: minAvailableExtent (thermodynamic max) = ${minAvailableExtentForThisReaction}`);

        // Apply kinetic adjustment factor
        const k_secondary = reaction.rateConstantAtT; // Forward rate constant of this secondary reaction
        let kineticAdjustmentFactor = 0.5; // Default fallback

        if (k_secondary !== undefined && k_secondary >= 0 && volumetricFlowRate > 0 && reactorVolume >= 0) {
          const tau = reactorVolume / volumetricFlowRate;
          
          if (!isFinite(tau) && tau > 0 && k_secondary > 0) { // Effectively infinite time (e.g. v0 -> 0, V > 0)
            kineticAdjustmentFactor = 1.0;
          } else if (tau === 0 || k_secondary === 0) { // No time to react or no rate constant
            kineticAdjustmentFactor = 0.0;
          } else if (isFinite(tau)) {
             if (reactorType === 'CSTR') {
                kineticAdjustmentFactor = (k_secondary * tau) / (1 + k_secondary * tau);
             } else if (reactorType === 'PFR') {
                kineticAdjustmentFactor = 1 - Math.exp(-k_secondary * tau);
             }
          } else {
             kineticAdjustmentFactor = 0.0; // Should not happen if v0 > 0, V >= 0
          }
          kineticAdjustmentFactor = Math.max(0, Math.min(kineticAdjustmentFactor, 1)); // Clamp to [0,1]
          console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1}: k_secondary=${k_secondary}, tau=${tau.toFixed(3)}, kineticFactor=${kineticAdjustmentFactor.toFixed(3)}`);
        } else {
            console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1}: Using fallback kineticFactor=0.5 (k_secondary=${k_secondary}, v0=${volumetricFlowRate}, V=${reactorVolume})`);
        }
        
        currentReactionExtent = Math.max(0, minAvailableExtentForThisReaction * kineticAdjustmentFactor);
      }
    }
    reactionExtents[reactionIndex] = currentReactionExtent;
    console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1} final extent: ${currentReactionExtent}`);
  });

  // Pass 2: Calculate outlet flow rates using the calculated extents
  components.forEach(c => {
    let F_i_outlet = parseFloat(c.initialFlowRate) || 0;
    parsedParallelReactions.reactions.forEach((reaction, reactionIndex) => {
      const compInReaction = reaction.allInvolved.find(comp => comp.name === c.name);
      if (compInReaction) {
        const nu_i_signed = compInReaction.role === 'reactant' 
          ? -compInReaction.stoichiometricCoefficient 
          : compInReaction.stoichiometricCoefficient;
        F_i_outlet += nu_i_signed * reactionExtents[reactionIndex];
      }
    });
    outletFlowRates[c.name] = Math.max(0, F_i_outlet); // Ensure non-negative flow rate
    console.log(`[DEBUG OutletCalc] Component ${c.name}: initial_flow=${c.initialFlowRate || '0'}, final_outlet_flow=${outletFlowRates[c.name].toFixed(4)}`);
  });

  console.log("[DEBUG calculateOutletFlowRatesParallel] Finished.", outletFlowRates);
  return outletFlowRates;
};

// Function to determine the limiting reactant for the primary reaction
// For parallel reactions, we consider the first reaction as the primary one
export const findLimitingReactant = (
  parsedParallelReactions: ParsedParallelReactions,
  components: Component[]
): { name: string; initialFlow: number } | null => {
  if (parsedParallelReactions.reactions.length === 0) return null;
  
  // Use the first reaction to determine limiting reactant
  const primaryReaction = parsedParallelReactions.reactions[0];
  
  if (primaryReaction.reactants.length === 0) return null;
  
  // Find the reactant with the smallest (initial flow / stoichiometric coefficient) ratio
  let limitingReactant: { name: string; initialFlow: number } | null = null;
  let minRatio = Infinity;
  
  primaryReaction.reactants.forEach(reactant => {
    const component = components.find(c => c.name === reactant.name);
    const initialFlow = parseFloat(component?.initialFlowRate || '0');
    
    if (initialFlow > 0 && reactant.stoichiometricCoefficient > 0) {
      const ratio = initialFlow / reactant.stoichiometricCoefficient;
      if (ratio < minRatio) {
        minRatio = ratio;
        limitingReactant = { name: reactant.name, initialFlow };
      }
    }
  });
  
  return limitingReactant;
};

// Keep original functions for backward compatibility
export const solveCSTR = (
  rateLaw: (concentrations: { [name: string]: number }) => number,
  F_A0: number,
  C_A0: number,
  V: number,
  initialConcentrations: { [name: string]: number },
  keyReactant: ParsedComponent,
  allInvolvedComponents: ParsedComponent[]
): number => {
  // This is kept for backward compatibility, but should be replaced with parallel version
  let X_guess = 0.5;
  const maxIter = 100;
  let iter = 0;
  const tol = 1e-6;

  const C_key0 = initialConcentrations[keyReactant.name];
  const nu_key_abs = keyReactant.stoichiometricCoefficient;

  while (iter < maxIter) {
    const currentConcentrations: { [name: string]: number } = {};
    
    allInvolvedComponents.forEach(comp => {
      const C_i0 = initialConcentrations[comp.name];
      let currentConc: number;

      if (comp.name === keyReactant.name) {
        currentConc = C_key0 * (1 - X_guess);
      } else if (comp.role === 'reactant') {
        const nu_i_signed = -comp.stoichiometricCoefficient;
        currentConc = C_i0 + (nu_i_signed / nu_key_abs) * C_key0 * X_guess;
      } else if (comp.role === 'product') {
        const nu_i_signed = comp.stoichiometricCoefficient;
        currentConc = C_i0 + (nu_i_signed / nu_key_abs) * C_key0 * X_guess;
      } else {
        currentConc = C_i0 !== undefined ? C_i0 : 0;
      }
      currentConcentrations[comp.name] = Math.max(0, currentConc);
    });

    const r_A = rateLaw(currentConcentrations);
    if (r_A <= 0) break;

    const X_calc = (r_A * V) / F_A0;
    
    if (Math.abs(X_calc - X_guess) < tol) {
      X_guess = X_calc;
      break;
    }
    X_guess = (X_calc + X_guess) / 2;
    iter++;
  }
  
  return Math.max(0, Math.min(X_guess, 1.0));
};

export const solvePFR = (
  rateLaw: (concentrations: { [name: string]: number }) => number,
  F_A0: number,
  C_A0: number,
  V: number,
  initialConcentrations: { [name: string]: number },
  keyReactant: ParsedComponent,
  allInvolvedComponents: ParsedComponent[]
): number => {
  // This is kept for backward compatibility, but should be replaced with parallel version
  const steps = 100;
  const dV = V / steps;
  let X_current = 0;

  const C_key0 = initialConcentrations[keyReactant.name];
  const nu_key_abs = keyReactant.stoichiometricCoefficient;

  for (let i = 0; i < steps; i++) {
    const currentConcentrations: { [name: string]: number } = {};
    
    allInvolvedComponents.forEach(comp => {
      const C_i0 = initialConcentrations[comp.name];
      let currentConc: number;

      if (comp.name === keyReactant.name) {
        currentConc = C_key0 * (1 - X_current);
      } else if (comp.role === 'reactant') {
        const nu_i_signed = -comp.stoichiometricCoefficient;
        currentConc = C_i0 + (nu_i_signed / nu_key_abs) * C_key0 * X_current;
      } else if (comp.role === 'product') {
        const nu_i_signed = comp.stoichiometricCoefficient;
        currentConc = C_i0 + (nu_i_signed / nu_key_abs) * C_key0 * X_current;
      } else {
        currentConc = C_i0 !== undefined ? C_i0 : 0;
      }
      currentConcentrations[comp.name] = Math.max(0, currentConc);
    });

    const r_A = rateLaw(currentConcentrations);
    if (r_A <= 0) break;
    
    const dX = (r_A / F_A0) * dV;
    X_current += dX;
    
    if (X_current >= 1.0) {
      X_current = 1.0;
      break;
    }
  }
  
  return Math.max(0, Math.min(X_current, 1.0));
};
