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
  F_key0: number,
  conversion: number
): { [componentName: string]: number } => {
  const outletFlowRates: { [componentName: string]: number } = {};

  components.forEach(c => {
    const F_i0 = parseFloat(c.initialFlowRate) || 0;
    let totalChange = 0;

    // Sum changes from all parallel reactions
    parsedParallelReactions.reactions.forEach(reaction => {
      const compInReaction = reaction.allInvolved.find(comp => comp.name === c.name);
      const keyInReaction = reaction.allInvolved.find(comp => comp.name === keyReactantName);

      if (compInReaction && keyInReaction && keyInReaction.role === 'reactant') {
        const nu_i = compInReaction.role === 'reactant' 
          ? -compInReaction.stoichiometricCoefficient 
          : compInReaction.stoichiometricCoefficient;
        const nu_key_abs = keyInReaction.stoichiometricCoefficient;
        
        // Assume equal conversion across all reactions for simplification
        totalChange += (nu_i / nu_key_abs) * F_key0 * conversion;
      }
    });

    outletFlowRates[c.name] = Math.max(0, F_i0 + totalChange);
  });

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
