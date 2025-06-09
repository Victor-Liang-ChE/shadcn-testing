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
  let totalRate = 0;
  parsedParallelReactions.reactions.forEach((reaction, reactionIndex) => { // Added reactionIndex here
    // Find the component in this reaction
    const componentInReaction = reaction.allInvolved.find(comp => comp.name === componentName);

    if (!componentInReaction) {
      return; // Component not involved in this reaction
    }

    // Calculate rate for this specific reaction
    const { reactants, products, rateConstantAtT, rateConstantBackwardAtT, isEquilibriumReaction } = reaction;

    if (rateConstantAtT === undefined) {
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex}: Skipped due to undefined forward rate constant.`);
      return; // Skip if no forward rate constant
    }

    // Calculate forward rate
    let forwardRate = rateConstantAtT;
    console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex}: Initial forwardRate (k_fwd) = ${forwardRate.toExponential(3)}`);
    reactants.forEach(reactant => {
      const conc = concentrations[reactant.name] || 0;
      const order = reactant.reactionOrderNum !== undefined ? reactant.reactionOrderNum : 1;
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} - Reactant ${reactant.name}: Conc=${conc.toFixed(4)}, Order_fwd=${order}, forwardRate_before_mult=${forwardRate.toExponential(3)}`);
      forwardRate *= Math.pow(Math.max(0, conc), order);
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} - Reactant ${reactant.name}: forwardRate_after_mult=${forwardRate.toExponential(3)}`);
    });
    console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex}: Final calculated forwardRate = ${forwardRate.toExponential(3)}`);

    // Calculate backward rate (if equilibrium)
    let backwardRate = 0;
    if (isEquilibriumReaction && rateConstantBackwardAtT !== undefined && rateConstantBackwardAtT >= 0) {
      backwardRate = rateConstantBackwardAtT; // Initialize with k_bwd
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} (Equilibrium): Initial backwardRate (k_bwd) = ${backwardRate.toExponential(3)}`);
      products.forEach(product => {
        const conc = concentrations[product.name] || 0;
        const order = product.reactionOrderNumBackward !== undefined ? product.reactionOrderNumBackward : 1;
        // Ensure product.reactionOrderNumBackward is logged if it exists
        const orderSource = product.reactionOrderNumBackward !== undefined ? `explicit (${product.reactionOrderNumBackward})` : `default (1)`;
        console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} - Product ${product.name} (for backward rate): Conc=${conc.toFixed(4)}, Order_bwd=${order} (source: ${orderSource}), k_bwd_term_before_mult=${backwardRate.toExponential(3)}`);
        if (conc === 0 && order > 0) { // if concentration is zero, this term will make backwardRate zero unless order is 0
             backwardRate = 0; // Optimization: if any product conc is 0 and its order > 0, that term is 0, so overall product of terms is 0.
             console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} - Product ${product.name}: Conc is 0 and order > 0, setting backwardRate to 0.`);
        } else {
            backwardRate *= Math.pow(Math.max(0, conc), order);
        }
        console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} - Product ${product.name}: backwardRate_after_product_term=${backwardRate.toExponential(3)}`);
      });
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} (Equilibrium): Final calculated backwardRate = ${backwardRate.toExponential(3)}`);
    } else if (isEquilibriumReaction) {
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex} (Equilibrium): Backward rate not calculated. k_bwd=${rateConstantBackwardAtT}, isEquilibriumReaction=${isEquilibriumReaction}`);
    } else {
      console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex}: Not an equilibrium reaction. Backward rate is 0.`);
    }

    // Net rate for this reaction
    const netReactionRate = forwardRate - backwardRate;
    console.log(`[DEBUG Solver RateCalc] Reaction Index ${reactionIndex}: NetReactionRate (forward - backward) = ${netReactionRate.toExponential(3)} (Forward: ${forwardRate.toExponential(3)}, Backward: ${backwardRate.toExponential(3)})`);

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
    } else {
      // Secondary reaction (does not involve the keyReactantName directly)
      // Calculate extent based on available reactants from prior reactions and its own kinetics.
      let minAvailableExtentForThisReaction = Infinity;

      if (reaction.reactants.length === 0) {
        currentReactionExtent = 0; // No reactants, extent is 0
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
            }
          }

          availableFlowOfReactant = Math.max(0, availableFlowOfReactant); // Ensure non-negative

          if (reactant.stoichiometricCoefficient > 0) {
            const maxExtentFromThisReactant = availableFlowOfReactant / reactant.stoichiometricCoefficient;
            minAvailableExtentForThisReaction = Math.min(minAvailableExtentForThisReaction, maxExtentFromThisReactant);
          } else {
            // Stoichiometric coefficient is zero or negative, this reactant cannot limit extent.
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
          let effective_k_for_adjustment = k_secondary;

          // If the secondary reaction has reactants
          if (reaction.reactants.length > 0) {
            // Consider the first reactant of the secondary reaction for order adjustment
            // This is a simplification, especially if the secondary reaction has multiple reactants.
            const mainReactantOfSecondary = reaction.reactants[0];
            const order_n = mainReactantOfSecondary.reactionOrderNum !== undefined ? mainReactantOfSecondary.reactionOrderNum : 1;

            if (order_n !== 1) {
              // Estimate the concentration of this mainReactantOfSecondary at the "inlet" of this reaction step
              const reactantCompInfoGlobal = components.find(comp => comp.name === mainReactantOfSecondary.name);
              let F_reactant_in_step = parseFloat(reactantCompInfoGlobal?.initialFlowRate || '0');

              // Sum changes from all previously calculated extents in this pass
              for (let prevReactionIdx = 0; prevReactionIdx < reactionIndex; prevReactionIdx++) {
                const prevReaction = parsedParallelReactions.reactions[prevReactionIdx];
                const compInPrevReaction = prevReaction.allInvolved.find(comp => comp.name === mainReactantOfSecondary.name);
                if (compInPrevReaction) {
                  const nu_signed = compInPrevReaction.role === 'reactant' 
                    ? -compInPrevReaction.stoichiometricCoefficient 
                    : compInPrevReaction.stoichiometricCoefficient;
                  F_reactant_in_step += nu_signed * reactionExtents[prevReactionIdx]; // Use finalized extents from this pass
                }
              }
              F_reactant_in_step = Math.max(0, F_reactant_in_step);
              const C_reactant_in_step = F_reactant_in_step / volumetricFlowRate;

              if (C_reactant_in_step > 1e-9) { // Avoid issues with zero concentration
                effective_k_for_adjustment = k_secondary * Math.pow(C_reactant_in_step, order_n - 1);
                console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1} (${mainReactantOfSecondary.name}): order_n=${order_n}, C_in_step=${C_reactant_in_step.toFixed(4)}, k_orig=${k_secondary.toExponential(3)}, effective_k=${effective_k_for_adjustment.toExponential(3)}`);
              } else if (order_n > 1) { // If C_in is zero and order > 1, rate is effectively zero
                effective_k_for_adjustment = 0;
                console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1} (${mainReactantOfSecondary.name}): order_n=${order_n}, C_in_step=${C_reactant_in_step.toFixed(4)}, effective_k set to 0`);
              } else { // order_n = 1 or (order_n < 1 and C_in_step is zero - potentially problematic)
                // Default to k_secondary if C_in_step is zero and order < 1 to avoid Math.pow errors, or if order is 1.
                console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1} (${mainReactantOfSecondary.name}): order_n=${order_n}, C_in_step=${C_reactant_in_step.toFixed(4)}, using k_orig=${k_secondary.toExponential(3)} as effective_k`);
              }
            }
          }

          const dimensionlessGroup = effective_k_for_adjustment * tau;
          // Simplified CSTR-like model for adjustment factor
          kineticAdjustmentFactor = dimensionlessGroup / (1 + dimensionlessGroup); 
          if (reactorType === 'PFR') {
             // Simplified PFR-like model for adjustment factor (for 1st order effective kinetics)
             kineticAdjustmentFactor = 1 - Math.exp(-dimensionlessGroup);
          }
          kineticAdjustmentFactor = Math.max(0, Math.min(kineticAdjustmentFactor, 1)); // Clamp between 0 and 1
          console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1}: k_secondary=${k_secondary.toExponential(3)}, tau=${tau.toFixed(2)}, effective_k_adj=${effective_k_for_adjustment.toExponential(3)}, dimGroup=${dimensionlessGroup.toFixed(3)}, kineticFactor=${kineticAdjustmentFactor.toFixed(3)}`);
        } else {
            // Fallback if k_secondary is not available or other params are invalid
            kineticAdjustmentFactor = 0.0; 
            console.log(`[DEBUG Extent Pass] Reaction ${reactionIndex + 1}: Fallback kineticFactor = 0.0 (k_secondary issue or invalid params)`);
        }
        
        currentReactionExtent = Math.max(0, minAvailableExtentForThisReaction * kineticAdjustmentFactor);
      }
    }
    reactionExtents[reactionIndex] = currentReactionExtent;
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
