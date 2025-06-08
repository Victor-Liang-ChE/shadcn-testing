// Helper types
type RateConstantInputMode = 'directK' | 'arrhenius';

export interface Component {
  id: string;
  name: string;
  initialFlowRate: string; // F_i0 (e.g., mol/s)
  // For stoichiometry parsing results:
  isReactant?: boolean;
  isProduct?: boolean;
  stoichiometricCoefficient?: number;
  reactionOrders?: { [reactionId: string]: string }; // Reaction-specific orders, e.g., "1-fwd": "1", "1-rev": "1"
}

export interface Kinetics {
  rateConstantInputMode: RateConstantInputMode;
  kValue: string; // Rate constant k
  AValue: string; // Pre-exponential factor A
  EaValue: string; // Activation energy Ea (e.g., J/mol)
  reactionTempK: string; // Reaction temperature in Kelvin
}

export interface ParsedComponent extends Component {
  stoichiometricCoefficient: number; // Now mandatory after parsing
  role: 'reactant' | 'product' | 'inert'; // Inert if present in feed but not in reaction
  reactionOrderNum?: number; // Parsed reaction order for forward reaction if component is a reactant
  reactionOrderNumBackward?: number; // Parsed reaction order for backward reaction if component is a product (acting as reactant in reverse)
}

export interface ParsedReaction {
  reactants: ParsedComponent[];
  products: ParsedComponent[];
  allInvolved: ParsedComponent[]; // All components from stoichiometry
  rateConstantAtT?: number; // Calculated k at reaction temperature for forward
  rateConstantBackwardAtT?: number; // Calculated k for backward reaction
  isEquilibriumReaction?: boolean; // Flag if the reaction is equilibrium
}

export interface ParsedParallelReactions {
  reactions: ParsedReaction[];
  allInvolvedComponents: ParsedComponent[]; // All unique components across all reactions
  error?: string;
}

export interface ReactionData {
  id: string;
  reactants: string;
  products: string;
  AValue: string;
  EaValue: string;
  AValueBackward?: string;
  EaValueBackward?: string;
  isEquilibrium?: boolean;
}

const R_GAS_CONSTANT = 8.314; // J/(mol*K) or L*kPa/(mol*K) depending on units

// Arrhenius equation: k = A * exp(-Ea / (R*T))
export const calculateRateConstant = (A: number, Ea: number, T: number): number => {
  return A * Math.exp(-Ea / (R_GAS_CONSTANT * T));
};

// Function to parse compounds from reaction stoichiometry
export const parseCompoundsFromReaction = (reaction: string): string[] => {
  // Remove spaces and split by -> to get reactants and products
  const cleanReaction = reaction.replace(/\s/g, '');
  const [reactants, products] = cleanReaction.split('->');
  
  if (!reactants || !products) return [];
  
  // Extract compound names (letters only, ignoring coefficients)
  const reactantCompounds = reactants.split('+').map(r => r.replace(/\d/g, ''));
  const productCompounds = products.split('+').map(p => p.replace(/\d/g, ''));
  
  // Combine and remove duplicates
  return [...new Set([...reactantCompounds, ...productCompounds])].filter(c => c.length > 0);
};

// Stoichiometry Parsing Logic
export const parseStoichiometry = (
  reactionString: string,
  components: Component[],
  kinetics: Kinetics,
  reactions: ReactionData[],
  forCalculation = false
): { parsedReaction: ParsedReaction | null; error: string | null } => {
  if (!reactionString.trim()) {
    return {
      parsedReaction: null,
      error: forCalculation ? "Reaction stoichiometry cannot be empty." : null
    };
  }

  const arrow = reactionString.includes('->') ? '->' : reactionString.includes('<=>') ? '<=>' : null;
  if (!arrow) {
    return {
      parsedReaction: null,
      error: forCalculation ? "Invalid reaction string: missing '->' or '<=>'." : null
    };
  }

  const [reactantsStr, productsStr] = reactionString.split(arrow).map(s => s.trim());
  if (!reactantsStr || !productsStr) {
    return {
      parsedReaction: null,
      error: forCalculation ? "Invalid reaction string: reactants or products side is empty." : null
    };
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
        // Skip components not in feed list instead of throwing error
        return null;
      }
      
      const firstReactionDetails = reactions[0]; // Assuming parsing is based on the first reaction
      let orderNum: number | undefined = undefined;
      let orderNumBackward: number | undefined = undefined;

      if (role === 'reactant') {
        const orderKeyForward = firstReactionDetails?.id ? `${firstReactionDetails.id}-fwd` : undefined;
        const orderStr = (orderKeyForward && feedComponent.reactionOrders?.[orderKeyForward] !== undefined)
          ? feedComponent.reactionOrders[orderKeyForward]
          : (firstReactionDetails?.id && feedComponent.reactionOrders?.[firstReactionDetails.id] !== undefined)
            ? feedComponent.reactionOrders[firstReactionDetails.id] // Fallback to just reactionId for non-equilibrium or older format
            : '1'; // Default to 1
        orderNum = parseFloat(orderStr);
        if (isNaN(orderNum) || orderNum < 0) {
            throw new Error(`Invalid forward reaction order for ${name} ('${orderStr}'). Must be a non-negative number.`);
        }
      } else if (role === 'product' && firstReactionDetails?.isEquilibrium) {
        const orderKeyBackward = firstReactionDetails?.id ? `${firstReactionDetails.id}-rev` : undefined;
        const orderStr = (orderKeyBackward && feedComponent.reactionOrders?.[orderKeyBackward] !== undefined)
          ? feedComponent.reactionOrders[orderKeyBackward]
          : '1'; // Default to 1 for reverse if not specified
        orderNumBackward = parseFloat(orderStr);
        if (isNaN(orderNumBackward) || orderNumBackward < 0) {
            throw new Error(`Invalid backward reaction order for ${name} ('${orderStr}'). Must be a non-negative number.`);
        }
      }

      return {
        ...feedComponent,
        stoichiometricCoefficient: coeff,
        role,
        reactionOrderNum: orderNum,
        reactionOrderNumBackward: orderNumBackward
      };
    }).filter(Boolean) as ParsedComponent[]; // Filter out null values
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
    if (isNaN(T_K) || T_K <= 0) {
      return {
        parsedReaction: null,
        error: forCalculation ? "Reaction temperature must be a positive number." : null
      };
    }

    // Use A and Ea from the first reaction for now (could be extended for multiple reactions)
    const firstReaction = reactions[0];
    const A = parseFloat(firstReaction?.AValue || kinetics.AValue);
    const Ea = parseFloat(firstReaction?.EaValue || kinetics.EaValue);
    
    if (isNaN(A) || isNaN(Ea) || A <= 0 || Ea < 0) {
      return {
        parsedReaction: null,
        error: forCalculation ? "Invalid A or Ea value for Arrhenius equation. A must be positive, Ea must be non-negative." : null
      };
    }
    
    kAtT = calculateRateConstant(A, Ea, T_K);
    
    if (kAtT <= 0 && A > 0) { // k can be 0 if A is 0
      return {
        parsedReaction: null,
        error: forCalculation ? "Forward rate constant (k) must be positive if A > 0." : null
      };
    }

    let kBackwardAtT: number | undefined;
    const isEquilibrium = firstReaction?.isEquilibrium || false;

    if (isEquilibrium) {
      const A_bwd_str = firstReaction?.AValueBackward;
      const Ea_bwd_str = firstReaction?.EaValueBackward;

      // Use "0" as default if string is empty or undefined, consistent with user snippet
      const A_bwd = parseFloat(A_bwd_str || "0");
      const Ea_bwd = parseFloat(Ea_bwd_str || "0");

      if (isNaN(A_bwd) || isNaN(Ea_bwd) || A_bwd < 0 || Ea_bwd < 0) {
          return {
            parsedReaction: null,
            error: forCalculation ? "Invalid A or Ea value for backward Arrhenius equation. A_bwd and Ea_bwd must be non-negative." : null
          };
      } else if (A_bwd > 0) { // Only calculate if A_bwd is positive, otherwise k_bwd is 0
          kBackwardAtT = calculateRateConstant(A_bwd, Ea_bwd, T_K);
          if (kBackwardAtT < 0 && forCalculation) { // Should not happen if A_bwd > 0
               return {
                 parsedReaction: null,
                 error: "Calculated backward rate constant (k_bwd) is negative."
               };
          }
      } else {
          kBackwardAtT = 0; // If A_bwd is 0, k_bwd is 0
      }
    }

    const newParsedData: ParsedReaction = {
      reactants: parsedReactants,
      products: parsedProducts,
      allInvolved: allInvolvedFromFeed,
      rateConstantAtT: kAtT,
      rateConstantBackwardAtT: kBackwardAtT,
      isEquilibriumReaction: isEquilibrium
    };
    
    return { parsedReaction: newParsedData, error: null };
  } catch (e: any) {
    return {
      parsedReaction: null,
      error: forCalculation ? `Stoichiometry/Kinetics Error: ${e.message}` : null
    };
  }
};

// Parse multiple parallel reactions
export const parseParallelReactions = (
  reactionString: string,
  components: Component[],
  kinetics: Kinetics,
  reactions: ReactionData[],
  forCalculation = false
): ParsedParallelReactions => {
  // Check if we have reactions to process
  if (reactions.length === 0) {
    return {
      reactions: [],
      allInvolvedComponents: [],
      error: forCalculation ? "No reactions defined." : undefined
    };
  }

  const parsedReactions: ParsedReaction[] = [];
  const allComponentsMap = new Map<string, ParsedComponent>();

  // Process each reaction in the reactions array
  for (let i = 0; i < reactions.length; i++) {
    const reactionData = reactions[i];
    
    if (!reactionData.reactants || !reactionData.products) {
      if (forCalculation) {
        return {
          reactions: [],
          allInvolvedComponents: [],
          error: `Reaction ${i + 1} is missing reactants or products.`
        };
      }
      continue; // Skip incomplete reactions
    }

    // Build reaction string from the reaction data
    const arrow = reactionData.isEquilibrium ? '<=>' : '->';
    const rxnString = `${reactionData.reactants} ${arrow} ${reactionData.products}`;
    
    // Parse this single reaction with its specific kinetic data
    const { parsedReaction, error } = parseSingleReaction(
      rxnString,
      components,
      kinetics,
      reactionData,
      forCalculation
    );

    if (error && forCalculation) {
      return {
        reactions: [],
        allInvolvedComponents: [],
        error: `Error in reaction ${i + 1}: ${error}`
      };
    }

    if (parsedReaction) {
      parsedReactions.push(parsedReaction);
      
      // Add all involved components to the global map
      parsedReaction.allInvolved.forEach(comp => {
        if (!allComponentsMap.has(comp.name)) {
          allComponentsMap.set(comp.name, comp);
        }
      });
    }
  }

  return {
    reactions: parsedReactions,
    allInvolvedComponents: Array.from(allComponentsMap.values()),
    error: undefined
  };
};

// Parse a single reaction (refactored from original parseStoichiometry)
export const parseSingleReaction = (
  reactionString: string,
  components: Component[],
  kinetics: Kinetics,
  reactionData: ReactionData,
  forCalculation = false
): { parsedReaction: ParsedReaction | null; error: string | null } => {
  if (!reactionString.trim()) {
    return {
      parsedReaction: null,
      error: forCalculation ? "Reaction stoichiometry cannot be empty." : null
    };
  }

  const arrow = reactionString.includes('->') ? '->' : reactionString.includes('<=>') ? '<=>' : null;
  if (!arrow) {
    return {
      parsedReaction: null,
      error: forCalculation ? "Invalid reaction string: missing '->' or '<=>'." : null
    };
  }

  const [reactantsStr, productsStr] = reactionString.split(arrow).map(s => s.trim());
  if (!reactantsStr || !productsStr) {
    return {
      parsedReaction: null,
      error: forCalculation ? "Invalid reaction string: reactants or products side is empty." : null
    };
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
        // Skip components not in feed list instead of throwing error
        return null;
      }
      
      let orderNum: number | undefined = undefined;
      let orderNumBackward: number | undefined = undefined;

      if (role === 'reactant') {
        const orderKeyForward = reactionData?.id ? `${reactionData.id}-fwd` : undefined;
        const orderStr = (orderKeyForward && feedComponent.reactionOrders?.[orderKeyForward] !== undefined)
          ? feedComponent.reactionOrders[orderKeyForward]
          : (reactionData?.id && feedComponent.reactionOrders?.[reactionData.id] !== undefined)
            ? feedComponent.reactionOrders[reactionData.id] // Fallback to just reactionId for non-equilibrium or older format
            : '1'; // Default to 1
        orderNum = parseFloat(orderStr);
        if (isNaN(orderNum) || orderNum < 0) {
            throw new Error(`Invalid forward reaction order for ${name} ('${orderStr}'). Must be a non-negative number.`);
        }
      } else if (role === 'product' && reactionData?.isEquilibrium) {
        const orderKeyBackward = reactionData?.id ? `${reactionData.id}-rev` : undefined;
        const orderStr = (orderKeyBackward && feedComponent.reactionOrders?.[orderKeyBackward] !== undefined)
          ? feedComponent.reactionOrders[orderKeyBackward]
          : '1'; // Default to 1 for reverse if not specified
        orderNumBackward = parseFloat(orderStr);
        if (isNaN(orderNumBackward) || orderNumBackward < 0) {
            throw new Error(`Invalid backward reaction order for ${name} ('${orderStr}'). Must be a non-negative number.`);
        }
      }

      return {
        ...feedComponent,
        stoichiometricCoefficient: coeff,
        role,
        reactionOrderNum: orderNum,
        reactionOrderNumBackward: orderNumBackward
      };
    }).filter(Boolean) as ParsedComponent[]; // Filter out null values
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
    if (isNaN(T_K) || T_K <= 0) {
      return {
        parsedReaction: null,
        error: forCalculation ? "Reaction temperature must be a positive number." : null
      };
    }

    // Use A and Ea from the specific reaction data
    const A = parseFloat(reactionData?.AValue || kinetics.AValue);
    const Ea = parseFloat(reactionData?.EaValue || kinetics.EaValue);
    
    if (isNaN(A) || isNaN(Ea) || A <= 0 || Ea < 0) {
      return {
        parsedReaction: null,
        error: forCalculation ? "Invalid A or Ea value for Arrhenius equation. A must be positive, Ea must be non-negative." : null
      };
    }
    
    kAtT = calculateRateConstant(A, Ea, T_K);
    
    if (kAtT <= 0 && A > 0) { // k can be 0 if A is 0
      return {
        parsedReaction: null,
        error: forCalculation ? "Forward rate constant (k) must be positive if A > 0." : null
      };
    }

    let kBackwardAtT: number | undefined;
    const isEquilibrium = reactionData?.isEquilibrium || false;

    if (isEquilibrium) {
      const A_bwd_str = reactionData?.AValueBackward;
      const Ea_bwd_str = reactionData?.EaValueBackward;

      // Use "0" as default if string is empty or undefined, consistent with user snippet
      const A_bwd = parseFloat(A_bwd_str || "0");
      const Ea_bwd = parseFloat(Ea_bwd_str || "0");

      if (isNaN(A_bwd) || isNaN(Ea_bwd) || A_bwd < 0 || Ea_bwd < 0) {
          return {
            parsedReaction: null,
            error: forCalculation ? "Invalid A or Ea value for backward Arrhenius equation. A_bwd and Ea_bwd must be non-negative." : null
          };
      } else if (A_bwd > 0) { // Only calculate if A_bwd is positive, otherwise k_bwd is 0
          kBackwardAtT = calculateRateConstant(A_bwd, Ea_bwd, T_K);
          if (kBackwardAtT < 0 && forCalculation) { // Should not happen if A_bwd > 0
               return {
                 parsedReaction: null,
                 error: "Calculated backward rate constant (k_bwd) is negative."
               };
          }
      } else {
          kBackwardAtT = 0; // If A_bwd is 0, k_bwd is 0
      }
    }

    const newParsedData: ParsedReaction = {
      reactants: parsedReactants,
      products: parsedProducts,
      allInvolved: allInvolvedFromFeed,
      rateConstantAtT: kAtT,
      rateConstantBackwardAtT: kBackwardAtT,
      isEquilibriumReaction: isEquilibrium
    };
    
    return { parsedReaction: newParsedData, error: null };
  } catch (e: any) {
    return {
      parsedReaction: null,
      error: forCalculation ? `Stoichiometry/Kinetics Error: ${e.message}` : null
    };
  }
};
