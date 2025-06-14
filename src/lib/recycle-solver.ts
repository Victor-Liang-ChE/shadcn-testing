import { ParsedReaction as Reaction, Component } from './reaction-parser';

export interface RecycleCalculationResults {
  freshFeedFlows: { [key: string]: number };
  outletFlows: { [key: string]: number };
  recycleFlows: { [key: string]: number };
}

/**
 * Solves the overall material balances for a process with recycle.
 * Based on the principles from Doherty & Malone, "Conceptual Process Design".
 *
 * @param reactions - The array of parsed chemical reactions.
 * @param components - The array of all components in the system.
 * @param reactorOutletFlows - The molar flow rates of each component LEAVING the reactor.
 * @param reactorConversion - The per-pass conversion of the limiting reactant in the reactor (a value between 0 and 1).
 * @param desiredProduct - An object containing the name and desired production rate (outlet flow) of the main product.
 * @param limitingReactantName - The name of the limiting reactant for the overall process.
 * @returns An object containing the calculated fresh feed, final outlet, and recycle flow rates.
 */
export function solveRecycleSystem(
  reactions: Reaction[],
  components: Component[],
  reactorOutletFlows: { [key: string]: number },
  reactorConversion: number,
  desiredProduct: { name: string; flowRate: number },
  limitingReactantName: string
): RecycleCalculationResults {
  
  if (!limitingReactantName || !desiredProduct.name) {
    throw new Error("Limiting reactant and desired product must be defined for recycle calculations.");
  }

  // Helper: Create an augmented reaction object that includes the stoichiometry map
  // that the rest of the function's logic expects. This adapts the ParsedReaction
  // structure (with 'reactants' and 'products' arrays) to the expected format.
  const reactionsWithStoichiometry = reactions.map(r => {
    const stoichiometry: { [key: string]: number } = {};
    
    // Reactants have negative coefficients
    r.reactants.forEach(reactant => {
      stoichiometry[reactant.name] = -reactant.stoichiometricCoefficient;
    });

    // Products have positive coefficients
    r.products.forEach(product => {
      stoichiometry[product.name] = product.stoichiometricCoefficient;
    });

    return { ...r, stoichiometry };
  });

  const freshFeedFlows: { [key: string]: number } = {};
  const outletFlows: { [key: string]: number } = {};
  const recycleFlows: { [key: string]: number } = {};

  // --- Step 1: Determine total consumption of Limiting Reactant (LR) based on reactor performance ---
  const F_LR_out_reactor = reactorOutletFlows[limitingReactantName] || 0;
  const F_LR_in_reactor = F_LR_out_reactor / (1 - reactorConversion);
  
  const moles_LR_consumed = F_LR_in_reactor - F_LR_out_reactor;

  if (moles_LR_consumed < 1e-9) {
    // If no reaction is happening, there's no recycle needed.
    components.forEach(c => {
        freshFeedFlows[c.name] = reactorOutletFlows[c.name] || 0;
        outletFlows[c.name] = reactorOutletFlows[c.name] || 0;
        recycleFlows[c.name] = 0;
    });
    return { freshFeedFlows, outletFlows, recycleFlows };
  }
  
  // --- Step 2: Calculate Overall Process Selectivity for all products ---
  const productSelectivities: { [key: string]: number } = {};
  components.forEach(c => {
    const isProduct = reactionsWithStoichiometry.some(r => (r.stoichiometry[c.name] || 0) > 0);
    if (isProduct) {
        // Assuming no product in feed to reactor.
        // For products, the amount formed is simply their outlet flow rate from the reactor.
        const product_formed = (reactorOutletFlows[c.name] || 0); 
        productSelectivities[c.name] = product_formed / moles_LR_consumed;
    }
  });

  // --- Step 3: Calculate Global IO (Input-Output) Balances ---
  const mainProductSelectivity = productSelectivities[desiredProduct.name];
  if (mainProductSelectivity === undefined || mainProductSelectivity < 1e-9) {
      throw new Error(`Selectivity for desired product "${desiredProduct.name}" is zero or undefined. Cannot calculate recycle system.`);
  }
  const F_LR_fresh = desiredProduct.flowRate / mainProductSelectivity;
  freshFeedFlows[limitingReactantName] = F_LR_fresh;
  
  // Calculate outlet flows for all byproducts and fresh feeds for other reactants
  components.forEach(c => {
    if (c.name === limitingReactantName) return;

    const isReactant = reactionsWithStoichiometry.some(r => (r.stoichiometry[c.name] || 0) < 0);
    const isProduct = reactionsWithStoichiometry.some(r => (r.stoichiometry[c.name] || 0) > 0);

    if (isProduct) {
        // Byproduct outlet flow is based on its selectivity relative to the LR fresh feed
        outletFlows[c.name] = F_LR_fresh * (productSelectivities[c.name] || 0);
        freshFeedFlows[c.name] = 0; // Products are not fed
    } else if (isReactant) {
        // Assume other reactants are fed stoichiometrically relative to the limiting reactant
        const totalStoichRatio = reactionsWithStoichiometry.reduce((sum, r) => {
            const reactantStoich = r.stoichiometry[c.name] || 0;
            const limitingStoich = r.stoichiometry[limitingReactantName] || -1;
            // Ensure we don't divide by zero if the limiting reactant isn't in this specific reaction
            if (limitingStoich === 0) return sum;
            return sum + (reactantStoich / limitingStoich);
        }, 0) / reactionsWithStoichiometry.length; // Averaging for simplicity across parallel reactions

        freshFeedFlows[c.name] = F_LR_fresh * Math.abs(totalStoichRatio);
        outletFlows[c.name] = 0; // Assume complete recycle for other reactants for simplicity
    } else {
        // Inert components are neither fed nor produced in the outlet (pass through)
        freshFeedFlows[c.name] = 0;
        outletFlows[c.name] = 0;
    }
  });
  outletFlows[desiredProduct.name] = desiredProduct.flowRate;

  // --- Step 4: Calculate Recycle Balances ---
  // Recycle = (Reactor Outlet) - (Process Outlet)
  // This assumes perfect separation where only unreacted reactants and desired product are recycled.
  components.forEach(c => {
      const reactorOut = reactorOutletFlows[c.name] || 0;
      const processOut = outletFlows[c.name] || 0;
      recycleFlows[c.name] = Math.max(0, reactorOut - processOut);
  });

  return { freshFeedFlows, outletFlows, recycleFlows };
}
