// reactor-solver.ts

// --- TYPE DEFINITIONS (Unchanged) ---
export interface Component {
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
export interface ParsedParallelReactions {
  reactions: ParsedReaction[];
  allInvolvedComponents: ParsedComponent[];
  error?: string;
}
export interface CalculationResults {
  conversion?: { reactantName: string; value: number };
  outletFlowRates?: { [componentName: string]: number };
  error?: string;
}
export interface PFR_ODE_ProfilePoint {
  volume: number;
  flowRates: { [componentName: string]: number };
}

// --- RATE CALCULATION (Unchanged) ---
export const calculateCombinedRate = (
  componentName: string,
  concentrations: { [name: string]: number },
  parsedParallelReactions: ParsedParallelReactions
): number => {
  let totalRate = 0;
  const MIN_CONC = 1e-12; // Numerical floor to prevent issues near zero

  parsedParallelReactions.reactions.forEach((reaction) => {
    const componentInReaction = reaction.allInvolved.find(comp => comp.name === componentName);
    if (!componentInReaction) return;

    const { reactants, products, rateConstantAtT, rateConstantBackwardAtT, isEquilibriumReaction } = reaction;
    if (rateConstantAtT === undefined) return;

    let forwardRate = rateConstantAtT;
    reactants.forEach(reactant => {
      const conc = Math.max(MIN_CONC, concentrations[reactant.name] || 0);
      forwardRate *= Math.pow(conc, reactant.reactionOrderNum ?? 1);
    });

    let backwardRate = 0;
    if (isEquilibriumReaction && rateConstantBackwardAtT !== undefined && rateConstantBackwardAtT >= 0) {
      backwardRate = rateConstantBackwardAtT;
      for (const product of products) {
        const conc = Math.max(MIN_CONC, concentrations[product.name] || 0);
        backwardRate *= Math.pow(conc, product.reactionOrderNumBackward ?? 1);
      }
    }

    const netReactionRate = forwardRate - backwardRate;
    let stoichCoeff = 0;
    if (componentInReaction.role === 'reactant') stoichCoeff = -componentInReaction.stoichiometricCoefficient;
    else if (componentInReaction.role === 'product') stoichCoeff = componentInReaction.stoichiometricCoefficient;
    
    // The calculated `netReactionRate` is the rate of consumption based on the rate law (e.g., -r_A).
    // For parallel reactions, we must normalize each reaction's rate by its key reactant's stoichiometry
    // to get the true "extent of reaction" before applying it to other components.
    const keyReactant = reaction.reactants[0];
    if (!keyReactant) return; // Should not happen in valid reactions
    
    // Apply the normalized rate to the current component using its stoichiometric coefficient
    totalRate += stoichCoeff * (netReactionRate / Math.abs(keyReactant.stoichiometricCoefficient));
  });
  return totalRate;
};

/**
 * Solves a linear system of equations Ax = b using Gaussian elimination with partial pivoting.
 */
function solveLinearSystem(A: number[][], b: number[]): number[] | null {
    const n = A.length;
    const b_copy = [...b];

    for (let i = 0; i < n; i++) {
        let maxRow = i;
        for (let k = i + 1; k < n; k++) {
            if (Math.abs(A[k][i]) > Math.abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        [A[i], A[maxRow]] = [A[maxRow], A[i]];
        [b_copy[i], b_copy[maxRow]] = [b_copy[maxRow], b_copy[i]];
        
        if (Math.abs(A[i][i]) < 1e-12) return null; // Singular matrix

        for (let k = i + 1; k < n; k++) {
            const factor = A[k][i] / A[i][i];
            for (let j = i; j < n; j++) {
                A[k][j] -= factor * A[i][j];
            }
            b_copy[k] -= factor * b_copy[i];
        }
    }

    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
        let sum = 0;
        for (let j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b_copy[i] - sum) / A[i][i];
    }
    return x;
}

/**
 * Approximates the Jacobian of the ODE system dy/dt = f(t, y) using finite differences.
 */
function calculateNumericalJacobian(
    derivatives: (t: number, y: number[]) => number[],
    t: number,
    y: number[]
): number[][] {
    const n = y.length;
    const J: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
    const fx = derivatives(t, y);
    const h_eps = 1e-8;

    for (let j = 0; j < n; j++) {
        const y_h = [...y];
        const original_yj = y[j];
        let h = h_eps * (Math.abs(original_yj) + 1);
        y_h[j] = original_yj + h;
        h = y_h[j] - original_yj;
        
        const fx_h = derivatives(t, y_h);
        for (let i = 0; i < n; i++) {
            J[i][j] = (fx_h[i] - fx[i]) / h;
        }
    }
    return J;
}


/**
 * Solves a system of stiff ODEs using a variable-step BDF method with a Newton-Raphson solver.
 */
function solveODE_BDF(
  derivatives: (t: number, y: number[]) => number[],
  y0: number[],
  tSpan: [number, number],
  componentNames: string[],
  stoppingCondition?: (t: number, y: number[], dy: number[]) => boolean
): { t: number[], y: number[][] } {
  const [t0, tf] = tSpan;
  let t = t0;
  let y = [...y0];
  const t_out = [t0];
  const y_out = [y0];

  const min_h = 1e-12;
  let h = Math.min(1e-4, tf / 1000); // Smaller initial step size for smoother curves

  const newton_tol = 1e-10; // Tighter tolerance
  const newton_max_iter = 10; // More iterations
  const n = y0.length;
  
  // Maximum number of output points to prevent excessive memory usage
  const max_output_points = 2000;
  let output_interval = (tf - t0) / max_output_points;
  let next_output_time = t0 + output_interval;
  
  while (t < tf && t_out.length < max_output_points) {
    if (t + h > tf) h = tf - t;
    if (h < min_h) break;

    let y_next = [...y];
    let converged = false;

    for (let iter = 0; iter < newton_max_iter; iter++) {
        const f_next = derivatives(t + h, y_next);
        const G = y_next.map((val, i) => val - y[i] - h * f_next[i]);
        
        const normG = Math.sqrt(G.reduce((sum, val) => sum + val*val, 0) / n);
        if (normG < newton_tol) {
            converged = true;
            break;
        }

        const J_f = calculateNumericalJacobian(derivatives, t + h, y_next);
        const J_G = Array.from({length: n}, (_, i) => 
            Array.from({length: n}, (_, j) => (i === j ? 1 : 0) - h * J_f[i][j])
        );

        const neg_G = G.map(val => -val);
        const delta_y = solveLinearSystem(J_G, neg_G);

        if (!delta_y) { break; } 
        
        y_next = y_next.map((val, i) => val + delta_y[i]);
    }
    
    if (converged) {
        y = y_next.map(v => Math.max(0, v));
        t += h;
        
        // Only store output at specified intervals to keep smooth curves
        if (t >= next_output_time || t >= tf) {
          t_out.push(t);
          y_out.push([...y]);
          next_output_time += output_interval;
        }

        if (stoppingCondition) {
            const dy = derivatives(t, y);
            if (stoppingCondition(t, y, dy)) {
                break; 
            }
        }
        
        // Adaptive step size control for smoothness
        h = Math.min(h * 1.5, tf - t, output_interval / 2);
    } else {
        h *= 0.5; // Step failed, reduce step size
    }
  }
  
  // Ensure we have the final point
  if (t_out[t_out.length - 1] < tf) {
    t_out.push(tf);
    y_out.push([...y]);
  }

  // Transpose results for plotting
  const y_transposed: number[][] = Array.from({ length: n }, () => []);
  for (let i = 0; i < t_out.length; i++) {
    for (let j = 0; j < n; j++) {
      y_transposed[j][i] = y_out[i][j];
    }
  }

  return { t: t_out, y: y_transposed };
}

// --- PFR SOLVER WRAPPER ---
export const solvePFR_ODE_System = (
  parsedParallelReactions: ParsedParallelReactions,
  V_total: number,
  initialFlowRates: { [componentName: string]: number },
  volumetricFlowRate: number,
  componentNames: string[]
): PFR_ODE_ProfilePoint[] => {
  if (V_total <= 0) return [{ volume: 0, flowRates: initialFlowRates }];

  const y0 = componentNames.map(name => initialFlowRates[name] || 0);
  const V_span: [number, number] = [0, V_total];

  const derivatives = (V: number, F: number[]): number[] => {
    const currentFlowRates: { [key: string]: number } = {};
    componentNames.forEach((name, i) => { 
        currentFlowRates[name] = Math.max(0, F[i]); 
    });

    const concentrations: { [key: string]: number } = {};
    componentNames.forEach(name => { 
        concentrations[name] = currentFlowRates[name] / volumetricFlowRate; 
    });
    
    return componentNames.map(name => 
        calculateCombinedRate(name, concentrations, parsedParallelReactions)
    );
  };

  const rawSolution = solveODE_BDF(derivatives, y0, V_span, componentNames);

  // Improved interpolation logic to create smoother plots
  const plotPoints = 250; // Reduced from 500 for better performance while maintaining smoothness
  const smoothProfile: PFR_ODE_ProfilePoint[] = [];

  // Add the initial point
  smoothProfile.push({ volume: 0, flowRates: initialFlowRates });

  for (let i = 1; i <= plotPoints; i++) {
    const targetVolume = (i / plotPoints) * V_total;
    
    // Find the appropriate indices for interpolation
    let leftIndex = 0;
    let rightIndex = 1;
    
    for (let j = 0; j < rawSolution.t.length - 1; j++) {
      if (rawSolution.t[j] <= targetVolume && rawSolution.t[j + 1] >= targetVolume) {
        leftIndex = j;
        rightIndex = j + 1;
        break;
      }
    }
    
    // Handle edge cases
    if (rightIndex >= rawSolution.t.length) {
      rightIndex = rawSolution.t.length - 1;
      leftIndex = Math.max(0, rightIndex - 1);
    }

    const V1 = rawSolution.t[leftIndex];
    const V2 = rawSolution.t[rightIndex];
    
    const interpolatedFlowRates: { [key: string]: number } = {};

    if (Math.abs(V2 - V1) < 1e-12 || leftIndex === rightIndex) {
      // No interpolation needed - use the exact value
      componentNames.forEach((name, compIndex) => {
        interpolatedFlowRates[name] = Math.max(0, rawSolution.y[compIndex][leftIndex]);
      });
    } else {
      // Linear interpolation
      const fraction = (targetVolume - V1) / (V2 - V1);
      componentNames.forEach((name, compIndex) => {
        const F1 = rawSolution.y[compIndex][leftIndex];
        const F2 = rawSolution.y[compIndex][rightIndex];
        const F_interp = F1 + fraction * (F2 - F1);
        interpolatedFlowRates[name] = Math.max(0, F_interp);
      });
    }
    
    smoothProfile.push({ volume: targetVolume, flowRates: interpolatedFlowRates });
  }

  return smoothProfile;
};


// --- CSTR SOLVER (Unchanged from last step) ---
export const solveCSTRParallel = (
    parsedParallelReactions: ParsedParallelReactions,
    V: number,
    initialFlowRates: { [componentName: string]: number },
    volumetricFlowRate: number
): { [componentName: string]: number } => {
    if (V <= 0 || volumetricFlowRate <= 0) {
        return { ...initialFlowRates };
    }

    const componentNames = Object.keys(initialFlowRates);
    const tau = V / volumetricFlowRate; 

    const inletConcentrations: { [key: string]: number } = {};
    componentNames.forEach(name => {
        inletConcentrations[name] = initialFlowRates[name] / volumetricFlowRate;
    });

    const derivatives = (t: number, C_vector: number[]): number[] => {
        const currentConcentrations: { [key: string]: number } = {};
        componentNames.forEach((name, i) => {
            currentConcentrations[name] = Math.max(0, C_vector[i]);
        });

        const dCdt = componentNames.map((name, i) => {
            const rate_i = calculateCombinedRate(name, currentConcentrations, parsedParallelReactions);
            const C_in_i = inletConcentrations[name] || 0;
            const C_i = currentConcentrations[name];
            return (C_in_i - C_i) / tau + rate_i;
        });
        return dCdt;
    };

    const C0_vector = componentNames.map(name => inletConcentrations[name]);

    const integrationTime = Math.max(tau * 200, 10);
    const t_span: [number, number] = [0, integrationTime];

    const steadyStateTolerance = 1e-7;
    const stoppingCondition = (t: number, y: number[], dy: number[]): boolean => {
        if (t < Math.min(tau, 2.0)) {
            return false;
        }
        const norm_dy_sq = dy.reduce((sum, val) => sum + val * val, 0);
        const rms = Math.sqrt(norm_dy_sq / dy.length);
        return rms < steadyStateTolerance;
    };

    const solution = solveODE_BDF(
        derivatives, 
        C0_vector, 
        t_span, 
        componentNames, 
        stoppingCondition
    );

    const steadyStateConcentrations: { [key: string]: number } = {};
    componentNames.forEach((name, i) => {
        const concentrationProfile = solution.y[i];
        const finalConcentration = concentrationProfile?.[concentrationProfile.length - 1] ?? 0;
        steadyStateConcentrations[name] = finalConcentration;
    });

    const outletFlowRates: { [key: string]: number } = {};
    componentNames.forEach(name => {
        outletFlowRates[name] = Math.max(0, steadyStateConcentrations[name] * volumetricFlowRate);
    });

    return outletFlowRates;
};


// --- LIMITING REACTANT (Unchanged) ---
export const findLimitingReactant = (
  parsedParallelReactions: ParsedParallelReactions,
  components: Component[]
): { name: string; initialFlow: number } | null => {
  if (parsedParallelReactions.reactions.length === 0) return null;
  const primaryReaction = parsedParallelReactions.reactions[0];
  if (primaryReaction.reactants.length === 0) return null;

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
  
  if (!limitingReactant && primaryReaction.reactants.length > 0) {
    const firstReactant = primaryReaction.reactants[0];
    const component = components.find(c => c.name === firstReactant.name);
    return { name: firstReactant.name, initialFlow: parseFloat(component?.initialFlowRate || '0') };
  }
  return limitingReactant;
};