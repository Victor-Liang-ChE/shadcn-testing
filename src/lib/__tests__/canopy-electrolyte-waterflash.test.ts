import { describe, expect, it } from 'vitest';

import { PRESET_COMPOUNDS, useCanopyStore } from '../canopy/store';
import { solveElectrolyteEquilibrium, solveFlashDrum, solveHeater, solveValve } from '../canopy/solver';
import {
  bubblePointTWaterSaturated,
  dewPointTWaterSaturated,
  flashPT,
  flashPTWaterSaturated,
  waterSaturatedComposition,
} from '../canopy/thermo';
import type { CanopyCompound, MaterialStream } from '../canopy/types';

const benzene = PRESET_COMPOUNDS.find(c => c.name === 'BENZENE')!;

const water: CanopyCompound = {
  ...benzene,
  name: 'WATER',
  displayName: 'Water',
  molecularWeight: 18.015,
  Tc_K: 647.096,
  Pc_Pa: 22064000,
  omega: 0.3443,
  Tb_K: 373.15,
  specificGravity: 1,
  isElectrolyteSolvent: true,
  dielectricConst: { e1: 78.36, e2: 0, e3_K: 298.15 },
  antoine: {
    C1: 73.649,
    C2: -7258.2,
    C3: 0,
    C4: 0,
    C5: -7.3037,
    C6: 4.1653e-6,
    C7: 2,
    Tmin_K: 273.16,
    Tmax_K: 647.096,
  },
  waterSolubility: undefined,
  hcSolubility: undefined,
};

const sodium: CanopyCompound = {
  ...benzene,
  name: 'NA+',
  displayName: 'Na+',
  chargeNumber: 1,
};

const chloride: CanopyCompound = {
  ...benzene,
  name: 'CL-',
  displayName: 'Cl-',
  chargeNumber: -1,
};

const apparentSalt: CanopyCompound = {
  ...benzene,
  name: 'NACL',
  displayName: 'NaCl',
  chargeNumber: 0,
};

function makeStream(
  compounds: CanopyCompound[],
  totalFlow_molps: number,
  moleFractions: number[],
  T_K = 298.15,
  P_Pa = 101325,
): MaterialStream {
  const flash = flashPT(compounds, moleFractions, T_K, P_Pa);
  return {
    id: 'S1',
    name: 'S1',
    T_K,
    P_Pa,
    totalFlow_molps,
    moleFractions,
    phase: flash.phase,
    vaporFraction: flash.vaporFraction,
    H_Jpmol: flash.H_Jpmol,
    x_liquid: flash.x,
    y_vapor: flash.y,
    sourceUnitId: null,
    sourcePortId: null,
    targetUnitId: null,
    targetPortId: null,
    solved: true,
  };
}

describe('electrolyte equilibrium and water-saturated flash behavior', () => {
  it('builds a water-saturated composition and enriched flash result', () => {
    const compounds = [water, benzene];
    const saturated = waterSaturatedComposition(compounds, [0, 1], 298.15);
    const flash = flashPTWaterSaturated(compounds, [0, 1], 298.15, 101325);

    expect(saturated).not.toBeNull();
    expect((saturated?.requiredWaterFraction ?? 0)).toBeGreaterThan(0);
    expect((saturated?.zSaturated[0] ?? 0)).toBeGreaterThan(0);
    expect(flash.waterAddedFraction).toBeGreaterThan(0);
  });

  it('computes finite water-saturated bubble and dew points', () => {
    const compounds = [water, benzene];
    const z = [0.01, 0.99];
    const Tbub = bubblePointTWaterSaturated(compounds, z, 101325);
    const Tdew = dewPointTWaterSaturated(compounds, z, 101325);

    expect(Tbub).not.toBeNull();
    expect(Tdew).not.toBeNull();
    expect((Tdew ?? 0)).toBeGreaterThan(Tbub ?? 0);
  });

  it('flash drum can enforce water-saturated feed handling', () => {
    const compounds = [water, benzene];
    const feed = makeStream(compounds, 100, [0.005, 0.995], 330, 101325);
    const result = solveFlashDrum(compounds, feed, {
      T_K: 330,
      P_Pa: 101325,
      waterSaturatedFeed: true,
    });

    expect(result.waterAddedForSaturation_molps).toBeGreaterThan(0);
    expect(result.vapor.totalFlow_molps! + result.liquid.totalFlow_molps!).toBeGreaterThan(feed.totalFlow_molps);
  });

  it('propagates water-saturated handling into heater and valve units', () => {
    const compounds = [water, benzene];
    const feed = makeStream(compounds, 25, [0.005, 0.995], 320, 250000);
    const heated = solveHeater(compounds, feed, {
      targetT_K: 360,
      outletP_Pa: 220000,
      waterSaturatedFeed: true,
    });
    const letdown = solveValve(compounds, materializeStreamLike(feed, heated.outlet), {
      outletP_Pa: 101325,
      waterSaturatedFeed: true,
    });

    expect(heated.waterAddedForSaturation_molps).toBeGreaterThan(0);
    expect(heated.outlet.totalFlow_molps!).toBeGreaterThan(feed.totalFlow_molps);
    expect(letdown.waterAddedForSaturation_molps).toBeGreaterThanOrEqual(0);
    expect(letdown.outlet.totalFlow_molps!).toBeGreaterThanOrEqual(heated.outlet.totalFlow_molps!);
  });

  it('electrolyte equilibrium dissociates an apparent salt into ions', () => {
    const compounds = [water, apparentSalt, sodium, chloride];
    const inlet = makeStream(compounds, 10, [0.8, 0.2, 0, 0], 298.15, 101325);

    const result = solveElectrolyteEquilibrium(compounds, inlet, {
      apparentSpeciesIndex: 1,
      cationIndex: 2,
      anionIndex: 3,
      cationStoich: 1,
      anionStoich: 1,
      log10K: 4,
      maxDissociationFraction: 0.99,
    }, {
      fluidPackage: 'Electrolyte-NRTL',
      interactionParams: {
        elecNrtl: {
          nrtl: {
            A: Array.from({ length: 4 }, () => Array(4).fill(0)),
            B: Array.from({ length: 4 }, () => Array(4).fill(0)),
            alpha: Array.from({ length: 4 }, () => Array(4).fill(0.2)),
          },
          charges: [0, 0, 1, -1],
          epsilon_r: 78.36,
          solventMW: 18.015,
        },
      },
    });

    expect(result.dissociationFraction).toBeGreaterThan(0.5);
    expect(result.ionicStrength).toBeGreaterThan(0);
    expect(result.outlet.moleFractions?.[2] ?? 0).toBeGreaterThan(0);
    expect(result.outlet.moleFractions?.[3] ?? 0).toBeGreaterThan(0);
  });

  it('supports multi-reaction electrolyte sets with sequential speciation passes', () => {
    const potassium: CanopyCompound = { ...benzene, name: 'K+', displayName: 'K+', chargeNumber: 1 };
    const bisulfate: CanopyCompound = { ...benzene, name: 'HSO4-', displayName: 'HSO4-', chargeNumber: -1 };
    const sulfate: CanopyCompound = { ...benzene, name: 'SO4--', displayName: 'SO4--', chargeNumber: -2 };
    const acidSalt: CanopyCompound = { ...benzene, name: 'KHSO4', displayName: 'KHSO4', chargeNumber: 0 };
    const compounds = [water, acidSalt, potassium, bisulfate, sodium, chloride, apparentSalt, sulfate];
    const inlet = makeStream(compounds, 12, [0.78, 0.08, 0, 0, 0, 0, 0.14, 0], 298.15, 101325);

    const result = solveElectrolyteEquilibrium(compounds, inlet, {
      apparentSpeciesIndex: 6,
      cationIndex: 4,
      anionIndex: 5,
      reactionSet: [
        { name: 'NaCl', apparentSpeciesIndex: 6, cationIndex: 4, anionIndex: 5, log10K: 4.5, active: true },
        { name: 'KHSO4', apparentSpeciesIndex: 1, cationIndex: 2, anionIndex: 3, log10K: 3.5, active: true },
      ],
      maxPasses: 6,
    }, {
      fluidPackage: 'Electrolyte-NRTL',
      interactionParams: {
        elecNrtl: {
          nrtl: {
            A: Array.from({ length: compounds.length }, () => Array(compounds.length).fill(0)),
            B: Array.from({ length: compounds.length }, () => Array(compounds.length).fill(0)),
            alpha: Array.from({ length: compounds.length }, () => Array(compounds.length).fill(0.2)),
          },
          charges: compounds.map(compound => compound.chargeNumber ?? 0),
          epsilon_r: 78.36,
          solventMW: 18.015,
        },
      },
    });

    expect(result.activeReactionCount).toBe(2);
    expect(result.equilibriumPasses).toBeGreaterThanOrEqual(1);
    expect(result.dissociationFraction).toBeGreaterThan(0.2);
    expect(result.outlet.moleFractions?.[2] ?? 0).toBeGreaterThan(0);
    expect(result.outlet.moleFractions?.[4] ?? 0).toBeGreaterThan(0);
  });

  it('supports generic stoichiometric association reactions with temperature-dependent equilibrium constants', () => {
    const monomer: CanopyCompound = { ...benzene, name: 'HA', displayName: 'HA' };
    const dimer: CanopyCompound = { ...benzene, name: 'H2A2', displayName: '(HA)2', molecularWeight: monomer.molecularWeight * 2 };
    const compounds = [water, monomer, dimer];
    const coldInlet = makeStream(compounds, 10, [0.75, 0.25, 0], 285.15, 101325);
    const hotInlet = makeStream(compounds, 10, [0.75, 0.25, 0], 335.15, 101325);
    const reaction = {
      name: 'HA dimerization',
      stoichiometry: [0, -2, 1],
      log10K: 1.1,
      referenceT_K: 298.15,
      deltaH_Jmol: -12000,
      active: true,
    };

    const coldResult = solveElectrolyteEquilibrium(compounds, coldInlet, {
      apparentSpeciesIndex: 1,
      cationIndex: 1,
      anionIndex: 2,
      reactionSet: [reaction],
      maxPasses: 4,
    });
    const hotResult = solveElectrolyteEquilibrium(compounds, hotInlet, {
      apparentSpeciesIndex: 1,
      cationIndex: 1,
      anionIndex: 2,
      reactionSet: [reaction],
      maxPasses: 4,
    });

    expect(coldResult.reactionExtents_molps).toHaveLength(1);
    expect(coldResult.reactionExtents_molps[0]).toBeGreaterThan(0);
    expect(coldResult.outlet.moleFractions?.[2] ?? 0).toBeGreaterThan(0);
    expect(coldResult.equilibriumConstant).toBeGreaterThan(hotResult.equilibriumConstant);
    expect(coldResult.reactionExtents_molps[0]).toBeGreaterThan(hotResult.reactionExtents_molps[0]);
  });

  it('store can add the electrolyte-equilibrium unit family', () => {
    useCanopyStore.setState(state => ({ ...state, units: {}, streams: {}, nodePositions: {} }));
    useCanopyStore.getState().addUnitFromPalette('ElectrolyteEquilibrium', { x: 160, y: 120 });

    const units = Object.values(useCanopyStore.getState().units);
    expect(units).toHaveLength(1);
    expect(units[0].type).toBe('ElectrolyteEquilibrium');
    expect(units[0].ports.map(port => port.label)).toEqual(['Feed', 'Equilibrated Outlet']);
  });
});

function materializeStreamLike(fallback: MaterialStream, partial: Partial<MaterialStream>): MaterialStream {
  return {
    ...fallback,
    ...partial,
    totalFlow_molps: partial.totalFlow_molps ?? fallback.totalFlow_molps,
    moleFractions: partial.moleFractions ?? fallback.moleFractions,
    x_liquid: partial.x_liquid ?? fallback.x_liquid,
    y_vapor: partial.y_vapor ?? fallback.y_vapor,
    phase: partial.phase ?? fallback.phase,
    vaporFraction: partial.vaporFraction ?? fallback.vaporFraction,
    H_Jpmol: partial.H_Jpmol ?? fallback.H_Jpmol,
    T_K: partial.T_K ?? fallback.T_K,
    P_Pa: partial.P_Pa ?? fallback.P_Pa,
    solved: partial.solved ?? true,
  };
}
