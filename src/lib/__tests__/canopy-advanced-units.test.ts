import { describe, expect, it } from 'vitest';
import type { MaterialStream } from '../canopy/types';
import { PRESET_COMPOUNDS } from '../canopy/store';
import {
  solveCentrifuge,
  solveCyclone,
  solveFilter,
  solveMembrane,
  solvePIDController,
  solvePIDControllerSignal,
  solveColumnRigorous,
  solveSteamHeater,
  solveSteamTurbine,
} from '../canopy/solver';
import { flashPT } from '../canopy/thermo';

const benzene = PRESET_COMPOUNDS.find(c => c.name === 'BENZENE')!;
const cumene = PRESET_COMPOUNDS.find(c => c.name === 'CUMENE')!;
const dipb = PRESET_COMPOUNDS.find(c => c.name === 'P-DIISOPROPYLBENZENE')!;
const compounds = [benzene, cumene, dipb];

function makeStream(name: string, totalFlow_molps: number, moleFractions: number[], T_K = 360, P_Pa = 101325): MaterialStream {
  const flash = flashPT(compounds, moleFractions, T_K, P_Pa);
  return {
    id: name,
    name,
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

describe('advanced canopy unit models', () => {
  it('membrane conserves mass and enriches the more permeable component', () => {
    const feed = makeStream('mem-feed', 100, [0.6, 0.3, 0.1], 370);
    const result = solveMembrane(compounds, feed, {
      stageCut: 0.25,
      selectivities: [3.0, 0.7, 0.3],
      permeateP_Pa: 60000,
      retentateP_Pa: 101325,
    });

    expect(result.permeate.totalFlow_molps + result.retentate.totalFlow_molps).toBeCloseTo(100, 8);
    expect(result.permeate.moleFractions[0]).toBeGreaterThan(feed.moleFractions[0]);
    expect(result.retentate.moleFractions[2]).toBeGreaterThan(result.permeate.moleFractions[2]);
  });

  it('membrane solution-diffusion and rejection modes remain mass conservative', () => {
    const feed = makeStream('mem-rigorous', 90, [0.55, 0.30, 0.15], 365, 202650);

    const diffusion = solveMembrane(compounds, feed, {
      model: 'solution-diffusion',
      stageCut: 0.2,
      membraneArea_m2: 80,
      permeanceGPU: [500, 160, 40],
      foulingFactor: 0.85,
      permeateP_Pa: 70000,
      retentateP_Pa: 180000,
    });
    const rejection = solveMembrane(compounds, feed, {
      model: 'rejection',
      stageCut: 0.2,
      rejectionCoefficients: [0.15, 0.55, 0.85],
      permeateP_Pa: 70000,
      retentateP_Pa: 180000,
    });

    expect(diffusion.permeate.totalFlow_molps + diffusion.retentate.totalFlow_molps).toBeCloseTo(feed.totalFlow_molps, 8);
    expect(rejection.permeate.totalFlow_molps + rejection.retentate.totalFlow_molps).toBeCloseTo(feed.totalFlow_molps, 8);
    expect(diffusion.permeate.moleFractions[0]).toBeGreaterThan(diffusion.retentate.moleFractions[0]);
    expect(rejection.retentate.moleFractions[2]).toBeGreaterThan(rejection.permeate.moleFractions[2]);
  });

  it('cyclone uses a cut curve to send larger particles to underflow', () => {
    const feed = makeStream('cyclone-feed', 80, [0.45, 0.35, 0.20]);
    const result = solveCyclone(compounds, feed, {
      cutSize_um: 8,
      particleSizes_um: [18, 6, 2],
      pressureDrop_Pa: 2000,
    });

    expect(result.overflow.totalFlow_molps + result.underflow.totalFlow_molps).toBeCloseTo(80, 8);
    expect(result.underflow.moleFractions[0]).toBeGreaterThan(result.overflow.moleFractions[0]);
    expect(result.overflow.P_Pa).toBeCloseTo(feed.P_Pa - 2000, 8);
  });

  it('filter retains specified components into the cake stream', () => {
    const feed = makeStream('filter-feed', 50, [0.2, 0.5, 0.3]);
    const result = solveFilter(compounds, feed, {
      retainedFractions: [0.95, 0.10, 0.80],
      pressureDrop_Pa: 4000,
    });

    expect(result.filtrate.totalFlow_molps + result.cake.totalFlow_molps).toBeCloseTo(50, 8);
    expect(result.cake.moleFractions[0]).toBeGreaterThan(result.filtrate.moleFractions[0]);
    expect(result.cake.moleFractions[2]).toBeGreaterThan(result.filtrate.moleFractions[2]);
  });

  it('centrifuge captures coarse particles while carrying some liquid with the cake', () => {
    const feed = makeStream('centrifuge-feed', 60, [0.25, 0.50, 0.25], 330, 150000);
    const result = solveCentrifuge(compounds, feed, {
      cutSize_um: 4,
      sharpness: 8,
      particleSizes_um: [12, 1, 0],
      liquidSplitToCake: 0.12,
      pressureDrop_Pa: 8000,
    });

    expect(result.centrate.totalFlow_molps + result.cake.totalFlow_molps).toBeCloseTo(60, 8);
    expect(result.cake.moleFractions[0]).toBeGreaterThan(result.centrate.moleFractions[0]);
    expect(result.cake.totalFlow_molps).toBeGreaterThan(0);
    expect(result.centrate.P_Pa).toBeCloseTo(feed.P_Pa - 8000, 8);
  });

  it('steam heater computes duty and implied steam consumption', () => {
    const feed = makeStream('steam-feed', 40, [0.5, 0.3, 0.2], 320);
    const result = solveSteamHeater(compounds, feed, {
      targetT_K: 440,
      steamLevel: 'Auto',
    });

    expect(result.outlet.T_K!).toBeCloseTo(440, 1);
    expect(result.duty_W).toBeGreaterThan(0);
    expect(result.utilityType).toBe('mpSteam');
    expect(result.steamFlow_kgps).toBeGreaterThan(0);
  });

  it('steady-state PID controller forwards the stream and computes bounded output', () => {
    const feed = makeStream('ctrl-feed', 25, [0.4, 0.4, 0.2], 330);
    const result = solvePIDController(feed, {
      measuredProperty: 'T_K',
      setpoint: 350,
      gain: 2,
      bias: 10,
      lowLimit: 0,
      highLimit: 100,
    });

    expect(result.measuredValue).toBe(330);
    expect(result.error).toBe(20);
    expect(result.controllerOutput).toBe(50);
    expect(result.outlet.totalFlow_molps).toBe(feed.totalFlow_molps);
    expect(result.outlet.T_K).toBe(feed.T_K);
  });

  it('dynamic PID signal applies measurement lag and output slew limits', () => {
    const result = solvePIDControllerSignal(360, {
      measuredProperty: 'T_K',
      setpoint: 420,
      gain: 1.5,
      timestep_s: 1,
      measurementLag_s: 2,
      outputRateLimit_per_s: 12,
      lowLimit: 0,
      highLimit: 500,
      filteredMeasuredValue: 300,
      previousMeasuredValue: 300,
      previousError: 90,
      previousPreviousError: 80,
      lastOutput: 200,
    });

    expect(result.filteredMeasuredValue).toBeCloseTo(320, 8);
    expect(result.controllerOutput).toBeCloseTo(212, 8);
    expect(result.error).toBeCloseTo(100, 8);
  });

  it('steam turbine recovers power as pressure is let down', () => {
    const feed = makeStream('steam-turbine-feed', 45, [1, 0, 0], 560, 900000);
    const result = solveSteamTurbine(compounds, feed, {
      outletP_Pa: 250000,
      efficiency: 0.78,
      generatorEfficiency: 0.95,
    });

    expect(result.outlet.P_Pa).toBe(250000);
    expect(result.outlet.T_K!).toBeLessThan(feed.T_K);
    expect(result.shaftWork_W).toBeGreaterThan(0);
    expect(result.electricPower_W).toBeGreaterThan(0);
  });

  it('rigorous column side draws preserve overall material balance and expose stage duties', () => {
    const feed = makeStream('radfrac-feed', 100, [0.5, 0.35, 0.15], 390);
    const result = solveColumnRigorous(compounds, feed, {
      nStages: 18,
      feedStage: 9,
      refluxRatio: 2.4,
      distillateRate_molps: 40,
      liquidSideDrawStage: 11,
      liquidSideDrawFraction: 0.08,
      vaporSideDrawStage: 5,
      vaporSideDrawFraction: 0.03,
      condenserP_Pa: 101325,
      reboilerP_Pa: 140000,
      condenserType: 'Total',
    });

    const totalProducts =
      result.distillate.totalFlow_molps!
      + result.bottoms.totalFlow_molps!
      + (result.liquidSideDraw?.totalFlow_molps ?? 0)
      + (result.vaporSideDraw?.totalFlow_molps ?? 0);

    expect(totalProducts).toBeCloseTo(feed.totalFlow_molps, 6);
    expect(result.liquidSideDraw?.totalFlow_molps ?? 0).toBeGreaterThan(0);
    expect(result.vaporSideDraw?.totalFlow_molps ?? 0).toBeGreaterThan(0);
    expect(result.stageDuties_W).toHaveLength(18);
    expect(result.stageDuties_W.every(value => Number.isFinite(value))).toBe(true);
  });
});
