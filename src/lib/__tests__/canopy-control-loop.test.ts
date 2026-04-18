import { beforeEach, describe, expect, it } from 'vitest';
import type { MaterialStream } from '../canopy/types';
import { PRESET_COMPOUNDS, useCanopyStore } from '../canopy/store';
import { flashPT } from '../canopy/thermo';

function resetState() {
  useCanopyStore.setState(state => ({
    ...state,
    compounds: PRESET_COMPOUNDS.slice(0, 3),
    fluidPackage: 'Ideal',
    interactionParams: {},
    units: {},
    streams: {},
    nodePositions: {},
    designSpecs: [],
    solverErrors: [],
    solved: false,
  }));
}

function makeFeedStream(): MaterialStream {
  const compounds = PRESET_COMPOUNDS.slice(0, 3);
  const moleFractions = [0.7, 0.2, 0.1];
  const T_K = 315;
  const P_Pa = 101325;
  const flash = flashPT(compounds, moleFractions, T_K, P_Pa);
  return {
    id: 'S1',
    name: 'Feed',
    T_K,
    P_Pa,
    totalFlow_molps: 50,
    moleFractions,
    phase: flash.phase,
    vaporFraction: flash.vaporFraction,
    H_Jpmol: flash.H_Jpmol,
    x_liquid: flash.x,
    y_vapor: flash.y,
    sourceUnitId: null,
    sourcePortId: null,
    targetUnitId: 'U1',
    targetPortId: 'U1-in',
    solved: true,
  };
}

describe('canopy control loop', () => {
  beforeEach(() => {
    resetState();
  });

  it('PID controllers manipulate unit parameters from measured stream values', () => {
    const store = useCanopyStore.getState();
    store.addUnitFromPalette('Heater', { x: 100, y: 100 });
    store.addUnitFromPalette('PIDController', { x: 260, y: 100 });

    const { units } = useCanopyStore.getState();
    const heater = units.U1;
    const controller = units.U2;
    const productStream: MaterialStream = {
      ...makeFeedStream(),
      id: 'S2',
      name: 'Product',
      solved: false,
      sourceUnitId: 'U1',
      sourcePortId: 'U1-out',
      targetUnitId: null,
      targetPortId: null,
    };

    useCanopyStore.setState(state => ({
      ...state,
      streams: {
        S1: makeFeedStream(),
        S2: productStream,
      },
      units: {
        U1: {
          ...heater,
          params: { ...heater.params, targetT_K: 330, outletP_Pa: 101325 },
          ports: heater.ports.map(port => ({
            ...port,
            streamId: port.label === 'Inlet' ? 'S1' : 'S2',
          })),
        },
        U2: {
          ...controller,
          params: {
            ...controller.params,
            setpoint: 385,
            gain: 0.8,
            controllerTimestep_s: 1,
            measurementLag_s: 1,
            outputRateLimit_per_s: 20,
            lowLimit: 320,
            highLimit: 420,
            measuredObjectType: 'stream',
            measuredObjectId: 'S2',
            measuredProperty: 'T_K',
            manipulatedUnitId: 'U1',
            manipulatedParamName: 'targetT_K',
            bias: 330,
          },
        },
      },
    }));

    useCanopyStore.getState().solveAll();
    const solved = useCanopyStore.getState();

    expect(solved.units.U1.params.targetT_K as number).toBeGreaterThan(330);
    expect(solved.streams.S2.T_K).toBeGreaterThan(330);
    expect((solved.units.U2.params.controllerOutput as number) ?? 0).toBeGreaterThan(330);
    expect((solved.units.U2.params.filteredMeasuredValue as number) ?? 0).toBeGreaterThan(315);
    expect(solved.solverErrors).toEqual([]);
  });
});
