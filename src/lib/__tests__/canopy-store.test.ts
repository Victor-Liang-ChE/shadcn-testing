import { beforeEach, describe, expect, it } from 'vitest';
import { useCanopyStore } from '../canopy/store';

function resetFlowsheetState() {
  useCanopyStore.setState(state => ({
    ...state,
    units: {},
    streams: {},
    nodePositions: {},
    openTabs: [],
    activeTab: null,
    solved: false,
    solverErrors: [],
  }));
}

describe('canopy store unit palette integration', () => {
  beforeEach(() => {
    resetFlowsheetState();
  });

  it('adds RadFrac as a first-class rigorous column', () => {
    useCanopyStore.getState().addUnitFromPalette('RadFrac', { x: 120, y: 80 });
    const units = Object.values(useCanopyStore.getState().units);
    expect(units).toHaveLength(1);
    expect(units[0].type).toBe('RadFrac');
    expect(units[0].params.rigorousMode).toBe(true);
    expect(units[0].ports.map(port => port.label)).toEqual(['Feed', 'Distillate', 'Side Vapor', 'Side Liquid', 'Bottoms']);
  });

  it('adds Extractor with feed, solvent, extract, and raffinate ports', () => {
    useCanopyStore.getState().addUnitFromPalette('Extractor', { x: 160, y: 120 });
    const units = Object.values(useCanopyStore.getState().units);
    expect(units).toHaveLength(1);
    expect(units[0].type).toBe('Extractor');
    expect(units[0].params.nStages).toBe(6);
    expect(units[0].ports.map(port => port.label)).toEqual([
      'Feed Inlet',
      'Solvent Inlet',
      'Extract',
      'Raffinate',
    ]);
  });

  it('adds membrane, steam-network, and controller unit families to the palette', () => {
    useCanopyStore.getState().addUnitFromPalette('Membrane', { x: 200, y: 120 });
    useCanopyStore.getState().addUnitFromPalette('Centrifuge', { x: 230, y: 120 });
    useCanopyStore.getState().addUnitFromPalette('SteamTurbine', { x: 245, y: 130 });
    useCanopyStore.getState().addUnitFromPalette('SteamTrap', { x: 260, y: 140 });
    useCanopyStore.getState().addUnitFromPalette('PIDController', { x: 290, y: 160 });
    const units = Object.values(useCanopyStore.getState().units);
    expect(units.map(unit => unit.type)).toEqual(['Membrane', 'Centrifuge', 'SteamTurbine', 'SteamTrap', 'PIDController']);
    expect(units[0].ports.map(port => port.label)).toEqual(['Feed', 'Permeate', 'Retentate']);
    expect(units[1].ports.map(port => port.label)).toEqual(['Feed', 'Centrate', 'Cake']);
    expect(units[2].ports.map(port => port.label)).toEqual(['Inlet', 'Outlet']);
    expect(units[3].ports.map(port => port.label)).toEqual(['Inlet', 'Flash Vapor', 'Condensate']);
    expect(units[4].ports.map(port => port.label)).toEqual(['Inlet', 'Outlet']);
  });
});
