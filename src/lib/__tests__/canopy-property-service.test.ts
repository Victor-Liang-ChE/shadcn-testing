import { describe, expect, it } from 'vitest';

import {
  buildCanopyCompoundFromRemoteRecords,
  loadRemoteInteractionParams,
  loadRemotePetroAssay,
  searchRemotePetroAssays,
} from '../canopy/property-service';
import type { CanopyCompound } from '../canopy/types';

function createMockClient(responseMap: Record<string, any[]>) {
  return {
    from(table: string) {
      return {
        select() {
          return this;
        },
        ilike() {
          return this;
        },
        in() {
          return this;
        },
        order() {
          return this;
        },
        limit() {
          return Promise.resolve({ data: responseMap[table] ?? [] });
        },
      };
    },
  } as any;
}

describe('canopy property-service', () => {
  it('builds a Canopy compound from remote DB records', () => {
    const compound = buildCanopyCompoundFromRemoteRecords(
      'WATER',
      {
        MW: 18.015,
        TC: 647.1,
        PC: 22064000,
        OMEGA: 0.344,
        GMUQR: 0.92,
        GMUQQ: 1.4,
        CPIGDP_A: 33363,
        CPIGDP_B: 26790,
        CPIGDP_C: 2610.5,
        CPIGDP_D: 8896,
        CPIGDP_E: 1169,
      },
      {
        C1: 73.649,
        C2: -7258.2,
        C3: 0,
        C4: 0,
        C5: -7.3037,
        C6: 4.1653e-6,
        C7: 2,
        Tmin: 273.16,
        Tmax: 647.1,
      },
      {
        MolecularWeight: { value: '18.015' },
        HeatOfFormation: { value: '-241820000' },
        GibbsEnergyOfFormation: { value: '-228570000' },
        DielectricConstant: { value: '78.36' },
        SpecificGravity: { value: '1.0' },
        FlashPoint: { value: '400' },
        AIT: { value: '860' },
        RefractiveIndex: { value: '1.333' },
        ANILPT: { value: '350' },
        CETANE: { value: '5' },
        HigherHeatingValue: { value: '-285830' },
        NetHeatingValue: { value: '-241820000' },
        UnifacVLE: { group: [{ id: '16', count: '1' }] },
        LiquidViscosity: {
          A: { value: '-52.843' },
          B: { value: '3703.6' },
          C: { value: '5.866' },
          D: { value: '-5.879e-29' },
        },
      },
    );

    expect(compound.displayName).toBe('Water');
    expect(compound.molecularWeight).toBeCloseTo(18.015, 6);
    expect(compound.uniquac_r).toBeCloseTo(0.92, 6);
    expect(compound.dielectricConst?.e1).toBeCloseTo(78.36, 6);
    expect(compound.isElectrolyteSolvent).toBe(true);
    expect(compound.unifacGroups).toEqual({ 16: 1 });
    expect(compound.Hf_298_Jmol).toBeCloseTo(-241820, 6);
    expect(compound.specificGravity).toBeCloseTo(1, 6);
    expect(compound.flashPoint_K).toBeCloseTo(400, 6);
    expect(compound.autoignitionTemp_K).toBeCloseTo(860, 6);
    expect(compound.refractiveIndex).toBeCloseTo(1.333, 6);
    expect(compound.anilinePoint_K).toBeCloseTo(350, 6);
    expect(compound.cetaneNumber).toBeCloseTo(5, 6);
    expect(compound.Hcomb_Jmol).toBeCloseTo(-285830, 6);
    expect(compound.hCombustion_Jpkmol).toBeCloseTo(-241820000, 6);
  });

  it('derives Electrolyte-NRTL parameters from remote NRTL rows and solvent metadata', async () => {
    const compounds: CanopyCompound[] = [
      {
        name: 'WATER',
        displayName: 'Water',
        molecularWeight: 18.015,
        Tc_K: 647.1,
        Pc_Pa: 22064000,
        omega: 0.344,
        dipprCoeffs: {},
        antoine: null,
        chargeNumber: 0,
        isElectrolyteSolvent: true,
        dielectricConst: { e1: 78.36, e2: 0, e3_K: 298.15 },
      },
      {
        name: 'NA+',
        displayName: 'Na+',
        molecularWeight: 22.99,
        Tc_K: 0,
        Pc_Pa: 0,
        omega: 0,
        dipprCoeffs: {},
        antoine: null,
        chargeNumber: 1,
      },
    ];

    const client = createMockClient({
      apv140_binary_params: [
        { Compound1Name: 'WATER', Compound2Name: 'NA+', ElementLabel: 'aij', Value: 1.2 },
        { Compound1Name: 'WATER', Compound2Name: 'NA+', ElementLabel: 'bij', Value: 250.0 },
        { Compound1Name: 'WATER', Compound2Name: 'NA+', ElementLabel: 'cij', Value: 0.2 },
      ],
    });

    const result = await loadRemoteInteractionParams('Electrolyte-NRTL', compounds, client);

    expect(result.params?.elecNrtl).toBeDefined();
    expect(result.params?.elecNrtl?.charges).toEqual([0, 1]);
    expect(result.params?.elecNrtl?.epsilon_r).toBeCloseTo(78.36, 6);
    expect(result.params?.elecNrtl?.nrtl.a_ext?.[0][1]).toBeCloseTo(1.2, 6);
  });

  it('searches and loads remote petroleum assays into pseudo-components', async () => {
    const client = {
      from(table: string) {
        return {
          select() { return this; },
          ilike() { return this; },
          eq() { return this; },
          order() {
            if (table === 'petro_assay_cuts') {
              return Promise.resolve({
                data: [
                  { assay_id: 'A1', cut_order: 1, cut_name: 'Naphtha', tbp_start_k: 320, tbp_end_k: 420, liquid_volume_fraction: 0.3, specific_gravity: 0.72 },
                  { assay_id: 'A1', cut_order: 2, cut_name: 'Diesel', tbp_start_k: 520, tbp_end_k: 650, liquid_volume_fraction: 0.7, specific_gravity: 0.85 },
                ],
              });
            }
            if (table === 'petro_assay_light_ends') {
              return Promise.resolve({
                data: [
                  { assay_id: 'A1', component_name: 'WATER', fraction_value: 0.02, fraction_basis: 'mole', normal_boiling_point_k: 373.15 },
                ],
              });
            }
            return Promise.resolve({ data: [] });
          },
          limit() {
            return Promise.resolve({
              data: table === 'petro_assays'
                ? [{ assay_id: 'A1', assay_name: 'Example Crude', source: 'Aspen Extract', location: 'Field 1', cut_count: 2 }]
                : table === 'apv140_pure_props_wide'
                  ? [{ NAME: 'WATER', MW: 18.015, TC: 647.1, PC: 22064000, OMEGA: 0.344 }]
                  : table === 'apv140_plxant_wide'
                    ? [{ NAME: 'WATER', C1: 73.649, C2: -7258.2, C3: 0, C4: 0, C5: -7.3037, C6: 4.1653e-6, C7: 2, Tmin: 273.16, Tmax: 647.1 }]
                    : table === 'compound_properties'
                      ? [{ Name: 'WATER', MolecularWeight: { value: '18.015' }, ElectrolyteSolvent: { value: 'true' } }]
                : [],
            });
          },
        };
      },
    } as any;

    const search = await searchRemotePetroAssays('Example', client);
    const loaded = await loadRemotePetroAssay('A1', client);

    expect(search).toHaveLength(1);
    expect(search[0].assayName).toBe('Example Crude');
    expect(loaded.assay?.cuts).toHaveLength(2);
    expect(loaded.pseudoComponents).toHaveLength(2);
    expect(loaded.lightEndCompounds).toHaveLength(1);
    expect(loaded.lightEndCompounds[0].name).toBe('WATER');
    expect(loaded.pseudoComponents[0].molecularWeight).toBeLessThan(loaded.pseudoComponents[1].molecularWeight);
  });
});
