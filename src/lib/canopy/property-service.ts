import type { SupabaseClient } from '@supabase/supabase-js';

import { formatCompoundName } from '@/lib/antoine-utils';
import { supabase } from '@/lib/supabaseClient';

import { buildInteractionParamsForPackage, type AspenBinaryParamRow } from './interaction-params';
import { generatePseudoComponentsFromAssay, type PetroAssay } from './petro';
import { UNIFAC_DORTMUND_PARAMS } from './store';
import type { CanopyCompound, FluidPackage } from './types';
import type {
  ElecNRTLParams,
  InteractionParams,
  UNIFACData,
  UNIFACSubgroup,
} from './thermo';

export interface PropertyProvenanceRecord {
  label: string;
  source: string;
  detail: string;
}

export interface RemoteCompoundLoadResult {
  compound: CanopyCompound | null;
  provenance: PropertyProvenanceRecord[];
  warnings: string[];
}

export interface RemoteInteractionLoadResult {
  params: Partial<InteractionParams> | null;
  provenance: PropertyProvenanceRecord[];
  warnings: string[];
}

export interface RemotePetroAssaySummary {
  assayId: string;
  assayName: string;
  source: string;
  location?: string;
  cutCount: number;
}

export interface RemotePetroAssayLoadResult {
  assay: PetroAssay | null;
  pseudoComponents: CanopyCompound[];
  lightEndCompounds: CanopyCompound[];
  provenance: PropertyProvenanceRecord[];
  warnings: string[];
}

type JsonLike = Record<string, any>;

function getNumericValue(...values: unknown[]): number | undefined {
  for (const value of values) {
    if (value == null) continue;
    if (typeof value === 'number' && Number.isFinite(value)) return value;
    if (typeof value === 'string') {
      const parsed = parseFloat(value);
      if (Number.isFinite(parsed)) return parsed;
    }
    if (typeof value === 'object' && value && 'value' in (value as JsonLike)) {
      const nested = getNumericValue((value as JsonLike).value);
      if (nested !== undefined) return nested;
    }
  }
  return undefined;
}

function getBooleanValue(...values: unknown[]): boolean | undefined {
  for (const value of values) {
    if (typeof value === 'boolean') return value;
    if (typeof value === 'number') return value !== 0;
    if (typeof value === 'string') {
      const normalized = value.trim().toLowerCase();
      if (['true', 'yes', 'y', '1'].includes(normalized)) return true;
      if (['false', 'no', 'n', '0'].includes(normalized)) return false;
    }
    if (typeof value === 'object' && value && 'value' in (value as JsonLike)) {
      const nested = getBooleanValue((value as JsonLike).value);
      if (nested !== undefined) return nested;
    }
  }
  return undefined;
}

function buildDipprCoefficients(properties: JsonLike | undefined): Record<string, { A: number; B: number; C: number; D: number; E?: number }> {
  const dipprCoeffs: Record<string, { A: number; B: number; C: number; D: number; E?: number }> = {};
  if (!properties) return dipprCoeffs;

  for (const key of ['IdealGasHeatCapacityCp', 'HeatOfVaporization', 'LiquidHeatCapacityCp']) {
    const propData = properties[key];
    if (propData?.A?.value !== undefined) {
      dipprCoeffs[key] = {
        A: parseFloat(propData.A.value),
        B: parseFloat(propData.B?.value ?? '0'),
        C: parseFloat(propData.C?.value ?? '0'),
        D: parseFloat(propData.D?.value ?? '0'),
        E: propData.E?.value ? parseFloat(propData.E.value) : undefined,
      };
    }
  }

  return dipprCoeffs;
}

function buildTransportCoefficients(pure: JsonLike, properties: JsonLike | undefined): CanopyCompound['transportCoeffs'] {
  const transportCoeffs: CanopyCompound['transportCoeffs'] = {};
  const transportMap: Record<string, keyof NonNullable<CanopyCompound['transportCoeffs']>> = {
    LiquidViscosity: 'LiquidViscosity',
    VaporViscosity: 'VaporViscosity',
    LiquidThermalConductivity: 'LiquidThermalConductivity',
    VaporThermalConductivity: 'VaporThermalConductivity',
    SurfaceTension: 'SurfaceTension',
  };

  if (properties) {
    for (const [jsonKey, tsKey] of Object.entries(transportMap)) {
      const propData = properties[jsonKey];
      if (propData?.A?.value !== undefined) {
        (transportCoeffs as any)[tsKey] = {
          A: parseFloat(propData.A.value),
          B: parseFloat(propData.B?.value ?? '0'),
          C: parseFloat(propData.C?.value ?? '0'),
          D: parseFloat(propData.D?.value ?? '0'),
          E: propData.E?.value ? parseFloat(propData.E.value) : undefined,
        };
      }
    }
  }

  if (!transportCoeffs.LiquidViscosity && pure.MULDIP_A != null) {
    transportCoeffs.LiquidViscosity = {
      A: Number(pure.MULDIP_A),
      B: Number(pure.MULDIP_B ?? 0),
      C: Number(pure.MULDIP_C ?? 0),
      D: Number(pure.MULDIP_D ?? 0),
      E: pure.MULDIP_E != null ? Number(pure.MULDIP_E) : undefined,
    };
  }
  if (!transportCoeffs.SurfaceTension && pure.SIGDIP_A != null) {
    transportCoeffs.SurfaceTension = {
      A: Number(pure.SIGDIP_A),
      B: Number(pure.SIGDIP_B ?? 0),
      C: Number(pure.SIGDIP_C ?? 0),
      D: Number(pure.SIGDIP_D ?? 0),
      E: pure.SIGDIP_E != null ? Number(pure.SIGDIP_E) : undefined,
    };
  }
  if (!transportCoeffs.LiquidThermalConductivity && pure.KLDIP_A != null) {
    transportCoeffs.LiquidThermalConductivity = {
      A: Number(pure.KLDIP_A),
      B: Number(pure.KLDIP_B ?? 0),
      C: Number(pure.KLDIP_C ?? 0),
      D: Number(pure.KLDIP_D ?? 0),
      E: pure.KLDIP_E != null ? Number(pure.KLDIP_E) : undefined,
    };
  }
  if (!transportCoeffs.VaporViscosity && pure.MUVDIP_A != null) {
    transportCoeffs.VaporViscosity = {
      A: Number(pure.MUVDIP_A),
      B: Number(pure.MUVDIP_B ?? 0),
      C: Number(pure.MUVDIP_C ?? 0),
      D: Number(pure.MUVDIP_D ?? 0),
      E: pure.MUVDIP_E != null ? Number(pure.MUVDIP_E) : undefined,
    };
  }

  return Object.keys(transportCoeffs).length > 0 ? transportCoeffs : undefined;
}

function buildUnifacGroups(properties: JsonLike | undefined): Record<number, number> | undefined {
  if (!properties?.UnifacVLE?.group) return undefined;
  const unifacGroups: Record<number, number> = {};
  const groups = Array.isArray(properties.UnifacVLE.group)
    ? properties.UnifacVLE.group
    : [properties.UnifacVLE.group];

  for (const group of groups) {
    const id = parseInt(group.id ?? group['@id'] ?? group.groupId, 10);
    const count = parseInt(group.count ?? group['@count'] ?? group.value ?? '1', 10);
    if (!Number.isNaN(id) && !Number.isNaN(count)) {
      unifacGroups[id] = count;
    }
  }

  return Object.keys(unifacGroups).length > 0 ? unifacGroups : undefined;
}

export function buildCanopyCompoundFromRemoteRecords(
  name: string,
  pure: JsonLike,
  ant: JsonLike | null | undefined,
  properties: JsonLike | undefined,
): CanopyCompound {
  const dipprCoeffs = buildDipprCoefficients(properties);
  const molecularWeight = getNumericValue(properties?.MolecularWeight, pure.MW) ?? 0;
  const uniquac_r = getNumericValue(pure.GMUQR);
  const uniquac_q = getNumericValue(pure.GMUQQ);
  const Hf_298_Jmol = getNumericValue(properties?.HeatOfFormation) != null
    ? getNumericValue(properties?.HeatOfFormation)! / 1000
    : undefined;
  const Gf_298_Jmol = getNumericValue(properties?.GibbsEnergyOfFormation) != null
    ? getNumericValue(properties?.GibbsEnergyOfFormation)! / 1000
    : undefined;
  const cpigdp = pure.CPIGDP_A != null
    ? {
        A: Number(pure.CPIGDP_A),
        B: Number(pure.CPIGDP_B ?? 0),
        C: Number(pure.CPIGDP_C ?? 0),
        D: Number(pure.CPIGDP_D ?? 0),
        E: Number(pure.CPIGDP_E ?? 0),
      }
    : undefined;
  const transportCoeffs = buildTransportCoefficients(pure, properties);
  const unifacGroups = buildUnifacGroups(properties);
  const chargeNumber = getNumericValue(
    properties?.ChargeNumber,
    properties?.IonicCharge,
    properties?.SpeciesCharge,
    properties?.FormalCharge,
    properties?.NetCharge,
  );
  const dielectricE1 = getNumericValue(
    properties?.DielectricConstant,
    properties?.RelativePermittivity,
    properties?.DielectricConstant25C,
  );
  const dielectricE2 = getNumericValue(properties?.DielectricConstantB);
  const dielectricE3 = getNumericValue(properties?.DielectricConstantRefTemp, 298.15);
  const isElectrolyteSolvent = getBooleanValue(properties?.ElectrolyteSolvent, properties?.SolventFlag)
    ?? name.toUpperCase() === 'WATER';
  const specificGravity = getNumericValue(
    properties?.SpecificGravity,
    properties?.SG,
    properties?.SG60_60,
    properties?.SGAIR,
    pure.SG,
    pure.SG60_60,
  );
  const flashPoint_K = getNumericValue(
    properties?.FlashPoint,
    properties?.FLASHPOINT,
    properties?.FlashPointClosedCup,
    properties?.ClosedCupFlashPoint,
    pure.FLASHPOINT,
  );
  const autoignitionTemp_K = getNumericValue(
    properties?.AIT,
    properties?.AutoIgnitionTemperature,
    properties?.AutoignitionTemperature,
    properties?.AUTOIGNITIONTEMP,
    pure.AIT,
  );
  const refractiveIndex = getNumericValue(
    properties?.REFINDEX,
    properties?.REFI,
    properties?.RefractiveIndex,
    properties?.NDEX,
    pure.REFI,
  );
  const anilinePoint_K = getNumericValue(
    properties?.ANILPT,
    properties?.ANILINEPT,
    properties?.AnilinePoint,
    pure.ANILPT,
  );
  const cetaneNumber = getNumericValue(
    properties?.CETANE,
    properties?.CetaneNumber,
    pure.CETANE,
  );
  const Hcomb_Jmol = getNumericValue(
    properties?.HigherHeatingValue,
    properties?.GrossHeatingValue,
    properties?.HeatOfCombustion,
    properties?.HHV,
    properties?.GHV,
    pure.HCOMB,
  );
  const hCombustion_Jpkmol = getNumericValue(
    properties?.NetHeatingValue,
    properties?.LowerHeatingValue,
    properties?.HeatOfCombustionNet,
    properties?.NHV,
    properties?.LHV,
    pure.NCOMB,
  );

  return {
    name,
    displayName: formatCompoundName(name),
    molecularWeight,
    Tc_K: Number(pure.TC ?? 0),
    Pc_Pa: Number(pure.PC ?? 0),
    omega: Number(pure.OMEGA ?? 0),
    dipprCoeffs,
    antoine: ant
      ? {
          C1: ant.C1 ?? 0,
          C2: ant.C2 ?? 0,
          C3: ant.C3 ?? 0,
          C4: ant.C4 ?? 0,
          C5: ant.C5 ?? 0,
          C6: ant.C6 ?? 0,
          C7: ant.C7 ?? 1,
          Tmin_K: ant.Tmin ?? 0,
          Tmax_K: ant.Tmax ?? 10000,
        }
      : null,
    ...(uniquac_r !== undefined ? { uniquac_r } : {}),
    ...(uniquac_q !== undefined ? { uniquac_q } : {}),
    ...(unifacGroups ? { unifacGroups } : {}),
    ...(transportCoeffs ? { transportCoeffs } : {}),
    ...(cpigdp ? { cpigdp } : {}),
    ...(Hf_298_Jmol !== undefined ? { Hf_298_Jmol } : {}),
    ...(Gf_298_Jmol !== undefined ? { Gf_298_Jmol } : {}),
    ...(chargeNumber !== undefined ? { chargeNumber } : {}),
    ...(dielectricE1 !== undefined ? { dielectricConst: { e1: dielectricE1, e2: dielectricE2 ?? 0, e3_K: dielectricE3 ?? 298.15 } } : {}),
    ...(isElectrolyteSolvent ? { isElectrolyteSolvent } : {}),
    ...(specificGravity !== undefined ? { specificGravity } : {}),
    ...(flashPoint_K !== undefined ? { flashPoint_K } : {}),
    ...(autoignitionTemp_K !== undefined ? { autoignitionTemp_K } : {}),
    ...(refractiveIndex !== undefined ? { refractiveIndex } : {}),
    ...(anilinePoint_K !== undefined ? { anilinePoint_K } : {}),
    ...(cetaneNumber !== undefined ? { cetaneNumber } : {}),
    ...(Hcomb_Jmol !== undefined ? { Hcomb_Jmol } : {}),
    ...(hCombustion_Jpkmol !== undefined ? { hCombustion_Jpkmol } : {}),
  } as CanopyCompound;
}

function deriveElectrolyteParams(
  compounds: CanopyCompound[],
  params: Partial<InteractionParams>,
): { elecNrtl?: ElecNRTLParams; warnings: string[] } {
  const warnings: string[] = [];
  if (!params.nrtl) return { warnings: ['Electrolyte-NRTL requires an NRTL short-range parameter set.'] };

  const solvent = compounds.find(c => c.isElectrolyteSolvent) ?? compounds.find(c => c.name.toUpperCase() === 'WATER');
  if (!solvent) {
    warnings.push('Electrolyte-NRTL could not identify a solvent compound. Mark a solvent compound before using this package.');
    return { warnings };
  }

  const epsilon_r = solvent.dielectricConst?.e1 ?? (solvent.name.toUpperCase() === 'WATER' ? 78.36 : undefined);
  if (epsilon_r === undefined) {
    warnings.push(`Electrolyte-NRTL is missing dielectric constant data for solvent ${solvent.displayName}.`);
    return { warnings };
  }

  const missingCharges = compounds.filter(c => c.chargeNumber == null).map(c => c.displayName);
  if (missingCharges.length > 0) {
    warnings.push(`Electrolyte-NRTL has no charge metadata for: ${missingCharges.join(', ')}. Unspecified species are treated as neutral.`);
  }

  return {
    elecNrtl: {
      nrtl: params.nrtl,
      charges: compounds.map(c => c.chargeNumber ?? 0),
      epsilon_r,
      solventMW: solvent.molecularWeight,
    },
    warnings,
  };
}

export async function searchRemoteCompounds(
  query: string,
  selectedNames: string[] = [],
  client: SupabaseClient = supabase,
): Promise<string[]> {
  if (query.trim().length < 2) return [];

  const { data } = await client
    .from('compounds_master')
    .select('"Name"')
    .ilike('Name', `${query.trim()}%`)
    .limit(40);

  if (!data) return [];

  const selected = new Set(selectedNames);
  const seen = new Set<string>();
  const names: string[] = [];
  for (const row of data as any[]) {
    const name = row.Name as string;
    if (!seen.has(name) && !selected.has(name)) {
      seen.add(name);
      names.push(name);
    }
    if (names.length >= 8) break;
  }
  return names;
}

export async function loadRemoteCompound(
  name: string,
  client: SupabaseClient = supabase,
): Promise<RemoteCompoundLoadResult> {
  const provenance: PropertyProvenanceRecord[] = [
    { label: 'Component Search', source: 'compounds_master', detail: name },
    { label: 'Pure Properties', source: 'apv140_pure_props_wide', detail: name },
    { label: 'Vapor Pressure', source: 'apv140_plxant_wide', detail: name },
    { label: 'Compound Properties', source: 'compound_properties', detail: name },
  ];

  const [pureRes, antoineRes, propsRes] = await Promise.all([
    client.from('apv140_pure_props_wide').select('*').ilike('CompoundName', name).limit(1),
    client.from('apv140_plxant_wide').select('*').ilike('CompoundName', name).limit(1),
    client.from('compound_properties').select('properties').ilike('name', name).limit(1),
  ]);

  const pure = pureRes.data?.[0] as any;
  const ant = antoineRes.data?.[0] as any;
  const props = (propsRes.data?.[0] as any)?.properties;

  if (!pure) {
    return {
      compound: null,
      provenance,
      warnings: [`Compound "${name}" was not found in the remote property tables.`],
    };
  }

  return {
    compound: buildCanopyCompoundFromRemoteRecords(name, pure, ant, props),
    provenance,
    warnings: [],
  };
}

export async function loadRemoteInteractionParams(
  pkg: FluidPackage,
  compounds: CanopyCompound[],
  client: SupabaseClient = supabase,
): Promise<RemoteInteractionLoadResult> {
  const warnings: string[] = [];
  const provenance: PropertyProvenanceRecord[] = [];
  if (compounds.length < 2 || pkg === 'Ideal') {
    return { params: null, provenance, warnings };
  }

  if (pkg === 'UNIFAC' || pkg === 'UNIFAC-DMD') {
    const allSubgroupIds = new Set<number>();
    const compGroups: Record<number, number>[] = [];
    for (const compound of compounds) {
      const groups = compound.unifacGroups ?? {};
      compGroups.push(groups);
      for (const subgroupId of Object.keys(groups)) allSubgroupIds.add(Number(subgroupId));
    }
    if (allSubgroupIds.size === 0) {
      warnings.push(`${pkg} needs UNIFAC subgroup definitions on each selected compound.`);
      return { params: null, provenance, warnings };
    }

    provenance.push({
      label: 'UNIFAC Subgroups',
      source: 'UNIFAC - Rk and Qk',
      detail: `${allSubgroupIds.size} subgroup records requested`,
    });

    const { data: subgroupRows } = await client
      .from('UNIFAC - Rk and Qk')
      .select('"Subgroup #", "Main Group #", "Rk", "Qk"')
      .in('"Subgroup #"', Array.from(allSubgroupIds));

    if (!subgroupRows || subgroupRows.length === 0) {
      warnings.push(`${pkg} subgroup definitions could not be loaded from the remote databank.`);
      return { params: null, provenance, warnings };
    }

    const subgroups: UNIFACSubgroup[] = [];
    const mainIds = new Set<number>();
    for (const row of subgroupRows as any[]) {
      const subgroupId = Number(row['Subgroup #']);
      const mainGroupId = Number(row['Main Group #']);
      subgroups.push({ subgroupId, mainGroupId, R: Number(row.Rk), Q: Number(row.Qk) });
      mainIds.add(mainGroupId);
    }

    if (pkg === 'UNIFAC') {
      const mainArr = Array.from(mainIds);
      provenance.push({
        label: 'UNIFAC Interactions',
        source: 'UNIFAC - a(ij)',
        detail: `${mainArr.length} main groups`,
      });

      const { data: intRows } = await client
        .from('UNIFAC - a(ij)')
        .select('i,j,"a(ij)"')
        .in('i', mainArr)
        .in('j', mainArr);

      const interactions: Record<string, number> = {};
      if (intRows) {
        for (const row of intRows as any[]) {
          interactions[`${row.i}-${row.j}`] = Number(row['a(ij)']) || 0;
        }
      }
      const data: UNIFACData = { subgroups, interactions };
      return { params: { unifac: { compGroups, data } }, provenance, warnings };
    }

    provenance.push({
      label: 'Dortmund UNIFAC Interactions',
      source: 'UNIFAC_DORTMUND_PARAMS',
      detail: `${Object.keys(UNIFAC_DORTMUND_PARAMS).length} built-in interaction pairs`,
    });

    const interactions: Record<string, { a: number; b: number; c: number }> = {};
    for (const [key, value] of Object.entries(UNIFAC_DORTMUND_PARAMS)) {
      interactions[key.replace('|', '-')] = value;
    }
    const coveredMainGroups = new Set<number>();
    for (const key of Object.keys(UNIFAC_DORTMUND_PARAMS)) {
      const [i, j] = key.split('|').map(Number);
      coveredMainGroups.add(i);
      coveredMainGroups.add(j);
    }
    const unsupported = Array.from(mainIds).filter(id => !coveredMainGroups.has(id));
    if (unsupported.length > 0) {
      warnings.push(`Dortmund UNIFAC will fall back to zero interaction for main groups ${unsupported.join(', ')}.`);
    }

    return { params: { unifacDmd: { compGroups, data: { subgroups, interactions } } }, provenance, warnings };
  }

  const propMap: Partial<Record<FluidPackage, string>> = {
    NRTL: 'NRTL',
    Wilson: 'WILSON',
    UNIQUAC: 'UNIQ',
    'Electrolyte-NRTL': 'NRTL',
    'Peng-Robinson': 'PRKA',
    SRK: 'SRKA',
  };
  const dbProp = propMap[pkg];
  if (!dbProp) return { params: null, provenance, warnings };

  const names = compounds.map(c => c.name);
  provenance.push({
    label: 'Binary Interaction Parameters',
    source: 'apv140_binary_params',
    detail: `${dbProp} for ${names.join(', ')}`,
  });

  const { data } = await client
    .from('apv140_binary_params')
    .select('Compound1Name,Compound2Name,ElementLabel,Value')
    .in('Compound1Name', names)
    .in('Compound2Name', names)
    .ilike('Property', `${dbProp}%`)
    .order('Compound1Name')
    .limit(1000);

  if (!data || data.length === 0) {
    warnings.push(`No ${pkg} binary interaction parameters were found for the selected compounds.`);
    return { params: null, provenance, warnings };
  }

  const parsed = buildInteractionParamsForPackage(
    pkg === 'Electrolyte-NRTL' ? 'NRTL' : pkg,
    compounds,
    data as AspenBinaryParamRow[],
  );
  if (!parsed) return { params: null, provenance, warnings };

  if (pkg === 'Electrolyte-NRTL') {
    const elec = deriveElectrolyteParams(compounds, parsed);
    warnings.push(...elec.warnings);
    if (!elec.elecNrtl) return { params: null, provenance, warnings };
    return { params: { nrtl: parsed.nrtl, elecNrtl: elec.elecNrtl }, provenance, warnings };
  }

  return { params: parsed, provenance, warnings };
}

export async function searchRemotePetroAssays(
  query: string,
  client: SupabaseClient = supabase,
): Promise<RemotePetroAssaySummary[]> {
  if (query.trim().length < 2) return [];
  const { data } = await client
    .from('petro_assays')
    .select('assay_id, assay_name, source, location, cut_count')
    .ilike('assay_name', `%${query.trim()}%`)
    .limit(20);

  return (data ?? []).map((row: any) => ({
    assayId: String(row.assay_id),
    assayName: String(row.assay_name),
    source: String(row.source ?? 'unknown'),
    location: row.location ? String(row.location) : undefined,
    cutCount: Number(row.cut_count ?? 0),
  }));
}

export async function loadRemotePetroAssay(
  assayId: string,
  client: SupabaseClient = supabase,
): Promise<RemotePetroAssayLoadResult> {
  const provenance: PropertyProvenanceRecord[] = [
    { label: 'Petroleum Assay Header', source: 'petro_assays', detail: assayId },
    { label: 'Petroleum Assay Cuts', source: 'petro_assay_cuts', detail: assayId },
  ];
  const warnings: string[] = [];

  const [assayRes, cutRes, lightEndRes] = await Promise.all([
    client.from('petro_assays').select('*').eq('assay_id', assayId).limit(1),
    client.from('petro_assay_cuts').select('*').eq('assay_id', assayId).order('cut_order'),
    client.from('petro_assay_light_ends').select('*').eq('assay_id', assayId).order('component_name'),
  ]);

  const assayRow = assayRes.data?.[0] as any;
  const cutRows = (cutRes.data ?? []) as any[];
  const lightEndRows = (lightEndRes.data ?? []) as any[];
  if (!assayRow || cutRows.length === 0) {
    return {
      assay: null,
      pseudoComponents: [],
      lightEndCompounds: [],
      provenance,
      warnings: [`Petroleum assay "${assayId}" was not found in the remote petro tables.`],
    };
  }

  const assay: PetroAssay = {
    assayName: String(assayRow.assay_name),
    cuts: cutRows.map(row => ({
      cutName: String(row.cut_name ?? `Cut ${row.cut_order ?? ''}`),
      tbpStart_K: Number(row.tbp_start_k),
      tbpEnd_K: Number(row.tbp_end_k),
      liquidVolumeFraction: Number(row.liquid_volume_fraction),
      specificGravity: Number(row.specific_gravity),
    })),
    lightEnds: lightEndRows.map(row => ({
      componentName: String(row.component_name),
      fractionValue: Number(row.fraction_value ?? 0),
      fractionBasis: row.fraction_basis ? String(row.fraction_basis) : undefined,
      normalBoilingPoint_K: row.normal_boiling_point_k != null ? Number(row.normal_boiling_point_k) : undefined,
    })),
  };

  let lightEndCompounds: CanopyCompound[] = [];
  if (lightEndRows.length > 0) {
    provenance.push({ label: 'Petroleum Light Ends', source: 'petro_assay_light_ends', detail: assayId });
    const uniqueLightEndNames = Array.from(
      new Set(
        lightEndRows
          .map(row => row.component_name)
          .filter((value): value is string => typeof value === 'string' && value.trim().length > 0)
          .map(value => value.trim().toUpperCase()),
      ),
    );

    const lightEndResults = await Promise.all(
      uniqueLightEndNames.map(name => loadRemoteCompound(name, client)),
    );

    lightEndCompounds = lightEndResults
      .map(result => result.compound)
      .filter((compound): compound is CanopyCompound => Boolean(compound));

    for (const result of lightEndResults) {
      provenance.push(...result.provenance);
      warnings.push(...result.warnings);
    }

    const missingLightEnds = uniqueLightEndNames.filter(name => !lightEndCompounds.some(compound => compound.name === name));
    if (missingLightEnds.length > 0) {
      warnings.push(`Light-end compounds were listed for assay "${assayId}" but not found in the remote compound tables: ${missingLightEnds.join(', ')}.`);
    }
  } else {
    warnings.push('No explicit light-end records were found for this assay; only pseudo-components were generated.');
  }

  return {
    assay,
    pseudoComponents: generatePseudoComponentsFromAssay(assay),
    lightEndCompounds,
    provenance,
    warnings,
  };
}
