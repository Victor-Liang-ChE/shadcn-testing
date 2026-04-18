import type { PetroAssay } from './petro';

export interface PetroSupabaseHeaderRow {
  assay_id: string;
  assay_name: string;
  source: string;
  location: string;
  cut_count: number;
}

export interface PetroSupabaseCutRow {
  assay_id: string;
  cut_order: number;
  cut_name: string;
  tbp_start_k: number;
  tbp_end_k: number;
  liquid_volume_fraction: number;
  specific_gravity: number;
}

export interface PetroSupabaseLightEndRow {
  assay_id: string;
  component_name: string;
  fraction_value: number;
  fraction_basis: string;
  normal_boiling_point_k: number | '';
}

export interface PetroSupabaseSeedBundle {
  assays: PetroSupabaseHeaderRow[];
  cuts: PetroSupabaseCutRow[];
  lightEnds: PetroSupabaseLightEndRow[];
}

export interface PetroSupabaseSeedOptions {
  assayId: string;
  source: string;
  location?: string;
}

function escapeCsv(value: string | number | ''): string {
  if (typeof value === 'number') return Number.isFinite(value) ? String(value) : '';
  const text = String(value);
  if (text.includes(',') || text.includes('"') || text.includes('\n')) {
    return `"${text.replace(/"/g, '""')}"`;
  }
  return text;
}

export function buildPetroSupabaseSeedBundle(
  assay: PetroAssay,
  options: PetroSupabaseSeedOptions,
): PetroSupabaseSeedBundle {
  return {
    assays: [{
      assay_id: options.assayId,
      assay_name: assay.assayName,
      source: options.source,
      location: options.location ?? '',
      cut_count: assay.cuts.length,
    }],
    cuts: assay.cuts.map((cut, index) => ({
      assay_id: options.assayId,
      cut_order: index + 1,
      cut_name: cut.cutName?.trim() || `Cut ${index + 1}`,
      tbp_start_k: cut.tbpStart_K,
      tbp_end_k: cut.tbpEnd_K,
      liquid_volume_fraction: cut.liquidVolumeFraction,
      specific_gravity: cut.specificGravity,
    })),
    lightEnds: (assay.lightEnds ?? []).map(lightEnd => ({
      assay_id: options.assayId,
      component_name: lightEnd.componentName,
      fraction_value: lightEnd.fractionValue,
      fraction_basis: lightEnd.fractionBasis ?? 'mole',
      normal_boiling_point_k: lightEnd.normalBoilingPoint_K ?? '',
    })),
  };
}

export function serializeRowsToCsv<Row extends Record<string, string | number | ''>>(
  rows: Row[],
  headers: (keyof Row)[],
): string {
  const headerLine = headers.map(header => String(header)).join(',');
  const body = rows.map(row => headers.map(header => escapeCsv(row[header])).join(','));
  return [headerLine, ...body].join('\n');
}
