import { describe, expect, it } from 'vitest';

import { createExampleCrudeAssay } from '../canopy/petro';
import { buildPetroSupabaseSeedBundle, serializeRowsToCsv } from '../canopy/petro-supabase';

describe('canopy petro supabase seed helpers', () => {
  it('builds schema-ready assay and cut rows from a petro assay', () => {
    const assay = createExampleCrudeAssay();
    const bundle = buildPetroSupabaseSeedBundle(assay, {
      assayId: 'EXAMPLE_CRUDE',
      source: 'Canopy Built-In Example',
      location: 'Seed Data',
    });

    expect(bundle.assays).toEqual([{
      assay_id: 'EXAMPLE_CRUDE',
      assay_name: 'Example Crude',
      source: 'Canopy Built-In Example',
      location: 'Seed Data',
      cut_count: 5,
    }]);
    expect(bundle.cuts).toHaveLength(5);
    expect(bundle.cuts[0].cut_name).toBe('Naphtha');
    expect(bundle.cuts[4].specific_gravity).toBeCloseTo(0.99, 8);
    expect(bundle.lightEnds).toEqual([]);
  });

  it('serializes headers even when the optional light-end file is empty', () => {
    const csv = serializeRowsToCsv([], [
      'assay_id',
      'component_name',
      'fraction_value',
      'fraction_basis',
      'normal_boiling_point_k',
    ]);

    expect(csv).toBe('assay_id,component_name,fraction_value,fraction_basis,normal_boiling_point_k');
  });
});
