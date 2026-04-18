import { describe, expect, it } from 'vitest';

import {
  buildFlowsheetDownloadFilename,
  normalizeCanopyRedirectPath,
  parseCanopyAuthCallback,
  parseCanopyFlowsheetJson,
  sanitizeCloudSaveName,
} from '../canopy/save-utils';

describe('canopy save utilities', () => {
  it('sanitizes cloud save names consistently', () => {
    expect(sanitizeCloudSaveName('   My    Column   Case   ')).toBe('My Column Case');
    expect(sanitizeCloudSaveName('')).toBe('Canopy Flowsheet');
  });

  it('builds deterministic download filenames', () => {
    const fileName = buildFlowsheetDownloadFilename('Cumene Case 01', new Date('2026-04-17T15:04:00Z'));
    expect(fileName).toBe('cumene-case-01-20260417-0804.json');
  });

  it('parses canopy flowsheet JSON shape', () => {
    expect(parseCanopyFlowsheetJson('{"version":1,"compounds":[],"streams":{},"units":{}}')).toEqual({
      ok: true,
      data: { version: 1, compounds: [], streams: {}, units: {} },
    });
    expect(parseCanopyFlowsheetJson('{"hello":"world"}')).toEqual({
      ok: false,
      error: 'File is not a valid Canopy flowsheet export.',
    });
  });

  it('normalizes post-auth redirect paths defensively', () => {
    expect(normalizeCanopyRedirectPath('/canopy')).toBe('/canopy');
    expect(normalizeCanopyRedirectPath('//evil.example')).toBe('/canopy');
    expect(normalizeCanopyRedirectPath('https://evil.example')).toBe('/canopy');
    expect(normalizeCanopyRedirectPath(null)).toBe('/canopy');
  });

  it('parses an implicit-flow supabase callback fragment', () => {
    expect(
      parseCanopyAuthCallback(
        new URLSearchParams('next=%2Fcanopy'),
        '#access_token=test-token&token_type=bearer',
      ),
    ).toEqual({
      accessToken: 'test-token',
      errorMessage: null,
      nextPath: '/canopy',
    });
  });

  it('surfaces callback errors from either query or fragment params', () => {
    expect(
      parseCanopyAuthCallback(
        new URLSearchParams('next=%2Fcanopy&error_description=Denied'),
        '',
      ),
    ).toEqual({
      accessToken: null,
      errorMessage: 'Denied',
      nextPath: '/canopy',
    });

    expect(
      parseCanopyAuthCallback(
        new URLSearchParams('next=%2Fcanopy'),
        '#error=access_denied',
      ),
    ).toEqual({
      accessToken: null,
      errorMessage: 'access_denied',
      nextPath: '/canopy',
    });
  });
});
