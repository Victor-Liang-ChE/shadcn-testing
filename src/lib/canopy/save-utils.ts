export const CANOPY_BROWSER_DRAFT_KEY = 'canopy-browser-draft';
export const CANOPY_BROWSER_DRAFT_META_KEY = 'canopy-browser-draft-meta';

export interface CanopyDraftMetadata {
  name: string;
  savedAt: string;
}

export interface CanopyAuthCallbackData {
  accessToken: string | null;
  errorMessage: string | null;
  nextPath: string;
}

export function sanitizeCloudSaveName(name: string): string {
  const trimmed = name.replace(/\s+/g, ' ').trim();
  return trimmed.length > 0 ? trimmed.slice(0, 120) : 'Canopy Flowsheet';
}

export function sanitizeDownloadStem(name: string): string {
  const safe = sanitizeCloudSaveName(name)
    .toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-+|-+$/g, '');
  return safe || 'canopy-flowsheet';
}

export function buildFlowsheetDownloadFilename(name: string, now = new Date()): string {
  const stamp = [
    now.getFullYear(),
    String(now.getMonth() + 1).padStart(2, '0'),
    String(now.getDate()).padStart(2, '0'),
    '-',
    String(now.getHours()).padStart(2, '0'),
    String(now.getMinutes()).padStart(2, '0'),
  ].join('');
  return `${sanitizeDownloadStem(name)}-${stamp}.json`;
}

export function parseCanopyFlowsheetJson(json: string): { ok: true; data: any } | { ok: false; error: string } {
  try {
    const data = JSON.parse(json);
    if (!data || typeof data !== 'object') {
      return { ok: false, error: 'File does not contain a JSON object.' };
    }
    if (!('version' in data) || !('compounds' in data) || !('streams' in data) || !('units' in data)) {
      return { ok: false, error: 'File is not a valid Canopy flowsheet export.' };
    }
    return { ok: true, data };
  } catch {
    return { ok: false, error: 'File is not valid JSON.' };
  }
}

export function normalizeCanopyRedirectPath(nextPath: string | null | undefined): string {
  if (!nextPath || !nextPath.startsWith('/')) return '/canopy';
  if (nextPath.startsWith('//')) return '/canopy';
  return nextPath;
}

export function parseCanopyAuthCallback(searchParams: URLSearchParams, hashFragment: string): CanopyAuthCallbackData {
  const nextPath = normalizeCanopyRedirectPath(searchParams.get('next'));
  const queryError = searchParams.get('error_description') ?? searchParams.get('error');
  if (queryError) {
    return {
      accessToken: null,
      errorMessage: queryError,
      nextPath,
    };
  }

  const hashParams = new URLSearchParams(hashFragment.replace(/^#/, ''));
  const hashError = hashParams.get('error_description') ?? hashParams.get('error');
  if (hashError) {
    return {
      accessToken: null,
      errorMessage: hashError,
      nextPath,
    };
  }

  return {
    accessToken: hashParams.get('access_token'),
    errorMessage: null,
    nextPath,
  };
}
