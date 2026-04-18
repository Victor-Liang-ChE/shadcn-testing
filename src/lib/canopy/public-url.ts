export interface ResolveCanopyPublicOriginOptions {
  configuredOrigin?: string | null;
  requestOrigin?: string | null;
  forwardedHost?: string | null;
  forwardedProto?: string | null;
  host?: string | null;
}

function normalizeOrigin(input: string | null | undefined): string | null {
  if (!input) return null;
  const trimmed = input.trim();
  if (!trimmed) return null;

  const withProtocol = /^https?:\/\//i.test(trimmed)
    ? trimmed
    : trimmed.includes('localhost')
      ? `http://${trimmed}`
      : `https://${trimmed}`;

  try {
    const url = new URL(withProtocol);
    return url.origin;
  } catch {
    return null;
  }
}

export function resolveCanopyPublicOrigin(options: ResolveCanopyPublicOriginOptions): string {
  const configured = normalizeOrigin(options.configuredOrigin);
  if (configured) return configured;

  const forwardedHost = options.forwardedHost?.trim() || options.host?.trim();
  if (forwardedHost) {
    const proto = options.forwardedProto?.trim() || (forwardedHost.includes('localhost') ? 'http' : 'https');
    const fromForwarded = normalizeOrigin(`${proto}://${forwardedHost}`);
    if (fromForwarded) return fromForwarded;
  }

  const requestOrigin = normalizeOrigin(options.requestOrigin);
  if (requestOrigin) return requestOrigin;

  return 'http://localhost:3000';
}
