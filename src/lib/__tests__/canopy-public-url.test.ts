import { describe, expect, it } from 'vitest';

import { resolveCanopyPublicOrigin } from '../canopy/public-url';

describe('canopy public origin resolution', () => {
  it('prefers the configured canonical site URL when present', () => {
    expect(resolveCanopyPublicOrigin({
      configuredOrigin: 'https://www.victorliang.com',
      requestOrigin: 'https://preview.vercel.app',
      forwardedHost: 'preview.vercel.app',
      forwardedProto: 'https',
    })).toBe('https://www.victorliang.com');
  });

  it('falls back to forwarded production headers and keeps localhost http', () => {
    expect(resolveCanopyPublicOrigin({
      configuredOrigin: null,
      requestOrigin: 'https://www.victorliang.com',
      forwardedHost: 'www.victorliang.com',
      forwardedProto: 'https',
    })).toBe('https://www.victorliang.com');

    expect(resolveCanopyPublicOrigin({
      configuredOrigin: null,
      requestOrigin: 'http://localhost:3000',
      forwardedHost: 'localhost:3000',
      forwardedProto: null,
    })).toBe('http://localhost:3000');
  });
});
