import { describe, expect, it } from 'vitest';

import { createCanopyAuthCookie, verifyCanopyAuthCookie } from '../canopy/server-auth';

describe('canopy server auth cookie', () => {
  it('creates and verifies a signed canopy auth cookie', () => {
    const cookie = createCanopyAuthCookie({
      id: 'user-123',
      email: 'user@example.com',
      fullName: 'Example User',
    });

    const verified = verifyCanopyAuthCookie(cookie);
    expect(verified).toEqual({
      id: 'user-123',
      email: 'user@example.com',
      fullName: 'Example User',
    });
  });

  it('rejects malformed cookies', () => {
    expect(verifyCanopyAuthCookie('bad-cookie')).toBeNull();
    expect(verifyCanopyAuthCookie('a.b')).toBeNull();
  });
});
