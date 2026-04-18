'use client';

import { useEffect, useRef, useState } from 'react';
import { useRouter } from 'next/navigation';

import { parseCanopyAuthCallback } from '@/lib/canopy/save-utils';

export default function CanopyAuthCallbackPage() {
  const router = useRouter();
  const attemptedRef = useRef(false);
  const [errorMessage, setErrorMessage] = useState<string | null>(null);

  useEffect(() => {
    if (attemptedRef.current) return;
    attemptedRef.current = true;

    const { accessToken, errorMessage: callbackError, nextPath } = parseCanopyAuthCallback(
      new URLSearchParams(typeof window !== 'undefined' ? window.location.search : ''),
      typeof window !== 'undefined' ? window.location.hash : '',
    );

    if (callbackError) {
      setErrorMessage(callbackError);
      return;
    }

    if (!accessToken) {
      setErrorMessage('Missing access token in callback URL fragment.');
      return;
    }

    if (typeof window !== 'undefined' && window.location.hash) {
      window.history.replaceState(null, '', window.location.pathname + window.location.search);
    }

    fetch('/api/canopy/auth/session', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ access_token: accessToken }),
    })
      .then(async (response) => {
        const payload = await response.json().catch(() => ({}));
        if (!response.ok) {
          throw new Error(payload.error || 'Failed to establish canopy session.');
        }
        router.replace(nextPath);
      })
      .catch((error: unknown) => {
        setErrorMessage(error instanceof Error ? error.message : 'OAuth callback failed.');
      });
  }, [router]);

  return (
    <div className="min-h-screen flex items-center justify-center bg-background px-6">
      <div className="w-full max-w-md rounded-lg border border-border bg-card p-6 shadow-sm">
        <h1 className="text-lg font-semibold mb-2">Canopy Sign-In</h1>
        {errorMessage ? (
          <div className="space-y-2">
            <p className="text-sm text-red-600 dark:text-red-400">{errorMessage}</p>
            <button
              type="button"
              onClick={() => router.replace('/canopy')}
              className="text-sm text-blue-600 hover:underline"
            >
              Return to Canopy
            </button>
          </div>
        ) : (
          <p className="text-sm text-muted-foreground">
            Finalizing sign-in and creating a secure canopy session.
          </p>
        )}
      </div>
    </div>
  );
}
