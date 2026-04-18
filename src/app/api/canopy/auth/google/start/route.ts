import { NextRequest, NextResponse } from 'next/server';

function normalizeNextPath(input: string | null): string {
  if (!input || !input.startsWith('/') || input.startsWith('//')) return '/canopy';
  return input;
}

export async function GET(request: NextRequest) {
  const supabaseUrl = process.env.SUPABASE_URL;
  if (!supabaseUrl) {
    return NextResponse.json({ error: 'SUPABASE_URL is missing.' }, { status: 500 });
  }

  const next = normalizeNextPath(request.nextUrl.searchParams.get('next'));
  const callbackUrl = new URL('/canopy/auth/callback', request.url);
  callbackUrl.searchParams.set('next', next);

  const redirectUrl = new URL('/auth/v1/authorize', supabaseUrl);
  redirectUrl.searchParams.set('provider', 'google');
  redirectUrl.searchParams.set('redirect_to', callbackUrl.toString());
  redirectUrl.searchParams.set('scopes', 'openid email profile');

  return NextResponse.redirect(redirectUrl);
}
