import { NextRequest, NextResponse } from 'next/server';

import {
  clearCanopyAuthCookie,
  fetchSupabaseUserFromAccessToken,
  getCanopyUserFromRequest,
  setCanopyAuthCookie,
} from '@/lib/canopy/server-auth';

export async function GET(request: NextRequest) {
  const user = getCanopyUserFromRequest(request);
  return NextResponse.json({ authenticated: Boolean(user), user });
}

export async function POST(request: NextRequest) {
  const body = await request.json().catch(() => null) as { access_token?: string } | null;
  const accessToken = body?.access_token;
  if (!accessToken) {
    return NextResponse.json({ error: 'Missing access token.' }, { status: 400 });
  }

  const user = await fetchSupabaseUserFromAccessToken(accessToken);
  if (!user) {
    return NextResponse.json({ error: 'Supabase token validation failed.' }, { status: 401 });
  }

  const response = NextResponse.json({ authenticated: true, user });
  setCanopyAuthCookie(response, user);
  return response;
}

export async function DELETE() {
  const response = NextResponse.json({ authenticated: false });
  clearCanopyAuthCookie(response);
  return response;
}
