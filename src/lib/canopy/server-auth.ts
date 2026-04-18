import crypto from 'crypto';

import { createClient } from '@supabase/supabase-js';
import { NextRequest, NextResponse } from 'next/server';

export const CANOPY_AUTH_COOKIE = 'canopy_auth';

export interface CanopyServerUser {
  id: string;
  email: string | null;
  fullName: string | null;
}

interface CanopyAuthCookiePayload {
  sub: string;
  email: string | null;
  fullName: string | null;
  exp: number;
}

function getSupabaseUrl(): string {
  const value = process.env.SUPABASE_URL;
  if (!value) throw new Error('SUPABASE_URL is missing.');
  return value;
}

function getSupabaseAnonKey(): string {
  const value = process.env.SUPABASE_ANON_KEY;
  if (!value) throw new Error('SUPABASE_ANON_KEY is missing.');
  return value;
}

function getAuthSecret(): string {
  return process.env.CANOPY_AUTH_SECRET || process.env.SUPABASE_SERVICE_KEY || process.env.SUPABASE_ANON_KEY || 'canopy-dev-secret';
}

function base64UrlEncode(input: string | Buffer): string {
  return Buffer.from(input).toString('base64url');
}

function base64UrlDecode(input: string): string {
  return Buffer.from(input, 'base64url').toString('utf8');
}

function signPayload(encodedPayload: string): string {
  return crypto.createHmac('sha256', getAuthSecret()).update(encodedPayload).digest('base64url');
}

export function createCanopyAuthCookie(user: CanopyServerUser, lifetimeDays = 30): string {
  const payload: CanopyAuthCookiePayload = {
    sub: user.id,
    email: user.email,
    fullName: user.fullName,
    exp: Math.floor(Date.now() / 1000) + lifetimeDays * 24 * 60 * 60,
  };
  const encodedPayload = base64UrlEncode(JSON.stringify(payload));
  const signature = signPayload(encodedPayload);
  return `${encodedPayload}.${signature}`;
}

export function verifyCanopyAuthCookie(rawCookie: string | undefined): CanopyServerUser | null {
  if (!rawCookie) return null;
  const [encodedPayload, signature] = rawCookie.split('.');
  if (!encodedPayload || !signature) return null;
  const expectedSignature = signPayload(encodedPayload);
  const signatureBuffer = Buffer.from(signature);
  const expectedBuffer = Buffer.from(expectedSignature);
  if (signatureBuffer.length !== expectedBuffer.length) return null;
  if (!crypto.timingSafeEqual(signatureBuffer, expectedBuffer)) return null;

  try {
    const payload = JSON.parse(base64UrlDecode(encodedPayload)) as CanopyAuthCookiePayload;
    if (!payload.sub || !payload.exp || payload.exp <= Math.floor(Date.now() / 1000)) return null;
    return {
      id: payload.sub,
      email: payload.email ?? null,
      fullName: payload.fullName ?? null,
    };
  } catch {
    return null;
  }
}

export function setCanopyAuthCookie(response: NextResponse, user: CanopyServerUser): void {
  response.cookies.set({
    name: CANOPY_AUTH_COOKIE,
    value: createCanopyAuthCookie(user),
    httpOnly: true,
    sameSite: 'lax',
    secure: process.env.NODE_ENV === 'production',
    path: '/',
    maxAge: 30 * 24 * 60 * 60,
  });
}

export function clearCanopyAuthCookie(response: NextResponse): void {
  response.cookies.set({
    name: CANOPY_AUTH_COOKIE,
    value: '',
    httpOnly: true,
    sameSite: 'lax',
    secure: process.env.NODE_ENV === 'production',
    path: '/',
    maxAge: 0,
  });
}

export function getCanopyUserFromRequest(request: NextRequest): CanopyServerUser | null {
  return verifyCanopyAuthCookie(request.cookies.get(CANOPY_AUTH_COOKIE)?.value);
}

export async function fetchSupabaseUserFromAccessToken(accessToken: string): Promise<CanopyServerUser | null> {
  const supabase = createClient(getSupabaseUrl(), getSupabaseAnonKey(), {
    auth: {
      persistSession: false,
      autoRefreshToken: false,
      detectSessionInUrl: false,
    },
  });
  const { data, error } = await supabase.auth.getUser(accessToken);
  if (error || !data.user) return null;
  return {
    id: data.user.id,
    email: data.user.email ?? null,
    fullName: (data.user.user_metadata?.full_name as string | undefined) ?? (data.user.user_metadata?.name as string | undefined) ?? null,
  };
}

export function createCanopyServiceClient() {
  const supabaseUrl = process.env.SUPABASE_URL;
  const supabaseServiceKey = process.env.SUPABASE_SERVICE_KEY;
  if (!supabaseUrl || !supabaseServiceKey) {
    throw new Error('SUPABASE_URL or SUPABASE_SERVICE_KEY is missing.');
  }
  return createClient(supabaseUrl, supabaseServiceKey, {
    auth: {
      persistSession: false,
      autoRefreshToken: false,
      detectSessionInUrl: false,
    },
  });
}
