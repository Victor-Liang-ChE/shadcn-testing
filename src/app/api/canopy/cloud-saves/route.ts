import { NextRequest, NextResponse } from 'next/server';

import { parseCanopyFlowsheetJson, sanitizeCloudSaveName } from '@/lib/canopy/save-utils';
import { createCanopyServiceClient, getCanopyUserFromRequest } from '@/lib/canopy/server-auth';

function unauthorizedResponse() {
  return NextResponse.json({ error: 'Authentication required.' }, { status: 401 });
}

export async function GET(request: NextRequest) {
  const user = getCanopyUserFromRequest(request);
  if (!user) return unauthorizedResponse();

  const supabase = createCanopyServiceClient();
  const { data, error } = await supabase
    .from('canopy_flowsheets')
    .select('id, name, created_at, updated_at')
    .eq('user_id', user.id)
    .order('updated_at', { ascending: false });

  if (error) {
    return NextResponse.json({ error: error.message }, { status: 500 });
  }
  return NextResponse.json({ saves: data ?? [] });
}

export async function POST(request: NextRequest) {
  const user = getCanopyUserFromRequest(request);
  if (!user) return unauthorizedResponse();

  const body = await request.json().catch(() => null) as { name?: string; flowsheetJson?: string } | null;
  const parsed = parseCanopyFlowsheetJson(body?.flowsheetJson ?? '');
  if (!parsed.ok) {
    return NextResponse.json({ error: parsed.error }, { status: 400 });
  }

  const supabase = createCanopyServiceClient();
  const { data, error } = await supabase
    .from('canopy_flowsheets')
    .insert({
      user_id: user.id,
      name: sanitizeCloudSaveName(body?.name ?? 'Canopy Flowsheet'),
      flowsheet_json: parsed.data,
    })
    .select('id, name, created_at, updated_at')
    .single();

  if (error) {
    return NextResponse.json({ error: error.message }, { status: 500 });
  }
  return NextResponse.json({ save: data });
}
