import { NextRequest, NextResponse } from 'next/server';

import { parseCanopyFlowsheetJson, sanitizeCloudSaveName } from '@/lib/canopy/save-utils';
import { createCanopyServiceClient, getCanopyUserFromRequest } from '@/lib/canopy/server-auth';

function unauthorizedResponse() {
  return NextResponse.json({ error: 'Authentication required.' }, { status: 401 });
}

export async function GET(request: NextRequest, context: { params: Promise<{ id: string }> }) {
  const user = getCanopyUserFromRequest(request);
  if (!user) return unauthorizedResponse();

  const { id } = await context.params;
  const supabase = createCanopyServiceClient();
  const { data, error } = await supabase
    .from('canopy_flowsheets')
    .select('id, name, created_at, updated_at, flowsheet_json')
    .eq('id', id)
    .eq('user_id', user.id)
    .single();

  if (error) {
    return NextResponse.json({ error: error.message }, { status: error.code === 'PGRST116' ? 404 : 500 });
  }
  return NextResponse.json({ save: data });
}

export async function PUT(request: NextRequest, context: { params: Promise<{ id: string }> }) {
  const user = getCanopyUserFromRequest(request);
  if (!user) return unauthorizedResponse();

  const { id } = await context.params;
  const body = await request.json().catch(() => null) as { name?: string; flowsheetJson?: string } | null;
  const parsed = parseCanopyFlowsheetJson(body?.flowsheetJson ?? '');
  if (!parsed.ok) {
    return NextResponse.json({ error: parsed.error }, { status: 400 });
  }

  const supabase = createCanopyServiceClient();
  const { data, error } = await supabase
    .from('canopy_flowsheets')
    .update({
      name: sanitizeCloudSaveName(body?.name ?? 'Canopy Flowsheet'),
      flowsheet_json: parsed.data,
    })
    .eq('id', id)
    .eq('user_id', user.id)
    .select('id, name, created_at, updated_at')
    .single();

  if (error) {
    return NextResponse.json({ error: error.message }, { status: error.code === 'PGRST116' ? 404 : 500 });
  }
  return NextResponse.json({ save: data });
}

export async function DELETE(request: NextRequest, context: { params: Promise<{ id: string }> }) {
  const user = getCanopyUserFromRequest(request);
  if (!user) return unauthorizedResponse();

  const { id } = await context.params;
  const supabase = createCanopyServiceClient();
  const { error } = await supabase
    .from('canopy_flowsheets')
    .delete()
    .eq('id', id)
    .eq('user_id', user.id);

  if (error) {
    return NextResponse.json({ error: error.message }, { status: 500 });
  }
  return NextResponse.json({ success: true });
}
