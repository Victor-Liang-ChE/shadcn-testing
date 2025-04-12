import { NextResponse } from 'next/server';

// Cloud Run URL - store this in an environment variable
const CLOUD_RUN_URL = process.env.CLOUD_RUN_URL || 'https://nextjsbackend-23789472506.us-west1.run.app/calculate';

export async function POST(request: Request) {
  console.log("DEBUG: McCabe-Thiele API route called");
  try {
    const body = await request.json();
    console.log("DEBUG: Request body received:", body);
    const { comp1, comp2, temperature, pressure } = body;
    console.log(`DEBUG: Extracted parameters - comp1: ${comp1}, comp2: ${comp2}, temperature: ${temperature}, pressure: ${pressure}`);
    
    // Validate inputs
    if (!comp1 || !comp2) {
      console.error("DEBUG ERROR: Missing component names");
      return NextResponse.json(
        { error: 'Component names are required' },
        { status: 400 }
      );
    }

    console.log("DEBUG: Forwarding request to Cloud Run");
    
    // Forward the request to Cloud Run
    const response = await fetch(CLOUD_RUN_URL, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
      },
      body: JSON.stringify({ comp1, comp2, temperature, pressure })
    });

    if (!response.ok) {
      console.error(`DEBUG ERROR: Cloud Run API responded with status: ${response.status}`);
      return NextResponse.json(
        { error: `Backend calculation failed with status: ${response.status}` },
        { status: response.status }
      );
    }

    const processResult = await response.json();
    console.log("DEBUG: Successfully processed McCabe-Thiele data");
    return NextResponse.json(processResult);
    
  } catch (error: any) {
    console.error('DEBUG ERROR: API route error:', error);
    return NextResponse.json(
      { error: error.message || 'Failed to process request' },
      { status: 500 }
    );
  }
}
