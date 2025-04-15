import { NextResponse } from 'next/server';

// Define the Cloud Run URL
const CLOUD_RUN_LYRICS_URL = 'https://nextjsbackend-23789472506.us-west1.run.app/process-lyrics';

// Determine the backend URL
// 1. Check environment variable (recommended for flexibility)
// 2. Use the Cloud Run URL if env var is not set
// 3. Default to localhost as a last resort (for local Python testing)
const PYTHON_BACKEND_URL = process.env.PYTHON_LYRICS_BACKEND_URL || CLOUD_RUN_LYRICS_URL || 'http://127.0.0.1:8080/process-lyrics';


export async function POST(request: Request) {
  console.log(`DEBUG: /api/process-lyrics route called. Target backend: ${PYTHON_BACKEND_URL}`);

  try {
    const body = await request.json();
    const { text } = body;

    if (text === undefined || text === null) {
      console.error("DEBUG ERROR: 'text' field missing in request body");
      return NextResponse.json({ error: "'text' field is required" }, { status: 400 });
    }

    console.log(`DEBUG: Forwarding request to Python backend. Text length: ${text.length}`);

    // Forward the request to the Python backend
    const response = await fetch(PYTHON_BACKEND_URL, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        // Add any other headers if required by your backend
      },
      body: JSON.stringify({ text }),
      // Add cache: 'no-store' if you don't want Next.js to cache the API route response
      cache: 'no-store',
    });

    // Check if the backend response is successful
    if (!response.ok) {
      const errorBody = await response.text();
      console.error(`DEBUG ERROR: Python backend responded with status ${response.status}. Body: ${errorBody}`);
      // Try to parse error JSON from backend if possible
      let errorJson = { error: `Backend failed with status ${response.status}` };
      try {
        errorJson = JSON.parse(errorBody);
      } catch (e) { /* Ignore parsing error, use default message */ }

      return NextResponse.json(errorJson, { status: response.status });
    }

    // Parse the JSON response from the backend
    const result = await response.json();
    console.log("DEBUG: Successfully received response from Python backend.");

    // Return the result from the backend to the frontend client
    return NextResponse.json(result);

  } catch (error: any) {
    console.error('DEBUG ERROR: Error in /api/process-lyrics route:', error);

    // Handle potential fetch errors (e.g., backend not reachable)
    if (error.cause?.code === 'ECONNREFUSED') {
         console.error(`DEBUG ERROR: Connection refused. Is the Python backend running at ${PYTHON_BACKEND_URL}?`);
         return NextResponse.json(
             { error: `Could not connect to the backend service at ${PYTHON_BACKEND_URL}. Please ensure it's running.` },
             { status: 503 } // Service Unavailable
         );
    }

    return NextResponse.json(
      { error: error.message || 'Failed to process lyrics request' },
      { status: 500 } // Internal Server Error
    );
  }
}
