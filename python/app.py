import os
from flask import Flask, request, jsonify
import json
import sys # Import sys for stderr

# Import lyrics processor functions
try:
    from japanese_lyrics_processor import process_lyrics
    print("DEBUG: Successfully imported japanese_lyrics_processor", file=sys.stderr)
except ImportError as e:
    print(f"DEBUG ERROR: Failed to import japanese_lyrics_processor: {e}", file=sys.stderr)
    # Define a dummy function if import fails, so the app can still start (maybe)
    def process_lyrics(text):
        return {"error": "Lyrics processor module failed to load."}


# Import chemistry tools functions if still needed
# from chemtools import (
#     parse_reaction_equation,
#     handle_parse_reaction,
#     handle_molecule_visualization,
#     get_molar_mass
# )

app = Flask(__name__)

# Keep existing /calculate route if needed
# @app.route('/calculate', methods=['POST'])
# def calculate():
#     # ... existing calculate code ...

# Keep existing chemistry endpoints if needed
# @app.route('/chemistry/parse-reaction', methods=['POST'])
# def parse_reaction():
#     # ... existing parse_reaction code ...

# @app.route('/chemistry/visualize-molecule', methods=['POST'])
# def visualize_molecule():
#     # ... existing visualize_molecule code ...

# New route for Japanese lyrics processing
@app.route('/process-lyrics', methods=['POST'])
def handle_process_lyrics():
    print("DEBUG: /process-lyrics endpoint hit", file=sys.stderr)
    if not request.is_json:
        print("DEBUG ERROR: Request is not JSON", file=sys.stderr)
        return jsonify({"error": "Request must be JSON"}), 400

    data = request.get_json()
    text = data.get('text')

    if text is None: # Check for None explicitly, empty string is allowed
        print("DEBUG ERROR: 'text' field missing in JSON payload", file=sys.stderr)
        return jsonify({"error": "'text' field is required"}), 400

    print(f"DEBUG: Received text length: {len(text)}", file=sys.stderr)
    try:
        result = process_lyrics(text)
        print("DEBUG: Lyrics processed successfully", file=sys.stderr)
        # Check if the processor returned an error internally
        if 'error' in result:
             return jsonify(result), 500
        return jsonify(result)
    except Exception as e:
        print(f"DEBUG ERROR: Exception during lyrics processing: {e}", file=sys.stderr)
        return jsonify({"error": f"An internal error occurred: {e}"}), 500


@app.route('/health', methods=['GET'])
def health():
    return jsonify({"status": "OK"}), 200

if __name__ == "__main__":
    # Default port 8080 if not specified by environment
    port = int(os.environ.get('PORT', 8080))
    print(f"DEBUG: Starting Flask server on host 0.0.0.0, port {port}", file=sys.stderr)
    # Use debug=False for production/deployment
    # Use debug=True for local development if you want auto-reloading and more detailed errors
    app.run(debug=False, host='0.0.0.0', port=port)