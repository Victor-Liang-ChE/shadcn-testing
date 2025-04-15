import os
from flask import Flask, request, jsonify
# Enable CORS if requests are coming directly from browser in some scenarios
# from flask_cors import CORS
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

# Import McCabe-Thiele calculator function
try:
    # Assuming mccabe_thiele_calculator.py is in the same directory
    from mccabe_thiele_calculator import calculate_vle
    print("DEBUG: Successfully imported mccabe_thiele_calculator", file=sys.stderr)
except ImportError as e:
    print(f"DEBUG ERROR: Failed to import mccabe_thiele_calculator: {e}", file=sys.stderr)
    def calculate_vle(comp1, comp2, temperature, pressure):
        return {"error": "McCabe-Thiele calculator module failed to load."}

# Import chemistry tools functions
try:
    # Assuming chemtools.py is in the same directory
    from chemtools import (
        handle_parse_reaction,
        handle_molecule_visualization
        # get_molar_mass # Import if needed directly by routes, otherwise handled within chemtools
    )
    print("DEBUG: Successfully imported chemtools functions", file=sys.stderr)
except ImportError as e:
    print(f"DEBUG ERROR: Failed to import chemtools: {e}", file=sys.stderr)
    def handle_parse_reaction(equation):
        return {"error": "ChemTools module failed to load."}
    def handle_molecule_visualization(name):
        return {"error": "ChemTools module failed to load."}


app = Flask(__name__)
# If needed, enable CORS for specific origins or all:
# CORS(app) # Allows all origins
# Or more specific:
# CORS(app, resources={r"/*": {"origins": "YOUR_FRONTEND_URL"}}) # Example for all routes


# McCabe-Thiele route
@app.route('/calculate', methods=['POST'])
def calculate():
    print("DEBUG: /calculate endpoint hit", file=sys.stderr)
    if not request.is_json:
        print("DEBUG ERROR: Request is not JSON", file=sys.stderr)
        return jsonify({"error": "Request must be JSON"}), 400
    data = request.json
    comp1 = data.get('comp1')
    comp2 = data.get('comp2')
    temperature = data.get('temperature')
    pressure = data.get('pressure')

    if not comp1 or not comp2:
         print("DEBUG ERROR: Missing component names in /calculate", file=sys.stderr)
         return jsonify({"error": "Component names (comp1, comp2) are required"}), 400

    try:
        # Call your existing calculation function
        result = calculate_vle(comp1, comp2, temperature, pressure)
        print("DEBUG: /calculate processed successfully", file=sys.stderr)
        if 'error' in result:
             print(f"DEBUG WARN: /calculate processor returned an error: {result['error']}", file=sys.stderr)
             return jsonify(result), 500
        return jsonify(result)
    except Exception as e:
        print(f"DEBUG ERROR: Exception during /calculate processing: {e}", file=sys.stderr)
        return jsonify({"error": f"An internal server error occurred: {e}"}), 500

# Chemistry routes
@app.route('/chemistry/parse-reaction', methods=['POST'])
def parse_reaction():
    print("DEBUG: /chemistry/parse-reaction endpoint hit", file=sys.stderr)
    if not request.is_json:
        print("DEBUG ERROR: Request is not JSON", file=sys.stderr)
        return jsonify({"error": "Request must be JSON"}), 400
    data = request.json
    equation = data.get('equation')

    if not equation:
        print("DEBUG ERROR: No equation provided in /chemistry/parse-reaction", file=sys.stderr)
        return jsonify({"success": False, "message": "No equation provided"}), 400

    try:
        result = handle_parse_reaction(equation)
        print("DEBUG: /chemistry/parse-reaction processed successfully", file=sys.stderr)
        return jsonify(result)
    except Exception as e:
        print(f"DEBUG ERROR: Exception during /chemistry/parse-reaction processing: {e}", file=sys.stderr)
        return jsonify({"success": False, "message": f"An internal server error occurred: {e}"}), 500


@app.route('/chemistry/visualize-molecule', methods=['POST'])
def visualize_molecule():
    print("DEBUG: /chemistry/visualize-molecule endpoint hit", file=sys.stderr)
    if not request.is_json:
        print("DEBUG ERROR: Request is not JSON", file=sys.stderr)
        return jsonify({"error": "Request must be JSON"}), 400
    data = request.json
    chemical_name = data.get('name')

    if not chemical_name:
        print("DEBUG ERROR: No chemical name provided in /chemistry/visualize-molecule", file=sys.stderr)
        return jsonify({"success": False, "message": "No chemical name provided"}), 400

    try:
        result = handle_molecule_visualization(chemical_name)
        print("DEBUG: /chemistry/visualize-molecule processed successfully", file=sys.stderr)
        return jsonify(result)
    except Exception as e:
        print(f"DEBUG ERROR: Exception during /chemistry/visualize-molecule processing: {e}", file=sys.stderr)
        return jsonify({"success": False, "message": f"An internal server error occurred: {e}"}), 500


# Japanese lyrics route (existing)
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
             print(f"DEBUG WARN: Processor returned an error: {result['error']}", file=sys.stderr)
             # Decide if internal processing errors should be 500 or maybe 4xx
             return jsonify(result), 500 # Internal Server Error seems appropriate
        return jsonify(result)
    except Exception as e:
        print(f"DEBUG ERROR: Exception during lyrics processing: {e}", file=sys.stderr)
        return jsonify({"error": f"An internal server error occurred: {e}"}), 500


# Health check route (existing)
@app.route('/health', methods=['GET'])
def health():
    return jsonify({"status": "OK"}), 200

# Main execution block (existing)
if __name__ == "__main__":
    # Default port 8080 if not specified by environment
    port = int(os.environ.get('PORT', 8080))
    print(f"DEBUG: Starting Flask server on host 0.0.0.0, port {port}", file=sys.stderr)
    # Use debug=False for production/deployment
    # Use debug=True for local development if you want auto-reloading and more detailed errors
    app.run(debug=False, host='0.0.0.0', port=port)