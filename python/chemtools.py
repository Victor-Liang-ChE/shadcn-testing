import re
from chempy import Substance
import json
from typing import Dict, List, Tuple, Any, Union, Optional

# Try to import optional dependencies
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

try:
    import pubchempy as pcp
    PUBCHEM_AVAILABLE = True
except ImportError:
    PUBCHEM_AVAILABLE = False

# Function to format chemical formulas with subscripts for display
def format_chemical_formula_components(formula):
    """Convert a chemical formula into a list of components (regular text and subscripts)
    Returns a list of tuples where each tuple is ('text', is_subscript)
    For example, H2O would return [('H', False), ('2', True), ('O', False)]
    """
    # First, handle the case where there's a leading coefficient (like 2H2O)
    leading_coef_match = re.match(r'^(\d+)([A-Za-z(].*)', formula)
    coefficient = ""
    chemical_part = formula
    
    if leading_coef_match:
        coefficient = leading_coef_match.group(1)
        chemical_part = leading_coef_match.group(2)
    
    components = []
    # Add the coefficient if it exists
    if (coefficient):
        components.append((coefficient, False))
    
    # Improved parsing that properly handles two-letter element symbols
    i = 0
    while i < len(chemical_part):
        char = chemical_part[i]
        
        # Case 1: Opening parenthesis
        if (char == '('):
            # Find the closing parenthesis
            paren_depth = 1
            j = i + 1
            while j < len(chemical_part) and paren_depth > 0:
                if chemical_part[j] == '(':
                    paren_depth += 1
                elif chemical_part[j] == ')':
                    paren_depth -= 1
                j += 1
            
            # Get the content inside parentheses
            paren_content = chemical_part[i:j]
            components.append((paren_content, False))
            i = j
            
            # Check if there's a subscript after the parenthesis
            if i < len(chemical_part) and chemical_part[i].isdigit():
                subscript_start = i
                while i < len(chemical_part) and chemical_part[i].isdigit():
                    i += 1
                subscript = chemical_part[subscript_start:i]
                components.append((subscript, True))
        
        # Case 2: Capital letter (start of an element symbol)
        elif char.isupper():
            element = char
            i += 1
            
            # Check for lowercase letters that are part of the same element symbol (e.g., Na, Mg)
            while i < len(chemical_part) and chemical_part[i].islower():
                element += chemical_part[i]
                i += 1
                
            components.append((element, False))
            
            # Check for subscripts (digits after an element)
            if i < len(chemical_part) and chemical_part[i].isdigit():
                subscript_start = i
                while i < len(chemical_part) and chemical_part[i].isdigit():
                    i += 1
                subscript = chemical_part[subscript_start:i]
                components.append((subscript, True))
        
        # Case 3: Any other character (should not normally occur in valid formulas)
        else:
            components.append((char, False))
            i += 1
            
    return components

# Function to parse chemical reaction equation
def parse_reaction_equation(equation):
    # Remove spaces and split by arrow
    parts = equation.replace(" ", "").split("->")
    if len(parts) != 2:
        return None, "Invalid reaction equation. Use format like 'A + B -> C + D'"
    
    reactants_str, products_str = parts
    
    # Split reactants and products by plus sign
    reactants = reactants_str.split("+")
    products = products_str.split("+")
    
    # Parse coefficients and compounds for reactants
    parsed_reactants = []
    for reactant in reactants:
        # Match pattern for coefficient and compound, handling complex formulas with parentheses
        match = re.match(r"^(\d*)([A-Za-z0-9()]+)$", reactant)
        if not match:
            return None, f"Invalid reactant format: {reactant}"
        coef, compound = match.groups()
        coef = int(coef) if coef else 1
        parsed_reactants.append({"compound": compound, "coefficient": coef, 
                               "display": format_chemical_formula_components(reactant)})
    
    # Parse coefficients and compounds for products
    parsed_products = []
    for product in products:
        match = re.match(r"^(\d*)([A-Za-z0-9()]+)$", product)
        if not match:
            return None, f"Invalid product format: {product}"
        coef, compound = match.groups()
        coef = int(coef) if coef else 1
        parsed_products.append({"compound": compound, "coefficient": coef, 
                              "display": format_chemical_formula_components(product)})
    
    return {
        "reactants": parsed_reactants,
        "products": parsed_products
    }, "Success"

# Function to get molar mass for a compound
def get_molar_mass(compound):
    try:
        # Add debugging to see what elements are detected in the compound
        print(f"Parsing compound: {compound}")
        
        # Debug parsed elements using our format_chemical_formula_components function
        parsed_components = format_chemical_formula_components(compound)
        elements = []
        current_element = None
        current_count = 1
        
        # Go through the parsed components and extract element-count pairs
        for text, is_subscript in parsed_components:
            if not is_subscript and text != '(' and text != ')':
                if current_element and current_count:
                    elements.append((current_element, current_count))
                current_element = text
                current_count = 1
            elif is_subscript and current_element:
                try:
                    current_count = int(text)
                except ValueError:
                    print(f"  Warning: Could not parse subscript {text} as integer")
        
        # Don't forget the last element
        if current_element and current_count:
            elements.append((current_element, current_count))
        
        print(f"  Detected elements: {elements}")
        
        # Calculate molar mass using ChemPy
        try:
            substance = Substance.from_formula(compound)
            molar_mass = substance.molar_mass()
            print(f"  Molar mass from ChemPy: {molar_mass:.4f} g/mol")
            
            # If we got here, ChemPy successfully calculated a molar mass
            # Store this value in case subsequent operations fail
            successful_molar_mass = molar_mass
            
            # Verify by calculating element-by-element
            calc_molar_mass = 0
            for element, count in elements:
                try:
                    # Get molar mass for individual element using ChemPy
                    element_substance = Substance.from_formula(element)
                    element_mw = element_substance.molar_mass()
                    calc_molar_mass += element_mw * count
                    print(f"  {element}: {count} × {element_mw:.4f} g/mol = {element_mw * count:.4f} g/mol")
                except Exception as elem_error:
                    print(f"  Warning: Could not get molar mass for {element}: {str(elem_error)}")
            
            print(f"  Calculated molar mass: {calc_molar_mass:.2f} g/mol ({molar_mass:.2f} g/mol)")
            
            return molar_mass
            
        except Exception as chempy_error:
            # Important fix: Check if we already got a valid molar mass value before the error
            if 'successful_molar_mass' in locals() and successful_molar_mass > 0:
                print(f"  ChemPy calculated mass ({successful_molar_mass:.4f} g/mol) but then encountered an error: {str(chempy_error)}")
                print(f"  Using the successfully calculated mass anyway")
                return successful_molar_mass
            
            print(f"  ChemPy error: {str(chempy_error)}")
            
            # Try element-by-element calculation as fallback
            if 'calc_molar_mass' in locals() and calc_molar_mass > 0:
                print(f"  Using our manual element-by-element calculation as fallback: {calc_molar_mass:.4f} g/mol")
                return calc_molar_mass
            
            # If we couldn't calculate the mass, return None
            print(f"  ERROR: Could not calculate molar mass for {compound}")
            return None
        
    except Exception as e:
        print(f"  ERROR getting molar mass for {compound}: {str(e)}")
        return None

# Functions for 3D molecular visualization
def get_smiles_from_name(name):
    """Get SMILES string directly from chemical name using PubChem"""
    if not PUBCHEM_AVAILABLE:
        return None, "PubChem library not available"
        
    try:
        compounds = pcp.get_compounds(name, 'name')
        if not compounds:
            return None, "Chemical not found in database"
        return compounds[0].canonical_smiles, compounds[0].iupac_name
    except Exception as e:
        # Check if it's a network error
        if 'urlopen error' in str(e) or 'getaddrinfo failed' in str(e):
            return None, "PubChem servers appear to be unavailable. Please try again later."
        return None, f"Error: {str(e)}"

def create_molecule_from_smiles(smiles):
    if not RDKIT_AVAILABLE:
        return None, "RDKit library not available"
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv2()
        params.randomSeed = 42
        params.clearConfs = True
        params.useSmallRingTorsions = True
        params.useBasicKnowledge = True
        
        status = AllChem.EmbedMolecule(mol, params=params)
        if (status == -1):
            status = AllChem.EmbedMolecule(mol, useRandomCoords=True, 
                                         randomSeed=42,
                                         clearConfs=True)
            
        if (status == -1):
            raise ValueError("Could not generate 3D coordinates")
            
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception as force_field_error:
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=200)
            except Exception:
                # Continue without optimization
                pass
        
        return mol, None
    except Exception as e:
        return None, str(e)

def get_element_list(mol):
    """Extract a set of unique elements in the molecule"""
    elements = set()
    if mol:
        for atom in mol.GetAtoms():
            elements.add(atom.GetSymbol())
    return sorted(list(elements))

def create_viewer_html(mol):
    if not RDKIT_AVAILABLE or mol is None:
        return ""
        
    pdb_data = Chem.MolToPDBBlock(mol)
    pdb_data = pdb_data.replace('\\', '\\\\').replace('\n', '\\n').replace("'", "\\'")
    
    # Get list of elements in the molecule
    element_list = get_element_list(mol)
    
    # Complete JMol element color mapping (abbreviated for brevity)
    element_colors = {
        'H': '#FFFFFF',   # White
        'He': '#D9FFFF',  # Light cyan
        'Li': '#CC80FF',  # Light purple
        'B': '#FFB5B5',   # Light red
        'C': '#909090',   # Gray
        'N': '#3050F8',   # Blue
        'O': '#FF0D0D',   # Red
        'F': '#90E050',   # Light green
        'P': '#FF8000',   # Orange
        'S': '#FFFF30',   # Yellow
        'Cl': '#1FF01F',  # Bright green
        'Br': '#A62929',  # Brown
        'I': '#940094',   # Purple
        'default': '#1FF01F'  # Bright green for unknown elements
    }   
    
    # Create legend HTML based on actual elements in the molecule
    legend_html = '<div style="display: flex; justify-content: center; align-items: center; flex-wrap: wrap;">\n'
    for element in element_list:
        color = element_colors.get(element, element_colors['default'])
        legend_html += f'<div style="margin: 0 10px;"><span class="elementColor" style="background-color: {color};"></span>{element}</div>\n'
    legend_html += '</div>'
    
    # Use the exact format for JSmol viewer
    jsmol_html = '''
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="UTF-8">
        <style>
            html, body {
                height: 100%;
                margin: 0;
                padding: 0;
                overflow: hidden;
            }
            #loading {
                position: absolute;
                top: 50%;
                left: 50%;
                transform: translate(-50%, -50%);
                text-align: center;
                z-index: 1000;
                background-color: rgba(255,255,255,0.7);
                padding: 20px;
                border-radius: 10px;
            }
            #jsmolDiv {
                width: 100%;
                height: 100%;
                position: relative;
            }
            #elementLegend {
                position: absolute;
                bottom: 0;
                left: 0;
                right: 0;
                padding: 5px;
                background-color: rgba(250,250,250,0.9);
                border-radius: 5px 5px 0 0;
                z-index: 100;
                font-size: 12px;
                font-family: Arial, sans-serif;
                border-top: 1px solid #ccc;
                text-align: center;
            }
            .elementColor {
                display: inline-block;
                width: 12px;
                height: 12px;
                margin-right: 5px;
                border: 1px solid #666;
                vertical-align: middle;
            }
        </style>
        <script src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
    </head>
    <body>
        <div id="loading">Loading molecule viewer...</div>
        <div id="jsmolDiv"></div>
        
        <!-- Element color legend - positioned at bottom -->
        <div id="elementLegend">
            {1}
        </div>
        
        <script>
            var pdbData = '{0}';
            
            // Load JSmol immediately
            document.addEventListener("DOMContentLoaded", function() {{
                loadJSmol();
            }});
            
            function loadJSmol() {{
                // Create a script element to load JSmol
                var script = document.createElement('script');
                script.src = "https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js";
                script.onload = initJSmol;
                document.head.appendChild(script);
            }}
            
            function initJSmol() {{
                var Info = {{ 
                    width: "100%", 
                    height: "100%", 
                    use: "HTML5", 
                    j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s", 
                    addSelectionOptions: false, 
                    debug: false, 
                    color: "white", 
                    disableInitialConsole: true, 
                    disableJ2SLoadMonitor: true, 
                    allowJavaScript: true 
                }};
                
                $("#jsmolDiv").html(Jmol.getAppletHtml("jmolApplet0", Info));
                setTimeout(loadMolecule, 500);
            }}
            
            function loadMolecule() {{
                try {{
                    Jmol.script(jmolApplet0, "zap");
                    var loadCmd = 'load DATA "pdb"\\n' + pdbData + '\\nEND "pdb"';
                    Jmol.script(jmolApplet0, loadCmd);

                    // Styling for the molecule
                    var styling = [
                        "select all;",
                        "wireframe 0.15;", 
                        "spacefill 25%;",
                        "select hydrogen; color white;",
                        "select carbon; color [144,144,144];", 
                        "select nitrogen; color [48,80,248];", 
                        "select oxygen; color [255,13,13];",   
                        "select sulfur; color [255,255,48];", 
                        "set showHydrogens TRUE;",
                        "zoom 80;", 
                        "center;",
                        "frank off;"
                    ].join('\\n');
                    
                    Jmol.script(jmolApplet0, styling);
                    $("#loading").hide();
                }} catch (e) {{
                    console.error("Error in loadMolecule:", e);
                    $("#loading").html("<p style='color:red'>Error loading molecule: " + e.message + "</p>");
                }}
            }}
        </script>
    </body>
    </html>
    '''.format(pdb_data, legend_html)
    
    return jsmol_html

def calculate_stoichiometry(reaction_data, input_data, conversion_percentage=100):
    # Extract reactants and products data
    reactants = reaction_data["reactants"]
    products = reaction_data["products"]
    
    # Get input amounts
    reactant_amounts = {}
    for reactant in reactants:
        compound = reactant["compound"]
        if compound in input_data:
            moles = input_data[compound].get("moles")
            grams = input_data[compound].get("grams")
            molar_mass = input_data[compound].get("molar_mass")
            
            # Convert grams to moles if needed
            if moles is None and grams is not None and molar_mass is not None:
                moles = grams / molar_mass
            
            if moles is not None:
                reactant_amounts[compound] = {
                    "moles": moles,
                    "coefficient": reactant["coefficient"],
                    "moles_per_coefficient": moles / reactant["coefficient"]
                }
    
    # If fewer than one reactant has amounts, return an error
    if len(reactant_amounts) < 1:
        return None, "Please provide amount data for at least one reactant"
    
    # Find limiting reactant
    limiting_reactant = min(reactant_amounts.items(), key=lambda x: x[1]["moles_per_coefficient"])
    limiting_compound = limiting_reactant[0]
    limiting_data = limiting_reactant[1]
    
    # Apply conversion percentage to limit how much of the limiting reactant is consumed
    conversion_factor = conversion_percentage / 100.0
    
    # Calculate product amounts
    results = {
        "limiting_reactant": limiting_compound,
        "reactants": {},
        "products": {},
        "conversion_percentage": conversion_percentage
    }
    
    # Calculate amounts for all reactants
    for reactant in reactants:
        compound = reactant["compound"]
        coef = reactant["coefficient"]
        molar_mass = input_data.get(compound, {}).get("molar_mass", None)
        
        if compound in reactant_amounts:
            moles = reactant_amounts[compound]["moles"]
            # Apply conversion factor to used moles
            used_moles = limiting_data["moles_per_coefficient"] * coef * conversion_factor
            excess_moles = moles - used_moles if moles > used_moles else 0
            
            results["reactants"][compound] = {
                "initial_moles": moles,
                "used_moles": used_moles,
                "excess_moles": excess_moles,
            }
            
            if molar_mass:
                results["reactants"][compound].update({
                    "initial_grams": moles * molar_mass,
                    "used_grams": used_moles * molar_mass,
                    "excess_grams": excess_moles * molar_mass
                })
    
    # Calculate product amounts
    for product in products:
        compound = product["compound"]
        coef = product["coefficient"]
        molar_mass = input_data.get(compound, {}).get("molar_mass", None)
        
        # Apply conversion factor to produced moles
        produced_moles = limiting_data["moles_per_coefficient"] * coef * conversion_factor
        
        results["products"][compound] = {
            "produced_moles": produced_moles
        }
        
        if molar_mass:
            results["products"][compound]["produced_grams"] = produced_moles * molar_mass
    
    return results, "Success"

# API handlers for Next.js
def handle_parse_reaction(equation):
    reaction_data, message = parse_reaction_equation(equation)
    
    if reaction_data is None:
        return {"success": False, "message": message}
    
    # Get molar masses for all compounds
    molar_masses = {}
    
    for reactant in reaction_data["reactants"]:
        compound = reactant["compound"]
        molar_masses[compound] = get_molar_mass(compound)
    
    for product in reaction_data["products"]:
        compound = product["compound"]
        molar_masses[product["compound"]] = get_molar_mass(product["compound"])
    
    # Return both the reaction data and molar masses
    return {
        "success": True,
        "reaction_data": reaction_data,
        "molar_masses": {k: float(v) if v is not None else None for k, v in molar_masses.items()}
    }

def handle_molecule_visualization(name):
    smiles, name_or_error = get_smiles_from_name(name)
    if smiles is None:
        return {"success": False, "message": name_or_error}
    
    mol, error = create_molecule_from_smiles(smiles)
    if error:
        return {"success": False, "message": f"Error creating 3D model: {error}"}
    
    viewer_html = create_viewer_html(mol)
    return {
        "success": True,
        "chemical_name": name_or_error or name,
        "viewer_html": viewer_html
    }

def handle_stoichiometry_calculation(reaction_data, input_data, conversion_percentage=100):
    results, message = calculate_stoichiometry(reaction_data, input_data, conversion_percentage)
    
    if results is None:
        return {"success": False, "message": message}
    
    return {"success": True, "results": results}

# Add to the end of the file to support command-line interface

# Main entry point for command-line operations
if __name__ == "__main__":
    import sys
    import json

    # Command line arguments:
    # sys.argv[1]: operation (parse-reaction, visualize-molecule)
    # sys.argv[2:]: arguments for the operation

    if len(sys.argv) < 2:
        print(json.dumps({"success": False, "message": "No operation specified"}))
        sys.exit(1)

    operation = sys.argv[1]

    if operation == "parse-reaction":
        if len(sys.argv) < 3:
            print(json.dumps({"success": False, "message": "No reaction equation provided"}))
            sys.exit(1)

        equation = sys.argv[2]
        result = handle_parse_reaction(equation)
        print(json.dumps(result))

    elif operation == "visualize-molecule":
        if len(sys.argv) < 3:
            print(json.dumps({"success": False, "message": "No chemical name provided"}))
            sys.exit(1)

        name = sys.argv[2]
        result = handle_molecule_visualization(name)
        print(json.dumps(result))

    else:
        print(json.dumps({"success": False, "message": f"Unknown operation: {operation}"}))
        sys.exit(1)
