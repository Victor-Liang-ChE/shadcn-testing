import sys
import json
import numpy as np
import os
print("DEBUG: Starting McCabe-Thiele calculator", file=sys.stderr)
print(f"DEBUG: Python version: {sys.version}", file=sys.stderr)
print(f"DEBUG: NumPy version: {np.__version__}", file=sys.stderr)
print(f"DEBUG: Current working directory: {os.getcwd()}", file=sys.stderr)

# Import necessary thermodynamic libraries
try:
    print("DEBUG: Attempting to import thermodynamic libraries", file=sys.stderr)
    from thermo.chemical import Chemical
    from thermo.chemical_package import ChemicalConstantsPackage
    from thermo.eos_mix import PRMIX, VDWMIX, SRKMIX, RKMIX
    from thermo.phases import CEOSGas, GibbsExcessLiquid
    from thermo.unifac import UNIFAC, DOUFIP2016, DOUFSG
    from thermo import FlashVL
    print("DEBUG: Successfully imported thermodynamic libraries", file=sys.stderr)
except ImportError as e:
    print(f"DEBUG ERROR: Failed to import thermodynamic libraries: {str(e)}", file=sys.stderr)
    print("Error: Required thermodynamic packages not installed", file=sys.stderr)
    sys.exit(1)

def xy(comp1, comp2, T=None, P=None, values=False, show=True):
    """Calculate vapor-liquid equilibrium data using rigorous thermodynamic models"""
    print(f"DEBUG: xy function called with comp1={comp1}, comp2={comp2}, T={T}, P={P}", file=sys.stderr)
    try:
        # Load constants and properties
        print(f"DEBUG: Loading constants and properties for {comp1} and {comp2}", file=sys.stderr)
        constants, properties = ChemicalConstantsPackage.from_IDs([comp1, comp2])
        print("DEBUG: Successfully loaded constants and properties", file=sys.stderr)
        
        if P is not None:
            T = 273.15 if T is None else T  # K
            P = P*1e5  # bar to Pa
            Pgiven = True
            print(f"DEBUG: Working with fixed pressure P={P} Pa, T={T} K", file=sys.stderr)
        elif T is not None:
            P = 1e5  # Pa
            Pgiven = False
            print(f"DEBUG: Working with fixed temperature T={T} K, P={P} Pa", file=sys.stderr)
        else:
            T = 298.15  # Default temperature
            P = 1e5    # Default pressure
            Pgiven = False
            print(f"DEBUG: Using default conditions T={T} K, P={P} Pa", file=sys.stderr)
            
        zs = [.5, .5]  # initial mole fraction of comp1 and comp2
        
        # Use Peng-Robinson for the vapor phase
        print("DEBUG: Setting up gas phase model with Peng-Robinson", file=sys.stderr)
        eos_kwargs = {'Pcs': constants.Pcs, 'Tcs': constants.Tcs, 'omegas': constants.omegas}
        gas = CEOSGas(PRMIX, HeatCapacityGases=properties.HeatCapacityGases, eos_kwargs=eos_kwargs)
        print("DEBUG: Gas phase model setup complete", file=sys.stderr)
        
        # Configure the activity model
        print("DEBUG: Setting up UNIFAC activity model", file=sys.stderr)
        GE = UNIFAC.from_subgroups(
            chemgroups=constants.UNIFAC_Dortmund_groups,
            version=1,
            T=T,
            xs=zs,
            interaction_data=DOUFIP2016, 
            subgroups=DOUFSG
        )
        print("DEBUG: UNIFAC activity model setup complete", file=sys.stderr)
        
        # Configure the liquid model with activity coefficients
        print("DEBUG: Setting up liquid phase model", file=sys.stderr)
        liquid = GibbsExcessLiquid(
            VaporPressures=properties.VaporPressures,
            HeatCapacityGases=properties.HeatCapacityGases,
            VolumeLiquids=properties.VolumeLiquids,
            GibbsExcessModel=GE,
            equilibrium_basis='Psat', 
            caloric_basis='Psat',
            T=T, 
            P=P, 
            zs=zs
        )
        print("DEBUG: Liquid phase model setup complete", file=sys.stderr)
        
        # Create flash calculator with the required parameters
        print("DEBUG: Creating flash calculator", file=sys.stderr)
        flasher = FlashVL(constants, properties, liquid=liquid, gas=gas)
        print("DEBUG: Flash calculator created successfully", file=sys.stderr)
        
        # Use the built-in plot_xy method to get VLE data
        print("DEBUG: Calculating VLE data using plot_xy method", file=sys.stderr)
        if Pgiven:
            print(f"DEBUG: Generating xy diagram at constant P={P/1e5:.2f} bar", file=sys.stderr)
            z1, z2, x1_bubble, y1_bubble = flasher.plot_xy(P=P, pts=100, values=True, show=False)
        else:
            print(f"DEBUG: Generating xy diagram at constant T={T:.2f} K", file=sys.stderr)
            z1, z2, x1_bubble, y1_bubble = flasher.plot_xy(T=T, pts=100, values=True, show=False)
        
        # Print a sample of points for debugging
        print("DEBUG: VLE calculation successful. Sample of x,y points:", file=sys.stderr)
        for i in range(0, len(x1_bubble), 10):
            print(f"DEBUG: Point {i}: x={x1_bubble[i]:.3f}, y={y1_bubble[i]:.3f}", file=sys.stderr)
        
        if values:
            return x1_bubble, y1_bubble
        
        return x1_bubble, y1_bubble
        
    except Exception as e:
        print(f"DEBUG ERROR: VLE calculation error: {str(e)}", file=sys.stderr)
        raise

def calculate_vle(comp1, comp2, temperature=None, pressure=None):
    """Calculate VLE data for the given components and conditions"""
    print(f"DEBUG: calculate_vle called with comp1={comp1}, comp2={comp2}, temperature={temperature}, pressure={pressure}", file=sys.stderr)
    try:
        # Convert inputs to appropriate types
        if temperature is not None and temperature != 'null':
            T = float(temperature) + 273.15  # Convert Celsius to Kelvin
            P = None
            print(f"DEBUG: Using specified temperature T={T} K ({temperature}°C)", file=sys.stderr)
        elif pressure is not None and pressure != 'null':
            P = float(pressure)  # Keep as bar for the xy function
            T = None
            print(f"DEBUG: Using specified pressure P={P} bar", file=sys.stderr)
        else:
            T = 298.15  # Default temperature in K
            P = None
            print(f"DEBUG: Using default temperature T={T} K (25°C)", file=sys.stderr)
        
        # Get VLE data with robust thermodynamic models
        print("DEBUG: Calling xy function to get VLE data", file=sys.stderr)
        x_values, y_values = xy(comp1, comp2, T=T, P=P, values=True, show=False)
        print(f"DEBUG: xy function returned {len(x_values)} data points", file=sys.stderr)
        
        # Convert to list for JSON serialization
        x_values_list = x_values.tolist() if hasattr(x_values, 'tolist') else list(x_values)
        y_values_list = y_values.tolist() if hasattr(y_values, 'tolist') else list(y_values)
        
        # Calculate polynomial fit coefficients
        print("DEBUG: Calculating polynomial fit coefficients", file=sys.stderr)
        z = np.polyfit(x_values, y_values, 20)
        z_list = z.tolist()
        print(f"DEBUG: Polynomial fit complete", file=sys.stderr)
        
        # Calculate volatility information
        print(f"DEBUG: Getting volatility information for {comp1} and {comp2}", file=sys.stderr)
        volatility_info = get_component_volatility(comp1, comp2)
        print(f"DEBUG: Volatility info: {volatility_info}", file=sys.stderr)
        
        # Calculate average volatility
        avg_volatility = calculate_average_volatility(x_values, y_values)
        print(f"DEBUG: Average relative volatility: {avg_volatility}", file=sys.stderr)
        
        result = {
            'x_values': x_values_list,
            'y_values': y_values_list,
            'poly_coeffs': z_list,
            'volatility': volatility_info,
            'temperature': T,
            'pressure': P,
            'comp1': comp1,
            'comp2': comp2,
            'avg_volatility': avg_volatility
        }
        print("DEBUG: VLE calculation successful", file=sys.stderr)
        return result
    
    except Exception as e:
        print(f"DEBUG ERROR: Error in calculate_vle: {str(e)}", file=sys.stderr)
        return {
            'error': str(e),
            'x_values': [],
            'y_values': [],
            'poly_coeffs': []
        }

def get_component_volatility(comp1, comp2):
    """Get boiling points and determine which component is more volatile"""
    print(f"DEBUG: get_component_volatility called for {comp1} and {comp2}", file=sys.stderr)
    boiling_points = {
        # Common solvents and chemicals with boiling points in Kelvin
        "methanol": 337.8,
        "ethanol": 351.4,
        "water": 373.15,
        "acetone": 329.2,
        "benzene": 353.2,
        "toluene": 383.8,
        "chloroform": 334.0,
        "hexane": 342.0,
        "heptane": 371.6,
        "octane": 398.8,
        "propanol": 370.4,
        "butanol": 390.9,
        "acetic acid": 391.2,
        "acetonitrile": 354.8,
        "carbon tetrachloride": 349.9,
        "diethyl ether": 307.6,
        "dmf": 426.0,
        "dmso": 462.0,
        "ethyl acetate": 350.3,
        "isopropanol": 355.4,
        "methyl ethyl ketone": 352.8,
        "pentane": 309.2,
    }
    
    comp1_lower = comp1.lower()
    comp2_lower = comp2.lower()
    
    bp1 = boiling_points.get(comp1_lower, None)
    bp2 = boiling_points.get(comp2_lower, None)
    
    print(f"DEBUG: Boiling point of {comp1}: {bp1} K, {comp2}: {bp2} K", file=sys.stderr)
    
    if bp1 is None or bp2 is None:
        print(f"DEBUG: Could not find boiling point data for {comp1} or {comp2}", file=sys.stderr)
        return {
            'message': f"Could not find boiling point data for {comp1} or {comp2}",
            'more_volatile': None,
            'bp1': bp1,
            'bp2': bp2
        }
    
    more_volatile = comp1 if bp1 < bp2 else comp2
    less_volatile = comp2 if bp1 < bp2 else comp1
    
    print(f"DEBUG: {more_volatile} is more volatile (lower boiling point)", file=sys.stderr)
    
    return {
        'message': f"{more_volatile} is more volatile (lower boiling point)",
        'more_volatile': more_volatile,
        'less_volatile': less_volatile,
        'bp1': bp1,
        'bp2': bp2
    }

def calculate_average_volatility(x_values, y_values):
    """Calculate average relative volatility from equilibrium data"""
    print("DEBUG: Calculating average volatility from VLE data", file=sys.stderr)
    volatilities = []
    for x, y in zip(x_values, y_values):
        if 0 < x < 1 and 0 < y < 1:  # Avoid division by zero
            volatility = y*(1-x)/(x*(1-y))
            volatilities.append(volatility)
    
    if not volatilities:
        print("DEBUG: No valid points for volatility calculation", file=sys.stderr)
        return None
    
    avg = sum(volatilities) / len(volatilities)
    print(f"DEBUG: Average volatility calculated: {avg} from {len(volatilities)} points", file=sys.stderr)
    return avg

if __name__ == "__main__":
    print("DEBUG: Script executed as main", file=sys.stderr)
    # Get arguments from command line
    if len(sys.argv) < 5:
        print(f"DEBUG: Invalid number of arguments: {len(sys.argv)}", file=sys.stderr)
        print(f"DEBUG: Arguments received: {sys.argv}", file=sys.stderr)
        result = {'error': 'Invalid number of arguments'}
    else:
        comp1 = sys.argv[1]
        comp2 = sys.argv[2]
        temperature = None if sys.argv[3] == 'null' else float(sys.argv[3])
        pressure = None if sys.argv[4] == 'null' else float(sys.argv[4])
        
        print(f"DEBUG: Processing with args: comp1={comp1}, comp2={comp2}, T={temperature}, P={pressure}", file=sys.stderr)
        result = calculate_vle(comp1, comp2, temperature, pressure)
    
    print(f"DEBUG: Final result structure: keys={list(result.keys())}", file=sys.stderr)
    print(f"DEBUG: Result contains {len(result.get('x_values', []))} data points", file=sys.stderr)
    
    # Print ONLY the JSON result to stdout, with no other text
    print(json.dumps(result))
    print("DEBUG: Script execution complete", file=sys.stderr)
