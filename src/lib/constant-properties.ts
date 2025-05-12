import { baseColors } from './property-equations'; // Added import

export interface ConstantPropertyDefinition {
  displayName: string;
  jsonKey: string; // Key used in the JSON blob from the database
  targetUnitName: string; // The unit the value is typically stored in or should be treated as base
  color?: string; // Optional: for bar chart consistency across different constants if needed
  symbol?: string; // Optional
  availableUnits?: Array<{
    unit: string;
    conversionFactorFromBase: number; // Factor to convert from targetUnitName to this display unit
    displayName?: string;
  }>;
}

// Helper to assign colors if needed, or can be done ad-hoc
const constantColors = ['#2ECC71', '#3498DB', '#9B59B6', '#F1C40F', '#E67E22', '#E74C3C', '#1ABC9C', '#34495E'];

export const constantPropertiesConfig: ConstantPropertyDefinition[] = [
  {
    displayName: "Molecular Weight",
    jsonKey: "Molecular weight",
    targetUnitName: "kg/kmol", // Assuming the base value from DB/calculation is in kg/kmol (or g/mol, which is numerically the same)
    availableUnits: [
      { unit: "kg/kmol", conversionFactorFromBase: 1 },
      { unit: "g/mol", conversionFactorFromBase: 1 },
      { unit: "amu", conversionFactorFromBase: 1 }, // amu is numerically same as g/mol
      { unit: "lb/lbmol", conversionFactorFromBase: 1 }, // Numerically same as g/mol or kg/kmol
    ],
    color: constantColors[0 % constantColors.length],
    symbol: "MW"
  },
  {
    displayName: "Critical Temperature",
    jsonKey: "Critical temperature",
    targetUnitName: "K",
    availableUnits: [
      { unit: "°C", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "°F", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "K", conversionFactorFromBase: 1 },
    ],
    color: constantColors[1 % constantColors.length],
    symbol: "T_c"
  },
  {
    displayName: "Critical Pressure",
    jsonKey: "Critical pressure",
    targetUnitName: "Pa",
    availableUnits: [
      { unit: "Pa", conversionFactorFromBase: 1 },
      { unit: "kPa", conversionFactorFromBase: 1e-3 },
      { unit: "bar", conversionFactorFromBase: 1e-5 },
      { unit: "atm", conversionFactorFromBase: 1 / 101325 },
    ],
    symbol: "P_c",
    color: constantColors[2 % constantColors.length]
  },
  {
    displayName: "Critical Volume",
    jsonKey: "Critical volume",
    targetUnitName: "m³/kmol",
    availableUnits: [
      { unit: "m³/kmol", conversionFactorFromBase: 1 },
      { unit: "L/mol", conversionFactorFromBase: 1 }, // m³/kmol is numerically equal to L/mol
      { unit: "cm³/mol", conversionFactorFromBase: 1000 },
    ],
    symbol: "V_c",
    color: constantColors[3 % constantColors.length]
  },
  {
    displayName: "Critical Compressibility Factor",
    jsonKey: "Critical compressibility factor",
    targetUnitName: "-", // Dimensionless
    availableUnits: [{ unit: "-", conversionFactorFromBase: 1 }],
    symbol: "Z_c",
    color: constantColors[4 % constantColors.length]
  },
  {
    displayName: "Normal Boiling Point",
    jsonKey: "Normal boiling point",
    targetUnitName: "K",
    availableUnits: [
      { unit: "°C", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "°F", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "K", conversionFactorFromBase: 1 },
    ],
    color: constantColors[5 % constantColors.length],
    symbol: "T_b"
  },
  {
    displayName: "Triple Point Temperature",
    jsonKey: "Triple point temperature",
    targetUnitName: "K",
    availableUnits: [
      { unit: "°C", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "°F", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "K", conversionFactorFromBase: 1 },
    ],
    color: constantColors[6 % constantColors.length],
    symbol: "T_tp"
  },
  {
    displayName: "Triple Point Pressure",
    jsonKey: "Triple point pressure",
    targetUnitName: "Pa",
    availableUnits: [
      { unit: "Pa", conversionFactorFromBase: 1 },
      { unit: "kPa", conversionFactorFromBase: 1e-3 },
      { unit: "bar", conversionFactorFromBase: 1e-5 },
    ],
    symbol: "P_tp",
    color: constantColors[7 % constantColors.length]
  },
  {
    displayName: "Flash Point",
    jsonKey: "Flash point",
    targetUnitName: "K",
    availableUnits: [
      { unit: "°C", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "°F", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "K", conversionFactorFromBase: 1 },
    ],
    color: constantColors[8 % constantColors.length],
    symbol: "T_f"
  },
  {
    displayName: "Autoignition Temperature",
    jsonKey: "Autoignition temperature",
    targetUnitName: "K",
    availableUnits: [
      { unit: "°C", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "°F", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "K", conversionFactorFromBase: 1 },
    ],
    color: constantColors[9 % constantColors.length],
    symbol: "T_ai"
  },
  {
    displayName: "Melting Point",
    jsonKey: "Melting point",
    targetUnitName: "K",
    availableUnits: [
      { unit: "°C", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "°F", conversionFactorFromBase: 1 }, // Actual conversion handled in plotting logic
      { unit: "K", conversionFactorFromBase: 1 },
    ],
    color: constantColors[10 % constantColors.length],
    symbol: "T_m"
  },
  {
    displayName: "Critical Density",
    jsonKey: "Critical density",
    targetUnitName: "kg/m³",
    availableUnits: [
      { unit: "kg/m³", conversionFactorFromBase: 1 },
      { unit: "g/cm³", conversionFactorFromBase: 1e3 },
      { unit: "lb/ft³", conversionFactorFromBase: 1 / 16.0185 },
    ],
    symbol: "ρ_c",
    color: constantColors[11 % constantColors.length]
  },
  {
    displayName: "Ideal Gas Constant",
    jsonKey: "Ideal gas constant",
    targetUnitName: "J/(mol·K)",
    availableUnits: [
      { unit: "J/(mol·K)", conversionFactorFromBase: 1 },
      { unit: "kJ/(mol·K)", conversionFactorFromBase: 1e-3 },
      { unit: "cal/(mol·K)", conversionFactorFromBase: 1 / 4184 },
      { unit: "Btu/(lb·°R)", conversionFactorFromBase: 1 / 778.169 },
    ],
    symbol: "R",
    color: constantColors[12 % constantColors.length]
  },
  {
    displayName: "Heat of Vaporization",
    jsonKey: "Heat of vaporization",
    targetUnitName: "J/kmol",
    availableUnits: [
      { unit: "J/kmol", conversionFactorFromBase: 1 },
      { unit: "kJ/kmol", conversionFactorFromBase: 1e-3 },
      { unit: "J/mol", conversionFactorFromBase: 1e-3 },
      { unit: "kJ/mol", conversionFactorFromBase: 1e-6 },
    ],
    symbol: "ΔH_vap",
    color: constantColors[0 % constantColors.length]
  },
  {
    displayName: "Heat of Fusion (Melting Point)",
    jsonKey: "Heat of fusion at melting point",
    targetUnitName: "J/kmol",
    availableUnits: [
      { unit: "J/kmol", conversionFactorFromBase: 1 },
      { unit: "kJ/kmol", conversionFactorFromBase: 1e-3 },
      { unit: "J/mol", conversionFactorFromBase: 1e-3 },
      { unit: "kJ/mol", conversionFactorFromBase: 1e-6 },
    ],
    symbol: "ΔH_fus",
    color: constantColors[1 % constantColors.length]
  },
  {
    displayName: "Std Net Heat of Combustion (LHV)",
    jsonKey: "Standard net heat of combustion LHV",
    targetUnitName: "J/kmol",
    availableUnits: [
      { unit: "J/kmol", conversionFactorFromBase: 1 },
      { unit: "kJ/kmol", conversionFactorFromBase: 1e-3 },
      { unit: "MJ/kmol", conversionFactorFromBase: 1e-6 },
      { unit: "J/mol", conversionFactorFromBase: 1e-3 },
      { unit: "kJ/mol", conversionFactorFromBase: 1e-6 },
      { unit: "MJ/mol", conversionFactorFromBase: 1e-9 },
    ],
    symbol: "ΔH_c°(net)",
    color: constantColors[2 % constantColors.length]
  },
  {
    displayName: "Lennard-Jones Diameter",
    jsonKey: "Lennard Jones diameter",
    targetUnitName: "m",
    availableUnits: [
      { unit: "m", conversionFactorFromBase: 1 },
      { unit: "Å", conversionFactorFromBase: 1e10 }, // Angstrom
      { unit: "nm", conversionFactorFromBase: 1e9 },  // Nanometer
    ],
    symbol: "σ_LJ",
    color: constantColors[3 % constantColors.length]
  },
  {
    displayName: "Lennard-Jones Energy",
    jsonKey: "Lennard Jones energy", // Given as K (epsilon/k_B)
    targetUnitName: "K",
    availableUnits: [{ unit: "K", conversionFactorFromBase: 1 }],
    symbol: "ε/k_B",
    color: constantColors[4 % constantColors.length]
  },
  {
    displayName: "Fuller et al. Diffusion Volume",
    jsonKey: "Fuller et al. diffusion volume",
    targetUnitName: "-", // Dimensionless
    availableUnits: [{ unit: "-", conversionFactorFromBase: 1 }],
    symbol: "Σv_Fuller",
    color: constantColors[5 % constantColors.length]
  },
  {
    displayName: "Parachor",
    jsonKey: "Parachor", // Units are complex: (density_g_cm3 * MW_g_mol) * (surface_tension_dyn_cm)^0.25 is not standard.
                          // DWSIM database seems to store it in (kg^0.25 * m^3) / (s^0.5 * kmol)
                          // Or mN/m * (cm3/mol)^4 = 10^-3 N/m * (10^-6 m3/mol)^4 ... this is not it.
                          // Standard definition: P = M/ρ * γ^(1/4) where M is molar mass, ρ density, γ surface tension.
                          // Units: (g/mol) / (g/cm³) * (dyn/cm)^(1/4) = cm³/mol * (dyn/cm)^(1/4)
                          // Let's assume the database unit is the one provided: kg⁰·²⁵·m³/s⁰·⁵/kmol
    targetUnitName: "kg∜·m³/(s√·kmol)", // Updated to use root symbols
    availableUnits: [{ unit: "kg∜·m³/(s√·kmol)", conversionFactorFromBase: 1 }], // Updated
    color: baseColors[16],
    symbol: "P_chor"
  },
  {
    displayName: "Solubility Parameter",
    jsonKey: "Solubility parameter",
    targetUnitName: "Pa√", // Updated to use root symbol
    availableUnits: [
      { unit: "Pa√", conversionFactorFromBase: 1 }, // (J/m³)^0.5
      { unit: "(MPa)√", displayName: "MPa√", conversionFactorFromBase: 1e-3 }, // (MPa)^0.5
      { unit: "(cal/cm³)√", displayName: "(cal/cm³)√", conversionFactorFromBase: 1/2.045484e3 } // (cal/cm³)^0.5, using the factor from previous comments
    ],
    color: baseColors[0], // Re-using a color, consider a unique one if available
    symbol: "δ"
  },
  {
    displayName: "UNIQUAC r (Van der Waals volume)",
    jsonKey: "UNIQUAC r (Van der Waals volume)",
    targetUnitName: "m³/kmol",
    availableUnits: [
      { unit: "m³/kmol", conversionFactorFromBase: 1 },
      { unit: "cm³/mol", conversionFactorFromBase: 1 }, // m³/kmol is numerically equal to cm³/mol
      { unit: "Å³/molecule", conversionFactorFromBase: 1 / (6.022e23 * 1e-30) }, // 1 m³/kmol * (1 kmol / 6.022e26 molecules) * (1e30 Å³/m³)
                                                                              // = (1e30 / 6.022e26) Å³/molecule = (10000 / 6.022) = 1660.5 Å³/molecule
                                                                              // So 1 Å³/molecule = 1 / 1660.5 m³/kmol
    ],
    symbol: "r_uniquac",
    color: baseColors[1],
  },
  {
    displayName: "UNIQUAC q (Van der Waals area)",
    jsonKey: "UNIQUAC q (Van der Waals area)",
    targetUnitName: "m²/kmol",
    availableUnits: [
      { unit: "m²/kmol", conversionFactorFromBase: 1 },
      { unit: "Å²/molecule", conversionFactorFromBase: 1 / (6.022e23 * 1e-20) }, // 1 m²/kmol * (1 kmol / 6.022e26 molecules) * (1e20 Å²/m²)
                                                                              // = (1e20 / 6.022e26) Å²/molecule = (1 / 6.022e6) Å²/molecule
                                                                              // So 1 Å²/molecule = 6.022e6 m²/kmol
    ],
    symbol: "q_uniquac",
    color: baseColors[2],
  },
];

