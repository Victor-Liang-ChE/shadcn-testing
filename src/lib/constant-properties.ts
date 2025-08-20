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
    jsonKey: "MolecularWeight",
    targetUnitName: "kg/kmol", // Assuming the base value from DB/calculation is in kg/kmol (or g/mol, which is numerically the same)
    availableUnits: [
      { unit: "g/mol", conversionFactorFromBase: 1, displayName: "g/mol" }, // Default
      { unit: "kg/kmol", conversionFactorFromBase: 1, displayName: "kg/kmol" },
      { unit: "kg/mol", conversionFactorFromBase: 0.001, displayName: "kg/mol" }, // 1 kg/kmol = 0.001 kg/mol
      { unit: "Da", conversionFactorFromBase: 1, displayName: "Da" }, // Changed displayName
      { unit: "amu", conversionFactorFromBase: 1, displayName: "amu" } // Added amu
    ],
    color: constantColors[0 % constantColors.length],
    symbol: "MW"
  },
  {
    displayName: "Critical Temperature",
    jsonKey: "CriticalTemperature",
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
    jsonKey: "CriticalPressure",
    targetUnitName: "Pa", // Changed from bar
    availableUnits: [
      { unit: "bar", conversionFactorFromBase: 1e-5 },     // 1 Pa = 0.00001 bar
      { unit: "Pa", conversionFactorFromBase: 1 },
      { unit: "kPa", conversionFactorFromBase: 1e-3 },     // 1 Pa = 0.001 kPa
      { unit: "atm", conversionFactorFromBase: 1 / 101325 }  // 1 Pa = 1/101325 atm
    ],
    symbol: "P_c",
    color: constantColors[2 % constantColors.length]
  },
  {
    displayName: "Critical Volume",
    jsonKey: "CriticalVolume",
    targetUnitName: "m³/kmol", // Changed from cm³/mol
    availableUnits: [
      { unit: "cm³/mol", conversionFactorFromBase: 1000 },    // 1 m³/kmol = 1000 cm³/mol
      { unit: "m³/kmol", conversionFactorFromBase: 1 },
      { unit: "L/mol", conversionFactorFromBase: 1 },        // 1 m³/kmol = 1 L/mol
    ],
    symbol: "V_c",
    color: constantColors[3 % constantColors.length]
  },
  {
    displayName: "Critical Compressibility Factor",
    jsonKey: "CriticalCompressibility",
    targetUnitName: "-", // Dimensionless
    availableUnits: [{ unit: "-", conversionFactorFromBase: 1 }],
    symbol: "Z_c",
    color: constantColors[4 % constantColors.length]
  },
  {
    displayName: "Normal Boiling Point",
    jsonKey: "NormalBoilingPointTemperature",
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
    jsonKey: "TriplePointTemperature",
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
    jsonKey: "TriplePointPressure",
    targetUnitName: "Pa", // Changed from bar
    availableUnits: [
      { unit: "Pa", conversionFactorFromBase: 1 },
      { unit: "kPa", conversionFactorFromBase: 1e-3 },    // 1 Pa = 0.001 kPa
      { unit: "bar", conversionFactorFromBase: 1e-5 }     // 1 Pa = 0.00001 bar
    ],
    symbol: "P_tp",
    color: constantColors[7 % constantColors.length]
  },
  {
    displayName: "Normal Melting Point",
    jsonKey: "NormalMeltingPointTemperature",
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
    displayName: "Heat of Formation",
    jsonKey: "HeatOfFormation",
    targetUnitName: "J/kmol",
    availableUnits: [
      { unit: "J/kmol", conversionFactorFromBase: 1 },
      { unit: "kJ/kmol", conversionFactorFromBase: 1e-3 },
      { unit: "MJ/kmol", conversionFactorFromBase: 1e-6 },
      { unit: "J/mol", conversionFactorFromBase: 1e-3 },
      { unit: "kJ/mol", conversionFactorFromBase: 1e-6 },
      { unit: "MJ/mol", conversionFactorFromBase: 1e-9 },
    ],
    symbol: "ΔH_f°",
    color: constantColors[1 % constantColors.length]
  },
  {
    displayName: "Melting Point Heat of Fusion",
    jsonKey: "HeatOfFusionAtMeltingPoint",
    targetUnitName: "J/kmol", // Changed from J/mol
    availableUnits: [
      { unit: "J/mol", conversionFactorFromBase: 1e-3 },     // 1 J/kmol = 0.001 J/mol
      { unit: "J/kmol", conversionFactorFromBase: 1 },
      { unit: "kJ/kmol", conversionFactorFromBase: 1e-3 },   // 1 J/kmol = 0.001 kJ/kmol
      { unit: "kJ/mol", conversionFactorFromBase: 1e-6 }     // 1 J/kmol = 0.000001 kJ/mol
    ],
    symbol: "ΔH_fus",
    color: constantColors[1 % constantColors.length]
  },
  {
    displayName: "Heat of Combustion (LHV)",
    jsonKey: "HeatOfCombustion",
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
    displayName: "Acentric Factor",
    jsonKey: "AcentricityFactor",
    targetUnitName: "-", // Dimensionless
    availableUnits: [{ unit: "-", conversionFactorFromBase: 1 }],
    symbol: "ω",
    color: constantColors[3 % constantColors.length]
  },
  {
    displayName: "Dipole Moment",
    jsonKey: "DipoleMoment",
    targetUnitName: "C·m",
    availableUnits: [
      { unit: "C·m", conversionFactorFromBase: 1 },
      { unit: "D", conversionFactorFromBase: 1/(3.33564e-30) }, // 1 Debye = 3.33564e-30 C·m
    ],
    symbol: "μ",
    color: constantColors[4 % constantColors.length]
  },
  {
    displayName: "Lennard-Jones Diameter",
    jsonKey: "DiameterLJ",
    targetUnitName: "m", // Changed target unit to meters
    availableUnits: [
      { unit: "Å", conversionFactorFromBase: 1e10, displayName: "Å" }, // 1 m = 1e10 Å
      { unit: "m", conversionFactorFromBase: 1 },  // Base unit
      { unit: "nm", conversionFactorFromBase: 1e9 },   // 1 m = 1e9 nm
      { unit: "pm", conversionFactorFromBase: 1e12 },  // 1 m = 1e12 pm
    ],
    symbol: "σ_LJ", // Symbol was σ_LJ, kept it consistent with the other Lennard-Jones sigma
    color: constantColors[3 % constantColors.length]
  },
  {
    displayName: "Lennard-Jones Energy",
    jsonKey: "EnergyLJ", // Given as K (epsilon/k_B)
    targetUnitName: "K",
    availableUnits: [{ unit: "K", conversionFactorFromBase: 1 }],
    symbol: "ε/k_B",
    color: constantColors[4 % constantColors.length]
  },
  {
    displayName: "Fuller Diffusion Volume",
    jsonKey: "FullerVolume",
    targetUnitName: "-", // Dimensionless
    availableUnits: [{ unit: "-", conversionFactorFromBase: 1 }],
    symbol: "Σv_Fuller",
    color: constantColors[5 % constantColors.length]
  },
  {
    displayName: "Parachor",
    jsonKey: "Parachor",
    targetUnitName: "∜kg·m³/(√s·kmol)", // Corrected unit string
    availableUnits: [
      { unit: "∜kg·m³/(√s·kmol)", conversionFactorFromBase: 1, displayName: "∜kg·m³/(√s·kmol)" }, // Base unit
      { 
        unit: "∜g·cm³/(√s·mol)", 
        conversionFactorFromBase: 10**3.75, // (1000 g/kg)^(1/4) * (1e6 cm³/m³) / (1000 mol/kmol) = 10^(3/4) * 1e6 / 1e3 = 10^(0.75) * 1e3 = 5.623 * 1000 = 5623 approx.
                                            // More precisely: (1000)^(1/4) * 100^3 / 1000 = 10^(3/4) * 10^6 / 10^3 = 10^(0.75) * 10^3.
                                            // (kg/1000)^0.25 * (m*100)^3 / (s^0.5 * kmol*1000)
                                            // kg^0.25 * m^3 / (s^0.5 * kmol)  to g^0.25 * cm^3 / (s^0.5 * mol)
                                            // (1000g)^0.25 * (100cm)^3 / (s^0.5 * 1000mol)
                                            // = 1000^0.25 * 100^3 / 1000  * [g^0.25 cm^3 / (s^0.5 mol)]
                                            // = 10^(3*0.25) * 10^6 / 10^3 = 10^0.75 * 10^3 = 10^3.75
        displayName: "∜g·cm³/(√s·mol)" 
      },
      { 
        unit: "∜kg·cm³/(√s·kmol)", 
        conversionFactorFromBase: 1e6, // m³ to cm³
        displayName: "∜kg·cm³/(√s·kmol)" 
      }
    ], 
    color: baseColors[16],
    symbol: "P_chor"
  },
  {
    displayName: "Solubility Parameter",
    jsonKey: "SolubilityParameter",
    targetUnitName: "√Pa", // Updated to √Pa
    availableUnits: [
      { unit: "√Pa", conversionFactorFromBase: 1, displayName: "√Pa" }, // Updated unit and displayName
      { unit: "√MPa", displayName: "√MPa", conversionFactorFromBase: 1e-3 }, // Updated unit and displayName
      { unit: "cal/cm³√", displayName: "√(cal/cm³)", conversionFactorFromBase: 1/2.045484e3 } // Kept as is, already in √Unit form
    ],
    color: baseColors[0], // Re-using a color, consider a unique one if available
    symbol: "δ"
  },
  {
    displayName: "UNIQUAC r",
    jsonKey: "UniquacR",
    targetUnitName: "-",
    availableUnits: [
      { unit: "-", conversionFactorFromBase: 1 }
    ],
    symbol: "r_UNIQUAC",
    color: baseColors[1],
  },
  {
    displayName: "UNIQUAC q",
    jsonKey: "UniquacQ",
    targetUnitName: "-",
    availableUnits: [
      { unit: "-", conversionFactorFromBase: 1 }
    ],
    symbol: "q_UNIQUAC",
    color: baseColors[2],
  },
  {
    displayName: "Van der Waals Volume",
    jsonKey: "VanDerWaalsVolume",
    targetUnitName: "m³/kmol",
    availableUnits: [
      { unit: "m³/kmol", conversionFactorFromBase: 1 },
      { unit: "cm³/mol", conversionFactorFromBase: 1000 },
      { unit: "L/mol", conversionFactorFromBase: 1 },
    ],
    symbol: "V_vdW",
    color: constantColors[1],
  },
  {
    displayName: "Van der Waals Area",
    jsonKey: "VanDerWaalsArea",
    targetUnitName: "m²/kmol",
    availableUnits: [
      { unit: "m²/kmol", conversionFactorFromBase: 1 },
      { unit: "cm²/mol", conversionFactorFromBase: 1e7 }, // 1 m²/kmol = 10^7 cm²/mol
    ],
    symbol: "A_vdW",
    color: constantColors[2],
  }
];

