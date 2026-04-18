import type { CanopyCompound } from './types';

export interface AssayCutInput {
  cutName?: string;
  tbpStart_K: number;
  tbpEnd_K: number;
  liquidVolumeFraction: number;
  specificGravity: number;
}

export interface PetroLightEndInput {
  componentName: string;
  fractionValue: number;
  fractionBasis?: string;
  normalBoilingPoint_K?: number;
}

export interface PetroAssay {
  assayName: string;
  cuts: AssayCutInput[];
  lightEnds?: PetroLightEndInput[];
}

function clamp(value: number, lower: number, upper: number): number {
  return Math.max(lower, Math.min(upper, value));
}

function pseudoCriticalTemperature_K(Tb_K: number, sg: number): number {
  return Math.max(Tb_K * (1.55 - 0.2 * clamp(sg, 0.6, 1.1)), Tb_K + 80);
}

function pseudoCriticalPressure_Pa(Tb_K: number, sg: number): number {
  const baseBar = 60 - 0.06 * (Tb_K - 300) - 20 * (sg - 0.75);
  return clamp(baseBar, 8, 45) * 1e5;
}

function pseudoAcentricFactor(Tb_K: number, Tc_K: number): number {
  const Trb = clamp(Tb_K / Math.max(Tc_K, Tb_K + 1), 0.45, 0.95);
  return clamp((3 / 7) * (Math.log(1 / Math.max(Trb, 1e-6)) - 1), 0.15, 1.2);
}

function pseudoMolecularWeight(Tb_K: number, sg: number): number {
  return clamp(0.012 * Tb_K * (10 + 8 * sg), 40, 1200);
}

function pseudoLiquidVolumeAtNBP(MW_gpmol: number, sg: number): number {
  const rho = Math.max(sg * 1000, 500);
  return (MW_gpmol / 1000) / rho;
}

export function generatePseudoComponentsFromAssay(assay: PetroAssay): CanopyCompound[] {
  const totalVol = assay.cuts.reduce((sum, cut) => sum + Math.max(cut.liquidVolumeFraction, 0), 0) || 1;
  return assay.cuts.map((cut, index) => {
    const sg = clamp(cut.specificGravity, 0.55, 1.2);
    const tbMid_K = (cut.tbpStart_K + cut.tbpEnd_K) / 2;
    const Tc_K = pseudoCriticalTemperature_K(tbMid_K, sg);
    const Pc_Pa = pseudoCriticalPressure_Pa(tbMid_K, sg);
    const omega = pseudoAcentricFactor(tbMid_K, Tc_K);
    const MW = pseudoMolecularWeight(tbMid_K, sg);
    const volFrac = Math.max(cut.liquidVolumeFraction, 0) / totalVol;
    const name = `${assay.assayName.toUpperCase().replace(/[^A-Z0-9]+/g, '-')}-PC${String(index + 1).padStart(2, '0')}`;
    const displayName = cut.cutName?.trim() || `${assay.assayName} Cut ${index + 1}`;
    const Vb_m3pmol = pseudoLiquidVolumeAtNBP(MW, sg);

    return {
      name,
      displayName,
      molecularWeight: MW,
      Tc_K,
      Pc_Pa,
      omega,
      dipprCoeffs: {},
      antoine: null,
      Tb_K: tbMid_K,
      Vb_m3pmol,
      specificGravity: sg,
      VLSTD_m3pkmol: Vb_m3pmol * 1000,
      anilinePoint_K: clamp(250 + 80 * (1 - sg), 220, 420),
      rackett_ZRA: clamp(0.27 + 0.04 * (sg - 0.75), 0.2, 0.35),
      omegaCostald: clamp(omega + 0.05, 0.15, 1.25),
      cpigPoly: {
        A: (120 + 220 * volFrac) * 1000,
        B: 180,
        C: 0.25,
        D: 0,
        E: 0,
      },
      dipoleMoment_Cm: 0,
    };
  });
}

export function createExampleCrudeAssay(): PetroAssay {
  return {
    assayName: 'Example Crude',
    cuts: [
      { cutName: 'Naphtha', tbpStart_K: 320, tbpEnd_K: 430, liquidVolumeFraction: 0.22, specificGravity: 0.70 },
      { cutName: 'Kerosene', tbpStart_K: 430, tbpEnd_K: 520, liquidVolumeFraction: 0.18, specificGravity: 0.78 },
      { cutName: 'Diesel', tbpStart_K: 520, tbpEnd_K: 620, liquidVolumeFraction: 0.27, specificGravity: 0.84 },
      { cutName: 'VGO', tbpStart_K: 620, tbpEnd_K: 760, liquidVolumeFraction: 0.21, specificGravity: 0.91 },
      { cutName: 'Resid', tbpStart_K: 760, tbpEnd_K: 900, liquidVolumeFraction: 0.12, specificGravity: 0.99 },
    ],
  };
}
