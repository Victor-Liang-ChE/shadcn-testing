import type { CanopyCompound, FluidPackage } from './types';
import type { InteractionParams } from './thermo';

export interface PackageMatrixDescriptor {
  field: string;
  label: string;
  matrix: number[][];
}

export interface PackageDiagnostic {
  level: 'ok' | 'warn' | 'error';
  message: string;
}

function createMatrix(size: number, fill = 0): number[][] {
  return Array.from({ length: size }, () => new Array(size).fill(fill));
}

function hasNonZeroMatrix(matrix?: number[][]): boolean {
  if (!matrix) return false;
  for (let i = 0; i < matrix.length; i++) {
    for (let j = 0; j < matrix[i].length; j++) {
      if (Math.abs(matrix[i][j] ?? 0) > 1e-18) return true;
    }
  }
  return false;
}

function countBinaryPairs(count: number): number {
  return count > 1 ? (count * (count - 1)) / 2 : 0;
}

export function getEditableInteractionMatrices(
  fluidPackage: FluidPackage,
  interactionParams: InteractionParams,
  compoundCount: number,
): PackageMatrixDescriptor[] {
  if (fluidPackage === 'NRTL' && interactionParams.nrtl) {
    const { nrtl } = interactionParams;
    if (
      hasNonZeroMatrix(nrtl.a_ext) ||
      hasNonZeroMatrix(nrtl.b_ext) ||
      hasNonZeroMatrix(nrtl.e_ext) ||
      hasNonZeroMatrix(nrtl.f_ext) ||
      hasNonZeroMatrix(nrtl.d_alpha)
    ) {
      return [
        { label: 'aij', matrix: nrtl.a_ext ?? createMatrix(compoundCount), field: 'nrtl.a_ext' },
        { label: 'bij', matrix: nrtl.b_ext ?? createMatrix(compoundCount), field: 'nrtl.b_ext' },
        { label: 'alpha (cij)', matrix: nrtl.alpha, field: 'nrtl.alpha' },
        { label: 'dij', matrix: nrtl.d_alpha ?? createMatrix(compoundCount), field: 'nrtl.d_alpha' },
        { label: 'eij', matrix: nrtl.e_ext ?? createMatrix(compoundCount), field: 'nrtl.e_ext' },
        { label: 'fij', matrix: nrtl.f_ext ?? createMatrix(compoundCount), field: 'nrtl.f_ext' },
      ];
    }
    return [
      { label: 'A', matrix: nrtl.A ?? createMatrix(compoundCount), field: 'nrtl.A' },
      { label: 'B', matrix: nrtl.B ?? createMatrix(compoundCount), field: 'nrtl.B' },
      { label: 'alpha', matrix: nrtl.alpha, field: 'nrtl.alpha' },
    ];
  }

  if (fluidPackage === 'Electrolyte-NRTL' && interactionParams.elecNrtl) {
    return getEditableInteractionMatrices('NRTL', { nrtl: interactionParams.elecNrtl.nrtl }, compoundCount);
  }

  if (fluidPackage === 'UNIQUAC' && interactionParams.uniquac) {
    const { uniquac } = interactionParams;
    if (
      hasNonZeroMatrix(uniquac.a_ext) ||
      hasNonZeroMatrix(uniquac.b_ext) ||
      hasNonZeroMatrix(uniquac.c_ext) ||
      hasNonZeroMatrix(uniquac.d_ext) ||
      hasNonZeroMatrix(uniquac.e_ext)
    ) {
      return [
        { label: 'aij', matrix: uniquac.a_ext ?? createMatrix(compoundCount), field: 'uniquac.a_ext' },
        { label: 'bij', matrix: uniquac.b_ext ?? createMatrix(compoundCount), field: 'uniquac.b_ext' },
        { label: 'cij', matrix: uniquac.c_ext ?? createMatrix(compoundCount), field: 'uniquac.c_ext' },
        { label: 'dij', matrix: uniquac.d_ext ?? createMatrix(compoundCount), field: 'uniquac.d_ext' },
        { label: 'eij', matrix: uniquac.e_ext ?? createMatrix(compoundCount), field: 'uniquac.e_ext' },
      ];
    }
    return [
      { label: 'a', matrix: uniquac.a, field: 'uniquac.a' },
      { label: 'b', matrix: uniquac.b, field: 'uniquac.b' },
    ];
  }

  if (fluidPackage === 'Wilson' && interactionParams.wilson) {
    const { wilson } = interactionParams;
    if (
      hasNonZeroMatrix(wilson.a_ext) ||
      hasNonZeroMatrix(wilson.b_ext) ||
      hasNonZeroMatrix(wilson.c_ext) ||
      hasNonZeroMatrix(wilson.d_ext) ||
      hasNonZeroMatrix(wilson.e_ext)
    ) {
      return [
        { label: 'aij', matrix: wilson.a_ext ?? createMatrix(compoundCount), field: 'wilson.a_ext' },
        { label: 'bij', matrix: wilson.b_ext ?? createMatrix(compoundCount), field: 'wilson.b_ext' },
        { label: 'cij', matrix: wilson.c_ext ?? createMatrix(compoundCount), field: 'wilson.c_ext' },
        { label: 'dij', matrix: wilson.d_ext ?? createMatrix(compoundCount), field: 'wilson.d_ext' },
        { label: 'eij', matrix: wilson.e_ext ?? createMatrix(compoundCount), field: 'wilson.e_ext' },
      ];
    }
    return [{ label: 'A', matrix: wilson.A, field: 'wilson.A' }];
  }

  if (fluidPackage === 'Peng-Robinson' || fluidPackage === 'SRK') {
    const matrices: PackageMatrixDescriptor[] = [
      {
        label: 'kij',
        matrix: interactionParams.kij ?? createMatrix(compoundCount),
        field: 'kij',
      },
    ];
    if (interactionParams.kij_a) {
      matrices.push({ label: 'aij', matrix: interactionParams.kij_a, field: 'kij_a' });
    }
    if (interactionParams.kij_b) {
      matrices.push({ label: 'bij', matrix: interactionParams.kij_b, field: 'kij_b' });
    }
    if (interactionParams.kij_c) {
      matrices.push({ label: 'cij', matrix: interactionParams.kij_c, field: 'kij_c' });
    }
    return matrices;
  }

  return [];
}

export function getFluidPackageDiagnostics(
  compounds: CanopyCompound[],
  fluidPackage: FluidPackage,
  interactionParams: InteractionParams,
): PackageDiagnostic[] {
  if (compounds.length === 0) {
    return [{ level: 'warn', message: 'Add compounds to evaluate property-method readiness.' }];
  }

  const diagnostics: PackageDiagnostic[] = [];
  const pairCount = countBinaryPairs(compounds.length);
  const namesMissingAntoine = compounds.filter(c => !c.antoine).map(c => c.displayName);
  const namesMissingCritical = compounds
    .filter(c => !(c.Tc_K > 0 && c.Pc_Pa > 0 && Number.isFinite(c.omega)))
    .map(c => c.displayName);

  if (fluidPackage === 'Ideal') {
    if (namesMissingAntoine.length > 0) {
      diagnostics.push({
        level: 'error',
        message: `Missing vapor-pressure data for: ${namesMissingAntoine.join(', ')}.`,
      });
    } else {
      diagnostics.push({
        level: 'ok',
        message: 'Ideal package is ready. Vapor-pressure data is available for all selected compounds.',
      });
    }
    return diagnostics;
  }

  if (fluidPackage === 'NRTL') {
    if (!interactionParams.nrtl) {
      diagnostics.push({ level: 'error', message: 'NRTL binary parameters have not been loaded.' });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `NRTL parameters loaded for ${pairCount} binary pair${pairCount === 1 ? '' : 's'}.`,
      });
      if (!hasNonZeroMatrix(interactionParams.nrtl.a_ext) && !hasNonZeroMatrix(interactionParams.nrtl.b_ext)) {
        diagnostics.push({
          level: 'warn',
          message: 'NRTL is using zero interaction energy terms. Results will trend toward ideal behavior.',
        });
      }
    }
  } else if (fluidPackage === 'Electrolyte-NRTL') {
    if (!interactionParams.elecNrtl) {
      diagnostics.push({ level: 'error', message: 'Electrolyte-NRTL parameters have not been loaded.' });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `Electrolyte-NRTL short-range parameters loaded for ${pairCount} binary pair${pairCount === 1 ? '' : 's'}.`,
      });
      const chargedSpecies = compounds.filter(c => (c.chargeNumber ?? 0) !== 0).map(c => c.displayName);
      const solvent = compounds.find(c => c.isElectrolyteSolvent) ?? compounds.find(c => c.name.toUpperCase() === 'WATER');
      if (chargedSpecies.length === 0) {
        diagnostics.push({
          level: 'warn',
          message: 'No charged species are defined, so Electrolyte-NRTL will reduce to its short-range NRTL contribution.',
        });
      }
      if (!solvent) {
        diagnostics.push({
          level: 'error',
          message: 'No electrolyte solvent has been identified for Electrolyte-NRTL.',
        });
      } else if (!(interactionParams.elecNrtl.epsilon_r > 0)) {
        diagnostics.push({
          level: 'error',
          message: `Dielectric constant data is missing for solvent ${solvent.displayName}.`,
        });
      }
    }
  } else if (fluidPackage === 'Wilson') {
    const missingVolumes = compounds
      .filter((_, i) => !(interactionParams.wilson?.molarVolumes?.[i] && interactionParams.wilson.molarVolumes[i] > 0))
      .map(c => c.displayName);
    if (!interactionParams.wilson) {
      diagnostics.push({ level: 'error', message: 'Wilson binary parameters have not been loaded.' });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `Wilson parameters loaded for ${pairCount} binary pair${pairCount === 1 ? '' : 's'}.`,
      });
      if (missingVolumes.length > 0) {
        diagnostics.push({
          level: 'error',
          message: `Wilson liquid molar volumes are missing for: ${missingVolumes.join(', ')}.`,
        });
      }
    }
  } else if (fluidPackage === 'UNIQUAC') {
    const missingRQ = compounds.filter(c => c.uniquac_r == null || c.uniquac_q == null).map(c => c.displayName);
    if (!interactionParams.uniquac) {
      diagnostics.push({ level: 'error', message: 'UNIQUAC binary parameters have not been loaded.' });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `UNIQUAC parameters loaded for ${pairCount} binary pair${pairCount === 1 ? '' : 's'}.`,
      });
    }
    if (missingRQ.length > 0) {
      diagnostics.push({
        level: 'error',
        message: `UNIQUAC structural parameters r/q are missing for: ${missingRQ.join(', ')}.`,
      });
    }
  } else if (fluidPackage === 'UNIFAC' || fluidPackage === 'UNIFAC-DMD') {
    const missingGroups = compounds.filter(c => !c.unifacGroups || Object.keys(c.unifacGroups).length === 0).map(c => c.displayName);
    const dataLoaded = fluidPackage === 'UNIFAC' ? !!interactionParams.unifac : !!interactionParams.unifacDmd;
    if (missingGroups.length > 0) {
      diagnostics.push({
        level: 'error',
        message: `UNIFAC subgroup definitions are missing for: ${missingGroups.join(', ')}.`,
      });
    }
    if (!dataLoaded) {
      diagnostics.push({
        level: 'error',
        message: `${fluidPackage} group-interaction data has not been loaded.`,
      });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `${fluidPackage} subgroup definitions and interaction data are loaded.`,
      });
    }
  } else if (fluidPackage === 'Peng-Robinson' || fluidPackage === 'SRK') {
    if (namesMissingCritical.length > 0) {
      diagnostics.push({
        level: 'error',
        message: `Critical property data is missing for: ${namesMissingCritical.join(', ')}.`,
      });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `${fluidPackage} has Tc/Pc/omega data for all selected compounds.`,
      });
    }
    if (!interactionParams.kij && !interactionParams.kij_a) {
      diagnostics.push({
        level: 'warn',
        message: `${fluidPackage} is using zero binary interaction parameters.`,
      });
    } else if (interactionParams.kij_a) {
      diagnostics.push({
        level: 'ok',
        message: `${fluidPackage} temperature-dependent binary interaction terms are loaded.`,
      });
    } else {
      diagnostics.push({
        level: 'ok',
        message: `${fluidPackage} constant binary interaction parameters are loaded.`,
      });
    }
  }

  if (namesMissingAntoine.length > 0 && fluidPackage !== 'Peng-Robinson' && fluidPackage !== 'SRK') {
    diagnostics.push({
      level: 'warn',
      message: `Missing vapor-pressure data will limit phase-equilibrium calculations for: ${namesMissingAntoine.join(', ')}.`,
    });
  }

  if (diagnostics.length === 0) {
    diagnostics.push({ level: 'warn', message: 'No diagnostics are available for this package yet.' });
  }

  return diagnostics;
}
