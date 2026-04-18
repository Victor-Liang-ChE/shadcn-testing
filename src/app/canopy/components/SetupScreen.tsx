'use client';

import { useState, useCallback } from 'react';
import { useCanopyStore, type UnitSystem } from '@/lib/canopy/store';
import type { FluidPackage, CanopyCompound } from '@/lib/canopy/types';
import type { InteractionParams, UNIFACSubgroup, UNIFACData } from '@/lib/canopy/thermo';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import {
  Sheet,
  SheetContent,
  SheetDescription,
  SheetHeader,
  SheetTitle,
} from '@/components/ui/sheet';
import { ArrowRight, X, FlaskConical, Beaker, Search, Plus, Ruler, AlertTriangle } from 'lucide-react';
import { formatCompoundName } from '@/lib/antoine-utils';
import { supabase } from '@/lib/supabaseClient';
import { UNIFAC_DORTMUND_PARAMS } from '@/lib/canopy/store';
import {
  buildInteractionParamsForPackage,
  type AspenBinaryParamRow,
} from '@/lib/canopy/interaction-params';
import {
  getEditableInteractionMatrices,
  getFluidPackageDiagnostics,
} from '@/lib/canopy/property-package';
import {
  loadRemoteCompound,
  loadRemoteInteractionParams,
  loadRemotePetroAssay,
  searchRemoteCompounds,
  searchRemotePetroAssays,
  type PropertyProvenanceRecord,
  type RemotePetroAssaySummary,
} from '@/lib/canopy/property-service';
import {
  createExampleCrudeAssay,
  generatePseudoComponentsFromAssay,
} from '@/lib/canopy/petro';

const FLUID_PACKAGES: { value: FluidPackage; label: string; description: string; implemented: boolean }[] = [
  { value: 'Ideal', label: 'Ideal (Raoult\'s Law)', description: 'Best for ideal mixtures of similar compounds (e.g. benzene + toluene).', implemented: true },
  { value: 'NRTL', label: 'NRTL', description: 'Non-random two-liquid model. Excellent for polar/non-ideal systems.', implemented: true },
  { value: 'Electrolyte-NRTL', label: 'Electrolyte-NRTL', description: 'Aspen-style electrolyte package with Debye-Huckel long-range plus NRTL short-range terms.', implemented: true },
  { value: 'Wilson', label: 'Wilson', description: 'Activity coefficient model for miscible systems.', implemented: true },
  { value: 'UNIQUAC', label: 'UNIQUAC', description: 'Universal quasi-chemical model for polar & non-polar mixtures.', implemented: true },
  { value: 'UNIFAC', label: 'UNIFAC', description: 'Group contribution method — no binary parameters needed.', implemented: true },
  { value: 'UNIFAC-DMD', label: 'Dortmund UNIFAC', description: 'Modified UNIFAC with temperature-dependent interactions for stronger non-ideal VLE work.', implemented: true },
  { value: 'Peng-Robinson', label: 'Peng-Robinson', description: 'Cubic EOS for high-pressure and gas-phase systems.', implemented: true },
  { value: 'SRK', label: 'SRK', description: 'Soave-Redlich-Kwong EOS for hydrocarbon systems.', implemented: true },
];

export default function SetupScreen() {
  const {
    compounds, fluidPackage, setFluidPackage,
    addCompound, removeCompound, setScreen,
    unitSystem, setUnitSystem, setInteractionParams,
    dbWarnings, addDbWarning, clearDbWarnings,
  } = useCanopyStore();

  const [searchQuery, setSearchQuery] = useState('');
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [searching, setSearching] = useState(false);
  const [petroSheetOpen, setPetroSheetOpen] = useState(false);
  const [petroQuery, setPetroQuery] = useState('');
  const [petroLoading, setPetroLoading] = useState(false);
  const [petroSuggestions, setPetroSuggestions] = useState<RemotePetroAssaySummary[]>([]);
  const [compoundProvenance, setCompoundProvenance] = useState<Record<string, PropertyProvenanceRecord[]>>({});
  const [packageProvenance, setPackageProvenance] = useState<PropertyProvenanceRecord[]>([]);

  /** Fetch binary interaction params from DB for current compounds + fluid package */
  const fetchBinaryParamsLegacy = useCallback(async (pkg: FluidPackage, comps: CanopyCompound[]) => {
    if (comps.length < 2 || pkg === 'Ideal') return;

    // Group-contribution packages: fetch subgroup definitions, then wire package-specific interactions
    if (pkg === 'UNIFAC' || pkg === 'UNIFAC-DMD') {
      // Collect all subgroup IDs from compounds' unifacGroups
      const allSubgroupIds = new Set<number>();
      const compGroups: Record<number, number>[] = [];
      for (const c of comps) {
        const groups = c.unifacGroups ?? {};
        compGroups.push(groups);
        for (const sgId of Object.keys(groups)) allSubgroupIds.add(Number(sgId));
      }
      if (allSubgroupIds.size === 0) return;

      try {
        // Fetch Rk, Qk, Main Group # for each subgroup
        const { data: rkData } = await supabase
          .from('UNIFAC - Rk and Qk')
          .select('"Subgroup #", "Main Group #", "Rk", "Qk"')
          .in('"Subgroup #"', Array.from(allSubgroupIds));

        if (!rkData || rkData.length === 0) return;

        const subgroups: UNIFACSubgroup[] = [];
        const mainIds = new Set<number>();
        for (const row of rkData as any[]) {
          const sg = Number(row['Subgroup #']);
          const mg = Number(row['Main Group #']);
          subgroups.push({ subgroupId: sg, mainGroupId: mg, R: Number(row.Rk), Q: Number(row.Qk) });
          mainIds.add(mg);
        }

        if (pkg === 'UNIFAC') {
          // Fetch interaction params a(ij) for all relevant main groups
          const mainArr = Array.from(mainIds);
          const { data: intData } = await supabase
            .from('UNIFAC - a(ij)')
            .select('i,j,"a(ij)"')
            .in('i', mainArr)
            .in('j', mainArr);

          const interactions: Record<string, number> = {};
          if (intData) {
            for (const row of intData as any[]) {
              interactions[`${row.i}-${row.j}`] = Number(row['a(ij)']) || 0;
            }
          }

          const unifacData: UNIFACData = { subgroups, interactions };
          setInteractionParams({ unifac: { compGroups, data: unifacData } } as InteractionParams);
          return;
        }

        const dmdInteractions: Record<string, { a: number; b: number; c: number }> = {};
        for (const [key, val] of Object.entries(UNIFAC_DORTMUND_PARAMS)) {
          dmdInteractions[key.replace('|', '-')] = val;
        }
        const coveredMainGroups = new Set<number>();
        for (const key of Object.keys(UNIFAC_DORTMUND_PARAMS)) {
          const [i, j] = key.split('|').map(Number);
          coveredMainGroups.add(i);
          coveredMainGroups.add(j);
        }
        const unsupportedMainGroups = Array.from(mainIds).filter(id => !coveredMainGroups.has(id));
        if (unsupportedMainGroups.length > 0) {
          addDbWarning(
            `Dortmund UNIFAC is using the built-in Aspen interaction subset. Main groups ${unsupportedMainGroups.join(', ')} will fall back to zero interaction where no parameters are available.`
          );
        }
        setInteractionParams({
          unifacDmd: {
            compGroups,
            data: { subgroups, interactions: dmdInteractions },
          },
        } as InteractionParams);
      } catch (e) {
        console.error(`Error fetching ${pkg} params:`, e);
        addDbWarning(`${pkg} group parameters could not be loaded from the database. Default values will be used.`);
      }
      return;
    }

    // Map fluid package to DB Property filter
    const propMap: Record<string, string> = {
      NRTL: 'NRTL', Wilson: 'WILSON', UNIQUAC: 'UNIQ',
      'Peng-Robinson': 'PRKA', SRK: 'SRKA',
    };
    const dbProp = propMap[pkg];
    if (!dbProp) return;

    const names = comps.map(c => c.name);

    try {
      // Fetch all binary pair rows matching our compounds & property
      const { data } = await supabase
        .from('apv140_binary_params')
        .select('Compound1Name,Compound2Name,ElementLabel,Value')
        .in('Compound1Name', names)
        .in('Compound2Name', names)
        .ilike('Property', `${dbProp}%`)
        .order('Compound1Name')
        .limit(1000);

      if (!data || data.length === 0) {
        addDbWarning(`No ${pkg} binary interaction parameters found for the selected compounds. Zero interaction parameters will be used.`);
        return;
      }

      const parsedParams = buildInteractionParamsForPackage(pkg, comps, data as AspenBinaryParamRow[]);
      if (parsedParams) {
        setInteractionParams(parsedParams);
      }
    } catch (e) {
      console.error('Error fetching binary params:', e);
      addDbWarning(`${pkg} binary interaction parameters could not be loaded. The table may not exist yet — results will use zero interaction parameters.`);
    }
  }, [setInteractionParams, addDbWarning]);

  const fetchBinaryParams = useCallback(async (pkg: FluidPackage, comps: CanopyCompound[]) => {
    try {
      const result = await loadRemoteInteractionParams(pkg, comps);
      setPackageProvenance(result.provenance);
      for (const warning of result.warnings) addDbWarning(warning);
      if (result.params) setInteractionParams(result.params);
    } catch (e) {
      console.error('Error fetching binary params via property service:', e);
      addDbWarning(`${pkg} interaction parameters could not be loaded from the remote property service.`);
    }
  }, [setInteractionParams, addDbWarning]);

  const handleSearch = async (query: string) => {
    setSearchQuery(query);
    if (query.trim().length < 2) {
      setSuggestions([]);
      return;
    }
    setSearching(true);
    try {
      setSuggestions(await searchRemoteCompounds(query, compounds.map(c => c.name)));
    } catch (e) {
      console.error('Search error:', e);
    }
    setSearching(false);
  };

  const handleAddCompound = async (name: string) => {
    // Fetch therm data from DB
    try {
      const result = await loadRemoteCompound(name);
      for (const warning of result.warnings) addDbWarning(warning);
      if (!result.compound) {
        alert(`Compound "${name}" not found in database.`);
        return;
      }
      addCompound(result.compound);
      setCompoundProvenance(prev => ({ ...prev, [name]: result.provenance }));
      setSearchQuery('');
      setSuggestions([]);
      fetchBinaryParams(fluidPackage, [...compounds, result.compound]);
      return;

      // Fetch pure props, Antoine, and compound_properties in parallel
      const [pureRes, antoineRes, propsRes] = await Promise.all([
        supabase.from('apv140_pure_props_wide').select('*').ilike('CompoundName', name).limit(1),
        supabase.from('apv140_plxant_wide').select('*').ilike('CompoundName', name).limit(1),
        supabase.from('compound_properties').select('properties').ilike('name', name).limit(1),
      ]);

      const pure = pureRes.data?.[0] as any;
      const ant = antoineRes.data?.[0] as any;
      const props = (propsRes.data?.[0] as any)?.properties;

      if (!pure) {
        alert(`Compound "${name}" not found in database.`);
        return;
      }

      // Extract DIPPR coefficients from compound_properties JSON
      const dipprCoeffs: Record<string, { A: number; B: number; C: number; D: number; E?: number }> = {};
      if (props) {
        for (const key of ['IdealGasHeatCapacityCp', 'HeatOfVaporization', 'LiquidHeatCapacityCp']) {
          const propData = props[key];
          if (propData?.A?.value !== undefined) {
            dipprCoeffs[key] = {
              A: parseFloat(propData.A.value),
              B: parseFloat(propData.B?.value ?? '0'),
              C: parseFloat(propData.C?.value ?? '0'),
              D: parseFloat(propData.D?.value ?? '0'),
              E: propData.E?.value ? parseFloat(propData.E.value) : undefined,
            };
          }
        }
      }

      const mw = props?.MolecularWeight?.value
        ? parseFloat(props.MolecularWeight.value)
        : pure.MW ?? 0;

      // DIPPR 107 (Aly-Lee) ideal gas Cp from pure props wide (CPIGDP_A..E)
      let cpigdp: CanopyCompound['cpigdp'];
      if (pure.CPIGDP_A != null) {
        cpigdp = {
          A: Number(pure.CPIGDP_A),
          B: Number(pure.CPIGDP_B ?? 0),
          C: Number(pure.CPIGDP_C ?? 0),
          D: Number(pure.CPIGDP_D ?? 0),
          E: Number(pure.CPIGDP_E ?? 0),
        };
      }

      // UNIQUAC r/q from pure props wide (GMUQR, GMUQQ columns)
      const uniquac_r = pure.GMUQR !== null && pure.GMUQR !== undefined ? Number(pure.GMUQR) : undefined;
      const uniquac_q = pure.GMUQQ !== null && pure.GMUQQ !== undefined ? Number(pure.GMUQQ) : undefined;

      // UNIFAC groups from compound_properties JSON
      let unifacGroups: Record<number, number> | undefined;
      if (props?.UnifacVLE?.group) {
        const parsedUnifacGroups: Record<number, number> = {};
        const groups = Array.isArray(props.UnifacVLE.group)
          ? props.UnifacVLE.group
          : [props.UnifacVLE.group];
        for (const g of groups) {
          const id = parseInt(g.id ?? g['@id'] ?? g.groupId, 10);
          const count = parseInt(g.count ?? g['@count'] ?? g.value ?? '1', 10);
          if (!isNaN(id) && !isNaN(count)) parsedUnifacGroups[id] = count;
        }
        unifacGroups = Object.keys(parsedUnifacGroups).length > 0 ? parsedUnifacGroups : undefined;
      }

      // Transport DIPPR coefficients from compound_properties JSON
      const transportCoeffs: NonNullable<CanopyCompound['transportCoeffs']> = {};
      const transportMap: Record<string, keyof NonNullable<CanopyCompound['transportCoeffs']>> = {
        LiquidViscosity: 'LiquidViscosity',
        VaporViscosity: 'VaporViscosity',
        LiquidThermalConductivity: 'LiquidThermalConductivity',
        VaporThermalConductivity: 'VaporThermalConductivity',
        SurfaceTension: 'SurfaceTension',
      };
      if (props) {
        for (const [jsonKey, tsKey] of Object.entries(transportMap)) {
          const propData = props[jsonKey];
          if (propData?.A?.value !== undefined) {
            (transportCoeffs as any)[tsKey] = {
              A: parseFloat(propData.A.value),
              B: parseFloat(propData.B?.value ?? '0'),
              C: parseFloat(propData.C?.value ?? '0'),
              D: parseFloat(propData.D?.value ?? '0'),
              E: propData.E?.value ? parseFloat(propData.E.value) : undefined,
            };
          }
        }
      }
      // Also try pure_props_wide DIPPR columns as fallback
      if (!transportCoeffs.LiquidViscosity && pure.MULDIP_A != null) {
        transportCoeffs.LiquidViscosity = {
          A: Number(pure.MULDIP_A), B: Number(pure.MULDIP_B ?? 0),
          C: Number(pure.MULDIP_C ?? 0), D: Number(pure.MULDIP_D ?? 0),
          E: pure.MULDIP_E != null ? Number(pure.MULDIP_E) : undefined,
        };
      }
      if (!transportCoeffs.SurfaceTension && pure.SIGDIP_A != null) {
        transportCoeffs.SurfaceTension = {
          A: Number(pure.SIGDIP_A), B: Number(pure.SIGDIP_B ?? 0),
          C: Number(pure.SIGDIP_C ?? 0), D: Number(pure.SIGDIP_D ?? 0),
          E: pure.SIGDIP_E != null ? Number(pure.SIGDIP_E) : undefined,
        };
      }
      if (!transportCoeffs.LiquidThermalConductivity && pure.KLDIP_A != null) {
        transportCoeffs.LiquidThermalConductivity = {
          A: Number(pure.KLDIP_A), B: Number(pure.KLDIP_B ?? 0),
          C: Number(pure.KLDIP_C ?? 0), D: Number(pure.KLDIP_D ?? 0),
          E: pure.KLDIP_E != null ? Number(pure.KLDIP_E) : undefined,
        };
      }
      if (!transportCoeffs.VaporViscosity && pure.MUVDIP_A != null) {
        transportCoeffs.VaporViscosity = {
          A: Number(pure.MUVDIP_A), B: Number(pure.MUVDIP_B ?? 0),
          C: Number(pure.MUVDIP_C ?? 0), D: Number(pure.MUVDIP_D ?? 0),
          E: pure.MUVDIP_E != null ? Number(pure.MUVDIP_E) : undefined,
        };
      }

      const hasTransport = Object.keys(transportCoeffs).length > 0;

      // Heat of formation & Gibbs energy of formation from compound_properties
      const Hf_raw = props?.HeatOfFormation?.value;
      const Gf_raw = props?.GibbsEnergyOfFormation?.value;
      // compound_properties stores in J/kmol → convert to J/mol
      const Hf_298_Jmol = Hf_raw ? parseFloat(Hf_raw) / 1000 : undefined;
      const Gf_298_Jmol = Gf_raw ? parseFloat(Gf_raw) / 1000 : undefined;

      const compound: CanopyCompound = {
        name,
        displayName: formatCompoundName(name),
        molecularWeight: mw,
        Tc_K: pure.TC ?? 0,
        Pc_Pa: pure.PC ?? 0,
        omega: pure.OMEGA ?? 0,
        dipprCoeffs,
        antoine: ant ? {
          C1: ant.C1 ?? 0, C2: ant.C2 ?? 0, C3: ant.C3 ?? 0,
          C4: ant.C4 ?? 0, C5: ant.C5 ?? 0, C6: ant.C6 ?? 0,
          C7: ant.C7 ?? 1, Tmin_K: ant.Tmin ?? 0, Tmax_K: ant.Tmax ?? 10000,
        } : null,
        ...(uniquac_r !== undefined ? { uniquac_r } : {}),
        ...(uniquac_q !== undefined ? { uniquac_q } : {}),
        ...(unifacGroups ? { unifacGroups } : {}),
        ...(hasTransport ? { transportCoeffs } : {}),
        ...(cpigdp ? { cpigdp } : {}),
        ...(Hf_298_Jmol !== undefined ? { Hf_298_Jmol } : {}),
        ...(Gf_298_Jmol !== undefined ? { Gf_298_Jmol } : {}),
      };

      addCompound(compound);
      setSearchQuery('');
      setSuggestions([]);
      // Re-fetch binary params with updated compound list
      fetchBinaryParams(fluidPackage, [...compounds, compound]);
    } catch (e) {
      console.error('Error adding compound:', e);
      addDbWarning(`Failed to fetch data for "${name}". Check your database connection.`);
    }
  };

  const registerPseudoComponents = useCallback((
    pseudoComponents: CanopyCompound[],
    provenance: PropertyProvenanceRecord[],
  ) => {
    const existing = new Set(compounds.map(compound => compound.name));
    const newPseudoComponents = pseudoComponents.filter(compound => !existing.has(compound.name));
    if (newPseudoComponents.length === 0) {
      addDbWarning('All pseudo-components from the selected petroleum assay are already loaded.');
      return;
    }

    for (const compound of newPseudoComponents) {
      addCompound(compound);
    }

    setCompoundProvenance(previous => {
      const next = { ...previous };
      for (const compound of newPseudoComponents) {
        next[compound.name] = provenance;
      }
      return next;
    });

    if (fluidPackage === 'Ideal' || fluidPackage === 'UNIFAC' || fluidPackage === 'UNIFAC-DMD') {
      addDbWarning('Petroleum pseudo-components are typically used with Peng-Robinson or SRK. Current package remains unchanged.');
    }
    if (fluidPackage !== 'Ideal') {
      addDbWarning('Generated petroleum pseudo-components do not yet have remote petro binary parameters. EOS methods will use zero interaction terms unless petro tables are uploaded.');
    }
  }, [addCompound, addDbWarning, compounds, fluidPackage]);

  const registerResolvedCompounds = useCallback((
    loadedCompounds: CanopyCompound[],
    provenance: PropertyProvenanceRecord[],
    emptyWarning: string,
  ) => {
    const existing = new Set(compounds.map(compound => compound.name));
    const newCompounds = loadedCompounds.filter(compound => !existing.has(compound.name));
    if (newCompounds.length === 0) {
      addDbWarning(emptyWarning);
      return;
    }

    for (const compound of newCompounds) {
      addCompound(compound);
    }

    setCompoundProvenance(previous => {
      const next = { ...previous };
      for (const compound of newCompounds) {
        next[compound.name] = provenance;
      }
      return next;
    });
  }, [addCompound, addDbWarning, compounds]);

  const handlePetroSearch = async (query: string) => {
    setPetroQuery(query);
    if (query.trim().length < 2) {
      setPetroSuggestions([]);
      return;
    }
    setPetroLoading(true);
    try {
      setPetroSuggestions(await searchRemotePetroAssays(query));
    } catch (error) {
      console.error('Petro assay search error:', error);
      addDbWarning('Remote petroleum assay search failed.');
    } finally {
      setPetroLoading(false);
    }
  };

  const handleLoadExampleAssay = () => {
    const assay = createExampleCrudeAssay();
    const provenance: PropertyProvenanceRecord[] = [
      { label: 'Petroleum Assay', source: 'canopy built-in example', detail: assay.assayName },
      { label: 'Pseudo-components', source: 'canopy petro generator', detail: `${assay.cuts.length} TBP cuts` },
    ];
    registerPseudoComponents(generatePseudoComponentsFromAssay(assay), provenance);
    setPetroSheetOpen(false);
  };

  const handleLoadRemoteAssay = useCallback(async (assayId: string) => {
    setPetroLoading(true);
    try {
      const result = await loadRemotePetroAssay(assayId);
      for (const warning of result.warnings) addDbWarning(warning);
      if (!result.assay || (result.pseudoComponents.length === 0 && result.lightEndCompounds.length === 0)) {
        addDbWarning(`Remote petroleum assay "${assayId}" could not be loaded.`);
        return;
      }
      if (result.pseudoComponents.length > 0) {
        registerPseudoComponents(result.pseudoComponents, result.provenance);
      }
      if (result.lightEndCompounds.length > 0) {
        registerResolvedCompounds(
          result.lightEndCompounds,
          result.provenance,
          `All explicit light-end compounds for assay "${assayId}" are already loaded.`,
        );
      }
      setPetroQuery('');
      setPetroSuggestions([]);
      setPetroSheetOpen(false);
    } catch (error) {
      console.error('Petro assay load error:', error);
      addDbWarning(`Remote petroleum assay "${assayId}" failed to load.`);
    } finally {
      setPetroLoading(false);
    }
  }, [addDbWarning, registerPseudoComponents, registerResolvedCompounds]);

  return (
    <div className="flex-1 overflow-auto p-6">
      <div className="max-w-5xl mx-auto">
        {/* DB Warnings Banner */}
        {dbWarnings.length > 0 && (
          <div className="mb-4 rounded-md border border-amber-500/40 bg-amber-500/10 px-4 py-3">
            <div className="flex items-start justify-between gap-2">
              <div className="flex items-start gap-2">
                <AlertTriangle className="w-4 h-4 mt-0.5 text-amber-600 dark:text-amber-400 shrink-0" />
                <div className="space-y-1">
                  {dbWarnings.map((w, i) => (
                    <p key={i} className="text-xs text-amber-700 dark:text-amber-300">{w}</p>
                  ))}
                </div>
              </div>
              <button
                onClick={clearDbWarnings}
                className="text-amber-600 dark:text-amber-400 hover:text-amber-800 dark:hover:text-amber-200 p-0.5"
              >
                <X className="w-3.5 h-3.5" />
              </button>
            </div>
          </div>
        )}
      </div>
      <div className="max-w-5xl mx-auto grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Compound Selection */}
        <Card className="bg-card border-border">
          <CardHeader>
            <div className="flex items-center justify-between gap-3">
              <CardTitle className="flex items-center gap-2">
                <Beaker className="w-5 h-5" />
                Compounds
              </CardTitle>
              <Button
                type="button"
                variant="outline"
                size="sm"
                onClick={() => setPetroSheetOpen(true)}
              >
                Petroleum Assays
              </Button>
            </div>
            <CardDescription>
              Search and add compounds to your simulation.
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-4">
            {/* Search */}
            <div className="relative">
              <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
              <input
                type="text"
                value={searchQuery}
                onChange={(e) => handleSearch(e.target.value)}
                placeholder="Search compounds (e.g. benzene, water, ethanol...)"
                className="w-full pl-9 pr-3 py-2 rounded-md border border-input bg-background text-sm focus:outline-none focus:ring-2 focus:ring-ring"
              />
              {suggestions.length > 0 && (
                <div className="absolute z-50 top-full left-0 right-0 mt-1 bg-popover border border-border rounded-md shadow-lg max-h-48 overflow-auto">
                  {suggestions.map(name => (
                    <button
                      key={name}
                      onClick={() => handleAddCompound(name)}
                      className="w-full text-left px-3 py-1.5 text-sm hover:bg-accent transition-colors flex items-center gap-2"
                    >
                      <Plus className="w-3 h-3 text-green-500" />
                      {formatCompoundName(name)}
                    </button>
                  ))}
                </div>
              )}
            </div>

            {/* Selected compounds */}
            <div className="space-y-1">
              {compounds.length === 0 && (
                <p className="text-sm text-muted-foreground italic">No compounds added yet.</p>
              )}
              {compounds.map((c, i) => (
                <div key={c.name}
                  className="flex items-center justify-between px-3 py-1.5 rounded-md bg-muted/50 border border-border"
                >
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono text-muted-foreground w-5">{i + 1}.</span>
                    <span className="text-sm font-medium">{c.displayName}</span>
                    <span className="text-xs text-muted-foreground">
                      MW: {c.molecularWeight.toFixed(1)} | T<sub>c</sub>: {c.Tc_K.toFixed(1)} K
                    </span>
                  </div>
                  <button onClick={() => removeCompound(c.name)}
                    className="text-red-400 hover:text-red-600 transition-colors p-0.5">
                    <X className="w-4 h-4" />
                  </button>
                </div>
              ))}
            </div>

            {/* Properties table for selected compounds */}
            {compounds.length > 0 && (
              <div className="mt-3 pt-3 border-t border-border">
                <h4 className="text-xs font-semibold text-muted-foreground mb-2">Physical Properties</h4>
                <div className="overflow-x-auto">
                  <table className="text-xs w-full border-collapse">
                    <thead>
                      <tr className="border-b border-border">
                        <th className="text-left p-1 font-semibold text-muted-foreground">Property</th>
                        {compounds.map(c => (
                          <th key={c.name} className="text-right p-1 font-semibold min-w-[80px]">{c.displayName}</th>
                        ))}
                      </tr>
                    </thead>
                    <tbody>
                      <tr className="border-b border-border/50">
                        <td className="p-1 text-muted-foreground">MW (g/mol)</td>
                        {compounds.map(c => <td key={c.name} className="p-1 text-right font-mono">{c.molecularWeight.toFixed(2)}</td>)}
                      </tr>
                      <tr className="border-b border-border/50">
                        <td className="p-1 text-muted-foreground">T<sub>c</sub> (K)</td>
                        {compounds.map(c => <td key={c.name} className="p-1 text-right font-mono">{c.Tc_K.toFixed(2)}</td>)}
                      </tr>
                      <tr className="border-b border-border/50">
                        <td className="p-1 text-muted-foreground">P<sub>c</sub> (Pa)</td>
                        {compounds.map(c => <td key={c.name} className="p-1 text-right font-mono">{(c.Pc_Pa / 1e5).toFixed(2)} bar</td>)}
                      </tr>
                      <tr className="border-b border-border/50">
                        <td className="p-1 text-muted-foreground">ω (acentric)</td>
                        {compounds.map(c => <td key={c.name} className="p-1 text-right font-mono">{c.omega.toFixed(4)}</td>)}
                      </tr>
                      {compounds.some(c => c.uniquac_r) && (
                        <tr className="border-b border-border/50">
                          <td className="p-1 text-muted-foreground">UNIQUAC r</td>
                          {compounds.map(c => <td key={c.name} className="p-1 text-right font-mono">{c.uniquac_r?.toFixed(3) ?? '—'}</td>)}
                        </tr>
                      )}
                      {compounds.some(c => c.uniquac_q) && (
                        <tr className="border-b border-border/50">
                          <td className="p-1 text-muted-foreground">UNIQUAC q</td>
                          {compounds.map(c => <td key={c.name} className="p-1 text-right font-mono">{c.uniquac_q?.toFixed(3) ?? '—'}</td>)}
                        </tr>
                      )}
                      <tr className="border-b border-border/50">
                        <td className="p-1 text-muted-foreground">Antoine</td>
                        {compounds.map(c => (
                          <td key={c.name} className="p-1 text-right">
                            {c.antoine ? <span className="text-green-500">✓</span> : <span className="text-red-400">✗</span>}
                          </td>
                        ))}
                      </tr>
                      <tr className="border-b border-border/50">
                        <td className="p-1 text-muted-foreground">DIPPR Cp<sub>ig</sub></td>
                        {compounds.map(c => (
                          <td key={c.name} className="p-1 text-right">
                            {c.cpigdp ? <span className="text-green-500">✓</span> : c.dipprCoeffs?.IdealGasHeatCapacityCp ? <span className="text-green-500">✓</span> : <span className="text-red-400">✗</span>}
                          </td>
                        ))}
                      </tr>
                    </tbody>
                  </table>
                </div>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Fluid Package Selection */}
        <Card className="bg-card border-border">
          <CardHeader>
            <CardTitle className="flex items-center gap-2">
              <FlaskConical className="w-5 h-5" />
              Fluid Package
            </CardTitle>
            <CardDescription>
              Select the thermodynamic model for your simulation.
              All compounds use the same fluid package.
            </CardDescription>
          </CardHeader>
          <CardContent className="space-y-2">
            {FLUID_PACKAGES.map(pkg => (
              <button
                key={pkg.value}
                onClick={() => { setFluidPackage(pkg.value); fetchBinaryParams(pkg.value, compounds); }}
                className={`w-full text-left px-4 py-3 rounded-md border transition-all
                  ${fluidPackage === pkg.value
                    ? 'border-blue-500 bg-blue-500/10 shadow-sm'
                    : 'border-border hover:border-muted-foreground/30 hover:bg-accent/50'
                  }`}
              >
                <div className="flex items-center gap-2">
                  <span className="font-medium text-sm">{pkg.label}</span>
                  {!pkg.implemented && (
                    <span className="text-[10px] font-medium px-1.5 py-0.5 rounded bg-amber-500/15 text-amber-600 dark:text-amber-400 border border-amber-500/30">
                      coming soon
                    </span>
                  )}
                </div>
                <div className="text-xs text-muted-foreground mt-0.5">{pkg.description}</div>
              </button>
            ))}
            {fluidPackage !== 'Ideal' && (
              <div className="flex items-start gap-2 mt-3 px-3 py-2.5 rounded-md bg-blue-500/10 border border-blue-500/30 text-blue-700 dark:text-blue-400">
                <FlaskConical className="w-4 h-4 mt-0.5 shrink-0" />
                <p className="text-xs leading-relaxed">
                  <span className="font-semibold">{fluidPackage}</span>{' '}
                  {fluidPackage === 'UNIFAC'
                    ? 'group contribution parameters (Rk, Qk, a_mn) will be auto-loaded from the database when compounds are selected.'
                    : fluidPackage === 'UNIFAC-DMD'
                    ? 'Dortmund UNIFAC subgroup definitions will be auto-loaded, and temperature-dependent Aspen interaction parameters will be applied from the local simulator dataset.'
                    : fluidPackage === 'Electrolyte-NRTL'
                    ? 'short-range NRTL parameters will be loaded remotely, and electrolyte metadata will be checked for charges, solvent identity, and dielectric constant.'
                    : 'binary interaction parameters will be auto-loaded from the Aspen database when compounds are selected.'}
                </p>
              </div>
            )}
          </CardContent>
        </Card>

        {/* Binary Interaction Parameter Matrix */}
        {compounds.length > 0 && fluidPackage !== 'Ideal' && (
          <PropertyMethodStatusCard />
        )}

        {(Object.keys(compoundProvenance).length > 0 || packageProvenance.length > 0) && (
          <PropertyProvenanceCard
            compoundProvenance={compoundProvenance}
            packageProvenance={packageProvenance}
          />
        )}

        {/* Binary Interaction Parameter Matrix */}
        {compounds.length >= 2 && fluidPackage !== 'Ideal' && fluidPackage !== 'UNIFAC' && fluidPackage !== 'UNIFAC-DMD' && (
          <BIPMatrixPanel />
        )}

        {/* Unit System */}
        <Card className="bg-card border-border lg:col-span-2">
          <CardHeader>
            <CardTitle className="flex items-center gap-2">
              <Ruler className="w-5 h-5" />
              Unit System
            </CardTitle>
            <CardDescription>
              Display units for the flowsheet. Internal calculations always use SI.
            </CardDescription>
          </CardHeader>
          <CardContent>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
              <UnitSelect
                label="Temperature"
                value={unitSystem.temperature}
                options={[
                  { value: 'C', label: '°C' },
                  { value: 'K', label: 'K' },
                  { value: 'F', label: '°F' },
                ]}
                onChange={(v) => setUnitSystem({ temperature: v as UnitSystem['temperature'] })}
              />
              <UnitSelect
                label="Pressure"
                value={unitSystem.pressure}
                options={[
                  { value: 'bar', label: 'bar' },
                  { value: 'atm', label: 'atm' },
                  { value: 'kPa', label: 'kPa' },
                  { value: 'Pa', label: 'Pa' },
                  { value: 'psi', label: 'psi' },
                ]}
                onChange={(v) => setUnitSystem({ pressure: v as UnitSystem['pressure'] })}
              />
              <UnitSelect
                label="Flow"
                value={unitSystem.flow}
                options={[
                  { value: 'mol/s', label: 'mol/s' },
                  { value: 'kmol/h', label: 'kmol/h' },
                  { value: 'kg/s', label: 'kg/s' },
                  { value: 'kg/h', label: 'kg/h' },
                ]}
                onChange={(v) => setUnitSystem({ flow: v as UnitSystem['flow'] })}
              />
              <UnitSelect
                label="Energy"
                value={unitSystem.energy}
                options={[
                  { value: 'kW', label: 'kW' },
                  { value: 'W', label: 'W' },
                  { value: 'BTU/h', label: 'BTU/h' },
                  { value: 'kcal/h', label: 'kcal/h' },
                ]}
                onChange={(v) => setUnitSystem({ energy: v as UnitSystem['energy'] })}
              />
            </div>
          </CardContent>
        </Card>
      </div>

      {/* Next button */}
      <div className="max-w-5xl mx-auto mt-6 flex justify-end">
        <Button
          onClick={() => setScreen('flowsheet')}
          disabled={compounds.length === 0}
          size="lg"
          className="bg-green-600 hover:bg-green-700 text-white"
        >
          Next: Flowsheet
          <ArrowRight className="w-4 h-4 ml-2" />
        </Button>
      </div>

      <Sheet open={petroSheetOpen} onOpenChange={setPetroSheetOpen}>
        <SheetContent side="right" className="sm:max-w-lg">
          <SheetHeader>
            <SheetTitle>Petroleum Assays</SheetTitle>
            <SheetDescription>
              Generate pseudo-components from a crude assay locally, or load a remote assay when petro tables are available in Supabase.
            </SheetDescription>
          </SheetHeader>
          <div className="space-y-4 px-4 pb-4">
            <div className="flex flex-col gap-3">
              <Button
                type="button"
                variant="outline"
                onClick={handleLoadExampleAssay}
              >
                Load Example Crude Assay
              </Button>
              <div className="relative flex-1">
                <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-muted-foreground" />
                <input
                  type="text"
                  value={petroQuery}
                  onChange={(e) => handlePetroSearch(e.target.value)}
                  placeholder="Search remote assays (e.g. crude, brent, vacuum gas oil...)"
                  className="w-full pl-9 pr-3 py-2 rounded-md border border-input bg-background text-sm focus:outline-none focus:ring-2 focus:ring-ring"
                />
                {petroSuggestions.length > 0 && (
                  <div className="absolute z-50 top-full left-0 right-0 mt-1 bg-popover border border-border rounded-md shadow-lg max-h-56 overflow-auto">
                    {petroSuggestions.map(assay => (
                      <button
                        key={assay.assayId}
                        onClick={() => handleLoadRemoteAssay(assay.assayId)}
                        className="w-full text-left px-3 py-2 text-sm hover:bg-accent transition-colors"
                      >
                        <div className="font-medium">{assay.assayName}</div>
                        <div className="text-xs text-muted-foreground">
                          {assay.source}
                          {assay.location ? ` | ${assay.location}` : ''}
                          {` | ${assay.cutCount} cuts`}
                        </div>
                      </button>
                    ))}
                  </div>
                )}
              </div>
            </div>

            <div className="rounded-md border border-border/60 px-3 py-2 text-xs text-muted-foreground">
              Built-in assay generation needs no database upload. Remote petro loading expects `petro_assays` and `petro_assay_cuts`, with optional `petro_assay_light_ends`.
            </div>

            {petroLoading && (
              <p className="text-xs text-muted-foreground">Loading petroleum assay data...</p>
            )}
          </div>
        </SheetContent>
      </Sheet>
    </div>
  );
}

function UnitSelect({
  label,
  value,
  options,
  onChange,
}: {
  label: string;
  value: string;
  options: { value: string; label: string }[];
  onChange: (v: string) => void;
}) {
  return (
    <div className="space-y-1">
      <label className="text-xs font-medium text-muted-foreground">{label}</label>
      <select
        value={value}
        onChange={(e) => onChange(e.target.value)}
        className="w-full px-2 py-1.5 text-sm border border-input rounded bg-background focus:outline-none focus:ring-2 focus:ring-ring"
      >
        {options.map(opt => (
          <option key={opt.value} value={opt.value}>{opt.label}</option>
        ))}
      </select>
    </div>
  );
}

function PropertyMethodStatusCard() {
  const { compounds, fluidPackage, interactionParams } = useCanopyStore();
  const diagnostics = getFluidPackageDiagnostics(compounds, fluidPackage, interactionParams);

  return (
    <Card className="bg-card border-border lg:col-span-2">
      <CardHeader>
        <CardTitle className="flex items-center gap-2 text-base">
          <Beaker className="w-5 h-5" />
          Property Method Status
        </CardTitle>
        <CardDescription>
          Aspen-style readiness checks for the selected property package and current component set.
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-2">
        {diagnostics.map((diagnostic, index) => {
          const tone =
            diagnostic.level === 'ok'
              ? 'border-emerald-500/30 bg-emerald-500/10 text-emerald-700 dark:text-emerald-300'
              : diagnostic.level === 'warn'
              ? 'border-amber-500/30 bg-amber-500/10 text-amber-700 dark:text-amber-300'
              : 'border-red-500/30 bg-red-500/10 text-red-700 dark:text-red-300';

          return (
            <div key={`${diagnostic.level}-${index}`} className={`rounded-md border px-3 py-2 text-xs ${tone}`}>
              {diagnostic.message}
            </div>
          );
        })}
      </CardContent>
    </Card>
  );
}

// ─── Binary Interaction Parameter Matrix Editor ────────────────

function BIPMatrixCard() {
  const { compounds, fluidPackage, interactionParams, setInteractionParams } = useCanopyStore();
  const nc = compounds.length;

  // Extract current matrix based on fluid package
  const getMatrix = (): { label: string; matrix: number[][]; field: string }[] => {
    if (fluidPackage === 'NRTL' && interactionParams.nrtl) {
      const zeroMatrix = Array.from({ length: nc }, () => new Array(nc).fill(0));
      return [
        { label: 'A (energy param)', matrix: interactionParams.nrtl.A ?? zeroMatrix, field: 'A' },
        { label: 'B (temperature dep)', matrix: interactionParams.nrtl.B ?? zeroMatrix, field: 'B' },
        { label: 'α (non-randomness)', matrix: interactionParams.nrtl.alpha, field: 'alpha' },
      ];
    }
    if (fluidPackage === 'UNIQUAC' && interactionParams.uniquac) {
      return [
        { label: 'a (energy param)', matrix: interactionParams.uniquac.a, field: 'a' },
        { label: 'b (temperature dep)', matrix: interactionParams.uniquac.b, field: 'b' },
      ];
    }
    if (fluidPackage === 'Wilson' && interactionParams.wilson) {
      return [
        { label: 'A (energy param)', matrix: interactionParams.wilson.A, field: 'A' },
      ];
    }
    if ((fluidPackage === 'Peng-Robinson' || fluidPackage === 'SRK') && interactionParams.kij) {
      return [
        { label: 'kij', matrix: interactionParams.kij, field: 'kij' },
      ];
    }
    // Default: show empty kij matrix
    return [{ label: 'kij', matrix: Array.from({ length: nc }, () => new Array(nc).fill(0)), field: 'kij' }];
  };

  const matrices = getMatrix();

  const handleCellChange = (matIdx: number, i: number, j: number, value: string) => {
    const v = parseFloat(value);
    if (isNaN(v)) return;
    const m = matrices[matIdx];
    const newMatrix = m.matrix.map(row => [...row]);
    newMatrix[i][j] = v;

    // Update the appropriate params
    if (fluidPackage === 'NRTL') {
      const nrtl = { ...interactionParams.nrtl! };
      (nrtl as any)[m.field] = newMatrix;
      setInteractionParams({ ...interactionParams, nrtl });
    } else if (fluidPackage === 'UNIQUAC') {
      const uniquac = { ...interactionParams.uniquac! };
      (uniquac as any)[m.field] = newMatrix;
      setInteractionParams({ ...interactionParams, uniquac });
    } else if (fluidPackage === 'Wilson') {
      const wilson = { ...interactionParams.wilson! };
      (wilson as any)[m.field] = newMatrix;
      setInteractionParams({ ...interactionParams, wilson });
    } else {
      setInteractionParams({ ...interactionParams, kij: newMatrix });
    }
  };

  return (
    <Card className="bg-card border-border lg:col-span-2">
      <CardHeader>
        <CardTitle className="flex items-center gap-2 text-base">
          <FlaskConical className="w-5 h-5" />
          Binary Interaction Parameters — {fluidPackage}
        </CardTitle>
        <CardDescription>
          Loaded from the Aspen database. Edit cells to override values.
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-4">
        {matrices.map((m, mi) => (
          <div key={mi}>
            <h4 className="text-xs font-semibold text-muted-foreground mb-2">{m.label}</h4>
            <div className="overflow-x-auto">
              <table className="text-xs border-collapse">
                <thead>
                  <tr>
                    <th className="p-1.5 text-left text-muted-foreground">i \ j</th>
                    {compounds.map(c => (
                      <th key={c.name} className="p-1.5 text-center font-medium min-w-[70px]">{c.displayName}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {compounds.map((ci, i) => (
                    <tr key={ci.name} className="border-t border-border/50">
                      <td className="p-1.5 font-medium text-muted-foreground">{ci.displayName}</td>
                      {compounds.map((cj, j) => (
                        <td key={cj.name} className="p-0.5">
                          {i === j ? (
                            <div className="w-full text-center text-muted-foreground/50 py-1">—</div>
                          ) : (
                            <input
                              type="number"
                              step="any"
                              value={m.matrix[i]?.[j] ?? 0}
                              onChange={e => handleCellChange(mi, i, j, e.target.value)}
                              className="w-full px-1.5 py-1 text-center text-xs border border-input rounded bg-background focus:outline-none focus:ring-1 focus:ring-ring font-mono"
                            />
                          )}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        ))}
        {matrices[0]?.matrix.length === 0 && (
          <p className="text-xs text-muted-foreground italic">
            No parameters loaded yet. Add compounds and select a fluid package to load binary parameters from the database.
          </p>
        )}
      </CardContent>
    </Card>
  );
}

function PropertyProvenanceCard({
  compoundProvenance,
  packageProvenance,
}: {
  compoundProvenance: Record<string, PropertyProvenanceRecord[]>;
  packageProvenance: PropertyProvenanceRecord[];
}) {
  return (
    <Card className="bg-card border-border lg:col-span-2">
      <CardHeader>
        <CardTitle className="flex items-center gap-2 text-base">
          <Ruler className="w-5 h-5" />
          Property Data Sources
        </CardTitle>
        <CardDescription>
          Remote property-service provenance for loaded compounds and the active package.
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-4 text-xs">
        {Object.entries(compoundProvenance).map(([compoundName, records]) => (
          <div key={compoundName}>
            <h4 className="font-semibold text-foreground mb-1">{formatCompoundName(compoundName)}</h4>
            <div className="space-y-1">
              {records.map((record, index) => (
                <div key={`${compoundName}-${record.source}-${index}`} className="rounded border border-border/60 px-3 py-2">
                  <span className="font-medium">{record.label}:</span> {record.source}
                  <span className="text-muted-foreground"> — {record.detail}</span>
                </div>
              ))}
            </div>
          </div>
        ))}
        {packageProvenance.length > 0 && (
          <div>
            <h4 className="font-semibold text-foreground mb-1">Active Package</h4>
            <div className="space-y-1">
              {packageProvenance.map((record, index) => (
                <div key={`${record.source}-${index}`} className="rounded border border-border/60 px-3 py-2">
                  <span className="font-medium">{record.label}:</span> {record.source}
                  <span className="text-muted-foreground"> — {record.detail}</span>
                </div>
              ))}
            </div>
          </div>
        )}
      </CardContent>
    </Card>
  );
}

function BIPMatrixPanel() {
  const { compounds, fluidPackage, interactionParams, setInteractionParams } = useCanopyStore();
  const matrices = getEditableInteractionMatrices(fluidPackage, interactionParams, compounds.length);

  const setMatrixAtPath = (field: string, matrix: number[][]) => {
    if (field === 'kij' || field === 'kij_a' || field === 'kij_b' || field === 'kij_c') {
      setInteractionParams({ [field]: matrix });
      return;
    }

    const [section, key] = field.split('.');
    if (!section || !key) return;

    if (section === 'nrtl' && interactionParams.nrtl) {
      setInteractionParams({ nrtl: { ...interactionParams.nrtl, [key]: matrix } });
    } else if (section === 'uniquac' && interactionParams.uniquac) {
      setInteractionParams({ uniquac: { ...interactionParams.uniquac, [key]: matrix } });
    } else if (section === 'wilson' && interactionParams.wilson) {
      setInteractionParams({ wilson: { ...interactionParams.wilson, [key]: matrix } });
    }
  };

  const handleCellChange = (field: string, matrix: number[][], i: number, j: number, value: string) => {
    const numeric = parseFloat(value);
    if (Number.isNaN(numeric)) return;
    const updated = matrix.map(row => [...row]);
    updated[i][j] = numeric;
    setMatrixAtPath(field, updated);
  };

  return (
    <Card className="bg-card border-border lg:col-span-2">
      <CardHeader>
        <CardTitle className="flex items-center gap-2 text-base">
          <FlaskConical className="w-5 h-5" />
          Binary Interaction Parameters — {fluidPackage}
        </CardTitle>
        <CardDescription>
          Loaded from the Aspen database using the active canopy package convention. Edit cells to override values.
        </CardDescription>
      </CardHeader>
      <CardContent className="space-y-4">
        {matrices.length === 0 && (
          <p className="text-xs text-muted-foreground italic">
            No editable binary matrices are available for this package yet.
          </p>
        )}
        {matrices.map((matrixDescriptor) => (
          <div key={matrixDescriptor.field}>
            <h4 className="text-xs font-semibold text-muted-foreground mb-2">{matrixDescriptor.label}</h4>
            <div className="overflow-x-auto">
              <table className="text-xs border-collapse">
                <thead>
                  <tr>
                    <th className="p-1.5 text-left text-muted-foreground">i \ j</th>
                    {compounds.map(c => (
                      <th key={c.name} className="p-1.5 text-center font-medium min-w-[70px]">{c.displayName}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {compounds.map((ci, i) => (
                    <tr key={ci.name} className="border-t border-border/50">
                      <td className="p-1.5 font-medium text-muted-foreground">{ci.displayName}</td>
                      {compounds.map((cj, j) => (
                        <td key={cj.name} className="p-0.5">
                          {i === j ? (
                            <div className="w-full text-center text-muted-foreground/50 py-1">—</div>
                          ) : (
                            <input
                              type="number"
                              step="any"
                              value={matrixDescriptor.matrix[i]?.[j] ?? 0}
                              onChange={e => handleCellChange(matrixDescriptor.field, matrixDescriptor.matrix, i, j, e.target.value)}
                              className="w-full px-1.5 py-1 text-center text-xs border border-input rounded bg-background focus:outline-none focus:ring-1 focus:ring-ring font-mono"
                            />
                          )}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        ))}
      </CardContent>
    </Card>
  );
}
