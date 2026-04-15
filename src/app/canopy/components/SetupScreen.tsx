'use client';

import { useState } from 'react';
import { useCanopyStore, type UnitSystem } from '@/lib/canopy/store';
import type { FluidPackage, CanopyCompound } from '@/lib/canopy/types';
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from '@/components/ui/card';
import { Button } from '@/components/ui/button';
import { ArrowRight, X, FlaskConical, Beaker, Search, Plus, Ruler, AlertTriangle } from 'lucide-react';
import { supabase } from '@/lib/supabaseClient';
import { formatCompoundName } from '@/lib/antoine-utils';

const FLUID_PACKAGES: { value: FluidPackage; label: string; description: string; implemented: boolean }[] = [
  { value: 'Ideal', label: 'Ideal (Raoult\'s Law)', description: 'Best for ideal mixtures of similar compounds (e.g. benzene + toluene).', implemented: true },
  { value: 'NRTL', label: 'NRTL', description: 'Non-random two-liquid model. Excellent for polar/non-ideal systems.', implemented: false },
  { value: 'Wilson', label: 'Wilson', description: 'Activity coefficient model for miscible systems.', implemented: false },
  { value: 'UNIQUAC', label: 'UNIQUAC', description: 'Universal quasi-chemical model for polar & non-polar mixtures.', implemented: false },
  { value: 'UNIFAC', label: 'UNIFAC', description: 'Group contribution method — no binary parameters needed.', implemented: false },
  { value: 'Peng-Robinson', label: 'Peng-Robinson', description: 'Cubic EOS for high-pressure and gas-phase systems.', implemented: false },
  { value: 'SRK', label: 'SRK', description: 'Soave-Redlich-Kwong EOS for hydrocarbon systems.', implemented: false },
];

export default function SetupScreen() {
  const {
    compounds, fluidPackage, setFluidPackage,
    addCompound, removeCompound, setScreen,
    unitSystem, setUnitSystem,
  } = useCanopyStore();

  const [searchQuery, setSearchQuery] = useState('');
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [searching, setSearching] = useState(false);

  const handleSearch = async (query: string) => {
    setSearchQuery(query);
    if (query.trim().length < 2) {
      setSuggestions([]);
      return;
    }
    setSearching(true);
    try {
      const { data } = await supabase
        .from('compounds_master')
        .select('"Name"')
        .ilike('Name', `${query.trim()}%`)
        .limit(40);
      if (data) {
        const seen = new Set<string>();
        const names: string[] = [];
        for (const row of data as any[]) {
          const name = row.Name as string;
          if (!seen.has(name) && !compounds.some(c => c.name === name)) {
            seen.add(name);
            names.push(name);
          }
          if (names.length >= 8) break;
        }
        setSuggestions(names);
      }
    } catch (e) {
      console.error('Search error:', e);
    }
    setSearching(false);
  };

  const handleAddCompound = async (name: string) => {
    // Fetch therm data from DB
    try {
      // Fetch pure props
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
      // DB format: { IdealGasHeatCapacityCp: { A: {value: "34010.24"}, B: {value: "-588"}, eqno: {value: "16"}, ... } }
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
      };

      addCompound(compound);
      setSearchQuery('');
      setSuggestions([]);
    } catch (e) {
      console.error('Error adding compound:', e);
    }
  };

  return (
    <div className="flex-1 overflow-auto p-6">
      <div className="max-w-5xl mx-auto grid grid-cols-1 lg:grid-cols-2 gap-6">
        {/* Compound Selection */}
        <Card className="bg-card border-border">
          <CardHeader>
            <CardTitle className="flex items-center gap-2">
              <Beaker className="w-5 h-5" />
              Compounds
            </CardTitle>
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
                onClick={() => setFluidPackage(pkg.value)}
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
              <div className="flex items-start gap-2 mt-3 px-3 py-2.5 rounded-md bg-amber-500/10 border border-amber-500/30 text-amber-700 dark:text-amber-400">
                <AlertTriangle className="w-4 h-4 mt-0.5 shrink-0" />
                <p className="text-xs leading-relaxed">
                  <span className="font-semibold">{fluidPackage}</span> is not yet implemented.
                  The simulator will fall back to <span className="font-semibold">Ideal (Raoult&apos;s Law)</span> for all calculations.
                </p>
              </div>
            )}
          </CardContent>
        </Card>

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
