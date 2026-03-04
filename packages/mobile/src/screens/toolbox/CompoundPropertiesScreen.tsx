import React, { useState, useCallback, useMemo, useRef } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity,
  FlatList, ActivityIndicator, useWindowDimensions,
} from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import { useTheme } from '../../theme/ThemeContext';
import SectionCard from '../../components/SectionCard';
import SliderRow from '../../components/SliderRow';
import { supabase } from '../../lib/supabaseClient';

// ─── Property definitions ────────────────────────────────────────────────────
type PropKey = 'VaporPressure' | 'LiquidDensity' | 'LiquidViscosity' |
  'IdealGasHeatCapacityCp' | 'LiquidHeatCapacityCp' | 'HeatOfVaporization' |
  'VaporThermalConductivity' | 'LiquidThermalConductivity';

interface PropDef {
  label: string;
  key: PropKey;
  unit: string;
  color: string;
  eval: (coeffs: Record<string, number>, T: number, Tc: number) => number;
}

const PROPS: PropDef[] = [
  {
    label: 'Vapor Pressure', key: 'VaporPressure', unit: 'Pa', color: '#60a5fa',
    eval: ({ A, B, C, D, E }, T) => Math.exp(A + B / T + C * Math.log(T) + D * Math.pow(T, E)),
  },
  {
    label: 'Liquid Density', key: 'LiquidDensity', unit: 'kmol/m³', color: '#22c55e',
    eval: ({ A, B, C, D }, T, Tc) => A / Math.pow(B, 1 + Math.pow(1 - T / Tc, D)),
  },
  {
    label: 'Liquid Viscosity', key: 'LiquidViscosity', unit: 'Pa·s', color: '#f97316',
    eval: ({ A, B, C, D, E }, T) => Math.exp(A + B / T + C * Math.log(T) + D * Math.pow(T, E)),
  },
  {
    label: 'Liq. Heat Capacity', key: 'LiquidHeatCapacityCp', unit: 'J/kmol/K', color: '#a855f7',
    eval: ({ A, B, C, D, E }, T) => A + Math.exp(B / T + C + D * T + E * T * T),
  },
  {
    label: 'Heat of Vaporization', key: 'HeatOfVaporization', unit: 'J/kmol', color: '#ec4899',
    eval: ({ A, B, C, D, E }, T, Tc) => {
      const tr = T / Tc;
      return A * Math.pow(1 - tr, B + C * tr + D * tr * tr + E * tr * tr * tr);
    },
  },
  {
    label: 'IG Heat Capacity', key: 'IdealGasHeatCapacityCp', unit: 'J/kmol/K', color: '#facc15',
    eval: ({ A, B, C, D, E }, T) => A + B * Math.pow(C / T / Math.sinh(C / T), 2) + D * Math.pow(E / T / Math.cosh(E / T), 2),
  },
];

const N_PTS = 80;

interface CompoundData {
  name: string;
  properties: Record<string, any>;
}

function getCoeffs(props: Record<string, any>, key: PropKey): Record<string, number> {
  const p = props[key];
  if (!p || !p.coefficients) return {};
  const c = p.coefficients;
  const toN = (v: any) => Number(v);
  return { A: toN(c[0]), B: toN(c[1]), C: toN(c[2]), D: toN(c[3]), E: toN(c[4]) };
}

function getTc(props: Record<string, any>): number {
  return Number(props?.CriticalTemperature?.value) || 500;
}

function getTRange(props: Record<string, any>, key: PropKey): [number, number] {
  const p = props[key];
  const lo = Number(p?.Tmin);
  const hi = Number(p?.Tmax);
  return [isFinite(lo) && lo > 0 ? lo : 200, isFinite(hi) && hi > lo ? hi : 600];
}

export default function CompoundPropertiesScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const CHART_W = width - 64;

  const [query, setQuery] = useState('');
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [compound, setCompound] = useState<CompoundData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [selectedProp, setSelectedProp] = useState<PropKey>('VaporPressure');
  const [tMin, setTMin] = useState(250);
  const [tMax, setTMax] = useState(500);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

  const propDef = PROPS.find(p => p.key === selectedProp)!;

  async function fetchSuggestions(text: string) {
    if (text.length < 2) { setSuggestions([]); return; }
    const { data } = await supabase
      .from('compound_properties')
      .select('name')
      .ilike('name', `${text}%`)
      .limit(6);
    setSuggestions(data?.map((d: any) => d.name) ?? []);
  }

  function onChangeText(text: string) {
    setQuery(text);
    if (debounceRef.current !== null) clearTimeout(debounceRef.current);
    debounceRef.current = setTimeout(() => fetchSuggestions(text), 300);
  }

  async function selectCompound(name: string) {
    setQuery(name);
    setSuggestions([]);
    setLoading(true);
    setError(null);
    try {
      const { data, error: err } = await supabase
        .from('compound_properties')
        .select('name, properties')
        .ilike('name', name)
        .limit(1)
        .single();
      if (err || !data) throw new Error(err?.message ?? 'Not found');
      setCompound(data as CompoundData);
      const range = getTRange(data.properties, selectedProp);
      setTMin(Math.round(range[0]));
      setTMax(Math.round(range[1]));
    } catch (e: any) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }

  const chartData = useMemo(() => {
    if (!compound) return [];
    const coeffs = getCoeffs(compound.properties, selectedProp);
    const Tc = getTc(compound.properties);
    if (!Object.keys(coeffs).length) return [];
    const step = (tMax - tMin) / N_PTS;
    const pts = [];
    const stride = Math.max(1, Math.floor(N_PTS / 5));
    for (let i = 0; i <= N_PTS; i++) {
      const T = tMin + i * step;
      try {
        const val = propDef.eval(coeffs, T, Tc);
        pts.push({ value: isFinite(val) ? val : 0, label: i % stride === 0 ? T.toFixed(0) : '' });
      } catch {
        pts.push({ value: 0, label: '' });
      }
    }
    return pts;
  }, [compound, selectedProp, tMin, tMax, propDef]);

  // Summary constants
  const constants = useMemo(() => {
    if (!compound) return [];
    const p = compound.properties;
    const items: { label: string; value: string }[] = [];
    const mw = Number(p.MolecularWeight?.value);
    const tc = Number(p.CriticalTemperature?.value);
    const pc = Number(p.CriticalPressure?.value);
    const omega = Number(p.AcentricFactor?.value);
    const tb = Number(p.NormalBoilingPoint?.value);
    if (isFinite(mw) && mw > 0) items.push({ label: 'MW', value: `${mw.toFixed(2)} g/mol` });
    if (isFinite(tc) && tc > 0) items.push({ label: 'Tc', value: `${tc.toFixed(1)} K` });
    if (isFinite(pc) && pc > 0) items.push({ label: 'Pc', value: `${(pc / 1e5).toFixed(2)} bar` });
    if (isFinite(omega)) items.push({ label: 'ω', value: omega.toFixed(3) });
    if (isFinite(tb) && tb > 0) items.push({ label: 'Tb', value: `${tb.toFixed(1)} K` });
    return items;
  }, [compound]);

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}
      keyboardShouldPersistTaps="handled">
      <SectionCard title="Compound Search">
        <TextInput
          value={query}
          onChangeText={onChangeText}
          placeholder="Type compound name..."
          placeholderTextColor={colors.mutedForeground}
          style={[styles.input, { color: colors.foreground, borderColor: colors.border, backgroundColor: colors.card }]}
        />
        {suggestions.length > 0 && (
          <View style={[styles.dropdown, { backgroundColor: colors.card, borderColor: colors.border }]}>
            {suggestions.map(s => (
              <TouchableOpacity key={s} onPress={() => selectCompound(s)} style={styles.suggestion}>
                <Text style={{ color: colors.foreground, fontSize: 14 }}>{s}</Text>
              </TouchableOpacity>
            ))}
          </View>
        )}
        {loading && <ActivityIndicator style={{ marginTop: 8 }} color={colors.blue} />}
        {error && <Text style={styles.errorText}>{error}</Text>}
      </SectionCard>

      {compound && (
        <>
          {constants.length > 0 && (
            <SectionCard title={`${compound.name} — Constants`}>
              <View style={styles.constantsGrid}>
                {constants.map(c => (
                  <View key={c.label} style={[styles.constItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
                    <Text style={[styles.constLabel, { color: colors.mutedForeground }]}>{c.label}</Text>
                    <Text style={[styles.constValue, { color: colors.foreground }]}>{c.value}</Text>
                  </View>
                ))}
              </View>
            </SectionCard>
          )}

          <SectionCard title="Property">
            <ScrollView horizontal showsHorizontalScrollIndicator={false}>
              <View style={styles.chipRow}>
                {PROPS.map(p => (
                  <TouchableOpacity key={p.key} onPress={() => setSelectedProp(p.key)}
                    style={[styles.chip, { backgroundColor: selectedProp === p.key ? p.color : colors.card, borderColor: selectedProp === p.key ? p.color : colors.border }]}>
                    <Text style={[styles.chipText, { color: selectedProp === p.key ? '#fff' : colors.foreground }]}>{p.label}</Text>
                  </TouchableOpacity>
                ))}
              </View>
            </ScrollView>
          </SectionCard>

          <SectionCard title="Temperature Range (K)">
            <SliderRow label="T min (K)" value={tMin} min={100} max={tMax - 5} step={5} decimals={0} onChange={setTMin} />
            <SliderRow label="T max (K)" value={tMax} min={tMin + 5} max={1500} step={5} decimals={0} onChange={setTMax} />
          </SectionCard>

          {chartData.length > 0 && (
            <SectionCard title={`${propDef.label} vs Temperature`}>
              <Text style={[styles.unitLabel, { color: colors.mutedForeground }]}>Units: {propDef.unit}</Text>
              <LineChart
                data={chartData}
                width={CHART_W}
                height={220}
                color={propDef.color}
                thickness={2.5}
                xAxisColor={colors.border}
                yAxisColor={colors.border}
                xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
                yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
                dataPointsRadius={0}
                hideDataPoints
                initialSpacing={0}
                endSpacing={4}
                noOfSections={4}
                curved
              />
              <Text style={[styles.xLabel, { color: colors.mutedForeground }]}>Temperature (K)</Text>
            </SectionCard>
          )}
        </>
      )}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  input: { borderWidth: 1, borderRadius: 10, padding: 12, fontSize: 14 },
  dropdown: { borderWidth: 1, borderRadius: 10, marginTop: 4, overflow: 'hidden' },
  suggestion: { paddingVertical: 10, paddingHorizontal: 14, borderBottomWidth: StyleSheet.hairlineWidth },
  errorText: { color: '#ef4444', fontSize: 13, marginTop: 6 },
  constantsGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 8 },
  constItem: { borderWidth: 1, borderRadius: 8, padding: 10, minWidth: 80, alignItems: 'center' },
  constLabel: { fontSize: 11, marginBottom: 2 },
  constValue: { fontSize: 14, fontWeight: '700' },
  chipRow: { flexDirection: 'row', gap: 8, paddingVertical: 2 },
  chip: { paddingHorizontal: 12, paddingVertical: 7, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 12, fontWeight: '600' },
  unitLabel: { fontSize: 12, marginBottom: 6 },
  xLabel: { fontSize: 11, textAlign: 'center', marginTop: 4 },
});
