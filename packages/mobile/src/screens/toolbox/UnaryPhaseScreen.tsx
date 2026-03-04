import React, { useState, useMemo, useRef } from 'react';
import { ScrollView, View, Text, StyleSheet, TextInput, TouchableOpacity, ActivityIndicator, useWindowDimensions } from 'react-native';
import Svg, { Path, Line, Text as SvgText, Circle } from 'react-native-svg';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';
import { supabase } from '../../lib/supabaseClient';

// ─── Triple-point and phase boundary approach ──────────────────────────────
// We use the Antoine-style Psat for vaporization, plus approximations for fusion.
// Compound data comes from compound_properties Supabase table.

interface CompoundInfo {
  name: string;
  Tc: number;   // K
  Pc: number;   // Pa
  Tb: number;   // K (normal boiling point)
  Tm: number;   // K (normal melting point) — from DB or default
  omega: number;
  vp_coeffs: number[];  // PLXANT-style: ln(P[Pa]) = C1 + C2/(T+C3) + C4*T + C5*ln(T) + C6*T^C7
}

function psatPa(c: CompoundInfo, T: number) {
  const [C1, C2, C3, C4, C5, C6, C7] = c.vp_coeffs;
  return Math.exp(C1 + C2 / (T + C3) + C4 * T + C5 * Math.log(T) + C6 * Math.pow(T, C7));
}

// Fusion curve: Clausius-Clapeyron simplified, P = Ptp + (ΔHfus/ΔVfus)*(T - Ttp)
// ΔHfus ≈ Trouton estimate, ΔVfus: negative slope approx
// For mobile we use the standard approximation: P_fusion scales very steeply with T
// Since we don't have exact fusion data, use Waals-Platteeuw-like approximation

const SVG_W = 300;
const SVG_H = 250;
const PAD = { l: 50, r: 16, t: 16, b: 36 };

export default function UnaryPhaseScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const svgWidth = Math.min(width - 32, 380);

  const [query, setQuery] = useState('');
  const [suggestions, setSuggestions] = useState<string[]>([]);
  const [compound, setCompound] = useState<CompoundInfo | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const debounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);

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
        .limit(1).single();
      if (err || !data) throw new Error(err?.message ?? 'Not found');
      const p = (data as any).properties;
      const Tc = Number(p?.CriticalTemperature?.value);
      const Pc = Number(p?.CriticalPressure?.value);
      const Tb = Number(p?.NormalBoilingPoint?.value);
      const TmRaw = Number(p?.NormalMeltingPoint?.value);
      const Tm = isFinite(TmRaw) && TmRaw > 0 ? TmRaw : (isFinite(Tb) && Tb > 0 ? Tb * 0.6 : 273);
      const omega = Number(p?.AcentricFactor?.value) || 0.3;
      const vp = p?.VaporPressure?.coefficients ?? [];
      if (!isFinite(Tc) || Tc <= 0 || !isFinite(Pc) || Pc <= 0 || !isFinite(Tb) || Tb <= 0 || vp.length === 0)
        throw new Error('Insufficient data for this compound');
      setCompound({ name: data.name, Tc, Pc, Tb, Tm, omega, vp_coeffs: vp });
    } catch (e: any) {
      setError(e.message);
    } finally {
      setLoading(false);
    }
  }

  const svgData = useMemo(() => {
    if (!compound) return null;
    const { Tc, Pc, Tb, Tm } = compound;

    // Vapor-liquid curve: from Tm to Tc
    const Ttp = Tm * 0.98; // approximate triple point temperature
    const Ptp = psatPa(compound, Ttp);
    const vapPoints: [number, number][] = [];
    const N = 80;
    for (let i = 0; i <= N; i++) {
      const T = Ttp + (Tc - Ttp) * i / N;
      try {
        const P = psatPa(compound, T);
        if (isFinite(P) && P > 0) vapPoints.push([T, P]);
      } catch {}
    }

    // Sublimation curve: below Ttp, approximate using Clausius-Clapeyron
    // ln(P/Ptp) = -(ΔHsub/R) * (1/T - 1/Ttp)
    // ΔHsub ≈ 1.3 * ΔHvap(Tb)  Trouton: ΔHvap ≈ 10.5*R*Tb
    const dHvap = 10.5 * 8.314 * Tb; // J/mol
    const dHsub = 1.3 * dHvap;
    const subPoints: [number, number][] = [];
    const Nub = 40;
    for (let i = 0; i <= Nub; i++) {
      const T = Ttp * (1 - 0.5 * i / Nub);
      const P = Ptp * Math.exp(-(dHsub / 8.314) * (1 / T - 1 / Ttp));
      if (isFinite(P) && P > 0) subPoints.push([T, P]);
    }

    // Fusion curve: nearly vertical, slight positive slope
    const fusPoints: [number, number][] = [[Ttp, Ptp], [Tm * 1.15, Ptp * 5000]];

    // Determine plot ranges
    const allT = [...vapPoints, ...subPoints, ...fusPoints].map(p => p[0]);
    const allP = [...vapPoints, ...subPoints].map(p => p[1]);
    const Tmin = Math.min(...allT), Tmax = Tc * 1.05;
    const Pmin = Math.max(1, Math.min(...allP) * 0.01), Pmax = Pc * 1.2;

    const plotW = SVG_W - PAD.l - PAD.r;
    const plotH = SVG_H - PAD.t - PAD.b;

    function txf(T: number) { return PAD.l + (T - Tmin) / (Tmax - Tmin) * plotW; }
    function tyf(P: number) {
      const logP = Math.log10(P), logPmin = Math.log10(Pmin), logPmax = Math.log10(Pmax);
      return PAD.t + plotH - (logP - logPmin) / (logPmax - logPmin) * plotH;
    }

    function toPath(pts: [number, number][]) {
      if (pts.length === 0) return '';
      return pts.map((p, i) => `${i === 0 ? 'M' : 'L'} ${txf(p[0]).toFixed(1)} ${tyf(p[1]).toFixed(1)}`).join(' ');
    }

    const yTicks: number[] = [];
    for (let e = Math.floor(Math.log10(Pmin)); e <= Math.ceil(Math.log10(Pmax)); e++) {
      yTicks.push(Math.pow(10, e));
    }
    const xTicks = [Tmin, Ttp, Tb, Tc].filter(t => t >= Tmin && t <= Tmax);

    return {
      vapPath: toPath(vapPoints),
      subPath: toPath(subPoints),
      fusPts: fusPoints.map(p => [txf(p[0]), tyf(p[1])] as [number, number]),
      critPt: [txf(Tc), tyf(Pc)] as [number, number],
      tripPt: [txf(Ttp), tyf(Ptp)] as [number, number],
      yTicks: yTicks.map(P => ({ P, y: tyf(P) })),
      xTicks: xTicks.map(T => ({ T, x: txf(T), label: T.toFixed(0) })),
    };
  }, [compound]);

  const axisColor = colors.border;
  const textColor = colors.mutedForeground;

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

      {compound && svgData && (
        <SectionCard title={`${compound.name} — Phase Diagram`}>
          <Svg width={svgWidth} height={SVG_H * (svgWidth / SVG_W)} viewBox={`0 0 ${SVG_W} ${SVG_H}`}>
            {/* Y-axis ticks */}
            {svgData.yTicks.map(t => (
              <React.Fragment key={t.P}>
                <Line x1={PAD.l} y1={t.y} x2={SVG_W - PAD.r} y2={t.y} stroke={axisColor} strokeWidth={0.5} />
                <SvgText x={PAD.l - 3} y={t.y + 3} fontSize={7} fill={textColor} textAnchor="end">
                  {t.P >= 1e5 ? `${(t.P / 1e5).toFixed(0)}b` : t.P >= 1000 ? `${(t.P / 1000).toFixed(0)}k` : t.P.toFixed(0)}
                </SvgText>
              </React.Fragment>
            ))}
            {/* X-axis ticks */}
            {svgData.xTicks.map(t => (
              <React.Fragment key={t.T}>
                <Line x1={t.x} y1={PAD.t} x2={t.x} y2={SVG_H - PAD.b} stroke={axisColor} strokeWidth={0.5} />
                <SvgText x={t.x} y={SVG_H - PAD.b + 12} fontSize={7} fill={textColor} textAnchor="middle">{t.label}K</SvgText>
              </React.Fragment>
            ))}

            {/* Sublimation */}
            <Path d={svgData.subPath} stroke="#60a5fa" strokeWidth={2} fill="none" />
            {/* Fusion */}
            {svgData.fusPts.length >= 2 && (
              <Line x1={svgData.fusPts[0][0]} y1={svgData.fusPts[0][1]}
                x2={svgData.fusPts[1][0]} y2={svgData.fusPts[1][1]}
                stroke="#f97316" strokeWidth={2} />
            )}
            {/* Vaporization */}
            <Path d={svgData.vapPath} stroke="#22c55e" strokeWidth={2} fill="none" />

            {/* Critical point */}
            <Circle cx={svgData.critPt[0]} cy={svgData.critPt[1]} r={4} fill="#a855f7" />
            <SvgText x={svgData.critPt[0] - 4} y={svgData.critPt[1] - 6} fontSize={7} fill="#a855f7">Tc/Pc</SvgText>

            {/* Triple point */}
            <Circle cx={svgData.tripPt[0]} cy={svgData.tripPt[1]} r={3} fill="#facc15" />
            <SvgText x={svgData.tripPt[0] + 4} y={svgData.tripPt[1] - 4} fontSize={7} fill="#facc15">TP</SvgText>

            {/* Axis labels */}
            <SvgText x={SVG_W / 2} y={SVG_H - 2} fontSize={8} fill={textColor} textAnchor="middle">Temperature (K)</SvgText>
            <SvgText x={10} y={SVG_H / 2} fontSize={8} fill={textColor} textAnchor="middle" transform={`rotate(-90,10,${SVG_H / 2})`}>Pressure (log)</SvgText>
          </Svg>
          <View style={styles.legend}>
            {[
              { color: '#60a5fa', label: 'Sublimation' },
              { color: '#f97316', label: 'Fusion' },
              { color: '#22c55e', label: 'Vaporization' },
              { color: '#a855f7', label: 'Critical point' },
              { color: '#facc15', label: 'Triple point' },
            ].map(l => (
              <View key={l.label} style={styles.legendItem}>
                <View style={[styles.dot, { backgroundColor: l.color }]} />
                <Text style={[styles.legendText, { color: colors.mutedForeground }]}>{l.label}</Text>
              </View>
            ))}
          </View>
          <Text style={[styles.note, { color: colors.mutedForeground }]}>
            Tc = {compound.Tc.toFixed(1)} K · Pc = {(compound.Pc / 1e5).toFixed(2)} bar · Tb = {compound.Tb.toFixed(1)} K
          </Text>
        </SectionCard>
      )}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  input: { borderWidth: 1, borderRadius: 10, padding: 12, fontSize: 14 },
  dropdown: { borderWidth: 1, borderRadius: 10, marginTop: 4 },
  suggestion: { paddingVertical: 10, paddingHorizontal: 14, borderBottomWidth: StyleSheet.hairlineWidth },
  errorText: { color: '#ef4444', fontSize: 13, marginTop: 6 },
  legend: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginTop: 10 },
  legendItem: { flexDirection: 'row', alignItems: 'center', gap: 4 },
  dot: { width: 8, height: 8, borderRadius: 4 },
  legendText: { fontSize: 11 },
  note: { fontSize: 11, marginTop: 8, fontStyle: 'italic' },
});
