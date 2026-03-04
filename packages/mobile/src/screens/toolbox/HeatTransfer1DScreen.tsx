import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions, TextInput } from 'react-native';
import Svg, { Rect, Line, Text as SvgText, Defs, LinearGradient, Stop } from 'react-native-svg';
import { useTheme } from '../../theme/ThemeContext';
import SectionCard from '../../components/SectionCard';
import SliderRow from '../../components/SliderRow';

const LAYER_COLORS = ['#60a5fa', '#f97316', '#22c55e', '#a855f7', '#ec4899'];

interface Layer {
  id: number;
  k: number;
  L: number;
  label: string;
}

let nextId = 4;

export default function HeatTransfer1DScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [Thot, setThot] = useState(200);
  const [Tcold, setTcold] = useState(20);
  const [area, setArea] = useState(1);
  const [layers, setLayers] = useState<Layer[]>([
    { id: 1, k: 50, L: 0.05, label: 'Steel' },
    { id: 2, k: 0.04, L: 0.1, label: 'Insulation' },
    { id: 3, k: 1.5, L: 0.03, label: 'Concrete' },
  ]);

  const SVG_W = width - 48;
  const SVG_H = 120;

  const results = useMemo(() => {
    const Rs = layers.map(l => l.L / (l.k * area));
    const Rtotal = Rs.reduce((a, b) => a + b, 0);
    const Q = Rtotal > 0 ? (Thot - Tcold) / Rtotal : 0;
    const temps: number[] = [];
    let T = Thot;
    temps.push(T);
    for (const r of Rs) {
      T -= Q * r;
      temps.push(T);
    }
    return { Rs, Rtotal, Q, temps };
  }, [layers, Thot, Tcold, area]);

  const totalL = layers.reduce((s, l) => s + l.L, 0) || 1;

  function addLayer() {
    if (layers.length >= 5) return;
    setLayers(prev => [...prev, { id: nextId++, k: 1, L: 0.05, label: `Layer ${prev.length + 1}` }]);
  }

  function removeLayer(id: number) {
    setLayers(prev => prev.filter(l => l.id !== id));
  }

  function updateLayer(id: number, field: keyof Layer, val: number | string) {
    setLayers(prev => prev.map(l => l.id === id ? { ...l, [field]: val } : l));
  }

  // SVG positions
  const PAD_L = 10, PAD_R = 10, BAR_H = 70, BAR_Y = 25;
  const barW = SVG_W - PAD_L - PAD_R;

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <SectionCard title="Boundary Conditions">
        <SliderRow label="T hot (°C)" value={Thot} min={50} max={500} step={10} decimals={0} onChange={setThot} />
        <SliderRow label="T cold (°C)" value={Tcold} min={-20} max={100} step={5} decimals={0} onChange={setTcold} />
        <SliderRow label="Area (m²)" value={area} min={0.1} max={10} step={0.1} decimals={1} onChange={setArea} />
      </SectionCard>

      <SectionCard title="Wall Layers">
        {layers.map((layer, idx) => (
          <View key={layer.id} style={[styles.layerBox, { borderColor: LAYER_COLORS[idx % LAYER_COLORS.length], backgroundColor: colors.card }]}>
            <View style={styles.layerHeader}>
              <View style={[styles.colorDot, { backgroundColor: LAYER_COLORS[idx % LAYER_COLORS.length] }]} />
              <TextInput
                value={layer.label}
                onChangeText={v => updateLayer(layer.id, 'label', v)}
                style={[styles.labelInput, { color: colors.foreground, borderColor: colors.border }]}
              />
              <TouchableOpacity onPress={() => removeLayer(layer.id)}>
                <Text style={styles.removeBtn}>✕</Text>
              </TouchableOpacity>
            </View>
            <SliderRow label="Conductivity k (W/m·K)" value={layer.k} min={0.01} max={400} step={0.01} decimals={2} onChange={v => updateLayer(layer.id, 'k', v)} />
            <SliderRow label="Thickness L (m)" value={layer.L} min={0.001} max={1} step={0.001} decimals={3} onChange={v => updateLayer(layer.id, 'L', v)} />
            <Text style={[styles.rText, { color: colors.mutedForeground }]}>
              R = {results.Rs[idx]?.toFixed(4)} K/W
            </Text>
          </View>
        ))}
        {layers.length < 5 && (
          <TouchableOpacity onPress={addLayer}
            style={[styles.addBtn, { backgroundColor: colors.blue, borderRadius: 10 }]}>
            <Text style={{ color: '#fff', fontWeight: '700', textAlign: 'center' }}>+ Add Layer</Text>
          </TouchableOpacity>
        )}
      </SectionCard>

      <SectionCard title="Wall Cross-Section">
        <Svg width={SVG_W} height={SVG_H}>
          <Defs>
            <LinearGradient id="grad" x1="0" y1="0" x2="1" y2="0">
              <Stop offset="0" stopColor="#ef4444" stopOpacity="0.7" />
              <Stop offset="1" stopColor="#60a5fa" stopOpacity="0.7" />
            </LinearGradient>
          </Defs>
          {layers.map((layer, idx) => {
            const x = PAD_L + layers.slice(0, idx).reduce((s, l) => s + (l.L / totalL) * barW, 0);
            const w = (layer.L / totalL) * barW;
            const color = LAYER_COLORS[idx % LAYER_COLORS.length];
            return (
              <React.Fragment key={layer.id}>
                <Rect x={x} y={BAR_Y} width={w} height={BAR_H} fill={color} opacity={0.55} />
                {w > 28 && (
                  <SvgText x={x + w / 2} y={BAR_Y + BAR_H / 2 + 5} fontSize={10} fill="#fff" textAnchor="middle">{layer.label}</SvgText>
                )}
              </React.Fragment>
            );
          })}
          {/* Temperature labels */}
          <SvgText x={PAD_L} y={20} fontSize={11} fill="#ef4444" textAnchor="middle">{Thot}°C</SvgText>
          <SvgText x={SVG_W - PAD_R} y={20} fontSize={11} fill="#60a5fa" textAnchor="end">{Tcold}°C</SvgText>
        </Svg>
      </SectionCard>

      <SectionCard title="Results">
        <View style={styles.resultGrid}>
          <View style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.resultLabel, { color: colors.mutedForeground }]}>Total Resistance</Text>
            <Text style={[styles.resultValue, { color: colors.foreground }]}>{results.Rtotal.toFixed(4)} K/W</Text>
          </View>
          <View style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.resultLabel, { color: colors.mutedForeground }]}>Heat Flux Q</Text>
            <Text style={[styles.resultValue, { color: colors.blue }]}>{results.Q.toFixed(2)} W</Text>
          </View>
        </View>
        <Text style={[styles.tempTitle, { color: colors.foreground }]}>Interface Temperatures:</Text>
        {results.temps.map((T, i) => (
          <Text key={i} style={[styles.tempRow, { color: colors.mutedForeground }]}>
            {i === 0 ? 'Hot surface' : i === layers.length ? 'Cold surface' : `Interface ${i}`}: {T.toFixed(1)} °C
          </Text>
        ))}
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  layerBox: { borderWidth: 1.5, borderRadius: 10, padding: 10, marginBottom: 10 },
  layerHeader: { flexDirection: 'row', alignItems: 'center', marginBottom: 6, gap: 8 },
  colorDot: { width: 12, height: 12, borderRadius: 6 },
  labelInput: { flex: 1, fontSize: 14, fontWeight: '600', borderBottomWidth: 1, paddingBottom: 2 },
  removeBtn: { color: '#ef4444', fontSize: 18, fontWeight: '700', paddingHorizontal: 4 },
  rText: { fontSize: 12, marginTop: 2 },
  addBtn: { padding: 12, marginTop: 4 },
  resultGrid: { flexDirection: 'row', gap: 10, marginBottom: 12 },
  resultItem: { flex: 1, borderRadius: 10, borderWidth: 1, padding: 12, alignItems: 'center' },
  resultLabel: { fontSize: 12, marginBottom: 4 },
  resultValue: { fontSize: 18, fontWeight: '700' },
  tempTitle: { fontSize: 13, fontWeight: '600', marginBottom: 4 },
  tempRow: { fontSize: 13, marginBottom: 2 },
});
