import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import Svg, { Path, Circle, Line, Text as SvgText, Polygon } from 'react-native-svg';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ─── Ternary coordinate transform ────────────────────────────────────────────
// Equilateral triangle: vertices at (0,0), (1,0), (0.5, sqrt(3)/2)
const SQ3_2 = Math.sqrt(3) / 2;
const SV = 320; // SVG size
const PAD = 24;
const SIDE = SV - 2 * PAD;

function ternaryToXY(a: number, b: number, c: number): [number, number] {
  // a=bottom-left, b=bottom-right, c=top
  const norm = a + b + c || 1;
  const an = a / norm, bn = b / norm, cn = c / norm;
  const x = PAD + (bn + cn / 2) * SIDE;
  const y = PAD + (1 - cn * SQ3_2) * SIDE;
  return [x, y];
}

// ─── Modified Raoult + Antoine for 3 components ──────────────────────────────
interface Antoine { A: number; B: number; C: number }

function psat(a: Antoine, T: number) {
  return Math.pow(10, a.A - a.B / (T + a.C)); // mmHg
}

const SYSTEMS: {
  label: string; comp: [string, string, string];
  ant: [Antoine, Antoine, Antoine];
  alpha: [number, number, number]; // relative to comp 2 (middle boiler)
}[] = [
  {
    label: 'Acetone / Methanol / Water',
    comp: ['Acetone', 'MeOH', 'Water'],
    ant: [
      { A: 7.11714, B: 1210.595, C: 229.664 },
      { A: 7.87863, B: 1473.11, C: 230.0 },
      { A: 8.07131, B: 1730.63, C: 233.426 },
    ],
    alpha: [3.2, 1.0, 0.45],
  },
  {
    label: 'Benzene / Toluene / o-Xylene',
    comp: ['Benzene', 'Toluene', 'o-Xylene'],
    ant: [
      { A: 6.87987, B: 1196.760, C: 219.161 },
      { A: 6.95087, B: 1342.310, C: 219.187 },
      { A: 6.99891, B: 1474.679, C: 213.686 },
    ],
    alpha: [4.0, 1.0, 0.26],
  },
  {
    label: 'Ethanol / Water / EtAc',
    comp: ['Ethanol', 'Water', 'EtAc'],
    ant: [
      { A: 8.11220, B: 1592.864, C: 226.184 },
      { A: 8.07131, B: 1730.63, C: 233.426 },
      { A: 7.10179, B: 1244.951, C: 217.881 },
    ],
    alpha: [1.2, 0.8, 2.8],
  },
];

// Bubble temperature using relative volatility (simplified)
function bubbleT_3(x: [number, number, number], P: number, ant: [Antoine, Antoine, Antoine]): number {
  let lo = 40, hi = 120;
  for (let i = 0; i < 60; i++) {
    const T = (lo + hi) / 2;
    const sum = x[0] * psat(ant[0], T) + x[1] * psat(ant[1], T) + x[2] * psat(ant[2], T);
    if (sum > P) hi = T; else lo = T;
  }
  return (lo + hi) / 2;
}

// Vapor compositions at bubble T
function vapComps(x: [number, number, number], T: number, P: number, ant: [Antoine, Antoine, Antoine]): [number, number, number] {
  const y0 = x[0] * psat(ant[0], T) / P;
  const y1 = x[1] * psat(ant[1], T) / P;
  const y2 = x[2] * psat(ant[2], T) / P;
  const sum = y0 + y1 + y2 || 1;
  return [y0 / sum, y1 / sum, y2 / sum];
}

// Residue curve ODE: dx/dt = x - y
function residueCurve(
  x0: [number, number, number],
  P: number,
  ant: [Antoine, Antoine, Antoine],
  nSteps = 80,
  dt = 0.15,
): [number, number, number][] {
  let x = [...x0] as [number, number, number];
  const pts: [number, number, number][] = [[...x] as [number, number, number]];
  for (let i = 0; i < nSteps; i++) {
    const T = bubbleT_3(x, P, ant);
    const y = vapComps(x, T, P, ant);
    const dx0 = x[0] - y[0];
    const dx1 = x[1] - y[1];
    x = [
      Math.max(0.001, Math.min(0.998, x[0] + dx0 * dt)),
      Math.max(0.001, Math.min(0.998, x[1] + dx1 * dt)),
      0,
    ];
    x[2] = Math.max(0.001, 1 - x[0] - x[1]);
    pts.push([...x] as [number, number, number]);
  }
  return pts;
}

// Generate grid of starting points
function generateStartPoints(n: number): [number, number, number][] {
  const pts: [number, number, number][] = [];
  const step = 1 / (n + 1);
  for (let i = 1; i <= n; i++) {
    for (let j = 1; j <= n - i + 1; j++) {
      const a = i * step, b = j * step;
      if (a + b < 1 - step / 2) {
        pts.push([a, b, 1 - a - b]);
      }
    }
  }
  return pts;
}

const CURVE_COLORS = ['#60a5fa', '#22c55e', '#f97316', '#a855f7', '#ec4899', '#facc15', '#06b6d4'];

export default function ResidueCurveMapScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const svgSize = Math.min(width - 32, 360);

  const [sysIdx, setSysIdx] = useState(0);
  const [P_mmHg, setP] = useState(760);
  const [nCurves, setNCurves] = useState(6);

  const sys = SYSTEMS[sysIdx];
  const isDark = colors.background.startsWith('#0') || colors.background === '#09090b';

  const { curves, azeoPts } = useMemo(() => {
    const starts = generateStartPoints(nCurves);
    const curves = starts.map(s => residueCurve(s, P_mmHg, sys.ant, 100, 0.12));
    // Detect azeotropes: points where |x - y| is small across all three
    // For now just note if diagram shows crossing — skip complex detection
    return { curves, azeoPts: [] as [number, number, number][] };
  }, [sys, P_mmHg, nCurves]);

  function curveToPath(pts: [number, number, number][]): string {
    return pts.map((p, i) => {
      const [cx, cy] = ternaryToXY(p[0], p[2], p[1]);
      return `${i === 0 ? 'M' : 'L'} ${cx.toFixed(1)} ${cy.toFixed(1)}`;
    }).join(' ');
  }

  const triPts = [
    ...ternaryToXY(1, 0, 0),
    ...ternaryToXY(0, 1, 0),
    ...ternaryToXY(0, 0, 1),
  ];
  const triStr = `${triPts[0]},${triPts[1]} ${triPts[2]},${triPts[3]} ${triPts[4]},${triPts[5]}`;
  const v0 = ternaryToXY(1, 0, 0);
  const v1 = ternaryToXY(0, 1, 0);
  const v2 = ternaryToXY(0, 0, 1);

  const axisColor = isDark ? '#444' : '#ccc';
  const textColor = isDark ? '#aaa' : '#666';

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <View style={styles.sysRow}>
        {SYSTEMS.map((s, i) => (
          <TouchableOpacity key={s.label} onPress={() => setSysIdx(i)}
            style={[styles.chip, { backgroundColor: sysIdx === i ? colors.blue : colors.card, borderColor: sysIdx === i ? colors.blue : colors.border }]}>
            <Text style={[styles.chipText, { color: sysIdx === i ? '#fff' : colors.foreground }]}>{s.label}</Text>
          </TouchableOpacity>
        ))}
      </View>

      <SectionCard title="Conditions">
        <SliderRow label="Pressure (mmHg)" value={P_mmHg} min={200} max={3000} step={50} decimals={0} onChange={setP} />
        <SliderRow label="Number of Curves" value={nCurves} min={3} max={8} step={1} decimals={0} onChange={setNCurves} />
      </SectionCard>

      <SectionCard title="Residue Curve Map">
        <Svg width={svgSize} height={svgSize * 0.9} viewBox={`0 0 ${SV} ${SV * 0.9}`}>
          {/* Triangle outline */}
          <Polygon points={triStr} fill="none" stroke={axisColor} strokeWidth={1.5} />

          {/* Residue curves */}
          {curves.map((pts, i) => (
            <Path key={i} d={curveToPath(pts)}
              stroke={CURVE_COLORS[i % CURVE_COLORS.length]}
              strokeWidth={1.5} fill="none" opacity={0.85} />
          ))}

          {/* Vertex labels */}
          <SvgText x={v0[0] - 4} y={v0[1] + 14} fontSize={9} fill={textColor} textAnchor="middle">{sys.comp[0]}</SvgText>
          <SvgText x={v1[0] + 4} y={v1[1] + 14} fontSize={9} fill={textColor} textAnchor="middle">{sys.comp[1]}</SvgText>
          <SvgText x={v2[0]} y={v2[1] - 6} fontSize={9} fill={textColor} textAnchor="middle">{sys.comp[2]}</SvgText>

          {/* Grid lines at 0.25, 0.5, 0.75 */}
          {[0.25, 0.5, 0.75].map(f => {
            const pts = [
              ternaryToXY(f, 1 - f, 0),
              ternaryToXY(f, 0, 1 - f),
              ternaryToXY(0, f, 1 - f),
              ternaryToXY(1 - f, f, 0),
              ternaryToXY(0, 1 - f, f),
              ternaryToXY(1 - f, 0, f),
            ];
            return (
              <React.Fragment key={f}>
                <Line x1={pts[0][0]} y1={pts[0][1]} x2={pts[1][0]} y2={pts[1][1]} stroke={axisColor} strokeWidth={0.4} />
                <Line x1={pts[2][0]} y1={pts[2][1]} x2={pts[3][0]} y2={pts[3][1]} stroke={axisColor} strokeWidth={0.4} />
                <Line x1={pts[4][0]} y1={pts[4][1]} x2={pts[5][0]} y2={pts[5][1]} stroke={axisColor} strokeWidth={0.4} />
              </React.Fragment>
            );
          })}
        </Svg>
        <Text style={[styles.note, { color: colors.mutedForeground }]}>
          Ideal Raoult's Law · {P_mmHg} mmHg · Arrows follow residue curve direction (increasing T)
        </Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  sysRow: { flexDirection: 'row', flexWrap: 'wrap', gap: 8, marginBottom: 14 },
  chip: { paddingHorizontal: 12, paddingVertical: 8, borderRadius: 20, borderWidth: 1 },
  chipText: { fontSize: 12, fontWeight: '600' },
  note: { fontSize: 11, textAlign: 'center', marginTop: 6, fontStyle: 'italic' },
});
