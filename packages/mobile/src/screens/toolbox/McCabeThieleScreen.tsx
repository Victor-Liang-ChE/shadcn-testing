import React, { useState, useMemo } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import Svg, { Path, Line, Circle, Text as SvgText, G } from 'react-native-svg';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ─── Equilibrium curve: constant relative volatility ─────────────────────────
function yEq(x: number, alpha: number) {
  return (alpha * x) / (1 + (alpha - 1) * x);
}

function xFromY(y: number, alpha: number) {
  return y / (alpha - (alpha - 1) * y);
}

// ─── Operating line intersections ────────────────────────────────────────────
function calcLines(xF: number, xD: number, xB: number, q: number, R: number) {
  // Rectifying: y = (R/(R+1))*x + xD/(R+1)
  const mRect = R / (R + 1);
  const bRect = xD / (R + 1);
  // q-line: y = (q/(q-1))*x - xF/(q-1)   (q != 1)
  // Intersection x_int
  let xInt: number, yInt: number;
  if (Math.abs(q - 1) < 0.0001) {
    xInt = xF;
    yInt = mRect * xF + bRect;
  } else {
    const mQ = q / (q - 1);
    const bQ = -xF / (q - 1);
    // R: y = mRect*x + bRect; Q: y = mQ*x + bQ
    xInt = (bQ - bRect) / (mRect - mQ);
    yInt = mRect * xInt + bRect;
  }
  // Stripping: through (xB, xB) and (xInt, yInt)
  const mStrip = (yInt - xB) / (xInt - xB || 1e-9);
  const bStrip = xB - mStrip * xB;
  return { mRect, bRect, mStrip, bStrip, xInt, yInt };
}

// ─── Stage stepping ──────────────────────────────────────────────────────────
function countStages(
  alpha: number, xD: number, xB: number,
  mRect: number, bRect: number,
  mStrip: number, bStrip: number,
  xInt: number,
): { stages: number; steps: { x1: number; y1: number; x2: number; y2: number }[] } {
  const steps: { x1: number; y1: number; x2: number; y2: number }[] = [];
  let x = xD;
  let y = xD;
  let stages = 0;

  for (let i = 0; i < 100; i++) {
    // Step horizontally to equilibrium curve
    const xNext = xFromY(y, alpha);
    steps.push({ x1: x, y1: y, x2: xNext, y2: y }); // horizontal
    x = xNext;
    stages++;

    if (x <= xB + 1e-6) break;

    // Switch operating line if needed
    const useStrip = x < xInt;
    const yNext = useStrip ? mStrip * x + bStrip : mRect * x + bRect;
    steps.push({ x1: x, y1: y, x2: x, y2: yNext }); // vertical
    y = yNext;
  }

  return { stages, steps };
}

// ─── SVG coordinate transform ────────────────────────────────────────────────
const SV = 300; // SVG viewBox size
const PAD = 28;
const PLOT = SV - 2 * PAD;

function tx(x: number) { return PAD + x * PLOT; }
function ty(y: number) { return SV - PAD - y * PLOT; }

export default function McCabeThieleScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const svgSize = Math.min(width - 32, 360);

  const [alpha, setAlpha] = useState(2.5);
  const [xF, setXF] = useState(0.4);
  const [xD, setXD] = useState(0.9);
  const [xB, setXB] = useState(0.1);
  const [q, setQ] = useState(1.0);
  const [R, setR] = useState(2.0);

  const { mRect, bRect, mStrip, bStrip, xInt, yInt, stages, steps, eqPath, diagPath, rectPath, stripPath, feedX, feedY } = useMemo(() => {
    const { mRect, bRect, mStrip, bStrip, xInt, yInt } = calcLines(xF, xD, xB, q, R);
    const { stages, steps } = countStages(alpha, xD, xB, mRect, bRect, mStrip, bStrip, xInt);

    // Generate equilibrium curve path
    const N = 60;
    let eqPath = `M ${tx(0)} ${ty(0)}`;
    for (let i = 1; i <= N; i++) {
      const x = i / N;
      const y = yEq(x, alpha);
      eqPath += ` L ${tx(x)} ${ty(y)}`;
    }

    // Diagonal y=x
    const diagPath = `M ${tx(0)} ${ty(0)} L ${tx(1)} ${ty(1)}`;

    // Rectifying operating line: xInt → xD
    const rectPath = `M ${tx(xInt)} ${ty(yInt)} L ${tx(xD)} ${ty(xD)}`;

    // Stripping: xB → xInt
    const stripPath = `M ${tx(xB)} ${ty(xB)} L ${tx(xInt)} ${ty(yInt)}`;

    // Feed line (q-line) from xF to xInt
    const feedX: [number, number] = [tx(xF), tx(xInt)];
    const feedY2: [number, number] = [ty(xF), ty(yInt)];

    return { mRect, bRect, mStrip, bStrip, xInt, yInt, stages, steps, eqPath, diagPath, rectPath, stripPath, feedX, feedY: feedY2 };
  }, [alpha, xF, xD, xB, q, R]);

  const isDark = colors.background === '#09090b' || colors.background.startsWith('#0');
  const axisColor = isDark ? '#444' : '#ddd';
  const textColor = isDark ? '#aaa' : '#666';

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <SectionCard title="Parameters">
        <SliderRow label="Relative Volatility α" value={alpha} min={1.1} max={10} step={0.1} decimals={1} onChange={setAlpha} />
        <SliderRow label="Feed xF" value={xF} min={0.05} max={0.95} step={0.01} decimals={2} onChange={v => setXF(Math.min(v, xD - 0.05))} />
        <SliderRow label="Distillate xD" value={xD} min={0.5} max={0.99} step={0.01} decimals={2} onChange={v => setXD(Math.max(v, xF + 0.05))} />
        <SliderRow label="Bottoms xB" value={xB} min={0.01} max={0.49} step={0.01} decimals={2} onChange={v => setXB(Math.min(v, xF - 0.05))} />
        <SliderRow label="q (feed quality)" value={q} min={0} max={2} step={0.05} decimals={2} onChange={setQ} />
        <SliderRow label="Reflux Ratio R" value={R} min={0.1} max={20} step={0.1} decimals={1} onChange={setR} />
      </SectionCard>

      <SectionCard title="McCabe-Thiele Diagram">
        <Svg width={svgSize} height={svgSize} viewBox={`0 0 ${SV} ${SV}`}>
          {/* Grid */}
          {[0, 0.25, 0.5, 0.75, 1].map(v => (
            <React.Fragment key={v}>
              <Line x1={tx(v)} y1={ty(0)} x2={tx(v)} y2={ty(1)} stroke={axisColor} strokeWidth={0.5} />
              <Line x1={tx(0)} y1={ty(v)} x2={tx(1)} y2={ty(v)} stroke={axisColor} strokeWidth={0.5} />
              <SvgText x={tx(v)} y={ty(0) + 14} fontSize={8} fill={textColor} textAnchor="middle">{v.toFixed(2)}</SvgText>
              <SvgText x={PAD - 4} y={ty(v) + 3} fontSize={8} fill={textColor} textAnchor="end">{v.toFixed(2)}</SvgText>
            </React.Fragment>
          ))}

          {/* Diagonal */}
          <Path d={diagPath} stroke={axisColor} strokeWidth={1} fill="none" />

          {/* Equilibrium curve */}
          <Path d={eqPath} stroke="#60a5fa" strokeWidth={2} fill="none" />

          {/* q-line */}
          <Line x1={feedX[0]} y1={feedY[0]} x2={feedX[1]} y2={feedY[1]} stroke="#facc15" strokeWidth={1.5} strokeDasharray="5,3" />

          {/* Rectifying OL */}
          <Path d={rectPath} stroke="#22c55e" strokeWidth={1.8} fill="none" />

          {/* Stripping OL */}
          <Path d={stripPath} stroke="#f97316" strokeWidth={1.8} fill="none" />

          {/* Stages */}
          {steps.map((s, i) => (
            <Line key={i} x1={tx(s.x1)} y1={ty(s.y1)} x2={tx(s.x2)} y2={ty(s.y2)}
              stroke="#a855f7" strokeWidth={1.5} />
          ))}

          {/* Key points */}
          {[{ x: xD, y: xD, c: '#22c55e', l: 'xD' }, { x: xB, y: xB, c: '#f97316', l: 'xB' }, { x: xF, y: xF, c: '#facc15', l: 'xF' }].map(pt => (
            <React.Fragment key={pt.l}>
              <Circle cx={tx(pt.x)} cy={ty(pt.y)} r={4} fill={pt.c} />
              <SvgText x={tx(pt.x) + 6} y={ty(pt.y) - 3} fontSize={8} fill={pt.c}>{pt.l}</SvgText>
            </React.Fragment>
          ))}

          {/* Axis labels */}
          <SvgText x={SV / 2} y={SV - 4} fontSize={9} fill={textColor} textAnchor="middle">x (liquid mole fraction)</SvgText>
          <SvgText x={8} y={SV / 2} fontSize={9} fill={textColor} textAnchor="middle" transform={`rotate(-90,8,${SV / 2})`}>y (vapor)</SvgText>
        </Svg>
      </SectionCard>

      <SectionCard title="Results">
        <View style={styles.resultGrid}>
          <View style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.label, { color: colors.mutedForeground }]}>Theoretical Stages</Text>
            <Text style={[styles.value, { color: colors.blue }]}>{stages}</Text>
          </View>
          <View style={[styles.resultItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.label, { color: colors.mutedForeground }]}>Reflux Ratio R</Text>
            <Text style={[styles.value, { color: colors.foreground }]}>{R.toFixed(1)}</Text>
          </View>
        </View>
        <View style={styles.legend}>
          {[
            { color: '#60a5fa', label: 'Equilibrium curve' },
            { color: '#22c55e', label: 'Rectifying line' },
            { color: '#f97316', label: 'Stripping line' },
            { color: '#facc15', label: 'q-line (feed)' },
            { color: '#a855f7', label: 'Stage steps' },
          ].map(l => (
            <View key={l.label} style={styles.legendItem}>
              <View style={[styles.legendDot, { backgroundColor: l.color }]} />
              <Text style={[styles.legendText, { color: colors.mutedForeground }]}>{l.label}</Text>
            </View>
          ))}
        </View>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  resultGrid: { flexDirection: 'row', gap: 10, marginBottom: 12 },
  resultItem: { flex: 1, borderRadius: 10, borderWidth: 1, padding: 12, alignItems: 'center' },
  label: { fontSize: 12, marginBottom: 4 },
  value: { fontSize: 28, fontWeight: '700' },
  legend: { flexDirection: 'row', flexWrap: 'wrap', gap: 10 },
  legendItem: { flexDirection: 'row', alignItems: 'center', gap: 5 },
  legendDot: { width: 10, height: 10, borderRadius: 5 },
  legendText: { fontSize: 12 },
});
