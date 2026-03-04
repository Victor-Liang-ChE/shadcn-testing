import React, { useState, useMemo, useCallback, useRef } from 'react';
import { ScrollView, View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { LineChart } from 'react-native-gifted-charts';
import SliderRow from '../../components/SliderRow';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

// ─── Tuning methods ──────────────────────────────────────────────────────────
type Method = 'ZN-OL' | 'ITAE-D' | 'AMIGO' | 'IMC';

interface PIDGains { Kc: number; tauI: number; tauD: number }

function zieglerNicholsOL(Kp: number, tau: number, theta: number): PIDGains {
  // Open-loop (reaction curve) Ziegler-Nichols
  return {
    Kc: 1.2 * tau / (Kp * theta),
    tauI: 2 * theta,
    tauD: 0.5 * theta,
  };
}

function itaeDisturbance(Kp: number, tau: number, theta: number): PIDGains {
  // ITAE for disturbance rejection
  const r = theta / tau;
  return {
    Kc: (0.859 * Math.pow(r, -0.977)) / Kp,
    tauI: tau / (0.674 * Math.pow(r, -0.680)),
    tauD: 0.381 * tau * Math.pow(r, 0.995),
  };
}

function amigo(Kp: number, tau: number, theta: number): PIDGains {
  const r = theta / tau;
  return {
    Kc: (0.2 + 0.45 * (tau / theta)) / Kp,
    tauI: (0.4 * theta + 0.8 * tau) / (theta + 0.1 * tau) * theta,
    tauD: 0.5 * theta * tau / (0.3 * theta + tau),
  };
}

function imc(Kp: number, tau: number, theta: number, lambda: number): PIDGains {
  return {
    Kc: tau / (Kp * (lambda + theta)),
    tauI: tau,
    tauD: theta / 2,
  };
}

// ─── FOPTD Step-response simulator (Euler) ──────────────────────────────────
const SIM_DT = 0.05;
const SIM_N = 600; // 30 s at dt=0.05

function simulatePID(
  Kp: number, tau: number, theta: number,
  gains: PIDGains,
  setpoint: number,
): { t: number[]; pv: number[]; sp: number[] } {
  const { Kc, tauI, tauD } = gains;
  const nDelay = Math.max(1, Math.round(theta / SIM_DT));
  const buffer: number[] = new Array(nDelay).fill(0);
  let pv = 0, integral = 0, prevErr = 0, bufIdx = 0;
  const ts: number[] = [], pvs: number[] = [], sps: number[] = [];

  for (let i = 0; i < SIM_N; i++) {
    const t = i * SIM_DT;
    const sp = t >= 2 ? setpoint : 0;
    const err = sp - pv;
    integral += err * SIM_DT;
    const deriv = (err - prevErr) / SIM_DT;
    const u = Kc * (err + integral / (tauI || 1e-9) + tauD * deriv);

    // Apply u to FOPTD plant via delay buffer
    const uDelayed = buffer[bufIdx];
    buffer[bufIdx] = u;
    bufIdx = (bufIdx + 1) % nDelay;

    const dpv = (dt: number) => (-pv + Kp * uDelayed) / tau * dt;
    pv += dpv(SIM_DT);
    prevErr = err;

    if (i % 4 === 0) { ts.push(t); pvs.push(pv); sps.push(sp); }
  }
  return { t: ts, pv: pvs, sp: sps };
}

const METHODS: { key: Method; label: string; desc: string }[] = [
  { key: 'ZN-OL', label: 'Ziegler-Nichols', desc: 'Reaction-curve method' },
  { key: 'ITAE-D', label: 'ITAE', desc: 'Optimized for disturbance rejection' },
  { key: 'AMIGO', label: 'AMIGO', desc: 'Approximately optimal M-constrained' },
  { key: 'IMC', label: 'IMC', desc: 'Internal Model Control (λ tunable)' },
];

export default function PidTuningScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [method, setMethod] = useState<Method>('IMC');
  const [Kp, setKp] = useState(2);
  const [tau, setTau] = useState(10);
  const [theta, setTheta] = useState(2);
  const [lambda, setLambda] = useState(3);
  const [setpoint, setSetpoint] = useState(1);
  const CHART_W = width - 64;

  const gains = useMemo((): PIDGains => {
    switch (method) {
      case 'ZN-OL': return zieglerNicholsOL(Kp, tau, theta);
      case 'ITAE-D': return itaeDisturbance(Kp, tau, theta);
      case 'AMIGO': return amigo(Kp, tau, theta);
      case 'IMC': return imc(Kp, tau, theta, lambda);
    }
  }, [method, Kp, tau, theta, lambda]);

  const simResult = useMemo(() => simulatePID(Kp, tau, theta, gains, setpoint), [Kp, tau, theta, gains, setpoint]);

  const pvData = simResult.pv.map((v, i) => ({
    value: v,
    label: i % Math.floor(simResult.t.length / 5) === 0 ? simResult.t[i].toFixed(0) : '',
  }));
  const spData = simResult.sp.map(v => ({ value: v }));

  // Performance
  const perf = useMemo(() => {
    const { pv, t, sp } = simResult;
    const ss = sp[sp.length - 1];
    if (ss === 0) return null;
    const peak = Math.max(...pv);
    const OS = ((peak - ss) / ss) * 100;
    let idx95 = pv.findIndex((v, i) => i > 0 && sp[i] > 0 && Math.abs(v - ss) / ss < 0.05);
    const ts = idx95 >= 0 ? t[idx95] : null;
    return { OS, ts };
  }, [simResult]);

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      {/* Method selector */}
      {METHODS.map(m => (
        <TouchableOpacity key={m.key} onPress={() => setMethod(m.key)}
          style={[styles.methodCard, { backgroundColor: method === m.key ? colors.blue : colors.card, borderColor: method === m.key ? colors.blue : colors.border }]}>
          <Text style={[styles.methodLabel, { color: method === m.key ? '#fff' : colors.foreground }]}>{m.label}</Text>
          <Text style={[styles.methodDesc, { color: method === m.key ? 'rgba(255,255,255,0.75)' : colors.mutedForeground }]}>{m.desc}</Text>
        </TouchableOpacity>
      ))}

      <SectionCard title="Process Model (FOPTD)">
        <SliderRow label="Process Gain Kp" value={Kp} min={0.1} max={10} step={0.1} decimals={1} onChange={setKp} />
        <SliderRow label="Time Constant τ (s)" value={tau} min={1} max={50} step={1} decimals={0} onChange={setTau} />
        <SliderRow label="Dead Time θ (s)" value={theta} min={0.1} max={20} step={0.1} decimals={1} onChange={setTheta} />
        {method === 'IMC' && (
          <SliderRow label="IMC Filter λ (s)" value={lambda} min={0.1} max={30} step={0.5} decimals={1} onChange={setLambda} />
        )}
        <SliderRow label="Setpoint" value={setpoint} min={0.1} max={5} step={0.1} decimals={1} onChange={setSetpoint} />
      </SectionCard>

      <SectionCard title="PID Parameters">
        <View style={styles.gainsRow}>
          {[
            { label: 'Kc', value: gains.Kc },
            { label: 'τI (s)', value: gains.tauI },
            { label: 'τD (s)', value: gains.tauD },
          ].map(g => (
            <View key={g.label} style={[styles.gainItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
              <Text style={[styles.gainLabel, { color: colors.mutedForeground }]}>{g.label}</Text>
              <Text style={[styles.gainValue, { color: colors.blue }]}>{g.value.toFixed(3)}</Text>
            </View>
          ))}
        </View>
      </SectionCard>

      <SectionCard title="Closed-Loop Step Response">
        <LineChart
          data={pvData}
          width={CHART_W}
          height={220}
          color={colors.blue}
          thickness={2.5}
          secondaryData={spData}
          secondaryLineConfig={{ color: '#ef4444', thickness: 1.5 }}
          xAxisColor={colors.border}
          yAxisColor={colors.border}
          xAxisLabelTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          yAxisTextStyle={{ color: colors.mutedForeground, fontSize: 10 }}
          dataPointsRadius={0}
          hideDataPoints
          initialSpacing={0}
          endSpacing={4}
          noOfSections={4}
        />
        <View style={styles.legendRow}>
          <View style={[styles.legendDot, { backgroundColor: colors.blue }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>PV (Process Variable)</Text>
          <View style={[styles.legendDot, { backgroundColor: '#ef4444' }]} />
          <Text style={[styles.legendText, { color: colors.mutedForeground }]}>SP (Setpoint)</Text>
        </View>
      </SectionCard>

      {perf && (
        <SectionCard title="Performance">
          <Text style={[styles.perfText, { color: colors.foreground }]}>
            Overshoot: {perf.OS.toFixed(1)}%
          </Text>
          {perf.ts !== null && (
            <Text style={[styles.perfText, { color: colors.foreground }]}>
              Settling Time (5%): {perf.ts.toFixed(1)} s
            </Text>
          )}
        </SectionCard>
      )}
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40, gap: 0 },
  methodCard: { borderWidth: 1.5, borderRadius: 12, padding: 12, marginBottom: 8, flexDirection: 'row', alignItems: 'center', justifyContent: 'space-between' },
  methodLabel: { fontSize: 15, fontWeight: '700' },
  methodDesc: { fontSize: 12 },
  gainsRow: { flexDirection: 'row', gap: 8 },
  gainItem: { flex: 1, borderWidth: 1, borderRadius: 10, padding: 10, alignItems: 'center' },
  gainLabel: { fontSize: 12, marginBottom: 4 },
  gainValue: { fontSize: 18, fontWeight: '700' },
  legendRow: { flexDirection: 'row', alignItems: 'center', gap: 6, marginTop: 8, justifyContent: 'center' },
  legendDot: { width: 10, height: 10, borderRadius: 5 },
  legendText: { fontSize: 12 },
  perfText: { fontSize: 14, marginBottom: 4 },
});
