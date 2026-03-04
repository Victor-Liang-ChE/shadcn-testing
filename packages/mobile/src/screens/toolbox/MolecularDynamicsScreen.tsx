import React, { useState, useRef, useCallback, useEffect } from 'react';
import { View, Text, StyleSheet, TouchableOpacity, useWindowDimensions } from 'react-native';
import { ScrollView } from 'react-native';
import Svg, { Circle as SvgCircle, Line as SvgLine } from 'react-native-svg';
import SectionCard from '../../components/SectionCard';
import SliderRow from '../../components/SliderRow';
import { useTheme } from '../../theme/ThemeContext';

// ─── 2D Lennard-Jones molecular dynamics ──────────────────────────────────────
const EPSILON = 1.0; // energy scale
const SIGMA = 1.0;   // length scale
const MASS = 1.0;

interface Particle {
  x: number; y: number;
  vx: number; vy: number;
  fx: number; fy: number;
}

function lj(r: number): { fMag: number; pe: number } {
  if (r < 0.7) return { fMag: 0, pe: 0 }; // hard cutoff
  const sr = SIGMA / r;
  const sr6 = Math.pow(sr, 6);
  const sr12 = sr6 * sr6;
  const pe = 4 * EPSILON * (sr12 - sr6);
  const fMag = 24 * EPSILON * (2 * sr12 - sr6) / (r * r);
  return { fMag, pe };
}

function initParticles(N: number, L: number): Particle[] {
  const particles: Particle[] = [];
  const cols = Math.ceil(Math.sqrt(N));
  const spacing = L / (cols + 1);
  for (let i = 0; i < N; i++) {
    const col = i % cols, row = Math.floor(i / cols);
    const x = spacing * (col + 1);
    const y = spacing * (row + 1);
    const theta = Math.random() * 2 * Math.PI;
    const speed = 0.5 + Math.random() * 0.5;
    particles.push({ x, y, vx: speed * Math.cos(theta), vy: speed * Math.sin(theta), fx: 0, fy: 0 });
  }
  return particles;
}

function computeForces(particles: Particle[], L: number): { pe: number } {
  let pe = 0;
  for (const p of particles) { p.fx = 0; p.fy = 0; }
  for (let i = 0; i < particles.length; i++) {
    for (let j = i + 1; j < particles.length; j++) {
      let dx = particles[i].x - particles[j].x;
      let dy = particles[i].y - particles[j].y;
      // Minimum image convention
      if (dx > L / 2) dx -= L;
      if (dx < -L / 2) dx += L;
      if (dy > L / 2) dy -= L;
      if (dy < -L / 2) dy += L;
      const r = Math.sqrt(dx * dx + dy * dy) || 1e-9;
      if (r > 3 * SIGMA) continue; // cutoff
      const { fMag, pe: pij } = lj(r);
      pe += pij;
      const fx = fMag * dx, fy = fMag * dy;
      particles[i].fx += fx; particles[i].fy += fy;
      particles[j].fx -= fx; particles[j].fy -= fy;
    }
  }
  return { pe };
}

function velVerlet(particles: Particle[], L: number, dt: number): { ke: number; pe: number } {
  // Update positions
  for (const p of particles) {
    p.x += p.vx * dt + 0.5 * p.fx / MASS * dt * dt;
    p.y += p.vy * dt + 0.5 * p.fy / MASS * dt * dt;
    // Periodic boundary
    p.x = ((p.x % L) + L) % L;
    p.y = ((p.y % L) + L) % L;
  }
  const oldFx = particles.map(p => p.fx);
  const oldFy = particles.map(p => p.fy);
  const { pe } = computeForces(particles, L);
  let ke = 0;
  for (let i = 0; i < particles.length; i++) {
    particles[i].vx += 0.5 * (oldFx[i] + particles[i].fx) / MASS * dt;
    particles[i].vy += 0.5 * (oldFy[i] + particles[i].fy) / MASS * dt;
    ke += 0.5 * MASS * (particles[i].vx ** 2 + particles[i].vy ** 2);
  }
  return { ke, pe };
}

// Velocity rescaling thermostat
function thermostat(particles: Particle[], targetT: number): void {
  const ke = particles.reduce((s, p) => s + 0.5 * MASS * (p.vx ** 2 + p.vy ** 2), 0);
  const Tcur = (2 * ke) / (2 * particles.length);
  if (Tcur < 1e-9) return;
  const scale = Math.sqrt(targetT / Tcur);
  for (const p of particles) { p.vx *= scale; p.vy *= scale; }
}

export default function MolecularDynamicsScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const svgSize = Math.min(width - 32, 360);
  const BOX = 12; // simulation box size in LJ units

  const [N, setN] = useState(16);
  const [T, setT] = useState(1.0);
  const [running, setRunning] = useState(false);
  const [particles, setParticles] = useState<Particle[]>(() => initParticles(16, BOX));
  const [energies, setEnergies] = useState<{ ke: number; pe: number; step: number }>({ ke: 0, pe: 0, step: 0 });
  const rafRef = useRef<number | null>(null);
  const particlesRef = useRef<Particle[]>(particles);
  const stepRef = useRef(0);
  const lastFrameRef = useRef<number>(0);
  const DT = 0.015;

  function reset() {
    setRunning(false);
    if (rafRef.current) cancelAnimationFrame(rafRef.current);
    const p = initParticles(N, BOX);
    computeForces(p, BOX);
    particlesRef.current = p;
    stepRef.current = 0;
    setParticles([...p]);
    setEnergies({ ke: 0, pe: 0, step: 0 });
  }

  useEffect(() => {
    reset();
  }, [N]);

  const animate = useCallback((now: number) => {
    if (now - lastFrameRef.current < 33) { // ~30fps
      rafRef.current = requestAnimationFrame(animate);
      return;
    }
    lastFrameRef.current = now;

    // 5 steps per frame
    let ke = 0, pe = 0;
    for (let i = 0; i < 5; i++) {
      const result = velVerlet(particlesRef.current, BOX, DT);
      ke = result.ke; pe = result.pe;
    }
    thermostat(particlesRef.current, T);
    stepRef.current += 5;
    setParticles([...particlesRef.current]);
    setEnergies({ ke, pe, step: stepRef.current });
    rafRef.current = requestAnimationFrame(animate);
  }, [T]);

  useEffect(() => {
    if (running) {
      rafRef.current = requestAnimationFrame(animate);
    } else {
      if (rafRef.current) cancelAnimationFrame(rafRef.current);
    }
    return () => { if (rafRef.current) cancelAnimationFrame(rafRef.current); };
  }, [running, animate]);

  const scale = svgSize / BOX;
  const r_px = Math.max(3, scale * 0.45);

  // Temperature indicator color
  const TColor = T < 0.5 ? '#60a5fa' : T < 1.5 ? '#22c55e' : T < 2.5 ? '#f97316' : '#ef4444';

  return (
    <ScrollView style={{ backgroundColor: colors.background }} contentContainerStyle={styles.container}>
      <SectionCard title="Simulation Box">
        <View style={[styles.boxContainer, { width: svgSize, height: svgSize, backgroundColor: colors.card, borderColor: colors.border }]}>
          <Svg width={svgSize} height={svgSize}>
            {particles.map((p, i) => (
              <SvgCircle
                key={i}
                cx={p.x * scale}
                cy={p.y * scale}
                r={r_px}
                fill={TColor}
                opacity={0.85}
              />
            ))}
          </Svg>
        </View>
        <Text style={[styles.stepLabel, { color: colors.mutedForeground }]}>Step: {energies.step}</Text>
      </SectionCard>

      <View style={styles.controls}>
        <TouchableOpacity
          onPress={() => setRunning(r => !r)}
          style={[styles.btn, { backgroundColor: running ? '#ef4444' : '#22c55e' }]}>
          <Text style={styles.btnText}>{running ? '⏸ Pause' : '▶ Run'}</Text>
        </TouchableOpacity>
        <TouchableOpacity onPress={reset}
          style={[styles.btn, { backgroundColor: colors.card, borderWidth: 1, borderColor: colors.border }]}>
          <Text style={[styles.btnText, { color: colors.foreground }]}>↺ Reset</Text>
        </TouchableOpacity>
      </View>

      <SectionCard title="Parameters">
        <SliderRow label="N Particles" value={N} min={4} max={36} step={1} decimals={0} onChange={v => { setN(v); }} />
        <SliderRow label="Target Temperature T*" value={T} min={0.1} max={4} step={0.1} decimals={1} onChange={setT} />
      </SectionCard>

      <SectionCard title="Energy">
        <View style={styles.energyRow}>
          <View style={[styles.energyItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.energyLabel, { color: colors.mutedForeground }]}>KE</Text>
            <Text style={[styles.energyValue, { color: '#22c55e' }]}>{energies.ke.toFixed(2)}</Text>
          </View>
          <View style={[styles.energyItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.energyLabel, { color: colors.mutedForeground }]}>PE</Text>
            <Text style={[styles.energyValue, { color: '#60a5fa' }]}>{energies.pe.toFixed(2)}</Text>
          </View>
          <View style={[styles.energyItem, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <Text style={[styles.energyLabel, { color: colors.mutedForeground }]}>Total</Text>
            <Text style={[styles.energyValue, { color: colors.foreground }]}>{(energies.ke + energies.pe).toFixed(2)}</Text>
          </View>
        </View>
      </SectionCard>

      <SectionCard title="About">
        <Text style={[styles.about, { color: colors.mutedForeground }]}>
          2D Lennard-Jones simulation with {N} particles in a periodic box.{'\n'}
          Forces: U(r) = 4ε[(σ/r)¹² - (σ/r)⁶]{'\n'}
          Integration: Velocity Verlet · Thermostat: Velocity rescaling{'\n'}
          Color scale: Blue (cold) → Orange → Red (hot)
        </Text>
      </SectionCard>
    </ScrollView>
  );
}

const styles = StyleSheet.create({
  container: { padding: 16, paddingBottom: 40 },
  boxContainer: { borderWidth: 1.5, borderRadius: 8, overflow: 'hidden' },
  stepLabel: { fontSize: 11, textAlign: 'center', marginTop: 4 },
  controls: { flexDirection: 'row', gap: 10, marginVertical: 12 },
  btn: { flex: 1, paddingVertical: 12, borderRadius: 10, alignItems: 'center' },
  btnText: { fontSize: 15, fontWeight: '700', color: '#fff' },
  energyRow: { flexDirection: 'row', gap: 8 },
  energyItem: { flex: 1, borderWidth: 1, borderRadius: 8, padding: 10, alignItems: 'center' },
  energyLabel: { fontSize: 11, marginBottom: 2 },
  energyValue: { fontSize: 18, fontWeight: '700' },
  about: { fontSize: 12, lineHeight: 20 },
});
