import React, { useState, useEffect, useCallback } from 'react';
import {
  ScrollView, View, Text, StyleSheet, TouchableOpacity, ActivityIndicator,
  Modal, useWindowDimensions, Image,
} from 'react-native';
import { createClient } from '@supabase/supabase-js';
import SectionCard from '../../components/SectionCard';
import { useTheme } from '../../theme/ThemeContext';

const supabase = createClient(
  process.env.EXPO_PUBLIC_SUPABASE_URL!,
  process.env.EXPO_PUBLIC_SUPABASE_ANON_KEY!
);

interface Champion {
  fighter: string;
  weight_class: string;
  start_date: string;
  end_date: string | null;
  days?: number;
}

interface FighterStats {
  height?: string;
  weight?: string;
  reach?: string;
  nationality?: string;
}

function getDays(start: string, end: string | null): number {
  const s = new Date(start).getTime();
  const e = end ? new Date(end).getTime() : Date.now();
  return Math.floor((e - s) / 86400000);
}

function getReignColor(days: number): string {
  const t = Math.min(days / 2550, 1);
  const h = Math.round(200 + t * 40);
  const l = Math.round(30 + t * 30);
  return `hsl(${h}, 70%, ${l}%)`;
}

const GENDER_MAP: Record<string, 'M' | 'F'> = {
  "Strawweight Women's": 'F', "Flyweight Women's": 'F',
  "Bantamweight Women's": 'F', "Featherweight Women's": 'F',
};

export default function UfcChampionsScreen() {
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [champions, setChampions] = useState<Champion[]>([]);
  const [loading, setLoading] = useState(true);
  const [genderFilter, setGenderFilter] = useState<'all' | 'men' | 'women'>('all');
  const [selected, setSelected] = useState<Champion | null>(null);
  const [stats, setStats] = useState<FighterStats | null>(null);
  const [statsLoading, setStatsLoading] = useState(false);
  const [weightClass, setWeightClass] = useState<string>('all');

  useEffect(() => {
    (async () => {
      const { data } = await supabase
        .from('ufc_champion_history')
        .select('fighter, weight_class, start_date, end_date')
        .order('weight_class')
        .order('start_date');
      if (data) {
        setChampions(data.map(c => ({ ...c, days: getDays(c.start_date, c.end_date) })));
      }
      setLoading(false);
    })();
  }, []);

  const weightClasses = ['all', ...Array.from(new Set(champions.map(c => c.weight_class)))];

  const filtered = champions.filter(c => {
    const isFemale = !!GENDER_MAP[c.weight_class];
    if (genderFilter === 'men' && isFemale) return false;
    if (genderFilter === 'women' && !isFemale) return false;
    if (weightClass !== 'all' && c.weight_class !== weightClass) return false;
    return true;
  });

  // Group by weight class
  const groups: Record<string, Champion[]> = {};
  for (const c of filtered) {
    if (!groups[c.weight_class]) groups[c.weight_class] = [];
    groups[c.weight_class].push(c);
  }

  const openFighter = useCallback(async (c: Champion) => {
    setSelected(c);
    setStats(null);
    setStatsLoading(true);
    try {
      const { data } = await supabase
        .from('ufc_fighter_tott')
        .select('height, weight, reach, nationality')
        .ilike('name', c.fighter)
        .single();
      setStats(data || null);
    } catch { setStats(null); }
    setStatsLoading(false);
  }, []);

  return (
    <View style={{ flex: 1, backgroundColor: colors.background }}>
      {/* Gender filter */}
      <View style={[styles.filterBar, { backgroundColor: colors.card, borderBottomColor: colors.border }]}>
        {(['all', 'men', 'women'] as const).map(g => (
          <TouchableOpacity key={g} onPress={() => setGenderFilter(g)}
            style={[styles.filterChip, genderFilter === g && { backgroundColor: colors.blue }]}>
            <Text style={[styles.filterText, { color: genderFilter === g ? '#fff' : colors.foreground }]}>
              {g.charAt(0).toUpperCase() + g.slice(1)}
            </Text>
          </TouchableOpacity>
        ))}
      </View>

      {loading ? (
        <View style={styles.center}><ActivityIndicator size="large" color={colors.blue} /></View>
      ) : (
        <ScrollView contentContainerStyle={styles.container}>
          {Object.entries(groups).map(([wc, champs]) => (
            <SectionCard key={wc} title={wc}>
              {champs.map((c, i) => {
                const days = c.days || 0;
                const barWidth = Math.min((days / 2550) * (width - 80), width - 80);
                return (
                  <TouchableOpacity key={i} onPress={() => openFighter(c)}
                    style={[styles.champRow, { borderBottomColor: colors.border }]}>
                    <View style={[styles.reignBar, { width: barWidth, backgroundColor: getReignColor(days) }]} />
                    <View style={styles.champInfo}>
                      <Text style={[styles.champName, { color: colors.foreground }]}>{c.fighter}</Text>
                      <Text style={[styles.champDates, { color: colors.mutedForeground }]}>
                        {c.start_date.slice(0, 10)} — {c.end_date ? c.end_date.slice(0, 10) : 'current'} ({days} days)
                      </Text>
                    </View>
                  </TouchableOpacity>
                );
              })}
            </SectionCard>
          ))}
        </ScrollView>
      )}

      {/* Fighter detail modal */}
      <Modal visible={!!selected} transparent animationType="slide" onRequestClose={() => setSelected(null)}>
        <View style={styles.modalOverlay}>
          <View style={[styles.modalCard, { backgroundColor: colors.card, borderColor: colors.border }]}>
            <TouchableOpacity style={styles.closeBtn} onPress={() => setSelected(null)}>
              <Text style={{ color: colors.mutedForeground, fontSize: 18 }}>✕</Text>
            </TouchableOpacity>
            {selected && (
              <>
                <Text style={[styles.modalName, { color: colors.foreground }]}>{selected.fighter}</Text>
                <Text style={[styles.modalWC, { color: colors.blue }]}>{selected.weight_class}</Text>
                <Text style={[styles.modalDates, { color: colors.mutedForeground }]}>
                  {selected.start_date.slice(0, 10)} — {selected.end_date ? selected.end_date.slice(0, 10) : 'Ongoing'}
                  {' '}({selected.days} days)
                </Text>
                {statsLoading && <ActivityIndicator color={colors.blue} style={{ marginTop: 10 }} />}
                {!statsLoading && stats && (
                  <View style={styles.statsGrid}>
                    {[
                      { label: 'Height', val: stats.height },
                      { label: 'Weight', val: stats.weight },
                      { label: 'Reach', val: stats.reach },
                      { label: 'Nationality', val: stats.nationality },
                    ].map(s => s.val ? (
                      <View key={s.label} style={styles.statItem}>
                        <Text style={[styles.statLabel, { color: colors.mutedForeground }]}>{s.label}</Text>
                        <Text style={[styles.statVal, { color: colors.foreground }]}>{s.val}</Text>
                      </View>
                    ) : null)}
                  </View>
                )}
              </>
            )}
          </View>
        </View>
      </Modal>
    </View>
  );
}

const styles = StyleSheet.create({
  filterBar: { flexDirection: 'row', padding: 10, borderBottomWidth: 1, gap: 8 },
  filterChip: { paddingHorizontal: 14, paddingVertical: 6, borderRadius: 16 },
  filterText: { fontSize: 13, fontWeight: '600' },
  center: { flex: 1, justifyContent: 'center', alignItems: 'center' },
  container: { padding: 12, paddingBottom: 40 },
  champRow: { paddingVertical: 10, borderBottomWidth: StyleSheet.hairlineWidth, minHeight: 50 },
  reignBar: { height: 4, borderRadius: 2, marginBottom: 6, opacity: 0.7 },
  champInfo: {},
  champName: { fontSize: 14, fontWeight: '600', marginBottom: 2 },
  champDates: { fontSize: 11 },
  modalOverlay: { flex: 1, backgroundColor: 'rgba(0,0,0,0.5)', justifyContent: 'flex-end' },
  modalCard: { borderTopLeftRadius: 20, borderTopRightRadius: 20, borderWidth: 1, padding: 24, paddingBottom: 40 },
  closeBtn: { position: 'absolute', top: 16, right: 20 },
  modalName: { fontSize: 22, fontWeight: '800', marginBottom: 4 },
  modalWC: { fontSize: 15, fontWeight: '600', marginBottom: 4 },
  modalDates: { fontSize: 13, marginBottom: 16 },
  statsGrid: { flexDirection: 'row', flexWrap: 'wrap', gap: 16 },
  statItem: {},
  statLabel: { fontSize: 11 },
  statVal: { fontSize: 14, fontWeight: '600' },
});
