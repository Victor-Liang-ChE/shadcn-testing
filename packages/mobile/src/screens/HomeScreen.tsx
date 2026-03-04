import React, { useState, useMemo } from 'react';
import {
  View,
  Text,
  ScrollView,
  StyleSheet,
  TextInput,
  TouchableOpacity,
  Image,
  useWindowDimensions,
  Linking,
} from 'react-native';
import { useNavigation } from '@react-navigation/native';
import { NativeStackNavigationProp } from '@react-navigation/native-stack';
import { Ionicons } from '@expo/vector-icons';
import { useTheme } from '../theme/ThemeContext';
import { FEATURED_SIMULATIONS, ToolItem } from '../constants/data';
import type { RootStackParamList } from '../navigation/types';

type NavProp = NativeStackNavigationProp<RootStackParamList>;

function SimCard({ item, cardWidth }: { item: ToolItem; cardWidth: number }) {
  const navigation = useNavigation<NavProp>();
  const { colors, colorScheme } = useTheme();
  const thumbUri = colorScheme === 'light' ? item.thumbnailLightPath : item.thumbnailPath;

  return (
    <TouchableOpacity
      style={[
        styles.simCard,
        {
          backgroundColor: colors.card,
          borderColor: colors.border,
          width: cardWidth,
        },
      ]}
      activeOpacity={0.75}
      onPress={() =>
        item.screenKey === 'WebView'
          ? navigation.navigate('WebView', { url: `https://victorliang.com${item.path}`, title: item.name })
          : navigation.navigate(item.screenKey as any)
      }
    >
      <View style={styles.thumbContainer}>
        <Image
          source={{ uri: thumbUri }}
          style={styles.thumbnail}
          resizeMode="contain"
        />
      </View>
      <Text style={[styles.simName, { color: colors.foreground }]}>{item.name}</Text>
      <Text
        style={[styles.simDesc, { color: colors.mutedForeground }]}
        numberOfLines={3}
      >
        {item.description}
      </Text>
    </TouchableOpacity>
  );
}

export default function HomeScreen() {
  const navigation = useNavigation<NavProp>();
  const { colors } = useTheme();
  const { width } = useWindowDimensions();
  const [searchQuery, setSearchQuery] = useState('');

  const PADDING = 12;
  const GAP = 10;
  const cardWidth = (width - PADDING * 2 - GAP) / 2;

  const filteredSims = useMemo(() => {
    if (!searchQuery.trim()) return FEATURED_SIMULATIONS;
    const q = searchQuery.toLowerCase();
    return FEATURED_SIMULATIONS.filter(
      s =>
        s.name.toLowerCase().includes(q) ||
        s.description.toLowerCase().includes(q),
    );
  }, [searchQuery]);

  // Build rows of 2
  const rows: ToolItem[][] = [];
  for (let i = 0; i < filteredSims.length; i += 2) {
    rows.push(filteredSims.slice(i, i + 2));
  }

  return (
    <ScrollView
      style={{ backgroundColor: colors.background }}
      contentContainerStyle={[styles.container, { paddingHorizontal: PADDING }]}
      keyboardShouldPersistTaps="handled"
    >
      {/* ── BIO CARD ── */}
      <View style={[styles.bioCard, { backgroundColor: colors.card, borderColor: colors.border }]}>
        <Text style={[styles.bioTitle, { color: colors.foreground }]}>
          Hi, I'm{' '}
          <Text style={{ color: colors.blue }}>Victor</Text>
        </Text>
        <Text style={[styles.bioSubtitle, { color: colors.foreground }]}>
          Chemical Engineering Optimization{'&'} Predictive Modeling Enthusiast
        </Text>

        <View style={[styles.divider, { backgroundColor: colors.border }]} />

        <BioRow label="Software" color={colors} text="SolidWorks, Aspen HYSYS, Aspen Plus, MATLAB, VS Code 💻" />
        <BioRow label="Python" color={colors} text="NumPy, SciPy, Pandas, Matplotlib, Scikit-learn, Plotly 📦" />
        <BioRow label="Location" color={colors} text="Based in SF — MS Materials Science @ UCSB" />
        <BioRow label="Degree" color={colors} text="BS Chemical Engineering" />

        <Text style={[styles.bioText, { color: colors.mutedForeground }]}>
          📧{' '}
          <Text
            style={{ color: colors.blue }}
            onPress={() => Linking.openURL('mailto:victorl1725@gmail.com')}
          >
            victorl1725@gmail.com
          </Text>
        </Text>
        <Text style={[styles.bioText, { color: colors.mutedForeground }]}>
          Under heavy construction, but take a look around! 😄
        </Text>

        <View style={styles.socialRow}>
          <TouchableOpacity
            style={[styles.socialBtn, { borderColor: colors.border, backgroundColor: colors.background }]}
            onPress={() =>
              navigation.navigate('WebView', {
                url: 'https://github.com/Victor-Liang-ChE',
                title: 'GitHub',
              })
            }
          >
            <Ionicons name="logo-github" size={18} color={colors.foreground} />
            <Text style={[styles.socialBtnText, { color: colors.foreground }]}>GitHub</Text>
          </TouchableOpacity>
          <TouchableOpacity
            style={[styles.socialBtn, { borderColor: colors.border, backgroundColor: colors.background }]}
            onPress={() =>
              navigation.navigate('WebView', {
                url: 'https://www.linkedin.com/in/victor-liang-567238231/',
                title: 'LinkedIn',
              })
            }
          >
            <Ionicons name="logo-linkedin" size={18} color="#0a66c2" />
            <Text style={[styles.socialBtnText, { color: colors.foreground }]}>LinkedIn</Text>
          </TouchableOpacity>
        </View>
      </View>

      {/* ── FEATURED SIMULATIONS HEADER ── */}
      <Text style={[styles.sectionHeader, { color: colors.foreground }]}>
        Featured Simulations and Tools
      </Text>

      {/* ── SEARCH BAR ── */}
      <View style={[styles.searchBar, { backgroundColor: colors.searchBackground, borderColor: colors.border }]}>
        <Ionicons name="search" size={16} color={colors.mutedForeground} style={{ marginRight: 8 }} />
        <TextInput
          style={[styles.searchInput, { color: colors.foreground }]}
          placeholder="Search tools..."
          placeholderTextColor={colors.mutedForeground}
          value={searchQuery}
          onChangeText={setSearchQuery}
          returnKeyType="search"
          clearButtonMode="while-editing"
        />
      </View>

      {/* ── SIMULATION GRID ── */}
      {rows.map((row, idx) => (
        <View key={idx} style={[styles.row, { gap: GAP }]}>
          {row.map(item => (
            <SimCard key={item.path} item={item} cardWidth={cardWidth} />
          ))}
          {/* Spacer if odd row */}
          {row.length === 1 && <View style={{ width: cardWidth }} />}
        </View>
      ))}

      {filteredSims.length === 0 && (
        <Text style={[styles.emptyText, { color: colors.mutedForeground }]}>
          No tools match "{searchQuery}"
        </Text>
      )}
    </ScrollView>
  );
}

function BioRow({
  label,
  text,
  color,
}: {
  label: string;
  text: string;
  color: ReturnType<typeof useTheme>['colors'];
}) {
  return (
    <Text style={[styles.bioText, { color: color.mutedForeground }]}>
      <Text style={{ fontWeight: '600', color: color.foreground }}>{label}: </Text>
      {text}
    </Text>
  );
}

const styles = StyleSheet.create({
  container: {
    paddingTop: 16,
    paddingBottom: 32,
  },
  bioCard: {
    borderRadius: 16,
    borderWidth: 1,
    padding: 20,
    marginBottom: 20,
  },
  bioTitle: {
    fontSize: 22,
    fontWeight: '700',
    marginBottom: 6,
  },
  bioSubtitle: {
    fontSize: 16,
    fontWeight: '700',
    marginBottom: 12,
    lineHeight: 22,
  },
  divider: {
    height: 1,
    marginBottom: 12,
  },
  bioText: {
    fontSize: 13,
    lineHeight: 20,
    marginBottom: 6,
  },
  socialRow: {
    flexDirection: 'row',
    gap: 10,
    marginTop: 12,
  },
  socialBtn: {
    flexDirection: 'row',
    alignItems: 'center',
    gap: 6,
    borderWidth: 1,
    borderRadius: 8,
    paddingHorizontal: 14,
    paddingVertical: 8,
  },
  socialBtnText: {
    fontSize: 14,
    fontWeight: '500',
  },
  sectionHeader: {
    fontSize: 20,
    fontWeight: '700',
    textAlign: 'center',
    marginBottom: 12,
  },
  searchBar: {
    flexDirection: 'row',
    alignItems: 'center',
    borderWidth: 1,
    borderRadius: 10,
    paddingHorizontal: 12,
    paddingVertical: 8,
    marginBottom: 14,
  },
  searchInput: {
    flex: 1,
    fontSize: 14,
    padding: 0,
  },
  row: {
    flexDirection: 'row',
    marginBottom: 10,
  },
  simCard: {
    borderRadius: 14,
    borderWidth: 1,
    padding: 12,
    overflow: 'hidden',
  },
  thumbContainer: {
    width: '100%',
    aspectRatio: 1,
    borderRadius: 14,
    overflow: 'hidden',
    marginBottom: 10,
    backgroundColor: 'transparent',
  },
  thumbnail: {
    width: '100%',
    height: '100%',
  },
  simName: {
    fontSize: 14,
    fontWeight: '700',
    marginBottom: 4,
  },
  simDesc: {
    fontSize: 11,
    lineHeight: 16,
  },
  emptyText: {
    textAlign: 'center',
    marginTop: 32,
    fontSize: 15,
  },
});
