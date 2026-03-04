import React from 'react';
import { View, Text, StyleSheet, TouchableOpacity, Linking } from 'react-native';
import { useTheme } from '../theme/ThemeContext';
import type { WebViewScreenProps } from '../navigation/types';

// Web platform: react-native-webview has no web implementation.
// Show a simple "open in browser" card instead.
export default function WebViewScreen({ route }: WebViewScreenProps) {
  const { url, title } = route.params;
  const { colors } = useTheme();

  return (
    <View style={[styles.container, { backgroundColor: colors.background }]}>
      <Text style={[styles.title, { color: colors.foreground }]}>
        {title ?? 'Open in Browser'}
      </Text>
      <Text style={[styles.sub, { color: colors.mutedForeground }]} numberOfLines={2}>
        {url}
      </Text>
      <TouchableOpacity
        style={[styles.btn, { backgroundColor: colors.blue }]}
        onPress={() => Linking.openURL(url)}
      >
        <Text style={styles.btnText}>Open in Browser</Text>
      </TouchableOpacity>
    </View>
  );
}

const styles = StyleSheet.create({
  container: {
    flex: 1,
    justifyContent: 'center',
    alignItems: 'center',
    padding: 32,
  },
  title: {
    fontSize: 20,
    fontWeight: '700',
    marginBottom: 8,
    textAlign: 'center',
  },
  sub: {
    fontSize: 12,
    marginBottom: 28,
    textAlign: 'center',
  },
  btn: {
    paddingHorizontal: 28,
    paddingVertical: 12,
    borderRadius: 10,
  },
  btnText: {
    color: '#ffffff',
    fontSize: 15,
    fontWeight: '600',
  },
});
