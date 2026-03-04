import React, { useRef, useState } from 'react';
import { View, StyleSheet, ActivityIndicator, TouchableOpacity, Text } from 'react-native';
import { WebView } from 'react-native-webview';
import { useTheme } from '../theme/ThemeContext';
import type { WebViewScreenProps } from '../navigation/types';

export default function WebViewScreen({ route }: WebViewScreenProps) {
  const { url } = route.params;
  const { colors } = useTheme();
  const webViewRef = useRef<WebView>(null);
  const [loading, setLoading] = useState(true);
  const [errorOccurred, setErrorOccurred] = useState(false);

  return (
    <View style={[styles.container, { backgroundColor: colors.background }]}>
      <WebView
        ref={webViewRef}
        source={{ uri: url }}
        style={{ flex: 1, backgroundColor: colors.background }}
        onLoadStart={() => { setLoading(true); setErrorOccurred(false); }}
        onLoadEnd={() => setLoading(false)}
        onError={() => { setLoading(false); setErrorOccurred(true); }}
        // Allow cookies and JS so the Next.js app works properly
        javaScriptEnabled
        domStorageEnabled
        sharedCookiesEnabled
        // Pass dark/light preference via header userAgent isn't ideal;
        // the site will use its own stored preference or system default
        renderLoading={() => (
          <View style={[styles.loading, { backgroundColor: colors.background }]}>
            <ActivityIndicator size="large" color={colors.blue} />
          </View>
        )}
        startInLoadingState
      />

      {errorOccurred && (
        <View style={[styles.errorOverlay, { backgroundColor: colors.background }]}>
          <Text style={[styles.errorTitle, { color: colors.foreground }]}>
            Couldn't load page
          </Text>
          <Text style={[styles.errorSub, { color: colors.mutedForeground }]}>
            {url}
          </Text>
          <TouchableOpacity
            style={[styles.retryBtn, { backgroundColor: colors.blue }]}
            onPress={() => {
              setErrorOccurred(false);
              webViewRef.current?.reload();
            }}
          >
            <Text style={styles.retryText}>Try Again</Text>
          </TouchableOpacity>
        </View>
      )}

      {loading && !errorOccurred && (
        <View
          pointerEvents="none"
          style={[styles.loadingOverlay, { backgroundColor: colors.background }]}
        >
          <ActivityIndicator size="large" color={colors.blue} />
        </View>
      )}
    </View>
  );
}

const styles = StyleSheet.create({
  container: {
    flex: 1,
  },
  loading: {
    flex: 1,
    justifyContent: 'center',
    alignItems: 'center',
  },
  loadingOverlay: {
    ...StyleSheet.absoluteFillObject,
    justifyContent: 'center',
    alignItems: 'center',
  },
  errorOverlay: {
    ...StyleSheet.absoluteFillObject,
    justifyContent: 'center',
    alignItems: 'center',
    padding: 24,
  },
  errorTitle: {
    fontSize: 18,
    fontWeight: '700',
    marginBottom: 8,
  },
  errorSub: {
    fontSize: 12,
    marginBottom: 24,
    textAlign: 'center',
  },
  retryBtn: {
    paddingHorizontal: 28,
    paddingVertical: 12,
    borderRadius: 10,
  },
  retryText: {
    color: '#ffffff',
    fontSize: 15,
    fontWeight: '600',
  },
});
