import React from 'react';
import { TouchableOpacity } from 'react-native';
import { createNativeStackNavigator } from '@react-navigation/native-stack';
import { createBottomTabNavigator } from '@react-navigation/bottom-tabs';
import { Ionicons } from '@expo/vector-icons';
import { useTheme } from '../theme/ThemeContext';
import HomeScreen from '../screens/HomeScreen';
import ToolboxScreen from '../screens/ToolboxScreen';
import MiscScreen from '../screens/MiscScreen';
import WebViewScreen from '../screens/WebViewScreen';
// Misc screens
import DropChanceScreen from '../screens/misc/DropChanceScreen';
import PortolaMenuScreen from '../screens/misc/PortolaMenuScreen';
import JapaneseLyricsScreen from '../screens/misc/JapaneseLyricsScreen';
import ChemistryToolsScreen from '../screens/misc/ChemistryToolsScreen';
import ChemeEconScreen from '../screens/misc/ChemeEconScreen';
import UfcChampionsScreen from '../screens/misc/UfcChampionsScreen';
import LatexConverterScreen from '../screens/misc/LatexConverterScreen';
import CsvToPlotScreen from '../screens/misc/CsvToPlotScreen';
// Toolbox screens
import IsothermsScreen from '../screens/toolbox/IsothermsScreen';
import ProcessDynamicsScreen from '../screens/toolbox/ProcessDynamicsScreen';
import HeatTransfer1DScreen from '../screens/toolbox/HeatTransfer1DScreen';
import KineticsScreen from '../screens/toolbox/KineticsScreen';
import CompoundPropertiesScreen from '../screens/toolbox/CompoundPropertiesScreen';
import PidTuningScreen from '../screens/toolbox/PidTuningScreen';
import McCabeThieleScreen from '../screens/toolbox/McCabeThieleScreen';
import AzeotropeFinderScreen from '../screens/toolbox/AzeotropeFinderScreen';
import BinaryPhaseScreen from '../screens/toolbox/BinaryPhaseScreen';
import UnaryPhaseScreen from '../screens/toolbox/UnaryPhaseScreen';
import ReactorDesignScreen from '../screens/toolbox/ReactorDesignScreen';
import FUGScreen from '../screens/toolbox/FUGScreen';
import ResidueCurveMapScreen from '../screens/toolbox/ResidueCurveMapScreen';
import MolecularDynamicsScreen from '../screens/toolbox/MolecularDynamicsScreen';
import type { RootStackParamList, TabParamList } from './types';

const Stack = createNativeStackNavigator<RootStackParamList>();
const Tab = createBottomTabNavigator<TabParamList>();

function TabNavigator() {
  const { colors, colorScheme, toggleTheme } = useTheme();

  const ThemeToggleButton = () => (
    <TouchableOpacity onPress={toggleTheme} style={{ marginRight: 14 }}>
      <Ionicons
        name={colorScheme === 'dark' ? 'sunny-outline' : 'moon-outline'}
        size={22}
        color={colors.foreground}
      />
    </TouchableOpacity>
  );

  return (
    <Tab.Navigator
      screenOptions={({ route }) => ({
        headerStyle: { backgroundColor: colors.navbarBackground },
        headerTintColor: colors.foreground,
        headerShadowVisible: false,
        headerRight: () => <ThemeToggleButton />,
        tabBarIcon: ({ focused, color, size }) => {
          let iconName: React.ComponentProps<typeof Ionicons>['name'];
          if (route.name === 'Home') {
            iconName = focused ? 'home' : 'home-outline';
          } else if (route.name === 'Toolbox') {
            iconName = focused ? 'construct' : 'construct-outline';
          } else {
            iconName = focused ? 'grid' : 'grid-outline';
          }
          return <Ionicons name={iconName} size={size} color={color} />;
        },
        tabBarActiveTintColor: colors.tabBarActive,
        tabBarInactiveTintColor: colors.tabBarInactive,
        tabBarStyle: {
          backgroundColor: colors.tabBar,
          borderTopColor: colors.border,
          borderTopWidth: 1,
        },
        tabBarLabelStyle: { fontSize: 12 },
      })}
    >
      <Tab.Screen
        name="Home"
        component={HomeScreen}
        options={{ title: 'Victor Liang' }}
      />
      <Tab.Screen
        name="Toolbox"
        component={ToolboxScreen}
        options={{ title: 'Toolbox' }}
      />
      <Tab.Screen
        name="Misc"
        component={MiscScreen}
        options={{ title: 'Misc' }}
      />
    </Tab.Navigator>
  );
}

export default function RootNavigator() {
  const { colors } = useTheme();

  return (
    <Stack.Navigator
      screenOptions={{
        headerStyle: { backgroundColor: colors.navbarBackground },
        headerTintColor: colors.foreground,
        headerShadowVisible: false,
        headerBackTitle: 'Back',
      }}
    >
      <Stack.Screen
        name="Main"
        component={TabNavigator}
        options={{ headerShown: false }}
      />
      <Stack.Screen
        name="WebView"
        component={WebViewScreen}
        options={({ route }) => ({ title: route.params.title })}
      />
      {/* Misc screens */}
      <Stack.Screen name="DropChance" component={DropChanceScreen} options={{ title: 'Drop Chance Calculator' }} />
      <Stack.Screen name="PortolaMenu" component={PortolaMenuScreen} options={{ title: 'Portola Menu' }} />
      <Stack.Screen name="JapaneseLyrics" component={JapaneseLyricsScreen} options={{ title: 'Japanese Lyrics Analyzer' }} />
      <Stack.Screen name="ChemistryTools" component={ChemistryToolsScreen} options={{ title: 'Chemistry Tools' }} />
      <Stack.Screen name="ChemeEcon" component={ChemeEconScreen} options={{ title: 'ChemE Economics' }} />
      <Stack.Screen name="UfcChampions" component={UfcChampionsScreen} options={{ title: 'UFC Championship Lineage' }} />
      <Stack.Screen name="LatexConverter" component={LatexConverterScreen} options={{ title: 'LaTeX Converter' }} />
      <Stack.Screen name="CsvToPlot" component={CsvToPlotScreen} options={{ title: 'CSV Plotter' }} />
      {/* Toolbox screens */}
      <Stack.Screen name="Isotherms" component={IsothermsScreen} options={{ title: 'Isotherms' }} />
      <Stack.Screen name="ProcessDynamics" component={ProcessDynamicsScreen} options={{ title: 'Process Dynamics' }} />
      <Stack.Screen name="HeatTransfer1D" component={HeatTransfer1DScreen} options={{ title: '1D Heat Transfer' }} />
      <Stack.Screen name="Kinetics" component={KineticsScreen} options={{ title: 'Reaction Kinetics' }} />
      <Stack.Screen name="CompoundProperties" component={CompoundPropertiesScreen} options={{ title: 'Compound Properties' }} />
      <Stack.Screen name="PidTuning" component={PidTuningScreen} options={{ title: 'PID Tuning' }} />
      <Stack.Screen name="McCabeThiele" component={McCabeThieleScreen} options={{ title: 'McCabe-Thiele' }} />
      <Stack.Screen name="AzeotropeFinder" component={AzeotropeFinderScreen} options={{ title: 'Azeotrope Finder' }} />
      <Stack.Screen name="BinaryPhase" component={BinaryPhaseScreen} options={{ title: 'Binary Phase Diagram' }} />
      <Stack.Screen name="UnaryPhase" component={UnaryPhaseScreen} options={{ title: 'Unary Phase Diagram' }} />
      <Stack.Screen name="ReactorDesign" component={ReactorDesignScreen} options={{ title: 'Reactor Design' }} />
      <Stack.Screen name="FUG" component={FUGScreen} options={{ title: 'FUG Distillation' }} />
      <Stack.Screen name="ResidueCurveMap" component={ResidueCurveMapScreen} options={{ title: 'Residue Curve Map' }} />
      <Stack.Screen name="MolecularDynamics" component={MolecularDynamicsScreen} options={{ title: 'Molecular Dynamics' }} />
    </Stack.Navigator>
  );
}
