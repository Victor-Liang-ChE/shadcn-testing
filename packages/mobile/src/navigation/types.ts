import type { NativeStackScreenProps } from '@react-navigation/native-stack';
import type { BottomTabScreenProps } from '@react-navigation/bottom-tabs';

export type RootStackParamList = {
  Main: undefined;
  WebView: { url: string; title: string };
  // Misc screens
  DropChance: undefined;
  PortolaMenu: undefined;
  JapaneseLyrics: undefined;
  ChemistryTools: undefined;
  ChemeEcon: undefined;
  UfcChampions: undefined;
  LatexConverter: undefined;
  CsvToPlot: undefined;
  // Toolbox screens
  Isotherms: undefined;
  ProcessDynamics: undefined;
  HeatTransfer1D: undefined;
  Kinetics: undefined;
  CompoundProperties: undefined;
  PidTuning: undefined;
  McCabeThiele: undefined;
  AzeotropeFinder: undefined;
  BinaryPhase: undefined;
  UnaryPhase: undefined;
  ReactorDesign: undefined;
  FUG: undefined;
  ResidueCurveMap: undefined;
  MolecularDynamics: undefined;
};

export type TabParamList = {
  Home: undefined;
  Toolbox: undefined;
  Misc: undefined;
};

export type WebViewScreenProps = NativeStackScreenProps<RootStackParamList, 'WebView'>;
export type HomeScreenProps = BottomTabScreenProps<TabParamList, 'Home'>;
export type ToolboxScreenProps = BottomTabScreenProps<TabParamList, 'Toolbox'>;
export type MiscScreenProps = BottomTabScreenProps<TabParamList, 'Misc'>;
