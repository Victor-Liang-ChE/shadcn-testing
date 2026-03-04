import type { RootStackParamList } from '../navigation/types';

const BASE = 'https://victorliang.com';

export interface ToolItem {
  name: string;
  path: string;
  description: string;
  thumbnailPath: string;
  thumbnailLightPath: string;
  category: 'cheme' | 'matsci';
  screenKey: keyof RootStackParamList;
}

export interface MiscProject {
  name: string;
  path: string;
  screenKey: keyof RootStackParamList;
}

export const FEATURED_SIMULATIONS: ToolItem[] = [
  {
    name: 'Reactor Design',
    path: '/toolbox/reactor-design',
    description: 'Design and analyze chemical reactors. Calculate conversions and outlet flow rates for CSTR and PFR configurations.',
    thumbnailPath: `${BASE}/thumbnails/reactor-design-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/reactor-design-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'ReactorDesign',
  },
  {
    name: 'Residue Curve Map',
    path: '/toolbox/residue-curve-map',
    description: 'Visualize and analyze residue curve maps for ternary mixtures, aiding in distillation sequence design.',
    thumbnailPath: `${BASE}/thumbnails/residue-curve-map-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/residue-curve-map-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'ResidueCurveMap',
  },
  {
    name: 'PID Tuning',
    path: '/toolbox/pid-tuning',
    description: 'Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time.',
    thumbnailPath: `${BASE}/thumbnails/pid-tuning-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/pid-tuning-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'PidTuning',
  },
  {
    name: 'McCabe-Thiele',
    path: '/toolbox/mccabe-thiele',
    description: 'Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams.',
    thumbnailPath: `${BASE}/thumbnails/mccabe-thiele-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/mccabe-thiele-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'McCabeThiele',
  },
  {
    name: 'Molecular Dynamics',
    path: '/toolbox/molecular-dynamics',
    description: 'Interactive molecular dynamics simulator for studying particle interactions, energy evolution, and radial distribution functions.',
    thumbnailPath: `${BASE}/thumbnails/molecular-dynamics-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/molecular-dynamics-thumbnail-light.png`,
    category: 'matsci',
    screenKey: 'MolecularDynamics',
  },
  {
    name: 'Unary Phase Diagram',
    path: '/toolbox/unary-phase-diagrams',
    description: 'Interactive phase diagrams for pure compounds showing vaporization, fusion, and sublimation curves with key thermodynamic points.',
    thumbnailPath: `${BASE}/thumbnails/unary-phase-diagrams-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/unary-phase-diagrams-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'UnaryPhase',
  },
  {
    name: 'Binary Phase Diagram',
    path: '/toolbox/binary-phase-diagrams',
    description: 'Generate and visualize binary phase diagrams using various thermodynamic models.',
    thumbnailPath: `${BASE}/thumbnails/binary-phase-diagrams-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/binary-phase-diagrams-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'BinaryPhase',
  },
  {
    name: 'Compound Properties',
    path: '/toolbox/compound-properties',
    description: 'Fetch, plot, and compare various physical and thermodynamic properties of chemical compounds.',
    thumbnailPath: `${BASE}/thumbnails/compound-properties-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/compound-properties-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'CompoundProperties',
  },
  {
    name: 'Azeotrope Finder',
    path: '/toolbox/azeotrope-finder',
    description: 'Predict and visualize azeotropic behavior of binary mixtures using various thermodynamic models.',
    thumbnailPath: `${BASE}/thumbnails/azeotrope-finder-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/azeotrope-finder-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'AzeotropeFinder',
  },
  {
    name: 'Isotherms',
    path: '/toolbox/isotherms',
    description: 'Explore and visualize adsorption isotherms including Langmuir, Freundlich, and Temkin models for surface chemistry analysis.',
    thumbnailPath: `${BASE}/thumbnails/isotherms-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/isotherms-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'Isotherms',
  },
  {
    name: 'Reaction Kinetics',
    path: '/toolbox/kinetics',
    description: 'Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time.',
    thumbnailPath: `${BASE}/thumbnails/kinetics-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/kinetics-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'Kinetics',
  },
  {
    name: 'Process Dynamics',
    path: '/toolbox/process-dynamics',
    description: 'Simulate process dynamics with various inputs and understand system behavior in chemical processes.',
    thumbnailPath: `${BASE}/thumbnails/process-dynamics-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/process-dynamics-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'ProcessDynamics',
  },
  {
    name: '1D Heat Transfer',
    path: '/toolbox/1d-heat-transfer',
    description: 'Interactive heat transfer visualization through multiple layers with temperature controls.',
    thumbnailPath: `${BASE}/thumbnails/1d-heat-transfer-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/1d-heat-transfer-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'HeatTransfer1D',
  },
  {
    name: 'FUG',
    path: '/toolbox/FUG',
    description: 'Multicomponent flash and distillation design using Fenske, Underwood, and Gilliland relations.',
    thumbnailPath: `${BASE}/thumbnails/FUG-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/FUG-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'FUG',
  },
];

export const EXTRA_TOOLS: ToolItem[] = [
  {
    name: 'Biot Number',
    path: '/toolbox/biot',
    description: 'Calculate and visualize the Biot number for heat transfer problems.',
    thumbnailPath: `${BASE}/thumbnails/biot-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/biot-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'WebView',
  },
  {
    name: 'Laplace Transform',
    path: '/toolbox/laplace',
    description: 'Compute and visualize Laplace transforms for process control applications.',
    thumbnailPath: `${BASE}/thumbnails/laplace-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/laplace-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'WebView',
  },
  {
    name: 'Monte Carlo',
    path: '/toolbox/monte-carlo',
    description: 'Monte Carlo simulation tools for stochastic analysis.',
    thumbnailPath: `${BASE}/thumbnails/monte-carlo-thumbnail.png`,
    thumbnailLightPath: `${BASE}/thumbnails/monte-carlo-thumbnail-light.png`,
    category: 'cheme',
    screenKey: 'WebView',
  },
];

export const ALL_TOOLS: ToolItem[] = [...FEATURED_SIMULATIONS, ...EXTRA_TOOLS];

export const MISC_PROJECTS: MiscProject[] = [
  { name: 'Drop Chance Calculator', path: '/misc/drop-chance', screenKey: 'DropChance' },
  { name: 'Japanese Lyrics Analyzer', path: '/misc/japanese-lyrics', screenKey: 'JapaneseLyrics' },
  { name: 'Portola Menu', path: '/misc/portola-menu', screenKey: 'PortolaMenu' },
  { name: 'ChemE Economics Calculator', path: '/misc/cheme-econ', screenKey: 'ChemeEcon' },
  { name: 'LaTeX Constructor & Converter', path: '/misc/latex', screenKey: 'LatexConverter' },
  { name: 'Chemistry Tools', path: '/misc/chemistry-tools', screenKey: 'ChemistryTools' },
  { name: 'UFC Championship Lineage', path: '/misc/ufc-champions', screenKey: 'UfcChampions' },
  { name: 'Fast .csv Plotter', path: '/misc/csv-to-plot', screenKey: 'CsvToPlot' },
];
