'use client';

import React, { useState, useEffect, useRef } from 'react';
import { Card, CardContent, CardHeader, CardTitle, CardDescription, CardFooter } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Button } from "@/components/ui/button";
import { Slider } from "@/components/ui/slider";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Checkbox } from "@/components/ui/checkbox";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import {
  Tooltip,
  TooltipContent,
  TooltipProvider,
  TooltipTrigger,
} from "@/components/ui/tooltip";
import { Terminal, Info } from "lucide-react";

// Global type declaration for JSMol
declare global {
  interface Window {
    Jmol: any;
    jQuery: any;
  }
}

// --- Interfaces ---
interface ChemicalFormulaPart { text: string; isSubscript: boolean; }
interface Compound { compound: string; coefficient: number; }
interface ReactionData { reactants: Compound[]; products: Compound[]; }
interface MolarMasses { [compound: string]: number; }
interface InputData {
  [compound: string]: {
    moles?: number; grams?: number; molar_mass?: number;
    moles_display?: string; grams_display?: string;
  };
}
interface StoichiometryResults {
  limiting_reactants: string[];
  reactants: { [compound: string]: { initial_moles: number; used_moles: number; excess_moles: number; initial_grams?: number; used_grams?: number; excess_grams?: number; }; };
  products: { [compound: string]: { produced_moles: number; produced_grams?: number; }; };
  conversion_percentage: number;
}

// --- Helper Functions ---
const formatNumber = (num: number | undefined, options: { type?: 'default' | 'molarMass', sigFigs?: number } = {}): string => { /* ... (no change) ... */
  const { type = 'default', sigFigs = 3 } = options;
  if (num === undefined || isNaN(num)) return "0";
  if (num === 0) { return type === 'molarMass' ? num.toFixed(3) : "0"; }
  if (type === 'molarMass') { return num.toFixed(3); }
  const scientificStr = num.toExponential(sigFigs - 1);
  const [coefficient, exponentPart] = scientificStr.split('e');
  const exponent = parseInt(exponentPart);
  if (exponent >= -3 && exponent <= 5) {
    const decimalPlaces = Math.max(0, sigFigs - Math.floor(Math.log10(Math.abs(num))) - 1);
    if (Math.abs(num) < Math.pow(10, -decimalPlaces)) { return num.toFixed(decimalPlaces); }
    return Number(num.toFixed(decimalPlaces)).toString();
  }
  return scientificStr;
};
const ChemicalFormula: React.FC<{ formula: string }> = ({ formula }) => { /* ... (no change) ... */
    const formatFormula = (formula: string): ChemicalFormulaPart[] => {
    const parts: ChemicalFormulaPart[] = []; let i = 0;
    const coeffMatch = formula.match(/^(\d+)/); let displayFormula = formula;
    if (coeffMatch) { parts.push({ text: coeffMatch[1], isSubscript: false }); displayFormula = formula.substring(coeffMatch[1].length).trim(); }
    i = 0;
    while (i < displayFormula.length) {
        const char = displayFormula[i];
        if (char === '(') { parts.push({ text: char, isSubscript: false }); i++; }
        else if (char === ')') { parts.push({ text: char, isSubscript: false }); i++; let subscript = ''; while (i < displayFormula.length && /\d/.test(displayFormula[i])) { subscript += displayFormula[i]; i++; } if (subscript) { parts.push({ text: subscript, isSubscript: true }); } }
        else if (/[A-Z]/.test(char)) { let element = char; i++; while (i < displayFormula.length && /[a-z]/.test(displayFormula[i])) { element += displayFormula[i]; i++; } parts.push({ text: element, isSubscript: false }); let subscript = ''; while (i < displayFormula.length && /\d/.test(displayFormula[i])) { subscript += displayFormula[i]; i++; } if (subscript) { parts.push({ text: subscript, isSubscript: true }); } }
        else { console.warn("Unexpected character in formula:", char); parts.push({ text: char, isSubscript: false }); i++; }
    } return parts;
  };
  const parts = formatFormula(formula);
  return (<span className="inline-block whitespace-nowrap">{parts.map((part, index) => part.isSubscript ? <sub key={index} className="text-xs relative -bottom-1">{part.text}</sub> : <span key={index}>{part.text}</span>)}</span>);
};

// --- Main Component ---
export default function ChemistryTools() { /* ... (no change) ... */
  const [activeTab, setActiveTab] = useState<string>('stoichiometry');
  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 md:p-8">
        <Tabs value={activeTab} onValueChange={setActiveTab} className="w-full mb-6">
          <TabsList className="grid w-full grid-cols-2"><TabsTrigger value="stoichiometry">Stoichiometry Calculator</TabsTrigger><TabsTrigger value="molecular">Molecular Visualization</TabsTrigger></TabsList>
          <TabsContent value="stoichiometry"><StoichiometryCalculator /></TabsContent>
          <TabsContent value="molecular"><MolecularVisualization /></TabsContent>
        </Tabs>
      </div>
    </TooltipProvider>
  );
}

// --- Stoichiometry Calculator Component ---
const StoichiometryCalculator: React.FC = () => {
  const [equation, setEquation] = useState<string>('');
  const [reactionData, setReactionData] = useState<ReactionData | null>(null);
  const [molarMasses, setMolarMasses] = useState<MolarMasses>({});
  const [inputData, setInputData] = useState<InputData>({});
  const [conversionPercentage, setConversionPercentage] = useState<number>(100);
  const [showConversion, setShowConversion] = useState<boolean>(false);
  const [results, setResults] = useState<StoichiometryResults | null>(null);
  const [errorMessage, setErrorMessage] = useState<string>('');
  const [isLoading, setIsLoading] = useState<boolean>(false);
  const [inputModes, setInputModes] = useState<Record<string, 'moles' | 'grams'>>({});
  const [hasCalculated, setHasCalculated] = useState<boolean>(false);
  const [reactants, setReactants] = useState<string>('2H2 + O2');
  const [products, setProducts] = useState<string>('2H2O');

  // Initialize with default reaction on component mount
  useEffect(() => {
    const initializeDefaultReaction = async () => {
      const fullEquation = `${reactants} → ${products}`;
      setEquation(fullEquation);
      setErrorMessage('');
      setIsLoading(true);
      
      try {
        const parsed = parseReactionEquation(fullEquation);
        
        if (!parsed) {
          console.warn('Failed to parse default reaction:', fullEquation);
          return; // Silently fail for default reaction
        }
        
        const compoundsList = [...parsed.reactants, ...parsed.products];
        const molarMassesData: MolarMasses = {};
        const initialInputData: InputData = {};
        const initialModes: Record<string, 'moles' | 'grams'> = {};
        
        for (const compound of compoundsList) {
          try {
            const molarMass = calculateMolarMassLocal(compound.compound, elementMasses);
            molarMassesData[compound.compound] = molarMass;
            initialInputData[compound.compound] = { molar_mass: molarMass };
            initialModes[compound.compound] = 'moles';
          } catch (error) {
            console.error(`Molar mass calculation error for ${compound.compound}:`, error);
            return; // Silently fail for default reaction
          }
        }
        
        // Set default input values: 2 mol H2 and 3 mol O2
        initialInputData['H2'] = {
          ...initialInputData['H2'],
          moles: 2,
          moles_display: '2'
        };
        initialInputData['O2'] = {
          ...initialInputData['O2'],
          moles: 3,
          moles_display: '3'
        };
        
        setMolarMasses(molarMassesData);
        setInputData(initialInputData);
        setInputModes(initialModes);
        setReactionData(parsed);
        
        // Auto-calculate with default values
        setTimeout(() => {
          try {
            const newResults = calculateStoichiometry(parsed, initialInputData, 100);
            setResults(newResults);
            setHasCalculated(true);
          } catch (calcError) {
            console.error('Default calculation error:', calcError);
            // Don't show error for default calculation
          }
        }, 100);
        
      } catch (error) {
        console.error('Default reaction initialization error:', error);
        // Don't show error message for default reaction failures
      } finally {
        setIsLoading(false);
      }
    };
    
    initializeDefaultReaction();
  }, []); // Only run on mount

  // --- Calculation Logic (Unchanged) ---
   const calculateStoichiometry = (reaction: ReactionData, currentInputData: InputData, convPercentage: number): StoichiometryResults => {
        const reactants = reaction.reactants; const products = reaction.products; const reactantData: Record<string, { moles: number; coefficient: number; moles_per_coefficient: number }> = {}; let hasPositiveInput = false; let hasZeroInputReactant = false; const epsilon = 1e-9;
        for (const reactant of reactants) {
            const compound = reactant.compound; const data = currentInputData[compound] || {}; const molarMass = data.molar_mass; let molesFromInput: number = 0;
            if (data.moles !== undefined && data.moles !== null) { molesFromInput = Math.max(0, data.moles); } else if (data.grams !== undefined && data.grams !== null && molarMass !== undefined && molarMass > 0) { molesFromInput = Math.max(0, data.grams / molarMass); }
            if (molesFromInput > epsilon) { hasPositiveInput = true; } else { molesFromInput = 0; if (reactant.coefficient > 0) { hasZeroInputReactant = true; } }
            reactantData[compound] = { moles: molesFromInput, coefficient: reactant.coefficient, moles_per_coefficient: molesFromInput > 0 && reactant.coefficient > 0 ? molesFromInput / reactant.coefficient : Infinity };
        }
        const calcResults: StoichiometryResults = { limiting_reactants: [], reactants: {}, products: {}, conversion_percentage: convPercentage };
        let minPositiveMolesPerCoef = Infinity; let calculationBasis = 0;
        if (hasPositiveInput && !hasZeroInputReactant) { for (const compound in reactantData) { if (reactantData[compound].moles_per_coefficient > 0 && reactantData[compound].moles_per_coefficient < minPositiveMolesPerCoef) { minPositiveMolesPerCoef = reactantData[compound].moles_per_coefficient; } } if (isFinite(minPositiveMolesPerCoef)) { calculationBasis = minPositiveMolesPerCoef; } }
        const limitingSet = new Set<string>();
        if (hasPositiveInput) { let theoreticalMinRatio = Infinity; for (const compound in reactantData) { if (reactantData[compound].moles_per_coefficient > 0 && reactantData[compound].moles_per_coefficient < theoreticalMinRatio) { theoreticalMinRatio = reactantData[compound].moles_per_coefficient; } } if (!isFinite(theoreticalMinRatio)) theoreticalMinRatio = Infinity; for (const compound in reactantData) { const data = reactantData[compound]; if (data.moles === 0 && data.coefficient > 0) { limitingSet.add(compound); } else if (data.moles > 0 && data.moles_per_coefficient > 0 && Math.abs(data.moles_per_coefficient - theoreticalMinRatio) < epsilon) { limitingSet.add(compound); } } } else { reactants.forEach(r => limitingSet.add(r.compound)); }
        calcResults.limiting_reactants = Array.from(limitingSet); const conversion = convPercentage / 100.0;
        for (const reactant of reactants) { const compound = reactant.compound; const coef = reactant.coefficient; const molarMass = currentInputData[compound]?.molar_mass; const initialMoles = reactantData[compound]?.moles ?? 0; const usedMoles = calculationBasis * coef * conversion; const excessMoles = Math.max(0, initialMoles - usedMoles); calcResults.reactants[compound] = { initial_moles: initialMoles, used_moles: usedMoles, excess_moles: excessMoles }; if (molarMass) { calcResults.reactants[compound].initial_grams = initialMoles * molarMass; calcResults.reactants[compound].used_grams = usedMoles * molarMass; calcResults.reactants[compound].excess_grams = excessMoles * molarMass; } }
        for (const product of products) { const compound = product.compound; const coef = product.coefficient; const molarMass = currentInputData[compound]?.molar_mass; const producedMoles = calculationBasis * coef * conversion; calcResults.products[compound] = { produced_moles: producedMoles }; if (molarMass) { calcResults.products[compound].produced_grams = producedMoles * molarMass; } }
        return calcResults;
    };

  // --- Event Handlers & Effects ---
  const handleCalculate = () => { /* ... (no change) ... */
    setErrorMessage(''); if (!reactionData) { setErrorMessage("Please set a valid reaction equation first."); return; }
    const hasValidInput = reactionData.reactants.some(reactant => { const data = inputData[reactant.compound]; return (data?.moles !== undefined && data.moles > 0) || (data?.grams !== undefined && data.grams > 0); });
    if (!hasCalculated && !hasValidInput) { setErrorMessage("Please provide a positive amount for at least one reactant."); setResults(null); return; }
    setIsLoading(true);
    try {
        const newResults = calculateStoichiometry(reactionData, inputData, showConversion ? conversionPercentage : 100);
        setResults(newResults);
        setHasCalculated(true);
    } catch (error) { setErrorMessage(`Calculation Error: ${error instanceof Error ? error.message : String(error)}`); setResults(null); }
    finally { setIsLoading(false); }
  };
  useEffect(() => { /* ... (no change) ... */
    if (reactionData && hasCalculated) { handleCalculate(); }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [conversionPercentage, showConversion, reactionData, hasCalculated]);
  const toggleInputMode = (compound: string, newMode: 'moles' | 'grams') => { /* ... (no change) ... */
    if (inputModes[compound] !== newMode) {
      setInputModes(prev => ({ ...prev, [compound]: newMode }));
      setInputData(prev => ({ ...prev, [compound]: { ...prev[compound], moles: undefined, grams: undefined, moles_display: '', grams_display: '', } }));
    }
  };
  const parseReactionEquation = (equation: string): ReactionData | null => { /* ... (no change) ... */
    const sides = equation.split(/→|->|=|ArrowRight/); if (sides.length !== 2) return null; const [reactantsSide, productsSide] = sides;
    const parseCompounds = (side: string): Compound[] => {
      return side.split('+').map(part => part.trim()).filter(part => part).map(part => {
        const match = part.match(/^(\d*)\s*([A-Za-z0-9()]+)$/); if (!match) return null; const [, coeffStr, compoundFormula] = match; if (!/^[A-Za-z0-9()]+$/.test(compoundFormula)) return null; const coefficient = coeffStr ? parseInt(coeffStr, 10) : 1; if (isNaN(coefficient) || coefficient < 0) return null; return { compound: compoundFormula.trim(), coefficient };
      }).filter((item): item is Compound => item !== null);
    };
    const reactants = parseCompounds(reactantsSide); const products = parseCompounds(productsSide); if (reactants.length === 0 || products.length === 0) return null;
    const reactantNames = reactants.map(r => r.compound); const productNames = products.map(p => p.compound); if (new Set(reactantNames).size !== reactantNames.length || new Set(productNames).size !== productNames.length) { throw new Error("Duplicate compounds found on the same side of the equation."); } return { reactants, products };
  };
  const handleEquationSubmit = async () => { /* ... (no change) ... */
    const fullEquation = `${reactants} → ${products}`; setEquation(fullEquation); setErrorMessage(''); setIsLoading(true); setReactionData(null); setResults(null); setInputData({}); setMolarMasses({}); setInputModes({}); setHasCalculated(false); setShowConversion(false); setConversionPercentage(100);
    if (!reactants.trim() && !products.trim()) { setErrorMessage('Please enter reactants and products.'); setIsLoading(false); return; } if (!reactants.trim()) { setErrorMessage('Please enter reactants.'); setIsLoading(false); return; } if (!products.trim()) { setErrorMessage('Please enter products.'); setIsLoading(false); return; }
    try {
      const parsed = parseReactionEquation(fullEquation); if (!parsed) { if (fullEquation.split(/→|->|=/ ).length !== 2) { throw new Error('Invalid format. Use "+" between compounds and "→" (or "->" or "=") between reactants and products.'); } else { throw new Error('Invalid format. Check compound formulas (e.g., Na2SO4, H2O, Fe(OH)3) and coefficients.'); } }
      const compoundsList = [...parsed.reactants, ...parsed.products]; const molarMassesData: MolarMasses = {}; const initialInputData: InputData = {}; const initialModes: Record<string, 'moles' | 'grams'> = {};
      for (const compound of compoundsList) { try { const molarMass = calculateMolarMassLocal(compound.compound, elementMasses); molarMassesData[compound.compound] = molarMass; initialInputData[compound.compound] = { molar_mass: molarMass }; initialModes[compound.compound] = 'moles'; } catch (error) { console.error(`Molar mass calculation error for ${compound.compound}:`, error); throw new Error(`Invalid formula or unknown element in "${compound.compound}". ${error instanceof Error ? error.message : ''}`); } }
      setMolarMasses(molarMassesData); setInputData(initialInputData); setInputModes(initialModes); setReactionData(parsed);
    } catch (error) { setErrorMessage(`Equation Error: ${error instanceof Error ? error.message : String(error)}`); setReactionData(null); setMolarMasses({}); setInputData({}); setInputModes({}); }
    finally { setIsLoading(false); }
  };
  const calculateMolarMassLocal = (formula: string, masses: Record<string, number>): number => { /* ... (no change) ... */
    let i = 0; let totalMass = 0;
    const parseSegment = (multiplier: number = 1): number => {
        let segmentMass = 0; let currentElement = ''; let currentCountStr = '';
        const closeSegment = () => { if (currentElement) { if (typeof masses[currentElement] !== 'number') { console.error("Element Masses Data:", masses); throw new Error(`Unknown or invalid element mass for: ${currentElement}`); } const count = currentCountStr ? parseInt(currentCountStr, 10) : 1; if (isNaN(count)) throw new Error(`Invalid count for element ${currentElement}`); segmentMass += masses[currentElement] * count; currentElement = ''; currentCountStr = ''; } };
        while (i < formula.length) {
            const char = formula[i];
            if (char === '(') { closeSegment(); i++; const subMass = parseSegment(); let groupMultiplierStr = ''; while (i < formula.length && /\d/.test(formula[i])) { groupMultiplierStr += formula[i]; i++; } const groupMultiplier = groupMultiplierStr ? parseInt(groupMultiplierStr, 10) : 1; if (isNaN(groupMultiplier)) throw new Error(`Invalid multiplier for group ()`); segmentMass += subMass * groupMultiplier; continue; }
            if (char === ')') { closeSegment(); i++; return segmentMass * multiplier; } if (/[A-Z]/.test(char)) { closeSegment(); currentElement = char; i++; } else if (/[a-z]/.test(char)) { if (!currentElement) throw new Error(`Lowercase letter '${char}' must follow an uppercase element letter.`); currentElement += char; i++; } else if (/\d/.test(char)) { if (!currentElement && (i === 0 || formula[i-1] !== ')')) { throw new Error(`Number '${char}' must follow an element or parenthesis.`); } if (currentElement) { currentCountStr += char; } else { throw new Error(`Unexpected number '${char}'.`); } i++; } else { throw new Error(`Invalid character in formula: '${char}'`); }
        } closeSegment(); return segmentMass * multiplier;
    }; totalMass = parseSegment(); if (totalMass <= 0) { throw new Error("Calculated molar mass is zero or negative, check formula structure."); } return totalMass;
  };
  const handleInputChange = (compound: string, inputType: 'moles' | 'grams', valueStr: string | undefined) => { /* ... (no change) ... */
    let numericValue: number | undefined = undefined; let displayValue = valueStr ?? '';
    if (valueStr !== undefined && valueStr.trim() !== '') {
        const cleanedValue = valueStr.replace(/[^0-9.]/g, ''); if ((cleanedValue.match(/\./g) || []).length > 1) { displayValue = inputData[compound]?.[`${inputType}_display`] ?? ''; } else {
            displayValue = cleanedValue; const parsed = parseFloat(cleanedValue); if (!isNaN(parsed) && parsed >= 0) { numericValue = parsed; } else if (cleanedValue === '.' || cleanedValue === '0.' || cleanedValue === '0' || cleanedValue === '') { numericValue = undefined; } else { numericValue = undefined; }
        }
    } else { displayValue = ''; numericValue = undefined; }
    setInputData(prev => {
        const currentCompoundData = prev[compound] || {}; const otherInputType = inputType === 'moles' ? 'grams' : 'moles'; const otherDisplayKey = `${otherInputType}_display` as keyof InputData[string]; const currentDisplayKey = `${inputType}_display` as keyof InputData[string];
        return { ...prev, [compound]: { ...currentCompoundData, [inputType]: numericValue, [currentDisplayKey]: displayValue, [otherInputType]: (numericValue !== undefined && numericValue > 0) ? undefined : currentCompoundData[otherInputType], [otherDisplayKey]: (numericValue !== undefined && numericValue > 0) ? '' : currentCompoundData[otherDisplayKey], } };
    });
  };

  // --- JSX Rendering ---
  return (
    <div className="space-y-6">
      {/* Reaction Input */}
      <Card>
        <CardContent className="pt-2 pb-2 space-y-4">
          <div className="flex flex-col sm:flex-row items-stretch sm:items-center gap-2 sm:gap-4">
            <div className="flex-1 flex items-center gap-2"><Input id="reactants-input" type="text" value={reactants} onChange={(e) => setReactants(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && handleEquationSubmit()} placeholder="Reactants (e.g., 2H2 + O2)" className="flex-grow" /></div>
            <span className="text-xl font-semibold mx-2 text-center sm:text-left">→</span>
            <div className="flex-1 flex items-center gap-2"><Input id="products-input" type="text" value={products} onChange={(e) => setProducts(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && handleEquationSubmit()} placeholder="Products (e.g., 2H2O)" className="flex-grow" /></div>
            <Button onClick={handleEquationSubmit} disabled={isLoading} className="w-full sm:w-auto mt-2 sm:mt-0">{isLoading ? 'Processing...' : 'Set Reaction'}</Button>
          </div>
          {reactionData && (<div className="text-center text-lg mt-4 p-3 bg-muted rounded-md shadow-sm">{reactionData.reactants.map((r, i) => (<React.Fragment key={`r-${i}`}>{i > 0 && ' + '}<ChemicalFormula formula={r.coefficient > 1 ? `${r.coefficient}${r.compound}` : r.compound} /></React.Fragment>))}<span className="mx-2 font-semibold">→</span>{reactionData.products.map((p, i) => (<React.Fragment key={`p-${i}`}>{i > 0 && ' + '}<ChemicalFormula formula={p.coefficient > 1 ? `${p.coefficient}${p.compound}` : p.compound} /></React.Fragment>))}</div>)}
          {errorMessage && (<Alert variant="destructive" className="mt-4"><Terminal className="h-4 w-4" /><AlertTitle>Error</AlertTitle><AlertDescription>{errorMessage}</AlertDescription></Alert>)}
        </CardContent>
      </Card>

      {/* Input/Results Area */}
      {reactionData && (
        <div className="space-y-6">
          {/* Reactant Cards */}
          <Card>
            <CardContent className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 pt-2">
              {reactionData.reactants.map((reactant) => {
                const compound = reactant.compound; const molarMass = molarMasses[compound]; const mode = inputModes[compound] || 'moles'; const res = results?.reactants[compound]; const isLimiting = results?.limiting_reactants?.includes(compound) ?? false; const currentInputValue = inputData[compound]?.[`${mode}_display`] ?? '';
                return (
                  <Card key={`reactant-${compound}`} className="flex flex-col overflow-hidden">
                    <CardHeader className="space-y-2">
                       <div className="flex justify-between items-baseline gap-2">
                            <div className="flex items-baseline gap-2">
                                <Tooltip delayDuration={150}>
                                    <TooltipTrigger asChild><span className="inline-flex items-baseline gap-1 cursor-help"><CardTitle className="text-lg"><ChemicalFormula formula={compound} /></CardTitle>{molarMass && <Info className="h-3 w-3 text-muted-foreground shrink-0" />}</span></TooltipTrigger>
                                    {molarMass && (<TooltipContent><p>Molar Mass: {formatNumber(molarMass, { type: 'molarMass' })} g/mol</p></TooltipContent>)}
                                </Tooltip>
                                {isLimiting && (<span className="text-xs font-semibold text-red-500">(Limiting)</span>)}
                            </div>
                            <Tabs value={mode} onValueChange={(newMode) => toggleInputMode(compound, newMode as 'moles' | 'grams')} className="w-auto shrink-0">
                                <TabsList className="grid grid-cols-2 h-8 text-xs px-0"><TabsTrigger value="moles" className="text-xs h-full px-3">Moles</TabsTrigger><TabsTrigger value="grams" className="text-xs h-full px-3">Grams</TabsTrigger></TabsList>
                            </Tabs>
                       </div>
                    </CardHeader>
                    <CardContent className="flex-grow">
                      {mode === 'moles' ? (<div className="relative"><Label htmlFor={`moles-${compound}`} className="sr-only">Initial Moles for {compound}</Label><Input id={`moles-${compound}`} type="text" inputMode="decimal" placeholder="Initial Moles" value={currentInputValue} onChange={(e) => handleInputChange(compound, 'moles', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && handleCalculate()} className="pr-12 hide-number-arrows" aria-label={`Initial moles for ${compound}`} /><span className="absolute right-3 top-1/2 -translate-y-1/2 text-xs text-muted-foreground pointer-events-none">mol</span></div>)
                       : (<div className="relative"><Label htmlFor={`grams-${compound}`} className="sr-only">Initial Grams for {compound}</Label><Input id={`grams-${compound}`} type="text" inputMode="decimal" placeholder="Initial Grams" value={currentInputValue} onChange={(e) => handleInputChange(compound, 'grams', e.target.value)} onKeyDown={(e) => e.key === 'Enter' && handleCalculate()} className="pr-8 hide-number-arrows" aria-label={`Initial grams for ${compound}`} /><span className="absolute right-3 top-1/2 -translate-y-1/2 text-xs text-muted-foreground pointer-events-none">g</span></div>)}
                    </CardContent>
                    {results && res && (<CardFooter className="flex justify-between text-sm pt-3 pb-3 border-t"><span>Used: {mode === 'moles' ? `${formatNumber(res.used_moles)} mol` : `${formatNumber(res.used_grams)} g`}</span><span>Excess: {mode === 'moles' ? `${formatNumber(res.excess_moles)} mol` : `${formatNumber(res.excess_grams)} g`}</span></CardFooter>)}
                  </Card>
                );
              })}
            </CardContent>
            <CardFooter className="pt-0 flex justify-center"><Button onClick={handleCalculate} disabled={isLoading} size="lg" className="w-full sm:w-auto">{isLoading ? 'Calculating...' : 'Calculate Stoichiometry'}</Button></CardFooter>
          </Card>

          {/* Product Cards */}
          <Card>
            <CardHeader>
                {hasCalculated && (
                    <div className="flex flex-col sm:flex-row items-start sm:items-center gap-4">
                        <div className="flex items-center space-x-2">
                            <Checkbox id="conversion-check" checked={showConversion} onCheckedChange={(checked) => { setShowConversion(Boolean(checked)); }} aria-label="Enable reaction conversion percentage" />
                            <Label htmlFor="conversion-check" className="text-sm font-medium whitespace-nowrap cursor-pointer">Set Conversion ({conversionPercentage}%)</Label>
                        </div>
                        {showConversion && (
                            <div className="flex-1">
                                <Label htmlFor="conversion-slider" className="sr-only">Conversion Percentage Slider</Label>
                                <Slider id="conversion-slider" min={0} max={100} step={1} value={[conversionPercentage]} onValueChange={(value) => setConversionPercentage(value[0])} aria-label={`Conversion percentage: ${conversionPercentage}%`} />
                            </div>
                        )}
                    </div>
                )}
            </CardHeader>
            <CardContent className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4">
              {reactionData.products.map((product) => {
                const compound = product.compound; const molarMass = molarMasses[compound]; const res = results?.products[compound];
                return (
                  <Card key={`product-${compound}`} className="overflow-hidden">
                    <CardHeader className="pb-0 flex flex-col items-start">
                        <Tooltip delayDuration={150}>
                           <TooltipTrigger asChild><span className="inline-flex items-baseline gap-1 cursor-help"><CardTitle className="text-lg"><ChemicalFormula formula={compound} /></CardTitle>{molarMass && <Info className="h-3 w-3 text-muted-foreground shrink-0" />}</span></TooltipTrigger>
                            {molarMass && (<TooltipContent><p>Molar Mass: {formatNumber(molarMass, { type: 'molarMass' })} g/mol</p></TooltipContent>)}
                        </Tooltip>
                    </CardHeader>
                    <CardContent className="pt-0 pb-0">
                       {results && res ? (<p className="text-sm">Produced: <span className="font-medium">{formatNumber(res.produced_grams)} g</span><span className="text-muted-foreground"> ({formatNumber(res.produced_moles)} mol)</span></p>)
                        : (<p className="text-sm text-muted-foreground">Calculate to see results.</p>)}
                    </CardContent>
                  </Card>
                );
              })}
            </CardContent>
          </Card>
        </div>
      )}
      <style jsx global>{`.hide-number-arrows::-webkit-outer-spin-button,.hide-number-arrows::-webkit-inner-spin-button {-webkit-appearance: none; margin: 0;} .hide-number-arrows {-moz-appearance: textfield;}`}</style>
    </div>
  );
};

// --- Molecular Visualization Component ---
const MolecularVisualization: React.FC = () => {
  const [chemicalName, setChemicalName] = useState<string>('aspirin');
  const [loading, setLoading] = useState<boolean>(false);
  const [error, setError] = useState<string>('');
  const [moleculeData, setMoleculeData] = useState<string>('');
  const jsmolContainerRef = useRef<HTMLDivElement>(null);
  const [jsmolLoaded, setJsmolLoaded] = useState<boolean>(false);
  const [jsmolInitialized, setJsmolInitialized] = useState<boolean>(false);

  // Load JSMol library
  useEffect(() => {
    const loadJSMol = () => {
      if (typeof window !== 'undefined' && !window.Jmol) {
        // First load the jQuery dependency if needed
        if (!(window as any).jQuery) {
          const jqueryScript = document.createElement('script');
          jqueryScript.src = 'https://code.jquery.com/jquery-3.6.0.min.js';
          jqueryScript.onload = () => {
            // After jQuery loads, load JSMol
            const jmolScript = document.createElement('script');
            jmolScript.src = 'https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js';
            jmolScript.onload = () => {
              console.log('JSMol library loaded successfully');
              setJsmolLoaded(true);
            };
            jmolScript.onerror = () => {
              console.error('Failed to load JSMol library');
              setError('Failed to load molecular viewer library');
            };
            document.head.appendChild(jmolScript);
          };
          jqueryScript.onerror = () => {
            setError('Failed to load required dependencies');
          };
          document.head.appendChild(jqueryScript);
        } else {
          // jQuery already exists, just load JSMol
          const jmolScript = document.createElement('script');
          jmolScript.src = 'https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js';
          jmolScript.onload = () => {
            console.log('JSMol library loaded successfully');
            setJsmolLoaded(true);
          };
          jmolScript.onerror = () => {
            console.error('Failed to load JSMol library');
            setError('Failed to load molecular viewer library');
          };
          document.head.appendChild(jmolScript);
        }
      } else if (window.Jmol) {
        console.log('JSMol already available');
        setJsmolLoaded(true);
      }
    };
    
    loadJSMol();
  }, []);

  // Auto-load aspirin when JSMol is ready
  useEffect(() => {
    if (jsmolLoaded && !moleculeData && chemicalName === 'aspirin') {
      handleVisualize();
    }
  }, [jsmolLoaded]);

  const handleVisualize = async () => {
    if (!chemicalName.trim()) { 
      setError("Please enter a chemical name."); 
      return; 
    }
    
    setLoading(true); 
    setError(''); 
    setMoleculeData(''); 
    setJsmolInitialized(false);
    
    // Try to fetch from PubChem for any molecule
    try {
      // First, try to get the compound CID by name
      const cidResponse = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${encodeURIComponent(chemicalName)}/cids/JSON`);
      
      if (!cidResponse.ok) {
        throw new Error(`Sorry, "${chemicalName}" couldn't be found in our database. Please try common chemical names like: water, caffeine, aspirin, glucose, benzene, ethanol, or acetone.`);
      }
      
      const cidData = await cidResponse.json();
      const cid = cidData.IdentifierList.CID[0];
      
      // Now get the 3D SDF structure using the CID
      const sdfResponse = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF?record_type=3d`);
      
      if (!sdfResponse.ok) {
        // Try 2D if 3D is not available
        const sdf2dResponse = await fetch(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/${cid}/SDF`);
        
        if (!sdf2dResponse.ok) {
          throw new Error(`No molecular structure data available for "${chemicalName}".`);
        }
        
        const sdf2dData = await sdf2dResponse.text();
        setMoleculeData(sdf2dData);
      } else {
        const sdfData = await sdfResponse.text();
        setMoleculeData(sdfData);
      }
      
    } catch (err) {
      const errorMsg = err instanceof Error ? err.message : String(err);
      setError(`${errorMsg}`);
      console.error("Visualization fetch error:", err);
    } finally {
      setLoading(false);
    }
  };

  // Initialize JSMol viewer when molecule data is available
  useEffect(() => {
    if (jsmolLoaded && moleculeData && jsmolContainerRef.current && typeof window !== 'undefined' && window.Jmol && !jsmolInitialized) {
      try {
        console.log('Initializing JSMol with data:', moleculeData.substring(0, 50) + '...');
        
        // Clear previous content
        jsmolContainerRef.current.innerHTML = '';
        
        // Create a unique container div for this applet
        const appletContainer = document.createElement('div');
        appletContainer.id = 'jmol-container-' + Date.now();
        appletContainer.style.width = '100%';
        appletContainer.style.height = '600px';
        jsmolContainerRef.current.appendChild(appletContainer);
        
        // JSMol configuration
        const jmolInfo = {
          width: '100%',
          height: 600,
          color: '#F0F8FF',
          use: 'HTML5',
          j2sPath: 'https://chemapps.stolaf.edu/jmol/jsmol/j2s',
          serverURL: 'https://chemapps.stolaf.edu/jmol/jsmol/php/jsmol.php',
          readyFunction: () => {
            console.log('JSMol applet ready');
            setJsmolInitialized(true);
          },
          script: `
            load data "molecule"
${moleculeData}
end "molecule";
            wireframe 0.15;
            spacefill 20%;
            color cpk;
            spin on;
            background white;
            zoom 80;
            set echo top left;
            echo "";
          `,
          disableJ2SLoadMonitor: true,
          disableInitialConsole: true,
          addSelectionOptions: false
        };
        
        // Create unique applet ID
        const appletId = 'jmolApplet_' + Date.now();
        
        // Use JSMol's direct HTML insertion method
        window.Jmol.setDocument(document);
        const htmlContent = window.Jmol.getAppletHtml(appletId, jmolInfo);
        
        // Insert the HTML content directly
        appletContainer.innerHTML = htmlContent;
        
        // Execute any scripts that JSMol needs
        const scripts = appletContainer.querySelectorAll('script');
        scripts.forEach(script => {
          if (script.innerHTML) {
            try {
              eval(script.innerHTML);
            } catch (e) {
              console.warn('Script execution warning:', e);
            }
          }
        });
        
        console.log('JSMol applet created and added to DOM');
        
      } catch (err) {
        console.error('JSMol initialization error:', err);
        setError('Failed to initialize molecular viewer: ' + (err instanceof Error ? err.message : String(err)));
        setJsmolInitialized(false);
      }
    }
  }, [jsmolLoaded, moleculeData, jsmolInitialized]);

  return (
    <div className="space-y-6">
      <Card>
        <CardContent className="pt-6 pb-6 space-y-4">
          <div className="flex flex-col sm:flex-row items-stretch sm:items-center gap-2">
             <Label htmlFor="chemical-name-vis" className="sr-only sm:not-sr-only sm:w-auto mb-1 sm:mb-0">Chemical Name:</Label>
            <Input id="chemical-name-vis" type="text" value={chemicalName} onChange={(e) => setChemicalName(e.target.value)} onKeyDown={(e) => e.key === 'Enter' && handleVisualize()} placeholder="e.g., caffeine, aspirin, glucose, benzene" className="flex-1" aria-label="Chemical name for visualization" />
            <Button onClick={handleVisualize} disabled={loading} className="w-full sm:w-auto mt-2 sm:mt-0">{loading ? 'Loading...' : 'Visualize'}</Button>
          </div>
          {error && (<Alert variant="destructive" className="mt-4"><Terminal className="h-4 w-4" /><AlertTitle>Visualization Error</AlertTitle><AlertDescription style={{whiteSpace: 'pre-line'}}>{error}</AlertDescription></Alert>)}
        </CardContent>
      </Card>
      {(loading || moleculeData || (error && !moleculeData)) && (
         <Card className="min-h-[650px] flex flex-col">
            <CardContent className="flex-grow flex items-center justify-center p-6 overflow-hidden relative">
              {loading && (<div className="absolute inset-0 flex items-center justify-center bg-background/80 z-10"><p>Loading molecule...</p></div>)}
              {moleculeData ? (
                <div className="w-full h-full min-h-[600px] relative rounded-lg overflow-hidden">
                  <div ref={jsmolContainerRef} className="w-full h-full" />
                  {!jsmolInitialized && jsmolLoaded && (
                    <div className="absolute inset-0 flex items-center justify-center bg-background/80">
                      <div className="text-center">
                        <p className="mb-2">Initializing 3D viewer...</p>
                        <p className="text-sm text-muted-foreground">If this takes too long, try refreshing the page</p>
                      </div>
                    </div>
                  )}
                </div>
              ) : !loading && !error ? (
                <p className="text-muted-foreground">Enter a chemical name above to visualize.</p>
              ) : null }
            </CardContent>
         </Card>
      )}
    </div>
  );
};

// --- Element Masses Data ---
const elementMasses: Record<string, number> = {
  'H': 1.008, 'He': 4.0026, 'Li': 6.94, 'Be': 9.0122, 'B': 10.81, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998, 'Ne': 20.180, 'Na': 22.990, 'Mg': 24.305, 'Al': 26.982, 'Si': 28.085, 'P': 30.974, 'S': 32.06, 'Cl': 35.45, 'Ar': 39.948, 'K': 39.098, 'Ca': 40.078, 'Sc': 44.956, 'Ti': 47.867, 'V': 50.942, 'Cr': 51.996, 'Mn': 54.938, 'Fe': 55.845, 'Co': 58.933, 'Ni': 58.693, 'Cu': 63.546, 'Zn': 65.38, 'Ga': 69.723, 'Ge': 72.630, 'As': 74.922, 'Se': 78.971, 'Br': 79.904, 'Kr': 83.798, 'Rb': 85.468, 'Sr': 87.62, 'Y': 88.906, 'Zr': 91.224, 'Nb': 92.906, 'Mo': 95.95, 'Tc': 98, 'Ru': 101.07, 'Rh': 102.91, 'Pd': 106.42, 'Ag': 107.87, 'Cd': 112.41, 'In': 114.82, 'Sn': 118.71, 'Sb': 121.76, 'Te': 127.60, 'I': 126.90, 'Xe': 131.29, 'Cs': 132.91, 'Ba': 137.33, 'La': 138.91, 'Ce': 140.12, 'Pr': 140.91, 'Nd': 144.24, 'Pm': 145, 'Sm': 150.36, 'Eu': 151.96, 'Gd': 157.25, 'Tb': 158.93, 'Dy': 162.50, 'Ho': 164.93, 'Er': 167.26, 'Tm': 168.93, 'Yb': 173.05, 'Lu': 174.97, 'Hf': 178.49, 'Ta': 180.95, 'W': 183.84, 'Re': 186.21, 'Os': 190.23, 'Ir': 192.22, 'Pt': 195.08, 'Au': 196.97, 'Hg': 200.59, 'Tl': 204.38, 'Pb': 207.2, 'Bi': 208.98, 'Po': 209, 'At': 210, 'Rn': 222, 'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232.04, 'Pa': 231.04, 'U': 238.03, 'Np': 237, 'Pu': 244, 'Am': 243, 'Cm': 247, 'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 262, 'Rf': 267, 'Db': 268, 'Sg': 271, 'Bh': 272, 'Hs': 270, 'Mt': 276, 'Ds': 281, 'Rg': 280, 'Cn': 285, 'Nh': 286, 'Fl': 289, 'Mc': 290, 'Lv': 293, 'Ts': 294, 'Og': 294
};