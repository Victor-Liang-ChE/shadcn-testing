'use client';

import React, { useState, useEffect, useCallback, useMemo } from 'react';
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Select, SelectContent, SelectItem, SelectTrigger, SelectValue } from "@/components/ui/select";
import { Alert, AlertDescription, AlertTitle } from "@/components/ui/alert";
import { TooltipProvider, Tooltip, TooltipTrigger, TooltipContent } from "@/components/ui/tooltip";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Slider } from "@/components/ui/slider";
import { AlertCircle, PlusCircle, RefreshCw, Trash2 } from 'lucide-react';
import { RadioGroup, RadioGroupItem } from "@/components/ui/radio-group";

// Import from lib files
import { 
  Component, 
  Kinetics, 
  ParsedComponent, 
  ParsedParallelReactions,
  ReactionData,
  calculateRateConstant,
  parseCompoundsFromReaction,
  parseParallelReactions
} from '@/lib/reaction-parser';
import { 
  CalculationResults, 
  solveCSTR, 
  solvePFR,
  solveCSTRParallel,
  solvePFRParallel,
  calculateOutletFlowRatesParallel
} from '@/lib/reactor-solver';

// Helper types
type ReactorType = 'PFR' | 'CSTR';
type ReactionPhase = 'Liquid' | 'Gas';

const R_GAS_CONSTANT = 8.314; // J/(mol*K) or L*kPa/(mol*K) depending on units

export default function ReactorDesignPage() {
  const [reactorType, setReactorType] = useState<ReactorType>('CSTR');
  const [reactionPhase, setReactionPhase] = useState<ReactionPhase>('Liquid');
  const [reactionString, setReactionString] = useState<string>('A + B -> C');
  const [reactions, setReactions] = useState<ReactionData[]>([
    { id: '1', reactants: 'A + B', products: 'C', AValue: '1e6', EaValue: '50000', AValueBackward: '1e4', EaValueBackward: '60000', isEquilibrium: false },
    { id: '2', reactants: 'A', products: 'D', AValue: '5e5', EaValue: '45000', AValueBackward: '1e3', EaValueBackward: '55000', isEquilibrium: false }
  ]); // Parallel reactions: A+B->C and A->D
  
  const [kinetics, setKinetics] = useState<Kinetics>({
    rateConstantInputMode: 'directK',
    kValue: '0.1', // Units depend on reaction order and concentration units
    AValue: '1e6',
    EaValue: '50000', // e.g., J/mol
    reactionTempK: '300', // Kelvin
  });

  const [reactorVolume, setReactorVolume] = useState<string>('100'); // e.g., Liters
  const [totalPressure, setTotalPressure] = useState<string>('1'); // e.g., bar (for gas phase)
  const [volumetricFlowRate, setVolumetricFlowRate] = useState<string>('1'); // v0 (for liquid phase, L/s)
  const [maxVolumeSlider, setMaxVolumeSlider] = useState<string>('1000'); // User-defined max for volume slider
  const [maxTemperatureSlider, setMaxTemperatureSlider] = useState<string>('350'); // User-defined max for temperature slider

  const [calculationResults, setCalculationResults] = useState<CalculationResults | null>(null);
  const [isLoading, setIsLoading] = useState(false);
  const [calculationError, setCalculationError] = useState<string | null>(null);

  // Handle Enter key press to trigger calculation
  const handleKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      handleCalculate();
    }
  };

  // Component Management
  const addReaction = () => {
    const newId = (reactions.length + 1).toString();
    setReactions([...reactions, { id: newId, reactants: '', products: '', AValue: '1e6', EaValue: '50000', AValueBackward: '1e4', EaValueBackward: '60000', isEquilibrium: false }]);
  };

  const removeReaction = (id: string) => {
    const updatedReactions = reactions.filter(reaction => reaction.id !== id);
    setReactions(updatedReactions);
    
    // Update the main reaction string to reflect the first remaining reaction
    if (updatedReactions.length > 0) {
      const firstReaction = updatedReactions[0];
      const arrow = firstReaction.isEquilibrium ? '<=>' : '->';
      setReactionString(`${firstReaction.reactants} ${arrow} ${firstReaction.products}`);
    } else {
      // If no reactions left, set a default
      setReactionString('A + B -> C');
    }
  };

  const updateReaction = (id: string, field: 'reactants' | 'products' | 'AValue' | 'EaValue' | 'AValueBackward' | 'EaValueBackward' | 'isEquilibrium', value: string | boolean) => {
    setReactions(reactions.map(reaction => 
      reaction.id === id ? { ...reaction, [field]: value } : reaction
    ));
    
    // Update the main reaction string based on the updated reaction
    if (field === 'reactants' || field === 'products' || field === 'isEquilibrium') {
      const updatedReaction = reactions.find(r => r.id === id);
      if (updatedReaction) {
        const newReactants = field === 'reactants' ? value as string : updatedReaction.reactants;
        const newProducts = field === 'products' ? value as string : updatedReaction.products;
        const reactionIsEquilibrium = field === 'isEquilibrium' ? value as boolean : updatedReaction.isEquilibrium;
        const arrow = reactionIsEquilibrium ? '<=>' : '->';
        setReactionString(`${newReactants} ${arrow} ${newProducts}`);
      }
    }
  };

  const handleComponentChange = (id: string, field: keyof Component, value: string | boolean | { [reactionId: string]: string }) => {
    // Components are now auto-generated from reactions, so this function is no longer needed
  };

  const handleKineticsChange = (field: keyof Kinetics, value: string) => {
    setKinetics(prev => ({ ...prev, [field]: value }));
  };

  // Auto-detect components from reactions using useMemo to prevent infinite loops
  const components = useMemo(() => {
    const allDetectedNames = new Set<string>();
    const reactantNames = new Set<string>();
    const productNames = new Set<string>();
    
    reactions.forEach(reaction => {
      const reactionStr = `${reaction.reactants} -> ${reaction.products}`;
      const names = parseCompoundsFromReaction(reactionStr);
      names.forEach(name => allDetectedNames.add(name));
      
      // Track which components are reactants vs products
      if (reaction.reactants) {
        const reactants = parseCompoundsFromReaction(reaction.reactants + ' -> dummy');
        reactants.forEach(name => reactantNames.add(name));
      }
      if (reaction.products) {
        const products = parseCompoundsFromReaction('dummy -> ' + reaction.products);
        products.forEach(name => productNames.add(name));
      }
    });
    
    const detectedNames = Array.from(allDetectedNames);
    return detectedNames.map(name => {
      // Set initial flow rate: 1 if reactant in any reaction, 0 if only product
      const initialFlowRate = reactantNames.has(name) ? '1' : '0';
      
      return {
        id: `comp-${name}`,
        name,
        initialFlowRate,
        reactionOrders: {}
      };
    });
  }, [reactions]);

  // Memoize parsed parallel reactions to prevent unnecessary recalculations
  const parsedParallelReactions = useMemo(() => {
    const parallelResult = parseParallelReactions(reactionString, components, kinetics, reactions, false);
    return parallelResult;
  }, [reactionString, components, kinetics, reactions]);

  // Update reaction string when reactions change (but don't trigger infinite loops)
  useEffect(() => {
    if (reactions.length > 0) {
      const firstReaction = reactions[0];
      const arrow = firstReaction.isEquilibrium ? '<=>' : '->';
      const newReactionString = `${firstReaction.reactants} ${arrow} ${firstReaction.products}`;
      
      // Only update if different to prevent loops
      if (newReactionString !== reactionString) {
        setReactionString(newReactionString);
      }
    }
  }, [reactions, reactionString]);

  // Main Calculation Logic - Updated for parallel reactions
  const handleCalculate = useCallback(async () => {
    setIsLoading(true);
    setCalculationResults(null);
    setCalculationError(null); // Clear previous errors

    // Check if parallel reactions were parsed successfully
    if (!parsedParallelReactions || parsedParallelReactions.reactions.length === 0) {
      if (!parsedParallelReactions?.error) { 
         setCalculationError("No valid reactions found. Please check reaction stoichiometry and kinetics.");
      } else {
         setCalculationError(parsedParallelReactions.error);
      }
      setIsLoading(false);
      return;
    }

    // Check if any reaction has invalid kinetics
    const invalidReaction = parsedParallelReactions.reactions.find(r => r.rateConstantAtT === undefined);
    if (invalidReaction) {
      if (!calculationError) {
         setCalculationError("Some reactions have invalid kinetics. Please check A/Ea values and temperature.");
      }
      setIsLoading(false);
      return;
    }

    // Check backward kinetics for equilibrium reactions
    const invalidEquilibriumReaction = parsedParallelReactions.reactions.find(r => 
      r.isEquilibriumReaction && r.rateConstantBackwardAtT === undefined
    );
    if (invalidEquilibriumReaction) {
      if (!calculationError) {
        setCalculationError("Some equilibrium reactions have invalid backward kinetics. Please check A_backward and Ea_backward values.");
      }
      setIsLoading(false);
      return;
    }

    const V = parseFloat(reactorVolume);
    const T_K = parseFloat(kinetics.reactionTempK);

    if (isNaN(V) || V <= 0) {
      setCalculationError("Reactor volume must be a positive number.");
      setIsLoading(false);
      return;
    }

    // Find all reactants across all reactions to determine key reactant
    const allReactants = parsedParallelReactions.reactions.flatMap(r => r.reactants);
    if (allReactants.length === 0) {
      setCalculationError("No reactants identified in any reaction.");
      setIsLoading(false);
      return;
    }

    // For simplicity, assume the first reactant in the first reaction is the key reactant
    const keyReactant = allReactants[0];
    const keyReactantName = keyReactant.name;
    const F_key0 = parseFloat(keyReactant.initialFlowRate);
    if (isNaN(F_key0) || F_key0 <= 0) {
        setCalculationError(`Initial flow rate for key reactant ${keyReactantName} must be positive.`);
        setIsLoading(false);
        return;
    }

    // Initial concentrations and total volumetric flow rate
    let C_key0: number; // Concentration of key reactant at inlet
    let v0_calc: number; // Total volumetric flow rate at inlet (L/s or m3/s)

    const totalPressureInBar = parseFloat(totalPressure); // For gas phase, input is in bar
    const v0_input_Ls = parseFloat(volumetricFlowRate); // For liquid phase if provided

    // Calculate initial concentrations based on phase
    const initialConcentrations: { [name: string]: number } = {};
    const F_total0 = components.reduce((sum, c) => sum + (parseFloat(c.initialFlowRate) || 0), 0);

    if (reactionPhase === 'Liquid') {
        if (!isNaN(v0_input_Ls) && v0_input_Ls > 0) {
            v0_calc = v0_input_Ls;
        } else {
            setCalculationError("For liquid phase, please provide total volumetric flow rate (v0).");
            setIsLoading(false); return;
        }
        components.forEach(c => {
            const Fi0 = parseFloat(c.initialFlowRate) || 0;
            initialConcentrations[c.name] = Fi0 / v0_calc;
        });
        C_key0 = initialConcentrations[keyReactantName];

    } else { // Gas Phase
        const P_total_kPa = totalPressureInBar * 100; // Convert bar to kPa

        if (isNaN(P_total_kPa) || P_total_kPa <= 0) {
            setCalculationError("For gas phase, total pressure must be a positive number.");
            setIsLoading(false); return;
        }
        // Ideal Gas Law: C_i = P_i / (R*T) = (y_i * P_total) / (R*T)
        // y_i = F_i0 / F_total0
        // R in L*kPa/(mol*K) is 8.314
        components.forEach(c => {
            const Fi0 = parseFloat(c.initialFlowRate) || 0;
            const yi0 = F_total0 > 0 ? Fi0 / F_total0 : 0;
            initialConcentrations[c.name] = (yi0 * P_total_kPa) / (R_GAS_CONSTANT * T_K);
        });
        C_key0 = initialConcentrations[keyReactantName];
        // v0 for gas phase: F_total0 * R * T_K / P_total_kPa
        v0_calc = (F_total0 * R_GAS_CONSTANT * T_K) / P_total_kPa; 
    }
    
    if (C_key0 <= 0 && F_key0 > 0) { // Check if key reactant concentration is valid
        setCalculationError(`Initial concentration of key reactant ${keyReactantName} is zero or negative. Check feed conditions.`);
        setIsLoading(false); return;
    }

    try {
        let conversion = 0;

        if (reactorType === 'CSTR') {
            // Use the parallel CSTR solver
            conversion = solveCSTRParallel(
                parsedParallelReactions,
                F_key0,
                C_key0,
                V,
                initialConcentrations,
                keyReactantName,
                parsedParallelReactions.allInvolvedComponents
            );

        } else { // PFR
            // Use the parallel PFR solver
            conversion = solvePFRParallel(
                parsedParallelReactions,
                F_key0,
                C_key0,
                V,
                initialConcentrations,
                keyReactantName,
                parsedParallelReactions.allInvolvedComponents
            );
        }

        // Calculate outlet flow rates for all components using parallel reactions
        const outletFlowRates = calculateOutletFlowRatesParallel(
            components,
            parsedParallelReactions,
            keyReactantName,
            F_key0,
            conversion
        );
        
        setCalculationResults({
            conversion: { reactantName: keyReactantName, value: conversion },
            outletFlowRates,
        });

    } catch (e: any) {
        setCalculationError(`Calculation Error: ${e.message}`);
    } finally {
        setIsLoading(false);
    }
  }, [parsedParallelReactions, reactorType, reactorVolume, kinetics.reactionTempK, components, reactionPhase, totalPressure, volumetricFlowRate, reactions]);

  // Auto-calculate when key dependencies change
  useEffect(() => {
    if (parsedParallelReactions && parsedParallelReactions.reactions.length > 0 && !parsedParallelReactions.error) {
      handleCalculate();
    }
  }, [handleCalculate, parsedParallelReactions]);

  // Debounced calculation for slider changes
  useEffect(() => {
    if (parsedParallelReactions && parsedParallelReactions.reactions.length > 0 && !parsedParallelReactions.error) {
      const timer = setTimeout(() => {
        handleCalculate();
      }, 300); // Debounce for smooth slider interaction
      return () => clearTimeout(timer);
    }
  }, [reactorVolume, kinetics.reactionTempK, handleCalculate]);

  // Force recalculation when A or Ea changes (for manual button press)
  // This is now handled automatically by the parseReactionData useCallback dependencies


  // Reactor SVG Visualization
  const renderReactorSVG = () => {
    if (reactorType === 'PFR') {
      return (
        <div className="w-full h-64 bg-card rounded flex items-center justify-center">
          <svg viewBox="0 0 400 200" className="w-full h-full">
            {/* PFR Tube - slightly smaller */}
            <rect x="60" y="70" width="280" height="60" fill="lightblue" stroke="currentColor" strokeWidth="2" rx="30"/>
            
            {/* Feed Arrow - longer tail */}
            <defs>
              <marker id="arrowhead" markerWidth="7" markerHeight="5" refX="6" refY="2.5" orient="auto">
                <polygon points="0 0, 7 2.5, 0 5" fill="currentColor"/>
              </marker>
            </defs>
            <line x1="10" y1="100" x2="55" y2="100" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="32" y="85" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Feed</text>
            
            {/* Product Arrow - longer tail */}
            <line x1="345" y1="100" x2="390" y2="100" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="367" y="85" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Product</text>
            
            {/* PFR Label */}
            <text x="200" y="110" fontSize="20" fill="currentColor" className="text-foreground" textAnchor="middle" fontWeight="bold">PFR</text>
          </svg>
        </div>
      );
    } else {
      return (
        <div className="w-full h-64 bg-card rounded flex items-center justify-center">
          <svg viewBox="0 0 400 200" className="w-full h-full">
            {/* CSTR Tank - much larger rectangular shape */}
            <rect x="120" y="40" width="160" height="120" fill="lightblue" stroke="currentColor" strokeWidth="2" rx="15"/>
            
            {/* Feed Arrow - longer tail */}
            <defs>
              <marker id="arrowhead" markerWidth="7" markerHeight="5" refX="6" refY="2.5" orient="auto">
                <polygon points="0 0, 7 2.5, 0 5" fill="currentColor"/>
              </marker>
            </defs>
            <line x1="60" y1="100" x2="115" y2="100" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="87" y="85" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Feed</text>
            
            {/* Product Arrow - longer tail */}
            <line x1="285" y1="100" x2="340" y2="100" stroke="currentColor" strokeWidth="3" markerEnd="url(#arrowhead)"/>
            <text x="312" y="85" fontSize="14" fill="currentColor" className="text-foreground font-medium" textAnchor="middle">Product</text>
            
            {/* Stirrer shaft */}
            <line x1="200" y1="20" x2="200" y2="140" stroke="currentColor" strokeWidth="4"/>
            
            {/* Stirrer blades with spinning animation - larger */}
            <g className="animate-spin" style={{transformOrigin: '200px 100px'}}>
              <line x1="160" y1="100" x2="240" y2="100" stroke="currentColor" strokeWidth="5"/>
              <line x1="200" y1="70" x2="200" y2="130" stroke="currentColor" strokeWidth="5"/>
            </g>
            
            {/* CSTR Label */}
            <text x="200" y="180" fontSize="20" fill="currentColor" className="text-foreground" textAnchor="middle" fontWeight="bold">CSTR</text>
          </svg>
        </div>
      );
    }
  };




  return (
    <TooltipProvider>
      <div className="container mx-auto p-4 space-y-6">
        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          {/* Left Column: Inputs */}
          <div className="space-y-6">
            {/* Reactor Volume/Flow Rate */}
            <Card>
              <CardContent className="space-y-6 pt-6">
                {/* Reactor Volume Slider */}
                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <div className="flex items-center gap-1">
                      <Label htmlFor="reactorVolumeSlider" className="text-sm font-medium">Reactor Volume:</Label>
                      <span className="text-sm font-medium">{reactorVolume}</span>
                      <span className="text-xs text-muted-foreground">L</span>
                    </div>
                  </div>
                  <div className="flex items-center gap-3">
                    <Slider
                      id="reactorVolumeSlider"
                      min={1}
                      max={parseFloat(maxVolumeSlider)}
                      step={1}
                      value={[parseFloat(reactorVolume)]}
                      onValueChange={(value) => setReactorVolume(value[0].toString())}
                      className="flex-1"
                    />
                    <div className="flex items-center gap-1">
                      <Label htmlFor="maxVolumeInput" className="text-xs text-muted-foreground">Max:</Label>
                      <Input
                        id="maxVolumeInput"
                        type="number"
                        value={maxVolumeSlider}
                        onChange={(e) => setMaxVolumeSlider(e.target.value)}
                        className="w-20 h-8 text-xs"
                        min="1"
                      />
                    </div>
                  </div>
                </div>

                {/* Temperature Slider */}
                <div className="space-y-3">
                  <div className="flex items-center justify-between">
                    <div className="flex items-center gap-1">
                      <Label htmlFor="temperatureSlider" className="text-sm font-medium">Temperature:</Label>
                      <span className="text-sm font-medium">{kinetics.reactionTempK}</span>
                      <span className="text-xs text-muted-foreground">K</span>
                    </div>
                  </div>
                  <div className="flex items-center gap-3">
                    <Slider
                      id="temperatureSlider"
                      min={250}
                      max={parseFloat(maxTemperatureSlider)}
                      step={5}
                      value={[parseFloat(kinetics.reactionTempK)]}
                      onValueChange={(value) => handleKineticsChange('reactionTempK', value[0].toString())}
                      className="flex-1"
                    />
                    <div className="flex items-center gap-1">
                      <Label htmlFor="maxTempInput" className="text-xs text-muted-foreground">Max:</Label>
                      <Input
                        id="maxTempInput"
                        type="number"
                        value={maxTemperatureSlider}
                        onChange={(e) => setMaxTemperatureSlider(e.target.value)}
                        className="w-20 h-8 text-xs"
                        min="250"
                      />
                    </div>
                  </div>
                </div>

                {/* Volumetric Flow Rate / Pressure */}
                <div className="space-y-3">
                  {reactionPhase === 'Gas' && (
                    <div className="flex items-center gap-2">
                      <Label htmlFor="totalPressure" className="text-sm whitespace-nowrap">Total Pressure:</Label>
                      <div className="flex items-center gap-1 flex-1 min-w-0">
                        <Input
                          id="totalPressure"
                          type="number"
                          value={totalPressure}
                          onChange={(e) => setTotalPressure(e.target.value)}
                          onKeyDown={handleKeyDown}
                          placeholder="1"
                        />
                        <span className="text-xs text-muted-foreground">bar</span>
                      </div>
                    </div>
                  )}
                  {reactionPhase === 'Liquid' && (
                    <div className="flex items-center gap-2">
                      <Label htmlFor="volumetricFlowRate" className="text-sm whitespace-nowrap">Vol Flow Rate:</Label>
                      <div className="flex items-center gap-1 flex-1 min-w-0">
                        <Input
                          id="volumetricFlowRate"
                          type="number"
                          value={volumetricFlowRate}
                          onChange={(e) => setVolumetricFlowRate(e.target.value)}
                          onKeyDown={handleKeyDown}
                          placeholder="1.0"
                        />
                        <span className="text-xs text-muted-foreground">L/s</span>
                      </div>
                    </div>
                  )}
                </div>
              </CardContent>
            </Card>

            {/* Reaction Stoichiometry */}
            <Card>
              <CardHeader>
                <CardTitle className="flex items-center justify-between">
                  Reaction Stoichiometry
                  <Button variant="outline" size="sm" onClick={addReaction}>
                    <PlusCircle className="mr-2 h-4 w-4" />Add Reaction
                  </Button>
                </CardTitle>
              </CardHeader>
              <CardContent className="space-y-4">
                {reactions.map((reaction, index) => (
                  <div key={reaction.id}>
                    <div className="grid grid-cols-12 gap-3 items-center">
                      <div className="col-span-2 flex items-center justify-start">
                        <Label className="text-sm font-medium whitespace-nowrap">Reaction {index + 1}:</Label>
                      </div>
                      <div className="col-span-3">
                        <Input
                          placeholder="A + B"
                          value={reaction.reactants}
                          onChange={(e) => updateReaction(reaction.id, 'reactants', e.target.value)}
                        />
                      </div>
                      <div className="col-span-2 text-center flex items-center justify-center">
                        <div className="flex flex-col items-center gap-2">
                          <Tabs value={reaction.isEquilibrium ? "equilibrium" : "forward"} onValueChange={(value) => updateReaction(reaction.id, 'isEquilibrium', value === "equilibrium")}>
                            <TabsList className="grid grid-cols-2 h-8">
                              <TabsTrigger value="forward" className="px-3">→</TabsTrigger>
                              <TabsTrigger value="equilibrium" className="px-3">⇌</TabsTrigger>
                            </TabsList>
                          </Tabs>
                        </div>
                      </div>
                      <div className="col-span-3">
                        <Input
                          placeholder="C"
                          value={reaction.products}
                          onChange={(e) => updateReaction(reaction.id, 'products', e.target.value)}
                        />
                      </div>
                      <div className="col-span-1"></div>
                      <div className="col-span-1 flex justify-start">
                        {index > 0 && (
                          <button
                            onClick={() => removeReaction(reaction.id)}
                            className="text-red-500 hover:text-red-700"
                          >
                            <Trash2 className="h-4 w-4" />
                          </button>
                        )}
                      </div>
                    </div>
                    
                    {/* A and Ea values for this reaction */}
                    <div className="mt-2 pl-4 space-y-3">
                      {/* Forward reaction parameters */}
                      <div className="grid grid-cols-2 gap-4">
                        <div className="flex items-center gap-2">
                          <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                            {reaction.isEquilibrium ? (
                              <span>A<sub style={{fontSize: '0.75em', lineHeight: '1'}}>forward</sub>:</span>
                            ) : (
                              'A:'
                            )}
                          </Label>
                          <Input
                            type="number"
                            value={reaction.AValue}
                            onChange={(e) => updateReaction(reaction.id, 'AValue', e.target.value)}
                            placeholder="1e6"
                            className="flex-1"
                          />
                        </div>
                        <div className="flex items-center gap-2">
                          <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                            {reaction.isEquilibrium ? (
                              <span>Ea<sub style={{fontSize: '0.75em', lineHeight: '1'}}>forward</sub>:</span>
                            ) : (
                              'Ea:'
                            )}
                          </Label>
                          <div className="flex items-center gap-1 flex-1">
                            <Input
                              type="number"
                              value={reaction.EaValue}
                              onChange={(e) => updateReaction(reaction.id, 'EaValue', e.target.value)}
                              placeholder="50000"
                            />
                            <span className="text-xs text-muted-foreground">J/mol</span>
                          </div>
                        </div>
                      </div>
                      
                      {/* Backward reaction parameters (only shown for equilibrium) */}
                      {reaction.isEquilibrium && (
                        <div className="grid grid-cols-2 gap-4">
                          <div className="flex items-center gap-2">
                            <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                              <span>A<sub style={{fontSize: '0.75em', lineHeight: '1'}}>backward</sub>:</span>
                            </Label>
                            <Input
                              type="number"
                              value={reaction.AValueBackward || ''}
                              onChange={(e) => updateReaction(reaction.id, 'AValueBackward', e.target.value)}
                              placeholder="1e4"
                              className="flex-1"
                            />
                          </div>
                          <div className="flex items-center gap-2">
                            <Label className="text-sm whitespace-nowrap" style={{fontFeatureSettings: '"subs" 1'}}>
                              <span>Ea<sub style={{fontSize: '0.75em', lineHeight: '1'}}>backward</sub>:</span>
                            </Label>
                            <div className="flex items-center gap-1 flex-1">
                              <Input
                                type="number"
                                value={reaction.EaValueBackward || ''}
                                onChange={(e) => updateReaction(reaction.id, 'EaValueBackward', e.target.value)}
                                placeholder="60000"
                              />
                              <span className="text-xs text-muted-foreground">J/mol</span>
                            </div>
                          </div>
                        </div>
                      )}
                    </div>
                    
                    {index < reactions.length - 1 && <hr className="my-4" />}
                  </div>
                ))}
              </CardContent>
            </Card>

            {/* Component Feed Composition */}
            <Card>
              <CardHeader>
                <CardTitle className="font-bold">Components</CardTitle>
                {/* Remove the Plus button - components are auto-detected from reaction stoichiometry */}
              </CardHeader>
              <CardContent className="space-y-4">
                <div className="mb-2">
                  <div className={`grid gap-2 items-center text-xs font-medium text-muted-foreground`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                    <div></div>
                    <div></div>
                    <div className="text-center" style={{gridColumn: `span ${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}`}}>
                      <Label className="text-sm font-medium text-center">Reaction Order</Label>
                    </div>
                  </div>
                </div>
                <div className={`grid gap-2 items-center text-xs font-medium text-muted-foreground border-b pb-1`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                  <div className="text-center">
                    <Label className="text-sm font-medium text-center">Name</Label>
                  </div>
                  <div className="text-center">
                    <Label className="text-sm font-medium text-center">Initial Flow Rate (mol/s)</Label>
                  </div>
                  {reactions.map((reaction, index) => (
                    reaction.isEquilibrium ? (
                      <React.Fragment key={reaction.id}>
                        <div className="text-center">
                          <Label className="text-xs font-medium text-center">Rxn {index + 1} Fwd</Label>
                        </div>
                        <div className="text-center">
                          <Label className="text-xs font-medium text-center">Rxn {index + 1} Rev</Label>
                        </div>
                      </React.Fragment>
                    ) : (
                      <div key={reaction.id} className="text-center">
                        <Label className="text-xs font-medium text-center">Rxn {index + 1}</Label>
                      </div>
                    )
                  ))}
                </div>
                {components.map((comp, index) => {
                  return (
                    <div key={comp.id} className={`grid gap-2 items-center`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                      <div>
                        <div className="text-sm font-medium p-2 text-left">
                          {comp.name}
                        </div>
                      </div>
                      <div>
                        <Input
                          type="number"
                          value={comp.initialFlowRate}
                          onChange={(e) => handleComponentChange(comp.id, 'initialFlowRate', e.target.value)}
                          onKeyDown={handleKeyDown}
                          placeholder="1.0"
                        />
                      </div>
                      {reactions.map((reaction, reactionIndex) => {
                        // Simplified: just show that the reaction order is assumed to be 1
                        // for all components participating in the reaction
                        const isReactantInReaction = reaction.reactants && reaction.reactants.toLowerCase().includes(comp.name.toLowerCase());
                        const isProductInReaction = reaction.products && reaction.products.toLowerCase().includes(comp.name.toLowerCase());
                        
                        if (reaction.isEquilibrium) {
                          return (
                            <React.Fragment key={reaction.id}>
                              {/* Forward reaction order - simplified to "1" */}
                              <div className="flex items-center justify-center text-sm">
                                {isReactantInReaction ? "1" : "-"}
                              </div>
                              {/* Reverse reaction order - simplified to "1" */}
                              <div className="flex items-center justify-center text-sm">
                                {isProductInReaction ? "1" : "-"}
                              </div>
                            </React.Fragment>
                          );
                        } else {
                          return (
                            <div key={reaction.id} className="flex items-center justify-center text-sm">
                              {isReactantInReaction ? "1" : "-"}
                            </div>
                          );
                        }
                      })}
                    </div>
                  );
                })}
              </CardContent>
            </Card>
          </div>

          {/* Right Column: Visualization and Results */}
          <div className="space-y-6">
            {/* Reactor Visualization and Type Selection */}
            <Card>
              <CardContent className="p-0">
                <div className="flex justify-center pt-6">
                  <Tabs value={reactorType} onValueChange={(value) => setReactorType(value as 'PFR' | 'CSTR')}>
                    <TabsList className="grid grid-cols-2 h-8">
                      <TabsTrigger value="CSTR" className="px-3">CSTR</TabsTrigger>
                      <TabsTrigger value="PFR" className="px-3">PFR</TabsTrigger>
                    </TabsList>
                  </Tabs>
                </div>
                <div className="p-6 pb-2">
                  {renderReactorSVG()}
                </div>
                {/* <div className="px-6 pb-1">
                  <Button 
                    onClick={handleCalculate} 
                    disabled={isLoading}
                    className="w-full bg-secondary text-secondary-foreground hover:bg-secondary/80"
                  >
                    {isLoading ? 'Calculating...' : ''}
                  </div> */}
              </CardContent>
            </Card>

            {/* Results Display */}
            {calculationResults && (
              <Card>
                <CardHeader>
                  <CardTitle>Conversion and Outlet Flow Rates</CardTitle>
                </CardHeader>
                <CardContent>
                  <div className="overflow-x-auto">
                    <table className="w-full border-collapse border border-border">
                      <thead>
                        <tr className="bg-muted/50">
                          <th className="border border-border px-3 py-2 text-center font-medium">Comp.</th>
                          <th className="border border-border px-3 py-2 text-center font-medium">Conversion (%)</th>
                          <th className="border border-border px-3 py-2 text-center font-medium">Selectivity (%)</th>
                          <th className="border border-border px-3 py-2 text-center font-medium">Outlet Flow Rate (mol/s)</th>
                        </tr>
                      </thead>
                      <tbody>
                        {components.map((comp) => {
                          // Calculate conversion for this component
                          let conversionValue = '';
                          if (calculationResults.conversion && parsedParallelReactions) {
                            // Check if this component is a reactant across ALL reactions
                            const isReactantInAnyReaction = reactions.some(reaction => {
                              const reactantsString = reaction.reactants || '';
                              return reactantsString.toLowerCase().includes(comp.name.toLowerCase());
                            });
                            
                            console.log(`${comp.name} isReactantInAnyReaction:`, isReactantInAnyReaction);
                            
                            if (isReactantInAnyReaction) {
                              // This component is a reactant, calculate conversion
                              if (comp.name === calculationResults.conversion.reactantName) {
                                conversionValue = (calculationResults.conversion.value * 100).toPrecision(3);
                              } else if (calculationResults.conversion.value > 0) {
                                // Calculate conversion for non-limiting reactants by looking across all reactions
                                const limitingReactantName = calculationResults.conversion.reactantName;
                                const limitingConversion = calculationResults.conversion.value;
                                
                                // Find this component in any of the parallel reactions
                                let thisReactantInfo = null;
                                let limitingReactantInfo = null;
                                
                                for (const reaction of parsedParallelReactions.reactions) {
                                  if (!thisReactantInfo) {
                                    thisReactantInfo = reaction.reactants.find(r => r.name === comp.name);
                                  }
                                  if (!limitingReactantInfo) {
                                    limitingReactantInfo = reaction.reactants.find(r => r.name === limitingReactantName);
                                  }
                                }
                                
                                if (thisReactantInfo && limitingReactantInfo) {
                                  const stoichRatio = thisReactantInfo.stoichiometricCoefficient / limitingReactantInfo.stoichiometricCoefficient;
                                  const thisReactantConversion = limitingConversion * stoichRatio;
                                  
                                  const initialFlow = parseFloat(comp.initialFlowRate || '0');
                                  const limitingInitialFlow = parseFloat(components.find(c => c.name === limitingReactantName)?.initialFlowRate || '0');
                                  const maxPossibleConversion = limitingInitialFlow > 0 ? (stoichRatio * limitingInitialFlow) / initialFlow : 0;
                                  
                                  const actualConversion = Math.min(thisReactantConversion, maxPossibleConversion);
                                  conversionValue = Math.max(0, actualConversion * 100).toPrecision(3);
                                } else {
                                  // Simple calculation if not found in parallel reactions
                                  conversionValue = (limitingConversion * 100).toPrecision(3);
                                }
                              }
                            } else {
                              // This component is not a reactant in any reaction, show "-"
                              conversionValue = '-';
                            }
                          }
                          
                          // Get outlet flow rate with 3 sig figs
                          const outletFlow = calculationResults.outletFlowRates?.[comp.name];
                          const outletFlowDisplay = outletFlow !== undefined ? outletFlow.toPrecision(3) : '-';
                          
                          // Calculate selectivity for products
                          let selectivityValue = '';
                          if (calculationResults.conversion && parsedParallelReactions && outletFlow !== undefined) {
                            // Check if this component is a product across ALL reactions
                            const isProductInAnyReaction = reactions.some(reaction => {
                              const productsString = reaction.products || '';
                              return productsString.toLowerCase().includes(comp.name.toLowerCase());
                            });
                            
                            console.log(`${comp.name} isProductInAnyReaction:`, isProductInAnyReaction);
                            
                            if (isProductInAnyReaction) {
                              // This component is a product, calculate selectivity
                              // Find this component in any of the parallel reactions
                              let productInfo = null;
                              for (const reaction of parsedParallelReactions.reactions) {
                                productInfo = reaction.products.find(p => p.name === comp.name);
                                if (productInfo) break;
                              }
                              
                              if (calculationResults.conversion.value > 0) {
                                // Selectivity = moles of product formed / moles of limiting reactant consumed
                                const limitingReactantName = calculationResults.conversion.reactantName;
                                const limitingReactantInitialFlow = components.find(c => c.name === limitingReactantName)?.initialFlowRate;
                                const limitingReactantInitial = parseFloat(limitingReactantInitialFlow || '0');
                                
                                if (limitingReactantInitial > 0) {
                                  const limitingReactantConsumed = limitingReactantInitial * calculationResults.conversion.value;
                                  const productInitialFlow = parseFloat(comp.initialFlowRate || '0');
                                  const productFormed = outletFlow - productInitialFlow;
                                  
                                  if (limitingReactantConsumed > 0 && productFormed >= 0) {
                                    if (productInfo) {
                                      // Get stoichiometric coefficients to calculate theoretical selectivity
                                      const productStoich = productInfo.stoichiometricCoefficient || 1;
                                      
                                      // Find limiting reactant info across all reactions
                                      let limitingReactantInfo = null;
                                      for (const reaction of parsedParallelReactions.reactions) {
                                        limitingReactantInfo = reaction.reactants.find(r => r.name === limitingReactantName);
                                        if (limitingReactantInfo) break;
                                      }
                                      
                                      const limitingReactantStoich = limitingReactantInfo?.stoichiometricCoefficient || 1;
                                      
                                      // Theoretical moles of product = (productStoich/limitingReactantStoich) * moles of limiting reactant consumed
                                      const theoreticalProductFormed = (productStoich / limitingReactantStoich) * limitingReactantConsumed;
                                      
                                      if (theoreticalProductFormed > 0) {
                                        const selectivity = (productFormed / theoreticalProductFormed) * 100;
                                        selectivityValue = Math.min(100, Math.max(0, selectivity)).toPrecision(3);
                                      }
                                    } else {
                                      // Simple selectivity calculation if not found in parallel reactions
                                      const selectivity = (productFormed / limitingReactantConsumed) * 100;
                                      selectivityValue = Math.min(100, Math.max(0, selectivity)).toPrecision(3);
                                    }
                                  }
                                }
                              }
                            } else {
                              // This component is not a product in any reaction, show "-"
                              selectivityValue = '-';
                            }
                          }
                          
                          // Check if this component is the limiting reactant (simplified check)
                          const isLimitingReactant = calculationResults.conversion && 
                            calculationResults.conversion.reactantName === comp.name;
                          
                          return (
                            <tr key={comp.id} className={`hover:bg-muted/25 ${isLimitingReactant ? 'bg-yellow-100 dark:bg-yellow-900/30' : ''}`}>
                              <td className="border border-border px-3 py-2 font-medium text-center">{comp.name}</td>
                              <td className="border border-border px-3 py-2 text-center">
                                {conversionValue}
                              </td>
                              <td className="border border-border px-3 py-2 text-center">{selectivityValue}</td>
                              <td className="border border-border px-3 py-2 text-center">{outletFlowDisplay}</td>
                            </tr>
                          );
                        })}
                      </tbody>
                    </table>
                  </div>
                </CardContent>
              </Card>
            )}

            {/* Error Display */}
            {calculationError && (
              <Alert variant="destructive">
                <AlertCircle className="h-4 w-4" />
                <AlertTitle>Error</AlertTitle>
                <AlertDescription>{calculationError}</AlertDescription>
              </Alert>
            )}
          </div>
        </div>
      </div>
    </TooltipProvider>
  );
}
