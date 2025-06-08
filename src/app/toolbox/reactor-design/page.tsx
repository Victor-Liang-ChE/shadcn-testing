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
  calculateOutletFlowRatesParallel,
  findLimitingReactant
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
  
  // State for components (auto-generated but with editable reaction orders)
  const [components, setComponents] = useState<Component[]>([]);
  
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
  const [maxFlowRateSlider, setMaxFlowRateSlider] = useState<string>('10'); // User-defined max for flow rate slider

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
    // Cap at 6 reactions maximum
    if (reactions.length >= 6) {
      return;
    }
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
    if (field === 'reactionOrders') {
      // Update the reaction orders for a specific component
      setComponents(prevComponents => 
        prevComponents.map(comp => 
          comp.id === id ? { ...comp, reactionOrders: value as { [reactionId: string]: string } } : comp
        )
      );
    } else if (field === 'initialFlowRate') {
      // Update the initial flow rate for a specific component
      setComponents(prevComponents => 
        prevComponents.map(comp => 
          comp.id === id ? { ...comp, initialFlowRate: value as string } : comp
        )
      );
    }
  };

  const handleKineticsChange = (field: keyof Kinetics, value: string) => {
    setKinetics(prev => ({ ...prev, [field]: value }));
  };

  // Auto-detect components from reactions and update state
  useEffect(() => {
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
    const newComponents = detectedNames.map(name => {
      // Find existing component to preserve user input
      const existingComponent = components.find(comp => comp.name === name);
      
      // Set initial flow rate: 1 if reactant in any reaction, 0 if only product
      const initialFlowRate = reactantNames.has(name) ? '1' : '0';
      
      // Initialize reaction orders for each reaction
      const reactionOrders: { [reactionId: string]: string } = {};
      reactions.forEach(reaction => {
        const isReactantInReaction = reaction.reactants && reaction.reactants.toLowerCase().includes(name.toLowerCase());
        const isProductInReaction = reaction.products && reaction.products.toLowerCase().includes(name.toLowerCase());
        
        if (reaction.isEquilibrium) {
          // For equilibrium reactions, add both forward and backward keys
          const forwardKey = `${reaction.id}-fwd`;
          const backwardKey = `${reaction.id}-rev`;
          
          if (isReactantInReaction) {
            reactionOrders[forwardKey] = existingComponent?.reactionOrders?.[forwardKey] || '1';
          }
          if (isProductInReaction) {
            reactionOrders[backwardKey] = existingComponent?.reactionOrders?.[backwardKey] || '1';
          }
        } else {
          // For forward-only reactions
          if (isReactantInReaction) {
            reactionOrders[reaction.id] = existingComponent?.reactionOrders?.[reaction.id] || '1';
          }
        }
      });
      
      return {
        id: `comp-${name}`,
        name,
        initialFlowRate: existingComponent?.initialFlowRate || initialFlowRate,
        reactionOrders
      };
    });
    
    // Only update if components actually changed to prevent loops
    if (JSON.stringify(newComponents) !== JSON.stringify(components)) {
      setComponents(newComponents);
    }
  }, [reactions]); // Only depend on reactions, not components to avoid loops

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

    // Find the limiting reactant based on stoichiometry and initial flow rates
    const limitingReactantInfo = findLimitingReactant(parsedParallelReactions, components);
    if (!limitingReactantInfo) {
      setCalculationError("No valid reactants found to determine limiting reactant.");
      setIsLoading(false);
      return;
    }
    
    const keyReactantName = limitingReactantInfo.name;
    const F_key0 = limitingReactantInfo.initialFlow;
    if (isNaN(F_key0) || F_key0 <= 0) {
        setCalculationError(`Initial flow rate for limiting reactant ${keyReactantName} must be positive.`);
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
  }, [reactorVolume, kinetics.reactionTempK, volumetricFlowRate, handleCalculate]);

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
                  <div className="flex items-center gap-3">
                    <div className="flex items-center gap-1 min-w-0">
                      <Label htmlFor="reactorVolumeSlider" className="text-sm font-medium whitespace-nowrap">Reactor Volume:</Label>
                      <span className="text-sm font-medium">{reactorVolume}</span>
                      <span className="text-xs text-muted-foreground">L</span>
                    </div>
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
                  <div className="flex items-center gap-3">
                    <div className="flex items-center gap-1 min-w-0">
                      <Label htmlFor="temperatureSlider" className="text-sm font-medium whitespace-nowrap">Temperature:</Label>
                      <span className="text-sm font-medium">{kinetics.reactionTempK}</span>
                      <span className="text-xs text-muted-foreground">K</span>
                    </div>
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
                    <div className="flex items-center gap-3">
                      <div className="flex items-center gap-1 min-w-0">
                        <Label htmlFor="volumetricFlowRateSlider" className="text-sm font-medium whitespace-nowrap">Vol Flow Rate:</Label>
                        <span className="text-sm font-medium">{volumetricFlowRate}</span>
                        <span className="text-xs text-muted-foreground">L/s</span>
                      </div>
                      <Slider
                        id="volumetricFlowRateSlider"
                        min={0.1}
                        max={parseFloat(maxFlowRateSlider)}
                        step={0.1}
                        value={[parseFloat(volumetricFlowRate)]}
                        onValueChange={(value) => setVolumetricFlowRate(value[0].toString())}
                        className="flex-1"
                      />
                      <div className="flex items-center gap-1">
                        <Label htmlFor="maxFlowRateInput" className="text-xs text-muted-foreground">Max:</Label>
                        <Input
                          id="maxFlowRateInput"
                          type="number"
                          value={maxFlowRateSlider}
                          onChange={(e) => setMaxFlowRateSlider(e.target.value)}
                          className="w-20 h-8 text-xs"
                          min="0.1"
                          step="0.1"
                        />
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
                  <Button 
                    variant="outline" 
                    size="sm" 
                    onClick={addReaction}
                    disabled={reactions.length >= 6}
                  >
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
                    <div className="text-center flex justify-center" style={{gridColumn: `span ${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}`}}>
                      <Label className="text-sm font-medium">Reaction Order</Label>
                    </div>
                  </div>
                </div>
                <div className={`grid gap-2 items-center text-xs font-medium text-muted-foreground border-b pb-1`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                  <div className="text-center flex justify-center">
                    <Label className="text-sm font-medium">Name</Label>
                  </div>
                  <div className="text-center flex justify-center">
                    <Label className="text-sm font-medium">Initial Flow Rate (mol/s)</Label>
                  </div>
                  {reactions.map((reaction, index) => (
                    reaction.isEquilibrium ? (
                      <React.Fragment key={reaction.id}>
                        <div className="text-center flex justify-center">
                          <Label className="text-xs font-medium">Rxn {index + 1} Fwd</Label>
                        </div>
                        <div className="text-center flex justify-center">
                          <Label className="text-xs font-medium">Rxn {index + 1} Rev</Label>
                        </div>
                      </React.Fragment>
                    ) : (
                      <div key={reaction.id} className="text-center flex justify-center">
                        <Label className="text-xs font-medium">Rxn {index + 1}</Label>
                      </div>
                    )
                  ))}
                </div>
                {components.map((comp, index) => {
                  return (
                    <div key={comp.id} className={`grid gap-2 items-center`} style={{gridTemplateColumns: `1fr 1fr repeat(${reactions.reduce((acc, r) => acc + (r.isEquilibrium ? 2 : 1), 0)}, 1fr)`}}>
                      <div className="flex justify-center">
                        <div className="text-sm font-medium p-2">
                          {comp.name}
                        </div>
                      </div>
                      <div className="flex justify-center">
                        <Input
                          type="number"
                          value={comp.initialFlowRate}
                          onChange={(e) => handleComponentChange(comp.id, 'initialFlowRate', e.target.value)}
                          onKeyDown={handleKeyDown}
                          placeholder="1.0"
                          className="w-full text-center"
                        />
                      </div>
                      {reactions.map((reaction, reactionIndex) => {
                        const isReactantInReaction = reaction.reactants && reaction.reactants.toLowerCase().includes(comp.name.toLowerCase());
                        const isProductInReaction = reaction.products && reaction.products.toLowerCase().includes(comp.name.toLowerCase());
                        
                        if (reaction.isEquilibrium) {
                          const forwardKey = `${reaction.id}-fwd`;
                          const backwardKey = `${reaction.id}-rev`;
                          
                          return (
                            <React.Fragment key={reaction.id}>
                              {/* Forward reaction order */}
                              <div className="flex items-center justify-center">
                                {isReactantInReaction ? (
                                  <Input
                                    type="number"
                                    value={comp.reactionOrders?.[forwardKey] || '1'}
                                    onChange={(e) => {
                                      const newReactionOrders = { ...comp.reactionOrders, [forwardKey]: e.target.value };
                                      handleComponentChange(comp.id, 'reactionOrders', newReactionOrders);
                                    }}
                                    className="w-16 h-8 text-center text-xs"
                                    step="0.1"
                                    min="0"
                                  />
                                ) : (
                                  <Input
                                    type="text"
                                    value="-"
                                    disabled
                                    className="w-16 h-8 text-center text-xs bg-muted text-muted-foreground cursor-not-allowed"
                                  />
                                )}
                              </div>
                              {/* Reverse reaction order */}
                              <div className="flex items-center justify-center">
                                {isProductInReaction ? (
                                  <Input
                                    type="number"
                                    value={comp.reactionOrders?.[backwardKey] || '1'}
                                    onChange={(e) => {
                                      const newReactionOrders = { ...comp.reactionOrders, [backwardKey]: e.target.value };
                                      handleComponentChange(comp.id, 'reactionOrders', newReactionOrders);
                                    }}
                                    className="w-16 h-8 text-center text-xs"
                                    step="0.1"
                                    min="0"
                                  />
                                ) : (
                                  <Input
                                    type="text"
                                    value="-"
                                    disabled
                                    className="w-16 h-8 text-center text-xs bg-muted text-muted-foreground cursor-not-allowed"
                                  />
                                )}
                              </div>
                            </React.Fragment>
                          );
                        } else {
                          return (
                            <div key={reaction.id} className="flex items-center justify-center">
                              {isReactantInReaction ? (
                                <Input
                                  type="number"
                                  value={comp.reactionOrders?.[reaction.id] || '1'}
                                  onChange={(e) => {
                                    const newReactionOrders = { ...comp.reactionOrders, [reaction.id]: e.target.value };
                                    handleComponentChange(comp.id, 'reactionOrders', newReactionOrders);
                                  }}
                                  className="w-16 h-8 text-center text-xs"
                                  step="0.1"
                                  min="0"
                                />
                              ) : (
                                <Input
                                  type="text"
                                  value="-"
                                  disabled
                                  className="w-16 h-8 text-center text-xs bg-muted text-muted-foreground cursor-not-allowed"
                                />
                              )}
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
                  <div className="flex items-center justify-between">
                    <CardTitle>Conversion and Outlet Flow Rates</CardTitle>
                    <div className="flex items-center gap-2 text-xs text-muted-foreground">
                      <div className="flex items-center gap-1">
                        <div className="px-2 py-1 bg-yellow-100 dark:bg-yellow-900/50 border-l-4 border-yellow-500 rounded-sm text-[10px] font-medium text-gray-700 dark:text-gray-300">
                          Limiting Reactant
                        </div>
                      </div>
                    </div>
                  </div>
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
                                // UNIVERSAL SELECTIVITY = (F_product,out - F_product,in) / (F_limitingReactant,in - F_limitingReactant,out)
                                // This definition works for both intermediates and final products
                                const limitingReactantName = calculationResults.conversion.reactantName;
                                
                                // Get limiting reactant flows
                                const limitingReactantInitial = parseFloat(components.find(c => c.name === limitingReactantName)?.initialFlowRate || '0');
                                const limitingReactantFinal = calculationResults.outletFlowRates?.[limitingReactantName] || 0;
                                const limitingReactantConsumed = limitingReactantInitial - limitingReactantFinal;
                                
                                // Get product flows
                                const productInitial = parseFloat(comp.initialFlowRate || '0');
                                const productFinal = outletFlow;
                                const netProductProduced = productFinal - productInitial;
                                
                                if (limitingReactantConsumed > 0 && netProductProduced >= 0) {
                                  // Universal selectivity: Net moles of desired product produced / Moles of limiting reactant consumed
                                  const selectivity = (netProductProduced / limitingReactantConsumed) * 100;
                                  selectivityValue = Math.max(0, selectivity).toPrecision(3);
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
                            <tr key={comp.id} className={`hover:bg-muted/25 ${isLimitingReactant ? 'bg-yellow-100 dark:bg-yellow-900/50 border-l-4 border-yellow-500' : ''}`}>
                              <td className="border border-border px-3 py-2 font-medium text-center">
                                {comp.name}
                              </td>
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
