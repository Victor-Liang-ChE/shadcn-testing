"use client";

import { useState, useMemo, FC, useEffect, useLayoutEffect, useRef } from "react";
import { motion } from "framer-motion";

import { Button } from "@/components/ui/button";
import { Input } from "@/components/ui/input";
import { Slider } from "@/components/ui/slider";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Label } from "@/components/ui/label";
import { ArrowUp, ArrowDown, Trash2 } from "lucide-react";

// Define the structure for a single layer
interface Layer {
  id: string;
  type: "solid" | "fluid";
  thickness: number; // in meters
  kValue: number; // for solids (W/m·K)
  hValue: number; // for fluids (W/m²·K)
}

// Function to generate a unique ID for new layers
const generateId = () => `layer_${Date.now()}_${Math.random()}`;

// A custom vertical switch for Solid/Fluid selection
const VerticalSwitch = ({
  isSolid,
  onTypeChange,
}: {
  isSolid: boolean;
  onTypeChange: (type: "solid" | "fluid") => void;
}) => (
  <div className="flex flex-col items-center justify-between self-stretch bg-muted rounded-lg p-1">
    <Button
      variant={!isSolid ? 'default' : 'ghost'}
      size="sm"
      onClick={() => onTypeChange('fluid')}
      className="flex-1 w-full"
    >
      Fluid
    </Button>
    <Button
      variant={isSolid ? 'default' : 'ghost'}
      size="sm"
      onClick={() => onTypeChange('solid')}
      className="flex-1 w-full"
    >
      Solid
    </Button>
  </div>
);

// Animated heat wave component
const HeatWaveAnimation = ({ totalWidth }: { totalWidth: number }) => {
  const pathRef = useRef<SVGPathElement>(null);
  const [pathLength, setPathLength] = useState(0);
  const [animationKey, setAnimationKey] = useState(0);

  useLayoutEffect(() => {
    if (pathRef.current) {
      const newPathLength = pathRef.current.getTotalLength();
      setPathLength(newPathLength);
      // Only update the key (and thus restart the animation) after the new length is set
      setAnimationKey(prevKey => prevKey + 1);
    }
  }, [totalWidth]);

  const svgWidth = totalWidth + 40;
  const pathStartX = 10;
  const pathEndX = pathStartX + totalWidth + 20;
  const pathWidth = pathEndX - pathStartX;

  const segmentLength = 60;
  const yOffset = 8;
  const mainPathD = `M ${pathStartX} 25 Q ${pathStartX + pathWidth / 4} 0, ${pathStartX + pathWidth / 2} 25 T ${pathEndX} 25`;
  const topPathD = `M ${pathStartX} ${25 - yOffset} Q ${pathStartX + pathWidth / 4} ${0 - yOffset}, ${pathStartX + pathWidth / 2} ${25 - yOffset} T ${pathEndX} ${25 - yOffset}`;
  const bottomPathD = `M ${pathStartX} ${25 + yOffset} Q ${pathStartX + pathWidth / 4} ${0 + yOffset}, ${pathStartX + pathWidth / 2} ${25 + yOffset} T ${pathEndX} ${25 + yOffset}`;

  return (
    <svg 
      width={svgWidth} 
      height="50" 
      viewBox={`0 0 ${svgWidth} 50`} 
      className="absolute inset-0 m-auto overflow-visible z-20"
    >
        <path
          ref={pathRef}
          d={mainPathD}
          stroke="hsl(var(--muted-foreground))"
          strokeOpacity="0.3"
          strokeWidth="3"
          fill="transparent"
          strokeLinecap="round"
        />
        <path d={topPathD} stroke="hsl(var(--muted-foreground))" strokeOpacity="0.1" strokeWidth="3" fill="transparent" strokeLinecap="round" />
        <path d={bottomPathD} stroke="hsl(var(--muted-foreground))" strokeOpacity="0.1" strokeWidth="3" fill="transparent" strokeLinecap="round" />
        {pathLength > 0 && (
          <>
            {/* Main Center Path */}
            <motion.path
              key={animationKey}
              d={mainPathD}
              stroke="url(#gradient)"
              strokeWidth="3"
              fill="transparent"
              strokeLinecap="round"
              strokeDasharray={`${segmentLength} ${pathLength}`}
              animate={{ strokeDashoffset: [pathLength + segmentLength, 0] }}
              transition={{
                duration: 3,
                repeat: Infinity,
                ease: "linear",
              }}
            />
            {/* Top Path */}
            <motion.path
              key={`${animationKey}-top`}
              d={topPathD}
              stroke="url(#gradient)"
              strokeOpacity={0.3}
              strokeWidth="3"
              fill="transparent"
              strokeLinecap="round"
              strokeDasharray={`${segmentLength} ${pathLength}`}
              animate={{ strokeDashoffset: [pathLength + segmentLength, 0] }}
              transition={{
                duration: 3,
                repeat: Infinity,
                ease: "linear",
              }}
            />
            {/* Bottom Path */}
            <motion.path
              key={`${animationKey}-bottom`}
              d={bottomPathD}
              stroke="url(#gradient)"
              strokeOpacity={0.3}
              strokeWidth="3"
              fill="transparent"
              strokeLinecap="round"
              strokeDasharray={`${segmentLength} ${pathLength}`}
              animate={{ strokeDashoffset: [pathLength + segmentLength, 0] }}
              transition={{
                duration: 3,
                repeat: Infinity,
                ease: "linear",
              }}
            />
          </>
        )}
        <defs>
          <linearGradient id="gradient" x1="0%" y1="0%" x2="100%" y2="0%">
            <stop offset="0%" stopColor="#ef4444" />
            <stop offset="100%" stopColor="#3b82f6" />
          </linearGradient>
        </defs>
    </svg>
  );
};


export default function HeatTransferPage() {
  const [isClient, setIsClient] = useState(false);
  const [layers, setLayers] = useState<Layer[]>([]);

  useEffect(() => {
      setIsClient(true);
      setLayers([{ id: generateId(), type: "solid", thickness: 0.1, kValue: 50, hValue: 10 }]);
  }, []);

  const [area, setArea] = useState<number | string>(1);
  // Toggles to determine which side(s) have user-defined temperatures
  const [hotActive, setHotActive] = useState<boolean>(true);
  const [coldActive, setColdActive] = useState<boolean>(false);

  // Temperature slider bounds
  const [tempMin, setTempMin] = useState<string>('0');
  const [tempMax, setTempMax] = useState<string>('200');

  const [hotTemp, setHotTemp] = useState<number | string>(100);
  const [coldTemp, setColdTemp] = useState<number | string>(0);
  const [heatFlow, setHeatFlow] = useState<number | string>(5000);
  const [heatFlowMax, setHeatFlowMax] = useState<string>('10000');
  const [areaMax, setAreaMax] = useState<string>('10');

  const addLayer = () => {
    if (layers.length < 9) {
      setLayers([...layers, { id: generateId(), type: "solid", thickness: 0.1, kValue: 50, hValue: 10 }]);
    }
  };

  const updateLayer = (id: string, newValues: Partial<Layer>) => {
    setLayers(layers.map((layer) => (layer.id === id ? { ...layer, ...newValues } : layer)));
  };
    
  const deleteLayer = (id: string) => {
    if (layers.length > 1) {
      setLayers(layers.filter(layer => layer.id !== id));
    }
  };

  const moveLayer = (index: number, direction: 'up' | 'down') => {
    const newLayers = [...layers];
    const targetIndex = direction === 'up' ? index - 1 : index + 1;

    if (targetIndex < 0 || targetIndex >= newLayers.length) {
      return; // Prevent moving out of bounds
    }

    // Swap the elements
    const temp = newLayers[index];
    newLayers[index] = newLayers[targetIndex];
    newLayers[targetIndex] = temp;

    setLayers(newLayers);
  };

  const gridClassName = useMemo(() => {
    const itemCount = layers.length;
    if (itemCount <= 3) return "grid grid-cols-1 gap-4";
    if (itemCount <= 6) return "grid grid-cols-1 md:grid-cols-2 gap-4";
    return "grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4";
  }, [layers.length]);

  const layersContainerRef = useRef<HTMLDivElement>(null);
  const [visualWidth, setVisualWidth] = useState(0);

  const simulationResults = useMemo(() => {
    const numArea = Number(area) || 0;
    const numHotTempInput = Number(hotTemp) || 0;
    const numColdTempInput = Number(coldTemp) || 0;
    const numHeatFlowInput = Number(heatFlow) || 0;

    if (numArea === 0) {
      return {
        heatFlow: 0,
        temps: [numHotTempInput, numColdTempInput],
        hotTemp: numHotTempInput,
        coldTemp: numColdTempInput,
      } as const;
    }

    // Calculate thermal resistance of each layer
    const resistances = layers.map((layer) =>
      layer.type === "solid"
        ? layer.thickness / (layer.kValue * numArea)
        : 1 / (layer.hValue * numArea)
    );
    const totalResistance = resistances.reduce((acc, r) => acc + r, 0);

    const bothActive = hotActive && coldActive;

    let q: number; // heat flow (W)
    let hotT: number; // calculated hot side temperature (°C)
    let coldT: number; // calculated cold side temperature (°C)

    if (bothActive) {
      hotT = numHotTempInput;
      coldT = numColdTempInput;
      q = (hotT - coldT) / totalResistance;
    } else if (hotActive) {
      hotT = numHotTempInput;
      q = numHeatFlowInput;
      coldT = Math.max(hotT - q * totalResistance, -273.15);
    } else {
      // only cold active
      coldT = numColdTempInput;
      q = numHeatFlowInput;
      hotT = Math.max(coldT + q * totalResistance, -273.15);
    }

    // Temperature profile across layers
    const temps: number[] = [hotT];
    let current = hotT;
    resistances.forEach((r) => {
      current = Math.max(current - q * r, -273.15);
      temps.push(current);
    });

    return {
      heatFlow: q,
      temps,
      hotTemp: hotT,
      coldTemp: coldT,
    } as const;
  }, [layers, area, hotTemp, coldTemp, heatFlow, hotActive, coldActive]);
  
  useLayoutEffect(() => {
    if (layersContainerRef.current) {
        setVisualWidth(layersContainerRef.current.offsetWidth);
    }
  }, [layers, simulationResults.temps]);

  const handleNumberChange = (setter: React.Dispatch<React.SetStateAction<number | string>>) => (e: React.ChangeEvent<HTMLInputElement>) => {
    const value = e.target.value;
    if (value === "") {
      setter("");
    } else {
      const numValue = Number(value);
      // Clamp temperature inputs to absolute zero
      if (setter === setHotTemp || setter === setColdTemp) {
        setter(Math.max(numValue, -273.15));
      } else {
        setter(numValue);
      }
    }
  };

  // Helper function to get the number of decimal places from a string or number
  const getDecimalPlaces = (value: number | string): number => {
    if (typeof value === 'string') {
      const parts = value.split('.');
      return parts.length > 1 ? parts[1].length : 0;
    }
    const str = value.toString();
    const parts = str.split('.');
    return parts.length > 1 ? parts[1].length : 0;
  };

  // Helper to compute a 'nice' slider step given the maximum value
  const computeStep = (maxVal: number): number => {
    const desiredSteps = 100;
    const raw = maxVal / desiredSteps;
    const exponent = Math.floor(Math.log10(raw));
    const pow10 = Math.pow(10, exponent);
    const mant = raw / pow10;
    let niceMant: number;
    if (mant <= 1) niceMant = 1;
    else if (mant <= 2) niceMant = 2;
    else if (mant <= 5) niceMant = 5;
    else niceMant = 10;
    return niceMant * pow10;
  };

  // Constrain hot temperature to never drop below cold, and vice-versa
  const updateTemperature = (
    type: 'hot' | 'cold',
    value: number,
  ) => {
    if (type === 'hot') {
      const clamped = Math.max(value, Number(coldTemp));
      setHotTemp(clamped);
    } else {
      const clamped = Math.min(value, Number(hotTemp));
      setColdTemp(clamped);
    }
  };

  // Helper function to format number to 3 significant figures
  const to3SigFigs = (value: number): string => {
    if (value === 0) return '0';
    const magnitude = Math.floor(Math.log10(Math.abs(value)));
    const precision = Math.max(0, 2 - magnitude);
    return value.toFixed(precision);
  };

  // Helper function to convert heat flow to appropriate unit with 3 sig figs
  const formatHeatFlow = (value: number): string => {
    if (value === 0) return '0 W';
    
    const absValue = Math.abs(value);
    let unit = 'W';
    let convertedValue = value;
    
    if (absValue >= 1e15) {
      unit = 'PW';
      convertedValue = value / 1e15;
    } else if (absValue >= 1e12) {
      unit = 'TW';
      convertedValue = value / 1e12;
    } else if (absValue >= 1e9) {
      unit = 'GW';
      convertedValue = value / 1e9;
    } else if (absValue >= 1e6) {
      unit = 'MW';
      convertedValue = value / 1e6;
    } else if (absValue >= 1e3) {
      unit = 'kW';
      convertedValue = value / 1e3;
    }
    
    return `${to3SigFigs(convertedValue)} ${unit}`;
  };

  const toggleHot = () => {
    if (hotActive && !coldActive) return; // ensure at least one remains active
    setHotActive(!hotActive);
  };

  const toggleCold = () => {
    if (coldActive && !hotActive) return; // ensure at least one remains active
    setColdActive(!coldActive);
  };

  const bothActive = hotActive && coldActive;
  const isHotLeft = hotActive; // left slider shows hot when hot is active, else cold

  if (!isClient) {
      return null;
  }

  return (
    <div className="container mx-auto p-4 space-y-8">
      {/* Visualization Section */}
      <Card>
        <CardContent>
          <div className="relative flex items-center h-48 rounded-lg p-4">
            <HeatWaveAnimation totalWidth={visualWidth} />
            <div className="flex-1 z-10 text-center flex flex-col items-center justify-center">
              <div className="font-bold text-lg text-red-500">Hot Side</div>
              <div className="text-foreground">{simulationResults.hotTemp.toFixed(1)}°C</div>
            </div>
            <div ref={layersContainerRef} className="z-10 flex h-full items-center">
              {layers.map((layer, index) => (
                <div key={layer.id} className="flex items-center h-full">
                  <div
                    className="h-full flex items-center justify-center text-white px-2 rounded-sm"
                    style={{
                      width: `${Math.max(layer.thickness * 400, 20)}px`,
                      backgroundColor: layer.type === "solid" ? "#22c55e" : "#F5A623",
                    }}
                  >
                    <span className="font-bold text-xl">{index + 1}</span>
                  </div>
                  {index < layers.length - 1 && (
                     <div className="text-center h-full flex flex-col justify-center px-2">
                       <div className="font-bold text-foreground">{simulationResults.temps[index + 1].toFixed(1)}°C</div>
                    </div>
                  )}
                </div>
              ))}
            </div>
            <div className="flex-1 z-10 text-center flex flex-col items-center justify-center">
              <div className="font-bold text-lg text-blue-500">Cold Side</div>
              <div className="text-foreground">{simulationResults.coldTemp.toFixed(1)}°C</div>
            </div>
          </div>

          <div className="text-center mt-4">
            <p className="text-foreground">Heat Flow: <span className="font-bold">{formatHeatFlow(simulationResults.heatFlow)}</span></p>
          </div>
        </CardContent>
      </Card>

      {/* Controls Section */}
      <Card>
        <CardContent className="space-y-6 pt-6">
          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            {/* Left Column: Hot/Cold Switch and Temperature Slider */}
            <div className="space-y-6">
              {/* Hot/Cold Switch */}
              <div className="flex items-center justify-center gap-2 bg-muted rounded-lg p-1 w-full mt-1.5">
                <Button variant={hotActive ? 'default' : 'ghost'} size="sm" onClick={toggleHot} className="flex-1">Hot</Button>
                <Button variant={coldActive ? 'default' : 'ghost'} size="sm" onClick={toggleCold} className="flex-1">Cold</Button>
              </div>

              {/* Temperature Slider(s) */}
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label htmlFor="leftTempSlider">{isHotLeft ? 'Hot Temp (°C):' : 'Cold Temp (°C):'} {Number(isHotLeft ? hotTemp : coldTemp).toFixed(2)}</Label>
                  <div className="flex items-center gap-1">
                    <Label className="text-xs text-muted-foreground">Min:</Label>
                    <Input type="text" value={tempMin} onChange={(e)=>setTempMin(e.target.value)} className="w-16 h-8 text-xs" />
                    <Label className="text-xs text-muted-foreground">Max:</Label>
                    <Input type="text" value={tempMax} onChange={(e)=>setTempMax(e.target.value)} className="w-16 h-8 text-xs" />
                  </div>
                </div>
                <Slider
                  id="leftTempSlider"
                  min={parseFloat(tempMin)}
                  max={parseFloat(tempMax)}
                  step={computeStep(parseFloat(tempMax))}
                  value={[Number(isHotLeft ? hotTemp : coldTemp)]}
                  onValueChange={([v]) => (isHotLeft ? updateTemperature('hot', v) : updateTemperature('cold', v))}
                  className="w-full"
                />
              </div>
            </div>

            {/* Right Column: Conditional Slider(s) */}
            <div className="space-y-6">
              {bothActive ? (
                <div className="space-y-2">
                   <div className="flex items-center justify-between">
                     <Label htmlFor="rightTempSlider">Cold Temp (°C): {Number(coldTemp).toFixed(2)}</Label>
                     <div className="flex items-center gap-1">
                       <Label className="text-xs text-muted-foreground">Min:</Label>
                       <Input type="text" value={tempMin} onChange={(e)=>setTempMin(e.target.value)} className="w-16 h-8 text-xs" />
                       <Label className="text-xs text-muted-foreground">Max:</Label>
                       <Input type="text" value={tempMax} onChange={(e)=>setTempMax(e.target.value)} className="w-16 h-8 text-xs" />
                     </div>
                   </div>
                   <Slider
                     id="rightTempSlider"
                     min={parseFloat(tempMin)}
                     max={parseFloat(tempMax)}
                     step={computeStep(parseFloat(tempMax))}
                     value={[Number(coldTemp)]}
                     onValueChange={([v]) => updateTemperature('cold', v)}
                     className="w-full"
                   />
                 </div>
               ) : (
                <div className="space-y-2">
                   <div className="flex items-center justify-between">
                     <Label htmlFor="heatFlowSlider">Heat Flow (W): {Number(heatFlow).toFixed(getDecimalPlaces(heatFlow))}</Label>
                     <div className="flex items-center gap-1">
                       <Label className="text-xs text-muted-foreground">Max:</Label>
                       <Input type="text" value={heatFlowMax} onChange={(e) => setHeatFlowMax(e.target.value)} className="w-20 h-8 text-xs" />
                     </div>
                   </div>
                   <Slider id="heatFlowSlider" min={0} max={parseFloat(heatFlowMax)} step={computeStep(parseFloat(heatFlowMax))} value={[Number(heatFlow)]} onValueChange={([v]) => setHeatFlow(v)} className="w-full" />
                 </div>
               )}

              {/* Area Slider */}
              <div className="space-y-2">
                <div className="flex items-center justify-between">
                  <Label htmlFor="areaSlider">Area (m²): {Number(area).toFixed(getDecimalPlaces(area))}</Label>
                  <div className="flex items-center gap-1">
                    <Label className="text-xs text-muted-foreground">Max:</Label>
                    <Input type="text" value={areaMax} onChange={(e) => setAreaMax(e.target.value)} className="w-16 h-8 text-xs" />
                  </div>
                </div>
                <Slider id="areaSlider" min={0} max={parseFloat(areaMax)} step={computeStep(parseFloat(areaMax))} value={[Number(area)]} onValueChange={([v]) => setArea(v)} className="w-full" />
              </div>
            </div>
          </div>

          <div className={gridClassName}>
            {/* Column 1 */}
            <div className="flex flex-col gap-4">
              {layers.slice(0, 3).map((layer, localIndex) => {
                const globalIndex = localIndex;
                return (
                  <LayerItem
                    key={layer.id}
                    layer={layer}
                    index={globalIndex}
                    isFirst={globalIndex === 0}
                    isLast={globalIndex === layers.length - 1}
                    updateLayer={updateLayer}
                    deleteLayer={deleteLayer}
                    moveLayer={moveLayer}
                  />
                );
              })}
            </div>

            {/* Column 2 */}
            {layers.length > 3 && (
              <div className="flex flex-col gap-4">
                {layers.slice(3, 6).map((layer, localIndex) => {
                  const globalIndex = localIndex + 3;
                  return (
                    <LayerItem
                      key={layer.id}
                      layer={layer}
                      index={globalIndex}
                      isFirst={globalIndex === 0}
                      isLast={globalIndex === layers.length - 1}
                      updateLayer={updateLayer}
                      deleteLayer={deleteLayer}
                      moveLayer={moveLayer}
                    />
                  );
                })}
              </div>
            )}

            {/* Column 3 */}
            {layers.length > 6 && (
              <div className="flex flex-col gap-4">
                {layers.slice(6, 9).map((layer, localIndex) => {
                  const globalIndex = localIndex + 6;
                  return (
                    <LayerItem
                      key={layer.id}
                      layer={layer}
                      index={globalIndex}
                      isFirst={globalIndex === 0}
                      isLast={globalIndex === layers.length - 1}
                      updateLayer={updateLayer}
                      deleteLayer={deleteLayer}
                      moveLayer={moveLayer}
                    />
                  );
                })}
              </div>
            )}
          </div>

          <Button onClick={addLayer} disabled={layers.length >= 9}>Add Layer</Button>
        </CardContent>
      </Card>
    </div>
  );
}

const LayerItem: FC<{
  layer: Layer;
  index: number;
  isFirst: boolean;
  isLast: boolean;
  updateLayer: (id: string, newValues: Partial<Layer>) => void;
  deleteLayer: (id: string) => void;
  moveLayer: (index: number, direction: 'up' | 'down') => void;
}> = ({ layer, index, isFirst, isLast, updateLayer, deleteLayer, moveLayer }) => {
  return (
     <motion.div layout>
        <div className="flex items-center space-x-4 p-4 border rounded-lg bg-background">
          <div className="flex flex-col items-center space-y-2">
            <Button
              variant="ghost"
              size="icon"
              disabled={isFirst}
              onClick={() => moveLayer(index, 'up')}
              className="h-6 w-6"
            >
              <ArrowUp className={`h-4 w-4 ${isFirst ? 'text-muted-foreground' : ''}`} />
            </Button>
            <div className="font-bold text-lg">{index + 1}</div>
            <Button
              variant="ghost"
              size="icon"
              disabled={isLast}
              onClick={() => moveLayer(index, 'down')}
              className="h-6 w-6"
            >
              <ArrowDown className={`h-4 w-4 ${isLast ? 'text-muted-foreground' : ''}`} />
            </Button>
          </div>
          <div className="flex-grow flex items-stretch space-x-4">
             <VerticalSwitch 
                isSolid={layer.type === 'solid'}
                onTypeChange={(type) => updateLayer(layer.id, { type })}
             />
             <div className="w-full space-y-3 flex flex-col justify-center">
               <div className="w-full space-y-2">
                 <Label>Thickness (m): {layer.thickness.toFixed(3)}</Label>
                 <Slider
                   min={0.001} max={0.5} step={0.001}
                   value={[layer.thickness]}
                   onValueChange={([val]) => updateLayer(layer.id, { thickness: val })}
                 />
               </div>
                {layer.type === 'solid' ? (
                   <div className="w-full space-y-2">
                     <Label>Thermal Conductivity (k): {layer.kValue} W/m·K</Label>
                     <Slider
                       min={1} max={500} step={1}
                       value={[layer.kValue]}
                       onValueChange={([val]) => updateLayer(layer.id, { kValue: val })}
                     />
                   </div>
                 ) : (
                   <div className="w-full space-y-2">
                     <Label>Convection Coefficient (h): {layer.hValue} W/m²·K</Label>
                     <Slider
                       min={1} max={1000} step={1}
                       value={[layer.hValue]}
                       onValueChange={([val]) => updateLayer(layer.id, { hValue: val })}
                     />
                   </div>
                 )}
             </div>
          </div>
          <div className="ml-auto">
            <Button
              variant="ghost"
              size="icon"
              disabled={isFirst}
              onClick={() => !isFirst && deleteLayer(layer.id)}
              className={isFirst ? "cursor-not-allowed" : ""}
            >
              <Trash2 className={`h-4 w-4 ${isFirst ? 'text-muted-foreground' : 'text-destructive'}`} />
            </Button>
          </div>
        </div>
      </motion.div>
  )
}
