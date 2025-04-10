"use client";

import { useState, useEffect } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Slider } from "@/components/ui/slider"; // Import Slider

const DropChanceCalculator = () => {
  const [percentValue, setPercentValue] = useState<number>(50); // Default to 50 for slider
  const [attemptsValue, setAttemptsValue] = useState<number>(10); // Default to 10 for slider
  const [maxAttemptsValue, setMaxAttemptsValue] = useState<number>(100);
  const [desiredDropsValue, setDesiredDropsValue] = useState<number>(1);
  const [desiredDropsRangeValue, setDesiredDropsRangeValue] = useState<[number, number]>([1, 2]);
  const [toggleState, setToggleState] = useState<number>(0);
  const [result, setResult] = useState<string>('');

  // Binomial coefficient function with additional safeguards for large numbers
  const binomialCoefficient = (n: number, k: number): number => {
    if (k < 0 || k > n) return 0;
    if (k === 0 || k === n) return 1;
    if (n <= 0) return 0;

    if (n > 1000 || k > 100) {
      let logResult = 0;
      for (let i = n - k + 1; i <= n; i++) {
        logResult += Math.log(i);
      }
      for (let i = 1; i <= k; i++) {
        logResult -= Math.log(i);
      }
      return Math.exp(logResult);
    }
    
    let coeff = 1;
    for (let x = n - k + 1; x <= n; x++) coeff *= x;
    for (let x = 1; x <= k; x++) coeff /= x;
    return coeff;
  };

  const binomialProbability = (n: number, k: number, p: number): number => {
    return binomialCoefficient(n, k) * Math.pow(p, k) * Math.pow(1 - p, n - k);
  };

  const calculateCumulativeProbability = (n: number, k: number, p: number, isAtLeast: boolean): number => {
    if (n > 1000) {
      const mean = n * p;
      const stdDev = Math.sqrt(n * p * (1 - p));
      
      if (isAtLeast) {
        const z = (k - 0.5 - mean) / stdDev;
        return 1 - normalCDF(z);
      } else {
        const z1 = (k - 0.5 - mean) / stdDev;
        const z2 = (k + 0.5 - mean) / stdDev;
        return normalCDF(z2) - normalCDF(z1);
      }
    }
    
    if (isAtLeast) {
      let cumProb = 0;
      for (let i = 0; i < k; i++) {
        cumProb += binomialProbability(n, i, p);
      }
      return 1 - cumProb;
    } else {
      return binomialProbability(n, k, p);
    }
  };
  
  const normalCDF = (z: number): number => {
    if (z < -6) return 0;
    if (z > 6) return 1;
    
    let sum = 0;
    let term = z;
    for (let i = 3; sum + term !== sum; i += 2) {
      sum += term;
      term = term * z * z / i;
    }
    
    return 0.5 + sum * Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI);
  };

  useEffect(() => {
    // Cap desiredDropsValue based on attemptsValue
    if (desiredDropsValue > attemptsValue) {
      setDesiredDropsValue(attemptsValue);
    }
    // Cap desiredDropsRangeValue based on attemptsValue
    if (desiredDropsRangeValue[0] > attemptsValue) {
        setDesiredDropsRangeValue([attemptsValue, Math.max(attemptsValue, desiredDropsRangeValue[1])]);
    }
     if (desiredDropsRangeValue[1] > attemptsValue) {
        setDesiredDropsRangeValue([desiredDropsRangeValue[0], attemptsValue]);
    }
  }, [attemptsValue, desiredDropsValue, desiredDropsRangeValue]);


  useEffect(() => {
    if (isNaN(percentValue) || isNaN(attemptsValue) || 
        isNaN(desiredDropsValue) || 
        isNaN(desiredDropsRangeValue[0]) || isNaN(desiredDropsRangeValue[1])) {
      setResult("Please enter valid values for all fields.");
      return;
    }

    const probability = percentValue / 100;
    let resultText = '';

    if (toggleState === 0) {
       // Ensure desiredDropsValue is not greater than attemptsValue before calculation
       const currentDesiredDrops = Math.min(desiredDropsValue, attemptsValue);
       if (currentDesiredDrops < desiredDropsValue) {
           setDesiredDropsValue(currentDesiredDrops); // Adjust state if needed
       }
      const exactProb = calculateCumulativeProbability(attemptsValue, currentDesiredDrops, probability, false);
      const atLeastProb = calculateCumulativeProbability(attemptsValue, currentDesiredDrops, probability, true);

      resultText = `Probability of getting exactly ${currentDesiredDrops} drops: ${(exactProb * 100).toFixed(4)}%\n`;
      resultText += `Probability of getting at least ${currentDesiredDrops} drops: ${(atLeastProb * 100).toFixed(4)}%`;
    } else {
       // Ensure range values are valid and capped
       const minDrops = Math.min(desiredDropsRangeValue[0], attemptsValue);
       const maxDrops = Math.min(desiredDropsRangeValue[1], attemptsValue);
       if (minDrops !== desiredDropsRangeValue[0] || maxDrops !== desiredDropsRangeValue[1]) {
           setDesiredDropsRangeValue([minDrops, maxDrops]); // Adjust state if needed
       }

       // Ensure minDrops is not greater than maxDrops after capping
       const finalMinDrops = Math.min(minDrops, maxDrops);

       // Calculate probability for the valid range [finalMinDrops, maxDrops]
       // P(min <= X <= max) = P(X <= max) - P(X < min) = P(X <= max) - P(X <= min - 1)
       // Using the 'at least' logic: P(min <= X <= max) = P(X >= min) - P(X > max) = P(X >= min) - P(X >= max + 1)
       const probAtLeastMin = calculateCumulativeProbability(attemptsValue, finalMinDrops, probability, true);
       const probAtLeastMaxPlusOne = calculateCumulativeProbability(attemptsValue, maxDrops + 1, probability, true);
       const rangeProb = probAtLeastMin - probAtLeastMaxPlusOne;


      resultText = `Probability of getting between ${finalMinDrops} and ${maxDrops} drops: ${(rangeProb * 100).toFixed(4)}%`;
    }

    setResult(resultText);
  }, [percentValue, attemptsValue, desiredDropsValue, desiredDropsRangeValue, toggleState]);

   // Handle max attempts input change with upper limit
  const handleMaxAttemptsChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const newMax = Number(e.target.value);
    if (newMax > 0 && newMax <= 1000000) { // Set a reasonable upper limit
      setMaxAttemptsValue(newMax);
      // If current attempts value is larger than newMax, adjust it
      if (attemptsValue > newMax) {
        setAttemptsValue(newMax);
      }
    }
  };

  // Handle range slider changes, ensuring min <= max
  const handleRangeChange = (index: number, value: number) => {
    const newRange = [...desiredDropsRangeValue] as [number, number];
    newRange[index] = value;

    if (index === 0 && newRange[0] > newRange[1]) {
      newRange[1] = newRange[0]; // If min exceeds max, set max to min
    } else if (index === 1 && newRange[1] < newRange[0]) {
      newRange[0] = newRange[1]; // If max goes below min, set min to max
    }

    setDesiredDropsRangeValue(newRange);
  };


  return (
    <main className="flex flex-col items-center justify-between pt-12 md:pt-24 px-6 md:px-24"> {/* Reduced top padding */}
      <Card className="w-full max-w-2xl">
        <CardHeader>
          <CardTitle>Drop Chance Calculator</CardTitle>
        </CardHeader>
        <CardContent className="space-y-6">
          <Tabs defaultValue="exact" onValueChange={(value) => setToggleState(value === "exact" ? 0 : 1)}>
            <TabsList className="grid w-full grid-cols-2">
              <TabsTrigger value="exact">Exact/At Least</TabsTrigger>
              <TabsTrigger value="range">Range</TabsTrigger>
            </TabsList>

            {/* Common Inputs */}
            <div className="space-y-4 pt-4">
               <div className="space-y-2">
                 <Label htmlFor="percent">Drop Chance: {percentValue}%</Label>
                 <Slider
                   id="percent"
                   min={0}
                   max={100}
                   step={0.1} // Finer step for percentage
                   value={[percentValue]}
                   onValueChange={(value) => setPercentValue(value[0])}
                 />
               </div>
               <div className="space-y-2">
                 <div className="flex justify-between items-center">
                    <Label htmlFor="attempts">Number of Attempts: {attemptsValue}</Label>
                    <div className="flex items-center gap-2">
                       <Label htmlFor="max-attempts" className="text-sm text-muted-foreground whitespace-nowrap">Max:</Label>
                       <Input
                         id="max-attempts"
                         type="number"
                         min="1"
                         max="1000000" // Set a reasonable upper limit
                         value={maxAttemptsValue}
                         onChange={handleMaxAttemptsChange}
                         className="h-8 w-24" // Smaller input for max
                       />
                    </div>
                 </div>
                 <Slider
                   id="attempts"
                   min={0}
                   max={maxAttemptsValue}
                   step={1}
                   value={[attemptsValue]}
                   onValueChange={(value) => setAttemptsValue(value[0])}
                 />
               </div>
            </div>

            <TabsContent value="exact" className="space-y-4 pt-4 border-t mt-4">
               <div className="space-y-2">
                 <Label htmlFor="drops">Desired Drops (At Least): {desiredDropsValue}</Label>
                 <Slider
                   id="drops"
                   min={0}
                   max={attemptsValue} // Dynamic max based on attempts
                   step={1}
                   value={[desiredDropsValue]}
                   onValueChange={(value) => setDesiredDropsValue(value[0])}
                   disabled={attemptsValue === 0} // Disable if 0 attempts
                 />
               </div>
            </TabsContent>
            <TabsContent value="range" className="space-y-4 pt-4 border-t mt-4">
               <div className="space-y-2">
                 <Label>Desired Drops Range: {desiredDropsRangeValue[0]} to {desiredDropsRangeValue[1]}</Label>
                 <div className="grid grid-cols-2 gap-4">
                    <div className="space-y-1">
                       <Label htmlFor="min-drops" className="text-sm">Min: {desiredDropsRangeValue[0]}</Label>
                       <Slider
                         id="min-drops"
                         min={0}
                         max={attemptsValue} // Dynamic max
                         step={1}
                         value={[desiredDropsRangeValue[0]]}
                         onValueChange={(value) => handleRangeChange(0, value[0])}
                         disabled={attemptsValue === 0}
                       />
                    </div>
                    <div className="space-y-1">
                       <Label htmlFor="max-drops" className="text-sm">Max: {desiredDropsRangeValue[1]}</Label>
                       <Slider
                         id="max-drops"
                         min={0}
                         max={attemptsValue} // Dynamic max
                         step={1}
                         value={[desiredDropsRangeValue[1]]}
                         onValueChange={(value) => handleRangeChange(1, value[0])}
                         disabled={attemptsValue === 0}
                       />
                    </div>
                 </div>
               </div>
            </TabsContent>
          </Tabs>

          <div className="mt-6 p-4 bg-muted rounded-md whitespace-pre-line text-sm"> {/* Smaller text for result */}
            {result || "Adjust sliders to calculate probability"}
          </div>
        </CardContent>
      </Card>
    </main>
  );
}

export default DropChanceCalculator;
