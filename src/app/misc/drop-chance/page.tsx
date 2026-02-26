"use client";

import { useState, useEffect } from 'react';
import { Card, CardContent, CardHeader } from "@/components/ui/card";
import { Input } from "@/components/ui/input";
import { Label } from "@/components/ui/label";
import { Tabs, TabsContent, TabsList, TabsTrigger } from "@/components/ui/tabs";
import { Slider } from "@/components/ui/slider";

const DropChanceCalculator = () => {
  const [percentValue, setPercentValue] = useState<number>(50);
  const [attemptsValue, setAttemptsValue] = useState<number>(10);
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
      // Exact probability calculation
      const currentDesiredDrops = Math.min(desiredDropsValue, attemptsValue);
      if (currentDesiredDrops < desiredDropsValue) { setDesiredDropsValue(currentDesiredDrops); }
      const exactProb = calculateCumulativeProbability(attemptsValue, currentDesiredDrops, probability, false);
      resultText = `Probability of getting exactly ${currentDesiredDrops} ${currentDesiredDrops === 1 ? 'drop' : 'drops'}: <b>${(exactProb * 100).toPrecision(3)}%</b>`;
    } else if (toggleState === 1) {
      // At least probability calculation
      const currentDesiredDrops = Math.min(desiredDropsValue, attemptsValue);
      if (currentDesiredDrops < desiredDropsValue) { setDesiredDropsValue(currentDesiredDrops); }
      const atLeastProb = calculateCumulativeProbability(attemptsValue, currentDesiredDrops, probability, true);
      resultText = `Probability of getting at least ${currentDesiredDrops} ${currentDesiredDrops === 1 ? 'drop' : 'drops'}: <b>${(atLeastProb * 100).toPrecision(3)}%</b>`;
    } else {
      // Range probability calculation
      const minDrops = Math.min(desiredDropsRangeValue[0], attemptsValue);
      const maxDrops = Math.min(desiredDropsRangeValue[1], attemptsValue);
      if (minDrops !== desiredDropsRangeValue[0] || maxDrops !== desiredDropsRangeValue[1]) { setDesiredDropsRangeValue([minDrops, maxDrops]); }
      const finalMinDrops = Math.min(minDrops, maxDrops);
      const probAtLeastMin = calculateCumulativeProbability(attemptsValue, finalMinDrops, probability, true);
      const probAtLeastMaxPlusOne = calculateCumulativeProbability(attemptsValue, maxDrops + 1, probability, true);
      const rangeProb = probAtLeastMin - probAtLeastMaxPlusOne;
      resultText = `Probability of getting between ${finalMinDrops} and ${maxDrops} ${maxDrops === 1 ? 'drop' : 'drops'}: <b>${(rangeProb * 100).toPrecision(3)}%</b>`;
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
    <main className="flex flex-col items-center justify-between pt-12 md:pt-24 px-8 md:px-32 pb-12"> {/* Increased horizontal padding: px-8 md:px-32 */}
      <Card className="w-full max-w-4xl">
        <CardHeader>
        </CardHeader>
        <CardContent className="space-y-8">
          <Tabs defaultValue="exact" onValueChange={(value) => {
            if (value === "exact") setToggleState(0);
            else if (value === "atleast") setToggleState(1);
            else if (value === "range") setToggleState(2);
          }}>
            <TabsList className="grid w-full grid-cols-3">
              <TabsTrigger value="exact">Exact</TabsTrigger>
              <TabsTrigger value="atleast">At Least</TabsTrigger>
              <TabsTrigger value="range">Range</TabsTrigger>
            </TabsList>

            {/* Common Inputs */}
            <div className="space-y-6 pt-6">
               <div className="space-y-3">
                 <Label htmlFor="percent" className="font-bold">Drop Chance: {percentValue}%</Label>
                 <Slider
                   id="percent"
                   min={0}
                   max={100}
                   step={0.1}
                   value={[percentValue]}
                   onValueChange={(value) => setPercentValue(value[0])}
                 />
               </div>
               <div className="space-y-3">
                 <div className="flex justify-between items-center">
                    <Label htmlFor="attempts" className="font-bold">Number of Attempts: {attemptsValue}</Label>
                    <div className="flex items-center gap-2">
                       <Label htmlFor="max-attempts" className="text-sm text-muted-foreground whitespace-nowrap font-bold">Max:</Label>
                       <Input
                         id="max-attempts"
                         type="number"
                         min="1"
                         max="1000000"
                         value={maxAttemptsValue}
                         onChange={handleMaxAttemptsChange}
                         className="h-8 w-24"
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

            <TabsContent value="exact" className="space-y-6 pt-6 mt-4">
               <div className="space-y-3">
                 <Label htmlFor="drops" className="font-bold">Desired {desiredDropsValue === 1 ? 'Drop' : 'Drops'}: {desiredDropsValue}</Label>
                 <Slider
                   id="drops"
                   min={0}
                   max={attemptsValue}
                   step={1}
                   value={[desiredDropsValue]}
                   onValueChange={(value) => setDesiredDropsValue(value[0])}
                   disabled={attemptsValue === 0}
                 />
               </div>
            </TabsContent>
            <TabsContent value="atleast" className="space-y-6 pt-6 mt-4">
               <div className="space-y-3">
                 <Label htmlFor="drops-atleast" className="font-bold">Desired {desiredDropsValue === 1 ? 'Drop' : 'Drops'}: {desiredDropsValue}</Label>
                 <Slider
                   id="drops-atleast"
                   min={0}
                   max={attemptsValue}
                   step={1}
                   value={[desiredDropsValue]}
                   onValueChange={(value) => setDesiredDropsValue(value[0])}
                   disabled={attemptsValue === 0}
                 />
               </div>
            </TabsContent>
            <TabsContent value="range" className="space-y-6 pt-6 mt-4">
               <div className="space-y-3">
                 <Label className="font-bold">Desired Range: {desiredDropsRangeValue[0]} to {desiredDropsRangeValue[1]} {desiredDropsRangeValue[1] === 1 ? 'Drop' : 'Drops'}</Label>
                 <div className="grid grid-cols-2 gap-6">
                    <div className="space-y-2">
                       <Label htmlFor="min-drops" className="text-sm font-bold">Min: {desiredDropsRangeValue[0]}</Label>
                       <Slider
                         id="min-drops"
                         min={0}
                         max={attemptsValue}
                         step={1}
                         value={[desiredDropsRangeValue[0]]}
                         onValueChange={(value) => handleRangeChange(0, value[0])}
                         disabled={attemptsValue === 0}
                       />
                    </div>
                    <div className="space-y-2">
                       <Label htmlFor="max-drops" className="text-sm font-bold">Max: {desiredDropsRangeValue[1]}</Label>
                       <Slider
                         id="max-drops"
                         min={0}
                         max={attemptsValue}
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

          <div className="mt-8 p-6 bg-muted rounded-md whitespace-pre-line">
            <div dangerouslySetInnerHTML={{ __html: result || "Adjust sliders to calculate probability" }} />
          </div>
        </CardContent>
      </Card>
    </main>
  );
}

export default DropChanceCalculator;
