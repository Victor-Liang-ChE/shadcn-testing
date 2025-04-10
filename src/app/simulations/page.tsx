'use client';

import Link from "next/link";
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";

export default function Page() {
  const simulations = [
    { 
      name: "Reaction Kinetics", 
      path: "/simulations/kinetics",
      description: "Interactive simulator for chemical reaction kinetics. Model various reaction types and visualize concentration profiles over time."
    },
    { 
      name: "McCabe-Thiele", 
      path: "/simulations/mccabe-thiele",
      description: "Select components and specify operating conditions to visualize distillation processes with accurate equilibrium diagrams."
    },
    { 
      name: "Process Control", 
      path: "/simulations/process-control",
      description: "Simulate process control systems with various inputs and understand system dynamics in chemical processes."
    },
    { 
      name: "PID Tuning", 
      path: "/simulations/pid-tuning",
      description: "Interactive PID controller tuning simulation. Adjust parameters and observe system response in real-time."
    }
  ];

  return (
    <main className="flex flex-col items-center justify-between pt-24 px-6 md:px-24">
      <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
        {simulations.map((simulation, index) => (
          <Link key={index} href={simulation.path} className="block h-full">
            <Card className="h-full hover:shadow-lg transition-shadow bg-card border border-border">
              <CardHeader>
                <div className="h-32 w-full bg-muted rounded-md mb-4 flex items-center justify-center text-muted-foreground">{simulation.name}</div>
                <CardTitle>{simulation.name}</CardTitle>
              </CardHeader>
              <CardContent>
                <CardDescription>{simulation.description}</CardDescription>
              </CardContent>
            </Card>
          </Link>
        ))}
      </div>
    </main>
  );
}
