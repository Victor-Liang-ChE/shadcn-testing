'use client';

import Link from "next/link";
import Image from "next/image"; // Import Image component
import { Card, CardContent, CardHeader, CardTitle, CardDescription } from "@/components/ui/card";

// Define the path for the thumbnail
const mccabeThumbnailPath = "/thumbnails/mccabe-thiele-thumbnail.png";

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
    <main className="flex flex-col items-center justify-between pt-24 px-6">
      <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 w-full"> 
        {simulations.map((simulation, index) => (
          <Link key={index} href={simulation.path} className="block h-full">
            <Card className="h-full hover:shadow-lg transition-shadow bg-card border border-border flex flex-col"> 
              <CardHeader>
                {/* Outer container for centering the block */}
                <div className="relative mx-auto mb-4"> 
                  {/* Inner container for aspect ratio, size, rounding, clipping, and relative positioning */}
                  <div className="aspect-square w-64 h-64 rounded-3xl overflow-hidden relative"> 
                    {simulation.name === "McCabe-Thiele" ? (
                      <Image 
                        src={mccabeThumbnailPath} 
                        alt={`${simulation.name} Thumbnail`} 
                        layout="fill" 
                        objectFit="contain" 
                      />
                    ) : (
                      // Placeholder text - absolutely positioned to center within the inner container
                      <div className="absolute inset-0 flex items-center justify-center text-muted-foreground bg-muted"> {/* Optional: bg-muted for placeholder */}
                        <span>{simulation.name}</span>
                      </div>
                    )}
                  </div>
                </div>
                <CardTitle>{simulation.name}</CardTitle>
              </CardHeader>
              <CardContent className="flex-grow"> 
                <CardDescription>{simulation.description}</CardDescription>
              </CardContent>
            </Card>
          </Link>
        ))}
      </div>
    </main>
  );
}
