'use client';

import React, { useState } from 'react';
import Link from 'next/link';
import { ScrollArea } from "@/components/ui/scroll-area";
import { Button } from "@/components/ui/button";
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { cn } from "@/lib/utils";
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faBarsStaggered, faBars } from '@fortawesome/free-solid-svg-icons';

type Category =
  | "Fluid Dynamics"
  | "Heat Transfer"
  | "Interactive Diagram Builders"
  | "Kinetics and Reactor Design"
  | "Material and Energy Balances"
  | "Materials Science"
  | "Physical Chemistry"
  | "Process Control"
  | "Separations / Mass Transfer"
  | "Statics"
  | "Statistics"
  | "Thermodynamics";

const categories: Category[] = [
  "Fluid Dynamics",
  "Heat Transfer",
  "Interactive Diagram Builders",
  "Kinetics and Reactor Design",
  "Material and Energy Balances",
  "Materials Science",
  "Physical Chemistry",
  "Process Control",
  "Separations / Mass Transfer",
  "Statics",
  "Statistics",
  "Thermodynamics",
];

// Placeholder component for categories without specific content yet
const PlaceholderContent = ({ category }: { category: Category }) => (
  <div className="p-6">
    <h2 className="text-2xl font-semibold mb-4">{category}</h2>
    <p>Content for {category} will be available soon.</p>
  </div>
);

// Component for Fluid Dynamics content
const FluidDynamicsContent = () => (
  <div className="p-6">
    <h2 className="text-2xl font-semibold mb-4">Fluid Dynamics Simulations</h2>
    <div className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-3 gap-4">
       <Link href="/LearnChemE/fluid-dynamics/bernoulli-equation" className="block">
         <Card className="hover:shadow-md transition-shadow h-full flex flex-col">
           <CardHeader>
             <CardTitle>Bernoulli Equation</CardTitle>
           </CardHeader>
           <CardContent className="flex-grow">
             <p className="text-sm text-muted-foreground">
               Interactive simulation demonstrating the principles of the Bernoulli equation for fluid flow.
             </p>
           </CardContent>
         </Card>
       </Link>
       <Link href="/LearnChemE/fluid-dynamics/buoyancy" className="block">
         <Card className="hover:shadow-md transition-shadow h-full flex flex-col">
           <CardHeader>
             <CardTitle>Buoyancy of a Floating Cube</CardTitle>
           </CardHeader>
           <CardContent className="flex-grow">
             <p className="text-sm text-muted-foreground">
               Explore the principles of buoyancy with an interactive floating cube simulation.
             </p>
           </CardContent>
         </Card>
       </Link>
       {/* Add more Fluid Dynamics simulations here as cards */}
    </div>
  </div>
);


export default function LearnChemEPage() {
  const [selectedCategory, setSelectedCategory] = useState<Category>("Fluid Dynamics");
  const [isSidebarCollapsed, setIsSidebarCollapsed] = useState(false); // State for collapse

  const renderContent = () => {
    switch (selectedCategory) {
      case "Fluid Dynamics":
        return <FluidDynamicsContent />;
      // Add cases for other categories here when content is ready
      // case "Heat Transfer":
      //   return <HeatTransferContent />;
      default:
        return <PlaceholderContent category={selectedCategory} />;
    }
  };

  return (
    <main className="flex h-[calc(100vh-var(--navbar-height,4rem))]"> {/* Adjust height calculation */}
      {/* Collapsible Sidebar */}
      <div className={cn(
        "border-r bg-background transition-all duration-300 ease-in-out", // Removed relative positioning here
        isSidebarCollapsed ? "w-16" : "w-64"
      )}>
        {/* Wrap ScrollArea content to manage header separately */}
        <div className="flex flex-col h-full">
          {/* Sidebar Header */}
          <div className="flex items-center justify-between p-4 border-b">
            {/* Conditionally render title */}
            {!isSidebarCollapsed && (
              <h2 className="text-lg font-semibold">Categories</h2>
            )}
            {/* Collapse Toggle Button - Use Correct FontAwesome Icons */}
            <Button
              variant="ghost"
              size="icon"
              className={cn(
                  "h-8 w-8 border rounded-full",
                  isSidebarCollapsed && "mx-auto" // Center button when collapsed
              )}
              onClick={() => setIsSidebarCollapsed(!isSidebarCollapsed)}
            >
              {/* Use faBarsStaggered when open (to collapse), faBars when closed (to expand) */}
              <FontAwesomeIcon icon={isSidebarCollapsed ? faBars : faBarsStaggered} className="h-4 w-4" />
            </Button>
          </div>

          {/* Sidebar Content */}
          <ScrollArea className="flex-grow p-4">
            <div className="flex flex-col space-y-1">
              {categories.map((category) => (
                <Button
                  key={category}
                  variant="ghost"
                  className={cn(
                    "w-full justify-start text-left h-auto py-2 px-2",
                    // Apply selected styles only if not collapsed
                    selectedCategory === category && !isSidebarCollapsed && "bg-accent text-accent-foreground",
                    isSidebarCollapsed && "justify-center pointer-events-none", // Disable pointer events and center
                    // Remove hover effect when collapsed
                    isSidebarCollapsed ? "hover:bg-transparent" : ""
                  )}
                  // Keep onClick, but pointer-events: none will prevent it
                  onClick={() => !isSidebarCollapsed && setSelectedCategory(category)}
                  title={isSidebarCollapsed ? category : ""} // Only show tooltip when collapsed
                  // Optionally disable the button visually when collapsed
                  // disabled={isSidebarCollapsed}
                >
                  {/* Conditionally render full text or nothing */}
                  {!isSidebarCollapsed && category}
                  {/* If using icons, render them here conditionally */}
                  {/* {isSidebarCollapsed && <IconComponent size={16} />} */}
                </Button>
              ))}
            </div>
          </ScrollArea>
        </div>
      </div>

      {/* Main Content Area */}
      <div className="flex-grow"> {/* Takes remaining width */}
        <ScrollArea className="h-full">
          {renderContent()}
        </ScrollArea>
      </div>
    </main>
  );
}
