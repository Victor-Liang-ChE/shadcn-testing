"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";
import { ThemeToggle } from "@/components/ThemeToggle";
import { Button } from "@/components/ui/button";
import { useState } from "react";
import { Menu, X } from "lucide-react";

function Navbar() {
  const [isMenuOpen, setIsMenuOpen] = useState(false);
  const pathname = usePathname();
  
  const getPageTitle = () => {
    // Use startsWith for nested routes
    if (pathname.startsWith("/LearnChemE/fluid-dynamics/bernoulli-equation")) {
      return "Bernoulli Equation";
    }
    if (pathname.startsWith("/LearnChemE/fluid-dynamics/couette-flow")) {
      return "Couette Flow";
    }
    if (pathname.startsWith("/LearnChemE/fluid-dynamics/buoyancy")) {
      return "Buoyancy";
    }
    if (pathname.startsWith("/LearnChemE")) {
      return "LearnChemE 2.0";
    }

    switch (pathname) {      
      
      case "/toolbox":
        return "Chemical Engineering Toolbox";

      case "/toolbox/mccabe-thiele":
        return "McCabe-Thiele Graphical Method";

      case "/toolbox/kinetics":
        return "Reaction Kinetics Simulator";

      case "/toolbox/process-dynamics":
        return "Process Dynamics Simulator";

      case "/toolbox/pid-tuning":
        return "PID Tuning Simulator";

      case "/toolbox/azeotrope-finder":
        return "Azeotrope Finder";

      case "/toolbox/compound-properties":
        return "Compound Properties";

      case "/toolbox/residue-curve-map":
        return "Residue Curve Map";

      case "/toolbox/reactor-design":
        return "Reactor Design Simulator";

      case "/toolbox/isotherms":
        return "Isotherm Simulator";

      case "/toolbox/binary-phase-diagrams":
        return "Binary Phase Diagram";

      case "/toolbox/unary-phase-diagrams":
        return "Unary Phase Diagram";

      case "/toolbox/1d-heat-transfer":
        return "1D Heat Transfer";

      case "/toolbox/molecular-dynamics":
        return "Molecular Dynamics Simulator";

      case "/toolbox/FUG":
        return "FUG";

      case "/misc":
          return "Miscellaneous Projects";
  
      case "/misc/drop-chance":
          return "Drop Chance Calculator";

      case "/misc/portola-menu":
          return "Portola Menu";

      case "/misc/latex":
          return "LaTeX Constructor & Converter";
        
      case "/misc/chemistry-tools":
          return "Chemistry Tools";

      case "/misc/cheme-econ":
          return "Chemical Engineering Economics Calculator";

      case "/misc/japanese-lyrics":
          return "Japanese Lyrics Analyzer";

      case "/misc/ufc-champions":
          return "UFC Championship Lineage";

      case "/misc/csv-to-plot":
          return "Fast .csv Plotter";

      default:
        // Check for base paths if no specific match
        if (pathname === "/LearnChemE") return "LearnChemE 2.0";
        if (pathname === "/toolbox") return "Chemical Engineering Toolbox";
        if (pathname === "/misc") return "Miscellaneous Projects";
        return ""; // Default empty title
    }
  };

  // Consistent hover style for links and button - Added padding and background hover
  const linkHoverStyle = "hover:text-accent-foreground hover:bg-accent rounded-md px-3 py-2 transition-colors duration-200";
  // Style specifically for the brand link - Removed background hover, kept text hover
  const brandLinkHoverStyle = "hover:text-accent-foreground rounded-md px-2 py-1 transition-colors duration-200";


  return (
    // Use inline style to directly reference CSS variables
    <nav
      className="flex items-center justify-between p-4 border-b min-w-full text-navbar-foreground"
      style={{ backgroundColor: 'var(--navbar-background)' }}
    >
      {/* Left: Brand */}
      <div>
        <Link href="/" className={`${brandLinkHoverStyle} inline-block`}>
          <img 
            src="/favicon.ico" 
            alt="Victor Liang" 
            className="w-8 h-8"
          />
        </Link>
      </div>
      {/* Center: Page Title */}
      <div className="absolute left-1/2 transform -translate-x-1/2">
        {/* Increased font size */}
        <span className="text-xl font-semibold">{getPageTitle()}</span>
      </div>
      {/* Mobile Menu Button */}
      <div className="md:hidden">
        <Button
          variant="ghost"
          size="icon"
          onClick={() => setIsMenuOpen(!isMenuOpen)}
          // Use navbar-foreground for icon color, accent for hover
          className={`text-navbar-foreground ${linkHoverStyle}`}
        >
          {isMenuOpen ? <X size={24} /> : <Menu size={24} />}
        </Button>
      </div>
      {/* Desktop Menu */}
      <div className="hidden md:flex items-center space-x-1"> {/* Reduced space for tighter fit with padding */}
        {/* Remove legacyBehavior, passHref, and nested <a>. Apply styles directly to Link. */}
        <Link href="/toolbox" className={`text-xl font-semibold ${linkHoverStyle}`}>
          Toolbox
        </Link>
        {/* Remove legacyBehavior, passHref, and nested <a>. Apply styles directly to Link. */}
        <Link href="/misc" className={`text-xl font-semibold ${linkHoverStyle}`}>
          Misc
        </Link>
        {/* Apply linkHoverStyle to the wrapper div */}
        <div className={linkHoverStyle}>
           <ThemeToggle />
        </div>
      </div>
      {/* Mobile Menu Dropdown */}
      <div
        className={`absolute top-16 right-0 w-48 z-50 md:hidden border-border overflow-hidden transition-all duration-200 ease-in-out ${
          isMenuOpen ? 'opacity-100 translate-y-0' : 'opacity-0 -translate-y-2 pointer-events-none'
        }`}
        style={{ backgroundColor: 'var(--navbar-background)' }}
      >
        <div className="p-4 flex flex-col space-y-2"> {/* Reduced space-y */}
           {/* Remove legacyBehavior, passHref, and nested <a>. Apply styles directly to Link. */}
           <Link href="/toolbox" className={`block text-xl font-semibold ${linkHoverStyle}`}>
             Toolbox
           </Link>
           {/* Remove legacyBehavior, passHref, and nested <a>. Apply styles directly to Link. */}
           <Link href="/misc" className={`block text-xl font-semibold ${linkHoverStyle}`}>
             Misc
           </Link>
          {/* Center the toggle button and apply linkHoverStyle to the wrapper */}
          <div className={`flex justify-center pt-2 ${linkHoverStyle}`}> {/* Added padding top and hover style */}
            <ThemeToggle />
          </div>
        </div>
      </div>
    </nav>
  );
}

export default Navbar;
