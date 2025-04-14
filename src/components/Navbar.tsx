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
    switch (pathname) {      
      
      case "/simulations":
        return "Chemical Engineering Simulations";

      case "/simulations/mccabe-thiele":
        return "McCabe-Thiele Graphical Method";

      case "/simulations/kinetics":
        return "Reaction Kinetics Simulator";

      case "/simulations/process-control":
        return "Process Control Simulator";

      case "/simulations/pid-tuning":
        return "PID Tuning Simulator";

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

      default:
        return "Welcome";
    }
  };

  // Consistent hover style for links and button
  const linkHoverStyle = "hover:text-accent-foreground transition-colors duration-200";

  return (
    // Use inline style to directly reference CSS variables
    <nav
      className="flex items-center justify-between p-4 border-b min-w-full text-navbar-foreground"
      style={{ backgroundColor: 'var(--navbar-background)' }}
    >
      {/* Left: Brand */}
      <div>
        <Link href="/" legacyBehavior passHref>
          {/* Increased font size */}
          <a className={`text-xl font-semibold ${linkHoverStyle} inline-block`}>
            Victor Liang
          </a>
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
      <div className="hidden md:flex items-center space-x-4">
        <Link href="/simulations" legacyBehavior passHref>
          {/* Increased font size */}
          <a className={`text-xl font-semibold ${linkHoverStyle}`}>Simulations</a>
        </Link>
        <Link href="/misc" legacyBehavior passHref>
          {/* Increased font size */}
          <a className={`text-xl font-semibold ${linkHoverStyle}`}>Misc</a>
        </Link>
        {/* Apply consistent hover to ThemeToggle's button */}
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
        <div className="p-4 flex flex-col space-y-4 border-x border-b">
          <Link href="/simulations" legacyBehavior passHref>
            {/* Increased font size */}
            <a className={`text-xl font-semibold ${linkHoverStyle}`}>Simulations</a>
          </Link>
          <Link href="/misc" legacyBehavior passHref>
            {/* Increased font size */}
            <a className={`text-xl font-semibold ${linkHoverStyle}`}>Misc</a>
          </Link>
          <div className={`flex justify-center ${linkHoverStyle}`}>
            <ThemeToggle />
          </div>
        </div>
      </div>
    </nav>
  );
}

export default Navbar;
