"use client";

import Link from "next/link";
import { ThemeToggle } from "@/components/ThemeToggle";
import { Button } from "@/components/ui/button";
import { useState } from "react";
import { Menu, X } from "lucide-react";

function Navbar() {
  const [isMenuOpen, setIsMenuOpen] = useState(false);
  const pageTitle = "Home"; // Placeholder

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
          {/* Removed font-serif, apply consistent hover */}
          <a className={`text-lg font-semibold ${linkHoverStyle} inline-block`}>
            Victor Liang
          </a>
        </Link>
      </div>

      {/* Center: Page Title */}
      <div className="absolute left-1/2 transform -translate-x-1/2">
        {/* Removed font-semibold, rely on base styles */}
        <span className="text-lg">{pageTitle}</span>
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
          {/* Removed font-serif, apply consistent hover */}
          <a className={`text-lg font-semibold ${linkHoverStyle}`}>Simulations</a>
        </Link>
        <Link href="/misc" legacyBehavior passHref>
          {/* Removed font-serif, apply consistent hover */}
          <a className={`text-lg font-semibold ${linkHoverStyle}`}>Misc</a>
        </Link>
        {/* Apply consistent hover to ThemeToggle's button */}
        <div className={linkHoverStyle}>
           <ThemeToggle />
        </div>
      </div>

      {/* Mobile Menu Dropdown */}
      {isMenuOpen && (
        // Apply background color directly with CSS variable
        <div 
          className="absolute top-16 right-0 left-0 z-50 md:hidden p-4 flex flex-col space-y-4 border-b border-border"
          style={{ backgroundColor: 'var(--navbar-background)' }}
        >
          <Link href="/simulations" legacyBehavior passHref>
            <a className={`text-lg font-semibold ${linkHoverStyle}`}>Simulations</a>
          </Link>
          <Link href="/misc" legacyBehavior passHref>
            <a className={`text-lg font-semibold ${linkHoverStyle}`}>Misc</a>
          </Link>
          <div className={`flex justify-center ${linkHoverStyle}`}>
            <ThemeToggle />
          </div>
        </div>
      )}
    </nav>
  );
}

export default Navbar;
