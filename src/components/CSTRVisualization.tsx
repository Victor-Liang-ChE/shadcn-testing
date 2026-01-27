import React from 'react';

interface CSTRVisualizationProps {
  className?: string;
  showLabel?: boolean;
}

export const CSTRVisualization: React.FC<CSTRVisualizationProps> = ({ 
  className = "w-full h-full", 
  showLabel = true 
}) => {
  return (
    <div className={`bg-card rounded flex items-center justify-center ${className}`} style={{perspective: '1000px'}}>
      <svg viewBox="0 0 400 200" className="w-full h-full">
        {/* CSTR Tank - even larger rectangular shape */}
        <rect x="80" y="5" width="240" height="220" fill="lightblue" stroke="currentColor" strokeWidth="4" rx="25"/>
        
        {/* Remove Feed and Product arrows and labels */}
        {/* ... */}
        
        {/* Stirrer shaft */}
        <line x1="200" y1="-40" x2="200" y2="150" stroke="currentColor" strokeWidth="8"/>
        
        {/* Stirrer blades with spinning animation - smaller */}
        <g style={{transformOrigin: '200px 150px', animation: 'spinVertical 2s linear infinite'}}>
          {/* A smaller horizontal blade centered at (200, 150) */}
          <line x1="160" y1="150" x2="240" y2="150" stroke="currentColor" strokeWidth="8"/>
          {/* A vertical blade of the same length, also centered at (200, 150) */}
          <line x1="200" y1="110" x2="200" y2="150" stroke="currentColor" strokeWidth="8"/>
        </g>
        
        {/* CSTR Label */}
        {showLabel && (
          <text x="200" y="215" fontSize="28" fill="currentColor" className="text-foreground" textAnchor="middle" fontWeight="bold">CSTR</text>
        )}
      </svg>
    </div>
  );
};

export default CSTRVisualization;
