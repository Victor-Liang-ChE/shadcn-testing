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
    <div className={`bg-card rounded flex items-center justify-center ${className}`}>
      <svg viewBox="0 0 400 200" className="w-full h-full">
        {/* CSTR Tank - even larger rectangular shape */}
        <rect x="80" y="5" width="240" height="220" fill="lightblue" stroke="currentColor" strokeWidth="4" rx="25"/>
        
        {/* Remove Feed and Product arrows and labels */}
        {/* ... */}
        
        {/* Stirrer shaft */}
        <line x1="200" y1="-40" x2="200" y2="150" stroke="currentColor" strokeWidth="8"/>
        
        {/* Stirrer blades with spinning animation - smaller */}
        <g className="animate-spin" style={{transformOrigin: '200px 90px'}}>
          <line x1="140" y1="90" x2="260" y2="90" stroke="currentColor" strokeWidth="8"/>
          <line x1="200" y1="50" x2="200" y2="130" stroke="currentColor" strokeWidth="8"/>
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
