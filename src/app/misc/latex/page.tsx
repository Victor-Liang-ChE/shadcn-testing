'use client';

import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Copy, Check } from 'lucide-react'; // Icons for copy button

// --- Conversion Logic (Ported from Python) ---

function replaceFracManually(expr: string): string {
    const fracIndex = expr.indexOf('\\frac');
    if (fracIndex === -1) return expr;

    let braceLevel = 0;
    let numStart = -1, numEnd = -1, denStart = -1, denEnd = -1;

    // Find numerator braces
    for (let i = fracIndex + 5; i < expr.length; i++) {
        if (expr[i] === '{') {
            if (braceLevel === 0) numStart = i;
            braceLevel++;
        } else if (expr[i] === '}') {
            braceLevel--;
            if (braceLevel === 0 && numStart !== -1) {
                numEnd = i;
                break;
            }
        }
    }

    if (numStart === -1 || numEnd === -1) return expr; // Malformed

    // Find denominator braces
    braceLevel = 0;
    for (let i = numEnd + 1; i < expr.length; i++) {
        if (expr[i] === '{') {
            if (braceLevel === 0) denStart = i;
            braceLevel++;
        } else if (expr[i] === '}') {
            braceLevel--;
            if (braceLevel === 0 && denStart !== -1) {
                denEnd = i;
                break;
            }
        }
    }

    if (denStart === -1 || denEnd === -1) return expr; // Malformed

    const numerator = expr.substring(numStart + 1, numEnd);
    const denominator = expr.substring(denStart + 1, denEnd);

    // Recursively process numerator and denominator
    const newNum = replaceFracManually(numerator);
    const newDen = replaceFracManually(denominator);

    const replacement = `(${newNum})/(${newDen})`;
    const newExpr = expr.substring(0, fracIndex) + replacement + expr.substring(denEnd + 1);

    // Process the rest of the string
    return replaceFracManually(newExpr);
}

function fullyReplaceFrac(expr: string): string {
    return replaceFracManually(expr);
}

function latex2python(latexStr: string): string {
    let expr = latexStr.replace(/^\$|\$$/g, ''); // Remove starting/ending $
    expr = expr.replace(/\\left|\\right/g, '');
    expr = fullyReplaceFrac(expr);
    expr = expr.replace(/\\operatorname\{([^{}]+)\}/g, '$1');
    expr = expr.replace(/\\?exp\(/g, 'e**(');
    expr = expr.replace(/\\(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh|log|ln)\b/g, '$1'); // Use \b for word boundary
    expr = expr.replace(/\\cdot\s*/g, '*');
    // Handle powers like x^y and x^{...} carefully
    // This needs refinement to handle nested braces correctly, but basic cases:
    expr = expr.replace(/([a-zA-Z0-9_]+)\s*\^\s*([a-zA-Z0-9_]+)/g, '$1**($2)'); // Simple power x^y
    expr = expr.replace(/([a-zA-Z0-9_]+)\s*\^\s*\{([^}]+?)\}/g, '$1**($2)'); // Power with braces x^{y}
    expr = expr.replace(/\{\s*([^}]+?)\s*\}/g, '($1)'); // Replace remaining braces with parentheses
    expr = expr.replace(/\^/g, '**'); // Replace any remaining ^

    // Final cleanup for potential double parentheses or misplaced operators might be needed
    expr = expr.replace(/\(\(([^)]+)\)\)/g, '($1)'); // Simplify double parentheses

    return expr;
}

function python_to_numpy(expr: string): string {
    let numpyExpr = expr;
    // Handle exp first
    numpyExpr = numpyExpr.replace(/\be\s*\*\*\s*\(([^)]+)\)/g, 'np.exp($1)');
    // Handle power recursively/iteratively - simple approach for now
    // A more robust parser would be better here.
    while (numpyExpr.includes('**')) {
         // Prioritize powers within parentheses if possible, or simple variable/number powers
         numpyExpr = numpyExpr.replace(/([a-zA-Z0-9_.]+)\s*\*\*\s*([a-zA-Z0-9_.]+|\([^)]+\))/g, 'np.power($1, $2)');
         // Fallback for more complex bases (might need refinement)
         numpyExpr = numpyExpr.replace(/\(([^)]+)\)\s*\*\*\s*([a-zA-Z0-9_.]+|\([^)]+\))/g, 'np.power($1, $2)');
         // Break if no more replacements are made to avoid infinite loops on complex cases
         if (!numpyExpr.includes('**')) break;
         // Add a safeguard if replacement isn't working
         if (numpyExpr.split('**').length === expr.split('**').length && numpyExpr.includes('**')) {
             console.warn("Potential infinite loop in numpy power conversion for:", expr);
             break; // Prevent infinite loop
         }
         expr = numpyExpr; // Update expr for loop check
    }

    // Handle log/ln
    numpyExpr = numpyExpr.replace(/(?<!np\.)\bln\s*\(([^)]+)\)/g, 'np.log($1)');
    numpyExpr = numpyExpr.replace(/(?<!np\.)\blog\s*\(([^)]+)\)/g, 'np.log10($1)'); // Assuming log is log10

    // Handle trig functions
    numpyExpr = numpyExpr.replace(/\b(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)\s*\(([^)]+)\)/g, 'np.$1($2)');

    return numpyExpr;
}

// --- React Component ---

// Simple Copy Button Component
const CopyButton = ({ textToCopy }: { textToCopy: string }) => {
    const [copied, setCopied] = useState(false);

    const handleCopy = async () => {
        try {
            await navigator.clipboard.writeText(textToCopy);
            setCopied(true);
            setTimeout(() => setCopied(false), 1500); // Reset after 1.5s
        } catch (err) {
            console.error('Failed to copy text: ', err);
            // Optionally show an error state
        }
    };

    return (
        <Button variant="ghost" size="icon" onClick={handleCopy} className="h-6 w-6 ml-2">
            {copied ? <Check className="h-4 w-4 text-green-500" /> : <Copy className="h-4 w-4" />}
        </Button>
    );
};


export default function LatexConverterPage() {
    const mathQuillEl = useRef<HTMLDivElement>(null);
    const mathField = useRef<any>(null); // To store the MathQuill instance

    const [rawLatex, setRawLatex] = useState('');
    const [pythonExpr, setPythonExpr] = useState('');
    const [numpyExpr, setNumpyExpr] = useState('');

    // Initialize MathQuill
    useEffect(() => {
        // Ensure this runs only on the client
        const initializeMathQuill = () => {
            if (typeof window !== 'undefined' && (window as any).jQuery && (window as any).MathQuill && mathQuillEl.current) {
                console.log("jQuery and MathQuill found, initializing...");
                const MQ = (window as any).MathQuill.getInterface(2); // Use API version 2
                // Prevent re-initialization
                if (!mathField.current && mathQuillEl.current) {
                     // Check if the element already has MathQuill data to avoid re-init
                     if (!(window as any).jQuery(mathQuillEl.current).data('[[mathquill internal data]]')) {
                         mathField.current = MQ.MathField(mathQuillEl.current, {
                             spaceBehavesLikeTab: true,
                             handlers: {
                                 edit: () => {
                                     // Optional: update state on edit
                                 }
                             }
                         });
                         console.log("MathQuill Initialized successfully.");
                     } else {
                         console.log("MathQuill already initialized on this element.");
                         // Attempt to retrieve existing instance if needed, though direct re-assignment might be complex
                     }
                } else {
                    console.log("MathField ref already exists or element ref is missing.");
                }
            } else {
                console.log("MathQuill or jQuery not ready or element not found. Retrying...");
                // Retry initialization after a short delay
                setTimeout(initializeMathQuill, 100); // Retry after 100ms
            }
        };

        initializeMathQuill(); // Start the initialization attempt

        // Cleanup function (optional, as MathField lacks a destroy method)
        // return () => { ... };
    }, []); // Empty dependency array ensures this runs once on mount

    const handleGetExpression = useCallback(() => {
        if (mathField.current) {
            const latex = mathField.current.latex();
            setRawLatex(latex);
            try {
                const pyExpr = latex2python(latex);
                setPythonExpr(pyExpr);
                const npExpr = python_to_numpy(pyExpr);
                setNumpyExpr(npExpr);
            } catch (error) {
                console.error("Conversion Error:", error);
                setPythonExpr("Error during conversion.");
                setNumpyExpr("Error during conversion.");
            }
        } else {
            console.error("MathField not initialized.");
            setRawLatex("Error: MathQuill not ready.");
            setPythonExpr("");
            setNumpyExpr("");
        }
    }, []); // Dependencies remain empty

    // Add KeyDown handler for the input area
    const handleInputKeyDown = (event: React.KeyboardEvent<HTMLDivElement>) => {
        if (event.key === 'Enter') {
            event.preventDefault(); // Prevent default Enter behavior (like adding a newline)
            handleGetExpression(); // Trigger conversion
        }
    };

    return (
        <div className="container mx-auto p-4 md:p-8">
            <Card className="mb-4">
                {/* Reduce top padding */}
                <CardContent className="pt-2"> {/* Changed pt-6 to pt-4 */}
                    <Label htmlFor="mathquill-input-area">Enter LaTeX Equation:</Label>
                    {/* MathQuill Input Area */}
                    <div
                        id="mathquill-input-area"
                        ref={mathQuillEl}
                        className="mathquill-input-field mt-4 p-3 border rounded-md min-h-[80px] bg-background text-foreground focus-within:ring-2 focus-within:ring-ring w-full"
                        style={{ fontFamily: 'mathquill-font, serif' }}
                        onKeyDown={handleInputKeyDown} // Add keydown listener
                        tabIndex={0} // Make div focusable
                    >
                        {/* Placeholder or initial content */}
                    </div>
                    {/* Button below the input */}
                    <Button onClick={handleGetExpression} className="mt-6 w-full md:w-auto">
                        Convert Expression
                    </Button>
                </CardContent>
            </Card>

            {rawLatex && (
                <div className="space-y-4">
                    <Card className="gap-2">
                        <CardHeader className="flex flex-row items-center justify-between pb-0">
                            <CardTitle className="text-sm font-medium">Raw LaTeX:</CardTitle>
                            <CopyButton textToCopy={rawLatex} />
                        </CardHeader>
                        <CardContent className="pt-0">
                            <pre className="text-xs p-3 bg-muted rounded-md overflow-x-auto">
                                <code>{rawLatex}</code>
                            </pre>
                        </CardContent>
                    </Card>

                    <Card className="gap-2">
                        <CardHeader className="flex flex-row items-center justify-between pb-0">
                            <CardTitle className="text-sm font-medium">Python Expression:</CardTitle>
                            <CopyButton textToCopy={pythonExpr} />
                        </CardHeader>
                        <CardContent className="pt-0">
                            <pre className="text-xs p-3 bg-muted rounded-md overflow-x-auto">
                                <code>{pythonExpr}</code>
                            </pre>
                        </CardContent>
                    </Card>

                    <Card className="gap-2">
                        <CardHeader className="flex flex-row items-center justify-between pb-0">
                            <CardTitle className="text-sm font-medium">NumPy Expression:</CardTitle>
                            <CopyButton textToCopy={numpyExpr} />
                        </CardHeader>
                        <CardContent className="pt-0">
                            <pre className="text-xs p-3 bg-muted rounded-md overflow-x-auto">
                                <code>{numpyExpr}</code>
                            </pre>
                        </CardContent>
                    </Card>
                </div>
            )}

            {/* Add specific styles for MathQuill if needed */}
            <style jsx global>{`
                .mathquill-input-field {
                    /* Add any specific overrides for the input field appearance */
                    font-size: 1.1rem; /* Adjust font size as needed */
                }
                .mq-editable-field .mq-cursor {
                    border-left: 1px solid currentColor; /* Use theme color for cursor */
                }
                /* Add more overrides as needed */
            `}</style>
        </div>
    );
}
