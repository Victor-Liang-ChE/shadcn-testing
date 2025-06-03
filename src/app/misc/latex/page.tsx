'use client';

import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Copy, Check } from 'lucide-react'; // Icons for copy button

// --- Conversion Logic (Ported from Python) ---

// NOTE: Conversion functions (replaceFracManually, fullyReplaceFrac,
// latex2python, python_to_numpy v4) are included here as they were
// in the previous "final" version.

function replaceFracManually(expr: string): string {
    const fracIndex = expr.indexOf('\\frac');
    if (fracIndex === -1) return expr;

    let braceLevel = 0;
    let numStart = -1, numEnd = -1, denStart = -1, denEnd = -1;

    // Find numerator braces {}
    for (let i = fracIndex + 5; i < expr.length; i++) {
        if (expr[i] === '{') {
            if (braceLevel === 0) numStart = i;
            braceLevel++;
        } else if (expr[i] === '}') {
            braceLevel--;
            if (braceLevel === 0 && numStart !== -1) {
                numEnd = i;
                break; // Found numerator
            }
        }
        // Handle simple frac like \frac12 early exit if applicable
        if (braceLevel === 0 && numStart === -1 && i === fracIndex + 5) {
             if (fracIndex + 6 < expr.length) {
               const numChar = expr[fracIndex + 5];
               const denChar = expr[fracIndex + 6];
               if (!['{', '}', '\\', ' '].includes(numChar) && !['{', '}', '\\', ' '].includes(denChar)) {
                   const replacement = `(${numChar})/(${denChar})`;
                   const before = expr.substring(0, fracIndex);
                   const after = expr.substring(fracIndex + 7);
                   return replaceFracManually(before + replacement + after);
               }
            }
        }
        if (braceLevel < 0) return expr; // Malformed
    }

    if (numStart === -1 || numEnd === -1) return expr; // Malformed or simple frac handled above

    // Find denominator braces {} starting after numerator
    braceLevel = 0;
    for (let i = numEnd + 1; i < expr.length; i++) {
        if (expr[i] === '{') {
            if (braceLevel === 0) denStart = i;
            braceLevel++;
        } else if (expr[i] === '}') {
            braceLevel--;
            if (braceLevel === 0 && denStart !== -1) {
                denEnd = i;
                break; // Found denominator
            }
        }
        if (braceLevel < 0) return expr; // Malformed
    }

    if (denStart === -1 || denEnd === -1) return expr; // Malformed

    // Extract content inside braces
    const numerator = expr.substring(numStart + 1, numEnd);
    const denominator = expr.substring(denStart + 1, denEnd);

    // Recursively process the extracted numerator and denominator first
    const newNum = replaceFracManually(numerator);
    const newDen = replaceFracManually(denominator);

    // Construct the replacement: (numerator)/(denominator)
    const replacement = `(${newNum})/(${newDen})`;

    // Assemble the new expression string
    const beforeFrac = expr.substring(0, fracIndex);
    const afterFrac = expr.substring(denEnd + 1);
    const newExpr = beforeFrac + replacement + afterFrac;

    // Recursively process the *entire new string* to find subsequent \frac instances
    return replaceFracManually(newExpr);
}

function fullyReplaceFrac(expr: string): string {
    return replaceFracManually(expr);
}

function latex2python(latexStr: string): string {
    let expr = latexStr.replace(/^\$|\$$/g, '');
    expr = expr.replace(/\\left|\\right/g, '');
    expr = fullyReplaceFrac(expr);
    expr = expr.replace(/\\operatorname\{([^{}]+)\}/g, '$1');
    expr = expr.replace(/\\?exp\(/g, 'e**(');
    expr = expr.replace(/\\(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh|log|ln)\b/g, '$1');
    expr = expr.replace(/\\cdot\s*/g, '*');
    expr = expr.replace(/\\pi\b/g, 'pi');
    expr = expr.replace(/\\sqrt\{([^}]+?)\}/g, 'sqrt($1)');

    expr = expr.replace(/([a-zA-Z0-9_]+|\([^)]+\))\s*\^\s*\{([^}]+?)\}/g, '$1**($2)');
    expr = expr.replace(/([a-zA-Z0-9_]+|\([^)]+\))\s*\^\s*([a-zA-Z0-9_]+)/g, '$1**($2)');
    expr = expr.replace(/\{\s*([^}]+?)\s*\}/g, '($1)');
    expr = expr.replace(/\^/g, '**');

    while (/\(\(([^()]*)\)\)/.test(expr)) {
      expr = expr.replace(/\(\(([^()]*)\)\)/g, '($1)');
    }
    expr = expr.replace(/\s+/g, ' ');
    return expr.trim();
}

function python_to_numpy(expr: string): string {
    let numpyExpr = expr;
    let initialExpr;
    let maxIterations = 100;

    // --- STEP 1: Convert e**(...) to np.exp(...) iteratively ---
    let safeguard = 0;
    const expRegex = /\be\s*\*\*\s*(\((?:[^()]|\([^()]*\))*\)|[a-zA-Z0-9_.]+)/g;
    while (/\be\s*\*\*/.test(numpyExpr) && safeguard < maxIterations) {
        initialExpr = numpyExpr;
        numpyExpr = numpyExpr.replace(expRegex, 'np.exp($1)');
        if (numpyExpr === initialExpr) {
            numpyExpr = numpyExpr.replace(/\be\s*\*\*\s*\(([^)]+)\)/g, 'np.exp($1)');
        }
        if (numpyExpr === initialExpr) {
             // console.warn("e** conversion stagnated for:", initialExpr);
             break;
        }
        safeguard++;
    }
     if (safeguard >= maxIterations) console.error("Max iterations reached during e** conversion for:", expr);

    // --- STEP 2: Convert remaining ** to np.power ---
    let iterations = 0;
    while (numpyExpr.includes('**') && iterations < maxIterations) {
        initialExpr = numpyExpr;
        const powerBaseRegex = /(?:[a-df-zA-Z0-9_][a-zA-Z0-9_.]*|\d+|\((?:[^()]|\([^()]*\))*\))/;
        const powerExpRegex = /(?:[a-zA-Z0-9_.]+|\((?:[^()]|\([^()]*\))*\))/;
        const powerRegexStr = `(${powerBaseRegex.source})\\s*\\*\\*\\s*(${powerExpRegex.source})`;
        const powerRegex = new RegExp(powerRegexStr, 'g');
        numpyExpr = numpyExpr.replace(powerRegex, 'np.power($1, $2)');
        if (numpyExpr === initialExpr) {
             numpyExpr = numpyExpr.replace(/\(((?:[^()]|\([^()]*\))*)\)\s*\*\*\s*([a-zA-Z0-9_.]+)/g, 'np.power($1, $2)');
             if (numpyExpr === initialExpr) {
                // console.warn("Power conversion loop stagnated for:", initialExpr);
                 break;
             }
        }
        iterations++;
    }
     if (iterations >= maxIterations) console.error("Max iterations reached during power conversion for:", expr);

    // --- STEP 3: Handle sqrt() ---
    numpyExpr = numpyExpr.replace(/(?<!np\.)\bsqrt\s*\(([^)]+)\)/g, 'np.sqrt($1)');
    // --- STEP 4: Handle log/ln ---
    numpyExpr = numpyExpr.replace(/(?<!np\.)\bln\s*\(([^)]+)\)/g, 'np.log($1)');
    numpyExpr = numpyExpr.replace(/(?<!np\.)\blog\s*\(([^)]+)\)/g, 'np.log10($1)');
    // --- STEP 5: Handle trig functions ---
    numpyExpr = numpyExpr.replace(/(?<!np\.)\b(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)\s*\(([^)]+)\)/g, 'np.$1($2)');
    // --- STEP 6: Handle 'pi' constant ---
    numpyExpr = numpyExpr.replace(/(?<![a-zA-Z0-9_.])\bpi\b(?![a-zA-Z0-9_.])/g, 'np.pi');
    // --- STEP 7: REMOVED Flawed Parenthesis Cleanup & standalone 'e' conversion ---

    // console.log("Final NumPy Expr (Attempt 4):", numpyExpr);
    return numpyExpr;
}


// --- React Component ---

const CopyButton = ({ textToCopy }: { textToCopy: string }) => {
    const [copied, setCopied] = useState(false);
    const handleCopy = async () => {
        if (!textToCopy) return;
        try {
            await navigator.clipboard.writeText(textToCopy);
            setCopied(true);
            setTimeout(() => setCopied(false), 1500);
        } catch (err) {
            console.error('Failed to copy text: ', err);
            alert('Failed to copy text to clipboard.');
        }
    };
    return (
        <Button variant="ghost" size="icon" onClick={handleCopy} className="h-6 w-6 ml-2" aria-label="Copy to clipboard">
            {copied ? <Check className="h-4 w-4 text-green-500" /> : <Copy className="h-4 w-4" />}
        </Button>
    );
};


export default function LatexConverterPage() {
    const mathQuillEl = useRef<HTMLDivElement>(null);
    const mathField = useRef<any>(null);
    const [rawLatex, setRawLatex] = useState('');
    const [pythonExpr, setPythonExpr] = useState('');
    const [numpyExpr, setNumpyExpr] = useState('');
    const [errorMsg, setErrorMsg] = useState('');

    // ============================================================
    // START OF ROBUST useEffect FOR MATHQUILL INITIALIZATION
    // ============================================================
    useEffect(() => {
        let attempts = 0;
        const maxAttempts = 30; // Try for ~6 seconds (30 * 200ms)
        let timeoutId: NodeJS.Timeout | null = null;

        const tryInitialize = () => {
            attempts++;
            // console.log(`MathQuill init attempt ${attempts}`); // Optional: for debugging

            // Check more thoroughly if libraries are ready and element exists
            const jQueryReady = typeof window !== 'undefined' && typeof (window as any).jQuery === 'function';
            // Check if MathQuill object exists AND its core function exists
            const MQReady = jQueryReady && typeof (window as any).MathQuill?.getInterface === 'function';
            const elementReady = mathQuillEl.current !== null;

            if (MQReady && elementReady) {
                // console.log("Dependencies and element ready."); // Optional: for debugging
                const MQ = (window as any).MathQuill.getInterface(2);

                // Prevent re-initialization checks (important!)
                if (mathField.current || (window as any).jQuery(mathQuillEl.current).data('[[mathquill internal data]]')) {
                   // console.log("Skipping: Already initialized."); // Optional: for debugging
                   return;
                }

                // Initialize
                try {
                    mathField.current = MQ.MathField(mathQuillEl.current, {
                         spaceBehavesLikeTab: true,
                         handlers: {
                             edit: () => {
                                 // Clear errors/results on edit if desired by user
                                 // setErrorMsg('');
                                 // setRawLatex(''); setPythonExpr(''); setNumpyExpr('');
                             }
                         }
                    });
                    // console.log("MathQuill Initialized successfully."); // Optional: for debugging
                } catch (initError) {
                    console.error("MathQuill initialization failed:", initError);
                    setErrorMsg("Failed to initialize math editor."); // Show error to user
                }

            } else {
                 // console.log(`Not ready: jQuery=${jQueryReady}, MQ=${MQReady}, Element=${elementReady}`); // Optional: for debugging
                if (attempts < maxAttempts) {
                    // Retry if not ready and under max attempts
                    timeoutId = setTimeout(tryInitialize, 200);
                } else {
                    // Stop retrying and inform user
                    console.error(`MathQuill initialization failed after ${maxAttempts} attempts.`);
                    setErrorMsg("Math editor failed to load. Please try reloading the page.");
                }
            }
        };

        // Start the first attempt slightly delayed
        timeoutId = setTimeout(tryInitialize, 150);

        // Cleanup function to clear timeout if component unmounts
        return () => {
            if (timeoutId) {
                clearTimeout(timeoutId);
            }
            // Note: MathQuill lacks a standard 'destroy' method.
        };
    }, []); // Empty dependency array ensures this runs once on mount
    // ============================================================
    // END OF ROBUST useEffect FOR MATHQUILL INITIALIZATION
    // ============================================================


    // Callback for triggering the conversion (Unchanged)
    const handleGetExpression = useCallback(() => {
        setErrorMsg(''); // Clear previous errors
        if (mathField.current) {
            const latex = mathField.current.latex();
            if (!latex.trim()) {
                 setRawLatex(''); setPythonExpr(''); setNumpyExpr(''); return;
            }
            setRawLatex(latex);
            try {
                const pyExpr = latex2python(latex);
                setPythonExpr(pyExpr);
                const npExpr = python_to_numpy(pyExpr);
                setNumpyExpr(npExpr);
            } catch (error: any) {
                console.error("Conversion Error:", error);
                setErrorMsg(`Conversion Error: ${error.message || 'Unknown error'}`);
                setPythonExpr("Error"); setNumpyExpr("Error");
            }
        } else {
            console.error("MathField not initialized.");
            setErrorMsg("Error: Math editor not ready.");
            setRawLatex(""); setPythonExpr(""); setNumpyExpr("");
        }
    }, []);


    // KeyDown handler (Unchanged)
    const handleInputKeyDown = (event: React.KeyboardEvent<HTMLDivElement>) => {
        if (event.key === 'Enter' && !event.shiftKey && !event.ctrlKey && !event.altKey) {
            event.preventDefault();
            handleGetExpression();
        }
    };

    // --- UI Section --- (Restored to original structure/styles)
    return (
        <div className="container mx-auto p-4 md:p-8">
            {/* Title Removed */}
            <Card className="mb-4">
                {/* Original Padding */}
                <CardContent className="pt-2">
                    <Label htmlFor="mathquill-input-area">Enter LaTeX Equation:</Label>
                    {/* MathQuill Input Area - Original Classes */}
                    <div
                        id="mathquill-input-area"
                        ref={mathQuillEl}
                        className="mathquill-input-field mt-4 p-3 border rounded-md min-h-[80px] bg-background text-foreground focus-within:ring-2 focus-within:ring-ring w-full"
                        style={{ fontFamily: 'mathquill-font, serif' }}
                        onKeyDown={handleInputKeyDown}
                        tabIndex={0}
                        aria-label="LaTeX Input Area"
                    >
                        {/* MathQuill fills this */}
                    </div>
                    {/* Error Message Area */}
                    {errorMsg && (
                       <p className="mt-2 text-sm text-red-600">{errorMsg}</p>
                    )}
                    {/* Button - Original Classes */}
                    <Button onClick={handleGetExpression} className="mt-6 w-full md:w-auto">
                        Convert Expression
                    </Button>
                </CardContent>
            </Card>

             {/* Results Section - Original Structure/Styles */}
            {rawLatex && !errorMsg && (
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

            {/* Original Style block */}
            <style jsx global>{`
                .mathquill-input-field {
                    font-size: 1.1rem;
                }
                .mq-editable-field .mq-cursor {
                    border-left: 1px solid currentColor;
                }
            `}</style>
        </div>
    );
}