'use client';

import React, { useState, useEffect, useRef, useCallback } from 'react';
import { Card, CardContent, CardHeader, CardTitle } from "@/components/ui/card";
import { Button } from "@/components/ui/button";
import { Label } from "@/components/ui/label";
import { Copy, Check } from 'lucide-react'; // Icons for copy button

// --- Conversion Logic ---

function replaceFracManually(expr: string): string {
    const fracIndex = expr.indexOf('\\frac');
    if (fracIndex === -1) return expr;
    let braceLevel=0,numStart=-1,numEnd=-1,denStart=-1,denEnd=-1;
    for(let i=fracIndex+5;i<expr.length;i++){if(expr[i]==='{'){if(braceLevel===0)numStart=i;braceLevel++;}else if(expr[i]==='}'){braceLevel--;if(braceLevel===0&&numStart!==-1){numEnd=i;break;}}if(braceLevel===0&&numStart===-1&&i===fracIndex+5){if(fracIndex+6<expr.length){const n=expr[fracIndex+5],d=expr[fracIndex+6];if(!['{','}','\\',' '].includes(n)&&!['{','}','\\',' '].includes(d)){const r=`(${n})/(${d})`;return replaceFracManually(expr.substring(0,fracIndex)+r+expr.substring(fracIndex+7));}}}if(braceLevel<0)return expr;}
    if(numStart===-1||numEnd===-1)return expr;braceLevel=0;for(let i=numEnd+1;i<expr.length;i++){if(expr[i]==='{'){if(braceLevel===0)denStart=i;braceLevel++;}else if(expr[i]==='}'){braceLevel--;if(braceLevel===0&&denStart!==-1){denEnd=i;break;}}if(braceLevel<0)return expr;}
    if(denStart===-1||denEnd===-1)return expr;const num=replaceFracManually(expr.substring(numStart+1,numEnd)),den=replaceFracManually(expr.substring(denStart+1,denEnd)),rep=`(${num})/(${den})`,newExpr=expr.substring(0,fracIndex)+rep+expr.substring(denEnd+1);return replaceFracManually(newExpr);
}

function latex2python(latexStr: string): string {
    let expr = latexStr.replace(/^\$|\$$/g, '');
    expr = expr.replace(/\\left|\\right/g, '');
    expr = replaceFracManually(expr);
    expr = expr.replace(/\\operatorname\{([^{}]+)\}/g, '$1');
    expr = expr.replace(/\\?exp\(/g, 'e**(');
    expr = expr.replace(/\\(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh|log|ln)\b/g, '$1');
    expr = expr.replace(/\\cdot\s*/g, '*');
    expr = expr.replace(/\\pi\b/g, 'pi');
    expr = expr.replace(/\\sqrt\{([^}]+?)\}/g, 'sqrt($1)');
    expr = expr.replace(/([a-zA-Z0-9_]+|\([^)]+\))\s*\^\s*\{([^}]+?)\}/g, '$1**($2)');
    expr = expr.replace(/([a-zA-Z0-9_]+|\([^)]+\))\s*\^\s*([a-zA-Z0-9_.]+)/g, '$1**$2');
    expr = expr.replace(/\{\s*([^}]+?)\s*\}/g, '($1)');
    expr = expr.replace(/\^/g, '**');
    while (/\(\(([^()]*)\)\)/.test(expr)) {
      expr = expr.replace(/\(\(([^()]*)\)\)/g, '($1)');
    }
    expr = expr.replace(/\s+/g, ' ');
    return expr.trim();
}

const convertPower = (expr: string, powFunc: string): string => {
    let resultExpr = expr;
    while (resultExpr.includes('**')) {
        const baseRegex = /(?:[a-df-zA-Z0-9_][a-zA-Z0-9_.]*|\d+|\((?:[^()]|\((?:[^()]|\([^()]*\))*\))*\))/;
        const expRegex = /(?:[a-zA-Z0-9_.]+|\((?:[^()]|\((?:[^()]|\([^()]*\))*\))*\))/;
        const powerRegex = new RegExp(`(${baseRegex.source})\\s*\\*\\*\\s*(${expRegex.source})`);
        const match = resultExpr.match(powerRegex);
        if (!match) break;
        resultExpr = resultExpr.replace(powerRegex, `${powFunc}(${match[1]}, ${match[2]})`);
    }
    return resultExpr;
}

const convertExp = (expr: string, expFunc: string): string => {
    return expr.replace(/\be\s*\*\*\s*(\((?:[^()]|\((?:[^()]|\([^()]*\))*\))*\)|[a-zA-Z0-9_.]+)/g, `${expFunc}($1)`);
}

function python_to_numpy(expr: string): string {
    let numpyExpr = convertExp(expr, 'np.exp');
    numpyExpr = convertPower(numpyExpr, 'np.power');
    numpyExpr = numpyExpr.replace(/(?<!np\.)\bsqrt\s*\(([^)]+)\)/g, 'np.sqrt($1)');
    numpyExpr = numpyExpr.replace(/(?<!np\.)\bln\s*\(([^)]+)\)/g, 'np.log($1)');
    numpyExpr = numpyExpr.replace(/(?<!np\.)\blog\s*\(([^)]+)\)/g, 'np.log10($1)');
    numpyExpr = numpyExpr.replace(/(?<!np\.)\b(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)\s*\(([^)]+)\)/g, 'np.$1($2)');
    numpyExpr = numpyExpr.replace(/(?<![a-zA-Z0-9_.])\bpi\b/g, 'np.pi');
    return numpyExpr;
}

function python_to_javascript(expr: string): string {
    let jsExpr = convertExp(expr, 'Math.exp');
    jsExpr = convertPower(jsExpr, 'Math.pow');
    jsExpr = jsExpr.replace(/(?<!Math\.)\bsqrt\s*\(([^)]+)\)/g, 'Math.sqrt($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\bln\s*\(([^)]+)\)/g, 'Math.log($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\blog\s*\(([^)]+)\)/g, 'Math.log10($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\b(sinh|cosh|tanh|sin|cos|tan)\s*\(([^)]+)\)/g, 'Math.$1($2)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\barcsinh\s*\(([^)]+)\)/g, 'Math.asinh($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\barccosh\s*\(([^)]+)\)/g, 'Math.acosh($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\barctanh\s*\(([^)]+)\)/g, 'Math.atanh($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\barcsin\s*\(([^)]+)\)/g, 'Math.asin($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\barccos\s*\(([^)]+)\)/g, 'Math.acos($1)');
    jsExpr = jsExpr.replace(/(?<!Math\.)\barctan\s*\(([^)]+)\)/g, 'Math.atan($1)');
    jsExpr = jsExpr.replace(/(?<![a-zA-Z0-9_.])\bpi\b/g, 'Math.PI');
    jsExpr = jsExpr.replace(/(?<![a-zA-Z0-9_.])\be\b/g, 'Math.E');
    return jsExpr;
}

function python_to_r(expr: string): string {
    let rExpr = convertExp(expr, 'exp');
    rExpr = rExpr.replace(/\*\*/g, '^');
    rExpr = rExpr.replace(/\bsqrt\s*\(([^)]+)\)/g, 'sqrt($1)');
    // FIX: Swapped order to prevent 'ln' from becoming 'log10'
    rExpr = rExpr.replace(/\blog\s*\(([^)]+)\)/g, 'log10($1)');
    rExpr = rExpr.replace(/\bln\s*\(([^)]+)\)/g, 'log($1)');
    rExpr = rExpr.replace(/\b(sin|cos|tan|sinh|cosh|tanh)\s*\(([^)]+)\)/g, '$1($2)');
    rExpr = rExpr.replace(/\barcsin/g, 'asin').replace(/\barccos/g, 'acos').replace(/\barctan/g, 'atan');
    rExpr = rExpr.replace(/(?<![a-zA-Z0-9_.])\bpi\b/g, 'pi');
    rExpr = rExpr.replace(/(?<![a-zA-Z0-9_.])\be\b/g, 'exp(1)');
    return rExpr;
}

function python_to_matlab(expr: string): string {
    let matlabExpr = convertExp(expr, 'exp');
    matlabExpr = matlabExpr.replace(/\*\*/g, '^');
    matlabExpr = matlabExpr.replace(/\bsqrt\s*\(([^)]+)\)/g, 'sqrt($1)');
    // FIX: Swapped order to prevent 'ln' from becoming 'log10'
    matlabExpr = matlabExpr.replace(/\blog\s*\(([^)]+)\)/g, 'log10($1)');
    matlabExpr = matlabExpr.replace(/\bln\s*\(([^)]+)\)/g, 'log($1)');
    matlabExpr = matlabExpr.replace(/\barcsinh/g, 'asinh').replace(/\barccosh/g, 'acosh').replace(/\barctanh/g, 'atanh');
    matlabExpr = matlabExpr.replace(/\barcsin/g, 'asin').replace(/\barccos/g, 'acos').replace(/\barctan/g, 'atan');
    matlabExpr = matlabExpr.replace(/\b(sin|cos|tan|sinh|cosh|tanh)\s*\(([^)]+)\)/g, '$1($2)');
    matlabExpr = matlabExpr.replace(/(?<![a-zA-Z0-9_.])\bpi\b/g, 'pi');
    matlabExpr = matlabExpr.replace(/(?<![a-zA-Z0-9_.])\be\b/g, 'exp(1)');
    return matlabExpr;
}

function python_to_excel(expr: string): string {
    let excelExpr = convertExp(expr, 'EXP');
    excelExpr = excelExpr.replace(/\*\*/g, '^');
    const funcs = ['sqrt', 'ln', 'log', 'sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh', 'asin', 'acos', 'atan', 'asinh', 'acosh', 'atanh', 'pi'];
    funcs.forEach(f => {
        excelExpr = excelExpr.replace(new RegExp(`\\b${f}\\b`, 'g'), f.toUpperCase());
    });
    excelExpr = excelExpr.replace(/\bPI\b/g, 'PI()');
    excelExpr = excelExpr.replace(/(?<![A-Z])\bE\b/g, 'EXP(1)');
    excelExpr = excelExpr.replace(/\bx\b/g, 'A1');
    return `= ${excelExpr}`;
}

function python_to_wolfram(expr: string): string {
    let wolframExpr = expr.replace(/\*\*/g, '^');
    wolframExpr = wolframExpr.replace(/\be\s*\^\s*/g, 'Exp');
    const funcs = {'sqrt':'Sqrt','ln':'Log','log':'Log10','sin':'Sin','cos':'Cos','tan':'Tan','sinh':'Sinh','cosh':'Cosh','tanh':'Tanh','arcsin':'ArcSin','arccos':'ArcCos','arctan':'ArcTan','arcsinh':'ArcSinh','arccosh':'ArcCosh','arctanh':'ArcTanh','pi':'Pi'};
    for (const [key, val] of Object.entries(funcs)) {
        wolframExpr = wolframExpr.replace(new RegExp(`\\b${key}\\b`, 'g'), val);
    }
    let previousExpr = "";
    while(previousExpr !== wolframExpr) {
        previousExpr = wolframExpr;
        wolframExpr = wolframExpr.replace(/([A-Z][a-zA-Z0-9]*)\(([^()]*)\)/g, '$1[$2]');
    }
    wolframExpr = wolframExpr.replace(/(?<![a-zA-Z0-9_.])\be\b/g, 'E');
    return wolframExpr;
}

function python_to_csharp(expr: string): string {
    let csExpr = convertExp(expr, 'Math.Exp');
    csExpr = convertPower(csExpr, 'Math.Pow');
    csExpr = csExpr.replace(/(?<!Math\.)\bsqrt\s*\(([^)]+)\)/g, 'Math.Sqrt($1)');
    csExpr = csExpr.replace(/(?<!Math\.)\bln\s*\(([^)]+)\)/g, 'Math.Log($1)');
    csExpr = csExpr.replace(/(?<!Math\.)\blog\s*\(([^)]+)\)/g, 'Math.Log10($1)');
    const arcFuncs = ['Asinh', 'Acosh', 'Atanh', 'Asin', 'Acos', 'Atan'];
    arcFuncs.forEach(f => {
        csExpr = csExpr.replace(new RegExp(`(?<!Math\\.)\\barc${f.substring(1).toLowerCase()}\\b`, 'g'), `Math.${f}`);
    });
    const funcs = ['Sin', 'Cos', 'Tan', 'Sinh', 'Cosh', 'Tanh'];
    funcs.forEach(f => {
        csExpr = csExpr.replace(new RegExp(`(?<!Math\\.)\\b${f.toLowerCase()}\\b`, 'g'), `Math.${f}`);
    });
    csExpr = csExpr.replace(/(?<![a-zA-Z0-9_.])\bpi\b/g, 'Math.PI');
    csExpr = csExpr.replace(/(?<![a-zA-Z0-9_.])\be\b/g, 'Math.E');
    return csExpr;
}

function python_to_cpp(expr: string): string {
    let cppExpr = convertExp(expr, 'std::exp');
    cppExpr = convertPower(cppExpr, 'std::pow');
    cppExpr = cppExpr.replace(/(?<!std::)\bsqrt\s*\(([^)]+)\)/g, 'std::sqrt($1)');
    cppExpr = cppExpr.replace(/(?<!std::)\bln\s*\(([^)]+)\)/g, 'std::log($1)');
    cppExpr = cppExpr.replace(/(?<!std::)\blog\s*\(([^)]+)\)/g, 'std::log10($1)');
    const arcFuncs = ['asinh', 'acosh', 'atanh', 'asin', 'acos', 'atan'];
    arcFuncs.forEach(f => {
        cppExpr = cppExpr.replace(new RegExp(`(?<!std::)\\barc${f.substring(1)}\\b`, 'g'), `std::${f}`);
    });
    const funcs = ['sin', 'cos', 'tan', 'sinh', 'cosh', 'tanh'];
    funcs.forEach(f => {
        cppExpr = cppExpr.replace(new RegExp(`(?<!std::)\\b${f}\\b`, 'g'), `std::${f}`);
    });
    cppExpr = cppExpr.replace(/(?<![a-zA-Z0-9_:.])\bpi\b/g, 'std::acos(-1.0)');
    cppExpr = cppExpr.replace(/(?<![a-zA-Z0-9_:.])\be\b/g, 'std::exp(1.0)');
    return cppExpr;
}

// --- React Components ---
const CopyButton = ({ textToCopy }: { textToCopy: string }) => {
    const [copied, setCopied] = useState(false);
    const handleCopy = async () => {
        if (!textToCopy) return;
        try { await navigator.clipboard.writeText(textToCopy); setCopied(true); setTimeout(() => setCopied(false), 1500); } catch (err) { console.error('Failed to copy text: ', err); alert('Failed to copy text.'); }
    };
    return (<Button variant="ghost" size="icon" onClick={handleCopy} className="h-6 w-6 ml-2" aria-label="Copy">{copied ? <Check className="h-4 w-4 text-green-500" /> : <Copy className="h-4 w-4" />}</Button>);
};
const ResultCard = ({ title, expression }: { title: string; expression: string }) => {
    if (!expression) return null;
    return (<Card className="gap-2"><CardHeader className="flex flex-row items-center justify-between pb-0"><CardTitle className="text-sm font-medium">{title}</CardTitle><CopyButton textToCopy={expression} /></CardHeader><CardContent className="pt-0"><pre className="text-xs p-3 bg-muted rounded-md overflow-x-auto"><code>{expression}</code></pre></CardContent></Card>);
};

export default function LatexConverterPage() {
    const mathQuillEl = useRef<HTMLDivElement>(null);
    const mathField = useRef<any>(null);
    const [rawLatex, setRawLatex] = useState('');
    const [pythonExpr, setPythonExpr] = useState('');
    const [numpyExpr, setNumpyExpr] = useState('');
    const [jsExpr, setJsExpr] = useState('');
    const [rExpr, setRExpr] = useState('');
    const [matlabExpr, setMatlabExpr] = useState('');
    const [excelExpr, setExcelExpr] = useState('');
    const [wolframExpr, setWolframExpr] = useState('');
    const [csharpExpr, setCsharpExpr] = useState('');
    const [javaExpr, setJavaExpr] = useState('');
    const [cppExpr, setCppExpr] = useState('');
    const [errorMsg, setErrorMsg] = useState('');
    
    useEffect(() => {
        let timeoutId: NodeJS.Timeout | null = null;
        const tryInitialize = (attempts = 0) => {
            const maxAttempts = 30;
            if (typeof window === 'undefined') return;
            const jQueryReady = typeof (window as any).jQuery === 'function';
            const MQReady = jQueryReady && typeof (window as any).MathQuill?.getInterface === 'function';
            if (MQReady && mathQuillEl.current && !mathField.current) {
                try {
                    const MQ = (window as any).MathQuill.getInterface(2);
                    mathField.current = MQ.MathField(mathQuillEl.current, { spaceBehavesLikeTab: true });
                    const defaultLatex = '42\\cdot \\frac{e^{\\frac{2}{\\ln\\left(2\\right)}}}{0.5}\\cdot \\sinh\\left(x\\right)';
                    mathField.current.latex(defaultLatex);
                    setTimeout(() => handleGetExpression(), 100);
                } catch (initError) { console.error("MathQuill init failed:", initError); setErrorMsg("Failed to init math editor."); }
            } else if (attempts < maxAttempts) {
                timeoutId = setTimeout(() => tryInitialize(attempts + 1), 200);
            } else { setErrorMsg("Math editor failed to load. Please reload."); }
        };
        tryInitialize();
        return () => { if (timeoutId) clearTimeout(timeoutId); };
    }, []);

    const handleGetExpression = useCallback(() => {
        setErrorMsg('');
        if (!mathField.current) { setErrorMsg("Editor not ready."); return; }
        const latex = mathField.current.latex();
        const statesToClear: React.Dispatch<React.SetStateAction<string>>[]=[setRawLatex,setPythonExpr,setNumpyExpr,setJsExpr,setRExpr,setMatlabExpr,setExcelExpr,setWolframExpr,setCsharpExpr,setJavaExpr,setCppExpr];
        if (!latex.trim()) { statesToClear.forEach(s => s('')); return; }
        setRawLatex(latex);
        try {
            const pyExpr = latex2python(latex);
            setPythonExpr(pyExpr);
            setNumpyExpr(python_to_numpy(pyExpr));
            setJsExpr(python_to_javascript(pyExpr));
            setRExpr(python_to_r(pyExpr));
            setMatlabExpr(python_to_matlab(pyExpr));
            setExcelExpr(python_to_excel(pyExpr));
            setWolframExpr(python_to_wolfram(pyExpr));
            setCsharpExpr(python_to_csharp(pyExpr));
            setJavaExpr(python_to_javascript(pyExpr)); 
            setCppExpr(python_to_cpp(pyExpr));
        } catch (error: any) {
            console.error("Conversion Error:", error);
            setErrorMsg(`Conversion Error: ${error.message || 'Unknown'}`);
            statesToClear.slice(1).forEach(s => s('Error'));
            setRawLatex(latex); 
        }
    }, []);

    const handleInputKeyDown = (event: React.KeyboardEvent<HTMLDivElement>) => {
        if (event.key === 'Enter' && !event.shiftKey) { event.preventDefault(); handleGetExpression(); }
    };

    return (
        <div className="container mx-auto p-4 md:p-8">
            <Card className="mb-4">
                <CardContent className="pt-6">
                    <Label htmlFor="mathquill-input-area">Enter LaTeX Equation:</Label>
                    <div id="mathquill-input-area" ref={mathQuillEl} className="mathquill-input-field mt-4 p-3 border rounded-md min-h-[80px] bg-background text-foreground focus-within:ring-2 focus-within:ring-ring w-full" onKeyDown={handleInputKeyDown} tabIndex={0} aria-label="LaTeX Input Area"></div>
                    {errorMsg && <p className="mt-2 text-sm text-red-600">{errorMsg}</p>}
                    <Button onClick={handleGetExpression} className="mt-6 w-full md:w-auto">Convert Expression</Button>
                </CardContent>
            </Card>
            {rawLatex && !errorMsg && (
                <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                    <ResultCard title="Raw LaTeX:" expression={rawLatex} />
                    <ResultCard title="TypeScript/JavaScript Expression:" expression={jsExpr} />
                    <ResultCard title="Python Expression:" expression={pythonExpr} />
                    <ResultCard title="NumPy Expression:" expression={numpyExpr} />
                    <ResultCard title="Wolfram Language (Mathematica):" expression={wolframExpr} />
                    <ResultCard title="Excel / Google Sheets:" expression={excelExpr} />
                    <ResultCard title="C++ Expression:" expression={cppExpr} />
                    <ResultCard title="C# Expression:" expression={csharpExpr} />
                    <ResultCard title="R Expression:" expression={rExpr} />
                    <ResultCard title="MATLAB / Octave Expression:" expression={matlabExpr} />
                    <ResultCard title="Java Expression:" expression={javaExpr} />
                </div>
            )}
            <style jsx global>{`.mathquill-input-field{font-size:1.1rem;}.mq-editable-field .mq-cursor{border-left:1px solid currentColor;}`}</style>
        </div>
    );
}