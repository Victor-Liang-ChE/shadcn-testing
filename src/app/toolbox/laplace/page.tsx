'use client'

import React, { useState, useEffect, useRef } from 'react'

interface ComponentProps {
  className?: string;
  [key: string]: any;
}

interface IconProps {
  [key: string]: any;
}

interface MathQuillField {
  latex: (value?: string) => string;
}

interface MathQuillInterface {
  MathField: (element: HTMLElement, config?: any) => MathQuillField;
}

declare global {
  interface Window {
    MathQuill?: {
      getInterface: (version: number) => MathQuillInterface;
    };
  }
}

// --- shadcn/ui Component Definitions (self-contained) ---

const Card = ({ className = '', children, ...props }: ComponentProps) => (
  <div
    className={`rounded-xl border bg-card text-card-foreground shadow ${className}`}
    {...props}
  >
    {children}
  </div>
);

// Minimal local placeholders for shadcn/ui components used in this page
const CardHeader = ({ children, className = '', ...props }: ComponentProps) => (
  <div className={`px-4 pt-4 ${className}`} {...props}>{children}</div>
);
const CardTitle = ({ children, className = '', ...props }: ComponentProps) => (
  <h3 className={`text-lg font-semibold ${className}`} {...props}>{children}</h3>
);
const CardDescription = ({ children, className = '', ...props }: ComponentProps) => (
  <p className={`text-sm text-muted ${className}`} {...props}>{children}</p>
);
const CardContent = ({ children, className = '', ...props }: ComponentProps) => (
  <div className={`p-4 ${className}`} {...props}>{children}</div>
);
const Label = ({ children, className = '', ...props }: ComponentProps) => (
  <label className={`block text-sm font-medium ${className}`} {...props}>{children}</label>
);

const ArrowRightLeft = (props: IconProps) => (
    <svg
      {...props}
      xmlns="http://www.w3.org/2000/svg"
      width="24"
      height="24"
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth="2"
      strokeLinecap="round"
      strokeLinejoin="round"
    >
      <path d="M8 3L4 7l4 4" />
      <path d="M4 7h16" />
      <path d="M16 21l4-4-4-4" />
      <path d="M20 17H4" />
    </svg>
)


// --- AST Definition for Mathematical Expressions ---

// Define the building blocks of our math expressions
type ASTNode = BinaryOperation | Constant | Variable | FunctionCall | UnaryMinus;

interface BinaryOperation {
  type: 'BinaryOperation';
  operator: '+' | '-' | '*' | '/' | '^';
  left: ASTNode;
  right: ASTNode;
}

interface Constant {
  type: 'Constant';
  value: number;
}

interface Variable {
  type: 'Variable';
  name: string;
}

interface FunctionCall {
  type: 'FunctionCall';
  name: string;
  argument: ASTNode;
}

interface UnaryMinus {
  type: 'UnaryMinus';
  operand: ASTNode;
}

// --- LaTeX Preprocessing (simplified) ---
const preprocessLatex = (latex: string): string => {
  console.log("=== preprocessLatex START ===");
  console.log("Raw input:", JSON.stringify(latex));

  let str = latex;

  const knownTokens = ['sin', 'cos', 'delta', 'pi'];
  const placeholders = new Map<string, string>();

  // Step 1: Protect all known tokens (functions and constants) with numeric placeholders
  knownTokens.forEach((token, index) => {
    const placeholder = `@${index};`;
    placeholders.set(placeholder, token);
    
    // Functions can have an optional backslash, but pi must have one
    const pattern = token === 'pi' ? `\\\\${token}` : `\\\\?${token}`;
    const regex = new RegExp(pattern, 'g');
    str = str.replace(regex, placeholder);
  });
  console.log("After func placeholders:", str);

  // CORRECT ORDER: Handle fractions before exponent braces
  // Step 2: Normalize fractions
  str = str.replace(/\\frac\s*\{([^}]*)\}\s*\{([^}]*)\}/g, '($1)/($2)');
  
  // Step 3: Normalize exponent braces
  str = str.replace(/([a-zA-Z0-9])\^\{([^}]*)\}/g, '$1^($2)');

  // Step 4: General cleanup
  str = str.replace(/\\left|\\right/g, '');
  str = str.replace(/\\cdot|⋅/g, '*');
  str = str.replace(/\\/g, ''); // remove leftover backslashes
  str = str.replace(/\s+/g, '');
  console.log("After cleanup:", str);

  // Step 5: Braces to parentheses
  str = str.replace(/\{/g, '(').replace(/\}/g, ')');

  // Step 6: Insert explicit multiplication
  str = str.replace(/(\))([a-zA-Z0-9@\(])/g, '$1*$2');
  str = str.replace(/(\d)([a-zA-Z@\(])/g, '$1*$2');
  str = str.replace(/([a-zA-Z])([a-zA-Z0-9@\(])/g, '$1*$2');

  // Step 7: Restore placeholders
  placeholders.forEach((token, placeholder) => {
    str = str.replace(new RegExp(placeholder, 'g'), token);
  });

  console.log("=== preprocessLatex END === Final:", str);
  return str;
};

// --- Simple Parser for Mathematical Expressions ---
function parse(processedInput: string): ASTNode | null {
  processedInput = processedInput.trim();
  
  if (!processedInput) return null;

  // Remove outer parentheses if they wrap the entire expression
  if (processedInput.startsWith('(') && processedInput.endsWith(')')) {
    let parenCount = 0;
    let canRemove = true;
    for (let i = 0; i < processedInput.length; i++) {
      if (processedInput[i] === '(') parenCount++;
      if (processedInput[i] === ')') parenCount--;
      if (parenCount === 0 && i < processedInput.length - 1) {
        canRemove = false;
        break;
      }
    }
    if (canRemove) {
      processedInput = processedInput.slice(1, -1);
    }
  }

  // Handle unary minus at the beginning
  if (processedInput.startsWith('-')) {
    const operand = parse(processedInput.substring(1));
    if (operand) {
      return { type: 'UnaryMinus', operand };
    }
  }

  // Look for top-level + or - (lowest precedence)
  let parenDepth = 0;
  for (let i = processedInput.length - 1; i >= 0; i--) {
    const char = processedInput[i];
    if (char === ')') parenDepth++;
    if (char === '(') parenDepth--;
    if ((char === '+' || char === '-') && parenDepth === 0 && i > 0) {
      const left = parse(processedInput.substring(0, i));
      const right = parse(processedInput.substring(i + 1));
      if (left && right) {
        return {
          type: 'BinaryOperation',
          operator: char as '+' | '-',
          left,
          right,
        };
      }
    }
  }
  
  // Then look for top-level * or / (medium precedence)
  parenDepth = 0;
  for (let i = processedInput.length - 1; i >= 0; i--) {
    const char = processedInput[i];
    if (char === ')') parenDepth++;
    if (char === '(') parenDepth--;
    if ((char === '*' || char === '/') && parenDepth === 0) {
      const left = parse(processedInput.substring(0, i));
      const right = parse(processedInput.substring(i + 1));
      if (left && right) {
        return {
          type: 'BinaryOperation',
          operator: char as '*' | '/',
          left,
          right,
        };
      }
    }
  }

  // Then look for ^ (highest precedence)
  parenDepth = 0;
  for (let i = processedInput.length - 1; i >= 0; i--) {
    const char = processedInput[i];
    if (char === ')') parenDepth++;
    if (char === '(') parenDepth--;
    if (char === '^' && parenDepth === 0) {
      const left = parse(processedInput.substring(0, i));
      const right = parse(processedInput.substring(i + 1));
      if (left && right) {
        return {
          type: 'BinaryOperation',
          operator: '^',
          left,
          right,
        };
      }
    }
  }

  // Handle function calls like sin(t), e^(at)
  const funcMatch = processedInput.match(/^(\w+)\((.+)\)$/);
  if (funcMatch) {
    const argument = parse(funcMatch[2]);
    if (argument) {
      return {
        type: 'FunctionCall',
        name: funcMatch[1],
        argument
      };
    }
  }

  // Handle numbers (including decimals)
  const numMatch = processedInput.match(/^-?\d+(?:\.\d+)?$/);
  if (numMatch) {
    return { type: 'Constant', value: parseFloat(processedInput) };
  }
  
  // Handle variables
  if (/^[a-zA-Z]+$/.test(processedInput)) {
    return { type: 'Variable', name: processedInput };
  }
  
  // If all else fails, return null
  console.log(`Failed to parse: "${processedInput}"`);
  return null;
}

// --- AST to LaTeX Conversion ---
function astToLatex(ast: ASTNode | null): string {
  if (!ast) return '';

  if (ast.type === "BinaryOperation") {
    const left = astToLatex(ast.left);
    const right = astToLatex(ast.right);

    switch (ast.operator) {
      case "+":
      case "-":
        return `${left} ${ast.operator} ${right}`;

      case "*": {
        const leftStr = (ast.left.type === "BinaryOperation" &&
                         (ast.left.operator === "+" || ast.left.operator === "-"))
                        ? `(${left})` : left;
        const rightStr = (ast.right.type === "BinaryOperation" &&
                          (ast.right.operator === "+" || ast.right.operator === "-"))
                        ? `(${right})` : right;
        return `${leftStr}${rightStr}`; // implicit multiplication
      }

      case "/":
        return `\\frac{${left}}{${right}}`;

      case "^": {
        const baseStr = (ast.left.type === "BinaryOperation" &&
                         (ast.left.operator === "+" || ast.left.operator === "-"))
                        ? `(${left})` : left;
        return `${baseStr}^{${right}}`;
      }
    }
  }

  if (ast.type === "Variable") {
    if (ast.name === "pi") return "\\pi";
    return ast.name;
  }
  if (ast.type === "Constant") return ast.value.toString();
  if (ast.type === "FunctionCall")
    return `\\${ast.name}(${astToLatex(ast.argument)})`;
  if (ast.type === "UnaryMinus")
    return `-${astToLatex(ast.operand)}`;

  return "";
}

// Helper to check if two AST nodes are structurally identical
const areNodesEqual = (node1: ASTNode | null, node2: ASTNode | null): boolean => {
    if (!node1 || !node2) return node1 === node2;
    return JSON.stringify(node1) === JSON.stringify(node2);
};

// Helper to find all non-integer constants in an expression tree
const findNonIntegers = (node: ASTNode | null): number[] => {
    if (!node) return [];
    if (isConstant(node) && node.value % 1 !== 0) {
        return [node.value];
    }
    if (isBinaryOp(node)) {
        return [...findNonIntegers(node.left), ...findNonIntegers(node.right)];
    }
    if (isFunction(node)) return findNonIntegers(node.argument);
    if (isUnaryMinus(node)) return findNonIntegers(node.operand);
    return [];
}

// A smarter multiplication helper that knows how to distribute
const multiplyAndSimplify = (node: ASTNode, factor: ASTNode): ASTNode => {
    // If node is A+B or A-B, distribute: factor*(A±B) -> factor*A ± factor*B
    if (isBinaryOp(node) && (node.operator === '+' || node.operator === '-')) {
        const newLeft = multiplyAndSimplify(node.left, factor);
        const newRight = multiplyAndSimplify(node.right, factor);
        const result = simplifyAst(BinOp(node.operator, newLeft, newRight));
        return result || BinOp(node.operator, newLeft, newRight);
    }
    // If node is A*B, try to multiply factor into the term with decimals
    if (isBinaryOp(node) && node.operator === '*') {
        const leftHasDecimal = findNonIntegers(node.left).length > 0;
        if (leftHasDecimal) {
            const newLeft = multiplyAndSimplify(node.left, factor);
            const result = simplifyAst(BinOp('*', newLeft, node.right));
            return result || BinOp('*', newLeft, node.right);
        } else {
            const newRight = multiplyAndSimplify(node.right, factor);
            const result = simplifyAst(BinOp('*', node.left, newRight));
            return result || BinOp('*', node.left, newRight);
        }
    }
    // Default: just multiply the whole expression
    const result = simplifyAst(BinOp('*', factor, node));
    return result || BinOp('*', factor, node);
}

function simplifyAst(node: ASTNode | null): ASTNode | null {
  if (!node) return null;

  // 1. Simplify children first
  if (node.type === 'BinaryOperation') {
    node.left = simplifyAst(node.left) as ASTNode;
    node.right = simplifyAst(node.right) as ASTNode;
  } else if (node.type === 'FunctionCall') {
    node.argument = simplifyAst(node.argument) as ASTNode;
  } else if (node.type === 'UnaryMinus') {
    node.operand = simplifyAst(node.operand) as ASTNode;
    if (node.operand.type === 'UnaryMinus') {
      return node.operand.operand;
    }
  }

  // 2. Simplify the current node with fundamental rules
  if (isBinaryOp(node)) {
    const left = node.left;
    const right = node.right;

    // --- MODIFIED Constant Folding ---
    if (isConstant(left) && isConstant(right)) {
      switch (node.operator) {
        case '+': return Const(left.value + right.value);
        case '-': return Const(left.value - right.value);
        case '*': return Const(left.value * right.value);
        case '/':
          if (right.value !== 0) {
            const result = left.value / right.value;
            // ONLY perform division if the result is a whole number
            if (Number.isInteger(result)) {
              return Const(result);
            }
          }
          return node; // Otherwise, keep the fraction symbolic
      }
    }
    
    // (Other simplification rules remain the same)
    if (node.operator === '*') {
      if (isConstant(left, 1)) return right;
      if (isConstant(right, 1)) return left;
      if (isConstant(left, 0) || isConstant(right, 0)) return Const(0);
      // Cancel constants like (1/2) * (2 / expr) -> 1 / expr
      if (isBinaryOp(node.left, '/') && isConstant(node.left.left, 1) && isConstant(node.left.right, 2)) {
        if (isBinaryOp(node.right, '/') && isConstant(node.right.left, 2)) {
          return BinOp("/", Const(1), node.right.right);
        }
      }
      if (isBinaryOp(node.right, '/') && isConstant(node.right.left, 1) && isConstant(node.right.right, 2)) {
        if (isBinaryOp(node.left, '/') && isConstant(node.left.left, 2)) {
          return BinOp("/", Const(1), node.left.right);
        }
      }
      // Correctly flatten (1/3) * (1/(s+1/2)) → 2/(3(2s+1))
      if (
        node.operator === "*" &&
        isBinaryOp(node.left, "/") &&
        isConstant(node.left.left, 1) &&
        isConstant(node.left.right) &&
        isBinaryOp(node.right, "/") &&
        isConstant(node.right.left, 1) &&
        isBinaryOp(node.right.right, "+") &&
        isVariable(node.right.right.left, "s") &&
        isBinaryOp(node.right.right.right, "/") &&
        isConstant(node.right.right.right.left, 1) &&
        isConstant(node.right.right.right.right, 2)
      ) {
        // matches (1/3) * (1/(s+1/2))
        return BinOp(
          "/",
          Const(2),
          BinOp("*", node.left.right, BinOp("+", BinOp("*", Const(2), Var("s")), Const(1)))
        );
      }
      // Flatten (a/b) * (c/d) -> (a*c)/(b*d)
      if (isBinaryOp(node.left, '/') && isBinaryOp(node.right, '/')) {
        return simplifyAst(BinOp(
          '/',
          BinOp('*', node.left.left, node.right.left),   // numerator = a*c
          BinOp('*', node.left.right, node.right.right) // denominator = b*d
        ));
      }
      // Flatten (a/b) * c -> (a*c)/b
      if (isBinaryOp(node.left, '/') && isConstant(node.right)) {
        return simplifyAst(BinOp('/', BinOp('*', node.left.left, node.right), node.left.right));
      }
      // Flatten c * (a/b) -> (c*a)/b
      if (isConstant(node.left) && isBinaryOp(node.right, '/')) {
        return simplifyAst(BinOp('/', BinOp('*', node.left, node.right.left), node.right.right));
      }
    }
    if (node.operator === '+' || node.operator === '-') {
      if (isConstant(right, 0)) return left;
      if (node.operator === '+' && isConstant(left, 0)) return right;
      if (node.operator === '-' && areNodesEqual(left, right)) return Const(0);
      if (node.operator === '-' && isBinaryOp(left, '+')) {
        if (areNodesEqual(left.left, right)) return left.right;
        if (areNodesEqual(left.right, right)) return left.left;
      }
      if (isBinaryOp(left, '/') && isBinaryOp(right, '/')) {
        const A = left.left, B = left.right, C = right.left, D = right.right;
        const newNumerator = BinOp(node.operator, BinOp('*', A, D), BinOp('*', C, B));
        const newDenominator = BinOp('*', B, D);
        return simplifyAst(BinOp('/', newNumerator, newDenominator));
      }
    }
    // Flatten (A/s)/B -> A/(B*s)
    if (node.operator === '/' && isBinaryOp(node.left, '/') && isVariable(node.left.right, 's') && isConstant(node.right)) {
      return simplifyAst(BinOp("/", node.left.left, BinOp("*", node.right, Var("s"))));
    }
    // --- NEW: Simplify numeric coefficients in fractions (FIXED) ---
    // Rule: C1 / (C2 * expr) -> (C1/gcd) / ((C2/gcd) * expr)
    if (node.operator === '/') {
      let num: Constant | null = isConstant(node.left) ? node.left : null;
      let den: ASTNode | null = node.right;

      if (num && isBinaryOp(den, '*')) {
        let denConstNode: Constant | null = null;
        let denExprNode: ASTNode | null = null;

        if (isConstant(den.left)) {
          denConstNode = den.left;
          denExprNode = den.right;
        } else if (isConstant(den.right)) {
          // handles the case where the expression is (expr * C2)
          denConstNode = den.right;
          denExprNode = den.left;
        }

        if (denConstNode && denExprNode) {
          const numVal = num.value;
          const denConstVal = denConstNode.value;

          if (denConstVal !== 0 && Number.isInteger(numVal) && Number.isInteger(denConstVal)) {
            const commonDivisor = gcd(numVal, denConstVal);

            // GUARD CONDITION: If already simplified, stop to prevent infinite loop.
            if (commonDivisor <= 1) {
                return node;
            }

            const newNumVal = numVal / commonDivisor;
            const newDenConstVal = denConstVal / commonDivisor;

            const newNumNode = Const(newNumVal);

            if (newDenConstVal === 1) {
              // Denominator's constant factor is 1, so it disappears
              // e.g., 2 / (2 * (s+1)) becomes 1 / (s+1)
              return simplifyAst(BinOp('/', newNumNode, denExprNode));
            } else {
              // Denominator's constant is just reduced
              // e.g., 2 / (4 * (s+1)) becomes 1 / (2 * (s+1))
              const newDenNode = BinOp('*', Const(newDenConstVal), denExprNode);
              return simplifyAst(BinOp('/', newNumNode, newDenNode));
            }
          }
        }
      }
    }
  }

  // 3. Perform cosmetic cleanup on fractions containing 1/2
  if (isBinaryOp(node, '/')) {
    const containsHalf = (n: ASTNode | null): boolean => {
      if (!n) return false;
      if (isBinaryOp(n, '/') && isConstant(n.left, 1) && isConstant(n.right, 2)) return true;
      if (isBinaryOp(n)) return containsHalf(n.left) || containsHalf(n.right);
      if (n.type === 'FunctionCall') return containsHalf(n.argument);
      if (n.type === 'UnaryMinus') return containsHalf(n.operand);
      return false;
    };

    if (containsHalf(node)) {
      const multiplyByTwo = (n: ASTNode): ASTNode => {
        if (isBinaryOp(n, '/') && isConstant(n.left, 1) && isConstant(n.right, 2)) {
            return Const(1); // 2 * (1/2) = 1
        }
        if (n.type === 'BinaryOperation' && (n.operator === '+' || n.operator === '-')) {
          const simplifiedResult = simplifyAst(BinOp(n.operator, multiplyByTwo(n.left), multiplyByTwo(n.right)));
          return simplifiedResult || BinOp(n.operator, multiplyByTwo(n.left), multiplyByTwo(n.right));
        }
        if (n.type === 'BinaryOperation' && n.operator === '*') {
          if (containsHalf(n.right)) {
            const simplifiedResult = simplifyAst(BinOp('*', n.left, multiplyByTwo(n.right)));
            return simplifiedResult || BinOp('*', n.left, multiplyByTwo(n.right));
          }
          const simplifiedResult = simplifyAst(BinOp('*', multiplyByTwo(n.left), n.right));
          return simplifiedResult || BinOp('*', multiplyByTwo(n.left), n.right);
        }
        const simplifiedResult = simplifyAst(BinOp('*', Const(2), n));
        return simplifiedResult || BinOp('*', Const(2), n);
      };
      const newNum = multiplyByTwo(node.left);
      const newDenom = multiplyByTwo(node.right);
      return simplifyAst(BinOp('/', newNum, newDenom));
    }
  }

  return node;
}


// --- Utility Functions ---

const factorial = (n: string | number): number => {
  let num = typeof n === 'string' ? parseInt(n) : n;
  if (num < 0) return NaN;
  if (num === 0 || num === 1) return 1;
  for (let i = num - 1; i >= 1; i--) {
    num *= i;
  }
  return num;
};

const gcd = (a: number, b: number): number => {
  return b === 0 ? Math.abs(a) : gcd(b, a % b);
};

// --- Utility: substitute variable occurrences in an AST ---
function substitute(node: ASTNode | null, findName: string, replaceNode: ASTNode): ASTNode | null {
  if (!node) return null;

  // If the current node is the variable we're looking for, return a deep copy of the replacement
  if (node.type === 'Variable' && node.name === findName) {
    return JSON.parse(JSON.stringify(replaceNode));
  }

  // Recursively traverse and replace inside BinaryOperation
  if (node.type === 'BinaryOperation') {
    return {
      ...node,
      left: substitute(node.left, findName, replaceNode) as ASTNode,
      right: substitute(node.right, findName, replaceNode) as ASTNode,
    };
  }

  // Recursively traverse and replace inside FunctionCall
  if (node.type === 'FunctionCall') {
    return {
      ...node,
      argument: substitute(node.argument, findName, replaceNode) as ASTNode,
    };
  }

  // Recursively handle UnaryMinus
  if (node.type === 'UnaryMinus') {
    return {
      ...node,
      operand: substitute(node.operand, findName, replaceNode) as ASTNode,
    };
  }

  // Constants and other variables remain unchanged
  return node;
}

// --- Rule-Based Laplace Transform System ---

// Rule interface for pattern matching and transformation
interface LaplaceRule {
  name: string; // descriptive label
  match: (ast: ASTNode) => boolean; // pattern matcher
  apply: (ast: ASTNode) => ASTNode; // builds the transformed AST
}

// Utility AST Builders
const Var = (name: string): ASTNode => ({ type: 'Variable', name });
const Const = (value: number): ASTNode => ({ type: 'Constant', value });
const BinOp = (op: '+' | '-' | '*' | '/' | '^', left: ASTNode, right: ASTNode): ASTNode =>
  ({ type: 'BinaryOperation', operator: op, left, right });
const UnaryMinus = (node: ASTNode): ASTNode =>
  ({ type: 'UnaryMinus', operand: node });

// --- Pattern Matching for Laplace Transform Rules ---

// Helper functions for pattern matching
function isConstant(node: ASTNode | null, value?: number): node is Constant {
  return node?.type === 'Constant' && (value === undefined || node.value === value);
}

function isVariable(node: ASTNode | null, name?: string): node is Variable {
  return node?.type === 'Variable' && (name === undefined || node.name === name);
}

function isFunction(node: ASTNode | null, name?: string): node is FunctionCall {
  return node?.type === 'FunctionCall' && (name === undefined || node.name === name);
}

function isBinaryOp(node: ASTNode | null, operator?: string): node is BinaryOperation {
  return node?.type === 'BinaryOperation' && (operator === undefined || node.operator === operator);
}

function isUnaryMinus(node: ASTNode | null): node is UnaryMinus {
  return node?.type === 'UnaryMinus';
}

// Pattern matching for exponential expressions like e^(a*t)
function matchExponential(node: ASTNode): { a: ASTNode } | null {
  if (isBinaryOp(node, '^') && isVariable(node.left, 'e')) {
    const exponent = node.right;

    // Case: e^(a*t) like e^(4*t)
    if (isBinaryOp(exponent, '*') && isVariable(exponent.right, 't')) {
      return { a: exponent.left };
    }

    // NEW Case: e^(-a*t) like e^(-4*t)
    if (
      isUnaryMinus(exponent) &&
      isBinaryOp(exponent.operand, '*') &&
      isVariable(exponent.operand.right, 't')
    ) {
      // The operand is 'a*t', so we need '-a'
      return {
        a: {
          type: 'UnaryMinus',
          operand: exponent.operand.left,
        },
      };
    }

    // Case: e^t
    if (isVariable(exponent, 't')) {
      return { a: { type: 'Constant', value: 1 } };
    }

    // Case: e^(-t)
    if (isUnaryMinus(exponent) && isVariable(exponent.operand, 't')) {
      return { a: { type: 'Constant', value: -1 } };
    }
  }
  return null;
}

// Pattern matching for power of t like t^n
function matchPowerOfT(node: ASTNode): { n: number } | null {
  if (isBinaryOp(node, '^') && isVariable(node.left, 't') && isConstant(node.right)) {
    return { n: node.right.value };
  }
  if (isVariable(node, 't')) {
    return { n: 1 };
  }
  return null;
}

// Pattern matching for trigonometric functions like sin(a*t)
function matchTrigFunction(node: ASTNode, funcName: string): { a: ASTNode } | null {
  if (isFunction(node, funcName)) {
    if (isBinaryOp(node.argument, '*') && isVariable(node.argument.right, 't')) {
      return { a: node.argument.left };
    }
    if (isVariable(node.argument, 't')) {
      return { a: { type: 'Constant', value: 1 } };
    }
  }
  return null;
}

// --- Helper Functions for Pattern Matching ---

// Check if node is exp(-b t)
function isExpMinusBt(node: ASTNode): boolean {
  return node.type === "BinaryOperation" &&
         node.operator === "^" &&
         node.left.type === "Variable" &&
         node.left.name === "e" &&
         node.right.type === "UnaryMinus" &&
         node.right.operand.type === "BinaryOperation" &&
         node.right.operand.operator === "*" &&
         node.right.operand.right.type === "Variable" &&
         node.right.operand.right.name === "t" &&
         node.right.operand.left.type === "Constant";
}

// Extract b from exp(-b t)
function extractB(node: ASTNode): ASTNode {
  if (node.type === "UnaryMinus" &&
      node.operand.type === "BinaryOperation" &&
      node.operand.operator === "*" &&
      node.operand.right.type === "Variable" &&
      node.operand.right.name === "t") {
    return node.operand.left; // the constant b
  }
  throw new Error("Not an exp(-b t) form");
}

// Check if node is exp(-t/τ)
function isExpMinusTOverTau(node: ASTNode): boolean {
  return node.type === "BinaryOperation" &&
         node.operator === "^" &&
         node.left.type === "Variable" &&
         node.left.name === "e" &&
         node.right.type === "UnaryMinus" &&
         node.right.operand.type === "BinaryOperation" &&
         node.right.operand.operator === "/" &&
         node.right.operand.left.type === "Variable" &&
         node.right.operand.left.name === "t" &&
         node.right.operand.right.type === "Constant";
}

// Extract τ from exp(-t/τ)
function extractTau(node: ASTNode): ASTNode {
  if (node.type === "UnaryMinus" &&
      node.operand.type === "BinaryOperation" &&
      node.operand.operator === "/" &&
      node.operand.left.type === "Variable" &&
      node.operand.left.name === "t") {
    return node.operand.right; // denominator is τ
  }
  throw new Error("Not an exp(-t/τ) form");
}

// Pattern matcher for (e^{-b1 t} - e^{-b2 t}) / (b2 - b1)
function isDiffOfExponentialsWithBt(ast: ASTNode): boolean {
  return ast.type === "BinaryOperation" && ast.operator === "/" &&
         ast.left.type === "BinaryOperation" && ast.left.operator === "-" &&
         isExpMinusBt(ast.left.left) && isExpMinusBt(ast.left.right) &&
         ast.right.type === "BinaryOperation" && ast.right.operator === "-";
}

// Pattern matcher for (e^{-t/τ1} - e^{-t/τ2}) / (τ2 - τ1)
function isDiffOfExponentialsWithTau(ast: ASTNode): boolean {
  return ast.type === "BinaryOperation" && ast.operator === "/" &&
         ast.left.type === "BinaryOperation" && ast.left.operator === "-" &&
         isExpMinusTOverTau(ast.left.left) && isExpMinusTOverTau(ast.left.right) &&
         ast.right.type === "BinaryOperation" && ast.right.operator === "-";
}

// Core Laplace Transform Rules
const laplaceRules: LaplaceRule[] = [
  {
    name: "Constant or Fraction C → C/s",
    match: ast =>
      isConstant(ast) ||
      (isBinaryOp(ast, '/') && isConstant(ast.left) && isConstant(ast.right)),
    apply: ast => {
      if (isConstant(ast)) {
        // e.g. 2 → 2/s
        return BinOp("/", Const(ast.value), Var("s"));
      }
      if (isBinaryOp(ast, '/') && isConstant(ast.left) && isConstant(ast.right)) {
        // e.g. (5)/(2) → 5/(2*s)
        return BinOp("/", ast.left, BinOp("*", ast.right, Var("s")));
      }
      return null as any;
    }
  },
  {
    name: "Ramp: t → 1/s^2",
    match: ast => ast.type === "Variable" && ast.name === "t",
    apply: _ => BinOp("/", Const(1), BinOp("^", Var("s"), Const(2)))
  },
  {
    name: "Power: t^n → n!/s^(n+1)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "^" &&
                  ast.left.type === "Variable" && ast.left.name === "t" &&
                  ast.right.type === "Constant",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const n = (binOp.right as Constant).value;
      // factorial(n)
      const fact = (k: number): number => k <= 1 ? 1 : k * fact(k - 1);
      return BinOp("/",
        Const(fact(n)),
        BinOp("^", Var("s"), Const(n + 1))
      );
    }
  },
  {
    name: "Exponential: e^(-bt) or e^(-t/τ) → 1/(s+b)",
    match: ast => {
      // Match e^(-b t) where b could be a constant or 1 (implicit)
      if (ast.type === "BinaryOperation" &&
          ast.operator === "^" &&
          ast.left.type === "Variable" && ast.left.name === "e" &&
          ast.right.type === "UnaryMinus") {
        
        const exponent = ast.right.operand;
        
        // Case 1: e^(-b t) where b is explicit
        if (exponent.type === "BinaryOperation" &&
            exponent.operator === "*" &&
            exponent.right.type === "Variable" &&
            exponent.right.name === "t") {
          return true;
        }
        
        // Case 2: e^(-t) where b is implicit (1)
        if (exponent.type === "Variable" &&
            exponent.name === "t") {
          return true;
        }
        
        // Case 3: e^(-t/τ)
        if (exponent.type === "BinaryOperation" &&
            exponent.operator === "/" &&
            exponent.left.type === "Variable" &&
            exponent.left.name === "t") {
          return true;
        }
      }
      return false;
    },
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const unaryMinus = binOp.right as UnaryMinus;
      const exponent = unaryMinus.operand;
      
      let b: ASTNode;
      
      // Case 1: e^(-b t) - extract b
      if (exponent.type === "BinaryOperation" &&
          exponent.operator === "*" &&
          exponent.right.type === "Variable" &&
          exponent.right.name === "t") {
        b = exponent.left;
      }
      // Case 2: e^(-t) - b is 1
      else if (exponent.type === "Variable" &&
               exponent.name === "t") {
        b = Const(1);
      }
      // Case 3: e^(-t/τ) - the transform constant is 1/τ
      else if (exponent.type === "BinaryOperation" &&
               exponent.operator === "/" &&
               exponent.left.type === "Variable" &&
               exponent.left.name === "t") {
        const tau = exponent.right;
        b = BinOp("/", Const(1), tau);
      }
      else {
        throw new Error("Invalid exponential form for apply function");
      }
      
      return BinOp("/", Const(1), BinOp("+", Var("s"), b));
    }
  },
  {
    name: "Sin: sin(εt) → ε / (s^2 + ε^2)",
    match: ast => ast.type === "FunctionCall" &&
                  ast.name === "sin" &&
                  ast.argument.type === "BinaryOperation" &&
                  ast.argument.operator === "*" &&
                  ast.argument.right.type === "Variable" &&
                  ast.argument.right.name === "t",
    apply: ast => {
      const funcCall = ast as FunctionCall;
      const mulOp = funcCall.argument as BinaryOperation;
      const ε = mulOp.left;
      return BinOp("/",
        ε,
        BinOp("+", BinOp("^", Var("s"), Const(2)), BinOp("^", ε, Const(2)))
      );
    }
  },
  {
    name: "Cos: cos(εt) → s / (s^2 + ε^2)",
    match: ast => ast.type === "FunctionCall" &&
                  ast.name === "cos" &&
                  ast.argument.type === "BinaryOperation" &&
                  ast.argument.operator === "*" &&
                  ast.argument.right.type === "Variable" &&
                  ast.argument.right.name === "t",
    apply: ast => {
      const funcCall = ast as FunctionCall;
      const mulOp = funcCall.argument as BinaryOperation;
      const ε = mulOp.left;
      return BinOp("/",
        Var("s"),
        BinOp("+", BinOp("^", Var("s"), Const(2)), BinOp("^", ε, Const(2)))
      );
    }
  },
  {
    name: "Exp * Sin: e^(-bt) sin(εt) → ε / ((s+b)^2 + ε^2)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "*" &&
                  ast.left.type === "BinaryOperation" &&
                  ast.left.operator === "^" &&
                  ast.left.left.type === "Variable" && ast.left.left.name === "e" &&
                  ast.left.right.type === "UnaryMinus" &&
                  ast.left.right.operand.type === "BinaryOperation" &&
                  ast.left.right.operand.operator === "*" &&
                  ast.left.right.operand.right.type === "Variable" &&
                  ast.left.right.operand.right.name === "t" &&
                  ast.right.type === "FunctionCall" && ast.right.name === "sin",
    apply: ast => {
      const mainOp = ast as BinaryOperation;
      const expOp = mainOp.left as BinaryOperation;
      const unaryMinus = expOp.right as UnaryMinus;
      const mulOp = unaryMinus.operand as BinaryOperation;
      const b = mulOp.left;
      const sinFunc = mainOp.right as FunctionCall;
      const sinArg = sinFunc.argument as BinaryOperation;
      const ε = sinArg.left;
      return BinOp("/",
        ε,
        BinOp("+",
          BinOp("^", BinOp("+", Var("s"), b), Const(2)),
          BinOp("^", ε, Const(2))
        )
      );
    }
  },
  {
    name: "Exp * Cos: e^(-bt) cos(εt) → (s+b)/((s+b)^2 + ε^2)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "*" &&
                  ast.left.type === "BinaryOperation" &&
                  ast.left.operator === "^" &&
                  ast.left.left.type === "Variable" && ast.left.left.name === "e" &&
                  ast.left.right.type === "UnaryMinus" &&
                  ast.left.right.operand.type === "BinaryOperation" &&
                  ast.left.right.operand.operator === "*" &&
                  ast.left.right.operand.right.type === "Variable" &&
                  ast.left.right.operand.right.name === "t" &&
                  ast.right.type === "FunctionCall" && ast.right.name === "cos",
    apply: ast => {
      const mainOp = ast as BinaryOperation;
      const expOp = mainOp.left as BinaryOperation;
      const unaryMinus = expOp.right as UnaryMinus;
      const mulOp = unaryMinus.operand as BinaryOperation;
      const b = mulOp.left;
      const cosFunc = mainOp.right as FunctionCall;
      const cosArg = cosFunc.argument as BinaryOperation;
      const ε = cosArg.left;
      return BinOp("/",
        BinOp("+", Var("s"), b),
        BinOp("+",
          BinOp("^", BinOp("+", Var("s"), b), Const(2)),
          BinOp("^", ε, Const(2))
        )
      );
    }
  },
  {
    name: "Impulse δ(t) → 1",
    match: ast => ast.type === "FunctionCall" &&
                  ast.name === "delta" &&
                  ast.argument.type === "Variable" &&
                  ast.argument.name === "t",
    apply: _ => Const(1)
  },
  {
    name: "Diff of exponentials with time constants",
    match: isDiffOfExponentialsWithTau,
    apply: ast => {
      const mainOp = ast as BinaryOperation;
      const leftOp = mainOp.left as BinaryOperation;
      const tau1 = extractTau((leftOp.left as BinaryOperation).right);  // from first exponential
      const tau2 = extractTau((leftOp.right as BinaryOperation).right); // from second exponential

      return BinOp("/",
        Const(1),
        BinOp("*",
          BinOp("+", BinOp("*", tau1, Var("s")), Const(1)),
          BinOp("+", BinOp("*", tau2, Var("s")), Const(1))
        )
      );
    }
  },
  // Case: (e^{-b1 t} - e^{-b2 t}) / (b2 - b1)
  {
    name: "Diff of exponentials (rate constants)",
    match: ast =>
      ast.type === "BinaryOperation" && ast.operator === "/" &&
      ast.left.type === "BinaryOperation" && ast.left.operator === "-" &&
      isExpMinusBt(ast.left.left) && isExpMinusBt(ast.left.right),
    apply: ast => {
      const mainOp = ast as BinaryOperation;
      const leftOp = mainOp.left as BinaryOperation;
      const b1 = (((leftOp.left as BinaryOperation).right as UnaryMinus).operand as BinaryOperation).left;  // constant
      const b2 = (((leftOp.right as BinaryOperation).right as UnaryMinus).operand as BinaryOperation).left; // constant
      return BinOp("/",
        Const(1),
        BinOp("*",
          BinOp("+", Var("s"), b1),
          BinOp("+", Var("s"), b2)
        )
      );
    }
  },
  // Case: (e^{-t/τ1} - e^{-t/τ2}) / (τ2 - τ1)
  {
    name: "Diff of exponentials (time constants)",
    match: ast =>
      ast.type === "BinaryOperation" && ast.operator === "/" &&
      ast.left.type === "BinaryOperation" && ast.left.operator === "-" &&
      isExpMinusTOverTau(ast.left.left) && isExpMinusTOverTau(ast.left.right),
    apply: ast => {
      const mainOp = ast as BinaryOperation;
      const leftOp = mainOp.left as BinaryOperation;
      const tau1 = extractTau((leftOp.left as BinaryOperation).right);
      const tau2 = extractTau((leftOp.right as BinaryOperation).right);
      return BinOp("/",
        Const(1),
        BinOp("*",
          BinOp("+", BinOp("*", tau1, Var("s")), Const(1)),
          BinOp("+", BinOp("*", tau2, Var("s")), Const(1))
        )
      );
    }
  },
  {
    name: "sin(εt+η)",
    match: ast => ast.type === "FunctionCall" &&
                  ast.name === "sin" &&
                  ast.argument.type === "BinaryOperation" &&
                  ast.argument.operator === "+" &&
                  ast.argument.left.type === "BinaryOperation" &&
                  ast.argument.left.operator === "*" &&
                  ast.argument.left.right.type === "Variable" &&
                  ast.argument.left.right.name === "t",
    apply: ast => {
      const funcCall = ast as FunctionCall;
      const addOp = funcCall.argument as BinaryOperation;
      const mulOp = addOp.left as BinaryOperation;
      const ε = mulOp.left;
      const η = addOp.right;
      return BinOp("/",
        BinOp("+",
          BinOp("*", ε, { type: "FunctionCall", name: "cos", argument: η } as FunctionCall),
          BinOp("*", Var("s"), { type: "FunctionCall", name: "sin", argument: η } as FunctionCall)
        ),
        BinOp("+", BinOp("^", Var("s"), Const(2)), BinOp("^", ε, Const(2)))
      );
    }
  }
];

// Inverse Laplace Transform Rules
const inverseLaplaceRules: LaplaceRule[] = [
  {
    name: "1/s → 1 (unit step)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "Constant" && ast.left.value === 1 &&
                  ast.right.type === "Variable" && ast.right.name === "s",
    apply: _ => Const(1)
  },
  {
    name: "1/s^2 → t (ramp)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "Constant" && ast.left.value === 1 &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "^" &&
                  ast.right.left.type === "Variable" && ast.right.left.name === "s" &&
                  ast.right.right.type === "Constant" && ast.right.right.value === 2,
    apply: _ => Var("t")
  },
  {
    name: "n!/s^(n+1) → t^n",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "Constant" &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "^" &&
                  ast.right.left.type === "Variable" && ast.right.left.name === "s" &&
                  ast.right.right.type === "Constant",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const denomOp = binOp.right as BinaryOperation;
      const exp = (denomOp.right as Constant).value;
      const n = exp - 1; // since denom is s^(n+1)
      // check factorial consistency (optional)
      return BinOp("^", Var("t"), Const(n));
    }
  },
  {
    name: "1/(s+b) → e^(-bt)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "Constant" && ast.left.value === 1 &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "+" &&
                  ast.right.left.type === "Variable" && ast.right.left.name === "s",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const denomOp = binOp.right as BinaryOperation;
      const b = denomOp.right;
      return BinOp("^", Var("e"),
        UnaryMinus(BinOp("*", b, Var("t")))
      );
    }
  },
  {
    name: "ε/(s^2 + ε^2) → sin(εt)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "+" &&
                  ast.right.left.type === "BinaryOperation" &&
                  ast.right.left.operator === "^" &&
                  ast.right.left.left.type === "Variable" &&
                  ast.right.left.left.name === "s",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const ε = binOp.left;
      return { type: "FunctionCall", name: "sin",
               argument: BinOp("*", ε, Var("t")) } as FunctionCall;
    }
  },
  {
    name: "s/(s^2 + ε^2) → cos(εt)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "Variable" && ast.left.name === "s" &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "+" &&
                  ast.right.left.type === "BinaryOperation" &&
                  ast.right.left.operator === "^" &&
                  ast.right.left.left.type === "Variable" &&
                  ast.right.left.left.name === "s",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const denomOp = binOp.right as BinaryOperation;
      const rightTerm = denomOp.right as BinaryOperation;
      const ε = rightTerm.left; // ε from ε^2
      return { type: "FunctionCall", name: "cos",
               argument: BinOp("*", ε, Var("t")) } as FunctionCall;
    }
  },
  {
    name: "ε/((s+b)^2 + ε^2) → e^(-bt) sin(εt)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "+" &&
                  ast.right.left.type === "BinaryOperation" &&
                  ast.right.left.operator === "^" &&
                  ast.right.left.left.type === "BinaryOperation" &&
                  ast.right.left.left.operator === "+" &&
                  ast.right.left.left.left.type === "Variable" &&
                  ast.right.left.left.left.name === "s",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const ε = binOp.left;
      const denomOp = binOp.right as BinaryOperation;
      const leftTerm = denomOp.left as BinaryOperation;
      const sPlusB = leftTerm.left as BinaryOperation;
      const b = sPlusB.right;
      return BinOp("*",
        BinOp("^", Var("e"), UnaryMinus(BinOp("*", b, Var("t")))),
        { type: "FunctionCall", name: "sin",
          argument: BinOp("*", ε, Var("t")) } as FunctionCall
      );
    }
  },
  {
    name: "(s+b)/((s+b)^2 + ε^2) → e^(-bt) cos(εt)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "BinaryOperation" &&
                  ast.left.operator === "+" &&
                  ast.left.left.type === "Variable" && ast.left.left.name === "s" &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "+" &&
                  ast.right.left.type === "BinaryOperation" &&
                  ast.right.left.operator === "^" &&
                  ast.right.left.left.type === "BinaryOperation" &&
                  ast.right.left.left.operator === "+" &&
                  ast.right.left.left.left.type === "Variable" &&
                  ast.right.left.left.left.name === "s",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const numerator = binOp.left as BinaryOperation;
      const b = numerator.right; // s+b
      const denomOp = binOp.right as BinaryOperation;
      const rightTerm = denomOp.right as BinaryOperation;
      const ε = rightTerm.left; // ε from ε^2
      return BinOp("*",
        BinOp("^", Var("e"), UnaryMinus(BinOp("*", b, Var("t")))),
        { type: "FunctionCall", name: "cos",
          argument: BinOp("*", ε, Var("t")) } as FunctionCall
      );
    }
  },
  {
    name: "1 → δ(t)",
    match: ast => ast.type === "Constant" && ast.value === 1,
    apply: _ => ({ type: "FunctionCall", name: "delta", argument: Var("t") } as FunctionCall)
  },
  {
    name: "1/((s+b1)(s+b2)) → (e^{-b2t} - e^{-b1t})/(b1-b2)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "Constant" && ast.left.value === 1 &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "*" &&
                  ast.right.left.type === "BinaryOperation" &&
                  ast.right.left.operator === "+" &&
                  ast.right.right.type === "BinaryOperation" &&
                  ast.right.right.operator === "+",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const denomOp = binOp.right as BinaryOperation;
      const leftFactor = denomOp.left as BinaryOperation;
      const rightFactor = denomOp.right as BinaryOperation;
      const b1 = leftFactor.right;
      const b2 = rightFactor.right;
      return BinOp("/",
        BinOp("-",
          BinOp("^", Var("e"), UnaryMinus(BinOp("*", b2, Var("t")))),
          BinOp("^", Var("e"), UnaryMinus(BinOp("*", b1, Var("t"))))
        ),
        BinOp("-", b1, b2)
      );
    }
  },
  {
    name: "(ε cos η + s sin η)/(s^2+ε^2) → sin(εt+η)",
    match: ast => ast.type === "BinaryOperation" &&
                  ast.operator === "/" &&
                  ast.left.type === "BinaryOperation" &&
                  ast.left.operator === "+" &&
                  ast.right.type === "BinaryOperation" &&
                  ast.right.operator === "+" &&
                  ast.right.left.type === "BinaryOperation" &&
                  ast.right.left.operator === "^" &&
                  ast.right.left.left.type === "Variable" &&
                  ast.right.left.left.name === "s",
    apply: ast => {
      const binOp = ast as BinaryOperation;
      const numerator = binOp.left as BinaryOperation; // ε cos η + s sin η
      const leftTerm = numerator.left as BinaryOperation; // ε * cos(η)
      const rightTerm = numerator.right as BinaryOperation; // s * sin(η)
      
      // Extract ε from the first term
      const ε = leftTerm.left;
      // Extract η from cos(η) in the first term
      const cosFunc = leftTerm.right as FunctionCall;
      const η = cosFunc.argument;
      
      return {
        type: "FunctionCall",
        name: "sin",
        argument: BinOp("+",
          BinOp("*", ε, Var("t")),
          η
        )
      } as FunctionCall;
    }
  }
];

// --- Core Transform Functions ---

// Forward transform: f(t) -> F(s)
// --- Core Transform Functions ---

// Forward transform: f(t) -> F(s)
function transformForward(node: ASTNode | null): ASTNode | null {
  if (!node) return null;

  // --- START: Linearity and High-Level Rules ---

  // Handle linearity: L{a*f(t) + b*g(t)} = a*L{f(t)} + b*L{g(t)}
  if (isBinaryOp(node) && (node.operator === '+' || node.operator === '-')) {
    const leftTransform = transformForward(node.left);
    const rightTransform = transformForward(node.right);
    if (leftTransform && rightTransform) {
      return {
        type: 'BinaryOperation',
        operator: node.operator,
        left: leftTransform,
        right: rightTransform
      };
    }
    return null;
  }

  // Handle constant multiplication: L{c*f(t)} = c*L{f(t)}
  if (isBinaryOp(node) && node.operator === '*') {
    // Case: left is a pure fraction constant
    if (isBinaryOp(node.left, '/') && isConstant(node.left.left) && isConstant(node.left.right)) {
      const funcTransform = transformForward(node.right);
      if (funcTransform) {
        return BinOp('*', node.left, funcTransform);
      }
    }
    // Case: right is a pure fraction constant
    if (isBinaryOp(node.right, '/') && isConstant(node.right.left) && isConstant(node.right.right)) {
      const funcTransform = transformForward(node.left);
      if (funcTransform) {
        return BinOp('*', node.right, funcTransform);
      }
    }
    // Existing integer constant cases
    if (isConstant(node.left)) {
      const funcTransform = transformForward(node.right);
      if (funcTransform) {
        return BinOp('*', node.left, funcTransform);
      }
    } else if (isConstant(node.right)) {
      const funcTransform = transformForward(node.left);
      if (funcTransform) {
        return BinOp('*', node.right, funcTransform);
      }
    }
  }

  // Handle division by a constant: L{f(t)/c} = (1/c) * L{f(t)}
  if (isBinaryOp(node, '/') && isConstant(node.right)) {
    const funcTransform = transformForward(node.left);
    if (funcTransform) {
        return {
            type: 'BinaryOperation',
            operator: '/',
            left: funcTransform,
            right: node.right,
        };
    }
  }

  // Handle unary minus: L{-f(t)} = -L{f(t)}
  if (isUnaryMinus(node)) {
    const operandTransform = transformForward(node.operand);
    if (operandTransform) {
      return { type: 'UnaryMinus', operand: operandTransform };
    }
  }

  // Handle frequency shifting: L{e^(a*t) * f(t)} = F(s-a)
  if (isBinaryOp(node, '*')) {
      let expNode: ASTNode | null = null;
      let otherNode: ASTNode | null = null;

      const expMatchLeft = matchExponential(node.left);
      if (expMatchLeft) {
          expNode = node.left;
          otherNode = node.right;
      } else {
          const expMatchRight = matchExponential(node.right);
          if (expMatchRight) {
              expNode = node.right;
              otherNode = node.left;
          }
      }

      if (expNode && otherNode) {
          const expMatch = matchExponential(expNode)!;
          const a = expMatch.a;
          
          const fOfS = transformForward(otherNode);
          if (fOfS) {
              const sMinusA: BinaryOperation = {
                  type: 'BinaryOperation',
                  operator: '-',
                  left: { type: 'Variable', name: 's' },
                  right: a
              };
              return substitute(fOfS, 's', sMinusA);
          }
      }
  }

  // --- END: Linearity and High-Level Rules ---

  // --- START: Rule-Based Transform ---

  // Apply core transform rules
  for (const rule of laplaceRules) {
    if (rule.match(node)) {
      console.log("Matched rule:", rule.name);
      return rule.apply(node);
    }
  }

  // Additional rule: δ(t) -> 1
  if (isFunction(node, 'delta') && isVariable(node.argument, 't')) {
    return { type: 'Constant', value: 1 };
  }

  // If no rule matches, return null
  return null;
}

// Inverse transform: F(s) -> f(t)
function transformInverse(node: ASTNode | null): ASTNode | null {
  if (!node) return null;

  // Handle linearity for inverse transforms
  if (isBinaryOp(node) && (node.operator === '+' || node.operator === '-')) {
    const leftTransform = transformInverse(node.left);
    const rightTransform = transformInverse(node.right);
    if (leftTransform && rightTransform) {
      return {
        type: 'BinaryOperation',
        operator: node.operator,
        left: leftTransform,
        right: rightTransform
      };
    }
    return null;
  }

  // Handle constant multiplication
  if (isBinaryOp(node) && node.operator === '*') {
    if (isConstant(node.left)) {
      const funcTransform = transformInverse(node.right);
      if (funcTransform) {
        return {
          type: 'BinaryOperation',
          operator: '*',
          left: node.left,
          right: funcTransform
        };
      }
    }
  }

  // Apply inverse transform rules
  for (const rule of inverseLaplaceRules) {
    if (rule.match(node)) {
      console.log("Matched inverse rule:", rule.name);
      return rule.apply(node);
    }
  }

  // Additional rule: 1 -> δ(t)
  if (isConstant(node, 1)) {
    return {
      type: 'FunctionCall',
      name: 'delta',
      argument: { type: 'Variable', name: 't' }
    };
  }

  return null;
}

const simplifyLatex = (latex: string): string => {
  if (!latex) return '';

  let simplified = latex;
  
  // NEW: Replace "--" with "+"
  // Example: s - -1  => s + 1
  simplified = simplified.replace(/-\s*-/g, '+ ');

  // Rule 1: Replace "+ (-number)" with "- number"
  // Example: s + (-2)  =>  s - 2
  simplified = simplified.replace(/\+\s*\((-[^)]+)\)/g, (match, capturedGroup) => {
    return `- ${capturedGroup.substring(1)}`;
  });

  // Rule 2: Replace "- (-number)" with "+ number"
  // Example: s - (-2)  =>  s + 2
  simplified = simplified.replace(/-\s*\((-[^)]+)\)/g, (match, capturedGroup) => {
    return `+ ${capturedGroup.substring(1)}`;
  });
  
  // Rule 3: Remove parentheses around single terms like (3) or (a) but not (a+b)
  simplified = simplified.replace(/\(([^()]+)\)/g, (match, content) => {
    // If content has addition or a binary subtraction, keep the parentheses
    if (content.includes('+') || content.includes(' - ')) {
        return match; // Keep parens for expressions like (a+b) or (a - b)
    }
    // Otherwise, it's a single term, so remove the parens
    return content;
  });

  return simplified.trim(); // Added trim() for good measure
};

const transform = (input: string, from: string): string => {
  console.log("=== transform START ===");
  console.log(`Input raw: ${JSON.stringify(input)}, from: ${from}`);

  if (!input || input.trim() === '') {
    console.log("Empty input, returning ''");
    return '';
  }

  try {
    const processedInput = preprocessLatex(input);
    console.log("Processed input for parse:", processedInput);

    const ast = parse(processedInput);
    console.log("Parsed AST:", JSON.stringify(ast, null, 2));

    if (!ast) {
      console.log("Failed to parse input into AST");
      return '';
    }

    let transformedAst: ASTNode | null = null;
    if (from === 't') {
      transformedAst = transformForward(ast);
    } else if (from === 's') {
      transformedAst = transformInverse(ast);
    }
    console.log("Transformed AST:", JSON.stringify(transformedAst, null, 2));

    if (!transformedAst) {
      console.log("No transform rule matched");
      return '';
    }
    
    // --- NEW SIMPLIFICATION STEP ---
    const simplifiedAst = simplifyAst(transformedAst);
    console.log("Simplified AST:", JSON.stringify(simplifiedAst, null, 2));

    const result = astToLatex(simplifiedAst); // Use the simplified AST
    console.log("Raw LaTeX result:", result);

    const simplified = simplifyLatex(result);
    console.log("Simplified result:", simplified);

    console.log("=== transform END ===");
    return simplified;
  } catch (error) {
    console.error("Error in transform:", error);
    return '';
  }
};


// --- Main Application Component ---
export default function LaplaceTransformPage() {
    const [isUpdating, setIsUpdating] = useState(false);
    const isUpdatingRef = useRef(false);
    const tMathQuillEl = useRef<HTMLDivElement>(null);
    const sMathQuillEl = useRef<HTMLDivElement>(null);
    const tMathField = useRef<MathQuillField | null>(null);
    const sMathField = useRef<MathQuillField | null>(null);

    useEffect(() => {
        const loadScript = (id: string, src: string, onLoad: () => void) => {
            if (document.getElementById(id)) {
                if(onLoad) onLoad();
                return;
            }
            const script = document.createElement('script');
            script.id = id;
            script.src = src;
            script.async = false;
            script.onload = onLoad;
            document.head.appendChild(script);
        };
        
        const loadCss = (id: string, href: string) => {
             if (document.getElementById(id)) return;
             const link = document.createElement('link');
             link.id = id;
             link.rel = 'stylesheet';
             link.href = href;
             document.head.appendChild(link);
        }

        const initializeMathQuill = () => {
            if (window.MathQuill && tMathQuillEl.current && sMathQuillEl.current && !tMathField.current) {
                const MQ = window.MathQuill.getInterface(2);
                
                tMathField.current = MQ.MathField(tMathQuillEl.current, {
                    handlers: {
                        edit: function() {
                            if (isUpdatingRef.current) return;
                            
                            isUpdatingRef.current = true;
                            setIsUpdating(true);
                            
                            setTimeout(() => {
                                if (tMathField.current && sMathField.current) {
                                    const latex = tMathField.current.latex();
                                    console.log('t field latex:', latex);
                                    const result = transform(latex, 't');
                                    console.log('t field transform result:', result);
                                    sMathField.current.latex(result || '');
                                }
                                isUpdatingRef.current = false;
                                setIsUpdating(false);
                            }, 0);
                        }
                    }
                });
                
                sMathField.current = MQ.MathField(sMathQuillEl.current, {
                    handlers: {
                        edit: function() {
                            if (isUpdatingRef.current) return;
                            
                            isUpdatingRef.current = true;
                            setIsUpdating(true);
                            
                            setTimeout(() => {
                                if (tMathField.current && sMathField.current) {
                                    const latex = sMathField.current.latex();
                                    console.log('s field latex:', latex);
                                    const result = transform(latex, 's');
                                    console.log('s field transform result:', result);
                                    tMathField.current.latex(result || '');
                                }
                                isUpdatingRef.current = false;
                                setIsUpdating(false);
                            }, 0);
                        }
                    }
                });
                
                // Set initial value
                if (tMathField.current) {
                    tMathField.current.latex('t');
                }
            }
        }
        
        loadCss('mathquill-css', "https://cdnjs.cloudflare.com/ajax/libs/mathquill/0.10.1/mathquill.css");
        loadScript('jquery-script', "https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js", () => {
             loadScript('mathquill-script', "https://cdnjs.cloudflare.com/ajax/libs/mathquill/0.10.1/mathquill.min.js", initializeMathQuill);
        });

    }, []); // Empty dependency array - only run once on mount

  return (
    <div className="flex min-h-screen w-full items-center justify-center bg-background text-foreground font-sans p-4">
      <Card className="w-full max-w-3xl">
        <CardHeader>
          <CardTitle className="text-2xl">
            Laplace Transform Calculator
          </CardTitle>
          {/* Instruction removed as requested */}
        </CardHeader>
        <CardContent>
          <div className="flex flex-col md:flex-row items-center justify-center space-y-4 md:space-y-0 md:space-x-4">
            {/* T-Space Input */}
            <div className="w-full space-y-2">
              <Label className="text-lg">
                f(t) - Time Domain
              </Label>
               <div ref={tMathQuillEl} className="math-input w-full p-3 border rounded-md min-h-[60px] bg-background text-foreground text-lg font-mono focus-within:ring-2 focus-within:ring-ring" />
            </div>

            {/* Icon */}
            <div className="px-2 pt-8">
                <ArrowRightLeft className="h-8 w-8 text-muted-foreground" />
            </div>
            
            {/* S-Space Input */}
            <div className="w-full space-y-2">
              <Label className="text-lg">
                F(s) - Frequency Domain
              </Label>
              <div ref={sMathQuillEl} className="math-input w-full p-3 border rounded-md min-h-[60px] bg-background text-foreground text-lg font-mono focus-within:ring-2 focus-within:ring-ring" />
            </div>
          </div>
        </CardContent>
      </Card>
      <style jsx global>{`
        .math-input {
            font-size: 1.2rem;
        }
        .mq-editable-field .mq-cursor {
            border-left: 2px solid hsl(var(--primary));
        }
        .mq-editable-field {
            color: hsl(var(--foreground));
        }
        /* Fix cursor visibility in dark mode */
        [data-theme="dark"] .mq-editable-field .mq-cursor,
        .dark .mq-editable-field .mq-cursor {
            border-left: 2px solid #ffffff !important;
        }
        /* Ensure cursor is visible in light mode too */
        [data-theme="light"] .mq-editable-field .mq-cursor,
        .light .mq-editable-field .mq-cursor {
            border-left: 2px solid #000000 !important;
        }
        /* Fallback for systems without theme attribute */
        @media (prefers-color-scheme: dark) {
            .mq-editable-field .mq-cursor {
                border-left: 2px solid #ffffff !important;
            }
        }
        @media (prefers-color-scheme: light) {
            .mq-editable-field .mq-cursor {
                border-left: 2px solid #000000 !important;
            }
        }
      `}</style>
    </div>
  )
}

