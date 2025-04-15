import re
import sys

def replace_frac_manually(expr: str) -> str:
    # Find the first occurrence of "\frac"
    idx = expr.find(r'\frac')
    if idx == -1:
        return expr  # No fraction found, return as-is.
    
    # Find the numerator: first '{' after "\frac"
    start_num = expr.find('{', idx)
    if start_num == -1:
        return expr  # Malformed fraction.
    stack = []
    i = start_num
    while i < len(expr):
        if expr[i] == '{':
            stack.append(i)
        elif expr[i] == '}':
            stack.pop()
            if not stack:
                end_num = i
                break
        i += 1
    else:
        return expr  # Unbalanced braces.
    
    numerator = expr[start_num+1:end_num]
    
    # Find the denominator: first '{' after numerator.
    start_den = expr.find('{', end_num)
    if start_den == -1:
        return expr  # Malformed fraction.
    stack = []
    i = start_den
    while i < len(expr):
        if expr[i] == '{':
            stack.append(i)
        elif expr[i] == '}':
            stack.pop()
            if not stack:
                end_den = i
                break
        i += 1
    else:
        return expr
    
    denominator = expr[start_den+1:end_den]
    
    # Recursively process numerator and denominator in case they contain nested fractions.
    new_num = replace_frac_manually(numerator)
    new_den = replace_frac_manually(denominator)
    
    # Build the replacement string.
    replacement = f"({new_num})/({new_den})"
    
    # Replace the whole "\frac{...}{...}" with the replacement.
    new_expr = expr[:idx] + replacement + expr[end_den+1:]
    
    # Process the new expression in case more fractions exist.
    return replace_frac_manually(new_expr)

def fully_replace_frac(expr: str) -> str:
    return replace_frac_manually(expr)

def latex2python(latex_str: str) -> str:
    expr = latex_str.strip('$')
    expr = re.sub(r'\\left', '', expr)
    expr = re.sub(r'\\right', '', expr)
    expr = fully_replace_frac(expr)
    expr = re.sub(r'\\operatorname\{([^{}]+)\}', r'\1', expr)
    expr = re.sub(r'\\?exp\(', 'e**(', expr)
    expr = re.sub(r'\\(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh|log|ln)', r'\1', expr)
    expr = re.sub(r'([a-zA-Z0-9_]+)\s*\^([a-zA-Z0-9_]+)', r'\1**(\2)', expr)
    expr = re.sub(r'\\cdot\s*', '*', expr)
    expr = re.sub(r'(\S+)\s*\^\{\s*([^}]+?)\s*\}', r'\1**(\2)', expr)
    
    # Simple substitution: replace { with (, } with ), and ^ with **
    expr = expr.replace('{', '(').replace('}', ')').replace('^', '**')
    return expr

def python_to_numpy(expr: str) -> str:
    # First handle all power operations recursively
    expr = re.sub(r'\be\*\*\(([^)]+)\)', r'np.exp(\1)', expr)
    while '**' in expr:
        expr = re.sub(r'([a-zA-Z0-9_.()]+)\*\*([a-zA-Z0-9_.()]+|\([^)]+\))', r'np.power(\1, \2)', expr)
    
    # Then handle other numpy conversions
    expr = re.sub(r'(?<!np\.)\bln\(([^)]+)\)', r'np.log(\1)', expr)
    expr = re.sub(r'(?<!np\.)\blog\(([^)]+)\)', r'np.log10(\1)', expr)
    expr = re.sub(r'\b(sin|cos|tan|sinh|cosh|tanh|arcsin|arccos|arctan|arcsinh|arccosh|arctanh)\(([^)]+)\)', r'np.\1(\2)', expr)
    
    return expr

def main():
    if len(sys.argv) > 1:
        latex_str = sys.argv[1]
        python_expr = latex2python(latex_str)
        numpy_expr = python_to_numpy(python_expr)
        
        print(f"Raw LaTeX: {latex_str}")
        print(f"Python Expression: {python_expr}")
        print(f"NumPy Expression: {numpy_expr}")
        
        # Generate full Python code that uses the equation
        full_code = f"""
# Python code that uses your equation
import numpy as np
from matplotlib import pyplot as plt

# Your equation converted to Python
equation = lambda x: {python_expr}

# Example usage
x = np.linspace(-10, 10, 1000)
y = equation(x)

# Plot
plt.figure(figsize=(10, 6))
plt.plot(x, y)
plt.grid(True)
plt.title("Plot of your equation")
plt.xlabel("x")
plt.ylabel("y")
plt.axhline(y=0, color='k', linestyle='--', alpha=0.3)
plt.axvline(x=0, color='k', linestyle='--', alpha=0.3)
plt.show()
"""
        print("\nFull Python code:")
        print(full_code)
    else:
        print("No LaTeX string provided.")

if __name__ == "__main__":
    main()
