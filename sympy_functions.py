# -*- coding: utf-8 -*-
"""
Started on Wed Apr  9 13:39:41 2025

@author: Vincent MAHE

Sympy functions useful to the use of MMS functions.
"""

#%% Imports
import copy
from sympy import sqrt, solve

#%% Functions
def sub_deep(expr, sub):
    """
    Performs deep substitutions of an expression. This is needed when a substitution involves terms that can still be substituted.
    For instance, one wants to substitute a1 and a2 by expressions, but a1 is actually a function of a2, so at least 2 substitutions are required.

    Parameters
    ----------
    expr : expr
        Expression on which substitutions are to be performed.
    sub : list of tuples
        The substitutions to perform.

    Returns
    -------
    expr_sub : expr
        The expression with substitutions performed.

    """
    expr_sub  = copy.copy(expr) 
    expr_init = 0
    while expr_init != expr_sub: # Check if the 2 are the same
        expr_init = copy.copy(expr_sub)       # Update expr_init -> now the 2 are the same
        expr_sub  = expr_sub.subs(sub).doit() # Update expr_sub  -> now the two are different if substitutions were performed
    
    return expr_sub

def solve_poly2(poly,x):
    """
    Finds the roots of a polynomial of degree 2 of the form
    p(x) = a*x**2 + b*x + c
    Note that b can be null but not a or c
    
    It is a workaround to using solve() or solveset(). 
    These two work but can be very long when coefficients a, b, c 
    are expressions involving many parameters.
    Note that solve() is significantly slower than solveset() .
    
    Parameters
    ----------
    poly: Expr
        polynomial whose roots are to be computed
    x: Symbol
        Variable of the polynomial
        
    Returns
    -------
    x_sol: list of Expr
        list containing the two roots of the polynomial
    """
    
    # Polynomial terms
    dic_x = poly.expand().collect(x, evaluate=False)
    keys = set(dic_x.keys())
    
    # Ensure the polynomial is written only with positive powers of x
    min_power = min(list(keys), key=lambda expr: get_exponent(expr, x))
    min_expo = get_exponent(min_power, x)
    if min_expo<=0:
        poly = (poly/min_power).expand()
    
    dic_x = poly.expand().collect(x, evaluate=False)
    keys = set(dic_x.keys())
    
    # Solve
    if keys == set([x**2, x, 1]):
        a     = dic_x[x**2].factor()
        b     = dic_x[x].factor()
        c     = dic_x[1]
        D     = (b**2 - 4*a*c).factor()
        T1    = (-b/(2*a)).factor()
        T2    = sqrt(D)/(2*a)
        x1    = T1 - T2
        x2    = T1 + T2
        x_sol = [x1,x2]
        
    elif keys == set([x**2, 1]):
        a        = dic_x[x**2].factor()
        c        = dic_x[1]
        T2       = sqrt(- 4*a*c)/(2*a)
        x_sol    = [-T2, T2]

    elif keys == set([x, 1]) or keys == set([x**2, x]):
        x_sol = solve(poly.expand(), x)

    else:
        x_sol = []
        print('Trying to use solve_poly2() with a polynomial different from p(x) = a*x**2 + b*x + c')
        
    return x_sol

def get_exponent(expr, x):
    """
    Get the exponent of x in an expression
    """
    # This assumes expr is a power of x
    if expr.is_Number:
        return 0
    elif expr == x:
        return 1
    elif expr.is_Pow and expr.base == x:
        return expr.exp
    
    else:
        return float('inf')  # Handle unexpected expressions
    

# %%
