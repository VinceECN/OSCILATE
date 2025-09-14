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

def solve_poly2(poly, x):
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

    # Check the solvability
    if not check_solvability(poly, x):
        print("The polynomial cannot be solved")
        return False

    # Polynomial terms
    dic_x = poly_positive_pow(poly, x)
    keys  = set(dic_x.keys())

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

def poly_positive_pow(poly, x):
    """
    Identify the terms of a polynomial. If the expression given for poly is of the form
    p(x) = q(x)*x**(-n), where the powers of x in q(x) are all superior or equal to 0,
    then an auxiliary polynomial P(x) = p(x)/x**-n is constructed. 

    Parameters
    ----------
    poly : Expr
        The polynomial considered.
    x : Symbol
        The variable to solve for.

    Returns
    -------
    dic_x: dict
        The polynomial terms.
    """
    # Polynomial terms
    dic_x = poly.expand().collect(x, evaluate=False)
    keys = set(dic_x.keys())

    # Increase the polynomial order if it contains negative powers of x so the lowest possible order is x**0=1
    min_power = min(list(keys), key=lambda expr: get_exponent(expr, x))
    min_expo = get_exponent(min_power, x)
    if min_expo<=0:
        poly = (poly/min_power).expand()
    
    # Terms of the increased-order polynomial
    dic_x = poly.expand().collect(x, evaluate=False)

    return dic_x

def check_solvability(poly, x):
    """
    Check the solvability of a polynomial.

    Parameters
    ----------
    poly : Expr
        The polynomial considered.
    x : Symbol
        The variable to solve for.

    Returns
    -------
    bool : True is solvable, False otherwise.
    """
    dic_x = poly_positive_pow(poly, x)
    poly_terms = set(dic_x.keys())
    min_power  = min(poly_terms, key=lambda expr: get_exponent(expr, x))
    poly_terms = set([poly_term/min_power for poly_term in poly_terms])

    if poly_terms in [set([x**2, x, 1]), set([x**2, 1]), set([x, 1])]:
        return True
    else:
        return False
    
def get_exponent(expr, x):
    """
    Get the exponent of :math:`x` in an expression
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

def get_block_diagonal_indices(matrix, block_sizes):
    """
    Generate a list of :math:`(i, j)` indices for all elements in the diagonal blocks of a block-diagonal matrix.

    Parameters
    ----------
    matrix: Matrix.
        The matrix to check for block diagonality.
    block_sizes: int or list of int
        Size(s) of the diagonal blocks.

    Returns
    -------
    indices : list
        A list of tuples (i, j) representing the indices of elements in the diagonal blocks.
    """

    if isinstance(block_sizes, int):
        block_sizes = [block_sizes]*(matrix.rows // block_sizes)

    indices = []
    start = 0

    for size in block_sizes:
        end = start + size
        # Iterate over the current block
        for ii in range(start, end):
            for jj in range(start, end):
                indices.append((ii, jj))
        start = end

    return indices

def is_block_diagonal(matrix, block_sizes):
    """
    Check if a SymPy matrix is block-diagonal given block sizes.

    Parameters
    ----------
    matrix: Matrix.
        The matrix to check for block diagonality.
    block_sizes: int or list of int
        Size(s) of the diagonal blocks.

    Returns
    -------
        bool: True if the matrix is block-diagonal, False otherwise.
    """
    n = matrix.rows
    if isinstance(block_sizes, int):
        block_sizes = [block_sizes]*(n // block_sizes)
    if sum(block_sizes) != n:
        return False

    # Get the block-diagonal elements indices
    indices_diag_blocks = get_block_diagonal_indices(matrix, block_sizes)

    # Initialize the starting index of the current block
    start = 0

    # Iterate over the matrix elements
    for i in range(n):
        for j in range(n):
            # Skip elements inside the diagonal blocks
            if (i, j) not in indices_diag_blocks:
                if matrix[i, j] != 0:
                    return False

    return True