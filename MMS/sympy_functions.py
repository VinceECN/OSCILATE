# -*- coding: utf-8 -*-
"""
Started on Wed Apr  9 13:39:41 2025

@author: Vincent MAHE

Sympy functions useful to the use of MMS functions.
"""

#%% Imports
import copy
from sympy import sqrt, solve, lambdify
import warnings

#%% Functions
def sub_deep(expr, sub):
    r"""
    Performs deep substitutions of an expression. 
    
    Parameters
    ----------
    expr : sympy.Expr
        Expression on which substitutions are to be performed.
    sub : list of tuples
        The substitutions to perform.

    Returns
    -------
    expr_sub : sympy.Expr
        The expression with substitutions performed.

    Notes
    -----
    Deep substitutions are needed when a substitution involves terms that can still be substituted.
    For instance, one wants to substitute :math:`a_1` and :math:`a_2` by expressions, but :math:`a_1` is actually a function of :math:`a_2`, so at least 2 substitutions are required.

    """
    expr_sub  = copy.copy(expr) 
    expr_init = 0
    while expr_init != expr_sub: # Check if the 2 are the same
        expr_init = copy.copy(expr_sub)       # Update expr_init -> now the 2 are the same
        expr_sub  = expr_sub.subs(sub).doit() # Update expr_sub  -> now the two are different if substitutions were performed
    
    return expr_sub

def solve_poly2(poly, x):
    r"""
    Finds the roots of a polynomial of degree 2. 

    Parameters
    ----------
    poly: sympy.Expr
        polynomial whose roots are to be computed
    x: sympy.Symbol
        Variable of the polynomial
        
    Returns
    -------
    x_sol: list of sympy.Expr
        list containing the two roots of the polynomial

    Notes
    -----
    The polynomial of degree 2 takes the form

    .. math::
        p(x) = a x^2 + bx + c.
    
    Note that :math:`b` can be null but not :math:`a` nor :math:`c`.
    
    It is a workaround to using :func:`~sympy.solvers.solvers.solve` or :func:`~sympy.solvers.solveset.solveset`. 
    These two work but can be very long when coefficients :math:`a,\; b,\; c` are expressions involving many parameters.
    Note that :func:`~sympy.solvers.solvers.solve` is significantly slower than :func:`~sympy.solvers.solveset.solveset`.
    
    """

    # Check the solvability
    if not check_solvability(poly, x):
        print("The polynomial cannot be solved")
        return False

    # Polynomial terms
    dic_x = polynomial_terms(poly, x)
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

def polynomial_terms(poly, x):
    r"""
    Identify the terms of a polynomial. 
    
    Parameters
    ----------
    poly : sympy.Expr
        The polynomial considered.
    x : sympy.Symbol
        The variable to solve for.

    Returns
    -------
    dic_x: dict
        The polynomial terms.

    Notes
    -----
    If the expression given for poly is of the form
    
    .. math::
        p(x) = q(x) x^{-n},

    where the powers of :math:`x` in :math:`q(x)` are all superior or equal to :math:`0`,
    then an auxiliary polynomial 
    
    .. math::
        P(x) = \dfrac{p(x)}{x^{-n}} 
    
    is constructed. It is the terms of that positive powers polynomial that are returned.
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
    r"""
    Check the solvability of a polynomial :math:`p(x)`.

    Parameters
    ----------
    poly : sympy.Expr
        The polynomial considered.
    x : sympy.Symbol
        The variable to solve for.

    Returns
    -------
    bool : bool,
        True is solvable, False otherwise.
    """
    dic_x = polynomial_terms(poly, x)
    poly_terms = set(dic_x.keys())
    min_power  = min(poly_terms, key=lambda expr: get_exponent(expr, x))
    poly_terms = set([poly_term/min_power for poly_term in poly_terms])

    if poly_terms in [set([x**2, x, 1]), set([x**2, 1]), set([x, 1])]:
        return True
    else:
        return False
    
def get_exponent(expr, x):
    r"""
    Get the exponent of :math:`x` in an expression of the type :math:`\lambda x^n` where :math:`\lambda` is a constant while :math:`n` is an integer or rational.

    Parameters
    ----------
    expr: sympy.Expr
        The expression in which one wants to identify the exponent of x.
    x: sympy.Symbol
        The variable whose exponent is to be known.
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
    r"""
    Generate a list of :math:`(i, j)` indices for all elements in the diagonal blocks of a block-diagonal matrix.

    Parameters
    ----------
    matrix: sympy.Matrix
        The matrix to check for block diagonality.
    block_sizes: int or list of int
        Size(s) of the diagonal blocks.

    Returns
    -------
    indices : list
        A list of tuples `(i, j)` representing the indices of elements in the diagonal blocks.
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
    Check if a matrix is block-diagonal given block sizes.

    Parameters
    ----------
    matrix: sympy.Matrix
        The matrix to check for block diagonality.
    block_sizes: int or list of int
        Size(s) of the diagonal blocks.

    Returns
    -------
        bool: `True` if the matrix is block-diagonal, `False` otherwise.
    """
    n = matrix.rows
    if isinstance(block_sizes, int):
        block_sizes = [block_sizes]*(n // block_sizes)
    if sum(block_sizes) != n:
        return False

    # Get the block-diagonal elements indices
    indices_diag_blocks = get_block_diagonal_indices(matrix, block_sizes)

    # Iterate over the matrix elements
    for i in range(n):
        for j in range(n):
            # Skip elements inside the diagonal blocks
            if (i, j) not in indices_diag_blocks:
                if matrix[i, j] != 0:
                    return False

    return True

def sympy_to_numpy(expr_sy, param):
    """
    Transform a sympy expression into a numpy array.

    Parameters
    ----------
    expr_sy : sympy.Expr
        A sympy expression.
    param : dict
        A dictionnary whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter
        
        2. The numerical value(s) taken by that parameter

    Returns
    -------
    expr_np : numpy.ndarray
        The numerical values taken by the sympy expression evaluated.
    """
    
    args, values = zip(*param.values())
    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")
        expr_np = lambdify(args, expr_sy, modules="numpy")(*values)

    return expr_np