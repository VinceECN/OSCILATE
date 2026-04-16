# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the dynamical system.
"""

#%% Imports and initialisation
from sympy import sympify, Symbol, Function, Expr
from typing import Union, TYPE_CHECKING


#%% Classes and functions
class Forcing:
    r"""
    Define the forcing on the system as

    - A forcing amplitude `F`,
    
    - Forcing coefficients `fF`, used to introduce parametric forcing or simply weight the harmonic forcing.
    
    For the :math:`i^\textrm{th}` oscillator, denoting `fF[i]` as :math:`f_{F,i}(\boldsymbol{x}(t), \dot{\boldsymbol{x}}(t), \ddot{\boldsymbol{x}}(t))`, 
    the forcing term on that oscillator is :math:`f_{F,i} F \cos(\omega t)`.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        F : Symbol
        fF: list[Union[Expr, int]]
    
    def __init__(self, F, fF):
        self.F       = F
        self.fF = fF

class Dynamical_system:
    r"""
    The dynamical system studied. 
    See :ref:`dyn_sys` for a detailed description of the dynamical system.

    Parameters
    ----------
    t : sympy.Symbol
        time :math:`t`.
    x : sympy.Function or list of sympy.Function
        Unknown(s) of the problem.
    Eq : sympy.Expr or list of sympy.Expr
        System's equations without forcing, which can be defined separately (see parameters `F` and `fF`).
        Eq is the unforced system of equations describing the system's dynamics. 
    omegas : sympy.Symbol or list of sympy.Symbol
        The natural frequency of each oscillator.
    F : sympy.Symbol or 0, optional
        Forcing amplitude :math:`F`. 
        Default is 0.
    fF : sympy.Expr or list of sympy.Expr, optional
        For each oscillator, specify the coefficient multiplying the forcing terms in the equation.
        It can be used to define parametric forcing. Typically, if the forcing is :math:`x F \cos(\omega t)`, then ``fF = x``.
        Default is a list of 1, so the forcing is direct. 
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        Eq:      list[Expr]
        forcing: Forcing
        ndof:    int
        omegas:  list[Symbol]
        t:       Symbol
        x:       list[Function]
        
    def __init__(self, t, x, Eq, omegas, F = 0, fF = None):
        r"""
        Initialisation of the dynamical system.
        """
        
        # Information
        print('Creation of the dynamical system')
        
        # Time
        self.t = t
        
        # Variables and equations
        if isinstance(x, list):
            self.ndof = len(x)
            self.x  = x
            self.Eq = Eq
            self.omegas = omegas
        else:
            self.ndof = 1
            self.x    = [x]
            self.Eq   = [Eq]
            self.omegas = [omegas]
            
        # Forcing
        F  = sympify(F)
        if fF == None:
            fF = [1]*self.ndof
        
        if not isinstance(fF, list): 
            fF = [fF]
        for ix, coeff in enumerate(fF):
            if isinstance(coeff, int):
                fF[ix] = sympify(coeff)
        self.forcing = Forcing(F, fF)
        
