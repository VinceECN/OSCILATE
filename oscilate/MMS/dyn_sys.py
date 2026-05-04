# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the dynamical system.
"""

#%% Imports and initialisation
from sympy import sympify, Symbol, Function, Expr, I
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
        Eqz:     list[Expr]
        forcing: Forcing
        form:    Union[str, list[str]]
        ndof:    int
        omegas:  list[Symbol]
        sub_z:   list[tuple]
        t:       Symbol
        x:       list[Function]
        z:       list[Function]
        
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

        # System form
        self.form = "oscillator"
        
    def complex_form(self):
        r"""
        Rewrite the dynamical system in complex form, resulting in a system of complex 1st order coupled ODEs from the initial real 2nd order ODEs. 
        See :ref:`dyn_sys` for details.
        """

        # Create the complex coordinates
        self.complex_coordinates()
            
        # Rewrite the equations in complex form
        self.Eqz = []
        for Eqi, omegai in zip(self.Eq, self.omegas):
            self.Eqz.append( (-I/(2*omegai) * Eqi.subs(self.sub_z)).simplify() )

        # Rewrite the forcing terms
        self.forcing.fFz= []
        for fFi, omegai in zip(self.forcing.fF, self.omegas):
            self.forcing.fFz.append( (-I/(2*omegai) * fFi.subs(self.sub_z)).simplify() )
        
        # System form
        self.form = [self.form, "complex"]

    def complex_coordinates(self):
        r"""
        Introduce the complex coordinates, denoted :math:`z_i` for oscillator :math:`i`, and defined as
        
        .. math::
            \begin{cases}
            x_i(t)       & = z_i(t) + \bar{z}_i(t), \\
            \dot{x}_i(t) & = \textrm{j} \omega_{i} (z_i(t) - \bar{z}_i(t)),
            \end{cases}

        where :math:`\bar{\bullet}` denotes the transpose. 
        """

        # Create the complex coordinates
        self.z     = [] # Complex coordinates
        self.sub_z = [] # Substitutions from x (real coordinates) to z
        for ix, (xi, omegai) in enumerate(zip(self.x, self.omegas)):
            zi = Function(r'z_{}'.format(ix), complex=True)(self.t)
            self.z.append(zi)
            self.sub_z += [(xi.diff(self.t, 2), 2*I*omegai*zi.diff(self.t) + omegai**2*(zi - zi.conjugate())),
                           (xi.diff(self.t)   , I*omegai*(zi - zi.conjugate())),
                           (xi                , zi + zi.conjugate())] 
            