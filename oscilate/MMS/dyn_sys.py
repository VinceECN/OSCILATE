# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the dynamical system.
"""

#%% Imports and initialisation
from sympy import sympify

#%% Classes and functions
class Dynamical_system:
    r"""
    The dynamical system studied.

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

    Notes
    -----
    Systems considered are typically composed of :math:`N` coupled nonlinear equations of the form

    .. math::
        \begin{cases}
        \ddot{x}_0 + \omega_0^2 x_0 & = f_0(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t), \\
        & \vdots \\
        \ddot{x}_{N-1} + \omega_{N-1}^2 x_{N-1} & = f_{N-1}(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t).
        \end{cases}
    
    The :math:`x_i(t)` (:math:`i=0,...,N-1`) are the oscillators' coordinates (dof for degrees of freedom), 

    .. math::

        \boldsymbol{x}(t)^\intercal = [x_0(t), x_1(t), \cdots, x_{N-1}(t)]
         
    is the vector containing all the oscillators' coordinates (:math:`^\intercal` denotes the transpose), 
    :math:`\omega_i` are their natural frequencies, 
    :math:`t` is the time, 
    :math:`\dot{(\bullet)} = \textrm{d}(\bullet)/\textrm{d}t` denotes a time-derivative. 
    :math:`f_i` are functions which can contain:

    - **Weak linear terms** in :math:`x_i`, :math:`\dot{x}_i`, or :math:`\ddot{x}_i`.
    
    - **Weak linear coupling terms** involving :math:`x_j`, :math:`\dot{x}_j`, or :math:`\ddot{x}_j` with :math:`j \neq i`.
    
    - **Weak nonlinear terms**. Taylor expansions are performed to approximate nonlinear terms as polynomial nonlinearities.
    
    - **Forcing terms**, which can be:
    
        - *Hard* (appearing at leading order) or *weak* (small).
        
        - Primarily harmonic, e.g., :math:`F \cos(\omega t)`, where :math:`F` and :math:`\omega` are the forcing amplitude and frequency, respectively.
        
        - Modulated by any function (constant, linear, or nonlinear) to model parametric forcing (e.g., :math:`x_i(t) F \cos(\omega t)`).

    Internal resonance relations among oscillators can be specified in a second step by expressing the :math:`\omega_i` as a function of a reference frequency. 
    Detuning can also be introduced during this step.
    """
    
    def __init__(self, t, x, Eq, omegas, **kwargs):
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
        F  = kwargs.get("F", sympify(0))
        fF = kwargs.get("fF", [1]*self.ndof)
        if not isinstance(fF, list): 
            fF = [fF]
        for ix, coeff in enumerate(fF):
            if isinstance(coeff, int):
                fF[ix] = sympify(coeff)
        self.forcing = Forcing(F, fF)
        
class Forcing:
    r"""
    Define the forcing on the system as

    - A forcing amplitude `F`,
    
    - Forcing coefficients `fF`, used to introduce parametric forcing or simply weight the harmonic forcing.
    
    For the :math:`i^\textrm{th}` oscillator, denoting `fF[i]` as :math:`f_{F,i}(\boldsymbol{x}(t), \dot{\boldsymbol{x}}(t), \ddot{\boldsymbol{x}}(t))`, 
    the forcing term on that oscillator is :math:`f_{F,i} F \cos(\omega t)`.
    """
    
    def __init__(self, F, fF):
        self.F       = F
        self.fF = fF
