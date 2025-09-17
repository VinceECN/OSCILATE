# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
"""

#%% Imports and initialisation
from sympy import (exp, I, conjugate, re, im, Rational, 
                   symbols, Symbol, Function, solve, dsolve,
                   cos, sin, tan, srepr, sympify, simplify, 
                   zeros, det, trace, eye, Mod, lambdify)
from sympy.simplify.fu import TR5, TR8, TR10
from . import sympy_functions as sfun
import numpy as np
import itertools
import warnings
import matplotlib.pyplot as plt

#%% Classes and functions
class Dynamical_system:
    r"""
    The dynamical system studied.

    Systems considered are typically composed of :math:`N` coupled nonlinear equations of the form

    .. math::
        \begin{cases}
        \ddot{x}_0 + \omega_0^2 x_0 & = f_0(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t), \\
        & \vdots \\
        \ddot{x}_{N-1} + \omega_{N-1}^2 x_{N-1} & = f_{N-1}(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t).
        \end{cases}
    
    The :math:`x_i(t)` (:math:`i=0,...,N-1`) are the oscillators' coordinates, 
    :math:`\omega_i` are their natural frequencies, 
    :math:`\boldsymbol{x}` is the vector containing all the oscillators' coordinates, 
    :math:`t` is the time, 
    :math:`\dot{(\bullet)}` denotes a time-derivative :math:`d(\bullet)/dt`, 
    :math:`f_i` is a function which can contain:

    - Linear terms in :math:`x_i`, :math:`\dot{x}_i` or :math:`\ddot{x}_i`, typically those that will be considered small in the MMS,
    
    - Weak coupling terms in :math:`x_j`, :math:`\dot{x}_j` or :math:`\ddot{x}_j`, :math:`j\neq i`,
    
    - Weak nonlinear terms. Only polynomial nonlinearities are supported. Taylor expansions are performed if nonlinearities are not polynomial,
    
    - Forcing, which can be hard (at first order) or weak (small). Harmonic and parametric forcing are supported.

    Internal resonance relations among oscillators can be specified in a second step by expressing the :math:`\omega_i` as a function of a reference frequency. Detuning can also be introduced during this step.
    """
    
    def __init__(self, t, x, Eq, omegas, **kwargs):
        r"""
        Initialisation of the dynamical system.

        Parameters
        ----------
        t : Symbol
            time :math:`t`.
        x : Function or list of Function
            Unknown(s) of the problem.
        Eq : Expr or list of Expr
            System's equations without forcing, which can be defined separately (see `F` and `f_coeff`).
            Eq is the unforced system of equations describing the system's dynamics. 
        omegas : Symbol or list of Symbol
            The natural frequency of each oscillator.
        F : Symbol or 0, optional
            Forcing amplitude :math:`F`. 
            Default is 0.
        f_coeff : Expr of list of Expr, optional
            For each dof, specify the coefficient multiplying the forcing terms in the equation.
            It can be used to define parametric forcing. Typically, if the forcing is :math:`x F \cos(\omega t)`, then f_coeff = x.
            Default is a list of 1, so the forcing is direct. 
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
        F       = kwargs.get("F", sympify(0))
        f_coeff = kwargs.get("f_coeff", [1]*self.ndof)
        if not isinstance(f_coeff, list): 
            f_coeff = [f_coeff]
        for ix, coeff in enumerate(f_coeff):
            if isinstance(coeff, int):
                f_coeff[ix] = sympify(coeff)
        self.forcing = Forcing(F, f_coeff)
        
class Forcing:
    r"""
    Define the forcing on the system as

    - A forcing amplitude `F`,
    
    - Forcing coefficients `f_coeff`, used to introduce parametric forcing.
    
    For the :math:`i^\textrm{th}` oscillator, denoting `f_coeff[i]` as :math:`\Gamma_i(\boldsymbol{x}(t), \dot{\boldsymbol{x}}(t), \ddot{\boldsymbol{x}}(t))`, 
    the forcing term on that oscillator is :math:`\Gamma_i F \cos(\omega t)`.
    """
    
    def __init__(self, F, f_coeff):
        self.F       = F
        self.f_coeff = f_coeff

        
def scale_parameters(param, scaling, eps):
    r"""
    Scale parameters with the scaling parameter :math:`\epsilon`.
    For a given parameter :math:`p` and a scaling order :math:`\lambda`, the associated scaled parameter :math:`\tilde{p}` is 

    .. math::
        p = \epsilon^{\lambda} \tilde{p} .
    

    Parameters
    ----------
    param : list of Symbol and/or Function
        Unscaled parameters.
    scaling : list of int or float
        The scaling for each parameter.
    eps : Symbol
        Small parameter :math:`\epsilon`.

    Returns
    -------
    param_scaled: list of Symbol and/or Function
        Scaled parameters.
    sub_scaling: list of 2 lists of tuple
        Substitutions from scaled to unscaled parameters and vice-versa. 

        - :math:`1^{\text{st}}` list: The substitutions to do to introduce the scaled parameters in an expression.
        
        - :math:`2^{\text{nd}}` list: The substitutions to do to reintroduce the unscaled parameters in a scaled expression.
    """
    
    param_scaled     = []
    sub_scaling      = [[], []]
    
    for ii, (p, pow_p) in enumerate(zip(param, scaling)):
        if isinstance(p, Symbol):
            param_scaled.append(symbols(r"\tilde{{{}}}".format(p.name), **p.assumptions0))
        elif isinstance(p, Function):
            param_scaled.append(Function(r"\tilde{{{}}}".format(p.name), **p.assumptions0)(*p.args))
            
        sub_scaling[0].append( (p, param_scaled[ii] * eps**pow_p) )
        sub_scaling[1].append( (param_scaled[ii], p / eps**pow_p) )
        
    return param_scaled, sub_scaling
        
        
class Multiple_scales_system:
    r"""
    The multiple scales system.

    The starting point is to introduce asymptotic series and multiple time scales in the initial dynamical system. 
    The solution for oscillator :math:`i` is sought as a series expansion up to order :math:`N_e` (for a leading order term :math:`\epsilon^0 = 1`). This expansion takes the form
        
    .. math::
        x_i(t) = x_{i,0}(t) + \epsilon x_{i,1}(t) + \epsilon^2 x_{i,2}(t) + \cdots + \epsilon^{N_e} x_{i,N_e}(t).

    Time scales are introduced as follows:
    
    .. math::
        t_0 = t, \; t_1 = \epsilon t, \; t_2 = \epsilon^2 t, \cdots, t_{N_e} = \epsilon^{N_e} t,
    
    where :math:`t_0` is the fast time, i.e. the time used to describe the oscillations, 
    while :math:`t_1, \; t_2,\; t_{N_e}` are slow times, associated to amplitude and phase variations of the solutions in time.
    
    Introducing the asymptotic series and time scales in the initial dynamical system (see :class:`~MMS.MMS.Dynamical_system`) results in :math:`N_e+1` dynamical systems, each one appearing at different orders of :math:`\epsilon`.
    Denoting time scales derivatives as :math:`\textrm{D}_i(\bullet) = \partial (\bullet) / \partial t_i`, introducing the vectors of asymptotic coordinates 
    
    .. math::
        \boldsymbol{x}_i^\intercal = [x_{0,i}, x_{1,i}, \cdots, x_{N-1, i}],

    containing all the asymptotic terms of order :math:`i`, with :math:`\intercal` denoting the transpose, and defining the stackings

    .. math::
        \begin{split}
        \boldsymbol{x}_{(p)}^\intercal                 = & [\boldsymbol{x}_0, \boldsymbol{x}_1, \cdots, \boldsymbol{x}_{p}], \\
        \textrm{D}_j\boldsymbol{x}_{(p)}^\intercal     = & [\textrm{D}_j\boldsymbol{x}_0, \textrm{D}_j\boldsymbol{x}_1, \cdots, \textrm{D}_j\boldsymbol{x}_{p}], \\
        \textrm{D}_{(q)}\boldsymbol{x}_i^\intercal     = & [\textrm{D}_0\boldsymbol{x}_i, \textrm{D}_1\boldsymbol{x}_i, \cdots, \textrm{D}_q\boldsymbol{x}_i], \\
        \textrm{D}_{(q)}\boldsymbol{x}_{(p)}^\intercal = & [\textrm{D}_0\boldsymbol{x}_{(p)}, \textrm{D}_1\boldsymbol{x}_{(p)}, \cdots, \textrm{D}_q\boldsymbol{x}_{(p)}], 
        \end{split}
        
    the MMS equations can be written as

    .. math::
        \begin{aligned}
        \begin{cases}
        \textrm{D}_0 x_{0,0}   + \omega_0^2 x_{0,0}   & = f_{0,0}(t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,0} + \omega_0^2 x_{N-1,0} & = f_{N-1,0}(t_0, t_1),
        \end{cases} \\
        \begin{cases}
        \textrm{D}_0 x_{0,1}   + \omega_0^2 x_{0,1}   & = f_{0,1}(\boldsymbol{x}_0, \textrm{D}_{(1)} \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,1} + \omega_0^2 x_{N-1,1} & = f_{N-1,1}(\boldsymbol{x}_0, \textrm{D}_{(1)} \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, t_0, t_1),
        \end{cases} \\
        \begin{cases}
        \textrm{D}_0 x_{0,2}   + \omega_0^2 x_{0,2}   & = f_{0,2}  (\boldsymbol{x}_{(1)}, \textrm{D}_{(1)} \boldsymbol{x}_{(1)}, \textrm{D}_0^2 \boldsymbol{x}_{(1)}, \textrm{D}_1 \boldsymbol{x}_0, t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,2} + \omega_0^2 x_{N-1,2} & = f_{N-1,2}(\boldsymbol{x}_{(1)}, \textrm{D}_{(1)} \boldsymbol{x}_{(1)}, \textrm{D}_0^2 \boldsymbol{x}_{(1)}, \textrm{D}_1 \boldsymbol{x}_0, t_0, t_1),
        \end{cases} \\
        \vdots \\
        \begin{cases}
        \textrm{D}_0 x_{0,N_e}   + \omega_0^2 x_{0,N_e}   & = f_{0,N_e}(t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,N_e} + \omega_0^2 x_{N-1,N_e} & = f_{N-1,N_e}(t_0, t_1),
        \end{cases}
        \end{aligned}


    The solutions for the :math:`x_i` will be sought around at frequency :math:`\omega`, defined as 
    
    .. math::
        \omega = \omega_{\textrm{MMS}} + \epsilon \sigma,
    
    where
    :math:`\omega_{\textrm{MMS}}` is the *central* MMS frequency and :math:`\sigma` is a detuning about that frequency. 


    """
    
    def __init__(self, dynamical_system, eps, Ne, omega_ref, sub_scaling, 
                 ratio_omegaMMS=1, eps_pow_0=0, **kwargs):
        r"""
        Transform the dynamical system introducing asymptotic series and multiple time scales. 

        Parameters
        ----------
        dynamical_system : Dynamical_system
            The dynamical system. 
        eps : Symbol
            Small perturbation parameter :math:`\epsilon`.
        Ne : int
            Truncation order of the asymptotic and order of the slowest time scale.
        omega_ref : Symbol
            Reference frequency :math:`\omega_{\textrm{ref}}` of the MMS. 
            Not necessarily the frequency around which the MMS is going to be applied, see `ratio_omegaMMS`.
        sub_scaling : list of tuples
            Substitutions to do to scale the equations. 
            Links small parameters to their scaled counterpart through :math:`\epsilon`.
        ratio_omegaMMS : int or Rational, optional
            Specify the frequency `omegaMMS` around which the MMS is going to be applied in terms of :math:`\omega_{\textrm{ref}}`.
            Denoting `ratio_omegaMMS` as :math:`r_{\textrm{MMS}}`, this means that
            
            .. math::
                \omega_{\textrm{MMS}} = r_{\textrm{MMS}} \omega_{\textrm{ref}}.
            
            Use `Rational(p,q)` for 
            
            .. math::
                q \omega_{\textrm{MMS}} = p \omega_{\textrm{ref}}
            
            to get better-looking results than the float :math:`p/q`.
            Default is 1.
        eps_pow_0 : int, optional
            Order of the leading-order term in the asymptotic series of each oscillators' response.
            For the :math:`i^{\textrm{th}}` dof and denoting `eps_pow_0` as :math:`\lambda_0`, this means that
            .. math::
                x_i = \epsilon^{\lambda_0} x_{i,0} + \epsilon^{\lambda_0+1} x_{i,1} + \cdots.
            
            Default is 0.
        ratio_omega_osc : list of int or Rational or None, optional
            Specify the natural frequencies of the oscillators :math:`\omega_i` in terms of the reference frequency :math:`\omega_{\textrm{ref}}`. 
            Denoting `ratio_omega_osc[i]` as :math:`r_i`, this means that 
            
            .. math::
                \omega_i \approx r_i \omega_{\textrm{ref}}.
            
            Default is `None` for each oscillator, so the :math:`\omega_i` are arbitrary and there are no internal resonances.
            Detuning can be introduced through the `detunings` keyword argument. 
        detunings : list of Symbol or int, optional
            The detuning of each oscillator. Denoting `detunings[i]` as :math:`\delta_i`, this means that 
            
            .. math::
            \omega_i = r_i \omega_{\textrm{ref}} + \delta_i.
            
            Defaut is 0 for each oscillator. 
        """
        
        # Information
        print('Initialisation of the multiple scales sytem')
        
        # Order of the method
        self.eps = eps
        self.Ne  = Ne
        self.eps_pow_0 = eps_pow_0
        
        # MMS reference frequency
        self.omega_ref = omega_ref
        
        # MMS frequencies of interest
        self.ratio_omegaMMS = ratio_omegaMMS
        self.omegaMMS       = self.ratio_omegaMMS * self.omega_ref
        self.sigma          = symbols(r'\sigma',real=True) # Detuning to investigate the response around omegaMMS
        self.omega          = symbols(r'\omega',real=True) # Frequencies investigated 
        sub_omega           = (self.omega, self.omegaMMS + self.eps*self.sigma) # Definition of omega
        sub_sigma           = (self.sigma, (self.omega - self.omegaMMS)/self.eps) # Definition of sigma (equivalent to that of omega)
        
        # Number of dof
        self.ndof = dynamical_system.ndof
        
        # Multiple time scales
        self.t = dynamical_system.t
        self.tS, sub_t = self._time_scales()
        
        # Asymptotic series of x
        self.xMMS, sub_xMMS_t, sub_x = self._asymptotic_series(dynamical_system, eps_pow_0=self.eps_pow_0)
        
        # Substitutions required
        self.sub = Substitutions_MMS(sub_t, sub_xMMS_t, sub_x, sub_scaling, sub_omega, sub_sigma)
    
        # Forcing
        self.forcing = self._forcing_MMS(dynamical_system)
        
        # Oscillators' frequencies (internal resonances and detuning)
        self.omegas = dynamical_system.omegas
        self.ratio_omega_osc = kwargs.get("ratio_omega_osc", [None]   *self.ndof)
        self.detunings       = kwargs.get("detunings",       [0]*self.ndof)
        self._oscillators_frequencies()        
        
        # Coordinates
        self.coord = Coord_MMS(self)
        self.polar_coordinates()
        
        # Solutions
        self.sol = Sol_MMS()
        
        # Compute the MMS equations
        self.compute_EqMMS(dynamical_system)
        
    def _time_scales(self):
        r"""
        Create the time scales 
        
        .. math::
            t_i = \epsilon^i t, \quad i=0, 1, ...
        
        and prepare substitutions from the physical time :math:`t` to these time scales :math:`t_i`.
        """
        
        tS  = []
        sub_t = []
        for it in range(self.Ne+1):
            tS.append(symbols(r't_{}'.format(it), real=True, positive=True))
            sub_t.append((self.eps**it * self.t, tS[it]))
            
        sub_t.reverse() # Start substitutions with the slowest time scale
        
        return tS, sub_t
    
    def _asymptotic_series(self, dynamical_system, eps_pow_0=0):
        r"""
        Define series expansions for each dof and prepare substitutions from 

        1. dof `x` to temporary t-dependent asymptotic terms `xMMS_t`

        2. Temporary `xMMS_t` to the time scales-dependent `xMMS` 

        For the :math:`i^\textrm{th}` oscillator and with `eps_pow_0=0`, stage 1 introduces temporary :math:`t`-dependent terms such that 
        
        .. math::
            x_i(t) = x_{i,0}(t) + \epsilon x_{i,1}(t) + \epsilon^2 x_{i,2}(t) + \cdots, 
        
        while stage 2 introduces the :math:`t_i`-dependent terms :math:`x_{i0}(\boldsymbol{t}),\; x_{i,1}(\boldsymbol{t}),\; x_{i,2}(\boldsymbol{t}),\; \cdots`, where 
        :math:`\boldsymbol{t} = [t_0, t_1, t_2, ...]` is the vector containing the time scales, such that
        
        .. math::
            x_i(\boldsymbol{t}) = x_{i,0}(\boldsymbol{t}) + \epsilon x_{i,1}(\boldsymbol{t}) + \epsilon^2 x_{i,2}(\boldsymbol{t}) + \cdots. 
        
        """
        
        # Initialisation
        xMMS         = [] # Terms x00, x01, ..., x10, x11, ... of the asymptotic series of the xi
        sub_xMMS_t   = [] # Substitutions from xMMS(t) to xMMS(*tS)
        x_expanded   = [] # x in terms of xMMS(t)
        sub_x        = [] # Substitutions from x to xMMS(t)
        
        for ix in range(self.ndof):
            
            # Initialisations 
            xMMS.append([])      # A list that will contain the different expansion orders of the current x
            xMMS_t = []          # Temporary xMMS(t) -> depend on the physical time t
            x_expanded.append(0) # Initialise the current x to 0
            
            for it in range(self.Ne+1):
            
                # Define time-dependent asymptotic terms
                xMMS_t.append(Function(r'x_{{{},{}}}'.format(ix,it), real=True)(self.t))
                x_expanded[ix] += self.eps**(it+eps_pow_0) * xMMS_t[it]
                
                # Define time scales-dependent asymptotic terms
                xMMS[ix].append(Function(xMMS_t[it].name, real=True)(*self.tS))
                
                # Substitutions from xMMS(t) and its time derivatives to xMMS(*tS) and its time scales derivatives
                sub_xMMS_t.extend( [(xMMS_t[it].diff(self.t,2), Chain_rule_d2fdt2(xMMS[ix][it], self.tS, self.eps)), 
                                    (xMMS_t[it].diff(self.t,1), Chain_rule_dfdt  (xMMS[ix][it], self.tS, self.eps)), 
                                    (xMMS_t[it]               , xMMS[ix][it])] )
            
            # Substitutions from x to xMMS(t)
            sub_x.append((dynamical_system.x[ix], x_expanded[ix]))
        
        return xMMS, sub_xMMS_t, sub_x
        
    def _forcing_MMS(self, dynamical_system):
        r"""
        Rewrite the forcing terms :math:`\Gamma_i(\boldsymbol{x}(t), \dot{\boldsymbol{x}}(t), \ddot{\boldsymbol{x}}(t)) F \cos(\omega t)`, :math:`i=1,...,N` where :math:`N` is the number of oscillators.
        This involves

        1. Replacing the :math:`x_i(t)` by their series expansions written in terms of time scales,
        
        2. Scaling the forcing and the parameters in :math:`\Gamma_i` if any,
        
        3. Truncating terms whose order is larger than the largest order retained in the MMS,
        
        4. Rewrite the :math:`\cos(\omega t)` as 
        
        .. math::
            \cos(\omega t) = \frac{1}{2} e^{\mathrm{i}(\omega_{\textrm{MMS}} + \epsilon \sigma)t} + cc = \frac{1}{2} e^{\mathrm{i}(\omega_{\textrm{MMS}}t_0 + \sigma t_1)} + cc.
        

        Parameters
        ----------
        dynamical_system : Dynamical_system
            The dynamical system.
        """
        
        # Get the expression of the forcing frequency
        omega = self.sub.sub_omega[1]
        
        # Get the forcing amplitude and order for the substitutions
        forcing = False
        for item in self.sub.sub_scaling:
            if dynamical_system.forcing.F in item:
                dic_F_MMS = item[1].collect(self.eps, evaluate=False)
                scaling_coeff = list(dic_F_MMS.keys())[0]
                forcing = True
                if scaling_coeff == 1:
                    f_order = 0
                elif scaling_coeff == self.eps:
                    f_order = 1 
                else:
                    f_order = scaling_coeff.args[1]
                    
                F       = list(dic_F_MMS.values())[0]
                forcing = True
                
        if not forcing:
            F       = 0
            f_order = self.eps_pow_0+self.Ne
            
        # Get the forcing term for each dof
        forcing_term = []
        f_coeff      = []
        sub_t, sub_x, sub_xMMS_t, sub_scaling = list(map(self.sub.__dict__.get,["sub_t", "sub_x", "sub_xMMS_t", "sub_scaling"]))

        for ix in range(self.ndof):
            f_coeff.append( (dynamical_system.forcing.f_coeff[ix].subs(self.sub.sub_scaling)
                             .subs(sub_x).doit().subs(sub_xMMS_t).expand().subs(sub_t).doit())
                            .series(self.eps, n=self.eps_pow_0+self.Ne+1).removeO())
            
            forcing_term.append( (f_coeff[ix] * Rational(1,2)*F*self.eps**f_order)
                                .series(self.eps, n=self.eps_pow_0+self.Ne+1).removeO() * 
                                (exp( I*((omega*self.t).expand().subs(sub_t).doit())) + 
                                 exp(-I*((omega*self.t).expand().subs(sub_t).doit())) ) 
                                )
        
        forcing = Forcing_MMS(F, f_order, f_coeff, forcing_term)
    
        return forcing
    
    def _oscillators_frequencies(self):
        r"""
        Gives the expression of every oscillator frequency in terms of the reference frequency, possibly with a detuning.
        For the :math:`i^\textrm{th}` oscillator, this corresponds to expression its frequency as 
        
        .. math::
            \omega_i = r_i \omega_{\textrm{ref}} + \delta_i .
        
        An associated first-order natural frequency :math:`\omega_{i,0}` is defined by neglecting the detuning :math:`\delta_i`, which is at least of order :math:`\epsilon`, resulting in
        
        .. math::
            \omega_{i,0} = r_i \omega_{\textrm{ref}}.
        
        """
        
        self.sub.sub_omegas = [] # Substitutions from the omegas to their expression in terms of omega_ref
        self.omegas_O0      = [] # Leading order oscillators' natural frequencies
        for ix in range(self.ndof):
            
            # Check if ratio_omega_osc should be modified
            if self.ratio_omega_osc[ix] is None:
                if Mod(self.omegas[ix], self.omega_ref)==0:
                    self.ratio_omega_osc[ix] = self.omegas[ix] // self.omega_ref
                elif Mod(self.omega_ref, self.omegas[ix])==0:
                    self.ratio_omega_osc[ix] = Rational(1, self.omega_ref // self.omegas[ix])
            
            if self.ratio_omega_osc[ix] is not None:
                self.sub.sub_omegas.append( (self.omegas[ix], self.ratio_omega_osc[ix]*self.omega_ref + self.detunings[ix]) )
                self.omegas_O0.append( self.ratio_omega_osc[ix] * self.omega_ref )
            else:
                self.sub.sub_omegas.append( (self.omegas[ix], self.omegas[ix]) )
                self.omegas_O0.append( self.omegas[ix] )
        
        
    def compute_EqMMS(self, dynamical_system):
        r"""
        Compute the system of equations at each order of :math:`\epsilon` for each dof.
        The output `EqMMS` is a list of lists:

        - The :math:`1^{\text{st}}` level lists are associated to the orders of :math:`\epsilon` from the lowest to the highest order,
        
        - The :math:`2^{\text{nd}}` level lists are associated to the equations for each oscillator.

        Parameters
        ----------
        dynamical_system : Dynamical_system
            The dynamical system.
        """
    
        # Equations with every epsilon appearing
        sub_t, sub_x, sub_xMMS_t, sub_scaling, sub_omegas = list(map(self.sub.__dict__.get,["sub_t", "sub_x", "sub_xMMS_t", "sub_scaling", "sub_omegas"]))
    
        Eq_eps = []
        for ix in range(self.ndof):
            Eq_eps.append( ((dynamical_system.Eq[ix].expand().subs(sub_omegas).doit().subs(sub_scaling).doit()
                             .subs(sub_x).doit().subs(sub_xMMS_t).doit().expand().subs(sub_t).doit())
                          .series(self.eps, n=self.eps_pow_0+self.Ne+1).removeO() 
                          - self.forcing.forcing_term[ix]).expand())
            
            if self.eps_pow_0 != 0: # Set the leading order to eps**0 = 1
                Eq_eps[-1] = (Eq_eps[-1] / self.eps**(self.eps_pow_0)).expand()
                
        # MMS equations system
        EqMMS = []
        for ix in range(self.ndof):
            
            # Initialise a list of the equations at each order. Start with the lowest order
            EqMMSO = [Eq_eps[ix].series(self.eps, n=1).removeO()] 
            
            # What has to be substracted to keep only the terms of order eps**io in equation at order io.
            retrieve_EqMMSO = EqMMSO[0] 
            
            # Feed EqMMSO with higher orders of epsilon
            for io in range(1, self.Ne+1):
                EqMMSO.append( ((Eq_eps[ix].series(self.eps, n=io+1).removeO() - retrieve_EqMMSO) / self.eps**io).simplify().expand() )
                
                # Update the terms that are to be substracted at order io+1
                retrieve_EqMMSO += self.eps**io * EqMMSO[io]
                
            EqMMS.append(EqMMSO)
            
        self.EqMMS = EqMMS
        
    def apply_MMS(self, rewrite_polar=0):
        r"""
        Apply the MMS. This is operated as follows:

        #. An equivalent system is written in terms of the fast scale :math:`t_0`. This introduces the temporary unknowns :math:`\tilde{x}_{i,j}(t_0)`, which allows the use of ``dsolve()``.

        #. Leading order solutions are introduced.

        #. The leading order solutions are introduced in the equations and the secular terms at each order are identified. 
           Cancelling those secular terms is a condition for bounded solutions and leads to a system of equations governing the slow evolution of the leading order solutions.
           After cancelling the secular terms the equations are solved to express the higher order solutions in terms of the leading order ones.

        #. The phase coordinates are changed to cancel the slow time in the secular terms. This will be used afterwards to obtain an autonomous system.

        #. The secular conditions are split into real and imaginary parts, polar coordinates are introduced and the autonomous phases are introduced,
           resulting in a system of evolution equations. This is the key result of the MMS.

        #. The leading and higher order solutions are rewritten in terms of polar coordinates.

        Parameters
        ----------
        rewrite_polar : str, int or list of int, optional
            The orders at which the solutions will be rewritten in polar form.
            See :func:`sol_xMMS_polar`.
        """
        
        # Write a temporary equivalent system depending only on t0
        self.system_t0()
        
        # Compute the solutions at order 0
        self.sol_order_0()

        # Analyse the secular terms
        self.secular_analysis()

        # Change the phase coordinates for autonomous purposes
        self.autonomous_phases()

        # Derive the evolution equations
        self.evolution_equations()
        
        # Write the x solutions in terms of polar coordinates
        self.sol_xMMS_polar(rewrite_polar=rewrite_polar)


    def system_t0(self):
        r"""
        Rewrite the equations in terms of new coordinates :math:`\tilde{x}_{i,j}(t_0)`, with :math:`i,j` denoting the oscillator number and :math:`\epsilon` order, respectively. 
        This is a trick to use ``dsolve()``, which only accepts functions of 1 variable. 
        """
        
        xMMS_t0  = [] # t0-dependent variables xij(t0). Higher time scales dependency is ignored.
        EqMMS_t0 = [] # Equations at each order with only t0 as an explicit variable. Leads to a harmonic oscillator at each order with a t0-periodic forcing coming from lower order solutions.
        
        for ix in range(self.ndof):
            xMMS_t0 .append([ Function(r'\tilde{x_'+'{{{},{}}}'.format(ix,io)+'}', real=True)(self.tS[0]) for io in range(0, 1+self.Ne) ]) # XXX : should it be real or complex?
            EqMMS_t0.append([ self.EqMMS[ix][0].subs(self.xMMS[ix][0], xMMS_t0[ix][0]).doit() ])
            
        self.EqMMS_t0 = EqMMS_t0
        self.xMMS_t0  = xMMS_t0
        
        
    def sol_order_0(self):
        r"""
        Compute the leading-order solutions for each oscillator. 

        For oscillator :math:`i`, the homogeneous solution takes the general form
        
        .. math::
            x_{i,0}^{(\textrm{h})}(\boldsymbol{t}) = A_i(t_1, t_2, ...) e^{\textrm{j} \omega_{i,0} t_0} + cc,
        
        where :math:`A_i` is a complex amplitudes to be determined.
        
        If the oscillator is subject to hard forcing (i.e. forcing appears at leading order), then the particular solution
        
        .. math::
            x_{i,0}^{(\textrm{p})}(\boldsymbol{t}) = B_i e^{\textrm{j} \omega t} + cc = B_i e^{\textrm{j} (\omega_{\textrm{MMS}} t_0 + \sigma t_1)} + cc,
        
        is also taken into account. :math:`B_i` is a time-independent function of the forcing parameters.
        """
        
        # Information
        print('Definition of leading order multiple scales solutions')
        
        # Initialisation
        xMMS0    = [] # leading order solutions
        sub_xMMS = [] # Substitutions from xij to its solution
        sub_B    = [] # Substitutions from the particular solution amplitude Bi to its expression
        
        # Compute the solutions
        for ix in range(self.ndof):
            
            # Homogeneous leading order solution 
            xMMS0_h_ix = (            self.coord.A[ix]*exp(I*self.omegas_O0[ix]*self.tS[0]) 
                          + conjugate(self.coord.A[ix]*exp(I*self.omegas_O0[ix]*self.tS[0])) )
            
            # Particular leading order solution - if the equation is not homogeneous (due to hard forcing)
            if not self.EqMMS[ix][0] == self.xMMS[ix][0].diff(self.tS[0],2) + (self.omegas_O0[ix])**2 * self.xMMS[ix][0]:
                hint="nth_linear_constant_coeff_undetermined_coefficients"
                
                # General solution, containing both homogeneous and particular solutions
                xMMS0_sol_general = ( dsolve(self.EqMMS_t0[ix][0], self.xMMS_t0[ix][0], hint=hint) ).rhs
                
                # Cancel the homogeneous solutions
                C      = list(xMMS0_sol_general.atoms(Symbol).difference(self.EqMMS[ix][0].atoms(Symbol)))
                sub_IC = [(Ci, 0) for Ci in C]
                xMMS0_p_ix = xMMS0_sol_general.subs(sub_IC).doit()
                
                # Get the real amplitude of the particular solution
                exp_keys = list(xMMS0_p_ix.atoms(exp))
                if exp_keys:
                    sub_B.append( (self.coord.B[ix], xMMS0_p_ix.coeff(exp_keys[0])) )
                else:
                    print("Static hard forcing is currently not handled")
                
                # Rewrite the particular solution in terms of B for the sake of readability and computational efficiency
                xMMS0_p_ix = (          self.coord.B[ix]*exp(I*self.omega*self.t).subs([self.sub.sub_omega]).expand().subs(self.sub.sub_t).expand() + 
                              conjugate(self.coord.B[ix]*exp(I*self.omega*self.t).subs([self.sub.sub_omega]).expand().subs(self.sub.sub_t).expand()))
                    
            else:
                xMMS0_p_ix = sympify(0)
                
            # Total leading order solution
            xMMS0.append( xMMS0_h_ix + xMMS0_p_ix ) 
            sub_xMMS.append( ( self.xMMS[ix][0], xMMS0[ix] ) )
        
        # Store the solutions
        self.sol.xMMS = [[xMMS0_dof] for xMMS0_dof in xMMS0]
        self.sub.sub_xMMS = sub_xMMS
        self.sub.sub_B    = sub_B
        
    def secular_analysis(self):
        r"""
        Identification of the secular terms in the equations. 
        This allows to:

        1. Compute the slow-times evolution of the complex amplitudes :math:`A_i` that cancel the secular terms. 
           These slow time evolutions are written :math:`D_j A_i`, with the notation :math:`D_j(\bullet) = \partial (\bullet)/\partial t_j`.
        
        2. Write the equations with the secular terms cancelled (nonsecular equations).
        
        3. Compute the higher order solutions :math:`x_{i,j},\; j>0`, in terms of the :math:`A_k` from the nonsecular equations.
        """
        
        # Information
        print("Secular analysis")
        
        # Initialisations - secular analysis
        DA_sol     = [] # Solutions Di(Aj) cancelling the secular terms for each oscillator j, in terms of Aj 
        sub_DA_sol = [] # Substitutions from DiAj to its solution 
        sec        = [] # The ith secular term in the equations of the jth oscillator is written only in terms of Di(Aj) and Aj (i.e. Dk(Aj) with k<i are substituted for their solution)
        
        for ix in range(self.ndof):
            DA_sol    .append([ 0 ]) # dAi/dt0 = 0 
            sub_DA_sol.append([ (self.coord.A[ix].diff(self.tS[0]), 0)] )
            sec       .append([ 0])
        
        E = symbols('E') # Symbol to substitue exponentials and use collect() in the following
        
        # Computation of the secular terms, DA solutions, equations with the secular terms cancelled and x solutions in terms of A
        for io in range(1,self.Ne+1):
            
            print('   Analysing the secular terms at order {}'.format(io))
            
            # Substitutions from x(t0, t1, ...) to x(t0) at order io to use sy.dsolve() in the following
            sub_xMMS_t0 = [ (self.xMMS[ix][io], self.xMMS_t0[ix][io]) for ix in range(self.ndof) ]
            
            # Substitute the solutions at previous orders in the MMS equations and make it t0-dependent. Contains the secular terms.
            EqMMS_t0_sec = [ self.EqMMS[ix][io].subs(self.sub.sub_xMMS).subs(sub_xMMS_t0).doit() for ix in range(self.ndof) ] 
            
            # Find the secular terms and deduce the D(A) that cancel them
            dicE = [] 
            for ix in range(self.ndof):
                
                # Define the exponential corresponding to secular terms
                sub_exp = [(exp(I*self.omegas_O0[ix]*self.tS[0]), E)] # Substitute exp(I*omegas_O0*t0) by E to use sy.collect() in the following
                
                # Substitute the low order DA to get rid of all A derivatives except the current one
                EqMMS_t0_sec[ix] = sfun.sub_deep(EqMMS_t0_sec[ix], sub_DA_sol[ix])
                
                # Identify the secular term
                dicE_ix = EqMMS_t0_sec[ix].expand().subs(sub_exp).doit().expand().collect(E, evaluate=False)
                if E in dicE_ix.keys():
                    sec_ix  = dicE_ix[E]
                else:
                    sec_ix = sympify(0)
                dicE.append(dicE_ix)
                
                # Solve D(A) such that the secular term is cancelled
                DA_sol[ix].append( solve(sec_ix, self.coord.A[ix].diff(self.tS[io]))[0] ) # Solution for the current D(A) cancelling the secular terms
                sub_DA_sol[ix].append( (self.coord.A[ix].diff(self.tS[io]), DA_sol[ix][io].expand()) )
            
                # Store the current secular term
                sec[ix].append(sec_ix)
                
            # Substitute the expression of the just computed DA in EqMMS_t0_sec to obtain nonsecular equations governing xMMS_t0 at the current order
            for ix in range(self.ndof):
                self.EqMMS_t0[ix].append(EqMMS_t0_sec[ix].subs(sub_DA_sol[ix]).doit().simplify())
            
            # Compute the x solution at order io in terms of the amplitudes A
            print('   Computing the higher order solutions at order {}'.format(io))
            for ix in range(self.ndof): 
                self.sol_higher_order(self.EqMMS_t0, self.xMMS_t0, io, ix)
            
        # Store the solutions
        self.sol.sec  = sec      # Secular terms
        self.sol.DA   = DA_sol   # Solutions that cancel the secular terms
    
    def sol_higher_order(self, EqMMS_t0, xMMS_t0, io, ix):
        r"""
        Compute higher order solutions :math:`x_{i,j}, j>0`.

        Parameters
        ----------
        EqMMS_t0 : list of list of expr
            The MMS equations at each order and for each oscillator written with :math:`t_0` as the only independent variable. 
        xMMS_t0 : list of list of Function
            Oscillators' solutions at each order, :math:`\tilde{x}_{i,j}(t_0)`.
        io : int
            The current order of :math:`\epsilon`.
        ix : int
            The current oscillator number.
        """
        
        # Hint for dsolve()
        if not EqMMS_t0[ix][io] == xMMS_t0[ix][io].diff(self.tS[0],2) + (self.omegas_O0[ix])**2 * xMMS_t0[ix][io]:
            hint="nth_linear_constant_coeff_undetermined_coefficients"
        else:
            hint="default"
        
        # General solution, containing both homogeneous and particular solutions
        xMMS_sol_general = ( dsolve(EqMMS_t0[ix][io], xMMS_t0[ix][io], hint=hint) ).rhs
        
        # Cancel the homogeneous solutions
        C      = list(xMMS_sol_general.atoms(Symbol).difference(EqMMS_t0[ix][-1].atoms(Symbol)))
        sub_IC = [(Ci, 0) for Ci in C]
        
        # Append the solution for dof ix at order io
        self.sol.xMMS[ix].append(xMMS_sol_general.subs(sub_IC).doit())
        
        # Update the list of substitutions from the x to their expression
        self.sub.sub_xMMS.append( (self.xMMS[ix][io], self.sol.xMMS[ix][io]) )
        
        
    def polar_coordinates(self):
        r"""
        Introduce polar coordinates such that, for oscillator :math:`i`, 
        
        .. math::
            A_i = \frac{1}{2} a_i e^{\textrm{j} \phi_i}.
        
        """
        
        self.coord.a   = [ Function(r'a_{}'.format(ix)   , real=True, positive=True)(*self.tS[1:]) for ix in range(self.ndof) ]
        self.coord.phi = [ Function(r'\phi_{}'.format(ix), real=True)               (*self.tS[1:]) for ix in range(self.ndof) ]
        self.sub.sub_A = [ ( self.coord.A[ix], Rational(1/2)*self.coord.a[ix]*exp(I*self.coord.phi[ix]) ) for ix in range(self.ndof)]
        
    
    def autonomous_phases(self):
        r"""
        Introduce new phase coordinates :math:`\beta_i` to transform nonautonomous equations into autonomous ones. 
        The :math:`\beta_i` are defined as
        
        .. math::
            \beta_i = - \frac{r_{\textrm{MMS}}}{r_i} \phi_i + \sigma t_1,
        
        where we recall that
        :math:`\omega = r_{\textrm{MMS}} \omega_{\textrm{ref}} + \epsilon \sigma` and
        :math:`\omega_{i,0} = r_i \omega_{\textrm{ref}}`. 
        """
        
        self.coord.beta   = [ Function(r'\beta_{}'.format(ix), real=True)(*self.tS[1:])                                            for ix in range(self.ndof) ]
        def_beta          = [ - Rational(self.ratio_omegaMMS, self.ratio_omega_osc[ix])*self.coord.phi[ix] + self.sigma*self.tS[1] for ix in range(self.ndof) ]
        def_phi           = [ solve(def_beta[ix]-self.coord.beta[ix], self.coord.phi[ix])[0]                                       for ix in range(self.ndof) ]
        self.sub.sub_phi  = [ (self.coord.phi[ix], def_phi[ix])                                                                    for ix in range(self.ndof) ]
        self.sub.sub_beta = [ (self.coord.beta[ix], def_beta[ix])                                                                  for ix in range(self.ndof) ]

    def evolution_equations(self):
        r"""
        Derive the evolution equations from the secular conditions. For oscillator :math:`i`, these are defined as
        
        .. math::
            \begin{cases}
            \frac{\textrm{d} a_i}{\textrm{d} t}         & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \frac{\textrm{d} \beta_i}{\textrm{d} t} & = f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}),
            \end{cases}
        
        where :math:`\boldsymbol{a}` and :math:`\boldsymbol{\beta}` are vectors containing the polar amplitudes and phases.

        The aim here is to compute all the :math:`f_{a_i}` and :math:`f_{\beta_i}`.
        This is done by:
        
        - Introducing polar coordinates in the secular terms
        
        - Splitting the real and imaginary parts
        
        - Using the autonomous phase coordinates
        
        - Collecting the terms governing the slow amplitude and phase dynamics.
        """
        
        # Information
        print('Computing the evolution equations')

        # Initialisation
        sec_re    = [[] for dummy in range(self.ndof)] # Real part of the secular terms
        sec_im    = [[] for dummy in range(self.ndof)] # Imaginary part of the secular terms
        sub_re_im = [] # To overcome a sympy limitation: derivatives of real functions w.r.t. a real variable are not recognised as real
        
        faO    = [[] for dummy in range(self.ndof)] # Defined as      Di(a) = faO[i](a,beta)
        fbetaO = [[] for dummy in range(self.ndof)] # Defined as a*Di(beta) = fbetaO[i](a,beta)
        fa     = [0 for dummy in range(self.ndof)]  # Defined as      da/dt = fa(a,beta)
        fbeta  = [0 for dummy in range(self.ndof)]  # Defined as a*dbeta/dt = fbeta(a,beta)
        
        for io in range(0, self.Ne+1):
            
            # print("    Evolution equations at order {}".format(io))
            
            for ix in range(self.ndof):

                # Order 0 -> there are no secular terms
                if io==0:
                    sec_re[ix].append(0)
                    sec_im[ix].append(0)
                    faO[ix].append(sympify(0))
                    fbetaO[ix].append(sympify(0))
                
                # Deal with the secular terms at order eps**io
                else:        
                    
                    # State that da/dti and dphi/dti are real
                    sub_re_im.extend( [(re(self.coord.a[ix].diff(self.tS[io]))  , self.coord.a[ix].diff(self.tS[io])  ),
                                       (im(self.coord.a[ix].diff(self.tS[io]))  , 0                                   ),
                                       (re(self.coord.phi[ix].diff(self.tS[io])), self.coord.phi[ix].diff(self.tS[io])),
                                       (im(self.coord.phi[ix].diff(self.tS[io])), 0                                   )] )
                    
                    # Split sec into real and imaginary parts and change phase coordinates from phi to beta
                    sec_re[ix].append( re( (self.sol.sec[ix][io]*4*exp(-I*self.coord.phi[ix]))
                                          .expand().subs(self.sub.sub_A).doit().expand() )
                                      .simplify().subs(sub_re_im).subs(self.sub.sub_phi).doit().expand() )
                    
                    sec_im[ix].append( im( (self.sol.sec[ix][io]*4*exp(-I*self.coord.phi[ix]))
                                          .expand().subs(self.sub.sub_A).doit().expand() )
                                      .simplify().subs(sub_re_im).subs(self.sub.sub_phi).doit().expand() )
                    
                    # Derive the evolution equations at each order
                    faO[ix]   .append( solve(sec_im[ix][io], self.coord.a[ix]   .diff(self.tS[io]))[0]                  )
                    fbetaO[ix].append( solve(sec_re[ix][io], self.coord.beta[ix].diff(self.tS[io]))[0]*self.coord.a[ix] )
                    
                    # Global evolution equations
                    fa[ix]    += self.eps**io * faO[ix][io]
                    fbeta[ix] += self.eps**io * fbetaO[ix][io]
        
        # Store the results
        self.sol.faO    = faO
        self.sol.fbetaO = fbetaO
        self.sol.fa     = fa
        self.sol.fbeta  = fbeta
                    

    def sol_xMMS_polar(self, rewrite_polar=0):
        r"""
        Write the solutions using the polar coordinates.

        Parameters
        ----------
        rewrite_polar : str or int or list of int or optional
            The orders at which the solutions will be rewritten in polar form.
            If ``"all"``, then all solution orders will be rewritten.
            If `int`, then only a single order will be rewritten.
            If `list` of `int`, then the listed orders will be rewritten.
            Default is 0, so only the leading order solution will be rewritten.
        """
        
        # Information
        print("Rewritting the solutions in polar form")
        
        # Orders to rewrite
        if rewrite_polar=="all":
            rewrite_polar = range(self.Ne+1)
        elif not isinstance(rewrite_polar, list):
            rewrite_polar = [rewrite_polar]
        if max(rewrite_polar)>self.Ne:
            print("Trying to rewrite a solution order that exceeds the maximum order computed.")
            return

        # Prepare substitutions
        sub_t_back = [ (item[1], item[0]) for item in self.sub.sub_t]
        sub_sigma  = [ (self.eps*self.sigma, self.omega-self.omegaMMS)]
        
        # Prepare the collection of sin and cos terms
        harmonics = self._find_harmonics()
        collect_omega = [sin(h*self.omega*self.t) for h in harmonics] + [cos(h*self.omega*self.t) for h in harmonics]
        
        # Rewrite the solutions
        xMMS_polar = []
        x          = [0 for dummy in range(self.ndof)]
        for ix in range(self.ndof):
            xMMS_polar.append([])
            for io in rewrite_polar:
                xMMS_polar[ix].append( TR10(TR8((self.sol.xMMS[ix][io]
                                        .subs(self.sub.sub_A).doit().expand()
                                        .subs(self.sub.sub_phi).doit()) 
                                      .rewrite(cos).simplify()) 
                                      .subs(sub_t_back).subs(sub_sigma).simplify()) 
                                      .expand()
                                      .collect(collect_omega)
                                      )
            if rewrite_polar == range(self.Ne+1): # Construct the full response if relevant
                x[ix] = sum([self.eps**(io+self.eps_pow_0) * xMMS_polar[ix][io] for io in range(self.Ne+1)]).simplify()
            else:
                x[ix] = "all solution orders were not rewritten in polar form"
        # Store
        self.sol.xMMS_polar = xMMS_polar
        self.sol.x          = x

    
    def _find_harmonics(self):
        list_xMMS = list(itertools.chain.from_iterable(self.sol.xMMS))
        harmonics = []
        for xMMS_ix in list_xMMS:
            exponents = [exp_term.args[0].subs(self.tS[1],0) for exp_term in xMMS_ix.atoms(exp)]
            for exponent in exponents:
                if self.tS[0] in exponent.atoms() and im(exponent)>0:
                    harmonics.append( exponent/(I*self.omega_ref*self.tS[0]) / self.ratio_omegaMMS )
            
        harmonics = list(dict.fromkeys(harmonics))
        harmonics.sort()
        return harmonics
    
class Substitutions_MMS:
    """
    Substitutions used in the MMS.
    """
    
    def __init__(self, sub_t, sub_xMMS_t, sub_x, sub_scaling, sub_omega, sub_sigma): 
        self.sub_t            = sub_t
        self.sub_xMMS_t       = sub_xMMS_t
        self.sub_x            = sub_x
        self.sub_scaling      = sub_scaling[0]
        self.sub_scaling_back = sub_scaling[1]
        self.sub_omega        = sub_omega
        self.sub_sigma        = sub_sigma
        
class Forcing_MMS:
    r"""
    Define the forcing on the system as
    
    - A forcing amplitude `F`
    
    - A scaling order `f_order` for the forcing
    
    - Forcing coefficients `f_coeff`
    
    - Forcing terms (direct or parametric) `forcing_term`
    """
    
    def __init__(self, F, f_order, f_coeff, forcing_term):
        self.F            = F
        self.f_order      = f_order
        self.f_coeff      = f_coeff
        self.forcing_term = forcing_term
        
class Coord_MMS:
    """
    The coordinates used in the MMS.
    """      
    
    def __init__(self, mms):
    
        self.A = [] # Complex amplitudes of the homogeneous leading order solutions
        self.B = [] # Real amplitudes of the particular leading order solutions (nonzero only if the forcing is hard)
        
        for ix in range(mms.ndof):
            self.A.append( Function(r'A_{}'.format(ix), complex=True)(*mms.tS[1:]) ) 
            
            if mms.forcing.f_order == 0: # Condition for hard forcing
                self.B.append( symbols(r'B_{}'.format(ix), real=True) ) 

class Sol_MMS:
    """
    Solutions obtained when applying the MMS.
    """                
    
    def __init__(self):
        pass


class Steady_state:
    """
    Steady-state analysis of previous computations. 
    Steady-state means that amplitudes and phases are time-independent. 
    """
    
    def __init__(self, mms):
        """
        Evaluate the MMS quantities at steady-state. 
        """
        
        # Information
        print('Initialisation of the steady state analysis')
        
        # Small parameter
        self.eps = mms.eps
        
        # MMS frequencies of interest
        self.omega_ref      = mms.omega_ref
        self.ratio_omegaMMS = mms.ratio_omegaMMS
        self.sigma          = mms.sigma
        self.omegaMMS       = mms.omegaMMS
        
        # Oscillators' internal resonances relations
        self.ratio_omega_osc = mms.ratio_omega_osc
        
        # Number of dof
        self.ndof = mms.ndof
        
        # Substitutions (initialisation)
        self.sub = Substitutions_SS(mms)
    
        # Forcing
        self.forcing = Forcing_SS(mms)
        
        # Coordinates
        self.coord = Coord_SS()
        self.SS_coord(mms)
        
        # Solutions
        self.sol = Sol_SS(self, mms)
        
        # Stability
        self.stab = Stab_SS()
        
        # Evolution equations at steady state
        self.SS_evolution_equations(mms)
        
    
    def SS_coord(self, mms):
        """
        Create time-independent amplitudes and phases.
        """
        
        a, beta, sub_SS = [], [], []
        for ix in range(self.ndof):
            a   .append( symbols(r'a_{}'.format(ix),positive=True))
            beta.append( symbols(r'\beta_{}'.format(ix),real=True) )
            sub_SS .extend( [(mms.coord.a[ix] , a[ix]), (mms.coord.beta[ix], beta[ix])] )
    
        self.coord.a    = a
        self.coord.beta = beta
        self.sub.sub_SS = sub_SS
        
    def SS_evolution_equations(self, mms):
        """
        Evaluate the evolution equations at steady state.
        """
        
        fa, fbeta, faO, fbetaO = [], [], [], []
        for ix in range(self.ndof):
            fa    .append( mms.sol.fa[ix]   .subs(self.sub.sub_SS).doit().expand() .collect([cos(self.coord.beta[ix]), sin(self.coord.beta[ix])]) )
            fbeta .append( mms.sol.fbeta[ix].subs(self.sub.sub_SS).doit().expand() .collect([cos(self.coord.beta[ix]), sin(self.coord.beta[ix])]) )
        
            faO    .append( [mms.sol.faO[ix][io]   .subs(self.sub.sub_SS).doit().expand() .collect([cos(self.coord.beta[ix]), sin(self.coord.beta[ix])]) for io in range(mms.Ne+1)] )
            fbetaO .append( [mms.sol.fbetaO[ix][io].subs(self.sub.sub_SS).doit().expand() .collect([cos(self.coord.beta[ix]), sin(self.coord.beta[ix])]) for io in range(mms.Ne+1)] )
        
        self.sol.fa     = fa
        self.sol.fbeta  = fbeta
        self.sol.faO    = faO
        self.sol.fbetaO = fbetaO
        
        # Check if the evolution equations are autonomous
        if 't_1' in srepr(fa) or 't_1' in srepr(fbeta):
            print("The evolution equations do not form an autonomous system")

    def solve_forced(self, solve_dof=None):
        """
        Find the steady state solution for a given oscillator.
        The response of other oscillators is set to 0.
                
        Parameters
        ----------
        solve_dof: None or int, oprtional
            The dof number to solve for. Start from 0. 
            If `None`, no dof is solved for.
            Default is `None`.
        """
        
        # Conditions for not solving the forced response
        if solve_dof==None or self.forcing.F==0:
            return
        
        # Information
        print('Computing the forced response for dof {}'.format(solve_dof))
        
        # Store the dof that is solved for
        self.sol.solve_dof = solve_dof
        
        # Set the other oscillator's amplitudes to zero
        self.substitution_solve_dof(solve_dof)
        
        # Phase response
        self.solve_phase()
        
        # Frequency (detuning) response
        self.solve_sigma()
        
        # Frequency response (in terms of oscillator's amplitude)
        self.solve_a()
        
        # Amplitude (forcing) respose
        self.solve_F()
        
    def substitution_solve_dof(self, solve_dof):
        """
        Set every dof amplitude to 0 except the one to solve for.
        """
        sub_solve = []
        for ix in range(self.ndof):
            if ix != solve_dof:
                sub_solve.append( (self.coord.a[ix], 0) )
                
        self.sub.sub_solve = sub_solve    
        
    def solve_phase(self):
        r"""
        Find implicit solutions for the oscillator's phase :math:`\beta_i`. 
        The solutions returned are
        :math:`\sin(k \beta_i)` and :math:`\cos(k \beta_i)` where :math:`k` is an integer or rational.
        """
        
        # Evaluate the evolution equations for a single oscillator responding
        fa_dof    = self.sol.fa[self.sol.solve_dof]   .expand().subs(self.sub.sub_solve)
        fbeta_dof = self.sol.fbeta[self.sol.solve_dof].expand().subs(self.sub.sub_solve)
        
        # Collect sin and cos terms in the evolution equations
        collect_sin_cos = list(fa_dof.atoms(cos, sin)) + list(fbeta_dof.atoms(cos, sin))
        collect_sin_cos = [item for item in collect_sin_cos if item.has(self.coord.beta[self.sol.solve_dof])]
    
        def sort_key(expr):
            """
            Sorting function. Assign a lower value to sine terms and a higher value to cosine terms
            """
            if expr.func == sin:
                return 0
            elif expr.func == cos:
                return 1
            else:
                return 2  

        collect_sin_cos = sorted(collect_sin_cos, key=sort_key) # sin terms first, cos terms then

        dic_fa    = fa_dof   .collect(collect_sin_cos, evaluate=False)
        dic_fbeta = fbeta_dof.collect(collect_sin_cos, evaluate=False)
    
        # Check the possibility to solve using standard procedure (quadratic sum) -> enforce the presence of 3 keys : {1, sin(phase), cos(phase)}
        if ( (len(list(set(list(dic_fa.keys()) + list(dic_fbeta.keys())))) != 3) or # cos and sin terms both appear in the same expression 
             (collect_sin_cos[0].args != collect_sin_cos[1].args) ): # Too many phases involved
            print('    No implemented analytical solution')
            return
    
        # Compute the expression of sin/cos as a function of the amplitude
        print('   Computing the phase response')
        
        if collect_sin_cos[0] in dic_fa: # sin in fa
            if 1 in dic_fa.keys():
                sin_phase = (dic_fa[1] / (-dic_fa[collect_sin_cos[0]])).simplify()
            else:
                sin_phase = sympify(0)
                
            if 1 in dic_fbeta.keys():
                cos_phase = (dic_fbeta[1] / (-dic_fbeta[collect_sin_cos[1]])).simplify()
            else:
                cos_phase = sympify(0)
                
        elif collect_sin_cos[1] in dic_fa: # cos in fa
            if 1 in dic_fa.keys():
                cos_phase = (dic_fa[1] / (-dic_fa[collect_sin_cos[1]])).simplify()
            else:
                cos_phase = sympify(0)
            if 1 in dic_fbeta.keys():
                sin_phase = (dic_fbeta[1] / (-dic_fbeta[collect_sin_cos[0]])).simplify()
            else:
                sin_phase = sympify(0)
        
        else:
            print("   dof {} is not forced".format(self.sol.solve_dof))
            return
    
        # Store the solutions
        self.sol.sin_phase = (collect_sin_cos[0], sin_phase)
        self.sol.cos_phase = (collect_sin_cos[1], cos_phase)
        self.sub.sub_phase = [self.sol.sin_phase, self.sol.cos_phase]
    
    def solve_sigma(self):
        r"""
        Solve the forced response in terms of the detuning :math:`\sigma`. 
        We recall that :math:`\omega = \omega_{\textrm{MMS}} + \epsilon \sigma`. 
        """
        
        sin_phase = self.sol.sin_phase[1]
        cos_phase = self.sol.cos_phase[1]
        
        Eq_sig = (sin_phase**2).expand() + (cos_phase**2).expand() - 1
    
        print('   Computing the frequency response')
        sol_sigma = sfun.solve_poly2(Eq_sig, self.sigma)
        sol_sigma = [sol_sigma_i.simplify() for sol_sigma_i in sol_sigma]
        
        self.sol.sigma = sol_sigma
    
    def solve_a(self):
        r"""
        Solve the forced response in terms of the oscillator's amplitude.
        For readability, the output actually returned in :math:`a^2`. 
        """
        
        sin_phase = self.sol.sin_phase[1]
        cos_phase = self.sol.cos_phase[1]
        a         = self.coord.a[self.sol.solve_dof]
        
        # Equation on a
        Eq_a = (sin_phase**2).expand() + (cos_phase**2).expand() - 1
        keys = Eq_a.expand().collect(a, evaluate=False)
        min_power = min(list(keys), key=lambda expr: sfun.get_exponent(expr, a))
        Eq_a = (Eq_a/min_power).expand()
        
        # Solve
        if set(Eq_a.collect(a, evaluate=False).keys()) in [set([1, a**4]), set([1, a**2, a**4])]:
            print("   Computing the response wrt the oscillator's amplitude")
            sol_a2 = sfun.solve_poly2(Eq_a, a**2)
            sol_a2 = [sol_a2_i.simplify() for sol_a2_i in sol_a2]
        else:
            print("   Not computing the response wrt the oscillator's amplitude as the equation to solve is not of 2nd degree")
            sol_a2 = None
        self.sol.sol_a2 = sol_a2
    
    def solve_F(self):
        r"""
        Solve the forced response in terms of the forcing amplitude :math:`F`.
        """
        
        sin_phase = self.sol.sin_phase[1]
        cos_phase = self.sol.cos_phase[1]
        F         = self.forcing.F
        
        # Equation on F
        Eq_F   = (((sin_phase*F).simplify()**2).expand() 
               + ( (cos_phase*F).simplify()**2).expand() 
               -    self.forcing.F**2).subs(self.sub.sub_B)
        keys = Eq_F.expand().collect(F, evaluate=False)
        min_power = min(list(keys), key=lambda expr: sfun.get_exponent(expr, F))
        Eq_F = (Eq_F/min_power).expand()
        
        # Solve
        if set(Eq_F.collect(F, evaluate=False).keys()) in [set([1, F**2]), set([1, F, F**2])]:
            print('   Computing the response wrt the forcing amplitude')
            sol_F = abs(sfun.solve_poly2(Eq_F, F)[1])
        else:
            print('   Not computing the response wrt the forcing amplitude as the equation to solve is not of 2nd degree')
            sol_F = None
        self.sol.F = sol_F    
    
    def solve_bbc(self, c=[], solve_dof=None):
        """
        Find the backbone curve of a given oscillator.
        The response of other oscillators is set to 0.
        
        Parameters
        ----------
        c: list, oprtional
            Damping terms. They will be set to 0 to compute the backbone curve.
            Note that these are the scaled damping terms.
            Default is `[]`.
        solve_dof: None or int, oprtional
            The dof number to solve for. Start from 0. If `None`, no dof is solved for.
            Default is `None`.
        """
        
        if solve_dof==None:
            return
        
        # Information
        print('Computing the backbone curve for dof {}'.format(solve_dof))
        
        # Set every dof amplitude to 0 except the one to solve for
        self.substitution_solve_dof(solve_dof)
        
        # Substitutions for the free response
        if not isinstance(c, list): 
            c = [c]
        sub_free = [(self.forcing.F,0), *[(ci, 0) for ci in c]]
        
        # Establish the backbone curve equation
        Eq_bbc        = self.sol.fbeta[solve_dof].subs(self.sub.sub_solve).subs(self.sub.sub_B).subs(sub_free)
        
        # Compute the backbone curve
        self.sol.sigma_bbc = solve(Eq_bbc, self.sigma)[0].simplify()
        self.sol.omega_bbc = self.omegaMMS + self.eps*self.sol.sigma_bbc
        
        self.sub.sub_free = sub_free


    def Jacobian_polar(self):
        r"""
        Compute the Jacobian of the evolution equations systems using polar coordinates (see :func:`polar_coordinates` and :func:`evolution_equations`).
        
        Returns
        -------
        J : Matrix
            Jacobian of the polar system, defined as
            
            .. math::
                \textrm{J} = 
                \begin{bmatrix}
                \textrm{J}_0 \\
                \vdots \\
                \textrm{J}_{N-1} 
                \end{bmatrix}
            
            where :math:`\textrm{J}_i` is the :math:`2 \times N` matrix associated to the evolution of dof :math:`i` and defined as
            
            .. math::
                \textrm{J}_i = 
                \begin{bmatrix}
                \frac{\partial f_{a_i}}      {a_0} & \frac{\partial f_{a_i}}      {\beta_0} & \cdots & \frac{\partial f_{a_i}}      {a_{N-1}} & \frac{\partial f_{a_i}}      {\beta_{N-1}} \\
                \frac{\partial f_{\beta_i}^*}{a_0} & \frac{\partial f_{\beta_i}^*}{\beta_0} & \cdots & \frac{\partial f_{\beta_i}^*}{a_{N-1}} & \frac{\partial f_{\beta_i}^*}{\beta_{N-1}}
                \end{bmatrix}.
            
            Functions :math:`f_{\beta_i}^*` are defined as
            
            .. math::
                f_{\beta_i}^* = \frac{f_{\beta_i}}{a_i}.
            
        """
        
        J = zeros(2*self.ndof,2*self.ndof)
        
        for ii, (fai, fbetai, ai) in enumerate(zip(self.sol.fa, self.sol.fbeta, self.coord.a)):
            for jj, (aj, betaj) in enumerate(zip(self.coord.a, self.coord.beta)):
                J[2*ii,2*jj]     = fai.diff(aj)
                J[2*ii,2*jj+1]   = fai.diff(betaj)
                J[2*ii+1,2*jj]   = (fbetai/ai).simplify().diff(aj)
                J[2*ii+1,2*jj+1] = (fbetai/ai).simplify().diff(betaj)
        
        return J
        
        
    def cartesian_coordinates(self):
        r"""
        Define cartesian coordinates from the polar ones.
        The leading order solution for oscillator :math:`i` expressed in polar coordinates takes the form

        .. math::
            \begin{split}
            x_{i,0}(t) & = a_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}(\omega t - \beta_i)\right), \\
                    & = a_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\beta_i\right)\cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) 
                    + a_i \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\beta_i\right)\sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right)
            \end{split}

        The polar coordinates are defined as

        .. math::
            \begin{cases}
            p_i = a_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\beta_i\right), \\
            q_i = a_i \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\beta_i\right),
            \end{cases}

        Such that the leading order solution can be written as

        .. math::
            x_{i,0}(t) = p_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) + q_i \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right).
        """

        
        # Define the cartesian coordinates
        p = [symbols(r'p_{}'.format(ix), real=True) for ix in range(0, self.ndof)]
        q = [symbols(r'q_{}'.format(ix), real=True) for ix in range(0, self.ndof)]
        
        # Define relations between the polar and cartesian coordinates
        a, beta = list(map(self.coord.__dict__.get, ["a", "beta"]))
        sub_cart  = []
        sub_polar = []
        for ix in range(self.ndof):
            ratio_omega = Rational(self.ratio_omega_osc[ix], self.ratio_omegaMMS)
            
            sub_cart.append( (a[ix]*cos(beta[ix]*ratio_omega)     , p[ix]) )
            sub_cart.append( (a[ix]*sin(beta[ix]*ratio_omega)     , q[ix]) )
            sub_cart.append( (a[ix]**2*cos(2*beta[ix]*ratio_omega), p[ix]**2 - q[ix]**2) )
            sub_cart.append( (a[ix]**2*sin(2*beta[ix]*ratio_omega), 2*p[ix]*q[ix]) )
            sub_cart.append( (a[ix]**2                            , p[ix]**2 + q[ix]**2) )
            
            sub_polar.append( (p[ix], a[ix]*cos(beta[ix]*ratio_omega)) )
            sub_polar.append( (q[ix], a[ix]*sin(beta[ix]*ratio_omega)) )
    
        # Store the results
        self.coord.p = p
        self.coord.q = q
        self.sub.sub_cart  = sub_cart
        self.sub.sub_polar = sub_polar
        
    
    def evolution_equations_cartesian(self):
        r"""
        Write the evolution equations using the cartesian coordinates. 
        For oscillator :math:`i`, this results in
        
        .. math::
            \begin{cases}
            \frac{\textrm{d} p_i}{\textrm{d} t} & = f_{p_i}(\boldsymbol{p}, \boldsymbol{q}), \\
            \frac{\textrm{d} q_i}{\textrm{d} t} & = f_{q_i}(\boldsymbol{p}, \boldsymbol{q}),
            \end{cases}
        
        where :math:`\boldsymbol{p}` and :math:`\boldsymbol{q}` are vectors containing the cartesian coordinates.
        """
        
        # Compute the functions fp(p,q) and fq(p,q)
        fp = []
        fq = []
        
        a, beta   = list(map(self.coord.__dict__.get, ["a", "beta"]))
        fa, fbeta = list(map(self.sol.__dict__.get, ["fa", "fbeta"]))
        
        for ix in range(self.ndof):
            ratio_omega = Rational(self.ratio_omega_osc[ix], self.ratio_omegaMMS)
            
            fp.append( TR10( ( fa[ix]*cos(beta[ix]*ratio_omega) - 
                               ratio_omega*fbeta[ix]*sin(beta[ix]*ratio_omega)
                              ).expand().simplify()
                            ).expand().subs(self.sub.sub_cart) )
            
            fq.append( TR10( ( fa[ix]*sin(beta[ix]*ratio_omega) + 
                               ratio_omega*fbeta[ix]*cos(beta[ix]*ratio_omega)
                              ).expand().simplify()
                            ).expand().subs(self.sub.sub_cart) )
        
        # Check if the a and beta have all been substituted
        substitution_OK = self._check_cartesian_substitutions(a, beta, fp, fq)
        
        # Try additional substitutions if the change of coordinates is incomplete
        if not substitution_OK:
            fp, fq = self._additional_cartesian_substitutions(fp, fq)
            
            # Check if the a and beta have all been substituted
            substitution_OK = self._check_cartesian_substitutions(a, beta, fp, fq)
            
        
        if not substitution_OK:
            print("   The substitution from polar to cartesian coordinates is incomplete")
        
        # Store the evolution equations
        self.sol.fp = fp
        self.sol.fq = fq
        
    def _check_cartesian_substitutions(self, a, beta, fp, fq):
        r"""
        Check if substitutions from polar to cartesian coordinates are complete.
        
        Parameters
        ----------
        a: list of Symbol
            Amplitudes of the leading order solutions.
        beta: list of Symbol
            Phases of the leading order solutions.
        fp: list of expr
            Evolution functions for the cartesian coordinates :math:`p_i`.
        fq: list of expr
            Evolution functions for the cartesian coordinates :math:`q_i`.

        Returns
        -------
        substitution_OK : bool
            `True` if substitutions are complete.
            `False` otherwise.
        """
        polar_coordinates = a + beta
        substitution_OK   = True
        count             = 0
        while substitution_OK and count<self.ndof:
            for ix in range(self.ndof):
                symbols_fpq = list(fp[ix].atoms(Symbol)) + list(fq[ix].atoms(Symbol))

                for polar_coordinate in polar_coordinates:
                    if polar_coordinate in symbols_fpq:
                        substitution_OK = False
                        
                count += 1
        
        return substitution_OK
        
    def _additional_cartesian_substitutions(self, fp, fq):
        r"""
        Reformulate the already-existing substitutions from polar to cartesian to try and substitute leftover polar terms.
        
        Parameters
        ----------
        fp: list of expr
            Evolution functions for the cartesian coordinates :math:`p_i`.
            There are polar coordinates remaining.
        fq: list of expr
            Evolution functions for the cartesian coordinates :math:`q_i`.
            There are polar coordinates remaining.

        Returns
        -------
        fp: list of expr
            Evolution functions for the cartesian coordinates :math:`p_i`.
            Additional substitutions were performed to get rid of polar coordinates.
        fq: list of expr
            Evolution functions for the cartesian coordinates :math:`q_i`.
            Additional substitutions were performed to get rid of polar coordinates.
        """
        
        sub_cart_add = [] # Additional substitutions required
        
        for f in fp+fq:
            # Terms that were not properly substituted
            left_polar_terms = list(f.atoms(sin, cos))
            
            # Check if unsubstituted terms were identified
            if left_polar_terms: 
                
                # Loop over the unsubstituted terms
                for term in left_polar_terms:
                    
                    # Loop over already-defined substitutions and look for the current unsubstituted term
                    for sub_cart_item in self.sub.sub_cart:
                        dic = sub_cart_item[0].collect(term, evaluate=False)
                        
                        # Identify possible additional substitutions
                        if term in dic.keys():
                            sub_cart_add.append( (term, solve(sub_cart_item[0]-sub_cart_item[1], term)[0].subs(self.sub.sub_cart)) )
        
        # Apply these new substitutions                                
        for ix in range(self.ndof):
            fp[ix] = fp[ix].subs(sub_cart_add)
            fq[ix] = fq[ix].subs(sub_cart_add)
            
        return fp, fq
    
    
    def Jacobian_cartesian(self):
        r"""
        Compute the Jacobian of the evolution equations expressed in cartesian coordinates (see :func:`cartesian_coordinates` and :func:`evolution_equations_cartesian`).
        
        Returns
        -------
        J : Matrix
            Jacobian of the cartesian system, defined as
            
            .. math::
                \textrm{J} = 
                \begin{bmatrix}
                \textrm{J}_0 \\
                \vdots \\
                \textrm{J}_{N-1} 
                \end{bmatrix}
            
            where :math:`\textrm{J}_i` is the :math:`2 \times N` matrix associated to the evolution of dof :math:`i` and defined as
            
            .. math::
                \textrm{J}_i = 
                \begin{bmatrix}
                \frac{\partial f_{p_i}}{p_0} & \frac{\partial f_{p_i}}{q_0} & \cdots & \frac{\partial f_{p_i}}{p_{N-1}} & \frac{\partial f_{p_i}}{q_{N-1}} \\
                \frac{\partial f_{q_i}}{p_0} & \frac{\partial f_{q_i}}{q_0} & \cdots & \frac{\partial f_{q_i}}{p_{N-1}} & \frac{\partial f_{q_i}}{q_{N-1}} 
                \end{bmatrix}.
            
        """
        
        J = zeros(2*self.ndof,2*self.ndof)
        
        for ii, (fpi, fqi) in enumerate(zip(self.sol.fp, self.sol.fq)):
            for jj, (pj, qj) in enumerate(zip(self.coord.p, self.coord.q)):
                J[2*ii,2*jj]     = fpi.diff(pj)
                J[2*ii,2*jj+1]   = fpi.diff(qj)
                J[2*ii+1,2*jj]   = fqi.diff(pj)
                J[2*ii+1,2*jj+1] = fqi.diff(qj)
        
        return J

    def eval_sol_stability(self, coord="cartesian", rewrite_polar=False, eigenvalues=False, bifurcation_curves=False, analyse_blocks=False, kwargs_bif=dict()):
        r"""
        Evaluate the stability of a solution. 
        
        Parameters
        ----------
        coord: str, optional
            Either ``"cartesian"`` or ``"polar"``. 
            Specifies the coordinates to use for the stability analysis.
            ``"cartesian"`` is recommended as it prevents divisions by 0, which occur when at least one of the dof has a null ampliutude.
            Default is ``"cartesian"``.
        rewrite_polar: bool, optional
            Rewrite the Jacobian's determinant and trace in polar coordinates (if computed using cartesian ones).
            This is time consuming and the current back substitutions from cartesian to polar coordinates are not always sufficient.
            Default is `False`.
        eigenvalues: bool, optional
            Compute the eigenvalues of the Jacobian.
            Default is `False`.
        bifurcation_curves: bool, optional
            Compute the bifurcation curves.
            Default is `False`.
        analyse_blocks: bool, optional
            Analyse the diagonal blocks of the Jacobian rather than the Jacobian itself. This is relevant if the Jacobian is block-diagonal.
        kwargs_bif: dict, optional
            Passed to :func:`bifurcation_curves`
            Default is `dict()`.
        """
        
        # Check if a solution has been computed
        if not "sigma" in self.sol.__dict__.keys():
            print("There is no solution to evaluate the stability of.")
            return
        
        # Information
        print("Evaluating the stability of the solution of dof {}".format(self.sol.solve_dof))
        
        # Introduce the cartesian coordinates and evolution equations
        if coord == "cartesian":
            print("   Rewritting the system in cartesian coordinates")

            self.stab.analysis_coord = "cartesian"
            self.cartesian_coordinates()
            self.evolution_equations_cartesian()

        else:
            self.stab.analysis_coord = "polar"
        
        # Compute the Jacobian
        print("   Computing the Jacobian")
        if coord=="cartesian":
            J = self.Jacobian_cartesian()
        else:
            J = self.Jacobian_polar()
        
        # Set every dof amplitude to 0 except the one solved for
        if coord=="cartesian":
            for ix in range(self.ndof):
                if ix != self.sol.solve_dof:
                    self.sub.sub_solve.extend( [(self.coord.p[ix], 0), (self.coord.q[ix], 0)] )
         
        # Use the steady-state solutions to perform substitutions
        self.sub.sub_solve.extend( [(self.forcing.F*self.sol.sin_phase[0], self.forcing.F*self.sol.sin_phase[1]),
                                    (self.forcing.F*self.sol.cos_phase[0], self.forcing.F*self.sol.cos_phase[1])] )
        
        if coord=="cartesian": 
            if self.forcing.F in self.sol.fp[self.sol.solve_dof].atoms(Symbol): 
                self.sub.sub_solve.append( (self.forcing.F, solve(self.sol.fp[self.sol.solve_dof].subs(self.sub.sub_solve), self.forcing.F)[0]) )
            else:
                self.sub.sub_solve.append( (self.forcing.F, solve(self.sol.fq[self.sol.solve_dof].subs(self.sub.sub_solve), self.forcing.F)[0]) )
        
        # Evaluate the Jacobian on the solution
        Jsol = simplify(J.subs(self.sub.sub_solve)) 
        
        # Analyse the Jacobian
        tr_Jsol  = trace(Jsol).simplify()
        det_Jsol = det(Jsol).simplify()
        
        # Rewrite the results in polar form if cartesian coordinates were used (time consuming)
        if coord=="cartesian": 
            # Save cartesian results
            self.stab.Jsolc     = Jsol
            self.stab.tr_Jsolc  = tr_Jsol
            self.stab.det_Jsolc = det_Jsol
            
            # Write the results in polar form
            if rewrite_polar:
                print("   Expressing the stability results in polar coordinates")
                Jsol     = cartesian_to_polar(Jsol, self.sub.sub_polar, sub_phase=self.sub.sub_phase)
                tr_Jsol  = cartesian_to_polar(tr_Jsol, self.sub.sub_polar, sub_phase=self.sub.sub_phase)
                det_Jsol = cartesian_to_polar(det_Jsol, self.sub.sub_polar, sub_phase=self.sub.sub_phase)

        # Store results
        self.stab.Jsol     = Jsol
        self.stab.tr_Jsol  = tr_Jsol
        self.stab.det_Jsol = det_Jsol
        
        # Compute eigenvalues and bifurcation curves from the analysis of Jsol
        if not analyse_blocks:
            if eigenvalues:
                self.stab.eigvals = self.eigenvalues(Jsol)
            if bifurcation_curves:
                self.stab.bif_a, self.stab.bif_sigma = self.bifurcation_curves(det_Jsol, tr_Jsol, **kwargs_bif)

        # Analyse the blocks of Jsol
        if analyse_blocks:
            print("   Block analysis")

            if coord == "cartesian":
                Jsol = self.stab.Jsolc

            if sfun.is_block_diagonal(Jsol, 2):
                self.stab.blocks         = []
                self.stab.blocks_det     = []
                self.stab.blocks_tr      = []
                self.stab.blocks_eigvals = []
                self.stab.blocks_bif_a   = []
                self.stab.blocks_bif_sig = []

                for idx in range(0, Jsol.rows, 2):
                    A = Jsol[idx:idx+2, idx:idx+2] 
                    self.stab.blocks.append(A)
                    detA = det(A)
                    trA  = trace(A) 
                    if coord=="cartesian":
                        detA = cartesian_to_polar(detA, self.sub.sub_polar, sub_phase=self.sub.sub_phase).factor()
                        trA  = cartesian_to_polar(trA, self.sub.sub_polar, sub_phase=self.sub.sub_phase).factor()
                    
                    self.stab.blocks_det.append(detA)
                    self.stab.blocks_tr.append(trA)

                    if eigenvalues:
                        eigvalsA = self.eigenvalues(A)
                        if coord=="cartesian":
                            eigvalsA = [cartesian_to_polar(eigval, self.sub.sub_polar, sub_phase=self.sub.sub_phase) for eigval in eigvalsA]
                        self.stab.blocks_eigvals.append(eigvalsA)

                    if bifurcation_curves:
                        bif_aA, bif_sigA = self.bifurcation_curves(detA, trA, **kwargs_bif)
                        self.stab.blocks_bif_a.append(bif_aA)
                        self.stab.blocks_bif_sig.append(bif_sigA)

            else:
                print("Trying to perform a block analysis while the Jacobian is not block-diagonal")

    def eigenvalues(self, J):
        r"""
        Computes the eigenvalues of a matrix :math:`\textrm{J}`.

        Parameters
        ----------
        J: Matrix
            The matrix whose eigenvalues are to be computed.

        Returns
        -------
        eigvals: list
            The eigenvalues of :math:`\textrm{J}`.
        """

        print("   Computing eigenvalues")
        
        lamb        = symbols(r"\lambda")
        eig_problem = J - lamb * eye(*J.shape)
        detEP       = eig_problem.det()
        eigvals     = solve(detEP, lamb)
        
        return eigvals
            
    def bifurcation_curves(self, detJ, trJ, var_a=False, var_sig=True, solver=sfun.solve_poly2):
        r"""
        Compute bifurcation curves, i.e. the curves defining the bifurcation points. 

        Parameters
        ----------
        detJ: Expr
            The determinant of the matrix.
        trJ: Expr
            The trace of the matrix.
        var_a: bool, optional
            Consider the oscillator's amplitude :math:`a` as the variable and find the bifurcation curve as an expression for :math:`a`.
            `detJ` is rarely a quadratic polynomial in :math:`a`, so this can rarely be computed easily.
            Default is `False`.
        var_sig: bool, optional
            Consider the detuning :math:`\sigma` as the variable and find the bifurcation curve as an expression for :math:`\sigma`.
            `detJ` is often a quadratic polynomial in :math:`\sigma`, so this can often be computed.
            Default is `True`.
        solver: function, optional
            The solver to use to compute the bifurcation curves.
            Available are solver called as `solve(expr, x)`, which solve `expr=0` for `x`.
            `sy.solve()` can be used but is sometimes slow.
            Default is :func:`sfun.solve_poly2`.

        Returns
        -------
        bif_a : list
            The bifurcation curves expressed in terms of :math:`a^2`.
        bif_sig : list
            The bifurcation curves expressed in terms of :math:`\sigma`. 
        """
        
        print("   Computing bifurcation curves")

        # Check if a stability analysis was performed
        if not "Jsol" in self.stab.__dict__.keys():
            print("There was no stability analysis performed.")
            return

        # Check if the stability analysis is expressed in polar coordinates
        if "p" in self.coord.__dict__.keys():
            cartesian_coordinates = self.coord.p + self.coord.q
            symbols_det = list(detJ.atoms(Symbol)) 
            for cartesian_coordinate in cartesian_coordinates:
                if cartesian_coordinate in symbols_det:
                    print("Substitutions from cartesian back to polar coordinates were incomplete. \n ",
                          "Try other substitutions manually or compute the Jacobian's determinant using block partitions if possible")

        # Compute the bifurcation curves from the determinant of the Jacobian
        if var_a and sfun.check_solvability(detJ, self.coord.a[self.sol.solve_dof]**2):
            bif_a = solver(detJ, self.coord.a[self.sol.solve_dof]**2)
        else:
            bif_a = []

        if var_sig and sfun.check_solvability(detJ, self.sigma):
            bif_sig = solver(detJ, self.sigma)
        else:
            bif_sig = []
        
        # Add bifurcation curves related to the trace of the Jacobian if it is not a constant
        if self.coord.a[self.sol.solve_dof] in trJ.atoms(Symbol):
            bif_a   += solver(trJ, self.coord.a[self.sol.solve_dof]**2)
        if self.sigma in list(trJ.atoms(Symbol)):
            bif_sig += solver(trJ, self.sigma)
        
        # Return
        return bif_a, bif_sig
    
    @staticmethod
    def plot_FRC(FRC, **kwargs):
        r"""
        Plots the frequency response curves (FRC), both frequency-amplitude and frequency-phase.
        Also includes the stability information if given.

        Parameters
        ----------
        FRC : dict
            Dictionary containing the frequency response curves and the bifurcation curves.
        
        Returns
        -------
        fig1 : Figure
            The amplitude plot :math:`a(\omega)`.
        fig2 : Figure
            The phase plot :math:`\beta(\omega)`.
        """

        # Extract the FRC data
        a         = FRC.get("a", np.full(10, np.nan))
        omega_bbc = FRC.get("omega_bbc", np.full_like(a, np.nan))
        omega     = FRC.get("omega", [np.full_like(a, np.nan)])
        phase     = FRC.get("phase", [np.full_like(a, np.nan)])
        omega_bif = FRC.get("omega_bif", [np.full_like(a, np.nan)])
        phase_bif = FRC.get("phase_bif", [np.full_like(a, np.nan)])
        
        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        amp_name   = kwargs.get("amp_name", "amplitude")
        phase_name = kwargs.get("phase_name", "phase")
        xlim       = kwargs.get("xlim", [coeff*np.min(omega_bbc) for coeff in (0.9, 1.1)])
        if np.isnan(xlim).any():
            xlim = [None, None]
            
        # FRC - amplitude 
        fig1, ax = plt.subplots(**fig_param)
        ax.plot(omega_bbc, a, c="tab:grey", lw=0.7)
        ax.axvline(np.min(omega_bbc), c="k")
        [ax.plot(omegai, a, c="tab:blue") for omegai in omega]
        [ax.plot(omegai, a, c="tab:red", lw=0.7) for omegai in omega_bif]
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"${}$".format(amp_name))
        ax.margins(y=0)
        plt.show(block=False)
        
        # FRC - phase
        fig2, ax = plt.subplots(**fig_param)
        ax.axvline(np.min(omega_bbc), c="k")
        ax.axhline(np.pi/2, c="k", lw=0.7)
        [ax.plot(omegai, phasei, c="tab:blue") for (omegai, phasei) in zip(omega, phase)]
        [ax.plot(omegai, phasei, c="tab:red", lw=0.7) for (omegai, phasei) in zip(omega_bif, phase_bif)]
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"${}$".format(phase_name))
        plt.show(block=False)
        
        # Return
        return fig1, fig2
    
    @staticmethod
    def plot_ARC(ARC, **kwargs):
        r"""
        Plots the amplitude-response curves (ARC), both forcing amplitude-amplitude and forcing amplitude-phase.

        Parameters
        ----------
        ARC : dict
            Dictionary containing the amplitude response curves.

        Returns
        -------
        fig1 : Figure
            The amplitude plot :math:`a(\omega)`.
        fig2 : Figure
            The phase plot :math:`\beta(\omega)`.
        """
    
        # Extract the FRC data and keyword arguments
        a     = ARC.get("a", np.full(10, np.nan))
        F     = ARC.get("F", np.full_like(a, np.nan))
        phase = ARC.get("phase", np.full_like(a, np.nan))
        
        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        amp_name   = kwargs.get("amp_name", "amplitude")
        phase_name = kwargs.get("phase_name", "phase")
        xlim       = kwargs.get("xlim", [0, np.max(F)])

        # ARC - amplitude 
        fig1, ax = plt.subplots(**fig_param)
        ax.plot(F, a, c="tab:blue")
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$F$")
        ax.set_ylabel(r"${}$".format(amp_name))
        ax.margins(x=0, y=0)
        plt.show(block=False)

        # ARC - phase
        fig2, ax = plt.subplots(**fig_param)
        ax.axhline(np.pi/2, c="k", lw=0.7)
        ax.plot(F, phase, c="tab:blue")
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$F$")
        ax.set_ylabel(r"${}$".format(phase_name))
        ax.margins(x=0)
        plt.show(block=False)
        
        # Return
        return fig1, fig2


class Substitutions_SS:
    """
    Substitutions used in the steady state evaluations.
    """
    
    def __init__(self, mms):
        
        self.sub_scaling_back = mms.sub.sub_scaling_back
        self.sub_B            = mms.sub.sub_B
        pass
        
    
class Forcing_SS:
    """
    Define the forcing on the system.
    """
    
    def __init__(self, mms):
        self.F            = mms.forcing.F
        self.f_order      = mms.forcing.f_order
        
class Coord_SS:
    """
    The steady-state coordinates.
    """      
    
    def __init__(self):
        pass
    
class Sol_SS:
    """
    Solutions obtained when evaluating at steady state.
    """                
    
    def __init__(self, ss, mms):
        
        self.x = []
        for ix in range(ss.ndof):
            self.x.append( [xio.subs(mms.sub.sub_t[:-1]+ss.sub.sub_SS) for xio in mms.sol.xMMS_polar[ix]] )
        
class Stab_SS:
    """
    Stability analysis.
    """                
    
    def __init__(self):
        pass  
    
#%% Chain rule functions written for the MMS
def Chain_rule_dfdt(f, tS, eps):
    r"""
    Consider a function :math:`f_t(t)` and its expression in terms of the time scales :math:`f(t_0, t_1, ...)`, 
    where :math:`t_0` is the fast time and :math:`t_1, ...` are the slow times. 
    The Chain Rule is applied to give the expression of :math:`\mathrm{d}f_t/ \mathrm{d}t` in terms of the time scales.
    
    Parameters
    ----------
    f: Function  
        Function :math:`f(t_0,t_1,...)`, i.e. :math:`f_t(t)` expressed in terms of the time scales.
    tS: list 
        Time scales.
    eps: Symbol
        Small parameter :math:`\epsilon`.
    
    Returns
    -------
    dfdt: Function
        :math:`\mathrm{d}f_t/ \mathrm{d}t` expressed in terms of the time scales.
    """
    
    Nt = len(tS)
    dfdt = 0
    for ii in range(Nt):
        dfdt += eps**ii * f.diff(tS[ii])
    
    return dfdt

def Chain_rule_d2fdt2(f, tS, eps):
    r"""
    Consider a function :math:`f_t(t)` and its expression in terms of the time scales :math:`f(t_0, t_1, ...)`, 
    where :math:`t_0` is the fast time and :math:`t_1, ...` are the slow times. 
    The Chain Rule is applied to give the expression of :math:`\mathrm{d}^2f_t/ \mathrm{d}t^2` in terms of the time scales.
    
    Parameters
    ----------
    f: Function  
        Function :math:`f(t_0,t_1,...)`, i.e. :math:`f_t(t)` expressed in terms of the time scales.
    tS: list 
        Time scales.
    eps: Symbol
        Small parameter :math:`\epsilon`.
    
    Returns
    -------
    d2fdt2: Function
        :math:`\mathrm{d}^2f_t/ \mathrm{d}t^2` expressed in terms of the time scales.
    """
    
    Nt = len(tS)
    d2fdt2 = 0
    for jj in range(Nt):
        for ii in range(Nt):
            d2fdt2 += eps**(jj+ii) * f.diff(tS[ii]).diff(tS[jj])
    
    return d2fdt2

def cartesian_to_polar(y, sub_polar, sub_phase=None):
    r"""
    Rewrites an expression or a Matrix `y` from cartesian to polar coordinates.

    Parameters
    ----------
    y : Expr or Matrix
        A sympy expression or matrix written in cartesian coordinates.
    sub_polar: list
        A list of substitutions to perform to go from cartesian to polar coordinates. 
    sub_phase: list, optional
        Additional substitutions to try and get rid of phases, so that only the amplitude remains in the expression.
        Default is `None`.

    Returns
    -------
    yp: Expr or Matrix
        The initial expression or matrix written in polar coordinates.
    """

    if y.is_Matrix:
        yp = simplify(y.subs(sub_polar))
    else:
        if sub_phase == None:
            yp = y.subs(sub_polar).expand().simplify()
        else:
            phase     = sub_phase[0][0].args[0]
            sub_tan   = [(tan(phase/2), sin(phase)/(cos(phase)+1))]
            yp = TR8((TR5(y.subs(sub_polar)).expand())).subs(sub_phase).simplify().subs(sub_tan).subs(sub_phase).expand().simplify()
    return yp

#%% Numpy transforms and plot functions
def sympy_to_numpy(expr_sy, param):
    """
    Transform a sympy expression into a numpy array.

    Parameters
    ----------
    expr_sy : Expr
        A sympy expression.
    param : dict
        A dictionnary whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter
        
        2. The numerical value(s) taken by that parameter

    Returns
    -------
    expr_np : ndarray
        The numerical values taken by the sympy expression evaluated.
    """
    
    args, values = zip(*param.values())
    with warnings.catch_warnings(): 
        warnings.filterwarnings("ignore", message="invalid value encountered in sqrt")
        expr_np = lambdify(args, expr_sy, modules="numpy")(*values)

    return expr_np

def rescale(expr, mms):
    expr_rescaled = expr.subs(*mms.sub.sub_sigma).subs(mms.sub.sub_scaling_back).simplify()
    return expr_rescaled

def numpise_FRC(mms, ss, dyn, param, bbc=True, forced=True, bif=True):
    r"""
    Evaluate the frequency-response and bifurcation curves at given numerical values.
    This transforms the sympy expressions to numpy arrays.

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The MMS results evaluated at steady state.
    dyn : Dynamical_system
        The initial dynamical system.
    param : dict
        A dictionary whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.

        The key of the amplitude vector must be ``"a"``.
        The key of the forcing amplitude must be ``"F"``.
    bbc : bool, optional
        Evaluate the backbone curve. 
        Default is `True`.
    forced : bool, optional
        Evaluate the forced response. 
        Default is `True`.
    bif : bool, optional
        Evaluate the bifurcation curves. 
        Default is `True`.

    Returns
    -------
    FRC : dict
        The frequency-response curves data.
    """

    
    # Information
    print("Converting sympy FRC expressions to numpy")

    # Initialisation
    a     = param.get("a")[1]
    F_val = param.get("F")[1]
    FRC   = {"a": a}

    # Evaluation of the FRC
    if bbc:
        FRC["omega_bbc"] = numpise_omega_bbc(mms, ss, param)
    if forced:
        FRC["omega"] = numpise_omega_FRC(mms, ss, param)
        FRC["phase"] = numpise_phase(mms, ss, dyn, param, FRC["omega"], F_val)
    if bif:
        FRC["omega_bif"] = numpise_omega_bif(mms, ss, param)
        FRC["phase_bif"] = numpise_phase(mms, ss, dyn, param, FRC["omega_bif"], F_val)

    return FRC

def numpise_ARC(mms, ss, dyn, param):
    r"""
    Evaluate the amplitude-response curves at given numerical values. 
    This transforms the sympy expressions to numpy arrays. 

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The MMS results evaluated at steady state.
    dyn : Dynamical_system
        The initial dynamical system.
    param : dict
        A dictionnary whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.

        The key of the amplitude vector must be ``"a"``.
        The key of the angular frequency must be ``"omega"``.

    Returns
    -------
    ARC: dict
        The amplitude-response curves data.
    """
    
    # Information
    print("Converting sympy ARC expressions to numpy")

    # Initialisation
    a         = param.get("a")[1]
    omega_val = param.get("omega")[1]
    ARC       = {"a": a}

    # Evaluation of the FRC
    ARC["F"]     = numpise_F_ARC(mms, ss, param)
    ARC["phase"] = numpise_phase(mms, ss, dyn, param, omega_val, ARC["F"])[0]

    return ARC

def numpise_omega_bbc(mms, ss, param):
    omega_bbc  = sympy_to_numpy(rescale(ss.sol.omega_bbc, mms), param)
    return omega_bbc

def numpise_omega_FRC(mms, ss, param):
    omega = [np.real(sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.sol.sigma]
    return omega

def numpise_omega_bif(mms, ss, param):
    omega_bif = [np.real(sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.stab.bif_sigma]
    return omega_bif

def numpise_phase(mms, ss, dyn, param, omega, F):
    
    if not isinstance(omega,list):
        omega = [omega]
        
    phase = []
    
    for omegai in omega:
        param_phase = param | dict(omega=(mms.omega, omegai), F=(dyn.forcing.F, F))
        sin_phase = sympy_to_numpy( rescale(ss.sol.sin_phase[1], mms), param_phase )
        cos_phase = sympy_to_numpy( rescale(ss.sol.cos_phase[1], mms), param_phase )
        phase.append(np.arctan2(sin_phase, cos_phase))

    return phase

def numpise_F_ARC(mms, ss, param):
    F = sympy_to_numpy(rescale(mms.eps**mms.forcing.f_order * ss.sol.F, mms), param)
    return F