# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the multiple scales system from the dynamical one, and the application of the MMS.
"""

#%% Imports and initialisation
from sympy import (exp, I, re, im, Rational, 
                   symbols, Symbol, Function, Expr, sympify, simplify, 
                   solve, dsolve, cos, sin, tan, sympify, Mod)
from sympy.simplify.fu import TR5, TR8, TR10
import itertools
from typing import Union, TYPE_CHECKING

#%% Classes and functions
def scale_parameters(param, scaling, eps):
    r"""
    Scale parameters with the scaling parameter :math:`\epsilon`.

    Parameters
    ----------
    param : list of sympy.Symbol and/or sympy.Function
        Unscaled parameters.
    scaling : list of int or float
        The scaling for each parameter.
    eps : sympy.Symbol
        Small parameter :math:`\epsilon`.

    Returns
    -------
    param_scaled: list of sympy.Symbol and/or sympy.Function
        Scaled parameters.
    sub_scaling: list of 2 lists of tuple
        Substitutions from scaled to unscaled parameters and vice-versa. 

        - :math:`1^{\text{st}}` list: The substitutions to do to introduce the scaled parameters in an expression.
        
        - :math:`2^{\text{nd}}` list: The substitutions to do to reintroduce the unscaled parameters in a scaled expression.
    
    Notes
    -----
    For a given parameter :math:`p` and a scaling order :math:`\lambda`, the associated scaled parameter :math:`\tilde{p}` is 

    .. math::
        p = \epsilon^{\lambda} \tilde{p} .
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
            
class Substitutions_MMS:
    """
    Substitutions used in the MMS.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        sub_A: list[tuple[Expr]]
        sub_B: list 
        sub_beta: list[tuple[Expr]]
        sub_omega: tuple[Expr]
        sub_omegas: list[tuple[Expr]] 
        sub_phi: list[tuple[Expr]] 
        sub_scaling: list[tuple[Expr]] 
        sub_scaling_back: list[tuple[Expr]] 
        sub_sigma: tuple[Expr] 
        sub_t: list[tuple[Expr]]  
        sub_tS_to_t_func: list[tuple[Function]]   
        sub_x: list[tuple[Expr]]   
        sub_xO : list[tuple[Expr]]   
        sub_xO_t: list[tuple[Expr]]   
    
    def __init__(self, sub_t, sub_scaling, sub_omega, sub_sigma): 
        self.sub_t            = sub_t
        self.sub_scaling      = sub_scaling[0]
        self.sub_scaling_back = sub_scaling[1]
        self.sub_omega        = sub_omega
        self.sub_sigma        = sub_sigma
        
class Forcing_MMS:
    r"""
    Define the forcing on the system as
    
    - A forcing amplitude `F`
    
    - A scaling order `f_order` for the forcing
    
    - Forcing coefficients `fF`
    
    - Forcing terms (direct or parametric) `forcing_term`
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        F           : Symbol
        fF          : list[Union[Expr, int]]
        f_order     : list[int]
        forcing_term: list[Expr]
    
    def __init__(self, F, f_order, fF, forcing_term):
        self.F       = F
        self.f_order = f_order
        self.fF      = fF
        self.forcing_term = forcing_term
        
class Coord_MMS:
    """
    The coordinates used in the MMS.
    """      

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        A:     list[Symbol]
        B:     list[Symbol]
        a:     list[Symbol]
        at:    list[Symbol]
        beta:  list[Symbol]
        betat: list[Symbol]
        phi:   list[Symbol]
    
    def __init__(self, mms):
    
        self.A = [] # Complex amplitudes of the homogeneous leading order solutions
        self.B = [] # Real amplitudes of the particular leading order solutions (nonzero only if the forcing is hard)
        
        for ix in range(mms.ndof):
            self.A.append( Function(r'A_{}'.format(ix), complex=True)(*mms.tS[1:]) ) 
            self.B.append( symbols(r'B_{}'.format(ix), real=True) ) 

class Sol_MMS:
    """
    Solutions obtained when applying the MMS.
    """             

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        DA      : list[list[Expr]]
        fa      : list[Expr]
        faO     : list[list[Expr]]
        fbeta   : list[Expr]
        fbetaO  : list[list[Expr]]
        sec     : list[list[Expr]]
        x       : Union[list[str], list[Expr]]
        xO      : list[list[Expr]]
        xO_polar: list[list[Expr]]
    
    def __init__(self):
        pass


class Sol_transient:
    """
    Solutions obtained when applying the MMS and solving for the transient response.
    """             

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        solve_dof   : int
        IC          : dict
        a           : Expr
        beta        : Expr
        psi         : Expr
        omega       : Expr
        x           : Expr
    
    def __init__(self):
        pass

class Multiple_scales_system:
    r"""
    The multiple scales system.
    See :ref:`mms` for a detailed description of the dynamical system.

    Parameters
    ----------
    dynamical_system : Dynamical_system
        The dynamical system.

    eps : sympy.Symbol
        Small perturbation parameter :math:`\epsilon`.

    Ne : int
        Truncation order of the asymptotic series and order of the slowest time scale.

    omega_ref : sympy.Symbol
        Reference frequency :math:`\omega_{\textrm{ref}}` of the MMS.
        Not necessarily the frequency around which the MMS is going to be applied, see `ratio_omegaMMS`.

    sub_scaling : list of tuples
        Substitutions to do to scale the equations.
        Links small parameters to their scaled counterpart through :math:`\epsilon`.

    ratio_omegaMMS : int or sympy.Rational, optional
        Specify the frequency `omegaMMS` around which the MMS is going to be applied in terms of :math:`\omega_{\textrm{ref}}`.
        Denoting `ratio_omegaMMS` as :math:`r_{\textrm{MMS}}`, this means that

        .. math::
            \omega_{\textrm{MMS}} = r_{\textrm{MMS}} \omega_{\textrm{ref}}.

        Use ``ratio_omegaMMS=Rational(p,q)`` for

        .. math::
            q \omega_{\textrm{MMS}} = p \omega_{\textrm{ref}}

        to get better-looking results than the float :math:`p/q`.
        Default is 1.

    eps_pow_0 : int, optional
        Order of the leading order term in the asymptotic series of each oscillators' response.
        For the :math:`i^{\textrm{th}}` oscillator and denoting `eps_pow_0` as :math:`\lambda_0`, this means that

        .. math::
            x_i = \epsilon^{\lambda_0} x_{i,0} + \epsilon^{\lambda_0+1} x_{i,1} + \cdots.

        Default is 0.

    ratio_omega_osc : list of int or sympy.Rational or None, optional
        Specify the natural frequencies of the oscillators :math:`\omega_i` in terms of the reference frequency :math:`\omega_{\textrm{ref}}`.
        Denoting ``ratio_omega_osc[i]`` as :math:`r_i`, this means that

        .. math::
            \omega_i \approx r_i \omega_{\textrm{ref}}.

        Use ``ratio_omega_osc[i]=Rational(p,q)`` for

        .. math::
            q \omega_{i} \approx p \omega_{\textrm{ref}}

        to get better-looking results than the float :math:`p/q`.
        Default is `None` for each oscillator, so the :math:`\omega_i` are arbitrary and there are no internal resonances.
        Detuning can be introduced through the `detunings` keyword argument.

    detunings : list of sympy.Symbol or int, optional
        The detuning of each oscillator. Denoting ``detunings[i]`` as :math:`\delta_i`, this means that

        .. math::
            \omega_i = r_i \omega_{\textrm{ref}} + \delta_i.

        Default is 0 for each oscillator.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        Ne:              int
        coord:           Coord_MMS
        detunings:       list
        eps:             Symbol
        eps_pow_0:       int
        ndof:            int
        omega:           Symbol
        omegaMMS:        Expr
        omega_ref:       Symbol
        omegas:          list[Symbol]
        omegas_O0:       list[Expr]
        ratio_omegaMMS:  Union[int, Rational]
        ratio_omega_osc: list[Union[int, Rational]]
        sigma:           Symbol
        sol:             Sol_MMS
        sol_transient:   Sol_transient
        sub:             Substitutions_MMS
        t:               Symbol
        tS:              list[Symbol]
    
    def __init__(self,
             dynamical_system, eps, Ne, omega_ref, sub_scaling,
             ratio_omegaMMS = 1, eps_pow_0 = 0, ratio_omega_osc = None,
             detunings = 0
             ):
        """
        Transform the dynamical system introducing multiple time scales and frequency relations. For a complete transformation, asymptotic series must also be introduced. This is done in the children classes :class:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator` and :class:`~oscilate.MMS.mms_complex.Multiple_scales_complex`, which make use of an oscillator or complex form of the equations, respectively.
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
        self.tS, sub_t = self.time_scales()
        
        # Substitutions
        self.sub = Substitutions_MMS(sub_t, sub_scaling, sub_omega, sub_sigma)
        
        # Oscillators' frequencies (internal resonances and detuning)
        self.omegas          = dynamical_system.omegas
        self.ratio_omega_osc = ratio_omega_osc
        self.detunings       = detunings
        if self.ratio_omega_osc == None:
            self.ratio_omega_osc = [None]*self.ndof
        if self.detunings == 0:
            self.detunings = [0]*self.ndof
        self.oscillators_frequencies()        
        
        # Coordinates
        self.coord = Coord_MMS(self)
        self.polar_coordinates()
        
        # Solutions
        self.sol           = Sol_MMS()
        self.sol_transient = Sol_transient()
        
        
    def time_scales(self):
        r"""
        Define the time scales.

        Notes
        -----
        The time scales are defined as (see :ref:`mms`)
        
        .. math::
            t_0 = t, \quad t_1 = \epsilon t, \quad \cdots, \quad t_{N_e} = \epsilon^{N_e} t.
        
        Substitutions from the physical time :math:`t` to the time scales :math:`t_i, \; i=0, ..., N_e` are also prepared.
        Note that :math:`t_0` is refered-to as the fast time as it captures the oscillations. 
        Other time scales are refered-to as slow time as they capture the modulations.
        """
        
        tS  = []
        sub_t = []
        for it in range(self.Ne+1):
            tS.append(symbols(r't_{}'.format(it), real=True, positive=True))
            sub_t.append((self.eps**it * self.t, tS[it]))
            
        sub_t.reverse() # Start substitutions with the slowest time scale
        
        return tS, sub_t
    
    
    def oscillators_frequencies(self):
        r"""
        Gives the expression of every oscillator frequency in terms of the reference frequency, possibly with a detuning.
        
        Notes
        -----
        For the :math:`i^\textrm{th}` oscillator, this leads to 
        
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
        
        
    def _apply_MMS_shared(self):
        r"""
        Apply the MMS. 
        This method contains the operations that are common to the oscillator and complex form. These include a change of phase variables to make the system autonomous, a separation of solvability conditions into real and imaginary parts, and their reconstitution, resulting in the modulation equations.
        See :func:`~oscilate.MMS.mms_oscillator.apply_func` for details.
        """
        
        # Change the phase coordinates for autonomous purposes
        self.autonomous_phases()

        # Derive the modulation equations
        self.modulation_equations()

        # Reconstitution
        self.reconstitution() 
        

    def polar_coordinates(self):
        r"""
        Define the polar coordinates.

        Notes
        -----
        Define polar coordinates such that, for oscillator :math:`i`, the complex amplitude of the homogeneous leading order solution is defined as
        
        .. math::
            A_i(\boldsymbol{t}_s) = \frac{1}{2} a_i(\boldsymbol{t}_s) e^{\textrm{j} \phi_i(\boldsymbol{t}_s)},

        where :math:`a_i(\boldsymbol{t}_s)` and :math:`\phi_i(\boldsymbol{t}_s)` are the solution's amplitude and phase, respectively. 
        """
        
        self.coord.a   = [ Function(r'a_{}'.format(ix)   , real=True, positive=True)(*self.tS[1:]) for ix in range(self.ndof) ]
        self.coord.phi = [ Function(r'\phi_{}'.format(ix), real=True)               (*self.tS[1:]) for ix in range(self.ndof) ]
        self.sub.sub_A = [ ( self.coord.A[ix], Rational(1/2)*self.coord.a[ix]*exp(I*self.coord.phi[ix]) ) for ix in range(self.ndof)]
        
    def autonomous_phases(self):
        r"""
        Define phase coordinates that render an autonomous system.

        Notes
        -----
        Define new phase coordinates :math:`\beta_i` to transform nonautonomous equations into autonomous ones. 
        The :math:`\beta_i` are defined as
        
        .. math::
            \beta_i = \frac{r_i}{r_{\textrm{MMS}}} \sigma t_1 - \phi_i,
        
        where we recall that
        :math:`\omega = r_{\textrm{MMS}} \omega_{\textrm{ref}} + \epsilon \sigma` and :math:`\omega_{i,0} = r_i \omega_{\textrm{ref}}`. 
        See details on this choice in :ref:`mms`.
        """
        
        self.coord.beta   = [ Function(r'\beta_{}'.format(ix), real=True)(*self.tS[1:])                                            for ix in range(self.ndof) ]
        def_beta          = [ Rational(self.ratio_omega_osc[ix], self.ratio_omegaMMS) * self.sigma*self.tS[1] - self.coord.phi[ix] for ix in range(self.ndof) ]
        def_phi           = [ solve(def_beta[ix]-self.coord.beta[ix], self.coord.phi[ix])[0]                                       for ix in range(self.ndof) ]
        self.sub.sub_phi  = [ (self.coord.phi[ix], def_phi[ix])                                                                    for ix in range(self.ndof) ]
        self.sub.sub_beta = [ (self.coord.beta[ix], def_beta[ix])                                                                  for ix in range(self.ndof) ]

    def modulation_equations(self):
        r"""
        Derive the modulation equations of the polar coordinates system.
        
        Notes
        -----
        Derive the modulation equations of the polar coordinates system (defined in :func:`polar_coordinates` and :func:`autonomous_phases`) from the secular conditions. 
        For oscillator :math:`i` and at order :math:`j`, these are defined as
        
        .. math::
            \begin{cases}
            \textrm{D}_{j} a_i        & = f_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \textrm{D}_{j}\beta_i & = f_{\beta_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}),
            \end{cases}
        
        where :math:`\boldsymbol{a}` and :math:`\boldsymbol{\beta}` are vectors containing the polar amplitudes and phases.

        The aim here is to compute all the :math:`f_{a_i}^{(j)}` and :math:`f_{\beta_i}^{(j)}`.
        This is done by:
        
        #. Introducing polar coordinates in the secular terms
        
        #. Splitting the real and imaginary parts of the (complex) secular terms
        
        #. Using the autonomous phase coordinates
        
        #. Collecting the terms governing the slow amplitude and phase dynamics.
        """
        
        # Information
        print('Computing the modulation equations')

        # Initialisation
        sec_re    = [[] for dummy in range(self.ndof)] # Real part of the secular terms
        sec_im    = [[] for dummy in range(self.ndof)] # Imaginary part of the secular terms
        sub_re_im = [] # To overcome a sympy limitation: derivatives of real functions w.r.t. a real variable are not recognised as real
        
        faO    = [[] for dummy in range(self.ndof)] # Defined as      Di(a) = faO[i](a,beta)
        fbetaO = [[] for dummy in range(self.ndof)] # Defined as a*Di(beta) = fbetaO[i](a,beta)
        
        for io in range(0, self.Ne+1):
            
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
                    
                    # Derive the modulation equations at each order
                    faO[ix]   .append( solve(sec_im[ix][io].subs(self.sub.sub_B), self.coord.a[ix]   .diff(self.tS[io]))[0]                  )
                    fbetaO[ix].append( solve(sec_re[ix][io].subs(self.sub.sub_B), self.coord.beta[ix].diff(self.tS[io]))[0]*self.coord.a[ix] )
                    
        # Store the results
        self.sol.faO    = faO
        self.sol.fbetaO = fbetaO
                    
    def reconstitution(self):
        r"""
        Use the reconstitution method to combine the modulation equations at each order. This reconstitution is based on the chain rule relation

        .. math::
            \dfrac{\textrm{d}(\bullet)}{\textrm{d}t} = \sum_{i=0}^{N_e} \epsilon^{i} \mathrm{D}_i (\bullet) + \mathcal{O}(\epsilon^{N_e+1}).

        For oscillator :math:`i`, this results in the reconstituted modulation equations

        .. math::
            \begin{cases}
            \dfrac{\textrm{d} a_i}{\textrm{d} t}         & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \dfrac{\textrm{d} \beta_i}{\textrm{d} t} & = f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}),
            \end{cases}

        Note that some MMS approaches do not apply this reconstitution step.
        """

        # Introduce amplitudes and phases as functions of the physical time
        self.coord.at    = [ Function(r'a_{}'.format(ix), real=True, positive=True)(self.t) for ix in range(self.ndof) ]
        self.coord.betat = [ Function(r'\beta_{}'.format(ix), real=True, positive=True)(self.t) for ix in range(self.ndof) ]

        # Substitutions from functions of the tS (time scales) to t (physical time)
        self.sub.sub_tS_to_t_func = sum([[(a, at), (beta, betat)] for (a, at, beta, betat) in zip(*list(map(self.coord.__dict__.get, ["a","at","beta","betat"])))], [])

        # Reconstitute the modulation equations
        fa    = [0 for dummy in range(self.ndof)]  # Defined as      da/dt = fa(a,beta)
        fbeta = [0 for dummy in range(self.ndof)]  # Defined as a*dbeta/dt = fbeta(a,beta)
        for ix in range(self.ndof):
            for io in range(self.Ne+1):
                fa[ix]    += self.eps**io * self.sol.faO[ix][io]
                fbeta[ix] += self.eps**io * self.sol.fbetaO[ix][io]
            
        self.sol.fa    = [ fai.subs(self.sub.sub_tS_to_t_func) for fai in fa ]
        self.sol.fbeta = [ fbetai.subs(self.sub.sub_tS_to_t_func) for fbetai in fbeta ]

    def sol_x_polar(self, rewrite_polar=0):
        r"""
        Write the solutions using the polar coordinates and :math:`\cos` and :math:`\sin` functions.

        Parameters
        ----------
        rewrite_polar : str or int or list of int, optional
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
        harmonics = self.find_harmonics()
        collect_omega = [sin(h*self.omega*self.t) for h in harmonics] + [cos(h*self.omega*self.t) for h in harmonics]
        
        # Rewrite the solutions
        xO_polar = []
        x        = [0 for dummy in range(self.ndof)]
        for ix in range(self.ndof):
            xO_polar.append([])
            for io in rewrite_polar:
                xO_polar[ix].append( TR10(TR8((self.sol.xO[ix][io]
                                        .subs(self.sub.sub_A).doit().expand()
                                        .subs(self.sub.sub_phi).doit()) 
                                      .rewrite(cos).simplify()) 
                                      .subs(self.sub.sub_tS_to_t_func).subs(sub_t_back).subs(sub_sigma).simplify()) 
                                      .expand()
                                      .collect(collect_omega)
                                      )
            if rewrite_polar == range(self.Ne+1): # Construct the full response if relevant
                x[ix] = sum([self.eps**(io+self.eps_pow_0) * xO_polar[ix][io] for io in range(self.Ne+1)])
                x[ix] = x[ix].expand().collect(collect_omega) # Factor by the cos and sin terms
            else:
                x[ix] = "all solution orders were not rewritten in polar form"
        # Store
        self.sol.xO_polar = xO_polar
        self.sol.x        = x

    
    def find_harmonics(self):
        """
        Determine the harmonics contained in the MMS solutions. 

        Returns
        -------
        harmonics: list
            list of the harmonics appearing in the MMS solutions.
        """
        list_xO = list(itertools.chain.from_iterable(self.sol.xO))
        harmonics = []
        for xO_ix in list_xO:
            exponents = [exp_term.args[0].subs(self.tS[1],0) for exp_term in xO_ix.atoms(exp)]
            for exponent in exponents:
                if self.tS[0] in exponent.atoms() and im(exponent)>0:
                    harmonics.append( exponent/(I*self.omega_ref*self.tS[0]) / self.ratio_omegaMMS )
            
        harmonics = list(dict.fromkeys(harmonics))
        harmonics.sort()
        return harmonics
    
    def solve_transient(self, solve_dof=None, IC=dict()):
        r"""
        Solve the transient response of an oscillator.

        Parameters
        ----------
        solve_dof: None or int, optional
            The oscillator to solve for. 
            If `None`, no oscillator is solved for.
            Default is `None`.
        IC: dict(), optional
            See :func:`solve_slow_time`.
            
        Notes
        -----
        Find the transient solution for a given oscillator with the other oscillators' amplitude set to 0. Setting the other oscillators' amplitude to zero is done using the method :func:`substitution_solve_dof`.
        """

        # Conditions for not solving the forced response
        if solve_dof==None:
            return
        
        # Information
        print('Computing the transient response for oscillator {}'.format(solve_dof))
        
        # Store the oscillator that is solved for
        self.sol_transient.solve_dof = solve_dof
        
        # Set the other oscillator's amplitudes to zero
        self.substitution_solve_dof(solve_dof)

        # Compute the slow time evolutions
        self.solve_slow_time(IC=IC)

        # Introduce an absolute phase
        self.absolute_phase()

        # Compute the instantaneous frequency
        self.solve_instantaneous_frequency()

        #Compute the oscillator's motion
        self.solve_x_transient()

    def substitution_solve_dof(self, solve_dof):
        r"""
        Set every oscillator amplitude to 0 except the one to solve for.

        Notes
        -----
        If one wants to solve for :math:`a_i`, then the system is evaluated for :math:`a_j=0, \; \forall j \neq i`.
        """
        sub_solve = []
        for ix in range(self.ndof):
            if ix != solve_dof:
                sub_solve.append( (self.coord.a[ix], 0) )
                
        self.sub.sub_solve = sub_solve   

    def solve_slow_time(self, IC):
        r"""
        Compute the slow time evolution of the amplitude and phase of an oscillator.

        Parameters
        ----------
        IC: dict(), optional
            Initial conditions on :math:`a` and :math:`\beta`, given as a `dict` with keys `a` and `beta`.
            The values associated to these keys are `dict` of the form :math:`\{\dot{a}(0), a_{\text{IC}} \}` and :math:`\{\dot{\beta}(0), \beta_{\text{IC}} \}`
            Default is an empty `dict`, subsequently filled with some symbols for the :math:`a_{\text{IC}}` and :math:`\beta_{\text{IC}}`.
        """

        # Construct the equations to solve
        solve_dof = self.sol_transient.solve_dof
        Eqa    = self.coord.at[solve_dof].diff(self.t)                             - self.sol.fa[solve_dof]    # Equation on a
        Eqbeta = self.coord.at[solve_dof]*self.coord.betat[solve_dof].diff(self.t) - self.sol.fbeta[solve_dof] # Equation on beta

        # Get the initial conditions
        if not IC:
            ai    = symbols(r"a_i", real=True, positive=True) # Initial amplitude
            betai = symbols(r"\beta_i", real=True)            # Initial phase
            IC["a"]    = {self.coord.at[solve_dof].subs(self.t,0) : ai}         # Initial condition on a
            IC["beta"] = {self.coord.betat[solve_dof].subs(self.t,0) : betai}   # Initial condition on beta

        # Solve the transient
        a_sol    = dsolve(Eqa, self.coord.at[solve_dof], ics=IC["a"]).rhs 
        beta_sol = dsolve(Eqbeta, self.coord.betat[solve_dof], ics=IC["beta"]).rhs.subs(self.coord.at[solve_dof], a_sol).doit().expand().simplify()

        # Store results
        self.sol_transient.IC   = IC
        self.sol_transient.a    = a_sol
        self.sol_transient.beta = beta_sol
        self.sub.sub_solve += [(self.coord.a[solve_dof], self.sol_transient.a),
                               (self.coord.beta[solve_dof], self.sol_transient.beta)]

    def absolute_phase(self):
        r"""
        Introduce an absolute phase :math:`\psi` defined as (for oscillator :math:`i`)

        .. math::
            \psi = \frac{r_i}{r_{\textrm{MMS}}} \omega t - \beta_i,

        such that the leading order, homogeneous solution for oscillator :math:`i` takes the form 

        .. math::
            x^{\textrm{h}}_{i,0}(t) = a_i(t) \cos(\psi).

        Then, deduce the solution on :math:`\psi` from that on :math:`\beta`.
        The introduction of that absolute phase is useful to write the time response in a more compact way, and to derive the instantaneous frequency of oscillation.
        """
        solve_dof = self.sol_transient.solve_dof

        # Introduce psi and its relation with beta and omega
        psi = Function(r"\psi", real=True, positive=True)(self.t) # Absolute phase
        sub_psi = [(Rational(self.ratio_omega_osc[solve_dof], self.ratio_omegaMMS)*self.omega*self.t, psi+self.coord.betat[solve_dof])]                    # Substitution from the relative phase beta to the absolute one psi
        psi_sol = (Rational(self.ratio_omega_osc[solve_dof], self.ratio_omegaMMS)*self.omega*self.t - self.sol_transient.beta).subs([self.sub.sub_omega]).expand().simplify() # From the beta solution to the psi one

        # Store the results
        self.coord.psi   = psi
        self.sub.sub_psi = sub_psi
        self.sol_transient.psi = psi_sol
        self.sub.sub_solve.append( (self.coord.psi, self.sol_transient.psi) )

    def solve_instantaneous_frequency(self):
        r"""
        Compute the instantaneous frequency from the absolute phase through :math:`\omega_\text{NL} = \dot{\psi}`.
        """
        self.sol_transient.omega = self.sol_transient.psi.diff(self.t).simplify().factor() # Instantaneous oscillation frequency
        self.sub.sub_solve.append( (self.omega, self.sol_transient.omega) )

    def solve_x_transient(self):
        r"""
        Compute the displacement :math:`x` associated to a transient trajectory.
        """
        self.sol_transient.x = self.sol.x[self.sol_transient.solve_dof].subs(self.sub.sub_psi).simplify() # Adding .subs(self.sub.sub_solve) would result in too complex expressions

def Chain_rule_dfdt(f, tS, eps):
    r"""
    Apply the chain rule to express first order time derivatives in terms of the time scales' derivatives.
    
    Parameters
    ----------
    f: sympy.Function  
        Time scales-dependent function :math:`f(\boldsymbol{t})`.
    tS: list 
        Time scales :math:`\boldsymbol{t}^\intercal = [t_0, \cdots, t_{N_e}]`.
    eps: sympy.Symbol
        Small parameter :math:`\epsilon`.
    
    Returns
    -------
    dfdt: sympy.Function
        :math:`\mathrm{d} f/ \mathrm{d}t` expressed in terms of the time scales.
    
    Notes
    -----
    Consider a time scales-dependent function :math:`f(t_0, t_1, ...)`, where :math:`t_0` is the fast time and :math:`t_1, ...` are the slow times. 
    The Chain Rule is applied to give the expression of :math:`\mathrm{d} f/ \mathrm{d}t` in terms of the time scales.
    """
    
    Nt = len(tS)
    dfdt = 0
    for ii in range(Nt):
        dfdt += eps**ii * f.diff(tS[ii])
    
    return dfdt

def Chain_rule_d2fdt2(f, tS, eps):
    r"""
    Apply the chain rule to express second order time derivatives in terms of the time scales' derivatives.

    Parameters
    ----------
    f: sympy.Function  
        Time scales-dependent function :math:`f(\boldsymbol{t})`.
    tS: list 
        Time scales :math:`\boldsymbol{t}^\intercal = [t_0, \cdots, t_{N_e}]`.
    eps: sympy.Symbol
        Small parameter :math:`\epsilon`.
    
    Returns
    -------
    d2fdt2: sympy.Function
        :math:`\mathrm{d}^2 f/ \mathrm{d}t^2` expressed in terms of the time scales.
    
    Notes
    -----
    Consider a time scales-dependent function :math:`f(t_0, t_1, ...)`, where :math:`t_0` is the fast time and :math:`t_1, ...` are the slow times. 
    The Chain Rule is applied to give the expression of :math:`\mathrm{d}^2 f/ \mathrm{d}t^2` in terms of the time scales.
    """
    
    Nt = len(tS)
    d2fdt2 = 0
    for jj in range(Nt):
        for ii in range(Nt):
            d2fdt2 += eps**(jj+ii) * f.diff(tS[ii]).diff(tS[jj])
    
    return d2fdt2

def cartesian_to_polar(y, sub_polar, sub_phase=None):
    r"""
    Rewrites an expression or a matrix `y` from cartesian to polar coordinates.

    Parameters
    ----------
    y : sympy.Expr or sympy.Matrix
        A sympy expression or matrix written in cartesian coordinates.
    sub_polar: list
        A list of substitutions to perform to go from cartesian to polar coordinates. 
    sub_phase: list, optional
        Additional substitutions to try and get rid of phases, so that only the amplitude remains in the expression.
        Default is `None`.

    Returns
    -------
    yp: sympy.Expr or sympy.Matrix
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

def rescale(expr, mms):
    r"""
    Rescales a scaled expression.
    
    Parameters
    ----------
    expr: sympy.Expr
        An unscaled expression, i.e. an expression appearing at some order of :math:`\epsilon`.
    mms: Multiple_scales_system
        The mms system, containing substitutions to scale an expression.

    Returns
    -------
    expr_scaled: sympy.Expr
        The scaled expression.
    """
    expr_rescaled = expr.subs(*mms.sub.sub_sigma).subs(mms.sub.sub_scaling_back).simplify()
    return expr_rescaled
# %%
