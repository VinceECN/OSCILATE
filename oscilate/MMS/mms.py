# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the multiple scales system from the dynamical one, and the application of the MMS.
"""

#%% Imports and initialisation
from sympy import (exp, I, conjugate, re, im, Rational, 
                   symbols, Symbol, Function, sympify, simplify, 
                   solve, dsolve, cos, sin, tan, sympify, Mod)
from sympy.simplify.fu import TR5, TR8, TR10
from .. import sympy_functions as sfun
import itertools

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

    Notes
    -----
    Description of the method of multiple scales.

    
    """
    
    def __init__(self, dynamical_system, eps, Ne, omega_ref, sub_scaling, 
                 ratio_omegaMMS=1, eps_pow_0=0, **kwargs):
        """
        Transform the dynamical system introducing asymptotic series and multiple time scales. 
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
        
        # Asymptotic series of x
        self.xO, sub_xO_t, sub_x = self.asymptotic_series(dynamical_system, eps_pow_0=self.eps_pow_0)
        
        # Substitutions required
        self.sub = Substitutions_MMS(sub_t, sub_xO_t, sub_x, sub_scaling, sub_omega, sub_sigma)
    
        # Forcing
        self.forcing = self.forcing_MMS(dynamical_system)
        
        # Oscillators' frequencies (internal resonances and detuning)
        self.omegas = dynamical_system.omegas
        self.ratio_omega_osc = kwargs.get("ratio_omega_osc", [None]   *self.ndof)
        self.detunings       = kwargs.get("detunings",       [0]*self.ndof)
        self.oscillators_frequencies()        
        
        # Coordinates
        self.coord = Coord_MMS(self)
        self.polar_coordinates()
        
        # Solutions
        self.sol = Sol_MMS()
        
        # Compute the MMS equations
        self.compute_EqO(dynamical_system)
        
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
    
    def asymptotic_series(self, dynamical_system, eps_pow_0=0):
        r"""
        Define the asymptotic series.

        Notes
        -----
        The series expansion for oscillator :math:`i` (and for a leading order term :math:`\epsilon^0 = 1`) takes the form (see :ref:`mms`)

        .. math::
            x_i(t) = x_{i,0}(t) + \epsilon x_{i,1}(t) + \epsilon^2 x_{i,2}(t) + \cdots + \epsilon^{N_e} x_{i,N_e}(t) + \mathcal{O}(\epsilon^{N_e+1}).

        On top of introducing the terms of the asymptotic series, this function prepares substitutions from

        1. dof :math:`x_i(t)` to temporary :math:`t`-dependent asymptotic terms :math:`x_{i,j}(t)`, such that

           .. math::
            x_i(t) = x_{i,0}(t) + \epsilon x_{i,1}(t) + \epsilon^2 x_{i,2}(t) + \cdots + \epsilon^{N_e} x_{i,N_e}(t), 

        2. Temporary :math:`x_{i,j}(t)` to the time scales-dependent terms :math:`x_{i,j}(\boldsymbol{t})`, such that
         
           .. math::
            x_i(\boldsymbol{t}) = x_{i,0}(\boldsymbol{t}) + \epsilon x_{i,1}(\boldsymbol{t}) + \epsilon^2 x_{i,2}(\boldsymbol{t}) + \cdots + \epsilon^{N_e} x_{i,N_e}(\boldsymbol{t}). 
        """
        
        # Initialisation
        xO         = [] # Terms x00, x01, ..., x10, x11, ... of the asymptotic series of the xi
        sub_xO_t   = [] # Substitutions from xO(t) to xO(*tS)
        x_expanded   = [] # x in terms of xO(t)
        sub_x        = [] # Substitutions from x to xO(t)
        
        for ix in range(self.ndof):
            
            # Initialisations 
            xO.append([])      # A list that will contain the different expansion orders of the current x
            xO_t = []          # Temporary xO(t) -> depend on the physical time t
            x_expanded.append(0) # Initialise the current x to 0
            
            for it in range(self.Ne+1):
            
                # Define time-dependent asymptotic terms
                xO_t.append(Function(r'x_{{{},{}}}'.format(ix,it), real=True)(self.t))
                x_expanded[ix] += self.eps**(it+eps_pow_0) * xO_t[it]
                
                # Define time scales-dependent asymptotic terms
                xO[ix].append(Function(xO_t[it].name, real=True)(*self.tS))
                
                # Substitutions from xO(t) and its time derivatives to xO(*tS) and its time scales derivatives
                sub_xO_t.extend( [(xO_t[it].diff(self.t,2), Chain_rule_d2fdt2(xO[ix][it], self.tS, self.eps)), 
                                  (xO_t[it].diff(self.t,1), Chain_rule_dfdt  (xO[ix][it], self.tS, self.eps)), 
                                  (xO_t[it]               , xO[ix][it])] )
            
            # Substitutions from x to xO(t)
            sub_x.append((dynamical_system.x[ix], x_expanded[ix]))
        
        return xO, sub_xO_t, sub_x
        
    def forcing_MMS(self, dynamical_system):
        r"""
        Rewrite the forcing terms.

        Parameters
        ----------
        dynamical_system : Dynamical_system
            The dynamical system.

        Notes
        -----
        The initial forcing terms are 
        
        .. math::
            f_{F,i}(\boldsymbol{x}(t), \dot{\boldsymbol{x}}(t), \ddot{\boldsymbol{x}}(t)) F \cos(\omega t), \quad i=1,...,N.
        
        Rewritting them involves

        1. Replacing the :math:`x_i(t)` by their series expansions written in terms of time scales,
        
        2. Scaling the forcing and the parameters in :math:`f_{F,i}` if any,
        
        3. Truncating terms whose order is larger than the largest order retained in the MMS,
        
        4. Rewrite :math:`\cos(\omega t)` as 
        
           .. math::
            \cos(\omega t) = \frac{1}{2} e^{\mathrm{j}(\omega_{\textrm{MMS}} + \epsilon \sigma)t} + cc = \frac{1}{2} e^{\mathrm{j}(\omega_{\textrm{MMS}}t_0 + \sigma t_1)} + cc.
        
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
            
        # Get the forcing term for each oscillator
        forcing_term = []
        fF           = []
        sub_t, sub_x, sub_xO_t = list(map(self.sub.__dict__.get,["sub_t", "sub_x", "sub_xO_t"]))

        for ix in range(self.ndof):
            fF.append( (dynamical_system.forcing.fF[ix].subs(self.sub.sub_scaling)
                        .subs(sub_x).doit().subs(sub_xO_t).expand().subs(sub_t).doit())
                        .series(self.eps, n=self.eps_pow_0+self.Ne+1).removeO())
            
            forcing_term.append( (fF[ix] * Rational(1,2)*F*self.eps**f_order)
                                .series(self.eps, n=self.eps_pow_0+self.Ne+1).removeO() * 
                                (exp( I*((omega*self.t).expand().subs(sub_t).doit())) + 
                                 exp(-I*((omega*self.t).expand().subs(sub_t).doit())) ) 
                                )
        
        forcing = Forcing_MMS(F, f_order, fF, forcing_term)
    
        return forcing
    
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
        
        
    def compute_EqO(self, dynamical_system):
        r"""
        Compute the system of equations for each oscillator at each order of :math:`\epsilon`. This system is described in :ref:`mms`.

        The output `EqO` is a list of lists:

        - The :math:`1^{\text{st}}` level lists are associated to the equations for each oscillator,
        
        - The :math:`2^{\text{nd}}` level lists are associated to the orders of :math:`\epsilon` from the lowest to the highest order.

        Parameters
        ----------
        dynamical_system : Dynamical_system
            The dynamical system.
        """
    
        # Equations with every epsilon appearing
        sub_t, sub_x, sub_xO_t, sub_scaling, sub_omegas = list(map(self.sub.__dict__.get,["sub_t", "sub_x", "sub_xO_t", "sub_scaling", "sub_omegas"]))
    
        Eq_eps = []
        for ix in range(self.ndof):
            Eq_eps.append( ((dynamical_system.Eq[ix].expand().subs(sub_omegas).doit().subs(sub_scaling).doit()
                             .subs(sub_x).doit().subs(sub_xO_t).doit().expand().subs(sub_t).doit())
                          .series(self.eps, n=self.eps_pow_0+self.Ne+1).removeO() 
                          - self.forcing.forcing_term[ix]).expand())
            
            if self.eps_pow_0 != 0: # Set the leading order to eps**0 = 1
                Eq_eps[-1] = (Eq_eps[-1] / self.eps**(self.eps_pow_0)).expand()
                
        # MMS equations system
        EqO = []
        for ix in range(self.ndof):
            
            # Initialise a list of the equations at each order. Start with the lowest order
            EqOx = [Eq_eps[ix].series(self.eps, n=1).removeO()] 
            
            # What has to be substracted to keep only the terms of order eps**io in equation at order io.
            retrieve_EqOx = EqOx[0] 
            
            # Feed EqOx with higher orders of epsilon
            for io in range(1, self.Ne+1):
                EqOx.append( ((Eq_eps[ix].series(self.eps, n=io+1).removeO() - retrieve_EqOx) / self.eps**io).simplify().expand() )
                
                # Update the terms that are to be substracted at order io+1
                retrieve_EqOx += self.eps**io * EqOx[io]
                
            EqO.append(EqOx)
            
        self.EqO = EqO
        
    def apply_MMS(self, rewrite_polar=0):
        r"""
        Apply the MMS. 

        Parameters
        ----------
        rewrite_polar : str, int or list of int, optional
            The orders at which the solutions will be rewritten in polar form.
            See :func:`sol_x_polar`.
        
        Notes
        -----
        The application of the MMS is operated by successively calling the following methods.

        #. :func:`system_t0`: An equivalent system is written in terms of the fast time scale :math:`t_0`. 
           This introduces the temporary unknowns :math:`\tilde{x}_{i,j}(t_0)` and allows the use of :func:`~sympy.solvers.ode.dsolve`.

        #. :func:`sol_order_0`: Leading order solutions :math:`x_{i,0}(\boldsymbol{t})` are defined.

        #. :func:`secular_analysis`: The leading order solutions are introduced in the equations and the secular terms at each order are identified. 
           Cancelling those secular terms is a condition for bounded solutions. 
           It leads to a system of modulation equations governing the slow time evolution of the complex amplitude of the homogeneous leading order solutions. 
           Each equation takes the form 
           
           .. math::
            \textrm{D}_{j} A_i(\boldsymbol{t}_s) = f_{A_i}^{(j)}(\boldsymbol{A}, t_1).

           After cancelling the secular terms the higher order equations are solved successively to express the higher order solutions :math:`x_{i,j}(\boldsymbol{t}),\; j>0` in terms of the leading order ones.

        #. :func:`autonomous_phases`: The phase coordinates are changed from :math:`\phi_i(\boldsymbol{t}_s)` to :math:`\beta_i(\boldsymbol{t}_s)` to cancel the slow time :math:`t_1` in the secular terms. This will be used afterwards to obtain an autonomous system.

        #. :func:`modulation_equations`: The secular conditions are split into real and imaginary parts, polar coordinates are used and the autonomous phases are introduced,
           resulting in an autonomous system of modulation equations on polar coordinates. 
           Equations come by two, one representing the amplitude modulation while the other represents the phase's, such that
           
           .. math::
            \begin{cases}
            \textrm{D}_{j} a_i(\boldsymbol{t}_s) & = f_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \textrm{D}_{j} \beta_i(\boldsymbol{t}_s) & = f_{\beta_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}).
            \end{cases}

           This is the key result of the MMS. The modulations on each time scale can be combined to reintroduce the physical time, resulting in a system of the form

           .. math::
            \begin{cases}
            \dfrac{\textrm{d}}{dt} a_i(t) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \dfrac{\textrm{d}}{dt} \beta_i(t) & = f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}).
            \end{cases}

           This is known as the reconstitution method.

        #. :func:`sol_x_polar`: The leading and higher order solutions are rewritten in terms of polar coordinates using :math:`\cos` and :math:`\sin` functions.
        """
        
        # Write a temporary equivalent system depending only on t0
        self.system_t0()
        
        # Compute the solutions at order 0
        self.sol_order_0()

        # Analyse the secular terms
        self.secular_analysis()

        # Change the phase coordinates for autonomous purposes
        self.autonomous_phases()

        # Derive the modulation equations
        self.modulation_equations()

        # Reconstitution
        self.reconstitution() 
        
        # Write the x solutions in terms of polar coordinates
        self.sol_x_polar(rewrite_polar=rewrite_polar)


    def system_t0(self):
        r"""
        Rewrite the system with the fast time scale as the only independent variable.

        Notes
        -----
        This is a trick to use :func:`~sympy.solvers.ode.dsolve`, which only accepts functions of 1 variable, to solve higher order equations. 
        The higher order equations are rewritten in terms of temporary coordinates :math:`\tilde{x}_{i,j}(t_0)` in place of :math:`x_{i,j}(\boldsymbol{t})`, with :math:`i,j` denoting the oscillator number and :math:`\epsilon` order, respectively. 
        This is equivalent to temporary considering that :math:`\boldsymbol{A}(\boldsymbol{t}_s)` does not depend on the slow times, which is of no consequence as there are no slow time derivatives appearing in the higher order equations at this stage. 
        Indeed, they were either substituted using the complex modulation equations, or they disappeared when eliminating the secular terms. 
        
        """
        
        xO_t0  = [] # t0-dependent variables xij(t0). Higher time scales dependency is ignored.
        EqO_t0 = [] # Equations at each order with only t0 as an explicit variable. Leads to a harmonic oscillator at each order with a t0-periodic forcing coming from lower order solutions.
        
        for ix in range(self.ndof):
            xO_t0 .append([ Function(r'\tilde{x_'+'{{{},{}}}'.format(ix,io)+'}', real=True)(self.tS[0]) for io in range(0, 1+self.Ne) ]) 
            EqO_t0.append([ self.EqO[ix][0].subs(self.xO[ix][0], xO_t0[ix][0]).doit() ])
            
        self.EqO_t0 = EqO_t0
        self.xO_t0  = xO_t0
        
        
    def sol_order_0(self):
        r"""
        Compute the leading-order solutions for each oscillator. 

        Notes
        -----
        For oscillator :math:`i`, the homogeneous solution takes the general form
        
        .. math::
            x_{i,0}^{\textrm{h}}(\boldsymbol{t}) = A_i(\boldsymbol{t}_s) e^{\textrm{j} \omega_{i,0} t_0} + cc,
        
        where :math:`A_i(\boldsymbol{t}_s)` is an unknown complex amplitude.
        
        If the oscillator is subject to hard forcing (i.e. forcing appears at leading order), then the particular solution
        
        .. math::
            x_{i,0}^{\textrm{p}}(t_0, t_1) = B_i e^{\textrm{j} \omega t} + cc = B_i e^{\textrm{j} (\omega_{\textrm{MMS}} t_0 + \sigma t_1)} + cc
        
        is also taken into account. :math:`B_i` is a time-independent function of the forcing parameters.
        """
        
        # Information
        print('Definition of leading order multiple scales solutions')
        
        # Initialisation
        xO0    = [] # leading order solutions
        sub_xO = [] # Substitutions from xij to its solution
        sub_B    = [] # Substitutions from the particular solution amplitude Bi to its expression
        
        # Compute the solutions
        for ix in range(self.ndof):
            
            # Homogeneous leading order solution 
            xO0_h_ix = (            self.coord.A[ix]*exp(I*self.omegas_O0[ix]*self.tS[0]) 
                          + conjugate(self.coord.A[ix]*exp(I*self.omegas_O0[ix]*self.tS[0])) )
            
            # Particular leading order solution - if the equation is not homogeneous (due to hard forcing)
            if not self.EqO[ix][0] == self.xO[ix][0].diff(self.tS[0],2) + (self.omegas_O0[ix])**2 * self.xO[ix][0]:
                hint="nth_linear_constant_coeff_undetermined_coefficients"
                
                # General solution, containing both homogeneous and particular solutions
                xO0_sol_general = ( dsolve(self.EqO_t0[ix][0], self.xO_t0[ix][0], hint=hint) ).rhs
                
                # Cancel the homogeneous solutions
                C      = list(xO0_sol_general.atoms(Symbol).difference(self.EqO[ix][0].atoms(Symbol)))
                sub_IC = [(Ci, 0) for Ci in C]
                xO0_p_ix = xO0_sol_general.subs(sub_IC).doit()
                
                # Get the real amplitude of the particular solution
                exp_keys = list(xO0_p_ix.atoms(exp))
                if exp_keys:
                    sub_B.append( (self.coord.B[ix], xO0_p_ix.coeff(exp_keys[0])) )
                else:
                    print("Static hard forcing is currently not handled")
                
                # Rewrite the particular solution in terms of B for the sake of readability and computational efficiency
                xO0_p_ix = (          self.coord.B[ix]*exp(I*self.omega*self.t).subs([self.sub.sub_omega]).expand().subs(self.sub.sub_t).expand() + 
                            conjugate(self.coord.B[ix]*exp(I*self.omega*self.t).subs([self.sub.sub_omega]).expand().subs(self.sub.sub_t).expand()))
                    
            else:
                xO0_p_ix = sympify(0)
                
            # Total leading order solution
            xO0.append( xO0_h_ix + xO0_p_ix ) 
            sub_xO.append( ( self.xO[ix][0], xO0[ix] ) )
        
        # Store the solutions
        self.sol.xO = [[xO0_dof] for xO0_dof in xO0]
        self.sub.sub_xO = sub_xO
        self.sub.sub_B    = sub_B
        
    def secular_analysis(self):
        r"""
        Identify the secular terms in the MMS equations. 
        This allows to:

        1. Compute the modulation equations of the complex amplitudes :math:`A_i(\boldsymbol{t}_s)`, coming from the elimination of the secular terms,
        
        2. Derive nonsecular MMS equations, i.e. MMS equations with the secular terms cancelled,
        
        3. Use the nonsecular equations to express the higher order solutions :math:`x_{i,j}(\boldsymbol{t}),\; j>0` in terms of the :math:`\boldsymbol{A}(\boldsymbol{t}_s)`.
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
            sub_xO_t0 = [ (self.xO[ix][io], self.xO_t0[ix][io]) for ix in range(self.ndof) ]
            
            # Substitute the solutions at previous orders in the MMS equations and make it t0-dependent. Contains the secular terms.
            EqO_t0_sec = [ self.EqO[ix][io].subs(self.sub.sub_xO).subs(sub_xO_t0).doit() for ix in range(self.ndof) ] 
            
            # Find the secular terms and deduce the D(A) that cancel them
            dicE = [] 
            for ix in range(self.ndof):
                
                # Define the exponential corresponding to secular terms
                sub_exp = [(exp(I*self.omegas_O0[ix]*self.tS[0]), E)] # Substitute exp(I*omegas_O0*t0) by E to use sy.collect() in the following
                
                # Substitute the low order DA to get rid of all A derivatives except the current one
                EqO_t0_sec[ix] = sfun.sub_deep(EqO_t0_sec[ix], sub_DA_sol[ix])
                
                # Identify the secular term
                dicE_ix = EqO_t0_sec[ix].expand().subs(sub_exp).doit().expand().collect(E, evaluate=False)
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
                
            # Substitute the expression of the just computed DA in EqO_t0_sec to obtain nonsecular equations governing xO_t0 at the current order
            for ix in range(self.ndof):
                self.EqO_t0[ix].append(EqO_t0_sec[ix].subs(sub_DA_sol[ix]).doit().simplify())
            
            # Compute the x solution at order io in terms of the amplitudes A
            print('   Computing the higher order solutions at order {}'.format(io))
            for ix in range(self.ndof): 
                self.sol_higher_order(self.EqO_t0, self.xO_t0, io, ix)
            
        # Store the solutions
        self.sol.sec  = sec      # Secular terms
        self.sol.DA   = DA_sol   # Solutions that cancel the secular terms
    
    def sol_higher_order(self, EqO_t0, xO_t0, io, ix):
        r"""
        Compute the higher order solutions :math:`x_{i,j}(\boldsymbol{t}_s),\; j>0`.

        Parameters
        ----------
        EqO_t0 : list of list of sympy.Expr
            The MMS equations at each order and for each oscillator written with :math:`t_0` as the only independent variable. 
        xO_t0 : list of list of sympy.Function
            Oscillators' response at each order written in terms of :math:`t_0` only, :math:`\tilde{x}_{i,j}(t_0)`.
        io : int
            The order of :math:`\epsilon`.
        ix : int
            The dof number.
        """
        
        # Hint for dsolve()
        if not EqO_t0[ix][io] == xO_t0[ix][io].diff(self.tS[0],2) + (self.omegas_O0[ix])**2 * xO_t0[ix][io]:
            hint="nth_linear_constant_coeff_undetermined_coefficients"
        else:
            hint="default"
        
        # General solution, containing both homogeneous and particular solutions
        xO_sol_general = ( dsolve(EqO_t0[ix][io], xO_t0[ix][io], hint=hint) ).rhs
        
        # Cancel the homogeneous solutions
        C      = list(xO_sol_general.atoms(Symbol).difference(EqO_t0[ix][-1].atoms(Symbol)))
        sub_IC = [(Ci, 0) for Ci in C]
        
        # Append the solution for dof ix at order io
        self.sol.xO[ix].append(xO_sol_general.subs(sub_IC).doit())
        
        # Update the list of substitutions from the x to their expression
        self.sub.sub_xO.append( (self.xO[ix][io], self.sol.xO[ix][io]) )
        
        
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
        Derive the modulation equations of the polar coordinates system (defined in :func:`polar_coordinates` and :func:`autonomous_phases`) from the secular conditions. For oscillator :math:`i`, these are defined as
        
        .. math::
            \begin{cases}
            \dfrac{\textrm{d} a_i}{\textrm{d} t}         & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \dfrac{\textrm{d} \beta_i}{\textrm{d} t} & = f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}),
            \end{cases}
        
        where :math:`\boldsymbol{a}` and :math:`\boldsymbol{\beta}` are vectors containing the polar amplitudes and phases.

        The aim here is to compute all the :math:`f_{a_i}` and :math:`f_{\beta_i}`.
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
                    
                    # Derive the modulation equations at each order
                    faO[ix]   .append( solve(sec_im[ix][io], self.coord.a[ix]   .diff(self.tS[io]))[0]                  )
                    fbetaO[ix].append( solve(sec_re[ix][io], self.coord.beta[ix].diff(self.tS[io]))[0]*self.coord.a[ix] )
                    
                    # Reconstituted modulation equations
                    fa[ix]    += self.eps**io * faO[ix][io]
                    fbeta[ix] += self.eps**io * fbetaO[ix][io]
        
        # Store the results
        self.sol.faO    = faO
        self.sol.fbetaO = fbetaO
        self.sol.fa     = fa    # Modified in reconstitution() to account for physical time variables
        self.sol.fbeta  = fbeta # Modified in reconstitution() to account for physical time variables
                    
    def reconstitution(self):
        r"""
        Use the reconstitution method to combine the modulation equations at each order. This reconstitution is based on the chain rule relation

        .. math::
            \dfrac{\textrm{d}(\bullet)}{\textrm{d}t} = \sum_{i=0}^{N_e} \epsilon^{i} \mathrm{D}_i (\bullet) + \mathcal{O}(\epsilon^{N_e+1}).

        Note that some MMS approaches do not apply this reconstitution step.
        """

        # Introduce amplitudes and phases as functions of the physical time
        self.coord.at = [ Function(r'a_{}'.format(ix)   , real=True, positive=True)(self.t) for ix in range(self.ndof) ]
        self.coord.betat = [ Function(r'\beta_{}'.format(ix)   , real=True, positive=True)(self.t) for ix in range(self.ndof) ]

        # Substitutions from functions of the tS (time scales) to t (physical time)
        self.sub.sub_tS_to_t_func = sum([[(a, at), (beta, betat)] for (a, at, beta, betat) in zip(*list(map(self.coord.__dict__.get, ["a","at","beta","betat"])))], [])

        # Rewrite the previsouly obtained 
        self.sol.fa    = [ fa.subs(self.sub.sub_tS_to_t_func) for fa in self.sol.fa ]
        self.sol.fbeta = [ fbeta.subs(self.sub.sub_tS_to_t_func) for fbeta in self.sol.fbeta ]

    def sol_x_polar(self, rewrite_polar=0):
        r"""
        Write the solutions using the polar coordinates and :math:`\cos` and :math:`\sin` functions.

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
        harmonics = self.find_harmonics()
        collect_omega = [sin(h*self.omega*self.t) for h in harmonics] + [cos(h*self.omega*self.t) for h in harmonics]
        
        # Rewrite the solutions
        xO_polar = []
        x          = [0 for dummy in range(self.ndof)]
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
                x[ix] = sum([self.eps**(io+self.eps_pow_0) * xO_polar[ix][io] for io in range(self.Ne+1)]).simplify()
                list_cos_sin = list(x[ix].atoms(cos, sin)) # Get the sin and cosine terms from x
                x[ix] = x[ix].expand().collect(list_cos_sin) # Factor by the cos and sin terms
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
    
class Substitutions_MMS:
    """
    Substitutions used in the MMS.
    """
    
    def __init__(self, sub_t, sub_xO_t, sub_x, sub_scaling, sub_omega, sub_sigma): 
        self.sub_t            = sub_t
        self.sub_xO_t       = sub_xO_t
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
    
    - Forcing coefficients `fF`
    
    - Forcing terms (direct or parametric) `forcing_term`
    """
    
    def __init__(self, F, f_order, fF, forcing_term):
        self.F       = F
        self.f_order = f_order
        self.fF      = fF
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