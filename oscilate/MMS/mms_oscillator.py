# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the multiple scales system from the dynamical one with an oscillator form, i.e. 2nd order real differential equations.
It also carries out the MMS derivations, still in oscillator form.
"""

#%% Imports and initialisation
from sympy import (exp, I, conjugate, Rational, 
                   symbols, Symbol, Function, Expr, sympify, 
                   solve, dsolve, sympify)
from .. import sympy_functions as sfun
from .mms import (Multiple_scales_system, Forcing_MMS, 
                  Chain_rule_dfdt, Chain_rule_d2fdt2)
from typing import TYPE_CHECKING

#%% Classes and functions
class Coord_MMS_oscillator:
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

class Multiple_scales_oscillator(Multiple_scales_system):
    r"""
    The multiple scales system in oscillator form, i.e. 2nd order differential equations.
    See the parent class :class:`~oscilate.MMS.mms.Multiple_scales_system` for a description of the input parameters.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        coord   : Coord_MMS_oscillator
        EqO     : list[list[Expr]]
        EqO_t0  : list[list[Expr]]
        forcing : Forcing_MMS
        xO      : list[list[Function]]
        xO_t0   : list[list[Function]]

    def __init__(self,
             dynamical_system, eps, Ne, omega_ref, sub_scaling,
             ratio_omegaMMS = 1, eps_pow_0 = 0, ratio_omega_osc = None, detunings = 0
             ):
        """
        Transform the dynamical system introducing asymptotic series of the oscillators' coordinates. 
        """

        # Initialise the parent class
        super().__init__(dynamical_system, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS, eps_pow_0=eps_pow_0, ratio_omega_osc=ratio_omega_osc, detunings=detunings)
        
        # Information
        print('   Oscillator form (2nd order differential equations)')

        # Coordinates
        self.coord = Coord_MMS_oscillator(self)
        self.polar_coordinates()
        
        # Asymptotic series of x
        self.xO, sub_xO_t, sub_x = self.asymptotic_series(dynamical_system, eps_pow_0=self.eps_pow_0)
        self.sub.sub_xO_t = sub_xO_t
        self.sub.sub_x    = sub_x

        # Forcing
        self.forcing = self.forcing_MMS(dynamical_system)
        
        # Compute the MMS equations
        self.compute_EqO(dynamical_system)
        
    def asymptotic_series(self, dynamical_system, eps_pow_0=0):
        r"""
        Define the asymptotic series on the oscillators' coordinates x.

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
            
            for io in range(self.Ne+1):
            
                # Define time-dependent asymptotic terms
                xO_t.append(Function(r'x_{{{},{}}}'.format(ix,io), real=True)(self.t))
                x_expanded[ix] += self.eps**(io+eps_pow_0) * xO_t[io]
                
                # Define time scales-dependent asymptotic terms
                xO[ix].append(Function(xO_t[io].name, real=True)(*self.tS))
                
                # Substitutions from xO(t) and its time derivatives to xO(*tS) and its time scales derivatives
                sub_xO_t.extend( [(xO_t[io].diff(self.t,2), Chain_rule_d2fdt2(xO[ix][io], self.tS, self.eps)), 
                                  (xO_t[io].diff(self.t,1), Chain_rule_dfdt  (xO[ix][io], self.tS, self.eps)), 
                                  (xO_t[io]               , xO[ix][io])] )
            
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
        
    def apply_MMS(self, orders_polar=0):
        r"""
        Apply the MMS. 

        Parameters
        ----------
        orders_polar : str, int or list of int, optional
            The orders at which the solutions will be rewritten in polar form.
            See :meth:`~oscilate.MMS.mms.Multiple_scales_system.sol_x_polar`.
        
        Notes
        -----
        The application of the MMS is operated by successively calling the following methods.

        #. :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.system_t0`: An equivalent system is written in terms of the fast time scale :math:`t_0`. 
           This introduces the temporary unknowns :math:`\tilde{x}_{i,j}(t_0)` and allows the use of :func:`~sympy.solvers.ode.dsolve`.

        #. :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.sol_order_0`: Leading order solutions :math:`x_{i,0}(\boldsymbol{t})` are defined.

        #. :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.secular_analysis`: The leading order solutions are introduced in the equations and the secular terms at each order are identified. 
           Cancelling those secular terms is a condition for bounded solutions. 
           It leads to a system of complex modulation equations governing the slow time evolution of the complex amplitude of the homogeneous leading order solutions. 
           Each equation takes the form 
           
           .. math::
            \textrm{D}_{j} A_i(\boldsymbol{t}_s) = f_{A_i}^{(j)}(\boldsymbol{A}, t_1).

           After cancelling the secular terms the higher order equations are solved successively to express the higher order solutions :math:`x_{i,j}(\boldsymbol{t}),\; j>0` in terms of the leading order ones.

        #. :meth:`~oscilate.MMS.mms.Multiple_scales_system.autonomous_phases`: The phase coordinates are changed from :math:`\phi_i(\boldsymbol{t}_s)` to :math:`\beta_i(\boldsymbol{t}_s)` to cancel the slow time :math:`t_1` in the secular terms. This will be used afterwards to obtain an autonomous system.

        #. :meth:`~oscilate.MMS.mms.Multiple_scales_system.apply_MMS_shared`: The functions that are shared among various MMS forms (oscillator and complex). This calls

            #. :meth:`~oscilate.MMS.mms.Multiple_scales_system.autonomous_phases`: The phase coordinates are changed from :math:`\phi_i(\boldsymbol{t}_s)` to :math:`\beta_i(\boldsymbol{t}_s)` to cancel the slow time :math:`t_1` in the secular terms. This will be used afterwards to obtain an autonomous system.

            #. :meth:`~oscilate.MMS.mms.Multiple_scales_system.modulation_equations`: The secular conditions are split into real and imaginary parts, polar coordinates are used and the autonomous phases are introduced,
            resulting in an autonomous system of modulation equations on polar coordinates. 
            Equations come by two, one representing the amplitude modulation while the other represents the phase's, such that
            
            .. math::
                \begin{cases}
                \textrm{D}_{j} a_i(\boldsymbol{t}_s) & = f_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
                a_i \textrm{D}_{j} \beta_i(\boldsymbol{t}_s) & = f_{\beta_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}).
                \end{cases}

            This is the key result of the MMS. 

            #. :meth:`~oscilate.MMS.mms.Multiple_scales_system.reconstitution`: The modulations on each time scale are now combined to reintroduce the physical time, resulting in a system of the form

            .. math::
                \begin{cases}
                \dfrac{\textrm{d}}{dt} a_i(t) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}), \\
                a_i \dfrac{\textrm{d}}{dt} \beta_i(t) & = f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}).
                \end{cases}

            This is known as the reconstitution method.

        #. :meth:`~oscilate.MMS.mms.Multiple_scales_system.sol_x_polar`: The leading and higher order solutions are rewritten in terms of polar coordinates using :math:`\cos` and :math:`\sin` functions.
        """
        
        # Write a temporary equivalent system depending only on t0
        self.system_t0()
        
        # Compute the solutions at order 0
        self.sol_order_0()

        # Analyse the secular terms
        self.secular_analysis()

        # Change phase coordinate and derive the modulation equations
        self.apply_MMS_shared() 
        
        # Write the x solutions in terms of polar coordinates
        self.sol_x_polar(orders_polar=orders_polar)


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
        sub_B  = [] # Substitutions from the particular solution amplitude Bi to its expression
        
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
        self.sub.sub_B  = sub_B
        
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
            sub_DA_sol.append([ (self.coord.A[ix].diff(self.tS[0]), 0) ] )
            sec       .append([ 0 ])
        
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
        self.sol.sec = sec      # Secular terms
        self.sol.DA  = DA_sol   # Solutions that cancel the secular terms
    
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
        
        
    
    