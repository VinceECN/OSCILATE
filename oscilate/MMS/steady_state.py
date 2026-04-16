# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the steady state system from the multiple scales one, and the functions for steady state evaluation.
"""

#%% Imports and initialisation
from sympy import (Rational, symbols, Symbol, Matrix, Expr, solve, 
                   cos, sin, srepr, sympify, simplify, 
                   zeros, det, trace, eye, sqrt, pi)
from sympy.simplify.fu import TR10
from .. import sympy_functions as sfun
from .mms import cartesian_to_polar
from typing import Union

#%% Classes and functions
class Substitutions_SS:
    """
    Substitutions used in the steady state evaluations.
    """

    # Class-level annotations for pyreverse
    sub_B               : list
    sub_SS              : list[tuple]
    sub_cart            : list[tuple[Expr]]
    sub_free            : list[tuple]
    sub_phase           : list[tuple[Expr]]
    sub_polar           : list[tuple[Expr]]
    sub_scaling_back    : list[tuple[Expr]]
    sub_solve_forced    : list[tuple[Expr]]
    sub_solve_bbc       : list[tuple[Expr]]
    sub_solve_LC        : list[tuple[Expr]]
    
    def __init__(self, mms):
        
        self.sub_scaling_back = mms.sub.sub_scaling_back
        self.sub_B            = mms.sub.sub_B
        pass
        
    
class Forcing_SS:
    """
    Define the forcing on the system.
    """

    # Class-level annotations for pyreverse
    F      : Symbol
    f_order: int
    
    def __init__(self, mms):
        self.F            = mms.forcing.F
        self.f_order      = mms.forcing.f_order
        
class Coord_SS:
    """
    The coordinates used in the steady state analysis.
    """      
    # Class-level annotations for pyreverse
    a   : list[Symbol]
    beta: list[Symbol]
    p   : list[Symbol]
    q   : list[Symbol]
    
    def __init__(self):
        pass

class Sol_SS:
    """
    Solutions obtained when evaluating at steady state.
    """        

    # Class-level annotations for pyreverse       
    fa          : list[Expr]
    faO         : list[list[Expr]]
    fbeta       : list[Expr]
    fbetaO      : list[list[Expr]]
    fp          : list[Expr]
    fq          : list[Expr]
    x           : list[Expr]
    xO          : list[list[Expr]]
    
    def __init__(self, ss, mms):
        
        self.xO = []
        self.x  = []
        for ix in range(ss.ndof):
            self.xO.append( [xio.subs(ss.sub.sub_SS) for xio in mms.sol.xO_polar[ix]] )
            if not isinstance(mms.sol.x[ix], str):
                self.x.append( mms.sol.x[ix].subs(ss.sub.sub_SS) )
            else:
                self.x.append( "all solution orders were not rewritten in polar form" )

class Sol_forced:
    """
    Solutions obtained when evaluating the forced response of the system.
    """        

    # Class-level annotations for pyreverse       
    a2          : Union[list[Expr], None]
    F           : Union[Expr, None]
    fa          : Expr
    fbeta       : Expr
    cos_phase   : tuple[Expr]
    sigma       : list[Expr]
    sin_phase   : tuple[Expr]
    solve_dof   : Union[int, None]
    stab        : Stability
    
    def __init__(self):
        pass

class Sol_bbc:
    """
    Solutions obtained when evaluating the backbone curve of the forced response.
    """        

    # Class-level annotations for pyreverse       
    beta      : Expr
    omega     : Expr
    sigma     : Expr
    solve_dof : Union[int, None]
    x         : Expr
    xO        : list[Expr]

    def __init__(self, ss, mms):
        
        self.xO = []
        self.x  = []
        for ix in range(ss.ndof):
            self.xO.append( [xio.subs(ss.sub.sub_SS) for xio in mms.sol.xO_polar[ix]] )
            if not isinstance(mms.sol.x[ix], str):
                self.x.append( mms.sol.x[ix].subs(ss.sub.sub_SS) )
            else:
                self.x.append( "all solution orders were not rewritten in polar form" )

class Sol_LC:
    """
    Solutions obtained when evaluating the limit cycle (LC) of the system.
    """        

    # Class-level annotations for pyreverse       
    a         : Expr
    beta      : Expr
    omega     : Expr
    sigma     : Expr
    solve_dof : Union[int, None]
    x         : Expr
    xO        : list[Expr]

    def __init__(self, ss, mms):
        
        self.xO = []
        self.x  = []
        for ix in range(ss.ndof):
            self.xO.append( [xio.subs(ss.sub.sub_SS) for xio in mms.sol.xO_polar[ix]] )
            if not isinstance(mms.sol.x[ix], str):
                self.x.append( mms.sol.x[ix].subs(ss.sub.sub_SS) )
            else:
                self.x.append( "all solution orders were not rewritten in polar form" )

class Stability:
    """
    Stability analysis parameters and outputs.
    """                

    # Class-level annotations for pyreverse      
    Jsol            : Matrix
    Jsolc           : Matrix
    analysis_coord  : str
    bif_a           : list[Expr]
    bif_sigma       : list[Expr]
    blocks          : list[Matrix]
    blocks_bif_a    : list[list[Expr]]
    blocks_bif_sig  : list[list[Expr]]
    blocks_det      : list[Expr]
    blocks_eigvals  : list[list[Expr]]
    blocks_tr       : list[Expr]
    blocks_tr_a     : list[list[Expr]]
    blocks_tr_sig   : list[list[Expr]]
    det_Jsol        : Expr
    det_Jsolc       : Expr
    eigvals         : list[Expr]
    tr_Jsol         : Expr
    tr_Jsolc        : Expr
    tr_a            : list[Expr]
    tr_sigma        : list[Expr]
    
    def __init__(self):
        pass  

class Steady_state:
    r"""
    Steady state analysis of the multiple scales system.
    See :ref:`steady_state` for a detailed description of the dynamical system.

    Parameters
    ----------
    mms : Multiple_scales_system
        The multiple scales system.
    """

    # Class-level annotations for pyreverse
    coord:           Coord_SS
    eps:             Symbol
    forcing:         Forcing_SS
    ndof:            int
    omega:           Symbol
    omegaMMS:        Expr
    omega_ref:       Symbol
    ratio_omegaMMS:  Union[int, Rational]
    ratio_omega_osc: list[Union[int, Rational]]
    sigma:           Symbol
    sol_forced:      Sol_forced
    sol_bbc:         Sol_bbc
    sol_LC:          Sol_LC
    sub:             Substitutions_SS
    
    def __init__(self, mms):
        """
        Evaluate the MMS quantities at steady state. 
        """
        
        # Information
        print('Initialisation of the steady state analysis')
        
        # Small parameter
        self.eps = mms.eps
        
        # MMS frequencies of interest
        self.omega_ref      = mms.omega_ref
        self.ratio_omegaMMS = mms.ratio_omegaMMS
        self.sigma          = mms.sigma
        self.omega          = mms.omega
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
        self.polar_coordinates_SS(mms)
        
        # Solutions
        self.sol        = Sol_SS(self, mms)
        self.sol_forced = Sol_forced()
        self.sol_bbc    = Sol_bbc(self, mms)
        self.sol_LC     = Sol_LC(self, mms)
                
        # Evolution equations at steady state
        self.modulation_equations_SS(mms)
        
    
    def polar_coordinates_SS(self, mms):
        """
        Introduce time-independent amplitudes and phases (polar coordinates).
        """
        
        a, beta, sub_SS = [], [], []
        for ix in range(self.ndof):
            a   .append( symbols(r'a_{}'.format(ix),positive=True))
            beta.append( symbols(r'\beta_{}'.format(ix),real=True) )
            sub_SS .extend( [(mms.coord.a[ix]  , a[ix]), (mms.coord.beta[ix], beta[ix]),
                             (mms.coord.at[ix] , a[ix]), (mms.coord.betat[ix], beta[ix])] )
    
        self.coord.a    = a
        self.coord.beta = beta
        self.sub.sub_SS = sub_SS
        
    def modulation_equations_SS(self, mms):
        """
        Evaluate the modulation equations at steady state (polar system).
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
        
        # Check if the modulation equations are autonomous
        if 't_1' in srepr(fa) or 't_1' in srepr(fbeta):
            print("The modulation equations do not form an autonomous system")

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
                
        return sub_solve    

    def solve_forced(self, solve_dof=None):
        r"""
        Solve the forced response of an oscillator.

        Parameters
        ----------
        solve_dof: None or int, optional
            The oscillator to solve for. 
            If `None`, no oscillator is solved for.
            Default is `None`.
        
        Notes
        -----
        Find the steady state solution for a given oscillator with the other oscillators' amplitude set to 0.

        To do so, one must choose an oscillator to chose for, say oscillator :math:`i`. Then, the following methods are called: 
        
        #. :func:`substitution_solve_dof`: Set the other oscillators' amplitude to 0, i.e. :math:`a_j = 0 \; \forall j \neq i`.

        #. :func:`solve_phase`: express the oscillator's phase :math:`\beta_i` as a function of its amplitude :math:`a_i`. 

        #. :func:`solve_sigma`: find the expression of :math:`\sigma(a_i)`.

        #. :func:`solve_a`: find the expression of :math:`a_i (\sigma, F)`.
                
        #. :func:`solve_F`: find the expression of :math:`F(a_i)`.

        """
        
        # Conditions for not solving the forced response
        if solve_dof==None or self.forcing.F==0:
            return
        
        # Information
        print('Computing the forced response for oscillator {}'.format(solve_dof))
        
        # Store the oscillator that is solved for
        self.sol_forced.solve_dof = solve_dof
        
        # Set the other oscillator's amplitudes to zero
        self.sub.sub_solve_forced = self.substitution_solve_dof(solve_dof)
        
        # Phase response
        self.solve_phase()
        
        # Frequency (detuning) response
        self.solve_sigma()
        
        # Frequency response (in terms of oscillator's amplitude)
        self.solve_a()
        
        # Amplitude (forcing) respose
        self.solve_F()
                
    def solve_phase(self):
        r"""
        Find solutions for the oscillator's phase :math:`\beta_i`. 
        The solutions actually returned are :math:`\sin(\beta_i)` and :math:`\cos(\beta_i)`.
        """
        
        # Evaluate the modulation equations for a single oscillator responding
        fa_dof    = self.sol.fa[self.sol_forced.solve_dof]   .expand().subs(self.sub.sub_solve_forced)
        fbeta_dof = self.sol.fbeta[self.sol_forced.solve_dof].expand().subs(self.sub.sub_solve_forced)
        
        # Collect sin and cos terms in the modulation equations
        collect_sin_cos = list(fa_dof.atoms(cos, sin)) + list(fbeta_dof.atoms(cos, sin))
        collect_sin_cos = [item for item in collect_sin_cos if item.has(self.coord.beta[self.sol_forced.solve_dof])]
    
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
        if ( (len(list(list(dic_fa.keys()) + list(dic_fbeta.keys()))) > 4) or # cos and sin terms both appear in the same expression 
             (collect_sin_cos[0].args != collect_sin_cos[1].args) ): # Too many phases involved
            print('    No analytical solution implemented for this problem')
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
            print("   oscillator {} is not forced".format(self.sol_forced.solve_dof))
            return
    
        # Store the solutions
        self.sol_forced.fa    = fa_dof
        self.sol_forced.fbeta = fbeta_dof
        self.sol_forced.sin_phase = (collect_sin_cos[0], sin_phase)
        self.sol_forced.cos_phase = (collect_sin_cos[1], cos_phase)
        self.sub.sub_phase        = [self.sol_forced.sin_phase, self.sol_forced.cos_phase]
    
    def solve_sigma(self):
        r"""
        Solve the forced response in terms of the detuning :math:`\sigma`. 
        Returns :math:`\sigma(a_i)`.
        It is recalled that :math:`\omega = \omega_{\textrm{MMS}} + \epsilon \sigma`. 
        """
        
        sin_phase = self.sol_forced.sin_phase[1]
        cos_phase = self.sol_forced.cos_phase[1]
        
        # Equation on sigma
        Eq_sig = (sin_phase**2).expand() + (cos_phase**2).expand() - 1
    
        # Solve
        print('   Computing the frequency response')
        sol_sigma = sfun.solve_poly2(Eq_sig, self.sigma)
        if sol_sigma != None:
            sol_sigma = [sol_sigma_i.simplify() for sol_sigma_i in sol_sigma]
        
        # Store the solution
        self.sol_forced.Eq_sig = Eq_sig
        self.sol_forced.sigma  = sol_sigma
    
    def solve_a(self):
        r"""
        Solve the forced response in terms of the oscillator's amplitude.
        For readability, the output actually returned is :math:`a_i^2(\sigma, F)`. 
        """
        
        sin_phase = self.sol_forced.sin_phase[1]
        cos_phase = self.sol_forced.cos_phase[1]
        a         = self.coord.a[self.sol_forced.solve_dof]
        
        # Equation on a
        Eq_a = (sin_phase**2).expand() + (cos_phase**2).expand() - 1
        keys = Eq_a.expand().collect(a, evaluate=False)
        min_power = min(list(keys), key=lambda expr: sfun.get_exponent(expr, a))
        Eq_a = (Eq_a/min_power).expand()
        
        # Solve
        if set(Eq_a.collect(a, evaluate=False).keys()) in [set([1, a**4]), set([1, a**2, a**4])]:
            print("   Computing the response with respect to the oscillator's amplitude")
            sol_a2 = sfun.solve_poly2(Eq_a, a**2)
            sol_a2 = [sol_a2_i.simplify() for sol_a2_i in sol_a2]
        else:
            print("   Not computing the response with respect to the oscillator's amplitude as the equation to solve is not of 2nd degree")
            sol_a2 = None
        
        # Store the solutions
        self.sol_forced.Eq_a = Eq_a
        self.sol_forced.a2   = sol_a2
    
    def solve_F(self):
        r"""
        Solve the forced response in terms of the forcing amplitude :math:`F`. Returns :math:`F(a_i)`
        """
        
        sin_phase = self.sol_forced.sin_phase[1]
        cos_phase = self.sol_forced.cos_phase[1]
        F         = self.forcing.F
        
        # Equation on F
        Eq_F   = (((sin_phase*F).simplify()**2).expand() 
               + ( (cos_phase*F).simplify()**2).expand() 
               -    F**2).subs(self.sub.sub_B)
        keys = Eq_F.expand().collect(F, evaluate=False)
        min_power = min(list(keys), key=lambda expr: sfun.get_exponent(expr, F))
        Eq_F = (Eq_F/min_power).expand()
        
        # Solve
        if set(Eq_F.collect(F, evaluate=False).keys()) in [set([1, F**2]), set([1, F, F**2])]:
            print('   Computing the response with respect to the forcing amplitude - 2nd order polynomial in F')
            sol_F = abs(sfun.solve_poly2(Eq_F, F)[1])
        elif set(Eq_F.collect(F, evaluate=False).keys()) == set([1, F**2, F**4]):
            print('   Computing the response with respect to the forcing amplitude - 2nd order polynomial in F**2')
            sol_F = [sqrt(abs(F2)) for F2 in sfun.solve_poly2(Eq_F, F**2)]
        elif set(Eq_F.collect(F, evaluate=False).keys()) == set([1, F**2, F**4, F**6]):
            print('   Computing the response with respect to the forcing amplitude - 3rd order polynomial in F**2')
            sol_F = [sqrt(abs(F2)) for F2 in solve(Eq_F, F**2)]
        else:
            print('   Not computing the response with respect to the forcing amplitude as the equation to solve is not of 2nd degree')
            sol_F = None

        # Store the solutions
        self.sol_forced.Eq_F = Eq_F
        self.sol_forced.F    = sol_F    
    
    def solve_bbc(self, c=[], solve_dof=None):
        r"""
        Find the backbone curve (bbc) of a given oscillator with the other oscillators' amplitude set to 0.
        
        Parameters
        ----------
        c: list, optional
            Damping terms. They will be set to 0 to compute the backbone curve.
            Note that these are the scaled damping terms.
            Default is `[]`.
        solve_dof: None or int, oprtional
            The oscillator number to solve for. 
            If `None`, no oscillator is solved for.
            Default is `None`.

        Notes
        -----
        The backbone curve describes the frequency of free oscillations as a function of the oscillator's amplitude.
        In the presence of small damping (as in the case in the MMS) the frequency of free oscillations is close from the resonance frequency. 
        The backbone curve can therefore be interpreted as the *backbone* of the forced response. 

        The backbone curve of oscillator :math:`i` typically takes the form

        .. math::
            \omega_{\textrm{bbc}}^{(i)} = k\omega_{i} + f_{\textrm{bbc}}^{(i)} (a_i),

        where :math:`k=1,\; k<1,\; k>1` are associated to direct, superharmonic and subharmonic responses, respectively. 
        """
        
        if solve_dof==None:
            return
        
        # Information
        print('Computing the backbone curve for oscillator {}'.format(solve_dof))
        
        # Store the oscillator that is solved for
        self.sol_bbc.solve_dof = solve_dof

        # Set every oscillator amplitude to 0 except the one to solve for
        self.sub.sub_solve_bbc = self.substitution_solve_dof(solve_dof)

        # Phase solution
        self.sol_bbc.beta = pi/2
        self.sub.sub_solve_bbc.append((self.coord.beta[solve_dof], self.sol_bbc.beta))
        
        # Substitutions for the free response
        if not isinstance(c, list): 
            c = [c]
        sub_free = [(self.forcing.F,0), *[(ci, 0) for ci in c]]
            
        # Compute the bbc equation
        if self.forcing.f_order != 0 and self.ratio_omegaMMS>=1: # using the free response
            self.sub.sub_solve_bbc += self.sub.sub_B + sub_free
            Eq_bbc = self.sol.fbeta[solve_dof].subs(self.sub.sub_solve_bbc)
        else:  # The backbone curve can be affected by the forcing (superharmonic resonance)
            self.sub.sub_solve_bbc += self.sub.sub_B 
            Eq_F   = self.sol.fa[solve_dof].subs(self.sub.sub_solve_bbc)
            sol_Fbbc = solve(Eq_F, self.forcing.F)
            if sol_Fbbc:
                Fbbc   = abs(sol_Fbbc[0])
                self.sub.sub_solve_bbc += [(self.forcing.F, Fbbc)]
            else:
                self.sub.sub_solve_bbc += sub_free
            Eq_bbc = self.sol.fbeta[solve_dof].subs(self.sub.sub_solve_bbc)

        # Compute the backbone curve
        self.sol_bbc.sigma = solve(Eq_bbc, self.sigma)[0].expand().collect(self.coord.a[solve_dof])
        self.sol_bbc.omega = self.omegaMMS + self.eps*self.sol_bbc.sigma
        
        # Store additional information
        self.sub.sub_free = sub_free
        self.sol_bbc.Eq_sig = Eq_bbc

        # Oscillator's motion
        self.solve_bbc_x()
        
    def solve_bbc_x(self):
        """
        Compute the displacement :math:`x` on the backbone curve.
        """
        self.sol_bbc.xO = [xio.subs(self.coord.beta[self.sol_bbc.solve_dof], self.sol_bbc.beta).simplify() for xio in self.sol.xO[self.sol_bbc.solve_dof]] 
        if not isinstance(self.sol.x[self.sol_bbc.solve_dof], str):
            self.sol_bbc.x  = self.sol.x[self.sol_bbc.solve_dof].subs(self.coord.beta[self.sol_bbc.solve_dof], self.sol_bbc.beta).simplify() 


    def solve_LC(self, solve_dof=None):
        r"""
        Find the limit cycle (LC) of a given oscillator with the other oscillators' amplitude set to 0.
        
        Parameters
        ----------
        solve_dof: None or int, oprtional
            The oscillator number to solve for. 
            If `None`, no oscillator is solved for.
            Default is `None`.
        """

        if solve_dof==None:
            return
        
        # Information
        print('Computing the limit cycle for oscillator {}'.format(solve_dof))
        
        # Store the oscillator that is solved for
        self.sol_LC.solve_dof = solve_dof

        # Set every oscillator amplitude to 0 except the one to solve for
        self.sub.sub_solve_LC = self.substitution_solve_dof(solve_dof)

        # Phase solution
        self.sol_LC.beta = 0
        self.sub.sub_solve_LC.append((self.coord.beta[solve_dof], self.sol_LC.beta))

        # Amplitude solution
        self.solve_LC_amplitude()

        # Frequency solution
        self.solve_LC_frequency()

        # Oscillator's motion
        self.solve_LC_x()

    def solve_LC_amplitude(self):
        """
        Compute the amplitude of the homogeneous, leading order solution on the limit cycle.
        """
        self.sol_LC.a   = solve(self.sol.fa[self.sol_LC.solve_dof], self.coord.a[self.sol_LC.solve_dof])[0] # Amplitude solution
        self.sub.sub_solve_LC.append((self.coord.a[self.sol_LC.solve_dof], self.sol_LC.a))

    def solve_LC_frequency(self):
        """
        Compute the oscillation frequency on the limit cycle.
        """
        self.sol_LC.sigma = solve(self.sol.fbeta[self.sol_LC.solve_dof], self.sigma)[0].subs(self.coord.a[0], self.sol_LC.a) # Detuning solution
        self.sol_LC.omega = (self.omegaMMS + self.eps*self.sol_LC.sigma) # Frequency solution
        self.sub.sub_solve_LC.append((self.omega, self.sol_LC.omega))

    def solve_LC_x(self):
        """
        Compute the displacement :math:`x` on the limit cycle.
        """
        self.sol_LC.xO = [xio.subs(self.coord.beta[self.sol_LC.solve_dof], self.sol_LC.beta).simplify() for xio in self.sol.xO[self.sol_LC.solve_dof]] 
        if not isinstance(self.sol.x[self.sol_LC.solve_dof], str):
            self.sol_LC.x  = self.sol.x[self.sol_LC.solve_dof].subs(self.coord.beta[self.sol_LC.solve_dof], self.sol_LC.beta).simplify() # Using .subs(self.sub.sub_solve_LC) would result in too long expressions

    def Jacobian_polar(self):
        r"""
        Compute the Jacobian of the modulation equations systems expressed in polar coordinates (see :func:`stability_analysis`).
        
        Returns
        -------
        J : sympy.Matrix
            Jacobian of the polar system.
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

        Notes
        -----
        The homogeneous leading order solution for oscillator :math:`i` expressed in polar coordinates takes the form

        .. math::
            \begin{split}
            x_{i,0}^{\textrm{h}}(t) & = a_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t - \beta_i \right), \\
                                    & = a_i \cos(\beta_i) \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) 
                                      + a_i \sin(\beta_i) \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right)
            \end{split}

        The polar coordinates are defined as

        .. math::
            \begin{cases}
            p_i = a_i \cos (\beta_i), \\
            q_i = a_i \sin (\beta_i),
            \end{cases}

        such that the leading order solution can be written as

        .. math::
            x_{i,0}^{\textrm{h}}(t) = p_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) + q_i \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right).
        """

        # Define the cartesian coordinates
        p = [symbols(r'p_{}'.format(ix), real=True) for ix in range(0, self.ndof)]
        q = [symbols(r'q_{}'.format(ix), real=True) for ix in range(0, self.ndof)]
        
        # Define relations between the polar and cartesian coordinates
        a, beta = list(map(self.coord.__dict__.get, ["a", "beta"]))
        sub_cart  = []
        sub_polar = []
        for ix in range(self.ndof):
            sub_cart.append( (a[ix]*cos(beta[ix])     , p[ix]) )
            sub_cart.append( (a[ix]*sin(beta[ix])     , q[ix]) )
            sub_cart.append( (a[ix]**2*cos(2*beta[ix]), p[ix]**2 - q[ix]**2) )
            sub_cart.append( (a[ix]**2*sin(2*beta[ix]), 2*p[ix]*q[ix]) )
            sub_cart.append( (a[ix]**2                , p[ix]**2 + q[ix]**2) )
            
            sub_polar.append( (p[ix], a[ix]*cos(beta[ix])) )
            sub_polar.append( (q[ix], a[ix]*sin(beta[ix])) )
    
        # Store the results
        self.coord.p = p
        self.coord.q = q
        self.sub.sub_cart  = sub_cart
        self.sub.sub_polar = sub_polar
        
    
    def modulation_equations_cartesian(self):
        r"""
        Compute the modulation equations of the cartesian coordinates system.

        Notes
        -----
        Write the modulation equations using the cartesian coordinates (defined in :func:`cartesian_coordinates`). 
        For oscillator :math:`i`, this results in
        
        .. math::
            \begin{cases}
            \dfrac{\textrm{d} p_i}{\textrm{d} t} & = f_{p_i}(\boldsymbol{p}, \boldsymbol{q}), \\
            \dfrac{\textrm{d} q_i}{\textrm{d} t} & = f_{q_i}(\boldsymbol{p}, \boldsymbol{q}),
            \end{cases}
        
        where :math:`\boldsymbol{p}(t)` and :math:`\boldsymbol{q}(t)` are vectors containing the cartesian coordinates.
        """
        
        # Compute the functions fp(p,q) and fq(p,q)
        fp = []
        fq = []
        
        a, beta   = list(map(self.coord.__dict__.get, ["a", "beta"]))
        fa, fbeta = list(map(self.sol.__dict__.get, ["fa", "fbeta"]))
        
        for ix in range(self.ndof):
            fp.append( TR10( ( fa[ix]*cos(beta[ix]) - fbeta[ix]*sin(beta[ix]) ).expand().simplify()
                            ).expand().subs(self.sub.sub_cart) )
            
            fq.append( TR10( ( fa[ix]*sin(beta[ix]) + fbeta[ix]*cos(beta[ix]) ).expand().simplify()
                            ).expand().subs(self.sub.sub_cart) )
        
        # Check if the a and beta have all been substituted
        substitution_OK = self.check_cartesian_substitutions(a, beta, fp, fq)
        
        # Try additional substitutions if the change of coordinates is incomplete
        if not substitution_OK:
            fp, fq = self.additional_cartesian_substitutions(fp, fq)
            
            # Check if the a and beta have all been substituted
            substitution_OK = self.check_cartesian_substitutions(a, beta, fp, fq)
        
        if not substitution_OK:
            print("   The substitution from polar to cartesian coordinates is incomplete")
        
        # Store the modulation equations
        self.sol.fp = fp
        self.sol.fq = fq
        
    def check_cartesian_substitutions(self, a, beta, fp, fq):
        r"""
        Check if substitutions from polar to cartesian coordinates are complete.
        
        Parameters
        ----------
        a: list of sympy.Symbol
            Amplitudes of the leading order solutions.
        beta: list of sympy.Symbol
            Phases of the leading order solutions.
        fp: list of sympy.Expr
            Evolution functions for the cartesian coordinates :math:`p_i`.
        fq: list of sympy.Expr
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
        
    def additional_cartesian_substitutions(self, fp, fq):
        r"""
        Reformulate the already-existing substitutions from polar to cartesian to try and substitute leftover polar terms.
        
        Parameters
        ----------
        fp: list of sympy.Expr
            Evolution functions for the cartesian coordinates :math:`p_i`.
            There are polar coordinates remaining.
        fq: list of sympy.Expr
            Evolution functions for the cartesian coordinates :math:`q_i`.
            There are polar coordinates remaining.

        Returns
        -------
        fp: list of sympy.Expr
            Evolution functions for the cartesian coordinates :math:`p_i`.
            Additional substitutions were performed to get rid of polar coordinates.
        fq: list of sympy.Expr
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
        Compute the Jacobian of the modulation equations systems expressed in cartesian coordinates (see :func:`stability_analysis`).
        
        Returns
        -------
        J : sympy.Matrix
            Jacobian of the cartesian system.
        """
        
        J = zeros(2*self.ndof,2*self.ndof)
        
        for ii, (fpi, fqi) in enumerate(zip(self.sol.fp, self.sol.fq)):
            for jj, (pj, qj) in enumerate(zip(self.coord.p, self.coord.q)):
                J[2*ii,2*jj]     = fpi.diff(pj)
                J[2*ii,2*jj+1]   = fpi.diff(qj)
                J[2*ii+1,2*jj]   = fqi.diff(pj)
                J[2*ii+1,2*jj+1] = fqi.diff(qj)
        
        return J

    def stability_analysis_forced(self, coord="cartesian", rewrite_polar=False, eigenvalues=False, bifurcation_curves=False, trace_curves=False, analyse_blocks=False, kwargs_bif=dict()):
        r"""
        Evaluate the stability of a steady state, forced solution. 
        See :ref:`stability` for a detailed description of the dynamical system.

        Parameters
        ----------
        coord: str, optional
            Either ``"cartesian"`` or ``"polar"``. 
            Specifies the coordinates to use for the stability analysis.
            ``"cartesian"`` is recommended as it prevents divisions by 0, which occur when at least one of the oscillator has a 0 ampliutude.
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
        trace_curves: bool, optional
            Compute the trace curves, i.e. the zeros of the Jacobian's trace.
            Default is `False`.
        analyse_blocks: bool, optional
            Analyse the diagonal blocks of the Jacobian rather than the Jacobian itself. This is relevant if the Jacobian is block-diagonal.
        kwargs_bif: dict, optional
            Passed to :func:`bifurcation_curves`
            Default is `dict()`.
        """
        
        # Check if a solution has been computed
        if not "sigma" in self.sol_forced.__dict__.keys():
            print("There is no solution to evaluate the stability of.")
            return
        
        # Information
        print("Evaluating the stability of the solution of oscillator {}".format(self.sol_forced.solve_dof))
        
        # Create a Stability instance
        self.sol_forced.stab = Stability()

        # Introduce the cartesian coordinates and modulation equations
        if coord == "cartesian":
            print("   Rewritting the system in cartesian coordinates")

            self.sol_forced.stab.analysis_coord = "cartesian"
            self.cartesian_coordinates()
            self.modulation_equations_cartesian()

        else:
            self.sol_forced.stab.analysis_coord = "polar"
        
        # Compute the Jacobian
        print("   Computing the Jacobian")
        if coord=="cartesian":
            J = self.Jacobian_cartesian()
        else:
            J = self.Jacobian_polar()
        
        # Set every oscillator amplitude to 0 except the one solved for
        if coord=="cartesian":
            for ix in range(self.ndof):
                if ix != self.sol_forced.solve_dof:
                    self.sub.sub_solve_forced.extend( [(self.coord.p[ix], 0), (self.coord.q[ix], 0)] )
         
        # Use the steady state solutions to perform substitutions
        self.sub.sub_solve_forced.extend( [(self.forcing.F*self.sol_forced.sin_phase[0], self.forcing.F*self.sol_forced.sin_phase[1]),
                                    (self.forcing.F*self.sol_forced.cos_phase[0], self.forcing.F*self.sol_forced.cos_phase[1])] )
        
        if coord=="cartesian": 
            if self.forcing.F in self.sol.fp[self.sol_forced.solve_dof].atoms(Symbol): 
                self.sub.sub_solve_forced.append( (self.forcing.F, solve(self.sol.fp[self.sol_forced.solve_dof].subs(self.sub.sub_solve_forced), self.forcing.F)[0]) )
            else:
                self.sub.sub_solve_forced.append( (self.forcing.F, solve(self.sol.fq[self.sol_forced.solve_dof].subs(self.sub.sub_solve_forced), self.forcing.F)[0]) )
        
        # Evaluate the Jacobian on the solution
        Jsol = simplify(J.subs(self.sub.sub_solve_forced)) 
        
        # Analyse the Jacobian
        tr_Jsol  = trace(Jsol).simplify()
        det_Jsol = det(Jsol).simplify()
        
        # Rewrite the results in polar form if cartesian coordinates were used (time consuming)
        if coord=="cartesian": 
            # Save cartesian results
            self.sol_forced.stab.Jsolc     = Jsol
            self.sol_forced.stab.tr_Jsolc  = tr_Jsol
            self.sol_forced.stab.det_Jsolc = det_Jsol
            
            # Write the results in polar form
            if rewrite_polar:
                print("   Expressing the stability results in polar coordinates")
                Jsol     = cartesian_to_polar(Jsol, self.sub.sub_polar, sub_phase=self.sub.sub_phase)
                tr_Jsol  = cartesian_to_polar(tr_Jsol, self.sub.sub_polar, sub_phase=self.sub.sub_phase)
                det_Jsol = cartesian_to_polar(det_Jsol, self.sub.sub_polar, sub_phase=self.sub.sub_phase)

        # Store results
        self.sol_forced.stab.Jsol     = Jsol
        self.sol_forced.stab.tr_Jsol  = tr_Jsol
        self.sol_forced.stab.det_Jsol = det_Jsol
        
        # Compute eigenvalues and bifurcation curves from the analysis of Jsol
        if not analyse_blocks:
            if eigenvalues:
                self.sol_forced.stab.eigvals = self.eigenvalues(Jsol)
            if bifurcation_curves:
                self.sol_forced.stab.bif_a, self.sol_forced.stab.bif_sigma = self.bifurcation_curves(det_Jsol, **kwargs_bif)
            if trace_curves:
                self.sol_forced.stab.tr_a, self.sol_forced.stab.tr_sigma = self.trace_curves(tr_Jsol, **kwargs_bif)

        # Analyse the blocks of Jsol
        if analyse_blocks:
            print("   Block analysis")

            if coord == "cartesian":
                Jsol = self.sol_forced.stab.Jsolc

            if sfun.is_block_diagonal(Jsol, 2):
                self.sol_forced.stab.blocks         = []
                self.sol_forced.stab.blocks_det     = []
                self.sol_forced.stab.blocks_tr      = []
                self.sol_forced.stab.blocks_eigvals = []
                self.sol_forced.stab.blocks_bif_a   = []
                self.sol_forced.stab.blocks_bif_sig = []
                self.sol_forced.stab.blocks_tr_a    = []
                self.sol_forced.stab.blocks_tr_sig  = []

                for idx in range(0, Jsol.rows, 2):
                    A = Jsol[idx:idx+2, idx:idx+2] 
                    self.sol_forced.stab.blocks.append(A)
                    detA = det(A)
                    trA  = trace(A) 
                    if coord=="cartesian":
                        detA = cartesian_to_polar(detA, self.sub.sub_polar, sub_phase=self.sub.sub_phase).factor()
                        trA  = cartesian_to_polar(trA, self.sub.sub_polar, sub_phase=self.sub.sub_phase).factor()
                    
                    self.sol_forced.stab.blocks_det.append(detA)
                    self.sol_forced.stab.blocks_tr.append(trA)

                    if eigenvalues:
                        eigvalsA = self.eigenvalues(A, detA=detA, trA=trA)
                        if coord=="cartesian":
                            eigvalsA = [cartesian_to_polar(eigval, self.sub.sub_polar, sub_phase=self.sub.sub_phase) for eigval in eigvalsA]
                        self.sol_forced.stab.blocks_eigvals.append(eigvalsA)

                    if bifurcation_curves:
                        bif_aA, bif_sigA = self.bifurcation_curves(detA, **kwargs_bif)
                        self.sol_forced.stab.blocks_bif_a.append(bif_aA)
                        self.sol_forced.stab.blocks_bif_sig.append(bif_sigA)

                    if trace_curves:
                        tr_aA, tr_sigA = self.bifurcation_curves(trA, **kwargs_bif)
                        self.sol_forced.stab.blocks_tr_a.append(tr_aA)
                        self.sol_forced.stab.blocks_tr_sig.append(tr_sigA)

            else:
                print("Trying to perform a block analysis while the Jacobian is not block-diagonal")

    def eigenvalues(self, A, detA=None, trA=None):
        r"""
        Computes the eigenvalues of a matrix :math:`\textrm{A}`.

        Parameters
        ----------
        A: sympy.Matrix
            The matrix whose eigenvalues are to be computed.
        detA: sympy.Expr
            Determinant of A.
            Default is `None`.
        trA: sympy.Expr
            Trace of A.
            Default is `None`.

        Returns
        -------
        eigvals: list
            The eigenvalues of :math:`\textrm{A}`.
        """

        print("   Computing eigenvalues")
        
        if A.shape == (2,2) and (detA, trA) != (None, None):
            eigvals = [Rational(1,2)* (trA - sqrt(trA**2 - 4*detA)), 
                       Rational(1,2)* (trA + sqrt(trA**2 - 4*detA))]
        else:
            lamb        = symbols(r"\lambda")
            eig_problem = A - lamb * eye(*A.shape)
            detEP       = eig_problem.det()
            eigvals     = solve(detEP, lamb)
        
        return eigvals
            
    def bifurcation_curves(self, detJ, sol_type="forced", var_a=False, var_sig=True, solver=sfun.solve_poly2):
        r"""
        Compute bifurcation curves, defined by the simple bifurcation points of the slow time system obtained for any forcing frequency and amplitude. 

        Parameters
        ----------
        detJ: sympy.Expr
            The determinant of the matrix.
        sol_type: str, optional
            The solution type, among {"forced", "bbc", "LC"}
            Default is forced.
        var_a: bool, optional
            Consider the :math:`i^{\textrm{th}}` oscillator's amplitude :math:`a_i` as the variable and find the bifurcation curve as an expression for :math:`a_i`.
            `detJ` is rarely a quadratic polynomial in :math:`a_i`, so this can rarely be computed easily.
            Default is `False`.
        var_sig: bool, optional
            Consider the detuning :math:`\sigma` as the variable and find the bifurcation curve as an expression for :math:`\sigma`.
            `detJ` is often a quadratic polynomial in :math:`\sigma`, so this can often be computed.
            Default is `True`.
        solver: function, optional
            The solver to use to compute the bifurcation curves.
            Available are solver called as `solve(expr, x)`, which solve `expr=0` for `x`.
            :func:`~sympy.solvers.solvers.solve` can be used but is sometimes slow.
            Default is :func:`~oscilate.sympy_functions.solve_poly2`.

        Returns
        -------
        bif_a : list
            The bifurcation curves for :math:`a_i^2`.
        bif_sig : list
            The bifurcation curves for :math:`\sigma`. 

        Notes
        -----
        Simple bifurcations are associated to real eigenvalues of the Jacobian matrix crossing the imaginary axis, hence passing through 0. As the determinant is the product of the eigenvalues, the zeros of the determinant give the bifurcation points of the system. 
        """
        
        print("   Computing bifurcation curves")

        # Get the solution to evaluate the stability of
        if sol_type=="forced": 
            sol = self.sol_forced
        elif sol_type=="bbc":
            sol = self.sol_bbc
        elif sol_type=="LC":
            sol = self.sol_LC
        else:
            "Wrong solution type"

        # Check if a stability analysis was performed
        if not "Jsol" in sol.stab.__dict__.keys():
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
        if var_a and sfun.check_solvability(detJ, self.coord.a[sol.solve_dof]**2):
            bif_a = solver(detJ, self.coord.a[self.sol_forced.solve_dof]**2)
        else:
            bif_a = []

        if var_sig and sfun.check_solvability(detJ, self.sigma):
            bif_sig = solver(detJ, self.sigma)
        else:
            bif_sig = []
        
        # Return
        return bif_a, bif_sig
    
    def trace_curves(self, trJ, sol_type="forced", var_a=False, var_sig=True, solver=sfun.solve_poly2):
        r"""
        Compute curves capturing the variations of the trace of the Jacobian. For 2 by 2 Jacobians, a negative trace with a positive determinant indicates a stable response.

        Parameters
        ----------
        trJ: sympy.Expr
            The trace of the matrix.
        sol_type: str, optional
            The solution type, among {"forced", "bbc", "LC"}
            Default is forced.
        var_a: bool, optional
            Consider the :math:`i^{\textrm{th}}` oscillator's amplitude :math:`a_i` as the variable and find the trace curve as an expression for :math:`a_i`.
            Default is `False`.
        var_sig: bool, optional
            Consider the detuning :math:`\sigma` as the variable and find the trace curve as an expression for :math:`\sigma`.
            Default is `True`.
        solver: function, optional
            The solver to use to compute the trace curves.
            Available are solver called as `solve(expr, x)`, which solve `expr=0` for `x`.
            :func:`~sympy.solvers.solvers.solve` can be used but is sometimes slow.
            Default is :func:`~oscilate.sympy_functions.solve_poly2`.

        Returns
        -------
        tr_a : list
            The trace curves for :math:`a_i^2`.
        tr_sig : list
            The trace curves for :math:`\sigma`. 
        """
        
        print("   Computing trace curves")

        # Get the solution to evaluate the stability of
        if sol_type=="forced": 
            sol = self.sol_forced
        elif sol_type=="bbc":
            sol = self.sol_bbc
        elif sol_type=="LC":
            sol = self.sol_LC
        else:
            "Wrong solution type"

        # Check if a stability analysis was performed
        if not "Jsol" in sol.stab.__dict__.keys():
            print("There was no stability analysis performed.")
            return

        # Check if the stability analysis is expressed in polar coordinates
        if "p" in self.coord.__dict__.keys():
            cartesian_coordinates = self.coord.p + self.coord.q
            symbols_tr = list(trJ.atoms(Symbol)) 
            for cartesian_coordinate in cartesian_coordinates:
                if cartesian_coordinate in symbols_tr:
                    print("Substitutions from cartesian back to polar coordinates were incomplete. \n ",
                          "Try other substitutions manually or compute the Jacobian's  trace using block partitions if possible")

        # Add curves related to the trace of the Jacobian if it is not a constant
        if var_a and self.coord.a[sol.solve_dof] in trJ.atoms(Symbol):
            tr_a   += solver(trJ, self.coord.a[sol.solve_dof]**2)
        else:
            tr_a = []
        
        if var_sig and self.sigma in list(trJ.atoms(Symbol)):
            tr_sig += solver(trJ, self.sigma)
        else:
            tr_sig = []
        
        # Return
        return tr_a, tr_sig