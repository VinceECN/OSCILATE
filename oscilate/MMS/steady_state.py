# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module defines the steady state system from the multiple scales one, and the functions for steady state evaluation.
"""

#%% Imports and initialisation
from sympy import (exp, I, conjugate, re, im, Rational, 
                   symbols, Symbol, Function, solve, dsolve,
                   cos, sin, tan, srepr, sympify, simplify, 
                   zeros, det, trace, eye, Mod, sqrt)
from sympy.simplify.fu import TR10
from .. import sympy_functions as sfun
from .mms import cartesian_to_polar

#%% Classes and functions
class Steady_state:
    r"""
    Steady state analysis of the multiple scales system.

    Parameters
    ----------
    mms : Multiple_scales_system
        The multiple scales system.

    Notes
    -----
    Description of the steady state analysis.

    ----------------------
    Steady state solutions
    ----------------------

    ^^^^^^^^^^^^^^^^^^^^^^^
    Steady state conditions
    ^^^^^^^^^^^^^^^^^^^^^^^

    At steady state, the solutions' amplitudes and phases are time-independent. One therefore has, for each oscillator :math:`i=1,...,N`,

    .. math::
        \begin{cases}
        \dfrac{\textrm{d}}{dt} a_i & = 0, \\
        \dfrac{\textrm{d}}{dt} \beta_i & = 0, \\
        \end{cases}
    
    and the homogeneous steady state solutions take the form

    .. math::
        x^{\textrm{h}}_{i,0}(t) = a_i \cos\left( \frac{r_i}{r_{\textrm{MMS}}} \omega t - \frac{r_i}{r_{\textrm{MMS}}} \beta_i \right).    
    
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    MMS evolution equations at steady state
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The steady state amplitudes :math:`\boldsymbol{a}` and phases :math:`\boldsymbol{\beta}` are governed by the evolution equations evaluated at steady state, which take the form

    .. math::
        \begin{cases}
        f_{a_0}(\boldsymbol{a}, \boldsymbol{\beta})     & = 0, \\
        f_{\beta_0}(\boldsymbol{a}, \boldsymbol{\beta}) & = 0, \\
        & \vdots \\
        f_{a_{N-1}}(\boldsymbol{a}, \boldsymbol{\beta})     & = 0, \\
        f_{\beta_{N-1}}(\boldsymbol{a}, \boldsymbol{\beta}) & = 0.
        \end{cases}

    This is now an algebraic system of equations. 

    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    MMS solutions at steady state
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Solving the evolution equations at steady state yields the steady state solutions :math:`\boldsymbol{a}, \boldsymbol{\beta}`. 
    The system can be solved directly for :math:`\boldsymbol{a}` and :math:`\boldsymbol{\beta}`, yielding explicit analytical solutions, but this is often complex as the system is nonlinear. 

    A possibility that is sometimes available to obtain **analytical solutions** is to rearrange the equations to isolate the phase terms :math:`\cos(f(\boldsymbol{\beta})),\; \sin(f(\boldsymbol{\beta}))` where :math:`f(\boldsymbol{\beta})` is a linear function of the :math:`\beta_i`. 
    Then the equations can be squared and summed up to obtain an equation on :math:`a_i` only and/or :math:`a_j` as a function of :math:`a_i`. The :math:`a_j` can be expressed as a function of :math:`a_i`, leading to a polynomial equation on :math:`a_i`.
    The resulting polynomial equation can rarely be solved directly as the polynomial involved are often of high order. 
    However, the polynomial is often quadratic in the detuning :math:`\sigma` and forcing amplitude :math:`F`. 
    It can therefore be solved for :math:`\sigma` and :math:`F` with :math:`a_i` seen as a parameter. 
    This yields an implicit solution for :math:`a_i(\sigma) \Rightarrow a_i(\omega)` and :math:`a_i(F)`, from which one can deduce the other amplitudes :math:`a_j` and phases :math:`\boldsymbol{\beta}`, thus reconstructing the oscillators' solutions.
    
    The processus described above is not always feasible. The blocking points are typically to 

    (i) get rid of phase terms :math:`\cos(f(\boldsymbol{\beta})),\; \sin(f(\boldsymbol{\beta}))` in the equations, 

    (ii) express every amplitude :math:`a_j` in terms of a single one, :math:`a_i`,

    (iii) end up with a polynomial of order 2 in :math:`\sigma` and :math:`F` 
    
    These difficulties become more pronounced when the system involves several oscillator, in which case the amplitudes and phases may only be **computed numerically**.

    To facilitate the derivation of an analytical solution, it is possible to consider the **backbone curve** (bbc) of the forced solution rather than the forced solution itself. 
    This bbc is computed in the absence of damping and forcing, therefore simplifying the system. Typically, this reduces the number of phase terms appearing. 
    The solving procedure is then the same as that described previously.

    ------------------
    Stability analysis
    ------------------
    The stability analysis is described in details in :func:`stability_analysis`.
    """
    
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
        self.sol = Sol_SS(self, mms)
        
        # Stability
        self.stab = Stab_SS()
        
        # Evolution equations at steady state
        self.evolution_equations_SS(mms)
        
    
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
        
    def evolution_equations_SS(self, mms):
        """
        Evaluate the evolution equations at steady state (polar system).
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
        
    def solve_phase(self):
        r"""
        Find solutions for the oscillator's phase :math:`\beta_i`. 
        The solutions actually returned are :math:`\sin(\beta_i)` and :math:`\cos(\beta_i)`.
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
            print("   oscillator {} is not forced".format(self.sol.solve_dof))
            return
    
        # Store the solutions
        self.sol.sin_phase = (collect_sin_cos[0], sin_phase)
        self.sol.cos_phase = (collect_sin_cos[1], cos_phase)
        self.sub.sub_phase = [self.sol.sin_phase, self.sol.cos_phase]
    
    def solve_sigma(self):
        r"""
        Solve the forced response in terms of the detuning :math:`\sigma`. 
        Returns :math:`\sigma(a_i)`.
        It is recalled that :math:`\omega = \omega_{\textrm{MMS}} + \epsilon \sigma`. 
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
        For readability, the output actually returned is :math:`a_i^2(\sigma, F)`. 
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
            print("   Computing the response with respect to the oscillator's amplitude")
            sol_a2 = sfun.solve_poly2(Eq_a, a**2)
            sol_a2 = [sol_a2_i.simplify() for sol_a2_i in sol_a2]
        else:
            print("   Not computing the response with respect to the oscillator's amplitude as the equation to solve is not of 2nd degree")
            sol_a2 = None
        self.sol.sol_a2 = sol_a2
    
    def solve_F(self):
        r"""
        Solve the forced response in terms of the forcing amplitude :math:`F`. Returns :math:`F(a_i)`
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
            print('   Computing the response with respect to the forcing amplitude')
            sol_F = abs(sfun.solve_poly2(Eq_F, F)[1])
        else:
            print('   Not computing the response with respect to the forcing amplitude as the equation to solve is not of 2nd degree')
            sol_F = None
        self.sol.F = sol_F    
    
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
        
        # Set every oscillator amplitude to 0 except the one to solve for
        self.substitution_solve_dof(solve_dof)
        
        # Substitutions for the free response
        if not isinstance(c, list): 
            c = [c]
        sub_free = [(self.forcing.F,0), *[(ci, 0) for ci in c]]
        
        # Establish the backbone curve equation
        Eq_bbc = self.sol.fbeta[solve_dof].subs(self.sub.sub_solve).subs(self.sub.sub_B).subs(sub_free)
        
        # Compute the backbone curve
        self.sol.sigma_bbc = solve(Eq_bbc, self.sigma)[0].simplify()
        self.sol.omega_bbc = self.omegaMMS + self.eps*self.sol.sigma_bbc
        
        self.sub.sub_free = sub_free


    def Jacobian_polar(self):
        r"""
        Compute the Jacobian of the evolution equations systems expressed in polar coordinates (see :func:`stability_analysis`).
        
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
            sub_cart.append( (a[ix]**2                            , p[ix]**2 + q[ix]**2) )
            
            sub_polar.append( (p[ix], a[ix]*cos(beta[ix])) )
            sub_polar.append( (q[ix], a[ix]*sin(beta[ix])) )
    
        # Store the results
        self.coord.p = p
        self.coord.q = q
        self.sub.sub_cart  = sub_cart
        self.sub.sub_polar = sub_polar
        
    
    def evolution_equations_cartesian(self):
        r"""
        Compute the evolution equations of the cartesian coordinates system.

        Notes
        -----
        Write the evolution equations using the cartesian coordinates (defined in :func:`cartesian_coordinates`). 
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
        
        # Store the evolution equations
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
        Compute the Jacobian of the evolution equations systems expressed in cartesian coordinates (see :func:`stability_analysis`).
        
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

    def stability_analysis(self, coord="cartesian", rewrite_polar=False, eigenvalues=False, bifurcation_curves=False, trace_curves=False, analyse_blocks=False, kwargs_bif=dict()):
        r"""
        Evaluate the stability of a steady state solution. 

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

        Notes
        -----

        ^^^^^^^^^^^^^^^^^^^^^
        Stability information
        ^^^^^^^^^^^^^^^^^^^^^

        Consider a steady state solution :math:`(\hat{\boldsymbol{a}} , \hat{\boldsymbol{\beta}})` such that, for :math:`i=1,...,N`,

        .. math::
            \begin{cases}
            f_{a_i}(\hat{\boldsymbol{a}}, \hat{\boldsymbol{\beta}})     & = 0, \\
            f_{\beta_i}(\hat{\boldsymbol{a}}, \hat{\boldsymbol{\beta}}) & = 0.
            \end{cases}

        The aim is to determine the stability state of that steady solution, which corresponds to a fixed point in the phase space. 

        ---------------
        Jacobian matrix
        ---------------
        Let's first introduce the vector of polar coordinates and polar evolution functions

        .. math::
            \boldsymbol{x}^{(\textrm{p})\intercal} & = [a_0, \beta_0, \cdots, a_{N-1}, \beta_{N-1}], \\
            \boldsymbol{f}^{(\textrm{p})\intercal} & = [f_{a_0}, f_{\beta_0}^*, \cdots, f_{a_{N-1}}, f_{\beta_{N-1}}^*].

        Note that the appearance of :math:`f_{\beta_i}^*` in :math:`\boldsymbol{f}^{(\textrm{p})}` requires :math:`a_i \neq 0`, which strongly constraints the type of steady state solutions that can be considered in the approach described below. 
        To relax this constraint, one can use a change of coordinates from polar to cartesian ones. 
        This will be discussed in following sections, after the description of this polar approach. 
        
        Using the vectors of polar coordinates and evolution functions, one can write the evolution equations system as

        .. math::
            \dfrac{\textrm{d} \boldsymbol{x}^{(\textrm{p})}}{\textrm{d}t} = \textrm{J}^{(\textrm{p})} \boldsymbol{x}^{(\textrm{p})},

        where we introduced the Jacobian matrix 

        .. math::
            \textrm{J}^{(\textrm{p})} 
            = \dfrac{\partial \boldsymbol{f}^{(\textrm{p})} }{ \partial \boldsymbol{x}^{(\textrm{p})} }
            = \begin{bmatrix}
            \frac{\partial f_{a_0}}{\partial a_0} & \frac{\partial f_{a_0}}{\partial \beta_0} & \cdots & \frac{\partial f_{a_0}}{\partial \beta_{N-1}} \\
            \frac{\partial f_{\beta_0}^*}{\partial a_0} & \frac{\partial f_{\beta_0}^*}{\partial \beta_0} & \cdots & \frac{\partial f_{\beta_0}^*}{\partial \beta_{N-1}} \\
            \vdots & \vdots & \ddots & \vdots \\
            \frac{\partial f_{\beta_{N-1}}^*}{\partial a_0} & \frac{\partial f_{\beta_{N-1}}^*}{\partial \beta_0} & \cdots & \frac{\partial f_{\beta_{N-1}}^*}{\partial \beta_{N-1}}
            \end{bmatrix}.

        -----------------------------------------
        Perturbation of the steady state solution
        -----------------------------------------
        Let us now consider a small perturbation :math:`\tilde{\boldsymbol{x}}^{(\textrm{p})}` of the steady state solution such that

        .. math::
            \boldsymbol{x}^{(\textrm{p})} = \hat{\boldsymbol{x}}^{(\textrm{p})} + \tilde{\boldsymbol{x}}^{(\textrm{p})} \quad \Leftrightarrow \quad \tilde{\boldsymbol{x}}^{(\textrm{p})} = \boldsymbol{x}^{(\textrm{p})} - \hat{\boldsymbol{x}}^{(\textrm{p})}.

        Using a first order Taylor expansion for :math:`\textrm{d} \tilde{\boldsymbol{x}}^{(\textrm{p})} / \textrm{d} t`, one can write

        .. math::
            \dfrac{\textrm{d} \tilde{\boldsymbol{x}}^{(\textrm{p})}}{\textrm{d}t} = \left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}} \tilde{\boldsymbol{x}}^{(\textrm{p})} + \mathcal{O}(||\tilde{\boldsymbol{x}}^{(\textrm{p})}||^2),

        where :math:`\left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}}` denotes the Jacobian matrix evaluated on the steady state solution. 
        The perturbation solution takes the form

        .. math::
            \tilde{\boldsymbol{x}}^{(\textrm{p})} = \sum_{i=1}^{2N} C_i \boldsymbol{\psi}_i e^{\lambda_i t},

        where :math:`(\lambda_i, \boldsymbol{\psi}_i),\; i=1, ..., 2N` are the eigensolutions of the Jacobian (evaluated on :math:`\hat{\boldsymbol{x}}^{(\textrm{p})}`).
        
        -------------------
        Stability condition
        -------------------
        The steady state solution :math:`\hat{\boldsymbol{x}}^{(\textrm{p})}` is considered stable if a small perturbation :math:`\tilde{\boldsymbol{x}}^{(\textrm{p})}` vanishes in time, 
        such that solutions close from :math:`\hat{\boldsymbol{x}}^{(\textrm{p})}` are converging towards it. This condition is fulfilled if

        .. math::
            \Re[\lambda_i] < 0, \quad \forall i,

        meaning that all eigenvalues of the Jacobian evaluated on the steady state solution must have negative real parts.
        If this condition is not met, the system is either quasi stable or unstable.

        ------------
        Bifurcations
        ------------
        A bifurcation occurs when the stability state of the system changes. 
        
        #. Simple bifurcations occur when at least one eigenvalue crosses the imaginary axis through 0.
           Such bifurcations include saddle node and pitchfork bifurcations, which cause jumps of the response and the appearance of lower symmetry solutions, respectively.

        #. Neimark-Sacker bifurcations occur when a pair of complex conjugate eigenvalues with nonzero real parts cross the imaginary axis.
           These bifurcations lead to non periodic solutions. 

        Simple bifurcations can be detected by evaluating the sign of the Jacobian, as

        .. math::
            \det \left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}} = \prod_{i=1}^{2N} \lambda_i,

        thereby making :math:`\det \left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}}` an important stability indicator.
        Neimark-Sacker bifurcations are more difficult to detect. Information from the trace of :math:`\left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}}` can be considered, or the Routh-Hurwitz criterion can be used. This is not detailed here.

        ------------------
        Bifurcation curves
        ------------------
        Bifurcation curves are curves constructed by evaluating the coordinates of bifurcation points when varying one or more parameters.
        The stability state of a solution changes when the response curve crosses a bifurcation curve.

        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        Stability analysis in cartesian coordinates
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

        --------------------------------
        Limitations of polar coordinates
        --------------------------------
        As mentioned previously, the approach described above fails when one of the oscillator's leading order amplitude is 0.
        Indeed, the Jacobian is constructed using the evolution functions

        .. math::
                f_{\beta_i}^*(\boldsymbol{a}, \boldsymbol{\beta}) = \frac{f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta})}{a_i},

        which are defined only if :math:`a_i=0`. 
        This prevents evaluating the stability of

        - Trivial solutions, for which no oscillator responds,

        - 1 mode solutions, whose stability can be affected under perturbation from another mode.

        These limitations can be overcome using a change of coordinates.

        ---------------------
        Cartesian coordinates
        ---------------------
        The leading order homogeneous solution for oscillator :math:`i` can be written as

        .. math::
            \begin{split}
            x_{i,0}^{\textrm{h}}(t) & = a_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t - \beta_i \right), \\
                                    & = a_i \cos(\beta_i) \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) 
                                      + a_i \sin(\beta_i) \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right).
            \end{split}

        It then appears natural to introduce the cartesian coordinates

        .. math::
            \begin{cases}
            p_i = a_i \cos(\beta_i), \\
            q_i = a_i \sin(\beta_i),
            \end{cases}

        in order to rewrite the solution as

        .. math::
            x_{i,0}^{\textrm{h}}(t) = p_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) + q_i \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right).

        In the following it will be convenient to use the cartesian coordinates vectors

        .. math::
            \begin{aligned}
            \boldsymbol{p}(t)^\intercal & = [p_0(t), p_1(t), \cdots, p_{N-1}(t)], \\
            \boldsymbol{q}(t)^\intercal & = [q_0(t), q_1(t), \cdots, q_{N-1}(t)].
            \end{aligned}

        -----------------------------
        Cartesian evolution equations
        -----------------------------
        The cartesian evolution equations can be obtained from the polar ones. To do so, one can write

        .. math::
            \begin{cases}
            \dfrac{\textrm{d} p_i}{\textrm{d}t} & = \dfrac{\textrm{d} a_i}{\textrm{d}t} \cos(\beta_i) - a_i \sin(\beta_i) \dfrac{\textrm{d} \beta_i}{\textrm{d}t}, \\
            \dfrac{\textrm{d} q_i}{\textrm{d}t} & = \dfrac{\textrm{d} a_i}{\textrm{d}t} \sin(\beta_i) + a_i \cos(\beta_i) \dfrac{\textrm{d} \beta_i}{\textrm{d}t}.
            \end{cases}

        Then, by identification, one necessarily has

        .. math::
            \begin{cases}
            f_{p_i}(\boldsymbol{p}, \boldsymbol{q}) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}) \cos(\beta_i) - f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}) \sin(\beta_i), \\
            f_{q_i}(\boldsymbol{p}, \boldsymbol{q}) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}) \sin(\beta_i) + f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}) \cos(\beta_i),
            \end{cases}

        in order to write the evolution equations in cartesian coordinates

        .. math::
            \begin{cases}
            \dfrac{\textrm{d}}{dt} p_0(t) & = f_{p_0}(\boldsymbol{p}, \boldsymbol{q}), \\
            \dfrac{\textrm{d}}{dt} q_0(t) & = f_{q_0}(\boldsymbol{p}, \boldsymbol{q}), \\
            & \vdots \\
            \dfrac{\textrm{d}}{dt} p_{N-1}(t) & = f_{p_{N-1}}(\boldsymbol{p}, \boldsymbol{q}), \\
            \dfrac{\textrm{d}}{dt} q_{N-1}(t) & = f_{q_{N-1}}(\boldsymbol{p}, \boldsymbol{q}).
            \end{cases}

        ---------------
        Jacobian matrix
        ---------------

        As done previously for polar coordinates, let's introduce the vectors of cartesian coordinates and evolution equations as

        .. math::
            \boldsymbol{x}^{(\textrm{c})\intercal} & = [p_0, q_0, \cdots, p_{N-1}, q_{N-1}], \\
            \boldsymbol{f}^{(\textrm{c})\intercal} & = [f_{p_0}, f_{q_0}, \cdots, f_{p_{N-1}}, f_{q_{N-1}}].

        Then one can write

        .. math::
            \dfrac{\textrm{d} \boldsymbol{x}^{(\textrm{c})}}{\textrm{d}t} = \textrm{J}^{(\textrm{c})} \boldsymbol{x}^{(\textrm{c})},

        where we introduced the Jacobian matrix 

        .. math::
            \textrm{J}^{(\textrm{c})}
            = \dfrac{\partial \boldsymbol{f}^{(\textrm{c})} }{ \partial \boldsymbol{x}^{(\textrm{c})} }
            = \begin{bmatrix}
            \frac{\partial f_{p_0}}{\partial p_0} & \frac{\partial f_{p_0}}{\partial q_0} & \cdots & \frac{\partial f_{p_0}}{\partial q_{N-1}} \\
            \frac{\partial f_{q_0}}{\partial p_0} & \frac{\partial f_{q_0}}{\partial q_0} & \cdots & \frac{\partial f_{q_0}}{\partial q_{N-1}} \\
            \vdots & \vdots & \ddots & \vdots \\
            \frac{\partial f_{q_{N-1}}}{\partial p_0} & \frac{\partial f_{q_{N-1}}}{\partial q_0} & \cdots & \frac{\partial f_{q_{N-1}}}{\partial q_{N-1}}
            \end{bmatrix}.

        Note that there are no constraints related to an oscillator's amplitude being 0 here. 
        This cartesian coordinates approach therefore allows to investigate how the stability of a steady state solution is affected by a perturbation from an oscillator who's amplitude is 0 in that steady state solution.

        The stability analysis with :math:`\textrm{J}^{(\textrm{c})}` is carried out as described previously with :math:`\textrm{J}^{(\textrm{p})}`.
        """
        
        # Check if a solution has been computed
        if not "sigma" in self.sol.__dict__.keys():
            print("There is no solution to evaluate the stability of.")
            return
        
        # Information
        print("Evaluating the stability of the solution of oscillator {}".format(self.sol.solve_dof))
        
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
        
        # Set every oscillator amplitude to 0 except the one solved for
        if coord=="cartesian":
            for ix in range(self.ndof):
                if ix != self.sol.solve_dof:
                    self.sub.sub_solve.extend( [(self.coord.p[ix], 0), (self.coord.q[ix], 0)] )
         
        # Use the steady state solutions to perform substitutions
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
                self.stab.bif_a, self.stab.bif_sigma = self.bifurcation_curves(det_Jsol, **kwargs_bif)
            if trace_curves:
                self.stab.tr_a, self.stab.tr_sigma = self.trace_curves(tr_Jsol, **kwargs_bif)

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
                self.stab.blocks_tr_a    = []
                self.stab.blocks_tr_sig  = []

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
                        eigvalsA = self.eigenvalues(A, detA=detA, trA=trA)
                        if coord=="cartesian":
                            eigvalsA = [cartesian_to_polar(eigval, self.sub.sub_polar, sub_phase=self.sub.sub_phase) for eigval in eigvalsA]
                        self.stab.blocks_eigvals.append(eigvalsA)

                    if bifurcation_curves:
                        bif_aA, bif_sigA = self.bifurcation_curves(detA, **kwargs_bif)
                        self.stab.blocks_bif_a.append(bif_aA)
                        self.stab.blocks_bif_sig.append(bif_sigA)

                    if trace_curves:
                        tr_aA, tr_sigA = self.bifurcation_curves(trA, **kwargs_bif)
                        self.stab.blocks_tr_a.append(tr_aA)
                        self.stab.blocks_tr_sig.append(tr_sigA)

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
            
    def bifurcation_curves(self, detJ, var_a=False, var_sig=True, solver=sfun.solve_poly2):
        r"""
        Compute bifurcation curves, defined by the simple bifurcation points of the slow time system obtained for any forcing frequency and amplitude. 

        Parameters
        ----------
        detJ: sympy.Expr
            The determinant of the matrix.
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
            Default is :func:`~MMS.sympy_functions.solve_poly2`.

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
        
        # Return
        return bif_a, bif_sig
    
    def trace_curves(self, trJ, var_a=False, var_sig=True, solver=sfun.solve_poly2):
        r"""
        Compute curves capturing the variations of the trace of the Jacobian. For 2 by 2 Jacobians, a negative trace with a positive determinant indicates a stable response.

        Parameters
        ----------
        trJ: sympy.Expr
            The trace of the matrix.
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
            Default is :func:`~MMS.sympy_functions.solve_poly2`.

        Returns
        -------
        tr_a : list
            The trace curves for :math:`a_i^2`.
        tr_sig : list
            The trace curves for :math:`\sigma`. 
        """
        
        print("   Computing trace curves")

        # Check if a stability analysis was performed
        if not "Jsol" in self.stab.__dict__.keys():
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
        if var_a and self.coord.a[self.sol.solve_dof] in trJ.atoms(Symbol):
            tr_a   += solver(trJ, self.coord.a[self.sol.solve_dof]**2)
        else:
            tr_a = []
        
        if var_sig and self.sigma in list(trJ.atoms(Symbol)):
            tr_sig += solver(trJ, self.sigma)
        else:
            tr_sig = []
        
        # Return
        return tr_a, tr_sig

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
    The coordinates used in the steady state analysis.
    """      
    
    def __init__(self):
        pass
    
class Sol_SS:
    """
    Solutions obtained when evaluating at steady state.
    """                
    
    def __init__(self, ss, mms):
        
        self.xO = []
        self.x  = []
        for ix in range(ss.ndof):
            self.xO.append( [xio.subs(ss.sub.sub_SS) for xio in mms.sol.xO_polar[ix]] )
            if not isinstance(mms.sol.x[ix], str):
                self.x.append( mms.sol.x[ix].subs(ss.sub.sub_SS) )
            else:
                self.x.append( "all solution orders were not rewritten in polar form" )


class Stab_SS:
    """
    Stability analysis parameters and outputs.
    """                
    
    def __init__(self):
        pass  