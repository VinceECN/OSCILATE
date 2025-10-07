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
                   zeros, det, trace, eye, Mod, sqrt)
from sympy.simplify.fu import TR5, TR8, TR10
from . import sympy_functions as sfun
import numpy as np
import itertools
import matplotlib.pyplot as plt

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

    ---------------------------------
    Asymptotic series and time scales
    ---------------------------------

    The starting point is to introduce asymptotic series and multiple time scales in the initial dynamical system.
    The solution for oscillator :math:`i` is sought as a series expansion up to order :math:`N_e` (for a leading order term :math:`\epsilon^0 = 1`). This expansion takes the form

    .. math::
        x_i(t) = x_{i,0}(t) + \epsilon x_{i,1}(t) + \epsilon^2 x_{i,2}(t) + \cdots + \epsilon^{N_e} x_{i,N_e}(t) + \mathcal{O}(\epsilon^{N_e+1}).

    Time scales are introduced as follows:

    .. math::
        t_0 = t, \; t_1 = \epsilon t, \; t_2 = \epsilon^2 t, \cdots, t_{N_e} = \epsilon^{N_e} t,

    where :math:`t_0` is the fast time, i.e. the time used to describe the oscillations,
    while :math:`t_1, \; t_2,\; \cdots,\; t_{N_e}` are slow times, associated to amplitude and phase variations of the solutions in time. In addition, the chain rule gives

    .. math::
        \begin{aligned}
        \dfrac{\textrm{d}(\bullet)}{\textrm{d}t}     & = \sum_{i=0}^{N_e} \epsilon^{i} \dfrac{\partial(\bullet)}{\partial t_i} + \mathcal{O}(\epsilon^{N_e+1}), \\
        \dfrac{\textrm{d}^2(\bullet)}{\textrm{d}t^2} & = \sum_{j=0}^{N_e}\sum_{i=0}^{N_e} \epsilon^{i+j} \dfrac{\partial}{\partial t_j}\dfrac{\partial(\bullet)}{\partial t_i} + \mathcal{O}(\epsilon^{N_e+1}).
        \end{aligned}

    The introduction of asymptotic series and time scales are performed using :func:`asymptotic_series` and :func:`time_scales`.

    -------
    Scaling
    -------

    The construction of the MMS system requires a scaling of the parameters. Most scalings are already passed to the MMS through the `sub_scaling` parameter. 
    However, the natural frequencies also need to be scaled as they can contain both a leading order term and a detuning term.
    Natural frequencies :math:`\omega_i` are defined as a function of the reference frequency :math:`\omega_{\textrm{ref}}` through the `ratio_omega_osc` optional parameter, which is then used in :func:`oscillators_frequencies`.
    This allows to define internal resonance relations among the oscillators. 
    If these internal resonances are not perfect, detunings can be introduced through the `detunings` optional parameter, which needs to be scaled and part of the `sub_scaling` parameter.

    To write the MMS system it is convenient to introduce the leading order natural frequencies

    .. math::
        \omega_{i,0} = r_i \omega_{\textrm{ref}},

    where :math:`r_i` stands for ``ratio_omega_osc[i]``.

    --------------------------
    The multiple scales system
    --------------------------

    Introducing the asymptotic series, the time scales and the scaled parameters in the initial dynamical system (see :class:`~MMS.MMS.Dynamical_system`) results in :math:`N_e+1` dynamical systems, each one appearing at different orders of :math:`\epsilon`.
    Denoting time scales derivatives as
    
    .. math::
        \textrm{D}_i(\bullet) = \partial (\bullet) / \partial t_i, 
        
    introducing the vector of time scales
     
    .. math::
        \boldsymbol{t}^\intercal = [t_0, t_1, \cdots, t_{N_e}],

    where :math:`^\intercal` denotes the transpose, and the vectors of asymptotic coordinates

    .. math::
        \boldsymbol{x}_i(\boldsymbol{t})^\intercal = [x_{0,i}(\boldsymbol{t}), x_{1,i}(\boldsymbol{t}), \cdots, x_{N-1, i}(\boldsymbol{t})],

    which contains all the asymptotic terms of order :math:`i`, the MMS equations can be written as

    .. math::
        \begin{aligned}
        & \epsilon^0 \rightarrow \;
        \begin{cases}
        \textrm{D}_0 x_{0,0} + \omega_{0,0}^2 x_{0,0} & = f_{0,0}(t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,0} + \omega_{N-1,0}^2 x_{N-1,0} & = f_{N-1,0}(t_0, t_1),
        \end{cases} \\[15pt]
        & \epsilon^1 \rightarrow \;
        \begin{cases}
        \textrm{D}_0 x_{0,1} + \omega_{0,0}^2 x_{0,1} & = f_{0,1}  (\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,1} + \omega_{N-1,0}^2 x_{N-1,1} & = f_{N-1,1}(\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1),
        \end{cases} \\[15pt]
        & \epsilon^2 \rightarrow \;
        \begin{cases}
        \textrm{D}_0 x_{0,2} + \omega_{0,0}^2 x_{0,2} & = f_{0,2}  (\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_2 \boldsymbol{x}_0, \boldsymbol{x}_1, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_1, t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,2} + \omega_{N-1,0}^2 x_{N-1,2} & = f_{N-1,2}(\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_2 \boldsymbol{x}_0, \boldsymbol{x}_1, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_1, t_0, t_1),
        \end{cases} \\[10pt]
        & \hspace{3cm} \vdots \\[10pt]
        & \epsilon^{N_e} \rightarrow \;
        \begin{cases}
        \textrm{D}_0 x_{0,N_e} + \omega_{0,0}^2 x_{0,N_e} & = f_{0,N_e}  (\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_{N_e} \boldsymbol{x}_0, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_{N_e-1}, t_0, t_1), \\
        & \vdots \\
        \textrm{D}_0 x_{N-1,N_e} + \omega_{N-1,0}^2 x_{N-1,N_e} & = f_{N-1,N_e}(\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_{N_e} \boldsymbol{x}_0, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_{N_e-1}, t_0, t_1).
        \end{cases}
        \end{aligned}

    Consider oscillator :math:`i` at order :math:`j`:

    - The left-hand side term represents a harmonic oscillator of frequency :math:`\omega_{i,0}` oscillating with respect to the fast time :math:`t_0`.
    - The right-hand side term :math:`f_{i,j}` is analogous to a forcing generated by all combinations of terms that appear on oscillator :math:`i`'s equation at order :math:`\epsilon^j`.
      This can involve lower order terms :math:`x_{i,\ell}, \; \ell \leq j`, coupling terms :math:`x_{k, \ell}, \; k \neq j,\; \ell \leq j`, their derivatives and cross-derivatives with respect to the time scales, and physical forcing terms.
      The later is responsible for the dependency on :math:`t_0,\; t_1`. The reason why slower time scales are not involved will be explained in the following.

    Function :math:`f_{i,j}` tends to get increasingly complex as the order increases because the initial equations generate more high order terms than low order ones.

    This operation is performed using :func:`compute_EqMMS`.

    Note that internal resonance relations can be given through the `ratio_omega_osc` optional parameter, which is then used in :func:`oscillators_frequencies`.

    ---------------------
    Frequency of interest
    ---------------------

    The response of :math:`x_i` will be analysed at a frequency :math:`\omega`, defined as

    .. math::
        \omega = \omega_{\textrm{MMS}} + \epsilon \sigma,

    where :math:`\omega_{\textrm{MMS}}` is the **central MMS frequency**, controlled through the `ratio_omegaMMS` optional parameter and expressed in terms of `omega_ref`, 
    and :math:`\sigma` is a detuning about that frequency.
    In case the forced response is studied, :math:`\omega` corresponds to the forcing frequency.
    In case the free response is studied, :math:`\omega` corresponds to the frequency of free oscillations, which generates the backbone curve of the forced response.
    Note that :math:`\omega t = \omega_{\textrm{MMS}} t_0 + \sigma t_1`. This is the reason why the forcing only involves these two time scales in the right-hand side functions of the MMS system.

    ----------------------------------
    Iteratively solving the MMS system
    ----------------------------------
    
    The multiple scales system can be solved iteratively by solving successively the systems of equations at each order.
    
    ^^^^^^^^^^^^^^^^^^^^^^
    Leading order solution
    ^^^^^^^^^^^^^^^^^^^^^^

    The leading order solution for oscillator :math:`i` must satisfy
    
    .. math::
        \textrm{D}_0 x_{i,0} + \omega_{i,0}^2 x_{i,0} = f_{i,0}(t_0, t_1).
    
    It is sought as

    .. math::
        x_{i,0}(\boldsymbol{t}) = x_{i,0}^\textrm{h}(\boldsymbol{t}) + x_{i,0}^\textrm{p}(t_0, t_1),

    where :math:`x_{i,0}^\textrm{h}(\boldsymbol{t})` and :math:`x_{i,0}^\textrm{p}(t_0, t_1)` are the leading order homogeneous and particular sollutions, respectively.

    It is now conveninent to introduce the slow times vector

    .. math::
        \boldsymbol{t}_s^\intercal = [t_1, \cdots, t_{N_e}].

    This way, one can express the leading order solutions as

    .. math::
        \begin{cases}
        x_{i,0}^\textrm{h}(\boldsymbol{t}) & = A_i(\boldsymbol{t}_s) e^{\textrm{j} \omega_{i,0} t_0} + cc = |A_i(\boldsymbol{t}_s)| \cos(\omega_{i,0} t_0 + \arg{A_i(\boldsymbol{t}_s)}), 
        \\
        x_{i,0}^\textrm{p}(t_0, t_1) & = B_i e^{\textrm{j} \omega t} + cc = B_i e^{\textrm{j} (\omega_{\textrm{MMS}} t_0 + \sigma t_1)} + cc = |B_i| \cos(\omega_{\textrm{MMS}} t_0 + \sigma t_1 + \arg{B_i}),
        \end{cases}

    where :math:`A_i` is a slow time-dependent complex amplitude to be determined while :math:`B_i` is a time-independent function of the forcing parameters. 
    :math:`cc` denotes the complex conjugate. 
    Note that in most situations, forcing does not appear at leading order (i.e. forcing is weak), so :math:`B_i=0`.

    In the following it will be convenient to use the notations

    .. math::
        \begin{split}
        \boldsymbol{A}(\boldsymbol{t}_s)^\intercal & = [A_0(\boldsymbol{t}_s), A_1(\boldsymbol{t}_s), \cdots, A_{N-1}(\boldsymbol{t}_s)], \\
        \boldsymbol{B}^\intercal & = [B_0, B_1, \cdots, B_{N-1}].
        \end{split}

    The leading order solutions are defined in :func:`sol_order_0`.

    ^^^^^^^^^^^^^^^^^^^^^^
    Higher order solutions
    ^^^^^^^^^^^^^^^^^^^^^^

    Once the leading order solutions are computed, they can be injected in the :math:`1^\textrm{st}` higher order equations, where they (and their derivatives) appear as *forcing terms*, potentially together with physical forcing.
    The :math:`1^\textrm{st}` higher order equation for oscillator :math:`i` is 
    
    .. math::
        \textrm{D}_0 x_{i,1} + \omega_{i,0}^2 x_{i,1} = f_{i,1}(\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1).
    
    The forcing terms that involve oscillations at :math:`\omega_{i,0}` would force the oscillator on its natural frequency. Moreover, damping is always weak in the MMS, so damping terms of the form 
    :math:`c \textrm{D}_0 x_{i,1}` do not appear at this order. 
    The aforementioned forcing terms would thus lead to unbounded solutions, which is unphysical. 
    These forcing terms, called **secular terms**, must therefore be eliminated. 
    For instance, the :math:`1^\textrm{st}` higher order equation for oscillator :math:`i` with the secular terms cancelled is
    
    .. math::
        \textrm{D}_0 x_{i,1} + \omega_{i,0}^2 x_{i,1} = \bar{f}_{i,1}(\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1).
    
    where :math:`\bar{f}_{i,1}` is :math:`f_{i,1}` with the secular terms cancelled, i.e. without terms oscillating as :math:`\omega_{i,0}`.
    After cancelation of the secular terms, each oscillator's equation can be solved as a forced harmonic oscillator with the independent variable :math:`t_0`.
    
    Note that only the particular solutions are considered when solving higher order terms, i.e.

    .. math::
        x_{i,1}(\boldsymbol{t}) = \underbrace{x_{i,1}^\textrm{h}(\boldsymbol{t})}_{=0} + x_{i,1}^\textrm{p}(t_0, t_1).

    This choice can be justified if one assumes that initial conditions are of leading order. Though this is questionable, it is assumed here.
    
    The higher order solutions :math:`x_{i,1}(\boldsymbol{t})` are expressed as a function of the leading order unknown amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)`, 
    their slow time derivatives :math:`\textrm{D}_i\boldsymbol{A}(\boldsymbol{t}_s), \; i=1, ..., N_e`, and forcing terms if any (including the hard forcing amplitudes :math:`\boldsymbol{B}`).  
    
    This process is repeated successively at each order, i.e. the computed solutions are introduced in the next higher order system of equations, 
    secular terms are cancelled and the next higher order solutions are computed. 

    
    The secular terms are identified in :func:`secular_analysis` and the leading order solutions are computed in :func:`sol_higher_order`. 
    Note that :func:`sol_higher_order` is applied on equations with only :math:`t_0` as the independent variable so as to allow the use of :func:`~sympy.solvers.ode.dsolve`. 
    This is enforced using :func:`system_t0`, which temporarily ignores the dependency of :math:`\boldsymbol{A}(\boldsymbol{t}_s)` on the slow time scales.

    ^^^^^^^^^^^^^^^^
    Secular analysis
    ^^^^^^^^^^^^^^^^

    At this stage, the solutions are all expressed in terms of the unknown amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)` and their slow time derivatives :math:`\textrm{D}_1 A_i(\boldsymbol{t}_s),\; \cdots,\; \textrm{D}_{N_e} A_i(\boldsymbol{t}_s)`. 
    These can be obtained from the elimination of the secular terms (called **secular conditions**), as described below. 
    
    The :math:`i^\textrm{th}` MMS equation at :math:`1^\textrm{st}` higher order involves the slow time derivative :math:`\textrm{D}_1 A_i(\boldsymbol{t}_s)`, which appears 
    in the secular term. It is coming from the chain rule 

    .. math::
        \dfrac{\textrm{d}^2 x_i(t)}{\textrm{d}t^2} = \dfrac{\partial^2 x_{i,0}(\boldsymbol{t})}{\partial t_0^2} + \epsilon \dfrac{\partial^2 x_{i,1}(\boldsymbol{t})}{\partial t_0^2} + 2 \epsilon \dfrac{\partial^2 x_{i,0}(\boldsymbol{t})}{\partial t_0 \partial t_1} + \mathcal{O}(\epsilon^2).

    In addition, the :math:`\textrm{D}_1 A_j(\boldsymbol{t}_s),\; j\neq i` do not appear in the :math:`i^\textrm{th}` MMS equation as couplings among oscillators are weak. 
    It is thus possible to use the secular conditions in order to express the :math:`\textrm{D}_1 A_i(\boldsymbol{t}_s)` as a function of :math:`\boldsymbol{A}(\boldsymbol{t}_s)`. 
    
    This process can be done successively at each order to obtain the system of complex evolution equations
    
    .. math::
        \begin{cases}
        \textrm{D}_1 A_i(\boldsymbol{t}_s) & = f_{A_i}^{(1)}(\boldsymbol{A}, t_1), \\
        & \vdots \\
        \textrm{D}_{N_e} A_i(\boldsymbol{t}_s) & = f_{A_i}^{(N_e)}(\boldsymbol{A}, t_1).
        \end{cases}

    :math:`f_{A_i}^{(j)}(\boldsymbol{A}, t_1)` are functions governing the evolution of :math:`A_i` with respect to the slow time :math:`t_j`. 
    Note the dependency of :math:`f_{A_i}^{(j)}` on :math:`t_1` due to the possible presence of forcing. The complex evolution equations are derived in :func:`secular_analysis`.
    
    The above system of :math:`1^\textrm{st}` order PDE can theoretically be solved to obtain the complex amplitudes :math:`\boldsymbol{A}`. 
    However, this approach is not the prefered one as 

    - It is more convenient to deal with real variables than complex ones to get a physical meaning from the analysis,

    - It is more convenient to deal with autonomous systems (without the explicit :math:`t_1`-dependency) than nonautonomous ones,

    - The PDEs are complex.

    The first two points can be achieved introducing new coordinates, as described thereafter.

    -------------------
    Evolution equations
    -------------------

    ^^^^^^^^^^^^^^^^^
    Polar coordinates
    ^^^^^^^^^^^^^^^^^

    As discussed previously, it is more convenient to deal with real variables than complex ones. This can be done introducing the polar coordinates :math:`a_i` and :math:`\phi_i` for oscillator :math:`i` such that

    .. math::
        A_i(\boldsymbol{t}_s) = \dfrac{1}{2} a_i(\boldsymbol{t}_s) e^{\textrm{j} \phi_i(\boldsymbol{t}_s)}.

    :math:`a_i` and :math:`\phi_i` correspond to the amplitude and phase of the leading solution for oscillator :math:`i`, respectively. 
    With these new coordinates and introducing
    
    .. math::
        \begin{aligned}
        \boldsymbol{a}(\boldsymbol{t}_s)^\intercal & = [a_0(\boldsymbol{t}_s), a_1(\boldsymbol{t}_s), \cdots, a_{N-1}(\boldsymbol{t}_s)], \\
        \boldsymbol{\phi}(\boldsymbol{t}_s)^\intercal & = [\phi_0(\boldsymbol{t}_s), \phi_1(\boldsymbol{t}_s), \cdots, \phi_{N-1}(\boldsymbol{t}_s)],
        \end{aligned}

    the evolution equations on the :math:`A_i(\boldsymbol{t}_s)` can be split into real and imaginary terms, leading to
    
    .. math::
        \begin{aligned}
        & \epsilon^1 \rightarrow \;
        \begin{cases}
        \textrm{D}_1 a_i & = \hat{f}_{a_i}^{(1)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1), \\
        a_i \textrm{D}_1 \phi_i & = \hat{f}_{\phi_i}^{(1)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1), 
        \end{cases}
        \\[5pt]
        & \hspace{3cm} \vdots \\[5pt]
        & \epsilon^{N_e} \rightarrow \;
        \begin{cases}
        \textrm{D}_{N_e} a_i & = \hat{f}_{a_i}^{(N_e)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1), \\
        a_i \textrm{D}_{N_e} \phi_i & = \hat{f}_{\phi_i}^{(N_e)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1).
        \end{cases}
        \end{aligned}
    
    The :math:`\epsilon^j \rightarrow` indicate a system of 2 equations originating from the secular analysis at order :math:`j`.
    As one has

    .. math::
        \textrm{D}_j A_i = \textrm{D}_j \dfrac{1}{2} a_i e^{\textrm{j} \phi_i} = \dfrac{1}{2} \textrm{D}_j \left( a_i \right) e^{\textrm{j} \phi_i} + \textrm{j} \dfrac{1}{2} a_i e^{\textrm{j} \phi_i} \textrm{D}_j \left( \phi_i \right),

    it is convenient to pre-multiply the evolution equations on the :math:`A_i(\boldsymbol{t}_s)` by :math:`e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)}` or even :math:`\gamma e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)}` with, for instance, :math:`\gamma = 2`. 
    This avoids the presence of :math:`\cos(\phi_i)`, :math:`\sin(\phi_i)` and many :math:`1/2` terms in the evolution equations. 
    The evolution functions on polar coordinates are therefore defined as

    .. math::
        \begin{cases}
        \hat{f}_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1) & = \Re\left[ 2 e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)} f_{A_i}^{(j)}(\boldsymbol{A}, t_1) \right], \\
        \hat{f}_{\phi_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1) & = \Im\left[ 2 e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)} f_{A_i}^{(j)}(\boldsymbol{A}, t_1) \right].
        \end{cases}

    The evolution equations system on the polar coordinates involves only real variables, but functions :math:`\hat{f}_{a_i}^{(j)}` and :math:`\hat{f}_{\phi_i}^{(j)}` are
    still nonautonomous due to the explicit dependency on :math:`t_1`. 

    The polar coordinates are introduced in :func:`polar_coordinates`. The real evolution equations are only computed for the autonomous system, as described below. 

    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Autonomous phase coordinates
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    The presence of nonautonomous terms stems from forcing terms, which involve :math:`\cos(\sigma t_1 - \phi_i), \; \sin(\sigma t_1 - \phi_i)` in the polar evolution functions of oscillator :math:`i`. 
    A change of phase coordinates is required to make this autonomous.
    Moreover, the change of phase coordinate is necessary even in the absence of forcing for a convenient representation of the leading order solution. 
    Indeed, the solution for :math:`x^{\textrm{h}}_{i,0}(\boldsymbol{t})` written in terms of the current polar coordinates is

    .. math::
        x^{\textrm{h}}_{i,0}(\boldsymbol{t}) = a_i(\boldsymbol{t}_s) \cos(\omega_{i,0} t_0 + \phi_i(\boldsymbol{t}_s)).

    However, one would eventually like to express the oscillations of oscillator :math:`i` in terms of the frequency :math:`\omega`. To force its appearance, we recall that

    .. math::
        \omega_{i,0} = r_i \omega_{\textrm{ref}}, \quad \omega_{\textrm{ref}} = \frac{1}{r_{\textrm{MMS}}} \omega_{\textrm{MMS}}, \quad \textrm{and} \quad \omega_{\textrm{MMS}} = \omega - \epsilon \sigma.  
    
    Introducing this in the leading order solution leads to

    .. math::
        x^{\textrm{h}}_{i,0}(\boldsymbol{t}) = a_i(\boldsymbol{t}_s) \cos\left( \frac{r_i}{r_{\textrm{MMS}}} \omega t_0 - \frac{r_i}{r_{\textrm{MMS}}} \sigma t_1 + \phi_i(\boldsymbol{t}_s)\right).

    It therefore appears convenient to introduce the new phase coordinate :math:`\beta_i(\boldsymbol{t}_s)` as

    .. math::
        \beta_i(\boldsymbol{t}_s) = \frac{r_i}{r_{\textrm{MMS}}} \sigma t_1 - \phi_i(\boldsymbol{t}_s),

    which allows to write the leading order solution as

    .. math::
        x^{\textrm{h}}_{i,0}(\boldsymbol{t}) = a_i(\boldsymbol{t}_s) \cos\left( \frac{r_i}{r_{\textrm{MMS}}} \omega t_0 - \beta_i(\boldsymbol{t}_s)\right).

    In addition, and as discussed previously, the introduction of these new phase coordinates removes the explicit dependency of the evolution functions on :math:`t_1`, which was due to terms :math:`\cos(\sigma t_1 - \phi_i), \; \sin(\sigma t_1 - \phi_i)`.
    The forcing phase being zero (i.e. reference phase), the :math:`\beta_i(\boldsymbol{t}_s)` can be seen as the phases relative to the forcing. 

    Introducing the notation
    
    .. math::
        \boldsymbol{\beta}(\boldsymbol{t}_s)^\intercal = [\beta_1(\boldsymbol{t}_s), \beta_2(\boldsymbol{t}_s), \dots, \beta_{N_e}(\boldsymbol{t}_s)],
    
    the evolution equations can be rewritten as

    .. math::
        \begin{aligned}
        & \epsilon^1 \rightarrow \;
        \begin{cases}
        \textrm{D}_1 a_0(\boldsymbol{t}_s) & = f_{a_0}^{(1)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        a_0 \textrm{D}_1 \beta_0(\boldsymbol{t}_s) & = f_{\beta_0}^{(1)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        & \vdots \\
        \textrm{D}_1 a_{N-1}(\boldsymbol{t}_s) & = f_{a_{N-1}}^{(1)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        a_{N-1} \textrm{D}_1 \beta_{N-1}(\boldsymbol{t}_s) & = f_{\beta_{N-1}}^{(1)}(\boldsymbol{a}, \boldsymbol{\beta}), 
        \end{cases} \\[5pt]
        & \hspace{3cm} \vdots \\[5pt]
        & \epsilon^{N_e} \rightarrow \;
        \begin{cases}
        \textrm{D}_{N_e} a_0(\boldsymbol{t}_s) & = f_{a_0}^{(N_e)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        a_0 \textrm{D}_{N_e} \beta_0(\boldsymbol{t}_s) & = f_{\beta_0}^{(N_e)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        & \vdots \\
        \textrm{D}_{N_e} a_{N-1}(\boldsymbol{t}_s) & = f_{a_{N-1}}^{(N_e)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        a_{N-1} \textrm{D}_{N_e} \beta_{N-1}(\boldsymbol{t}_s) & = f_{\beta_{N-1}}^{(N_e)}(\boldsymbol{a}, \boldsymbol{\beta}).
        \end{cases}
        \end{aligned}

    The above system is the key result of the application of the MMS as it governs the evolution of leading order amplitudes and phases, on which all higher order solutions depend. 
    It can be solved numerically or analytically, if analytical solutions exist, though it is generally not the case. 
    It can also be rewritten in a more compact form as discussed in the following.
    
    The autonomous phase coordinates are introduced in :func:`autonomous_phases` and the evolution equations are computed in :func:`evolution_equations`.
    
    All solutions previously computed using the complex amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)` can be rewritten in terms of the polar coordinates :math:`\boldsymbol{a}(\boldsymbol{t}_s),\; \boldsymbol{\beta}(\boldsymbol{t}_s)` using :func:`sol_xMMS_polar`. 

    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    Reintroduction of the physical time
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Recalling the chain rule

    .. math::
        \dfrac{\textrm{d}(\bullet)}{\textrm{d}t} = \sum_{i=0}^{N_e} \epsilon^{i} \textrm{D}_i(\bullet), 

    and reintroducing the physical time, the systems at each order can be summed up to write

    .. math::
        \begin{cases}
        \dfrac{\textrm{d}}{dt} a_0(t) & = f_{a_0}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        a_0 \dfrac{\textrm{d}}{dt} \beta_0(t) & = f_{\beta_0}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        & \vdots \\
        \dfrac{\textrm{d}}{dt} a_{N-1}(t) & = f_{a_{N-1}}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        a_{N-1} \dfrac{\textrm{d}}{dt} \beta_{N-1}(t) & = f_{\beta_{N-1}}(\boldsymbol{a}, \boldsymbol{\beta}), 
        \end{cases}

    where

    .. math::
        \begin{cases}
        f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}) & = \sum_{j=1}^{N_e} \epsilon^{j} f_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
        f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}) & = \sum_{j=1}^{N_e} \epsilon^{j} f_{\beta_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}).
        \end{cases}

    The MMS system then obtained represents a nonlinear autonomous system of :math:`2N` coupled :math:`1^{\textrm{st}}` order PDEs, with the physical time :math:`t` as the independent variable. 
    Like the time scales-dependent evolution equations, they can be solved numerically or analytically, if analytical solutions exist. 

    These physical time-dependent evolution equations are also computed in :func:`evolution_equations`.
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
        self.xMMS, sub_xMMS_t, sub_x = self.asymptotic_series(dynamical_system, eps_pow_0=self.eps_pow_0)
        
        # Substitutions required
        self.sub = Substitutions_MMS(sub_t, sub_xMMS_t, sub_x, sub_scaling, sub_omega, sub_sigma)
    
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
        self.compute_EqMMS(dynamical_system)
        
    def time_scales(self):
        r"""
        Define the time scales.

        Notes
        -----
        The time scales are defined as (see :class:`~MMS.MMS.Multiple_scales_system`)
        
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
        The series expansion for oscillator :math:`i` (and for a leading order term :math:`\epsilon^0 = 1`) takes the form (see :class:`~MMS.MMS.Multiple_scales_system`)

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
        sub_t, sub_x, sub_xMMS_t = list(map(self.sub.__dict__.get,["sub_t", "sub_x", "sub_xMMS_t"]))

        for ix in range(self.ndof):
            fF.append( (dynamical_system.forcing.fF[ix].subs(self.sub.sub_scaling)
                        .subs(sub_x).doit().subs(sub_xMMS_t).expand().subs(sub_t).doit())
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
        
        
    def compute_EqMMS(self, dynamical_system):
        r"""
        Compute the system of equations for each oscillator at each order of :math:`\epsilon`. This system is described in :class:`~MMS.MMS.Multiple_scales_system`.

        The output `EqMMS` is a list of lists:

        - The :math:`1^{\text{st}}` level lists are associated to the equations for each oscillator,
        
        - The :math:`2^{\text{nd}}` level lists are associated to the orders of :math:`\epsilon` from the lowest to the highest order.

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
        Apply the MMS. 

        Parameters
        ----------
        rewrite_polar : str, int or list of int, optional
            The orders at which the solutions will be rewritten in polar form.
            See :func:`sol_xMMS_polar`.
        
        Notes
        -----
        The application of the MMS is operated by successively calling the following methods.

        #. :func:`system_t0`: An equivalent system is written in terms of the fast time scale :math:`t_0`. 
           This introduces the temporary unknowns :math:`\tilde{x}_{i,j}(t_0)` and allows the use of :func:`~sympy.solvers.ode.dsolve`.

        #. :func:`sol_order_0`: Leading order solutions :math:`x_{i,0}(\boldsymbol{t})` are defined.

        #. :func:`secular_analysis`: The leading order solutions are introduced in the equations and the secular terms at each order are identified. 
           Cancelling those secular terms is a condition for bounded solutions. 
           It leads to a system of equations governing the slow evolution of the complex amplitude of the homogeneous leading order solutions. 
           Each equation takes the form 
           
           .. math::
            \textrm{D}_{j} A_i(\boldsymbol{t}_s) = f_{A_i}^{(j)}(\boldsymbol{A}, t_1).

           After cancelling the secular terms the higher order equations are solved successively to express the higher order solutions :math:`x_{i,j}(\boldsymbol{t}),\; j>0` in terms of the leading order ones.

        #. :func:`autonomous_phases`: The phase coordinates are changed from :math:`\phi_i(\boldsymbol{t}_s)` to :math:`\beta_i(\boldsymbol{t}_s)` to cancel the slow time :math:`t_1` in the secular terms. This will be used afterwards to obtain an autonomous system.

        #. :func:`evolution_equations`: The secular conditions are split into real and imaginary parts, polar coordinates are used and the autonomous phases are introduced,
           resulting in an autonomous system of evolution equations on polar coordinates. 
           Equations come by two, one representing the amplitude evolution while the other represents the phase's, such that
           
           .. math::
            \begin{cases}
            \textrm{D}_{j} a_i(\boldsymbol{t}_s) & = f_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \textrm{D}_{j} \beta_i(\boldsymbol{t}_s) & = f_{\beta_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\beta}).
            \end{cases}

           This is the key result of the MMS. The evolution on each time scale are combined to reintroduce the physical time, resulting in a system of the form

           .. math::
            \begin{cases}
            \dfrac{\textrm{d}}{dt} a_i(t) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}), \\
            a_i \dfrac{\textrm{d}}{dt} \beta_i(t) & = f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}).
            \end{cases}

        #. :func:`sol_xMMS_polar`: The leading and higher order solutions are rewritten in terms of polar coordinates using :math:`\cos` and :math:`\sin` functions.
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
        Rewrite the system with the fast time scale as the only independent variable.

        Notes
        -----
        This is a trick to use :func:`~sympy.solvers.ode.dsolve`, which only accepts functions of 1 variable, to solve higher order equations. 
        The higher order equations are rewritten in terms of temporary coordinates :math:`\tilde{x}_{i,j}(t_0)` in place of :math:`x_{i,j}(\boldsymbol{t})`, with :math:`i,j` denoting the oscillator number and :math:`\epsilon` order, respectively. 
        This is equivalent to temporary considering that :math:`\boldsymbol{A}(\boldsymbol{t}_s)` does not depend on the slow times, which is of no consequence as there are no slow time derivatives appearing in the higher order equations at this stage. 
        Indeed, they were either substituted using the complex evolution equations, or they disappeared when eliminating the secular terms. 
        
        """
        
        xMMS_t0  = [] # t0-dependent variables xij(t0). Higher time scales dependency is ignored.
        EqMMS_t0 = [] # Equations at each order with only t0 as an explicit variable. Leads to a harmonic oscillator at each order with a t0-periodic forcing coming from lower order solutions.
        
        for ix in range(self.ndof):
            xMMS_t0 .append([ Function(r'\tilde{x_'+'{{{},{}}}'.format(ix,io)+'}', real=True)(self.tS[0]) for io in range(0, 1+self.Ne) ]) 
            EqMMS_t0.append([ self.EqMMS[ix][0].subs(self.xMMS[ix][0], xMMS_t0[ix][0]).doit() ])
            
        self.EqMMS_t0 = EqMMS_t0
        self.xMMS_t0  = xMMS_t0
        
        
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
        Identify the secular terms in the MMS equations. 
        This allows to:

        1. Compute the evolution equations of the complex amplitudes :math:`A_i(\boldsymbol{t}_s)`, coming from the elimination of the secular terms,
        
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
        Compute the higher order solutions :math:`x_{i,j}(\boldsymbol{t}_s),\; j>0`.

        Parameters
        ----------
        EqMMS_t0 : list of list of sympy.Expr
            The MMS equations at each order and for each oscillator written with :math:`t_0` as the only independent variable. 
        xMMS_t0 : list of list of sympy.Function
            Oscillators' response at each order written in terms of :math:`t_0` only, :math:`\tilde{x}_{i,j}(t_0)`.
        io : int
            The order of :math:`\epsilon`.
        ix : int
            The dof number.
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
        See details on this choice in :class:`~MMS.MMS.Multiple_scales_system`.
        """
        
        self.coord.beta   = [ Function(r'\beta_{}'.format(ix), real=True)(*self.tS[1:])                                            for ix in range(self.ndof) ]
        def_beta          = [ Rational(self.ratio_omega_osc[ix], self.ratio_omegaMMS) * self.sigma*self.tS[1] - self.coord.phi[ix] for ix in range(self.ndof) ]
        def_phi           = [ solve(def_beta[ix]-self.coord.beta[ix], self.coord.phi[ix])[0]                                       for ix in range(self.ndof) ]
        self.sub.sub_phi  = [ (self.coord.phi[ix], def_phi[ix])                                                                    for ix in range(self.ndof) ]
        self.sub.sub_beta = [ (self.coord.beta[ix], def_beta[ix])                                                                  for ix in range(self.ndof) ]

    def evolution_equations(self):
        r"""
        Derive the evolution equations of the polar coordinates system.
        
        Notes
        -----
        Derive the evolution equations of the polar coordinates system (defined in :func:`polar_coordinates` and :func:`autonomous_phases`) from the secular conditions. For oscillator :math:`i`, these are defined as
        
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

    
    def find_harmonics(self):
        """
        Determine the harmonics contained in the MMS solutions. 

        Returns
        -------
        harmonics: list
            list of the harmonics appearing in the MMS solutions.
        """
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
            sub_SS .extend( [(mms.coord.a[ix] , a[ix]), (mms.coord.beta[ix], beta[ix])] )
    
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

    def stability_analysis(self, coord="cartesian", rewrite_polar=False, eigenvalues=False, bifurcation_curves=False, analyse_blocks=False, kwargs_bif=dict()):
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
                        eigvalsA = self.eigenvalues(A, detA=detA, trA=trA)
                        if coord=="cartesian":
                            eigvalsA = [cartesian_to_polar(eigval, self.sub.sub_polar, sub_phase=self.sub.sub_phase) for eigval in eigvalsA]
                        self.stab.blocks_eigvals.append(eigvalsA)

                    if bifurcation_curves:
                        bif_aA, bif_sigA = self.bifurcation_curves(detA, trA, **kwargs_bif)
                        self.stab.blocks_bif_a.append(bif_aA)
                        self.stab.blocks_bif_sig.append(bif_sigA)

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
            
    def bifurcation_curves(self, detJ, trJ, var_a=False, var_sig=True, solver=sfun.solve_poly2):
        r"""
        Compute bifurcation curves.

        Parameters
        ----------
        detJ: sympy.Expr
            The determinant of the matrix.
        trJ: sympy.Expr
            The trace of the matrix.
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
        The bifurcation curves computed here are the curves defined by the bifurcation points obtained for any forcing frequency and amplitude. 
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
            The amplitude plot :math:`a(F)`.
        fig2 : Figure
            The phase plot :math:`\beta(F)`.
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
    The coordinates used in the steady state analysis.
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
    Stability analysis parameters and outputs.
    """                
    
    def __init__(self):
        pass  
    
#%% Chain rule functions written for the MMS
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

#%% Numpy transforms and plot functions


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
    r"""
    Numpise the backbone curve's frequency :math:`\omega_{\textrm{bbc}}`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega_bbc: numpy.ndarray
        Numpised backbone curve's frequency.
    """
    omega_bbc  = sfun.sympy_to_numpy(rescale(ss.sol.omega_bbc, mms), param)
    return omega_bbc

def numpise_omega_FRC(mms, ss, param):
    r"""
    Numpise the forced response's frequency :math:`\omega`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega: numpy.ndarray
        Numpised forced response's frequency.
    """
    omega = [np.real(sfun.sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.sol.sigma]
    return omega

def numpise_omega_bif(mms, ss, param):
    r"""
    Numpise the bifurcation curves' frequency :math:`\omega_{\textrm{bif}}`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega_bif: list of numpy.ndarray
        Numpised bifurcation curves' frequency.
    """
    omega_bif = [np.real(sfun.sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.stab.bif_sigma]
    return omega_bif

def numpise_phase(mms, ss, dyn, param, omega, F):
    r"""
    Numpise the phase :math:`\beta_i`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    dyn: Dynamical_system
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.
    omega: numpy.ndarray or list of numpy.ndarray
        The frequency array.
    F: numpy.ndarray
        The forcing amplitude array.

    Returns
    -------
    phase: list of numpy.ndarray
        Numpised phase.
    """
    
    if not isinstance(omega,list):
        omega = [omega]
        
    phase = []
    
    for omegai in omega:
        param_phase = param | dict(omega=(mms.omega, omegai), F=(dyn.forcing.F, F))
        sin_phase = sfun.sympy_to_numpy( rescale(ss.sol.sin_phase[1], mms), param_phase )
        cos_phase = sfun.sympy_to_numpy( rescale(ss.sol.cos_phase[1], mms), param_phase )
        phase.append(np.arctan2(sin_phase, cos_phase))

    return phase

def numpise_F_ARC(mms, ss, param):
    r"""
    Numpise the forced response's forcing amplitude :math:`F`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    F: numpy.ndarray
        Numpised forced response's forcing amplitude.
    """
    F = sfun.sympy_to_numpy(rescale(mms.eps**mms.forcing.f_order * ss.sol.F, mms), param)
    return F