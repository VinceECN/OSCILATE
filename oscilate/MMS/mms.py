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

    This operation is performed using :func:`compute_EqO`.

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
    
    All solutions previously computed using the complex amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)` can be rewritten in terms of the polar coordinates :math:`\boldsymbol{a}(\boldsymbol{t}_s),\; \boldsymbol{\beta}(\boldsymbol{t}_s)` using :func:`sol_x_polar`. 

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
        Compute the system of equations for each oscillator at each order of :math:`\epsilon`. This system is described in :class:`~MMS.MMS.Multiple_scales_system`.

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

        # Derive the evolution equations
        self.evolution_equations()

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
        Indeed, they were either substituted using the complex evolution equations, or they disappeared when eliminating the secular terms. 
        
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