Application of the MMS
----------------------

Asymptotic series and time scales
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

The introduction of asymptotic series and time scales are performed using :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.asymptotic_series` and :meth:`~oscilate.MMS.mms.Multiple_scales_system.time_scales`.

Scaling
~~~~~~~

The construction of the MMS system requires a scaling of the parameters. Most scalings are already passed to the MMS through the `sub_scaling` parameter. 
However, the natural frequencies also need to be scaled as they can contain both a leading order term and a detuning term.
Natural frequencies :math:`\omega_i` are defined as a function of the reference frequency :math:`\omega_{\textrm{ref}}` through the `ratio_omega_osc` optional parameter, which is then used in :meth:`~oscilate.MMS.mms.Multiple_scales_system.oscillators_frequencies`.
This allows to define internal resonance relations among the oscillators. 
If these internal resonances are not perfect, detunings can be introduced through the `detunings` optional parameter, which needs to be scaled and part of the `sub_scaling` parameter.

To write the MMS system it is convenient to introduce the leading order natural frequencies

.. math::
    \omega_{i,0} = r_i \omega_{\textrm{ref}},

where :math:`r_i` stands for ``ratio_omega_osc[i]``.

The multiple scales system
~~~~~~~~~~~~~~~~~~~~~~~~~~

Introducing the asymptotic series, the time scales and the scaled parameters in the initial dynamical system (see :class:`~oscilate.MMS.dyn_sys.Dynamical_system`) results in :math:`N_e+1` dynamical systems, each one appearing at different orders of :math:`\epsilon`.
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
    \textrm{D}_0^2 x_{0,0} + \omega_{0,0}^2 x_{0,0} & = f_{0}^{(0)}(t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0^2 x_{N-1,0} + \omega_{N-1,0}^2 x_{N-1,0} & = f_{N-1}^{(0)}(t_0, t_1),
    \end{cases} \\[15pt]
    & \epsilon^1 \rightarrow \;
    \begin{cases}
    \textrm{D}_0^2 x_{0,1} + \omega_{0,0}^2 x_{0,1} & = f_{0}^{(1)}  (\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0^2 x_{N-1,1} + \omega_{N-1,0}^2 x_{N-1,1} & = f_{N-1}^{(1)}(\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1),
    \end{cases} \\[15pt]
    & \epsilon^2 \rightarrow \;
    \begin{cases}
    \textrm{D}_0^2 x_{0,2} + \omega_{0,0}^2 x_{0,2} & = f_{0}^{(2)}  (\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_2 \boldsymbol{x}_0, \boldsymbol{x}_1, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_1, t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0^2 x_{N-1,2} + \omega_{N-1,0}^2 x_{N-1,2} & = f_{N-1}^{(2)}(\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_2 \boldsymbol{x}_0, \boldsymbol{x}_1, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_1, t_0, t_1),
    \end{cases} \\[10pt]
    & \hspace{3cm} \vdots \\[10pt]
    & \epsilon^{N_e} \rightarrow \;
    \begin{cases}
    \textrm{D}_0^2 x_{0,N_e} + \omega_{0,0}^2 x_{0,N_e} & = f_{0}^{(N_e)}  (\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_{N_e} \boldsymbol{x}_0, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_{N_e-1}, t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0^2 x_{N-1,N_e} + \omega_{N-1,0}^2 x_{N-1,N_e} & = f_{N-1}^{(N_e)}(\boldsymbol{x}_0, \cdots, \textrm{D}_0 \textrm{D}_{N_e} \boldsymbol{x}_0, \cdots, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_{N_e-1}, t_0, t_1).
    \end{cases}
    \end{aligned}

Consider oscillator :math:`i` at order :math:`j`:

- The left-hand side term represents a harmonic oscillator of frequency :math:`\omega_{i,0}` oscillating with respect to the fast time :math:`t_0`.
- The right-hand side term :math:`f_{i}^{(j)}` is analogous to a forcing generated by all combinations of terms that appear on oscillator :math:`i`'s equation at order :math:`\epsilon^j`. This can involve lower order terms :math:`x_{i,\ell}, \; \ell \leq j`, coupling terms :math:`x_{k, \ell}, \; k \neq j,\; \ell \leq j`, their derivatives and cross-derivatives with respect to the time scales, and physical forcing terms. The later is responsible for the dependency on :math:`t_0,\; t_1`. The reason why slower time scales are not involved will be explained in the following.

Function :math:`f_{i}^{(j)}` tends to get increasingly complex as the order increases because the initial equations generate more high order terms than low order ones.

This operation is performed using :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.compute_EqO`.

Note that internal resonance relations can be given through the `ratio_omega_osc` optional parameter, which is then used in :meth:`~oscilate.MMS.mms.Multiple_scales_system.oscillators_frequencies`.

Frequency of interest
~~~~~~~~~~~~~~~~~~~~~

The response of :math:`x_i` will be analysed at a frequency :math:`\omega`, defined as

.. math::
    \omega = \omega_{\textrm{MMS}} + \epsilon \sigma,

where :math:`\omega_{\textrm{MMS}}` is the **central MMS frequency**, controlled through the `ratio_omegaMMS` optional parameter and expressed in terms of `omega_ref`, 
and :math:`\sigma` is a detuning about that frequency.
In case the forced response is studied, :math:`\omega` corresponds to the forcing frequency.
In case the free response is studied, :math:`\omega` corresponds to the frequency of free oscillations, which generates the backbone curve of the forced response.
Note that :math:`\omega t = \omega_{\textrm{MMS}} t_0 + \sigma t_1`. This is the reason why the forcing only involves these two time scales in the right-hand side functions of the MMS system.


Iteratively solving the MMS system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The multiple scales system can be solved iteratively by solving successively the systems of equations at each order.

Leading order solution
^^^^^^^^^^^^^^^^^^^^^^

The leading order solution for oscillator :math:`i` must satisfy

.. math::
    \textrm{D}_0^2 x_{i,0} + \omega_{i,0}^2 x_{i,0} = f_{i}^{(0)}(t_0, t_1).

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

The leading order solutions are defined in :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.sol_order_0`.

Higher order solutions
^^^^^^^^^^^^^^^^^^^^^^

Once the leading order solutions are computed, they can be injected in the :math:`1^\textrm{st}` higher order equations, where they (and their derivatives) appear as *forcing terms*, potentially together with physical forcing.
The :math:`1^\textrm{st}` higher order equation for oscillator :math:`i` is 

.. math::
    \textrm{D}_0 x_{i,1} + \omega_{i,0}^2 x_{i,1} = f_{i}^{(1)}(\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1).

The forcing terms that involve oscillations at :math:`\omega_{i,0}` would force the oscillator on its natural frequency. Moreover, damping is always weak in the MMS, so damping terms of the form 
:math:`c \textrm{D}_0 x_{i,1}` do not appear at this order. 
The aforementioned forcing terms would thus lead to unbounded solutions, which is unphysical. 
These forcing terms, called **secular terms**, must therefore be eliminated. 
For instance, the :math:`1^\textrm{st}` higher order equation for oscillator :math:`i` with the secular terms cancelled is

.. math::
    \textrm{D}_0 x_{i,1} + \omega_{i,0}^2 x_{i,1} = \breve{f}_{i}^{(1)}(\boldsymbol{x}_0, \textrm{D}_0 \boldsymbol{x}_0, \textrm{D}_0^2 \boldsymbol{x}_0, \textrm{D}_1 \boldsymbol{x}_0, \textrm{D}_0\textrm{D}_1 \boldsymbol{x}_0, t_0, t_1),

where :math:`\breve{f}_{i}^{(1)}` is :math:`f_{i}^{(1)}` with the secular terms cancelled, i.e. without terms oscillating as :math:`\omega_{i,0}`.
After cancelation of the secular terms, each oscillator's equation can be solved as a forced harmonic oscillator with the independent variable :math:`t_0`.

Note that only the particular solutions are considered when solving higher order terms, i.e.

.. math::
    x_{i,1}(\boldsymbol{t}) = \underbrace{x_{i,1}^\textrm{h}(\boldsymbol{t})}_{=0} + x_{i,1}^\textrm{p}(t_0, t_1).

This choice can be justified if one assumes that initial conditions are of leading order. Though this is questionable, it is assumed here. A more advanced treatment of the initial conditions of higher order solutions is detailed in :cite:`nayfehResolvingControversiesApplication2005` and section :ref:`mms_complex`.

The higher order solutions :math:`x_{i,1}(\boldsymbol{t})` are expressed as a function of the leading order unknown amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)`, 
their slow time derivatives :math:`\textrm{D}_i\boldsymbol{A}(\boldsymbol{t}_s), \; i=1, ..., N_e`, and forcing terms if any (including the hard forcing amplitudes :math:`\boldsymbol{B}`).  

This process is repeated successively at each order, i.e. the computed solutions are introduced in the next higher order system of equations, 
secular terms are cancelled and the next higher order solutions are computed. 


The secular terms are identified in :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.secular_analysis` and the leading order solutions are computed in :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.sol_higher_order`. 
Note that :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.sol_higher_order` is applied on equations with only :math:`t_0` as the independent variable so as to allow the use of :func:`~sympy.solvers.ode.dsolve`. 
This is enforced using :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.system_t0`, which temporarily ignores the dependency of :math:`\boldsymbol{A}(\boldsymbol{t}_s)` on the slow time scales.

Secular analysis
^^^^^^^^^^^^^^^^

At this stage, the solutions are all expressed in terms of the unknown amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)` and their slow time derivatives :math:`\textrm{D}_1 A_i(\boldsymbol{t}_s),\; \cdots,\; \textrm{D}_{N_e} A_i(\boldsymbol{t}_s)`. 
These can be obtained from the elimination of the secular terms (called **secular conditions**), as described below. 

The :math:`i^\textrm{th}` MMS equation at :math:`1^\textrm{st}` higher order involves the slow time derivative :math:`\textrm{D}_1 A_i(\boldsymbol{t}_s)`, which appears 
in the secular term. It is coming from the chain rule 

.. math::
    \begin{cases}
    \dfrac{\textrm{d} x_i(t)}{\textrm{d}t} & = \dfrac{\partial x_{i,0}(\boldsymbol{t})}{\partial t_0} + \epsilon \dfrac{\partial x_{i,1}(\boldsymbol{t})}{\partial t_0} + \epsilon \dfrac{\partial x_{i,0}(\boldsymbol{t})}{\partial t_1} + \mathcal{O}(\epsilon^2),
    \\
    \dfrac{\textrm{d}^2 x_i(t)}{\textrm{d}t^2} & = \dfrac{\partial^2 x_{i,0}(\boldsymbol{t})}{\partial t_0^2} + \epsilon \dfrac{\partial^2 x_{i,1}(\boldsymbol{t})}{\partial t_0^2} + 2 \epsilon \dfrac{\partial^2 x_{i,0}(\boldsymbol{t})}{\partial t_0 \partial t_1} + \mathcal{O}(\epsilon^2).
    \end{cases}

In addition, the :math:`\textrm{D}_1 A_j(\boldsymbol{t}_s),\; j\neq i` do not appear in the :math:`i^\textrm{th}` MMS equation as couplings among oscillators are weak. 
It is thus possible to use the secular conditions in order to express the :math:`\textrm{D}_1 A_i(\boldsymbol{t}_s)` as a function of :math:`\boldsymbol{A}(\boldsymbol{t}_s)`. 

This process can be done successively at each order to obtain the system of complex modulation (or slow time) equations

.. math::
    \begin{cases}
    \textrm{D}_1 A_i(\boldsymbol{t}_s) & = f_{A_i}^{(1)}(\boldsymbol{A}, t_1), \\
    & \vdots \\
    \textrm{D}_{N_e} A_i(\boldsymbol{t}_s) & = f_{A_i}^{(N_e)}(\boldsymbol{A}, t_1).
    \end{cases}

:math:`f_{A_i}^{(j)}(\boldsymbol{A}, t_1)` are functions governing the modulation of :math:`A_i` with respect to the slow time :math:`t_j`. 
Note the dependency of :math:`f_{A_i}^{(j)}` on :math:`t_1` due to the possible presence of forcing. The complex modulation equations are derived in :meth:`~oscilate.MMS.mms_oscillator.Multiple_scales_oscillator.secular_analysis`.

The above system of :math:`1^\textrm{st}` order PDE can theoretically be solved to obtain the complex amplitudes :math:`\boldsymbol{A}`. 
However, this approach is not the prefered one as 

- It is more convenient to deal with real variables than complex ones to get a physical meaning from the analysis,

- It is more convenient to deal with autonomous systems (without the explicit :math:`t_1`-dependency) than nonautonomous ones,

- The PDEs are complex.

The first two points can be achieved introducing new coordinates, as described thereafter.

Modulation equations
~~~~~~~~~~~~~~~~~~~~

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

the modulation equations on the :math:`A_i(\boldsymbol{t}_s)` can be split into real and imaginary terms, leading to

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

it is convenient to pre-multiply the modulation equations on the :math:`A_i(\boldsymbol{t}_s)` by :math:`e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)}` or even :math:`\gamma e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)}` with, for instance, :math:`\gamma = 2`. 
This avoids the presence of :math:`\cos(\phi_i)`, :math:`\sin(\phi_i)` and many :math:`1/2` terms in the modulation equations. 
The modulation functions on polar coordinates are therefore defined as

.. math::
    \begin{cases}
    \hat{f}_{a_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1) & = \Re\left[ 2 e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)} f_{A_i}^{(j)}(\boldsymbol{A}, t_1) \right], \\
    \hat{f}_{\phi_i}^{(j)}(\boldsymbol{a}, \boldsymbol{\phi}, t_1) & = \Im\left[ 2 e^{-\textrm{j} \phi_i(\boldsymbol{t}_s)} f_{A_i}^{(j)}(\boldsymbol{A}, t_1) \right].
    \end{cases}

The modulation equations system on the polar coordinates involves only real variables, but functions :math:`\hat{f}_{a_i}^{(j)}` and :math:`\hat{f}_{\phi_i}^{(j)}` are
still nonautonomous due to the explicit dependency on :math:`t_1`. 

The polar coordinates are introduced in :meth:`~oscilate.MMS.mms.Multiple_scales_system.polar_coordinates`. The real modulation equations are only computed for the autonomous system, as described below. 

Autonomous phase coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The presence of nonautonomous terms stems from forcing terms, which involve :math:`\cos(\sigma t_1 - \phi_i), \; \sin(\sigma t_1 - \phi_i)` in the polar modulation functions of oscillator :math:`i`. 
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

In addition, and as discussed previously, the introduction of these new phase coordinates removes the explicit dependency of the modulation functions on :math:`t_1`, which was due to terms :math:`\cos(\sigma t_1 - \phi_i), \; \sin(\sigma t_1 - \phi_i)`.
The forcing phase being zero (i.e. reference phase), the :math:`\beta_i(\boldsymbol{t}_s)` can be seen as the phases relative to the forcing. 

Introducing the notation

.. math::
    \boldsymbol{\beta}(\boldsymbol{t}_s)^\intercal = [\beta_1(\boldsymbol{t}_s), \beta_2(\boldsymbol{t}_s), \dots, \beta_{N_e}(\boldsymbol{t}_s)],

the modulation equations can be rewritten as

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

The above system is the key result of the application of the MMS as it governs the modulation of leading order amplitudes and phases, on which all higher order solutions depend. 
It can be solved numerically or analytically, if analytical solutions exist, though it is generally not the case. An analytical resolution is shown in :ref:`example_VdP`.
The modulation system can also be rewritten in a more compact form as discussed in the following.

The autonomous phase coordinates are introduced in :meth:`~oscilate.MMS.mms.Multiple_scales_system.autonomous_phases` and the modulation equations are computed in :meth:`~oscilate.MMS.mms.Multiple_scales_system.modulation_equations`.

All solutions previously computed using the complex amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)` can be rewritten in terms of the polar coordinates :math:`\boldsymbol{a}(\boldsymbol{t}_s),\; \boldsymbol{\beta}(\boldsymbol{t}_s)` using :meth:`~oscilate.MMS.mms.Multiple_scales_system.sol_x_polar`. 

Reconstituted equations
^^^^^^^^^^^^^^^^^^^^^^^

The reconstitution method consists in combining the modulation equations at each slow time scales in order to obtain modulation equations in the physical time. The approach used in the ``oscilate`` package corresponds to the *Complete Inconsistent Method* discussed in :cite:`luongoReconstitutionProblemMultiple1999a`, in opposition to consistent approaches which treat each order individually. 

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
Like the time scales-dependent modulation equations, they can be solved numerically or analytically, if analytical solutions exist. 

These physical time-dependent modulation equations are also computed in :meth:`~oscilate.MMS.mms.Multiple_scales_system.modulation_equations`.

Note that after reintroduction of the physical time, the leading order, homogeneous solution for oscillator :math:`i` takes the form 

.. math::
    x^{\textrm{h}}_{i,0}(t) = a_i(t) \cos\left( \frac{r_i}{r_{\textrm{MMS}}} \omega t - \beta_i(t)\right).


.. _mms_complex:

Complex form of the MMS
~~~~~~~~~~~~~~~~~~~~~~~

The application of the MMS to the system's equations in oscillator form was described above. Alternatively, the MMS can be applied to the system's equations in complex form (see section :ref:`dyn_sys_complex`). The procedure is very similar, but the fact that the complex variables :math:`\boldsymbol{z}` embed the state space :math:`(\boldsymbol{x}, \dot{\boldsymbol{x}})` allows to extract more information from the system's dynamics, resulting in more accurate solutions and modulation equations. The differences between the oscillator and complex forms of the equations are stressed thereafter.

The multiple scales system
^^^^^^^^^^^^^^^^^^^^^^^^^^

Series expansions of the :math:`z_i(t)` variables are performed such that

.. math::
    z_i(t) = z_{i,0}(t) + \epsilon z_{i,1}(t) + \epsilon^2 z_{i,2}(t) + \cdots + \epsilon^{N_e} z_{i,N_e}(t) + \mathcal{O}(\epsilon^{N_e+1}),

time scales are introduced in the same way as before, the chain rule is used (only the :math:`1^{\text{st}}` order time derivatives are needed), parameters are scaled, oscillators' natural frequencies and the frequency of interest are expressed in terms of a reference frequency, possibly with detunings. This leads to the complex MMS system

.. math::
    \begin{aligned}
    & \epsilon^0 \rightarrow \;
    \begin{cases}
    \textrm{D}_0 z_{0,0} & = \textrm{j} \omega_{0,0} z_{0,0} + g_{0}^{(0)}(t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0 z_{N-1,0} & = \textrm{j} \omega_{N-1,0} z_{N-1,0} + g_{N-1}^{(0)}(t_0, t_1),
    \end{cases} \\[15pt]
    & \epsilon^1 \rightarrow \;
    \begin{cases}
    \textrm{D}_0 z_{0,1} & = \textrm{j} \omega_{0,0} z_{0,1} + g_{0}^{(1)}  (\boldsymbol{z}_0, \bar{\boldsymbol{z}}_0, \textrm{D}_0 \boldsymbol{z}_0, \textrm{D}_1 \boldsymbol{z}_0, t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0 z_{N-1,1} & = \textrm{j} \omega_{0,0} z_{N-1,1} + g_{N-1}^{(1)}(\boldsymbol{z}_0, \bar{\boldsymbol{z}}_0, \textrm{D}_0 \boldsymbol{z}_0, \textrm{D}_1 \boldsymbol{z}_0, t_0, t_1),
    \end{cases} \\[10pt]
    & \hspace{3cm} \vdots \\[10pt]
    & \epsilon^{N_e} \rightarrow \;
    \begin{cases}
    \textrm{D}_0 z_{0,N_e} & = \textrm{j} \omega_{0,0} z_{0,N_e} + g_{0}^{(N_e)}  (\boldsymbol{z}_0, \cdots, \textrm{D}_{N_e} \boldsymbol{z}_0, \cdots, \textrm{D}_1 \boldsymbol{z}_{N_e-1}, t_0, t_1), \\
    & \vdots \\
    \textrm{D}_0 z_{N-1,N_e} & = \textrm{j} \omega_{N-1,0} z_{N-1,N_e} + g_{N-1}^{(N_e)}(\boldsymbol{z}_0, \cdots, \textrm{D}_{N_e} \boldsymbol{z}_0, \cdots, \textrm{D}_1 \boldsymbol{z}_{N_e-1}, t_0, t_1).
    \end{cases}
    \end{aligned}

Like with the oscillator form, these equations represent harmonic oscillators of natural frequencies :math:`+ \omega_{i,0}` in the fast time, and subjected to forcing-like terms. 

Leading order solution
^^^^^^^^^^^^^^^^^^^^^^

In the same way as with the oscillator form, the leading order solution for oscillator :math:`i` satisfies

.. math::
    \textrm{D}_0 z_{i,0} = \textrm{j} \omega_{i,0} z_{i,0} + g_{i}^{(0)}(t_0, t_1),

such that it takes the form

.. math::
    z_{i,0}(\boldsymbol{t}) = z_{i,0}^\textrm{h}(\boldsymbol{t}) + z_{i,0}^\textrm{p}(t_0, t_1).

These homogeneous and particular solutions can be expressed

.. math::
    \begin{cases}
    z_{i,0}^\textrm{h}(\boldsymbol{t}) & = A_i(\boldsymbol{t}_s) e^{\textrm{j} \omega_{i,0} t_0}, 
    \\
    z_{i,0}^\textrm{p}(t_0, t_1) & = B_{i,+} e^{\textrm{j} \omega t} + B_{i,-} e^{-\textrm{j} \omega t} = B_{i,+} e^{\textrm{j} (\omega_{\textrm{MMS}} t_0 + \sigma t_1)} + B_{i,-} e^{-\textrm{j} (\omega_{\textrm{MMS}} t_0 + \sigma t_1)},
    \end{cases}

where :math:`A_i` is a slow time-dependent complex amplitude to be determined while :math:`B_{i,\pm}` are time-independent functions of the forcing parameters. Note that in most situations, forcing does not appear at leading order (i.e. forcing is weak), so :math:`B_{i,\pm}=0`.

In the following it will be convenient to use the notations

.. math::
    \begin{split}
    \boldsymbol{A}(\boldsymbol{t}_s)^\intercal & = [A_0(\boldsymbol{t}_s), A_1(\boldsymbol{t}_s), \cdots, A_{N-1}(\boldsymbol{t}_s)], \\
    \boldsymbol{B}_\pm^\intercal & = [B_{0,\pm}, B_{1,\pm}, \cdots, B_{N-1,\pm}], \\
    \boldsymbol{B}^\intercal & = [\boldsymbol{B}_+^\intercal, \boldsymbol{B}_-^\intercal].
    \end{split}

The leading order solutions are defined in :meth:`~oscilate.MMS.mms_complex.Multiple_scales_complex.sol_order_0`.

Higher order solutions
^^^^^^^^^^^^^^^^^^^^^^

The higher order solutions can be derived in the same way as with the oscillator form, noticing that the terms in :math:`e^{+\textrm{j}\omega_{i,0} t_0}` in the equation of oscillator :math:`i` are secular and must therefore be eliminated. The :math:`1^{\textrm{st}}` higher order solution of oscillator :math:`i` is hence governed by the equation

.. math::
    \textrm{D}_0 z_{i,1} = \textrm{j} \omega_{i,0} z_{i,1} + \breve{g}_{i}^{(1)}(\boldsymbol{z}_0, \bar{\boldsymbol{z}}_0, \textrm{D}_0 \boldsymbol{z}_0, \textrm{D}_1 \boldsymbol{z}_0, t_0, t_1),

where :math:`\breve{g}_{i}^{(1)}` is :math:`g_{i}^{(1)}` with the secular terms cancelled. 

It is interesting to stress that this cancellation of secular terms is equivalent to fulfilling the condition that the particular solution of a linear problem must be orthogonal to the kernel of the adjoint problem. To apply this, it is convenient to introduce the differential operator :math:`\mathscr{L}_i` such that

.. math::
    \mathscr{L}_i z_{i,1} =  (\textrm{D}_0 - \textrm{j} \omega_{i,0}) z_{i,1}.

This way, the initial problem writes

.. math::
    \mathscr{L}_i z_{i,1} = \breve{g}_{i}^{(1)},

and the associated adjoint problem is

.. math::
    \mathscr{L}^{\dagger}_i \breve{g}_{i}^{(1)} = z_{i,1},

where the adjoint operator :math:`\mathscr{L}^{\dagger}_i` is defined such that

.. math::
    \langle \breve{g}_{i}^{(1)}, \mathscr{L}_i z_{i,1} \rangle = \langle \mathscr{L}^{\dagger}_i \breve{g}_{i}^{(1)}, z_{i,1} \rangle.

Given the :math:`\omega_{i,0}`-periodicity of the homogeneous solution, it is relevant to use the metric associated to the :math:`L^2` inner product

.. math::
    \langle f, g \rangle = \frac{\omega_{i,0}}{2\pi} \int_0^{2 \pi / \omega_{i,0}} \bar{f} g \textrm{d} t_0.

From there, one can show that the operator is anti self-adjoint, i.e. :math:`\mathscr{L}^{\dagger}_i = - \mathscr{L}_i`, and the adjoint kernel is 

.. math::
    \textrm{ker}(\mathscr{L}^{\dagger}_i) = \{ e^{\textrm{j} \omega_{i,0} t_0} \}.

The condition to fulfil is therefore

.. math::
    \langle e^{\textrm{j} \omega_{i,0} t_0}, {g}_{i}^{(1)} \rangle = \frac{\omega_{i,0}}{2\pi} \int_0^{2 \pi / \omega_{i,0}} e^{-\textrm{j} \omega_{i,0} t_0} {g}_{i}^{(1)} \textrm{d} t_0 = 0,

such that the adjoint kernel orthogonality condition is equivalent to imposing a zero :math:`e^{\textrm{j} \omega_{i,0} t_0}` component on :math:`{g}_{i}^{(1)}`, i.e. cancelling the secular terms. 

After the cancellation of secular terms, the higher order solutions :math:`z_{i,1}(\boldsymbol{t})` can be expressed as a function of the leading order unknown amplitudes :math:`\boldsymbol{A}(\boldsymbol{t}_s)`, 
their slow time derivatives :math:`\textrm{D}_i\boldsymbol{A}(\boldsymbol{t}_s), \; i=1, ..., N_e`, and forcing terms if any (including the hard forcing amplitudes :math:`\boldsymbol{B}`).  

At this stage, there is a key difference to highlight between the present approach and the application of the MMS to the system in oscillator form. In the complex form, the secular terms cancelled are those oscillating at :math:`e^{+\textrm{j} \omega_{i,0} t_0}`, but not :math:`e^{-\textrm{j} \omega_{i,0} t_0}`. Therefore, when reconstructing the oscillator's motion :math:`x_{i} = z_{i} + \bar{z}_{i}`, there can be higher order terms contributing to the term in :math:`\cos(\omega_{i,0} t_0), \; \sin(\omega_{i,0} t_0)`. In other words, the higher order particular solutions enrich the main oscillations term. On the contrary, with the oscillator form there are no contributions for higher order particular solutions. Indeed, the cancellation of secular terms is more stringent, imposing no terms in :math:`\cos(\omega_{i,0} t_0)` nor :math:`\sin(\omega_{i,0} t_0)`, therefore setting to zero both :math:`e^{+\textrm{j} \omega_{i,0} t_0}` and :math:`e^{-\textrm{j} \omega_{i,0} t_0}` terms. Contributions to the main oscillations term can therefore only arise from the homogeneous higher order solutions, but there are no conditions imposing a proper choice of higher order homogeneous amplitudes. These conditions can only be derived from a formulation of the system's dynamics that embeds its state space, as shown in :cite:`nayfehResolvingControversiesApplication2005`. 

There is a choice left for the higher order homogeneous solutions. With the present complex approach, it is relevant to set these to zero so as to fully capture oscillations at :math:`e^{+\textrm{j} \omega_{i,0} t_0}` through the leading order term. Their amplitude in the full solutions :math:`z_i` are thus simply the complex amplitudes :math:`A_i`.

As with the oscillator form, the process described above is repeated successively at each order, i.e. the computed solutions are introduced in the next higher order system of equations, 
secular terms are cancelled and the next higher order solutions are computed. 

The secular terms are identified in :meth:`~oscilate.MMS.mms_complex.Multiple_scales_complex.secular_analysis` and the leading order solutions are computed in :meth:`~oscilate.MMS.mms_complex.Multiple_scales_complex.sol_higher_order`. 

Secular analysis and modulation equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The secular analysis is carried out in the same way as with the oscillator form. A small, practical difference is that the secular terms are multiplied by :math:`\textrm{j}` to retrieve the system of complex modulation equations

.. math::
    \begin{cases}
    \textrm{D}_1 A_i(\boldsymbol{t}_s) & = f_{A_i}^{(1)}(\boldsymbol{A}, t_1), \\
    & \vdots \\
    \textrm{D}_{N_e} A_i(\boldsymbol{t}_s) & = f_{A_i}^{(N_e)}(\boldsymbol{A}, t_1),
    \end{cases}

with a structure similar to that of the oscillator form. From there on, the derivation of modulation equations is unchanged. Still, one must keep in mind that these modulation equations will govern the amplitude and phase evolution of the complex coordinates, only indirectly that of the oscillators' motion. 

The complex modulation equations are derived in :meth:`~oscilate.MMS.mms_complex.Multiple_scales_complex.secular_analysis`.


Reconstruction of the oscillators' motion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The oscillators' motion is simply reconstructed at each order from the definition of the complex coordinates. For oscillator :math:`i` and at order :math:`j`, this results in 

.. math::
    x_{i,j}(\boldsymbol{t}) = z_{i,j}(\boldsymbol{t}) + \bar{z}_{i,j}(\boldsymbol{t}).

At this stage, the complex form approach branches back with the oscillator form one. 

The oscillators' motion is retrieved using :meth:`~oscilate.MMS.mms_complex.Multiple_scales_complex.sol_xO`.

