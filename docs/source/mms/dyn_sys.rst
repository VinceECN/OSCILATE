Dynamical systems of interest
-----------------------------

Coupled oscillators system
~~~~~~~~~~~~~~~~~~~~~~~~~~

Systems considered are typically composed of :math:`N` coupled nonlinear equations of the form

.. math::
    \begin{cases}
    \ddot{x}_0 + \omega_0^2 x_0 & = f_0(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t), \\
    & \vdots \\
    \ddot{x}_{N-1} + \omega_{N-1}^2 x_{N-1} & = f_{N-1}(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t).
    \end{cases}

The :math:`x_i(t)` (:math:`i=0,...,N-1`) are the oscillators' coordinates, and represent the degrees of freedom (dof) of the system.

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

Internal resonance relations among oscillators can be specified in a second step by expressing the :math:`\omega_i` as a function of a reference frequency. See :class:`~oscilate.MMS.dyn_sys.Dynamical_system`. Detuning can also be introduced during this step.

.. figure:: /_static/coupled_NL_oscillators.svg
   :alt: Nonlinear system.
   :width: 80%
   :align: center

   Illustration of a system of coupled nonlinear oscillators that can be studied using :mod:`oscilate`. :math:`f_{\text{nl},i}(x_i, \dot{x}_i, \ddot{x}_i)`, :math:`g_{\text{nl},i}(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}})` and :math:`f_{\text{f},_i}(x_i, \dot{x}_i, \ddot{x}_i)` represent the intrinsic nonlinear of oscillator :math:`i`, the nonlinear couplings with other oscillators, and the parametric forcing coeffiicent. A unitary mass is considered for each oscillator, and :math:`\omega_i^2` and :math:`c_i` are the stiffness (natural frequencies) and the (weak) viscous damping coefficients. Note that this illustration does not represent weak linear terms nor complex coupled parametric forcing, which can be accounted for in the :mod:`oscilate` package.

.. _dyn_sys_complex:

Complex form
~~~~~~~~~~~~

It will later be shown that it may be convenient to rewrite the coupled :math:`2^{\text{nd}}` order equations as coupled :math:`1^{\text{st}}` order ones through the introduction of complex variables :math:`z_i(t)` (see section :ref:`mms_complex`). These are defined through the transformation

.. math::
    \begin{cases}
    x_i       & = z_i + \bar{z}_i, \\
    \dot{x}_i & = \textrm{j} \omega_{i} (z_i - \bar{z}_i),
    \end{cases}

where :math:`\bar{(\bullet)}` denotes the complex conjugate, or equivalently,  

.. math::
   \begin{cases}
   z_{i}       & = \frac{1}{2} (x_{i} - \frac{\textrm{j}}{\omega_i} \dot{x}_{i}), \\
   \bar{z}_{i} & = \frac{1}{2} (x_{i} + \frac{\textrm{j}}{\omega_i} \dot{x}_{i}).
   \end{cases}

The :math:`z_i(t)`, :math:`i=0,...,N-1`, complex variables therefore encapsulate the state space of the system. Now, noting that 

.. math::
    \ddot{x}_{i} = 2 \textrm{j} \omega_{i} \dot{z}_i + \omega_{i}^2 (z_{i} - \bar{z}_{i}),

one obtains the coupled, complex first order equations on the complex variables,

.. math::
    \begin{cases}
    \dot{z}_0 & = \textrm{j} \omega_0 z_0 + g_0(\boldsymbol{z}, \bar{\boldsymbol{z}}, \dot{\boldsymbol{z}}, t), \\
    & \vdots \\
    \dot{z}_{N-1} & = \textrm{j} \omega_{N-1} z_{N-1} + g_{N-1}(\boldsymbol{z}, \bar{\boldsymbol{z}}, \dot{\boldsymbol{z}}, t),
    \end{cases}

where

.. math::

    \boldsymbol{z}(t)^\intercal = [z_0(t), z_1(t), \cdots, z_{N-1}(t)]

is the vector of complex coordinates, and functions :math:`g_i` are related to :math:`f_i` such that

.. math::
    g_i(\boldsymbol{z}, \bar{\boldsymbol{z}}, \dot{\boldsymbol{z}}, t) = - \frac{\textrm{j}}{2\omega_i} f_i(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t).

The dynamical system is transformed from oscillator to complex form using :meth:`~oscilate.MMS.dyn_sys.Dynamical_system.complex_form`.
