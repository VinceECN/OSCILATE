Dynamical systems of interest
-----------------------------

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