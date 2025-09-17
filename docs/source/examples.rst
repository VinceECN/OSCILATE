Application Examples
====================

Here are some practical examples of using **MMS_solver**.

Example 1: Duffing oscillator
-----------------------------

MMS example on the Duffing oscillator subject to harmonic forcing. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), sections 4.1 and 4.1.1.

System description
^^^^^^^^^^^^^^^^^^

The system's equation is

.. math::
    \ddot{x} + c \dot{x} + \omega_0^2 x + \gamma x^3 = F \cos(\omega t),

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c` is the linear viscous damping coefficient,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\gamma` is the nonlinear coefficient,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A direct response is sought so the frequency (either backbone curve frequency or forcing frequency) is close from the oscilltor's such that

.. math::
    \omega = \omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the reference frequency :math:`\omega_0`.

The parameters are then scaled to indicate how weak they are:

- :math:`c = \epsilon^{N_e} \tilde{c}` indicates that damping is weak (of order :math:`N_e`),
- :math:`F = \epsilon^{N_e} \tilde{F}` indicates that forcing is weak (of order :math:`N_e`),
- :math:`\gamma = \epsilon \tilde{\gamma}` indicates that nonlinearities are weak (of order 1).

Note that :math:`N_e` is the order up to which the solutions are sought and the time scales are constructed. 
It is therefore chosen here to scale the damping and forcing at maximum order. 

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system up to order :math:`N_e = 3`. The MMS results can be visualized in LaTeX in an **IPython-based interactive environment** (e.g., VS Code's Python Interactive Window or Jupyter Notebook). For example, to visualize the evolution function :math:`f_\beta(a, \beta)`, use

  .. prompt:: python

     mms.sol.fbeta  # Display the symbolic expression for :math:`f_\beta(a, \beta)`

- Evaluate the MMS results at steady state. The solutions are stored in ``ss.sol`` (e.g., ``ss.sol.fa``, ``ss.sol.fbeta``).
- Compute the forced response and the backbone curve. These results are also stored in ``ss.sol`` (e.g., ``ss.sol.sigma``, ``ss.sol.F``, ``ss.sol.omega_bbc``).
- Evaluate the stability of the computed forced solution. The stability results are stored in ``ss.stab`` (e.g., ``ss.stab.eigvals``, ``ss.stab.bif_a``).
- Evaluate the steady-state results for given numerical values of the parameters and plot the results using methods like ``ss.plot_FRC()`` or ``ss.plot_ARC()``.

.. literalinclude:: ../../examples/Duffing_direct.py
   :language: python
   :linenos:



Example 2: Coupled Duffings in 1:3 internal resonance
-----------------------------------------------------

MMS example on coupled Duffing oscillators in 1:3 internal resonance subject to harmonic forcing. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), sections 6.6 and 6.6.2.

System description
^^^^^^^^^^^^^^^^^^

The system's equations are

.. math::
    
    \begin{cases}
        \ddot{x}_{0} + 2 \mu_{0} \dot{x}_{0} + \omega_{0}^{2} x_{0} + \alpha_{1} x_{0}^{3} + \alpha_{2} x_{0}^{2} x_{1} + \alpha_{3} x_{0} x_{1}^{2} + \alpha_{4} x_{1}^{3}    & = \Gamma_0 F \cos(\omega t), \\
        \ddot{x}_{1} + 2 \mu_{1} \dot{x}_{1} + \omega_{1}^{2} x_{1} + \alpha_{5} x_{0}^{3} + \alpha_{6} x_{0}^{2} x_{1} + \alpha_{7} x_{0} x_{1}^{2} + \alpha_{8} x_{1}^{3} & = \Gamma_1 F \cos(\omega t),
    \end{cases}
    

where 

- :math:`x_0,\; x_1` are the oscillators' coordinates,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`\mu_0,\; \mu_1` are the linear viscous damping coefficients,
- :math:`\omega_0,\; \omega_1` are the oscillator's natural frequencies,
- :math:`\alpha_i,\; i\in\{1,\cdots,8\}` are nonlinear coefficients,
- :math:`\Gamma_0,\; \Gamma_1` are forcing coefficients,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

The oscillators are in 1:3 internal resonance. Taking :math:`\omega_0` as the reference frequency, this implies

.. math::
    \omega_1 = 3\omega_0 + \sigma_1

where :math:`\sigma_1` is the detuning of oscillator 1 wrt the 1:3 internal resonance condition.

A response around :math:`3\omega_0 \approx \omega_1` is sought so the frequency is set to

.. math::
    \omega = 3\omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the frequency :math:`3\omega_0`.

The parameters are then scaled to indicate how weak they are:

- :math:`\mu_i = \epsilon \tilde{\mu}_i` indicates that damping is weak,
- :math:`F = \epsilon \tilde{F}` indicates that forcing is weak,
- :math:`\sigma_1 = \epsilon \tilde{\sigma}_1` indicates that detuning is small,
- :math:`\alpha_i = \epsilon \tilde{\alpha}_i` indicates that nonlinearities are weak.

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the backbone curve when only oscillator 1 responds. 

.. literalinclude:: ../../examples/Duffing_13.py
   :language: python
   :linenos:




Example 3: Coupled Duffings in 1:1 internal resonance
-----------------------------------------------------

MMS example on coupled Duffing oscillators in 1:1 internal resonance subject to harmonic forcing. 

System description
^^^^^^^^^^^^^^^^^^

The system's equations are

.. math::
    
    \begin{cases}
        \ddot{x}_{0} + \omega_{0}^{2} x_{0} + c_{0} \dot{x}_{0} + \gamma_{0} x_{0}^{3} + \gamma_{01} x_{0} x_{1}^{2} & = F \cos(\omega t), \\
        \ddot{x}_{1} + \omega_{0}^{2} x_{1} + c_{1} \dot{x}_{1} + \gamma_{1} x_{1}^{3} + \gamma_{01} x_{0}^{2} x_{1} & = F \cos(\omega t),
    \end{cases}
    

where 

- :math:`x_0,\; x_1` are the oscillators' coordinates,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c_0,\; c_1` are the linear viscous damping coefficients,
- :math:`\omega_0` is the oscillators' natural frequency,
- :math:`\gamma_0,\; \gamma_1,\; \gamma_{01}` are nonlinear coefficients,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A response around :math:`\omega_0` is sought so the frequency is set to

.. math::
    \omega = \omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the oscillators' frequency :math:`\omega_0`.

The parameters are then scaled to indicate how weak they are:

- :math:`c_i = \epsilon \tilde{c}_i` indicates that damping is weak,
- :math:`F = \epsilon \tilde{F}` indicates that forcing is weak,
- :math:`\gamma_i = \epsilon \tilde{\gamma}_i, \; i\in\{0, 1, 01\}` indicates that nonlinearities are weak.

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the backbone curve and forced response when only oscillator 1 responds,
- Evaluate the stability of the forced solution,
- Compute the coupled-mode backbone curve. 

.. literalinclude:: ../../examples/Duffing_11_direct.py
   :language: python
   :linenos:



Example 4: Coupled quadratic oscillators in 1:2 internal resonance
------------------------------------------------------------------

MMS example on coupled quadratic oscillators in 1:2 internal resonance subject to harmonic forcing. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), sections 6.5 and 6.5.1.


System description
^^^^^^^^^^^^^^^^^^

The system's equations are

.. math::
    
    \begin{cases}
        \ddot{x}_{0} + 2 \mu_{0} \dot{x}_{0} + \omega_{0}^{2} x_{0} + \alpha_{0} x_{0} x_{1} & = \eta_0 F \cos(\omega t), \\
        \ddot{x}_{1} + 2 \mu_{1} \dot{x}_{1} + \omega_{1}^{2} x_{1} + \alpha_{1} x_{0}^{2}   & = \eta_1 F \cos(\omega t),
    \end{cases}
    

where 

- :math:`x_0,\; x_1` are the oscillators' coordinates,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`\mu_0,\; \mu_1` are the linear viscous damping coefficients,
- :math:`\omega_1,\; \omega_1` are the oscillator's natural frequencies,
- :math:`\alpha_0,\; \alpha_1` are nonlinear coefficients,
- :math:`\eta_0,\; \eta_1` are forcing coefficients,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

The oscillators are in 1:2 internal resonance. Taking :math:`\omega_0` as the reference frequency, this implies

.. math::
    \omega_1 = 2\omega_0 + \sigma_1

where :math:`\sigma_1` is the detuning of oscillator 1 wrt the 1:2 internal resonance condition.

A response around :math:`2\omega_0 \approx \omega_1` is sought so the frequency is set to

.. math::
    \omega = 2\omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning.

The parameters are then scaled to indicate how weak they are:

- :math:`\mu_i = \epsilon \tilde{\mu}_i` indicates that damping is weak,
- :math:`\eta_0 = \epsilon \tilde{\eta}_0` indicates that forcing on oscillator 0 is weak,
- :math:`\eta_1 = \epsilon \tilde{\eta}_1` indicates that forcing on oscillator 1 is one order weaker than that on oscillator 0,
- :math:`\sigma_1 = \epsilon \tilde{\sigma}_1` indicates that detuning is small,
- :math:`\alpha_i = \epsilon \tilde{\alpha}_i` indicates that nonlinearities are weak.

In addition, the solutions are sought with leading order terms at order :math:`\epsilon` rather than :math:`\epsilon^0=1` which was used in previous examples.
This is controled through ``eps_pow_0=1``.

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the backbone curve when only oscillator 0 responds,
- Compute the coupled-mode forced response. 

.. literalinclude:: ../../examples/quadratic_12.py
   :language: python
   :linenos:



Example 5: Parametrically excited Duffing oscillator
----------------------------------------------------

MMS example on the Duffing oscillator subject to parametric forcing. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), section 5.7.3.

System description
^^^^^^^^^^^^^^^^^^

The system's equation is

.. math::
    \ddot{x} + c \dot{x} + \omega_0^2 x + \gamma x^3 = -2 x F \cos(\omega t),

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c` is the linear viscous damping coefficient,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\gamma` is the nonlinear coefficient,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A parametric response at twice the oscillator's frequency is sought so the frequency is set to

.. math::
    \omega = 2\omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the reference frequency :math:`\omega_0`.

The parameters are then scaled to indicate how weak they are:

- :math:`c = \epsilon \tilde{c}` indicates that damping is weak,
- :math:`F = \epsilon \tilde{F}` indicates that forcing is weak,
- :math:`\gamma = \epsilon \tilde{\gamma}` indicates that nonlinearities are weak.


Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the forced response and the backbone curve,
- Evaluate the stability of the computed forced solution.

.. literalinclude:: ../../examples/Duffing_parametric.py
   :language: python
   :linenos:



Example 6: Van der Pol oscillator
---------------------------------

MMS example on a Van der Pol oscillator. 

System description
^^^^^^^^^^^^^^^^^^

The system's equation is

.. math::
    \ddot{x} + \omega_{0}^{2} x + \mu \left(x^{2} - 1\right) \dot{x} = 0,

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\mu` is the linear and nonlinear damping coefficient.

A parametric response around the oscillator's frequency is sought so the frequency is set to

.. math::
    \omega = \omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the reference frequency :math:`\omega_0`.

The parameter :math:`\mu = \epsilon \tilde{\mu}` is then scaled such that 

.. math::
    \mu = \epsilon \tilde{\mu}

to indicate that linear and nonlinear dampings are weak. 

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state.

.. literalinclude:: ../../examples/Van_der_Pol.py
   :language: python
   :linenos:



Example 7: Parametrically excited Van der Pol oscillator
--------------------------------------------------------

MMS example on a Van der Pol oscillator subject to parametric forcing. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), section 5.7.2, where the Van der Pol oscillator is also called Rayleigh oscillator.

System description
^^^^^^^^^^^^^^^^^^

The system's equation is

.. math::
    \ddot{x} + c \dot{x} + \gamma \dot{x}^{3} + \omega_{0}^{2} x = -2x F \cos(\omega t),

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c` is the linear viscous damping coefficient,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\gamma` is the nonlinear damping coefficient,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A parametric response around twice the oscillator's frequency is sought so the frequency is set to

.. math::
    \omega = 2\omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning.

The parameters are then scaled to indicate how weak they are:

- :math:`c = \epsilon \tilde{c}` indicates that damping is weak,
- :math:`F = \epsilon \tilde{F}` indicates that forcing is weak,
- :math:`\gamma = \epsilon \tilde{\gamma}` indicates that nonlinear damping is weak.

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the forced response and the backbone curve,
- Evaluate the stability of the computed forced solution.

.. literalinclude:: ../../examples/Van_der_Pol_parametric.py
   :language: python
   :linenos:


Example 8: Superharmonic response of a Duffing oscillator
---------------------------------------------------------

MMS example on a Duffing oscillator subject to hard forcing triggering a superharmonic response. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), sections 4.1.2 and 4.1.3.

System description
^^^^^^^^^^^^^^^^^^

The system's equation is

.. math::
    \ddot{x} + c \dot{x} + \gamma \dot{x}^{3} + \omega_{0}^{2} x = -2x F \cos(\omega t),

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c` is the linear viscous damping coefficient,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\gamma` is the nonlinear coefficient,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A parametric response around :math:`1/3` times the oscillator's frequency is sought so the frequency is set to

.. math::
    \omega = \frac{1}{3}\omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning.

The parameters are then scaled to indicate how weak they are:

- :math:`c = \epsilon \tilde{c}` indicates that damping is weak,
- :math:`\gamma = \epsilon \tilde{\gamma}` indicates that nonlinear damping is weak.

Note that the forcing is not scaled, i.e. :math:`F` appears at leading order. This is called hard forcing. 

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the forced response and the backbone curve,

.. literalinclude:: ../../examples/Duffing_superharmonic.py
   :language: python
   :linenos:

Example 9: Subharmonic response of a Duffing oscillator
-------------------------------------------------------

MMS example on a Duffing oscillator subject to hard forcing triggering a subharmonic response. 
This configuration was studied by Nayfeh and Mook in *Nonlinear Oscillations* (1995), sections 4.1.2 and 4.1.4.

System description
^^^^^^^^^^^^^^^^^^

The system's equation is

.. math::
    \ddot{x} + c \dot{x} + \gamma \dot{x}^{3} + \omega_{0}^{2} x = -2x F \cos(\omega t),

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c` is the linear viscous damping coefficient,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\gamma` is the nonlinear coefficient,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A parametric response around 3 times the oscillator's frequency is sought so the frequency is set to

.. math::
    \omega = 3\omega_0 + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning.

The parameters are then scaled to indicate how weak they are:

- :math:`c = \epsilon \tilde{c}` indicates that damping is weak,
- :math:`\gamma = \epsilon \tilde{\gamma}` indicates that nonlinear damping is weak.

Note that the forcing is not scaled, i.e. :math:`F` appears at leading order. This is called hard forcing. 

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system.
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Compute the forced response and the backbone curve,

.. literalinclude:: ../../examples/Duffing_subharmonic.py
   :language: python
   :linenos:



Example 10: Subharmonic response of parametrically excited coupled centrifugal pendulums in 1:1 internal resonance
------------------------------------------------------------------------------------------------------------------

MMS example on two coupled centrifugal pendulums in 1:1 internal resonance subject to direct and parametric forcings, triggering a subharmonic response. 
The pendulums are indirectly and weakly coupled through the rotor that supports them.
This configuration was studied in V. Mahé et al., *Subharmonic centrifugal pendulum vibration absorbers allowing a rotational mobility*, MSSP (2022).

System description
^^^^^^^^^^^^^^^^^^

The system's equations in modal space are

.. math::
    
    \begin{cases}
        \ddot{\zeta}_1 + n_p^2 \zeta_1 & = f_1(\boldsymbol{\zeta}, \dot{\boldsymbol{\zeta}}, \ddot{\boldsymbol{\zeta}}, t) \\
        \ddot{\zeta}_2 + n_p^2 \zeta_2 & = f_2(\boldsymbol{\zeta}, \dot{\boldsymbol{\zeta}}, \ddot{\boldsymbol{\zeta}}, t),
    \end{cases}
    

where 

- :math:`\zeta_1,\; \zeta_2` are the (modal) oscillators' coordinates,
- :math:`\boldsymbol{\zeta}` is the vector of oscillator coordinates,
- :math:`t` is the angular position of the rotor supporting the pendulums, analogous to time
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a angle derivative,
- :math:`n_p` is the pendulums' tuning order and their natural frequency when they are uncoupled from the rotor,

- :math:`f_1(\boldsymbol{\zeta}, \dot{\boldsymbol{\zeta}}, \ddot{\boldsymbol{\zeta}}, t)` and :math:`f_2(\boldsymbol{\zeta}, \dot{\boldsymbol{\zeta}}, \ddot{\boldsymbol{\zeta}}, t)` are (small) functions involving linear damping, an additional linear stiffness on oscillator 2, nonlinear terms, and both direct and parametric forcings. They take the form:

  .. math::

      \begin{aligned}
      \begin{split}
      f_1(\boldsymbol{\zeta}, \dot{\boldsymbol{\zeta}}, \ddot{\boldsymbol{\zeta}}, t) = & -\Lambda_m^{-1} \Big[ (\Lambda_m \dot{\zeta}_1 - n_t^2(1+n_t^2)\zeta_1\zeta_2)(\mu n_p^2 \Lambda_c \zeta_2 + T_1\cos(\omega t)) \\
      & + 2 \mu n_t^2 \Lambda_m \dot{\zeta}_1 (\zeta_1\dot{\zeta}_1 + \zeta_2 \dot{\zeta}_2) \\
      & + 6 \eta \alpha_1 \alpha_3 (\zeta_1 \dot{\zeta}_1^2 + 2 \zeta_2 \dot{\zeta}_1 \dot{\zeta}_2 + \zeta_1 \dot{\zeta}_2^2 + \zeta_1^2 \ddot{\zeta}_1 + 2 \zeta_1 \zeta_2 \ddot{\zeta}_2 + \zeta_2^2 \ddot{\zeta}_1) \\
      & - 2 x_4 (\zeta_1^3 + 3\zeta_1\zeta_2^2) + b \dot{\zeta}_1 \Big],
      \end{split} \\[
      0.5em]
      \begin{split}
      f_2(\boldsymbol{\zeta}, \dot{\boldsymbol{\zeta}}, \ddot{\boldsymbol{\zeta}}, t) = & -\Lambda_m^{-1} \Bigg[ \left(\Lambda_c + \Lambda_m \dot{\zeta}_2 - \frac{n_t^2(1+n_t^2)}{2}(\zeta_1^2 + \zeta_2^2)\right)(\mu n_p^2 \Lambda_c \zeta_2 + T_1\cos(\omega t)) \\
      & + 2 \mu n_t^2 (\Lambda_c + \Lambda_m \dot{\zeta}_2) (\zeta_1 \dot{\zeta}_1 + \zeta_2 \dot{\zeta}_2) \\
      & - \Lambda_c \mu n_p^2\frac{n_t^2(1+n_t^2)}{2}(3\zeta_1^2\zeta_2 + \zeta_2^3) \\
      & + \Lambda_c \mu n_t^2(1+n_t^2)(2\zeta_1 \dot{\zeta}_1 \dot{\zeta}_2 + \zeta_2\dot{\zeta}_1^2 + \zeta_2 \dot{\zeta}_2^2) \\
      & + 6 \eta \alpha_1 \alpha_3 (2\zeta_1 \dot{\zeta}_1 \dot{\zeta}_2 + \zeta_2 \dot{\zeta}_1^2 + \zeta_2 \dot{\zeta}_2^2 + 2\zeta_1 \ddot{\zeta}_1 \zeta_2 + \ddot{\zeta}_2 \zeta_1^2 + \ddot{\zeta}_2 \zeta_2^2) \\
      & - 2\tilde{x}_{4}(3\zeta_1^2\zeta_2 + \zeta_2^3) + \tilde{b}\dot{\zeta}_2 \Bigg].
      \end{split}
      \end{aligned}

  where

  - :math:`\Lambda_m` is an equivalent pendulum mass,
  - :math:`\Lambda_c` is a linear coupling coefficient among the pendulums,
  - :math:`\mu` is a linear coupling coefficient between the pendulums and the rotor that supports them,
  - :math:`n_t` and :math:`x_4` are linear and nonlinear pendulums' path coefficients, respectively,
  - :math:`\alpha_1` and :math:`\alpha_3` are linear and nonlinear pendulums' rotation coefficients, respectively,
  - :math:`b` is the linear viscous damping coefficient,
  - :math:`T_1` is the forcing amplitude,
  - :math:`\omega` is the forcing frequency (often written :math:`n` as it is actually a forcing order).

A response around :math:`n_p` is sought so the frequency is set to

.. math::
    \omega = n_p + \epsilon \sigma

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the pendulums' tuning order.

The parameters are then scaled to indicate how weak they are:

- :math:`\mu = \epsilon \tilde{\mu}` indicates that pendulums are weakly coupled,
- :math:`b = \epsilon \tilde{b}` indicates that damping is weak,
- :math:`T_1 = \epsilon \tilde{T}_1` indicates that forcing is weak,
- :math:`x_4 = \epsilon \tilde{x}_4, \; \alpha_3 = \tilde{\alpha}_3` indicates that nonlinearities are weak.

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical (modal) system,
- Apply the MMS to the system,
- Evaluate the MMS results at steady state,
- Introduce compact notations,
- Compute the backbone curve and forced response when only oscillator (mode) 1 responds,
- Evaluate the stability of the forced solution using cartesian coordinates,

.. literalinclude:: ../../examples/CPVA_SH.py
   :language: python
   :linenos: