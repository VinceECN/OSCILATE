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

.. literalinclude:: ../../../examples/CPVA_SH.py
   :language: python
   :linenos: