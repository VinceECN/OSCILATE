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

.. literalinclude:: ../../../examples/Duffing_13.py
   :language: python
   :linenos:
