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

.. literalinclude:: ../../../examples/quadratic_12.py
   :language: python
   :linenos:
