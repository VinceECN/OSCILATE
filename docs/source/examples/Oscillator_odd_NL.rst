Oscillator with odd nonlinearities
----------------------------------

MMS example of an oscillator with odd (cubic, quintic and septic) nonlinearities subject to a direct harmonic forcing. 

System description
^^^^^^^^^^^^^^^^^^

.. figure:: /_static/examples/oscillator_odd_NL.svg
   :alt: Nonlinear system.
   :width: 70%
   :align: center

   Illustration of a forced nonlinear oscillator with odd nonlinearities of degree 3, 5 and 7.

The system's equation is

.. math::
    \ddot{x} + c \dot{x} + \omega_0^2 x + \gamma_3 x^3 + \gamma_5 x^5 + \gamma_7 x^7 = F \cos(\omega t),

where 

- :math:`x` is the oscillator's coordinate,
- :math:`t` is the time,
- :math:`\dot{(\bullet)} = \mathrm{d}(\bullet)/\mathrm{d}t` is a time derivative,
- :math:`c` is the linear viscous damping coefficient,
- :math:`\omega_0` is the oscillator's natural frequency,
- :math:`\gamma_3,\; \gamma_5, \; \gamma_7` are the nonlinear coefficients,
- :math:`F` is the forcing amplitude,
- :math:`\omega` is the forcing frequency.

A direct response is sought so the frequency (either backbone curve frequency or forcing frequency) is close from the oscillator's, such that

.. math::
    \omega = \omega_0 + \epsilon \sigma,

where

- :math:`\epsilon` is a small parameter involved in the MMS,
- :math:`\sigma` is the detuning wrt the reference frequency :math:`\omega_0`.

The parameters are then scaled to indicate how weak they are:

- :math:`c = \epsilon^{N_e} \tilde{c}` indicates that damping is weak (of order :math:`N_e`),
- :math:`F = \epsilon^{N_e} \tilde{F}` indicates that forcing is weak (of order :math:`N_e`),
- :math:`\gamma_i = \epsilon \tilde{\gamma}_i, \; i=3, 5, 7` indicates that nonlinearities are weak (of order 1).

Note that :math:`N_e` is the order up to which the solutions are sought and the time scales are constructed. 
It is therefore chosen here to scale the damping and forcing at maximum order. 

Code description
^^^^^^^^^^^^^^^^
The script below allows to

- Construct the dynamical system
- Apply the MMS to the system up to order :math:`N_e=3`
- Evaluate the MMS results at steady state
- Compute the forced response and the backbone curve
- Evaluate the stability of the computed forced solution. 
- Evaluate the steady-state results for given numerical values of the parameters and plot the results.

.. literalinclude:: ../../../examples/Oscillator_odd_NL.py
   :language: python
   :linenos:

