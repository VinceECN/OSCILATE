.. _example_VdP:

Van der Pol oscillator
----------------------

MMS example on a Van der Pol oscillator. 
This configuration is very close from that studied by Nayfeh and Mook :cite:`nayfehNonlinearOscillations1995`, section 3.3.4.


System description
^^^^^^^^^^^^^^^^^^

.. figure:: /_static/examples/VdP.svg
   :alt: Nonlinear system.
   :width: 50%
   :align: center

   Illustration of a Van der Pol oscillator. 

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
- Solve the modulation equations (yields the transient response),
- Evaluate the MMS results at steady state,
- Solve the steady state modulation equations (yields the limit cycle),
- Evaluate the symbolic expressions for given numerical parameters and plot the limit cycle and two trajectories in the phase portrait.

.. literalinclude:: ../../../examples/Van_der_Pol.py
   :language: python
   :linenos:

Plot outputs
^^^^^^^^^^^^
The plot output shown below is generated from the code above. 

The figure below displays the phase portrait of the Van der Pol oscillator. Three different trajectories are shown: the limit cycle (black), an external one (red) and an internal one (blue). The limit cycle is stable, such that the trajectories are converging towards it.

.. figure:: /_static/examples/VdP_plots/phase_portrait.svg
   :alt: Phase portrait
   :width: 80%
   :align: center

   Phase portrait of the Van der Pol oscillator.
