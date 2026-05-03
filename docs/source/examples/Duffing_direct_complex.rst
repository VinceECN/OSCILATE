Duffing oscillator - complex form
---------------------------------

MMS example on the Duffing oscillator subject to a direct harmonic forcing. 
This configuration was studied by Nayfeh and Mook :cite:`nayfehNonlinearOscillations1995`, sections 4.1 and 4.1.1, but here it is solved using a complex form of the dynamical system, as in :cite:`nayfehResolvingControversiesApplication2005`.

System description
^^^^^^^^^^^^^^^^^^

.. figure:: /_static/examples/Duffing_direct.svg
   :alt: Nonlinear system.
   :width: 70%
   :align: center

   Illustration of a forced Duffing oscillator.

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

This system of equations is transformed through the introduction of the complex coordinates

.. math::
   \begin{cases}
   x       & = z + \bar{z}, \\
   \dot{x} & = \textrm{j} \omega_0 (z - \bar{z}),
   \end{cases}

or equivalently

.. math::
   \begin{cases}
   z       & = \frac{1}{2} (x - \frac{\textrm{j}}{\omega_0} \dot{x}), \\
   \bar{z} & = \frac{1}{2} (x + \frac{\textrm{j}}{\omega_0} \dot{x}),
   \end{cases}

resulting in the complex first order dynamical system (on :math:`z(t)`)

.. math::
   \dot{z} = \textrm{j} \omega z + \frac{\textrm{j} \gamma \left(z + \overline{z}\right)^{3}}{2 \omega} - \frac{c z}{2} + \frac{c \overline{z}}{2}.

A direct response is sought so the frequency (either backbone curve frequency or forcing frequency) is close from the oscillator's, such that

.. math::
    \omega = \omega_0 + \epsilon \sigma,

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

- Construct the dynamical system in oscillator and complex forms
- Apply the MMS to the system up to order :math:`N_e=5` using the oscillator and complex forms
- Evaluate the MMS results at steady state
- Compute the backbone curve and the associated mapping from the normal coordinate's amplitude to the displacement (initial) coordinate. These results are in accordance with those from the normal form theory :cite:`defigueiredostabileNormalFormAnalysis2025`. 

.. literalinclude:: ../../../examples/Duffing_direct_complex.py
   :language: python
   :linenos:

Plot outputs
^^^^^^^^^^^^
The plot outputs shown below are generated from the code above. 

The figure below displays the backbone curve as the peak oscillator amplitude as a function of the oscillation frequency. Three curves are shown, associated to the results obtained with the complex (blue) and oscillator (orange) forms, and the later if pushing the expansions only up to :math:`N_e=1` (green). One can see how these results deviate at large amplitudes.

.. figure:: /_static/examples/Duffing_direct_plots/BBC_complex_form.svg
   :alt: Backbone curve
   :width: 80%
   :align: center

   Backbone curve of the Duffing oscillator.
