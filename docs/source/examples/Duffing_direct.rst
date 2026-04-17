Duffing oscillator
------------------

MMS example on the Duffing oscillator subject to a direct harmonic forcing. 
This configuration was studied by Nayfeh and Mook :cite:`nayfehNonlinearOscillations1995`, sections 4.1 and 4.1.1.

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

- Construct the dynamical system.
- Apply the MMS to the system up to order :math:`N_e = 3`. The MMS results can be visualized in LaTeX in an **IPython-based interactive environment** (e.g., VS Code's Python Interactive Window or Jupyter Notebook). For example, to visualize the modulation function :math:`f_\beta(a, \beta)`, run

  .. prompt:: python

     mms.sol.fbeta  # Display the symbolic expression for :math:`f_\beta(a, \beta)`

  in the Interactive Window.

- Evaluate the MMS results at steady state. The solutions are stored in ``ss.sol`` (e.g., ``ss.sol.fa``, ``ss.sol.fbeta``).
- Compute the forced response and the backbone curve. These results are stored in ``ss.sol_forced`` and ``ss.sol_bbc`` (e.g., ``ss.sol_forced.sigma``, ``ss.sol_forced.F``, ``ss.sol_bbc.omega``).
- Evaluate the stability of the computed forced solution. The stability results are stored in ``ss.sol_forced.stab`` (e.g., ``ss.sol_forced.stab.eigvals``, ``ss.sol_forced.stab.bif_a``).
- Evaluate the steady-state results for given numerical values of the parameters and plot the results.

.. literalinclude:: ../../../examples/Duffing_direct.py
   :language: python
   :linenos:

Plot outputs
^^^^^^^^^^^^
The plot outputs shown below are generated from the code above. 

Frequency response curve
~~~~~~~~~~~~~~~~~~~~~~~~

The two figures below display the amplitude and phase responses (blue) of the :math:`1^{\text{st}}` harmonic of the Duffing oscillator as a function of the excitation frequency and the associated bifurcation curves (red), delimiting unstable from stable zones. The backbone curve is shown in grey and the linear frequency in black.

.. figure:: /_static/examples/Duffing_direct_plots/FRC_a.svg
   :alt: Frequency response curve - amplitude
   :width: 80%
   :align: center

   Frequency response curve of the Duffing oscillator (amplitude).

.. figure:: /_static/examples/Duffing_direct_plots/FRC_beta.svg
   :alt: Frequency response curve - phase
   :width: 80%
   :align: center

   Frequency response curve of the Duffing oscillator (phase).

Amplitude response curve
~~~~~~~~~~~~~~~~~~~~~~~~

The two figures below display the amplitude and phase responses (blue) of the :math:`1^{\text{st}}` harmonic of the Duffing oscillator as a function of the excitation amplitude.

.. figure:: /_static/examples/Duffing_direct_plots/ARC_a.svg
   :alt: Amplitude response curve - amplitude
   :width: 80%
   :align: center

   Amplitude response curve of the Duffing oscillator (amplitude).

.. figure:: /_static/examples/Duffing_direct_plots/ARC_beta.svg
   :alt: Amplitude response curve - phase
   :width: 80%
   :align: center

   Amplitude-response curve of the Duffing oscillator (phase).

