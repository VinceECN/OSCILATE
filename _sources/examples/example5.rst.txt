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

.. literalinclude:: ../../../examples/Duffing_parametric.py
   :language: python
   :linenos:

Validation
^^^^^^^^^^
The notations used in the above code are related to those in *Nonlinear Oscillations* as follows:

.. table:: Link between notations from *Nonlinear Oscillations* and this document.

   ================================  =============================  ==================================
   Current Notation                  *Nonlinear Oscillations*       Description
   ================================  =============================  ==================================
   :math:`u`                         :math:`x`                      Oscillator's dof
   :math:`\omega_0`                  :math:`\omega`                 Natural frequency
   :math:`\omega`                    1                              Forcing frequency   
   :math:`2\omega_0+\epsilon\sigma`  :math:`\omega+\epsilon\sigma`  Forcing frequency definition
   :math:`F`                         :math:`1`                      Forcing amplitude
   :math:`\tilde{c}`                 :math:`2\mu`                   Damping coefficient (scaled)
   :math:`\tilde{\gamma}`            :math:`\alpha`                 Nonlinear coefficient (scaled)
   :math:`x_0^{\textrm{h}}`          [no name]                      Leading order homogeneous solution
   :math:`a_0`                       :math:`a`                      Oscillator's amplitude
   :math:`2\beta_0`                  :math:`\psi`                   Oscillator's autonomous phase
   :math:`\textrm{D}_1(\bullet)`     :math:`(\bullet)'`             Slow time derivative
   :math:`\lambda_i`                 :math:`\lambda_i`              Eigenvalues of the Jacobian matrix
   ================================  =============================  ==================================

Note that the above table implies that the detuning :math:`\sigma` in *Nonlinear Oscillations*'s notations is half that of the current notation.
Section 5.7.3 from *Nonlinear Oscillations* gives the following results

.. math::
    \begin{cases}
    u(t)          & = a \cos(t - \frac{1}{2}\psi) + \mathcal{O}(\epsilon), \\
    a'            & = - \dfrac{a}{2 \omega} \sin \psi - \mu a, \\
    a\psi'        & = 2\sigma a - \dfrac{a}{\omega} \cos \psi - \dfrac{3 \alpha}{4 \omega} a^3, \\
    \lambda_{1,2} & = - \mu \pm \sqrt{\mu^2 + \frac{3}{4} \alpha a^2 \cos \psi}, \quad \cos \psi = 2 \sigma \omega - \frac{3}{4} \alpha a^2.
    \end{cases}

In the current notations, this is equivalent to

.. math::
    \begin{cases}
    x_0^{\textrm{h}}(t)      & = a_0 \cos(\frac{\omega}{2} t  - \beta_0), \\
    \textrm{D}_1 a_0         & = - \dfrac{\tilde{c} a_{0}}{2} - \dfrac{\tilde{F} a_{0} \sin{\left(2 \beta_{0} \right)}}{2 \omega_{0}}, \\
    a_0 \textrm{D}_1 \beta_0 & = \dfrac{\sigma a_{0}}{2} - \dfrac{\tilde{F} a_{0} \cos{\left(2 \beta_{0} \right)}}{2 \omega_{0}} - \dfrac{3 \tilde{\gamma} a_{0}^{3}}{8 \omega_{0}}, \\
    \lambda_{1,2}            & = \epsilon \dfrac{1}{2} \dfrac{ \left(- \omega_{0} \tilde{c} \pm \sqrt{\omega_{0}^{2} \tilde{c}^{2} + 3 \omega_{0} \sigma \tilde{\gamma} a_{0}^{2} - \frac{9}{4} \tilde{\gamma}^{2} a_{0}^{4}}\right)}{ \omega_{0}},
    \end{cases}

which are the outputs for ``ss.sol.x[0][0].simplify()``, ``ss.sol.faO[0][1]``, ``ss.sol.fbetaO[0][1]`` and ``ss.stab.eigvals``, respectively. 
Note the :math:`\epsilon` factor for the eigenvalues. It is related to the eigenvalues being computed from the physical time evolution equations rather than just the :math:`1^{\textrm{st}` order slow time evolution equations.
