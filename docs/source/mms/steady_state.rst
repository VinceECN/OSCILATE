Evaluation at steady state
--------------------------

Steady state solutions
~~~~~~~~~~~~~~~~~~~~~~

Steady state conditions
^^^^^^^^^^^^^^^^^^^^^^^

At steady state, the solutions' amplitudes and phases are time-independent. One therefore has, for each oscillator :math:`i=1,...,N`,

.. math::
    \begin{cases}
    \dfrac{\textrm{d}}{dt} a_i & = 0, \\
    \dfrac{\textrm{d}}{dt} \beta_i & = 0, \\
    \end{cases}

and the homogeneous steady state solutions take the form

.. math::
    x^{\textrm{h}}_{i,0}(t) = a_i \cos\left( \frac{r_i}{r_{\textrm{MMS}}} \omega t - \beta_i \right).    

MMS modulation equations at steady state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The steady state amplitudes :math:`\boldsymbol{a}` and phases :math:`\boldsymbol{\beta}` are governed by the modulation equations evaluated at steady state, which take the form

.. math::
    \begin{cases}
    f_{a_0}(\boldsymbol{a}, \boldsymbol{\beta})     & = 0, \\
    f_{\beta_0}(\boldsymbol{a}, \boldsymbol{\beta}) & = 0, \\
    & \vdots \\
    f_{a_{N-1}}(\boldsymbol{a}, \boldsymbol{\beta})     & = 0, \\
    f_{\beta_{N-1}}(\boldsymbol{a}, \boldsymbol{\beta}) & = 0.
    \end{cases}

This is now an algebraic system of equations. 

MMS solutions at steady state
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Solving the modulation equations at steady state yields the steady state solutions :math:`\boldsymbol{a}, \boldsymbol{\beta}`. 
The system can be solved directly for :math:`\boldsymbol{a}` and :math:`\boldsymbol{\beta}`, yielding explicit analytical solutions, but this is often complex as the system is nonlinear. 

A possibility that is sometimes available to obtain **analytical solutions** is to rearrange the equations to isolate the phase terms :math:`\cos(f(\boldsymbol{\beta})),\; \sin(f(\boldsymbol{\beta}))` where :math:`f(\boldsymbol{\beta})` is a linear function of the :math:`\beta_i`. 
Then the equations can be squared and summed up to obtain an equation on :math:`a_i` only and/or :math:`a_j` as a function of :math:`a_i`. The :math:`a_j` can be expressed as a function of :math:`a_i`, leading to a polynomial equation on :math:`a_i`.
The resulting polynomial equation can rarely be solved directly as the polynomial involved are often of high order. 
However, the polynomial is often quadratic in the detuning :math:`\sigma` and forcing amplitude :math:`F`. 
It can therefore be solved for :math:`\sigma` and :math:`F` with :math:`a_i` seen as a parameter. 
This yields an implicit solution for :math:`a_i(\sigma) \Rightarrow a_i(\omega)` and :math:`a_i(F)`, from which one can deduce the other amplitudes :math:`a_j` and phases :math:`\boldsymbol{\beta}`, thus reconstructing the oscillators' solutions.

The processus described above is not always feasible. The blocking points are typically to 

(i) get rid of phase terms :math:`\cos(f(\boldsymbol{\beta})),\; \sin(f(\boldsymbol{\beta}))` in the equations, 

(ii) express every amplitude :math:`a_j` in terms of a single one, :math:`a_i`,

(iii) end up with a polynomial of order 2 in :math:`\sigma` and :math:`F` 

These difficulties become more pronounced when the system involves several oscillators, in which case the amplitudes and phases may only be **computed numerically**.

To facilitate the derivation of an analytical solution, it is possible to consider the **backbone curve** (bbc) of the forced solution rather than the forced solution itself. 
This bbc is computed in the absence of damping and forcing, therefore simplifying the system. Typically, this reduces the number of phase terms appearing. 
The solving procedure is then the same as that described previously.

Stability analysis
~~~~~~~~~~~~~~~~~~
The stability analysis is described in details in :ref:`stability`.