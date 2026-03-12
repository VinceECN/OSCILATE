Stability analysis
------------------

Stability information
~~~~~~~~~~~~~~~~~~~~~

Consider a steady state solution :math:`(\hat{\boldsymbol{a}} , \hat{\boldsymbol{\beta}})` such that, for :math:`i=1,...,N`,

.. math::
    \begin{cases}
    f_{a_i}(\hat{\boldsymbol{a}}, \hat{\boldsymbol{\beta}})     & = 0, \\
    f_{\beta_i}(\hat{\boldsymbol{a}}, \hat{\boldsymbol{\beta}}) & = 0.
    \end{cases}

The aim is to determine the stability state of that steady solution, which corresponds to a fixed point in the phase space. 


Jacobian matrix
^^^^^^^^^^^^^^^

Let's first introduce the vector of polar coordinates and polar modulation functions

.. math::
    \boldsymbol{x}^{(\textrm{p})\intercal} & = [a_0, \beta_0, \cdots, a_{N-1}, \beta_{N-1}], \\
    \boldsymbol{f}^{(\textrm{p})\intercal} & = [f_{a_0}, f_{\beta_0}^*, \cdots, f_{a_{N-1}}, f_{\beta_{N-1}}^*].

Note that the appearance of :math:`f_{\beta_i}^*` in :math:`\boldsymbol{f}^{(\textrm{p})}` requires :math:`a_i \neq 0`, which strongly constraints the type of steady state solutions that can be considered in the approach described below. 
To relax this constraint, one can use a change of coordinates from polar to cartesian ones. 
This will be discussed in following sections, after the description of this polar approach. 

Using the vectors of polar coordinates and modulation functions, one can write the modulation equations system as

.. math::
    \dfrac{\textrm{d} \boldsymbol{x}^{(\textrm{p})}}{\textrm{d}t} = \textrm{J}^{(\textrm{p})} \boldsymbol{x}^{(\textrm{p})},

where we introduced the Jacobian matrix 

.. math::
    \textrm{J}^{(\textrm{p})} 
    = \dfrac{\partial \boldsymbol{f}^{(\textrm{p})} }{ \partial \boldsymbol{x}^{(\textrm{p})} }
    = \begin{bmatrix}
    \frac{\partial f_{a_0}}{\partial a_0} & \frac{\partial f_{a_0}}{\partial \beta_0} & \cdots & \frac{\partial f_{a_0}}{\partial \beta_{N-1}} \\
    \frac{\partial f_{\beta_0}^*}{\partial a_0} & \frac{\partial f_{\beta_0}^*}{\partial \beta_0} & \cdots & \frac{\partial f_{\beta_0}^*}{\partial \beta_{N-1}} \\
    \vdots & \vdots & \ddots & \vdots \\
    \frac{\partial f_{\beta_{N-1}}^*}{\partial a_0} & \frac{\partial f_{\beta_{N-1}}^*}{\partial \beta_0} & \cdots & \frac{\partial f_{\beta_{N-1}}^*}{\partial \beta_{N-1}}
    \end{bmatrix}.


Perturbation of the steady state solution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let us now consider a small perturbation :math:`\tilde{\boldsymbol{x}}^{(\textrm{p})}` of the steady state solution such that

.. math::
    \boldsymbol{x}^{(\textrm{p})} = \hat{\boldsymbol{x}}^{(\textrm{p})} + \tilde{\boldsymbol{x}}^{(\textrm{p})} \quad \Leftrightarrow \quad \tilde{\boldsymbol{x}}^{(\textrm{p})} = \boldsymbol{x}^{(\textrm{p})} - \hat{\boldsymbol{x}}^{(\textrm{p})}.

Using a first order Taylor expansion for :math:`\textrm{d} \tilde{\boldsymbol{x}}^{(\textrm{p})} / \textrm{d} t`, one can write

.. math::
    \dfrac{\textrm{d} \tilde{\boldsymbol{x}}^{(\textrm{p})}}{\textrm{d}t} = \left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}} \tilde{\boldsymbol{x}}^{(\textrm{p})} + \mathcal{O}(||\tilde{\boldsymbol{x}}^{(\textrm{p})}||^2),

where :math:`\left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}}` denotes the Jacobian matrix evaluated on the steady state solution. 
The perturbation solution takes the form

.. math::
    \tilde{\boldsymbol{x}}^{(\textrm{p})} = \sum_{i=1}^{2N} C_i \boldsymbol{\psi}_i e^{\lambda_i t},

where :math:`(\lambda_i, \boldsymbol{\psi}_i),\; i=1, ..., 2N` are the eigensolutions of the Jacobian (evaluated on :math:`\hat{\boldsymbol{x}}^{(\textrm{p})}`).


Stability condition
^^^^^^^^^^^^^^^^^^^

The steady state solution :math:`\hat{\boldsymbol{x}}^{(\textrm{p})}` is considered stable if a small perturbation :math:`\tilde{\boldsymbol{x}}^{(\textrm{p})}` vanishes in time, 
such that solutions close from :math:`\hat{\boldsymbol{x}}^{(\textrm{p})}` are converging towards it. This condition is fulfilled if

.. math::
    \Re[\lambda_i] < 0, \quad \forall i,

meaning that all eigenvalues of the Jacobian evaluated on the steady state solution must have negative real parts.
If this condition is not met, the system is either quasi stable or unstable.


Bifurcations
^^^^^^^^^^^^

A bifurcation occurs when the stability state of the system changes. 

#. Simple bifurcations occur when at least one eigenvalue crosses the imaginary axis through 0.
    Such bifurcations include saddle node and pitchfork bifurcations, which cause jumps of the response and the appearance of lower symmetry solutions, respectively.

#. Neimark-Sacker bifurcations occur when a pair of complex conjugate eigenvalues with nonzero real parts cross the imaginary axis.
    These bifurcations lead to non periodic solutions. 

Simple bifurcations can be detected by evaluating the sign of the Jacobian, as

.. math::
    \det \left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}} = \prod_{i=1}^{2N} \lambda_i,

thereby making :math:`\det \left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}}` an important stability indicator.
Neimark-Sacker bifurcations are more difficult to detect. Information from the trace of :math:`\left.\textrm{J}^{(\textrm{p})}\right|_{\hat{\boldsymbol{x}}^{(\textrm{p})}}` can be considered, or the Routh-Hurwitz criterion can be used. This is not detailed here.


Bifurcation curves
^^^^^^^^^^^^^^^^^^

Bifurcation curves are curves constructed by evaluating the coordinates of bifurcation points when varying one or more parameters.
The stability state of a solution changes when the response curve crosses a bifurcation curve.

Stability analysis in cartesian coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Limitations of polar coordinates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned previously, the approach described above fails when one of the oscillator's leading order amplitude is 0.
Indeed, the Jacobian is constructed using the modulation functions

.. math::
        f_{\beta_i}^*(\boldsymbol{a}, \boldsymbol{\beta}) = \frac{f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta})}{a_i},

which are defined only if :math:`a_i=0`. 
This prevents evaluating the stability of

- Trivial solutions, for which no oscillator responds,

- 1 mode solutions, whose stability can be affected under perturbation from another mode.

These limitations can be overcome using a change of coordinates.

Cartesian coordinates
^^^^^^^^^^^^^^^^^^^^^

The leading order homogeneous solution for oscillator :math:`i` can be written as

.. math::
    \begin{split}
    x_{i,0}^{\textrm{h}}(t) & = a_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t - \beta_i \right), \\
                            & = a_i \cos(\beta_i) \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) 
                                + a_i \sin(\beta_i) \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right).
    \end{split}

It then appears natural to introduce the cartesian coordinates

.. math::
    \begin{cases}
    p_i = a_i \cos(\beta_i), \\
    q_i = a_i \sin(\beta_i),
    \end{cases}

in order to rewrite the solution as

.. math::
    x_{i,0}^{\textrm{h}}(t) = p_i \cos\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right) + q_i \sin\left(\frac{r_i}{r_{\textrm{MMS}}}\omega t\right).

In the following it will be convenient to use the cartesian coordinates vectors

.. math::
    \begin{aligned}
    \boldsymbol{p}(t)^\intercal & = [p_0(t), p_1(t), \cdots, p_{N-1}(t)], \\
    \boldsymbol{q}(t)^\intercal & = [q_0(t), q_1(t), \cdots, q_{N-1}(t)].
    \end{aligned}

Cartesian modulation equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cartesian modulation equations can be obtained from the polar ones. To do so, one can write

.. math::
    \begin{cases}
    \dfrac{\textrm{d} p_i}{\textrm{d}t} & = \dfrac{\textrm{d} a_i}{\textrm{d}t} \cos(\beta_i) - a_i \sin(\beta_i) \dfrac{\textrm{d} \beta_i}{\textrm{d}t}, \\
    \dfrac{\textrm{d} q_i}{\textrm{d}t} & = \dfrac{\textrm{d} a_i}{\textrm{d}t} \sin(\beta_i) + a_i \cos(\beta_i) \dfrac{\textrm{d} \beta_i}{\textrm{d}t}.
    \end{cases}

Then, by identification, one necessarily has

.. math::
    \begin{cases}
    f_{p_i}(\boldsymbol{p}, \boldsymbol{q}) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}) \cos(\beta_i) - f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}) \sin(\beta_i), \\
    f_{q_i}(\boldsymbol{p}, \boldsymbol{q}) & = f_{a_i}(\boldsymbol{a}, \boldsymbol{\beta}) \sin(\beta_i) + f_{\beta_i}(\boldsymbol{a}, \boldsymbol{\beta}) \cos(\beta_i),
    \end{cases}

in order to write the modulation equations in cartesian coordinates

.. math::
    \begin{cases}
    \dfrac{\textrm{d}}{dt} p_0(t) & = f_{p_0}(\boldsymbol{p}, \boldsymbol{q}), \\
    \dfrac{\textrm{d}}{dt} q_0(t) & = f_{q_0}(\boldsymbol{p}, \boldsymbol{q}), \\
    & \vdots \\
    \dfrac{\textrm{d}}{dt} p_{N-1}(t) & = f_{p_{N-1}}(\boldsymbol{p}, \boldsymbol{q}), \\
    \dfrac{\textrm{d}}{dt} q_{N-1}(t) & = f_{q_{N-1}}(\boldsymbol{p}, \boldsymbol{q}).
    \end{cases}

Jacobian matrix
^^^^^^^^^^^^^^^

As done previously for polar coordinates, let's introduce the vectors of cartesian coordinates and modulation equations as

.. math::
    \boldsymbol{x}^{(\textrm{c})\intercal} & = [p_0, q_0, \cdots, p_{N-1}, q_{N-1}], \\
    \boldsymbol{f}^{(\textrm{c})\intercal} & = [f_{p_0}, f_{q_0}, \cdots, f_{p_{N-1}}, f_{q_{N-1}}].

Then one can write

.. math::
    \dfrac{\textrm{d} \boldsymbol{x}^{(\textrm{c})}}{\textrm{d}t} = \textrm{J}^{(\textrm{c})} \boldsymbol{x}^{(\textrm{c})},

where we introduced the Jacobian matrix 

.. math::
    \textrm{J}^{(\textrm{c})}
    = \dfrac{\partial \boldsymbol{f}^{(\textrm{c})} }{ \partial \boldsymbol{x}^{(\textrm{c})} }
    = \begin{bmatrix}
    \frac{\partial f_{p_0}}{\partial p_0} & \frac{\partial f_{p_0}}{\partial q_0} & \cdots & \frac{\partial f_{p_0}}{\partial q_{N-1}} \\
    \frac{\partial f_{q_0}}{\partial p_0} & \frac{\partial f_{q_0}}{\partial q_0} & \cdots & \frac{\partial f_{q_0}}{\partial q_{N-1}} \\
    \vdots & \vdots & \ddots & \vdots \\
    \frac{\partial f_{q_{N-1}}}{\partial p_0} & \frac{\partial f_{q_{N-1}}}{\partial q_0} & \cdots & \frac{\partial f_{q_{N-1}}}{\partial q_{N-1}}
    \end{bmatrix}.

Note that there are no constraints related to an oscillator's amplitude being 0 here. 
This cartesian coordinates approach therefore allows to investigate how the stability of a steady state solution is affected by a perturbation from an oscillator who's amplitude is 0 in that steady state solution.

The stability analysis with :math:`\textrm{J}^{(\textrm{c})}` is carried out as described previously with :math:`\textrm{J}^{(\textrm{p})}`.
