Nonlinear systems considered
============================

The **OSCILATE** (Oscillators' nonlinear analysis through Symbolic ImplementATion of the mEthod of multiple scales) project allows the application of the **Method of Multiple Scales** (MMS) to a nonlinear equation or systems of :math:`N` coupled nonlinear equations of the form:

.. math::
    \begin{cases}
        \ddot{x}_0 + \omega_0^2 x_0 & = f_0(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t), \\
        & \vdots \\
        \ddot{x}_{N-1} + \omega_{N-1}^2 x_{N-1} & = f_{N-1}(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t).
    \end{cases}

The :math:`x_i(t)` (:math:`i=0,...,N-1`) are the oscillators' coordinates,

.. math::
    \boldsymbol{x}(t)^\intercal = [x_0(t), x_1(t), \cdots, x_{N-1}(t)]

is the vector containing all the oscillators' coordinates (the :math:`^\intercal` denotes the transpose),
:math:`\omega_i` are their natural frequencies,
:math:`t` is the time and
:math:`\dot{(\bullet)} = \textrm{d}(\bullet)/\textrm{d}t` denotes a time-derivative.

The :math:`f_i` are functions which can contain:

- **Weak linear terms** in :math:`x_i`, :math:`\dot{x}_i`, or :math:`\ddot{x}_i`.
- **Weak linear coupling terms** involving :math:`x_j`, :math:`\dot{x}_j`, or :math:`\ddot{x}_j` with :math:`j \neq i`.
- **Weak nonlinear terms**. Taylor expansions are performed to approximate nonlinear terms as *polynomial nonlinearities*.
- **Forcing terms**:

  - Can be hard (appearing at leading order) or weak (small).
  - Primarily *harmonic*, e.g., :math:`F \cos(\omega t)`, where :math:`F` and :math:`\omega` are the forcing amplitude and frequency, respectively.
  - Modulated by any function (constant, linear, or nonlinear) to model *parametric* forcing (e.g., :math:`x_i(t) F \cos(\omega t)`).

Internal resonance relations among oscillators can be specified in a second step by expressing the :math:`\omega_i` as a function of a reference frequency.
Detuning can also be introduced during this step.

Overview
========

The package associated with the **OSCILATE** project is called ``oscilate``.
It contains two modules:

- The ``MMS`` module is the MMS solver.
- The ``sympy_functions`` module contains additional functions that are not directly related to the MMS but which are used in ``MMS``.

Solver
------

``MMS`` contains 3 main classes:

- ``Dynamical_system``: the dynamical system considered.
- ``Multiple_scales_system``: the system obtained after applying the MMS to the dynamical system.
- ``Steady_state``: the MMS results evaluated at steady state and (if computed) the system's response and its stability.

These classes are described in detail in the `documentation <https://vinceECN.github.io/OSCILATE/>`_.

Examples
--------

Application examples are proposed in the `documentation <https://vinceECN.github.io/OSCILATE/>`_. They include:

- The Duffing oscillator.
- Coupled Duffing oscillators.
- Coupled nonlinear oscillators with quadratic nonlinearities.
- Parametrically excited oscillators.
- Hard forcing of a Duffing oscillator.
- Subharmonic response of 2 coupled centrifugal pendulum modes.

Outputs
-------

Results are returned as SymPy expressions.
They can be printed using :math:`\LaTeX` if the code is run in an appropriate interactive window. Here are possibilities:

* `VS Code's interactive Window <https://code.visualstudio.com/docs/python/jupyter-support-py>`_
* `Jupyter notebook <https://jupyter.org/>`_
* `Spyder's IPython console <https://docs.spyder-ide.org/current/panes/ipythonconsole.html>`_

SymPy expressions can also be printed as unformatted :math:`\LaTeX` using:

.. code-block:: python

    print(vlatex(the_expr))

Methods of ``Steady_state`` also allow evaluating SymPy results for given numerical values of system parameters and plotting them.

Documentation
=============

A full documentation is available `here <https://vinceECN.github.io/OSCILATE/>`_.

Citation
========

Please cite this package when using it. See the Citation section of the `documentation <https://vinceECN.github.io/OSCILATE/>`_ for details.
A regular entry and a LaTeX/BibTeX users entry are given.
A paper describing this work is currently in publication and will become the preferred citation once published. For now, please cite this repository.

Installation guide
==================

To install the ``oscilate`` package, refer to the Installation guide section of the `documentation <https://vinceECN.github.io/OSCILATE/>`_.

Disclaimer
==========

This code is provided as-is and has been tested on a limited number of nonlinear systems.
Other test cases might trigger bugs or unexpected behavior that I am not yet aware of.
If you encounter any issues, find a bug, or have suggestions for improvements, please feel free to:

- Open an issue on the GitHub repository (if applicable).
- Propose a solution.
- Contact me directly at vincent.mahe@ec-nantes.fr.

Your feedback is highly appreciated!

Vincent MAHÉ

License
=======

This project is licensed under the **Apache License 2.0** – see the LICENSE file for details.
