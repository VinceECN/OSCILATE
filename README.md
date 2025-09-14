# Nonlinear systems considered
This package allows the application of the **Method of Multiple Scales** (MMS) to nonlinear equations and systems of $N$ coupled nonlinear equations of the form

$$
\begin{cases}
        \ddot{x}_0 + \omega_0^2 x_0 & = f_0(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t), \\
        & \vdots \\
        \ddot{x}_{N-1} + \omega_{N-1}^2 x_{N-1} & = f_{N-1}(\boldsymbol{x}, \dot{\boldsymbol{x}}, \ddot{\boldsymbol{x}}, t).
        \end{cases}
$$

The $x_i$ ($i=0,...,N-1$) are the oscillators' coordinates, $\omega_i$ are their natural frequencies, $\boldsymbol{x}$ is the vector containing all the oscillators' coordinates, $t$ is the time, $\dot{(\bullet)}$ denotes a time-derivative $d(\bullet)/dt$, $f_i$ is a function which can contain:
- Linear terms in $x_i$, $\dot{x}_i$ or $\ddot{x}_i$, typically those that will be considered small in the MMS
- Weak coupling terms in $x_j$, $\dot{x}_j$ or $\ddot{x}_j$, $j\neq i$
- Weak nonlinear terms. Only polynomial nonlinearities are supported. Taylor expansions are performed if nonlinearities are not polynomial.
- Forcing, which can be hard (at first order) or weak (small). Harmonic and parametric forcing are supported.

Internal resonance relations among oscillators can be specified in a second step by expressing the $\omega_i$ as a function of a reference frequency. Detuning can also be introduced during this step.

# Solver
`MMS.py` is the MMS solver. It contains 3 main classes:
- `Dynamical_system` : the dynamical system considered
- `Multiple_scales_system` : the system obtained after applying the MMS to the dynamical system
- `Stead_state` : the MMS results evaluated at steady state and (if computed) the system's response and its stability. 

# Tests and examples
`test_MMS.py` is a script containing several application examples of the MMS solver. Among the examples are
- The Duffing oscillator
- Coupled Duffing oscillators
- Coupled nonlinear oscillators with quadratic nonlinearities
- Parametrically excited oscillators
- Hard forcing of a Duffing oscillator
- Subharmonic response of 2 coupled centrifugal pendulum modes

# Additional functions
`sympy_functions.py` contains additional functions that are not directly related to the MMS but which are used in `MMS.py`.

# Outputs
Results are returned as sympy expressions.
They can be printed using $\LaTeX$ if the code is run in an appropriate interactive Window. 
It is the case with VS Code's interactive Window or Spyder's IPython consol.

Methods of `Steady_state` also allow to evaluate sympy results for given numerical values of system parameters and plot them.

# Disclaimer
This code is provided as-is and has been tested on a limited number of nonlinear systems. 
Other test cases might trigger bugs or unexpected behavior that I am not yet aware of.
If you encounter any issues, find a bug, or have suggestions for improvements, please feel free to:
- Open an issue on the GitHub repository (if applicable).
- Contact me directly at [vincent.mahe@ec-nantes.fr].

Your feedback is highly appreciated!

Vincent MAHE

# License
This project is licensed under the **Apache License 2.0** – see the LICENSE file for details.

