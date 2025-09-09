The proposed code allows the application of the Method of Multiple Scales (MMS) to nonlinear equations and systems of coupled nonlinear equations.
MMS.py is the MMS solver. It contains 3 main classes:
- Dynamical_system : the dynamical system considered
- Multiple_scales_system : the system obtained after applying the MMS to the dynamical system
- Stead_state : the MMS results evaluated at steady state and (if computed) the system's response and its stability. Plotting functions 

test_MMS.py is a script containing several application examples of the MMS solver. Among the examples are
- The Duffing oscillator
- Coupled Duffing oscillators
- Coupled nonlinear oscillators with quadratic nonlinearities
- Parametrically excited oscillators
- Hard forcing of a Duffing oscillator
- Subharmonic response of 2 coupled centrifugal pendulum modes

sympy_functions.py contains additional functions that are not directly related to the MMS but which are used MMS.py.

Results are returned as sympy expressions.
They can be printed using LaTeX if the code is run in an appropriate interactive Window. 
It is the case VS Code's interactive Window or Spyder's IPython consol.

Methods of Steady_state also allow to evaluate sympy results for given numerical values of system parameters and plot them.


This code is provided as-is and has been tested on a limited number of nonlinear systems. 
Other test cases might trigger bugs or unexpected behavior that I am not yet aware of.
If you encounter any issues, find a bug, or have suggestions for improvements, please feel free to:
- Open an issue on the GitHub repository (if applicable).
- Contact me directly at [vincent.mahe@ec-nantes.fr].

Your feedback is highly appreciated!

Vincent MAHE

## License
This project is licensed under the **Apache License 2.0** – see the LICENSE file for details.

