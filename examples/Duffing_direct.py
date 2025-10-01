# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function
from sympy.physics.vector.printing import init_vprinting, vlatex
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from MMS import MMS

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True,positive=True)
gamma        = symbols(r'\gamma',real=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + gamma*x**3 + c*x.diff(t)
dyn = MMS.Dynamical_system(t, x, Eq, omega0, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 3      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 1      # Look for responses around the reference frequency (Default behaviour)

param_to_scale = (gamma, F , c )
scaling        = (1    , Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_forced(solve_dof=solve_dof)
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])

# Stability analysis
ss.stability_analysis(coord="polar", eigenvalues=True, bifurcation_curves=True)

# Plot the steady state results
# -----------------------------
import numpy as np

# Set parameters' numerical values
a0 = np.linspace(1e-10, 1.2, 1000)

dic_numpy = dict(
    omega0 = (omega0, 1),
    c      = (c, 1e-2),
    gamma  = (gamma, 0.2),
    a      = (ss.coord.a[0], a0),
    )

F_val     = 1e-2
omega_val = 1.05

# Compute and plot the frequency-response curves (FRC)
dic_FRC = dic_numpy | dict(F=(dyn.forcing.F, F_val)) # Parameters for the FRC
FRC     = MMS.numpise_FRC(mms, ss, dyn, dic_FRC)
kwargs  = dict(phase_name=vlatex(ss.sol.cos_phase[0].args[0]),  # Plot parameters
               amp_name=vlatex(ss.coord.a[0]))
ss.plot_FRC(FRC, **kwargs)

# Compute and plot the amplitude-response curves (ARC)
dic_ARC = dic_numpy | dict(omega=(mms.omega, omega_val)) # Parameters for the ARC
ARC     = MMS.numpise_ARC(mms, ss, dyn, dic_ARC)
kwargs  = dict(phase_name=vlatex(ss.sol.cos_phase[0].args[0]), # Plot parameters
               amp_name=vlatex(ss.coord.a[0]))
ss.plot_ARC(ARC, **kwargs)


# %%
