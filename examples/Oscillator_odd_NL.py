# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function
from sympy.physics.vector.printing import init_vprinting, vlatex
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from oscilate import MMS

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True,positive=True)
gamma3       = symbols(r'\gamma_3',real=True)
gamma5       = symbols(r'\gamma_5',real=True)
gamma7       = symbols(r'\gamma_7',real=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + gamma3*x**3 + gamma5*x**5 + gamma7*x**7 + c*x.diff(t)
dyn = MMS.Dynamical_system(t, x, Eq, omega0, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 3      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 1      # Look for responses around the reference frequency (Default behaviour)

param_to_scale = (gamma3, gamma5, gamma7, F , c )
scaling        = (1     , 1     , 1     , Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling)

# Application of the MMS
mms.apply_MMS(rewrite_polar=0)

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_forced(solve_dof=solve_dof)
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])

# Stability analysis
ss.stability_analysis_forced(coord="polar", eigenvalues=True, bifurcation_curves=True)

# Plot the steady state results
# -----------------------------

# Set parameters' numerical values
import numpy as np
param = [(omega0, 1),
         (c, 1e-2),
         (gamma3, 0.2),
         (gamma5, -0.3),
         (gamma7, 0.12),
         (ss.coord.a[0], np.linspace(1e-10, 2, 1000))]

# Frequency response
param_FRC = param + [(dyn.forcing.F, 1.4e-2)]
FRC = MMS.visualisation.Frequency_response_curve(mms, ss, dyn, param_FRC)
FRC.plot(ss=ss)

# Amplitude response
param_ARC = param + [(mms.omega, 1.02)]
ARC = MMS.visualisation.Amplitude_response_curve(mms, ss, dyn, param_ARC)
figs_ARC = ARC.plot(ss=ss)
[fig.get_axes()[0].set_xlim(0, 0.04) for fig in figs_ARC]


# import numpy as np
# import matplotlib.pyplot as plt

# # Set parameters' numerical values
# a0 = np.linspace(1e-10, 2, 1000)

# dic_numpy = dict(
#     omega0 = (omega0, 1),
#     c      = (c, 1e-2),
#     gamma3 = (gamma3, 0.2),
#     gamma5 = (gamma5, -0.3),
#     gamma7 = (gamma7, 0.12),
#     a      = (ss.coord.a[0], a0),
#     )

# F_val     = 1.4e-2
# omega_val = 1.02

# # Compute and plot the frequency-response curves (FRC)
# dic_FRC = dic_numpy | dict(F=(dyn.forcing.F, F_val)) # Parameters for the FRC
# FRC     = MMS.visualisation.numpise_FRC(mms, ss, dyn, dic_FRC)
# kwargs  = dict(phase_name=vlatex(ss.sol.cos_phase[0].args[0]),  # Plot parameters
#                amp_name=vlatex(ss.coord.a[0]))
# MMS.visualisation.plot_FRC(FRC, **kwargs)

# # Compute and plot the amplitude-response curves (ARC)
# dic_ARC = dic_numpy | dict(omega=(mms.omega, omega_val)) # Parameters for the ARC
# ARC     = MMS.visualisation.numpise_ARC(mms, ss, dyn, dic_ARC)
# kwargs  = dict(phase_name=vlatex(ss.sol.cos_phase[0].args[0]), # Plot parameters
#                amp_name=vlatex(ss.coord.a[0]))
# figs_ARC = MMS.visualisation.plot_ARC(ARC, **kwargs)
# [fig.get_axes()[0].set_xlim(0, 0.04) for fig in figs_ARC]

# %%
