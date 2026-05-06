# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function, Rational
from sympy.physics.vector.printing import init_vprinting, vlatex
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from oscilate import MMS

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True,positive=True)
gamma        = symbols(r'\gamma',real=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq  = x.diff(t,2) + omega0**2*x + gamma*x**3 + c*x.diff(t) # Equation (unforced)
dyn = MMS.Dynamical_system(t, x, Eq, omega0, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = Rational(1,3) # Look for a superharmonic resonance of order 3

param_to_scale = (gamma, F, c)
scaling        = (1    , 0, 1)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_oscillator(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS(orders_polar="all")

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])
ss.solve_forced(solve_dof=solve_dof)

# Plot the steady state results
# -----------------------------
import numpy as np

# Set parameters' numerical values
param_FRC = [(omega0, 1),
             (c, 1e-2),
             (gamma, 0.2),
             (ss.coord.a[0], np.linspace(1e-10, 1.2, 1000)),
             (dyn.forcing.F, 0.5)]

# Frequency response
BBC = MMS.visualisation.Backbone_curve(mms, ss, dyn, param_FRC)
FRC = MMS.visualisation.Frequency_response_curve(mms, ss, dyn, param_FRC, bif=False)
figs = FRC.plot(ss=ss, bbc=BBC)
[fig.get_axes()[0].set_xlim(0.32, 0.38) for fig in figs ]

# %%
