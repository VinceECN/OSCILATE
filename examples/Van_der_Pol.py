# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function
from sympy.physics.vector.printing import init_vprinting
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from MMS import MMS

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True, positive=True)
mu           = symbols(r'\mu', real=True, positive=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + mu*(x**2 - 1 )*x.diff(t) 
dyn = MMS.Dynamical_system(t, x, Eq, omega0, F=0)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 3      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 1      # Look for a solution around omega_ref

param_to_scale = (mu,)
scaling        = (1 ,)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)


# %%
