# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function
from sympy.physics.vector.printing import init_vprinting, vlatex
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from oscilate import MMS

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True,positive=True)
gamma        = symbols(r'\gamma',real=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + gamma*x**3 + c*x.diff(t)
dyn = MMS.Dynamical_system(t, x, Eq, omega0, F=F)
dyn.complex_form() # Construct the complex form of the system

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 5      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 1      # Look for responses around the reference frequency (Default behaviour)

param_to_scale = (gamma, F , c )
scaling        = (1    , Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

# Complex form of the MMS
mms = MMS.Multiple_scales_complex(dyn, eps, Ne, omega_ref, sub_scaling)
mms.apply_MMS(orders_polar="all") # Application of the MMS
ss = MMS.Steady_state(mms) # Steady state analysis
ss.solve_bbc(solve_dof=0, c=param_scaled[-1]) # Backbone curve computation

# Oscillator form of the MMS
mms_o = MMS.Multiple_scales_oscillator(dyn, eps, Ne, omega_ref, sub_scaling)
mms_o.apply_MMS(orders_polar="all") # Application of the MMS
ss_o = MMS.Steady_state(mms_o) # Steady state analysis
ss_o.solve_bbc(solve_dof=0, c=param_scaled[-1]) # Backbone curve computation

# Oscillator form of the MMS - order 1 
mms_o1 = MMS.Multiple_scales_oscillator(dyn, eps, 1, omega_ref, sub_scaling)
mms_o1.apply_MMS(orders_polar="all") # Application of the MMS
ss_o1 = MMS.Steady_state(mms_o1) # Steady state analysis
ss_o1.solve_bbc(solve_dof=0, c=param_scaled[-1]) # Backbone curve computation

# Plot the steady state results
# -----------------------------

# Set parameters' numerical values
import numpy as np
param = [(omega0, 1),
         (gamma, 0.5),
         (ss.coord.a[0], np.linspace(1e-10, 1.6, 1000))]

# Frequency response
param_FRC = param + [(dyn.forcing.F, 1e-2)]
BBC    = MMS.visualisation.Backbone_curve(mms, ss, dyn, param_FRC)
BBC_o  = MMS.visualisation.Backbone_curve(mms_o, ss_o, dyn, param_FRC)
BBC_o1 = MMS.visualisation.Backbone_curve(mms_o1, ss_o1, dyn, param_FRC)

# Plot
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(BBC.omega, BBC.xmax, label=f"complex form - Ne={Ne}")
ax.plot(BBC_o.omega, BBC_o.xmax, label=f"oscillator form - Ne={Ne}")
ax.plot(BBC_o1.omega, BBC_o1.xmax, label=f"oscillator form - Ne=1")
ax.axvline(BBC.omegaMMS, c="k")
ax.set_xlabel(r"$\omega_{\text{nl}}$")
ax.set_ylabel(r"$max(x(t))$")
ax.legend()
ax.set_ymargin(0)

# %%
