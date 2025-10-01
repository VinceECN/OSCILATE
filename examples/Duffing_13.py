# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function
from sympy.physics.vector.printing import init_vprinting
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from MMS import MMS

# Parameters and variables
omega0, omega1, F, mu0, ma0 = symbols(r'\omega_0, \omega_1, F, \mu_0, \mu_1', real=True, positive=True)
alphas = [symbols(r'\alpha_{{{}}}'.format(ii), real=True) for ii in range(1,9)] # Nonlinear coefficients
Gam0, Gam1 = symbols(r"\Gamma_0, \Gamma_1", real=True) # Forcing coefficients

t  = symbols('t')
x0 = Function(r'x_0', real=True)(t)
x1 = Function(r'x_1', real=True)(t)

# Dynamical system
Eq0 = x0.diff(t,2) + omega0**2*x0 + 2*mu0*x0.diff(t) + sum([alphas[ii]   * x0**(3-ii)*x1**ii for ii in range(4)])
Eq1 = x1.diff(t,2) + omega1**2*x1 + 2*ma0*x1.diff(t) + sum([alphas[4+ii] * x0**(3-ii)*x1**ii for ii in range(4)])
fF  = [Gam0, Gam1]
dyn = MMS.Dynamical_system(t, [x0, x1], [Eq0, Eq1], [omega0,omega1], fF=fF, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 3      # Look for a solution around omega1

ratio_omega_osc = [1, 3] # Ratio between the oscillators' frequencies and the reference frequency
sigma1          = symbols(r"\sigma_1", real=True) # Detuning of oscillator 1
detunings       = [0, sigma1] # Detuning of the oscillators' frequency

param_to_scale = (*alphas,          F, mu0, ma0, sigma1)
scaling        = (*[1]*len(alphas), 1, 1,   1,   1)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

kwargs_mms = dict(ratio_omegaMMS=ratio_omegaMMS, ratio_omega_osc=ratio_omega_osc, detunings=detunings)
mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, **kwargs_mms)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)
ss.solve_bbc(solve_dof=1, c=param_scaled[9:11])
