# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function, sympify
from sympy.physics.vector.printing import init_vprinting
init_vprinting(use_latex=True, forecolor='White') 
from oscilate import MMS

# Parameters and variables
omega0, c = symbols(r'\omega_0, c', real=True, positive=True)
mu        = symbols(r'\mu', real=True, positive=True)
t         = symbols('t')
x         = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + mu*(x**2 - 1 )*x.diff(t) 
dyn = MMS.dyn_sys.Dynamical_system(t, x, Eq, omega0, F=0)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 2      # Order of the expansions
omega_ref      = sympify(1) # Reference frequency
ratio_omegaMMS = 1      # Look for a solution around omega_ref

sigma0          = symbols(r"\sigma_0", real=True) # Detuning of omega0 wrt 1
detunings       = [sigma0] # Detuning of the oscillator's frequency
sub_sigma0      = [(sigma0, omega0-1)] # Only for display purposes

param_to_scale = (mu, sigma0)
scaling        = (1 , 1     )
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

kwargs_mms = dict(ratio_omegaMMS=ratio_omegaMMS, ratio_omega_osc=[1], detunings=detunings)
mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, **kwargs_mms)

# Application of the MMS
mms.apply_MMS(rewrite_polar="all")

# Transient analysis
solve_dof = 0
mms.solve_transient(solve_dof=solve_dof)

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Steady state analysis - amplitude and frequency
ss.solve_LC(solve_dof=solve_dof)

# Steady state analysis - Stability analysis
J   = ss.Jacobian_polar() # Jacobian of the slow time system (in polar coordinates)
JSS = J.applyfunc(lambda expr: expr.subs(ss.sub.sub_solve_LC).subs(ss.sub.sub_scaling_back)) # Jacobian evaluated on the steady state solution
eigvals = list(JSS.eigenvals().keys())                    # Eigenvalues
eigvecs = list([item[2][0] for item in JSS.eigenvects()]) # Eigenvectors

# Plot the steady state results
# -----------------------------

# Set parameters' numerical values
import numpy as np
omega0 = 0.8
param = [(sigma0, 1-omega0),
         (mu, 0.3)]

# Limit cycle
LC = MMS.visualisation.Limit_cycle(mms, ss, param)
fig = LC.plot_PP(c="k")
ax = fig.get_axes()[0]

# Transient response - internal trajectory
param_IC = [(list(mms.sol_transient.IC["a"].values())[0], 0.5),
            (list(mms.sol_transient.IC["beta"].values())[0], 0)]
param_transient = param + [(mms.t, np.linspace(0, 60, 10000))] + param_IC
TR_in = MMS.visualisation.Transient_response(mms, param_transient)
ax.plot(TR_in.x, TR_in.dxdt, c="tab:blue")
ax.plot(TR_in.x[0], TR_in.dxdt[0], marker="o", mfc="tab:blue", mec="none", ms=4)

# Transient response - external trajectory
param_IC = [(list(mms.sol_transient.IC["a"].values())[0], 3.5),
            (list(mms.sol_transient.IC["beta"].values())[0], 0)]
param_transient = param + [(mms.t, np.linspace(0, 60, 10000))] + param_IC
TR_ex = MMS.visualisation.Transient_response(mms, param_transient)
ax.plot(TR_ex.x, TR_ex.dxdt, c="tab:red")
ax.plot(TR_ex.x[0], TR_ex.dxdt[0], marker="o", mfc="tab:red", mec="none", ms=4)

# Time signal (internal trajectory and LC)
figt = TR_in.plot_time()
axt = figt.get_axes()[0]
axt.axhline(LC.a, c="k", lw=0.5)
axt.axhline(-LC.a, c="k", lw=0.5)

# %%
