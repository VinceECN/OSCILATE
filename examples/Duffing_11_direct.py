# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function, solve, cos
from sympy.physics.vector.printing import init_vprinting
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from oscilate import MMS

# Parameters and variables
omega0, F, c0, c1       = symbols(r'\omega_0, F, c_0, c_1', real=True, positive=True)
gamma0, gamma1, gamma01 = symbols(r'\gamma_0, \gamma_1, \gamma_{01}', real=True)
t  = symbols('t')
x0 = Function(r'x_0', real=True)(t)
x1 = Function(r'x_1', real=True)(t)

# Dynamical system
Eq0 = x0.diff(t,2) + omega0**2*x0 + gamma0*x0**3 + gamma01*x0*x1**2 + c0*x0.diff(t)
Eq1 = x1.diff(t,2) + omega0**2*x1 + gamma1*x1**3 + gamma01*x0**2*x1 + c1*x1.diff(t)
dyn = MMS.Dynamical_system(t, [x0, x1], [Eq0, Eq1], [omega0, omega0], fF=[0, 1], F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 1      # Look for a solution around omega_ref

param_to_scale = (gamma0, gamma1, gamma01, F , c0, c1)
scaling        = (1     , 1     , 1      , Ne, Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 1 # dof to solve for
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-2:])
ss.solve_forced(solve_dof=solve_dof)

# Stability analysis of the 1dof solution
ss.stability_analysis(coord="cartesian", eigenvalues=False, rewrite_polar=True)

# Computation of the coupled free solution
print("Manual computation of the coupled free solution")
Dbeta         = symbols(r"\Delta\beta", real=True) # Phase difference beta0-beta0
sub_phase     = [(ss.coord.beta[1], ss.coord.beta[0]+Dbeta)] # Substitute the phases by the phase difference
fa_free       = [fai.subs(ss.sub.sub_free+sub_phase) for fai in ss.sol.fa]
Dbeta_sol     = solve(fa_free[0], Dbeta) # Solution in terms of phase difference
chi           = symbols(r"\chi", real=True) # +/- symbol introduced to represent cos(2*Dbeta_sol)
sub_cos2Dbeta = [(cos(2*Dbeta), chi)]
fbeta_free    = [fbetai.subs(ss.sub.sub_free + sub_phase + sub_cos2Dbeta) for fbetai in ss.sol.fbeta]
a02_sol       = solve(fbeta_free[0], ss.coord.a[0]**2)[1]
sig_bbc       = solve(fbeta_free[1].subs(ss.coord.a[0]**2, a02_sol), mms.sigma)[0].subs(chi**2,1)

