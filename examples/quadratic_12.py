# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function, solve, sin, cos, Rational
from sympy.physics.vector.printing import init_vprinting
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from oscilate import MMS

# Parameters and variables
omega0, omega1, F, mu0, ma0 = symbols(r'\omega_0, \omega_1, F, \mu_0, \mu_1', real=True, positive=True)
alpha0, alpha1 = symbols(r"\alpha_0, \alpha_1", real=True) # Nonlinear coefficients
eta0, eta1 = symbols(r"\eta_0, \eta_1", real=True) # Forcing coefficients

t  = symbols('t')
x0 = Function(r'x_0', real=True)(t)
x1 = Function(r'x_1', real=True)(t)

# Dynamical system
Eq0 = x0.diff(t,2) + omega0**2*x0 + 2*mu0*x0.diff(t) + alpha0*x0*x1
Eq1 = x1.diff(t,2) + omega1**2*x1 + 2*ma0*x1.diff(t) + alpha1*x0**2
fF = [eta0, eta1]
dyn = MMS.Dynamical_system(t, [x0, x1], [Eq0, Eq1], [omega0, omega1], fF=fF, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 2      # Look for a solution around omega1

ratio_omega_osc = [1, 2] # Ratio between the oscillators' frequencies and the reference frequency
sigma1          = symbols(r"\sigma_1", real=True) # Detuning of oscillator 1
detunings       = [0, sigma1] # Detuning of the oscillators' frequency

param_to_scale = (F, eta0, eta1, mu0, ma0, sigma1, alpha0, alpha1)
scaling        = (0, 1,    2,    1,   1,   1,      0,      0)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

kwargs_mms = dict(ratio_omegaMMS=ratio_omegaMMS, ratio_omega_osc=ratio_omega_osc, 
                  detunings=detunings, eps_pow_0=1)
mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, **kwargs_mms)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)
ss.solve_bbc(solve_dof=0, c=param_scaled[3:5])

# Computation of the coupled forced solution
print("Manual computation of the coupled forced solution")

a, beta = ss.coord.a, ss.coord.beta

Dbeta     = symbols(r"\Delta\beta", real=True) # Phase difference Dbeta = 2beta0-beta1
sub_phase = [(beta[0], Rational(1,2)*(beta[1]-Dbeta))] # Substitute the phases by the phase difference
fa        = [fai.subs(sub_phase) for fai in ss.sol.fa]
fbeta     = [fbetai.subs(sub_phase) for fbetai in ss.sol.fbeta]

dic_fa0    = fa[0].collect(sin(Dbeta), evaluate=False)
dic_fbeta0 = fbeta[0].collect(cos(Dbeta), evaluate=False)
Eq_a1      = (dic_fa0[1]/dic_fa0[sin(Dbeta)])**2 + (dic_fbeta0[1]/dic_fbeta0[cos(Dbeta)])**2 - 1
a1_2_sol   = solve(Eq_a1, a[1]**2)[0] # Solution a1**2

Dbeta_sol = solve(fa[0], Dbeta) # Phase difference solution
chi       = symbols(r"\chi", real=True) # +/- symbol introduced to represent the 2 phase difference solutions
sub_Dbeta = [(sin(Dbeta), sin(Dbeta_sol[0])), (cos(Dbeta), chi*cos(Dbeta_sol[0]))]

dic_fa1    = fa[1].collect(sin(beta[1]), evaluate=False)
dic_fbeta0 = fbeta[1].collect(cos(beta[1]), evaluate=False)
Eq_a0      = (((dic_fa1[1]/dic_fa1[sin(beta[1])])**2 + (dic_fbeta0[1]/dic_fbeta0[cos(beta[1])])**2 - 1).subs(sub_Dbeta).subs(ss.coord.a[1]**2, a1_2_sol)).simplify()
a0_2_sol   = solve(Eq_a0, a[0]**2) # Solution a0**2

