# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function, Rational, sin, cos, sqrt, solve
from sympy.physics.vector.printing import init_vprinting
init_vprinting(use_latex=True, forecolor='White') # Initialise latex printing 
from MMS import MMS

# Parameters and variables
nP, nt, mu, eta, Lm, b, T1  = symbols(r"n_p, n_t, \mu, \eta, \Lambda_m, b, T_1", real=True, positive=True)
Lc, alpha1, alpha3, x4 = symbols(r"\Lambda_c, \alpha_1, \alpha_3, x_4", real=True)
t = symbols('t') # Replaces the rotor's angular position \vartheta.
z1 = Function(r'\zeta_1', real=True)(t) # Mode 1 coordinate (phase-opposition mode)
z2 = Function(r'\zeta_2', real=True)(t) # Mode 2 coordinate (unison mode)

# Dynamical system
dz1 = z1.diff(t)
dz2 = z2.diff(t)
ddz1 = z1.diff(t,2)
ddz2 = z2.diff(t,2)

f1 = -Lm**-1 * \
     (
         (Lm*dz1 - nt**2*(1+nt**2)*z1*z2) * (mu*nP**2*Lc*z2)
         + 2*mu*nt**2*Lm*dz1*(z1*dz1 + z2*dz2)
         + 6*eta*alpha1*alpha3*(z1*dz1**2 + 2*z2*dz1*dz2 + z1*dz2**2 + z1**2*ddz1
                                + 2*z1*z2*ddz2 + z2**2*ddz1)
         - 2*x4*(z1**3 + 3*z1*z2**2)
         + b*dz1
       )

f2 = -Lm**-1 * \
     (
         (Lc + Lm*dz2 - Rational(1,2)*nt**2*(1+nt**2)*(z1**2 + z2**2)) * (mu*nP**2*Lc*z2)
         + 2*mu*nt**2*(Lc + Lm*dz2)*(z1*dz1 + z2*dz2)
         - Lc*mu*nP**2*Rational(1,2)*nt**2*(1+nt**2)*(3*z1**2*z2 + z2**3)
         + Lc*mu*nt**2*(1+nt**2)*(2*z1*dz1*dz2 + z2*dz1**2 + z2*dz2**2)
         + 6*eta*alpha1*alpha3*(2*z1*dz1*dz2 + z2*dz1**2 + z2*dz2**2 
                                + 2*z1*ddz1*z2 + ddz2*z1**2 + ddz2*z2**2)
         - 2*x4*(3*z1**2 *z2 + z2**3)
         + b*dz2
       )

fF1 = -Lm**-1 * (Lm*dz1 - nt**2*(1+nt**2)*z1*z2)
fF2 = -Lm**-1 * (Lc + Lm*dz2 - Rational(1,2)*nt**2*(1+nt**2)*(z1**2 + z2**2))

Eq1 = ddz1 + nP**2*z1 - f1
Eq2 = ddz2 + nP**2*z2 - f2

kwargs_dyn = dict(F=T1, fF=[fF1, fF2])
dyn = MMS.Dynamical_system(t, [z1, z2], [Eq1, Eq2], [nP, nP], **kwargs_dyn)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1  # Order of the expansions
omega_ref      = nP # Reference frequency
ratio_omegaMMS = 2

param_to_scale = (mu, alpha3, x4, b, T1)
scaling        = (1 , 1,      1,  1, 1)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling,
                                 ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Substitute parameters
cc = symbols(r"c_c", real=True, positive=True)
ct, cp = symbols(r"c_t, c_p", real=True)
mut, alpha3t, x4t = param_scaled[0:3]
def_param = [(cp, 3*(x4t + 2*nP**2*eta*alpha1*alpha3t)),
             (ct, Rational(1,2) * Lc*mut*nP**2*nt**2*(1+nt**2)),
             (cc, mut*nt**4)
             ]
sub_nt = [(nt, nP*sqrt(Lm))]
sub_param = [(x4t, solve(def_param[0][1] - def_param[0][0], x4t)[0]),
             (Lc*mut*nt, solve(def_param[1][1] - def_param[1][0], Lc*mut*nt)[0]),
             (mut*nt, solve(def_param[2][1] - def_param[2][0], mut*nt)[0])]

a, beta = ss.coord.a, ss.coord.beta
collect_f = [a[0]*a[1]**2*sin(beta[0]-beta[1]),
             a[0]**2*a[1]*sin(beta[0]-beta[1]),
             a[0]*a[1]**2*cos(beta[0]-beta[1]),
             a[0]**2*a[1]*cos(beta[0]-beta[1]),
             a[0]**3, a[1]**3, a[0]*a[1]**2, a[0]**2*a[1]**2]
ss.sol.fa    = [fa_ix.subs(sub_param+sub_nt).simplify().expand().collect(collect_f) for fa_ix in ss.sol.fa]
ss.sol.fbeta = [fbeta_ix.subs(sub_param+sub_nt).simplify().expand().collect(collect_f) for fbeta_ix in ss.sol.fbeta]

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_forced(solve_dof=solve_dof)
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-2])

# Stability analysis 
kwargs_stab = dict(coord="cartesian", eigenvalues=True, bifurcation_curves=True, analyse_blocks=True, kwargs_bif=dict(var_a=True, var_sig=True, solver=solve))
ss.stability_analysis(**kwargs_stab)
