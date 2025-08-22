# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Application of the Method of Multiple Scales (MMS) to nonlinear equations and systems of coupled nonlinear equations.
Several examples are proposed below:
- The Duffing oscillator
- Coupled Duffing oscillators
- Coupled nonlinear oscillators with quadratic nonlinearities
- Parametrically excited oscillators
- Hard forcing of a Duffing oscillator
- Subharmonic response of 2 coupled centrifugal pendulum modes

Results are returned as sympy expressions.
They can be printed using LaTeX if the code is run in an appropriate interactive Window. 
It is the case VS Code's interactive Window or Spyder's IPython consol.

On top of applying the MMS, the code can
- Evaluate the results at steady-state
- Compute steady-state responses (if analytical solutions exist)
- Evaluate the stability of a computed steady-state response
- Change coordinates from polar to cartesian
- Evaluate steady-state responses for some given numerical parameter values and plot them
"""

#%% Imports and initialisation
import sympy as sy
from sympy import (symbols, Function, Symbol, cos, sin, exp, I, conjugate, re, im, Rational, fraction, solve, dsolve)
from sympy.physics.vector.printing import init_vprinting, vlatex
import MMS

# Initialise latex printing - only works in Interactive Windows backed by the Jupyter kernel (not in the regular Python REPL in VS Code)
init_vprinting(use_latex=True, forecolor='White') 

# Count time
import time
start = time.time() 

# Close open figures
import matplotlib.pyplot as plt
plt.close("all") 

# Force the reload of MMS - useful when developping in VS Code using the IPython consol
import importlib
importlib.reload(MMS) 

# display() function - useful when in debugger mode in Spyder
from IPython.display import display 

# print_py() function - useful to display outputs as python code
from sympy.printing import sstrrepr as print_py 

#%% MMS - Duffing (Nayfeh & Mook, 4.1, 4.1.1)

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
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency

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
ss.eval_sol_stability(coord="polar", eigenvalues=True, bifurcation_curves=True)

# Plot the steady state results
# -------------------
import numpy as np

# Set parameters
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
dic_FRC  = dic_numpy | dict(F=(dyn.forcing.F, F_val)) # Parameters for the FRC
FRC      = MMS.numpise_FRC(mms, ss, dyn, dic_FRC)
kwargs   = dict(phase_name=vlatex(ss.sol.cos_phase[0].args[0]), 
                amp_name=vlatex(ss.coord.a[0]))
ss.plot_FRC(FRC, **kwargs)

# Compute and plot the amplitude-response curves (ARC)
dic_ARC   = dic_numpy | dict(omega=(mms.omega, omega_val)) # Parameters for the ARC
ARC      = MMS.numpise_ARC(mms, ss, dyn, dic_ARC)
kwargs   = dict(phase_name=vlatex(ss.sol.cos_phase[0].args[0]), 
                amp_name=vlatex(ss.coord.a[0]))
ss.plot_ARC(ARC, **kwargs)

#%% MMS - coupled Duffing in 1:3 internal resonance (Nayfeh & Mook, 6.6, 6.6.2) 

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
fcoeff = [Gam0, Gam1]
dyn = MMS.Dynamical_system(t, [x0, x1], [Eq0, Eq1], [omega0,omega1], f_coeff=fcoeff, F=F)

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

# Computation of the coupled free solution -> not working as does not yield a polynomial in sigma. Could try the amplitude response instead
# Dbeta         = symbols(r"\Delta\beta", real=True) # Phase difference beta0-beta0
# sub_phase     = [(ss.coord.beta[1], ss.coord.beta[0]+Dbeta)] # Substitute the phases by the phase difference
# fa_free       = [fai.subs(ss.sub.sub_free+sub_phase) for fai in ss.sol.fa]
# Dbeta_sol     = solve(fa_free[0], Dbeta) # Solution in terms of phase difference
# chi           = symbols(r"\chi", real=True) # +/- symbol introduced to represent cos(2*Dbeta_sol)
# sub_cos2Dbeta = [(cos(Dbeta), chi)]
# fbeta_free    = [fbetai.subs(ss.sub.sub_free + sub_phase + sub_cos2Dbeta) for fbetai in ss.sol.fbeta]
# fbeta_free[0] = (fbeta_free[0]/ss.coord.a[0]).simplify().expand()
# a0_sols       = solve(fbeta_free[0], ss.coord.a[0]) # Choose solution [1] but keep in mind there is a +/- sign in front of the sqrt
# chi2          = symbols(r"\chi_2", real=True) # +/- symbol introduced to represent +/- sqrt()
# a0_sol        = Rational(1,2)*(sum(a0_sols) + chi2*(a0_sols[1] - a0_sols[0])).simplify()
# sig_bbc       = solve(fbeta_free[1].subs(ss.coord.a[0], a0_sol), mms.sigma)[0].subs([(chi**2,1), (chi2**2,1)]) # Not possible as it is not a polynomial in sigma 

#%% MMS - coupled Duffing in 1:1 internal resonance (not treated in Nayfeh & Mook)

# Parameters and variables
omega0, F, c0, c1       = symbols(r'\omega_0, F, c_0, c_1', real=True, positive=True)
gamma0, gamma1, gamma01 = symbols(r'\gamma_0, \gamma_1, \gamma_{01}', real=True)
t  = symbols('t')
x0 = Function(r'x_0', real=True)(t)
x1 = Function(r'x_1', real=True)(t)

# Dynamical system
Eq0 = x0.diff(t,2) + omega0**2*x0 + gamma0*x0**3 + gamma01*x0*x1**2 + c0*x0.diff(t)
Eq1 = x1.diff(t,2) + omega0**2*x1 + gamma1*x1**3 + gamma01*x0**2*x1 + c1*x1.diff(t)
dyn = MMS.Dynamical_system(t, [x0, x1], [Eq0, Eq1], [omega0, omega0], f_coeff=[0, 1], F=F)

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
ss.eval_sol_stability(coord="cartesian", eigenvalues=False)

# Computation of the coupled free solution
Dbeta         = symbols(r"\Delta\beta", real=True) # Phase difference beta0-beta0
sub_phase     = [(ss.coord.beta[1], ss.coord.beta[0]+Dbeta)] # Substitute the phases by the phase difference
fa_free       = [fai.subs(ss.sub.sub_free+sub_phase) for fai in ss.sol.fa]
Dbeta_sol     = solve(fa_free[0], Dbeta) # Solution in terms of phase difference
chi           = symbols(r"\chi", real=True) # +/- symbol introduced to represent cos(2*Dbeta_sol)
sub_cos2Dbeta = [(cos(2*Dbeta), chi)]
fbeta_free    = [fbetai.subs(ss.sub.sub_free + sub_phase + sub_cos2Dbeta) for fbetai in ss.sol.fbeta]
a02_sol       = solve(fbeta_free[0], ss.coord.a[0]**2)[1]
sig_bbc       = solve(fbeta_free[1].subs(ss.coord.a[0]**2, a02_sol), mms.sigma)[0].subs(chi**2,1)

#%% MMS - coupled quadratic oscillators in 1:2 internal resonance (Nayfeh & Mook, 6.5, 6.5.1)

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
fcoeff = [eta0, eta1]
dyn = MMS.Dynamical_system(t, [x0, x1], [Eq0, Eq1], [omega0, omega1], f_coeff=fcoeff, F=F)

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
a, beta = ss.coord.a, ss.coord.beta

Dbeta     = symbols(r"\Delta\beta", real=True) # Phase difference beta0-beta0
sub_phase = [(beta[0], beta[1]-Dbeta)] # Substitute the phases by the phase difference
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

#%% Parametrically excited Duffing (Nayfeh & Mook, 5.7.3)

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True,positive=True)
gamma        = symbols(r'\gamma',real=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + gamma*x**3 + c*x.diff(t)
f_coeff = -2*x # Parametric forcing
dyn = MMS.Dynamical_system(t, x, Eq, omega0, f_coeff=f_coeff, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 2      # Look for a solution around n*omega_ref

param_to_scale = (gamma, F , c )
scaling        = (1    , Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])
ss.solve_forced(solve_dof=solve_dof)

# Stability analysis
ss.eval_sol_stability(coord="polar", eigenvalues=True)

#%% Van der Pol's oscillator- (Not treated in Nayfeh and Mook, Rayleigh oscillator's is considered instead)

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

param_to_scale = (mu, F , c )
scaling        = (1 , Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Amplitude evolution
a_sol      = sy.dsolve(mms.sol.fa[0].expand() - mms.coord.a[0].diff(mms.tS[1]), mms.coord.a[0])[1].rhs
C1         = list(a_sol.atoms(Symbol).difference(mms.sol.fa[0].expand().atoms(Symbol)))[0]
a_sol_init = 1
C1_sol     = sy.solve(a_sol.subs(mms.tS[1],0) - a_sol_init, C1)[0]
a_sol      = a_sol.subs(C1, C1_sol)

#%% Rayleigh's oscillator (Nayfeh and Mook, 5.7.2)

# Parameters and variables
omega0, F, c = symbols(r'\omega_0, F, c', real=True,positive=True)
gamma        = symbols(r'\gamma',real=True)
t            = symbols('t')
x            = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + gamma*x.diff(t)**3 + c*x.diff(t) 
f_coeff = -2*x # Parametric forcing
dyn = MMS.Dynamical_system(t, x, Eq, omega0, f_coeff=f_coeff, F=F)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 1      # Order of the expansions
omega_ref      = omega0 # Reference frequency
ratio_omegaMMS = 2      # Look for a solution around 2*omega_ref

param_to_scale = (gamma, F , c )
scaling        = (1    , Ne, Ne)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])
ss.solve_forced(solve_dof=solve_dof)

# Stability analysis
ss.eval_sol_stability(coord="polar", eigenvalues=True)

#%% MMS - Duffing superharmonic resonance (Nayfeh and Mook, 4.1.2, 4.1.3)

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

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])
ss.solve_forced(solve_dof=solve_dof)

#%% MMS - Duffing subharmonic resonance (Nayfeh and Mook, 4.1.2, 4.1.4)

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
ratio_omegaMMS = 3      # Look for a subharmonic resonance of order 3

param_to_scale = (gamma, F, c)
scaling        = (1    , 0, 1)
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, ratio_omegaMMS=ratio_omegaMMS)

# Application of the MMS
mms.apply_MMS()

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Solve the evolution equations for a given dof
solve_dof = 0 # dof to solve for
ss.solve_bbc(solve_dof=solve_dof, c=param_scaled[-1])
ss.solve_forced(solve_dof=solve_dof)

#%% CPVA with 2 pendulum modes in 1:1 internal resonance (Mahe, 2022, Subharmonic CPVAs & On the dynamic stability and efficiency of CPVAs)

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

f_coeff1 = -Lm**-1 * (Lm*dz1 - nt**2*(1+nt**2)*z1*z2)
f_coeff2 = -Lm**-1 * (Lc + Lm*dz2 - Rational(1,2)*nt**2*(1+nt**2)*(z1**2 + z2**2))

Eq1 = ddz1 + nP**2*z1 - f1
Eq2 = ddz2 + nP**2*z2 - f2

kwargs_dyn = dict(F=T1, f_coeff=[f_coeff1, f_coeff2])
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
sub_nt = [(nt, nP*sy.sqrt(Lm))]
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
kwargs_stab = dict(coord="cartesian", eigenvalues=True, bifurcation_curves=True, analyse_blocks=True, kwargs_bif=dict(var_a=True, var_sig=True, solver=sy.solve))
ss.eval_sol_stability(**kwargs_stab)

#%% Display the execution time
end = time.time()
print('exécution : ', '%.1f'%(end - start),' secondes')