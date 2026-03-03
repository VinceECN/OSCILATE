# -*- coding: utf-8 -*-

#%% Imports and initialisation
from sympy import symbols, Function, sympify, solve, dsolve, cos, sin
from sympy.physics.vector.printing import init_vprinting, vlatex
init_vprinting(use_latex=True, forecolor='White') 
from oscilate import MMS

# Parameters and variables
omega0, c = symbols(r'\omega_0, c', real=True, positive=True)
mu        = symbols(r'\mu', real=True, positive=True)
t         = symbols('t')
x         = Function(r'x', real=True)(t)

# Dynamical system
Eq = x.diff(t,2) + omega0**2*x + mu*(x**2 - 1 )*x.diff(t) 
dyn = MMS.Dynamical_system(t, x, Eq, omega0, F=0)

# Initialisation of the MMS sytem
eps            = symbols(r"\epsilon", real=True, positive=True) # Small parameter epsilon
Ne             = 2      # Order of the expansions
omega_ref      = sympify(1) # Reference frequency
ratio_omegaMMS = 1      # Look for a solution around omega_ref

sigma0          = symbols(r"\sigma_v", real=True) # Detuning of omega0 wrt 1
detunings       = [sigma0] # Detuning of the oscillator's frequency
sub_sigma0      = [(sigma0, omega0-1)] # Only for display purposes

param_to_scale = (mu, sigma0)
scaling        = (2 , 2     )
param_scaled, sub_scaling = MMS.scale_parameters(param_to_scale, scaling, eps)

kwargs_mms = dict(ratio_omegaMMS=ratio_omegaMMS, ratio_omega_osc=[1], detunings=detunings)
mms = MMS.Multiple_scales_system(dyn, eps, Ne, omega_ref, sub_scaling, **kwargs_mms)

# Application of the MMS
mms.apply_MMS(rewrite_polar="all")

# Transient analysis - slow time solutions
Eqa    = mms.coord.at[0].diff(t)                   - mms.sol.fa[0]    # Equation on a
Eqbeta = mms.coord.at[0]*mms.coord.betat[0].diff(t) - mms.sol.fbeta[0] # Equation on beta
ai    = symbols(r"a_i", real=True, positive=True)     # Initial amplitude
betai = symbols(r"\beta_i", real=True, positive=True) # Initial phase
ICa = {mms.coord.at[0].subs(t,0) : ai}           # Initial condition on a
ICbeta = {mms.coord.betat[0].subs(t,0) : betai}  # Initial condition on beta
a_sol    = dsolve(Eqa, mms.coord.at[0], ics=ICa).rhs.subs(mms.sub.sub_scaling_back).subs(sigma0+1, omega0)
beta_sol = dsolve(Eqbeta, mms.coord.betat[0], ics=ICbeta).rhs.subs(mms.coord.at[0], a_sol).doit().expand().simplify().subs(mms.sub.sub_scaling_back).subs(sigma0+1, omega0)

# Transient analysis - instantaneous oscillation frequency and time response
psi = Function(r"\psi", real=True, positive=True)(t) # Full phase
sub_psi = [(mms.omega*t, psi+mms.coord.betat[0])]     # Substitution from the relative phase beta to the full one psi
psi_sol = (mms.omega*t - beta_sol).subs([mms.sub.sub_omega]).expand().simplify() # From the beta solution to the psi one
omegaNL = psi_sol.diff(t).simplify().factor() # Instantaneous oscillation frequency
list_cos_sin = list(mms.sol.x[0].atoms(cos, sin)) # List of harmonics in the response
x_sol = mms.sol.x[0].expand().collect(list_cos_sin).subs(mms.sub.sub_t[:-1]).subs(mms.sub.sub_scaling_back).subs(sub_psi) # Time signal

# Evaluation at steady state
ss = MMS.Steady_state(mms)

# Steady state analysis - amplitude and frequency
aSS_sol     = solve(ss.sol.fa[0], ss.coord.a[0])[0] # Amplitude solution
sigSS_sol   = solve(ss.sol.fbeta[0], mms.sigma)[0].subs(ss.coord.a[0], aSS_sol) # Detuning solution
omegaSS_sol = (1 + eps*sigSS_sol).subs(mms.sub.sub_scaling_back).simplify() # Frequency solution
sub_SS_sol  = [(ss.coord.a[0], aSS_sol), (mms.sigma, sigSS_sol)] # Substitutions to steady state solutions
xSS_sol = mms.sol.x[0].expand().collect(list_cos_sin).subs(mms.sub.sub_t[:-1]).subs(ss.sub.sub_SS).subs(mms.sub.sub_scaling_back) # Time signal

# Steady state analysis - Stability analysis
J   = ss.Jacobian_polar() # Jacobian of the slow time system (in polar coordinates)
JSS = J.applyfunc(lambda expr: expr.subs(sub_SS_sol).subs(mms.sub.sub_scaling_back)) # Jacobian evaluated on the steady state solution
eigvals = list(JSS.eigenvals().keys()) # Eigenvalues
eigvecs = list([item[2][0] for item in JSS.eigenvects()]) # Eigenvectors

# Phase portrait - Numerical evaluation of the symbolic solutions
from oscilate.sympy_functions import sympy_to_numpy
import numpy as np

def get_x_v(omega0_num, mu_num, SS=False, **kwargs):
    """
    Retrieve the numerical expressions of the displacement x and velocity v from the symbolic results, for a given set of numerical parameters
    """

    sigma0_num = 1 - omega0_num # Detuning

    dic_time   = {
        "epsilon" : (mu, mu_num),
        "sigma0": (sigma0, sigma0_num),
        "omega0": (omega0, omega0_num)
    }

    if SS: # Limit cycle
        a_num     = sympy_to_numpy(aSS_sol, dic_time) # Amplitude response at steady state
        omega_num = sympy_to_numpy(omegaSS_sol, dic_time) # Frequency of the response at steady state
        t_num = np.linspace(0, 2*np.pi/omega_num, 10000)

        dic_a_beta = {
            "t": (t, t_num), 
            "a": (ss.coord.a[0], a_num),
            "beta": (ss.coord.beta[0], 0),
            "omega": (mms.omega, omega_num)
        }

        vSS_sol  = xSS_sol.diff(t) # Symbolic expression of the velocity at steady state
        x = sympy_to_numpy(xSS_sol, dic_time | dic_a_beta) # Displacement at steady state
        v = sympy_to_numpy(vSS_sol, dic_time | dic_a_beta) # Velocity at steady state

    else: # Transient trajectories
        ai_num, t_num  = list(map(kwargs.get,["ai","t"]))
        dic_time["t"]     = (t, t_num)
        dic_time["ai"]    = (ai, ai_num)
        dic_time["betai"] = (betai, 0)

        dic_a_psi = {
        "a": (mms.coord.at[0], sympy_to_numpy(a_sol, dic_time)),
        "psi": (psi, sympy_to_numpy(psi_sol, dic_time)),
        "da": (mms.coord.at[0].diff(t), sympy_to_numpy(a_sol.diff(t).simplify(), dic_time)),
        "dpsi": (psi.diff(t), sympy_to_numpy(psi_sol.diff(t).simplify(), dic_time))
        }

        v_sol = x_sol.diff(t)
        x = sympy_to_numpy(x_sol, dic_time | dic_a_psi)
        v = sympy_to_numpy(v_sol, dic_time | dic_a_psi)

    return x, v
omega0_num = 0.8
mu_num     = 0.3
xLC, vLC = get_x_v(omega0_num, mu_num, SS=True)
xLT, vLT = get_x_v(omega0_num, mu_num, SS=False, **dict(ai=0.5, t=np.linspace(0, 60, 10000)))
xUT, vUT = get_x_v(omega0_num, mu_num, SS=False, **dict(ai=3.5, t=np.linspace(0, 60, 10000)))

# Plot the phase portrait
import matplotlib.pyplot as plt

fig, ax = plt.subplots()
ax.plot(xLC, vLC, c="k", lw=2, zorder=3)
ax.plot(xLT, vLT, c="tab:blue")
ax.plot(xLT[0], vLT[0], marker="o", mfc="tab:blue", mec="none", ms=4)
ax.plot(xUT, vUT, c="tab:red")
ax.plot(xUT[0], vUT[0], marker="o", mfc="tab:red", mec="none", ms=4)
ax.set_xlabel(r"${}$".format(vlatex(x)))
ax.set_ylabel(r"${}$".format(vlatex(x.diff(t))))

