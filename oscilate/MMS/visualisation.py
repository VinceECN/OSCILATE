# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module evaluates the symbolic expressions for given numerical parameters and allows to plot the resulting numerical results.
"""

#%% Imports and initialisation
from .. import sympy_functions as sfun
import numpy as np
import matplotlib.pyplot as plt
from .mms import rescale

#%% Classes and functions
def plot_FRC(FRC, **kwargs):
    r"""
    Plots the frequency response curves (FRC), both frequency-amplitude and frequency-phase.
    Also includes the stability information if given.

    Parameters
    ----------
    FRC : dict
        Dictionary containing the frequency response curves and the bifurcation curves.
    
    Returns
    -------
    fig1 : Figure
        The amplitude plot :math:`a(\omega)`.
    fig2 : Figure
        The phase plot :math:`\beta(\omega)`.
    """

    # Extract the FRC data
    a         = FRC.get("a", np.full(10, np.nan))
    omega_bbc = FRC.get("omega_bbc", np.full_like(a, np.nan))
    omega     = FRC.get("omega", [np.full_like(a, np.nan)])
    phase     = FRC.get("phase", [np.full_like(a, np.nan)])
    omega_bif = FRC.get("omega_bif", [np.full_like(a, np.nan)])
    phase_bif = FRC.get("phase_bif", [np.full_like(a, np.nan)])
    
    # Backbone for zero amplitude
    if isinstance(omega_bbc, np.ndarray):
        omega_bbc_0 = omega_bbc[0]
    elif isinstance(omega_bbc, float) or isinstance(omega_bbc, int):
        omega_bbc_0 = omega_bbc

    # Extract the keyword arguments
    fig_param  = kwargs.get("fig_param", dict())
    amp_name   = kwargs.get("amp_name", "amplitude")
    phase_name = kwargs.get("phase_name", "phase")
    xlim       = kwargs.get("xlim", [coeff*omega_bbc_0 for coeff in (0.9, 1.1)])
    if np.isnan(xlim).any():
        xlim = [None, None]

    # FRC - amplitude 
    fig1, ax = plt.subplots(**fig_param)
    if isinstance(omega_bbc, np.ndarray): 
        ax.plot(omega_bbc, a, c="tab:grey", lw=0.7)
    ax.axvline(omega_bbc_0, c="k")
    [ax.plot(omegai, a, c="tab:blue") for omegai in omega]
    [ax.plot(omegai, a, c="tab:red", lw=0.7) for omegai in omega_bif]
    
    ax.set_xlim(xlim)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"${}$".format(amp_name))
    ax.margins(y=0)

    if not isinstance(omega_bbc, np.ndarray):
        ax.set_ylim(0, 1.2*np.max(a))
    
    # FRC - phase
    fig2, ax = plt.subplots(**fig_param)
    ax.axvline(omega_bbc_0, c="k")
    ax.axhline(0.5, c="k", lw=0.7)
    [ax.plot(omegai, phasei/np.pi, c="tab:blue") for (omegai, phasei) in zip(omega, phase)]
    [ax.plot(omegai, phasei/np.pi, c="tab:red", lw=0.7) for (omegai, phasei) in zip(omega_bif, phase_bif)]
    
    ax.set_xlim(xlim)
    ax.set_xlabel(r"$\omega$")
    ax.set_ylabel(r"${} \; [\pi]$".format(phase_name))
    
    # Return
    return fig1, fig2

def plot_ARC(ARC, **kwargs):
    r"""
    Plots the amplitude-response curves (ARC), both forcing amplitude-amplitude and forcing amplitude-phase.

    Parameters
    ----------
    ARC : dict
        Dictionary containing the amplitude response curves.

    Returns
    -------
    fig1 : Figure
        The amplitude plot :math:`a(F)`.
    fig2 : Figure
        The phase plot :math:`\beta(F)`.
    """

    # Extract the FRC data and keyword arguments
    a     = ARC.get("a", np.full(10, np.nan))
    F     = ARC.get("F", np.full_like(a, np.nan))
    phase = ARC.get("phase", np.full_like(a, np.nan))
    
    # Extract the keyword arguments
    fig_param  = kwargs.get("fig_param", dict())
    amp_name   = kwargs.get("amp_name", "amplitude")
    phase_name = kwargs.get("phase_name", "phase")
    if isinstance(F, np.ndarray):
        xlim = kwargs.get("xlim", [0, np.nanmax(F)])
    elif isinstance(F, list):
        xlim = kwargs.get("xlim", [0, np.nanmax(np.hstack(F))])

    # ARC - amplitude 
    fig1, ax = plt.subplots(**fig_param)
    if isinstance(F, np.ndarray):
        ax.plot(F, a, c="tab:blue")
    elif isinstance(F, list):
        [ax.plot(Fi, a, c="tab:blue") for Fi in F]
    
    ax.set_xlim(xlim)
    ax.set_xlabel(r"$F$")
    ax.set_ylabel(r"${}$".format(amp_name))
    ax.margins(x=0, y=0)

    # ARC - phase
    fig2, ax = plt.subplots(**fig_param)
    ax.axhline(0.5, c="k", lw=0.7)
    if isinstance(F, np.ndarray):
        ax.plot(F, phase/np.pi, c="tab:blue")
    elif isinstance(F, list):
        [ax.plot(Fi, phasei/np.pi, c="tab:blue") for (Fi, phasei) in zip(F, phase)]
    
    ax.set_xlim(xlim)
    ax.set_xlabel(r"$F$")
    ax.set_ylabel(r"${} \; [\pi]$".format(phase_name))
    ax.margins(x=0)
    
    # Return
    return fig1, fig2

def numpise_FRC(mms, ss, dyn, param, bbc=True, forced=True, bif=True):
    r"""
    Evaluate the frequency-response and bifurcation curves at given numerical values.
    This transforms the sympy expressions to numpy arrays.

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The MMS results evaluated at steady state.
    dyn : Dynamical_system
        The initial dynamical system.
    param : dict
        A dictionary whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.

        The key of the amplitude vector must be ``"a"``.
        The key of the forcing amplitude must be ``"F"``.
    bbc : bool, optional
        Evaluate the backbone curve. 
        Default is `True`.
    forced : bool, optional
        Evaluate the forced response. 
        Default is `True`.
    bif : bool, optional
        Evaluate the bifurcation curves. 
        Default is `True`.

    Returns
    -------
    FRC : dict
        The frequency-response curves data.
    """

    # Information
    print("Converting sympy FRC expressions to numpy")

    # Initialisation
    a     = param.get("a")[1]
    F_val = param.get("F")[1]
    FRC   = {"a": a}

    # Evaluation of the FRC
    if bbc:
        FRC["omega_bbc"] = numpise_omega_bbc(mms, ss, param)
    if forced:
        FRC["omega"] = numpise_omega_FRC(mms, ss, param)
        FRC["phase"] = []
        for omegai in FRC["omega"]:
            FRC["phase"].append(numpise_phase(mms, ss, dyn, param, omegai, F_val))
    if bif:
        FRC["omega_bif"] = numpise_omega_bif(mms, ss, param)
        if isinstance(FRC["omega_bif"], np.ndarray):
            FRC["phase_bif"] = numpise_phase(mms, ss, dyn, param, FRC["omega_bif"], F_val)
        elif isinstance(FRC["omega_bif"], list):
            FRC["phase_bif"] = []
            for omegai in FRC["omega_bif"]:
                FRC["phase_bif"].append(numpise_phase(mms, ss, dyn, param, omegai, F_val))

    return FRC

def numpise_ARC(mms, ss, dyn, param):
    r"""
    Evaluate the amplitude-response curves at given numerical values. 
    This transforms the sympy expressions to numpy arrays. 

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The MMS results evaluated at steady state.
    dyn : Dynamical_system
        The initial dynamical system.
    param : dict
        A dictionnary whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.

        The key of the amplitude vector must be ``"a"``.
        The key of the angular frequency must be ``"omega"``.

    Returns
    -------
    ARC: dict
        The amplitude-response curves data.
    """
    
    # Information
    print("Converting sympy ARC expressions to numpy")

    # Initialisation
    a         = param.get("a")[1]
    omega_val = param.get("omega")[1]
    ARC       = {"a": a}

    # Evaluation of the FRC
    ARC["F"]     = numpise_F_ARC(mms, ss, param)
    if isinstance(ARC["F"], np.ndarray):
        ARC["phase"] = numpise_phase(mms, ss, dyn, param, omega_val, ARC["F"])
    elif isinstance(ARC["F"], list):
        ARC["phase"] = []
        for Fi in ARC["F"]:
            ARC["phase"].append(numpise_phase(mms, ss, dyn, param, omega_val, Fi))

    return ARC

def numpise_omega_bbc(mms, ss, param):
    r"""
    Numpise the backbone curve's frequency :math:`\omega_{\textrm{bbc}}`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega_bbc: numpy.ndarray
        Numpised backbone curve's frequency.
    """
    omega_bbc  = sfun.sympy_to_numpy(rescale(ss.sol.omega_bbc, mms), param)
    return omega_bbc

def numpise_omega_FRC(mms, ss, param):
    r"""
    Numpise the forced response's frequency :math:`\omega`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega: numpy.ndarray
        Numpised forced response's frequency.
    """
    omega = [np.real(sfun.sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.sol.sigma]
    return omega

def numpise_omega_bif(mms, ss, param):
    r"""
    Numpise the bifurcation curves' frequency :math:`\omega_{\textrm{bif}}`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega_bif: list of numpy.ndarray
        Numpised bifurcation curves' frequency.
    """
    omega_bif = [np.real(sfun.sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.stab.bif_sigma]
    return omega_bif

def numpise_phase(mms, ss, dyn, param, omega, F):
    r"""
    Numpise the phase :math:`\beta_i`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    dyn: Dynamical_system
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.
    omega: numpy.ndarray
        The frequency array.
    F: numpy.ndarray
        The forcing amplitude array.

    Returns
    -------
    phase: numpy.ndarray
        Numpised phase.
    """
    
    param_phase = param | dict(omega=(mms.omega, omega), F=(dyn.forcing.F, F))
    sin_phase = sfun.sympy_to_numpy( rescale(ss.sol.sin_phase[1], mms), param_phase)
    cos_phase = sfun.sympy_to_numpy( rescale(ss.sol.cos_phase[1], mms), param_phase)
    phase = np.arctan2(sin_phase, cos_phase)

    return phase

def numpise_F_ARC(mms, ss, param):
    r"""
    Numpise the forced response's forcing amplitude :math:`F`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    F: numpy.ndarray
        Numpised forced response's forcing amplitude.
    """
    if not isinstance(ss.sol.F, list):
        F = sfun.sympy_to_numpy(rescale(mms.eps**mms.forcing.f_order * ss.sol.F, mms), param)
    else:
        F = [sfun.sympy_to_numpy(rescale(mms.eps**mms.forcing.f_order * Fi, mms), param) for Fi in ss.sol.F]
    return F