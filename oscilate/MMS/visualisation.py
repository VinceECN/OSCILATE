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
from sympy.physics.vector.printing import vlatex

#%% Classes and functions
class Frequency_response_curve:
    """
    Evaluate the frequency response curves (FRC) and bifurcation curves (if computed) for given numerical parameters. This transforms the sympy expressions to numpy arrays.
    They can then be plotted.

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The MMS results evaluated at steady state.
    dyn : Dynamical_system
        The initial dynamical system.
    param : list[tuple]
        A list whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.

    bbc : bool, optional
        Evaluate the backbone curve. 
        Default is `True`.
    forced : bool, optional
        Evaluate the forced response. 
        Default is `True`.
    bif : bool, optional
        Evaluate the bifurcation curves. 
        Default is `True`.
    """

    def __init__(self, mms, ss, dyn, param, bbc=True, forced=True, bif=True):
        
        # Information
        print("Converting sympy FRC expressions to numpy")

        # Construct a dictionary of substitutions
        param_dic = {}
        for ii in range(len(param)):
            if param[ii][0] == ss.coord.a[ss.sol_forced.solve_dof]:
                param_dic["a"] = param[ii]
            elif param[ii][0] == dyn.forcing.F:
                param_dic["F"] = param[ii]
            else:
                param_dic[f"param_{ii}"] = param[ii]
                
        # Initialisation
        self.a        = param_dic["a"][1]
        self.param    = param_dic
        self.omegaMMS = numpise_omegaMMS(mms, param_dic)

        # Evaluation of the FRC
        if bbc:
            self.omega_bbc = numpise_omega_bbc(mms, ss, param_dic)
        if forced:
            self.omega = numpise_omega_FRC(mms, ss, param_dic)
            self.phase = []
            for omegai in self.omega:
                self.phase.append(numpise_phase(mms, ss, dyn, param_dic, omegai, self.param["F"][1]))
        if bif:
            self.omega_bif = numpise_omega_bif(mms, ss, param_dic)
            if isinstance(self.omega_bif, np.ndarray):
                self.phase_bif = numpise_phase(mms, ss, dyn, param_dic, self.omega_bif, self.param["F"][1])
            elif isinstance(self.omega_bif, list):
                self.phase_bif = []
                for omegai in self.omega_bif:
                    self.phase_bif.append(numpise_phase(mms, ss, dyn, param_dic, omegai, self.param["F"][1]))
    
    def plot(self, **kwargs):
        r"""
        Plots the frequency response curves (FRC), both frequency-amplitude and frequency-phase.
        Also includes the stability information if given.

        Parameters
        ----------
        ss : Steady_state, optional
            Steady state object. Used to name the axis labels.
        
        Returns
        -------
        fig1 : Figure
            The amplitude plot :math:`a(\omega)`.
        fig2 : Figure
            The phase plot :math:`\beta(\omega)`.
        """

        # Extract the FRC data
        a         = self.__dict__.get("a",          np.full(10, np.nan))
        omega_bbc = self.__dict__.get("omega_bbc",  np.full_like(a, np.nan))
        omega     = self.__dict__.get("omega",      [np.full_like(a, np.nan)])
        omegaMMS  = self.__dict__.get("omegaMMS",   np.nan)
        phase     = self.__dict__.get("phase",      [np.full_like(a, np.nan)])
        omega_bif = self.__dict__.get("omega_bif",  [np.full_like(a, np.nan)])
        phase_bif = self.__dict__.get("phase_bif",  [np.full_like(a, np.nan)])
        
        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        if "ss" in kwargs.keys():
            ss = kwargs.get("ss")
            amp_name   = kwargs.get("amp_name", vlatex(ss.coord.a[ss.sol_forced.solve_dof]))
            phase_name = kwargs.get("phase_name", vlatex(ss.sol_forced.cos_phase[0].args[0]))
        else:
            amp_name   = kwargs.get("amp_name", "amplitude")
            phase_name = kwargs.get("phase_name", "phase")

        xlim = kwargs.get("xlim", [coeff*omegaMMS for coeff in (0.9, 1.1)])
        if np.isnan(xlim).any():
            xlim = [None, None]

        # FRC - amplitude 
        fig1, ax = plt.subplots(**fig_param)
        if isinstance(omega_bbc, np.ndarray): 
            ax.plot(omega_bbc, a, c="tab:grey", lw=0.7)
        ax.axvline(omegaMMS, c="k")
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
        ax.axvline(omegaMMS, c="k")
        ax.axhline(0.5, c="k", lw=0.7)
        [ax.plot(omegai, phasei/np.pi, c="tab:blue") for (omegai, phasei) in zip(omega, phase)]
        [ax.plot(omegai, phasei/np.pi, c="tab:red", lw=0.7) for (omegai, phasei) in zip(omega_bif, phase_bif)]
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"${} \; [\pi]$".format(phase_name))
        
        # Return
        return fig1, fig2

class Amplitude_response_curve:
    """
    Evaluate the amplitude response curves (ARC) for given numerical parameters. 
    'Amplitude' refers to the forcing amplitude.
    This transforms the sympy expressions to numpy arrays.
    They can then be plotted.

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The MMS results evaluated at steady state.
    dyn : Dynamical_system
        The initial dynamical system.
    param : list[tuple]
        A list whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.
    """

    def __init__(self, mms, ss, dyn, param):
        
        # Information
        print("Converting sympy ARC expressions to numpy")

        # Construct a dictionary of substitutions
        param_dic = {}
        for ii in range(len(param)):
            if param[ii][0] == ss.coord.a[ss.sol_forced.solve_dof]:
                param_dic["a"] = param[ii]
            elif param[ii][0] == mms.omega:
                param_dic["omega"] = param[ii]
            else:
                param_dic[f"param_{ii}"] = param[ii]
                
        # Initialisation
        self.a        = param_dic["a"][1]
        self.param    = param_dic
        self.omegaMMS = numpise_omegaMMS(mms, param_dic)

        # Evaluation of the ARC
        self.F     = numpise_F_ARC(mms, ss, self.param)
        if isinstance(self.F, np.ndarray):
            self.phase = numpise_phase(mms, ss, dyn, self.param, self.param["omega"][1], self.F)
        elif isinstance(self.F, list):
            self.phase = []
            for Fi in self.F:
                self.phase.append(numpise_phase(mms, ss, dyn, self.param, self.param["omega"][1], Fi))


    def plot(self, **kwargs):
        r"""
        Plots the amplitude-response curves (ARC), both forcing amplitude-amplitude and forcing amplitude-phase.

        Parameters
        ----------
        ss : Steady_state, optional
            Steady state object. Used to name the axis labels.

        Returns
        -------
        fig1 : Figure
            The amplitude plot :math:`a(F)`.
        fig2 : Figure
            The phase plot :math:`\beta(F)`.
        """

        # Extract the FRC data and keyword arguments
        a     = self.__dict__.get("a", np.full(10, np.nan))
        F     = self.__dict__.get("F", np.full_like(a, np.nan))
        phase = self.__dict__.get("phase", np.full_like(a, np.nan))
        
        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        if "ss" in kwargs.keys():
            ss = kwargs.get("ss")
            amp_name   = kwargs.get("amp_name", vlatex(ss.coord.a[ss.sol_forced.solve_dof]))
            phase_name = kwargs.get("phase_name", vlatex(ss.sol_forced.cos_phase[0].args[0]))
        else:
            amp_name   = kwargs.get("amp_name", "amplitude")
            phase_name = kwargs.get("phase_name", "phase")

        if isinstance(F, np.ndarray):
            xlim = kwargs.get("xlim", [0, np.nanmax(F)])
        elif isinstance(F, list):
            xlim = kwargs.get("xlim", [0, np.nanmax(np.hstack(F))])
        if np.isinf(xlim).any():
            xlim = [None, None]

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

def numpise_omegaMMS(mms, param):
    r"""
    Numpise the frequency around which a solution is sought.

    Parameters
    ----------
    mms: Multiple_scales_system
    param: dict
        See :func:`~MMS.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omegaMMS: float
        Numpised MMS frequency.
    """
    omegaMMS  = sfun.sympy_to_numpy(mms.omegaMMS, param)
    return omegaMMS

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
    omega_bbc  = sfun.sympy_to_numpy(rescale(ss.sol_bbc.omega, mms), param)
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
    omega = [np.real(sfun.sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.sol_forced.sigma]
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
    omega_bif = [np.real(sfun.sympy_to_numpy(mms.omegaMMS + rescale(mms.eps*sigmai, mms), param)) for sigmai in ss.sol_forced.stab.bif_sigma]
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
    sin_phase = sfun.sympy_to_numpy( rescale(ss.sol_forced.sin_phase[1], mms), param_phase)
    cos_phase = sfun.sympy_to_numpy( rescale(ss.sol_forced.cos_phase[1], mms), param_phase)
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
    if not isinstance(ss.sol_forced.F, list):
        F = sfun.sympy_to_numpy(rescale(mms.eps**mms.forcing.f_order * ss.sol_forced.F, mms), param)
    else:
        F = [sfun.sympy_to_numpy(rescale(mms.eps**mms.forcing.f_order * Fi, mms), param) for Fi in ss.sol_forced.F]
    return F