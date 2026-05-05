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
from typing import Union, TYPE_CHECKING
from numpy import ndarray

#%% Classes and functions
class Backbone_curve:
    """
    Evaluate the backbone curve for given numerical parameters. This transforms the sympy expressions to numpy arrays. They can then be plotted.

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

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        a           : np.ndarray
        omegaMMS    : float
        omega       : np.ndarray
        param       : dict

    def __init__(self, mms, ss, dyn, param):
        
        # Information
        print("Converting sympy BBC expressions to numpy")

        # Construct a dictionary of substitutions
        param_dic = {}
        for ii in range(len(param)):
            if param[ii][0] == ss.coord.a[ss.sol_bbc.solve_dof]:
                param_dic["a"] = param[ii]
            elif param[ii][0] == dyn.forcing.F:
                param_dic["F"] = param[ii]
            else:
                param_dic[f"param_{ii}"] = param[ii]
                
        # Initialisation
        self.a        = param_dic["a"][1]
        self.param    = param_dic
        self.omegaMMS = numpise_omegaMMS(mms, param_dic)

        # Evaluation of the bbc
        self.omega = numpise_omega_bbc(mms, ss, param_dic)
        if ss.sol_bbc.xmax != None:
            self.xmax = numpise_xmax_bbc(mms, ss, param_dic)
    
    def plot(self, **kwargs):
        r"""
        Plots the backbone curve.

        Parameters
        ----------
        ss : Steady_state, optional
            Steady state object. Used to name the axis labels.
        
        Returns
        -------
        fig : Figure
            The amplitude plot :math:`a(\omega)`.
        """

        # Extract the bbc data
        a        = self.__dict__.get("a",        np.full(10, np.nan))
        omega    = self.__dict__.get("omega",    np.full_like(a, np.nan))
        omegaMMS = self.__dict__.get("omegaMMS", np.nan)
        
        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        if "ss" in kwargs.keys():
            ss = kwargs.get("ss")
            amp_name   = kwargs.get("amp_name", vlatex(ss.coord.a[ss.sol_forced.solve_dof]))
        else:
            amp_name   = kwargs.get("amp_name", "amplitude")

        xlim = kwargs.get("xlim", [coeff*omegaMMS for coeff in (0.9, 1.1)])
        if np.isnan(xlim).any():
            xlim = [None, None]

        # Backbone curve 
        fig, ax = plt.subplots(**fig_param)
        ax.plot(omega, a, c="tab:grey", lw=0.7)
        ax.axvline(omegaMMS, c="k")
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$\omega_{\textrm{nl}}$")
        ax.set_ylabel(r"${}$".format(amp_name))
        ax.margins(y=0)

        # Return
        return fig

class Frequency_response_curve:
    """
    Evaluate the frequency response curves (FRC) and bifurcation curves (if computed) for given numerical parameters. This transforms the sympy expressions to numpy arrays. They can then be plotted.

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

    bif : bool, optional
        Evaluate the bifurcation curves. 
        Default is `True`.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        a           : np.ndarray
        omega       : list[ndarray]
        omegaMMS    : float
        omega_bif   : list[ndarray]
        param       : dict
        phase       : list[ndarray]
        phase_bif   : list[ndarray]

    def __init__(self, mms, ss, dyn, param, bif=True):
        
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
    
    def plot(self, bbc=None, **kwargs):
        r"""
        Plots the frequency response curves (FRC), both frequency-amplitude and frequency-phase.
        Also includes the stability information if computed, and the backbone curve if given.

        Parameters
        ----------
        bbc : Backbone_curve, optional
            The evaluated backbone curve. 
            Default is None.
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
        omega     = self.__dict__.get("omega",      [np.full_like(a, np.nan)])
        omegaMMS  = self.__dict__.get("omegaMMS",   np.nan)
        phase     = self.__dict__.get("phase",      [np.full_like(a, np.nan)])
        omega_bif = self.__dict__.get("omega_bif",  [np.full_like(a, np.nan)])
        phase_bif = self.__dict__.get("phase_bif",  [np.full_like(a, np.nan)])

        # Extract the backbone curve data if any
        if bbc != None:
            a_bbc     = bbc.__dict__.get("a"    ,  np.full_like(a, np.nan))
            omega_bbc = bbc.__dict__.get("omega",  np.full_like(a, np.nan))
            if not isinstance(omega_bbc, np.ndarray):
                omega_bbc = np.full_like(a_bbc, omega_bbc)
        
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
        if bbc != None: 
            ax.plot(omega_bbc, a_bbc, c="tab:grey", lw=0.7)
        ax.axvline(omegaMMS, c="k")
        [ax.plot(omegai, a, c="tab:blue") for omegai in omega]
        [ax.plot(omegai, a, c="tab:red", lw=0.7) for omegai in omega_bif]
        
        ax.set_xlim(xlim)
        ax.set_xlabel(r"$\omega$")
        ax.set_ylabel(r"${}$".format(amp_name))
        ax.margins(y=0)

        if bbc == None:
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

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        F           : Union[list[ndarray], ndarray]
        a           : np.ndarray
        omegaMMS    : float
        param       : dict
        phase       : list[ndarray]

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
    
class Transient_response:
    """
    Evaluate the transient response for given numerical parameters. 
    This transforms the sympy expressions to numpy arrays.
    They can then be plotted in the phase portrait and as time signals.

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    param : list[tuple]
        A list whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        a           : np.ndarray
        dxdt        : np.ndarray
        param       : dict
        psi         : np.ndarray
        solve_dof   : int
        t           : np.ndarray
        x           : np.ndarray

    def __init__(self, mms, param):
        
        # Information
        print("Converting sympy transient response expressions to numpy")

        # Construct a dictionary of substitutions
        param_dic = {}
        for ii in range(len(param)):
            if param[ii][0] == mms.t:
                param_dic["t"] = param[ii]
            elif param[ii][0] in mms.sol_transient.IC["a"].values():
                param_dic["ai"] = param[ii]
            elif param[ii][0] in mms.sol_transient.IC["a"].values():
                param_dic["betai"] = param[ii]
            else:
                param_dic[f"param_{ii}"] = param[ii]

        # Compute the slow time solutions
        a, psi, dadt, dpsidt = numpise_transient_slow_time(mms, param_dic)

        slow_sol           = dict()
        slow_sol["a"]      = (mms.coord.at[mms.sol_transient.solve_dof], a)
        slow_sol["psi"]    = (mms.coord.psi, psi)
        slow_sol["dadt"]   = (mms.coord.at[mms.sol_transient.solve_dof].diff(mms.t), dadt)
        slow_sol["dpsidt"] = (mms.coord.psi.diff(mms.t), dpsidt)

        # Compute the time signals
        x, dxdt = numpise_transient_trajectory(mms, param_dic | slow_sol)

        # Store the results
        self.param = param_dic
        self.t     = param_dic["t"][1]
        self.a     = a
        self.psi   = psi
        self.x     = x
        self.dxdt  = dxdt
        self.solve_dof = mms.sol_transient.solve_dof

    def plot_PP(self, c="tab:blue", **kwargs):
        r"""
        Plots the transient trajectory in the phase portrait.

        Parameters
        ----------
        c : str, optional
            The color of the line.

        Returns
        -------
        fig : Figure
            The phase portrait plot.
        """

        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        
        # Trajectory plot
        fig, ax = plt.subplots(**fig_param)
        ax.plot(self.x, self.dxdt, c=c)
        ax.plot(self.x[0], self.dxdt[0], marker="o", mfc=c, mec="none", ms=4)
        
        # Labels
        ax.set_xlabel(r"$x_{}$".format(self.solve_dof))
        ax.set_ylabel(r"$\dot{{x}}_{}$".format(self.solve_dof))

        # Return
        return fig
    
    def plot_time(self, c="tab:blue", **kwargs):
        r"""
        Plots the transient time response.

        Parameters
        ----------
        c : str, optional
            The color of the line.

        Returns
        -------
        fig : Figure
            The time plot.
        """

        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        
        # Time plot
        fig, ax = plt.subplots(**fig_param)
        ax.plot(self.t, self.x, c=c)
        
        # Labels
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$x_{}$".format(self.solve_dof))

        # Return
        return fig
    

class Limit_cycle:
    """
    Evaluate the limit cycle for given numerical parameters. 
    This transforms the sympy expressions to numpy arrays.
    They can then be plotted in the phase portrait and as time signals.

    Parameters
    ----------
    mms : Multiple_scales_system
        The MMS object.
    ss : Steady_state
        The SS object.
    param : list[tuple]
        A list whose values are tuples with 2 elements:

        1. The sympy symbol of a parameter,

        2. The numerical value(s) taken by that parameter.
    Npts: int, optional
        Number of time points.
        Default is 1000.
    """

    # Class-level annotations for pyreverse
    if TYPE_CHECKING:
        a           : float
        beta        : float
        dxdt        : np.ndarray
        param       : dict
        solve_dof   : int
        t           : np.ndarray
        x           : np.ndarray

    def __init__(self, mms, ss, param, Npts=1000):
        
        # Information
        print("Converting sympy transient response expressions to numpy")

        # Construct a dictionary of substitutions
        param_dic = {}
        for ii in range(len(param)):
            param_dic[f"param_{ii}"] = param[ii]

        # Compute the LC amplitude, phase and frequency.
        a, beta, omega = numpise_LC(mms, ss, param_dic)
        LC_sol          = dict()
        LC_sol["a"]     = (ss.coord.a[ss.sol_LC.solve_dof], a)
        LC_sol["beta"]  = (ss.coord.beta[ss.sol_LC.solve_dof], beta)
        LC_sol["omega"] = (ss.omega, omega)

        # Compute the time signals
        param_dic["t"] = (mms.t, np.linspace(0, 2*np.pi/omega, Npts))
        x, dxdt = numpise_LC_trajectory(mms, ss, param_dic | LC_sol)

        # Store the results
        self.param = param_dic
        self.t     = param_dic["t"][1]
        self.a     = a
        self.beta  = beta
        self.x     = x
        self.dxdt  = dxdt
        self.solve_dof = ss.sol_LC.solve_dof

    def plot_PP(self, c="tab:blue", lw=2, **kwargs):
        r"""
        Plots the transient trajectory in the phase portrait.

        Parameters
        ----------
        c : str, optional
            The color of the line.
        lw : float, optional
            The linewidth
            Default is 2.

        Returns
        -------
        fig : Figure
            The phase portrait plot.
        """

        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        
        # Trajectory plot
        fig, ax = plt.subplots(**fig_param)
        ax.plot(self.x, self.dxdt, c=c, lw=lw, zorder=10)
        
        # Labels
        ax.set_xlabel(r"$x_{}$".format(self.solve_dof))
        ax.set_ylabel(r"$\dot{{x}}_{}$".format(self.solve_dof))

        # Return
        return fig
    
    def plot_time(self, c="tab:blue", **kwargs):
        r"""
        Plots the transient time response.

        Parameters
        ----------
        c : str, optional
            The color of the line.

        Returns
        -------
        fig : Figure
            The time plot.
        """

        # Extract the keyword arguments
        fig_param  = kwargs.get("fig_param", dict())
        
        # Time plot
        fig, ax = plt.subplots(**fig_param)
        ax.plot(self.t, self.x, c=c)
        
        # Labels
        ax.set_xlabel(r"$t$")
        ax.set_ylabel(r"$x_{}$".format(self.solve_dof))

        # Return
        return fig

def numpise_omegaMMS(mms, param):
    r"""
    Numpise the frequency around which a solution is sought.

    Parameters
    ----------
    mms: Multiple_scales_system
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

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
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    omega_bbc: numpy.ndarray
        Numpised backbone curve's frequency.
    """
    omega_bbc  = sfun.sympy_to_numpy(rescale(ss.sol_bbc.omega, mms), param)
    return omega_bbc

def numpise_xmax_bbc(mms, ss, param):
    r"""
    Numpise the peak oscillator's amplitude :math:`x_{\textrm{max}}` on the backbone curve.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    xmax: numpy.ndarray
        Numpised peak amplitude on the backbone curve.
    """
    xmax  = sfun.sympy_to_numpy(rescale(ss.sol_bbc.xmax, mms), param)
    return xmax

def numpise_omega_FRC(mms, ss, param):
    r"""
    Numpise the forced response's frequency :math:`\omega`.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

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
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

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
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.
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
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

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

def numpise_transient_slow_time(mms, param):
    r"""
    Numpise the slow time transient response.

    Parameters
    ----------
    mms: Multiple_scales_system
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    a: numpy.ndarray
        Numpised transient amplitude.
    psi: numpy.ndarray
        Numpised transient absolute phase.
    dadt: numpy.ndarray
        Numpised time derivative of the transient amplitude.
    dpsidt: numpy.ndarray
        Numpised time derivative of the transient absolute phase.
    """

    a      = sfun.sympy_to_numpy(rescale(mms.sol_transient.a, mms), param)
    psi    = sfun.sympy_to_numpy(rescale(mms.sol_transient.psi, mms), param)
    dadt   = sfun.sympy_to_numpy(rescale(mms.sol_transient.a.diff(mms.t), mms), param)
    dpsidt = sfun.sympy_to_numpy(rescale(mms.sol_transient.psi.diff(mms.t), mms), param)

    return a, psi, dadt, dpsidt

def numpise_transient_trajectory(mms, param):
    r"""
    Numpise the transient oscillator's trajectory

    Parameters
    ----------
    mms: Multiple_scales_system
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    x: numpy.ndarray
        Numpised transient motion.
    dxdt: numpy.ndarray
        Numpised transient velocity.
    """

    x    = sfun.sympy_to_numpy(rescale(mms.sol_transient.x, mms), param)
    dxdt = sfun.sympy_to_numpy(rescale(mms.sol_transient.x.diff(mms.t), mms), param)

    return x, dxdt

def numpise_LC(mms, ss, param):
    r"""
    Numpise the limit cycle solution.

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state_system
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    a: float
        Numpised LC amplitude.
    beta: float
        Numpised LC initial phase.
    omega: float
        Numpised LC frequency.
    """

    a     = sfun.sympy_to_numpy(rescale(ss.sol_LC.a, mms), param)
    beta  = sfun.sympy_to_numpy(rescale(ss.sol_LC.beta, mms), param)
    omega = sfun.sympy_to_numpy(rescale(ss.sol_LC.omega, mms), param)

    return a, beta, omega

def numpise_LC_trajectory(mms, ss, param):
    r"""
    Numpise the oscillator's LC trajectory

    Parameters
    ----------
    mms: Multiple_scales_system
    ss: Steady_state_system
    param: dict
        See :func:`~oscilate.sympy_functions.sympy_to_numpy`.

    Returns
    -------
    x: numpy.ndarray
        Numpised LC motion.
    dxdt: numpy.ndarray
        Numpised LC velocity.
    """

    x    = sfun.sympy_to_numpy(rescale(ss.sol_LC.x, mms), param)
    dxdt = sfun.sympy_to_numpy(rescale(ss.sol_LC.x.diff(mms.t), mms), param)

    return x, dxdt