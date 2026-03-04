# -*- coding: utf-8 -*-
"""
Started on Tue Feb 15 17:25:59 2022

@author: Vincent MAHE

Analyse systems of coupled nonlinear equations using the Method of Multiple Scales (MMS).
This sub-module evaluates the stability of a steady state response.
"""

#%% Imports and initialisation
from sympy import (exp, I, conjugate, re, im, Rational, 
                   symbols, Symbol, Function, solve, dsolve,
                   cos, sin, tan, srepr, sympify, simplify, 
                   zeros, det, trace, eye, Mod, sqrt)
from sympy.simplify.fu import TR5, TR8, TR10
from .. import sympy_functions as sfun
import numpy as np
import itertools
import matplotlib.pyplot as plt