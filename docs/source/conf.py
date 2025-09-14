# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# Project information
project = 'MMS_solver'
copyright = '2025, Vincent MAHE'
author = 'Vincent MAHE'
release = '1.0.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = []

# Enable math support in MyST
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
]

# Options for HTML output
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# MathJax configuration for bold symbols
# mathjax3_config = {
#     'loader': {
#         'load': ['[tex]/bm']
#     },
#     'tex': {
#         'packages': {'[+]': ['bm']},
#     },
# }