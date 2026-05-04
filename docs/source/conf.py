# Configuration file for the Sphinx documentation builder.
import os
import sys
sys.path.insert(0, os.path.abspath('../../'))

# Project information
project = 'OSCILATE'
copyright = '2025, Vincent MAHE'
author = 'Vincent MAHE'
release = '0.0.0'

# General configuration
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinx.ext.autosummary',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'sphinx_copybutton',
    'sphinx_prompt',
    'myst_parser',
    'sphinxcontrib.bibtex',
]

# Template configuration
templates_path = ['_templates']
exclude_patterns = []

# Generate autosummary
autosummary_generate = True
autosummary_imported_members = True

# Intersphinx configuration
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'sympy': ('https://docs.sympy.org/latest/', None),
}

# Enable math support in MyST
myst_enable_extensions = [
    "amsmath",
    "dollarmath",
]

# Options for HTML output
html_theme = "sphinx_rtd_theme"
html_static_path = ['_static']
html_theme_options = {
    'navigation_depth': 6,  # default is 4
}

# Bibliography configuration
bibtex_bibfiles         = ['refs_oscilate.bib']
bibtex_default_style    = 'unsrt'  
bibtex_reference_style  = 'label'

# Types
autodoc_typehints = "none"  # This hides type hints from the function signature

# Define standalone files
exclude_patterns = [
    'mms/dyn_sys.rst',
    'mms/mms_application.rst',
    'mms/steady_state.rst',
    'mms/stability.rst',
]