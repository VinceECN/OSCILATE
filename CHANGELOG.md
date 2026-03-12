# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.4] - 2025-10-11
### Added
- Links in the README
- Separate readme files (md for GitHub, 2 .rst for PyPI and for the doc itself)
- Make The install guide separate from readme
- Reduce the readme size and point to the doc

## [1.0.5] - 2026-13-03
### Added in the documentation
- Figures to illustrate the oscillators
- Minor corrections in the equations and text
- A doc page to describe the MMS rather than relying on the classes doc
- The Van der Pol example was extended with a transient solution and plots
- The package tree is now in the readme

### Changes in the code
- The organisation is different: the MMS.py file was replaced by an MMS/ folder, containing the sub-modules dyn_sys.py, mms.py, steady_state.py and visualisation.py
- The sol class attributes were modified to be clearer, typically using an O when refering to lists containing solutions at each orders. Similar changes were performed in the code.
- The rewritting of solutions in polar coordinates was improved for a better visualisation
- Distinguish between actual simple bifurcation curves and trace curves

### General changes
- Add an AUTHORS file for HAL
- Add a SWHID to point to SoftwareHeritage
