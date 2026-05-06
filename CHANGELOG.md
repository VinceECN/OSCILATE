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

## [1.0.6] - 2026-06-05
### Added in the documentation
- Added a doc link at the beginning of the README files
- Added the references section
- Detailed the complex form approach, with corrections to local hyperlinks and standalone reference directives
- Added a citation for the normal form approach with the Duffing in complex MMS form
- Clarified code output in the Duffing example
- Corrected the phase expression in the steady state part of the doc

### Changes in the code
- Implemented the complex form of the MMS
- Added `Multiple_scales_system` subclasses dedicated to the MMS form, either oscillator or complex
- Split `mms.sol` into `.sol` and `.sol_transient`; introduced `sol_forced`, `sol_bbc`, `sol_LC` in steady state
- Implemented computation of transient trajectories and limit cycles
- Added harmonic decomposition of solutions
- Split the backbone curve from the FRC in the visualisation module
- In the reconstruction of x, replace missing polar form terms from the asymptotic series by zero instead of skipping. Also, `order_polar` replaces `rewrite_polar`
- Changed the representation of x and z solutions
- Made visualisation classes more flexible and intuitive
- Removed the 3rd order polynomial solver from `steady_state.solve_F` (used for the superharmonic Duffing), as it was not working reliably
- Changed prints in `mms.py` and `visualisation.py`

### Changes in the examples
- Added a complex form example for the Duffing oscillator
- Extended examples with frequency plots for the superharmonic Duffing
- Added a time signal to the Van der Pol example
- Updated all examples to use `order_polar` (replacing `rewrite_polar`) and the split FRC/BBC visualisation
- Added SVG plots for application examples
- Reviewed and improved existing examples, added high-order nonlinearity example

### General changes
- Renamed `Duffing_SupH_plots` to `Duffing_supH_plots`
- Added more description to the GitHub repo