# Reachability Analysis for Switched Systems

This project implements reachability analysis techniques for pwa systems using matrix zonotopes and constraint zonotopes.

## Prerequisites

Before running the code in this project, you need to install the following packages:

1. **CORA 2020** - COntinuous Reachability Analyzer
   - Download from: [https://tumcps.github.io/CORA/](https://tumcps.github.io/CORA/)
   - CORA is a toolbox for reachability analysis of continuous and hybrid systems

2. **zonoLAB** - Zonotope Laboratory
   - Download from: [https://github.com/ESCL-at-UTD/zonoLAB](https://github.com/ESCL-at-UTD/zonoLAB)
   - zonoLAB is a MATLAB toolbox for set-based computation with zonotopes

## Installation

1. Download and install CORA 2020 and zonoLAB-main from the links provided above
2. Add both toolboxes to your MATLAB path:
   ```matlab
   addpath(genpath('path/to/CORA'));
   addpath(genpath('path/to/zonoLAB'));
   ```
3. Clone or download this project to your local machine

## License

This project is provided for academic and research purposes.

## Acknowledgements

This work builds upon the CORA and zonoLAB toolboxes for set-based analysis.
