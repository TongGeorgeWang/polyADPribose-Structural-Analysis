[![DOI](https://zenodo.org/badge/761409703.svg)](https://zenodo.org/doi/10.5281/zenodo.12823464)


![image](https://github.com/TongGeorgeWang/polyADPribose-Structural-Analysis/assets/160785251/db72b17f-fa75-41c8-8be3-dde28659a9a2)

**Overview** <br />
All analysis performed on structural ensembles of poly(ADP-ribose) that are included in: <br />
Length-dependent Intramolecular Coil-to-Globule Transition in Poly(ADP-ribose) Induced by Cations
Tong Wang, Kush Coshic, Mohsen Badiee, Aleksei Aksimentiev, Lois Pollack, Anthony K. L. Leung
BioRxiv, 2024
doi: https://doi.org/10.1101/2023.10.25.564012

Each folder represents a separate analysis step. <br />
EnsembleOptimizationMethod: EOM analysis to determine an ensemble of conformers that agrees with experimental small-angle scattering data  <br />
OrderParameters: computation of orientational correlation function, tortuosity (chain twist), and number of base stacking events in a structural ensemble <br />
HierarchicalClustering: dividing structural ensembles into ~3 distinct categories based on size parameters <br />
ClassAveraging_SpectralClustering: class averaging of disordered structural ensembles (identifying a handful of characteristic conformations) using graph theory and spectral clustering. <br />
MolecularDynamics: unique parameter set used to model PAR and example contact analysis outputs





**Hardware Requirements** <br />
This analysis should only require a standard computer that is able to run the commonly available software below, with enough memory to complete low to moderate intensity computing tasks. <br />
Please consult NAMD for specific requirements about running MD simulations, as those can deviate from this claim. https://www.ks.uiuc.edu/Research/namd/2.14/ug/node99.html

**Software Requrements** <br />
NAMD 2.14: To run MD simulations. <br />
VMD 1.9.4a43: To visualize MD simulation systems. <br />
Tcl 8.0 or Python 3.7.3: Custom code provided to analyze the MD trajectories. <br />
VMD (for tcl) or MDAnalysis 1.1.1 (for Python): For reading the trajectory (DCD) files. <br />
EOM in ATSAS v3: Structural ensemble analysis, using MD trajectories as input. <br />
PyMOL v2.5.5: To visualize structures and for structural alignment during class averaging.<br />
MATLAB R2021: Custom code provided for hierarchical clustering, computing structural order parameters, and class averaging via spectral clustering. <br />
  Special MATLAB dependency:  <br />
  SpectralClustering.m (Ingo Bürk, areslp, 2011) <br />
  After EOM, the CRYSOL executable (from ATSAS) is utilized to compute scattering profiles. <br /> 
The MATLAB code for running EOM may require a Windows OS, as it calls on remote DOS commands. <br />
The specific versions of each software that were used for analyses is included. Future versions will likely also work. Please see manuscript references for citations for these software entities. 

**Instructions for use** <br />
You need only install the above software and clone/download the folders in this repository. Separate, more detailed README files are provided in each folder where appropriate, wherein more detailed instructions are provided. <br />
For a computer with an i5 or later processor and 16gb or more of RAM, each script should take 30 minutes or less, with most processes being completed in under 2 minutes. Please consult NAMD for an idea of how long MD simulations will take for a given computer setup. https://www.ks.uiuc.edu/Research/namd/2.14/ug/node99.html    

**Demos** <br />
Example input/output data featuring 15mer poly(ADP-ribose) in 100mM NaCl are included in each corresponding folder. Please contact the corresponding authors of the manuscript if the full dataset on all examined constructs is desired. 

**License** <br />
This project is covered under the Apache 2.0 License.




