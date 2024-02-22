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





**Hardware Requirements** <br />
This analysis should only require a standard computer that is able to run the commonly available software below, with enough memory to complete low to moderate intensity computing tasks.

**Software Requrements** <br />
VMD 1.9.4a43: To visualize and construct MD simulation systems. <br />
Tcl 8.0 or Python 3.7.3: Custom code provided to analyze the MD trajectories. <br />
VMD (for tcl) or MDAnalysis 1.1.1 (for Python): For reading the trajectory (DCD) files. <br />
EOM v2.1 in ATSAS v3: Structural ensemble analysis, using MD trajectories as input. <br />
MATLAB R2021: Custom code provided for hierarchical clustering, computing structural order parameters, and class averaging via spectral clustering. <br />
PyMOL v2.5.5: To visualize structures and for structural alignment during class averaging.<br />

The specific versions of each software that were used for analyses is included. Future versions will likely also work. 

**Instructions for use** <br />
You need only install the above software and clone/download the folders in this repository. Separate, more detailed README files are provided in each folder where appropriate, wherein more detailed instructions are provided. 

**Demos** <br />
Example input/output data featuring 15mer poly(ADP-ribose) in 100mM NaCl are included in each corresponding folder. 

**License** <br />
This project is covered under the Apache 2.0 License.




