/////////////
This folder contains all the scripts used to class average disordered PAR ensembles in this study. An example dataset of PAR15 in 100mM NaCl is included.  
The analysis has already been performed and all outputs generated. The workflow for doing so is written below.
This analysis was performed using MATLAB R2021 and PyMol v1.5.5

An updated, more user-friendly class averaging workflow is available at: https://github.com/TongGeorgeWang/CASA-ToDiMo
/////////////

Copy all EOM PAR PDBs into folder
Use Power Rename function in Windows to give filenames systematic names:
	Shift-select all PDBs
	Hit F2
	The last file will be prompted to be renamed
	Rename the file to the basename (eg 'PAR15naPDB')
	All selected filenames will be renamed to '[basename] (i)' with i in ascending order
	The parentheses do not matter - when you import them into Pymol it will automatically remove them 

Load all systematically renamed files into Pymol 
	If ions are present, remove them in the command line:	
		remove resn SOD
		remove resn MGH
		remove resn CLA
Create alignAll pml script for pairwise RMSD matrix generation 
	Change basename variable and number of structures 
	Run (can go to File-->Run Script, sort by .pml or .txt)
		May need to change filetype to '.pml' if it's '.pml.txt'
	This should output the rmsd.txt file that contains all pairs of minimized RMSD values 

Input this into "PAR_classAverage.m" script, which is run to perform the class averaging (via spectral clustering). 
Save graph.fig and clusters.mat 

Run GroupClusters.m to automatically partition pdb files into cluster folders.

In each cluster folder, have an empty folder called Aligned and the script alignClusters.pml. 

Need to edit alignClusters.pml to have the right basename and also to list all
the PDB numbers in a comma separated array. To do this automatically: 
Copy all PDBs in each cluster folder to a throwaway directory for the purpose of renaming and 
automatically generating a comma-separated list of PDB numbers to put into alignCluster.pml script,
for clusters with many structures to avoid doing it by hand. 

Example for PAR15naPDB (running in Powershell): 

Truncate parts of filenames: 
get-childitem *.pdb | foreach { rename-item $_ $_.Name.Replace("PAR15naPDB (", "").Replace(")", "") }

Make comma-separated list of filenames:
(dir  | % { $_.basename }) -join ','


In Pymol, cd to a cluster folder, load all structures in a cluster, and run alignClusters.pml. Do this for all clusters.

You should now be ready to run sphagettiPlot.m, in the PDBs_SpectralClustered directory, to plot the final PAR backbone models.
