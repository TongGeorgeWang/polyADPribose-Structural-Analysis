basename = "PAR15naPDB_"


remove resn SOD
remove resn CLA
remove resn MGH

python 

import numpy as pl
PDBnumbers = pl.array([100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,126,127,128,129,130,163,265,266,289,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99])

Nstructures = len(PDBnumbers)
RMSDmatrix = pl.zeros((Nstructures,Nstructures))


## Compute all pairwise alignments
for i in range(0,Nstructures,1):
	print(PDBnumbers[i])
	print(i)
	for j in range(0,Nstructures,1):
		if j == i:
			print("special case: same structure. Not aligning.")
		else:
			print(PDBnumbers[j])
			data = cmd.align(basename+str(PDBnumbers[i]), basename+str(PDBnumbers[j]))
			print(data[0])
			RMSDmatrix[i,j] = data[0]
	print()


## Determine minimum RMSD values 
whereToAlign = pl.zeros((Nstructures,2))
whereToAlign[:,0] = pl.arange(0,Nstructures,1)
for k in range(0,Nstructures,1):
	currentRow = RMSDmatrix[k,:]
	currentRow[currentRow==0] = pl.nan
	print(currentRow)
	minValue = pl.nanmin(currentRow)
	print(minValue)
	minIndices = pl.where(currentRow==minValue)
	print(minIndices[0])
	whereToAlign[k,1] = minIndices[0]

print(whereToAlign)



## Align structures that have minimum resulting RMSDs
for i in range(0,Nstructures,1):
	structure1 = int(whereToAlign[i,0])
	structure2 = int(whereToAlign[i,1])
	print(structure1)
	print(structure2)
	data = cmd.align( basename+str(PDBnumbers[structure1]), basename+str(PDBnumbers[structure2]) )
	print(data[0])
	print()
	


## Save aligned PDBs
for i in range(0,Nstructures,1):
	cmd.save( 'Aligned/'+basename+str(PDBnumbers[i])+'.pdb' , basename+str(PDBnumbers[i]))


python end
