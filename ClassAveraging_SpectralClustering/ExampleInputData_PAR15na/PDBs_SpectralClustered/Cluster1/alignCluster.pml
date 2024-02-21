basename = "PAR15naPDB_"


remove resn SOD
remove resn CLA
remove resn MGH

python 

import numpy as pl
PDBnumbers = pl.array([136,162,188,19,191,192,2,216,217,218,220,228,23,238,239,24,245,246,247,248,249,25,250,26,260,263,264,267,268,269,27,270,271,272,277,278,279,28,280,281,282,283,284,285,29,293,296,297,31,32])

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
