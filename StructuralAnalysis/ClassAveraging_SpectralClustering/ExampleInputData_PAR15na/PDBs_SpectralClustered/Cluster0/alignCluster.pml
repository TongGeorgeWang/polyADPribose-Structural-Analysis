basename = "PAR15naPDB_"


remove resn SOD
remove resn CLA
remove resn MGH

python 

import numpy as pl
PDBnumbers = pl.array([1,10,11,12,125,13,131,132,133,134,135,137,138,139,14,140,141,142,143,144,145,146,147,148,149,15,150,151,152,153,154,155,156,157,158,159,16,160,161,17,18,190,237,261,262,3,30,4,5,6,7,8,9])

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
