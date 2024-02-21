basename = "PAR15naPDB_"


remove resn SOD
remove resn CLA
remove resn MGH

python 

import numpy as pl
PDBnumbers = pl.array([164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,189,193,194,195,196,197,198,199,20,200,201,202,203,204,205,206,207,208,209,21,210,211,212,213,214,215,219,22,221,222,223,224,225,226,227,229,230,231,232,233,234,235,236,240,241,242,243,244,251,252,253,254,255,256,257,258,259,273,274,275,276,286,287,288,290,291,292,294,295])

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
