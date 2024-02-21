output = open("rmsd.txt", "w")
basename = "PAR15naPDB_"
python 
Nstructures = 297
for i in range(1,Nstructures+1,1):
	print("alignment template: structure",str(i))
	for j in range(1,Nstructures+1,1):
		if j == i:
			print("special case: same structure")
			print(0)
			output.write("%s %s\n" % (str(i)+"0"+str(j), 0))
		else:
			print(basename+str(j)+")")
			data = cmd.align(basename+str(j), basename+str(i))
			print(data[0])
			output.write("%s %s\n" % (str(i)+"0"+str(j), data[0]))
python end 
output.close()