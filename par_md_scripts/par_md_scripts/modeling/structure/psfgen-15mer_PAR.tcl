package require psfgen
resetpsf

topology topology/custom_PAR.str
topology topology/custom_BONDpatch.str
topology topology/custom_TERpatch.str
topology topology/top_all36_carb.rtf   # no modifications
topology topology/top_all36_na.rtf     # no modifications

set pre "pdbs_tmp/par"
set files ""
set suf ".pdb"
for {set i 0} {$i<15} {incr i 1} {
    set in [join "$pre$i$suf"]
    segment P$i {
        pdb $in
    }
    coordpdb $in
}

puts "Applying patches now"

for {set i 0} {$i<14} {incr i 1} {
    set j [expr $i+1]
    set jj [expr $j+1]
    patch BND P$i:$j P$j:$jj
}

patch 2TER P0:1
patch 1TER P14:15

puts "applied patches"

pdbalias residue HOH TIP3
pdbalias atom HOH O OH2
pdbalias residue AR6 PAR
pdbalias residue HIS HSE
pdbalias residue MSE MET
pdbalias residue CSD CYS

guesscoord 	 
regenerate angles dihedrals

writepsf par_15mer.psf 	 
writepdb par_15mer.pdb

exit
