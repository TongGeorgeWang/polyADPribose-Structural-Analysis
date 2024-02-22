set skip 1

set ID [mol new ../subset_strucs/nuc+ions.psf ; mol addfile ../subset_strucs/nuc+ions.pdb waitfor all]

mol addfile ../dcds/eq1_nuc+ions.dcd  step $skip waitfor all
mol addfile ../dcds/eq-LR1_nuc+ions.dcd step $skip waitfor all
mol addfile ../dcds/eq-LR2_nuc+ions.dcd step $skip waitfor all
mol addfile ../dcds/eqLR3_nuc+ions.dcd step $skip waitfor all
mol addfile ../dcds/eqLR4_nuc+ions.dcd step $skip waitfor all

set n [as top "nucleic"]

set nf [molinfo top get numframes]

set frames [molinfo $ID get numframes]

set cut_m 10

for {set r 1} {$r<16} {incr r 1} {
    set mon [as top "nucleic and resid $r"]
    set sel [as top "nucleic and (within $cut_m of [$mon text])"]
    set out_m [open data_$cut_m/monomer$r.dat a]
    for {set f 0} {$f<$nf} {incr f 1} {
	animate goto $f
	$sel update
	set num_m [expr [llength [lsort -unique [$sel get resid]]]-1]
	puts $out_m "$num_m"
    }
    close $out_m
    $sel delete
    $mon delete
    puts "Done for monomer $r"
}

exit
