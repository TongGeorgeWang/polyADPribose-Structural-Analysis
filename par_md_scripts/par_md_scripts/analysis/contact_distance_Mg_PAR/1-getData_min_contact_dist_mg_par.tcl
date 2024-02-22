set skip 1
set ID [mol new ../subset_strucs/nuc+ions.psf ; mol addfile ../subset_strucs/nuc+ions.pdb waitfor all]

mol addfile unwrapped_eq1-LR4_nuc+mg_skip1.dcd waitfor all

set nf [molinfo top get numframes]

set mg [as top "name MG"]
set par [as top "nucleic"]

set frames_min_dist ""

for {set f 0} {$f<$nf} {incr f 1} {
    animate goto $f
    $mg update
    $par update

    # Choose a high value for min_dist, ideally infinity
    set min_dist 10000

    # Looping over all atom pairs based on selections
    foreach atom1 [$mg get index] {
        foreach atom2 [$par get index] {
            set list_atomPair "$atom1 $atom2"
            set dist [measure bond $list_atomPair]
            if {$dist < $min_dist} {
                set min_dist $dist
            }
        }
    }

    # Record the minimum distance for this frame
    lappend frames_min_dist $min_dist
    puts "Done frame $f"
}

set out [open data_min_dists.dat w]
foreach d $frames_min_dist {
    puts $out $d
}
close $out

exit
