#old_locs new_locs
biodiff.py $1 $2 -by ll | biointers.py $2 - -by cc | biodiff.py $2 - -by cc;
