
set -xe 

project_potential.jl ../potential_lmc_evolving.ini -T times.txt -k 1 --limits 500 -n 1001 -o projected_lmc.hdf5
