
pot_dir=../../potentials/vasiliev24/L2M11/
cp $pot_dir/potential.ini agama_potential.ini
cp $pot_dir/potential_lmc_init.ini .
cp $pot_dir/potential_mw_init.ini .
cp $pot_dir/trajlmc.txt .
cp $pot_dir/boundmass.txt .
cp $pot_dir/accel.txt .
ln -s $pot_dir/pot ./

