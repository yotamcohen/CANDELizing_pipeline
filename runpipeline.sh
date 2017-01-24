workingdir="/Users/yotam/CANDELizing_pipeline"
lcfile="$workingdir/example_lightcone.dat"
h5file="$workingdir/example_lightcone.h5"

python python_scripts/rawlightconefile_to_h5.py \
    lib/headerfile.txt $lcfile $h5file

python python_scripts/galaxy_size_corrections.py $h5file

python python_scripts/galaxy_profiles.py $h5file

python python_scripts/photometric_noise_new.py $h5file

python python_scripts/calculate_dp.py $h5file \
    ./lib/completeness_tables/goodss_expdisk.npz \
    ./lib/completeness_tables/goodss_devauc.npz

