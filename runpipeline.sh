#raiddir="/Volumes/G-RAID/yotam/SAM-CANDELS_pipeline/catalogs/SC"
raiddir="/Volumes/G-RAID/yotam/catalogs/SC/12.16"

#python python_scripts/rawlightconefile_to_h5.py \
#    lib/headerfile.txt $raiddir/goodsn.dat $raiddir/goodsn.h5
#python python_scripts/galaxy_size_corrections.py $raiddir/goodsn.h5
#python python_scripts/galaxy_profiles.py $raiddir/goodsn.h5
#python python_scripts/photometric_noise_new.py $raiddir/goodsn.h5

python python_scripts/calculate_dp.py $raiddir/cosmos.h5 \
    ./lib/completeness_tables/cosmos_expdisk.npz \
    ./lib/completeness_tables/cosmos_devauc.npz
python python_scripts/detection_scheme.py $raiddir/cosmos.h5

python python_scripts/calculate_dp.py $raiddir/uds.h5 \
    ./lib/completeness_tables/uds_expdisk.npz \
    ./lib/completeness_tables/uds_devauc.npz
python python_scripts/detection_scheme.py $raiddir/uds.h5

python python_scripts/calculate_dp.py $raiddir/egs.h5 \
    ./lib/completeness_tables/egs_expdisk.npz \
    ./lib/completeness_tables/egs_devauc.npz
python python_scripts/detection_scheme.py $raiddir/egs.h5

python python_scripts/calculate_dp.py $raiddir/goodsn.h5 \
    ./lib/completeness_tables/goodsn_expdisk.npz \
    ./lib/completeness_tables/goodsn_devauc.npz
python python_scripts/detection_scheme.py $raiddir/goodsn.h5

python python_scripts/calculate_dp.py $raiddir/goodss.h5 \
    ./lib/completeness_tables/goodss_expdisk.npz \
    ./lib/completeness_tables/goodss_devauc.npz
python python_scripts/detection_scheme.py $raiddir/goodss.h5





