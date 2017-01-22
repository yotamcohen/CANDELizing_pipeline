raiddir="/Volumes/G-RAID/yotam/catalogs/SC/12.16"

#python python_scripts/galaxy_size_corrections.py $raiddir/cosmos.h5
#python python_scripts/galaxy_size_corrections.py $raiddir/uds.h5
#python python_scripts/galaxy_size_corrections.py $raiddir/egs.h5
#python python_scripts/galaxy_size_corrections.py $raiddir/goodsn.h5
#python python_scripts/galaxy_size_corrections.py $raiddir/goodss.h5

python python_scripts/write_final_datfiles.py $raiddir/cosmos.h5 $raiddir/CANDELized_cosmos.dat
python python_scripts/write_final_datfiles.py $raiddir/uds.h5 $raiddir/CANDELized_uds.dat
python python_scripts/write_final_datfiles.py $raiddir/egs.h5 $raiddir/CANDELized_egs.dat
python python_scripts/write_final_datfiles.py $raiddir/goodsn.h5 $raiddir/CANDELized_goodsn.dat
python python_scripts/write_final_datfiles.py $raiddir/goodss.h5 $raiddir/CANDELized_goodss.dat

#python python_scripts/photoz_est.py $raiddir/cosmos.h5
#python python_scripts/photoz_est.py $raiddir/uds.h5
#python python_scripts/photoz_est.py $raiddir/egs.h5
#python python_scripts/photoz_est.py $raiddir/goodsn.h5
#python python_scripts/photoz_est.py $raiddir/goodss.h5


