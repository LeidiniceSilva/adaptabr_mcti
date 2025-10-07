#!/bin/bash

#SBATCH -N 2
#SBATCH -t 24:00:00
#SBATCH -J CMIP6
#SBATCH -p esp
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=mda_silv@ictp.it

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'May 21, 2023'
#__description__ = 'Posprocessing the observational database'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------" 

# Variables list
var_list=( 'q' 't' 'v' )     
#var_list=( 'ps' 'q' 'tp' 't' 'u' 'v' )     

# Database
type='ERA5'

# Date
dt='197901-201412'

for var in ${var_list[@]}; do

	path="/home/mda_silv/users/AdaptaBr_MCTI/database/paper_cmip6/obs"
	cd ${path}

	echo
	echo ${path}
	echo ${var}_${type}_mon_${dt}.nc

	echo
	echo "1. Select levels"
	if [ ${var} == 'q' ] || [ ${var} == 't' ] || [ ${var} == 'u' ] || [ ${var} == 'v' ]; then
	cdo sellevel,850,500,200 ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc

	else
	cp ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	fi
		
	echo
	echo "2. Conventing grade"
	/home/mda_silv/github_projects/shell/ufrn/regcm_post/./regrid ${var}_${type}_mon_${dt}_new.nc -60,15,1 -100,-20,1 bil

	echo
	echo "3. Conventing calendar"
	cdo setcalendar,standard ${var}_${type}_mon_${dt}_new_lonlat.nc ${var}_${type}_mon_${dt}_new_lonlat_std.nc

	echo 
	echo "4. Conventing unit"
	if [ ${var} == 'ps' ]; then
	cdo -b f32 divc,100 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc

	elif [ ${var} == 'q' ]; then
	cdo -b f32 mulc,1000 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc

	elif [ ${var} == 'tp' ]; then
	cdo -b f32 mulc,1000 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
		
	elif [ ${var} == 't' ];	then
	cdo -b f32 subc,273.15 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc

	else
	cp ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	fi
		
	echo 
	echo "5. Deleting file"
	rm ${var}_${type}_mon_${dt}_new.nc
	rm ${var}_${type}_mon_${dt}_new_lonlat.nc 
	rm ${var}_${type}_mon_${dt}_new_lonlat_std.nc

done

echo
echo "--------------- END POSPROCESSING OBS ----------------"















