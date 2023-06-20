#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'May 21, 2023'
#__description__ = 'Posprocessing the observational database'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------" 

# Variables list
var_list=( 'mslp' 'ps' 'q' 't2m' 't' 'u' 'v' 'z' )     

# Database
type='era5'

# Date
dt='197901-201412'

for var in ${var_list[@]}; do

	path="/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/database/obs/"
	cd ${path}

	echo
	echo ${path}

	echo ${var}_${type}_mon_${dt}.nc
	echo
		
	echo
	echo "1. Conventing grade"
	/home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_${type}_mon_${dt}.nc -60,15,0.25 -100,-20,0.25 bil

	echo
	echo "2. Conventing calendar"
	cdo setcalendar,standard ${var}_${type}_mon_${dt}_lonlat.nc ${var}_${type}_mon_${dt}_lonlat_std.nc

	echo 
	echo "3. Conventing unit"
	if [ ${var} == 'hus' ]
	then
	cdo mulc,1000 ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	elif [ ${var} == 'psl' ]
	then
	cdo divc,100 ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	elif [ ${var} == 'ps' ]
	then
	cdo divc,100 ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	elif [ ${var} == 'tas' ]
	then
	cdo subc,273.15 ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	elif [ ${var} == 'ta' ]
	then
	cdo subc,273.15 ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	elif [ ${var} == 'zg' ]
	then
	cdo mulc,10 ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	else
	mv ${var}_${type}_mon_${dt}_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	fi
		
	echo 
	echo "4. Deleting file"
	rm ${var}_${type}_mon_${dt}.nc
	rm ${var}_${type}_mon_${dt}_lonlat.nc 
	rm ${var}_${type}_mon_${dt}_lonlat_std.nc

done


echo
echo "--------------- END POSPROCESSING OBS ----------------"















