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
dt='1961-2014'

for var in ${var_list[@]}; do

	path="/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/database/obs"
	cd ${path}

	echo
	echo ${path}

	echo ${var}_${type}_mon_${dt}.nc
	echo
		
	echo
	echo "1. Conventing calendar"
	cdo setcalendar,standard ${var}_${type}_mon_${dt}.nc ${var}_sa_${type}_mon_${dt}.nc

	echo
	echo "2. Conventing grade"
	/home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_sa_${type}_mon_${dt}.nc -60,15,0.25 -100,-20,0.25 bil

	echo 
	echo "3. Deleting file"
	rm ${var}_${type}_mon_${dt}.nc
	rm ${var}_sa_${type}_mon_${dt}.nc

done

echo
echo "--------------- END POSPROCESSING OBS ----------------"














