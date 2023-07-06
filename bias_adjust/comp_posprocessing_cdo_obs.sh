#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'May 21, 2023'
#__description__ = 'Posprocessing the observational database'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------" 

# Variables list
var_list=( 'mslp' 'mtpr' 'ps' 'q' 't2m' 't' 'u' 'v' 'z' )     

# Database
type='era5'

# Date
dt='197901-201412'

for var in ${var_list[@]}; do

	path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/paper_cmip6/obs/"
	cd ${path}

	echo
	echo ${path}
	echo ${var}_${type}_mon_${dt}.nc

	echo
	echo "1. Select levels"
	if [ ${var} == 'q' ]
	then
	cdo sellevel,850,500,200 ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	
	elif [ ${var} == 't' ]
	then
	cdo sellevel,850,500,200 ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	
	elif [ ${var} == 'u' ]
	then
	cdo sellevel,850,500,200 ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	
	elif [ ${var} == 'v' ]
	then
	cdo sellevel,850,500,200 ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	
	elif [ ${var} == 'z' ]
	then
	cdo sellevel,850,500,200 ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	
	else
	cp ${var}_${type}_mon_${dt}.nc ${var}_${type}_mon_${dt}_new.nc
	fi
		
	echo
	echo "2. Conventing grade"
	/home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_${type}_mon_${dt}_new.nc -60,15,0.25 -100,-20,0.25 bil

	echo
	echo "3. Conventing calendar"
	cdo setcalendar,standard ${var}_${type}_mon_${dt}_new_lonlat.nc ${var}_${type}_mon_${dt}_new_lonlat_std.nc

	echo 
	echo "4. Conventing unit"
	if [ ${var} == 'mslp' ]
	then
	cdo -b f32 divc,100 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	
	elif [ ${var} == 'mtpr' ]
	then
	cdo -b f32 mulc,86400 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
		
	elif [ ${var} == 'ps' ]
	then
	cdo -b f32 divc,100 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	
	elif [ ${var} == 'q' ]
	then
	cdo -b f32 mulc,1000 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	
	elif [ ${var} == 't2m' ]
	then
	cdo -b f32 subc,273.15 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	
	elif [ ${var} == 't' ]
	then
	cdo -b f32 subc,273.15 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	
	elif [ ${var} == 'z' ]
	then
	cdo -b f32 divc,10 ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	
	else
	cp ${var}_${type}_mon_${dt}_new_lonlat_std.nc ${var}_sa_${type}_mon_${dt}_lonlat.nc
	fi
		
	echo 
	echo "5. Deleting file"
	rm ${var}_${type}_mon_${dt}.nc
	rm ${var}_${type}_mon_${dt}_new.nc
	rm ${var}_${type}_mon_${dt}_new_lonlat.nc 
	rm ${var}_${type}_mon_${dt}_new_lonlat_std.nc

done

echo
echo "--------------- END POSPROCESSING OBS ----------------"















