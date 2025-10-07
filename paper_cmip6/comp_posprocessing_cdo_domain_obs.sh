#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jun 10, 2024'
#__description__ = 'Posprocessing the observational database'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------" 

# Variables list
var_list=( 'msl' 'mtpr' 'sp' 'q' 't2m' 't' 'u' 'v' 'z' )     

# Database
dataset='era5'

# Date
dt='197901-201412'

for var in ${var_list[@]}; do

	path="/afs/ictp.it/home/m/mda_silv/Documents/CMIP6/database/obs/"
	cd ${path}
	echo ${path}

	echo
	echo "1. Select subregion"
	
	cdo sellonlatbox,-70,-45,-5,5 ${var}_sa_${dataset}_mon_${dt}_lonlat.nc ${var}_namz_${dataset}_mon_${dt}_lonlat.nc
	cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_sa_${dataset}_mon_${dt}_lonlat.nc ${var}_samz_${dataset}_mon_${dt}_lonlat.nc
	cdo sellonlatbox,-55,-45,-20,-10 ${var}_sa_${dataset}_mon_${dt}_lonlat.nc ${var}_sam_${dataset}_mon_${dt}_lonlat.nc
	cdo sellonlatbox,-45,-34,-15,-2 ${var}_sa_${dataset}_mon_${dt}_lonlat.nc ${var}_neb_${dataset}_mon_${dt}_lonlat.nc
	cdo sellonlatbox,-60,-45,-35,-20 ${var}_sa_${dataset}_mon_${dt}_lonlat.nc ${var}_lpb_${dataset}_mon_${dt}_lonlat.nc

done

echo
echo "--------------- END POSPROCESSING OBS ----------------"















