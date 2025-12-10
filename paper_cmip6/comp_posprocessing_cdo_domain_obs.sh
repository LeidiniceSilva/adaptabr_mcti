#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jun 10, 2024'
#__description__ = 'Posprocessing the observational database'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------" 

# Variables list
var_list=( 'tp' 'sp' 'msl' 'q' 't' 'u' 'v' )     

# Database
dataset='ERA5'

# Date
dt='197901-201412'

freq='mon'

for var in ${var_list[@]}; do

	path="/home/mda_silv/users/AdaptaBr_MCTI/database/paper_cmip6/obs"
	cd ${path}
	echo ${path}

	echo
	echo "1. Select subregion"
	
	cdo sellonlatbox,-70,-45,-5,5 ${var}_sa_${dataset}_${freq}_${dt}_lonlat.nc ${var}_namz_${dataset}_${freq}_${dt}_lonlat.nc
	cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_sa_${dataset}_${freq}_${dt}_lonlat.nc ${var}_samz_${dataset}_${freq}_${dt}_lonlat.nc
	cdo sellonlatbox,-55,-45,-20,-10 ${var}_sa_${dataset}_${freq}_${dt}_lonlat.nc ${var}_sam_${dataset}_${freq}_${dt}_lonlat.nc
	cdo sellonlatbox,-45,-34,-15,-2 ${var}_sa_${dataset}_${freq}_${dt}_lonlat.nc ${var}_neb_${dataset}_${freq}_${dt}_lonlat.nc
	cdo sellonlatbox,-60,-45,-35,-20 ${var}_sa_${dataset}_${freq}_${dt}_lonlat.nc ${var}_lpb_${dataset}_${freq}_${dt}_lonlat.nc

done

echo
echo "--------------- END POSPROCESSING OBS ----------------"















