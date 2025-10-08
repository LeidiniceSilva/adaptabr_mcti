#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'May 21, 2023'
#__description__ = 'Posprocessing the CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Models list
model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MIROC-ES2L' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' ) 
# Variables list
var_list=( 'hus' 'pr' 'ps' 'ta' 'ua' 'va' )        

for model in ${model_list[@]}; do

	for var in ${var_list[@]}; do

		path="/home/mda_silv/users/AdaptaBr_MCTI/database/paper_cmip6/cmip6/"${model}
		cd ${path}
		echo ${path}
		
		# Experiment name
		exp='historical'

		# Member name
		if [ ${model} == 'CNRM-CM6-1' ]; then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'CNRM-ESM2-1' ]; then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'GFDL-ESM4' ]; then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'INM-CM4-8' ]; then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'INM-CM5-0' ]; then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'KIOST-ESM' ]; then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'MIROC-ES2L' ]; then
		member='r1i1p1f2_gn'
		else
		member='r1i1p1f1_gn'
		fi

		dt='197901-201412'

		echo
		echo "1. Select subregion"
		cdo sellonlatbox,-70,-45,-5,5 ${var}_sa_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_namz_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_sa_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_samz_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_sa_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_sam_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc
		cdo sellonlatbox,-45,-34,-15,-2 ${var}_sa_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_neb_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc
		cdo sellonlatbox,-60,-45,-35,-20 ${var}_sa_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_lpb_Amon_${model}_${exp}_${member}_${dt}_lonlat.nc

	done
done

echo
echo "--------------- END POSPROCESSING CMIP6 MODELS ----------------"

