#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jun 10, 2024'
#__description__ = 'Posprocessing the CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Models list
model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MIROC-ES2L' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' ) 

# Variables list
var_list=( 'hus' 'pr' 'ps' 'ta' 'ua' 'va' )                

freq='day' 

for model in ${model_list[@]}; do

	for var in ${var_list[@]}; do

		path="/home/mda_silv/users/AdaptaBr_MCTI/database/paper_cmip6/cmip6/"${model}
		cd ${path}

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
		echo ${path}
		echo ${var}_${freq}_${model}_${exp}_${member}_${dt}.nc

		echo
		echo "1. Select levels"
		if [ ${var} == 'hus' ] || [ ${var} == 'ta' ] || [ ${var} == 'ua' ] || [ ${var} == 'va' ]; then
		cdo sellevel,85000,50000,20000 ${var}_${freq}_${model}_${exp}_${member}_${dt}.nc ${var}_${freq}_${model}_${exp}_${member}_${dt}_new.nc
		
		else
		cp ${var}_${freq}_${model}_${exp}_${member}_${dt}.nc ${var}_${freq}_${model}_${exp}_${member}_${dt}_new.nc
		fi

		echo
		echo "2. Conventing grade"
		/home/mda_silv/github_projects/shell/ufrn/regcm_post/./regrid ${var}_${freq}_${model}_${exp}_${member}_${dt}_new.nc -60,15,1 -100,-20,1 bil

		echo
		echo "3. Conventing calendar"
		cdo setcalendar,standard ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat.nc ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc

		echo 
		echo "4. Conventing unit"
		if [ ${var} == 'hus' ]; then
		cdo -b f32 mulc,1000 ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc ${var}_sa_${freq}_${model}_${exp}_${member}_${dt}_lonlat.nc
		
		elif [ ${var} == 'pr' ]; then
		cdo -b f32 mulc,86400 ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc ${var}_sa_${freq}_${model}_${exp}_${member}_${dt}_lonlat.nc

		elif [ ${var} == 'ps' ]; then
		cdo -b f32 divc,100 ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc ${var}_sa_${freq}_${model}_${exp}_${member}_${dt}_lonlat.nc
		
		elif [ ${var} == 'ta' ]; then
		cdo -b f32 subc,273.15 ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc ${var}_sa_${freq}_${model}_${exp}_${member}_${dt}_lonlat.nc
		
		else
		cp ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc ${var}_sa_${freq}_${model}_${exp}_${member}_${dt}_lonlat.nc
		fi

		echo 
		echo "5. Deleting file"
		rm ${var}_${freq}_${model}_${exp}_${member}_${dt}_new.nc
		rm ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat.nc 
		rm ${var}_${freq}_${model}_${exp}_${member}_${dt}_new_lonlat_std.nc
	
	done
done

echo
echo "--------------- END POSPROCESSING CMIP6 MODELS ----------------"
