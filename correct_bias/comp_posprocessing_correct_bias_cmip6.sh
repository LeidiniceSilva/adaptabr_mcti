#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jul 05, 2023'
#__description__ = 'Posprocessing correction bias from CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING ----------------"

# Model list
# model_list=( 'GFDL-ESM4' 'INM-CM5-0' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'NorESM2-MM' ) 
model_list=( 'NorESM2-MM' ) 

# Experiment list
# exp_list=( 'historical' 'ssp126' 'ssp245' 'ssp585' ) 
exp_list=( 'historical' ) 

# Variable list
var_list=( 'pr' 'tasmax' 'tasmin' )     

path_regrid="/home/nice/Documentos/github_projects/shell/regcm_pos"
path_mask="/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias/obs"

for model in ${model_list[@]}; do

	for exp in ${exp_list[@]}; do
	
		for var in ${var_list[@]}; do

			path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias/cmip6_correct/${model}/${exp}"
			cd ${path}

			# Member name
			if [ ${model} == 'GFDL-ESM4' ]
			then
			member='r1i1p1f1_gr1'
			elif [ ${model} == 'INM-CM5-0' ]
			then
			member='r1i1p1f1_gr1'
			else
			member='r1i1p1f1_gn'
			fi

			# Datetime
			if [ ${exp} == 'historical' ]
			then
			dt='19860101-20051231'
			else
			dt='20150101-21001231'
			fi
			
			echo
			echo ${path}
			echo ${var}_br_day_${model}_${exp}_${member}_${dt}_correct.nc

			echo
			echo "1. Merge date" 
			cdo mergetime ${var}_br_day_${model}_${exp}_${member}_*_correct.nc ${var}_br_day_${model}_${exp}_${member}_1986-2005_correct.nc

			echo
			echo "2. Creating mask"			
			if [ ${var} == 'pr' ]
			then
			cdo ifthen ${path_mask}/mask_br_lonlat.nc ${var}_br_day_${model}_${exp}_${member}_1986-2005_correct.nc ${var}_br_day_${model}_${exp}_${member}_${dt}_correct.nc
			else
			cp ${var}_br_day_${model}_${exp}_${member}_1986-2005_correct.nc ${var}_br_day_${model}_${exp}_${member}_${dt}_correct.nc
			fi
						
			echo 
			echo "3. Deleting file"
			rm ${var}_br_day_${model}_${exp}_${member}_1986-2005_correct.nc

		done
	done
done
			
echo
echo "--------------- END POSPROCESSING ----------------"
