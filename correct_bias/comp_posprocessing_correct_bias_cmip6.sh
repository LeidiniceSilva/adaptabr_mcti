#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jul 05, 2023'
#__description__ = 'Posprocessing correction bias from CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING ----------------"

# Model list
model_list=( 'GFDL-ESM4' 'INM-CM5-0' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'NorESM2-MM' ) 

# Experiment list
exp_list=( 'historical' 'ssp126' 'ssp245' 'ssp585' ) 

# Variable list
var_list=( 'pr' 'tasmax' 'tasmin' )     

path_regrid="/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shell/regcm_pos"
path_mask="/scratch/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias/obs"

for model in ${model_list[@]}; do

	for exp in ${exp_list[@]}; do
	
		for var in ${var_list[@]}; do

			path="/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias/cmip6_correct/${model}/${exp}"
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
			echo "1. Creating mask"			
			if [ ${var} == 'pr' ]
			then
			cdo ifthen ${path_mask}/mask_br_lonlat.nc ${var}_br_day_${model}_${exp}_${member}_${dt}_correct.nc ${var}_br_day_${model}_${exp}_${member}_${dt}_corrected.nc
			else
			cp ${var}_br_day_${model}_${exp}_${member}_${dt}_correct.nc ${var}_br_day_${model}_${exp}_${member}_${dt}_corrected.nc
			fi
						
			echo 
			echo "2. Deleting file"
			# rm ${var}_br_day_${model}_${exp}_${member}_${dt}_correct.nc

		done
	done
done
			
echo
echo "--------------- END POSPROCESSING ----------------"
