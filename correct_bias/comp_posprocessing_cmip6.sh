#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jul 05, 2023'
#__description__ = 'Posprocessing the CMIP6 models and obs database'

echo
echo "--------------- INIT POSPROCESSING ----------------"

# Model list
model_list=( 'GFDL-ESM4' 'INM-CM5-0' 'MPI-ESM1-2-HR' 'MRI-ESM2-0' 'NorESM2-MM' ) 

# Experiment list
# exp_list=( 'historical' 'ssp126' 'ssp245' 'ssp585' ) 
exp_list=( 'ssp585' ) 

# Variable list
var_list=( 'pr' 'tasmax' 'tasmin' )        

path_regrid="/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shell/regcm_pos"
path_mask="/scratch/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias/obs"

for model in ${model_list[@]}; do

	for exp in ${exp_list[@]}; do
	
		for var in ${var_list[@]}; do

			path="/scratch/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias/cmip6/${model}/${exp}"
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
			echo  ${var}_sa_day_${model}_${exp}_${member}_${dt}.nc

			echo
			echo "1. Conventing grade"
			${path_regrid}/./regrid ${var}_sa_day_${model}_${exp}_${member}_${dt}.nc -33.85,5.35,0.25 -73.85,-34.85,0.25 bil

			echo
			echo "2. Conventing calendar"
			cdo setcalendar,standard ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std.nc

			echo 
			echo "3. Conventing unit"
			if [ ${var} == 'pr' ]
			then
			cdo -b f32 mulc,86400 ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std.nc ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std_unit.nc
			else
			cdo -b f32 subc,273.15 ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std.nc ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std_unit.nc
			fi

			echo
			echo "4. Creating mask"
			cdo ifthen ${path_mask}/mask_br_lonlat.nc ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std_unit.nc ${var}_br_day_${model}_${exp}_${member}_${dt}_lonlat.nc

			#~ echo 
			#~ echo "5. Split years"
			#~ if [ ${exp} == 'historical' ]
			#~ then
			#~ for year in $(/usr/bin/seq -w 1986 2005); do
				#~ cdo seldate,${year}-01-01,${year}-12-31 ${var}_br_day_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_br_day_${model}_${exp}_${member}_${year}_lonlat.nc
			#~ done
			#~ else
			#~ for year in $(/usr/bin/seq -w 2015 2100); do
				#~ cdo seldate,${year}-01-01,${year}-12-31 ${var}_br_day_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_br_day_${model}_${exp}_${member}_${year}_lonlat.nc
			#~ done
			#~ fi
			
			echo 
			echo "6. Deleting file"
			rm ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat.nc 
			rm ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std.nc
			rm ${var}_sa_day_${model}_${exp}_${member}_${dt}_lonlat_std_unit.nc

		done
	done
done


echo
echo "--------------- END POSPROCESSING ----------------"
