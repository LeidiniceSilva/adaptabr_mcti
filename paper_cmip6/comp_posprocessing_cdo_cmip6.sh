#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 01, 2023'
#__description__ = 'Posprocessing the CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Models list
model_list=( 'MPI-ESM1-2-HR' ) 
#~ model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MIROC-ES2L' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' ) 

# Variables list
var_list=( 'hus' )     
#~ var_list=( 'hus' 'psl' 'ps' 'tas' 'ta' 'ua' 'va' 'zg' )     

for model in ${model_list[@]}; do

	for var in ${var_list[@]}; do

		path="/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/database/cmip6/"${model}
		cd ${path}

		echo
		echo ${path}

		# Experiment name
		exp='historical'

		# Member name
		if [ ${model} == 'CNRM-CM6-1' ]
		then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'CNRM-ESM2-1' ]
		then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'GFDL-ESM4' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'INM-CM4-8' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'INM-CM5-0' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'KIOST-ESM' ]
		then
		member='r1i1p1f1_gr1'
		elif [ ${model} == 'MIROC-ES2L' ]
		then
		member='r1i1p1f2_gn'
		else
		member='r1i1p1f1_gn'
		fi

		#~ # Datetime
		#~ if [ ${model} == 'ACCESS-CM2' ]
		#~ then
		#~ dt0='195001-201412'
		#~ elif [ ${model} == 'NorESM2-MM' ]
		#~ then
		#~ dt0='196001-201412'
		#~ else
		#~ dt0='185001-201412'
		#~ fi
		#~ dt1='197901-201412'

		dt0='197901-201412'
		dt1='1979-2014'

		echo ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc

		#~ echo
		#~ echo "1. Select date"
		#~ cdo seldate,1979-01-01T00:00:00,2014-12-31T00:00:00 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}.nc

		echo
		echo "2. Select levels"
		if [ ${var} == 'hus' ]
		then
		cdo sellevel,85000,70000,50000,20000 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		elif [ ${var} == 'ta' ]
		then
		cdo sellevel,85000,70000,50000,20000 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		elif [ ${var} == 'ua' ]
		then
		cdo sellevel,85000,70000,50000,20000 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		elif [ ${var} == 'va' ]
		then
		cdo sellevel,85000,70000,50000,20000 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		elif [ ${var} == 'zg' ]
		then
		cdo sellevel,85000,70000,50000,20000 ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		else
		mv ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		fi

		echo
		echo "3. Conventing grade"
		/home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc -60,15,0.25 -100,-20,0.25 bil

		echo
		echo "4. Conventing calendar"
		cdo setcalendar,standard ${var}_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc

		#~ echo 
		#~ echo "5. Conventing unit"
		#~ if [ ${var} == 'hus' ]
		#~ then
		#~ cdo mulc,1000 ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ elif [ ${var} == 'psl' ]
		#~ then
		#~ cdo divc,100 ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ elif [ ${var} == 'ps' ]
		#~ then
		#~ cdo divc,100 ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ elif [ ${var} == 'tas' ]
		#~ then
		#~ cdo subc,273.15 ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ elif [ ${var} == 'ta' ]
		#~ then
		#~ cdo subc,273.15 ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ elif [ ${var} == 'zg' ]
		#~ then
		#~ cdo mulc,10 ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ else
		#~ cp ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc ${var}_sa_Amon_${model}_${exp}_${member}_${dt1}_lonlat.nc
		#~ fi

		echo 
		echo "6. Deleting file"
		rm ${var}_Amon_${model}_${exp}_${member}_${dt0}.nc
		rm ${var}_Amon_${model}_${exp}_${member}_${dt1}_new.nc
		rm ${var}_Amon_${model}_${exp}_${member}_${dt1}_new_lonlat.nc
	
	done		
done

echo
echo "--------------- END POSPROCESSING CMIP6 MODELS ----------------"
