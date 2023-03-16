#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = '03/13/2023'
#__description__ = 'Posprocessing the CMIP6 models dataset'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Variables list
var_list=( 'pr' 'tas' 'tasmax' 'tasmin' )     

# Models list
model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-CM6-1-HR' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' )

for var in ${var_list[@]}; do

    for model in ${model_list[@]}; do

		path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6"
		cd ${path}
		
		echo
		echo ${path}
		
		# Experiment name
		exp = 'historical'

		# Member name
		if [ ${model} == 'CNRM-CM6-1' 'CNRM-CM6-1-HR' 'CNRM-ESM2-1' ]
		then
		member = 'r1i1p1f2_gr'
		else
		member = 'r1i1p1f1_gr1'
		fi
		
		# Datetime
		dt = '18500116-20141216'
		
		echo "1. Select date"
		cdo seldate,1981-01-01T00:00:00,2014-12-31T00:00:00 ${var}_Amon_${model}_${exp}_${member}_${dt}_*.nc ${var}_Amon_${model}_${exp}_${member}_19810101-20141231.nc

		echo
		echo "2. Convert calendar"
		cdo setcalendar,standard ${var}_Amon_${model}_${exp}_${member}_19810101-20141231.nc ${var}_Amon_${model}_${exp}_${member}_1981-2014.nc
		
		echo
		echo "3. Remapbil"
		/home/nice/Documents/github_projects/shell/regcm_pos/./regrid ${var}_Amon_${model}_${exp}_${member}_1981-2014.nc -60,15,0.5 -90,-30,0.5 bil

		echo 
		echo "4. Convention unit"
		cdo mulc,86400 ${var}_Amon_${model}_${exp}_${member}_1981-2014_lonlat.nc ${var}_SA_Amon_${model}_${exp}_${member}_1981-2014.nc
		
		echo
		echo "5. Select subdomain"
		cdo sellonlatbox,-68,-52,-12,-3  ${var}_SA_Amon_${model}_${exp}_${member}_1981-2014.nc ${var}_SAMZ_Amon_${model}_${exp}_${member}_1981-2014.nc
		cdo sellonlatbox,-61,-48,-37,-27 ${var}_SA_Amon_${model}_${exp}_${member}_1981-2014.nc ${var}_NAMZ_Amon_${model}_${exp}_${member}_1981-2014.nc
		cdo sellonlatbox,-47,-35,-15,-2  ${var}_SA_Amon_${model}_${exp}_${member}_1981-2014.nc ${var}_SAM_Amon_${model}_${exp}_${member}_1981-2014.nc
		cdo sellonlatbox,-47,-35,-15,-2  ${var}_SA_Amon_${model}_${exp}_${member}_1981-2014.nc ${var}_NEB_Amon_${model}_${exp}_${member}_1981-2014.nc
		cdo sellonlatbox,-47,-35,-15,-2  ${var}_SA_Amon_${model}_${exp}_${member}_1981-2014.nc ${var}_LPB_Amon_${model}_${exp}_${member}_1981-2014.nc

		echo 
		echo "6. Deleting file"
		rm ${var}_Amon_${model}_${exp}_${member}_19810101-20141231.nc
		rm ${var}_Amon_${model}_${exp}_${member}_1981-2014.nc

		rm ${var}_CRU_ts4_mon_1981-2014.nc
		rm ${var}_CRU_ts4_mon_1981-2014_lonlat.nc
	
	done
done

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"
