#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 01, 2023'
#__description__ = 'Posprocessing the CMIP6 models'

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"

# Variables list
var_list=( 'pr' 'tas' )     

# Models list
model_list=( 'ACCESS-CM2' 'BCC-CSM2-MR' 'CanESM5' 'CMCC-ESM2' 'CNRM-CM6-1' 'CNRM-CM6-1-HR' 'CNRM-ESM2-1' 'GFDL-ESM4' 'INM-CM4-8' 'INM-CM5-0' 'KIOST-ESM' 'MIROC6' 'MME-CMIP6' 'MPI-ESM1-2-HR' 'MPI-ESM1-2-LR' 'MRI-ESM2-0' 'NESM3' 'NorESM2-MM' )

for var in ${var_list[@]}; do

    for model in ${model_list[@]}; do

		path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6/"
		cd ${path}
		
		echo
		echo ${path}
		
		# Experiment name
		exp='historical'

		# Member name
		if [ ${model} == 'CNRM-CM6-1' ]
		then
		member='r1i1p1f2_gr'
		elif [ ${model} == 'CNRM-CM6-1-HR' ]
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
		else
		member='r1i1p1f1_gn'
		fi
		
		# Datetime
		dt='1980-2014'
		
		echo
		echo ${var}_Amon_${model}_${exp}_${member}_185001-201412.nc
		
		#~ echo "1. Select date"
		#~ cdo seldate,1980-01-01T00:00:00,2014-12-31T00:00:00 ${var}_Amon_${model}_${exp}_${member}_185001-201412.nc ${var}_Amon_${model}_${exp}_${member}_${dt}.nc

		#~ echo
		#~ echo "2. Conventing calendar"
		#~ cdo setcalendar,standard ${var}_Amon_${model}_${exp}_${member}_${dt}.nc ${var}_${model}_${exp}_${member}_${dt}.nc
		
		#~ echo
		#~ echo "3. Conventing grade"
		#~ /home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_${model}_${exp}_${member}_${dt}.nc -60,15,0.5 -85,-30,0.5 bil

		#~ echo 
		#~ echo "4. Conventing unit"
		#~ if [ ${var} == 'pr' ]
		#~ then
		#~ cdo mulc,86400 ${var}_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_${dt}_lonlat.nc
		#~ else
		#~ cdo subc,273.15 ${var}_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_${dt}_lonlat.nc
		#~ fi

		#~ echo
		#~ echo "5. Creating mask"
		#~ cdo -f nc -remapnn,${var}_SA_${model}_${exp}_${member}_${dt}_lonlat.nc -gtc,0 -topo ${var}_${model}_seamask.nc
		#~ cdo ifthen ${var}_${model}_seamask.nc ${var}_SA_${model}_${exp}_${member}_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc
		
		echo
		echo "6. Calculate periods"
		cdo -yearavg ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc
		cdo -r -timselavg,3 -selmon,1,2,12 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc
		cdo -r -timselavg,3 -selmon,3,4,5 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc
		cdo -r -timselavg,3 -selmon,6,7,8 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc
		cdo -r -timselavg,3 -selmon,9,10,11 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SA_${model}_${exp}_${member}_SON_${dt}_lonlat.nc

		echo
		echo "7. Select subregion"
		cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_NAMZ_${model}_${exp}_${member}_MON_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc ${var}_NAMZ_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc ${var}_NAMZ_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc ${var}_NAMZ_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc ${var}_NAMZ_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${model}_${exp}_${member}_SON_${dt}_lonlat.nc ${var}_NAMZ_${model}_${exp}_${member}_SON_${dt}_lonlat.nc

		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SAMZ_${model}_${exp}_${member}_MON_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc ${var}_SAMZ_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc ${var}_SAMZ_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc ${var}_SAMZ_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc ${var}_SAMZ_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc
		cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${model}_${exp}_${member}_SON_${dt}_lonlat.nc ${var}_SAMZ_${model}_${exp}_${member}_SON_${dt}_lonlat.nc

		cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_NEB_${model}_${exp}_${member}_MON_${dt}_lonlat.nc
		cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc ${var}_NEB_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc
		cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc ${var}_NEB_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc
		cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc ${var}_NEB_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc
		cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc ${var}_NEB_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc
		cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${model}_${exp}_${member}_SON_${dt}_lonlat.nc ${var}_NEB_${model}_${exp}_${member}_SON_${dt}_lonlat.nc
				
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_SAM_${model}_${exp}_${member}_MON_${dt}_lonlat.nc
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc ${var}_SAM_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc ${var}_SAM_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc ${var}_SAM_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc ${var}_SAM_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc
		cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${model}_${exp}_${member}_SON_${dt}_lonlat.nc ${var}_SAM_${model}_${exp}_${member}_SON_${dt}_lonlat.nc

		cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${model}_${exp}_${member}_MON_${dt}_lonlat.nc ${var}_LPB_${model}_${exp}_${member}_MON_${dt}_lonlat.nc
		cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc ${var}_LPB_${model}_${exp}_${member}_ANN_${dt}_lonlat.nc
		cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc ${var}_LPB_${model}_${exp}_${member}_DJF_${dt}_lonlat.nc
		cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc ${var}_LPB_${model}_${exp}_${member}_MAM_${dt}_lonlat.nc
		cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc ${var}_LPB_${model}_${exp}_${member}_JJA_${dt}_lonlat.nc
		cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${model}_${exp}_${member}_SON_${dt}_lonlat.nc ${var}_LPB_${model}_${exp}_${member}_SON_${dt}_lonlat.nc

		#~ echo 
		#~ echo "8. Deleting file"
		#~ rm ${var}_Amon_${model}_${exp}_${member}_${dt}.nc
		#~ rm ${var}_${model}_${exp}_${member}_${dt}.nc
		#~ rm ${var}_${model}_${exp}_${member}_${dt}_lonlat.nc
		#~ rm ${var}_SA_${model}_${exp}_${member}_${dt}_lonlat.nc
		#~ rm ${var}_${model}_seamask.nc
	
	done
done

echo
echo "--------------- INIT POSPROCESSING CMIP6 MODELS ----------------"
