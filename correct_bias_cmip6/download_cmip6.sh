#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jul 05, 2023'
#__description__ = 'Download CMIP6 simulations and projections'

# center
center=( 'NCC' )

# model 
model=( 'NorESM2-MM' )

# Variable list
var_list=( 'pr' 'tasmax' 'tasmin' )

for var in ${var_list[@]}; do

	# Path to save file
	PATH="/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/bias_adjust/cmip6/${model}/historical"    	    
	cd ${PATH}

	# Download file
	/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_historical_r1i1p1f1_gn_19800101-19891231.nc			
	/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_historical_r1i1p1f1_gn_19900101-19991231.nc		
	/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_historical_r1i1p1f1_gn_20000101-20091231.nc
	
done

# scenario list
scenario_list=( 'ssp126' 'ssp245' 'ssp585' )

for exp in ${scenario_list[@]}; do
	for var in ${var_list[@]}; do

		# Path to save file
		PATH="/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/bias_adjust/cmip6/${model}/${exp}"    	    
		cd ${PATH}

		# Download file
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20150101-20201231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20210101-20301231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20310101-20401231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20410101-20501231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20510101-20601231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20610101-20701231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20710101-20801231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20810101-20901231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20191108/${var}_day_${model}_${exp}_r1i1p1f1_gn_20910101-21001231.nc
		  
	done
done
