#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Jul 05, 2023'
#__description__ = 'Download CMIP6 simulations and projections'

# center
center=( 'DKRZ' )

# model 
model=( 'MPI-ESM1-2-HR' )

# Variable list
var_list=( 'pr' 'tasmax' 'tasmin' )

#~ for var in ${var_list[@]}; do

	#~ # Path to save file
	#~ PATH="/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/bias_adjust/cmip6/${model}/historical"    	    
	#~ cd ${PATH}

	#~ # Download file
	#~ /usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_historical_r1i1p1f1_gn_19850101-19891231.nc			
	#~ /usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_historical_r1i1p1f1_gn_19900101-19941231.nc		
	#~ /usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_historical_r1i1p1f1_gn_19950101-19991231.nc		
	#~ /usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_historical_r1i1p1f1_gn_20000101-20041231.nc		
	#~ /usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/CMIP/${center}/${model}/historical/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_historical_r1i1p1f1_gn_20050101-20091231.nc

#~ done

# scenario list
scenario_list=( 'ssp126' 'ssp245' 'ssp585' )

for exp in ${scenario_list[@]}; do
	for var in ${var_list[@]}; do

		# Path to save file
		PATH="/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/bias_adjust/cmip6/${model}/${exp}"    	    
		cd ${PATH}

		# Download file
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20150101-20191231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20200101-20241231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20250101-20291231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20300101-20341231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20350101-20391231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20400101-20441231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20450101-20491231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20500101-20541231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20550101-20591231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20600101-20641231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20650101-20691231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20700101-20741231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20750101-20791231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20800101-20841231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20850101-20891231.nc
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20900101-20941231.nc			
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_20950101-20991231.nc		
		/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/${center}/${model}/${exp}/r1i1p1f1/day/${var}/gn/files/d20190710/${var}_day_${model}_${exp}_r1i1p1f1_gn_21000101-21001231.nc
		  
	done
done
