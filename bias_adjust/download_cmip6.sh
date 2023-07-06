#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = '07/05/2023'
#__description__ = 'Download CMIP6 simulations and projections'

exp='ssp245'

PATH="/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/bias_adjust/cmip6/GFDL-ESM4/"${exp}    	    
cd ${PATH}


/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/pr/gr1/files/d20180701/pr_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20150101-20341231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/pr/gr1/files/d20180701/pr_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20350101-20541231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/pr/gr1/files/d20180701/pr_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20550101-20741231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/pr/gr1/files/d20180701/pr_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20750101-20941231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/pr/gr1/files/d20180701/pr_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20950101-21001231.nc

/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmax/gr1/files/d20180701/tasmax_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20150101-20341231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmax/gr1/files/d20180701/tasmax_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20350101-20541231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmax/gr1/files/d20180701/tasmax_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20550101-20741231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmax/gr1/files/d20180701/tasmax_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20750101-20941231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmax/gr1/files/d20180701/tasmax_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20950101-21001231.nc

/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmin/gr1/files/d20180701/tasmin_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20150101-20341231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmin/gr1/files/d20180701/tasmin_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20350101-20541231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmin/gr1/files/d20180701/tasmin_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20550101-20741231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmin/gr1/files/d20180701/tasmin_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20750101-20941231.nc	
/usr/bin/wget wget -N -c https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/NOAA-GFDL/GFDL-ESM4/${exp}/r1i1p1f1/day/tasmin/gr1/files/d20180701/tasmin_day_GFDL-ESM4_${exp}_r1i1p1f1_gr1_20950101-21001231.nc
