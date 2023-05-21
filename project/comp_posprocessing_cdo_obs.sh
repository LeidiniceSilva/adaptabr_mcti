#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 01, 2023'
#__description__ = 'Posprocessing the observational database'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------"

# Variable 
var='Tmin'  

# Database
type='BR-DWGD_UFES_UTEXAS_v_3.0'

# Date
dt='1986-2005'

path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs/"
cd ${path}

echo
echo ${path}

echo 
echo "1. Select date"
cdo seldate,1986-01-01T00:00:00,2005-12-31T00:00:00 ${var}_${type}_MON_198101-202007.nc ${var}_${type}_MON_${dt}.nc

echo
echo "2. Conventing calendar"
cdo setcalendar,standard ${var}_${type}_MON_${dt}.nc ${var}_SA_${type}_MON_${dt}.nc

echo
echo "3. Conventing grade"
/home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_SA_${type}_MON_${dt}.nc -60,15,0.25 -85,-30,0.25 bil

echo
echo "4. Calculate periods"
cdo -yearavg ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_SA_${type}_ANN_${dt}_lonlat.nc

echo
echo "5. Select subregion"
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_NAMZ_${type}_MON_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_${type}_ANN_${dt}_lonlat.nc ${var}_NAMZ_${type}_ANN_${dt}_lonlat.nc

cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_SAMZ_${type}_MON_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_${type}_ANN_${dt}_lonlat.nc ${var}_SAMZ_${type}_ANN_${dt}_lonlat.nc

cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_NEB_${type}_MON_${dt}_lonlat.nc
cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_${type}_ANN_${dt}_lonlat.nc ${var}_NEB_${type}_ANN_${dt}_lonlat.nc

cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_SAM_${type}_MON_${dt}_lonlat.nc
cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_${type}_ANN_${dt}_lonlat.nc ${var}_SAM_${type}_ANN_${dt}_lonlat.nc

cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_LPB_${type}_MON_${dt}_lonlat.nc
cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_${type}_ANN_${dt}_lonlat.nc ${var}_LPB_${type}_ANN_${dt}_lonlat.nc

cdo sellonlatbox,-73.85,-34.85,-33.85,5.85 ${var}_SA_${type}_MON_${dt}_lonlat.nc ${var}_BR_${type}_MON_${dt}_lonlat.nc
cdo sellonlatbox,-73.85,-34.85,-33.85,5.85 ${var}_SA_${type}_ANN_${dt}_lonlat.nc ${var}_BR_${type}_ANN_${dt}_lonlat.nc

echo 
echo "6. Deleting file"
rm ${var}_${type}_MON_${dt}.nc
rm ${var}_SA_${type}_MON_${dt}.nc

echo
echo "--------------- END POSPROCESSING OBS ----------------"















