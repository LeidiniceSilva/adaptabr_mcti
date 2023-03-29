#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Mar 01, 2023'
#__description__ = 'Posprocessing the CRU-ts4.06 obs dataset'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------"

# Variable 
var='tmn'  

# Date
dt='1980-2014'

path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs/"
cd ${path}

echo
echo ${path}

echo 
echo "1. Select date"
cdo seldate,1980-01-01T00:00:00,2014-12-31T00:00:00 cru_ts4.06.1901.2021.${var}.dat.nc ${var}_CRU_ts4_${dt}.nc

echo
echo "2. Conventing calendar"
cdo setcalendar,standard ${var}_CRU_ts4_${dt}.nc ${var}_CRU_ts4_mon_${dt}.nc

echo
echo "3. Conventing grade"
/home/nice/Documentos/github_projects/shell/regcm_pos/./regrid ${var}_CRU_ts4_mon_${dt}.nc -60,15,0.5 -85,-30,0.5 bil

echo 
echo "4. Conventing unit"
if [ ${var} == 'pre' ]
then
cdo divc,30.5 ${var}_CRU_ts4_mon_${dt}_lonlat.nc ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc
else
cp ${var}_CRU_ts4_mon_${dt}_lonlat.nc ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc
fi

echo
echo "5. Calculate periods"
cdo -yearavg ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SA_CRU_ts4_ANN_${dt}_lonlat.nc
cdo -r -timselavg,3 -selmon,1,2,12 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SA_CRU_ts4_DJF_${dt}_lonlat.nc
cdo -r -timselavg,3 -selmon,3,4,5 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SA_CRU_ts4_MAM_${dt}_lonlat.nc
cdo -r -timselavg,3 -selmon,6,7,8 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SA_CRU_ts4_JJA_${dt}_lonlat.nc
cdo -r -timselavg,3 -selmon,9,10,11 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SA_CRU_ts4_SON_${dt}_lonlat.nc

echo
echo "6. Select subregion"
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_NAMZ_CRU_ts4_MON_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_CRU_ts4_ANN_${dt}_lonlat.nc ${var}_NAMZ_CRU_ts4_ANN_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_CRU_ts4_DJF_${dt}_lonlat.nc ${var}_NAMZ_CRU_ts4_DJF_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_CRU_ts4_MAM_${dt}_lonlat.nc ${var}_NAMZ_CRU_ts4_MAM_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_CRU_ts4_JJA_${dt}_lonlat.nc ${var}_NAMZ_CRU_ts4_JJA_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-5,5 ${var}_SA_CRU_ts4_SON_${dt}_lonlat.nc ${var}_NAMZ_CRU_ts4_SON_${dt}_lonlat.nc

cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SAMZ_CRU_ts4_MON_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_CRU_ts4_ANN_${dt}_lonlat.nc ${var}_SAMZ_CRU_ts4_ANN_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_CRU_ts4_DJF_${dt}_lonlat.nc ${var}_SAMZ_CRU_ts4_DJF_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_CRU_ts4_MAM_${dt}_lonlat.nc ${var}_SAMZ_CRU_ts4_MAM_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_CRU_ts4_JJA_${dt}_lonlat.nc ${var}_SAMZ_CRU_ts4_JJA_${dt}_lonlat.nc
cdo sellonlatbox,-70,-45,-12.5,-5 ${var}_SA_CRU_ts4_SON_${dt}_lonlat.nc ${var}_SAMZ_CRU_ts4_SON_${dt}_lonlat.nc

cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_NEB_CRU_ts4_MON_${dt}_lonlat.nc
cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_CRU_ts4_ANN_${dt}_lonlat.nc ${var}_NEB_CRU_ts4_ANN_${dt}_lonlat.nc
cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_CRU_ts4_DJF_${dt}_lonlat.nc ${var}_NEB_CRU_ts4_DJF_${dt}_lonlat.nc
cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_CRU_ts4_MAM_${dt}_lonlat.nc ${var}_NEB_CRU_ts4_MAM_${dt}_lonlat.nc
cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_CRU_ts4_JJA_${dt}_lonlat.nc ${var}_NEB_CRU_ts4_JJA_${dt}_lonlat.nc
cdo sellonlatbox,-45,-34,-15,-2 ${var}_SA_CRU_ts4_SON_${dt}_lonlat.nc ${var}_NEB_CRU_ts4_SON_${dt}_lonlat.nc

cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_SAM_CRU_ts4_MON_${dt}_lonlat.nc
cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_CRU_ts4_ANN_${dt}_lonlat.nc ${var}_SAM_CRU_ts4_ANN_${dt}_lonlat.nc
cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_CRU_ts4_DJF_${dt}_lonlat.nc ${var}_SAM_CRU_ts4_DJF_${dt}_lonlat.nc
cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_CRU_ts4_MAM_${dt}_lonlat.nc ${var}_SAM_CRU_ts4_MAM_${dt}_lonlat.nc
cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_CRU_ts4_JJA_${dt}_lonlat.nc ${var}_SAM_CRU_ts4_JJA_${dt}_lonlat.nc
cdo sellonlatbox,-55,-45,-20,-10 ${var}_SA_CRU_ts4_SON_${dt}_lonlat.nc ${var}_SAM_CRU_ts4_SON_${dt}_lonlat.nc

cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_CRU_ts4_MON_${dt}_lonlat.nc ${var}_LPB_CRU_ts4_MON_${dt}_lonlat.nc
cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_CRU_ts4_ANN_${dt}_lonlat.nc ${var}_LPB_CRU_ts4_ANN_${dt}_lonlat.nc
cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_CRU_ts4_DJF_${dt}_lonlat.nc ${var}_LPB_CRU_ts4_DJF_${dt}_lonlat.nc
cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_CRU_ts4_MAM_${dt}_lonlat.nc ${var}_LPB_CRU_ts4_MAM_${dt}_lonlat.nc
cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_CRU_ts4_JJA_${dt}_lonlat.nc ${var}_LPB_CRU_ts4_JJA_${dt}_lonlat.nc
cdo sellonlatbox,-60,-45,-35,-20 ${var}_SA_CRU_ts4_SON_${dt}_lonlat.nc ${var}_LPB_CRU_ts4_SON_${dt}_lonlat.nc

echo 
echo "6. Deleting file"
rm ${var}_CRU_ts4_${dt}.nc
rm ${var}_CRU_ts4_mon_${dt}.nc
rm ${var}_CRU_ts4_mon_${dt}_lonlat.nc

echo
echo "--------------- END POSPROCESSING OBS ----------------"















