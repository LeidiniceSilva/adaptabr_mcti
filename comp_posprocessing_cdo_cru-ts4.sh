#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = '03/13/2023'
#__description__ = 'Posprocessing the CRU-ts4.06 dataset'

echo
echo "--------------- INIT POSPROCESSING OBS ----------------"

# Variable name
var = 'pr'     

path="/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs"
cd ${path}

echo
echo ${path}

echo 
echo "1. Select date"
cdo seldate,1981-01-01T00:00:00,2014-12-31T00:00:00 ${var}_ts4.06.1901.2021.pre.dat.nc ${var}_CRU_ts4_1981-2014.nc

echo
echo "2. Convert calendar"
cdo setcalendar,standard ${var}_CRU_ts4_1981-2014.nc ${var}_CRU_ts4_mon_1981-2014.nc

echo
echo "3. Remapbil"
/home/nice/Documents/github_projects/shell/regcm_pos/./regrid ${var}_CRU_ts4_mon_1981-2014.nc -20,10,0.25 -85,-15,0.25 bil

echo 
echo "4. Convention unit"

# Member name
if [ ${var} == 'pr' ]
then
cdo divc,30.5 ${var}_CRU_ts4_mon_1981-2014_lonlat.nc ${var}_SA_CRU_ts4_mon_1981-2014_lonlat.nc
else
cp ${var}_CRU_ts4_mon_1981-2014_lonlat.nc ${var}_SA_CRU_ts4_mon_1981-2014_lonlat.nc
fi

echo
echo "5. Select subdomain"
cdo sellonlatbox,-68,-52,-12,-3 pre_SA_CRU_ts4_mon_1981-2014_lonlat.nc pre_SAMZ_CRU_ts4_mon_1981-2014_lonlat.nc
cdo sellonlatbox,-40,-35,-16,-3 pre_SA_CRU_ts4_mon_1981-2014_lonlat.nc pre_NAMZ_CRU_ts4_mon_1981-2014_lonlat.nc
cdo sellonlatbox,-68,-52,-12,-3 pre_SA_CRU_ts4_mon_1981-2014_lonlat.nc pre_SAM_CRU_ts4_mon_1981-2014_lonlat.nc
cdo sellonlatbox,-40,-35,-16,-3 pre_SA_CRU_ts4_mon_1981-2014_lonlat.nc pre_NEB_CRU_ts4_mon_1981-2014_lonlat.nc
cdo sellonlatbox,-50.5,-42.5,-15,-2.5 pre_SA_CRU_ts4_mon_1981-2014_lonlat.nc pre_LPB_CRU_ts4_mon_1981-2014_lonlat.nc

echo 
echo "6. Deleting file"
rm ${var}_CRU_ts4_1981-2014.nc
rm ${var}_CRU_ts4_mon_1981-2014.nc
rm ${var}_CRU_ts4_mon_1981-2014_lonlat.nc

echo
echo "--------------- END POSPROCESSING OBS ----------------"















