# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot annual bias maps of cmip6 models"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from dict_cmip6_models_name import cmip6


def import_obs(param, date, min_lon, min_lat, max_lon, max_lat):

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/project/database/obs/' + '{0}_SA_BR-DWGD_UFES_UTEXAS_v_3.0_MON_{1}_lonlat.nc'.format(param, date))
	data = arq[param]
	data = data.sel(lat=slice(min_lat,max_lat), lon=slice(min_lon,max_lon))
	time = data.sel(time=slice('1986-01-01','2005-12-31'))
	var = time.groupby('time.month').mean('time')
	mean = np.nanmean(np.nanmean(var.values, axis=1), axis=1)
	
	return mean
	
	
def import_cmip(param, model, exp, date, min_lon, min_lat, max_lon, max_lat):

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/project/database/cmip6/' + '{0}_SA_{1}_historical_{2}_MON_{3}_lonlat.nc'.format(param, model, exp, date))
	data = arq[param]
	data = data.sel(lat=slice(min_lat,max_lat), lon=slice(min_lon,max_lon))
	time = data.sel(time=slice('1986-01-01','2005-12-31'))
	var = time.groupby('time.month').mean('time')
	mean = np.nanmean(np.nanmean(var.values, axis=1), axis=1)
	
	return mean
	

# Import cmip models and obs database 
var_obs = 'pr'
var_cmip6 = 'pr'
dt = '1986-2005'

namz_min_lon, namz_min_lat, namz_max_lon, namz_max_lat = -70., -5., -45., 5. 
samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat = -70., -12.5, -45., -5.
neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat = -45,-15,-34,-2
sam_min_lon, sam_min_lat, sam_max_lon, sam_max_lat = -55,-20,-45,-10
lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat = -60,-35,-45,-20
br_min_lon, br_min_lat, br_max_lon, br_max_lat = -73.85,-33.85,-34.85,5.85

clim_namz_obs = import_obs(var_obs, dt, namz_min_lon, namz_min_lat, namz_max_lon, namz_max_lat)
clim_samz_obs = import_obs(var_obs, dt, samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat)
clim_neb_obs = import_obs(var_obs, dt, neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat)
clim_sam_obs = import_obs(var_obs, dt, sam_min_lon, sam_min_lat, sam_max_lon, sam_max_lat)
clim_lpb_obs = import_obs(var_obs, dt, lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat)
clim_br_obs = import_obs(var_obs, dt, br_min_lon, br_min_lat, br_max_lon, br_max_lat)

clim_namz_cmip6 = []
clim_samz_cmip6 = []
clim_neb_cmip6 = []
clim_sam_cmip6 = []
clim_lpb_cmip6 = []
clim_br_cmip6 = []

best_models = [17, 7, 13, 9, 15]
for i in best_models:

	clim_namz_cmip6.append(import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt, namz_min_lon, namz_min_lat, namz_max_lon, namz_max_lat))	
	clim_samz_cmip6.append(import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt, samz_min_lon, samz_min_lat, samz_max_lon, samz_max_lat))
	clim_neb_cmip6.append(import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt, neb_min_lon, neb_min_lat, neb_max_lon, neb_max_lat))
	clim_sam_cmip6.append(import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt, sam_min_lon, sam_min_lat, sam_max_lon, sam_max_lat))
	clim_lpb_cmip6.append(import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt, lpb_min_lon, lpb_min_lat, lpb_max_lon, lpb_max_lat))
	clim_br_cmip6.append(import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt, br_min_lon, br_min_lat, br_max_lon, br_max_lat))
	
# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 7))
time = np.arange(0.5, 12 + 0.5)

ax = fig.add_subplot(3, 2, 1)  
annual_cycle = ax.plot(time, clim_namz_cmip6[0], time, clim_namz_cmip6[1], time, clim_namz_cmip6[2], 
time, clim_namz_cmip6[3], time, clim_namz_cmip6[4], time, clim_namz_obs)
plt.title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 45, 3), fontsize=8)
	plt.ylim(18, 42)
else:
	plt.yticks(np.arange(9, 36, 3), fontsize=8)
	plt.ylim(9, 33)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='green')
plt.setp(l4, color='orange')
plt.setp(l5, color='gray')
plt.setp(l6, color='black')

ax = fig.add_subplot(3, 2, 2)
annual_cycle = ax.plot(time, clim_samz_cmip6[0], time, clim_samz_cmip6[1], time, clim_samz_cmip6[2], 
time, clim_samz_cmip6[3], time, clim_samz_cmip6[4], time, clim_samz_obs)
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 45, 3), fontsize=8)
	plt.ylim(18, 42)
else:
	plt.yticks(np.arange(9, 36, 3), fontsize=8)
	plt.ylim(9, 33)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='green')
plt.setp(l4, color='orange')
plt.setp(l5, color='gray')
plt.setp(l6, color='black')

ax = fig.add_subplot(3, 2, 3)
annual_cycle = ax.plot(time, clim_neb_cmip6[0], time, clim_neb_cmip6[1], time, clim_neb_cmip6[2], 
time, clim_neb_cmip6[3], time, clim_neb_cmip6[4], time, clim_neb_obs)
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
if var_cmip6 == 'pr':
	plt.ylabel('Precipitação (mm d⁻¹)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.ylabel('Temperatura máxima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(18, 45, 3), fontsize=8)
	plt.ylim(18, 42)
else:
	plt.ylabel('Temperatura mínima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(9, 36, 3), fontsize=8)
	plt.ylim(9, 33)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='green')
plt.setp(l4, color='orange')
plt.setp(l5, color='gray')
plt.setp(l6, color='black')

ax = fig.add_subplot(3, 2, 4)
annual_cycle = ax.plot(time, clim_sam_cmip6[0], time, clim_sam_cmip6[1], time, clim_sam_cmip6[2], 
time, clim_sam_cmip6[3], time, clim_sam_cmip6[4], time, clim_sam_obs)
plt.title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
if var_cmip6 == 'pr':
	plt.ylabel('Precipitação (mm d⁻¹)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.ylabel('Temperatura máxima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(18, 45, 3), fontsize=8)
	plt.ylim(18, 42)
else:
	plt.ylabel('Temperatura mínima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(9, 36, 3), fontsize=8)
	plt.ylim(9, 33)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='green')
plt.setp(l4, color='orange')
plt.setp(l5, color='gray')
plt.setp(l6, color='black')

ax = fig.add_subplot(3, 2, 5)
annual_cycle = ax.plot(time, clim_lpb_cmip6[0], time, clim_lpb_cmip6[1], time, clim_lpb_cmip6[2], 
time, clim_lpb_cmip6[3], time, clim_lpb_cmip6[4], time, clim_lpb_obs)
plt.title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Meses', fontsize=8, fontweight='bold')
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 45, 3), fontsize=8)
	plt.ylim(18, 42)
else:
	plt.yticks(np.arange(9, 36, 3), fontsize=8)
	plt.ylim(9, 33)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='green')
plt.setp(l4, color='orange')
plt.setp(l5, color='gray')
plt.setp(l6, color='black')

ax = fig.add_subplot(3, 2, 6)
annual_cycle = ax.plot(time, clim_br_cmip6[0], time, clim_br_cmip6[1], time, clim_br_cmip6[2], 
time, clim_br_cmip6[3], time, clim_br_cmip6[4], time, clim_br_obs)
plt.title(u'(f) BR', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.xlabel('Meses', fontsize=8, fontweight='bold')
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2), fontsize=8)
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 45, 3), fontsize=8)
	plt.ylim(18, 42)
else:
	plt.yticks(np.arange(9, 36, 3), fontsize=8)
	plt.ylim(9, 33)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6 = annual_cycle
plt.setp(l1, color='blue')
plt.setp(l2, color='red')
plt.setp(l3, color='green')
plt.setp(l4, color='orange')
plt.setp(l5, color='gray')
plt.setp(l6, color='black')

legend = [cmip6[17][0],cmip6[7][0],cmip6[13][0],cmip6[9][0],cmip6[15][0], 'BR-DWGD']
plt.legend(annual_cycle, legend, ncol=6, loc=(-1.25, -0.5), fontsize=8)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.20, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/project/figs/figs_report-II'
name_out = 'pyplt_annual_cycle_best_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()
