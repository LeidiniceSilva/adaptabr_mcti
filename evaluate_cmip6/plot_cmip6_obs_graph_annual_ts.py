# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot annual cycle of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from dict_cmip6_models_name import cmip6


def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_ANN_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return obs_data

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_ANN_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return mdl_data
	              
               
# Import cmip models and obs database 
var_obs = 'Tmin'
var_cmip6 = 'tasmin'
dt = '1986-2005'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)
clim_br_obs = import_obs(var_obs, 'BR', dt)

clim_namz_cmip6 = []
clim_samz_cmip6 = []
clim_neb_cmip6 = []
clim_sam_cmip6 = []
clim_lpb_cmip6 = []
clim_br_cmip6 = []

for i in range(1, 18):

	clim_namz_cmip6.append(import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt))	
	clim_samz_cmip6.append(import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt))
	clim_neb_cmip6.append(import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt))
	clim_sam_cmip6.append(import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt))
	clim_lpb_cmip6.append(import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt))
	clim_br_cmip6.append(import_cmip(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], dt))

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 7))
time = np.arange(0.5, 20 + 0.5)

ax = fig.add_subplot(3, 2, 1)  
annual_cycle = ax.plot(time, clim_namz_cmip6[0], time, clim_namz_cmip6[1], time, clim_namz_cmip6[2], 
time, clim_namz_cmip6[3], time, clim_namz_cmip6[4], time, clim_namz_cmip6[5], time, clim_namz_cmip6[6], 
time, clim_namz_cmip6[7], time, clim_namz_cmip6[8], time, clim_namz_cmip6[9], time, clim_namz_cmip6[10], 
time, clim_namz_cmip6[11], time, clim_namz_cmip6[12], time, clim_namz_cmip6[13], time, clim_namz_cmip6[14],
time, clim_namz_cmip6[15], time, clim_namz_cmip6[16], time, clim_namz_obs)
plt.title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('1986', '', '1988', '', '1990', '', '1992', '', '1994', '', '1996', '', '1998', '', '2000', '', '2002', '', '2004', ''), fontsize=8)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 10, 1), fontsize=8)
	plt.ylim(0, 9)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(19, 41, 2), fontsize=8)
	plt.ylim(19, 39)
else:
	plt.yticks(np.arange(9, 31, 2), fontsize=8)
	plt.ylim(9, 29)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18 = annual_cycle
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18, color='black')

ax = fig.add_subplot(3, 2, 2)
annual_cycle = ax.plot(time, clim_samz_cmip6[0], time, clim_samz_cmip6[1], time, clim_samz_cmip6[2], 
time, clim_samz_cmip6[3], time, clim_samz_cmip6[4], time, clim_samz_cmip6[5], time, clim_samz_cmip6[6], 
time, clim_samz_cmip6[7], time, clim_samz_cmip6[8], time, clim_samz_cmip6[9], time, clim_samz_cmip6[10], 
time, clim_samz_cmip6[11], time, clim_samz_cmip6[12], time, clim_samz_cmip6[13], time, clim_samz_cmip6[14],
time, clim_samz_cmip6[15], time, clim_samz_cmip6[16], time, clim_samz_obs)
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('1986', '', '1988', '', '1990', '', '1992', '', '1994', '', '1996', '', '1998', '', '2000', '', '2002', '', '2004', ''), fontsize=8)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 10, 1), fontsize=8)
	plt.ylim(0, 9)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(19, 41, 2), fontsize=8)
	plt.ylim(19, 39)
else:
	plt.yticks(np.arange(9, 31, 2), fontsize=8)
	plt.ylim(9, 29)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18 = annual_cycle
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18, color='black')

ax = fig.add_subplot(3, 2, 3)
annual_cycle = ax.plot(time, clim_neb_cmip6[0], time, clim_neb_cmip6[1], time, clim_neb_cmip6[2], 
time, clim_neb_cmip6[3], time, clim_neb_cmip6[4], time, clim_neb_cmip6[5], time, clim_neb_cmip6[6], 
time, clim_neb_cmip6[7], time, clim_neb_cmip6[8], time, clim_neb_cmip6[9], time, clim_neb_cmip6[10], 
time, clim_neb_cmip6[11], time, clim_neb_cmip6[12], time, clim_neb_cmip6[13], time, clim_neb_cmip6[14],
time, clim_neb_cmip6[15], time, clim_neb_cmip6[16], time, clim_neb_obs)
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('1986', '', '1988', '', '1990', '', '1992', '', '1994', '', '1996', '', '1998', '', '2000', '', '2002', '', '2004', ''), fontsize=8)
if var_cmip6 == 'pr':
	plt.ylabel('Precipitação (mm d⁻¹)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 10, 1), fontsize=8)
	plt.ylim(0, 9)
elif var_cmip6 == 'tasmax':
	plt.ylabel('Temperatura máxima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(19, 41, 2), fontsize=8)
	plt.ylim(19, 39)
else:
	plt.ylabel('Temperatura mínima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(9, 31, 2), fontsize=8)
	plt.ylim(9, 29)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18 = annual_cycle
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18, color='black')

ax = fig.add_subplot(3, 2, 4)
annual_cycle = ax.plot(time, clim_sam_cmip6[0], time, clim_sam_cmip6[1], time, clim_sam_cmip6[2], 
time, clim_sam_cmip6[3], time, clim_sam_cmip6[4], time, clim_sam_cmip6[5], time, clim_sam_cmip6[6], 
time, clim_sam_cmip6[7], time, clim_sam_cmip6[8], time, clim_sam_cmip6[9], time, clim_sam_cmip6[10], 
time, clim_sam_cmip6[11], time, clim_sam_cmip6[12], time, clim_sam_cmip6[13], time, clim_sam_cmip6[14],
time, clim_sam_cmip6[15], time, clim_sam_cmip6[16], time, clim_sam_obs)
plt.title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('1986', '', '1988', '', '1990', '', '1992', '', '1994', '', '1996', '', '1998', '', '2000', '', '2002', '', '2004', ''), fontsize=8)
if var_cmip6 == 'pr':
	plt.ylabel('Precipitação (mm d⁻¹)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 10, 1), fontsize=8)
	plt.ylim(0, 9)
elif var_cmip6 == 'tasmax':
	plt.ylabel('Temperatura máxima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(19, 41, 2), fontsize=8)
	plt.ylim(19, 39)
else:
	plt.ylabel('Temperatura mínima (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(9, 31, 2), fontsize=8)
	plt.ylim(9, 29)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18 = annual_cycle
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18, color='black')

ax = fig.add_subplot(3, 2, 5)
annual_cycle = ax.plot(time, clim_lpb_cmip6[0], time, clim_lpb_cmip6[1], time, clim_lpb_cmip6[2], 
time, clim_lpb_cmip6[3], time, clim_lpb_cmip6[4], time, clim_lpb_cmip6[5], time, clim_lpb_cmip6[6], 
time, clim_lpb_cmip6[7], time, clim_lpb_cmip6[8], time, clim_lpb_cmip6[9], time, clim_lpb_cmip6[10], 
time, clim_lpb_cmip6[11], time, clim_lpb_cmip6[12], time, clim_lpb_cmip6[13], time, clim_lpb_cmip6[14],
time, clim_lpb_cmip6[15], time, clim_lpb_cmip6[16], time, clim_lpb_obs)
plt.title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('1986', '', '1988', '', '1990', '', '1992', '', '1994', '', '1996', '', '1998', '', '2000', '', '2002', '', '2004', ''), fontsize=8)
plt.xlabel('Years', fontsize=8, fontweight='bold')
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 10, 1), fontsize=8)
	plt.ylim(0, 9)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(19, 41, 2), fontsize=8)
	plt.ylim(19, 39)
else:
	plt.yticks(np.arange(9, 31, 2), fontsize=8)
	plt.ylim(9, 29)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18 = annual_cycle
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18, color='black')

ax = fig.add_subplot(3, 2, 6)
annual_cycle = ax.plot(time, clim_br_cmip6[0], time, clim_br_cmip6[1], time, clim_br_cmip6[2], 
time, clim_br_cmip6[3], time, clim_br_cmip6[4], time, clim_br_cmip6[5], time, clim_br_cmip6[6], 
time, clim_br_cmip6[7], time, clim_br_cmip6[8], time, clim_br_cmip6[9], time, clim_br_cmip6[10], 
time, clim_br_cmip6[11], time, clim_br_cmip6[12], time, clim_br_cmip6[13], time, clim_br_cmip6[14],
time, clim_br_cmip6[15], time, clim_br_cmip6[16], time, clim_br_obs)
plt.title(u'(f) BR', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('1986', '', '1988', '', '1990', '', '1992', '', '1994', '', '1996', '', '1998', '', '2000', '', '2002', '', '2004', ''), fontsize=8)
plt.xlabel('Years', fontsize=8, fontweight='bold')
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 10, 1), fontsize=8)
	plt.ylim(0, 9)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(19, 41, 2), fontsize=8)
	plt.ylim(19, 39)
else:
	plt.yticks(np.arange(9, 31, 2), fontsize=8)
	plt.ylim(9, 29)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18 = annual_cycle
plt.setp(l1)
plt.setp(l2)
plt.setp(l3)
plt.setp(l4)
plt.setp(l5)
plt.setp(l6)
plt.setp(l7)
plt.setp(l8)
plt.setp(l9)
plt.setp(l10)
plt.setp(l11)
plt.setp(l12)
plt.setp(l13)
plt.setp(l14)
plt.setp(l15)
plt.setp(l16)
plt.setp(l17)
plt.setp(l18, color='black')

legend = [cmip6[1][0],cmip6[2][0],cmip6[3][0],cmip6[4][0],cmip6[5][0],cmip6[6][0],
cmip6[7][0],cmip6[8][0],cmip6[9][0],cmip6[10][0],cmip6[11][0],cmip6[12][0],cmip6[13][0],
cmip6[14][0],cmip6[15][0],cmip6[16][0],cmip6[17][0],'BR-DWGD']
plt.legend(annual_cycle, legend, ncol=6, loc=(-1.35, -0.7), fontsize=8)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.20, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_annual_ts_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






