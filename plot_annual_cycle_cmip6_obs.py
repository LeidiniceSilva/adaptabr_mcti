# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 14, 2022"
__description__ = "This script plot annual cycle from cmip6"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from dict_name_cmip6 import cmip6


def import_obs(param, area):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_cru_ts4.05_obs_mon_1961-2014_lonlat.nc'.format(path, param, area)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	obs_data = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)
	obs_clim = []
	for mon in range(0, 11 + 1):
		obs = np.nanmean(obs_data[mon::12], axis=0)
		obs_clim.append(obs)
			
	return obs_clim
	

def import_cmip(path, param, area, cmip, exp, date):
	
	arq   = '{0}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat_mask.nc'.format(path, param, area, cmip, exp, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	sim_data = np.nanmean(np.nanmean(value[:,:,:], axis=1), axis=1)
	sim_clim = []
	for mon in range(0, 11 + 1):
		sim = np.nanmean(sim_data[mon::12], axis=0)
		sim_clim.append(sim)
	
	return sim_clim	
	              
               
# Import cmip models and obs database 
cmip5_path = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip/cmip5'
cmip6_path = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip/cmip6'

pre_clim_samz_cru = import_obs(u'pre', 'SAMZ')
pre_clim_slpb_cru = import_obs(u'pre', 'SLPB')
pre_clim_neb_cru  = import_obs(u'pre', 'NEB')

tas_clim_samz_cru = import_obs(u'tmp', 'SAMZ')
tas_clim_slpb_cru = import_obs(u'tmp', 'SLPB')
tas_clim_neb_cru  = import_obs(u'tmp', 'NEB')

pre_clim_samz_cmip5 = import_cmip(cmip5_path, u'pr', 'SAMZ', u'ensmean_cmip5', 'r1i1p1', '1961-2005')
pre_clim_slpb_cmip5 = import_cmip(cmip5_path, u'pr', 'SLPB', u'ensmean_cmip5', 'r1i1p1', '1961-2005')
pre_clim_neb_cmip5  = import_cmip(cmip5_path, u'pr', 'NEB', u'ensmean_cmip5', 'r1i1p1', '1961-2005')

tas_clim_samz_cmip5 = import_cmip(cmip5_path, u'tas', 'SAMZ', u'ensmean_cmip5', 'r1i1p1', '1961-2005')
tas_clim_slpb_cmip5 = import_cmip(cmip5_path, u'tas', 'SLPB', u'ensmean_cmip5', 'r1i1p1', '1961-2005')
tas_clim_neb_cmip5  = import_cmip(cmip5_path, u'tas', 'NEB', u'ensmean_cmip5', 'r1i1p1', '1961-2005')

pre_clim_samz_cmip6 = []
pre_clim_slpb_cmip6 = []
pre_clim_neb_cmip6 = []
tas_clim_samz_cmip6 = []
tas_clim_slpb_cmip6 = []
tas_clim_neb_cmip6 = []
legend = []

for i in range(1, 24):

	pre_clim_samz_cmip6.append(import_cmip(cmip6_path, u'pr', 'SAMZ', cmip6[i][0], cmip6[i][1], '1961-2014'))
	pre_clim_slpb_cmip6.append(import_cmip(cmip6_path, u'pr', 'SLPB', cmip6[i][0], cmip6[i][1], '1961-2014'))
	pre_clim_neb_cmip6.append(import_cmip(cmip6_path, u'pr', 'NEB', cmip6[i][0], cmip6[i][1], '1961-2014'))
	tas_clim_samz_cmip6.append(import_cmip(cmip6_path, u'tas', 'SAMZ', cmip6[i][0], cmip6[i][1], '1961-2014'))
	tas_clim_slpb_cmip6.append(import_cmip(cmip6_path, u'tas', 'SLPB', cmip6[i][0], cmip6[i][1], '1961-2014'))
	tas_clim_neb_cmip6.append(import_cmip(cmip6_path, u'tas', 'NEB', cmip6[i][0], cmip6[i][1], '1961-2014'))
	legend.append(cmip6[i][0])

# Plot regcm exps and obs database 
fig = plt.figure(figsize=(10, 8))
time = np.arange(0.5, 12 + 0.5)
colors = plt.matplotlib.cm.gist_rainbow(np.linspace(0.1, 0.6, 26))

ax = fig.add_subplot(3, 2, 1)  
annual_cycle = ax.plot(time, pre_clim_samz_cmip6[0], time, pre_clim_samz_cmip6[1], time, pre_clim_samz_cmip6[2], 
time, pre_clim_samz_cmip6[3], time, pre_clim_samz_cmip6[4], time, pre_clim_samz_cmip6[5], time, pre_clim_samz_cmip6[6], 
time, pre_clim_samz_cmip6[7], time, pre_clim_samz_cmip6[8], time, pre_clim_samz_cmip6[9], time, pre_clim_samz_cmip6[10], 
time, pre_clim_samz_cmip6[11], time, pre_clim_samz_cmip6[12], time, pre_clim_samz_cmip6[13], time, pre_clim_samz_cmip6[14],
time, pre_clim_samz_cmip6[15], time, pre_clim_samz_cmip6[16], time, pre_clim_samz_cmip6[17], time, pre_clim_samz_cmip6[18], 
time, pre_clim_samz_cmip6[19], time, pre_clim_samz_cmip6[20], time, pre_clim_samz_cmip6[21], time, pre_clim_samz_cmip6[22], 
time, pre_clim_samz_cmip5, time, pre_clim_samz_cru)
plt.title(u'(a) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 18, 2), fontsize=8)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25 = annual_cycle
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
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23, color='black', linestyle='--')
plt.setp(l24, color='black', linestyle='-.')
plt.setp(l25, color='black')

ax = fig.add_subplot(3, 2, 2)
annual_cycle = ax.plot(time, tas_clim_samz_cmip6[0], time, tas_clim_samz_cmip6[1], time, tas_clim_samz_cmip6[2], 
time, tas_clim_samz_cmip6[3], time, tas_clim_samz_cmip6[4], time, tas_clim_samz_cmip6[5], time, tas_clim_samz_cmip6[6], 
time, tas_clim_samz_cmip6[7], time, tas_clim_samz_cmip6[8], time, tas_clim_samz_cmip6[9], time, tas_clim_samz_cmip6[10], 
time, tas_clim_samz_cmip6[11], time, tas_clim_samz_cmip6[12], time, tas_clim_samz_cmip6[13], time, tas_clim_samz_cmip6[14],
time, tas_clim_samz_cmip6[15], time, tas_clim_samz_cmip6[16], time, tas_clim_samz_cmip6[17], time, tas_clim_samz_cmip6[18], 
time, tas_clim_samz_cmip6[19], time, tas_clim_samz_cmip6[20], time, tas_clim_samz_cmip6[21], time, tas_clim_samz_cmip6[22], 
time, tas_clim_samz_cmip5, time, tas_clim_samz_cru)
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(10, 36, 4), fontsize=8)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25 = annual_cycle
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
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23, color='black', linestyle='--')
plt.setp(l24, color='black', linestyle='-.')
plt.setp(l25, color='black')

ax = fig.add_subplot(3, 2, 3)
annual_cycle = ax.plot(time, pre_clim_neb_cmip6[0], time, pre_clim_neb_cmip6[1], time, pre_clim_neb_cmip6[2], 
time, pre_clim_neb_cmip6[3], time, pre_clim_neb_cmip6[4], time, pre_clim_neb_cmip6[5], time, pre_clim_neb_cmip6[6], 
time, pre_clim_neb_cmip6[7], time, pre_clim_neb_cmip6[8], time, pre_clim_neb_cmip6[9], time, pre_clim_neb_cmip6[10], 
time, pre_clim_neb_cmip6[11], time, pre_clim_neb_cmip6[12], time, pre_clim_neb_cmip6[13], time, pre_clim_neb_cmip6[14],
time, pre_clim_neb_cmip6[15], time, pre_clim_neb_cmip6[16], time, pre_clim_neb_cmip6[17], time, pre_clim_neb_cmip6[18], 
time, pre_clim_neb_cmip6[19], time, pre_clim_neb_cmip6[20], time, pre_clim_neb_cmip6[21], time, pre_clim_neb_cmip6[22], 
time, pre_clim_neb_cmip5, time, pre_clim_neb_cru)
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Precipitation (mm d$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 18, 2), fontsize=8)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25 = annual_cycle
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
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23, color='black', linestyle='--')
plt.setp(l24, color='black', linestyle='-.')
plt.setp(l25, color='black')

ax = fig.add_subplot(3, 2, 4)
annual_cycle = ax.plot(time, tas_clim_neb_cmip6[0], time, tas_clim_neb_cmip6[1], time, tas_clim_neb_cmip6[2], 
time, tas_clim_neb_cmip6[3], time, tas_clim_neb_cmip6[4], time, tas_clim_neb_cmip6[5], time, tas_clim_neb_cmip6[6], 
time, tas_clim_neb_cmip6[7], time, tas_clim_neb_cmip6[8], time, tas_clim_neb_cmip6[9], time, tas_clim_neb_cmip6[10], 
time, tas_clim_neb_cmip6[11], time, tas_clim_neb_cmip6[12], time, tas_clim_neb_cmip6[13], time, tas_clim_neb_cmip6[14],
time, tas_clim_neb_cmip6[15], time, tas_clim_neb_cmip6[16], time, tas_clim_neb_cmip6[17], time, tas_clim_neb_cmip6[18], 
time, tas_clim_neb_cmip6[19], time, tas_clim_neb_cmip6[20], time, tas_clim_neb_cmip6[21], time, tas_clim_neb_cmip6[22], 
time, tas_clim_neb_cmip5, time, tas_clim_neb_cru)
plt.title(u'(d) NEB', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Temperature (°C)', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(10, 36, 4), fontsize=8)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25 = annual_cycle
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
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23, color='black', linestyle='--')
plt.setp(l24, color='black', linestyle='-.')
plt.setp(l25, color='black')

ax = fig.add_subplot(3, 2, 5)
annual_cycle = ax.plot(time, pre_clim_slpb_cmip6[0], time, pre_clim_slpb_cmip6[1], time, pre_clim_slpb_cmip6[2], 
time, pre_clim_slpb_cmip6[3], time, pre_clim_slpb_cmip6[4], time, pre_clim_slpb_cmip6[5], time, pre_clim_slpb_cmip6[6], 
time, pre_clim_slpb_cmip6[7], time, pre_clim_slpb_cmip6[8], time, pre_clim_slpb_cmip6[9], time, pre_clim_slpb_cmip6[10], 
time, pre_clim_slpb_cmip6[11], time, pre_clim_slpb_cmip6[12], time, pre_clim_slpb_cmip6[13], time, pre_clim_slpb_cmip6[14],
time, pre_clim_slpb_cmip6[15], time, pre_clim_slpb_cmip6[16], time, pre_clim_slpb_cmip6[17], time, pre_clim_slpb_cmip6[18], 
time, pre_clim_slpb_cmip6[19], time, pre_clim_slpb_cmip6[20], time, pre_clim_slpb_cmip6[21], time, pre_clim_slpb_cmip6[22], 
time, pre_clim_slpb_cmip5, time, pre_clim_slpb_cru)
plt.title(u'(e) SLPB', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('Months', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 18, 2), fontsize=8)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25 = annual_cycle
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
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23, color='black', linestyle='--')
plt.setp(l24, color='black', linestyle='-.')
plt.setp(l25, color='black')

ax = fig.add_subplot(3, 2, 6)
annual_cycle = ax.plot(time, tas_clim_slpb_cmip6[0], time, tas_clim_slpb_cmip6[1], time, tas_clim_slpb_cmip6[2], 
time, tas_clim_slpb_cmip6[3], time, tas_clim_slpb_cmip6[4], time, tas_clim_slpb_cmip6[5], time, tas_clim_slpb_cmip6[6], 
time, tas_clim_slpb_cmip6[7], time, tas_clim_slpb_cmip6[8], time, tas_clim_slpb_cmip6[9], time, tas_clim_slpb_cmip6[10], 
time, tas_clim_slpb_cmip6[11], time, tas_clim_slpb_cmip6[12], time, tas_clim_slpb_cmip6[13], time, tas_clim_slpb_cmip6[14],
time, tas_clim_slpb_cmip6[15], time, tas_clim_slpb_cmip6[16], time, tas_clim_slpb_cmip6[17], time, tas_clim_slpb_cmip6[18], 
time, tas_clim_slpb_cmip6[19], time, tas_clim_slpb_cmip6[20], time, tas_clim_slpb_cmip6[21], time, tas_clim_slpb_cmip6[22], 
time, tas_clim_slpb_cmip5, time, tas_clim_slpb_cru)
plt.title(u'(f) SLPB', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('Months', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(10, 36, 4), fontsize=8)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20, l21, l22, l23, l24, l25 = annual_cycle
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
plt.setp(l18)
plt.setp(l19)
plt.setp(l20)
plt.setp(l21)
plt.setp(l22)
plt.setp(l23, color='black', linestyle='--')
plt.setp(l24, color='black', linestyle='-.')
plt.setp(l25, color='black')

legend = [cmip6[1][0],cmip6[2][0],cmip6[3][0],cmip6[4][0],cmip6[5][0],cmip6[6][0],
cmip6[7][0],cmip6[8][0],cmip6[9][0],cmip6[10][0],cmip6[11][0],cmip6[12][0],
cmip6[13][0],cmip6[14][0],cmip6[15][0],cmip6[16][0],cmip6[17][0],cmip6[18][0],
cmip6[19][0],cmip6[20][0],cmip6[21][0],cmip6[22][0],'CMIP6-MME','CMIP5-MME','CRU']
plt.legend(annual_cycle, legend, loc=(1.019, 0.5), fontsize=8)
plt.subplots_adjust(left=0.15, bottom=0.15, right=0.93, top=0.93, wspace=0.20, hspace=0.35)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'pyplt_clim_cmip6_1961-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()






