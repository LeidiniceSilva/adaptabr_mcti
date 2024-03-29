# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot ivs of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_stats_metrics import compute_ivs


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

ivs_namz_cmip6 = []
ivs_samz_cmip6 = []
ivs_neb_cmip6 = []
ivs_sam_cmip6 = []
ivs_lpb_cmip6 = []
ivs_br_cmip6 = []
legend = []

for i in range(1, 18):

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	ivs_namz_cmip6.append(compute_ivs(clim_namz_obs, clim_namz_cmip))
	
	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	ivs_samz_cmip6.append(compute_ivs(clim_samz_obs, clim_samz_cmip))
	
	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	ivs_neb_cmip6.append(compute_ivs(clim_neb_obs, clim_neb_cmip))
		
	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	ivs_sam_cmip6.append(compute_ivs(clim_sam_obs, clim_sam_cmip))
		
	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	ivs_lpb_cmip6.append(compute_ivs(clim_lpb_obs, clim_lpb_cmip))

	clim_br_cmip = import_cmip(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], dt)
	ivs_br_cmip6.append(compute_ivs(clim_br_obs, clim_br_cmip))
	
	legend.append(cmip6[i][0])

ivs_cmip6 = np.array([ivs_br_cmip6,ivs_lpb_cmip6,ivs_sam_cmip6,ivs_neb_cmip6,ivs_samz_cmip6,ivs_namz_cmip6])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 3))
norm = colors.BoundaryNorm(boundaries=np.arange(0, 5.5, 0.5), ncolors=256)
xlabels = legend
ylabels = ['BR', 'LPB', 'SAM', 'NEB', 'SAMZ', 'NAMZ']

if var_cmip6 == 'pr':
	color = cm.Blues
else:
	color = cm.Reds

ax = fig.add_subplot(1, 1, 1)  
pcm = ax.pcolormesh(ivs_cmip6, edgecolors ='white', linewidths = 2., norm=norm, cmap=color)
ax.set_title(u'(a) IVS', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
for y in range(ivs_cmip6.shape[0]):
    for x in range(ivs_cmip6.shape[1]):
        ax.text(x + 0.5, y + 0.5, '%.2f' % ivs_cmip6[y, x],
                 ha="center", va="center", color='k', size=8)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_portrait_ivs_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






