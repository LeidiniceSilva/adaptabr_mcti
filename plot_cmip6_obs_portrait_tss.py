# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot tss of cmip6 models"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_statistical_metrics import compute_tss


def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_ANN_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]	
	obs = np.nanmean(value, axis=0)
	
	latlon_obs = []
	for i in range(0, obs.shape[0]):
		for ii in obs[i]:
			latlon_obs.append(ii)
	latlon_obs = np.array(latlon_obs)
			
	return lat, lon, latlon_obs

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_ANN_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl = np.nanmean(value, axis=0)

	latlon_mdl = []
	for i in range(0, mdl.shape[0]):
		for ii in mdl[i]:
			latlon_mdl.append(ii)
	latlon_mdl = np.array(latlon_mdl)
	
	return lat, lon, latlon_mdl
	              
               
# Import cmip models and obs database 
var_obs = 'tmp'
var_cmip6 = 'tas'
dt = '1980-2014'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)

tss_namz_cmip6 = []
tss_samz_cmip6 = []
tss_neb_cmip6 = []
tss_sam_cmip6 = []
tss_lpb_cmip6 = []
legend = []

for i in range(1, 19):

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	tss_namz_cmip6.append(compute_tss(clim_namz_obs[2], clim_namz_cmip[2]))
	
	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	tss_samz_cmip6.append(compute_tss(clim_samz_obs[2], clim_samz_cmip[2]))
	
	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	tss_neb_cmip6.append(compute_tss(clim_neb_obs[2], clim_neb_cmip[2]))
		
	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	tss_sam_cmip6.append(compute_tss(clim_sam_obs[2], clim_sam_cmip[2]))
		
	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	tss_lpb_cmip6.append(compute_tss(clim_lpb_obs[2], clim_lpb_cmip[2]))
	
	legend.append(cmip6[i][0])

ivs_cmip6 = np.array([tss_lpb_cmip6,tss_sam_cmip6,tss_neb_cmip6,tss_samz_cmip6,tss_namz_cmip6])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 3))
norm = colors.BoundaryNorm(boundaries=np.arange(0, 1., 0.1), ncolors=256)
xlabels = legend
ylabels = [u'LPB', u'SAM', u'NEB', u'SAMZ', u'NAMZ']

if var_cmip6 == 'pr':
	color = cm.Blues
else:
	color = cm.Reds

ax = fig.add_subplot(1, 1, 1)  
pcm = ax.pcolormesh(ivs_cmip6, edgecolors ='white', linewidths = 2., norm=norm, cmap=color)
ax.set_title(u'(a) TSS', loc='left', fontweight='bold', fontsize=8)
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
name_out = 'pyplt_portrait_tss_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






