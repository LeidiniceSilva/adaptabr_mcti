# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot cri of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_statistical_metrics import compute_ivs


def import_obs(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
		
	return lat, lon, obs_data

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return lat, lon, mdl_data
	              

def compute_cri(p1, p2, p3, p4, p5):

	cri = (p1 +p2 +p3 + p4 + p5) / 90

	return cri
	
	     
# Import cmip models and obs database 
var_obs = 'pre'
var_cmip6 = 'pr'
dt = '1980-2014'

namz_obs_ann_latlon, namz_obs_ann_ts = import_obs(var_obs, 'NAMZ', 'ANN', dt)
namz_obs_ann_latlon, samz_obs_ann_ts = import_obs(var_obs, 'SAMZ', 'ANN', dt)
namz_obs_ann_latlon, neb_obs_ann_ts = import_obs(var_obs, 'NEB', 'ANN', dt)
namz_obs_ann_latlon, sam_obs_ann_ts = import_obs(var_obs, 'SAM', 'ANN', dt)
namz_obs_ann_latlon, lpb_obs_ann_ts = import_obs(var_obs, 'LPB', 'ANN', dt)

namz_obs_mon_latlon, namz_obs_mon_ts = import_obs(var_obs, 'NAMZ', 'MON', dt)
namz_obs_mon_latlon, samz_obs_mon_ts = import_obs(var_obs, 'SAMZ', 'MON', dt)
namz_obs_mon_latlon, neb_obs_mon_ts  = import_obs(var_obs, 'NEB', 'MON', dt)
namz_obs_mon_latlon, sam_obs_mon_ts = import_obs(var_obs, 'SAM', 'MON', dt)
namz_obs_mon_latlon, lpb_obs_mon_ts = import_obs(var_obs, 'LPB', 'MON', dt)

cri_namz_cmip6 = []
cri_samz_cmip6 = []
cri_neb_cmip6 = []
cri_sam_cmip6 = []
cri_lpb_cmip6 = []
legend = []

for i in range(1, 19):

	namz_cmip_ann_latlon, namz_cmip_ann_ts = import_cmip(var_cmip6, 'NAMZ', 'ANN', cmip6[i][0], cmip6[i][1], dt)
	namz_cmip_mon_latlon, namz_cmip_mon_ts = import_cmip(var_cmip6, 'NAMZ', 'MON', cmip6[i][0], cmip6[i][1], dt)
	mbe_namz_cmip = compute_mbe(namz_cmip_ann_ts[2], namz_obs_ann_ts[2])
	rmse_namz_cmip = compute_rsme(namz_cmip_ann_ts[2], namz_obs_ann_ts[2])
	tcc_namz_cmip = compute_tcc(namz_obs_ann_latlon[2], namz_cmip_ann_latlon[2])
	pcc_namz_cmip = compute_pcc(namz_obs_mon_ts[2], namz_cmip_mon_ts[2])
	ivs_namz_cmip = compute_ivs(namz_obs_ann_ts[2], namz_cmip_ann_ts[2])
	cri_namz_cmip = compute_cri(mbe_namz_cmip, rmse_namz_cmip, tcc_namz_cmip, pcc_namz_cmip, ivs_namz_cmip)
	cri_namz_cmip6.append(cri_namz_cmip)

	samz_cmip_ann_latlon, samz_cmip_ann_ts = import_cmip(var_cmip6, 'SAMZ', 'ANN', cmip6[i][0], cmip6[i][1], dt)
	samz_cmip_mon_latlon, samz_cmip_mon_ts = import_cmip(var_cmip6, 'SAMZ', 'MON', cmip6[i][0], cmip6[i][1], dt)
	mbe_samz_cmip = compute_mbe(samz_cmip_ann_ts[2], samz_obs_ann_ts[2])
	rmse_samz_cmip = compute_rsme(samz_cmip_ann_ts[2], samz_obs_ann_ts[2])
	tcc_samz_cmip = compute_tcc(samz_obs_ann_latlon[2], samz_cmip_ann_latlon[2])
	pcc_samz_cmip = compute_pcc(samz_obs_mon_ts[2], samz_cmip_mon_ts[2])
	ivs_samz_cmip = compute_ivs(samz_obs_ann_ts[2], samz_cmip_ann_ts[2])
	cri_samz_cmip = compute_cri(mbe_samz_cmip, rmse_samz_cmip, tcc_samz_cmip, pcc_samz_cmip, ivs_samz_cmip)
	cri_samz_cmip6.append(cri_samz_cmip)
	
	neb_cmip_ann_latlon, neb_cmip_ann_ts = import_cmip(var_cmip6, 'NEB', 'ANN', cmip6[i][0], cmip6[i][1], dt)
	neb_cmip_mon_latlon, neb_cmip_mon_ts = import_cmip(var_cmip6, 'NEB', 'MON', cmip6[i][0], cmip6[i][1], dt)
	mbe_neb_cmip = compute_mbe(neb_cmip_ann_ts[2], neb_obs_ann_ts[2])
	rmse_neb_cmip = compute_rsme(neb[2], neb_obs_ann_ts[2])
	tcc_neb_cmip = compute_tcc(neb_obs_ann_latlon[2], neb_cmip_ann_latlon[2])
	pcc_neb_cmip = compute_pcc(neb_obs_mon_ts[2], neb_cmip_mon_ts[2])
	ivs_neb_cmip = compute_ivs(neb_obs_ann_ts[2], neb_cmip_ann_ts[2])
	cri_neb_cmip = compute_cri(mbe_neb_cmip, rmse_neb_cmip, tcc_neb_cmip, pcc_neb_cmip, ivs_neb_cmip)
	cri_neb_cmip6.append(cri_neb_cmip)
	
	sam_cmip_ann_latlon, sam_cmip_ann_ts = import_cmip(var_cmip6, 'SAM', 'ANN', cmip6[i][0], cmip6[i][1], dt)
	sam_cmip_mon_latlon, sam_cmip_mon_ts = import_cmip(var_cmip6, 'SAM', 'MON', cmip6[i][0], cmip6[i][1], dt)
	mbe_sam_cmip = compute_mbe(sam_cmip_ann_ts[2], sam_obs_ann_ts[2])
	rmse_sam_cmip = compute_rsme(sam_cmip_ann_ts[2], sam_obs_ann_ts[2])
	tcc_sam_cmip = compute_tcc(sam_obs_ann_latlon[2], sam_cmip_ann_latlon[2])
	pcc_sam_cmip = compute_pcc(sam_obs_mon_ts[2], sam_cmip_mon_ts[2])
	ivs_sam_cmip = compute_ivs(sam_obs_ann_ts[2], sam_cmip_ann_ts[2])
	cri_sam_cmip = compute_cri(mbe_sam_cmip, rmse_sam_cmip, tcc_sam_cmip, pcc_sam_cmip, ivs_sam_cmip)
	cri_sam_cmip6.append(cri_sam_cmip)
	
	lpb_cmip_ann_latlon, lpb_cmip_ann_ts = import_cmip(var_cmip6, 'LPB', 'ANN', cmip6[i][0], cmip6[i][1], dt)
	lpb_cmip_mon_latlon, lpb_cmip_mon_ts = import_cmip(var_cmip6, 'LPB', 'MON', cmip6[i][0], cmip6[i][1], dt)
	mbe_lpb_cmip = compute_mbe(lpb_cmip_ann_ts[2], lpb_obs_ann_ts[2])
	rmse_lpb_cmip = compute_rsme(lpb_cmip_ann_ts[2], lpb_obs_ann_ts[2])
	tcc_lpb_cmip = compute_tcc(lpb_obs_ann_latlon[2], lpb_cmip_ann_latlon[2])
	pcc_lpb_cmip = compute_pcc(lpb_obs_mon_ts[2], lpb_cmip_mon_ts[2])
	ivs_lpb_cmip = compute_ivs(lpb_obs_ann_ts[2], lpb_cmip_ann_ts[2])
	cri_lpb_cmip = compute_cri(mbe_lpb_cmip, rmse_lpb_cmip, tcc_lpb_cmip, pcc_lpb_cmip, ivs_lpb_cmip)
	cri_lpb_cmip6.append(cri_lpb_cmip)

	legend.append(cmip6[i][0])

cri_cmip6 = np.array([cri_lpb_cmip6,cri_sam_cmip6,cri_neb_cmip6,cri_samz_cmip6,cri_namz_cmip6])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 3))
norm = colors.BoundaryNorm(boundaries=np.arange(0, 5.5, 0.5), ncolors=256)
xlabels = legend
ylabels = [u'LPB', u'SAM', u'NEB', u'SAMZ', u'NAMZ']

if var_cmip6 == 'pr':
	color = cm.Blues
else:
	color = cm.Reds

ax = fig.add_subplot(1, 1, 1)  
pcm = ax.pcolormesh(cri_cmip6, edgecolors ='white', linewidths = 2., norm=norm, cmap=color)
ax.set_title(u'(a) RCI', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(cri_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(cri_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
for y in range(cri_cmip6.shape[0]):
    for x in range(cri_cmip6.shape[1]):
        ax.text(x + 0.5, y + 0.5, '%.2f' % cri_cmip6[y, x],
                 ha="center", va="center", color='k', size=8)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_portrait_cri_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






