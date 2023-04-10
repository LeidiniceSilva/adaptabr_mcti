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
from comp_statistical_metrics import compute_mbe
from comp_statistical_metrics import compute_rmse
from comp_statistical_metrics import compute_tss
from comp_statistical_metrics import compute_pcc
from comp_statistical_metrics import compute_ivs


def import_obs_latlon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
	
	return ts_latlon
	
	
def import_obs_mon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mon_mean = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	ts_mon = []
	for i in range(0, 11 + 1):
		clim = np.nanmean(mon_mean[i::12], axis=0)
		ts_mon.append(clim)
	
	return ts_mon


def import_obs_ann(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	ts_ann = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return ts_ann
	
	
def import_cmip_latlon(param, area, model, exp, period, date):

	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
	
	return ts_latlon
	              

def import_cmip_mon(param, area, model, exp, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mon_mean = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	ts_mon = []
	for i in range(0, 11 + 1):
		clim = np.nanmean(mon_mean[i::12], axis=0)
		ts_mon.append(clim)
	
	return ts_mon
	

def import_cmip_ann(param, area, model, exp, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	ts_ann = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return ts_ann
	
	
def compute_cri(p1, p2, p3, p4, p5):

	cri = (p1 + p2 + p3 + p4 + p5) / 90

	return cri
	
	     
# Import cmip models and obs database 
var_obs = 'tmp'
var_cmip6 = 'tas'
dt = '1980-2014'

namz_obs_latlon = import_obs_latlon(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_latlon = import_obs_latlon(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_latlon = import_obs_latlon(var_obs, 'NEB', 'ANN', dt)
sam_obs_latlon = import_obs_latlon(var_obs, 'SAM', 'ANN', dt)
lpb_obs_latlon = import_obs_latlon(var_obs, 'LPB', 'ANN', dt)

namz_obs_mon_ts = import_obs_mon(var_obs, 'NAMZ', 'MON', dt)
samz_obs_mon_ts = import_obs_mon(var_obs, 'SAMZ', 'MON', dt)
neb_obs_mon_ts  = import_obs_mon(var_obs, 'NEB', 'MON', dt)
sam_obs_mon_ts = import_obs_mon(var_obs, 'SAM', 'MON', dt)
lpb_obs_mon_ts = import_obs_mon(var_obs, 'LPB', 'MON', dt)

namz_obs_ann_ts = import_obs_ann(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_ann_ts = import_obs_ann(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_ann_ts  = import_obs_ann(var_obs, 'NEB', 'ANN', dt)
sam_obs_ann_ts = import_obs_ann(var_obs, 'SAM', 'ANN', dt)
lpb_obs_ann_ts = import_obs_ann(var_obs, 'LPB', 'ANN', dt)

cri_namz_cmip6 = []
cri_samz_cmip6 = []
cri_neb_cmip6 = []
cri_sam_cmip6 = []
cri_lpb_cmip6 = []
legend = []

for i in range(1, 19):

	namz_cmip_latlon = import_cmip_latlon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_latlon = import_cmip_latlon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_latlon = import_cmip_latlon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	
	namz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	samz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	neb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	sam_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'MON', dt)
	lpb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	
	namz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
			
	mbe_namz_cmip = compute_mbe(namz_cmip_latlon, namz_obs_latlon)
	rmse_namz_cmip = compute_rmse(namz_cmip_latlon, namz_obs_latlon)	
	tcc_namz_cmip = compute_tss(namz_obs_latlon, namz_cmip_latlon)
	pcc_namz_cmip = compute_pcc(namz_obs_mon_ts, namz_cmip_mon_ts)
	ivs_namz_cmip = compute_ivs(namz_obs_ann_ts, namz_cmip_ann_ts)
	cri_namz_cmip = compute_cri(mbe_namz_cmip, rmse_namz_cmip, tcc_namz_cmip, pcc_namz_cmip, ivs_namz_cmip)
	cri_namz_cmip6.append(cri_namz_cmip)

	mbe_samz_cmip = compute_mbe(samz_cmip_latlon, samz_obs_latlon)
	rmse_samz_cmip = compute_rmse(samz_cmip_latlon, samz_obs_latlon)
	tcc_samz_cmip = compute_tss(samz_obs_latlon, samz_cmip_latlon)
	pcc_samz_cmip = compute_pcc(samz_obs_mon_ts, samz_cmip_mon_ts)
	ivs_samz_cmip = compute_ivs(samz_obs_ann_ts, samz_cmip_ann_ts)
	cri_samz_cmip = compute_cri(mbe_samz_cmip, rmse_samz_cmip, tcc_samz_cmip, pcc_samz_cmip, ivs_samz_cmip)
	cri_samz_cmip6.append(cri_samz_cmip)
	
	mbe_neb_cmip = compute_mbe(neb_cmip_latlon, neb_obs_latlon)
	rmse_neb_cmip = compute_rmse(neb_cmip_latlon, neb_obs_latlon)
	tcc_neb_cmip = compute_tss(neb_obs_latlon, neb_cmip_latlon)
	pcc_neb_cmip = compute_pcc(neb_obs_mon_ts, neb_cmip_mon_ts)
	ivs_neb_cmip = compute_ivs(neb_obs_ann_ts, neb_cmip_ann_ts)
	cri_neb_cmip = compute_cri(mbe_neb_cmip, rmse_neb_cmip, tcc_neb_cmip, pcc_neb_cmip, ivs_neb_cmip)
	cri_neb_cmip6.append(cri_neb_cmip)
	
	mbe_sam_cmip = compute_mbe(sam_cmip_latlon, sam_obs_latlon)
	rmse_sam_cmip = compute_rmse(sam_cmip_latlon, sam_obs_latlon)
	tcc_sam_cmip = compute_tss(sam_obs_latlon, sam_cmip_latlon)
	pcc_sam_cmip = compute_pcc(sam_obs_mon_ts, sam_cmip_mon_ts)
	ivs_sam_cmip = compute_ivs(sam_obs_ann_ts, sam_cmip_ann_ts)
	cri_sam_cmip = compute_cri(mbe_sam_cmip, rmse_sam_cmip, tcc_sam_cmip, pcc_sam_cmip, ivs_sam_cmip)
	cri_sam_cmip6.append(cri_sam_cmip)
	
	mbe_lpb_cmip = compute_mbe(lpb_cmip_latlon, lpb_obs_latlon)
	rmse_lpb_cmip = compute_rmse(lpb_cmip_latlon, lpb_obs_latlon)
	tcc_lpb_cmip = compute_tss(lpb_obs_latlon, lpb_cmip_latlon)
	pcc_lpb_cmip = compute_pcc(lpb_obs_mon_ts, lpb_cmip_mon_ts)
	ivs_lpb_cmip = compute_ivs(lpb_obs_ann_ts, lpb_cmip_ann_ts)
	cri_lpb_cmip = compute_cri(mbe_lpb_cmip, rmse_lpb_cmip, tcc_lpb_cmip, pcc_lpb_cmip, ivs_lpb_cmip)
	cri_lpb_cmip6.append(cri_lpb_cmip)

	legend.append(cmip6[i][0])

cri_cmip6 = np.array([cri_lpb_cmip6,cri_sam_cmip6,cri_neb_cmip6,cri_samz_cmip6,cri_namz_cmip6])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 3))
norm = colors.BoundaryNorm(boundaries=np.arange(0, 0.475, 0.025), ncolors=256)
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






