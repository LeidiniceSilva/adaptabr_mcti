# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot all statistical of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_stats_metrics import compute_mbe
from comp_stats_metrics import compute_rmse
from comp_stats_metrics import compute_tss
from comp_stats_metrics import compute_pcc
from comp_stats_metrics import compute_ivs


def import_obs_latlon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
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
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
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
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
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
	
		
# Import cmip models and obs database 
var_obs = 'Tmin'
var_cmip6 = 'tasmin'
dt = '1986-2005'

namz_obs_latlon = import_obs_latlon(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_latlon = import_obs_latlon(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_latlon = import_obs_latlon(var_obs, 'NEB', 'ANN', dt)
sam_obs_latlon = import_obs_latlon(var_obs, 'SAM', 'ANN', dt)
lpb_obs_latlon = import_obs_latlon(var_obs, 'LPB', 'ANN', dt)
br_obs_latlon = import_obs_latlon(var_obs, 'BR', 'ANN', dt)

namz_obs_mon_ts = import_obs_mon(var_obs, 'NAMZ', 'MON', dt)
samz_obs_mon_ts = import_obs_mon(var_obs, 'SAMZ', 'MON', dt)
neb_obs_mon_ts  = import_obs_mon(var_obs, 'NEB', 'MON', dt)
sam_obs_mon_ts = import_obs_mon(var_obs, 'SAM', 'MON', dt)
lpb_obs_mon_ts = import_obs_mon(var_obs, 'LPB', 'MON', dt)
br_obs_mon_ts = import_obs_mon(var_obs, 'BR', 'MON', dt)

namz_obs_ann_ts = import_obs_ann(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_ann_ts = import_obs_ann(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_ann_ts  = import_obs_ann(var_obs, 'NEB', 'ANN', dt)
sam_obs_ann_ts = import_obs_ann(var_obs, 'SAM', 'ANN', dt)
lpb_obs_ann_ts = import_obs_ann(var_obs, 'LPB', 'ANN', dt)
br_obs_ann_ts = import_obs_ann(var_obs, 'BR', 'ANN', dt)

mbe_namz_cmip6 = []
rmse_namz_cmip6 = []
tss_namz_cmip6 = []
pcc_namz_cmip6 = []
ivs_namz_cmip6 = []

mbe_samz_cmip6 = []
rmse_samz_cmip6 = []
tss_samz_cmip6 = []
pcc_samz_cmip6 = []
ivs_samz_cmip6 = []

mbe_neb_cmip6 = []
rmse_neb_cmip6 = []
tss_neb_cmip6 = []
pcc_neb_cmip6 = []
ivs_neb_cmip6 = []

mbe_sam_cmip6 = []
rmse_sam_cmip6 = []
tss_sam_cmip6 = []
pcc_sam_cmip6 = []
ivs_sam_cmip6 = []

mbe_lpb_cmip6 = []
rmse_lpb_cmip6 = []
tss_lpb_cmip6 = []
pcc_lpb_cmip6 = []
ivs_lpb_cmip6 = []

mbe_br_cmip6 = []
rmse_br_cmip6 = []
tss_br_cmip6 = []
pcc_br_cmip6 = []
ivs_br_cmip6 = []

legend = []

for i in range(1, 18):

	namz_cmip_latlon = import_cmip_latlon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_latlon = import_cmip_latlon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_latlon = import_cmip_latlon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	br_cmip_latlon = import_cmip_latlon(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	
	namz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	samz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	neb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	sam_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'MON', dt)
	lpb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	br_cmip_mon_ts = import_cmip_mon(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'MON', dt)
	
	namz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	br_cmip_ann_ts = import_cmip_ann(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'ANN', dt)
			
	mbe_namz_cmip6.append(compute_mbe(namz_cmip_latlon, namz_obs_latlon))
	rmse_namz_cmip6.append(compute_rmse(namz_cmip_latlon, namz_obs_latlon))	
	tss_namz_cmip6.append(compute_tss(namz_obs_latlon, namz_cmip_latlon))
	pcc_namz_cmip6.append(compute_pcc(namz_obs_mon_ts, namz_cmip_mon_ts))
	ivs_namz_cmip6.append(compute_ivs(namz_obs_ann_ts, namz_cmip_ann_ts))
	
	mbe_samz_cmip6.append(compute_mbe(samz_cmip_latlon, samz_obs_latlon))
	rmse_samz_cmip6.append(compute_rmse(samz_cmip_latlon, samz_obs_latlon))
	tss_samz_cmip6.append(compute_tss(samz_obs_latlon, samz_cmip_latlon))
	pcc_samz_cmip6.append(compute_pcc(samz_obs_mon_ts, samz_cmip_mon_ts))
	ivs_samz_cmip6.append(compute_ivs(samz_obs_ann_ts, samz_cmip_ann_ts))
		
	mbe_neb_cmip6.append(compute_mbe(neb_cmip_latlon, neb_obs_latlon))
	rmse_neb_cmip6.append(compute_rmse(neb_cmip_latlon, neb_obs_latlon))
	tss_neb_cmip6.append(compute_tss(neb_obs_latlon, neb_cmip_latlon))
	pcc_neb_cmip6.append(compute_pcc(neb_obs_mon_ts, neb_cmip_mon_ts))
	ivs_neb_cmip6.append(compute_ivs(neb_obs_ann_ts, neb_cmip_ann_ts))
		
	mbe_sam_cmip6.append(compute_mbe(sam_cmip_latlon, sam_obs_latlon))
	rmse_sam_cmip6.append(compute_rmse(sam_cmip_latlon, sam_obs_latlon))
	tss_sam_cmip6.append(compute_tss(sam_obs_latlon, sam_cmip_latlon))
	pcc_sam_cmip6.append(compute_pcc(sam_obs_mon_ts, sam_cmip_mon_ts))
	ivs_sam_cmip6.append(compute_ivs(sam_obs_ann_ts, sam_cmip_ann_ts))
		
	mbe_lpb_cmip6.append(compute_mbe(lpb_cmip_latlon, lpb_obs_latlon))
	rmse_lpb_cmip6.append(compute_rmse(lpb_cmip_latlon, lpb_obs_latlon))
	tss_lpb_cmip6.append(compute_tss(lpb_obs_latlon, lpb_cmip_latlon))
	pcc_lpb_cmip6.append(compute_pcc(lpb_obs_mon_ts, lpb_cmip_mon_ts))
	ivs_lpb_cmip6.append(compute_ivs(lpb_obs_ann_ts, lpb_cmip_ann_ts))
	
	mbe_br_cmip6.append(compute_mbe(br_cmip_latlon, br_obs_latlon))
	rmse_br_cmip6.append(compute_rmse(br_cmip_latlon, br_obs_latlon))
	tss_br_cmip6.append(compute_tss(br_obs_latlon, br_cmip_latlon))
	pcc_br_cmip6.append(compute_pcc(br_obs_mon_ts, br_cmip_mon_ts))
	ivs_br_cmip6.append(compute_ivs(br_obs_ann_ts, br_cmip_ann_ts))

	legend.append(cmip6[i][0])

mbe_cmip6 = np.array([mbe_br_cmip6,mbe_lpb_cmip6,mbe_sam_cmip6,mbe_neb_cmip6,mbe_samz_cmip6,mbe_namz_cmip6])
rmse_cmip6 = np.array([rmse_br_cmip6,rmse_lpb_cmip6,rmse_sam_cmip6,rmse_neb_cmip6,rmse_samz_cmip6,rmse_namz_cmip6])
tss_cmip6 = np.array([tss_br_cmip6,tss_lpb_cmip6,tss_sam_cmip6,tss_neb_cmip6,tss_samz_cmip6,tss_namz_cmip6])
pcc_cmip6 = np.array([pcc_br_cmip6,pcc_lpb_cmip6,pcc_sam_cmip6,pcc_neb_cmip6,pcc_samz_cmip6,pcc_namz_cmip6])
ivs_cmip6 = np.array([ivs_br_cmip6,ivs_lpb_cmip6,ivs_sam_cmip6,ivs_neb_cmip6,ivs_samz_cmip6,ivs_namz_cmip6])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 6))

xlabels = legend
ylabels = ['BR', 'LPB', 'SAM', 'NEB', 'SAMZ', 'NAMZ']

norm1 = colors.BoundaryNorm(boundaries=np.arange(-6, 7, 1), ncolors=256)
norm2 = colors.BoundaryNorm(boundaries=np.arange(0, 5.5, .5), ncolors=256)
norm3 = colors.BoundaryNorm(boundaries=np.arange(0, 1., 0.1), ncolors=256)
norm4 = colors.BoundaryNorm(boundaries=np.arange(-1, 1.1, 0.1), ncolors=256)

if var_cmip6 == 'pr':
	color1 = cm.BrBG
	color2 = cm.Blues
	color3 = cm.PiYG
else:
	color1 = cm.bwr
	color2 = cm.Reds
	color3 = cm.PiYG

ax = fig.add_subplot(3, 2, 1)  
pcm = ax.pcolormesh(mbe_cmip6, edgecolors='white', linewidths=2., norm=norm1, cmap=color1)
ax.set_title(u'(a) MBE', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(mbe_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(mbe_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(3, 2, 2)  
pcm = ax.pcolormesh(rmse_cmip6, edgecolors='white', linewidths=2., norm=norm2, cmap=color2)
ax.set_title(u'(b) RMSE', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rmse_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(rmse_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)                

ax = fig.add_subplot(3, 2, 3)  
pcm = ax.pcolormesh(tss_cmip6, edgecolors='white', linewidths=2., norm=norm3, cmap=color2)
ax.set_title(u'(c) TSS', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(tss_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(tss_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(3, 2, 4)  
pcm = ax.pcolormesh(pcc_cmip6, edgecolors='white', linewidths=2., norm=norm4, cmap=color3)
ax.set_title(u'(d) PCC', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(pcc_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(pcc_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)     

ax = fig.add_subplot(3, 2, 5)  
pcm = ax.pcolormesh(ivs_cmip6, edgecolors='white', linewidths=2., norm=norm2, cmap=color2)
ax.set_title(u'(e) IVS', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)   
     	                
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_portrait_indices_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

