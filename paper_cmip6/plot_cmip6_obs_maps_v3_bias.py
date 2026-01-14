# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "June 10, 2024"
__description__ = "This script plot portrait of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat

from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from dict_cmip6_models_name import cmip6
from comp_stats_metrics import compute_nrmse, compute_tss, compute_corr, compute_ivs

dt = '197901-201412'
latlon = [-100, -20, -60, 15]
font_size=9

path  = '/home/mda_silv/users/AdaptaBr_MCTI'
		    
			    
def import_obs_srf(param, area):

	arq  = '{0}/database/paper_cmip6/obs/{1}_{2}_ERA5_mon_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	
	time = data.variables['valid_time'][:]
	ntime = len(time)
	nyears = ntime // 12
	
	var_reshaped = var[:nyears*12].reshape(nyears, 12, var.shape[1], len(lat), len(lon))
	annual_mean = np.nanmean(var_reshaped, axis=1)
	mean_850 = np.nanmean(annual_mean[:, 0, :, :], axis=0)
	mean_500 = np.nanmean(annual_mean[:, 1, :, :], axis=0)
	mean_200 = np.nanmean(annual_mean[:, 2, :, :], axis=0)
		
	return lat, lon, mean_850, mean_500, mean_200


def import_cmip_srf(param, area, model, exp):
	
	arq   = '{0}/database/paper_cmip6/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	time = data.variables['time'][:]
	ntime = len(time)
	nyears = ntime // 12
	
	var_reshaped = var[:nyears*12].reshape(nyears, 12, var.shape[1], len(lat), len(lon))
	annual_mean = np.nanmean(var_reshaped, axis=1)
	mean_850 = annual_mean[:, 0, :, :]
	mean_500 = annual_mean[:, 1, :, :]
	mean_200 = annual_mean[:, 2, :, :]
		
	return lat, lon, mean_850, mean_500, mean_200
    

def configure_subplot(ax):

	ax.set_extent(latlon, crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(latlon[0], latlon[1], 20), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(latlon[2], latlon[3], 20), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(font_size)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)
	
		
# Import obs database and cmip models
lat, lon, q_850_obs, q_500_obs, q_200_obs = import_obs_srf('q', 'sa')
lat, lon, u_850_obs, u_500_obs, u_200_obs = import_obs_srf('u', 'sa')
lat, lon, v_850_obs, v_500_obs, v_200_obs = import_obs_srf('v', 'sa')

lat, lon, q_850_mdl1b, q_500_mdl1b, q_200_mdl1b = import_cmip_srf('hus', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, q_850_mdl2b, q_500_mdl2b, q_200_mdl2b = import_cmip_srf('hus', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, q_850_mdl3b, q_500_mdl3b, q_200_mdl3b = import_cmip_srf('hus', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, q_850_mdl1w, q_500_mdl1w, q_200_mdl1w = import_cmip_srf('hus', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, q_850_mdl2w, q_500_mdl2w, q_200_mdl2w = import_cmip_srf('hus', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, q_850_mdl3w, q_500_mdl3w, q_200_mdl3w = import_cmip_srf('hus', 'sa', 'NESM3', 'r1i1p1f1_gn')

lat, lon, u_850_mdl1b, u_500_mdl1b, u_200_mdl1b = import_cmip_srf('ua', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, u_850_mdl2b, u_500_mdl2b, u_200_mdl2b = import_cmip_srf('ua', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, u_850_mdl3b, u_500_mdl3b, u_200_mdl3b = import_cmip_srf('ua', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, u_850_mdl1w, u_500_mdl1w, u_200_mdl1w = import_cmip_srf('ua', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, u_850_mdl2w, u_500_mdl2w, u_200_mdl2w = import_cmip_srf('ua', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, u_850_mdl3w, u_500_mdl3w, u_200_mdl3w = import_cmip_srf('ua', 'sa', 'NESM3', 'r1i1p1f1_gn')

lat, lon, v_850_mdl1b, v_500_mdl1b, v_200_mdl1b = import_cmip_srf('va', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, v_850_mdl2b, v_500_mdl2b, v_200_mdl2b = import_cmip_srf('va', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, v_850_mdl3b, v_500_mdl3b, v_200_mdl3b = import_cmip_srf('va', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, v_850_mdl1w, v_500_mdl1w, v_200_mdl1w = import_cmip_srf('va', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, v_850_mdl2w, v_500_mdl2w, v_200_mdl2w = import_cmip_srf('va', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, v_850_mdl3w, v_500_mdl3w, v_200_mdl3w = import_cmip_srf('va', 'sa', 'NESM3', 'r1i1p1f1_gn')   
    
q_850_mme_best = np.nanmean(np.nanmean([q_850_mdl1b, q_850_mdl2b, q_850_mdl3b], axis=0), axis=0)
q_500_mme_best = np.nanmean(np.nanmean([q_500_mdl1b, q_500_mdl2b, q_500_mdl3b], axis=0), axis=0)
q_200_mme_best = np.nanmean(np.nanmean([q_200_mdl1b, q_200_mdl2b, q_200_mdl3b], axis=0), axis=0)

q_850_mme_worse = np.nanmean(np.nanmean([q_850_mdl1w, q_850_mdl2w, q_850_mdl3w], axis=0), axis=0)
q_500_mme_worse = np.nanmean(np.nanmean([q_500_mdl1w, q_500_mdl2w, q_500_mdl3w], axis=0), axis=0)
q_200_mme_worse = np.nanmean(np.nanmean([q_200_mdl1w, q_200_mdl2w, q_200_mdl3w], axis=0), axis=0)

u_850_mme_best = np.nanmean(np.nanmean([u_850_mdl1b, u_850_mdl2b, u_850_mdl3b], axis=0), axis=0)
u_500_mme_best = np.nanmean(np.nanmean([u_500_mdl1b, u_500_mdl2b, u_500_mdl3b], axis=0), axis=0)
u_200_mme_best = np.nanmean(np.nanmean([u_200_mdl1b, u_200_mdl2b, u_200_mdl3b], axis=0), axis=0)

u_850_mme_worse = np.nanmean(np.nanmean([u_850_mdl1w, u_850_mdl2w, u_850_mdl3w], axis=0), axis=0)
u_500_mme_worse = np.nanmean(np.nanmean([u_500_mdl1w, u_500_mdl2w, u_500_mdl3w], axis=0), axis=0)
u_200_mme_worse = np.nanmean(np.nanmean([u_200_mdl1w, u_200_mdl2w, u_200_mdl3w], axis=0), axis=0)

v_850_mme_best = np.nanmean(np.nanmean([v_850_mdl1b, v_850_mdl2b, v_850_mdl3b], axis=0), axis=0)
v_500_mme_best = np.nanmean(np.nanmean([v_500_mdl1b, v_500_mdl2b, v_500_mdl3b], axis=0), axis=0)
v_200_mme_best = np.nanmean(np.nanmean([v_200_mdl1b, v_200_mdl2b, v_200_mdl3b], axis=0), axis=0)

v_850_mme_worse = np.nanmean(np.nanmean([v_850_mdl1w, v_850_mdl2w, v_850_mdl3w], axis=0), axis=0)
v_500_mme_worse = np.nanmean(np.nanmean([v_500_mdl1w, v_500_mdl2w, v_500_mdl3w], axis=0), axis=0)
v_200_mme_worse = np.nanmean(np.nanmean([v_200_mdl1w, v_200_mdl2w, v_200_mdl3w], axis=0), axis=0)

q_850_clim_, q_500_clim_, q_200_clim_, u_850_clim_, u_500_clim_, u_200_clim_, v_850_clim_, v_500_clim_, v_200_clim_ = [], [], [], [], [], [], [], [], []
for i in range(1, 18):
	lat, lon, q_850_clim, q_500_clim, q_200_clim = import_cmip_srf('hus', 'sa', cmip6[i][0], cmip6[i][1])
	lat, lon, u_850_clim, u_500_clim, u_200_clim = import_cmip_srf('ua', 'sa', cmip6[i][0], cmip6[i][1])
	lat, lon, v_850_clim, v_500_clim, v_200_clim = import_cmip_srf('va', 'sa', cmip6[i][0], cmip6[i][1])

	q_850_clim_.append(q_850_clim)
	q_500_clim_.append(q_500_clim)
	q_200_clim_.append(q_200_clim)

	u_850_clim_.append(u_850_clim)
	u_500_clim_.append(u_500_clim)
	u_200_clim_.append(u_200_clim)
	
	v_850_clim_.append(v_850_clim)
	v_500_clim_.append(v_500_clim)
	v_200_clim_.append(v_200_clim)
	
q_850_mme = np.nanmean(np.nanmean([q_850_clim_[0], q_850_clim_[1], q_850_clim_[2], q_850_clim_[3], q_850_clim_[4],
q_850_clim_[5], q_850_clim_[6], q_850_clim_[7], q_850_clim_[8], q_850_clim_[9], q_850_clim_[10], q_850_clim_[11],
q_850_clim_[12], q_850_clim_[13], q_850_clim_[14], q_850_clim_[15], q_850_clim_[16]], axis=0), axis=0)

q_500_mme = np.nanmean(np.nanmean([q_500_clim_[0], q_500_clim_[1], q_500_clim_[2], q_500_clim_[3], q_500_clim_[4],
q_500_clim_[5], q_500_clim_[6], q_500_clim_[7], q_500_clim_[8], q_500_clim_[9], q_500_clim_[10], q_500_clim_[11],
q_500_clim_[12], q_500_clim_[13], q_500_clim_[14], q_500_clim_[15], q_500_clim_[16]], axis=0), axis=0)

q_200_mme = np.nanmean(np.nanmean([q_200_clim_[0], q_200_clim_[1], q_200_clim_[2], q_200_clim_[3], q_200_clim_[4],
q_200_clim_[5], q_200_clim_[6], q_200_clim_[7], q_200_clim_[8], q_200_clim_[9], q_200_clim_[10], q_200_clim_[11],
q_200_clim_[12], q_200_clim_[13], q_200_clim_[14], q_200_clim_[15], q_200_clim_[16]], axis=0), axis=0)

u_850_mme = np.nanmean(np.nanmean([u_850_clim_[0], u_850_clim_[1], u_850_clim_[2], u_850_clim_[3], u_850_clim_[4],
u_850_clim_[5], u_850_clim_[6], u_850_clim_[7], u_850_clim_[8], u_850_clim_[9], u_850_clim_[10], u_850_clim_[11],
u_850_clim_[12], u_850_clim_[13], u_850_clim_[14], u_850_clim_[15], u_850_clim_[16]], axis=0), axis=0)

u_500_mme = np.nanmean(np.nanmean([u_500_clim_[0], u_500_clim_[1], u_500_clim_[2], u_500_clim_[3], u_500_clim_[4],
u_500_clim_[5], u_500_clim_[6], u_500_clim_[7], u_500_clim_[8], u_500_clim_[9], u_500_clim_[10], u_500_clim_[11],
u_500_clim_[12], u_500_clim_[13], u_500_clim_[14], u_500_clim_[15], u_500_clim_[16]], axis=0), axis=0)

u_200_mme = np.nanmean(np.nanmean([u_200_clim_[0], u_200_clim_[1], u_200_clim_[2], u_200_clim_[3], u_200_clim_[4],
u_200_clim_[5], u_200_clim_[6], u_200_clim_[7], u_200_clim_[8], u_200_clim_[9], u_200_clim_[10], u_200_clim_[11],
u_200_clim_[12], u_200_clim_[13],u_200_clim_[14], u_200_clim_[15], u_200_clim_[16]], axis=0), axis=0)

v_850_mme = np.nanmean(np.nanmean([v_850_clim_[0], v_850_clim_[1], v_850_clim_[2], v_850_clim_[3], v_850_clim_[4],
v_850_clim_[5], v_850_clim_[6], v_850_clim_[7], v_850_clim_[8], v_850_clim_[9], v_850_clim_[10], v_850_clim_[11],
v_850_clim_[12], v_850_clim_[13], v_850_clim_[14], v_850_clim_[15], v_850_clim_[16]], axis=0), axis=0)

v_500_mme = np.nanmean(np.nanmean([v_500_clim_[0], v_500_clim_[1], v_500_clim_[2], v_500_clim_[3], v_500_clim_[4],
v_500_clim_[5], v_500_clim_[6], v_500_clim_[7], v_500_clim_[8], v_500_clim_[9], v_500_clim_[10], v_500_clim_[11],
v_500_clim_[12], v_500_clim_[13], v_500_clim_[14], v_500_clim_[15], v_500_clim_[16]], axis=0), axis=0)

v_200_mme = np.nanmean(np.nanmean([v_200_clim_[0], v_200_clim_[1], v_200_clim_[2], v_200_clim_[3], v_200_clim_[4],
v_200_clim_[5], v_200_clim_[6], v_200_clim_[7], v_200_clim_[8], v_200_clim_[9], v_200_clim_[10], v_200_clim_[11],
v_200_clim_[12], v_200_clim_[13], v_200_clim_[14], v_200_clim_[15], v_200_clim_[16]], axis=0), axis=0)

uv_850_obs = np.sqrt(u_850_obs**2+v_850_obs**2)
uv_500_obs = np.sqrt(u_500_obs**2+v_500_obs**2)
uv_200_obs = np.sqrt(u_200_obs**2+v_200_obs**2)

uv_850_mme = np.sqrt(u_850_mme**2+v_850_mme**2)
uv_500_mme = np.sqrt(u_500_mme**2+v_500_mme**2)
uv_200_mme = np.sqrt(u_200_mme**2+v_200_mme**2)

uv_850_mme_best = np.sqrt(u_850_mme_best**2+v_850_mme**2)
uv_500_mme_best = np.sqrt(u_500_mme_best**2+v_500_mme**2)
uv_200_mme_best = np.sqrt(u_200_mme_best**2+v_200_mme**2)

uv_850_mme_worse = np.sqrt(u_850_mme_worse**2+v_850_mme**2)
uv_500_mme_worse = np.sqrt(u_500_mme_worse**2+v_500_mme**2)
uv_200_mme_worse = np.sqrt(u_200_mme_worse**2+v_200_mme**2)

# Bias
bias_mme_obs_q_850 = q_850_mme - q_850_obs
bias_mme_obs_q_500 = q_500_mme - q_500_obs
bias_mme_obs_q_200 = q_200_mme - q_200_obs

bias_mme_best_obs_q_850 = q_850_mme_best - q_850_obs
bias_mme_best_obs_q_500 = q_500_mme_best - q_500_obs
bias_mme_best_obs_q_200 = q_200_mme_best - q_200_obs

bias_mme_worse_obs_q_850 = q_850_mme_worse - q_850_obs
bias_mme_worse_obs_q_500 = q_500_mme_worse - q_500_obs
bias_mme_worse_obs_q_200 = q_200_mme_worse - q_200_obs

bias_mme_obs_uv_850 = uv_850_mme - uv_850_obs
bias_mme_obs_uv_500 = uv_500_mme - uv_500_obs
bias_mme_obs_uv_200 = uv_200_mme - uv_200_obs

bias_mme_best_obs_uv_850 = uv_850_mme_best - uv_850_obs
bias_mme_best_obs_uv_500 = uv_500_mme_best - uv_500_obs
bias_mme_best_obs_uv_200 = uv_200_mme_best - uv_200_obs

bias_mme_worse_obs_uv_850 = uv_850_mme_worse - uv_850_obs
bias_mme_worse_obs_uv_500 = uv_500_mme_worse - uv_500_obs
bias_mme_worse_obs_uv_200 = uv_200_mme_worse - uv_200_obs

# Corr
corr_mme_obs_q_850 = np.corrcoef(q_850_mme.flatten(), q_850_obs.flatten())[0, 1]
corr_mme_obs_q_500 = np.corrcoef(q_500_mme.flatten(), q_500_obs.flatten())[0, 1]
corr_mme_obs_q_200 = np.corrcoef(q_200_mme.flatten(), q_200_obs.flatten())[0, 1]

corr_mme_best_obs_q_850 = np.corrcoef(q_850_mme_best.flatten(), q_850_obs.flatten())[0, 1]
corr_mme_best_obs_q_500 = np.corrcoef(q_500_mme_best.flatten(), q_500_obs.flatten())[0, 1]
corr_mme_best_obs_q_200 = np.corrcoef(q_200_mme_best.flatten(), q_200_obs.flatten())[0, 1]

corr_mme_worse_obs_q_850 = np.corrcoef(q_850_mme_worse.flatten(), q_850_obs.flatten())[0, 1]
corr_mme_worse_obs_q_500 = np.corrcoef(q_500_mme_worse.flatten(), q_500_obs.flatten())[0, 1]
corr_mme_worse_obs_q_200 = np.corrcoef(q_200_mme_worse.flatten(), q_200_obs.flatten())[0, 1]

corr_mme_obs_uv_850 = np.corrcoef(uv_850_mme.flatten(), uv_850_obs.flatten())[0, 1]
corr_mme_obs_uv_500 = np.corrcoef(uv_500_mme.flatten(), uv_500_obs.flatten())[0, 1]
corr_mme_obs_uv_200 = np.corrcoef(uv_200_mme.flatten(), uv_200_obs.flatten())[0, 1]

corr_mme_best_obs_uv_850 = np.corrcoef(uv_850_mme_best.flatten(), uv_850_obs.flatten())[0, 1]
corr_mme_best_obs_uv_500 = np.corrcoef(uv_500_mme_best.flatten(), uv_500_obs.flatten())[0, 1]
corr_mme_best_obs_uv_200 = np.corrcoef(uv_200_mme_best.flatten(), uv_200_obs.flatten())[0, 1]

corr_mme_worse_obs_uv_850 = np.corrcoef(uv_850_mme_worse.flatten(), uv_850_obs.flatten())[0, 1]
corr_mme_worse_obs_uv_500 = np.corrcoef(uv_500_mme_worse.flatten(), uv_500_obs.flatten())[0, 1]
corr_mme_worse_obs_uv_200 = np.corrcoef(uv_200_mme_worse.flatten(), uv_200_obs.flatten())[0, 1]

# Plot figure
fig, axes = plt.subplots(3, 3, figsize=(12, 13), subplot_kw={'projection': ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
cmap_color=cm.rainbow_r

cf1 = ax1.contourf(lon, lat, bias_mme_obs_q_200, levels=np.arange(-0.02, 0.021, 0.001), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct1 = ax1.contour(lon, lat, bias_mme_obs_uv_200, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax1.clabel(ct1, fontsize=font_size, colors='black')
ax1.set_title('(a) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('200 hPa', fontsize=font_size, fontweight='bold')
ax1.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_obs_q_200, 2), round(corr_mme_obs_uv_200, 2)), transform=ax1.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)
cbar = plt.colorbar(cf1, cax=fig.add_axes([0.25, 0.634, 0.5, 0.015]), orientation='horizontal')
cbar.set_label('Specific humidity 200 hPa (g kg$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf2 = ax2.contourf(lon, lat, bias_mme_best_obs_q_200, levels=np.arange(-0.02, 0.021, 0.001), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct2 = ax2.contour(lon, lat, bias_mme_best_obs_uv_200, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax2.clabel(ct2, fontsize=font_size, colors='black')
ax2.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_best_obs_q_200, 2), round(corr_mme_best_obs_uv_200, 2)), transform=ax2.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax2.set_title('(b) MME-best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, bias_mme_worse_obs_q_200, levels=np.arange(-0.02, 0.021, 0.001), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct3 = ax3.contour(lon, lat, bias_mme_worse_obs_uv_200, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax3.clabel(ct3, fontsize=font_size, colors='black')
ax3.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_worse_obs_q_200, 2), round(corr_mme_worse_obs_uv_200, 2)), transform=ax3.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax3.set_title('(c) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, bias_mme_obs_q_500, levels=np.arange(-2, 2.10, 0.10), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct4 = ax4.contour(lon, lat, bias_mme_obs_uv_500, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax4.clabel(ct4, fontsize=font_size, colors='black')
ax4.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_obs_q_500, 2), round(corr_mme_obs_uv_500, 2)), transform=ax4.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax4.set_title('(d) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel('500 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cf5 = ax5.contourf(lon, lat, bias_mme_best_obs_q_500, levels=np.arange(-2, 2.10, 0.10), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct5 = ax5.contour(lon, lat, bias_mme_best_obs_uv_500, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax5.clabel(ct5, fontsize=font_size, colors='black')
ax5.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_best_obs_q_500, 2), round(corr_mme_best_obs_uv_500, 2)), transform=ax5.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax5.set_title('(e) MME_best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)
cbar = plt.colorbar(cf5, cax=fig.add_axes([0.25, 0.363, 0.5, 0.015]), orientation='horizontal')
cbar.set_label('Specific humidity 500 hPa (g kg$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf6 = ax6.contourf(lon, lat, bias_mme_worse_obs_q_500, levels=np.arange(-2, 2.10, 0.10), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct6 = ax6.contour(lon, lat, bias_mme_best_obs_uv_500, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax6.clabel(ct6, fontsize=font_size, colors='black')
ax6.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_worse_obs_q_500, 2), round(corr_mme_worse_obs_uv_200, 2)), transform=ax6.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax6.set_title('(f) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

cf7 = ax7.contourf(lon, lat, bias_mme_obs_q_850, levels=np.arange(-4, 4.20, 0.20), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct7 = ax7.contour(lon, lat, bias_mme_obs_uv_850, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax7.clabel(ct7, fontsize=font_size, colors='black')
ax7.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_obs_q_850, 2), round(corr_mme_obs_uv_850, 2)), transform=ax7.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax7.set_title('(g) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_ylabel('850 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

cf8 = ax8.contourf(lon, lat, bias_mme_best_obs_q_850, levels=np.arange(-4, 4.20, 0.20), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct8 = ax8.contour(lon, lat, bias_mme_best_obs_uv_850, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax8.clabel(ct8, fontsize=font_size, colors='black')
ax8.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_best_obs_q_850, 2), round(corr_mme_best_obs_uv_850, 2)), transform=ax8.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax8.set_title('(h) MME-best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cf9 = ax9.contourf(lon, lat, bias_mme_worse_obs_q_850, levels=np.arange(-4, 4.20, 0.20), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct9 = ax9.contour(lon, lat, bias_mme_worse_obs_uv_850, levels=np.arange(-10, 10, 1), linewidths=0.6, colors='black')
ax9.clabel(ct9, fontsize=font_size, colors='black')
ax9.text(0.91, 0.04, '{0} ({1})'.format(round(corr_mme_worse_obs_q_850, 2), round(corr_mme_worse_obs_uv_850, 2)), transform=ax9.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax9.set_title('(i) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)
cbar = plt.colorbar(cf9, cax=fig.add_axes([0.25, 0.092, 0.5, 0.015]), orientation='horizontal')
cbar.set_label('Specific humidity 850 hPa (g kg$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

# Path out to save figure
path_out = '{0}/figs/paper_cmip6'.format(path)
name_out = 'pyplt_maps_atm_cmip6_obs_quv_{0}_bias.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

