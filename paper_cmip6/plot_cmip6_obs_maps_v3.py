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

# Plot figure
fig, axes = plt.subplots(3, 4, figsize=(14, 12), subplot_kw={'projection': ccrs.PlateCarree()})
(ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12) = axes
cmap_color=cm.hsv

cf1 = ax1.contourf(lon, lat, q_850_obs, levels=np.arange(0, 20.25, 0.25), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st1 = ax1.streamplot(lon, lat, u_850_obs[:,:], v_850_obs[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax1.set_title('(a) ', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('850 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)
cbar = plt.colorbar(cf1, cax=fig.add_axes([0.25, 0.634, 0.5, 0.02]), orientation='horizontal')
cbar.set_label('Specific humidity 850 hPa (g kg$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf2 = ax2.contourf(lon, lat, q_850_mme, levels=np.arange(0, 20.25, 0.25), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st2 = ax2.streamplot(lon, lat, u_850_mme[:,:], v_850_mme[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax2.set_title('(b) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, q_850_mme_best, levels=np.arange(0, 20.25, 0.25), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st3 = ax3.streamplot(lon, lat, u_850_mme_best[:,:], v_850_mme_best[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, q_850_mme_worse, levels=np.arange(0, 20.25, 0.25), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st4 = ax4.streamplot(lon, lat, u_850_mme_worse[:,:], v_850_mme_worse[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax4.set_title('(d) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cf5 = ax5.contourf(lon, lat, q_500_obs, levels=np.arange(0, 5.0625, 0.0625), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st5 = ax5.streamplot(lon, lat, u_500_obs[:,:], v_500_obs[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax5.set_title('(e) ', loc='left', fontsize=font_size, fontweight='bold')
ax5.set_ylabel('500 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)
cbar = plt.colorbar(cf5, cax=fig.add_axes([0.25, 0.363, 0.5, 0.02]), orientation='horizontal')
cbar.set_label('Specific humidity 500 hPa (g kg$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf6 = ax6.contourf(lon, lat, q_500_mme, levels=np.arange(0, 5.0625, 0.0625), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st6 = ax6.streamplot(lon, lat, u_500_mme[:,:], v_500_mme[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax6.set_title('(f) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

cf7 = ax7.contourf(lon, lat, q_500_mme_best, levels=np.arange(0, 5.0625, 0.0625), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st7 = ax7.streamplot(lon, lat, u_500_mme_best[:,:], v_500_mme_best[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax7.set_title('(g) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

cf8 = ax8.contourf(lon, lat, q_500_mme_worse, levels=np.arange(0, 5.0625, 0.0625), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st8 = ax8.streamplot(lon, lat, u_500_mme_worse[:,:], v_500_mme_worse[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax8.set_title('(h)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cf9 = ax9.contourf(lon, lat, q_200_obs, levels=np.arange(0, 0.1, 0.00125), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st9 = ax9.streamplot(lon, lat, u_200_obs[:,:], v_200_obs[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax9.set_title('(i) ', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_xlabel('ERA5', fontsize=font_size, fontweight='bold')
ax9.set_ylabel('200 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)
cbar = plt.colorbar(cf9, cax=fig.add_axes([0.25, 0.075, 0.5, 0.02]), orientation='horizontal')
cbar.set_label('Specific humidity 200 hPa (g kg$^-$$^1$)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf10 = ax10.contourf(lon, lat, q_200_mme, levels=np.arange(0, 0.1, 0.00125), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st10 = ax10.streamplot(lon, lat, u_200_mme[:,:], v_200_mme[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax10.set_title('(j) ', loc='left', fontsize=font_size, fontweight='bold')
ax10.set_xlabel('MME', fontsize=font_size, fontweight='bold')
configure_subplot(ax10)

cf11 = ax11.contourf(lon, lat, q_200_mme_best, levels=np.arange(0, 0.1, 0.00125), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st11 = ax11.streamplot(lon, lat, u_200_mme_best[:,:], v_200_mme_best[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax11.set_title('(k) ', loc='left', fontsize=font_size, fontweight='bold')
ax11.set_xlabel('MME-best', fontsize=font_size, fontweight='bold')
configure_subplot(ax11)

cf12 = ax12.contourf(lon, lat, q_200_mme_worse, levels=np.arange(0, 0.1, 0.00125), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
st12 = ax12.streamplot(lon, lat, u_200_mme_worse[:,:], v_200_mme_worse[:,:], arrowsize=1, arrowstyle='->', color='black', density=2, linewidth=0.5)
ax12.set_title('(l) ', loc='left', fontsize=font_size, fontweight='bold')
ax12.set_xlabel('MME-worst', fontsize=font_size, fontweight='bold')
configure_subplot(ax12)

# Path out to save figure
path_out = '{0}/figs/paper_cmip6'.format(path)
name_out = 'pyplt_maps_atm_cmip6_obs_quv_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

