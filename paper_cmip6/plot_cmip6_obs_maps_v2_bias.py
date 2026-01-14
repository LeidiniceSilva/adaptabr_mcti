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
	mean_850 = annual_mean[:, 0, :, :]
	mean_500 = annual_mean[:, 1, :, :]
	mean_200 = annual_mean[:, 2, :, :]
		
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
lat, lon, t_850_obs_, t_500_obs_, t_200_obs_ = import_obs_srf('t', 'sa')

t_850_obs = np.nanmean(t_850_obs_, axis=0)
t_500_obs = np.nanmean(t_500_obs_, axis=0)
t_200_obs = np.nanmean(t_200_obs_, axis=0)

lat, lon, t_850_mdl1b, t_500_mdl1b, t_200_mdl1b = import_cmip_srf('ta', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, t_850_mdl2b, t_500_mdl2b, t_200_mdl2b = import_cmip_srf('ta', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, t_850_mdl3b, t_500_mdl3b, t_200_mdl3b = import_cmip_srf('ta', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, t_850_mdl1w, t_500_mdl1w, t_200_mdl1w = import_cmip_srf('ta', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, t_850_mdl2w, t_500_mdl2w, t_200_mdl2w = import_cmip_srf('ta', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, t_850_mdl3w, t_500_mdl3w, t_200_mdl3w = import_cmip_srf('ta', 'sa', 'NESM3', 'r1i1p1f1_gn')

t_850_mme_best = np.nanmean(np.nanmean([t_850_mdl1b, t_850_mdl2b, t_850_mdl3b], axis=0), axis=0)
t_500_mme_best = np.nanmean(np.nanmean([t_500_mdl1b, t_500_mdl2b, t_500_mdl3b], axis=0), axis=0)
t_200_mme_best = np.nanmean(np.nanmean([t_200_mdl1b, t_200_mdl2b, t_200_mdl3b], axis=0), axis=0)

t_850_mme_worse = np.nanmean(np.nanmean([t_850_mdl1w, t_850_mdl2w, t_850_mdl3w], axis=0), axis=0)
t_500_mme_worse = np.nanmean(np.nanmean([t_500_mdl1w, t_500_mdl2w, t_500_mdl3w], axis=0), axis=0)
t_200_mme_worse = np.nanmean(np.nanmean([t_200_mdl1w, t_200_mdl2w, t_200_mdl3w], axis=0), axis=0)

t_850_clim_, t_500_clim_, t_200_clim_ = [], [], []
for i in range(1, 18):
	lat, lon, t_850_clim, t_500_clim, t_200_clim = import_cmip_srf('ta', 'sa', cmip6[i][0], cmip6[i][1])
	t_850_clim_.append(t_850_clim)
	t_500_clim_.append(t_500_clim)
	t_200_clim_.append(t_200_clim)

t_850_mme = np.nanmean(np.nanmean([t_850_clim_[0], t_850_clim_[1], t_850_clim_[2], t_850_clim_[3], t_850_clim_[4],
t_850_clim_[5], t_850_clim_[6], t_850_clim_[7], t_850_clim_[8], t_850_clim_[9], t_850_clim_[10], t_850_clim_[11],
t_850_clim_[12], t_850_clim_[13], t_850_clim_[14], t_850_clim_[15], t_850_clim_[16]], axis=0), axis=0)

t_500_mme = np.nanmean(np.nanmean([t_500_clim_[0], t_500_clim_[1], t_500_clim_[2], t_500_clim_[3], t_500_clim_[4],
t_500_clim_[5], t_500_clim_[6], t_500_clim_[7], t_500_clim_[8], t_500_clim_[9], t_500_clim_[10], t_500_clim_[11],
t_500_clim_[12], t_500_clim_[13], t_500_clim_[14], t_500_clim_[15], t_500_clim_[16]], axis=0), axis=0)

t_200_mme = np.nanmean(np.nanmean([t_200_clim_[0], t_200_clim_[1], t_200_clim_[2], t_200_clim_[3], t_200_clim_[4],
t_200_clim_[5], t_200_clim_[6], t_200_clim_[7], t_200_clim_[8], t_200_clim_[9], t_200_clim_[10], t_200_clim_[11],
t_200_clim_[12], t_200_clim_[13], t_200_clim_[14], t_200_clim_[15], t_200_clim_[16]], axis=0), axis=0)

# Bias
bias_mme_obs_t_850 = t_850_mme - t_850_obs
bias_mme_obs_t_500 = t_500_mme - t_500_obs
bias_mme_obs_t_200 = t_200_mme - t_200_obs

bias_mme_best_obs_t_850 = t_850_mme_best - t_850_obs
bias_mme_best_obs_t_500 = t_500_mme_best - t_500_obs
bias_mme_best_obs_t_200 = t_200_mme_best - t_200_obs

bias_mme_worse_obs_t_850 = t_850_mme_worse - t_850_obs
bias_mme_worse_obs_t_500 = t_500_mme_worse - t_500_obs
bias_mme_worse_obs_t_200 = t_200_mme_worse - t_200_obs

# Corr
corr_mme_obs_t_850 = np.corrcoef(t_850_mme.flatten(), t_850_obs.flatten())[0, 1]
corr_mme_obs_t_500 = np.corrcoef(t_500_mme.flatten(), t_500_obs.flatten())[0, 1]
corr_mme_obs_t_200 = np.corrcoef(t_200_mme.flatten(), t_200_obs.flatten())[0, 1]

# Corr
corr_mme_best_obs_t_850 = np.corrcoef(t_850_mme_best.flatten(), t_850_obs.flatten())[0, 1]
corr_mme_best_obs_t_500 = np.corrcoef(t_500_mme_best.flatten(), t_500_obs.flatten())[0, 1]
corr_mme_best_obs_t_200 = np.corrcoef(t_200_mme_best.flatten(), t_200_obs.flatten())[0, 1]

# Corr
corr_mme_worse_obs_t_850 = np.corrcoef(t_850_mme_worse.flatten(), t_850_obs.flatten())[0, 1]
corr_mme_worse_obs_t_500 = np.corrcoef(t_500_mme_worse.flatten(), t_500_obs.flatten())[0, 1]
corr_mme_worse_obs_t_200 = np.corrcoef(t_200_mme_worse.flatten(), t_200_obs.flatten())[0, 1]

# Plot figure
fig, axes = plt.subplots(3, 3, figsize=(12, 13), subplot_kw={'projection': ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
cmap_color=cm.bwr

cf1 = ax1.contourf(lon, lat, bias_mme_obs_t_200, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct1 = ax1.contour(lon, lat, bias_mme_obs_t_200, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax1.clabel(ct1, fontsize=font_size, colors='black')
ax1.set_title('(a) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('200 hPa', fontsize=font_size, fontweight='bold')
ax1.text(0.91, 0.04, '{0}'.format(round(corr_mme_obs_t_200, 2)), transform=ax1.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)
cbar = plt.colorbar(cf1, cax=fig.add_axes([0.25, 0.634, 0.5, 0.015]), orientation='horizontal')
cbar.set_label('Air temperature 200 hPa (°C)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf2 = ax2.contourf(lon, lat, bias_mme_best_obs_t_200, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct2 = ax2.contour(lon, lat, bias_mme_best_obs_t_200, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax2.clabel(ct2, fontsize=font_size, colors='black')
ax2.text(0.91, 0.04, '{0}'.format(round(corr_mme_best_obs_t_200, 2)), transform=ax2.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax2.set_title('(b) MME-best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, bias_mme_worse_obs_t_200, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct3 = ax3.contour(lon, lat, bias_mme_worse_obs_t_200, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax3.clabel(ct3, fontsize=font_size, colors='black')
ax3.text(0.91, 0.04, '{0}'.format(round(corr_mme_worse_obs_t_200, 2)), transform=ax3.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax3.set_title('(c) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, bias_mme_obs_t_500, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct4 = ax4.contour(lon, lat, bias_mme_obs_t_500, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax4.clabel(ct4, fontsize=font_size, colors='black')
ax4.text(0.91, 0.04, '0.9', transform=ax4.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax4.set_title('(d) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel('500 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cf5 = ax5.contourf(lon, lat, bias_mme_best_obs_t_500, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct5 = ax5.contour(lon, lat, bias_mme_best_obs_t_500, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax5.clabel(ct5, fontsize=font_size, colors='black')
ax5.text(0.91, 0.04, '0.92', transform=ax5.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax5.set_title('(e) MME_best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)
cbar = plt.colorbar(cf5, cax=fig.add_axes([0.25, 0.363, 0.5, 0.015]), orientation='horizontal')
cbar.set_label('Air temperature 500 hPa (°C)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf6 = ax6.contourf(lon, lat, bias_mme_worse_obs_t_500, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct6 = ax6.contour(lon, lat, bias_mme_worse_obs_t_500, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax6.clabel(ct6, fontsize=font_size, colors='black')
ax6.text(0.91, 0.04, '0.90', transform=ax6.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax6.set_title('(f) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

cf7 = ax7.contourf(lon, lat, bias_mme_obs_t_850, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct7 = ax7.contour(lon, lat, bias_mme_obs_t_850, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax7.clabel(ct7, fontsize=font_size, colors='black')
ax7.text(0.91, 0.04, '0.91', transform=ax7.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax7.set_title('(g) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_ylabel('850 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

cf8 = ax8.contourf(lon, lat, bias_mme_best_obs_t_850, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct8 = ax8.contour(lon, lat, bias_mme_best_obs_t_850, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax8.clabel(ct8, fontsize=font_size, colors='black')
ax8.text(0.91, 0.04, '0.91', transform=ax8.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax8.set_title('(h) MME-best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cf9 = ax9.contourf(lon, lat, bias_mme_worse_obs_t_850, levels=np.arange(-9, 9.5, 0.5), transform=ccrs.PlateCarree(), extend='both', cmap=cmap_color)
ct9 = ax9.contour(lon, lat, bias_mme_worse_obs_t_850, levels=np.arange(-9, 9.5, 0.5), linewidths=0.6, colors='black')
ax9.clabel(ct9, fontsize=font_size, colors='black')
ax9.text(0.91, 0.04, '0.91', transform=ax9.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
ax9.set_title('(i) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)
cbar = plt.colorbar(cf9, cax=fig.add_axes([0.25, 0.092, 0.5, 0.015]), orientation='horizontal')
cbar.set_label('Air temperature 850 hPa (°C)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

# Path out to save figure
path_out = '{0}/figs/paper_cmip6'.format(path)
name_out = 'pyplt_maps_atm_cmip6_obs_temp_{0}_bias.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

