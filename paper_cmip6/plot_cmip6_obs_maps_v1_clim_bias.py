# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "June 10, 2024"
__description__ = "This script plot portrait of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.colors
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
	var_reshaped = var[:nyears*12].reshape(nyears, 12, len(lat), len(lon))
	annual_mean = np.nanmean(var_reshaped, axis=1)
	climatology = np.nanmean(annual_mean, axis=0)
	data.close()
	
	return lat, lon, annual_mean


def import_cmip_srf(param, area, model, exp):
	
	arq   = '{0}/database/paper_cmip6/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	time = data.variables['time'][:]
	ntime = len(time)
	nyears = ntime // 12
	var_reshaped = var[:nyears*12].reshape(nyears, 12, len(lat), len(lon))
	annual_mean = np.nanmean(var_reshaped, axis=1)  
	climatology = np.nanmean(annual_mean, axis=0)
	data.close()
	
	return lat, lon, annual_mean
    
	
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
lat, lon, pr_obs_ = import_obs_srf('tp', 'sa')
lat, lon, ps_obs_ = import_obs_srf('msl', 'sa')

pr_obs = np.nanmean(pr_obs_, axis=0)
ps_obs = np.nanmean(ps_obs_, axis=0)

lat, lon, pr_mdl1b = import_cmip_srf('pr', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, pr_mdl2b = import_cmip_srf('pr', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, pr_mdl3b = import_cmip_srf('pr', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, ps_mdl1b = import_cmip_srf('psl', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, ps_mdl2b = import_cmip_srf('psl', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, ps_mdl3b = import_cmip_srf('psl', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, pr_mdl1w = import_cmip_srf('pr', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, pr_mdl2w = import_cmip_srf('pr', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, pr_mdl3w = import_cmip_srf('pr', 'sa', 'NESM3', 'r1i1p1f1_gn')

lat, lon, ps_mdl1w = import_cmip_srf('psl', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, ps_mdl2w = import_cmip_srf('psl', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, ps_mdl3w = import_cmip_srf('psl', 'sa', 'NESM3', 'r1i1p1f1_gn')

pr_mme_best = np.nanmean(np.nanmean([pr_mdl1b, pr_mdl2b, pr_mdl3b], axis=0), axis=0)
pr_mme_worse = np.nanmean(np.nanmean([pr_mdl1w, pr_mdl2w, pr_mdl3w], axis=0), axis=0)

ps_mme_best = np.nanmean(np.nanmean([ps_mdl1b, ps_mdl2b, ps_mdl3b], axis=0), axis=0)
ps_mme_worse = np.nanmean(np.nanmean([ps_mdl1w, ps_mdl2w, ps_mdl3w], axis=0), axis=0)

pr_clim_, ps_clim_  = [], []
for i in range(1, 18):
	lat, lon, pr_clim = import_cmip_srf('pr', 'sa', cmip6[i][0], cmip6[i][1])
	pr_clim_.append(pr_clim)
	
	lat, lon, ps_clim = import_cmip_srf('psl', 'sa', cmip6[i][0], cmip6[i][1])
	ps_clim_.append(ps_clim)

pr_mme = np.nanmean(np.nanmean([pr_clim_[0], pr_clim_[1], pr_clim_[2], pr_clim_[3], pr_clim_[4],
pr_clim_[5], pr_clim_[6], pr_clim_[7], pr_clim_[8], pr_clim_[9], pr_clim_[10], pr_clim_[11],
pr_clim_[12], pr_clim_[13], pr_clim_[14], pr_clim_[15], pr_clim_[16]], axis=0), axis=0)

ps_mme = np.nanmean(np.nanmean([ps_clim_[0], ps_clim_[1], ps_clim_[2], ps_clim_[3], ps_clim_[4],
ps_clim_[5], ps_clim_[6], ps_clim_[7], ps_clim_[8], ps_clim_[9], ps_clim_[10], ps_clim_[11],
ps_clim_[12], ps_clim_[13], ps_clim_[14], ps_clim_[15], ps_clim_[16]], axis=0), axis=0)

# Bias
bias_mme_obs_pr = pr_mme - pr_obs
bias_mme_best_obs_pr = pr_mme_best - pr_obs
bias_mme_worse_obs_pr = pr_mme_worse - pr_obs

bias_mme_obs_ps = ps_mme - ps_obs
bias_mme_best_obs_ps = ps_mme_best - ps_obs
bias_mme_worse_obs_ps = ps_mme_worse - ps_obs

# Corr
corr_mme_obs_pr = np.corrcoef(pr_mme.flatten(), pr_obs.flatten())[0, 1]
corr_mme_best_obs_pr = np.corrcoef(pr_mme_best.flatten(), pr_obs.flatten())[0, 1]
corr_mme_worse_obs_pr = np.corrcoef(pr_mme_worse.flatten(), pr_obs.flatten())[0, 1]

corr_mme_obs_ps = np.corrcoef(ps_mme.flatten(), ps_obs.flatten())[0, 1]
corr_mme_best_obs_ps = np.corrcoef(ps_mme_best.flatten(), ps_obs.flatten())[0, 1]
corr_mme_worse_obs_ps = np.corrcoef(ps_mme_worse.flatten(), ps_obs.flatten())[0, 1]

print(corr_mme_obs_pr)
print(corr_mme_best_obs_pr)
print(corr_mme_worse_obs_pr)

print(corr_mme_obs_ps)
print(corr_mme_best_obs_ps)
print(corr_mme_worse_obs_ps)

# Plot figure
fig, axes = plt.subplots(4, 4, figsize=(14, 16), subplot_kw={'projection': ccrs.PlateCarree()})
(ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12), (ax13, ax14, ax15, ax16) = axes
fig.delaxes(ax5)
fig.delaxes(ax13)

color=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']

cf1 = ax1.contourf(lon, lat, pr_obs, levels=np.arange(0, 18, 1), transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct1 = ax1.contour(lon, lat, pr_obs, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax1.clabel(ct1, fontsize=font_size, colors='black')
ax1.set_title('(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cf2 = ax2.contourf(lon, lat, pr_mme, levels=np.arange(0, 18, 1), transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct2 = ax2.contour(lon, lat, pr_mme, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax2.clabel(ct2, fontsize=font_size, colors='black')
ax2.set_title('(b) MME', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, pr_mme_best, levels=np.arange(0, 18, 1), transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct3 = ax3.contour(lon, lat, pr_mme_best, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax3.clabel(ct3, fontsize=font_size, colors='black')
ax3.set_title('(c) MME-best', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, pr_mme_worse, levels=np.arange(0, 18, 1), transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct4 = ax4.contour(lon, lat, pr_mme_worse, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax4.clabel(ct3, fontsize=font_size, colors='black')
ax4.set_title('(d) MME-worst', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

cbar = plt.colorbar(cf4, cax=fig.add_axes([0.15, 0.7, 0.7, 0.01]), orientation='horizontal')
cbar.set_label('Mean precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf6 = ax6.contourf(lon, lat, bias_mme_obs_pr, levels=np.arange(-9, 10, 1), transform=ccrs.PlateCarree(), extend='both', cmap=cm.BrBG)
ax6.set_title('(e) MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax6.text(0.96, 0.04, '0.82', transform=ax6.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

cf7 = ax7.contourf(lon, lat, bias_mme_best_obs_pr, levels=np.arange(-9, 10, 1), transform=ccrs.PlateCarree(), extend='both', cmap=cm.BrBG)
ax7.set_title('(f) MME-best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax7.text(0.96, 0.04, '0.83', transform=ax7.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

cf8 = ax8.contourf(lon, lat, bias_mme_worse_obs_pr, levels=np.arange(-9, 10, 1), transform=ccrs.PlateCarree(), extend='both', cmap=cm.BrBG)
ax8.set_title('(g) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax8.text(0.96, 0.04, '0.74', transform=ax8.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cbar = plt.colorbar(cf8, cax=fig.add_axes([0.15, 0.5, 0.7, 0.01]), orientation='horizontal')
cbar.set_label('Bias of precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf9 = ax9.contourf(lon, lat, ps_obs, levels=np.arange(990, 1026, 2), transform=ccrs.PlateCarree(), extend='max', cmap=cm.BuPu)
ct9 = ax9.contour(lon, lat, ps_obs, levels=np.arange(990, 1034, 4), linewidths=0.6, colors='black')
ax9.clabel(ct9, fontsize=font_size, colors='black')
ax9.set_title('(h) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

cf10 = ax10.contourf(lon, lat, ps_mme, levels=np.arange(990, 1026, 2), transform=ccrs.PlateCarree(), extend='max', cmap=cm.BuPu)
ct10 = ax10.contour(lon, lat, ps_mme, levels=np.arange(990, 1034, 4), linewidths=0.6, colors='black')
ax10.clabel(ct10, fontsize=font_size, colors='black')
ax10.set_title('(i) MME', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax10)

cf11 = ax11.contourf(lon, lat, ps_mme_best, levels=np.arange(990, 1026, 2), transform=ccrs.PlateCarree(), extend='max', cmap=cm.BuPu)
ct11 = ax11.contour(lon, lat, ps_mme_best, levels=np.arange(990, 1034, 4), linewidths=0.6, colors='black')
ax11.clabel(ct11, fontsize=font_size, colors='black')
ax11.set_title('(j) MME-best', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax11)

cf12 = ax12.contourf(lon, lat, ps_mme_worse, levels=np.arange(990, 1026, 2), transform=ccrs.PlateCarree(), extend='max', cmap=cm.BuPu)
ct12 = ax12.contour(lon, lat, ps_mme_worse, levels=np.arange(990, 1034, 4), linewidths=0.6, colors='black')
ax12.clabel(ct12, fontsize=font_size, colors='black')
ax12.set_title('(k) MME-worst', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax12)

cbar = plt.colorbar(cf12, cax=fig.add_axes([0.15, 0.3, 0.7, 0.01]), orientation='horizontal')
cbar.set_label('Mean sea level pressure (hPa)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

cf14 = ax14.contourf(lon, lat, bias_mme_obs_ps, levels=np.arange(-9, 10, 1), transform=ccrs.PlateCarree(), extend='both', cmap=cm.PuOr)
ax14.set_title('(l)  MME - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax14.text(0.96, 0.04, '0.9', transform=ax14.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax14)

cf15 = ax15.contourf(lon, lat, bias_mme_best_obs_ps, levels=np.arange(-9, 10, 1), transform=ccrs.PlateCarree(), extend='both', cmap=cm.PuOr)
ax15.set_title('(m) MME-best - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax15.text(0.89, 0.04, '0.91', transform=ax15.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax15)

cf16 = ax16.contourf(lon, lat, bias_mme_worse_obs_ps, levels=np.arange(-9, 10, 1), transform=ccrs.PlateCarree(), extend='both', cmap=cm.PuOr)
ax16.set_title('(n) MME-worst - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax16.text(0.91, 0.04, '0.88', transform=ax16.transAxes, ha='right', va='bottom', fontsize=font_size, fontweight='bold')
configure_subplot(ax16)

cbar = plt.colorbar(cf16, cax=fig.add_axes([0.15, 0.1, 0.7, 0.01]), orientation='horizontal')
cbar.set_label('Bias of sea level pressure (hPa)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size, direction='in')

# Path out to save figure
path_out = '{0}/figs/paper_cmip6'.format(path)
name_out = 'pyplt_maps_atm_cmip6_obs_pr_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


