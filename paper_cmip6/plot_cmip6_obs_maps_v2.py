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

path  = '/afs/ictp.it/home/m/mda_silv/Documents/AdaptaBr_MCTI'
		    
			    
def import_obs_srf(param, area):

	arq  = '{0}/database/obs/{1}_{2}_ERA5_mon_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]

	mean_850 = np.nanmean(var[:][:,0,:,:], axis=0)
	mean_500 = np.nanmean(var[:][:,1,:,:], axis=0)
	mean_200 = np.nanmean(var[:][:,2,:,:], axis=0)
			
	return lat, lon, mean_850, mean_500, mean_200


def import_cmip_srf(param, area, model, exp):
	
	arq   = '{0}/database/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]

	mean_850 = np.nanmean(var[:][:,0,:,:], axis=0)
	mean_500 = np.nanmean(var[:][:,1,:,:], axis=0)
	mean_200 = np.nanmean(var[:][:,2,:,:], axis=0)
			
	return lat, lon, mean_850, mean_500, mean_200
    

def configure_subplot(ax):

	ax.set_extent(latlon, crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(latlon[0], latlon[1], 20), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(latlon[2], latlon[3], 20), crs=ccrs.PlateCarree())

	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(6)

	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.4)
	ax.add_feature(cfeat.BORDERS, linewidth=0.5)
	ax.coastlines(linewidth=0.5)
	
		
# Import obs database and cmip models
lat, lon, t_850_obs, t_500_obs, t_200_obs = import_obs_srf('t', 'sa')

lat, lon, t_850_mdl1b, t_500_mdl1b, t_200_mdl1b = import_cmip_srf('ta', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, t_850_mdl2b, t_500_mdl2b, t_200_mdl2b = import_cmip_srf('ta', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, t_850_mdl3b, t_500_mdl3b, t_200_mdl3b = import_cmip_srf('ta', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, t_850_mdl1w, t_500_mdl1w, t_200_mdl1w = import_cmip_srf('ta', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, t_850_mdl2w, t_500_mdl2w, t_200_mdl2w = import_cmip_srf('ta', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, t_850_mdl3w, t_500_mdl3w, t_200_mdl3w = import_cmip_srf('ta', 'sa', 'NESM3', 'r1i1p1f1_gn')

t_850_mme_best = np.nanmean([t_850_mdl1b, t_850_mdl2b, t_850_mdl3b], axis=0)
t_500_mme_best = np.nanmean([t_500_mdl1b, t_500_mdl2b, t_500_mdl3b], axis=0)
t_200_mme_best = np.nanmean([t_200_mdl1b, t_200_mdl2b, t_200_mdl3b], axis=0)

t_850_mme_worse = np.nanmean([t_850_mdl1w, t_850_mdl2w, t_850_mdl3w], axis=0)
t_500_mme_worse = np.nanmean([t_500_mdl1w, t_500_mdl2w, t_500_mdl3w], axis=0)
t_200_mme_worse = np.nanmean([t_200_mdl1w, t_200_mdl2w, t_200_mdl3w], axis=0)

print(np.min(t_850_obs), np.max(t_850_obs))
print(np.min(t_500_obs), np.max(t_500_obs))
print(np.min(t_200_obs), np.max(t_200_obs))

# Plot figure
fig, axes = plt.subplots(3, 3, figsize=(12, 10), subplot_kw={'projection': ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
cmap_color=cm.jet
font_size=8

cf1 = ax1.contourf(lon, lat, t_850_obs, levels=np.arange(-10, 30, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct1 = ax1.contour(lon, lat, t_850_obs, linewidths=0.5, colors='black')
ax1.clabel(ct1, fontsize=font_size, colors='black')
ax1.set_title('(a) ', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel('850 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)
cbar = plt.colorbar(cf1, cax=fig.add_axes([0.89, 0.3, 0.015, 0.4]))
cbar.set_label('Air temperature 850 hPa (°C)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

cf2 = ax2.contourf(lon, lat, t_850_mme_best, levels=np.arange(-10, 30, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct2 = ax2.contour(lon, lat, t_850_mme_best, linewidths=0.5, colors='black')
ax2.clabel(ct2, fontsize=font_size, colors='black')
ax2.set_title('(b) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, t_850_mme_worse, levels=np.arange(-10, 30, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct3 = ax3.contour(lon, lat, t_850_mme_worse, linewidths=0.5, colors='black')
ax3.clabel(ct3, fontsize=font_size, colors='black')
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cf4 = ax4.contourf(lon, lat, t_500_obs, levels=np.arange(-35, 5, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct4 = ax4.contour(lon, lat, t_500_obs, linewidths=0.5, colors='black')
ax4.clabel(ct4, fontsize=font_size, colors='black')
ax4.set_title('(d) ', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel('500 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)
cbar = plt.colorbar(cf4, cax=fig.add_axes([0.96, 0.3, 0.015, 0.4]))
cbar.set_label('Air temperature 500 hPa (°C)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

cf5 = ax5.contourf(lon, lat, t_500_mme_best, levels=np.arange(-35, 5, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct5 = ax5.contour(lon, lat, t_500_mme_best, linewidths=0.5, colors='black')
ax5.clabel(ct5, fontsize=font_size, colors='black')
ax5.set_title('(e) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

cf6 = ax6.contourf(lon, lat, t_500_mme_worse, levels=np.arange(-35, 5, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct6 = ax6.contour(lon, lat, t_500_mme_worse, linewidths=0.5, colors='black')
ax6.clabel(ct6, fontsize=font_size, colors='black')
ax6.set_title('(f) ', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

cf7 = ax7.contourf(lon, lat, t_200_obs, levels=np.arange(-70, -30, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct7 = ax7.contour(lon, lat, t_200_obs, linewidths=0.5, colors='black')
ax7.clabel(ct7, fontsize=font_size, colors='black')
ax7.set_title('(g) ', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel('ERA5', fontsize=font_size, fontweight='bold')
ax7.set_ylabel('200 hPa', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)
cbar = plt.colorbar(cf7, cax=fig.add_axes([1.03, 0.3, 0.015, 0.4]))
cbar.set_label('Air temperature 200 hPa (°C)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

cf8 = ax8.contourf(lon, lat, t_200_mme_best, levels=np.arange(-70, -30, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct8 = ax8.contour(lon, lat, t_200_mme_best, linewidths=0.5, colors='black')
ax8.clabel(ct8, fontsize=font_size, colors='black')
ax8.set_title('(h)', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel('MME-best', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cf9 = ax9.contourf(lon, lat, t_200_mme_worse, levels=np.arange(-70, -30, 1), transform=ccrs.PlateCarree(), extend='max', cmap=cmap_color)
ct9 = ax9.contour(lon, lat, t_200_mme_worse, linewidths=0.5, colors='black')
ax9.clabel(ct9, fontsize=font_size, colors='black')
ax9.set_title('(i) ', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_xlabel('MME-worst', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_atm_cmip6_obs_{0}_temp.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

