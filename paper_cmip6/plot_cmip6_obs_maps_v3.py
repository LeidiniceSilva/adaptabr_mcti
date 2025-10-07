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

path  = '/afs/ictp.it/home/m/mda_silv/Documents/AdaptaBr_MCTI'
		    
			    
def import_obs_srf(param, area):

	arq  = '{0}/database/obs/{1}_{2}_ERA5_mon_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, mean


def import_cmip_srf(param, area, model, exp):
	
	arq   = '{0}/database/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, mean
    

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
lat, lon, pr_obs = import_obs_srf('tp', 'sa')
lat, lon, ps_obs = import_obs_srf('sp', 'sa')

lat, lon, pr_mdl1b = import_cmip_srf('pr', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, pr_mdl2b = import_cmip_srf('pr', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, pr_mdl3b = import_cmip_srf('pr', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, ps_mdl1b = import_cmip_srf('ps', 'sa', 'MRI-ESM2-0', 'r1i1p1f1_gn')
lat, lon, ps_mdl2b = import_cmip_srf('ps', 'sa', 'MPI-ESM1-2-HR', 'r1i1p1f1_gn')
lat, lon, ps_mdl3b = import_cmip_srf('ps', 'sa', 'ACCESS-CM2', 'r1i1p1f1_gn')

lat, lon, pr_mdl1w = import_cmip_srf('pr', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, pr_mdl2w = import_cmip_srf('pr', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, pr_mdl3w = import_cmip_srf('pr', 'sa', 'NESM3', 'r1i1p1f1_gn')

lat, lon, ps_mdl1w = import_cmip_srf('ps', 'sa', 'KIOST-ESM', 'r1i1p1f1_gr1')
lat, lon, ps_mdl2w = import_cmip_srf('ps', 'sa', 'CanESM5', 'r1i1p1f1_gn')
lat, lon, ps_mdl3w = import_cmip_srf('ps', 'sa', 'NESM3', 'r1i1p1f1_gn')

pr_mme_best = np.nanmean([pr_mdl1b, pr_mdl2b, pr_mdl3b], axis=0)
pr_mme_worse = np.nanmean([pr_mdl1w, pr_mdl2w, pr_mdl3w], axis=0)

ps_mme_best = np.nanmean([ps_mdl1b, ps_mdl2b, ps_mdl3b], axis=0)
ps_mme_worse = np.nanmean([ps_mdl1w, ps_mdl2w, ps_mdl3w], axis=0)

print(pr_obs.shape)
print(pr_mme_best.shape)
print(pr_mme_worse.shape)
print()

# Plot figure
fig, axes = plt.subplots(1, 3, figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
(ax1, ax2, ax3) = axes
levs=np.arange(0, 18.5, 0.5)
color=['#ffffffff','#d7f0fcff','#ade0f7ff','#86c4ebff','#60a5d6ff','#4794b3ff','#49a67cff','#55b848ff','#9ecf51ff','#ebe359ff','#f7be4aff','#f58433ff','#ed5a28ff','#de3728ff','#cc1f27ff','#b01a1fff','#911419ff']
font_size=8

cf1 = ax1.contourf(lon, lat, pr_obs, levels=levs, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct1 = ax1.contour(lon, lat, pr_obs, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax1.clabel(ct1, fontsize=font_size, colors='black')
ax1.set_title('(a) ', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_xlabel('ERA5', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

cf2 = ax2.contourf(lon, lat, pr_mme_best, levels=levs, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct2 = ax2.contour(lon, lat, pr_mme_best, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax2.clabel(ct2, fontsize=font_size, colors='black')
ax2.set_title('(b) ', loc='left', fontsize=font_size, fontweight='bold')
ax2.set_xlabel('MME-best', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

cf3 = ax3.contourf(lon, lat, pr_mme_worse, levels=levs, transform=ccrs.PlateCarree(), extend='max', cmap=matplotlib.colors.ListedColormap(color))
ct3 = ax3.contour(lon, lat, pr_mme_worse, levels=np.arange(0, 25, 5), linewidths=0.6, colors='black')
ax3.clabel(ct3, fontsize=font_size, colors='black')
ax3.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax3.set_xlabel('MME-worst', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

cbar = plt.colorbar(cf3, cax=fig.add_axes([0.91, 0.35, 0.015, 0.3]))
cbar.set_label('Precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_maps_atm_cmip6_obs_{0}_pr.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

