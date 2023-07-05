# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 14, 2023"
__description__ = "This script plot spatial climatology maps"

import os
import netCDF4
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from dict_cmip6_models_name import cmip6


def import_obs(param, date):

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/database/paper_cmip6/obs/' + '{0}_sa_era5_mon_{1}_lonlat.nc'.format(param, date))
	data = arq[param]
	time = data.sel(time=slice('1979-01-01','2014-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
		
	
def import_cmip(param, model, exp, date):

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/database/paper_cmip6/cmip6/{0}/'.format(model) + '{0}_sa_Amon_{1}_historical_{2}_{3}_lonlat.nc'.format(param, model, exp, date))
	data = arq[param]
	time = data.sel(time=slice('1979-01-01','2014-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean


def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]

	map = Basemap(projection='cyl', llcrnrlon=-100., llcrnrlat=-60., urcrnrlon=-20., urcrnrlat=15, resolution='c')
	map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=1.)
	map.drawmeridians(np.arange(-100., -20., 20.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-60., 15., 15.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	
	return map, xx, yy
	

# Import cmip models and obs database
var_cmip6 = 'tas'
dt = '197901-201412'

if var_cmip6 == 'psl':
	var_obs = 'msl'
elif var_cmip6 == 'ps':
	var_obs = 'sp'
else:
	var_obs = 't2m'
	
lat, lon, mean_obs = import_obs(var_obs, dt)

mean_cmip6 = []
for i in range(1, 18):
	print(cmip6[i][0])
	lat, lon, mean_cmip = import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt)
	mean_cmip6.append(mean_cmip)
			
# Plot cmip models and obs database 
fig = plt.figure(figsize=(7, 9))
font_size = 8

if var_cmip6 == 'psl':
	levs = [977, 980, 983, 986, 989, 992, 995, 998, 1001, 1004, 1007, 1010, 1013, 1016, 1019, 1022, 1025]
	color = cm.rainbow
	legend = 'Mean sea level pressure (hPa)'
elif var_cmip6 == 'ps':
	levs = [921, 928, 935, 942, 948, 955, 962, 966, 969, 976, 983, 990, 997, 1004, 1011, 1018, 1025]
	color = cm.rainbow
	legend = 'Surface pressure (hPa)'
else:
	levs = [-2, 0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
	color = cm.jet
	legend = 'Near-surface air temperature (°C)'
		
ax = fig.add_subplot(5, 4, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_obs, levels=levs, latlon=True, cmap=color, extend=None) 
plt.title(u'(a) ERA5', loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.28, 0.02, 0.43]))
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(5, 4, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[0], levels=levs, latlon=True, cmap=color) 
plt.title(u'(b) {0}'.format(cmip6[1][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[1], levels=levs, latlon=True, cmap=color) 
plt.title(u'(c) {0}'.format(cmip6[2][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[2], levels=levs, latlon=True, cmap=color) 
plt.title(u'(d) {0}'.format(cmip6[3][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[3], levels=levs, latlon=True, cmap=color) 
plt.title(u'(e) {0}'.format(cmip6[4][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[4], levels=levs, latlon=True, cmap=color) 
plt.title(u'(f) {0}'.format(cmip6[5][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 7)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[5], levels=levs, latlon=True, cmap=color) 
plt.title(u'(g) {0}'.format(cmip6[6][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[6], levels=levs, latlon=True, cmap=color) 
plt.title(u'(h) {0}'.format(cmip6[7][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[7], levels=levs, latlon=True, cmap=color) 
plt.title(u'(i) {0}'.format(cmip6[8][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[8], levels=levs, latlon=True, cmap=color) 
plt.title(u'(j) {0}'.format(cmip6[9][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[9], levels=levs, latlon=True, cmap=color) 
plt.title(u'(k) {0}'.format(cmip6[10][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[10], levels=levs, latlon=True, cmap=color) 
plt.title(u'(l) {0}'.format(cmip6[11][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 13)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[11], levels=levs, latlon=True, cmap=color) 
plt.title(u'(m) {0}'.format(cmip6[12][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 14)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[12], levels=levs, latlon=True, cmap=color) 
plt.title(u'(n) {0}'.format(cmip6[13][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 15)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[13], levels=levs, latlon=True, cmap=color) 
plt.title(u'(o) {0}'.format(cmip6[14][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 16)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[14], levels=levs, latlon=True, cmap=color) 
plt.title(u'(p) {0}'.format(cmip6[15][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 17)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[15], levels=levs, latlon=True, cmap=color) 
plt.title(u'(q) {0}'.format(cmip6[16][0]), loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(5, 4, 18)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[16], levels=levs, latlon=True, cmap=color) 
plt.title(u'(r) {0}'.format(cmip6[17][0]), loc='left', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/paper_cmip6'
name_out = 'pyplt_maps_cmip6_clim_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()



