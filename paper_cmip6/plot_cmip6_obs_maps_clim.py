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

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/database/obs/' + '{0}_sa_era5_mon_{1}_lonlat.nc'.format(param, date))
	data = arq[param]
	time = data.sel(time=slice('1979-01-01','2014-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
	
	
def import_cmip(param, model, exp, date):

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/database/cmip6/{0}/'.format(model) + '{0}_sa_Amon_{1}_historical_{2}_{3}_lonlat.nc'.format(param, model, exp, date))
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
	map.drawmeridians(np.arange(-100., -20., 10.), size=8, labels=[0,0,0,1], linewidth=0.4, color='black')
	map.drawparallels(np.arange(-60., 15., 10.), size=8, labels=[1,0,0,0], linewidth=0.4, color='black')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	# Import shapefile 
	map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=1.)
	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()
	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]
	map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul2', drawbounds=False)
	polys = polys + getattr(map, 'america_sul2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] # creating a PathPatch
	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)
	plt.gca().add_patch(patch)
	
	return map, xx, yy
	

# Import cmip models and obs database 
var_dict = {1 :['q',    'hus'],
			2 :['mslp', 'psl'],
			3 :['ps',   'ps'],
			4 :['t',    'ta'],
			5 :['t2m',  'tas'],
			6 :['u',    'ua'],
			7 :['v',    'va'],
			8 :['z',    'zg']}

var_obs = var_dict[5][0]
var_cmip6 = var_dict[5][1]
dt = '1979-2014'

lat, lon, mean_obs = import_obs(var_obs, dt)

mean_cmip6 = []

for i in range(1, 18):
	print(cmip6[i][0])
	lat, lon, mean_cmip = import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt)
	mean_cmip6.append(mean_cmip)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 3))

if var_cmip6 == 'pr':
	levs0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	color0 = cm.Blues
elif var_cmip6 == 'tasmax':
	levs0 = [18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38]
	color0 = cm.Reds
else:
	levs0 = [14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds

ax = fig.add_subplot(2, 6, 1)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_obs, levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(a) BR-DWGD', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 2)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[0], levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(b) NorESM2-MM \n Rank 1', loc='left', fontsize=7, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.91, 0.28, 0.02, 0.43]))
cbar.ax.tick_params(labelsize=7)

ax = fig.add_subplot(2, 6, 3)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[1], levels=levs0, latlon=True, cmap=color0, extend='both') 
plt.title(u'(c) GFDL-ESM4 \n Rank 2', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 4)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[2], levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(d) MPI-ESM1-2-HR \n Rank 3', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 5)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[3], levels=levs0, latlon=True, cmap=color0, extend='both') 
plt.title(u'(e) INM-CM5-0 \n Rank 4', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 6)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_cmip6[4], levels=levs0, latlon=True, cmap=color0, extend='max') 
plt.title(u'(f) MRI-ESM2-0 \n Rank 5', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 8)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[0], levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(g) NorESM2-MM \n Rank 1', loc='left', fontsize=7, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.97, 0.28, 0.02, 0.43]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(2, 6, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[1], levels=levs, latlon=True, cmap=color, extend='max') 
plt.title(u'(h) GFDL-ESM4 \n Rank 2', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[2], levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(i) MPI-ESM1-2-HR \n Rank 3', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[3], levels=levs, latlon=True, cmap=color, extend='max') 
plt.title(u'(j) INM-CM5-0 \n Rank 4', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[4], levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(l) MRI-ESM2-0 \n Rank 5', loc='left', fontsize=7, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/figs'
name_out = 'pyplt_maps_clim_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



