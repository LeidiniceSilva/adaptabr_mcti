# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot annual bias maps of cmip6 models"

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

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/project/database/obs/' + '{0}_SA_BR-DWGD_UFES_UTEXAS_v_3.0_MON_{1}_lonlat.nc'.format(param, date))
	data = arq[param]
	time = data.sel(time=slice('1986-01-01','2005-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean
	
	
def import_cmip(param, model, exp, date):

	arq = xr.open_dataset('/home/nice/Documentos/AdaptaBrasil_MCTI/project/database/cmip6/' + '{0}_SA_{1}_historical_{2}_MON_{3}_lonlat.nc'.format(param, model, exp, date))
	data = arq[param]
	time = data.sel(time=slice('1986-01-01','2005-12-31'))
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

	map = Basemap(projection='cyl', llcrnrlon=-76., llcrnrlat=-36., urcrnrlon=-32.,urcrnrlat=7., resolution='c')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	# Import shapefile 
	map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='white')
	map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black')
	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()
	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]
	map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed2', drawbounds=False)
	polys = polys + getattr(map, 'lim_unid_fed2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] # creating a PathPatch
	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)
	plt.gca().add_patch(patch)
	
	return map, xx, yy
	

# Import cmip models and obs database 
var_obs = 'pr'
var_cmip6 = 'pr'
dt = '1986-2005'

lat, lon, mean_obs = import_obs(var_obs, dt)

mean_cmip6 = []
bias_cmip6 = []

best_models = [17, 7, 13, 9, 15]
for i in best_models:
	print(cmip6[i][0])
	lat, lon, mean_cmip = import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt)
	mean_cmip6.append(mean_cmip)
	bias_cmip6.append(mean_cmip - mean_obs)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 3))

if var_cmip6 == 'pr':
	levs0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	color0 = cm.Blues
	levs = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
	color = cm.BrBG
elif var_cmip6 == 'tasmax':
	levs0 = [18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38]
	color0 = cm.Reds
	levs = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
	color = cm.bwr
else:
	levs0 = [14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34]
	color0 = cm.Reds
	levs = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
	color = cm.bwr

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
plt.title(u'(g) Rank 1', loc='left', fontsize=7, fontweight='bold')
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.97, 0.28, 0.02, 0.43]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(2, 6, 9)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[1], levels=levs, latlon=True, cmap=color, extend='max') 
plt.title(u'(h) Rank 2', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 10)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[2], levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(i) Rank 3', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 11)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[3], levels=levs, latlon=True, cmap=color, extend='max') 
plt.title(u'(j) Rank 4', loc='left', fontsize=7, fontweight='bold')

ax = fig.add_subplot(2, 6, 12)  
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[4], levels=levs, latlon=True, cmap=color, extend='both') 
plt.title(u'(l) Rank 5', loc='left', fontsize=7, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/project/figs/figs_report-II'
name_out = 'pyplt_best_cmip6_{0}_{1}_i.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



