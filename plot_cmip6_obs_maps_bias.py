# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot annual bias maps of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from dict_cmip6_models_name import cmip6


def import_obs(param, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_SA_CRU_ts4_ANN_{2}_lonlat.nc'.format(path, param, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(value[:,:,:], axis=0)

	return lat, lon, mean

	
def import_cmip(param, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_SA_{2}_historical_{3}_ANN_{4}_lonlat.nc'.format(path, param, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = np.nanmean(value[:,:,:], axis=0)
	
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
var_obs = 'pre'
var_cmip6 = 'pr'
dt = '1980-2014'

lat, lon, mean_obs = import_obs(var_obs, dt)

bias_cmip6 = []
for i in range(1, 19):
	print(cmip6[i][0])
	lat, lon, mean_cmip = import_cmip(var_cmip6, cmip6[i][0], cmip6[i][1], dt)
	bias_cmip = mean_cmip - mean_obs
	bias_cmip6.append(bias_cmip)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(7, 9))

if var_cmip6 == 'pr':
	levs0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
	color0 = cm.Blues
	levs = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
	color = cm.RdBu
else:
	levs0 = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]
	color0 = cm.Reds
	levs = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]
	color = cm.bwr
	
ax = fig.add_subplot(5, 4, 1)  
plt.title(u'(a) CRU', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean_obs, levels=levs0, latlon=True, cmap=color0, extend='max') 
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.28, 0.02, 0.43]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 4, 2)  
plt.title(u'(b) {0}'.format(cmip6[1][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[0], levels=levs, latlon=True, cmap=color, extend='both') 
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.98, 0.28, 0.02, 0.43]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 4, 3)  
plt.title(u'(c) {0}'.format(cmip6[2][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[1], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 4)  
plt.title(u'(d) {0}'.format(cmip6[3][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[2], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 5)  
plt.title(u'(e) {0}'.format(cmip6[4][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[3], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 6)  
plt.title(u'(f) {0}'.format(cmip6[5][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[4], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 7)  
plt.title(u'(g) {0}'.format(cmip6[6][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[5], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 8)  
plt.title(u'(h) {0}'.format(cmip6[7][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[6], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 9)  
plt.title(u'(i) {0}'.format(cmip6[8][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[7], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 10)  
plt.title(u'(j) {0}'.format(cmip6[9][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[8], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 11)  
plt.title(u'(k) {0}'.format(cmip6[10][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[9], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 12)  
plt.title(u'(l) {0}'.format(cmip6[11][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[10], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 13)  
plt.title(u'(m) {0}'.format(cmip6[12][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[11], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 14)  
plt.title(u'(n) {0}'.format(cmip6[13][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[12], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 15)  
plt.title(u'(o) {0}'.format(cmip6[14][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[13], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 16)  
plt.title(u'(p) {0}'.format(cmip6[15][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[14], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 17)  
plt.title(u'(q) {0}'.format(cmip6[16][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[15], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 18)  
plt.title(u'(e) {0}'.format(cmip6[17][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[16], levels=levs, latlon=True, cmap=color) 

ax = fig.add_subplot(5, 4, 19)  
plt.title(u'(s) {0}'.format(cmip6[18][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6[17], levels=levs, latlon=True, cmap=color) 

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_maps_bias_ann_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()




