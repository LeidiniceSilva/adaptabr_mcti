# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 14, 2022"
__description__ = "This script plot annual climatology maps from cmip6"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap

def import_obs(param):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_SA_cru_ts4.05_obs_yr_1961-2014_lonlat.nc'.format(path, param)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean1 = np.nanmean(value[0:44,:,:], axis=0)
	mean2 = np.nanmean(value[:,:,:], axis=0)

	return lat, lon, mean1, mean2

	
def import_cmip(param, cmip, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip/{0}'.format(cmip)
	arq   = '{0}/{1}_SA_Ayr_ensmean_{2}_historical_{3}_{4}_lonlat_mask.nc'.format(path, param, cmip, exp, date)	
		
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
	
	map = Basemap(projection='cyl', llcrnrlon=-90., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
	map.drawmapboundary(color='white')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	
	map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'world', drawbounds=True, color='black', linewidth=1.)
	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()
	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]

	map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul2', drawbounds=False)
	polys = polys + getattr(map, 'america_sul2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] 
	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)
	plt.gca().add_patch(patch) 
	
	return map, xx, yy
	
# Import regcm exps and obs database 
var = 'pre'

lat, lon, mean1_cru, mean2_cru = import_obs(var)
lat, lon, mean_cmip5 = import_cmip(u'pr', u'cmip5', 'r1i1p1', '1961-2005')
lat, lon, mean_cmip6 = import_cmip(u'pr', u'cmip6', 'r1i1p1f1_gn', '1961-2014')

# Compute and plot bias
bias_cmip5 = mean_cmip5 - mean1_cru
bias_cmip6 = mean_cmip6 - mean2_cru

fig = plt.figure(figsize=(7, 9))
levs1 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
levs2 = [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5]

ax = fig.add_subplot(5, 5, 1)  
plt.title(u'(a) CRU', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, mean1_cru, levels=levs1, latlon=True, cmap=cm.Blues, extend='max') 
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.25, 0.02, 0.47]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 5, 2)  
plt.title(u'(b) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both') 
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.96, 0.25, 0.02, 0.47]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 5, 3)  
plt.title(u'(c) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 4)  
plt.title(u'(d) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 5)  
plt.title(u'(e) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 6)  
plt.title(u'(f) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 7)  
plt.title(u'(g) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 8)  
plt.title(u'(h) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 9)  
plt.title(u'(i) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 10)  
plt.title(u'(j) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 11)  
plt.title(u'(k) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 12)  
plt.title(u'(l) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 13)  
plt.title(u'(m) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 14)  
plt.title(u'(n) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 15)  
plt.title(u'(o) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)

ax = fig.add_subplot(5, 5, 16)  
plt.title(u'(p) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 17)  
plt.title(u'(q) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 18)  
plt.title(u'(r) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 19)  
plt.title(u'(s) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 20)  
plt.title(u'(t) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 21)  
plt.title(u'(u) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 22)  
plt.title(u'(v) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 23)  
plt.title(u'(w) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 24)  
plt.title(u'(x) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

ax = fig.add_subplot(5, 5, 25)  
plt.title(u'(y) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG) 

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'pyplt_maps_bias_cmip6_1961-2014_{0}.png'.format(var)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



