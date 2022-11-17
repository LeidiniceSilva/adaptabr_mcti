# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Nov 14, 2022"
__description__ = "This script plot annual climatology maps from cmip6"

import os
import conda
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib as mpl 
import matplotlib.cm as cm
import scipy.stats as stats
import matplotlib.pyplot as plt

# mpl.use('Agg')

conda_file_dir = conda.__file__
conda_dir = conda_file_dir.split('lib')[0]
proj_lib = os.path.join(os.path.join(conda_dir, 'share'), 'proj')
os.environ["PROJ_LIB"] = proj_lib

from os.path import expanduser
from matplotlib.path import Path
from matplotlib import colors as c
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import BoundaryNorm

def import_obs(param):
	
	path  = '/home/nice/Documents/dataset/obs/cmip6'
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
	
	path  = '/home/nice/Documents/dataset/gcm/cmip6/{0}'.format(cmip)
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
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documents/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)

	x0, x1 = plt.xlim()
	y0, y1 = plt.ylim()
	map_edges = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
	polys = [map_edges]

	map.readshapefile('/home/nice/Documents/github_projects/shp/shp_america_sul/america_sul', 'america_sul2', drawbounds=False)

	polys = polys + getattr(map, 'america_sul2')
	codes = [[Path.MOVETO] + [Path.LINETO for p in p[1:]] for p in polys] 
	polys_lin = [v for p in polys for v in p]
	codes_lin = [cdg for cdgs in codes for cdg in cdgs]
	path  = Path(polys_lin, codes_lin)
	patch = PathPatch(path, facecolor='white', lw=0)

	plt.gca().add_patch(patch) 

	plt.text(-44, -55, u'\u25B2 \nN', fontsize=8)
	
	return map, xx, yy
	
	
# ~ def plot_maps_clim(mean1_cru, mean2_cru, mean_cmip5, mean_cmip6, bias_cmip5, bias_cmip6):
		
	#~ fig = plt.figure()
	#~ levs1 = [0, 2, 4, 6, 8, 10, 12]
	#~ levs2 = [-3, -2, -1, 1, 2, 3]

	#~ ax = fig.add_subplot(231)
	#~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plt_map = map.contourf(xx, yy, mean1_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	#~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')

	#~ y1,j1 = map(-68,-12)
	#~ y2,j2 = map(-68,-3)
	#~ y3,j3 = map(-52,-3)
	#~ y4,j4 = map(-52,-12)

	#~ poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.0)
	#~ plt.gca().add_patch(poly2)

	#~ x1,i1 = map(-61,-37)
	#~ x2,i2 = map(-61,-27)
	#~ x3,i3 = map(-48,-27)
	#~ x4,i4 = map(-48,-37)

	#~ poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.0)
	#~ plt.gca().add_patch(poly1)
	
	#~ z1,k1 = map(-47,-15)
	#~ z2,k2 = map(-47,-2)
	#~ z3,k3 = map(-35,-2)
	#~ z4,k4 = map(-35,-15)

	#~ poly2 = Polygon([(z1,k1),(z2,k2),(z3,k3),(z4,k4)], facecolor='none', edgecolor='black', linewidth=1.0)
	#~ plt.gca().add_patch(poly2)
		
	#~ ax = fig.add_subplot(232)
	#~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plt_maps_clim = map.contourf(xx, yy, mean_cmip5, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=8) 
	#~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	#~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	#~ ax = fig.add_subplot(233)
	#~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plt_maps_clim = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
	#~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	#~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	#~ cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=8) 

	#~ ax = fig.add_subplot(234)
	#~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	#~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plt_map = map.contourf(xx, yy, mean2_cru, levels=levs1, latlon=True, cmap=cm.YlGnBu) 
	#~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	#~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	#~ ax = fig.add_subplot(235)
	#~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plt_maps_clim = map.contourf(xx, yy, mean_cmip6, levels=levs1, latlon=True, cmap=cm.YlGnBu, extend='max')
	#~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	#~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	#~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=8) 
		
	#~ ax = fig.add_subplot(236)
	#~ plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
	#~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	#~ map, xx, yy = basemap(lat, lon)
	#~ plt_maps_clim = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.BrBG, extend='both')
	#~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	#~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	#~ cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	#~ cbar.ax.tick_params(labelsize=8) 

	# ~ fig = plt.figure()
	# ~ levs1 = [18, 20, 24, 26, 28, 30, 32]
	# ~ levs2 = [-3, -2, -1, 1, 2, 3]

	# ~ ax = fig.add_subplot(231)
	# ~ plt.title(u'A)', loc='left', fontsize=8, fontweight='bold')
	# ~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
	# ~ map, xx, yy = basemap(lat, lon)
	# ~ plt_map = map.contourf(xx, yy, mean1_cru, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	# ~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	# ~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	# ~ ax = fig.add_subplot(232)
	# ~ plt.title(u'B)', loc='left', fontsize=8, fontweight='bold')
	# ~ map, xx, yy = basemap(lat, lon)
	# ~ plt_maps_clim = map.contourf(xx, yy, mean_cmip5, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='max')
	# ~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	# ~ cbar.ax.tick_params(labelsize=8) 
	# ~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	# ~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
		
	# ~ ax = fig.add_subplot(233)
	# ~ plt.title(u'C)', loc='left', fontsize=8, fontweight='bold')
	# ~ map, xx, yy = basemap(lat, lon)
	# ~ plt_maps_clim = map.contourf(xx, yy, bias_cmip5, levels=levs2, latlon=True, cmap=cm.bwr, extend='both')
	# ~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	# ~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	# ~ cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	# ~ cbar.ax.tick_params(labelsize=8) 

	# ~ ax = fig.add_subplot(234)
	# ~ plt.title(u'D)', loc='left', fontsize=8, fontweight='bold')
	# ~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	# ~ plt.ylabel(u'Latitude', fontsize=8, labelpad=25, fontweight='bold')
	# ~ map, xx, yy = basemap(lat, lon)
	# ~ plt_map = map.contourf(xx, yy, mean2_cru, levels=levs1, latlon=True, cmap=cm.YlOrRd) 
	# ~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	# ~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
	
	# ~ ax = fig.add_subplot(235)
	# ~ plt.title(u'E)', loc='left', fontsize=8, fontweight='bold')
	# ~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	# ~ map, xx, yy = basemap(lat, lon)
	# ~ plt_maps_clim = map.contourf(xx, yy, mean_cmip6, levels=levs1, latlon=True, cmap=cm.YlOrRd, extend='max')
	# ~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	# ~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	# ~ cbar = map.colorbar(ticks=levs1, drawedges=True, ax=ax)
	# ~ cbar.ax.tick_params(labelsize=8) 
		
	# ~ ax = fig.add_subplot(236)
	# ~ plt.title(u'F)', loc='left', fontsize=8, fontweight='bold')
	# ~ plt.xlabel(u'Longitude', fontsize=8, labelpad=15, fontweight='bold')
	# ~ map, xx, yy = basemap(lat, lon)
	# ~ plt_maps_clim = map.contourf(xx, yy, bias_cmip6, levels=levs2, latlon=True, cmap=cm.bwr, extend='both')
	# ~ map.drawmeridians(np.arange(-90.,-30.,20.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')
	# ~ map.drawparallels(np.arange(-60.,15.,15.), size=8, labels=[0,0,0,0], linewidth=0.5, color='black')
	# ~ cbar = map.colorbar(ticks=levs2, drawedges=True, ax=ax)
	# ~ cbar.ax.tick_params(labelsize=8) 
	
	# ~ return plt_maps_clim


# Import regcm exps and obs database 
lat, lon, mean1_cru, mean2_cru = import_obs(u'pre')
lat, lon, mean_cmip5 = import_cmip(u'pr', u'cmip5', 'r1i1p1', '1961-2005')
lat, lon, mean_cmip6 = import_cmip(u'pr', u'cmip6', 'r1i1p1f1_gn', '1961-2014')

lat, lon, mean1_cru, mean2_cru = import_obs(u'tmp')
lat, lon, mean_cmip5 = import_cmip(u'tas', u'cmip5', 'r1i1p1', '1961-2005')
lat, lon, mean_cmip6 = import_cmip(u'tas', u'cmip6', 'r1i1p1f1_gn', '1961-2014')

# Compute and plot bias
bias_cmip5 = mean_cmip5 - mean1_cru
bias_cmip6 = mean_cmip6 - mean2_cru

# Plot maps with the function
plt_map = plot_maps_clim(mean1_cru, mean2_cru, mean_cmip5, mean_cmip6, bias_cmip5, bias_cmip6)

# Path out to save figure
path_out = '/home/nice/Documents/ufrn/papers/paper_cmip6/figs'
name_out = 'pyplt_maps_clim_cmip6_1961-2014_tas.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



