# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script plot skill of cmip6 correct models"

import os
import numpy as np
import xarray as xr
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = "/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias"


def import_observed(var_name, target_date):

    """ Import observed data.
    :param: model_name: str (BR-DWGD)
    :param: target_date: datetime
    :return: observed data (Prec, Tasmax and Tasmin)
    :return: flag that it refers to data control
    :rtype: 3D array
    """
    
    arq = xr.open_dataset('{0}/obs/'.format(dataset_dir) + '{0}_{1}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc'.format(var_name, target_date))
    data = arq[var_name]
    time = data.sel(time=slice('1986-01-01','2005-12-31'))
    var = time.groupby('time.year').mean('time')
    lat = var.lat
    lon = var.lon
    mean = np.nanmean(var.values, axis=0)
    
    return lat, lon, mean
   
   
def import_simulated(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	arq = xr.open_dataset('{0}/cmip6/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	time = data.sel(time=slice('1986-01-01','2005-12-31'))
	var = time.groupby('time.year').mean('time')
	lat = var.lat
	lon = var.lon
	mean = np.nanmean(var.values, axis=0)
	
	return lat, lon, mean


def import_simulated_correct(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	arq = xr.open_dataset('{0}/cmip6_correct/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_correct.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
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
	map.drawmeridians(np.arange(-76,-32,8), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	map.drawparallels(np.arange(-36,7,8), labels=[1,0,0,0], size=6, linewidth=0.5, color='black')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)

	# Import shapefile 
	map.readshapefile('/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shp/shp_world/world', 'world', drawbounds=True, color='white')
	map.readshapefile('/afs/ictp.it/home/m/mda_silv/Documents/github_projects/shp/lim_unid_fed/lim_unid_fed', 'lim_unid_fed', drawbounds=True, color='black')
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
	

# Best models list
best_models = [7, 9, 13, 15, 17]

# Variable dictionary
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}

experiment = 'historical'
dt = '19860101-20051231'

for i in best_models:
	for j in range(1, 4):	
		var_obs = var_dict[j][0]
		var_cmip6 = var_dict[j][1]
	
		print(cmip6[i][0])
		print(var_cmip6)

		# Import cmip models and obs database 
		lat, lon, obs = import_observed(var_obs, dt)
		lat, lon, sim = import_simulated(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
		lat, lon, sim_correct = import_simulated_correct(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)

		bias_sim = sim - obs
		bias_sim_correct = sim_correct - obs

		# Plot cmip models and obs database 
		fig = plt.figure(figsize=(8, 4))
		font_size = 7

		if var_cmip6 == 'pr':
			levs0 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
			color0 = cm.Blues
			levs = [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
			color = cm.BrBG
			levs1 = [-1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2]
			legend = 'Viés de precipitação (mm d⁻¹)'
		elif var_cmip6 == 'tasmax':
			levs0 = [26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
			color0 = cm.Reds
			levs = [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
			color = cm.bwr
			levs1 = [-1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2]
			legend = 'Viés de temperatura maxima (°C)'
		else:
			levs0 = [16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26]
			color0 = cm.Reds
			levs = [-3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3]
			color = cm.bwr
			levs1 = [-1.2, -1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1, 1.2]
			legend = 'Viés de temperatura mínima (°C)'

		ax = fig.add_subplot(2, 3, 1)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contourf(xx, yy, obs, levels=levs0, latlon=True, cmap=color0, extend='max') 
		plt.title(u'(a) BR-DWGD', loc='left', fontsize=font_size, fontweight='bold')
		plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')
		plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
		cbar = map.colorbar(plt_map, shrink=0.8, pad=0.1)
		cbar.ax.tick_params(labelsize=font_size) 

		ax = fig.add_subplot(2, 3, 2)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contourf(xx, yy, sim, levels=levs0, latlon=True, cmap=color0, extend='max') 
		plt.title(u'(b) {0}'.format(cmip6[i][0]), loc='left', fontsize=font_size, fontweight='bold')
		cbar = map.colorbar(plt_map, shrink=0.8, pad=0.1)
		cbar.ax.tick_params(labelsize=font_size) 

		ax = fig.add_subplot(2, 3, 3)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contourf(xx, yy, sim_correct, levels=levs0, latlon=True, cmap=color0, extend='max') 
		plt.title(u'(c) {0} correct'.format(cmip6[i][0]), loc='left', fontsize=font_size, fontweight='bold')
		cbar = map.colorbar(plt_map, shrink=0.8, pad=0.1)
		cbar.set_label('Climatologia', fontsize=font_size, fontweight='bold')
		cbar.ax.tick_params(labelsize=font_size)  

		ax = fig.add_subplot(2, 3, 5)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contourf(xx, yy, bias_sim, levels=levs, latlon=True, cmap=color, extend='both') 
		plt.title(u'(d) Bias {0}'.format(cmip6[i][0]), loc='left', fontsize=font_size, fontweight='bold')
		plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')
		plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
		cbar = map.colorbar(plt_map, shrink=0.8, pad=0.1)
		cbar.ax.tick_params(labelsize=font_size) 

		ax = fig.add_subplot(2, 3, 6)  
		map, xx, yy = basemap(lat, lon)
		plt_map = map.contourf(xx, yy, bias_sim_correct, levels=levs1, latlon=True, cmap=color, extend='both') 
		plt.title(u'(e) Bias {0} correct'.format(cmip6[i][0]), loc='left', fontsize=font_size, fontweight='bold')
		plt.xlabel(u'Longitude', labelpad=10, fontsize=font_size, fontweight='bold')
		cbar = map.colorbar(plt_map, shrink=0.8, pad=0.1)
		cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
		cbar.ax.tick_params(labelsize=font_size)  

		# Path out to save figure
		path_out = '/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/figs/correct_bias'
		name_out = 'pyplt_maps_clim_correct_bias_cmip6_{0}_{1}_{2}.png'.format(cmip6[i][0], var_cmip6, dt)
		plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
		plt.close('all')
		plt.cla()
exit()
