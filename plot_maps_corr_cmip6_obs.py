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

from dict_name_cmip6 import cmip6
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
	mean1 = value[0:45,:,:]
	mean2 = value[:,:,:]

	return lat, lon, mean1, mean2

	
def import_cmip(path, param, name, exp, date):
	
	arq   = '{0}/{1}_SA_Ayr_{2}_historical_{3}_{4}_lonlat_mask.nc'.format(path, param, name, exp, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mean = value[:,:,:]
	
	return lat, lon, mean
		       
  
def compute_corr(model, obs, **kwargs):

    """
	The input arrays must have the same dimensions
	:Param model: Numpy array with model data
	:Param obs: Numpy array with obs data
    Return: Pearson correlation
    """

    method = kwargs.pop('method', '3d')  # 1d ou 3d

    if method == '1d':

        corr, pvalue = st.pearsonr(model, obs)

        return corr

    elif method == '3d':

        timelen = float(obs.shape[0])

        obs_mean = np.nanmean(obs, axis=0)
        obs_std = np.nanstd(obs, axis=0)

        model_mean = np.nanmean(model, axis=0)
        model_std = np.nanstd(model, axis=0)

        x1 = (obs - obs_mean)/obs_std
        y1 = (model - model_mean)/model_std

        xymult = x1 * y1

        xysum = np.nansum(xymult, axis=0)

        corr = xysum/timelen

        return corr

    else:

        print('ClimateStats.compute_pearson function')
        print('Input data error')
        exit(1)
        
              
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
	
# Import cmip models and obs database 
var = 'tmp'
var_mdl = 'tas'
cmip6_dt = '1961-2014'
cmip5_path = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip/cmip5'
cmip6_path = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip/cmip6'
corr_cmip6 = []

lat, lon, mean1_cru, mean2_cru = import_obs(var)
lat, lon, mean_cmip5 = import_cmip(cmip5_path, var_mdl, u'ensmean_cmip5', 'r1i1p1', '1961-2005')
corr_cmip5 = compute_corr(mean_cmip5, mean1_cru)

for i in range(1, 24):
	lat, lon, mean_cmip = import_cmip(cmip6_path, var_mdl, cmip6[i][0], cmip6[i][1], cmip6_dt)
	corr_cmip = compute_corr(mean_cmip, mean2_cru)
	corr_cmip6.append(corr_cmip)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(7.5, 9.5))

if var == 'pre':
	levs2 = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.5]
	color2 = cm.PiYG
else:
	levs2 = [-1, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1]
	color2 = cm.PiYG
	
ax = fig.add_subplot(5, 5, 1)  
plt.title(u'(a) {0}'.format(cmip6[1][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[0], levels=levs2, latlon=True, cmap=color2, extend='both') 
cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.92, 0.25, 0.02, 0.47]))
cbar.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 5, 2)  
plt.title(u'(b) {0}'.format(cmip6[2][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[1], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 3)  
plt.title(u'(c) {0}'.format(cmip6[3][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[2], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 4)  
plt.title(u'(d) {0}'.format(cmip6[4][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[3], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 5)  
plt.title(u'(e) {0}'.format(cmip6[5][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[4], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 6)  
plt.title(u'(f) {0}'.format(cmip6[6][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[5], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 7)  
plt.title(u'(g) {0}'.format(cmip6[7][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[6], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 8)  
plt.title(u'(h) {0}'.format(cmip6[8][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[7], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 9)  
plt.title(u'(i) {0}'.format(cmip6[9][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[8], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 10)  
plt.title(u'(j) {0}'.format(cmip6[10][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[9], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 11)  
plt.title(u'(k) {0}'.format(cmip6[11][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[10], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 12)  
plt.title(u'(l) {0}'.format(cmip6[12][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[11], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 13)  
plt.title(u'(m) {0}'.format(cmip6[13][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[12], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 14)  
plt.title(u'(n) {0}'.format(cmip6[14][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[13], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 15)  
plt.title(u'(o) {0}'.format(cmip6[15][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[14], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 16)  
plt.title(u'(p) {0}'.format(cmip6[16][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[15], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 17)  
plt.title(u'(q) {0}'.format(cmip6[17][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[16], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 18)  
plt.title(u'(r) {0}'.format(cmip6[18][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[17], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 19)  
plt.title(u'(s) {0}'.format(cmip6[19][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[18], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 20)  
plt.title(u'(t) {0}'.format(cmip6[20][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[19], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 21)  
plt.title(u'(u) {0}'.format(cmip6[21][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[20], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 22)  
plt.title(u'(v) {0}'.format(cmip6[22][0]), loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[21], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 23)  
plt.title(u'(w) CMIP6-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip6[22], levels=levs2, latlon=True, cmap=color2) 

ax = fig.add_subplot(5, 5, 24)  
plt.title(u'(x) CMIP5-MME', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, corr_cmip5, levels=levs2, latlon=True, cmap=color2) 

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'pyplt_maps_corr_cmip6_1961-2014_{0}.png'.format(var)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



