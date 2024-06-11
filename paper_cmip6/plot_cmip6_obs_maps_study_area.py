# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "June 10, 2024"
__description__ = "This script plot study area"

import os
import numpy as np
import rockhound as rh
import matplotlib as mpl 
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap, maskoceans


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
	
	
# Import south america topography
grid = rh.fetch_etopo1(version="bedrock")
sa = grid.sel(latitude=slice(-60, 15), longitude=slice(-100, -20))
lat=sa.latitude.values
lon=sa.longitude.values

# Plot study area
fig = plt.figure(figsize=(6, 7))

ax = plt.subplot()
map, xx, yy = basemap(lat, lon)
plt_map = map.pcolormesh(xx, yy, sa.bedrock, cmap=cm.terrain, norm=mpl.colors.Normalize(vmin=0, vmax=3000))
plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')

x1,i1 = (-70,-5)
x2,i2 = (-70,5)
x3,i3 = (-45,5)
x4,i4 = (-45,-5)
poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

y1,j1 = (-70,-12.5)
y2,j2 = (-70,-5)
y3,j3 = (-45,-5)
y4,j4 = (-45,-12.5)
poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly2)

z1,k1 = (-45,-15)
z2,k2 = (-45,-2)
z3,k3 = (-34,-2)
z4,k4 = (-34,-15)
poly3 = Polygon([(z1,k1),(z2,k2),(z3,k3),(z4,k4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly3)

xx1,ii1 = (-55,-20)
xx2,ii2 = (-55,-10)
xx3,ii3 = (-45,-10)
xx4,ii4 = (-45,-20)
poly4 = Polygon([(xx1,ii1),(xx2,ii2),(xx3,ii3),(xx4,ii4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly4)

yy1,jj1 = (-60,-35)
yy2,jj2 = (-60,-20)
yy3,jj3 = (-45,-20)
yy4,jj4 = (-45,-35)
poly5 = Polygon([(yy1,jj1),(yy2,jj2),(yy3,jj3),(yy4,jj4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly5)

plt.text(-37, -53, u'\u25B2 \nN', color='black', fontweight='bold')
plt.text(-68, -4, u'NAMZ', color='gray', fontweight='bold')
plt.text(-55, -8, u'SAMZ', color='gray', fontweight='bold')
plt.text(-44, -6, u'NEB', color='gray', fontweight='bold')
plt.text(-54, -19, u'SAM', color='gray', fontweight='bold')
plt.text(-58, -23, u'LPB', color='gray', fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.15, 0.22, 0.02, 0.35]))
cbar.set_label('Topography (meters)', fontsize=8, fontweight='bold')
cbar.ax.tick_params(labelsize=8)

# Path out to save figure
path_out = '/afs/ictp.it/home/m/mda_silv/Documents/CMIP6/figs'
name_out = 'pyplt_maps_study_area.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()







