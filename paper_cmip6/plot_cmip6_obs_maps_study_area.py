# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 14, 2023"
__description__ = "This script plot study area"

import os
import cmocean
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

map, xx, yy = basemap(lat, lon)
# ~ sa_mask = maskoceans(xx, yy, sa.bedrock, resolution='c')
plt_map = map.pcolormesh(xx, yy, sa.bedrock, cmap=cm.terrain, norm=mpl.colors.Normalize(vmin=0, vmax=3000)) 
plt.xlabel(u'Longitude', labelpad=15, fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=25, fontsize=8, fontweight='bold')

# NWS
a1,b1 = (-83,-15)
a2,b2 = (-83,12)
a3,b3 = (-70,12)
a4,b4 = (-70,-15)
poly1 = Polygon([(a1,b1),(a2,b2),(a3,b3),(a4,b4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly1)

# NSA
c1,d1 = (-70,-7)
c2,d2 = (-70,12)
c3,d3 = (-50,12)
c4,d4 = (-50,-7)
poly2 = Polygon([(c1,d1),(c2,d2),(c3,d3),(c4,d4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly2)

# SAM
f1,g1 = (-70,-20)
f2,g2 = (-70,-7)
f3,g3 = (-50,-7)
f4,g4 = (-50,-20)
poly3 = Polygon([(f1,g1),(f2,g2),(f3,g3),(f4,g4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly3)

# SES
j1,k1 = (-68,-40)
j2,k2 = (-68,-20)
j3,k3 = (-40,-20)
j4,k4 = (-40,-40)
poly4 = Polygon([(j1,k1),(j2,k2),(j3,k3),(j4,k4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly4)

# NES
l1,m1 = (-50,-20)
l2,m2 = (-50,0)
l3,m3 = (-34,0)
l4,m4 = (-34,-20)
poly5 = Polygon([(l1,m1),(l2,m2),(l3,m3),(l4,m4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly5)

# SWS 
n1,p1 = (-76,-42)
n2,p2 = (-76,-15)
n3,p3 = (-68,-15)
n4,p4 = (-68,-42)
poly6 = Polygon([(n1,p1),(n2,p2),(n3,p3),(n4,p4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly6)

# SSA
q1,r1 = (-78,-56)
q2,r2 = (-78,-40)
q3,r3 = (-56,-40)
q4,r4 = (-56,-56)
poly7 = Polygon([(q1,r1),(q2,r2),(q3,r3),(q4,r4)], facecolor='none', edgecolor='black', linewidth=1.)
plt.gca().add_patch(poly7)

plt.text(-28, -54, u'\u25B2 \nN', color='black', fontsize=8, fontweight='bold')
plt.text(-79, -7, u'NWS', color='black', fontsize=8, fontweight='bold')
plt.text(-67, -3, u'NSA', color='black', fontsize=8, fontweight='bold')
plt.text(-58, -12, u'SAM', color='black', fontsize=8, fontweight='bold')
plt.text(-51, -24, u'SES', color='black', fontsize=8, fontweight='bold')
plt.text(-45, -9, u'NES', color='black', fontsize=8, fontweight='bold')
plt.text(-75, -24, u'SWS', color='black', fontsize=8, fontweight='bold')
plt.text(-71, -44, u'SSA', color='black', fontsize=8, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.15, 0.20, 0.02, 0.35]))
cbar.set_label('Topography (meters)', fontsize=8, fontweight='bold')
# ~ cbar.ax.tick_params(labelsize=8)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/figs'
name_out = 'pyplt_maps_study_area.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()







