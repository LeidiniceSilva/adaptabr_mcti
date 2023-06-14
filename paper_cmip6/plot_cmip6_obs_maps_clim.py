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
from mpl_toolkits.basemap import Basemap
	
# Import south america topography
grid = rh.fetch_etopo1(version="bedrock")
sa = grid.sel(latitude=slice(-60, 15), longitude=slice(-100, -20))
lat=sa.latitude.values
lon=sa.longitude.values

# Plot study area
fig = plt.figure()
ax = plt.subplot(111)
color = cmocean.cm.topo
norm = mpl.colors.Normalize(vmin=-8000, vmax=8000)
level = (-8000, -6000, -4000, -2000, 0, 2000, 4000, 6000, 8000)

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
map.drawmeridians(np.arange(-100., -20., 10.), size=6, labels=[0,0,0,1], linewidth=0.4, color='black')
map.drawparallels(np.arange(-60., 15., 10.), size=6, labels=[1,0,0,0], linewidth=0.4, color='black')
lons, lats = np.meshgrid(new_lon, new_lat)
xx, yy = map(lons,lats)

plt_map = map.pcolormesh(xx, yy, sa.bedrock, cmap=color, latlon=True) 
plt.xlabel(u'Longitude', labelpad=10, fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=8, fontweight='bold')

cbar = plt.colorbar(plt_map, cax=fig.add_axes([0.84, 0.28, 0.02, 0.43]))
cbar.set_label('Topography (meters)', fontsize=8, fontweight='bold')
cbar.ax.tick_params(labelsize=8)

plt.text(-37, -53, u'\u25B2 \nN', color='black', fontweight='bold')
plt.text(-68, -4, u'NAMZ', color='gray', fontweight='bold')
plt.text(-55, -8, u'SAMZ', color='gray', fontweight='bold')
plt.text(-44, -6, u'NEB', color='gray', fontweight='bold')
plt.text(-54, -19, u'SAM', color='gray', fontweight='bold')
plt.text(-58, -23, u'LPB', color='gray', fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/paper_cmip6/figs'
name_out = 'pyplt_maps_study_area.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()







