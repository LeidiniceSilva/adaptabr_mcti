# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot cri maps of cmip6 models"

import os
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import Polygon
from dict_cmip6_models_name import cmip6
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch

# Plot cmip models and obs database 
fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)  
map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
path = '/home/nice/Documentos/github_projects/shp'
map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='gray', linewidth=.8)
map.fillcontinents(color='beige', lake_color='wheat')
map.drawmapboundary(fill_color='skyblue')
map.drawparallels(np.arange(-60., 30., 15.), labels=[1,0,0,0], linewidth=0.5, color='black')
map.drawmeridians(np.arange(-85., -30., 15.), labels=[0,0,0,1], linewidth=0.5, color='black')

x1,i1 = map(-70,-5)
x2,i2 = map(-70,5)
x3,i3 = map(-45,5)
x4,i4 = map(-45,-5)
poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='gray', linewidth=1.)
plt.gca().add_patch(poly1)

y1,j1 = map(-70,-12.5)
y2,j2 = map(-70,-5)
y3,j3 = map(-45,-5)
y4,j4 = map(-45,-12.5)
poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='gray', linewidth=1.)
plt.gca().add_patch(poly2)

z1,k1 = map(-45,-15)
z2,k2 = map(-45,-2)
z3,k3 = map(-34,-2)
z4,k4 = map(-34,-15)
poly3 = Polygon([(z1,k1),(z2,k2),(z3,k3),(z4,k4)], facecolor='none', edgecolor='gray', linewidth=1.)
plt.gca().add_patch(poly3)

xx1,ii1 = map(-55,-20)
xx2,ii2 = map(-55,-10)
xx3,ii3 = map(-45,-10)
xx4,ii4 = map(-45,-20)
poly4 = Polygon([(xx1,ii1),(xx2,ii2),(xx3,ii3),(xx4,ii4)], facecolor='none', edgecolor='gray', linewidth=1.)
plt.gca().add_patch(poly4)

yy1,jj1 = map(-60,-35)
yy2,jj2 = map(-60,-20)
yy3,jj3 = map(-45,-20)
yy4,jj4 = map(-45,-35)
poly5 = Polygon([(yy1,jj1),(yy2,jj2),(yy3,jj3),(yy4,jj4)], facecolor='none', edgecolor='gray', linewidth=1.)
plt.gca().add_patch(poly5)

plt.text(-37, -53, u'\u25B2 \nN', color='gray', fontweight='bold')
plt.text(-68, -4, u'NAMZ', color='gray', fontweight='bold')
plt.text(-55, -8, u'SAMZ', color='gray', fontweight='bold')
plt.text(-44, -6, u'NEB', color='gray', fontweight='bold')
plt.text(-54, -19, u'SAM', color='gray', fontweight='bold')
plt.text(-52, -34, u'LPB', color='gray', fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_maps_cri_cmip6_pr_tas_1980-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



