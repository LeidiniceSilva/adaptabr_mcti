# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot cri maps of cmip6 models"

import os
import cmocean
import numpy as np
import rockhound as rh
import matplotlib.pyplot as plt

from matplotlib.path import Path
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from matplotlib.patches import PathPatch

grid = rh.fetch_etopo1(version="bedrock")
sa = grid.sel(latitude=slice(-60, 15), longitude=slice(-85, -30))

plt.figure(figsize=(6, 7))
ax = plt.subplot(111)

sa.bedrock.plot.pcolormesh(cmap=cmocean.cm.topo, cbar_kwargs=dict(pad=0.01, aspect=30, label='Topografia (m)'), ax=ax)
ax.set_xlabel(u'Longitude', fontweight='bold')
ax.set_ylabel(u'Latitude', fontweight='bold')

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
plt.text(-52, -34, u'LPB', color='gray', fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/project/figs/figs_report-II'
name_out = 'pyplt_maps_study_area.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

