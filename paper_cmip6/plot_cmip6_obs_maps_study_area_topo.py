# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "September 10, 2025"
__description__ = "This script plot map of genesis"

import os
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import matplotlib.patches as mpatches

from matplotlib.patches import Polygon
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# Open the file
file = '/afs/ictp.it/home/m/mda_silv/Downloads/orog_SAM-22_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-7_v0_fx.nc'
ds = xr.open_dataset(file)
alt = ds['orog']
lons = ds['lon']
lats = ds['lat']

# Create figure
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
lon1, lon2, lat1, lat2 = -90, -30, -60, 15

contour = ax.contourf(lons, lats, alt, levels=30, cmap='terrain', extend='max', transform=ccrs.PlateCarree())
ax.set_xlabel(u'Longitude', labelpad=15, fontweight='bold')
ax.set_ylabel(u'Latitude', labelpad=25, fontweight='bold')
ax.set_extent([lon1, lon2, lat1, lat2], crs=ccrs.PlateCarree())
ax.set_xticks(np.arange(lon1,lon2,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lat1,lat2,10), crs=ccrs.PlateCarree())
ax.xaxis.set_major_formatter(LongitudeFormatter())
ax.yaxis.set_major_formatter(LatitudeFormatter())
ax.grid(c='k', ls='--', alpha=0.5)       
ax.add_feature(cfeature.BORDERS, linestyle=":")
ax.coastlines(resolution='10m')
ax.add_feature(cfeature.OCEAN, facecolor='white', zorder=1)

cbar = plt.colorbar(contour, ax=ax, orientation='vertical', shrink=0.7, pad=0.05)
cbar.set_label('Topography (m)', fontweight='bold')

x1,i1 = (-70,-5)
x2,i2 = (-70,5)
x3,i3 = (-45,5)
x4,i4 = (-45,-5)
poly1 = Polygon([(x1,i1),(x2,i2),(x3,i3),(x4,i4)], facecolor='none', edgecolor='black', linewidth=1.5)
plt.gca().add_patch(poly1)

y1,j1 = (-70,-12.5)
y2,j2 = (-70,-5)
y3,j3 = (-45,-5)
y4,j4 = (-45,-12.5)
poly2 = Polygon([(y1,j1),(y2,j2),(y3,j3),(y4,j4)], facecolor='none', edgecolor='black', linewidth=1.5)
plt.gca().add_patch(poly2)

z1,k1 = (-45,-15)
z2,k2 = (-45,-2)
z3,k3 = (-34,-2)
z4,k4 = (-34,-15)
poly3 = Polygon([(z1,k1),(z2,k2),(z3,k3),(z4,k4)], facecolor='none', edgecolor='black', linewidth=1.5)
plt.gca().add_patch(poly3)

xx1,ii1 = (-55,-20)
xx2,ii2 = (-55,-10)
xx3,ii3 = (-45,-10)
xx4,ii4 = (-45,-20)
poly4 = Polygon([(xx1,ii1),(xx2,ii2),(xx3,ii3),(xx4,ii4)], facecolor='none', edgecolor='black', linewidth=1.5)
plt.gca().add_patch(poly4)

yy1,jj1 = (-60,-35)
yy2,jj2 = (-60,-20)
yy3,jj3 = (-45,-20)
yy4,jj4 = (-45,-35)
poly5 = Polygon([(yy1,jj1),(yy2,jj2),(yy3,jj3),(yy4,jj4)], facecolor='none', edgecolor='black', linewidth=1.5)
plt.gca().add_patch(poly5)

plt.text(-37, -53, u'\u25B2 \nN', color='white', fontweight='bold')
plt.text(-68, -4, u'NAMZ', color='white', fontweight='bold')
plt.text(-55, -8, u'SAMZ', color='white', fontweight='bold')
plt.text(-44, -6, u'NEB', color='white', fontweight='bold')
plt.text(-54, -19, u'SAM', color='white', fontweight='bold')
plt.text(-58, -23, u'LPB', color='white', fontweight='bold')

# Path out to save figure
path_out = '/afs/ictp.it/home/m/mda_silv/Downloads'
name_out = 'pyplt_maps_study_area_cordex_sam-22.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




