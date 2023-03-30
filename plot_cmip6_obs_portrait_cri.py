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
from matplotlib.patches import PathPatch
from mpl_toolkits.basemap import Basemap
from dict_cmip6_models_name import cmip6

# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 10))

ax = fig.add_subplot(1, 2, 1)  
map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
path = '/home/nice/Documentos/github_projects/shp'
map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)
map.fillcontinents(color='beige', lake_color='wheat')
map.drawmapboundary(fill_color='skyblue')
map.drawparallels(np.arange(-60., 30., 15.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
map.drawmeridians(np.arange(-85., -15., 15.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')

ax = fig.add_subplot(1, 2, 2)  
map = Basemap(projection='cyl', llcrnrlon=-85., llcrnrlat=-60., urcrnrlon=-30.,urcrnrlat=15., resolution='c')
path = '/home/nice/Documentos/github_projects/shp'
map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=.8)
map.fillcontinents(color='beige', lake_color='wheat')
map.drawmapboundary(fill_color='skyblue')
map.drawparallels(np.arange(-60., 30., 15.), size=8, labels=[1,0,0,0], linewidth=0.5, color='black')
map.drawmeridians(np.arange(-85., -15., 15.), size=8, labels=[0,0,0,1], linewidth=0.5, color='black')

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_maps_cri_cmip6_pr_tas_1980-2014.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()



