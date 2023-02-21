# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/13/2023"
__description__ = "This script plot figures to adapta_br-mcti report"

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties

print('Read data')
# Read data
scopus_cmip6 = np.array([0, 0, 0, 0, 0, 0, 51, 162, 299])
springer_cmip6 = np.array([0, 0, 0, 0, 0, 0, 89, 142, 270])
periodicos_cmip6 = np.array([0, 0, 0, 0, 0, 0, 57, 75, 106])

scopus_cmip5 = np.array([63, 142, 197, 181, 247, 275, 307, 444, 446])
springer_cmip5 = np.array([188, 210, 309, 350, 354, 471, 491, 499, 495])
periodicos_cmip5 = np.array([58, 101, 66, 83, 148, 111, 142, 90, 113])

scopus_cordex = np.array([0, 0, 39, 44, 80, 91, 107, 162, 154])
springer_cordex = np.array([28, 39, 68, 47, 110, 129, 140, 125, 166])
periodicos_cordex = np.array([14, 10, 20, 21, 57, 59, 48, 63, 58])

print('Plot figure')
# Plot figure 
fig = plt.figure()
x = np.arange(9)

ax = fig.add_subplot(3, 1, 1)
plot = plt.bar(x + 0.00, scopus_cmip6, color='orange', edgecolor='black', linewidth=1, width = 0.25, label=u'Scopus')
plot = plt.bar(x + 0.25, springer_cmip6, color='blue', edgecolor='black', linewidth=1, width = 0.25, label=u'Springer Link')
plot = plt.bar(x + 0.50, periodicos_cmip6, color='gray', edgecolor='black', linewidth=1, width = 0.25, label=u'Periódicos CAPES')
plt.title(u'(a)', loc='left', fontweight='bold')
plt.xticks(x + .25, ('2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.yticks(np.arange(0, 330, 30), fontsize=8)
plt.legend(loc=2, handlelength=0.75, handleheight=0.75, shadow=True, ncol=1, prop=FontProperties(size=8))

ax = fig.add_subplot(3, 1, 2)
plot = plt.bar(x + 0.00, scopus_cmip5, color='orange', edgecolor='black', linewidth=1, width = 0.25, label=u'Scopus')
plot = plt.bar(x + 0.25, springer_cmip5, color='blue', edgecolor='black', linewidth=1, width = 0.25, label=u'Springer Link')
plot = plt.bar(x + 0.50, periodicos_cmip5, color='gray', edgecolor='black', linewidth=1, width = 0.25, label=u'Periódicos CAPES')
plt.title(u'(b)', loc='left', fontweight='bold')
plt.xticks(x + .25, ('2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.yticks(np.arange(0, 550, 50), fontsize=8)
plt.ylabel('Número de publicações', fontweight='bold')

ax = fig.add_subplot(3, 1, 3)
plot = plt.bar(x + 0.00, scopus_cordex, color='orange', edgecolor='black', linewidth=1, width = 0.25, label=u'Scopus')
plot = plt.bar(x + 0.25, springer_cordex, color='blue', edgecolor='black', linewidth=1, width = 0.25, label=u'Springer Link')
plot = plt.bar(x + 0.50, periodicos_cordex, color='gray', edgecolor='black', linewidth=1, width = 0.25, label=u'Periódicos CAPES')
plt.title(u'(c)', loc='left', fontweight='bold')
plt.xticks(x + .25, ('2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'), fontsize=8)
plt.yticks(np.arange(0, 220, 20), fontsize=8)
plt.xlabel('Anos', fontweight='bold')

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'figure_1.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()


print('Read data')
# Read data
cmip6 = [48, 45, 33, 24, 19, 14, 13, 12, 10, 9]
cmip6_j = ['Science of the Total Environment', 'Atmospheric Research', 'Journal of Hydrology',
'Weather and Climate Extremes','Journal of Hydrology: Regional Studies','Advances in Climate Change Research',
'Journal of Cleaner Production','Global and Planetary Change','Ecological Indicators','Agricultural Water Management']

cmip5 = [1178, 1074, 773, 225, 213, 141, 46, 26, 12, 10]
cmip5_j = ['Environmental Science', 'Earth and Planetary Sciences', 'Agricultural and Biological Sciences',
'Social Sciences', 'Energy', 'Engineering', 'Computer Science', 'Economics, Econometrics and Finance',
'Business, Management and Accounting', 'Physics and Astronomy']

cordex = [87, 58, 42, 39, 27, 24, 18, 17, 15, 15]
cordex_j = ['Science of the Total Environment','Journal of Hydrology','Atmospheric Research',
'Journal of Hydrology: Regional Studies','Climate Services','Weather and Climate Extremes',
'Agricultural and Forest Meteorology','Renewable Energy','Applied Energy','Global and Planetary Change']

print('Plot figure')
# Plot figure 
fig = plt.figure()
y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

ax = fig.add_subplot(3, 1, 1)
plt.barh(y, cmip6, color='orange', edgecolor='black', linewidth=1)
for i, v in enumerate(cmip6):
    plt.text(v + 3, i + .50, str(v), color='black', fontsize=8)
plt.title(u'(a)', loc='left', fontweight='bold', fontsize=8)
plt.yticks(y, cmip6_j, fontsize=8)
plt.xticks(np.arange(0, 70, 10), fontsize=8)
plt.ylabel('Jornais', fontweight='bold', fontsize=8)

ax = fig.add_subplot(3, 1, 2)
plt.barh(y, cmip5, color='orange', edgecolor='black', linewidth=1)
for i, v in enumerate(cmip5):
    plt.text(v + 3, i + .50, str(v), color='black', fontsize=8)
plt.title(u'(b)', loc='left', fontweight='bold', fontsize=8)
plt.yticks(y, cmip5_j, fontsize=8)
plt.xticks(np.arange(0, 1400, 200), fontsize=8)
plt.ylabel('Jornais', fontweight='bold', fontsize=8)

ax = fig.add_subplot(3, 1, 3)
plt.barh(y, cordex, color='orange', edgecolor='black', linewidth=1)
for i, v in enumerate(cordex):
    plt.text(v + 3, i + .50, str(v), color='black', fontsize=8)
plt.title(u'(c)', loc='left', fontweight='bold', fontsize=8)
plt.yticks(y, cordex_j, fontsize=8)
plt.xticks(np.arange(0, 100, 10), fontsize=8)
plt.xlabel('Número de publicações', fontweight='bold', fontsize=8)
plt.ylabel('Jornais', fontweight='bold', fontsize=8)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'figure_2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit


