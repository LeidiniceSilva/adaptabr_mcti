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
fig = plt.figure(figsize=(6, 6))
x = np.arange(9)

ax = fig.add_subplot(3, 1, 1)
plot = plt.bar(x + 0.00, scopus_cmip6, color='orange', edgecolor='black', linewidth=1, width = 0.25, label=u'Scopus')
plot = plt.bar(x + 0.25, springer_cmip6, color='blue', edgecolor='black', linewidth=1, width = 0.25, label=u'Springer Link')
plot = plt.bar(x + 0.50, periodicos_cmip6, color='gray', edgecolor='black', linewidth=1, width = 0.25, label=u'Periódicos CAPES')
plt.title(u'(a) Resultados de publicações por ano (CMIP6)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(x + .25, ('2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.yticks(np.arange(0, 330, 30), fontsize=8)
plt.legend(loc=2, handlelength=0.75, handleheight=0.75, shadow=True, ncol=1, prop=FontProperties(size=8))

ax = fig.add_subplot(3, 1, 2)
plot = plt.bar(x + 0.00, scopus_cmip5, color='orange', edgecolor='black', linewidth=1, width = 0.25, label=u'Scopus')
plot = plt.bar(x + 0.25, springer_cmip5, color='blue', edgecolor='black', linewidth=1, width = 0.25, label=u'Springer Link')
plot = plt.bar(x + 0.50, periodicos_cmip5, color='gray', edgecolor='black', linewidth=1, width = 0.25, label=u'Periódicos CAPES')
plt.title(u'(b) Resultados de publicações por ano (CMIP5)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(x + .25, ('2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'), fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.yticks(np.arange(0, 550, 50), fontsize=8)
plt.ylabel('Número de publicações', fontweight='bold', fontsize=8)

ax = fig.add_subplot(3, 1, 3)
plot = plt.bar(x + 0.00, scopus_cordex, color='orange', edgecolor='black', linewidth=1, width = 0.25, label=u'Scopus')
plot = plt.bar(x + 0.25, springer_cordex, color='blue', edgecolor='black', linewidth=1, width = 0.25, label=u'Springer Link')
plot = plt.bar(x + 0.50, periodicos_cordex, color='gray', edgecolor='black', linewidth=1, width = 0.25, label=u'Periódicos CAPES')
plt.title(u'(c) Resultados de publicações por ano (CORDEX)', loc='left', fontweight='bold', fontsize=8)
plt.xticks(x + .25, ('2014', '2015', '2016', '2017', '2018', '2019', '2020', '2021', '2022'), fontsize=8)
plt.yticks(np.arange(0, 220, 20), fontsize=8)
plt.xlabel('Anos', fontweight='bold', fontsize=8)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
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
fig = plt.figure(figsize=(6, 8))
y = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

ax = fig.add_subplot(3, 1, 1)
plt.barh(y, cmip6, color='orange', edgecolor='black', linewidth=1)
for i, v in enumerate(cmip6):
    plt.text(v + 3, i + .50, str(v), color='black', fontsize=8)
plt.title(u'(a) Jornais que mais publicaram: CMIP6 (Top 10)', loc='left', fontweight='bold', fontsize=8)
plt.yticks(y, cmip6_j, fontsize=8)
plt.xticks(np.arange(0, 70, 10), fontsize=8)
plt.ylabel('Jornais', fontweight='bold', fontsize=8)

ax = fig.add_subplot(3, 1, 2)
plt.barh(y, cmip5, color='orange', edgecolor='black', linewidth=1)
for i, v in enumerate(cmip5):
    plt.text(v + 3, i + .50, str(v), color='black', fontsize=8)
plt.title(u'(b) Jornais que mais publicaram: CMIP5 (Top 10)', loc='left', fontweight='bold', fontsize=8)
plt.yticks(y, cmip5_j, fontsize=8)
plt.xticks(np.arange(0, 1400, 200), fontsize=8)
plt.ylabel('Jornais', fontweight='bold', fontsize=8)

ax = fig.add_subplot(3, 1, 3)
plt.barh(y, cordex, color='orange', edgecolor='black', linewidth=1)
for i, v in enumerate(cordex):
    plt.text(v + 3, i + .50, str(v), color='black', fontsize=8)
plt.title(u'(c) Jornais que mais publicaram: CORDEX (Top 10)', loc='left', fontweight='bold', fontsize=8)
plt.yticks(y, cordex_j, fontsize=8)
plt.xticks(np.arange(0, 100, 10), fontsize=8)
plt.xlabel('Número de publicações', fontweight='bold', fontsize=8)
plt.ylabel('Jornais', fontweight='bold', fontsize=8)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'figure_2.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()

print('Plot figure')
# Plot figure 
fig = plt.figure()

sizes = [3035, 197, 147, 38, 33]
labels = ['Articles', 'Review', 'Book chapter', 'Short communications', 'Encyclopedia']
colors = ['royalblue', 'yellow', 'limegreen', 'red', 'silver']
pie = plt.pie(sizes, wedgeprops={'linewidth': 1, 'edgecolor': 'black'}, colors=colors, startangle=90, shadow=True)
plt.title(u'(a) Tipos de documentos mais publicados (Top 5)', loc='left', fontweight='bold', fontsize=8)
labels = [f'{l} ({s:0.0f})' for l, s in zip(labels, sizes)]
plt.legend(bbox_to_anchor=(0.85, 1), loc='upper left', labels=labels, prop=FontProperties(size=8))

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'figure_3.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit

print('Plot figure')
# Plot figure 
fig = plt.figure(figsize=(8, 4))

ax = fig.add_subplot(1, 3, 1)
sizes = [429, 33, 29, 7, 3]
labels = ['Articles', 'Review', 'Book chapter', 'Short communications', 'Encyclopedia']
colors = ['royalblue', 'limegreen', 'red', 'yellow', 'silver']
pie = plt.pie(sizes, wedgeprops={'linewidth': 1, 'edgecolor': 'black'}, colors=colors, startangle=90, shadow=True)
plt.title(u'(a) Documentos mais publicados: \n CMIP6 (Top 5)', loc='left', fontweight='bold', fontsize=8)
labels = [f'{l} ({s:0.0f})' for l, s in zip(labels, sizes)]
ax.legend(bbox_to_anchor=(0.50, -0.50), loc=8, labels=labels, prop=FontProperties(size=8))

ax = fig.add_subplot(1, 3, 2)
sizes = [1982, 128, 100, 23, 27]
labels = ['Articles', 'Review', 'Book chapter', 'Short communications', 'Encyclopedia']
colors = ['royalblue', 'limegreen', 'red', 'yellow', 'silver']
pie = plt.pie(sizes, wedgeprops={'linewidth': 1, 'edgecolor': 'black'}, colors=colors, startangle=90, shadow=True)
plt.title(u'(b) Documentos mais publicados: \n CMIP5 (Top 5)', loc='left', fontweight='bold', fontsize=8)
labels = [f'{l} ({s:0.0f})' for l, s in zip(labels, sizes)]
ax.legend(bbox_to_anchor=(0.50, -0.50), loc=8, labels=labels, prop=FontProperties(size=8))

ax = fig.add_subplot(1, 3, 3)
sizes = [624, 36, 18, 8, 3]
labels = ['Articles', 'Review', 'Book chapter', 'Short communications', 'Encyclopedia']
colors = ['royalblue', 'limegreen', 'red', 'yellow', 'silver']
pie = plt.pie(sizes, wedgeprops={'linewidth': 1, 'edgecolor': 'black'}, colors=colors, startangle=90, shadow=True)
plt.title(u'(c) Documentos mais publicados: \n CORDEX (Top 5)', loc='left', fontweight='bold', fontsize=8)
labels = [f'{l} ({s:0.0f})' for l, s in zip(labels, sizes)]
ax.legend(bbox_to_anchor=(0.50, -0.50), loc=8, labels=labels, prop=FontProperties(size=8))

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'figure_3.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit

print('Plot figure')
# Plot figure 

fig = plt.figure(figsize=(6, 10))

ax = fig.add_subplot(3, 1, 1)
sizes = [291, 270, 126, 62, 43]
labels = ['Environmental \n Science (291)', 'Earth and \n Planetary \n Sciences (270)',
'Agricultural and \n Biological Sciences (126)', 'Energy (62)', 'Social Sciences (43)']
wedges, texts = ax.pie(sizes, wedgeprops=dict(width=0.5), startangle=-40)
bbox_props = dict(boxstyle='square,pad=0.3', fc='w', ec='k', lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
ax.set_title(u'(a) Áreas temáticas: CMIP6 (Top 5)', loc='left', fontweight='bold')
for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = f"angle,angleA=0,angleB={ang}"
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y), horizontalalignment=horizontalalignment, **kw)

ax = fig.add_subplot(3, 1, 2)
sizes = [1289, 1151, 789, 232, 220]
labels = ['Environmental \n Science (1289)', 'Earth and \n Planetary \n Sciences (1151)', 
'Agricultural and \n Biological Sciences (789)', 'Social Science (232)', 'Energy (220)']
wedges, texts = ax.pie(sizes, wedgeprops=dict(width=0.5), startangle=-40)
bbox_props = dict(boxstyle='square,pad=0.3', fc='w', ec='k', lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
ax.set_title(u'(b) Áreas temáticas: CMIP5 (Top 5)', loc='left', fontweight='bold')
for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = f"angle,angleA=0,angleB={ang}"
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y), horizontalalignment=horizontalalignment, **kw)

ax = fig.add_subplot(3, 1, 3)
sizes = [420, 328, 208, 110, 64]
labels = ['Environmental \n Science (420)', 'Earth and \n Planetary \n Sciences (328)', 
'Agricultural and \n Biological Sciences (208)', 'Energy (110)', 'Engineering (64)']
wedges, texts = ax.pie(sizes, wedgeprops=dict(width=0.5), startangle=-40)
bbox_props = dict(boxstyle='square,pad=0.3', fc='w', ec='k', lw=0.72)
kw = dict(arrowprops=dict(arrowstyle="-"), bbox=bbox_props, zorder=0, va="center")
ax.set_title(u'(c) Áreas temáticas: CORDEX (Top 5)', loc='left', fontweight='bold')
for i, p in enumerate(wedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = f"angle,angleA=0,angleB={ang}"
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax.annotate(labels[i], xy=(x, y), xytext=(1.35*np.sign(x), 1.4*y), horizontalalignment=horizontalalignment, **kw)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs'
name_out = 'figure_4.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit




