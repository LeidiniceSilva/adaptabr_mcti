# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot cri of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.pyplot as plt

from dict_cmip6_models_name import cmip6

# Import cmip models and obs database 
mdl1  = [14, 7, 2]
mdl2  = [11, 6, 11]
mdl3  = [15, 3, 13]
mdl4  = [9,  5, 14]
mdl5  = [13, 9, 7]
mdl6  = [12, 2, 4]
mdl7  = [2,  8, 1]
mdl8  = [10, 16, 6]
mdl9  = [4,  14, 9]
mdl10 = [17, 17, 12]
mdl11 = [7,  12, 15]
mdl12 = [8,  15, 17]
mdl13 = [3,  10, 8]
mdl14 = [6,  11, 5]
mdl15 = [5,  1, 3]
mdl16 = [16, 13, 16]
mdl17 = [1,  4, 10]
			 
# Plot cmip models and obs database 
fig = plt.figure(figsize=(10, 8))
xaxis = np.arange(1, 4, 1)
yaxis = np.arange(0, 19, 1)

ax = fig.add_subplot(1, 1, 1) 
annual_cycle = ax.plot(xaxis, mdl1, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl2, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl3, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl4, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl5, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl6, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl7, marker='o', ms=15, linewidth=3, color='red')
annual_cycle = ax.plot(xaxis, mdl8, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl9, marker='o', ms=15, linewidth=3, color='yellow')
annual_cycle = ax.plot(xaxis, mdl10, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl11, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl12, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl13, marker='o', ms=15, linewidth=3, color='green')
annual_cycle = ax.plot(xaxis, mdl14, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl15, marker='o', ms=15, linewidth=3, color='black')
annual_cycle = ax.plot(xaxis, mdl16, marker='o', ms=15, linewidth=3, color='lightgray')
annual_cycle = ax.plot(xaxis, mdl17, marker='o', ms=15, linewidth=3, color='blue')

plt.text(1.03, 16.8, 'KIOST-ESM')
plt.text(1.03, 15.8, 'NESM3')
plt.text(1.03, 14.8, 'CanESM5')
plt.text(1.03, 13.8, 'ACCESS-CM2')
plt.text(1.03, 12.8, 'CNRM-CM6-1')
plt.text(1.03, 11.8, 'CNRM-ESM2-1')
plt.text(1.03, 10.8, 'BCC-CSM2-MR')
plt.text(1.03, 9.8, 'INM-CM4-8')
plt.text(1.03, 8.8, 'CMCC-ESM2')
plt.text(1.03, 7.8, 'MIROC-ES2L')
plt.text(1.03, 6.8, 'MIROC6')
plt.text(1.03, 5.8, 'MPI-ESM1-2-LR')
plt.text(1.03, 4.8, 'MRI-ESM2-0')
plt.text(1.03, 3.8, 'INM-CM5-0')
plt.text(1.03, 2.8, 'MPI-ESM1-2-HR')
plt.text(1.03, 1.8, 'GFDL-ESM4')
plt.text(1.03, 0.8, 'NorESM2-MM')

plt.title(u'(a) BR', loc='left', fontweight='bold')
plt.xticks(xaxis, ('Precipitação', 'Temperatura máxima', 'Temperatura mínima'))
plt.yticks(yaxis, ('', 'Rank 1', 'Rank 2', 'Rank 3', 'Rank 4', 'Rank 5', 'Rank 6', 'Rank 7', 'Rank 8', 'Rank 9', 'Rank 10', 'Rank 11', 'Rank 12', 'Rank 13', 'Rank 14', 'Rank 15', 'Rank 16', 'Rank 17', ''))
plt.ylim(0, 18)
ax.invert_yaxis()
	
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_rank_cmip6.png'
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






