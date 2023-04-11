# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot ivs of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_statistical_metrics import compute_ivs


def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_ANN_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
		
	return obs_data

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_ANN_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	return mdl_data
	              

def sort_list(data_list):
	
	li = []
	for i in range(len(data_list)):
		li.append([data_list[i], i])
	  
	li.sort()
	sort_index = []
	for x in li:
		sort_index.append(x[1])
	
	return sort_index
	
	           
# Import cmip models and obs database 
var_obs = 'pre'
var_cmip6 = 'pr'
dt = '1980-2014'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)

ivs_namz_cmip6 = []
ivs_samz_cmip6 = []
ivs_neb_cmip6 = []
ivs_sam_cmip6 = []
ivs_lpb_cmip6 = []
legend = []

for i in range(1, 19):

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	ivs_namz_cmip6.append(compute_ivs(clim_namz_obs, clim_namz_cmip))
	
	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	ivs_samz_cmip6.append(compute_ivs(clim_samz_obs, clim_samz_cmip))
	
	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	ivs_neb_cmip6.append(compute_ivs(clim_neb_obs, clim_neb_cmip))
		
	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	ivs_sam_cmip6.append(compute_ivs(clim_sam_obs, clim_sam_cmip))
		
	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	ivs_lpb_cmip6.append(compute_ivs(clim_lpb_obs, clim_lpb_cmip))
	
	legend.append(cmip6[i][0])

ivs_cmip6 = np.array([ivs_lpb_cmip6,ivs_sam_cmip6,ivs_neb_cmip6,ivs_samz_cmip6,ivs_namz_cmip6])

sort_list_namz = sort_list(ivs_namz_cmip6)
model_list_namz = []
value_list_namz = []
for i in sort_list_namz:
	model_list_namz.append(cmip6[i+1][0])
	value_list_namz.append(ivs_namz_cmip6[i])

sort_list_samz = sort_list(ivs_samz_cmip6)
model_list_samz = []
value_list_samz = []
for ii in sort_list_samz:
	model_list_samz.append(cmip6[ii+1][0])
	value_list_samz.append(ivs_samz_cmip6[ii])

sort_list_neb = sort_list(ivs_neb_cmip6)
model_list_neb = []
value_list_neb = []
for iii in sort_list_neb:
	model_list_neb.append(cmip6[iii+1][0])
	value_list_neb.append(ivs_neb_cmip6[iii])

sort_list_sam = sort_list(ivs_sam_cmip6)
model_list_sam = []
value_list_sam = []
for iv in sort_list_sam:
	model_list_sam.append(cmip6[iv+1][0])
	value_list_sam.append(ivs_sam_cmip6[iv])

sort_list_lpb = sort_list(ivs_lpb_cmip6)
model_list_lpb = []
value_list_lpb = []
for v in sort_list_lpb:
	model_list_lpb.append(cmip6[v+1][0])
	value_list_lpb.append(ivs_lpb_cmip6[v])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(6, 10))

if var_cmip6 == 'pr':
	color = 'blue'
else:
	color = 'red'
	
ax = fig.add_subplot(5, 1, 1)  
ax.barh(model_list_namz, value_list_namz, color=color, edgecolor='white')
plt.title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
ax.invert_yaxis()
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
             
ax = fig.add_subplot(5, 1, 2)  
ax.barh(model_list_samz, value_list_samz, color=color, edgecolor='white')
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
ax.invert_yaxis()
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
        
ax = fig.add_subplot(5, 1, 3)  
ax.barh(model_list_neb, value_list_neb, color=color, edgecolor='white')
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
ax.invert_yaxis()
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
        
ax = fig.add_subplot(5, 1, 4)  
ax.barh(model_list_sam, value_list_sam, color=color, edgecolor='white')
plt.title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
ax.invert_yaxis()
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
        
ax = fig.add_subplot(5, 1, 5)  
ax.barh(model_list_lpb, value_list_lpb, color=color, edgecolor='white')
plt.title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
plt.yticks(fontsize=7)
plt.xticks(fontsize=7)
ax.xaxis.set_tick_params(pad=-5)
for s in ['top', 'bottom', 'left', 'right']:
    ax.spines[s].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')
ax.grid(b=True, color ='gray', linestyle='--', linewidth=0.5, alpha = 0.2)
ax.invert_yaxis()
for i in ax.patches:
    plt.text(i.get_width()+0.2, i.get_y()+0.5, str(round((i.get_width()), 2)), fontsize=7, fontweight='bold', color='gray')
        
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_rank_ivs_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 3))
norm = colors.BoundaryNorm(boundaries=np.arange(0, 5.5, 0.5), ncolors=256)
xlabels = legend
ylabels = [u'LPB', u'SAM', u'NEB', u'SAMZ', u'NAMZ']

if var_cmip6 == 'pr':
	color = cm.Blues
else:
	color = cm.Reds

ax = fig.add_subplot(1, 1, 1)  
pcm = ax.pcolormesh(ivs_cmip6, edgecolors ='white', linewidths = 2., norm=norm, cmap=color)
ax.set_title(u'(a) IVS', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(ivs_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(ivs_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
for y in range(ivs_cmip6.shape[0]):
    for x in range(ivs_cmip6.shape[1]):
        ax.text(x + 0.5, y + 0.5, '%.2f' % ivs_cmip6[y, x],
                 ha="center", va="center", color='k', size=8)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_portrait_ivs_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






