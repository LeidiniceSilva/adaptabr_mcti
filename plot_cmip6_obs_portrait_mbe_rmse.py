# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot mbe and rmse of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6
from comp_statistical_metrics import compute_mbe, compute_rmse


def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_ANN_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
			
	return ts_latlon

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_ANN_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
			
	return ts_latlon
	              

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
var_obs = 'tmp'
var_cmip6 = 'tas'
dt = '1980-2014'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)

mbe_namz_cmip6 = []
mbe_samz_cmip6 = []
mbe_neb_cmip6 = []
mbe_sam_cmip6 = []
mbe_lpb_cmip6 = []

rmse_namz_cmip6 = []
rmse_samz_cmip6 = []
rmse_neb_cmip6 = []
rmse_sam_cmip6 = []
rmse_lpb_cmip6 = []

legend = []

for i in range(1, 19):

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	mbe_namz_cmip6.append(compute_mbe(clim_namz_cmip, clim_namz_obs))
	rmse_namz_cmip6.append(compute_rmse(clim_namz_cmip, clim_namz_obs))

	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	mbe_samz_cmip6.append(compute_mbe(clim_samz_cmip, clim_samz_obs))
	rmse_samz_cmip6.append(compute_rmse(clim_samz_cmip, clim_samz_obs))

	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	mbe_neb_cmip6.append(compute_mbe(clim_neb_cmip, clim_neb_obs))
	rmse_neb_cmip6.append(compute_rmse(clim_neb_cmip, clim_neb_obs))

	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	mbe_sam_cmip6.append(compute_mbe(clim_sam_cmip, clim_sam_obs))
	rmse_sam_cmip6.append(compute_rmse(clim_sam_cmip, clim_sam_obs))

	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	mbe_lpb_cmip6.append(compute_mbe(clim_lpb_cmip, clim_lpb_obs))
	rmse_lpb_cmip6.append(compute_rmse(clim_lpb_cmip, clim_lpb_obs))

	legend.append(cmip6[i][0])

mbe_cmip6 = np.array([mbe_lpb_cmip6,mbe_sam_cmip6,mbe_neb_cmip6,mbe_samz_cmip6,mbe_namz_cmip6])
rmse_cmip6 = np.array([rmse_lpb_cmip6,rmse_sam_cmip6,rmse_neb_cmip6,rmse_samz_cmip6,rmse_namz_cmip6])

sort_list_namz = sort_list(rmse_namz_cmip6)
model_list_namz = []
value_list_namz = []
for i in sort_list_namz:
	model_list_namz.append(cmip6[i+1][0])
	value_list_namz.append(rmse_namz_cmip6[i])

sort_list_samz = sort_list(rmse_samz_cmip6)
model_list_samz = []
value_list_samz = []
for ii in sort_list_samz:
	model_list_samz.append(cmip6[ii+1][0])
	value_list_samz.append(rmse_samz_cmip6[ii])

sort_list_neb = sort_list(rmse_neb_cmip6)
model_list_neb = []
value_list_neb = []
for iii in sort_list_neb:
	model_list_neb.append(cmip6[iii+1][0])
	value_list_neb.append(rmse_neb_cmip6[iii])

sort_list_sam = sort_list(rmse_sam_cmip6)
model_list_sam = []
value_list_sam = []
for iv in sort_list_sam:
	model_list_sam.append(cmip6[iv+1][0])
	value_list_sam.append(rmse_sam_cmip6[iv])

sort_list_lpb = sort_list(rmse_lpb_cmip6)
model_list_lpb = []
value_list_lpb = []
for v in sort_list_lpb:
	model_list_lpb.append(cmip6[v+1][0])
	value_list_lpb.append(rmse_lpb_cmip6[v])

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
name_out = 'pyplt_rank_rmse_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 6))
norm1 = colors.BoundaryNorm(boundaries=np.arange(-6, 7, 1), ncolors=256)
norm2 = colors.BoundaryNorm(boundaries=np.arange(0, 5.5, .5), ncolors=256)

xlabels = legend
ylabels = [u'LPB', u'SAM', u'NEB', u'SAMZ', u'NAMZ']

if var_cmip6 == 'pr':
	color1 = cm.BrBG
	color2 = cm.Blues
else:
	color1 = cm.bwr
	color2 = cm.Reds
	
ax = fig.add_subplot(2, 1, 1)  
pcm = ax.pcolormesh(mbe_cmip6, edgecolors='white', linewidths=2., norm=norm1, cmap=color1)
ax.set_title(u'(a) MBE', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(mbe_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(mbe_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
for y in range(mbe_cmip6.shape[0]):
    for x in range(mbe_cmip6.shape[1]):
        ax.text(x + 0.5, y + 0.5, '%.2f' % mbe_cmip6[y, x],
                 ha="center", va="center", color='k', size=8)

ax = fig.add_subplot(2, 1, 2)  
pcm = ax.pcolormesh(rmse_cmip6, edgecolors='white', linewidths=2., norm=norm2, cmap=color2)
ax.set_title(u'(b) RMSE', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(rmse_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(rmse_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='max', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
for y in range(rmse_cmip6.shape[0]):
    for x in range(rmse_cmip6.shape[1]):
        ax.text(x + 0.5, y + 0.5, '%.2f' % rmse_cmip6[y, x],
                 ha="center", va="center", color='k', size=8)
                 
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_mbe_rmse_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






