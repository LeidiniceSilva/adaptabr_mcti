# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot annual cycle of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6

def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_MON_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	obs_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
	obs_clim = []
	for mon in range(0, 11 + 1):
		obs = np.nanmean(obs_data[mon::12], axis=0)
		obs_clim.append(obs)
		
	return obs_clim

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_MON_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	mdl_data = np.nanmean(np.nanmean(value, axis=1), axis=1)
	mdl_clim = []
	for mon in range(0, 11 + 1):
		mdl = np.nanmean(mdl_data[mon::12], axis=0)
		mdl_clim.append(mdl)
	
	return mdl_clim
	              
               
# Import cmip models and obs database 
var_obs = 'tmp'
var_cmip6 = 'tas'
dt = '1980-2014'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)

be_namz_cmip6 = []
be_samz_cmip6 = []
be_neb_cmip6 = []
be_sam_cmip6 = []
be_lpb_cmip6 = []
legend = []

for i in range(1, 19):

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	be_namz_cmip = np.array(clim_namz_cmip) - np.array(clim_namz_obs)
	be_namz_cmip6.append(be_namz_cmip)
		
	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	be_samz_cmip = np.array(clim_samz_cmip) - np.array(clim_samz_obs)
	be_samz_cmip6.append(be_samz_cmip)

	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	be_neb_cmip = np.array(clim_neb_cmip) - np.array(clim_neb_obs)
	be_neb_cmip6.append(be_neb_cmip)

	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	be_sam_cmip = np.array(clim_sam_cmip) - np.array(clim_sam_obs)
	be_sam_cmip6.append(be_sam_cmip)
	
	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	be_lpb_cmip = np.array(clim_lpb_cmip) - np.array(clim_lpb_obs)
	be_lpb_cmip6.append(be_lpb_cmip)

	legend.append(cmip6[i][0])

namz_cmip6 = np.array([be_namz_cmip6[0], be_namz_cmip6[1], be_namz_cmip6[2], be_namz_cmip6[3], be_namz_cmip6[4], be_namz_cmip6[5], 
be_namz_cmip6[6], be_namz_cmip6[7], be_namz_cmip6[8], be_namz_cmip6[9], be_namz_cmip6[10], be_namz_cmip6[11], be_namz_cmip6[12], 
be_namz_cmip6[13], be_namz_cmip6[14], be_namz_cmip6[15], be_namz_cmip6[16], be_namz_cmip6[17]])
 
samz_cmip6 = np.array([be_samz_cmip6[0], be_samz_cmip6[1], be_samz_cmip6[2], be_samz_cmip6[3], be_samz_cmip6[4], be_samz_cmip6[5], 
be_samz_cmip6[6], be_samz_cmip6[7], be_samz_cmip6[8], be_samz_cmip6[9], be_samz_cmip6[10], be_samz_cmip6[11], be_samz_cmip6[12], 
be_samz_cmip6[13], be_samz_cmip6[14], be_samz_cmip6[15], be_samz_cmip6[16], be_samz_cmip6[17]])

neb_cmip6 = np.array([be_neb_cmip6[0], be_neb_cmip6[1], be_neb_cmip6[2], be_neb_cmip6[3], be_neb_cmip6[4], be_neb_cmip6[5], 
be_neb_cmip6[6], be_neb_cmip6[7], be_neb_cmip6[8], be_neb_cmip6[9], be_neb_cmip6[10], be_neb_cmip6[11], be_neb_cmip6[12], 
be_neb_cmip6[13], be_neb_cmip6[14], be_neb_cmip6[15], be_neb_cmip6[16], be_neb_cmip6[17]])

sam_cmip6 = np.array([be_sam_cmip6[0], be_sam_cmip6[1], be_sam_cmip6[2], be_sam_cmip6[3], be_sam_cmip6[4], be_sam_cmip6[5], 
be_sam_cmip6[6], be_sam_cmip6[7], be_sam_cmip6[8], be_sam_cmip6[9], be_sam_cmip6[10], be_sam_cmip6[11], be_sam_cmip6[12], 
be_sam_cmip6[13], be_sam_cmip6[14], be_sam_cmip6[15], be_sam_cmip6[16], be_sam_cmip6[17]])

lpb_cmip6 = np.array([be_lpb_cmip6[0], be_lpb_cmip6[1], be_lpb_cmip6[2], be_lpb_cmip6[3], be_lpb_cmip6[4], be_lpb_cmip6[5], 
be_lpb_cmip6[6], be_lpb_cmip6[7], be_lpb_cmip6[8], be_lpb_cmip6[9], be_lpb_cmip6[10], be_lpb_cmip6[11], be_lpb_cmip6[12], 
be_lpb_cmip6[13], be_lpb_cmip6[14], be_lpb_cmip6[15], be_lpb_cmip6[16], be_lpb_cmip6[17]])

namz_cmip6 = np.transpose(namz_cmip6)
samz_cmip6 = np.transpose(samz_cmip6)
neb_cmip6 = np.transpose(neb_cmip6)
sam_cmip6 = np.transpose(sam_cmip6)
lpb_cmip6 = np.transpose(lpb_cmip6)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(6, 10))

xlabels = legend
ylabels = ['Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez']

if var_cmip6 == 'pr':
	norm = colors.BoundaryNorm(boundaries=np.arange(-5, 6, 1), ncolors=256)
	color = cm.BrBG
else:
	norm = colors.BoundaryNorm(boundaries=np.arange(-5, 6, 1), ncolors=256)
	color = cm.bwr	
	
ax = fig.add_subplot(5, 1, 1)  
pcm = ax.pcolormesh(namz_cmip6, edgecolors='white', linewidths=2., norm=norm, cmap=color)
ax.set_title(u'(a) NAMZ', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(lpb_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(lpb_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
if var_cmip6 == 'pr':
	clb.set_label(u'Precipitação (mm d⁻¹)', size=8, rotation=90)
else:
	clb.set_label(u'Temperatura (°C)', size=8, rotation=90)
	
ax = fig.add_subplot(5, 1, 2)  
pcm = ax.pcolormesh(samz_cmip6, edgecolors='white', linewidths=2., norm=norm, cmap=color)
ax.set_title(u'(b) SAMZ', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(lpb_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(lpb_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
if var_cmip6 == 'pr':
	clb.set_label(u'Precipitação (mm d⁻¹)', size=8, rotation=90)
else:
	clb.set_label(u'Temperatura (°C)', size=8, rotation=90)
	
ax = fig.add_subplot(5, 1, 3)  
pcm = ax.pcolormesh(neb_cmip6, edgecolors='white', linewidths=2., norm=norm, cmap=color)
ax.set_title(u'(c) NEB', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(lpb_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(lpb_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
if var_cmip6 == 'pr':
	clb.set_label(u'Precipitação (mm d⁻¹)', size=8, rotation=90)
else:
	clb.set_label(u'Temperatura (°C)', size=8, rotation=90)
	
ax = fig.add_subplot(5, 1, 4)  
pcm = ax.pcolormesh(sam_cmip6, edgecolors='white', linewidths=2., norm=norm, cmap=color)
ax.set_title(u'(d) SAM', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(lpb_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(lpb_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
if var_cmip6 == 'pr':
	clb.set_label(u'Precipitação (mm d⁻¹)', size=8, rotation=90)
else:
	clb.set_label(u'Temperatura (°C)', size=8, rotation=90)
	
ax = fig.add_subplot(5, 1, 5)  
pcm = ax.pcolormesh(lpb_cmip6, edgecolors='white', linewidths=2., norm=norm, cmap=color)
ax.set_title(u'(e) LPB', loc='left', fontweight='bold', fontsize=8)	
ax.set_xticks(np.arange(lpb_cmip6.shape[1]) + 0.5)
ax.set_yticks(np.arange(lpb_cmip6.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, ax=ax, extend='both', pad=0.01)
clb.ax.yaxis.set_label_position('right')
clb.ax.tick_params(labelsize=8)
if var_cmip6 == 'pr':
	clb.set_label(u'Precipitação (mm d⁻¹)', size=8, rotation=90)
else:
	clb.set_label(u'Temperatura (°C)', size=8, rotation=90)
	
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_portrait_annual_cycle_bias_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






