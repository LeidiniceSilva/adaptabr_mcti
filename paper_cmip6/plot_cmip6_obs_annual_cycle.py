# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "June 10, 2024"
__description__ = "This script plot portrait of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from matplotlib.patches import Patch
from dict_cmip6_models_name import cmip6
from comp_stats_metrics import compute_nrmse, compute_tss, compute_corr, compute_ivs

dt = '197901-201412'
path  = '/afs/ictp.it/home/m/mda_silv/Documents/AdaptaBr_MCTI'
		    
			    
def import_obs_srf(param, area):

	arq  = '{0}/database/obs/{1}_{2}_ERA5_mon_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]

	mean = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)	
	
	annual_cycle = []
	for mon in range(0, 11 + 1):
		mean_ = np.nanmean(mean[mon::12], axis=0)
		annual_cycle.append(mean_)
			
	return annual_cycle


def import_cmip_srf(param, area, model, exp):
	
	arq   = '{0}/database/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	mean = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)	
	
	annual_cycle = []
	for mon in range(0, 11 + 1):
		mean_ = np.nanmean(mean[mon::12], axis=0)
		annual_cycle.append(mean_)
			
	return annual_cycle



def comp_mme(subdomain, best_models):

	annual_cycles = []

	for model in best_models:
		member = next(v[1] for v in cmip6.values() if v[0] == model)
		ac = import_cmip_srf('pr', subdomain, model, member)
		annual_cycles.append(ac)

	annual_cycles = np.array(annual_cycles)
	ensemble_mean = np.nanmean(annual_cycles, axis=0)

	return annual_cycles, ensemble_mean
    
    
# Import obs database and cmip models
namz_obs_pr_ac = import_obs_srf('tp', 'namz')
samz_obs_pr_ac = import_obs_srf('tp', 'samz')
sam_obs_pr_ac = import_obs_srf('tp', 'sam')
neb_obs_pr_ac = import_obs_srf('tp', 'neb')
lpb_obs_pr_ac = import_obs_srf('tp', 'lpb')

best_models_by_domain = {'namz': ['MRI-ESM2-0', 'GFDL-ESM4', 'MPI-ESM1-2-HR'],
    'samz': ['MRI-ESM2-0', 'MPI-ESM1-2-HR', 'NorESM2-MM'],
    'neb': ['MRI-ESM2-0', 'ACCESS-CM2', 'MPI-ESM1-2-HR'],
    'sam': ['GFDL-ESM4', 'MRI-ESM2-0', 'ACCESS-CM2'],
    'lpb': ['ACCESS-CM2', 'MPI-ESM1-2-HR', 'GFDL-ESM4']}

worse_models_by_domain = {'namz': ['KIOST-ESM', 'CMCC-ESM2', 'NESM3'],
    'samz': ['CMCC-ESM2', 'KIOST-ESM', 'INM-CM4-8'],
    'neb': ['KIOST-ESM', 'MIROC6', 'CNRM-ESM2-1'],
    'sam': ['NESM3', 'KIOST-ESM', 'INM-CM4-8'],
    'lpb': ['KIOST-ESM', 'MIROC6', 'NESM3']}

results_best = {}
for subdomain, best_models in best_models_by_domain.items():
    ac, ens = comp_mme(subdomain, best_models)
    results_best[subdomain] = {
        'annual_cycles': ac,
        'ensemble_mean': ens
    }

results_worse = {}
for subdomain, worse_models in worse_models_by_domain.items():
    ac, ens = comp_mme(subdomain, worse_models)
    results_worse[subdomain] = {
        'annual_cycles': ac,
        'ensemble_mean': ens
    }
    
namz_mme_best = results_best['namz']['ensemble_mean']
samz_mme_best = results_best['samz']['ensemble_mean']
neb_mme_best = results_best['neb']['ensemble_mean']
sam_mme_best = results_best['sam']['ensemble_mean']
lpb_mme_best = results_best['lpb']['ensemble_mean']

namz_mme_worse = results_worse['namz']['ensemble_mean']
samz_mme_worse = results_worse['samz']['ensemble_mean']
neb_mme_worse = results_worse['neb']['ensemble_mean']
sam_mme_worse = results_worse['sam']['ensemble_mean']
lpb_mme_worse = results_worse['lpb']['ensemble_mean']

clim_namz_cmip6, clim_samz_cmip6, clim_neb_cmip6, clim_sam_cmip6, clim_lpb_cmip6 = [], [], [], [], []
for i in range(1, 18):
	clim_namz_cmip6.append(import_cmip_srf('pr', 'namz', cmip6[i][0], cmip6[i][1]))	
	clim_samz_cmip6.append(import_cmip_srf('pr', 'samz', cmip6[i][0], cmip6[i][1]))
	clim_neb_cmip6.append(import_cmip_srf('pr', 'neb', cmip6[i][0], cmip6[i][1]))
	clim_sam_cmip6.append(import_cmip_srf('pr', 'sam', cmip6[i][0], cmip6[i][1]))
	clim_lpb_cmip6.append(import_cmip_srf('pr', 'lpb', cmip6[i][0], cmip6[i][1]))

fig = plt.figure(figsize=(14, 6))
time = np.arange(0.5, 12 + 0.5)

ax = fig.add_subplot(2, 3, 1)  
annual_cycle = ax.plot(time, namz_obs_pr_ac, time, namz_mme_best, time, namz_mme_worse, time, clim_namz_cmip6[0], 
time, clim_namz_cmip6[1], time, clim_namz_cmip6[2], time, clim_namz_cmip6[3], time, clim_namz_cmip6[4], 
time, clim_namz_cmip6[5], time, clim_namz_cmip6[6], time, clim_namz_cmip6[7], time, clim_namz_cmip6[8], 
time, clim_namz_cmip6[9], time, clim_namz_cmip6[10], time, clim_namz_cmip6[11], time, clim_namz_cmip6[12], 
time, clim_namz_cmip6[13], time, clim_namz_cmip6[14], time, clim_namz_cmip6[15], time, clim_namz_cmip6[16])
plt.title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = annual_cycle
plt.setp(l1, color='black')
plt.setp(l2, color='blue')
plt.setp(l3, color='red')
plt.setp(l4, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l5, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l6, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l7, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l8, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l9, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l10, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l11, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l12, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l13, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l14, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l15, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l16, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l17, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l18, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l19, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l20, color='gray', linewidth=0.5, alpha=0.5)

legend = ['ERA5', 'MME-best', 'MME-worst']
plt.legend(annual_cycle, legend, ncol=3, loc=1, fontsize=8)

ax = fig.add_subplot(2, 3, 2)  
annual_cycle = ax.plot(time, samz_obs_pr_ac, time, samz_mme_best, time, samz_mme_worse, time, clim_samz_cmip6[0], 
time, clim_samz_cmip6[1], time, clim_samz_cmip6[2], time, clim_samz_cmip6[3], time, clim_samz_cmip6[4], 
time, clim_samz_cmip6[5], time, clim_samz_cmip6[6], time, clim_samz_cmip6[7], time, clim_samz_cmip6[8], 
time, clim_samz_cmip6[9], time, clim_samz_cmip6[10], time, clim_samz_cmip6[11], time, clim_samz_cmip6[12], 
time, clim_samz_cmip6[13], time, clim_samz_cmip6[14], time, clim_samz_cmip6[15], time, clim_samz_cmip6[16])
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = annual_cycle
plt.setp(l1, color='black')
plt.setp(l2, color='blue')
plt.setp(l3, color='red')
plt.setp(l4, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l5, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l6, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l7, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l8, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l9, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l10, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l11, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l12, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l13, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l14, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l15, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l16, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l17, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l18, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l19, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l20, color='gray', linewidth=0.5, alpha=0.5)

ax = fig.add_subplot(2, 3, 3)  
annual_cycle = ax.plot(time, neb_obs_pr_ac, time, neb_mme_best, time, neb_mme_worse, time, clim_neb_cmip6[0], 
time, clim_neb_cmip6[1], time, clim_neb_cmip6[2], time, clim_neb_cmip6[3], time, clim_neb_cmip6[4], 
time, clim_neb_cmip6[5], time, clim_neb_cmip6[6], time, clim_neb_cmip6[7], time, clim_neb_cmip6[8], 
time, clim_neb_cmip6[9], time, clim_neb_cmip6[10], time, clim_neb_cmip6[11], time, clim_neb_cmip6[12], 
time, clim_neb_cmip6[13], time, clim_neb_cmip6[14], time, clim_neb_cmip6[15], time, clim_neb_cmip6[16])
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.xlabel('Months', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = annual_cycle
plt.setp(l1, color='black')
plt.setp(l2, color='blue')
plt.setp(l3, color='red')
plt.setp(l4, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l5, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l6, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l7, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l8, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l9, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l10, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l11, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l12, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l13, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l14, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l15, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l16, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l17, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l18, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l19, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l20, color='gray', linewidth=0.5, alpha=0.5)

ax = fig.add_subplot(2, 3, 4)  
annual_cycle = ax.plot(time, sam_obs_pr_ac, time, sam_mme_worse, time, sam_mme_best, time, clim_sam_cmip6[0], 
time, clim_sam_cmip6[1], time, clim_sam_cmip6[2], time, clim_sam_cmip6[3], time, clim_sam_cmip6[4], 
time, clim_sam_cmip6[5], time, clim_sam_cmip6[6], time, clim_sam_cmip6[7], time, clim_sam_cmip6[8], 
time, clim_sam_cmip6[9], time, clim_sam_cmip6[10], time, clim_sam_cmip6[11], time, clim_sam_cmip6[12], 
time, clim_sam_cmip6[13], time, clim_sam_cmip6[14], time, clim_sam_cmip6[15], time, clim_sam_cmip6[16])
plt.title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.xlabel('Months', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = annual_cycle
plt.setp(l1, color='black')
plt.setp(l2, color='blue')
plt.setp(l3, color='red')
plt.setp(l4, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l5, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l6, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l7, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l8, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l9, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l10, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l11, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l12, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l13, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l14, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l15, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l16, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l17, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l18, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l19, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l20, color='gray', linewidth=0.5, alpha=0.5)

ax = fig.add_subplot(2, 3, 5)  
annual_cycle = ax.plot(time, lpb_obs_pr_ac, time, lpb_mme_best, time, lpb_mme_worse, time, clim_lpb_cmip6[0], 
time, clim_lpb_cmip6[1], time, clim_lpb_cmip6[2], time, clim_lpb_cmip6[3], time, clim_lpb_cmip6[4], 
time, clim_lpb_cmip6[5], time, clim_lpb_cmip6[6], time, clim_lpb_cmip6[7], time, clim_lpb_cmip6[8], 
time, clim_lpb_cmip6[9], time, clim_lpb_cmip6[10], time, clim_lpb_cmip6[11], time, clim_lpb_cmip6[12], 
time, clim_lpb_cmip6[13], time, clim_lpb_cmip6[14], time, clim_lpb_cmip6[15], time, clim_lpb_cmip6[16])
plt.title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
plt.yticks(np.arange(0, 16, 2), fontsize=8)
plt.ylim(0, 14)
plt.xlabel('Months', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15, l16, l17, l18, l19, l20 = annual_cycle
plt.setp(l1, color='black')
plt.setp(l2, color='blue')
plt.setp(l3, color='red')
plt.setp(l4, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l5, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l6, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l7, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l8, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l9, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l10, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l11, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l12, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l13, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l14, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l15, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l16, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l17, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l18, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l19, color='gray', linewidth=0.5, alpha=0.5)
plt.setp(l20, color='gray', linewidth=0.5, alpha=0.5)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_annual_cycle_cmip6_obs_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

