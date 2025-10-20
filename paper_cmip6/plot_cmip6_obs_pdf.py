# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "June 10, 2024"
__description__ = "This script plot annual cycle of cmip6 models"

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

dt = '197901-201412'
path  = '/home/mda_silv/users/AdaptaBr_MCTI'
		    
			    
def import_obs_srf(param, area):

	arq  = '{0}/database/paper_cmip6/obs/{1}_{2}_ERA5_day_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]
	value = var[:][:,:,:]

	return value


def import_cmip_srf(param, area, model):

	exp = next(v[1] for v in cmip6.values() if v[0] == model)
	
	arq   = '{0}/database/paper_cmip6/cmip6/{3}/{1}_{2}_day_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]

	return value
 

def comp_mme(value_i, value_ii, value_iii):

	arrays = [value_i, value_ii, value_iii]
	min_time = min(arr.shape[0] for arr in arrays)
	arrays_equal = [arr[:min_time, :, :] for arr in arrays]
	mme = np.nanmean(arrays_equal, axis=0)

	return mme


def comp_mme_(value_):

	arrays_ = [value_[0], value_[1], value_[2], value_[3], value_[4], value_[5], value_[6], value_[7], value_[8], value_[9], value_[10], value_[11], value_[12], value_[13], value_[14], value_[15], value_[16]]
	min_time_ = min(arr.shape[0] for arr in arrays_)
	arrays_equal_ = [arr[:min_time_, :, :] for arr in arrays_]
	mme_ = np.nanmean(arrays_equal_, axis=0)

	return mme_
		
  
def comp_pdf(value):

	step=1
	rain_min=0.1
	rain_max=500.0
	bins = np.arange(rain_min, rain_max + step, step)
	rain_hist = np.zeros(len(bins) - 1)
	for t in range(value.shape[0]):
		values_ = value[t, :, :]
		wet = values_[values_ >= rain_min]
		hist, _ = np.histogram(wet, bins=bins)
		rain_hist += hist

	total = np.sum(rain_hist)
	rain_hist[rain_hist < 1] = np.nan  
	pdf = rain_hist / (total * step)

	bin_centers = (bins[:-1] + bins[1:]) / 2

	return pdf
	
	
# Import obs database and cmip models
namz_obs_pr = import_obs_srf('tp', 'namz')
samz_obs_pr = import_obs_srf('tp', 'samz')
sam_obs_pr = import_obs_srf('tp', 'sam')
neb_obs_pr = import_obs_srf('tp', 'neb')
lpb_obs_pr = import_obs_srf('tp', 'lpb')

clim_namz_cmip6, clim_samz_cmip6, clim_neb_cmip6, clim_sam_cmip6, clim_lpb_cmip6 = [], [], [], [], []
for i in range(1, 18):
	clim_namz_cmip6.append(import_cmip_srf('pr', 'namz', cmip6[i][0]))	
	clim_samz_cmip6.append(import_cmip_srf('pr', 'samz', cmip6[i][0]))
	clim_neb_cmip6.append(import_cmip_srf('pr', 'neb', cmip6[i][0]))
	clim_sam_cmip6.append(import_cmip_srf('pr', 'sam', cmip6[i][0]))
	clim_lpb_cmip6.append(import_cmip_srf('pr', 'lpb', cmip6[i][0]))

clim_namz_cmip6_ = comp_mme_(clim_namz_cmip6)
clim_samz_cmip6_ = comp_mme_(clim_samz_cmip6)
clim_sam_cmip6_ = comp_mme_(clim_sam_cmip6)
clim_neb_cmip6_ = comp_mme_(clim_neb_cmip6)
clim_lpb_cmip6_ = comp_mme_(clim_lpb_cmip6)

# Best
namz_cmip_pr_best1 = import_cmip_srf('pr', 'namz', 'MRI-ESM2-0')
namz_cmip_pr_best2 = import_cmip_srf('pr', 'namz', 'GFDL-ESM4')
namz_cmip_pr_best3 = import_cmip_srf('pr', 'namz', 'MPI-ESM1-2-HR')
namz_cmip_pr_best_mme = comp_mme(namz_cmip_pr_best1, namz_cmip_pr_best2, namz_cmip_pr_best3)

samz_cmip_pr_best1 = import_cmip_srf('pr', 'samz', 'MRI-ESM2-0')
samz_cmip_pr_best2 = import_cmip_srf('pr', 'samz', 'MPI-ESM1-2-HR')
samz_cmip_pr_best3 = import_cmip_srf('pr', 'samz', 'NorESM2-MM')
samz_cmip_pr_best_mme = comp_mme(samz_cmip_pr_best1, samz_cmip_pr_best2, samz_cmip_pr_best3)

neb_cmip_pr_best1 = import_cmip_srf('pr', 'neb', 'MRI-ESM2-0')
neb_cmip_pr_best2 = import_cmip_srf('pr', 'neb', 'ACCESS-CM2')
neb_cmip_pr_best3 = import_cmip_srf('pr', 'neb', 'MPI-ESM1-2-HR')
neb_cmip_pr_best_mme = comp_mme(neb_cmip_pr_best1, neb_cmip_pr_best2, neb_cmip_pr_best3)

sam_cmip_pr_best1 = import_cmip_srf('pr', 'sam', 'GFDL-ESM4')
sam_cmip_pr_best2 = import_cmip_srf('pr', 'sam', 'MRI-ESM2-0')
sam_cmip_pr_best3 = import_cmip_srf('pr', 'sam', 'ACCESS-CM2')
sam_cmip_pr_best_mme = comp_mme(sam_cmip_pr_best1, sam_cmip_pr_best2, sam_cmip_pr_best3)

lpb_cmip_pr_best1 = import_cmip_srf('pr', 'lpb', 'ACCESS-CM2')
lpb_cmip_pr_best2 = import_cmip_srf('pr', 'lpb', 'MPI-ESM1-2-HR')
lpb_cmip_pr_best3 = import_cmip_srf('pr', 'lpb', 'GFDL-ESM4')
lpb_cmip_pr_best_mme = comp_mme(lpb_cmip_pr_best1, lpb_cmip_pr_best2, lpb_cmip_pr_best3)

# Worst
namz_cmip_pr_worst1 = import_cmip_srf('pr', 'namz', 'KIOST-ESM')
namz_cmip_pr_worst2 = import_cmip_srf('pr', 'namz', 'CMCC-ESM2')
namz_cmip_pr_worst3 = import_cmip_srf('pr', 'namz', 'NESM3')
namz_cmip_pr_worst_mme = comp_mme(namz_cmip_pr_worst1, namz_cmip_pr_worst2, namz_cmip_pr_worst3)

samz_cmip_pr_worst1 = import_cmip_srf('pr', 'samz', 'CMCC-ESM2')
samz_cmip_pr_worst2 = import_cmip_srf('pr', 'samz', 'KIOST-ESM')
samz_cmip_pr_worst3 = import_cmip_srf('pr', 'samz', 'INM-CM4-8')
samz_cmip_pr_worst_mme = comp_mme(samz_cmip_pr_worst1, samz_cmip_pr_worst2, samz_cmip_pr_worst3)

neb_cmip_pr_worst1 = import_cmip_srf('pr', 'neb', 'KIOST-ESM')
neb_cmip_pr_worst2 = import_cmip_srf('pr', 'neb', 'MIROC6')
neb_cmip_pr_worst3 = import_cmip_srf('pr', 'neb', 'CNRM-ESM2-1')
neb_cmip_pr_worst_mme = comp_mme(neb_cmip_pr_worst1, neb_cmip_pr_worst2, neb_cmip_pr_worst3)

sam_cmip_pr_worst1 = import_cmip_srf('pr', 'sam', 'NESM3')
sam_cmip_pr_worst2 = import_cmip_srf('pr', 'sam', 'KIOST-ESM')
sam_cmip_pr_worst3 = import_cmip_srf('pr', 'sam', 'INM-CM4-8')
sam_cmip_pr_worst_mme = comp_mme(sam_cmip_pr_worst1, sam_cmip_pr_worst2, sam_cmip_pr_worst3)

lpb_cmip_pr_worst1 = import_cmip_srf('pr', 'lpb', 'KIOST-ESM')
lpb_cmip_pr_worst2 = import_cmip_srf('pr', 'lpb', 'MIROC6')
lpb_cmip_pr_worst3 = import_cmip_srf('pr', 'lpb', 'NESM3')
lpb_cmip_pr_worst_mme = comp_mme(lpb_cmip_pr_worst1, lpb_cmip_pr_worst2, lpb_cmip_pr_worst3)

# Compute pdf
namz_obs_pr_pdf = comp_pdf(namz_obs_pr)
samz_obs_pr_pdf = comp_pdf(samz_obs_pr)
sam_obs_pr_pdf = comp_pdf(sam_obs_pr)
neb_obs_pr_pdf = comp_pdf(neb_obs_pr)
lpb_obs_pr_pdf = comp_pdf(lpb_obs_pr)

namz_cmip_pr_pdf = comp_pdf(clim_namz_cmip6_)
samz_cmip_pr_pdf = comp_pdf(clim_samz_cmip6_)
sam_cmip_pr_pdf = comp_pdf(clim_sam_cmip6_)
neb_cmip_pr_pdf = comp_pdf(clim_neb_cmip6_)
lpb_cmip_pr_pdf = comp_pdf(clim_lpb_cmip6_)

namz_cmip_pr_best_pdf = comp_pdf(namz_cmip_pr_best_mme)
samz_cmip_pr_best_pdf = comp_pdf(samz_cmip_pr_best_mme)
neb_cmip_pr_best_pdf = comp_pdf(neb_cmip_pr_best_mme)
sam_cmip_pr_best_pdf = comp_pdf(sam_cmip_pr_best_mme)
lpb_cmip_pr_best_pdf = comp_pdf(lpb_cmip_pr_best_mme)

namz_cmip_pr_worst_pdf = comp_pdf(namz_cmip_pr_worst_mme)
samz_cmip_pr_worst_pdf = comp_pdf(samz_cmip_pr_worst_mme)
neb_cmip_pr_worst_pdf = comp_pdf(neb_cmip_pr_worst_mme)
sam_cmip_pr_worst_pdf = comp_pdf(sam_cmip_pr_worst_mme)
lpb_cmip_pr_worst_pdf = comp_pdf(lpb_cmip_pr_worst_mme)

# Plot figure
fig = plt.figure(figsize=(14, 6))
time = np.arange(0.5, 12 + 0.5)

ax = fig.add_subplot(2, 3, 1)  
pdf_plot = ax.plot(namz_obs_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='ERA5')
pdf_plot = ax.plot(namz_cmip_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.5, linestyle='None', label='MME')
pdf_plot = ax.plot(namz_cmip_pr_best_pdf, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.75, linestyle='None', label='MME-best')
pdf_plot = ax.plot(namz_cmip_pr_worst_pdf, marker='o', markersize=4, mfc='red', mec='red', alpha=0.75, linestyle='None', label='MME-worst')
plt.title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
ax.set_yscale('log')
plt.legend(ncol=2, loc=1, fontsize=8)

ax = fig.add_subplot(2, 3, 2)  
pdf_plot = ax.plot(samz_obs_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='ERA5')
pdf_plot = ax.plot(samz_cmip_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.5, linestyle='None', label='MME')
pdf_plot = ax.plot(samz_cmip_pr_best_pdf, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.75, linestyle='None', label='MME-best')
pdf_plot = ax.plot(samz_cmip_pr_worst_pdf, marker='o', markersize=4, mfc='red', mec='red', alpha=0.75, linestyle='None', label='MME-worst')
plt.title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
ax.set_yscale('log')

ax = fig.add_subplot(2, 3, 3)  
pdf_plot = ax.plot(neb_obs_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='ERA5')
pdf_plot = ax.plot(neb_cmip_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.5, linestyle='None', label='MME')
pdf_plot = ax.plot(neb_cmip_pr_best_pdf, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.75, linestyle='None', label='MME-best')
pdf_plot = ax.plot(neb_cmip_pr_worst_pdf, marker='o', markersize=4, mfc='red', mec='red', alpha=0.75, linestyle='None', label='MME-worst')
plt.title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
ax.set_yscale('log')

ax = fig.add_subplot(2, 3, 4)  
pdf_plot = ax.plot(sam_obs_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='ERA5')
pdf_plot = ax.plot(sam_cmip_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.5, linestyle='None', label='MME')
pdf_plot = ax.plot(sam_cmip_pr_best_pdf, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.75, linestyle='None', label='MME-best')
pdf_plot = ax.plot(sam_cmip_pr_worst_pdf, marker='o', markersize=4, mfc='red', mec='red', alpha=0.75, linestyle='None', label='MME-worst')
plt.title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
plt.ylabel('Frequency (#)', fontsize=8, fontweight='bold')
plt.xlabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
ax.set_yscale('log')

ax = fig.add_subplot(2, 3, 5)  
pdf_plot = ax.plot(lpb_obs_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.75, linestyle='None', label='ERA5')
pdf_plot = ax.plot(lpb_cmip_pr_pdf, marker='o', markersize=4, mfc='black', mec='black', alpha=0.5, linestyle='None', label='MME')
pdf_plot = ax.plot(lpb_cmip_pr_best_pdf, marker='o', markersize=4, mfc='blue', mec='blue', alpha=0.75, linestyle='None', label='MME-best')
pdf_plot = ax.plot(lpb_cmip_pr_worst_pdf, marker='o', markersize=4, mfc='red', mec='red', alpha=0.75, linestyle='None', label='MME-worst')
plt.title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
plt.xlabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
plt.grid(linestyle='--')
ax.set_yscale('log')

# Path out to save figure
path_out = '{0}/figs/paper_cmip6'.format(path)
name_out = 'pyplt_pdf_cmip6_obs_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

