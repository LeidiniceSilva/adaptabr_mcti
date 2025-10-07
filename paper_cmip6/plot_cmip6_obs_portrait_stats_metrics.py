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

from dict_cmip6_models_name import cmip6
from comp_stats_metrics import compute_nrmse, compute_tss, compute_corr, compute_ivs

dt = '197901-201412'
path  = '/afs/ictp.it/home/m/mda_silv/Documents/AdaptaBr_MCTI'


def annual_cycle(dataset):

	ac_mean = []
	for mon in range(0, 11 + 1):
		ac = np.nanmean(dataset[mon::12], axis=0)
		ac_mean.append(ac)
	
	return ac_mean
	

def year_mean(dataset):
	
	years = dataset.shape[0] // 12	
	dataset_reshaped = dataset.reshape((years, 12))
	yr_mean = np.mean(dataset_reshaped, axis=1)
	
	return yr_mean
	

def latlon(dataset):

	latlon_ts = []
	for lat_idx in range(dataset.shape[0]):
		for lon_idx in range(dataset.shape[1]):
			value = dataset[lat_idx, lon_idx]
			latlon_ts.append(value)
				
	return latlon_ts
			    
			    
def import_obs_srf(param, area):

	arq  = '{0}/database/obs/{1}_{2}_ERA5_mon_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]

	mean = np.nanmean(var[:][:,:,:], axis=0)
	mean_ts = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)	
		
	xy_ts = latlon(mean)
	ac_ts = annual_cycle(mean_ts)
	yr_ts = year_mean(mean_ts)

	return xy_ts, ac_ts, yr_ts


def import_obs_atm(param, area):
	
	
	arq  = '{0}/database/obs/{1}_{2}_ERA5_mon_{3}_lonlat.nc'.format(path, param, area, dt)	
	data = netCDF4.Dataset(arq)
	var  = data.variables[param][:] 
	lat  = data.variables['lat'][:]
	lon  = data.variables['lon'][:]

	mean_850 = np.nanmean(var[:][:,0,:,:], axis=0)
	mean_500 = np.nanmean(var[:][:,1,:,:], axis=0)
	mean_200 = np.nanmean(var[:][:,2,:,:], axis=0)
	
	mean_850_ts = np.nanmean(np.nanmean(var[:][:,0,:,:], axis=1), axis=1)
	mean_500_ts = np.nanmean(np.nanmean(var[:][:,1,:,:], axis=1), axis=1)
	mean_200_ts = np.nanmean(np.nanmean(var[:][:,2,:,:], axis=1), axis=1)

	xy_ts_850 = latlon(mean_850)
	xy_ts_500 = latlon(mean_500)
	xy_ts_200 = latlon(mean_200)

	ac_ts_850 = annual_cycle(mean_850_ts)
	ac_ts_500 = annual_cycle(mean_500_ts)
	ac_ts_200 = annual_cycle(mean_200_ts)

	yr_ts_850 = year_mean(mean_850_ts)
	yr_ts_500 = year_mean(mean_500_ts)
	yr_ts_200 = year_mean(mean_200_ts)
		
	return xy_ts_850, xy_ts_500, xy_ts_200, ac_ts_850, ac_ts_500, ac_ts_200, yr_ts_850, yr_ts_500, yr_ts_200
	
		
def import_cmip_srf(param, area, model, exp):
	
	arq   = '{0}/database/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	mean = np.nanmean(var[:][:,:,:], axis=0)
	mean_ts = np.nanmean(np.nanmean(var[:][:,:,:], axis=1), axis=1)	
		
	xy_ts = latlon(mean)
	ac_ts = annual_cycle(mean_ts)
	yr_ts = year_mean(mean_ts)

	return xy_ts, ac_ts, yr_ts
	              
 
def import_cmip_atm(param, area, model, exp):
	
	
	arq   = '{0}/database/cmip6/{3}/{1}_{2}_Amon_{3}_historical_{4}_{5}_lonlat.nc'.format(path, param, area, model, exp, dt)			
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	
	mean_850 = np.nanmean(var[:][:,0,:,:], axis=0)
	mean_500 = np.nanmean(var[:][:,1,:,:], axis=0)
	mean_200 = np.nanmean(var[:][:,2,:,:], axis=0)
	
	mean_850_ts = np.nanmean(np.nanmean(var[:][:,0,:,:], axis=1), axis=1)
	mean_500_ts = np.nanmean(np.nanmean(var[:][:,1,:,:], axis=1), axis=1)
	mean_200_ts = np.nanmean(np.nanmean(var[:][:,2,:,:], axis=1), axis=1)

	xy_ts_850 = latlon(mean_850)
	xy_ts_500 = latlon(mean_500)
	xy_ts_200 = latlon(mean_200)

	ac_ts_850 = annual_cycle(mean_850_ts)
	ac_ts_500 = annual_cycle(mean_500_ts)
	ac_ts_200 = annual_cycle(mean_200_ts)

	yr_ts_850 = year_mean(mean_850_ts)
	yr_ts_500 = year_mean(mean_500_ts)
	yr_ts_200 = year_mean(mean_200_ts)
		
	return xy_ts_850, xy_ts_500, xy_ts_200, ac_ts_850, ac_ts_500, ac_ts_200, yr_ts_850, yr_ts_500, yr_ts_200
	
	      
# Import obs database and cmip models
namz_obs_pr_xy, namz_obs_pr_ac, namz_obs_pr_yr = import_obs_srf('tp', 'namz')
namz_obs_ps_xy, namz_obs_ps_ac, namz_obs_ps_yr = import_obs_srf('sp', 'namz')
namz_obs_t850_xy, namz_obs_t500_xy, namz_obs_t200_xy, namz_obs_t850_ac, namz_obs_t500_ac, namz_obs_t200_ac, namz_obs_t850_yr, namz_obs_t500_yr, namz_obs_t200_yr = import_obs_atm('t', 'namz')
namz_obs_u850_xy, namz_obs_u500_xy, namz_obs_u200_xy, namz_obs_u850_ac, namz_obs_u500_ac, namz_obs_u200_ac, namz_obs_u850_yr, namz_obs_u500_yr, namz_obs_u200_yr = import_obs_atm('u', 'namz')
namz_obs_v850_xy, namz_obs_v500_xy, namz_obs_v200_xy, namz_obs_v850_ac, namz_obs_v500_ac, namz_obs_v200_ac, namz_obs_v850_yr, namz_obs_v500_yr, namz_obs_v200_yr = import_obs_atm('v', 'namz')
namz_obs_q850_xy, namz_obs_q500_xy, namz_obs_q200_xy, namz_obs_q850_ac, namz_obs_q500_ac, namz_obs_q200_ac, namz_obs_q850_yr, namz_obs_q500_yr, namz_obs_q200_yr = import_obs_atm('q', 'namz')

samz_obs_pr_xy, samz_obs_pr_ac, samz_obs_pr_yr = import_obs_srf('tp', 'samz')
samz_obs_ps_xy, samz_obs_ps_ac, samz_obs_ps_yr = import_obs_srf('sp', 'samz')
samz_obs_t850_xy, samz_obs_t500_xy, samz_obs_t200_xy, samz_obs_t850_ac, samz_obs_t500_ac, samz_obs_t200_ac, samz_obs_t850_yr, samz_obs_t500_yr, samz_obs_t200_yr = import_obs_atm('t', 'samz')
samz_obs_u850_xy, samz_obs_u500_xy, samz_obs_u200_xy, samz_obs_u850_ac, samz_obs_u500_ac, samz_obs_u200_ac, samz_obs_u850_yr, samz_obs_u500_yr, samz_obs_u200_yr = import_obs_atm('u', 'samz')
samz_obs_v850_xy, samz_obs_v500_xy, samz_obs_v200_xy, samz_obs_v850_ac, samz_obs_v500_ac, samz_obs_v200_ac, samz_obs_v850_yr, samz_obs_v500_yr, samz_obs_v200_yr = import_obs_atm('v', 'samz')
samz_obs_q850_xy, samz_obs_q500_xy, samz_obs_q200_xy, samz_obs_q850_ac, samz_obs_q500_ac, samz_obs_q200_ac, samz_obs_q850_yr, samz_obs_q500_yr, samz_obs_q200_yr = import_obs_atm('q', 'samz')

sam_obs_pr_xy, sam_obs_pr_ac, sam_obs_pr_yr = import_obs_srf('tp', 'sam')
sam_obs_ps_xy, sam_obs_ps_ac, sam_obs_ps_yr = import_obs_srf('sp', 'sam')
sam_obs_t850_xy, sam_obs_t500_xy, sam_obs_t200_xy, sam_obs_t850_ac, sam_obs_t500_ac, sam_obs_t200_ac, sam_obs_t850_yr, sam_obs_t500_yr, sam_obs_t200_yr = import_obs_atm('t', 'sam')
sam_obs_u850_xy, sam_obs_u500_xy, sam_obs_u200_xy, sam_obs_u850_ac, sam_obs_u500_ac, sam_obs_u200_ac, sam_obs_u850_yr, sam_obs_u500_yr, sam_obs_u200_yr = import_obs_atm('u', 'sam')
sam_obs_v850_xy, sam_obs_v500_xy, sam_obs_v200_xy, sam_obs_v850_ac, sam_obs_v500_ac, sam_obs_v200_ac, sam_obs_v850_yr, sam_obs_v500_yr, sam_obs_v200_yr = import_obs_atm('v', 'sam')
sam_obs_q850_xy, sam_obs_q500_xy, sam_obs_q200_xy, sam_obs_q850_ac, sam_obs_q500_ac, sam_obs_q200_ac, sam_obs_q850_yr, sam_obs_q500_yr, sam_obs_q200_yr = import_obs_atm('q', 'sam')

neb_obs_pr_xy, neb_obs_pr_ac, neb_obs_pr_yr = import_obs_srf('tp', 'neb')
neb_obs_ps_xy, neb_obs_ps_ac, neb_obs_ps_yr = import_obs_srf('sp', 'neb')
neb_obs_t850_xy, neb_obs_t500_xy, neb_obs_t200_xy, neb_obs_t850_ac, neb_obs_t500_ac, neb_obs_t200_ac, neb_obs_t850_yr, neb_obs_t500_yr, neb_obs_t200_yr = import_obs_atm('t', 'neb')
neb_obs_u850_xy, neb_obs_u500_xy, neb_obs_u200_xy, neb_obs_u850_ac, neb_obs_u500_ac, neb_obs_u200_ac, neb_obs_u850_yr, neb_obs_u500_yr, neb_obs_u200_yr = import_obs_atm('u', 'neb')
neb_obs_v850_xy, neb_obs_v500_xy, neb_obs_v200_xy, neb_obs_v850_ac, neb_obs_v500_ac, neb_obs_v200_ac, neb_obs_v850_yr, neb_obs_v500_yr, neb_obs_v200_yr = import_obs_atm('v', 'neb')
neb_obs_q850_xy, neb_obs_q500_xy, neb_obs_q200_xy, neb_obs_q850_ac, neb_obs_q500_ac, neb_obs_q200_ac, neb_obs_q850_yr, neb_obs_q500_yr, neb_obs_q200_yr = import_obs_atm('q', 'neb')

lpb_obs_pr_xy, lpb_obs_pr_ac, lpb_obs_pr_yr = import_obs_srf('tp', 'lpb')
lpb_obs_ps_xy, lpb_obs_ps_ac, lpb_obs_ps_yr = import_obs_srf('sp', 'lpb')
lpb_obs_t850_xy, lpb_obs_t500_xy, lpb_obs_t200_xy, lpb_obs_t850_ac, lpb_obs_t500_ac, lpb_obs_t200_ac, lpb_obs_t850_yr, lpb_obs_t500_yr, lpb_obs_t200_yr = import_obs_atm('t', 'lpb')
lpb_obs_u850_xy, lpb_obs_u500_xy, lpb_obs_u200_xy, lpb_obs_u850_ac, lpb_obs_u500_ac, lpb_obs_u200_ac, lpb_obs_u850_yr, lpb_obs_u500_yr, lpb_obs_u200_yr = import_obs_atm('u', 'lpb')
lpb_obs_v850_xy, lpb_obs_v500_xy, lpb_obs_v200_xy, lpb_obs_v850_ac, lpb_obs_v500_ac, lpb_obs_v200_ac, lpb_obs_v850_yr, lpb_obs_v500_yr, lpb_obs_v200_yr = import_obs_atm('v', 'lpb')
lpb_obs_q850_xy, lpb_obs_q500_xy, lpb_obs_q200_xy, lpb_obs_q850_ac, lpb_obs_q500_ac, lpb_obs_q200_ac, lpb_obs_q850_yr, lpb_obs_q500_yr, lpb_obs_q200_yr = import_obs_atm('q', 'lpb')

# Define variables
m1_namz_pr, m1_namz_ps, m1_namz_t850, m1_namz_t500, m1_namz_t200, m1_namz_u850, m1_namz_u500, m1_namz_u200, m1_namz_v850, m1_namz_v500, m1_namz_v200, m1_namz_q850, m1_namz_q500, m1_namz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m2_namz_pr, m2_namz_ps, m2_namz_t850, m2_namz_t500, m2_namz_t200, m2_namz_u850, m2_namz_u500, m2_namz_u200, m2_namz_v850, m2_namz_v500, m2_namz_v200, m2_namz_q850, m2_namz_q500, m2_namz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m3_namz_pr, m3_namz_ps, m3_namz_t850, m3_namz_t500, m3_namz_t200, m3_namz_u850, m3_namz_u500, m3_namz_u200, m3_namz_v850, m3_namz_v500, m3_namz_v200, m3_namz_q850, m3_namz_q500, m3_namz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m4_namz_pr, m4_namz_ps, m4_namz_t850, m4_namz_t500, m4_namz_t200, m4_namz_u850, m4_namz_u500, m4_namz_u200, m4_namz_v850, m4_namz_v500, m4_namz_v200, m4_namz_q850, m4_namz_q500, m4_namz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []

m1_samz_pr, m1_samz_ps, m1_samz_t850, m1_samz_t500, m1_samz_t200, m1_samz_u850, m1_samz_u500, m1_samz_u200, m1_samz_v850, m1_samz_v500, m1_samz_v200, m1_samz_q850, m1_samz_q500, m1_samz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m2_samz_pr, m2_samz_ps, m2_samz_t850, m2_samz_t500, m2_samz_t200, m2_samz_u850, m2_samz_u500, m2_samz_u200, m2_samz_v850, m2_samz_v500, m2_samz_v200, m2_samz_q850, m2_samz_q500, m2_samz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m3_samz_pr, m3_samz_ps, m3_samz_t850, m3_samz_t500, m3_samz_t200, m3_samz_u850, m3_samz_u500, m3_samz_u200, m3_samz_v850, m3_samz_v500, m3_samz_v200, m3_samz_q850, m3_samz_q500, m3_samz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m4_samz_pr, m4_samz_ps, m4_samz_t850, m4_samz_t500, m4_samz_t200, m4_samz_u850, m4_samz_u500, m4_samz_u200, m4_samz_v850, m4_samz_v500, m4_samz_v200, m4_samz_q850, m4_samz_q500, m4_samz_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []

m1_sam_pr, m1_sam_ps, m1_sam_t850, m1_sam_t500, m1_sam_t200, m1_sam_u850, m1_sam_u500, m1_sam_u200, m1_sam_v850, m1_sam_v500, m1_sam_v200, m1_sam_q850, m1_sam_q500, m1_sam_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m2_sam_pr, m2_sam_ps, m2_sam_t850, m2_sam_t500, m2_sam_t200, m2_sam_u850, m2_sam_u500, m2_sam_u200, m2_sam_v850, m2_sam_v500, m2_sam_v200, m2_sam_q850, m2_sam_q500, m2_sam_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m3_sam_pr, m3_sam_ps, m3_sam_t850, m3_sam_t500, m3_sam_t200, m3_sam_u850, m3_sam_u500, m3_sam_u200, m3_sam_v850, m3_sam_v500, m3_sam_v200, m3_sam_q850, m3_sam_q500, m3_sam_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m4_sam_pr, m4_sam_ps, m4_sam_t850, m4_sam_t500, m4_sam_t200, m4_sam_u850, m4_sam_u500, m4_sam_u200, m4_sam_v850, m4_sam_v500, m4_sam_v200, m4_sam_q850, m4_sam_q500, m4_sam_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []

m1_neb_pr, m1_neb_ps, m1_neb_t850, m1_neb_t500, m1_neb_t200, m1_neb_u850, m1_neb_u500, m1_neb_u200, m1_neb_v850, m1_neb_v500, m1_neb_v200, m1_neb_q850, m1_neb_q500, m1_neb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m2_neb_pr, m2_neb_ps, m2_neb_t850, m2_neb_t500, m2_neb_t200, m2_neb_u850, m2_neb_u500, m2_neb_u200, m2_neb_v850, m2_neb_v500, m2_neb_v200, m2_neb_q850, m2_neb_q500, m2_neb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m3_neb_pr, m3_neb_ps, m3_neb_t850, m3_neb_t500, m3_neb_t200, m3_neb_u850, m3_neb_u500, m3_neb_u200, m3_neb_v850, m3_neb_v500, m3_neb_v200, m3_neb_q850, m3_neb_q500, m3_neb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m4_neb_pr, m4_neb_ps, m4_neb_t850, m4_neb_t500, m4_neb_t200, m4_neb_u850, m4_neb_u500, m4_neb_u200, m4_neb_v850, m4_neb_v500, m4_neb_v200, m4_neb_q850, m4_neb_q500, m4_neb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []

m1_lpb_pr, m1_lpb_ps, m1_lpb_t850, m1_lpb_t500, m1_lpb_t200, m1_lpb_u850, m1_lpb_u500, m1_lpb_u200, m1_lpb_v850, m1_lpb_v500, m1_lpb_v200, m1_lpb_q850, m1_lpb_q500, m1_lpb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m2_lpb_pr, m2_lpb_ps, m2_lpb_t850, m2_lpb_t500, m2_lpb_t200, m2_lpb_u850, m2_lpb_u500, m2_lpb_u200, m2_lpb_v850, m2_lpb_v500, m2_lpb_v200, m2_lpb_q850, m2_lpb_q500, m2_lpb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m3_lpb_pr, m3_lpb_ps, m3_lpb_t850, m3_lpb_t500, m3_lpb_t200, m3_lpb_u850, m3_lpb_u500, m3_lpb_u200, m3_lpb_v850, m3_lpb_v500, m3_lpb_v200, m3_lpb_q850, m3_lpb_q500, m3_lpb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []
m4_lpb_pr, m4_lpb_ps, m4_lpb_t850, m4_lpb_t500, m4_lpb_t200, m4_lpb_u850, m4_lpb_u500, m4_lpb_u200, m4_lpb_v850, m4_lpb_v500, m4_lpb_v200, m4_lpb_q850, m4_lpb_q500, m4_lpb_q200 = [], [], [], [], [], [], [], [], [], [], [], [], [], []

legend = []
for i in range(1, 18):
	print(cmip6[i][0])

	namz_cmip_pr_xy, namz_cmip_pr_ac, namz_cmip_pr_yr = import_cmip_srf('pr', 'namz', cmip6[i][0], cmip6[i][1])
	namz_cmip_ps_xy, namz_cmip_ps_ac, namz_cmip_ps_yr = import_cmip_srf('ps', 'namz', cmip6[i][0], cmip6[i][1])
	namz_cmip_t850_xy, namz_cmip_t500_xy, namz_cmip_t200_xy, namz_cmip_t850_ac, namz_cmip_t500_ac, namz_cmip_t200_ac, namz_cmip_t850_yr, namz_cmip_t500_yr, namz_cmip_t200_yr = import_cmip_atm('ta', 'namz', cmip6[i][0], cmip6[i][1])
	namz_cmip_u850_xy, namz_cmip_u500_xy, namz_cmip_u200_xy, namz_cmip_u850_ac, namz_cmip_u500_ac, namz_cmip_u200_ac, namz_cmip_u850_yr, namz_cmip_u500_yr, namz_cmip_u200_yr = import_cmip_atm('ua', 'namz', cmip6[i][0], cmip6[i][1])
	namz_cmip_v850_xy, namz_cmip_v500_xy, namz_cmip_v200_xy, namz_cmip_v850_ac, namz_cmip_v500_ac, namz_cmip_v200_ac, namz_cmip_v850_yr, namz_cmip_v500_yr, namz_cmip_v200_yr = import_cmip_atm('va', 'namz', cmip6[i][0], cmip6[i][1])
	namz_cmip_q850_xy, namz_cmip_q500_xy, namz_cmip_q200_xy, namz_cmip_q850_ac, namz_cmip_q500_ac, namz_cmip_q200_ac, namz_cmip_q850_yr, namz_cmip_q500_yr, namz_cmip_q200_yr = import_cmip_atm('hus', 'namz', cmip6[i][0], cmip6[i][1])

	samz_cmip_pr_xy, samz_cmip_pr_ac, samz_cmip_pr_yr = import_cmip_srf('pr', 'samz', cmip6[i][0], cmip6[i][1])
	samz_cmip_ps_xy, samz_cmip_ps_ac, samz_cmip_ps_yr = import_cmip_srf('ps', 'samz', cmip6[i][0], cmip6[i][1])
	samz_cmip_t850_xy, samz_cmip_t500_xy, samz_cmip_t200_xy, samz_cmip_t850_ac, samz_cmip_t500_ac, samz_cmip_t200_ac, samz_cmip_t850_yr, samz_cmip_t500_yr, samz_cmip_t200_yr = import_cmip_atm('ta', 'samz', cmip6[i][0], cmip6[i][1])
	samz_cmip_u850_xy, samz_cmip_u500_xy, samz_cmip_u200_xy, samz_cmip_u850_ac, samz_cmip_u500_ac, samz_cmip_u200_ac, samz_cmip_u850_yr, samz_cmip_u500_yr, samz_cmip_u200_yr = import_cmip_atm('ua', 'samz', cmip6[i][0], cmip6[i][1])
	samz_cmip_v850_xy, samz_cmip_v500_xy, samz_cmip_v200_xy, samz_cmip_v850_ac, samz_cmip_v500_ac, samz_cmip_v200_ac, samz_cmip_v850_yr, samz_cmip_v500_yr, samz_cmip_v200_yr = import_cmip_atm('va', 'samz', cmip6[i][0], cmip6[i][1])
	samz_cmip_q850_xy, samz_cmip_q500_xy, samz_cmip_q200_xy, samz_cmip_q850_ac, samz_cmip_q500_ac, samz_cmip_q200_ac, samz_cmip_q850_yr, samz_cmip_q500_yr, samz_cmip_q200_yr = import_cmip_atm('hus', 'samz', cmip6[i][0], cmip6[i][1])

	sam_cmip_pr_xy, sam_cmip_pr_ac, sam_cmip_pr_yr = import_cmip_srf('pr', 'sam', cmip6[i][0], cmip6[i][1])
	sam_cmip_ps_xy, sam_cmip_ps_ac, sam_cmip_ps_yr = import_cmip_srf('ps', 'sam', cmip6[i][0], cmip6[i][1])
	sam_cmip_t850_xy, sam_cmip_t500_xy, sam_cmip_t200_xy, sam_cmip_t850_ac, sam_cmip_t500_ac, sam_cmip_t200_ac, sam_cmip_t850_yr, sam_cmip_t500_yr, sam_cmip_t200_yr = import_cmip_atm('ta', 'sam', cmip6[i][0], cmip6[i][1])
	sam_cmip_u850_xy, sam_cmip_u500_xy, sam_cmip_u200_xy, sam_cmip_u850_ac, sam_cmip_u500_ac, sam_cmip_u200_ac, sam_cmip_u850_yr, sam_cmip_u500_yr, sam_cmip_u200_yr = import_cmip_atm('ua', 'sam', cmip6[i][0], cmip6[i][1])
	sam_cmip_v850_xy, sam_cmip_v500_xy, sam_cmip_v200_xy, sam_cmip_v850_ac, sam_cmip_v500_ac, sam_cmip_v200_ac, sam_cmip_v850_yr, sam_cmip_v500_yr, sam_cmip_v200_yr = import_cmip_atm('va', 'sam', cmip6[i][0], cmip6[i][1])
	sam_cmip_q850_xy, sam_cmip_q500_xy, sam_cmip_q200_xy, sam_cmip_q850_ac, sam_cmip_q500_ac, sam_cmip_q200_ac, sam_cmip_q850_yr, sam_cmip_q500_yr, sam_cmip_q200_yr = import_cmip_atm('hus', 'sam', cmip6[i][0], cmip6[i][1])

	neb_cmip_pr_xy, neb_cmip_pr_ac, neb_cmip_pr_yr = import_cmip_srf('pr', 'neb', cmip6[i][0], cmip6[i][1])
	neb_cmip_ps_xy, neb_cmip_ps_ac, neb_cmip_ps_yr = import_cmip_srf('ps', 'neb', cmip6[i][0], cmip6[i][1])
	neb_cmip_t850_xy, neb_cmip_t500_xy, neb_cmip_t200_xy, neb_cmip_t850_ac, neb_cmip_t500_ac, neb_cmip_t200_ac, neb_cmip_t850_yr, neb_cmip_t500_yr, neb_cmip_t200_yr = import_cmip_atm('ta', 'neb', cmip6[i][0], cmip6[i][1])
	neb_cmip_u850_xy, neb_cmip_u500_xy, neb_cmip_u200_xy, neb_cmip_u850_ac, neb_cmip_u500_ac, neb_cmip_u200_ac, neb_cmip_u850_yr, neb_cmip_u500_yr, neb_cmip_u200_yr = import_cmip_atm('ua', 'neb', cmip6[i][0], cmip6[i][1])
	neb_cmip_v850_xy, neb_cmip_v500_xy, neb_cmip_v200_xy, neb_cmip_v850_ac, neb_cmip_v500_ac, neb_cmip_v200_ac, neb_cmip_v850_yr, neb_cmip_v500_yr, neb_cmip_v200_yr = import_cmip_atm('va', 'neb', cmip6[i][0], cmip6[i][1])
	neb_cmip_q850_xy, neb_cmip_q500_xy, neb_cmip_q200_xy, neb_cmip_q850_ac, neb_cmip_q500_ac, neb_cmip_q200_ac, neb_cmip_q850_yr, neb_cmip_q500_yr, neb_cmip_q200_yr = import_cmip_atm('hus', 'neb', cmip6[i][0], cmip6[i][1])

	lpb_cmip_pr_xy, lpb_cmip_pr_ac, lpb_cmip_pr_yr = import_cmip_srf('pr', 'lpb', cmip6[i][0], cmip6[i][1])
	lpb_cmip_ps_xy, lpb_cmip_ps_ac, lpb_cmip_ps_yr = import_cmip_srf('ps', 'lpb', cmip6[i][0], cmip6[i][1])
	lpb_cmip_t850_xy, lpb_cmip_t500_xy, lpb_cmip_t200_xy, lpb_cmip_t850_ac, lpb_cmip_t500_ac, lpb_cmip_t200_ac, lpb_cmip_t850_yr, lpb_cmip_t500_yr, lpb_cmip_t200_yr = import_cmip_atm('ta', 'lpb', cmip6[i][0], cmip6[i][1])
	lpb_cmip_u850_xy, lpb_cmip_u500_xy, lpb_cmip_u200_xy, lpb_cmip_u850_ac, lpb_cmip_u500_ac, lpb_cmip_u200_ac, lpb_cmip_u850_yr, lpb_cmip_u500_yr, lpb_cmip_u200_yr = import_cmip_atm('ua', 'lpb', cmip6[i][0], cmip6[i][1])
	lpb_cmip_v850_xy, lpb_cmip_v500_xy, lpb_cmip_v200_xy, lpb_cmip_v850_ac, lpb_cmip_v500_ac, lpb_cmip_v200_ac, lpb_cmip_v850_yr, lpb_cmip_v500_yr, lpb_cmip_v200_yr = import_cmip_atm('va', 'lpb', cmip6[i][0], cmip6[i][1])
	lpb_cmip_q850_xy, lpb_cmip_q500_xy, lpb_cmip_q200_xy, lpb_cmip_q850_ac, lpb_cmip_q500_ac, lpb_cmip_q200_ac, lpb_cmip_q850_yr, lpb_cmip_q500_yr, lpb_cmip_q200_yr = import_cmip_atm('hus', 'lpb', cmip6[i][0], cmip6[i][1])

	# NAMZ
	m1_namz_pr.append(compute_nrmse(namz_cmip_pr_xy, namz_obs_pr_xy))
	m1_namz_ps.append(compute_nrmse(namz_cmip_ps_xy, namz_obs_ps_xy))
	m1_namz_t850.append(compute_nrmse(namz_cmip_t850_xy, namz_obs_t850_xy))
	m1_namz_t500.append(compute_nrmse(namz_cmip_t500_xy, namz_obs_t500_xy))
	m1_namz_t200.append(compute_nrmse(namz_cmip_t200_xy, namz_obs_t200_xy))
	m1_namz_u850.append(compute_nrmse(namz_cmip_u850_xy, namz_obs_u850_xy))
	m1_namz_u500.append(compute_nrmse(namz_cmip_u500_xy, namz_obs_u500_xy))
	m1_namz_u200.append(compute_nrmse(namz_cmip_u200_xy, namz_obs_u200_xy))
	m1_namz_v850.append(compute_nrmse(namz_cmip_v850_xy, namz_obs_v850_xy))
	m1_namz_v500.append(compute_nrmse(namz_cmip_v500_xy, namz_obs_v500_xy))
	m1_namz_v200.append(compute_nrmse(namz_cmip_v200_xy, namz_obs_v200_xy))
	m1_namz_q850.append(compute_nrmse(namz_cmip_q850_xy, namz_obs_q850_xy))
	m1_namz_q500.append(compute_nrmse(namz_cmip_q500_xy, namz_obs_q500_xy))
	m1_namz_q200.append(compute_nrmse(namz_cmip_q200_xy, namz_obs_q200_xy))

	m2_namz_pr.append(compute_tss(namz_cmip_pr_xy, namz_obs_pr_xy))
	m2_namz_ps.append(compute_tss(namz_cmip_ps_xy, namz_obs_ps_xy))
	m2_namz_t850.append(compute_tss(namz_cmip_t850_xy, namz_obs_t850_xy))
	m2_namz_t500.append(compute_tss(namz_cmip_t500_xy, namz_obs_t500_xy))
	m2_namz_t200.append(compute_tss(namz_cmip_t200_xy, namz_obs_t200_xy))
	m2_namz_u850.append(compute_tss(namz_cmip_u850_xy, namz_obs_u850_xy))
	m2_namz_u500.append(compute_tss(namz_cmip_u500_xy, namz_obs_u500_xy))
	m2_namz_u200.append(compute_tss(namz_cmip_u200_xy, namz_obs_u200_xy))
	m2_namz_v850.append(compute_tss(namz_cmip_v850_xy, namz_obs_v850_xy))
	m2_namz_v500.append(compute_tss(namz_cmip_v500_xy, namz_obs_v500_xy))
	m2_namz_v200.append(compute_tss(namz_cmip_v200_xy, namz_obs_v200_xy))
	m2_namz_q850.append(compute_tss(namz_cmip_q850_xy, namz_obs_q850_xy))
	m2_namz_q500.append(compute_tss(namz_cmip_q500_xy, namz_obs_q500_xy))
	m2_namz_q200.append(compute_tss(namz_cmip_q200_xy, namz_obs_q200_xy))

	m3_namz_pr.append(compute_corr(namz_cmip_pr_ac, namz_obs_pr_ac))
	m3_namz_ps.append(compute_corr(namz_cmip_ps_ac, namz_obs_ps_ac))
	m3_namz_t850.append(compute_corr(namz_cmip_t850_ac, namz_obs_t850_ac))
	m3_namz_t500.append(compute_corr(namz_cmip_t500_ac, namz_obs_t500_ac))
	m3_namz_t200.append(compute_corr(namz_cmip_t200_ac, namz_obs_t200_ac))
	m3_namz_u850.append(compute_corr(namz_cmip_u850_ac, namz_obs_u850_ac))
	m3_namz_u500.append(compute_corr(namz_cmip_u500_ac, namz_obs_u500_ac))
	m3_namz_u200.append(compute_corr(namz_cmip_u200_ac, namz_obs_u200_ac))
	m3_namz_v850.append(compute_corr(namz_cmip_v850_ac, namz_obs_v850_ac))
	m3_namz_v500.append(compute_corr(namz_cmip_v500_ac, namz_obs_v500_ac))
	m3_namz_v200.append(compute_corr(namz_cmip_v200_ac, namz_obs_v200_ac))
	m3_namz_q850.append(compute_corr(namz_cmip_q850_ac, namz_obs_q850_ac))
	m3_namz_q500.append(compute_corr(namz_cmip_q500_ac, namz_obs_q500_ac))
	m3_namz_q200.append(compute_corr(namz_cmip_q200_ac, namz_obs_q200_ac))

	m4_namz_pr.append(compute_ivs(namz_cmip_pr_yr, namz_obs_pr_yr))
	m4_namz_ps.append(compute_ivs(namz_cmip_ps_yr, namz_obs_ps_yr))
	m4_namz_t850.append(compute_ivs(namz_cmip_t850_yr, namz_obs_t850_yr))
	m4_namz_t500.append(compute_ivs(namz_cmip_t500_yr, namz_obs_t500_yr))
	m4_namz_t200.append(compute_ivs(namz_cmip_t200_yr, namz_obs_t200_yr))
	m4_namz_u850.append(compute_ivs(namz_cmip_u850_yr, namz_obs_u850_yr))
	m4_namz_u500.append(compute_ivs(namz_cmip_u500_yr, namz_obs_u500_yr))
	m4_namz_u200.append(compute_ivs(namz_cmip_u200_yr, namz_obs_u200_yr))
	m4_namz_v850.append(compute_ivs(namz_cmip_v850_yr, namz_obs_v850_yr))
	m4_namz_v500.append(compute_ivs(namz_cmip_v500_yr, namz_obs_v500_yr))
	m4_namz_v200.append(compute_ivs(namz_cmip_v200_yr, namz_obs_v200_yr))
	m4_namz_q850.append(compute_ivs(namz_cmip_q850_yr, namz_obs_q850_yr))
	m4_namz_q500.append(compute_ivs(namz_cmip_q500_yr, namz_obs_q500_yr))
	m4_namz_q200.append(compute_ivs(namz_cmip_q200_yr, namz_obs_q200_yr))

	# SAMZ
	m1_samz_pr.append(compute_nrmse(samz_cmip_pr_xy, samz_obs_pr_xy))
	m1_samz_ps.append(compute_nrmse(samz_cmip_ps_xy, samz_obs_ps_xy))
	m1_samz_t850.append(compute_nrmse(samz_cmip_t850_xy, samz_obs_t850_xy))
	m1_samz_t500.append(compute_nrmse(samz_cmip_t500_xy, samz_obs_t500_xy))
	m1_samz_t200.append(compute_nrmse(samz_cmip_t200_xy, samz_obs_t200_xy))
	m1_samz_u850.append(compute_nrmse(samz_cmip_u850_xy, samz_obs_u850_xy))
	m1_samz_u500.append(compute_nrmse(samz_cmip_u500_xy, samz_obs_u500_xy))
	m1_samz_u200.append(compute_nrmse(samz_cmip_u200_xy, samz_obs_u200_xy))
	m1_samz_v850.append(compute_nrmse(samz_cmip_v850_xy, samz_obs_v850_xy))
	m1_samz_v500.append(compute_nrmse(samz_cmip_v500_xy, samz_obs_v500_xy))
	m1_samz_v200.append(compute_nrmse(samz_cmip_v200_xy, samz_obs_v200_xy))
	m1_samz_q850.append(compute_nrmse(samz_cmip_q850_xy, samz_obs_q850_xy))
	m1_samz_q500.append(compute_nrmse(samz_cmip_q500_xy, samz_obs_q500_xy))
	m1_samz_q200.append(compute_nrmse(samz_cmip_q200_xy, samz_obs_q200_xy))

	m2_samz_pr.append(compute_tss(samz_cmip_pr_xy, samz_obs_pr_xy))
	m2_samz_ps.append(compute_tss(samz_cmip_ps_xy, samz_obs_ps_xy))
	m2_samz_t850.append(compute_tss(samz_cmip_t850_xy, samz_obs_t850_xy))
	m2_samz_t500.append(compute_tss(samz_cmip_t500_xy, samz_obs_t500_xy))
	m2_samz_t200.append(compute_tss(samz_cmip_t200_xy, samz_obs_t200_xy))
	m2_samz_u850.append(compute_tss(samz_cmip_u850_xy, samz_obs_u850_xy))
	m2_samz_u500.append(compute_tss(samz_cmip_u500_xy, samz_obs_u500_xy))
	m2_samz_u200.append(compute_tss(samz_cmip_u200_xy, samz_obs_u200_xy))
	m2_samz_v850.append(compute_tss(samz_cmip_v850_xy, samz_obs_v850_xy))
	m2_samz_v500.append(compute_tss(samz_cmip_v500_xy, samz_obs_v500_xy))
	m2_samz_v200.append(compute_tss(samz_cmip_v200_xy, samz_obs_v200_xy))
	m2_samz_q850.append(compute_tss(samz_cmip_q850_xy, samz_obs_q850_xy))
	m2_samz_q500.append(compute_tss(samz_cmip_q500_xy, samz_obs_q500_xy))
	m2_samz_q200.append(compute_tss(samz_cmip_q200_xy, samz_obs_q200_xy))

	m3_samz_pr.append(compute_corr(samz_cmip_pr_ac, samz_obs_pr_ac))
	m3_samz_ps.append(compute_corr(samz_cmip_ps_ac, samz_obs_ps_ac))
	m3_samz_t850.append(compute_corr(samz_cmip_t850_ac, samz_obs_t850_ac))
	m3_samz_t500.append(compute_corr(samz_cmip_t500_ac, samz_obs_t500_ac))
	m3_samz_t200.append(compute_corr(samz_cmip_t200_ac, samz_obs_t200_ac))
	m3_samz_u850.append(compute_corr(samz_cmip_u850_ac, samz_obs_u850_ac))
	m3_samz_u500.append(compute_corr(samz_cmip_u500_ac, samz_obs_u500_ac))
	m3_samz_u200.append(compute_corr(samz_cmip_u200_ac, samz_obs_u200_ac))
	m3_samz_v850.append(compute_corr(samz_cmip_v850_ac, samz_obs_v850_ac))
	m3_samz_v500.append(compute_corr(samz_cmip_v500_ac, samz_obs_v500_ac))
	m3_samz_v200.append(compute_corr(samz_cmip_v200_ac, samz_obs_v200_ac))
	m3_samz_q850.append(compute_corr(samz_cmip_q850_ac, samz_obs_q850_ac))
	m3_samz_q500.append(compute_corr(samz_cmip_q500_ac, samz_obs_q500_ac))
	m3_samz_q200.append(compute_corr(samz_cmip_q200_ac, samz_obs_q200_ac))

	m4_samz_pr.append(compute_ivs(samz_cmip_pr_yr, samz_obs_pr_yr))
	m4_samz_ps.append(compute_ivs(samz_cmip_ps_yr, samz_obs_ps_yr))
	m4_samz_t850.append(compute_ivs(samz_cmip_t850_yr, samz_obs_t850_yr))
	m4_samz_t500.append(compute_ivs(samz_cmip_t500_yr, samz_obs_t500_yr))
	m4_samz_t200.append(compute_ivs(samz_cmip_t200_yr, samz_obs_t200_yr))
	m4_samz_u850.append(compute_ivs(samz_cmip_u850_yr, samz_obs_u850_yr))
	m4_samz_u500.append(compute_ivs(samz_cmip_u500_yr, samz_obs_u500_yr))
	m4_samz_u200.append(compute_ivs(samz_cmip_u200_yr, samz_obs_u200_yr))
	m4_samz_v850.append(compute_ivs(samz_cmip_v850_yr, samz_obs_v850_yr))
	m4_samz_v500.append(compute_ivs(samz_cmip_v500_yr, samz_obs_v500_yr))
	m4_samz_v200.append(compute_ivs(samz_cmip_v200_yr, samz_obs_v200_yr))
	m4_samz_q850.append(compute_ivs(samz_cmip_q850_yr, samz_obs_q850_yr))
	m4_samz_q500.append(compute_ivs(samz_cmip_q500_yr, samz_obs_q500_yr))
	m4_samz_q200.append(compute_ivs(samz_cmip_q200_yr, samz_obs_q200_yr))

	# SAM
	m1_sam_pr.append(compute_nrmse(sam_cmip_pr_xy, sam_obs_pr_xy))
	m1_sam_ps.append(compute_nrmse(sam_cmip_ps_xy, sam_obs_ps_xy))
	m1_sam_t850.append(compute_nrmse(sam_cmip_t850_xy, sam_obs_t850_xy))
	m1_sam_t500.append(compute_nrmse(sam_cmip_t500_xy, sam_obs_t500_xy))
	m1_sam_t200.append(compute_nrmse(sam_cmip_t200_xy, sam_obs_t200_xy))
	m1_sam_u850.append(compute_nrmse(sam_cmip_u850_xy, sam_obs_u850_xy))
	m1_sam_u500.append(compute_nrmse(sam_cmip_u500_xy, sam_obs_u500_xy))
	m1_sam_u200.append(compute_nrmse(sam_cmip_u200_xy, sam_obs_u200_xy))
	m1_sam_v850.append(compute_nrmse(sam_cmip_v850_xy, sam_obs_v850_xy))
	m1_sam_v500.append(compute_nrmse(sam_cmip_v500_xy, sam_obs_v500_xy))
	m1_sam_v200.append(compute_nrmse(sam_cmip_v200_xy, sam_obs_v200_xy))
	m1_sam_q850.append(compute_nrmse(sam_cmip_q850_xy, sam_obs_q850_xy))
	m1_sam_q500.append(compute_nrmse(sam_cmip_q500_xy, sam_obs_q500_xy))
	m1_sam_q200.append(compute_nrmse(sam_cmip_q200_xy, sam_obs_q200_xy))

	m2_sam_pr.append(compute_tss(sam_cmip_pr_xy, sam_obs_pr_xy))
	m2_sam_ps.append(compute_tss(sam_cmip_ps_xy, sam_obs_ps_xy))
	m2_sam_t850.append(compute_tss(sam_cmip_t850_xy, sam_obs_t850_xy))
	m2_sam_t500.append(compute_tss(sam_cmip_t500_xy, sam_obs_t500_xy))
	m2_sam_t200.append(compute_tss(sam_cmip_t200_xy, sam_obs_t200_xy))
	m2_sam_u850.append(compute_tss(sam_cmip_u850_xy, sam_obs_u850_xy))
	m2_sam_u500.append(compute_tss(sam_cmip_u500_xy, sam_obs_u500_xy))
	m2_sam_u200.append(compute_tss(sam_cmip_u200_xy, sam_obs_u200_xy))
	m2_sam_v850.append(compute_tss(sam_cmip_v850_xy, sam_obs_v850_xy))
	m2_sam_v500.append(compute_tss(sam_cmip_v500_xy, sam_obs_v500_xy))
	m2_sam_v200.append(compute_tss(sam_cmip_v200_xy, sam_obs_v200_xy))
	m2_sam_q850.append(compute_tss(sam_cmip_q850_xy, sam_obs_q850_xy))
	m2_sam_q500.append(compute_tss(sam_cmip_q500_xy, sam_obs_q500_xy))
	m2_sam_q200.append(compute_tss(sam_cmip_q200_xy, sam_obs_q200_xy))
	
	m3_sam_pr.append(compute_corr(sam_cmip_pr_ac, sam_obs_pr_ac))
	m3_sam_ps.append(compute_corr(sam_cmip_ps_ac, sam_obs_ps_ac))
	m3_sam_t850.append(compute_corr(sam_cmip_t850_ac, sam_obs_t850_ac))
	m3_sam_t500.append(compute_corr(sam_cmip_t500_ac, sam_obs_t500_ac))
	m3_sam_t200.append(compute_corr(sam_cmip_t200_ac, sam_obs_t200_ac))
	m3_sam_u850.append(compute_corr(sam_cmip_u850_ac, sam_obs_u850_ac))
	m3_sam_u500.append(compute_corr(sam_cmip_u500_ac, sam_obs_u500_ac))
	m3_sam_u200.append(compute_corr(sam_cmip_u200_ac, sam_obs_u200_ac))
	m3_sam_v850.append(compute_corr(sam_cmip_v850_ac, sam_obs_v850_ac))
	m3_sam_v500.append(compute_corr(sam_cmip_v500_ac, sam_obs_v500_ac))
	m3_sam_v200.append(compute_corr(sam_cmip_v200_ac, sam_obs_v200_ac))
	m3_sam_q850.append(compute_corr(sam_cmip_q850_ac, sam_obs_q850_ac))
	m3_sam_q500.append(compute_corr(sam_cmip_q500_ac, sam_obs_q500_ac))
	m3_sam_q200.append(compute_corr(sam_cmip_q200_ac, sam_obs_q200_ac))
	
	m4_sam_pr.append(compute_ivs(sam_cmip_pr_yr, sam_obs_pr_yr))
	m4_sam_ps.append(compute_ivs(sam_cmip_ps_yr, sam_obs_ps_yr))
	m4_sam_t850.append(compute_ivs(sam_cmip_t850_yr, sam_obs_t850_yr))
	m4_sam_t500.append(compute_ivs(sam_cmip_t500_yr, sam_obs_t500_yr))
	m4_sam_t200.append(compute_ivs(sam_cmip_t200_yr, sam_obs_t200_yr))
	m4_sam_u850.append(compute_ivs(sam_cmip_u850_yr, sam_obs_u850_yr))
	m4_sam_u500.append(compute_ivs(sam_cmip_u500_yr, sam_obs_u500_yr))
	m4_sam_u200.append(compute_ivs(sam_cmip_u200_yr, sam_obs_u200_yr))
	m4_sam_v850.append(compute_ivs(sam_cmip_v850_yr, sam_obs_v850_yr))
	m4_sam_v500.append(compute_ivs(sam_cmip_v500_yr, sam_obs_v500_yr))
	m4_sam_v200.append(compute_ivs(sam_cmip_v200_yr, sam_obs_v200_yr))
	m4_sam_q850.append(compute_ivs(sam_cmip_q850_yr, sam_obs_q850_yr))
	m4_sam_q500.append(compute_ivs(sam_cmip_q500_yr, sam_obs_q500_yr))
	m4_sam_q200.append(compute_ivs(sam_cmip_q200_yr, sam_obs_q200_yr))		

	# NEB
	m1_neb_pr.append(compute_nrmse(neb_cmip_pr_xy, neb_obs_pr_xy))
	m1_neb_ps.append(compute_nrmse(neb_cmip_ps_xy, neb_obs_ps_xy))
	m1_neb_t850.append(compute_nrmse(neb_cmip_t850_xy, neb_obs_t850_xy))
	m1_neb_t500.append(compute_nrmse(neb_cmip_t500_xy, neb_obs_t500_xy))
	m1_neb_t200.append(compute_nrmse(neb_cmip_t200_xy, neb_obs_t200_xy))
	m1_neb_u850.append(compute_nrmse(neb_cmip_u850_xy, neb_obs_u850_xy))
	m1_neb_u500.append(compute_nrmse(neb_cmip_u500_xy, neb_obs_u500_xy))
	m1_neb_u200.append(compute_nrmse(neb_cmip_u200_xy, neb_obs_u200_xy))
	m1_neb_v850.append(compute_nrmse(neb_cmip_v850_xy, neb_obs_v850_xy))
	m1_neb_v500.append(compute_nrmse(neb_cmip_v500_xy, neb_obs_v500_xy))
	m1_neb_v200.append(compute_nrmse(neb_cmip_v200_xy, neb_obs_v200_xy))
	m1_neb_q850.append(compute_nrmse(neb_cmip_q850_xy, neb_obs_q850_xy))
	m1_neb_q500.append(compute_nrmse(neb_cmip_q500_xy, neb_obs_q500_xy))
	m1_neb_q200.append(compute_nrmse(neb_cmip_q200_xy, neb_obs_q200_xy))

	m2_neb_pr.append(compute_tss(neb_cmip_pr_xy, neb_obs_pr_xy))
	m2_neb_ps.append(compute_tss(neb_cmip_ps_xy, neb_obs_ps_xy))
	m2_neb_t850.append(compute_tss(neb_cmip_t850_xy, neb_obs_t850_xy))
	m2_neb_t500.append(compute_tss(neb_cmip_t500_xy, neb_obs_t500_xy))
	m2_neb_t200.append(compute_tss(neb_cmip_t200_xy, neb_obs_t200_xy))
	m2_neb_u850.append(compute_tss(neb_cmip_u850_xy, neb_obs_u850_xy))
	m2_neb_u500.append(compute_tss(neb_cmip_u500_xy, neb_obs_u500_xy))
	m2_neb_u200.append(compute_tss(neb_cmip_u200_xy, neb_obs_u200_xy))
	m2_neb_v850.append(compute_tss(neb_cmip_v850_xy, neb_obs_v850_xy))
	m2_neb_v500.append(compute_tss(neb_cmip_v500_xy, neb_obs_v500_xy))
	m2_neb_v200.append(compute_tss(neb_cmip_v200_xy, neb_obs_v200_xy))
	m2_neb_q850.append(compute_tss(neb_cmip_q850_xy, neb_obs_q850_xy))
	m2_neb_q500.append(compute_tss(neb_cmip_q500_xy, neb_obs_q500_xy))
	m2_neb_q200.append(compute_tss(neb_cmip_q200_xy, neb_obs_q200_xy))

	m3_neb_pr.append(compute_corr(neb_cmip_pr_ac, neb_obs_pr_ac))
	m3_neb_ps.append(compute_corr(neb_cmip_ps_ac, neb_obs_ps_ac))
	m3_neb_t850.append(compute_corr(neb_cmip_t850_ac, neb_obs_t850_ac))
	m3_neb_t500.append(compute_corr(neb_cmip_t500_ac, neb_obs_t500_ac))
	m3_neb_t200.append(compute_corr(neb_cmip_t200_ac, neb_obs_t200_ac))
	m3_neb_u850.append(compute_corr(neb_cmip_u850_ac, neb_obs_u850_ac))
	m3_neb_u500.append(compute_corr(neb_cmip_u500_ac, neb_obs_u500_ac))
	m3_neb_u200.append(compute_corr(neb_cmip_u200_ac, neb_obs_u200_ac))
	m3_neb_v850.append(compute_corr(neb_cmip_v850_ac, neb_obs_v850_ac))
	m3_neb_v500.append(compute_corr(neb_cmip_v500_ac, neb_obs_v500_ac))
	m3_neb_v200.append(compute_corr(neb_cmip_v200_ac, neb_obs_v200_ac))
	m3_neb_q850.append(compute_corr(neb_cmip_q850_ac, neb_obs_q850_ac))
	m3_neb_q500.append(compute_corr(neb_cmip_q500_ac, neb_obs_q500_ac))
	m3_neb_q200.append(compute_corr(neb_cmip_q200_ac, neb_obs_q200_ac))

	m4_neb_pr.append(compute_ivs(neb_cmip_pr_yr, neb_obs_pr_yr))
	m4_neb_ps.append(compute_ivs(neb_cmip_ps_yr, neb_obs_ps_yr))
	m4_neb_t850.append(compute_ivs(neb_cmip_t850_yr, neb_obs_t850_yr))
	m4_neb_t500.append(compute_ivs(neb_cmip_t500_yr, neb_obs_t500_yr))
	m4_neb_t200.append(compute_ivs(neb_cmip_t200_yr, neb_obs_t200_yr))
	m4_neb_u850.append(compute_ivs(neb_cmip_u850_yr, neb_obs_u850_yr))
	m4_neb_u500.append(compute_ivs(neb_cmip_u500_yr, neb_obs_u500_yr))
	m4_neb_u200.append(compute_ivs(neb_cmip_u200_yr, neb_obs_u200_yr))
	m4_neb_v850.append(compute_ivs(neb_cmip_v850_yr, neb_obs_v850_yr))
	m4_neb_v500.append(compute_ivs(neb_cmip_v500_yr, neb_obs_v500_yr))
	m4_neb_v200.append(compute_ivs(neb_cmip_v200_yr, neb_obs_v200_yr))
	m4_neb_q850.append(compute_ivs(neb_cmip_q850_yr, neb_obs_q850_yr))
	m4_neb_q500.append(compute_ivs(neb_cmip_q500_yr, neb_obs_q500_yr))
	m4_neb_q200.append(compute_ivs(neb_cmip_q200_yr, neb_obs_q200_yr))	

	# LPB
	m1_lpb_pr.append(compute_nrmse(lpb_cmip_pr_xy, lpb_obs_pr_xy))
	m1_lpb_ps.append(compute_nrmse(lpb_cmip_ps_xy, lpb_obs_ps_xy))
	m1_lpb_t850.append(compute_nrmse(lpb_cmip_t850_xy, lpb_obs_t850_xy))
	m1_lpb_t500.append(compute_nrmse(lpb_cmip_t500_xy, lpb_obs_t500_xy))
	m1_lpb_t200.append(compute_nrmse(lpb_cmip_t200_xy, lpb_obs_t200_xy))
	m1_lpb_u850.append(compute_nrmse(lpb_cmip_u850_xy, lpb_obs_u850_xy))
	m1_lpb_u500.append(compute_nrmse(lpb_cmip_u500_xy, lpb_obs_u500_xy))
	m1_lpb_u200.append(compute_nrmse(lpb_cmip_u200_xy, lpb_obs_u200_xy))
	m1_lpb_v850.append(compute_nrmse(lpb_cmip_v850_xy, lpb_obs_v850_xy))
	m1_lpb_v500.append(compute_nrmse(lpb_cmip_v500_xy, lpb_obs_v500_xy))
	m1_lpb_v200.append(compute_nrmse(lpb_cmip_v200_xy, lpb_obs_v200_xy))
	m1_lpb_q850.append(compute_nrmse(lpb_cmip_q850_xy, lpb_obs_q850_xy))
	m1_lpb_q500.append(compute_nrmse(lpb_cmip_q500_xy, lpb_obs_q500_xy))
	m1_lpb_q200.append(compute_nrmse(lpb_cmip_q200_xy, lpb_obs_q200_xy))

	m2_lpb_pr.append(compute_tss(lpb_cmip_pr_xy, lpb_obs_pr_xy))
	m2_lpb_ps.append(compute_tss(lpb_cmip_ps_xy, lpb_obs_ps_xy))
	m2_lpb_t850.append(compute_tss(lpb_cmip_t850_xy, lpb_obs_t850_xy))
	m2_lpb_t500.append(compute_tss(lpb_cmip_t500_xy, lpb_obs_t500_xy))
	m2_lpb_t200.append(compute_tss(lpb_cmip_t200_xy, lpb_obs_t200_xy))
	m2_lpb_u850.append(compute_tss(lpb_cmip_u850_xy, lpb_obs_u850_xy))
	m2_lpb_u500.append(compute_tss(lpb_cmip_u500_xy, lpb_obs_u500_xy))
	m2_lpb_u200.append(compute_tss(lpb_cmip_u200_xy, lpb_obs_u200_xy))
	m2_lpb_v850.append(compute_tss(lpb_cmip_v850_xy, lpb_obs_v850_xy))
	m2_lpb_v500.append(compute_tss(lpb_cmip_v500_xy, lpb_obs_v500_xy))
	m2_lpb_v200.append(compute_tss(lpb_cmip_v200_xy, lpb_obs_v200_xy))
	m2_lpb_q850.append(compute_tss(lpb_cmip_q850_xy, lpb_obs_q850_xy))
	m2_lpb_q500.append(compute_tss(lpb_cmip_q500_xy, lpb_obs_q500_xy))
	m2_lpb_q200.append(compute_tss(lpb_cmip_q200_xy, lpb_obs_q200_xy))
	
	m3_lpb_pr.append(compute_corr(lpb_cmip_pr_ac, lpb_obs_pr_ac))
	m3_lpb_ps.append(compute_corr(lpb_cmip_ps_ac, lpb_obs_ps_ac))
	m3_lpb_t850.append(compute_corr(lpb_cmip_t850_ac, lpb_obs_t850_ac))
	m3_lpb_t500.append(compute_corr(lpb_cmip_t500_ac, lpb_obs_t500_ac))
	m3_lpb_t200.append(compute_corr(lpb_cmip_t200_ac, lpb_obs_t200_ac))
	m3_lpb_u850.append(compute_corr(lpb_cmip_u850_ac, lpb_obs_u850_ac))
	m3_lpb_u500.append(compute_corr(lpb_cmip_u500_ac, lpb_obs_u500_ac))
	m3_lpb_u200.append(compute_corr(lpb_cmip_u200_ac, lpb_obs_u200_ac))
	m3_lpb_v850.append(compute_corr(lpb_cmip_v850_ac, lpb_obs_v850_ac))
	m3_lpb_v500.append(compute_corr(lpb_cmip_v500_ac, lpb_obs_v500_ac))
	m3_lpb_v200.append(compute_corr(lpb_cmip_v200_ac, lpb_obs_v200_ac))
	m3_lpb_q850.append(compute_corr(lpb_cmip_q850_ac, lpb_obs_q850_ac))
	m3_lpb_q500.append(compute_corr(lpb_cmip_q500_ac, lpb_obs_q500_ac))
	m3_lpb_q200.append(compute_corr(lpb_cmip_q200_ac, lpb_obs_q200_ac))
	
	m4_lpb_pr.append(compute_ivs(lpb_cmip_pr_yr, lpb_obs_pr_yr))
	m4_lpb_ps.append(compute_ivs(lpb_cmip_ps_yr, lpb_obs_ps_yr))
	m4_lpb_t850.append(compute_ivs(lpb_cmip_t850_yr, lpb_obs_t850_yr))
	m4_lpb_t500.append(compute_ivs(lpb_cmip_t500_yr, lpb_obs_t500_yr))
	m4_lpb_t200.append(compute_ivs(lpb_cmip_t200_yr, lpb_obs_t200_yr))
	m4_lpb_u850.append(compute_ivs(lpb_cmip_u850_yr, lpb_obs_u850_yr))
	m4_lpb_u500.append(compute_ivs(lpb_cmip_u500_yr, lpb_obs_u500_yr))
	m4_lpb_u200.append(compute_ivs(lpb_cmip_u200_yr, lpb_obs_u200_yr))
	m4_lpb_v850.append(compute_ivs(lpb_cmip_v850_yr, lpb_obs_v850_yr))
	m4_lpb_v500.append(compute_ivs(lpb_cmip_v500_yr, lpb_obs_v500_yr))
	m4_lpb_v200.append(compute_ivs(lpb_cmip_v200_yr, lpb_obs_v200_yr))
	m4_lpb_q850.append(compute_ivs(lpb_cmip_q850_yr, lpb_obs_q850_yr))
	m4_lpb_q500.append(compute_ivs(lpb_cmip_q500_yr, lpb_obs_q500_yr))
	m4_lpb_q200.append(compute_ivs(lpb_cmip_q200_yr, lpb_obs_q200_yr))
						
	legend.append(i)

m1_namz = np.array([m1_namz_q200, m1_namz_q500, m1_namz_q850, m1_namz_v200, m1_namz_v500, m1_namz_v850, m1_namz_u200, m1_namz_u500, m1_namz_u850, m1_namz_t200, m1_namz_t500, m1_namz_t850, m1_namz_ps, m1_namz_pr])
m2_namz = np.array([m2_namz_q200, m2_namz_q500, m2_namz_q850, m2_namz_v200, m2_namz_v500, m2_namz_v850, m2_namz_u200, m2_namz_u500, m2_namz_u850, m2_namz_t200, m2_namz_t500, m2_namz_t850, m2_namz_ps, m2_namz_pr])
m3_namz = np.array([m3_namz_q200, m3_namz_q500, m3_namz_q850, m3_namz_v200, m3_namz_v500, m3_namz_v850, m3_namz_u200, m3_namz_u500, m3_namz_u850, m3_namz_t200, m3_namz_t500, m3_namz_t850, m3_namz_ps, m3_namz_pr])
m4_namz = np.array([m4_namz_q200, m4_namz_q500, m4_namz_q850, m4_namz_v200, m4_namz_v500, m4_namz_v850, m4_namz_u200, m4_namz_u500, m4_namz_u850, m4_namz_t200, m4_namz_t500, m4_namz_t850, m4_namz_ps, m4_namz_pr])

m1_samz = np.array([m1_samz_q200, m1_samz_q500, m1_samz_q850, m1_samz_v200, m1_samz_v500, m1_samz_v850, m1_samz_u200, m1_samz_u500, m1_samz_u850, m1_samz_t200, m1_samz_t500, m1_samz_t850, m1_samz_ps, m1_samz_pr])
m2_samz = np.array([m2_samz_q200, m2_samz_q500, m2_samz_q850, m2_samz_v200, m2_samz_v500, m2_samz_v850, m2_samz_u200, m2_samz_u500, m2_samz_u850, m2_samz_t200, m2_samz_t500, m2_samz_t850, m2_samz_ps, m2_samz_pr])
m3_samz = np.array([m3_samz_q200, m3_samz_q500, m3_samz_q850, m3_samz_v200, m3_samz_v500, m3_samz_v850, m3_samz_u200, m3_samz_u500, m3_samz_u850, m3_samz_t200, m3_samz_t500, m3_samz_t850, m3_samz_ps, m3_samz_pr])
m4_samz = np.array([m4_samz_q200, m4_samz_q500, m4_samz_q850, m4_samz_v200, m4_samz_v500, m4_samz_v850, m4_samz_u200, m4_samz_u500, m4_samz_u850, m4_samz_t200, m4_samz_t500, m4_samz_t850, m4_samz_ps, m4_samz_pr])

m1_sam = np.array([m1_sam_q200, m1_sam_q500, m1_sam_q850, m1_sam_v200, m1_sam_v500, m1_sam_v850, m1_sam_u200, m1_sam_u500, m1_sam_u850, m1_sam_t200, m1_sam_t500, m1_sam_t850, m1_sam_ps, m1_sam_pr])
m2_sam = np.array([m2_sam_q200, m2_sam_q500, m2_sam_q850, m2_sam_v200, m2_sam_v500, m2_sam_v850, m2_sam_u200, m2_sam_u500, m2_sam_u850, m2_sam_t200, m2_sam_t500, m2_sam_t850, m2_sam_ps, m2_sam_pr])
m3_sam = np.array([m3_sam_q200, m3_sam_q500, m3_sam_q850, m3_sam_v200, m3_sam_v500, m3_sam_v850, m3_sam_u200, m3_sam_u500, m3_sam_u850, m3_sam_t200, m3_sam_t500, m3_sam_t850, m3_sam_ps, m3_sam_pr])
m4_sam = np.array([m4_sam_q200, m4_sam_q500, m4_sam_q850, m4_sam_v200, m4_sam_v500, m4_sam_v850, m4_sam_u200, m4_sam_u500, m4_sam_u850, m4_sam_t200, m4_sam_t500, m4_sam_t850, m4_sam_ps, m4_sam_pr])

m1_neb = np.array([m1_neb_q200, m1_neb_q500, m1_neb_q850, m1_neb_v200, m1_neb_v500, m1_neb_v850, m1_neb_u200, m1_neb_u500, m1_neb_u850, m1_neb_t200, m1_neb_t500, m1_neb_t850, m1_neb_ps, m1_neb_pr])
m2_neb = np.array([m2_neb_q200, m2_neb_q500, m2_neb_q850, m2_neb_v200, m2_neb_v500, m2_neb_v850, m2_neb_u200, m2_neb_u500, m2_neb_u850, m2_neb_t200, m2_neb_t500, m2_neb_t850, m2_neb_ps, m2_neb_pr])
m3_neb = np.array([m3_neb_q200, m3_neb_q500, m3_neb_q850, m3_neb_v200, m3_neb_v500, m3_neb_v850, m3_neb_u200, m3_neb_u500, m3_neb_u850, m3_neb_t200, m3_neb_t500, m3_neb_t850, m3_neb_ps, m3_neb_pr])
m4_neb = np.array([m4_neb_q200, m4_neb_q500, m4_neb_q850, m4_neb_v200, m4_neb_v500, m4_neb_v850, m4_neb_u200, m4_neb_u500, m4_neb_u850, m4_neb_t200, m4_neb_t500, m4_neb_t850, m4_neb_ps, m4_neb_pr])

m1_lpb = np.array([m1_lpb_q200, m1_lpb_q500, m1_lpb_q850, m1_lpb_v200, m1_lpb_v500, m1_lpb_v850, m1_lpb_u200, m1_lpb_u500, m1_lpb_u850, m1_lpb_t200, m1_lpb_t500, m1_lpb_t850, m1_lpb_ps, m1_lpb_pr])
m2_lpb = np.array([m2_lpb_q200, m2_lpb_q500, m2_lpb_q850, m2_lpb_v200, m2_lpb_v500, m2_lpb_v850, m2_lpb_u200, m2_lpb_u500, m2_lpb_u850, m2_lpb_t200, m2_lpb_t500, m2_lpb_t850, m2_lpb_ps, m2_lpb_pr])
m3_lpb = np.array([m3_lpb_q200, m3_lpb_q500, m3_lpb_q850, m3_lpb_v200, m3_lpb_v500, m3_lpb_v850, m3_lpb_u200, m3_lpb_u500, m3_lpb_u850, m3_lpb_t200, m3_lpb_t500, m3_lpb_t850, m3_lpb_ps, m3_lpb_pr])
m4_lpb = np.array([m4_lpb_q200, m4_lpb_q500, m4_lpb_q850, m4_lpb_v200, m4_lpb_v500, m4_lpb_v850, m4_lpb_u200, m4_lpb_u500, m4_lpb_u850, m4_lpb_t200, m4_lpb_t500, m4_lpb_t850, m4_lpb_ps, m4_lpb_pr])

# Plot cmip models and obs database 
fig = plt.figure(figsize=(10, 12))

norm_m1 = colors.BoundaryNorm(boundaries=np.arange(0, 1.1, 0.1), ncolors=256)
norm_m2 = colors.BoundaryNorm(boundaries=np.arange(0, 1.1, 0.1), ncolors=256)
norm_m3 = colors.BoundaryNorm(boundaries=np.arange(0, 1.1, 0.1), ncolors=256)
norm_m4 = colors.BoundaryNorm(boundaries=np.arange(0, 2.2, 0.2), ncolors=256)

color_m1 = cm.Blues
color_m2 = cm.Reds
color_m3 = cm.Greys
color_m4 = cm.Greens

xlabels = legend
ylabels = ['Q200', 'Q500', 'Q850', 'V200', 'V500', 'V850', 'U200', 'U500', 'U850', 'T200', 'T500', 'T850', 'SP', 'Pr']
	
ax = fig.add_subplot(5, 4, 1)  
pcm = ax.pcolormesh(m1_namz, edgecolors='white', linewidths=2., norm=norm_m1, cmap=color_m1)
ax.set_title(u'(a)', loc='left', fontweight='bold', fontsize=8)
ax.set_ylabel('NAMZ', fontsize=8)
ax.set_xticks(np.arange(m1_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m1_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 2)  
pcm = ax.pcolormesh(m2_namz, edgecolors='white', linewidths=2., norm=norm_m2, cmap=color_m2)
ax.set_title(u'(b)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m2_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m2_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 3)  
pcm = ax.pcolormesh(m3_namz, edgecolors='white', linewidths=2., norm=norm_m3, cmap=color_m3)
ax.set_title(u'(c)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m3_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m3_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 4)  
pcm = ax.pcolormesh(m4_namz, edgecolors='white', linewidths=2., norm=norm_m4, cmap=color_m4)
ax.set_title(u'(d)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m4_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m4_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 5)  
pcm = ax.pcolormesh(m1_samz, edgecolors='white', linewidths=2., norm=norm_m1, cmap=color_m1)
ax.set_title(u'(e)', loc='left', fontweight='bold', fontsize=8)
ax.set_ylabel('SAMZ', fontsize=8)
ax.set_xticks(np.arange(m1_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m1_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 6)  
pcm = ax.pcolormesh(m2_samz, edgecolors='white', linewidths=2., norm=norm_m2, cmap=color_m2)
ax.set_title(u'(f)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m2_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m2_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 7)  
pcm = ax.pcolormesh(m3_samz, edgecolors='white', linewidths=2., norm=norm_m3, cmap=color_m3)
ax.set_title(u'(g)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m3_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m3_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 8)  
pcm = ax.pcolormesh(m4_samz, edgecolors='white', linewidths=2., norm=norm_m4, cmap=color_m4)
ax.set_title(u'(h)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m4_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m4_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 9)  
pcm = ax.pcolormesh(m1_sam, edgecolors='white', linewidths=2., norm=norm_m1, cmap=color_m1)
ax.set_title(u'(i)', loc='left', fontweight='bold', fontsize=8)
ax.set_ylabel('SAM', fontsize=8)
ax.set_xticks(np.arange(m1_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m1_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 10)  
pcm = ax.pcolormesh(m2_sam, edgecolors='white', linewidths=2., norm=norm_m2, cmap=color_m2)
ax.set_title(u'(j)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m2_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m2_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 11)  
pcm = ax.pcolormesh(m3_sam, edgecolors='white', linewidths=2., norm=norm_m3, cmap=color_m3)
ax.set_title(u'(k)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m3_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m3_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 12)  
pcm = ax.pcolormesh(m4_sam, edgecolors='white', linewidths=2., norm=norm_m4, cmap=color_m4)
ax.set_title(u'(l)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m4_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m4_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 13)  
pcm = ax.pcolormesh(m1_neb, edgecolors='white', linewidths=2., norm=norm_m1, cmap=color_m1)
ax.set_title(u'(m)', loc='left', fontweight='bold', fontsize=8)
ax.set_ylabel('NEB', fontsize=8)
ax.set_xticks(np.arange(m1_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m1_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 14)  
pcm = ax.pcolormesh(m2_neb, edgecolors='white', linewidths=2., norm=norm_m2, cmap=color_m2)
ax.set_title(u'(n)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m2_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m2_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 15)  
pcm = ax.pcolormesh(m3_neb, edgecolors='white', linewidths=2., norm=norm_m3, cmap=color_m3)
ax.set_title(u'(o)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m3_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m3_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 16)  
pcm = ax.pcolormesh(m4_neb, edgecolors='white', linewidths=2., norm=norm_m4, cmap=color_m4)
ax.set_title(u'(p)', loc='left', fontweight='bold', fontsize=8)
ax.set_xticks(np.arange(m4_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m4_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_xticklabels(), visible=False)
plt.setp(ax.get_yticklabels(), visible=False)

ax = fig.add_subplot(5, 4, 17)  
pcm = ax.pcolormesh(m1_lpb, edgecolors='white', linewidths=2., norm=norm_m1, cmap=color_m1)
ax.set_title(u'(q)', loc='left', fontweight='bold', fontsize=8)
ax.set_xlabel('NRMSE', fontsize=8)
ax.set_ylabel('LPB', fontsize=8)
ax.set_xticks(np.arange(m1_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m1_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
clb = fig.colorbar(pcm, cax=fig.add_axes([0.12, 0.06, 0.18, 0.01]), orientation='horizontal')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 4, 18)  
pcm = ax.pcolormesh(m2_lpb, edgecolors='white', linewidths=2., norm=norm_m2, cmap=color_m2)
ax.set_title(u'(r)', loc='left', fontweight='bold', fontsize=8)
ax.set_xlabel('TSS', fontsize=8)
ax.set_xticks(np.arange(m2_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m2_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, cax=fig.add_axes([0.32, 0.06, 0.18, 0.01]), orientation='horizontal')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 4, 19)  
pcm = ax.pcolormesh(m3_lpb, edgecolors='white', linewidths=2., norm=norm_m3, cmap=color_m3)
ax.set_title(u'(s)', loc='left', fontweight='bold', fontsize=8)
ax.set_xlabel('R-squared', fontsize=8)
ax.set_xticks(np.arange(m3_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m3_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, cax=fig.add_axes([0.52, 0.06, 0.18, 0.01]), orientation='horizontal')
clb.ax.tick_params(labelsize=8)

ax = fig.add_subplot(5, 4, 20)  
pcm = ax.pcolormesh(m4_lpb, edgecolors='white', linewidths=2., norm=norm_m4, cmap=color_m4)
ax.set_title(u'(t)', loc='left', fontweight='bold', fontsize=8)
ax.set_xlabel('IVS', fontsize=8)
ax.set_xticks(np.arange(m4_namz.shape[1]) + 0.5)
ax.set_yticks(np.arange(m4_namz.shape[0]) + 0.5)
ax.set_xticklabels(xlabels, fontsize=8, rotation=90)
ax.set_yticklabels(ylabels, fontsize=8)
plt.setp(ax.get_yticklabels(), visible=False)
clb = fig.colorbar(pcm, cax=fig.add_axes([0.72, 0.06, 0.18, 0.01]), orientation='horizontal')
clb.ax.tick_params(labelsize=8)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_portrait_stats_metrics_cmip6_obs_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()






