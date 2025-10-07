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
	
	mean_850_ts = np.nanmean(np.nanmean(var[:][:,2,:,:], axis=1), axis=1)
	mean_500_ts = np.nanmean(np.nanmean(var[:][:,1,:,:], axis=1), axis=1)
	mean_200_ts = np.nanmean(np.nanmean(var[:][:,0,:,:], axis=1), axis=1)

	xy_ts_850 = latlon(mean_850)
	xy_ts_500 = latlon(mean_500)
	xy_ts_200 = latlon(mean_200)

	ac_ts_850 = annual_cycle(mean_850_ts)
	ac_ts_500 = annual_cycle(mean_500_ts)
	ac_ts_200 = annual_cycle(mean_200_ts)

	yr_ts_850 = year_mean(mean_850_ts)
	yr_ts_500 = year_mean(mean_850_ts)
	yr_ts_200 = year_mean(mean_850_ts)
		
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
	      

def compute_cri(rank1,rank2,rank3,rank4):
	
	p1 = (rank1+rank2+rank3+rank4)
	p2 = p1/68
	cri = 1 - p2
			
	return cri
	

def sort_list(data_list):
	
	li = []
	for i in range(len(data_list)):
		li.append([data_list[i], i])
	
	li.sort()
	sort_index = []
	for x in li:
		sort_index.append(x[1])

	model_list = []
	value_list = []
	for ii in sort_index:
		model_list.append(cmip6[ii+1][0])
		value_list.append(data_list[ii])
	
	data_argsort = np.argsort(model_list)

	data_argsort_i = []
	for idx in data_argsort:
		data_argsort_i.append(idx+1)

	return data_argsort_i
	

def sort_list_reverse(data_list):
	
	li = []
	for i in range(len(data_list)):
		li.append([data_list[i], i])
	  
	li.sort(reverse=True)
	sort_index = []
	for x in li:
		sort_index.append(x[1])

	model_list = []
	value_list = []
	for ii in sort_index:
		model_list.append(cmip6[ii+1][0])
		value_list.append(data_list[ii])

	data_argsort = np.argsort(model_list)
	
	data_argsort_ii = []
	for idx in data_argsort:
		data_argsort_ii.append(idx+1)
	
	return data_argsort_ii
	
	
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
m1_namz_pr, m1_namz_ps, m1_namz_t850, m1_namz_t500, m1_namz_t200, m1_namz_u850, m1_namz_u500, m1_namz_u200, m1_namz_v850, m1_namz_v500, m1_namz_v200, m1_namz_q850, m1_namz_q500, m1_namz_q200 = [[] for _ in range(14)]
m2_namz_pr, m2_namz_ps, m2_namz_t850, m2_namz_t500, m2_namz_t200, m2_namz_u850, m2_namz_u500, m2_namz_u200, m2_namz_v850, m2_namz_v500, m2_namz_v200, m2_namz_q850, m2_namz_q500, m2_namz_q200 = [[] for _ in range(14)]
m3_namz_pr, m3_namz_ps, m3_namz_t850, m3_namz_t500, m3_namz_t200, m3_namz_u850, m3_namz_u500, m3_namz_u200, m3_namz_v850, m3_namz_v500, m3_namz_v200, m3_namz_q850, m3_namz_q500, m3_namz_q200 = [[] for _ in range(14)]
m4_namz_pr, m4_namz_ps, m4_namz_t850, m4_namz_t500, m4_namz_t200, m4_namz_u850, m4_namz_u500, m4_namz_u200, m4_namz_v850, m4_namz_v500, m4_namz_v200, m4_namz_q850, m4_namz_q500, m4_namz_q200 = [[] for _ in range(14)]

m1_samz_pr, m1_samz_ps, m1_samz_t850, m1_samz_t500, m1_samz_t200, m1_samz_u850, m1_samz_u500, m1_samz_u200, m1_samz_v850, m1_samz_v500, m1_samz_v200, m1_samz_q850, m1_samz_q500, m1_samz_q200 = [[] for _ in range(14)]
m2_samz_pr, m2_samz_ps, m2_samz_t850, m2_samz_t500, m2_samz_t200, m2_samz_u850, m2_samz_u500, m2_samz_u200, m2_samz_v850, m2_samz_v500, m2_samz_v200, m2_samz_q850, m2_samz_q500, m2_samz_q200 = [[] for _ in range(14)]
m3_samz_pr, m3_samz_ps, m3_samz_t850, m3_samz_t500, m3_samz_t200, m3_samz_u850, m3_samz_u500, m3_samz_u200, m3_samz_v850, m3_samz_v500, m3_samz_v200, m3_samz_q850, m3_samz_q500, m3_samz_q200 = [[] for _ in range(14)]
m4_samz_pr, m4_samz_ps, m4_samz_t850, m4_samz_t500, m4_samz_t200, m4_samz_u850, m4_samz_u500, m4_samz_u200, m4_samz_v850, m4_samz_v500, m4_samz_v200, m4_samz_q850, m4_samz_q500, m4_samz_q200 = [[] for _ in range(14)]

m1_sam_pr, m1_sam_ps, m1_sam_t850, m1_sam_t500, m1_sam_t200, m1_sam_u850, m1_sam_u500, m1_sam_u200, m1_sam_v850, m1_sam_v500, m1_sam_v200, m1_sam_q850, m1_sam_q500, m1_sam_q200 = [[] for _ in range(14)]
m2_sam_pr, m2_sam_ps, m2_sam_t850, m2_sam_t500, m2_sam_t200, m2_sam_u850, m2_sam_u500, m2_sam_u200, m2_sam_v850, m2_sam_v500, m2_sam_v200, m2_sam_q850, m2_sam_q500, m2_sam_q200 = [[] for _ in range(14)]
m3_sam_pr, m3_sam_ps, m3_sam_t850, m3_sam_t500, m3_sam_t200, m3_sam_u850, m3_sam_u500, m3_sam_u200, m3_sam_v850, m3_sam_v500, m3_sam_v200, m3_sam_q850, m3_sam_q500, m3_sam_q200 = [[] for _ in range(14)]
m4_sam_pr, m4_sam_ps, m4_sam_t850, m4_sam_t500, m4_sam_t200, m4_sam_u850, m4_sam_u500, m4_sam_u200, m4_sam_v850, m4_sam_v500, m4_sam_v200, m4_sam_q850, m4_sam_q500, m4_sam_q200 = [[] for _ in range(14)]

m1_neb_pr, m1_neb_ps, m1_neb_t850, m1_neb_t500, m1_neb_t200, m1_neb_u850, m1_neb_u500, m1_neb_u200, m1_neb_v850, m1_neb_v500, m1_neb_v200, m1_neb_q850, m1_neb_q500, m1_neb_q200 = [[] for _ in range(14)]
m2_neb_pr, m2_neb_ps, m2_neb_t850, m2_neb_t500, m2_neb_t200, m2_neb_u850, m2_neb_u500, m2_neb_u200, m2_neb_v850, m2_neb_v500, m2_neb_v200, m2_neb_q850, m2_neb_q500, m2_neb_q200 = [[] for _ in range(14)]
m3_neb_pr, m3_neb_ps, m3_neb_t850, m3_neb_t500, m3_neb_t200, m3_neb_u850, m3_neb_u500, m3_neb_u200, m3_neb_v850, m3_neb_v500, m3_neb_v200, m3_neb_q850, m3_neb_q500, m3_neb_q200 = [[] for _ in range(14)]
m4_neb_pr, m4_neb_ps, m4_neb_t850, m4_neb_t500, m4_neb_t200, m4_neb_u850, m4_neb_u500, m4_neb_u200, m4_neb_v850, m4_neb_v500, m4_neb_v200, m4_neb_q850, m4_neb_q500, m4_neb_q200 = [[] for _ in range(14)]

m1_lpb_pr, m1_lpb_ps, m1_lpb_t850, m1_lpb_t500, m1_lpb_t200, m1_lpb_u850, m1_lpb_u500, m1_lpb_u200, m1_lpb_v850, m1_lpb_v500, m1_lpb_v200, m1_lpb_q850, m1_lpb_q500, m1_lpb_q200 = [[] for _ in range(14)]
m2_lpb_pr, m2_lpb_ps, m2_lpb_t850, m2_lpb_t500, m2_lpb_t200, m2_lpb_u850, m2_lpb_u500, m2_lpb_u200, m2_lpb_v850, m2_lpb_v500, m2_lpb_v200, m2_lpb_q850, m2_lpb_q500, m2_lpb_q200 = [[] for _ in range(14)]
m3_lpb_pr, m3_lpb_ps, m3_lpb_t850, m3_lpb_t500, m3_lpb_t200, m3_lpb_u850, m3_lpb_u500, m3_lpb_u200, m3_lpb_v850, m3_lpb_v500, m3_lpb_v200, m3_lpb_q850, m3_lpb_q500, m3_lpb_q200 = [[] for _ in range(14)]
m4_lpb_pr, m4_lpb_ps, m4_lpb_t850, m4_lpb_t500, m4_lpb_t200, m4_lpb_u850, m4_lpb_u500, m4_lpb_u200, m4_lpb_v850, m4_lpb_v500, m4_lpb_v200, m4_lpb_q850, m4_lpb_q500, m4_lpb_q200 = [[] for _ in range(14)]


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

# Sort NAMZ
sort_m1_namz_pr = sort_list(m1_namz_pr)
sort_m2_namz_pr = sort_list_reverse(m2_namz_pr)
sort_m3_namz_pr = sort_list_reverse(m3_namz_pr)
sort_m4_namz_pr = sort_list(m4_namz_pr)

sort_m1_namz_ps = sort_list(m1_namz_ps)
sort_m2_namz_ps = sort_list_reverse(m2_namz_ps)
sort_m3_namz_ps = sort_list_reverse(m3_namz_ps)
sort_m4_namz_ps = sort_list(m4_namz_ps)

sort_m1_namz_t850 = sort_list(m1_namz_t850)
sort_m2_namz_t850 = sort_list_reverse(m2_namz_t850)
sort_m3_namz_t850 = sort_list_reverse(m3_namz_t850)
sort_m4_namz_t850 = sort_list(m4_namz_t850)

sort_m1_namz_t500 = sort_list(m1_namz_t500)
sort_m2_namz_t500 = sort_list_reverse(m2_namz_t500)
sort_m3_namz_t500 = sort_list_reverse(m3_namz_t500)
sort_m4_namz_t500 = sort_list(m4_namz_t500)

sort_m1_namz_t200 = sort_list(m1_namz_t200)
sort_m2_namz_t200 = sort_list_reverse(m2_namz_t200)
sort_m3_namz_t200 = sort_list_reverse(m3_namz_t200)
sort_m4_namz_t200 = sort_list(m4_namz_t200)

sort_m1_namz_u850 = sort_list(m1_namz_u850)
sort_m2_namz_u850 = sort_list_reverse(m2_namz_u850)
sort_m3_namz_u850 = sort_list_reverse(m3_namz_u850)
sort_m4_namz_u850 = sort_list(m4_namz_u850)

sort_m1_namz_u500 = sort_list(m1_namz_u500)
sort_m2_namz_u500 = sort_list_reverse(m2_namz_u500)
sort_m3_namz_u500 = sort_list_reverse(m3_namz_u500)
sort_m4_namz_u500 = sort_list(m4_namz_u500)

sort_m1_namz_u200 = sort_list(m1_namz_u200)
sort_m2_namz_u200 = sort_list_reverse(m2_namz_u200)
sort_m3_namz_u200 = sort_list_reverse(m3_namz_u200)
sort_m4_namz_u200 = sort_list(m4_namz_u200)

sort_m1_namz_v850 = sort_list(m1_namz_v850)
sort_m2_namz_v850 = sort_list_reverse(m2_namz_v850)
sort_m3_namz_v850 = sort_list_reverse(m3_namz_v850)
sort_m4_namz_v850 = sort_list(m4_namz_v850)

sort_m1_namz_v500 = sort_list(m1_namz_v500)
sort_m2_namz_v500 = sort_list_reverse(m2_namz_v500)
sort_m3_namz_v500 = sort_list_reverse(m3_namz_v500)
sort_m4_namz_v500 = sort_list(m4_namz_v500)

sort_m1_namz_v200 = sort_list(m1_namz_v200)
sort_m2_namz_v200 = sort_list_reverse(m2_namz_v200)
sort_m3_namz_v200 = sort_list_reverse(m3_namz_v200)
sort_m4_namz_v200 = sort_list(m4_namz_v200)

sort_m1_namz_q850 = sort_list(m1_namz_q850)
sort_m2_namz_q850 = sort_list_reverse(m2_namz_q850)
sort_m3_namz_q850 = sort_list_reverse(m3_namz_q850)
sort_m4_namz_q850 = sort_list(m4_namz_q850)

sort_m1_namz_q500 = sort_list(m1_namz_q500)
sort_m2_namz_q500 = sort_list_reverse(m2_namz_q500)
sort_m3_namz_q500 = sort_list_reverse(m3_namz_q500)
sort_m4_namz_q500 = sort_list(m4_namz_q500)

sort_m1_namz_q200 = sort_list(m1_namz_q200)
sort_m2_namz_q200 = sort_list_reverse(m2_namz_q200)
sort_m3_namz_q200 = sort_list_reverse(m3_namz_q200)
sort_m4_namz_q200 = sort_list(m4_namz_q200)

# Sort SAMZ
sort_m1_samz_pr = sort_list(m1_samz_pr)
sort_m2_samz_pr = sort_list_reverse(m2_samz_pr)
sort_m3_samz_pr = sort_list_reverse(m3_samz_pr)
sort_m4_samz_pr = sort_list(m4_samz_pr)

sort_m1_samz_ps = sort_list(m1_samz_ps)
sort_m2_samz_ps = sort_list_reverse(m2_samz_ps)
sort_m3_samz_ps = sort_list_reverse(m3_samz_ps)
sort_m4_samz_ps = sort_list(m4_samz_ps)

sort_m1_samz_t850 = sort_list(m1_samz_t850)
sort_m2_samz_t850 = sort_list_reverse(m2_samz_t850)
sort_m3_samz_t850 = sort_list_reverse(m3_samz_t850)
sort_m4_samz_t850 = sort_list(m4_samz_t850)

sort_m1_samz_t500 = sort_list(m1_samz_t500)
sort_m2_samz_t500 = sort_list_reverse(m2_samz_t500)
sort_m3_samz_t500 = sort_list_reverse(m3_samz_t500)
sort_m4_samz_t500 = sort_list(m4_samz_t500)

sort_m1_samz_t200 = sort_list(m1_samz_t200)
sort_m2_samz_t200 = sort_list_reverse(m2_samz_t200)
sort_m3_samz_t200 = sort_list_reverse(m3_samz_t200)
sort_m4_samz_t200 = sort_list(m4_samz_t200)

sort_m1_samz_u850 = sort_list(m1_samz_u850)
sort_m2_samz_u850 = sort_list_reverse(m2_samz_u850)
sort_m3_samz_u850 = sort_list_reverse(m3_samz_u850)
sort_m4_samz_u850 = sort_list(m4_samz_u850)

sort_m1_samz_u500 = sort_list(m1_samz_u500)
sort_m2_samz_u500 = sort_list_reverse(m2_samz_u500)
sort_m3_samz_u500 = sort_list_reverse(m3_samz_u500)
sort_m4_samz_u500 = sort_list(m4_samz_u500)

sort_m1_samz_u200 = sort_list(m1_samz_u200)
sort_m2_samz_u200 = sort_list_reverse(m2_samz_u200)
sort_m3_samz_u200 = sort_list_reverse(m3_samz_u200)
sort_m4_samz_u200 = sort_list(m4_samz_u200)

sort_m1_samz_v850 = sort_list(m1_samz_v850)
sort_m2_samz_v850 = sort_list_reverse(m2_samz_v850)
sort_m3_samz_v850 = sort_list_reverse(m3_samz_v850)
sort_m4_samz_v850 = sort_list(m4_samz_v850)

sort_m1_samz_v500 = sort_list(m1_samz_v500)
sort_m2_samz_v500 = sort_list_reverse(m2_samz_v500)
sort_m3_samz_v500 = sort_list_reverse(m3_samz_v500)
sort_m4_samz_v500 = sort_list(m4_samz_v500)

sort_m1_samz_v200 = sort_list(m1_samz_v200)
sort_m2_samz_v200 = sort_list_reverse(m2_samz_v200)
sort_m3_samz_v200 = sort_list_reverse(m3_samz_v200)
sort_m4_samz_v200 = sort_list(m4_samz_v200)

sort_m1_samz_q850 = sort_list(m1_samz_q850)
sort_m2_samz_q850 = sort_list_reverse(m2_samz_q850)
sort_m3_samz_q850 = sort_list_reverse(m3_samz_q850)
sort_m4_samz_q850 = sort_list(m4_samz_q850)

sort_m1_samz_q500 = sort_list(m1_samz_q500)
sort_m2_samz_q500 = sort_list_reverse(m2_samz_q500)
sort_m3_samz_q500 = sort_list_reverse(m3_samz_q500)
sort_m4_samz_q500 = sort_list(m4_samz_q500)

sort_m1_samz_q200 = sort_list(m1_samz_q200)
sort_m2_samz_q200 = sort_list_reverse(m2_samz_q200)
sort_m3_samz_q200 = sort_list_reverse(m3_samz_q200)
sort_m4_samz_q200 = sort_list(m4_samz_q200)

# Sort SAM
sort_m1_sam_pr = sort_list(m1_sam_pr)
sort_m2_sam_pr = sort_list_reverse(m2_sam_pr)
sort_m3_sam_pr = sort_list_reverse(m3_sam_pr)
sort_m4_sam_pr = sort_list(m4_sam_pr)

sort_m1_sam_ps = sort_list(m1_sam_ps)
sort_m2_sam_ps = sort_list_reverse(m2_sam_ps)
sort_m3_sam_ps = sort_list_reverse(m3_sam_ps)
sort_m4_sam_ps = sort_list(m4_sam_ps)

sort_m1_sam_t850 = sort_list(m1_sam_t850)
sort_m2_sam_t850 = sort_list_reverse(m2_sam_t850)
sort_m3_sam_t850 = sort_list_reverse(m3_sam_t850)
sort_m4_sam_t850 = sort_list(m4_sam_t850)

sort_m1_sam_t500 = sort_list(m1_sam_t500)
sort_m2_sam_t500 = sort_list_reverse(m2_sam_t500)
sort_m3_sam_t500 = sort_list_reverse(m3_sam_t500)
sort_m4_sam_t500 = sort_list(m4_sam_t500)

sort_m1_sam_t200 = sort_list(m1_sam_t200)
sort_m2_sam_t200 = sort_list_reverse(m2_sam_t200)
sort_m3_sam_t200 = sort_list_reverse(m3_sam_t200)
sort_m4_sam_t200 = sort_list(m4_sam_t200)

sort_m1_sam_u850 = sort_list(m1_sam_u850)
sort_m2_sam_u850 = sort_list_reverse(m2_sam_u850)
sort_m3_sam_u850 = sort_list_reverse(m3_sam_u850)
sort_m4_sam_u850 = sort_list(m4_sam_u850)

sort_m1_sam_u500 = sort_list(m1_sam_u500)
sort_m2_sam_u500 = sort_list_reverse(m2_sam_u500)
sort_m3_sam_u500 = sort_list_reverse(m3_sam_u500)
sort_m4_sam_u500 = sort_list(m4_sam_u500)

sort_m1_sam_u200 = sort_list(m1_sam_u200)
sort_m2_sam_u200 = sort_list_reverse(m2_sam_u200)
sort_m3_sam_u200 = sort_list_reverse(m3_sam_u200)
sort_m4_sam_u200 = sort_list(m4_sam_u200)

sort_m1_sam_v850 = sort_list(m1_sam_v850)
sort_m2_sam_v850 = sort_list_reverse(m2_sam_v850)
sort_m3_sam_v850 = sort_list_reverse(m3_sam_v850)
sort_m4_sam_v850 = sort_list(m4_sam_v850)

sort_m1_sam_v500 = sort_list(m1_sam_v500)
sort_m2_sam_v500 = sort_list_reverse(m2_sam_v500)
sort_m3_sam_v500 = sort_list_reverse(m3_sam_v500)
sort_m4_sam_v500 = sort_list(m4_sam_v500)

sort_m1_sam_v200 = sort_list(m1_sam_v200)
sort_m2_sam_v200 = sort_list_reverse(m2_sam_v200)
sort_m3_sam_v200 = sort_list_reverse(m3_sam_v200)
sort_m4_sam_v200 = sort_list(m4_sam_v200)

sort_m1_sam_q850 = sort_list(m1_sam_q850)
sort_m2_sam_q850 = sort_list_reverse(m2_sam_q850)
sort_m3_sam_q850 = sort_list_reverse(m3_sam_q850)
sort_m4_sam_q850 = sort_list(m4_sam_q850)

sort_m1_sam_q500 = sort_list(m1_sam_q500)
sort_m2_sam_q500 = sort_list_reverse(m2_sam_q500)
sort_m3_sam_q500 = sort_list_reverse(m3_sam_q500)
sort_m4_sam_q500 = sort_list(m4_sam_q500)

sort_m1_sam_q200 = sort_list(m1_sam_q200)
sort_m2_sam_q200 = sort_list_reverse(m2_sam_q200)
sort_m3_sam_q200 = sort_list_reverse(m3_sam_q200)
sort_m4_sam_q200 = sort_list(m4_sam_q200)

# Sort NEB
sort_m1_neb_pr = sort_list(m1_neb_pr)
sort_m2_neb_pr = sort_list_reverse(m2_neb_pr)
sort_m3_neb_pr = sort_list_reverse(m3_neb_pr)
sort_m4_neb_pr = sort_list(m4_neb_pr)

sort_m1_neb_ps = sort_list(m1_neb_ps)
sort_m2_neb_ps = sort_list_reverse(m2_neb_ps)
sort_m3_neb_ps = sort_list_reverse(m3_neb_ps)
sort_m4_neb_ps = sort_list(m4_neb_ps)

sort_m1_neb_t850 = sort_list(m1_neb_t850)
sort_m2_neb_t850 = sort_list_reverse(m2_neb_t850)
sort_m3_neb_t850 = sort_list_reverse(m3_neb_t850)
sort_m4_neb_t850 = sort_list(m4_neb_t850)

sort_m1_neb_t500 = sort_list(m1_neb_t500)
sort_m2_neb_t500 = sort_list_reverse(m2_neb_t500)
sort_m3_neb_t500 = sort_list_reverse(m3_neb_t500)
sort_m4_neb_t500 = sort_list(m4_neb_t500)

sort_m1_neb_t200 = sort_list(m1_neb_t200)
sort_m2_neb_t200 = sort_list_reverse(m2_neb_t200)
sort_m3_neb_t200 = sort_list_reverse(m3_neb_t200)
sort_m4_neb_t200 = sort_list(m4_neb_t200)

sort_m1_neb_u850 = sort_list(m1_neb_u850)
sort_m2_neb_u850 = sort_list_reverse(m2_neb_u850)
sort_m3_neb_u850 = sort_list_reverse(m3_neb_u850)
sort_m4_neb_u850 = sort_list(m4_neb_u850)

sort_m1_neb_u500 = sort_list(m1_neb_u500)
sort_m2_neb_u500 = sort_list_reverse(m2_neb_u500)
sort_m3_neb_u500 = sort_list_reverse(m3_neb_u500)
sort_m4_neb_u500 = sort_list(m4_neb_u500)

sort_m1_neb_u200 = sort_list(m1_neb_u200)
sort_m2_neb_u200 = sort_list_reverse(m2_neb_u200)
sort_m3_neb_u200 = sort_list_reverse(m3_neb_u200)
sort_m4_neb_u200 = sort_list(m4_neb_u200)

sort_m1_neb_v850 = sort_list(m1_neb_v850)
sort_m2_neb_v850 = sort_list_reverse(m2_neb_v850)
sort_m3_neb_v850 = sort_list_reverse(m3_neb_v850)
sort_m4_neb_v850 = sort_list(m4_neb_v850)

sort_m1_neb_v500 = sort_list(m1_neb_v500)
sort_m2_neb_v500 = sort_list_reverse(m2_neb_v500)
sort_m3_neb_v500 = sort_list_reverse(m3_neb_v500)
sort_m4_neb_v500 = sort_list(m4_neb_v500)

sort_m1_neb_v200 = sort_list(m1_neb_v200)
sort_m2_neb_v200 = sort_list_reverse(m2_neb_v200)
sort_m3_neb_v200 = sort_list_reverse(m3_neb_v200)
sort_m4_neb_v200 = sort_list(m4_neb_v200)

sort_m1_neb_q850 = sort_list(m1_neb_q850)
sort_m2_neb_q850 = sort_list_reverse(m2_neb_q850)
sort_m3_neb_q850 = sort_list_reverse(m3_neb_q850)
sort_m4_neb_q850 = sort_list(m4_neb_q850)

sort_m1_neb_q500 = sort_list(m1_neb_q500)
sort_m2_neb_q500 = sort_list_reverse(m2_neb_q500)
sort_m3_neb_q500 = sort_list_reverse(m3_neb_q500)
sort_m4_neb_q500 = sort_list(m4_neb_q500)

sort_m1_neb_q200 = sort_list(m1_neb_q200)
sort_m2_neb_q200 = sort_list_reverse(m2_neb_q200)
sort_m3_neb_q200 = sort_list_reverse(m3_neb_q200)
sort_m4_neb_q200 = sort_list(m4_neb_q200)

# Sort LPB
sort_m1_lpb_pr = sort_list(m1_lpb_pr)
sort_m2_lpb_pr = sort_list_reverse(m2_lpb_pr)
sort_m3_lpb_pr = sort_list_reverse(m3_lpb_pr)
sort_m4_lpb_pr = sort_list(m4_lpb_pr)

sort_m1_lpb_ps = sort_list(m1_lpb_ps)
sort_m2_lpb_ps = sort_list_reverse(m2_lpb_ps)
sort_m3_lpb_ps = sort_list_reverse(m3_lpb_ps)
sort_m4_lpb_ps = sort_list(m4_lpb_ps)

sort_m1_lpb_t850 = sort_list(m1_lpb_t850)
sort_m2_lpb_t850 = sort_list_reverse(m2_lpb_t850)
sort_m3_lpb_t850 = sort_list_reverse(m3_lpb_t850)
sort_m4_lpb_t850 = sort_list(m4_lpb_t850)

sort_m1_lpb_t500 = sort_list(m1_lpb_t500)
sort_m2_lpb_t500 = sort_list_reverse(m2_lpb_t500)
sort_m3_lpb_t500 = sort_list_reverse(m3_lpb_t500)
sort_m4_lpb_t500 = sort_list(m4_lpb_t500)

sort_m1_lpb_t200 = sort_list(m1_lpb_t200)
sort_m2_lpb_t200 = sort_list_reverse(m2_lpb_t200)
sort_m3_lpb_t200 = sort_list_reverse(m3_lpb_t200)
sort_m4_lpb_t200 = sort_list(m4_lpb_t200)

sort_m1_lpb_u850 = sort_list(m1_lpb_u850)
sort_m2_lpb_u850 = sort_list_reverse(m2_lpb_u850)
sort_m3_lpb_u850 = sort_list_reverse(m3_lpb_u850)
sort_m4_lpb_u850 = sort_list(m4_lpb_u850)

sort_m1_lpb_u500 = sort_list(m1_lpb_u500)
sort_m2_lpb_u500 = sort_list_reverse(m2_lpb_u500)
sort_m3_lpb_u500 = sort_list_reverse(m3_lpb_u500)
sort_m4_lpb_u500 = sort_list(m4_lpb_u500)

sort_m1_lpb_u200 = sort_list(m1_lpb_u200)
sort_m2_lpb_u200 = sort_list_reverse(m2_lpb_u200)
sort_m3_lpb_u200 = sort_list_reverse(m3_lpb_u200)
sort_m4_lpb_u200 = sort_list(m4_lpb_u200)

sort_m1_lpb_v850 = sort_list(m1_lpb_v850)
sort_m2_lpb_v850 = sort_list_reverse(m2_lpb_v850)
sort_m3_lpb_v850 = sort_list_reverse(m3_lpb_v850)
sort_m4_lpb_v850 = sort_list(m4_lpb_v850)

sort_m1_lpb_v500 = sort_list(m1_lpb_v500)
sort_m2_lpb_v500 = sort_list_reverse(m2_lpb_v500)
sort_m3_lpb_v500 = sort_list_reverse(m3_lpb_v500)
sort_m4_lpb_v500 = sort_list(m4_lpb_v500)

sort_m1_lpb_v200 = sort_list(m1_lpb_v200)
sort_m2_lpb_v200 = sort_list_reverse(m2_lpb_v200)
sort_m3_lpb_v200 = sort_list_reverse(m3_lpb_v200)
sort_m4_lpb_v200 = sort_list(m4_lpb_v200)

sort_m1_lpb_q850 = sort_list(m1_lpb_q850)
sort_m2_lpb_q850 = sort_list_reverse(m2_lpb_q850)
sort_m3_lpb_q850 = sort_list_reverse(m3_lpb_q850)
sort_m4_lpb_q850 = sort_list(m4_lpb_q850)

sort_m1_lpb_q500 = sort_list(m1_lpb_q500)
sort_m2_lpb_q500 = sort_list_reverse(m2_lpb_q500)
sort_m3_lpb_q500 = sort_list_reverse(m3_lpb_q500)
sort_m4_lpb_q500 = sort_list(m4_lpb_q500)

sort_m1_lpb_q200 = sort_list(m1_lpb_q200)
sort_m2_lpb_q200 = sort_list_reverse(m2_lpb_q200)
sort_m3_lpb_q200 = sort_list_reverse(m3_lpb_q200)
sort_m4_lpb_q200 = sort_list(m4_lpb_q200)

cri_namz_pr, cri_namz_ps, cri_namz_t850, cri_namz_t500, cri_namz_t200, cri_namz_u850, cri_namz_u500, cri_namz_u200, cri_namz_v850, cri_namz_v500, cri_namz_v200, cri_namz_q850, cri_namz_q500, cri_namz_q200 = [[] for _ in range(14)]
cri_samz_pr, cri_samz_ps, cri_samz_t850, cri_samz_t500, cri_samz_t200, cri_samz_u850, cri_samz_u500, cri_samz_u200, cri_samz_v850, cri_samz_v500, cri_samz_v200, cri_samz_q850, cri_samz_q500, cri_samz_q200 = [[] for _ in range(14)]
cri_sam_pr, cri_sam_ps, cri_sam_t850, cri_sam_t500, cri_sam_t200, cri_sam_u850, cri_sam_u500, cri_sam_u200, cri_sam_v850, cri_sam_v500, cri_sam_v200, cri_sam_q850, cri_sam_q500, cri_sam_q200 = [[] for _ in range(14)]
cri_neb_pr, cri_neb_ps, cri_neb_t850, cri_neb_t500, cri_neb_t200, cri_neb_u850, cri_neb_u500, cri_neb_u200, cri_neb_v850, cri_neb_v500, cri_neb_v200, cri_neb_q850, cri_neb_q500, cri_neb_q200 = [[] for _ in range(14)]
cri_lpb_pr, cri_lpb_ps, cri_lpb_t850, cri_lpb_t500, cri_lpb_t200, cri_lpb_u850, cri_lpb_u500, cri_lpb_u200, cri_lpb_v850, cri_lpb_v500, cri_lpb_v200, cri_lpb_q850, cri_lpb_q500, cri_lpb_q200 = [[] for _ in range(14)]

for i in range(0, 17):

	cri_namz_pr.append(compute_cri(sort_m1_namz_pr[i], sort_m2_namz_pr[i], sort_m3_namz_pr[i], sort_m4_namz_pr[i]))
	cri_namz_ps.append(compute_cri(sort_m1_namz_ps[i], sort_m2_namz_ps[i], sort_m3_namz_ps[i], sort_m4_namz_ps[i]))
	cri_namz_t850.append(compute_cri(sort_m1_namz_t850[i], sort_m2_namz_t850[i], sort_m3_namz_t850[i], sort_m4_namz_t850[i]))
	cri_namz_t500.append(compute_cri(sort_m1_namz_t500[i], sort_m2_namz_t500[i], sort_m3_namz_t500[i], sort_m4_namz_t500[i]))
	cri_namz_t200.append(compute_cri(sort_m1_namz_t200[i], sort_m2_namz_t200[i], sort_m3_namz_t200[i], sort_m4_namz_t200[i]))
	cri_namz_u850.append(compute_cri(sort_m1_namz_u850[i], sort_m2_namz_u850[i], sort_m3_namz_u850[i], sort_m4_namz_u850[i]))
	cri_namz_u500.append(compute_cri(sort_m1_namz_u500[i], sort_m2_namz_u500[i], sort_m3_namz_u500[i], sort_m4_namz_u500[i]))
	cri_namz_u200.append(compute_cri(sort_m1_namz_u200[i], sort_m2_namz_u200[i], sort_m3_namz_u200[i], sort_m4_namz_u200[i]))
	cri_namz_v850.append(compute_cri(sort_m1_namz_v850[i], sort_m2_namz_v850[i], sort_m3_namz_v850[i], sort_m4_namz_v850[i]))
	cri_namz_v500.append(compute_cri(sort_m1_namz_v500[i], sort_m2_namz_v500[i], sort_m3_namz_v500[i], sort_m4_namz_v500[i]))
	cri_namz_v200.append(compute_cri(sort_m1_namz_v200[i], sort_m2_namz_v200[i], sort_m3_namz_v200[i], sort_m4_namz_v200[i]))
	cri_namz_q850.append(compute_cri(sort_m1_namz_q850[i], sort_m2_namz_q850[i], sort_m3_namz_q850[i], sort_m4_namz_q850[i]))
	cri_namz_q500.append(compute_cri(sort_m1_namz_q500[i], sort_m2_namz_q500[i], sort_m3_namz_q500[i], sort_m4_namz_q500[i]))
	cri_namz_q200.append(compute_cri(sort_m1_namz_q200[i], sort_m2_namz_q200[i], sort_m3_namz_q200[i], sort_m4_namz_q200[i]))
	
	cri_samz_pr.append(compute_cri(sort_m1_samz_pr[i], sort_m2_samz_pr[i], sort_m3_samz_pr[i], sort_m4_samz_pr[i]))
	cri_samz_ps.append(compute_cri(sort_m1_samz_ps[i], sort_m2_samz_ps[i], sort_m3_samz_ps[i], sort_m4_samz_ps[i]))
	cri_samz_t850.append(compute_cri(sort_m1_samz_t850[i], sort_m2_samz_t850[i], sort_m3_samz_t850[i], sort_m4_samz_t850[i]))
	cri_samz_t500.append(compute_cri(sort_m1_samz_t500[i], sort_m2_samz_t500[i], sort_m3_samz_t500[i], sort_m4_namz_t500[i]))
	cri_samz_t200.append(compute_cri(sort_m1_samz_t200[i], sort_m2_samz_t200[i], sort_m3_samz_t200[i], sort_m4_samz_t200[i]))
	cri_samz_u850.append(compute_cri(sort_m1_samz_u850[i], sort_m2_samz_u850[i], sort_m3_samz_u850[i], sort_m4_samz_u850[i]))
	cri_samz_u500.append(compute_cri(sort_m1_samz_u500[i], sort_m2_samz_u500[i], sort_m3_samz_u500[i], sort_m4_samz_u500[i]))
	cri_samz_u200.append(compute_cri(sort_m1_samz_u200[i], sort_m2_samz_u200[i], sort_m3_samz_u200[i], sort_m4_samz_u200[i]))
	cri_samz_v850.append(compute_cri(sort_m1_samz_v850[i], sort_m2_samz_v850[i], sort_m3_samz_v850[i], sort_m4_samz_v850[i]))
	cri_samz_v500.append(compute_cri(sort_m1_samz_v500[i], sort_m2_samz_v500[i], sort_m3_samz_v500[i], sort_m4_samz_v500[i]))
	cri_samz_v200.append(compute_cri(sort_m1_samz_v200[i], sort_m2_samz_v200[i], sort_m3_samz_v200[i], sort_m4_samz_v200[i]))
	cri_samz_q850.append(compute_cri(sort_m1_samz_q850[i], sort_m2_samz_q850[i], sort_m3_samz_q850[i], sort_m4_samz_q850[i]))
	cri_samz_q500.append(compute_cri(sort_m1_samz_q500[i], sort_m2_samz_q500[i], sort_m3_samz_q500[i], sort_m4_samz_q500[i]))
	cri_samz_q200.append(compute_cri(sort_m1_samz_q200[i], sort_m2_samz_q200[i], sort_m3_samz_q200[i], sort_m4_samz_q200[i]))
	
	cri_sam_pr.append(compute_cri(sort_m1_sam_pr[i], sort_m2_sam_pr[i], sort_m3_sam_pr[i], sort_m4_sam_pr[i]))
	cri_sam_ps.append(compute_cri(sort_m1_sam_ps[i], sort_m2_sam_ps[i], sort_m3_sam_ps[i], sort_m4_sam_ps[i]))
	cri_sam_t850.append(compute_cri(sort_m1_sam_t850[i], sort_m2_sam_t850[i], sort_m3_sam_t850[i], sort_m4_sam_t850[i]))
	cri_sam_t500.append(compute_cri(sort_m1_sam_t500[i], sort_m2_sam_t500[i], sort_m3_sam_t500[i], sort_m4_sam_t500[i]))
	cri_sam_t200.append(compute_cri(sort_m1_sam_t200[i], sort_m2_sam_t200[i], sort_m3_sam_t200[i], sort_m4_sam_t200[i]))
	cri_sam_u850.append(compute_cri(sort_m1_sam_u850[i], sort_m2_sam_u850[i], sort_m3_sam_u850[i], sort_m4_sam_u850[i]))
	cri_sam_u500.append(compute_cri(sort_m1_sam_u500[i], sort_m2_sam_u500[i], sort_m3_sam_u500[i], sort_m4_sam_u500[i]))
	cri_sam_u200.append(compute_cri(sort_m1_sam_u200[i], sort_m2_sam_u200[i], sort_m3_sam_u200[i], sort_m4_sam_u200[i]))
	cri_sam_v850.append(compute_cri(sort_m1_sam_v850[i], sort_m2_sam_v850[i], sort_m3_sam_v850[i], sort_m4_sam_v850[i]))
	cri_sam_v500.append(compute_cri(sort_m1_sam_v500[i], sort_m2_sam_v500[i], sort_m3_sam_v500[i], sort_m4_sam_v500[i]))
	cri_sam_v200.append(compute_cri(sort_m1_sam_v200[i], sort_m2_sam_v200[i], sort_m3_sam_v200[i], sort_m4_sam_v200[i]))
	cri_sam_q850.append(compute_cri(sort_m1_sam_q850[i], sort_m2_sam_q850[i], sort_m3_sam_q850[i], sort_m4_sam_q850[i]))
	cri_sam_q500.append(compute_cri(sort_m1_sam_q500[i], sort_m2_sam_q500[i], sort_m3_sam_q500[i], sort_m4_sam_q500[i]))
	cri_sam_q200.append(compute_cri(sort_m1_sam_q200[i], sort_m2_sam_q200[i], sort_m3_sam_q200[i], sort_m4_sam_q200[i]))
	
	cri_neb_pr.append(compute_cri(sort_m1_neb_pr[i], sort_m2_neb_pr[i], sort_m3_neb_pr[i], sort_m4_neb_pr[i]))
	cri_neb_ps.append(compute_cri(sort_m1_neb_ps[i], sort_m2_neb_ps[i], sort_m3_neb_ps[i], sort_m4_neb_ps[i]))
	cri_neb_t850.append(compute_cri(sort_m1_neb_t850[i], sort_m2_neb_t850[i], sort_m3_neb_t850[i], sort_m4_neb_t850[i]))
	cri_neb_t500.append(compute_cri(sort_m1_neb_t500[i], sort_m2_neb_t500[i], sort_m3_neb_t500[i], sort_m4_neb_t500[i]))
	cri_neb_t200.append(compute_cri(sort_m1_neb_t200[i], sort_m2_neb_t200[i], sort_m3_neb_t200[i], sort_m4_neb_t200[i]))
	cri_neb_u850.append(compute_cri(sort_m1_neb_u850[i], sort_m2_neb_u850[i], sort_m3_neb_u850[i], sort_m4_neb_u850[i]))
	cri_neb_u500.append(compute_cri(sort_m1_neb_u500[i], sort_m2_neb_u500[i], sort_m3_neb_u500[i], sort_m4_neb_u500[i]))
	cri_neb_u200.append(compute_cri(sort_m1_neb_u200[i], sort_m2_neb_u200[i], sort_m3_neb_u200[i], sort_m4_neb_u200[i]))
	cri_neb_v850.append(compute_cri(sort_m1_neb_v850[i], sort_m2_neb_v850[i], sort_m3_neb_v850[i], sort_m4_neb_v850[i]))
	cri_neb_v500.append(compute_cri(sort_m1_neb_v500[i], sort_m2_neb_v500[i], sort_m3_neb_v500[i], sort_m4_neb_v500[i]))
	cri_neb_v200.append(compute_cri(sort_m1_neb_v200[i], sort_m2_neb_v200[i], sort_m3_neb_v200[i], sort_m4_neb_v200[i]))
	cri_neb_q850.append(compute_cri(sort_m1_neb_q850[i], sort_m2_neb_q850[i], sort_m3_neb_q850[i], sort_m4_neb_q850[i]))
	cri_neb_q500.append(compute_cri(sort_m1_neb_q500[i], sort_m2_neb_q500[i], sort_m3_neb_q500[i], sort_m4_neb_q500[i]))
	cri_neb_q200.append(compute_cri(sort_m1_neb_q200[i], sort_m2_neb_q200[i], sort_m3_neb_q200[i], sort_m4_neb_q200[i]))

	cri_lpb_pr.append(compute_cri(sort_m1_lpb_pr[i], sort_m2_lpb_pr[i], sort_m3_lpb_pr[i], sort_m4_lpb_pr[i]))
	cri_lpb_ps.append(compute_cri(sort_m1_lpb_ps[i], sort_m2_lpb_ps[i], sort_m3_lpb_ps[i], sort_m4_lpb_ps[i]))
	cri_lpb_t850.append(compute_cri(sort_m1_lpb_t850[i], sort_m2_lpb_t850[i], sort_m3_lpb_t850[i], sort_m4_lpb_t850[i]))
	cri_lpb_t500.append(compute_cri(sort_m1_lpb_t500[i], sort_m2_lpb_t500[i], sort_m3_lpb_t500[i], sort_m4_lpb_t500[i]))
	cri_lpb_t200.append(compute_cri(sort_m1_lpb_t200[i], sort_m2_lpb_t200[i], sort_m3_lpb_t200[i], sort_m4_lpb_t200[i]))
	cri_lpb_u850.append(compute_cri(sort_m1_lpb_u850[i], sort_m2_lpb_u850[i], sort_m3_lpb_u850[i], sort_m4_lpb_u850[i]))
	cri_lpb_u500.append(compute_cri(sort_m1_lpb_u500[i], sort_m2_lpb_u500[i], sort_m3_lpb_u500[i], sort_m4_lpb_u500[i]))
	cri_lpb_u200.append(compute_cri(sort_m1_lpb_u200[i], sort_m2_lpb_u200[i], sort_m3_lpb_u200[i], sort_m4_lpb_u200[i]))
	cri_lpb_v850.append(compute_cri(sort_m1_lpb_v850[i], sort_m2_lpb_v850[i], sort_m3_lpb_v850[i], sort_m4_lpb_v850[i]))
	cri_lpb_v500.append(compute_cri(sort_m1_lpb_v500[i], sort_m2_lpb_v500[i], sort_m3_lpb_v500[i], sort_m4_lpb_v500[i]))
	cri_lpb_v200.append(compute_cri(sort_m1_lpb_v200[i], sort_m2_lpb_v200[i], sort_m3_lpb_v200[i], sort_m4_lpb_v200[i]))
	cri_lpb_q850.append(compute_cri(sort_m1_lpb_q850[i], sort_m2_lpb_q850[i], sort_m3_lpb_q850[i], sort_m4_lpb_q850[i]))
	cri_lpb_q500.append(compute_cri(sort_m1_lpb_q500[i], sort_m2_lpb_q500[i], sort_m3_lpb_q500[i], sort_m4_lpb_q500[i]))
	cri_lpb_q200.append(compute_cri(sort_m1_lpb_q200[i], sort_m2_lpb_q200[i], sort_m3_lpb_q200[i], sort_m4_lpb_q200[i]))
	
sort_cri_namz_pr = sort_list_reverse(cri_namz_pr)
sort_cri_namz_ps = sort_list_reverse(cri_namz_ps)
sort_cri_namz_t850 = sort_list_reverse(cri_namz_t850)
sort_cri_namz_t500 = sort_list_reverse(cri_namz_t500)
sort_cri_namz_t200 = sort_list_reverse(cri_namz_t200)
sort_cri_namz_u850 = sort_list_reverse(cri_namz_u850)
sort_cri_namz_u500 = sort_list_reverse(cri_namz_u500)
sort_cri_namz_u200 = sort_list_reverse(cri_namz_u200)
sort_cri_namz_v850 = sort_list_reverse(cri_namz_v850)
sort_cri_namz_v500 = sort_list_reverse(cri_namz_v500)
sort_cri_namz_v200 = sort_list_reverse(cri_namz_v200)
sort_cri_namz_q850 = sort_list_reverse(cri_namz_q850)
sort_cri_namz_q500 = sort_list_reverse(cri_namz_q500)
sort_cri_namz_q200 = sort_list_reverse(cri_namz_q200)

sort_cri_namz = np.array([sort_cri_namz_pr, sort_cri_namz_ps, sort_cri_namz_t850, sort_cri_namz_t500, sort_cri_namz_t200, sort_cri_namz_u850, sort_cri_namz_u500, sort_cri_namz_u200, sort_cri_namz_v850, sort_cri_namz_v500, sort_cri_namz_v200, sort_cri_namz_q850, sort_cri_namz_q500, sort_cri_namz_q200])

sort_cri_samz_pr = sort_list_reverse(cri_samz_pr)
sort_cri_samz_ps = sort_list_reverse(cri_samz_ps)
sort_cri_samz_t850 = sort_list_reverse(cri_samz_t850)
sort_cri_samz_t500 = sort_list_reverse(cri_samz_t500)
sort_cri_samz_t200 = sort_list_reverse(cri_samz_t200)
sort_cri_samz_u850 = sort_list_reverse(cri_samz_u850)
sort_cri_samz_u500 = sort_list_reverse(cri_samz_u500)
sort_cri_samz_u200 = sort_list_reverse(cri_samz_u200)
sort_cri_samz_v850 = sort_list_reverse(cri_samz_v850)
sort_cri_samz_v500 = sort_list_reverse(cri_samz_v500)
sort_cri_samz_v200 = sort_list_reverse(cri_samz_v200)
sort_cri_samz_q850 = sort_list_reverse(cri_samz_q850)
sort_cri_samz_q500 = sort_list_reverse(cri_samz_q500)
sort_cri_samz_q200 = sort_list_reverse(cri_samz_q200)

sort_cri_samz = np.array([sort_cri_samz_pr, sort_cri_samz_ps, sort_cri_samz_t850, sort_cri_samz_t500, sort_cri_samz_t200, sort_cri_samz_u850, sort_cri_samz_u500, sort_cri_samz_u200, sort_cri_samz_v850, sort_cri_samz_v500, sort_cri_samz_v200, sort_cri_samz_q850, sort_cri_samz_q500, sort_cri_samz_q200])

sort_cri_sam_pr = sort_list_reverse(cri_sam_pr)
sort_cri_sam_ps = sort_list_reverse(cri_sam_ps)
sort_cri_sam_t850 = sort_list_reverse(cri_sam_t850)
sort_cri_sam_t500 = sort_list_reverse(cri_sam_t500)
sort_cri_sam_t200 = sort_list_reverse(cri_sam_t200)
sort_cri_sam_u850 = sort_list_reverse(cri_sam_u850)
sort_cri_sam_u500 = sort_list_reverse(cri_sam_u500)
sort_cri_sam_u200 = sort_list_reverse(cri_sam_u200)
sort_cri_sam_v850 = sort_list_reverse(cri_sam_v850)
sort_cri_sam_v500 = sort_list_reverse(cri_sam_v500)
sort_cri_sam_v200 = sort_list_reverse(cri_sam_v200)
sort_cri_sam_q850 = sort_list_reverse(cri_sam_q850)
sort_cri_sam_q500 = sort_list_reverse(cri_sam_q500)
sort_cri_sam_q200 = sort_list_reverse(cri_sam_q200)

sort_cri_sam = np.array([sort_cri_sam_pr, sort_cri_sam_ps, sort_cri_sam_t850, sort_cri_sam_t500, sort_cri_sam_t200, sort_cri_sam_u850, sort_cri_sam_u500, sort_cri_sam_u200, sort_cri_sam_v850, sort_cri_sam_v500, sort_cri_sam_v200, sort_cri_sam_q850, sort_cri_sam_q500, sort_cri_sam_q200])

sort_cri_neb_pr = sort_list_reverse(cri_neb_pr)
sort_cri_neb_ps = sort_list_reverse(cri_neb_ps)
sort_cri_neb_t850 = sort_list_reverse(cri_neb_t850)
sort_cri_neb_t500 = sort_list_reverse(cri_neb_t500)
sort_cri_neb_t200 = sort_list_reverse(cri_neb_t200)
sort_cri_neb_u850 = sort_list_reverse(cri_neb_u850)
sort_cri_neb_u500 = sort_list_reverse(cri_neb_u500)
sort_cri_neb_u200 = sort_list_reverse(cri_neb_u200)
sort_cri_neb_v850 = sort_list_reverse(cri_neb_v850)
sort_cri_neb_v500 = sort_list_reverse(cri_neb_v500)
sort_cri_neb_v200 = sort_list_reverse(cri_neb_v200)
sort_cri_neb_q850 = sort_list_reverse(cri_neb_q850)
sort_cri_neb_q500 = sort_list_reverse(cri_neb_q500)
sort_cri_neb_q200 = sort_list_reverse(cri_neb_q200)

sort_cri_neb = np.array([sort_cri_neb_pr, sort_cri_neb_ps, sort_cri_neb_t850, sort_cri_neb_t500, sort_cri_neb_t200, sort_cri_neb_u850, sort_cri_neb_u500, sort_cri_neb_u200, sort_cri_neb_v850, sort_cri_neb_v500, sort_cri_neb_v200, sort_cri_neb_q850, sort_cri_neb_q500, sort_cri_neb_q200])

sort_cri_lpb_pr = sort_list_reverse(cri_lpb_pr)
sort_cri_lpb_ps = sort_list_reverse(cri_lpb_ps)
sort_cri_lpb_t850 = sort_list_reverse(cri_lpb_t850)
sort_cri_lpb_t500 = sort_list_reverse(cri_lpb_t500)
sort_cri_lpb_t200 = sort_list_reverse(cri_lpb_t200)
sort_cri_lpb_u850 = sort_list_reverse(cri_lpb_u850)
sort_cri_lpb_u500 = sort_list_reverse(cri_lpb_u500)
sort_cri_lpb_u200 = sort_list_reverse(cri_lpb_u200)
sort_cri_lpb_v850 = sort_list_reverse(cri_lpb_v850)
sort_cri_lpb_v500 = sort_list_reverse(cri_lpb_v500)
sort_cri_lpb_v200 = sort_list_reverse(cri_lpb_v200)
sort_cri_lpb_q850 = sort_list_reverse(cri_lpb_q850)
sort_cri_lpb_q500 = sort_list_reverse(cri_lpb_q500)
sort_cri_lpb_q200 = sort_list_reverse(cri_lpb_q200)

sort_cri_lpb = np.array([sort_cri_lpb_pr, sort_cri_lpb_ps, sort_cri_lpb_t850, sort_cri_lpb_t500, sort_cri_lpb_t200, sort_cri_lpb_u850, sort_cri_lpb_u500, sort_cri_lpb_u200, sort_cri_lpb_v850, sort_cri_lpb_v500, sort_cri_lpb_v200, sort_cri_lpb_q850, sort_cri_lpb_q500, sort_cri_lpb_q200])

# Sort CRI rank all variables
cri_namz_vars = np.nanmean(sort_cri_namz, axis=0)
cri_namz_vars_rank = sort_list(cri_namz_vars) 

cri_samz_vars = np.nanmean(sort_cri_samz, axis=0)
cri_samz_vars_rank = sort_list(cri_samz_vars) 

cri_neb_vars = np.nanmean(sort_cri_neb, axis=0)
cri_neb_vars_rank = sort_list(cri_neb_vars) 

cri_sam_vars = np.nanmean(sort_cri_sam, axis=0)
cri_sam_vars_rank = sort_list(cri_sam_vars) 

cri_lpb_vars = np.nanmean(sort_cri_lpb, axis=0)
cri_lpb_vars_rank = sort_list(cri_lpb_vars) 

sort_cri_br = np.array([cri_namz_vars_rank, cri_samz_vars_rank, cri_sam_vars_rank, cri_neb_vars_rank, cri_lpb_vars_rank])
sort_cri_br_vars = np.nanmean(sort_cri_br, axis=0)
sort_cri_br_rank = sort_list(sort_cri_br_vars) 

# Plot cmip models and obs database 
fig = plt.figure(figsize=(9, 7))
cm = plt.get_cmap('jet')
colors = cm(np.linspace(0.1, 0.9, 17))

ax = fig.add_subplot(231, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(cri_namz_vars_rank)) 
bars = ax.bar(theta, cri_namz_vars_rank, color=colors, width=0.4)
ax.set_title(u'(a) NAMZ', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(232, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(cri_samz_vars_rank)) 
bars = bars = ax.bar(theta, cri_samz_vars_rank, color=colors, width=0.4)
ax.set_title(u'(b) SAMZ', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(233, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(cri_sam_vars_rank)) 
bars = ax.bar(theta, cri_sam_vars_rank, color=colors, width=0.4)
ax.set_title(u'(c) NEB', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(234, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(cri_neb_vars_rank)) 
bars = ax.bar(theta, cri_neb_vars_rank, color=colors, width=0.4)
ax.set_title(u'(d) SAM', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(235, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(cri_lpb_vars_rank)) 
bars = ax.bar(theta, cri_lpb_vars_rank, color=colors, width=0.4)
ax.set_title(u'(e) LPB', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

ax = fig.add_subplot(236, polar=True)
theta = np.arange(0, 2*np.pi, 2*np.pi/len(sort_cri_br_rank)) 
bars = ax.bar(theta, sort_cri_br_rank, color=colors, width=0.4)
ax.set_title(u'(f) BR', loc='left', fontsize=8, fontweight='bold')
ax.set_xticks(theta)
ax.set_xticklabels(range(1, len(theta)+1), fontsize=8)
ax.yaxis.grid(True)

legend = ['ACCESS-CM2', 'BCC-CSM2-MR', 'CanESM5', 'CMCC-ESM2', 'CNRM-CM6-1', 
'CNRM-ESM2-1', 'GFDL-ESM4', 'INM-CM4-8', 'INM-CM5-0', 'KIOST-ESM', 'MIROC6', 'MIROC-ES2L',
'MPI-ESM1-2-HR', 'MPI-ESM1-2-LR', 'MRI-ESM2-0', 'NESM3', 'NorESM2-MM']
plt.legend(bars, legend, ncol=6, loc=(-2.53, -0.43), fontsize=8)

# Path out to save figure
path_out = '{0}/figs'.format(path)
name_out = 'pyplt_radar_rank_cmip6_obs_{0}.png'.format(dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

