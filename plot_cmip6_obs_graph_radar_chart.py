# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot cri of cmip6 models"

import os
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from dict_cmip6_models_name import cmip6, cmip6_i
from comp_statistical_metrics import compute_mbe
from comp_statistical_metrics import compute_rmse
from comp_statistical_metrics import compute_tss
from comp_statistical_metrics import compute_pcc
from comp_statistical_metrics import compute_ivs


def import_obs_latlon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
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
	
	
def import_obs_mon(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mon_mean = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	ts_mon = []
	for i in range(0, 11 + 1):
		clim = np.nanmean(mon_mean[i::12], axis=0)
		ts_mon.append(clim)
	
	return ts_mon


def import_obs_ann(param, area, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_{3}_{4}_lonlat.nc'.format(path, param, area, period, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	return value
	
	
def import_cmip_latlon(param, area, model, exp, period, date):

	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
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
	              

def import_cmip_mon(param, area, model, exp, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mon_mean = np.nanmean(np.nanmean(value, axis=1), axis=1)
	
	ts_mon = []
	for i in range(0, 11 + 1):
		clim = np.nanmean(mon_mean[i::12], axis=0)
		ts_mon.append(clim)
	
	return ts_mon
	

def import_cmip_ann(param, area, model, exp, period, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_{5}_{6}_lonlat.nc'.format(path, param, area, model, exp, period, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	return value
	

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
		model_list.append(cmip6_i[ii+1][0])
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
		model_list.append(cmip6_i[ii+1][0])
		value_list.append(data_list[ii])

	data_argsort = np.argsort(model_list)
	
	data_argsort_i = []
	for idx in data_argsort:
		data_argsort_i.append(idx+1)
	
	return data_argsort_i
		

def compute_cri(rank1,rank2,rank3,rank4):
	
	p1 = (rank1+rank2+rank3+rank4)
	p2 = p1/68
	cri = 1 - p2
	
	return cri
	
	
# Import cmip models and obs database 
area_cmip6 = 'namz'
var_obs = 'pr'
var_cmip6 = 'pr'
dt = '1986-2005'

namz_obs_latlon = import_obs_latlon(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_latlon = import_obs_latlon(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_latlon = import_obs_latlon(var_obs, 'NEB', 'ANN', dt)
sam_obs_latlon = import_obs_latlon(var_obs, 'SAM', 'ANN', dt)
lpb_obs_latlon = import_obs_latlon(var_obs, 'LPB', 'ANN', dt)
br_obs_latlon = import_obs_latlon(var_obs, 'BR', 'ANN', dt)

namz_obs_mon_ts = import_obs_mon(var_obs, 'NAMZ', 'MON', dt)
samz_obs_mon_ts = import_obs_mon(var_obs, 'SAMZ', 'MON', dt)
neb_obs_mon_ts  = import_obs_mon(var_obs, 'NEB', 'MON', dt)
sam_obs_mon_ts = import_obs_mon(var_obs, 'SAM', 'MON', dt)
lpb_obs_mon_ts = import_obs_mon(var_obs, 'LPB', 'MON', dt)
br_obs_mon_ts = import_obs_mon(var_obs, 'BR', 'MON', dt)

namz_obs_ann_ts = import_obs_ann(var_obs, 'NAMZ', 'ANN', dt)
samz_obs_ann_ts = import_obs_ann(var_obs, 'SAMZ', 'ANN', dt)
neb_obs_ann_ts  = import_obs_ann(var_obs, 'NEB', 'ANN', dt)
sam_obs_ann_ts = import_obs_ann(var_obs, 'SAM', 'ANN', dt)
lpb_obs_ann_ts = import_obs_ann(var_obs, 'LPB', 'ANN', dt)
br_obs_ann_ts = import_obs_ann(var_obs, 'BR', 'ANN', dt)

mbe_namz_cmip6 = []
rmse_namz_cmip6 = []
tss_namz_cmip6 = []
pcc_namz_cmip6 = []
ivs_namz_cmip6 = []

mbe_samz_cmip6 = []
rmse_samz_cmip6 = []
tss_samz_cmip6 = []
pcc_samz_cmip6 = []
ivs_samz_cmip6 = []

mbe_neb_cmip6 = []
rmse_neb_cmip6 = []
tss_neb_cmip6 = []
pcc_neb_cmip6 = []
ivs_neb_cmip6 = []

mbe_sam_cmip6 = []
rmse_sam_cmip6 = []
tss_sam_cmip6 = []
pcc_sam_cmip6 = []
ivs_sam_cmip6 = []

mbe_lpb_cmip6 = []
rmse_lpb_cmip6 = []
tss_lpb_cmip6 = []
pcc_lpb_cmip6 = []
ivs_lpb_cmip6 = []

mbe_br_cmip6 = []
rmse_br_cmip6 = []
tss_br_cmip6 = []
pcc_br_cmip6 = []
ivs_br_cmip6 = []

legend = []

for i in range(1, 18):

	namz_cmip_latlon = import_cmip_latlon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_latlon = import_cmip_latlon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_latlon = import_cmip_latlon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_latlon = import_cmip_latlon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	br_cmip_latlon = import_cmip_latlon(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	
	namz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	samz_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'MON', dt)
	neb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	sam_cmip_mon_ts = import_cmip_mon(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'MON', dt)
	lpb_cmip_mon_ts = import_cmip_mon(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'MON', dt)
	br_cmip_mon_ts = import_cmip_mon(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'MON', dt)
	
	namz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	samz_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	neb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	sam_cmip_ann_ts = import_cmip_ann(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	lpb_cmip_ann_ts = import_cmip_ann(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], 'ANN', dt)
	br_cmip_ann_ts = import_cmip_ann(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], 'ANN', dt)
			
	mbe_namz_cmip6.append(compute_mbe(namz_cmip_latlon, namz_obs_latlon))
	rmse_namz_cmip6.append(compute_rmse(namz_cmip_latlon, namz_obs_latlon))	
	tss_namz_cmip6.append(compute_tss(namz_obs_latlon, namz_cmip_latlon))
	pcc_namz_cmip6.append(compute_pcc(namz_obs_mon_ts, namz_cmip_mon_ts))
	ivs_namz_cmip = compute_ivs(namz_obs_ann_ts, namz_cmip_ann_ts)
	ivs_namz_cmip6.append(np.nanmean(ivs_namz_cmip))
	
	mbe_samz_cmip6.append(compute_mbe(samz_cmip_latlon, samz_obs_latlon))
	rmse_samz_cmip6.append(compute_rmse(samz_cmip_latlon, samz_obs_latlon))
	tss_samz_cmip6.append(compute_tss(samz_obs_latlon, samz_cmip_latlon))
	pcc_samz_cmip6.append(compute_pcc(samz_obs_mon_ts, samz_cmip_mon_ts))
	ivs_samz_cmip = compute_ivs(samz_obs_ann_ts, samz_cmip_ann_ts)
	ivs_samz_cmip6.append(np.nanmean(ivs_samz_cmip))
		
	mbe_neb_cmip6.append(compute_mbe(neb_cmip_latlon, neb_obs_latlon))
	rmse_neb_cmip6.append(compute_rmse(neb_cmip_latlon, neb_obs_latlon))
	tss_neb_cmip6.append(compute_tss(neb_obs_latlon, neb_cmip_latlon))
	pcc_neb_cmip6.append(compute_pcc(neb_obs_mon_ts, neb_cmip_mon_ts))
	ivs_neb_cmip = compute_ivs(neb_obs_ann_ts, neb_cmip_ann_ts)
	ivs_neb_cmip6.append(np.nanmean(ivs_neb_cmip))
		
	mbe_sam_cmip6.append(compute_mbe(sam_cmip_latlon, sam_obs_latlon))
	rmse_sam_cmip6.append(compute_rmse(sam_cmip_latlon, sam_obs_latlon))
	tss_sam_cmip6.append(compute_tss(sam_obs_latlon, sam_cmip_latlon))
	pcc_sam_cmip6.append(compute_pcc(sam_obs_mon_ts, sam_cmip_mon_ts))
	ivs_sam_cmip = compute_ivs(sam_obs_ann_ts, sam_cmip_ann_ts)
	ivs_sam_cmip6.append(np.nanmean(ivs_sam_cmip))
		
	mbe_lpb_cmip6.append(compute_mbe(lpb_cmip_latlon, lpb_obs_latlon))
	rmse_lpb_cmip6.append(compute_rmse(lpb_cmip_latlon, lpb_obs_latlon))
	tss_lpb_cmip6.append(compute_tss(lpb_obs_latlon, lpb_cmip_latlon))
	pcc_lpb_cmip6.append(compute_pcc(lpb_obs_mon_ts, lpb_cmip_mon_ts))
	ivs_lpb_cmip = compute_ivs(lpb_obs_ann_ts, lpb_cmip_ann_ts)
	ivs_lpb_cmip6.append(np.nanmean(ivs_lpb_cmip))
	
	mbe_br_cmip6.append(compute_mbe(br_cmip_latlon, br_obs_latlon))
	rmse_br_cmip6.append(compute_rmse(br_cmip_latlon, br_obs_latlon))
	tss_br_cmip6.append(compute_tss(br_obs_latlon, br_cmip_latlon))
	pcc_br_cmip6.append(compute_pcc(br_obs_mon_ts, br_cmip_mon_ts))
	ivs_br_cmip = compute_ivs(br_obs_ann_ts, br_cmip_ann_ts)
	ivs_br_cmip6.append(np.nanmean(ivs_br_cmip))

	legend.append(cmip6[i][0])

if area_cmip6 == 'namz':
	argsort_rmse = sort_list(rmse_namz_cmip6)
	argsort_tss = sort_list_reverse(tss_namz_cmip6)
	argsort_pcc = sort_list_reverse(pcc_namz_cmip6)
	argsort_ivs = sort_list(ivs_namz_cmip6)
elif area_cmip6 == 'samz':
	argsort_rmse = sort_list(rmse_samz_cmip6)
	argsort_tss = sort_list_reverse(tss_samz_cmip6)
	argsort_pcc = sort_list_reverse(pcc_samz_cmip6)
	argsort_ivs = sort_list(ivs_samz_cmip6)
elif area_cmip6 == 'neb':
	argsort_rmse = sort_list(rmse_neb_cmip6)
	argsort_tss = sort_list_reverse(tss_neb_cmip6)
	argsort_pcc = sort_list_reverse(pcc_neb_cmip6)
	argsort_ivs = sort_list(ivs_neb_cmip6)
elif area_cmip6 == 'sam':
	argsort_rmse = sort_list(rmse_sam_cmip6)
	argsort_tss = sort_list_reverse(tss_sam_cmip6)
	argsort_pcc = sort_list_reverse(pcc_sam_cmip6)
	argsort_ivs = sort_list(ivs_sam_cmip6)
elif area_cmip6 == 'lpb':
	argsort_rmse = sort_list(rmse_lpb_cmip6)
	argsort_tss = sort_list_reverse(tss_lpb_cmip6)
	argsort_pcc = sort_list_reverse(pcc_lpb_cmip6)
	argsort_ivs = sort_list(ivs_lpb_cmip6)
else:
	argsort_rmse = sort_list(rmse_br_cmip6)
	argsort_tss = sort_list_reverse(tss_br_cmip6)
	argsort_pcc = sort_list_reverse(pcc_br_cmip6)
	argsort_ivs = sort_list(ivs_br_cmip6)

# ~ print()
# ~ print('rmse_{0}_{1} ='.format(area_cmip6, var_cmip6), argsort_rmse)
# ~ print('tss_{0}_{1} ='.format(area_cmip6, var_cmip6), argsort_tss)
# ~ print('pcc_{0}_{1} ='.format(area_cmip6, var_cmip6), argsort_pcc)
# ~ print('ivs_{0}_{1} ='.format(area_cmip6, var_cmip6), argsort_ivs)
# ~ print()

rmse_namz_pr = [9, 16, 15, 7, 11, 10, 8, 4, 3, 17, 6, 1, 13, 14, 2, 12, 5]
tss_namz_pr = [4, 11, 1, 6, 17, 16, 10, 14, 15, 13, 7, 2, 12, 8, 5, 3, 9]
pcc_namz_pr = [16, 14, 12, 11, 8, 9, 5, 10, 4, 15, 6, 3, 7, 13, 2, 17, 1]
ivs_namz_pr = [4, 12, 8, 9, 11, 10, 3, 15, 13, 17, 2, 14, 16, 7, 1, 6, 5]

rmse_samz_pr = [9, 10, 16, 15, 11, 14, 6, 2, 8, 17, 3, 4, 1, 5, 12, 13, 7]
tss_samz_pr = [15, 10, 13, 17, 12, 14, 9, 5, 6, 7, 2, 11, 3, 4, 1, 16, 8]
pcc_samz_pr = [4, 8, 14, 7, 10, 11, 12, 13, 9, 16, 6, 15, 1, 5, 2, 17, 3]
ivs_samz_pr = [2, 12, 3, 15, 10, 11, 5, 16, 14, 17, 4, 13, 6, 8, 7, 1, 9]

rmse_neb_pr = [12, 15, 7, 16, 10, 8, 11, 13, 14, 1, 6, 9, 3, 5, 2, 17, 4]
tss_neb_pr = [9, 6, 16, 7, 14, 10, 4, 13, 11, 5, 15, 12, 3, 8, 1, 17, 2]
pcc_neb_pr = [7, 5, 13, 3, 6, 4, 8, 1, 2, 15, 12, 16, 9, 11, 14, 17, 10]
ivs_neb_pr = [16, 12, 6, 5, 8, 9, 14, 13, 10, 17, 11, 3, 4, 1, 7, 2, 15]

rmse_sam_pr = [16, 11, 12, 6, 9, 10, 7, 3, 1, 14, 13, 15, 4, 8, 5, 17, 2]
tss_sam_pr = [10, 7, 16, 11, 13, 14, 8, 6, 3, 9, 15, 17, 2, 5, 1, 12, 4]
pcc_sam_pr = [7, 1, 14, 9, 15, 17, 13, 16, 12, 11, 2, 4, 8, 6, 10, 5, 3]
ivs_sam_pr = [10, 15, 3, 12, 7, 1, 11, 16, 5, 17, 14, 13, 9, 6, 8, 2, 4]

rmse_lpb_pr = [17, 11, 2, 3, 14, 15, 8, 12, 10, 16, 6, 13, 1, 4, 5, 7, 9]
tss_lpb_pr = [12, 10, 4, 2, 6, 7, 3, 15, 14, 5, 9, 17, 1, 11, 8, 16, 13]
pcc_lpb_pr = [16, 7, 9, 8, 11, 5, 1, 2, 4, 13, 6, 12, 17, 14, 15, 10, 3]
ivs_lpb_pr = [15, 11, 5, 10, 2, 7, 3, 16, 13, 17, 4, 12, 1, 6, 14, 9, 8]

rmse_br_pr = [13, 16, 14, 11, 10, 12, 6, 5, 4, 17, 3, 8, 7, 9, 1, 15, 2]
tss_br_pr = [14, 15, 16, 13, 11, 10, 8, 4, 3, 12, 5, 7, 6, 9, 1, 17, 2]
pcc_br_pr = [7, 10, 9, 3, 4, 5, 8, 2, 1, 16, 13, 15, 11, 14, 12, 17, 6]
ivs_br_pr = [7, 11, 6, 12, 10, 9, 3, 16, 14, 17, 1, 13, 15, 8, 2, 5, 4]

rmse_namz_tasmax = [2, 11, 12, 5, 7, 3, 4, 16, 15, 17, 10, 13, 9, 8, 6, 14, 1]
tss_namz_tasmax = [7, 10, 12, 13, 1, 2, 11, 6, 5, 17, 15, 3, 16, 14, 4, 9, 8]
pcc_namz_tasmax = [16, 11, 15, 8, 7, 6, 2, 9, 4, 14, 3, 10, 5, 13, 12, 17, 1]
ivs_namz_tasmax = [4, 1, 5, 3, 8, 7, 6, 17, 15, 14, 16, 11, 12, 10, 13, 2, 9]

rmse_samz_tasmax = [1, 8, 5, 4, 7, 3, 6, 13, 14, 16, 9, 11, 12, 15, 10, 17, 2]
tss_samz_tasmax = [11, 13, 15, 16, 2, 1, 7, 5, 6, 17, 9, 10, 4, 8, 3, 12, 14]
pcc_samz_tasmax = [16, 5, 14, 2, 1, 3, 10, 7, 8, 4, 9, 6, 12, 15, 13, 17, 11]
ivs_samz_tasmax = [9, 3, 13, 4, 2, 5, 1, 14, 8, 15, 17, 16, 10, 6, 7, 11, 12]

rmse_neb_tasmax = [8, 2, 3, 5, 12, 6, 11, 15, 13, 14, 9, 4, 10, 16, 1, 17, 7]
tss_neb_tasmax = [8, 7, 10, 9, 1, 3, 4, 15, 16, 17, 14, 12, 6, 11, 2, 13, 5]
pcc_neb_tasmax = [9, 14, 10, 7, 13, 11, 8, 4, 6, 16, 17, 15, 12, 1, 5, 2, 3]
ivs_neb_tasmax = [1, 11, 7, 8, 10, 9, 5, 17, 14, 15, 16, 12, 2, 3, 4, 13, 6]

rmse_sam_tasmax = [8, 1, 2, 4, 7, 3, 12, 11, 10, 15, 5, 13, 14, 16, 6, 17, 9]
tss_sam_tasmax = [8, 11, 12, 10, 5, 7, 2, 14, 13, 15, 16, 17, 1, 3, 9, 4, 6]
pcc_sam_tasmax = [3, 13, 10, 5, 8, 9, 4, 7, 14, 11, 2, 1, 6, 15, 12, 17, 16]
ivs_sam_tasmax = [1, 3, 5, 4, 8, 12, 14, 7, 6, 17, 16, 15, 2, 13, 11, 9, 10]

rmse_lpb_tasmax = [11, 9, 4, 2, 1, 3, 10, 7, 8, 14, 12, 13, 15, 16, 5, 17, 6]
tss_lpb_tasmax = [8, 9, 4, 6, 2, 5, 3, 13, 12, 10, 14, 17, 7, 15, 1, 16, 11]
pcc_lpb_tasmax = [11, 12, 13, 10, 7, 1, 5, 16, 14, 17, 9, 2, 6, 4, 8, 15, 3]
ivs_lpb_tasmax = [6, 3, 7, 12, 10, 14, 13, 5, 8, 17, 11, 9, 16, 15, 1, 4, 2]

rmse_br_tasmax = [5, 7, 6, 1, 8, 2, 9, 14, 13, 16, 10, 11, 12, 15, 4, 17, 3]
tss_br_tasmax = [8, 15, 12, 11, 7, 10, 3, 4, 2, 14, 16, 17, 9, 6, 1, 13, 5]
pcc_br_tasmax = [15, 14, 7, 3, 5, 2, 1, 8, 9, 16, 6, 4, 10, 12, 13, 17, 11]
ivs_br_tasmax = [1, 2, 7, 3, 6, 8, 4, 15, 13, 17, 16, 14, 12, 10, 9, 5, 11]

rmse_namz_tasmin = [6, 7, 17, 9, 11, 4, 5, 3, 1, 10, 15, 13, 14, 12, 8, 16, 2]
tss_namz_tasmin = [6, 7, 13, 9, 2, 1, 10, 11, 8, 4, 16, 14, 17, 15, 3, 12, 5]
pcc_namz_tasmin = [9, 7, 6, 8, 3, 4, 1, 11, 12, 16, 15, 14, 10, 13, 5, 17, 2]
ivs_namz_tasmin = [1, 12, 8, 4, 13, 15, 6, 16, 17, 3, 14, 10, 5, 7, 2, 11, 9]

rmse_samz_tasmin = [8, 12, 17, 15, 6, 1, 2, 4, 5, 11, 13, 14, 9, 3, 10, 16, 7]
tss_samz_tasmin = [10, 3, 13, 17, 2, 6, 4, 12, 14, 9, 5, 16, 1, 7, 11, 8, 15]
pcc_samz_tasmin = [2, 13, 14, 9, 8, 12, 10, 4, 5, 17, 6, 16, 11, 7, 3, 15, 1]
ivs_samz_tasmin = [3, 9, 6, 1, 10, 13, 4, 16, 17, 5, 11, 14, 12, 15, 8, 2, 7]

rmse_neb_tasmin = [10, 8, 15, 17, 9, 1, 2, 4, 5, 6, 16, 12, 7, 3, 13, 14, 11]
tss_neb_tasmin = [7, 10, 9, 17, 1, 2, 3, 15, 16, 13, 6, 14, 4, 12, 11, 8, 5]
pcc_neb_tasmin = [7, 10, 4, 3, 2, 6, 1, 15, 13, 12, 17, 8, 14, 11, 5, 16, 9]
ivs_neb_tasmin = [3, 9, 11, 5, 15, 12, 7, 16, 17, 1, 13, 14, 4, 2, 6, 8, 10]

rmse_sam_tasmin = [7, 14, 15, 17, 8, 2, 1, 4, 3, 9, 13, 11, 6, 5, 12, 16, 10]
tss_sam_tasmin = [8, 15, 13, 10, 4, 3, 2, 16, 14, 11, 9, 17, 1, 6, 5, 7, 12]
pcc_sam_tasmin = [4, 11, 13, 14, 15, 16, 12, 2, 3, 17, 8, 5, 9, 1, 6, 10, 7]
ivs_sam_tasmin = [2, 1, 5, 12, 11, 7, 9, 6, 15, 10, 14, 17, 4, 8, 13, 3, 16]

rmse_lpb_tasmin = [4, 13, 15, 12, 5, 3, 2, 8, 9, 11, 14, 17, 1, 6, 7, 16, 10]
tss_lpb_tasmin = [10, 6, 16, 11, 2, 3, 5, 14, 15, 8, 9, 17, 1, 13, 4, 12, 7]
pcc_lpb_tasmin = [10, 5, 7, 3, 9, 4, 11, 2, 12, 17, 14, 13, 6, 8, 15, 16, 1]
ivs_lpb_tasmin = [9, 4, 7, 2, 10, 17, 6, 15, 8, 11, 13, 16, 14, 12, 3, 1, 5]

rmse_br_tasmin = [6, 11, 16, 14, 7, 1, 2, 3, 4, 12, 13, 15, 8, 5, 10, 17, 9]
tss_br_tasmin = [4, 9, 17, 16, 13, 15, 14, 6, 7, 10, 3, 2, 12, 11, 8, 1, 5]
pcc_br_tasmin = [9, 12, 14, 1, 8, 11, 3, 5, 13, 17, 16, 15, 4, 10, 6, 7, 2]
ivs_br_tasmin = [1, 10, 7, 2, 13, 15, 4, 16, 17, 3, 12, 14, 9, 11, 6, 5, 8]

cri_namz_pr = []
cri_samz_pr = []
cri_neb_pr = []
cri_sam_pr = []
cri_lpb_pr = []
cri_br_pr = []

cri_namz_tasmax = []
cri_samz_tasmax = []
cri_neb_tasmax = []
cri_sam_tasmax = []
cri_lpb_tasmax = []
cri_br_tasmax = []

cri_namz_tasmin = []
cri_samz_tasmin = []
cri_neb_tasmin = []
cri_sam_tasmin = []
cri_lpb_tasmin = []
cri_br_tasmin = []

for i in range(0, 17):

	cri_namz_pr.append(compute_cri(rmse_namz_pr[i],tss_namz_pr[i],pcc_namz_pr[i],ivs_namz_pr[i]))
	cri_samz_pr.append(compute_cri(rmse_samz_pr[i],tss_samz_pr[i],pcc_samz_pr[i],ivs_samz_pr[i]))
	cri_neb_pr.append(compute_cri(rmse_neb_pr[i],tss_neb_pr[i],pcc_neb_pr[i],ivs_neb_pr[i]))
	cri_sam_pr.append(compute_cri(rmse_sam_pr[i],tss_sam_pr[i],pcc_sam_pr[i],ivs_sam_pr[i]))
	cri_lpb_pr.append(compute_cri(rmse_lpb_pr[i],tss_lpb_pr[i],pcc_lpb_pr[i],ivs_lpb_pr[i]))
	cri_br_pr.append(compute_cri(rmse_br_pr[i],tss_br_pr[i],pcc_br_pr[i],ivs_br_pr[i]))

	cri_namz_tasmax.append(compute_cri(rmse_namz_tasmax[i],tss_namz_tasmax[i],pcc_namz_tasmax[i],ivs_namz_tasmax[i]))
	cri_samz_tasmax.append(compute_cri(rmse_samz_tasmax[i],tss_samz_tasmax[i],pcc_samz_tasmax[i],ivs_samz_tasmax[i]))
	cri_neb_tasmax.append(compute_cri(rmse_neb_tasmax[i],tss_neb_tasmax[i],pcc_neb_tasmax[i],ivs_neb_tasmax[i]))
	cri_sam_tasmax.append(compute_cri(rmse_sam_tasmax[i],tss_sam_tasmax[i],pcc_sam_tasmax[i],ivs_sam_tasmax[i]))
	cri_lpb_tasmax.append(compute_cri(rmse_lpb_tasmax[i],tss_lpb_tasmax[i],pcc_lpb_tasmax[i],ivs_lpb_tasmax[i]))
	cri_br_tasmax.append(compute_cri(rmse_br_tasmax[i],tss_br_tasmax[i],pcc_br_tasmax[i],ivs_br_tasmax[i]))

	cri_namz_tasmin.append(compute_cri(rmse_namz_tasmin[i],tss_namz_tasmin[i],pcc_namz_tasmin[i],ivs_namz_tasmin[i]))
	cri_samz_tasmin.append(compute_cri(rmse_samz_tasmin[i],tss_samz_tasmin[i],pcc_samz_tasmin[i],ivs_samz_tasmin[i]))
	cri_neb_tasmin.append(compute_cri(rmse_neb_tasmin[i],tss_neb_tasmin[i],pcc_neb_tasmin[i],ivs_neb_tasmin[i]))
	cri_sam_tasmin.append(compute_cri(rmse_sam_tasmin[i],tss_sam_tasmin[i],pcc_sam_tasmin[i],ivs_sam_tasmin[i]))
	cri_lpb_tasmin.append(compute_cri(rmse_lpb_tasmin[i],tss_lpb_tasmin[i],pcc_lpb_tasmin[i],ivs_lpb_tasmin[i]))
	cri_br_tasmin.append(compute_cri(rmse_br_tasmin[i],tss_br_tasmin[i],pcc_br_tasmin[i],ivs_br_tasmin[i]))

sort_cri_namz_pr = sort_list_reverse(cri_namz_pr)
sort_cri_samz_pr = sort_list_reverse(cri_samz_pr)
sort_cri_neb_pr = sort_list_reverse(cri_neb_pr)
sort_cri_sam_pr = sort_list_reverse(cri_sam_pr)
sort_cri_lpb_pr = sort_list_reverse(cri_lpb_pr)
sort_cri_br_pr = sort_list_reverse(cri_br_pr)

sort_cri_namz_tasmax = sort_list_reverse(cri_namz_tasmax)
sort_cri_samz_tasmax = sort_list_reverse(cri_samz_tasmax)
sort_cri_neb_tasmax = sort_list_reverse(cri_neb_tasmax)
sort_cri_sam_tasmax = sort_list_reverse(cri_sam_tasmax)
sort_cri_lpb_tasmax = sort_list_reverse(cri_lpb_tasmax)
sort_cri_br_tasmax = sort_list_reverse(cri_br_tasmax)

sort_cri_namz_tasmin = sort_list_reverse(cri_namz_tasmin)
sort_cri_samz_tasmin = sort_list_reverse(cri_samz_tasmin)
sort_cri_neb_tasmin = sort_list_reverse(cri_neb_tasmin)
sort_cri_sam_tasmin = sort_list_reverse(cri_sam_tasmin)
sort_cri_lpb_tasmin = sort_list_reverse(cri_lpb_tasmin)
sort_cri_br_tasmin = sort_list_reverse(cri_br_tasmin)

print(sort_cri_br_tasmin)
exit()














	


