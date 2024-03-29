# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script plot skill of cmip6 correct models"

import os
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

from pylab import setp
from netCDF4 import Dataset
from datetime import datetime
from scipy.stats import norm
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = "/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias"


def import_observed(var_name, target_date):

    """ Import observed data.
    :param: model_name: str (BR-DWGD)
    :param: target_date: datetime
    :return: observed data (Prec, Tasmax and Tasmin)
    :return: flag that it refers to data control
    :rtype: 3D array
    """
   
    arq = xr.open_dataset('{0}/obs/'.format(dataset_dir) + '{0}_{1}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc'.format(var_name, target_date))
    data = arq[var_name]
    var = data.sel(time=slice('{}-01-01'.format(year),'{}-12-31'.format(year)))
    day = var.values
    day_mean = np.nanmean(np.nanmean(day, axis=1), axis=1)

    return day_mean
   
   
def import_simulated(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	arq = xr.open_dataset('{0}/cmip6/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('{}-01-01'.format(year),'{}-12-31'.format(year)))
	day = var.values
	day_mean = np.nanmean(np.nanmean(day, axis=1), axis=1)

	return day_mean


def import_simulated_correct(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	arq = xr.open_dataset('{0}/cmip6_correct/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_corrected.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('{}-01-01'.format(year),'{}-12-31'.format(year)))
	day = var.values
	day_mean = np.nanmean(np.nanmean(day, axis=1), axis=1)

	return day_mean


def setBoxColors(bp):
    setp(bp['boxes'][0], color='black')
    setp(bp['medians'][0], color='black')
    setp(bp['caps'][0], color='black')
    setp(bp['caps'][1], color='black')
    setp(bp['whiskers'][0], color='black')
    setp(bp['whiskers'][1], color='black')
        
    setp(bp['boxes'][1], color='red')
    setp(bp['medians'][1], color='red')
    setp(bp['caps'][2], color='red')
    setp(bp['caps'][3], color='red')
    setp(bp['whiskers'][2], color='red')
    setp(bp['whiskers'][3], color='red')
       
    setp(bp['boxes'][2], color='blue')
    setp(bp['medians'][2], color='blue')
    setp(bp['caps'][4], color='blue')
    setp(bp['caps'][5], color='blue')
    setp(bp['whiskers'][4], color='blue')
    setp(bp['whiskers'][5], color='blue')
    

def compute_cdf(data):

	"""
	The input arrays must have the same dimensions
	:Param data: Numpy array with model or obs data
	:Return: Cumulative Density Function
	"""

	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	cdf = norm.cdf(x,y,z)

	return x, cdf
	

# Best models list
best_models = [7, 9, 13, 15, 17]

# Variable dictionary
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}

experiment = 'historical'
dt = '19860101-20051231'
year=2005 # 1986 1995 2005

for i in best_models:
	for j in range(1, 4):	
		var_obs = var_dict[j][0]
		var_cmip6 = var_dict[j][1]
	
		print(cmip6[i][0])
		print(var_cmip6)
		print(year)

		# Import cmip models and obs database 
		day_obs = import_observed(var_obs, dt)
		day_sim = import_simulated(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
		day_sim_correct = import_simulated_correct(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)

		day_boxplot = [day_obs, day_sim, day_sim_correct]

		# Import cdf function
		x_obs, cdf_obs = compute_cdf(day_obs)
		x_sim, cdf_sim = compute_cdf(day_sim)
		x_correct_bias, cdf_correct_bias = compute_cdf(day_sim_correct)

		# Plot figure
		fig = plt.figure(figsize=(12, 8))

		if var_cmip6 == 'pr':
			legend = 'Precipitação (mm d⁻¹)'
		elif var_cmip6 == 'tasmax':
			legend = 'Temperatura máxima (°C)'
		else:
			legend = 'Temperatura mínima(°C)'
			
		ax = fig.add_subplot(1, 2, 1)
		x = np.arange(1, 3 + 1)
		bp = plt.boxplot(day_boxplot, positions=[1, 2, 3], sym='.')
		setBoxColors(bp)
		plt.title('(a) Boxplot diário ({0}) {1}'.format(year, cmip6[i][0]), loc='left', fontweight='bold')
		plt.xticks(x, ('Observação','Simulação','Correção'))
		plt.xlabel('Conjunto de dados', fontweight='bold')
		plt.ylabel('{0}'.format(legend), fontweight='bold')
		plt.grid(linestyle='--')
		if var_cmip6 == 'pr':
			plt.yticks(np.arange(0, 18, 2))
			plt.ylim(0, 16)
		elif var_cmip6 == 'tasmax':
			plt.yticks(np.arange(18, 38, 2))
			plt.ylim(18, 36)
		else:
			plt.yticks(np.arange(14, 30, 2))
			plt.ylim(14, 28)
			
		ax = fig.add_subplot(1, 2, 2)
		plt.plot(x_obs, cdf_obs, color='black', label='Observação', linewidth=1.5)
		plt.plot(x_sim, cdf_sim,  color='red', label='Simulação', linewidth=1.5)
		plt.plot(x_correct_bias, cdf_correct_bias,  color='blue', label='Correção', linewidth=1.5)  
		plt.title('(b) CDF diário ({0}) {1}'.format(year, cmip6[i][0]), loc='left', fontweight='bold')
		plt.xlabel('{0}'.format(legend), fontweight='bold')
		plt.ylabel('CDF', fontweight='bold')
		if var_cmip6 == 'pr':
			plt.xticks(np.arange(0, 18, 2))
			plt.xlim(0, 16)
		elif var_cmip6 == 'tasmax':
			plt.xticks(np.arange(18, 38, 2))
			plt.xlim(18, 36)
		else:
			plt.xticks(np.arange(14, 30, 2))
			plt.xlim(14, 28)
		plt.grid(linestyle='--')
		plt.legend(loc=9, ncol=3, frameon=False)

		# Path out to save figure
		path_out = '/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/figs/correct_bias'
		name_out = 'pyplt_boxplot_cdf_correct_bias_cmip6_{0}_{1}_{2}.png'.format(cmip6[i][0], var_cmip6, year)
		plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
		plt.close('all')
		plt.cla()
exit()
