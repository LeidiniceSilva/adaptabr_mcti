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

from netCDF4 import Dataset
from datetime import datetime
from scipy.stats import norm
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = "/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias"

   
def import_simulated(model_name, var_name, member):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	arq = xr.open_dataset('{0}/cmip6/{1}/historical/'.format(dataset_dir, model_name) + '{0}_br_day_{1}_historical_{2}_19860101-20051231_lonlat.nc'.format(var_name, model_name, member))
	data = arq[var_name]
	var = data.sel(time=slice('1986-01-01','2005-12-31'))
	day = var.values
	ts = np.nanmean(np.nanmean(day, axis=1), axis=1)

	return ts


def import_projected(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	arq = xr.open_dataset('{0}/cmip6/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('2081-01-01','2100-12-31'))
	day = var.values
	ts = np.nanmean(np.nanmean(day, axis=1), axis=1)

	return ts


def import_projected_correct(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	arq = xr.open_dataset('{0}/cmip6_correct/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_corrected.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('2081-01-01','2100-12-31'))
	day = var.values
	ts = np.nanmean(np.nanmean(day, axis=1), axis=1)

	return ts
    

def compute_pdf(data):

	"""
	The input arrays must have the same dimensions
	:Param data: Numpy array with model or obs data
	:Return: Probability Density Function
	"""

	x = np.linspace(np.min(data), np.max(data))
	y = np.nanmean(x)
	z = np.nanstd(x)
	pdf = norm.pdf(x,y,z)

	return x, pdf
	

# Best models list
best_models = [7, 9, 13, 15, 17]

# Variable dictionary
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}

experiment = 'ssp585'
dt = '20150101-21001231'

for i in best_models:
	for j in range(1, 4):	
		var_cmip6 = var_dict[j][1]
	
		print(cmip6[i][0])
		print(var_cmip6)

		# Import cmip models and obs database 
		sim = import_simulated(cmip6[i][0], var_cmip6, cmip6[i][1])
		proj = import_projected(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
		proj_cor = import_projected_correct(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
		
		print(len(sim))
		print()
		print(len(proj))
		# ~ print(len(proj_cor))
		
		# Import cdf function
		x_sim, pdf_sim = compute_pdf(sim)
		x_proj, pdf_proj = compute_pdf(proj)
		x_proj_cor, pdf_proj_cor = compute_pdf(proj_cor)

		# Plot figure
		fig = plt.figure()

		if var_cmip6 == 'pr':
			legend = 'Precipitação (mm d⁻¹)'
		elif var_cmip6 == 'tasmax':
			legend = 'Temperatura máxima (°C)'
		else:
			legend = 'Temperatura mínima(°C)'
					
		ax = fig.add_subplot(1, 1, 1)
		plt.plot(x_sim, pdf_sim, color='black', label='Simulação', linewidth=1.5)
		plt.plot(x_proj, pdf_proj,  color='red', label='Projeção', linewidth=1.5)
		plt.plot(x_proj, pdf_proj,  color='red', label='Projeção corrigida', linestyle='--', linewidth=1.5)  
		plt.title('(a) PDF diário {0}'.format(cmip6[i][0]), loc='left', fontweight='bold')
		plt.xlabel('{0}'.format(legend), fontweight='bold')
		plt.ylabel('CDF', fontweight='bold')
		if var_cmip6 == 'pr':
			plt.xticks(np.arange(0, 32, 2))
			plt.xlim(0, 30)
		elif var_cmip6 == 'tasmax':
			plt.xticks(np.arange(20, 42, 2))
			plt.xlim(20, 40)
		else:
			plt.xticks(np.arange(10, 32, 2))
			plt.xlim(10, 30)
		plt.grid(linestyle='--')
		plt.legend(loc=9, ncol=3, frameon=False)

		# Path out to save figure
		path_out = '/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/figs/correct_bias'
		name_out = 'pyplt_boxplot_pdf_correct_bias_cmip6_{0}_{1}_{2}.png'.format(cmip6[i][0], var_cmip6, dt)
		plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
		plt.close('all')
		plt.cla()
exit()
