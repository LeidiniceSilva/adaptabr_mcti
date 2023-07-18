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

dataset_dir = "/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/correct_bias/"


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
    var = data.sel(time=slice('1986-01-01','1986-12-31'))
    
    clim_mean = var.groupby('time.month').mean('time')
    clim_mean = clim_mean.values
    clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)
    
    mon_mean = var.resample(time='1M').mean()
    mon_mean = mon_mean.values
    mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)
    
    day_mean = var.values
    day = np.nanmean(np.nanmean(day_mean, axis=1), axis=1)

    return clim, mon, day
   
   
def import_simulated(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	clim = []
	mon = []
	day = []

	arq = xr.open_dataset('{0}/cmip6/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('1986-01-01','1986-12-31'))

	clim_mean = var.groupby('time.month').mean('time')
	clim_mean = clim_mean.values
	clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)

	mon_mean = var.resample(time='1M').mean()
	mon_mean = mon_mean.values
	mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)

	day_mean = var.values
	day = np.nanmean(np.nanmean(day_mean, axis=1), axis=1)

	return clim, mon, day


def import_simulated_correct(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	clim = []
	mon = []
	day = []

	arq = xr.open_dataset('{0}/cmip6_correct/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_correct.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('1986-01-01','1986-12-31'))

	clim_mean = var.groupby('time.month').mean('time')
	clim_mean = clim_mean.values
	clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)

	mon_mean = var.resample(time='1M').mean()
	mon_mean = mon_mean.values
	mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)

	day_mean = var.values
	day = np.nanmean(np.nanmean(day_mean, axis=1), axis=1)

	return clim, mon, day


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
	

# Import cmip models and obs database 
best_models = [17, 7, 13, 9, 15]
i = 9

dt = '19860101-20051231'
experiment = 'historical'
var_obs = 'Tmin'
var_cmip6 = 'tasmin'

print(cmip6[i][0])
print(experiment)
print(var_cmip6)

clim_obs, mon_obs, day_obs = import_observed(var_obs, dt)
clim_sim, mon_sim, day_sim = import_simulated(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
clim_sim_correct, mon_sim_correct, day_sim_correct = import_simulated_correct(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)

day_boxplot = [day_obs, day_sim, day_sim_correct]

# Import cdf function
x_obs, cdf_obs = compute_cdf(day_obs)
x_sim, cdf_sim = compute_cdf(day_sim)
x_correct_bias, cdf_correct_bias = compute_cdf(day_sim_correct)

# Plot figure
fig = plt.figure(figsize=(10, 11))

if var_cmip6 == 'pr':
	legend = 'Precipitation (mm d⁻¹)'
elif var_cmip6 == 'tasmax':
	legend = 'Maximum Air Temperature (°C)'
else:
	legend = 'Minimum Air Temperature (°C)'
	
ax = fig.add_subplot(2, 2, 1)
x = np.arange(1, 3 + 1)
bp = plt.boxplot(day_boxplot, positions=[1, 2, 3], sym='.')
setBoxColors(bp)
plt.title('(a) Daily boxplot', loc='left', fontweight='bold')
plt.xticks(x, ('Observed','Simulated','Correction'))
plt.xlabel('Dataset', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.grid(linestyle='--')

ax = fig.add_subplot(2, 2, 2)
plt.plot(x_obs, cdf_obs, color='black', label='Observed', linewidth=1.5)
plt.plot(x_sim, cdf_sim,  color='red', label='Simulated', linewidth=1.5)
plt.plot(x_correct_bias, cdf_correct_bias,  color='blue', label='Correction', linewidth=1.5)
plt.xlabel('{0}'.format(legend), fontweight='bold')
plt.ylabel('CDF', fontweight='bold')
plt.title('(b) Daily CDF', loc='left', fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=1, fontsize=8, frameon=False)

ax = fig.add_subplot(2, 2, 3)
time = np.arange(0.5, 12 + 0.5)
plt.plot(time, clim_obs, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='black', label = 'Observed')
plt.plot(time, clim_sim, linewidth=1.5, linestyle='--', markersize=3, marker='*', markerfacecolor='white', color='red', label = 'Simulated')
plt.plot(time, clim_sim_correct, linewidth=1.5, linestyle='--', markersize=3, marker='^', markerfacecolor='white', color='blue', label = 'Correction')
plt.title('(c) Anual cycle', loc='left', fontweight='bold')
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.xlabel('Months', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=3, fontsize=8, frameon=False)

ax = fig.add_subplot(2, 2, 4)
mon_dt = pd.date_range(start="19860101", end="19861231", freq="M")
obs_dt = pd.Series(data=mon_obs, index=mon_dt)
sim_dt = pd.Series(data=mon_sim, index=mon_dt)
sim_correct_dt = pd.Series(data=mon_sim_correct, index=mon_dt)

plt.plot(obs_dt, linewidth=1.5, linestyle='--', color='black', label = 'Observed')
plt.plot(sim_dt, linewidth=1.5, linestyle='--', color='red', label = 'Simulated')
plt.plot(sim_correct_dt, linewidth=1.5, linestyle='--', color='blue', label = 'Correction')
plt.title('(d) Monthly time series', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.xlabel('Period', fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=3, fontsize=8, frameon=False)

# Path out to save figure
path_out = '/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/figs/correct_bias'
name_out = 'pyplt_skill_correct_bias_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
