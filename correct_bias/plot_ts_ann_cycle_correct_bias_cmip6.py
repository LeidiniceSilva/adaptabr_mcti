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

dataset_dir = "/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias"


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

    mon_mean = var.resample(time='1M').mean()
    mon_mean = mon_mean.values
    mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)

    clim_mean = var.groupby('time.month').mean('time')
    clim_mean = clim_mean.values
    clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)
    
    return mon, clim
   
   
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
	var = data.sel(time=slice('1986-01-01','1986-12-31'))
	
	mon_mean = var.resample(time='1M').mean()
	mon_mean = mon_mean.values
	mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)

	clim_mean = var.groupby('time.month').mean('time')
	clim_mean = clim_mean.values
	clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)

	return mon, clim


def import_simulated_correct(model_name, exp_name, var_name, member, target_date):
       
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	arq = xr.open_dataset('{0}/cmip6_correct/{1}/{2}/'.format(dataset_dir, model_name, exp_name) + '{0}_br_day_{1}_{2}_{3}_{4}_correct.nc'.format(var_name, model_name, exp_name, member, target_date))
	data = arq[var_name]
	var = data.sel(time=slice('1986-01-01','1986-12-31'))
	
	mon_mean = var.resample(time='1M').mean()
	mon_mean = mon_mean.values
	mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)

	clim_mean = var.groupby('time.month').mean('time')
	clim_mean = clim_mean.values
	clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)

	return mon, clim


# Import cmip models and obs database 
best_models = [17, 7, 13, 9, 15]
i = 7

target_dt = '1986'
dt = '19860101-20051231'
experiment = 'historical'
var_obs = 'pr'
var_cmip6 = 'pr'

print(cmip6[i][0])
print(experiment)
print(var_cmip6)

mon_obs, clim_obs = import_observed(var_obs, dt)
mon_sim, clim_sim = import_simulated(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
mon_sim_correct, clim_sim_correct = import_simulated_correct(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], target_dt)

# Plot figure
fig = plt.figure(figsize=(10, 8))

if var_cmip6 == 'pr':
	legend = 'Precipitação (mm d⁻¹)'
elif var_cmip6 == 'tasmax':
	legend = 'Temperatura máxima (°C)'
else:
	legend = 'Temperatura mínima(°C)'

ax = fig.add_subplot(2, 1, 1)
time = np.arange(0.5, 12 + 0.5)
plt.plot(time, clim_obs, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='black', label = 'Observação')
plt.plot(time, clim_sim, linewidth=1.5, linestyle='--', markersize=3, marker='s', markerfacecolor='white', color='red', label = 'Simulação')
plt.plot(time, clim_sim_correct, linewidth=1.5, linestyle='--', markersize=3, marker='^', markerfacecolor='white', color='blue', label = 'Correção')
plt.title('(a) Ciclo anual', loc='left', fontweight='bold')
plt.xticks(time, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'))
plt.xlabel('Meses', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=3, fontsize=8, frameon=False)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2))
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 45, 3))
	plt.ylim(18, 42)
else:
	plt.yticks(np.arange(9, 36, 3))
	plt.ylim(9, 33)
	
ax = fig.add_subplot(2, 1, 2)
mon_dt = pd.date_range(start="19860101", end="19861231", freq="M")
obs_dt = pd.Series(data=mon_obs, index=mon_dt)
sim_dt = pd.Series(data=mon_sim, index=mon_dt)
sim_correct_dt = pd.Series(data=mon_sim_correct, index=mon_dt)

plt.plot(obs_dt, linewidth=1.5, linestyle='--', color='black', label = 'Observação')
plt.plot(sim_dt, linewidth=1.5, linestyle='--', color='red', label = 'Simulação')
plt.plot(sim_correct_dt, linewidth=1.5, linestyle='--', color='blue', label = 'Correção')
plt.title('(b) Série temporal mensal', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.xlabel('Período', fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=3, frameon=False)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2))
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 45, 3))
	plt.ylim(18, 42)
else:
	plt.yticks(np.arange(9, 36, 3))
	plt.ylim(9, 33)
	
# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/correct_bias'
name_out = 'pyplt_ts_ann_cycle_correct_bias_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
