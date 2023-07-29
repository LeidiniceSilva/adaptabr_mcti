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
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = "/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/correct_bias"

# Best models list
best_models = [17, 7, 13, 9, 15]
i = 9

experiment = 'historical'
dt = '19860101-20051231'
var_obs = 'Tmin'
var_cmip6 = 'tasmin'

print(cmip6[i][0])
print(experiment)
print(var_cmip6)


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
    var = data.sel(time=slice('1986-01-01','2005-12-31'))
    
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
	var = data.sel(time=slice('1986-01-01','2005-12-31'))

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
	var = data.sel(time=slice('1986-01-01','2005-12-31'))

	mon_mean = var.resample(time='1M').mean()
	mon_mean = mon_mean.values
	mon = np.nanmean(np.nanmean(mon_mean, axis=1), axis=1)

	clim_mean = var.groupby('time.month').mean('time')
	clim_mean = clim_mean.values
	clim = np.nanmean(np.nanmean(clim_mean, axis=1), axis=1)

	return mon, clim


# Import cmip models and obs database 
mon_obs, clim_obs = import_observed(var_obs, dt)
mon_sim, clim_sim = import_simulated(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
mon_sim_correct, clim_sim_correct = import_simulated_correct(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)

# Plot figure
fig = plt.figure(figsize=(8, 8))

if var_cmip6 == 'pr':
	legend = 'Precipitação (mm d⁻¹)'
elif var_cmip6 == 'tasmax':
	legend = 'Temperatura máxima (°C)'
else:
	legend = 'Temperatura mínima(°C)'

ax = fig.add_subplot(2, 1, 1)
time = np.arange(0.5, 240 + 0.5)
plt.plot(time, mon_obs, linewidth=1.5, linestyle='--', color='black', label = 'Observação')
plt.plot(time, mon_sim, linewidth=1.5, linestyle='--', color='red', label = 'Simulação')
plt.plot(time, mon_sim_correct, linewidth=1.5, linestyle='--', color='blue', label = 'Correção')
plt.title('(a) Série temporal mensal - {0}'.format (cmip6[i][0]), loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.xlabel('Período 01/1986 - 12/2005', fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=3, frameon=False)
plt.xticks(np.arange(0, 252, 12))
plt.xlim(0, 240)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2))
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 38, 2))
	plt.ylim(18, 36)
else:
	plt.yticks(np.arange(14, 30, 2))
	plt.ylim(14, 28)
	
ax = fig.add_subplot(2, 1, 2)
time = np.arange(0.5, 12 + 0.5)
plt.plot(time, clim_obs, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='black', label = 'Observação')
plt.plot(time, clim_sim, linewidth=1.5, linestyle='--', markersize=3, marker='s', markerfacecolor='white', color='red', label = 'Simulação')
plt.plot(time, clim_sim_correct, linewidth=1.5, linestyle='--', markersize=3, marker='^', markerfacecolor='white', color='blue', label = 'Correção')
plt.title('(b) Ciclo anual', loc='left', fontweight='bold')
plt.xticks(time, ('Jan', 'Fev', 'Mar', 'Abr', 'Mai', 'Jun', 'Jul', 'Ago', 'Set', 'Out', 'Nov', 'Dez'))
plt.xlabel('Meses', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.grid(linestyle='--')
plt.legend(loc=9, ncol=3, frameon=False)
if var_cmip6 == 'pr':
	plt.yticks(np.arange(0, 18, 2))
	plt.ylim(0, 16)
elif var_cmip6 == 'tasmax':
	plt.yticks(np.arange(18, 38, 2))
	plt.ylim(18, 36)
else:
	plt.yticks(np.arange(14, 30, 2))
	plt.ylim(14, 28)	

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/correct_bias'
name_out = 'pyplt_ts_ann_cycle_correct_bias_cmip6_{0}_{1}_{2}.png'.format(cmip6[i][0], var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
