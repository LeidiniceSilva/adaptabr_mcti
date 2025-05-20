# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import os
import cftime
import netCDF4
import calendar
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import scipy.stats as stats

from netCDF4 import Dataset
from dict_cmip6_models_name import cmip6

warnings.filterwarnings("ignore")

# Dataset directory
dataset_dir = '/afs/ictp.it/home/m/mda_silv/Documents/projects/AdaptaBrasil_MCTI/database/correct_bias'

# Best models list
best_models = [7, 9, 13, 15, 17]
mdl = 9

# Variable dictionary
var = 1
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}

dt = '20150101-21001231'
experiment = 'ssp585'
var_obs = var_dict[var][0]
var_cmip6 = var_dict[var][1]

print(cmip6[mdl][0])
print(var_cmip6)


def import_observed(var_name):
	
    """ Import observed data.
    :param: model_name: str (BR-DWGD)
    :param: target_date: datetime
    :return: observed data (Prec, Tasmax and Tasmin)
    :return: flag that it refers to data control
    :rtype: 3D array
    """
    
    value = xr.open_dataset('{0}/obs/{1}_19860101-20051231_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc'.format(dataset_dir, var_name))
        
    return value	
    

def import_simulated(model_name, var_name, member):
        
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	dir_model = "{0}/cmip6/{1}/historical".format(dataset_dir, model_name)
	value = xr.open_dataset('{0}/{1}_br_day_{2}_historical_{3}_19860101-20051231_lonlat.nc'.format(dir_model, var_name, model_name, member))
	
	return value
	
	
def import_projected(model_name, exp_name, var_name, member, target_date):
        
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	dir_model = "{0}/cmip6/{1}/{2}".format(dataset_dir, model_name, exp_name)
	value = xr.open_dataset('{0}/{1}_br_day_{2}_{3}_{4}_{5}_lonlat.nc'.format(dir_model, var_name, model_name, exp_name, member, target_date))
	
	return value 
	

def quantile_mapping_i(observed, model, future, variable_obs, variable_mdl):
	
    # Calculate quantiles for observed and model data
    obs_quantiles = observed[variable_obs].quantile([0.1, 0.5, 0.9], dim='time')
    model_quantiles = model[variable_mdl].quantile([0.1, 0.5, 0.9], dim='time')

    # Calculate quantile mapping correction factors
    correction_factors = {}
    
    for quantile in [0.1, 0.5, 0.9]:
        obs_q = obs_quantiles.sel(quantile=quantile)
        model_q = model_quantiles.sel(quantile=quantile)
        correction_factors[quantile] = obs_q / model_q

    # Apply quantile mapping to correct the future projections
    corrected_future = future[variable_mdl] * correction_factors[0.5]

    return corrected_future


# ~ def quantile_mapping_ii(observed, model, variable_obs, variable_mdl):

	# ~ # Calculate quantiles for observed and model data
	# ~ obs_quantiles = observed[variable_obs].quantile([0.1, 0.5, 0.9], dim='time')
	# ~ model_quantiles = model[variable_mdl].quantile([0.1, 0.5, 0.9], dim='time')
	
	# ~ mapping_functions = {}
	# ~ for quantile in [0.1, 0.5, 0.9]:
		# ~ obs_values = obs_quantiles.sel(quantile=quantile)
		# ~ model_values = model_quantiles.sel(quantile=quantile)
		# ~ mapping_functions[quantile] = stats.linregress(model_values, obs_values)
	
	# ~ corrected_var = model_var.copy()
	# ~ for i in range(len(model_var.time)):
		# ~ for j in range(len(model_var.lat)):
			# ~ for k in range(len(model_var.lon)):
				# ~ for quantile in [0.1, 0.5, 0.9]:
					# ~ slope, intercept, r_value, _, _ = mapping_functions[quantile]
					# ~ corrected_var[i, j, k] = (model_var[i, j, k] - intercept) / slope

	# ~ # Save the corrected data 
	# ~ corrected_data = model_data.copy()
	# ~ corrected_data[variable_mdl] = corrected_var
	
	# ~ return corrected_data

    
# Import cmip models and obs database
observed = import_observed(var_obs)
simulated = import_simulated(cmip6[mdl][0], var_cmip6, cmip6[mdl][1])
projected = import_projected(cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt)

corrected_future = quantile_mapping_i(observed, simulated, projected, var_obs, var_cmip6)
print(np.min(corrected_future))
print(np.max(corrected_future))

corrected_future.to_netcdf('{0}/cmip6_correct/{1}/{2}/{3}_br_day_{1}_{2}_{4}_{5}_corrected.nc'.format(dataset_dir, cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt))
exit()
