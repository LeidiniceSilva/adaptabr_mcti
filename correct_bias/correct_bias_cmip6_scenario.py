# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import os
import netCDF4
import calendar
import numpy as np
import pandas as pd
import xarray as xr

from netCDF4 import Dataset
from scipy.signal import detrend
from scipy.stats import gamma, norm
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias'

# Best models list
best_models = [7, 9, 13, 15, 17]
mdl = 17

# Variable dictionary
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}
var = 1

experiment = 'ssp126'
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
    
    file_name = "{0}/obs/{1}_1986_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc".format(dataset_dir, var_name)
    data  = netCDF4.Dataset(file_name)
    var   = data.variables[var_name][:]
    lat   = data.variables['lat'][:]
    lon   = data.variables['lon'][:]
    value = var[:][:,:,:]
        
    return value	
    

def import_simulated(model_name, var_name):
        
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	dir_model = "{0}/cmip6/{1}/historical".format(dataset_dir, model_name)
	file_model = "{0}/{1}_br_day_NorESM2-MM_historical_r1i1p1f1_gn_1986_lonlat.nc".format(dir_model, var_name)
	data  = netCDF4.Dataset(file_model)
	var   = data.variables[var_name][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	return lat, lon, value 
	
	
def import_projected(model_name, exp_name, var_name, member, target_date):
        
	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""
	
	dir_model = "{0}/cmip6/{1}/{2}".format(dataset_dir, model_name, exp_name)
	file_model = "{0}/{1}_br_day_{2}_{3}_{4}_{5}_lonlat.nc".format(dir_model, var_name, model_name, exp_name, member, target_date)
	data  = netCDF4.Dataset(file_model)
	var   = data.variables[var_name][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	
	return lat, lon, value 
	

def remove_leap_days(array, indices):
	
    # Create a mask to identify the indices to delete
    mask = np.ones(array.shape[0], dtype=bool)
    mask[indices] = False

    # Use the mask to select the non-deleted indices
    new_array = array[mask]

    return new_array
    

def gamma_correction(obs_data, mod_data, sce_data, lower_limit=0.1, cdf_threshold=0.9999999):
	
    obs_raindays, mod_raindays, sce_raindays = [
        x[x >= lower_limit] for x in [obs_data, mod_data, sce_data]
    ]
    obs_gamma, mod_gamma, sce_gamma = [
        gamma.fit(x) for x in [obs_raindays, mod_raindays, sce_raindays]
    ]

    obs_cdf = gamma.cdf(np.sort(obs_raindays), *obs_gamma)
    mod_cdf = gamma.cdf(np.sort(mod_raindays), *mod_gamma)
    sce_cdf = gamma.cdf(np.sort(sce_raindays), *sce_gamma)

    obs_cdf[obs_cdf > cdf_threshold] = cdf_threshold
    mod_cdf[mod_cdf > cdf_threshold] = cdf_threshold
    sce_cdf[sce_cdf > cdf_threshold] = cdf_threshold

    obs_cdf_intpol = np.interp(
        np.linspace(1, len(obs_raindays), len(sce_raindays)),
        np.linspace(1, len(obs_raindays), len(obs_raindays)),
        obs_cdf,
    )

    mod_cdf_intpol = np.interp(
        np.linspace(1, len(mod_raindays), len(sce_raindays)),
        np.linspace(1, len(mod_raindays), len(mod_raindays)),
        mod_cdf,
    )

    obs_inverse, mod_inverse, sce_inverse = [
        1.0 / (1.0 - x) for x in [obs_cdf_intpol, mod_cdf_intpol, sce_cdf]
    ]

    adapted_cdf = 1 - 1.0 / (obs_inverse * sce_inverse / mod_inverse)
    adapted_cdf[adapted_cdf < 0.0] = 0.0

    initial = (
        gamma.ppf(np.sort(adapted_cdf), *obs_gamma)
        * gamma.ppf(sce_cdf, *sce_gamma)
        / gamma.ppf(sce_cdf, *mod_gamma)
    )

    mod_frequency = 1.0 * mod_raindays.shape[0] / mod_data.shape[0]
    sce_frequency = 1.0 * sce_raindays.shape[0] / sce_data.shape[0]

    days_min = len(sce_raindays) * sce_frequency / mod_frequency

    expected_sce_raindays = int(min(days_min, len(sce_data)))

    sce_argsort = np.argsort(sce_data)
    correction = np.zeros(len(sce_data))

    if len(sce_raindays) > expected_sce_raindays:
        initial = np.interp(
            np.linspace(1, len(sce_raindays), expected_sce_raindays),
            np.linspace(1, len(sce_raindays), len(sce_raindays)),
            initial,
        )
    else:
        initial = np.hstack(
            (np.zeros(expected_sce_raindays - len(sce_raindays)), initial)
        )

    correction[sce_argsort[:expected_sce_raindays]] = initial
    # correction = pd.Series(correction, index=sce_data.index)
    
    return correction
    
      
def write_3d_nc(ncname, var_array, time_array, lat_array, lon_array, var_units, var_shortname, var_longname, time_units, missing_value=-999.):
	
	"""Write 3 dimensional (time, lat, lon) netCDF4
	:param ncname: Name of the output netCDF
	:param var_array: 3D array (time, lat, lon)
	:param time_array: 1D array with unique value
	:param lat_array: 1D array which contains the lat values
	:param lon_array: 1D array which contains the lon values
	:param var_units: variable units
	:param var_shortname: variable short name
	:param var_longname: variable long name
	:param time_units: time units
	:param missing_value: float missing value
	"""

	if cmip6[mdl][0] == 'GFDL-ESM4':
		cmip6_inst = 'NOAA GFDL, Princeton, NJ 08540, USA'
	elif cmip6[mdl][0] == 'INM-CM5-0':
		cmip6_inst = 'INM, Russian Academy of Science, Moscow 119991, Russia'
	elif cmip6[mdl][0] == 'MPI-ESM1-2-HR':
		cmip6_inst = 'MPI for Meteorology, Hamburg 20146, Germany'
	elif cmip6[mdl][0] == 'MRI-ESM2-0':
		cmip6_inst = 'MRI, Tsukuba, Ibaraki 305-0052, Japan'
	else:
		cmip6_inst = 'NCC, c/o MET-Norway, Henrik Mohns plass 1, Oslo 0313, Norway'
		
	foo = Dataset(ncname, 'w', format='NETCDF4_CLASSIC')

	foo.Conventions = 'CF-1.7 CMIP-6.0 UGRID-1.0'
	foo.title 		= '{0} correct bias model output'.format(cmip6[mdl][0])
	foo.institution = '{0}'.format(cmip6_inst)
	foo.source 		= '{0}'.format(cmip6[mdl][0])
	foo.history 	= 'CMIP6 model with corrected bias'
	foo.references 	= 'https://esgf.ceda.ac.uk/thredds/catalog/esg_cmip6/catalog.html'
	foo.comment 	= 'This netCDF4 file was corrected using EQM function'

	foo.createDimension('time', None)
	foo.createDimension('lat', lat_array.shape[0])
	foo.createDimension('lon', lon_array.shape[0])

	times = foo.createVariable('time', float, ('time'))
	times.units = time_units
	times.calendar = '365_day'
	times[:] = range(len(time_array))

	laty = foo.createVariable('lat', 'f4', 'lat')
	laty.units = 'degrees_north'
	laty.long_name = 'latitude'
	laty[:] = lat_array

	lonx = foo.createVariable('lon', 'f4', 'lon')
	lonx.units = 'degrees_east'
	lonx.long_name = 'longitude'
	lonx[:] = lon_array

	var = foo.createVariable(var_shortname, float, ('time', 'lat', 'lon'))
	var.units = var_units
	var.long_name = var_longname
	var.missing_value = missing_value
	var[:] = var_array

	foo.close()


# Import cmip models and obs database 
obs = import_observed(var_obs)
sim = import_simulated(cmip6[mdl][0], var_cmip6)

for dt in range(2015, 2101):
	print(dt)

	time_i = pd.date_range('{0}-01-01'.format(dt),'{0}-12-31'.format(dt), freq='D').strftime('%Y-%m-%d').tolist()
	time_array = time_i
	print(len(time_array))
	
	lat_array, lon_array, proj_array = import_projected(cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt)

	# leap day indices to delete
	leap_day_indices = [59]

	# Remove leap day (Feb 29th) from datasets
	if proj_array.shape[0] == 366:
		projected = remove_leap_days(proj_array, leap_day_indices)
	else:
		projected = proj_array
		
	# Import correct bias function
	corrected_dataset = gamma_correction(obs, sim, projected)
	
	# ~ # To replace negative values with 0
	# ~ if var_cmip6 == 'pr':
		# ~ corrected_dataset = np.where(corrected_data<0, 0, corrected_data)
	# ~ else:
		# ~ corrected_dataset = corrected_data
		
	if var_cmip6 == 'pr':
		var_units = 'mm'
		var_shortname = 'pr'
		var_longname = 'Daily Total Precipitation'
		time_units = 'days since {}'.format(time_array[0])
	elif var_cmip6 == 'tasmax':
		var_units = 'Celcius degrees'
		var_shortname = 'tasmax'
		var_longname = 'Daily Maximum Near-Surface Air Temperature'
		time_units = 'days since {}'.format(time_array[0])
	else:
		var_units = 'Celcius degrees'
		var_shortname = 'tasmin'
		var_longname = 'Daily Minimum Near-Surface Air Temperature'
		time_units = 'days since {}'.format(time_array[0])

	# Path out to save netCDF4
	nc_output = '{0}/cmip6_correct/{1}/{2}/{3}_br_day_{1}_{2}_{4}_{5}_correct.nc'.format(dataset_dir, cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt)
	write_3d_nc(nc_output, corrected_dataset, time_array, lat_array, lon_array, var_units, var_shortname, var_longname, time_units, missing_value=-999.)
	exit()
