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

from netCDF4 import Dataset
from scipy.stats import gamma
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = "/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias"

# Best models list
best_models = [17, 7, 13, 9, 15]
mdl = 7

def import_observed(var_name, target_date):
	
    """ Import observed data.
    :param: model_name: str (BR-DWGD)
    :param: target_date: datetime
    :return: observed data (Prec, Tasmax and Tasmin)
    :return: flag that it refers to data control
    :rtype: 3D array
    """
    
    file_name = "{0}/obs/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc".format(dataset_dir, var_name, target_date)
    data  = netCDF4.Dataset(file_name)
    var   = data.variables[var_name][:]
    lat   = data.variables['lat'][:]
    lon   = data.variables['lon'][:]
    value = var[:][:,:,:]
        
    return value	
    
    
def import_simulated(model_name, exp_name, var_name, member, target_date):
        
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
    

def correct_bias(simulated_data, observed_data):

	# Get the shape of the simulated data array
	shape = simulated_data.shape

	# Correct each grid point separately
	corrected_data = np.zeros_like(simulated_data)
	for i in range(shape[1]):
		for j in range(shape[2]):
			corrected_data[:, i, j] = correct_bias_1d(simulated_data[:, i, j], observed_data[:, i, j])

	return corrected_data
    

def correct_bias_1d(simulated_values, observed_values):
	
    # Calculate the empirical quantiles of simulated and observed values
    simulated_quantiles = np.percentile(simulated_values, np.arange(0, 101))
    observed_quantiles = np.percentile(observed_values, np.arange(0, 101))
    
    # Find the mapping function by fitting a linear regression
    mapping_function = np.polyfit(simulated_quantiles, observed_quantiles, 1)
    
    # Apply the mapping function to the simulated values
    corrected_values = np.polyval(mapping_function, simulated_values)
    
    return corrected_values
    

def get_leap_day_indices(start_year, end_year):
	
    leap_day_indices = []
    for year in range(start_year, end_year + 1):
        if calendar.isleap(year):
            # Leap day (February 29) corresponds to index 59 in a non-leap year and index 60 in a leap year
            leap_day_indices.append((year - start_year) * 365 + 59 if calendar.monthrange(year, 2)[1] == 29 else (year - start_year) * 365 + 60)
    
    return leap_day_indices


def remove_leap_days(array, indices):
	
    # Create a mask to identify the indices to delete
    mask = np.ones(array.shape[0], dtype=bool)
    mask[indices] = False

    # Use the mask to select the non-deleted indices
    new_array = array[mask]

    return new_array
    

def group_elements(lst, group_size=365):
    
    grouped_list = [lst[i:i + group_size] for i in range(0, len(lst), group_size)]
    
    return grouped_list
    
      
def write_3d_nc(ncname, var_array, time_array, lat_array, lon_array, var_units,
                var_shortname, var_longname, time_units, missing_value=-999.):
	
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

	foo = Dataset(ncname, 'w', format='NETCDF4_CLASSIC')

	foo.Conventions = 'CF-1.7 CMIP-6.0 UGRID-1.0'
	foo.title 		= 'GFDL-ESM4 correct bias model output'
	foo.institution = 'NOAA GFDL, Princeton, NJ 08540, USA'
	foo.source 		= 'GFDL-ESM4'
	foo.history 	= 'CMIP6 model with corrected bias'
	foo.references 	= 'https://esgf.ceda.ac.uk/thredds/catalog/esg_cmip6/catalog.html'
	foo.comment 	= 'This netCDF4 file was corrected using GQM function'

	foo.createDimension('time', None)
	foo.createDimension('lat', lat_array.shape[0])
	foo.createDimension('lon', lon_array.shape[0])

	times = foo.createVariable('time', float, ('time'))
	times.units = time_units
	times.calendar = 'standard'
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


# Create date list non leap days
time = pd.date_range('1986','2006', freq='Y').strftime('%Y').tolist()
time_i = pd.date_range('1986-01-01','2005-12-31', freq='D').strftime('%Y-%m-%d').tolist()
time_i.remove('1988-02-29')
time_i.remove('1992-02-29')
time_i.remove('1996-02-29')
time_i.remove('2000-02-29')
time_i.remove('2004-02-29')
		
dt = '19860101-20051231'
experiment = 'historical'
var_obs = 'pr'
var_cmip6 = 'pr'

print(cmip6[mdl][0])
print(experiment)
print(var_cmip6)

# Import cmip models and obs database 
obs = import_observed(var_obs, dt)
lat_array, lon_array, sim_array = import_simulated(cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt)

print('first step')
print(obs.shape)
print(sim_array.shape)

# leap day indices to delete
# [789, 2249, 3709, 5169, 6629]
leap_day_indices = get_leap_day_indices(1986, 2005)
print(leap_day_indices)

# Remove leap day (Feb 29th) from datasets
if obs.shape[0] == 7305:
	observed = remove_leap_days(obs, leap_day_indices)
else:
	observed = obs
	
if sim_array.shape[0] == 7305:
	simulated = remove_leap_days(sim_array, leap_day_indices)
else:
	simulated = sim_array

print('second step')
print(observed.shape)
print(simulated.shape)

# Call the function with the sample_list and group size 
sample_list = list(range(0, 7300))
grouped_elements = group_elements(sample_list, group_size=365)

# Grouped elements
yr_i = []
yr_f = []
time_ii = []
for gp in grouped_elements:
	yr_i.append(gp[0])
	yr_f.append(gp[-1])
	time_ii.append(time_i[gp[0]:gp[-1]+1])

print('third step')
# Print the grouped elements
for idx in range(0, 20):
	time_array = time_ii[idx]
	print(time[idx], len(time_array), yr_i[idx], yr_f[idx])
		
	# Import correct bias function
	corrected_dataset = correct_bias(simulated[yr_i[idx]:yr_f[idx]+1,:,:], observed[yr_i[idx]:yr_f[idx]+1,:,:])

	# Print the shape of the corrected array
	print(corrected_dataset.shape)
	exit()
	
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

	print('fourth step')	
	# Path out to save netCDF4
	nc_output = '{0}/cmip6_correct/{1}/{2}/{3}_br_day_{1}_{2}_{4}_{5}_correct.nc'.format(dataset_dir, cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], time[idx])
	write_3d_nc(nc_output, corrected_dataset, time_array, lat_array, lon_array, var_units, var_shortname, var_longname, time_units, missing_value=-999.)


