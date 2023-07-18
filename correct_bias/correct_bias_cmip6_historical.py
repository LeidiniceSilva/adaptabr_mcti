# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import os
import netCDF4
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_cmip6_models_name import cmip6

dataset_dir = "/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/correct_bias"


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
        
    return lat, lon, value	
    
    
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


# Import cmip models and obs database 
best_models = [17, 7, 13, 9, 15]
i = 9

# create date list
time_array = pd.date_range('1986-01-01','1986-12-31', freq='D').strftime('%Y-%m-%d').tolist()

dt = '19860101-20051231'
experiment = 'historical'
var_obs = 'Tmin'
var_cmip6 = 'tasmin'

print(cmip6[i][0])
print(experiment)
print(var_cmip6)

lat, lon, observed = import_observed(var_obs, dt)
lat_array, lon_array, simulated = import_simulated(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)

# Import correct bias function
corrected_dataset = correct_bias(simulated[0:365, :, :], observed[0:365, :, :])

# Print the shape of the corrected array
print(corrected_dataset.shape)

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
nc_output = '{0}/cmip6_correct/{1}/{2}/{3}_br_day_{1}_{2}_{4}_{5}_correct.nc'.format(dataset_dir, cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)
write_3d_nc(nc_output, corrected_dataset, time_array, lat_array, lon_array, var_units, var_shortname, var_longname, time_units, missing_value=-999.)
exit()





