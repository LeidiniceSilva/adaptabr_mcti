# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import os
import netCDF4
import numpy as np

from dict_cmip6_models_name import cmip6

dataset_dir = "/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/correct_bias"


def import_observed_data(var_name, target_date):
	
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
    
    
def import_simulated_data(model_name, exp_name, var_name, member, target_date):
        
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
    :type ncname: string object
    :param var_array: 3d prec array (time, lat, lon)
    :type var_array: numpy array
    :param time_array: 1d array with unique value which is the prediction hour
    :type time_array: numpy array
    :param lat_array: 1d array which contains the lat values
    :type lat_array: numpy array
    :param lon_array: 1d array which contains the lon values
    :type lon_array: numpy array
    :param var_units: variable units
    :type var_units: string object
    :param var_shortname: variable short name
    :type var_shortname: string object
    :param var_longname: variable long name
    :type var_longname: string object
    :param time_units: time units
    :type time_units: string object
    :param missing_value: missing value
    :type missing_value: float
    """

    foo = Dataset(ncname, 'w', format='NETCDF4_CLASSIC')

    foo.createDimension('time', time_array.shape[0])
    foo.createDimension('lat', lat_array.shape[0])
    foo.createDimension('lon', lon_array.shape[0])

    times = foo.createVariable('time', 'f8', 'time')
    times.units = time_units
    times.calendar = 'standard'
    times[:] = time_array

    laty = foo.createVariable('lat', 'f4', 'lat')
    laty.units = 'degrees_north'
    laty.long_name = 'latitude'
    laty[:] = lat_array

    lonx = foo.createVariable('lon', 'f4', 'lon')
    lonx.units = 'degrees_east'
    lonx.long_name = 'longitude'
    lonx[:] = lon_array

    var = foo.createVariable(var_shortname, 'f', ('time', 'lat', 'lon'))
    var.units = var_units
    var.long_name = var_longname
    var.missing_value = missing_value
    var[:] = var_array

    foo.close()
	

# Import cmip models and obs database 
best_models = [17, 7, 13, 9, 15]
i = 9

experiment = 'historical'
if experiment == 'historical':
	dt = '19860101-20051231'
else:
	dt = '20060101-21001231'
	
var_obs = 'pr'
var_cmip6 = 'pr'

print(cmip6[i][0])
print(experiment)
print(var_cmip6)

lat, lon, observed = import_observed_data(var_obs, dt)
lat, lon, simulated = import_simulated_data(cmip6[i][0], experiment, var_cmip6, cmip6[i][1], dt)

print(simulated.shape)
exit()

# Import correct bias function
corrected_dataset = correct_bias(simulated, observed)

# Print the shape of the corrected array
print(corrected_dataset.shape)

# Path out to save figure
path_out = '{0}/cmip6_correct/{1}/{2}'.format(dataset_dir, cmip6[i][0], experiment)
name_out = '{0}_br_day_{1}_{2}_{3}_{4}_correct_bias.nc'.format(var_cmip6, cmip6[i][0], experiment, member, dt)

write_3d_nc(ncname, var_array, time_array, lat_array, lon_array, var_units,
                var_shortname, var_longname, time_units, missing_value=-999.)
                
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()





