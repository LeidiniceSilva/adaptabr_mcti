# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import os
import numpy as np

from netCDF4 import Dataset
from datetime import datetime

dataset_dir = "/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/correct_bias/"


def import_simulated_data(model_name, exp_name, var_name, member, target_date):
    
    """ Import simulated data.
    :param: model_name: str (GFDL, INM, MPI, MRI and Nor)
    :param: target_date: datetime
    :return: simulated data (Pr, Tasmax and Tasmin)
    :return: flag that it refers to data control
    :rtype: 3D array
    """
    
    dir_modelname = "{0}/cmip6/{1}/{2}".format(dataset_dir, model_name, exp_name)
    file_modelname = "{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc".format(var_name, model_name, exp_name, member, target_date)
    
    input_data = Dataset(os.path.join(dir_modelname, file_modelname))
    var = input_data.variables[var_name][:]
    input_data.close()

    return var


def import_observed_data(model_name, target_date, basin, dates_list):
	
    """ Import observed data.
    :param: model_name: str (BR-DWGD)
    :param: target_date: datetime
    :return: observed data (Pr, Tasmax and Tasmin)
    :return: flag that it refers to data control
    :rtype: 3D array
    """
    
    file_name = "{0}/obs/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc".format(dataset_dir, var_name, target_date)
    
    input_data = Dataset(file_name)
    var = input_data.variables[var_name][:]
    input_data.close()

    return var		   
    

def correct_bias(simulated, observed):
	
    # Get the shape of the simulated precipitation array
    shape = simulated.shape
    
    # Reshape the simulated precipitation array to 2D
    simulated_2d = simulated.reshape((shape[0], -1))
    
    # Correct each grid point separately
    corrected_2d = np.zeros_like(simulated_2d)
    for i in range(simulated_2d.shape[1]):
        corrected_2d[:, i] = correct_bias_1d(simulated_2d[:, i], observed[:, i])
    
    # Reshape the corrected precipitation array back to 3D
    corrected_dataset = corrected_2d.reshape(shape)
    
    return corrected_dataset


def correct_bias_1d(simulated, observed):
	
    # Calculate the empirical quantiles of simulated and observed precipitation
    simulated_quantiles = np.percentile(simulated, np.arange(0, 101))
    observed_quantiles = np.percentile(observed, np.arange(0, 101))
    
    # Find the mapping function by fitting a linear regression
    mapping_function = np.polyfit(simulated_quantiles, observed_quantiles, 1)
    
    # Apply the mapping function to the simulated precipitation
    corrected_dataset = np.polyval(mapping_function, simulated)

    return corrected_dataset
    
    
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


simulated = np.random.normal(10, 2, size=(365, 100, 100))  # Simulated data
observed = np.random.normal(12, 1, size=(365, 100, 100))  # Observed data

print(observed.shape)
exit()

corrected_dataset = correct_bias(simulated, observed)

# Print the shape of the corrected precipitation array
print(corrected_dataset.shape)
exit()


