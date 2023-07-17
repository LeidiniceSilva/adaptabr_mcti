# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import os
import numpy as np
import matplotlib.pyplot as plt

from pylab import setp
from netCDF4 import Dataset
from datetime import datetime
from scipy.stats import norm

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
	
	
observed = np.random.normal(12, 1, size=(365, 100, 100))  # Observed data
simulated = np.random.normal(10, 2, size=(365, 100, 100))  # Simulated data
corrected_dataset = correct_bias(simulated, observed)

obs = np.nanmean(np.nanmean(observed, axis=1), axis=1)
sim = np.nanmean(np.nanmean(simulated, axis=1), axis=1)
correct_bias = np.nanmean(np.nanmean(corrected_dataset, axis=1), axis=1)

# Print the shape of the corrected array
print(corrected_dataset.shape)

day_boxplot = [obs, sim, correct_bias]

# Import cdf function
x_obs, cdf_obs = compute_cdf(obs)
x_sim, cdf_sim = compute_cdf(sim)
x_correct_bias, cdf_correct_bias = compute_cdf(correct_bias)

# Plot figure
fig = plt.figure()

ax = fig.add_subplot(1, 2, 1)
x = np.arange(1, 3 + 1)
bp = plt.boxplot(day_boxplot, positions=[1, 2, 3], sym='.')
setBoxColors(bp)
plt.title('(a) Daily boxplot', loc='left', fontweight='bold')
plt.xticks(x, ('Observed','Simulated','Correct bias'))
plt.xlabel('Dataset', fontweight='bold')
plt.grid(linestyle='--')

ax = fig.add_subplot(1, 2, 2)
plt.plot(x_obs, cdf_obs, color='black', label='Observed', linewidth=1.5)
plt.plot(x_sim, cdf_sim,  color='red', label='Simulated', linewidth=1.5)
plt.plot(x_correct_bias, cdf_correct_bias,  color='blue', label='Correct bias', linewidth=1.5)
plt.title('(b) Daily CDF', loc='left', fontweight='bold')
plt.grid(linestyle='--')

plt.show()
exit()
