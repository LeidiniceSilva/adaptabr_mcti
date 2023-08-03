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
from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = "/media/nice/Nice/documentos/projetos/AdaptaBrasil_MCTI/database/correct_bias"

# Best models list
best_models = [7, 9, 13, 15, 17]
mdl = 7

# Scenario dictionary
ssp_dict = {1 :['ssp126'], 2 :['ssp245'], 3 :['ssp485']}
ssp = 1

# Variable dictionary
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}
var = 1

dt0 = '19860101-20051231'
dt1 = '20150101-21001231'
experiment = ssp_dict[ssp]
var_obs = var_dict[var][0]
var_cmip6 = var_dict[var][1]

print(cmip6[mdl][0])
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
    

def correct_bias(model_data, observed_data):
	
    # Calculate the CDF of model and observed data
    model_cdf = model_data.rank(dim='time') / float(model_data.time.size)
    observed_cdf = observed_data.rank(dim='time') / float(observed_data.time.size)

    # Compute the quantile mapping correction function
    correction_function = xr.interp(model_cdf, observed_cdf.quantile(dim='time', q=model_cdf), observed_data)

    # Apply the correction to the model data
    corrected_model_data = xr.where(model_data.notnull(), correction_function, np.nan)

    return corrected_model_data
    

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


# Create date list non leap days
time = pd.date_range('2015','2100', freq='Y').strftime('%Y').tolist()
time_i = pd.date_range('2015-01-01','2100-12-31', freq='D').strftime('%Y-%m-%d').tolist()
time_i.remove('2016-02-29')
time_i.remove('2020-02-29')
time_i.remove('2024-02-29')
time_i.remove('2028-02-29')
time_i.remove('2032-02-29')
time_i.remove('2036-02-29')
time_i.remove('2040-02-29')
time_i.remove('2044-02-29')
time_i.remove('2048-02-29')
time_i.remove('2052-02-29')
time_i.remove('2056-02-29')
time_i.remove('2060-02-29')
time_i.remove('2064-02-29')
time_i.remove('2068-02-29')
time_i.remove('2072-02-29')
time_i.remove('2076-02-29')
time_i.remove('2080-02-29')
time_i.remove('2084-02-29')
time_i.remove('2088-02-29')
time_i.remove('2092-02-29')
time_i.remove('2096-02-29')

# Import cmip models and obs database 
obs = import_observed(var_obs, dt0)
lat_array, lon_array, sim_array = import_simulated(cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt1)

print('first step')
print(obs.shape)
print(sim_array.shape)

# leap day indices to delete
leap_day = get_leap_day_indices(2015, 2100)
leap_day_indices = [424, 1884, 3344, 4804, 6264, 7724, 9184, 10644, 12104, 13564, 15024, 16484, 17944, 19404, 20864, 22324, 23784, 25244, 26704, 28164, 29624]

# Remove leap day (Feb 29th) from datasets
if obs.shape[0] == 31411:
	observed = remove_leap_days(obs, leap_day_indices)
else:
	observed = obs
	
if sim_array.shape[0] == 31411:
	simulated = remove_leap_days(sim_array, leap_day_indices)
else:
	simulated = sim_array

print('second step')
print(observed.shape)
print(simulated.shape)

# Call the function with the sample_list and group size 
sample_list = list(range(0, 31390))
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
for idx in range(0, 86):
	time_array = time_ii[idx]
	print(time[idx], len(time_array), yr_i[idx], yr_f[idx])
		
	# Import correct bias function
	corrected_data = correct_bias(simulated[yr_i[idx]:yr_f[idx]+1,:,:], observed[yr_i[idx]:yr_f[idx]+1,:,:])

	# Print the shape of the corrected array
	print(corrected_data.shape)
	exit()
	
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

	print('fourth step')	
	# Path out to save netCDF4
	nc_output = '{0}/cmip6_correct/{1}/{2}/{3}_br_day_{1}_{2}_{4}_{5}_correct.nc'.format(dataset_dir, cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], time[idx])
	write_3d_nc(nc_output, corrected_dataset, time_array, lat_array, lon_array, var_units, var_shortname, var_longname, time_units, missing_value=-999.)


