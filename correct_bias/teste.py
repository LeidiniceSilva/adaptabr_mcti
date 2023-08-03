# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script correct bias of cmip6 models"

import numpy as np
import xarray as xr

from dict_cmip6_models_name import cmip6

# Dataset directory
dataset_dir = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias'

# Best models list
best_models = [7, 9, 13, 15, 17]
mdl = 17

# Variable dictionary
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}
var = 1

dt = '19860101-20051231'
exp = 'historical'
model_name = cmip6[mdl][0]
member = cmip6[mdl][1]
var_obs = var_dict[var][0]
var_cmip6 = var_dict[var][1]

print(model_name)
print(var_cmip6)


import xarray as xr
import numpy as np

def quantile_delta_mapping(observed, model_data):
    # Flatten the 3D arrays to 1D arrays
    observed_flat = observed.flatten()
    model_flat = model_data.flatten()

    # Sort both observed and model data
    observed_sorted = np.sort(observed_flat)
    model_sorted = np.sort(model_flat)

    # Compute the quantiles for observed and model data
    observed_quantiles = np.linspace(0, 1, len(observed_flat))
    model_quantiles = np.linspace(0, 1, len(model_flat))

    # Calculate the quantile delta mapping
    quantile_mapping = np.interp(model_quantiles, observed_quantiles, observed_sorted)

    # Apply the quantile delta mapping to the model data
    corrected_data_flat = np.interp(model_sorted, model_sorted, quantile_mapping)

    # Reshape the corrected 1D array back to 3D
    corrected_data = corrected_data_flat.reshape(model_data.shape)

    return corrected_data

if __name__ == '__main__':
    # Replace with the paths to your observed and model data NetCDF files
    model_data_file = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias/cmip6/NorESM2-MM/historical/pr_br_day_NorESM2-MM_historical_r1i1p1f1_gn_19860101-20051231_lonlat.nc'
    observed_data_file = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/correct_bias/obs/pr_19860101-20051231_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc'

    # Read historical observations and model data using xarray
    observed_data = xr.open_dataset(observed_data_file)['pr']
    model_data = xr.open_dataset(model_data_file)['pr']

    # Perform quantile delta mapping correction
    corrected_data = quantile_delta_mapping(observed_data.values, model_data.values)

    # Save the corrected data to a new NetCDF file
    corrected_data_array = xr.DataArray(corrected_data, coords=model_data.coords, dims=model_data.dims)
    corrected_data_array.to_netcdf('corrected_model_data.nc')


exit()

def empirical_quantile_mapping(model_data, observed_data):
    # Calculate the empirical quantiles of model and observed data
    model_quantiles = model_data.quantile(dim='time', q=np.linspace(0, 1, 101))
    observed_quantiles = observed_data.quantile(dim='time', q=np.linspace(0, 1, 101))

    # Calculate the adjustment factors
    adjustment_factors = observed_quantiles / model_quantiles

    # Adjust the model data using the quantile mapping adjustment factors
    adjusted_model_data = model_data * adjustment_factors

    return adjusted_model_data


# Example usage
if __name__ == "__main__":
    # Load CMIP6 model data and observed data
    # Replace 'model_file_path' and 'observed_file_path' with the actual file paths
    model_data = xr.open_dataset('{0}/cmip6/{1}/{2}/'.format(dataset_dir, model_name, exp) + '{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc'.format(var_cmip6, model_name, exp, member, dt))
    print(model_data)
     
    observed_data = xr.open_dataset('{0}/obs/'.format(dataset_dir) + '{0}_{1}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc'.format(var_obs, dt))
    print(observed_data)
    
    # Select the variable of interest (e.g., 'tas' for temperature)
    variable = 'pr'

    # Extract the specific variable for historical simulations
    model_historical_data = model_data[variable].sel(time=slice('1986', '1986'))
    observed_historical_data = observed_data[variable].sel(time=slice('1986', '1986'))

    # Apply empirical quantile mapping bias correction to historical simulations
    corrected_historical_data = empirical_quantile_mapping(model_historical_data[0:365,:,:], observed_historical_data[0:365,:,:])

    # Update the model data with the corrected historical simulations
    model_data[variable].loc[{'time': slice('1986', '1986')}] = corrected_historical_data

    # Save the corrected data to a new NetCDF file
    model_data.to_netcdf('corrected_historical_simulations.nc')
  
exit()


# ~ def quantile_mapping_bias_correction(model_data, observed_data):
    # ~ # Calculate the CDF of model and observed data
    # ~ model_cdf = model_data.rank(dim='time') / float(model_data.time.size)
    # ~ observed_cdf = observed_data.rank(dim='time') / float(observed_data.time.size)

    # ~ # Compute the quantile mapping correction function
    # ~ correction_function = xr.interp(model_cdf, observed_cdf.quantile(dim='time', q=model_cdf), observed_data)

    # ~ # Apply the correction to the model data
    # ~ corrected_model_data = xr.where(model_data.notnull(), correction_function, np.nan)

    # ~ return corrected_model_data
      

# ~ def bias_correction_bcsd(model_data, observed_data, variable):
    # ~ # Create the BCSD mapping
    # ~ mapping = xc.core.utils.create_bcsd_mapping(
        # ~ model_data[variable], observed_data[variable]
    # ~ )

    # ~ # Apply BCSD bias correction
    # ~ corrected_data = xc.indices.bias_correct(
        # ~ model_data[variable], mapping=mapping
    # ~ )

    # ~ return corrected_data
    
# ~ # Example usage
# ~ if __name__ == "__main__":
    # ~ # Load CMIP6 model data and observed data
    # ~ # Replace 'model_file_path' and 'observed_file_path' with the actual file paths
    # ~ model_data = xr.open_dataset('{0}/cmip6/{1}/{2}/'.format(dataset_dir, model_name, exp) + '{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc'.format(var_cmip6, model_name, exp, member, dt))
    # ~ print(model_data)
     
    # ~ observed_data = xr.open_dataset('{0}/obs/'.format(dataset_dir) + '{0}_{1}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc'.format(var_obs, dt))
    # ~ print(observed_data)
    
    # ~ # Select the variable of interest (e.g., 'tas' for temperature)
    # ~ variable = 'pr'

    # ~ # Apply BCSD bias correction
    # ~ corrected_data = bias_correction_bcsd(model_data, observed_data, variable)

    # ~ # Save the corrected data to a new NetCDF file
    # ~ corrected_data.to_netcdf('corrected_ssp1-26_data.nc')
    
    # ~ # Extract the specific variable for historical simulations
    # ~ model_historical_data = model_data[variable].sel(time=slice('1986', '2005'))
    # ~ print(model_historical_data)

    # ~ # Apply quantile mapping bias correction to historical simulations
    # ~ corrected_historical_data = quantile_mapping_bias_correction(model_historical_data, observed_data[variable])
    # ~ print(corrected_historical_data)
    
    # ~ # Update the model data with the corrected historical simulations
    # ~ model_data[variable].loc[{'time': slice('1986', '2005')}] = corrected_historical_data

    # ~ # Save the corrected data to a new NetCDF file
    # ~ model_data.to_netcdf('corrected_historical_simulations.nc')
    
