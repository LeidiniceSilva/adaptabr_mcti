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
import scipy.stats as st
import geopandas as gpd
import matplotlib.pyplot as plt
import shapefile as shp

from netCDF4 import Dataset
from shapely.geometry import Point, Polygon

warnings.filterwarnings('ignore')

# Best models list
mdl = 17
best_models = [7, 9, 13, 15, 17]
cmip6 = {7 :['GFDL-ESM4', 'r1i1p1f1_gr1'], 9 :['INM-CM5-0', 'r1i1p1f1_gr1'], 13 :['MPI-ESM1-2-HR', 'r1i1p1f1_gn'], 15   :['MRI-ESM2-0', 'r1i1p1f1_gn'], 17 :['NorESM2-MM', 'r1i1p1f1_gn']}

# Variable dictionary
var = 1
var_dict = {1 :['pr', 'pr'], 2 :['Tmax', 'tasmax'], 3 :['Tmin', 'tasmin']}

dt = '19860101-20051231'
experiment = 'historical'
var_obs   = var_dict[var][0]
var_cmip6 = var_dict[var][1]

print()
print()
print('Model -->', cmip6[mdl][0])
print('Var   -->', var_cmip6)


def import_observed(var_name, target_date):

	""" Import observed data.
	:param: model_name: str (BR-DWGD)
	:param: target_date: datetime
	:return: observed data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	file_name = "./{0}_{1}_BR-DWGD_UFES_UTEXAS_v_3.2.2_lonlat.nc".format(var_name, target_date)
	data  = netCDF4.Dataset(file_name)
	var   = data.variables[var_name][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	data.close()

	return lat, lon, value


def import_simulated(model_name, exp_name, var_name, member, target_date):

	""" Import simulated data.
	:param: model_name: str (Nor, GFDL, MPI, INM and MRI)
	:param: target_date: datetime
	:return: simulated data (Prec, Tasmax and Tasmin)
	:return: flag that it refers to data control
	:rtype: 3D array
	"""

	file_model = "./{0}_br_day_{1}_{2}_{3}_{4}_lonlat.nc".format(var_name, model_name, exp_name, member, target_date)
	data  = netCDF4.Dataset(file_model)
	var   = data.variables[var_name][:]
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	data.close()
  
	return lat, lon, value


def mask_function_super_fast(data_, lat_, lon_, ind_):

	mask = np.full(data_[0,:,:].shape, np.nan)
	
	lats = []
	lons = []
	for i in range(0, mask.shape[0]):
		for j in range(0, mask.shape[1]):
			lats.append(lat_[i])
			lons.append(lon_[j])
			
	# Open the shapefile
	counties = gpd.GeoDataFrame.from_file('./BRA_adm0.shp')
	
	# Create geodataframe from numpy arrays
	df = pd.DataFrame({'lon':lons, 'lat':lats})
	df['coords'] = list(zip(df['lon'], df['lat']))
	df['coords'] = df['coords'].apply(Point)
	points = gpd.GeoDataFrame(df, geometry='coords', crs=counties.crs)

	# Perform spatial join to match points and polygons
	pointInPolys = gpd.tools.sjoin(points, counties, predicate="within", how='left')
	
	# Example use: get points in Los Angeles, CA.
	pnt_BR = points[pointInPolys.ISO=='BRA']

	if False: # Plot map with points in BRA in red

		base = counties.boundary.plot(linewidth=1, edgecolor="black")
		points.plot(ax=base, linewidth=1, color="blue", markersize=1)
		pnt_BR.plot(ax=base, linewidth=1, color="red", markersize=8)
		plt.show()
	
	# Applying the found points to the mask
		
	lon_point_inside = pnt_BR.iloc[:, 0]
	lat_point_inside = pnt_BR.iloc[:, 1]
	
	for la,lo in zip(lat_point_inside, lon_point_inside):

		la_index = np.where(la == lat_)[0][0]
		lo_index = np.where(lo == lon_)[0][0]
		
		mask[la_index,lo_index] = 1
	
	if False: # Plot to check if mask works
		
		for i in range(0, mask.shape[0],5):
			for j in range(0, mask.shape[1],5):
				if mask[i,j]==1:
					plt.scatter(lon_[j], lat_[i], marker='o', color='black')
				else:
					plt.scatter(lon_[j], lat_[i], marker='o', color='red')
		
		plt.show()
		plt.close()

	# Apply the mask to the data
	data_ = (mask*data_[ind_,:,:])
	print(data_.shape)
	print(data_[0,0::30,0::30])

	del lats
	del lons
	del counties
	del points
	del pointInPolys
	del pnt_BR
	del lon_point_inside
	del lat_point_inside
	del mask

	return data_.filled(np.nan)
	
	
def mask_function_slow(data_, lat_, lon_, ind_):

	mask = np.full(data_[0,:,:].shape, np.nan)
		
	# Shapefile of Brazil
	wc = gpd.read_file('./BRA_adm0.shp')
		
	if False: # Testing shape points

		sf = shp.Reader('./BRA_adm0.shp')

		plt.figure()
		for shape in sf.shapeRecords():
			x = [i[0] for i in shape.shape.points[:]]
			y = [i[1] for i in shape.shape.points[:]]
			plt.plot(x,y)
			
		for i in range(0, mask.shape[0], 20):
			for j in range(0, mask.shape[1], 20):
				p = Point(lon_[j], lat_[i])	
				print(p.within(wc.geometry[0]))
				if p.within(wc.geometry[0]):
					plt.scatter(lon_[j], lat_[i], marker='o', color='black')
					
		plt.show()
		plt.close()
		exit()

	# Mask to keep only lat and lon points  Brazil
	for i in range(0, mask.shape[0], 1):
		for j in range(0, mask.shape[1], 1):
			print('point -->', i,j)
			
			p = Point(lon_[j], lat_[i])
			
			# p = Point(round(lon_[j],2), round(lat_[i],2))
			check_inside = p.within(wc.geometry[0])
			
			if check_inside:
				mask[i, j] = 1
				# print('inside')
				
			if i==40  and j==40:  print(i,j)
			if i==80  and j==80:  print(i,j)
			if i==120 and j==120: print(i,j)
			if i==150 and j==150: print(i,j)
			
	# Apply the mask to the data
	data_ = (mask*data_[ind_,:,:])
	print(data_.shape)
	print(data_[0,0::30,0::30])
	
	return data_.filled(np.nan)
	
	
def mask_function(data_, lat_, lon_, ind_):

	# Shapefile of Brazil
	wc = gpd.read_file('/content/drive/My Drive/Shapefile/BRA_adm0.shp')

	# Mask to keep only lat and lon points  Brazil
	mask = np.zeros(np.shape(data_[0,:,:]))*np.nan
	for i in range(np.shape(data_[0,:,:])[0]):
		for j in range(np.shape(data_[0,:,:])[1]):
			la = lat_[i,j]
			lo = lon_[i,j]
			p = Point(la,lo)
			if p.within(wc.geometry[0]):
				mask[i,j] = 1

	# Apply the mask
	data_ = (mask*data_[ind_,:,:]).flatten()
	data_ = data_[~np.isnan(data_)]

	return data_


def quantile_mapping(dataset, fit_params_hist, th_distrib_hist, fit_params_obs, th_distrib_obs):

	k,  s,  loc,  scale  = fit_params_hist
	k_, s_, loc_, scale_ = fit_params_obs

	# Apply cumulative distribution function with historical params and distribution
	cdf_hist = th_distrib_hist.cdf(dataset, k=k, s=s, loc=loc, scale=scale)

	# Then apply inverse cdf (or ppf) with era5 fitted params and distribution
	correct_data = th_distrib_obs.ppf(cdf_hist, k=k_, s=s_, loc=loc_, scale=scale_)

	del cdf_hist
	
	return correct_data


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
		calendar = 'standard'
	elif cmip6[mdl][0] == 'INM-CM5-0':
		cmip6_inst = 'INM, Russian Academy of Science, Moscow 119991, Russia'
		calendar = '365_day'
	elif cmip6[mdl][0] == 'NorESM2-MM':
		cmip6_inst = 'NCC, c/o MET-Norway, Henrik Mohns plass 1, Oslo 0313, Norway'
		calendar = '365_day'
	elif cmip6[mdl][0] == 'MPI-ESM1-2-HR':
		cmip6_inst = 'MPI for Meteorology, Hamburg 20146, Germany'
		calendar = 'standard'
	elif cmip6[mdl][0] == 'MRI-ESM2-0':
		cmip6_inst = 'MRI, Tsukuba, Ibaraki 305-0052, Japan'
		calendar = 'standard'
	else:
		cmip6_inst = 'NCC, c/o MET-Norway, Henrik Mohns plass 1, Oslo 0313, Norway'
		calendar = 'standard'

	foo = Dataset(ncname, 'w', format='NETCDF4_CLASSIC')

	foo.Conventions = 'CF-1.7 CMIP-6.0 UGRID-1.0'
	foo.title = '{0} correct bias model output'.format(cmip6[mdl][0])
	foo.institution = '{0}'.format(cmip6_inst)
	foo.source = '{0}'.format(cmip6[mdl][0])
	foo.history = 'CMIP6 model with corrected bias'
	foo.references = 'https://esgf.ceda.ac.uk/thredds/catalog/esg_cmip6/catalog.html'
	foo.comment = 'This netCDF4 file was corrected using EQM function'

	foo.createDimension('time', None)
	foo.createDimension('lat', lat_array.shape[0])
	foo.createDimension('lon', lon_array.shape[0])

	times = foo.createVariable('time', float, ('time'))
	times.units = time_units
	times.calendar = '{0}'.format(calendar)
	times[:] = range(len(time_array))

	laty = foo.createVariable('lat', 'f4', 'lat')
	laty.units = 'degrees_north'
	laty.long_name = 'latitude'
	laty[:] = lat_array

	lonx = foo.createVariable('lon', 'f4', 'lon')
	lonx.units = 'degrees_east'
	lonx.long_name = 'longitude'
	lonx[:] = lon_array

	var = foo.createVariable(var_shortname, float, ('time', 'lat', 'lon'), zlib=True)
	var.units = var_units
	var.long_name = var_longname
	var.missing_value = missing_value
	var[:] = var_array

	foo.close()


print()
print()
print('# Import cmip models and obs database')

print()
lat_array, lon_array, obs = import_observed(var_obs, dt)
print('OBS -->', lat_array.shape, lon_array.shape, obs.shape)

print()
lat_array, lon_array, sim = import_simulated(cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt)
print('MOD -->', lat_array.shape, lon_array.shape, sim.shape)


print()
print()
print('# Date range of the historical period')

time_obs = xr.cftime_range(start=cftime.DatetimeGregorian(1986,1,1),end=cftime.DatetimeGregorian(2005,12,31),freq="D",calendar='standard')
ind_obs = range(len(time_obs))

# Select the appropraite time index
if np.shape(sim)[0]==7305:
	time_array = xr.cftime_range(start=cftime.DatetimeGregorian(1986,1,1),end=cftime.DatetimeGregorian(2005,12,31),freq="D",calendar='standard')
	ind_366 = range(len(time_array))
	ind_sim = ind_366
else:
	time_array = xr.cftime_range(start=cftime.DatetimeNoLeap(1986,1,1),end=cftime.DatetimeNoLeap(2005,12,31),freq="D",calendar='365_day')
	ind_365 = range(len(time_array))
	ind_sim = ind_365

print('ind_sim -->', ind_sim)


print()
print()
print('# Apply mask function')

# obs = mask_function         (obs, lat_array, lon_array, ind_obs) # original - not working
# obs = mask_function_slow    (obs, lat_array, lon_array, ind_obs) # fixed - work but super slow

obs_masked = mask_function_super_fast(obs, lat_array, lon_array, ind_obs) # new - fast
print()
print('obs masked -->', obs_masked.shape)

sim_masked = mask_function_super_fast(sim, lat_array, lon_array, ind_sim)
print()
print('sim masked -->', sim_masked.shape)

del obs
del sim

print()
print()
print('# Calculate st.mielke parameters')

print()
print('Fitting data')
print()

sim_corrected = sim_masked*np.nan

for i in range(0, obs_masked.shape[1], 1):
	
	for j in range(0, obs_masked.shape[2], 1):

		for k in range(0, obs_masked.shape[1], 10):
			if i==k and j==1: print('running --> i={} j={}'.format(i, j))
		
		if np.isnan(obs_masked[0,i,j]): continue
						
		obs_1d_params = st.mielke.fit(obs_masked[:,i,j]) # obs_k, obs_s, obs_loc, obs_scale
		sim_1d_params = st.mielke.fit(sim_masked[:,i,j]) # sim_k, sim_s, sim_loc, sim_scale

		sim_corrected[:,i,j] = quantile_mapping(sim_masked[:,i,j], sim_1d_params, st.mielke, obs_1d_params, st.mielke)
		
		"""
		if False: 
			plt.figure(figsize=(16, 8))
			plt.plot(ind_obs[0:360], obs_masked[0:360,i,j],    '-', color='green', lw=2, label='obs')
			plt.plot(ind_sim[0:360], sim_masked[0:360,i,j],    '-', color='blue',  lw=2, label='model')
			plt.plot(ind_sim[0:360], sim_corrected[0:360,i,j], '-', color='red',   lw=2, label='model corrected')
			plt.legend()
			# plt.show()
			plt.savefig('sim_1d_corrected.png')
			plt.close()
		"""
				
if var_cmip6 == 'pr':
	
	# To replace negative values with 0
	sim_corrected = np.where(sim_corrected<0, 0, sim_corrected)
	
	# Fix points with wrong high values to the maximum in obs
	for i in range(0, obs_masked.shape[1], 1):
		for j in range(0, obs_masked.shape[2], 1):
			sim_corrected[:,i,j] = np.where(sim_corrected[:,i,j]>np.nanmax(obs_masked[:,i,j]), np.nanmax(obs_masked[:,i,j]), obs_masked[:,i,j])

# Path out to save netCDF4
if  var_cmip6 == 'pr':
	var_units     = 'mm'
	var_shortname = 'pr'
	var_longname  = 'Daily Total Precipitation'
	time_units    = 'days since {}'.format(time_array[0])
elif var_cmip6 == 'tasmax':
	var_units     = 'Celcius degrees'
	var_shortname = 'tasmax'
	var_longname  = 'Daily Maximum Near-Surface Air Temperature'
	time_units    = 'days since {}'.format(time_array[0])
else:
	var_units     = 'Celcius degrees'
	var_shortname = 'tasmin'
	var_longname  = 'Daily Minimum Near-Surface Air Temperature'
	time_units    = 'days since {}'.format(time_array[0])

nc_output = './{0}_br_day_{1}_{2}_{3}_{4}_correct.nc'.format(cmip6[mdl][0], experiment, var_cmip6, cmip6[mdl][1], dt)
write_3d_nc(nc_output, sim_corrected, time_array, lat_array, lon_array, var_units, var_shortname, var_longname, time_units, missing_value=np.nan)
if os.path.exists(nc_output):
	print('done -->', nc_output)

# time execution ---> 3h 36min
# real	216m58.940s
# user	216m50.879s
# sys	0m3.329s



