# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jul 05, 2023"
__description__ = "This script adjust bias of cmip6 models"

import os
import numpy as np

from netCDF4 import Dataset
from datetime import datetime


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



import os
import glob
import argparse
import bisect
import calendar
import numpy as np
import numpy.ma as ma
import scipy.stats as ss

from netCDF4 import Dataset
from collections import Counter
from dateutil.relativedelta import *
from datetime import datetime, date
from hidropy.utils.hidropy_utils import date2index, basin_dict
from hidropy.utils.write_flow import write_flow

from os.path import expanduser

HIDROPY_DIR = os.environ['HIDROPY_DIR']
    
__author__ = 'Leidinice Silva'
__email__ = 'leidinice.silva@funceme.br'
__date__ = '11/08/2017'
__description__='This script remove daily mean flow bias from GFS05' \
                ' model. This script focus on operational' \
                ' period and creates files in the same flow pattern.'


def arguments():
    global args

    parser = argparse.ArgumentParser(description=__description__,
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
		
    parser.add_argument('--day_target', default='',
                        help='Insert initial forecast day.')
    parser.add_argument('--month_target', default='', 
                        help='Insert initial forecast month.')
    parser.add_argument('--year_target', default='', 
                        help='Insert initial forecast year.')
    parser.add_argument('--local_dir', default='/io',
                        help='Path to data directory.')
    parser.add_argument('--model_name', choices=['gfs05'],
                        type=str, nargs='?', help='Model option.')
    parser.add_argument('--hydro_model', choices=['hymod', 'smap', 'all'],
                        type=str, nargs='?', help='Hidrologicals models options.')
    parser.add_argument('--basin', help='Name of input basin.', default='')
    parser.add_argument('--macro', action='store_true',
                        help='Compute flow for all macro-basins.')
    parser.add_argument('--micro', action='store_true',
                        help='Compute flow for all micro-basins of '
                             'input basin.')
    parser.add_argument('--all_basins', action='store_true',
                        help='True to compute flow for all macro-basins'
                             ' and micro-basins.')
    args = parser.parse_args()

    return args


def define_initial_parameters(hydromodel):
    """ Define target data and models list.

    :param: modelname: str
    :return: target_date: datatime
    :return: model name list: list
    """

    # List Models
    if hydromodel  == 'all':
        list_model = ['hymod', 'smap']
    else:
        list_model = [hydromodel]

    if not day:
        target_date = datetime(date.today().year, date.today().month,
                               date.today().day)
    else:
        target_date = datetime(int(year), int(month), int(day))

    return target_date, list_model    
   
   
def define_dates(target_date):
    """ Define date start run, and dates that the forecast starts and ends.

    :param: target_date: datatime
    :return: start_rundate: str
    :return: start_fcstdate: str
    :return: end_fcstdate: str
    """    
    
    str_mon = target_date.strftime("%b").lower()

    target_ifcstdate = target_date + relativedelta(days=1)
    target_efcstdate = target_date + relativedelta(days=7)
    
    start_rundate = '{0}{1:02d}{2:02d}'.format(target_date.year,
                                          target_date.month, target_date.day)
    star_fcstdate = '{0}{1:02d}{2:02d}'.format(target_ifcstdate.year,
                                        target_ifcstdate.month, target_ifcstdate.day)
    end_fcstdate = '{0}{1:02d}{2:02d}'.format(target_efcstdate.year,
                                       target_efcstdate.month, target_efcstdate.day)

    return str_mon, start_rundate, star_fcstdate, end_fcstdate
    

def import_fcst_model_data(model_name, hydro_model, target_date, basin):
    """ Import fcst model data.

    :param: model_name: str (gfs05)
    :param: target_date: datetime
    :param: basin: str
    :return: flow time serie
    :return: flag that it refers to data control
    :rtype:list
    """

    str_mon, start_rundate, star_fcstdate, end_fcstdate = define_dates(target_date)
    basin_full_name = basin_dict(basin)[2]
    basin_name = basin_dict(basin)[1]
    
    dir_modelname   =  "{0}{1}/{2}/{3}_{4}/{5}/{6}/{7}/". \
		       format(HIDROPY_DIR, localdir, var_name, hydro_model, time_freq,
		       model_name, str_mon, basin_name)
		       
    file_modelname  = "{0}_{1}_{2}_fcst_{3}_{4}_{5}_{6}_{7}.nc". \
		       format(var_name, time_freq, model_name, start_rundate, star_fcstdate, end_fcstdate,
		       hydro_model, basin_full_name)
    try:
	input_data = Dataset(os.path.join(dir_modelname, file_modelname))
	var = input_data.variables[var_name][:]	
	input_data.close()
	flag = True
    except:
	print 'File: {0} is not available'.format(os.path.join(dir_modelname, file_modelname))
	print
	var = np.nan
	flag = False

    return var, flag


def import_hind_model_data(model_name, hydro_model, target_date, basin):
    """ Import fcst model data.

        :param: model_name: str ('gfs05')
        :param: target_date: datetime
        :param: basin: str
        :return: flow time serie
        :rtype: list
        """
	
    clim_list = []
    date_list = []
    ihind = 2009
    
    for y in range(ihind, 2012 + 1):
	ndays = calendar.monthrange(y, target_date.month)[1]
        for d in range(1, ndays + 1):
	
	    target_rundate = datetime(y, target_date.month, d)
	    str_mon, start_rundate, star_fcstdate, end_fcstdate = define_dates(target_rundate)
		    
	    basin_full_name = basin_dict(basin)[2]
	    basin_name = basin_dict(basin)[1]

	    dir_modelname   =  "{0}{1}/{2}/{3}_{4}/{5}/{6}/{7}/". \
			       format(HIDROPY_DIR, localdir, var_name, hydro_model, time_freq,
			       model_name, str_mon, basin_name)
			    
	    file_modelname  = "{0}_{1}_{2}_fcst_{3}_{4}_{5}_{6}_{7}.nc". \
			       format(var_name, time_freq, model_name, start_rundate, star_fcstdate, end_fcstdate,
			       hydro_model, basin_full_name)
			
	    try:
		input_data = Dataset(os.path.join(dir_modelname, file_modelname))
		var = input_data.variables[var_name][:]
		
		if ma.count_masked(var) == 0:
		    var = np.nansum(var)
		else:
		    var = np.nan
		    date_list.append(target_rundate)
		input_data.close()	
		    
	    except:
		print 'Hindcast File: {0} is not available'.format(os.path.join(dir_modelname, file_modelname))
		print
		var = np.nan
		date_list.append(target_rundate)
	    clim_list.append(var)

    clim = np.array(clim_list)
    clim = clim[~np.isnan(clim)]
    
    return clim, date_list


def import_hind_obs_data(model_name, target_date, basin, dates_list):
    """ Import obs data.

    :param: model_name: str (gfs05')
    :param: target_date: datetime
    :param: basin: str
    :return: flow time serie
    """
    
    obs = []
    ihind = 2009
    
    for y in range(ihind, 2012 + 1):
	ndays = calendar.monthrange(y, target_date.month)[1]
        for d in range(1, ndays + 1):
	    target_rundate = datetime(y, target_date.month, d)
	    
	    if not target_rundate in dates_list:
		str_mon, start_rundate, star_fcstdate, end_fcstdate = define_dates(target_rundate)
		basin_full_name = basin_dict(basin)[2]
		basin_name = basin_dict(basin)[1]
		
		file_name = '{0}{1}/{2}/ons_daily/1981-present/{3}/' \
			    '{2}_{4}_{5}_{6}_{7}.nc' \
			    .format(HIDROPY_DIR, localdir, var_name, basin_name, time_freq,
			     data_base, hind_obs_period, basin_full_name)
			     
		inc = Dataset(file_name)
		time_var = inc.variables['time']
		
		ifcstdate = datetime.strptime(star_fcstdate, "%Y%m%d")
		ffcstdate = datetime.strptime(end_fcstdate, "%Y%m%d")

		iidx = date2index(datetime(ifcstdate.year, ifcstdate.month, ifcstdate.day), time_var)
		fidx = date2index(datetime(ffcstdate.year, ffcstdate.month, ffcstdate.day), time_var)
		aux = inc.variables[var_name][iidx:fidx+1]
		obs.append(np.nansum(aux))
	    else:
                print "Model data for {0} does not exist.".format(target_rundate)    
		print
                continue

    inc.close()
    obs = np.array(obs)
    
    return obs


def closest(fcst_target, clim_list):
    aux = []    
    for ii in range(0, clim_list.shape[0]):
        aux.append(abs(fcst_target-clim_list[ii]))
    min_value = np.min(aux)    
    idx_min = np.where(aux == min_value)[-1][-1]

    return idx_min


def eqmres_bias_correction(hind, obs, fcst):    
    """ Remove precipitaiton bias via Empirical Quantile Mapping
    (eQM).

    :param: hind: list
    :param: obs: list
    :param: fcst: list
    :return: precipitation
    """

    xs_hind = np.sort(hind)
    #ys_hind = np.arange(1, len(xs_hind)+1)/float(len(xs_hind))
    xs_obs = np.sort(obs)
    #ys_obs = np.arange(1, len(xs_obs)+1)/float(len(xs_obs))
    
    idx_hind = closest(fcst, xs_hind)
    corr = xs_obs[idx_hind]
    
    return np.array(corr)


if __name__ == "__main__":

    arguments()
    global day, month, year, localdir, var_name, time_freq, data_base,  \
           hind_obs_period
	   
    # Target basin
    bname  = args.basin
    macros = args.macro
    micros = args.micro
    allbasins = args.all_basins

    # Local directory,variable name and Model
    localdir   = args.local_dir
    modelname  = args.model_name
    hydromodel = args.hydro_model
    
    type_corr  = '_cor_flow' 
    var_name   = 'flow'
    time_freq  = 'daily'
    data_base  = 'ons_obs'
    hind_obs_period = '19810101_20151231'

    # Target period
    day = args.day_target
    month = args.month_target
    year = args.year_target

    # Target basin
    if macros or micros or allbasins:
        basin_list = basin_dict(bname, macro=macros, micro=micros,
                                all_basins=allbasins)
    else:
        basin_list = [bname]
                
    basinlist = []
    for l in basin_list:
        if l == 'amazonas_guapore_inc' or l == 'paraguai_ponte_de_pedra'  \
	    or l == 'parana_ilha_solteira_equivalente' or l == 'sao_francisco_xingo_inc' \
	    or l == 'tiete_edgard_de_souza_inc' or l == 'tiete_traicao_inc':
            var_aux0 = []
        else:
            basinlist.append(l)
        
    target_date, list_model = define_initial_parameters(hydromodel)
    
    for hydroname in list_model:
        print "Run Remotion Bias for Model and Date: {0} - {1}/{2:02d}".format(modelname.upper(), target_date.year, target_date.month)
        print 
         
        for basin in basinlist:
            print 'Processing Basin: {0}'.format(basin)
	    print

            basin_full_name = basin_dict(basin)[2]
            basin_name = basin_dict(basin)[1]

	    str_mon, start_rundate, star_fcstdate, end_fcstdate = define_dates(target_date)
            
	    fcst, flag = import_fcst_model_data(modelname, hydroname, target_date, basin)
	    fcst_acc = np.array(np.nansum(fcst))
	    
	    if flag:
		
                model, dates_list = import_hind_model_data(modelname, hydroname, target_date, basin)
                hind = import_hind_obs_data(modelname, target_date, basin, dates_list)

                if  len(model) != len(hind):
                    print 'Error hind and obs climatology period1'
		    print
		
		# Removing Bias - Empirical Quantil Mapping Method
		fcst_corr = eqmres_bias_correction(model, hind, fcst_acc)

		corr = np.full((7), np.nan)
		for i in range(0, 7):
		    
		    if (float(fcst_acc) != 0.):
                        corr[i] = ((fcst[i] / float(fcst_acc)) * fcst_corr)
                    else:
                        corr[i] = fcst_corr / 7
                corr[np.isnan(corr)] = -999.    
		    			
		path_out =  "{0}{1}/{2}/{3}_{4}{5}/{6}/{7}/{8}/". \
			   format(HIDROPY_DIR, localdir, var_name, hydroname, time_freq,
				type_corr, modelname, str_mon, basin_name)

		if not os.path.exists(path_out):
		    os.makedirs(path_out)
			
		cor_name_nc = write_flow(np.array(corr), star_fcstdate, end_fcstdate, time_freq, 
			      var_name, modelname, 'fcst', '{0}{1}'.format(basin_full_name, type_corr), hydroname,
			      init_date=start_rundate, output_path=path_out)
	    else:
		continue


