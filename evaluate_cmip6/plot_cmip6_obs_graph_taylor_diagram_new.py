# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Aug 03, 2023"
__description__ = "This script plot taylor diagram of cmip6 models"

import os
import netCDF4
import numpy as np
import xarray as xr
import scipy.stats as st
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpl_toolkits.axisartist as axisartist


from comp_taylor_diagram import TaylorDiagram


def import_obs():

	arq = xr.open_dataset('/home/nice/Downloads/CMIP6_historical_pr_monthlytimeseries/' + 'pr_sam_monthly_199501-201412_ts.nc')
	data = arq['pr']
	ts = data.sel(time=slice('1995-01-01','2014-12-31'))
	value = ts.values
			
	return value

	
def import_cmip(model):

	arq = xr.open_dataset('/home/nice/Downloads/CMIP6_historical_pr_monthlytimeseries/' + 'pr_SAM_monthly_{0}_historical.nc'.format(model))
	data = arq['precipitation_flux']
	ts = data.sel(time=slice('1995','2014'))
	value = ts.values
			
	return value
	              
               
# Import cmip models and obs database 
clim_sam_obs = import_obs()
clim_sam_obs = np.squeeze(clim_sam_obs)
std_sam_obs = np.nanstd(clim_sam_obs, ddof=0)

std_sam_cmip6 = []
pcc_sam_cmip6 = []

model_list=['ACCESS-CM2','ACCESS-ESM1-5','AWI-CM-1-1-MR','BCC-CSM2-MR','CanESM5','CAS-ESM2-0','CESM2','CESM2-WACCM','CMCC-CM2-SR5',
'CMCC-ESM2','CNRM-CM6-1','CNRM-ESM2-1','E3SM-1-0','E3SM-1-1-ECA','E3SM-1-1','EC-Earth3-CC','EC-Earth3','EC-Earth3-Veg-LR',
'FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-H','HadGEM3-GC31-LL','HadGEM3-GC31-MM','IITM-ESM','INM-CM4-8',
'INM-CM5-0','IPSL-CM6A-LR','KACE-1-0-G','KIOST-ESM','MIROC6','MIROC-ES2L','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0',
'NESM3','NorESM2-MM','TaiESM1','UKESM1-0-LL']

for mdl in model_list:
	print(mdl)

	clim_sam_cmip = import_cmip(mdl)
	pcc_sam_cmip6.append(st.pearsonr(clim_sam_obs, clim_sam_cmip)[0])
	std_sam_cmip6.append(np.nanstd(clim_sam_cmip, ddof=0))
		
pcc_sam = np.array(pcc_sam_cmip6)
std_sam = np.array(std_sam_cmip6)
	
# Plot cmip models and obs database 
fig = plt.figure(figsize=(8, 8))

titleprops_dict = dict(loc='left', fontweight='bold', x=0.0, y=1.05)
fig, ax1 = TaylorDiagram(std_sam, pcc_sam, std_sam_obs, fig=fig, rect=111, title='(a) SAM', titleprops_dict=titleprops_dict, normalize=True, labels=model_list, ref_label='Reference')				
ax1.legend(bbox_to_anchor=(0.5, -0.4), loc='lower center', fontsize=10, ncol=5)

# Path out to save figure
path_out = '/home/nice/Downloads/CMIP6_historical_pr_monthlytimeseries'
name_out = 'pyplt_taylor_diagram_cmip6_pr_color2.png'
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

