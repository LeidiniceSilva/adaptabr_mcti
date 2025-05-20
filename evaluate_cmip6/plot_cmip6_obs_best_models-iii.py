# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot taylor diagram of cmip6 models"

import os
import netCDF4
import numpy as np
import numpy.ma as ma
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import mpl_toolkits.axisartist as axisartist

from dict_cmip6_models_name import cmip6
from comp_taylor_diagram import TaylorDiagram


def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/project/database/obs'
	arq   = '{0}/{1}_{2}_BR-DWGD_UFES_UTEXAS_v_3.0_ANN_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]	
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
			
	return ts_latlon

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/project/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_ANN_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	fld_mean = np.nanmean(value, axis=0)
	
	latlon = []
	for i in range(0, fld_mean.shape[0]):
		for ii in fld_mean[i]:
			latlon.append(ii)
	ts_latlon = np.array(latlon)
			
	return ts_latlon
	              
               
# Import cmip models and obs database 
var_obs = 'pr'
var_cmip6 = 'pr'
dt = '1986-2005'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)
clim_br_obs = import_obs(var_obs, 'BR', dt)

std_namz_obs = np.nanstd(clim_namz_obs, ddof=0)
std_samz_obs = np.nanstd(clim_samz_obs, ddof=0)
std_neb_obs = np.nanstd(clim_neb_obs, ddof=0)
std_sam_obs = np.nanstd(clim_sam_obs, ddof=0)
std_lpb_obs = np.nanstd(clim_lpb_obs, ddof=0)
std_br_obs = np.nanstd(clim_br_obs, ddof=0)

std_namz_cmip6 = []
std_samz_cmip6 = []
std_neb_cmip6 = []
std_sam_cmip6 = []
std_lpb_cmip6 = []
std_br_cmip6 = []

pcc_namz_cmip6 = []
pcc_samz_cmip6 = []
pcc_neb_cmip6 = []
pcc_sam_cmip6 = []
pcc_lpb_cmip6 = []
pcc_br_cmip6 = []

legend = []

best_models = [17, 7, 13, 9, 15]
for i in best_models:

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	pcc_namz_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_namz_obs), ma.masked_invalid(clim_namz_cmip))[0][1])
	std_namz_cmip6.append(np.nanstd(clim_namz_cmip, ddof=0))

	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	pcc_samz_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_samz_obs), ma.masked_invalid(clim_samz_cmip))[0][1])
	std_samz_cmip6.append(np.nanstd(clim_samz_cmip, ddof=0))
	
	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	pcc_neb_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_neb_obs), ma.masked_invalid(clim_neb_cmip))[0][1])
	std_neb_cmip6.append(np.nanstd(clim_neb_cmip, ddof=0))
		
	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	pcc_sam_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_sam_obs), ma.masked_invalid(clim_sam_cmip))[0][1])
	std_sam_cmip6.append(np.nanstd(clim_sam_cmip, ddof=0))
		
	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	pcc_lpb_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_lpb_obs), ma.masked_invalid(clim_lpb_cmip))[0][1])
	std_lpb_cmip6.append(np.nanstd(clim_lpb_cmip, ddof=0))

	clim_br_cmip = import_cmip(var_cmip6, 'BR', cmip6[i][0], cmip6[i][1], dt)
	pcc_br_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_br_obs), ma.masked_invalid(clim_br_cmip))[0][1])
	std_br_cmip6.append(np.nanstd(clim_br_cmip, ddof=0))
	
	legend.append(cmip6[i][0])

pcc_namz = np.array(pcc_namz_cmip6)
std_namz = np.array(std_namz_cmip6)

pcc_samz = np.array(pcc_samz_cmip6)
std_samz = np.array(std_samz_cmip6)

pcc_neb = np.array(pcc_neb_cmip6)
std_neb = np.array(std_neb_cmip6)

pcc_sam = np.array(pcc_sam_cmip6)
std_sam = np.array(std_sam_cmip6)

pcc_lpb = np.array(pcc_lpb_cmip6)
std_lpb = np.array(std_lpb_cmip6)

pcc_br = np.array(pcc_br_cmip6)
std_br = np.array(std_br_cmip6)

# Plot cmip models and obs database 
fig = plt.figure(figsize=(12,8))

titleprops_dict = dict(loc='left', fontweight='bold', x=0.0, y=1.05)

fig, ax1 = TaylorDiagram(std_namz, pcc_namz, std_namz_obs, fig=fig, rect=231, title='(a) NAMZ', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax2 = TaylorDiagram(std_samz, pcc_samz, std_samz_obs, fig=fig, rect=232, title='(b) SAMZ', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax3 = TaylorDiagram(std_neb, pcc_neb, std_neb_obs, fig=fig, rect=233, title='(c) NEB', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax4 = TaylorDiagram(std_sam, pcc_sam, std_sam_obs, fig=fig, rect=234, title='(d) SAM', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax5 = TaylorDiagram(std_lpb, pcc_lpb, std_lpb_obs, fig=fig, rect=235, title='(e) LPB', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax6 = TaylorDiagram(std_br, pcc_br, std_br_obs, fig=fig, rect=236, title='(e) BR', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
						
ax6.legend(bbox_to_anchor=(-2.5, -0.165), loc='upper left', ncol=6)
plt.subplots_adjust(hspace=0.9)
plt.subplots_adjust(bottom=0.15)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/project/figs/figs_report-II'
name_out = 'pyplt_taylor_diagram_best_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()

