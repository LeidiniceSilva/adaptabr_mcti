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

from taylor_diagram import TaylorDiagram
from dict_cmip6_models_name import cmip6
from comp_statistical_metrics import compute_tss


def import_obs(param, area, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/obs'
	arq   = '{0}/{1}_{2}_CRU_ts4_ANN_{3}_lonlat.nc'.format(path, param, area, date)	
		
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]	
	obs = np.nanmean(value, axis=0)
	
	latlon_obs = []
	for i in range(0, obs.shape[0]):
		for ii in obs[i]:
			latlon_obs.append(ii)
	latlon_obs = np.array(latlon_obs)
			
	return lat, lon, latlon_obs

	
def import_cmip(param, area, model, exp, date):
	
	path  = '/home/nice/Documentos/AdaptaBrasil_MCTI/database/cmip6'
	arq   = '{0}/{1}_{2}_{3}_historical_{4}_ANN_{5}_lonlat.nc'.format(path, param, area, model, exp, date)	
				
	data  = netCDF4.Dataset(arq)
	var   = data.variables[param][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	value = var[:][:,:,:]
	mdl = np.nanmean(value, axis=0)

	latlon_mdl = []
	for i in range(0, mdl.shape[0]):
		for ii in mdl[i]:
			latlon_mdl.append(ii)
	latlon_mdl = np.array(latlon_mdl)
	
	return lat, lon, latlon_mdl
	              
               
# Import cmip models and obs database 
var_obs = 'pre'
var_cmip6 = 'pr'
dt = '1980-2014'

clim_namz_obs = import_obs(var_obs, 'NAMZ', dt)
clim_samz_obs = import_obs(var_obs, 'SAMZ', dt)
clim_neb_obs  = import_obs(var_obs, 'NEB', dt)
clim_sam_obs = import_obs(var_obs, 'SAM', dt)
clim_lpb_obs = import_obs(var_obs, 'LPB', dt)

std_namz_obs = np.nanstd(clim_namz_obs[2], ddof=0)
std_samz_obs = np.nanstd(clim_samz_obs[2], ddof=0)
std_neb_obs = np.nanstd(clim_neb_obs[2], ddof=0)
std_sam_obs = np.nanstd(clim_sam_obs[2], ddof=0)
std_lpb_obs = np.nanstd(clim_lpb_obs[2], ddof=0)
	
std_namz_cmip6 = []
std_samz_cmip6 = []
std_neb_cmip6 = []
std_sam_cmip6 = []
std_lpb_cmip6 = []

pcc_namz_cmip6 = []
pcc_samz_cmip6 = []
pcc_neb_cmip6 = []
pcc_sam_cmip6 = []
pcc_lpb_cmip6 = []

legend = []

for i in range(1, 19):

	clim_namz_cmip = import_cmip(var_cmip6, 'NAMZ', cmip6[i][0], cmip6[i][1], dt)
	pcc_namz_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_namz_obs[2]), ma.masked_invalid(clim_namz_cmip[2]))[0][1])
	std_namz_cmip6.append(np.nanstd(clim_namz_cmip[2], ddof=0))

	clim_samz_cmip = import_cmip(var_cmip6, 'SAMZ', cmip6[i][0], cmip6[i][1], dt)
	pcc_samz_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_samz_obs[2]), ma.masked_invalid(clim_samz_cmip[2]))[0][1])
	std_samz_cmip6.append(np.nanstd(clim_samz_cmip[2], ddof=0))
	
	clim_neb_cmip = import_cmip(var_cmip6, 'NEB', cmip6[i][0], cmip6[i][1], dt)
	pcc_neb_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_neb_obs[2]), ma.masked_invalid(clim_neb_cmip[2]))[0][1])
	std_neb_cmip6.append(np.nanstd(clim_neb_cmip[2], ddof=0))
		
	clim_sam_cmip = import_cmip(var_cmip6, 'SAM', cmip6[i][0], cmip6[i][1], dt)
	pcc_sam_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_sam_obs[2]), ma.masked_invalid(clim_sam_cmip[2]))[0][1])
	std_sam_cmip6.append(np.nanstd(clim_sam_cmip[2], ddof=0))
		
	clim_lpb_cmip = import_cmip(var_cmip6, 'LPB', cmip6[i][0], cmip6[i][1], dt)
	pcc_lpb_cmip6.append(ma.corrcoef(ma.masked_invalid(clim_lpb_obs[2]), ma.masked_invalid(clim_lpb_cmip[2]))[0][1])
	std_lpb_cmip6.append(np.nanstd(clim_lpb_cmip[2], ddof=0))
	
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
	
# Plot cmip models and obs database 
fig = plt.figure(figsize=(12,8))
fig.tight_layout(h_pad=1)

titleprops_dict = dict(loc='left', fontweight='bold', x=0.0, y=1.05)

fig, ax1 = TaylorDiagram(std_namz, pcc_namz, std_namz_obs, fig=fig, rect=231, title='(a) NAMZ', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax2 = TaylorDiagram(std_samz, pcc_namz, std_samz_obs, fig=fig, rect=232, title='(b) SAMZ', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax3 = TaylorDiagram(std_neb, pcc_neb, std_neb_obs, fig=fig, rect=233, title='(c) NEB', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax4 = TaylorDiagram(std_sam, pcc_sam, std_sam_obs, fig=fig, rect=234, title='(d) SAM', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
fig, ax5 = TaylorDiagram(std_lpb, pcc_lpb, std_lpb_obs, fig=fig, rect=235, title='(e) LPB', titleprops_dict=titleprops_dict, normalize=True, labels=legend, ref_label='Reference')
						
ax5.legend(bbox_to_anchor=(1.1, 0.1), loc='lower left', ncol=2)
plt.subplots_adjust(top=0.95)

# Path out to save figure
path_out = '/home/nice/Documentos/AdaptaBrasil_MCTI/figs/figs_report-II'
name_out = 'pyplt_talyor_diagram_cmip6_{0}_{1}.png'.format(var_cmip6, dt)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300)
plt.show()
exit()

