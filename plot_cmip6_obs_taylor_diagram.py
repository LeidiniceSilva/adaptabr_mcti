# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "This script plot taylor diagram of cmip6 models"

import os
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist

from dict_cmip6_models_name import cmip6
from comp_taylor_diagram import TaylorDiagram
from comp_statistical_metrics import compute_pdf


print('Read std and corr')
# Read std and corr
var = 'uv10'

if var == 'pr':
	stddev1 = np.array([0.8854898761088505, 3.0378468, 1.8355739])
	corrcoeff1 = np.array([0.7198974667207679, 0.9840942804159893, 0.9802872855751843])
	refstd1 = 1.966482847259447 
	models1 = ['RegCM_INMET', 'RegCM_CMORPH', 'RegCM_ERA5']
	stddev2 = np.array([1.205941042365372, 2.2004964, 1.0286905])
	corrcoeff2 = np.array([0.6536889355370896, 0.5776006345736175, 0.8013250619890022])
	refstd2 = 1.172333241817898 
	models2 = ['RegCM_INMET', 'RegCM_CMORPH', 'RegCM_ERA5']
	stddev3 = np.array([2.873461344286867, 4.299233, 2.4813993])
	corrcoeff3 = np.array([0.9821012147374468, 0.9879343127825087, 0.9732309003963605])
	refstd3 = 2.7121768885670376
	models3 = ['RegCM_INMET', 'RegCM_CMORPH', 'RegCM_ERA5']
	stddev4 = np.array([1.552037245070002, 2.6403365, 1.5159656])
	corrcoeff4 = np.array([0.8839608606676465, 0.8417541081890552, 0.9245633960200523])
	refstd4 = 2.468589819352789 
	models4 = ['RegCM_INMET', 'RegCM_CMORPH', 'RegCM_ERA5']
	stddev5 = np.array([1.9228065483116121, 2.9381738, 1.64957064])
	corrcoeff5 = np.array([0.9604711456049633, 0.9601399548715892, 0.9638035353974743])
	refstd5 = 1.9303826437215652 
	models5 = ['RegCM_INMET', 'RegCM_CMORPH', 'RegCM_ERA5']

elif var == 't2m':
	stddev1 = np.array([2.2670769428163795, 2.2407956])
	corrcoeff1 = np.array([0.9777592933835016, 0.9990502042603551])
	refstd1 = 2.5691664 
	models1 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev2 = np.array([3.5095901461400203, 3.454323])
	corrcoeff2 = np.array([0.9985541951895013, 0.9990973047345347])
	refstd2 = 3.636886 
	models2 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev3 = np.array([1.9470102170043189, 1.9555959])
	corrcoeff3 = np.array([0.5687493763620622, 0.5721556928007774])
	refstd3 = 2.049306 
	models3 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev4 = np.array([2.765495991163395, 2.8351734])
	corrcoeff4 = np.array([0.9979425890852557, 0.9989362089049363])
	refstd4 = 3.074289 
	models4 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev5 = np.array([2.3732859018788663, 2.4083235])
	corrcoeff5 = np.array([0.9920793474963145, 0.99458438974648])
	refstd5 = 2.58902 
	models5 = ['RegCM_INMET', 'RegCM_ERA5']
else:
	stddev1 = np.array([0.2651691558540669, 0.25749257])
	corrcoeff1 = np.array([0.842388983778291, 0.902911888009654])
	refstd1 = 0.2872287 
	models1 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev2 = np.array([0.1606290011605906, 0.18485455])
	corrcoeff2 = np.array([0.9523818433899746, 0.6216208707906876])
	refstd2 = 0.24910957 
	models2 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev3 = np.array([0.22757475654989132, 0.32147142])
	corrcoeff3 = np.array([0.30609175324341387, 0.6895755523968694])
	refstd3 = 0.42634553  
	models3 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev4 = np.array([0.1909248302231144, 0.17310059])
	corrcoeff4 = np.array([0.776442835272769, 0.816551580794502])
	refstd4 = 0.1540892 
	models4 = ['RegCM_INMET', 'RegCM_ERA5']
	stddev5 = np.array([0.21959988872502148, 0.26360545])
	corrcoeff5 = np.array([0.9419243360286759, 0.9373718906351077])
	refstd5 = 0.32911777 
	models5 = ['RegCM_INMET', 'RegCM_ERA5']
	
print('Plot figure')
# Plot figure
fig = plt.figure(figsize=(12,8))
fig.tight_layout(h_pad=1)

titleprops_dict = dict(loc='left', fontweight='bold',
					   x=0.0, y=1.05)

fig, ax1 = TaylorDiagram(stddev1, corrcoeff1, refstd1, 
						 fig=fig, rect=231, 
						 title='(a) Cluster I', titleprops_dict=titleprops_dict,
						 normalize=True, labels=models1, ref_label='Reference')

fig, ax2 = TaylorDiagram(stddev2, corrcoeff2, refstd2, 
						 fig=fig, rect=232, 
						 title='(b) Cluster II', titleprops_dict=titleprops_dict,
						 normalize=True, labels=models2, ref_label='Reference')

fig, ax3 = TaylorDiagram(stddev3, corrcoeff3, refstd3, 
						 fig=fig, rect=233,  
						 title='(c) Cluster III', titleprops_dict=titleprops_dict,
						 normalize=True, labels=models3, ref_label='Reference')

fig, ax4 = TaylorDiagram(stddev4, corrcoeff4, refstd4, 
						 fig=fig, rect=234, 
						 title='(d) Cluster IV', titleprops_dict=titleprops_dict,
						 normalize=True, labels=models4, ref_label='Reference')

fig, ax5 = TaylorDiagram(stddev5, corrcoeff5, refstd5, 
						 fig=fig, rect=235,
						 title='(e) Cluster V', titleprops_dict=titleprops_dict,
						 normalize=True, labels=models5, ref_label='Reference')
						
ax5.legend(bbox_to_anchor=(1.2, 0.70), loc='lower left', ncol=1)
plt.subplots_adjust(top=0.95)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_stations_cluster_taylor_diagram_{0}_sesa.png'.format(var)
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600)
plt.close('all')
plt.cla()
exit()
