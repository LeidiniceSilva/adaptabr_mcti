# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "June 10, 2024"
__description__ = "Compute statistical metrics"

import math  
import numpy as np
import numpy.ma as ma
import scipy.stats as st

from scipy.stats import norm


def compute_nrmse(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Root Mean Square Error
    """
    
    model = np.array(model)
    obs = np.array(obs)
    
    p1 = np.nanmean((model - obs)**2)
    p2 = np.sqrt(p1)
    p3 = np.nanmax(obs) - np.nanmin(obs)
    nrmse = p2 / p3

    return nrmse
       

def compute_tss(obs, model):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Taylor Skill Score
    """
    
    p1 = ma.corrcoef(ma.masked_invalid(obs), ma.masked_invalid(model))[0][1]
    p2 = np.nanstd(obs, ddof=0)
    p3 = np.nanstd(model, ddof=0)
    p4 = p3/p2
    p5 = (1 + p1)**4
    p6 = 1/p4
    p7 = (p4 + p6)**2
    p8 = 4*p7
    tss = p5 / p8
    
    return tss
	
	
def compute_corr(obs, model):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Pearson Correlation Coefficient
    """
    
    pcc, pvalue = st.pearsonr(obs, model)
    r_squared = pcc ** 2

    return r_squared
           
	   
def compute_ivs(obs, model):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Interannual Variability Skill Score
    """
    
    p1 = np.nanstd(obs, ddof=0)
    p2 = np.nanstd(model, ddof=0)
    p3 = p2 / p1
    p4 = p1 / p2
    ivs = (p3 - p4)**2  
    
    return ivs    
    


	
	


