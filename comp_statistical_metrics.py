# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Mar 01, 2023"
__description__ = "Compute statistical metrics"

import math  
import numpy as np
import numpy.ma as ma
import scipy.stats as st

from scipy.stats import norm


def compute_pcc(obs, model):

    """
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
    Return: Pearson Correlation Coefficient
    """
    
    pcc, pvalue = st.pearsonr(obs, model)
    
    return pcc
           

def compute_r2(obs, model):

	"""
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
	Return: R-squared
	"""

	pcc, pvalue = st.pearsonr(obs, model)
	r2 = pcc ** 2

	return r2


def compute_mbe(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Mean Bias Error
    """

    mbe = np.nanmean(np.array(model) - np.array(obs))
    
    return mbe
    
    
def compute_mae(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Mean Absoluty Error
    """

    mae = np.nanmean(np.abs(np.array(model) - np.array(obs)))
    
    return mae
    

def compute_rmse(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Root Mean Square Error
    """
    
    mse = np.nanmean(np.square(np.subtract(obs, model)))
    rsme = math.sqrt(mse)
    
    return rsme
    
   
def compute_pbias(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Percentage Bias Error
    """

    pbias = 100.0 * sum(np.array(model) - np.array(obs)) / sum(np.array(obs))
    
    return pbias
        
    
def compute_apb(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Absolute Percent Bias Error
    """

    apb = 100.0 * sum(np.abs(np.array(model), np.array(obs))) / sum(np.array(obs))
    
    return apb


def compute_nse(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Nash–Sutcliffe Efficient Coefficient
    """
    
    p1 = sum((model - obs) ** 2)
    p2 = sum((obs - np.mean(obs)) ** 2)
    nse = 1 - p1 / p2
    
    return nse
    
        
def compute_diso(model, obs):

	"""
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
	Return: Distance between Indices of Simulation and Observation
	"""

	p0 = np.abs(np.array(model).mean())
	p1 = np.corrcoef(model, obs)[0][1]
	p2 = np.mean(np.abs(np.array(model) - np.array(obs))) / p0
	p3 = np.sqrt(((np.array(model) - np.array(obs)) ** 2).mean()) / p0
	diso = np.sqrt((p1)** 2 + (p2)** 2 + (p3)** 2)

	return diso
    
    
def compute_ioa(model, obs):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Index of Agreement
    """

    p1 = (model - obs)**2
    p2 = np.abs(model - np.mean(obs))
    p3 = np.abs(obs - np.mean(obs))
    ioa = 1 - sum(p1) / sum((p2 + p3)**2)
    
    return ioa
      
    
def compute_ivs(obs, model):

    """
    The input arrays must have the same dimensions
    Param model: Numpy array with model data
    Param obs: Numpy array with obs data
    Return: Interannual Variability Skill Score
    """

    p1 = np.nanstd(obs)
    p2 = np.nanstd(model)
    p3 = p2 / p1
    p4 = p1 / p2
    ivs = (p3 - p4)**2  
    
    return ivs    
    

def compute_kge(obs, model):

	"""
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
	Return: Kling-Gupta Efficiency
	"""

	p1 = np.corrcoef(obs, model)[0][1]
	p2 = np.nanmean(model) / np.nanmean(obs)
	p3 = np.nanstd(model, ddof=1) / np.nanmean(model) * 100 
	p4 = np.nanstd(obs, ddof=1) / np.nanmean(obs) * 100 
	p5 = p3/p4
	kge = 1 - np.sqrt((1 - p1)**2 + (1 - p2)**2 + (1 - p5)**2)

	return kge


def compute_tss(obs, model):

	"""
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
	Return: Taylor Skill Score
	"""

	p1 = ma.corrcoef(ma.masked_invalid(obs), ma.masked_invalid(model))[0][1]
	p2 = np.nanstd(model, ddof=1)
	p3 = (1 + p1)**4
	p4 = 1/p2
	p5 = (p2 + p4)**2
	p6 = 4*p5
	tss = p3 / p6

	return tss
	
	
def compute_cdf(dataset):

	"""
	Param dataset: Numpy array with model or obs dataset
	Return: Cumulative Distribution Function
	"""

	x = np.linspace(np.min(dataset), np.max(dataset))
	y = np.nanmean(x)
	z = np.nanstd(x)
	cdf = norm.cdf(x,y,z)

	return x, cdf


def compute_pdf(dataset):

	"""
	:Param dataset: Numpy array with model or obs dataset
	:Return: Probability Density Function
	"""

	x = np.linspace(np.min(dataset), np.max(dataset))
	y = np.nanmean(x)
	z = np.nanstd(x)
	pdf = norm.pdf(x,y,z)

	return x, pdf
    

