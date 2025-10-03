"""
Planet simulation for HWO SNR Calculator 
Author:    A. G. Sreejith
Version:  0.5 17.09.2025     Initial beta release    
"""

from astropy.io import ascii
from random import random
import numpy as np
import pollux_snr_util as psu 


cwd = "../"


def get_planet_radius(wave,flux,snr,planet_name,transit_duration,res,extime):
    """
    Parameters
    ----------
    wave:              wavelength of the source.
    flux:              flux of the source.
    snr:               snr of the source.
    planet_name:       name of the planet
    transit_duration:  transit duration in hours.
    res:               instrument resolution.
    extime:            exposure time in hours 

    Returns
    -------
    transmission spectrum with and without errors

    """

    tab    = ascii.read(cwd+'data/planets/'+planet_name+'_trans.txt', 
                        delimiter=',',comment='#', names=['wave','FinFout'])
    fin_fout  = tab["FinFout"]
    pwave     = tab["wave"]
    radius = (1-fin_fout)
    if np.all(np.diff(pwave) < 0): 
        # Sort the table by the wave 
        sort_ind = np.argsort(pwave)
        pwave    = pwave[sort_ind]
        radius   = radius[sort_ind]
    rprs_in = np.interp(wave, pwave, radius, left=0., right=0.)
    dw = wave/res
    st       = psu.find_nearest(wave, 1000)
    en       = psu.find_nearest(wave, 2100)
    tp_in = psu.gaussbroad(wave,rprs_in,np.mean(dw[st:en])/2.0)
    in_transit_flux = flux*tp_in    
    N1 = np.zeros(len(wave)) 
    N2 = np.zeros(len(wave))
    R1 = np.zeros(len(wave)) 
    R2 = np.zeros(len(wave))
    tp_unc = np.zeros(len(wave))

    #number of exposures in transit duration
    t_dur         = transit_duration*3600. #in seconds
    obs_time      = extime*3600.  #in seconds
    obs_transit   = t_dur/obs_time         #number of observations in transit
    
    #transit depth calculations per transit
    unc_transit_full=np.divide(np.sqrt(2.), (np.sqrt(obs_transit)*snr),  
                              where=(np.sqrt(obs_transit)*snr) != 0)
    #unc_radius_full    = (rprs_in/(np.sqrt(2.))*np.sqrt(obs_transit)*snr)
    
    for i in range((len(wave))):
        
        value1 = np.random.normal(0, 1, 1)
        R1[i] = value1
        N1[i]= np.divide(value1, snr[i], where= snr[i] !=0)
    
        value2 = np.random.normal(0, 1, 1)
        R2[i]= value2
        N2[i] = np.divide(value2,snr[i], where= snr[i] !=0)
        value3 = random()
        tp_unc[i] = abs(value3*unc_transit_full[i])
    
    out_tp = np.divide(((N1*in_transit_flux)+in_transit_flux),((N2*flux)+flux),
                       where= ((N2*flux)+flux) !=0) 
    
    return(out_tp,tp_in)


