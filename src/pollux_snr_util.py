"""
Functions for HWO SNR Calculator
Author:    A. G. Sreejith
Version:  0.5 17.09.2025     Initial beta release
"""

import numpy as np
import astropy.constants as const



def find_nearest(array, value):
  """
  Find the index of nearest value in an array
  """
  array = np.asarray(array)
  idx   = (np.abs(array - value)).argmin()
  
  return idx


def gaussbroad(w,s,hwhm):
    """
    Smooths a spectrum by convolution with a gaussian of specified hwhm.
    Parameters
    ----------
    w:    wavelength scale of spectrum to be smoothed
    s:    spectrum to be smoothed
    hwhm: half width at half maximum of smoothing gaussian.
    
    Returns
    --------
    a vector containing the gaussian-smoothed spectrum.
     
    Warn user if hwhm is negative.
      if hwhm lt 0.0 then $
      message,/info,'Warning! Forcing negative smoothing width to zero.'
    
    """
    #Return input argument if half-width is nonpositive.
    if hwhm <= 0.0:
        return(s)			                                                   # true: no broadening

    #Calculate (uniform) dispersion.
    dw = (w[-1] - w[0]) / len(w)		                                       
    for i in range(0, len(w)):
        #Make smoothing gaussian # extend to 4 sigma.
        if(hwhm > 5*(w[-1] - w[0])): 
            return np.full(len(w),np.sum(s)/len(w))
        nhalf = int(3.3972872*hwhm/dw)		                                   
        ng = 2 * nhalf + 1				
        wg = dw * (np.arange(ng) - (ng-1)/2.0)	                               # wavelength scale of gaussian
        xg = ( (0.83255461) / hwhm) * wg 		                               
        gpro = ( (0.46974832) * dw / hwhm) * np.exp(-xg*xg)                    # unit area gaussian w/ FWHM
        gpro=gpro/np.sum(gpro)

    #Pad spectrum ends to minimize impact of Fourier ringing.
    npad = nhalf + 2				                                           
    spad = np.concatenate((np.full(npad,s[0]),s,np.full(npad,s[-1])))
    #Convolve & trim.
    sout = np.convolve(spad,gpro,mode='full')			                       
    sout = sout[npad:npad+len(w)]			                                   
    return sout					                                               



def simulate_exposure(telescope, spectrograph, wave, flux, exptime):
    """
    Simulates the exposure and creates snr's and counts
    
    Parameters
    ----------
    telescope:    Class object with telescope information
    spectrograph: Class object with spectrograph paramaeters
    wave:         Wavelength from source
    flux:         Flux from source
    exptime:      Exposure time in hours

    Returns
    -------
    signal_to_noise:   Calculated SNR for the observation
    source_counts:     Calculated source counts at the detector
    bef_interp:        Background flux at instrument resolution
    background_counts: Backgrounbd counts at detector
    
    """     
    print("Attempting to create an exposure for Telescope: ", telescope.name, 
          telescope.aperture, ' m') 
    print("                              and Spectrograph: ", spectrograph.name
          , " in mode ", spectrograph.mode_name) 

    # obtain the interpolated effective areas for the input spectrum 
    aeff_interp = np.interp(wave, spectrograph.wave, spectrograph.aeff, 
                   left=0., right=0.) * np.pi*(telescope.aperture*1.0e2/2.0)**2                  # modified to keep Aeff free of telescope size
    bef_interp  = np.interp(wave, spectrograph.wave, spectrograph.bef, 
                   left=0., right=0.) # background to use 
    phot_energy = const.h.to('erg s').value*const.c.to('cm/s').value/(wave*1e-8)                 # convert from erg cm^-2 s^-1 A^-1  
    source_counts = (flux/phot_energy*aeff_interp*(exptime*3600.)
                     *(wave/spectrograph.R)) 
    source_counts[(wave < spectrograph.lambda_range[0])] = 0. 
    source_counts[(wave > spectrograph.lambda_range[1])] = 0. 
    background_counts = (bef_interp/phot_energy*aeff_interp*(exptime*3600.)
                         *(wave/spectrograph.R)) 
    signal_to_noise   = source_counts / (source_counts + background_counts)**0.5 
    
    return signal_to_noise,source_counts,bef_interp, background_counts


def text_head_generator(spectrograph,template,planet=False,distance=None):
    """
    Creates header information for SNR calculator output.

    Parameters
    ----------
    spectrograph: Calss object with spectrograph parameters.
    template:     Source name as string.
    planet:       Set to true if exoplanet transit is to be simulated, default 
                  is set to False.
    Returns
    -------
    hwo_uv_snr:  A header string for the SNR output file. 

    """
    hwo_uv_snr='--------------------------------------------------------------'
    hwo_uv_snr=hwo_uv_snr+'-----------------\n'
    hwo_uv_snr=hwo_uv_snr+'                      HWO UV SNR CALCULATIONS\n'
    hwo_uv_snr=hwo_uv_snr+'---------------------------------------------------'
    hwo_uv_snr=hwo_uv_snr+'----------------------------\n'
    hwo_uv_snr=hwo_uv_snr+'\n'
    hwo_uv_snr=hwo_uv_snr+'Input Parameters \n'
    hwo_uv_snr=hwo_uv_snr+'Source :'+(template)
    if distance != None:
        hwo_uv_snr=hwo_uv_snr+' @'+str(distance)+' parsec\n'
    hwo_uv_snr=hwo_uv_snr+'\n'
    hwo_uv_snr=hwo_uv_snr+'Instrument Parameters('+spectrograph.name+') \n'
    hwo_uv_snr=hwo_uv_snr+'Mode             : '+spectrograph.mode_name+'\n'
    hwo_uv_snr=hwo_uv_snr+'Resolution       : '+str(spectrograph.R)+'\n'
    hwo_uv_snr=hwo_uv_snr+'Wavelength Range : '+str(spectrograph.lambda_range[0])
    hwo_uv_snr=hwo_uv_snr+'-'+str(spectrograph.lambda_range[1])+'\n'
    hwo_uv_snr=hwo_uv_snr+'---------------------------------------------------'
    hwo_uv_snr=hwo_uv_snr+'----------------------------\n'
    hwo_uv_snr=hwo_uv_snr+'                        HWO SNR OUTPUTS\n'
    hwo_uv_snr=hwo_uv_snr+'---------------------------------------------------'
    hwo_uv_snr=hwo_uv_snr+'----------------------------\n'
    if planet == False:
        hwo_uv_snr=hwo_uv_snr+'Wavelength [A] Flux [ergs/cm2/s/A] Background' 
        hwo_uv_snr=hwo_uv_snr+'[ergs/cm2/s/res] Flux [counts] SNR \n' 
    else:
        hwo_uv_snr=hwo_uv_snr+'Wavelength [A] Input Radius [(Rp/Rs)2] Output'
        hwo_uv_snr=hwo_uv_snr+' Radius [(Rp/Rs)2] \n'
    return hwo_uv_snr

