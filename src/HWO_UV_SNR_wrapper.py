"""
HWO SNR Calculator wrapper

Calculates SNR and transit uncertanities for HWO instruments 
based on user inputs.

Produces SNR text files and plots. 


Author:    A. G. Sreejith
Version:  0.5 17.09.2025     Initial beta release    
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import get_pollux_spectra
import pollux_snr_util as psu 
import Telescope as T 
import os 
import planet_sim as pl
import pysynphot as S 

mpl.rcParams['agg.path.chunksize'] = 10000
distance = None

#Parameter set up
#**********************************************************************************
#setting path to current working directory
cwd = "../"
#Constants
parsec2AU = 206264.80624538

#Current avaliable modes for POLLUX: ['FUV','MUV','NUV','OPT','NIR']
luvoir = T.Telescope(8., 280., 500.)                                                # set up LUVOIR with 8 meters, T = 280, and diff limit at 500 nm 
pollux = T.Spectropolarimeter()                                                     # Selecting POLLUX as the HWO instrument
mode_selected = 'MUV'                                                               # Mode selection: FUVPOL,MUV,MUVPOL,NUV,NUVPOL,OPT,OPTPOL,NIR,NIRPOL 
extime        =  1.0                                                                # Exposure time in hours
template      = 'G2V'                                                               # Source name, refer readme for the list of sources

model         = True                                                               # set to True if the source specified is pectral model
if model != False : distance = 10                                                   # Provide distance (in parsec) for scaling flux if templeate is a model file with spectral type    

"""
# Exoplanet calculations
planet           = False                                                            # Set to true to include exoplanet calculations
planet_name      = "KELT9b"
transit_duration = 3.0                                                              # Transit duration in hours
"""

#************************************************************************************
#
#Reading inputs and SNR calculation 
#
#************************************************************************************
pollux.set_mode(mode_selected)                                                      # Setting instrument mode
spec_dict = get_pollux_spectra.add_spectrum_to_library(template)                    # Read all spectra
waves  = spec_dict[template].wave 
fluxes = spec_dict[template].flux

if model != False :
    fluxatE = fluxes/(distance*(parsec2AU)**2)                                      # Scaling flux to Erath (model flux at 1AU from the star)
    spec_dict[template].flux = fluxatE
    spec_dict[template].distance =distance
    

signal_to_noise,source_counts,bef_int,b_counts = psu.simulate_exposure(luvoir,      # Calculate SNR for selected parameters
                                                 pollux,spec_dict[template].wave
                                           ,spec_dict[template].flux, extime) 
flux_cut = spec_dict[template].flux 
flux_cut[spec_dict[template].wave < pollux.lambda_range[0]] = -999.  
flux_cut[spec_dict[template].wave > pollux.lambda_range[1]] = -999.  

#************************************************************************************
#
# Setting up plots
#
#************************************************************************************
fig = plt.figure(301,(15, 10))
gs  = GridSpec(2,2, figure=fig,wspace=0.2,hspace=0.3)

ax  = fig.add_subplot(gs[0,:-1])
plt.plot(spec_dict[template].wave,spec_dict[template].flux,color='red')
ax.set_ylabel(r"$\rm Flux\ [ergs~cm^{2}~s^{-1}~\AA^{-1}]$", fontsize=12)
ax.set_xlim(pollux.lambda_range[0], pollux.lambda_range[1])
st = psu.find_nearest(spec_dict[template].wave, pollux.lambda_range[0])
en = psu.find_nearest(spec_dict[template].wave, pollux.lambda_range[1])
ax.set_ylim(0,max(spec_dict[template].flux[st:en]))
ax.set_xlabel(r"$\rm Wavelength\ [\AA]$", fontsize=12)
plt.title('Flux vs Wavelength ('+template+')')

bx = fig.add_subplot(gs[1,:-1])
plt.plot(spec_dict[template].wave,signal_to_noise, color="black")
bx.set_xlim(pollux.lambda_range[0], pollux.lambda_range[1])
bx.set_xlabel(r"$\rm Wavelength\ [\AA]$", fontsize=12)
bx.set_ylabel(r"$\rm S/N\ per\ resel$", fontsize=12)
plt.title('SNR vs Wavelength ('+template+' in '+pollux.name+' mode '
          +pollux.mode_name+')')

cx = fig.add_subplot(gs[:,-1])
plt.plot(pollux.wave,pollux.aeff*np.pi*(luvoir.aperture*1.0e2/2.0)**2  , 
         color="black")
cx.set_xlabel(r"$\rm Wavelength\ [\AA]$", fontsize=12)
cx.set_ylabel(r"$\rm Effective\ area\ [cm^{2}]$", fontsize=12)

#************************************************************************************
#
#Saving output files
#
#************************************************************************************
if not os.path.exists(cwd+'output/'):
    os.makedirs(cwd+'output/')
snr_header = psu.text_head_generator(pollux,template,**({'distance':distance} 
                                                        if distance is not None 
                                                        else {}))
file_out   = cwd+'output/hwo_'+(template)+'_'+(mode_selected)+'_snr.txt'             # Output filenames have the format: source_mode_snr
figname    = cwd+'output/hwo_'+(template)+'_'+(mode_selected)+'_snr.png'

plt.savefig(figname,dpi=600)    
wave_array = np.asarray(spec_dict[template].wave) 
st         = psu.find_nearest(wave_array, pollux.lambda_range[0])
en         = psu.find_nearest(wave_array, pollux.lambda_range[1])
np.savetxt(file_out, np.c_[spec_dict[template].wave[st:en],
           spec_dict[template].flux[st:en],bef_int[st:en],source_counts[st:en],
           signal_to_noise[st:en]],header=snr_header,
           fmt=['% 9.3f','% 12.5e','% 12.5e','% 12.5e','% 12.5e'])
plt.close()
"""
if planet == True:
    snr_header_planet = psu.text_head_generator(pollux,template,planet)

    out_tp,tp_in      = pl.get_planet_radius(pollux,spec_dict[template].wave,
                                             spec_dict[template].flux,
                                             signal_to_noise,planet_name,
                                             transit_duration,extime)
      
    
    fig = plt.figure(302,(15, 10))
    ax  = plt.subplot(1,1,1)
    plt.plot(spec_dict[template].wave,out_tp,color='lightblue',label='output')
    plt.plot(spec_dict[template].wave,tp_in,color='pink',label='input')
    plt.xlabel(r"$\rm Wavelength\ [\AA]$", fontsize=12)
    plt.ylabel(r"$\rm (R_{p}/R_{s})^2 $", fontsize=12)
    plt.xlim(pollux.lambda_range[0],pollux.lambda_range[1])
    plt.ylim(0,0.1)
    plt.title('Wavelength vs Radius ('+template+')')
    plt.legend()
    plt.tight_layout()
    plfigname = (cwd+'output/hwo_'+template+'_'+(planet_name)+'_'+
                 (mode_selected)+'_snr.png')
    plt.savefig(plfigname,dpi=600)    
    plt.close()

    file_out = cwd+'output/hwo_'+(planet_name)+'_'+(mode_selected)+'_snr.txt'
    np.savetxt(file_out, np.c_[spec_dict[template].wave[st:en],tp_in[st:en],
                               out_tp[st:en]],header=snr_header_planet,
                             fmt=['% 9.3f','% 9.5f','% 9.5f'])
"""
#************************************************************************************
  
print("SNR plots and output file generated for HWO with "+pollux.name+
      " in mode "+pollux.mode_name+" for target "+template+ ".")    
    

