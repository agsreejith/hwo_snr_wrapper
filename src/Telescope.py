"""
Define a telescope class for HWO with instrument clases
Author:    A. G. Sreejith
Version:  0.5 17.09.2025     Initial beta release
Notes:    Modified from HWO tools for POLLUX and LUMOS    
"""
from __future__ import print_function
import numpy as np 
from astropy.table import Table 


cwd = "../"

class Telescope: 

    def __init__(self, aperture,temperature,diff_limit_wavelength):
        self.name = 'LUVOIR' 
        self.aperture = aperture # aperture in meters 
        self.temperature = temperature # temperature in Kelvin 
        self.ota_emissivity = 0.09 # emissivity factor for a TMA 
        self.diff_limit_wavelength = diff_limit_wavelength # in nanometers 
        diff_limit_in_arcsec = (1.22*(self.diff_limit_wavelength*0.000000001)*
                                206264.8062/self.aperture)

class Camera(): 

    def __init__(self): 

        self.name = 'HDI' 
        self.pivotwave = np.array([155., 228., 360., 440., 550., 640., 790., 
                                   1260., 1600., 2220.])
        self.bandnames = ['FUV', 'NUV', 'U','B','V','R','I', 'J', 'H', 'K'] 
        self.R_effective = np.array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5.])

        self.ab_zeropoint = np.array([35548., 24166., 15305., 12523., 10018., 
                                      8609., 6975., 4373., 3444., 2482.])

        self.total_qe = np.array([0.1, 0.1, 0.15, 0.45, 0.6, 0.6, 0.6, 0.6, 
                                  0.6, 0.6])
        self.aperture_correction = np.array([1., 1., 1., 1., 1., 1., 1., 1., 
                                             1., 1.])
        self.bandpass_r = np.array([5., 5., 5., 5., 5., 5., 5., 5., 5., 5.])
        self.derived_bandpass = self.pivotwave / self.bandpass_r
        self.dark_current = np.array([0.0005, 0.0005, 0.001, 0.001, 0.001, 
                                      0.001, 0.001, 0.002, 0.002, 0.002])
        self.detector_read_noise = np.array([3., 3., 3., 3., 3., 3., 3., 4., 
                                             4., 4.])

        self.pixel_size = np.array([0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 
                                    0.016, 0.04, 0.04, 0.04])

        
    def set_pixel_sizes(self, telescope): 

        self.pixel_size = (1.22*(self.pivotwave*0.000000001)*206264.8062/
                           telescope.aperture / 2.) 
        # this enforces the rule that the pixel sizes are set at the shortest wavelength in each channel 
        self.pixel_size[0:2] = (1.22*(self.pivotwave[2]*0.000000001)*
                                206264.8062/telescope.aperture / 2.)           # UV set at U 
        self.pixel_size[2:-3] = (1.22*(self.pivotwave[2]*0.000000001)*
                                 206264.8062/telescope.aperture / 2.)          # Opt set at U 
        self.pixel_size[-3:] = (1.22*(self.pivotwave[7]*0.000000001)*
                                206264.8062/telescope.aperture / 2.)           # NIR set at J 

class Spectrograph():
    """
    Class object for LUMOS spectrograph
    """
        
    def __init__(self): 

        self.name = 'LUMOS' 
        """
        lumos = Table.read(cwd+'/data/instrument/LUMOS_vals.dat', format='ascii') 
        self.wave = lumos['Wave']
        self.aeff = lumos['A_eff']
        self.bef = lumos['Med_Res_BEF'] 
        self.med_bef = lumos['Med_Res_BEF'] 
        self.low_bef = lumos['Low_Res_BEF'] 
        self.delta_lambda = self.wave / 30000. #  EXTREMELY ROUGH resel width 
        self.lumos_table = lumos 
        self.mode_name = 'G150M' 
        self.R = 30000. 
        """ 
    def set_mode(self, mode_name): 

        self.mode_names = mode_name 
        """
        old parameter reader         
        if 'G120M' in mode_name:
            print('Setting the spectrograph to mode: ', mode_name) 
            self.bef = self.lumos_table['Med_Res_BEF'] 
            self.delta_lambda = self.wave / 30000. 
            self.lambda_range = np.array([1000., 1425.]) 
            self.mode_name = 'G120M' 
            self.R = 30000. 
           
        if 'G150M' in mode_name: 
            print('Setting the spectrograph to mode: ', mode_name) 
            self.bef = self.lumos_table['Med_Res_BEF'] 
            self.delta_lambda = self.wave / 30000. 
            self.lambda_range = np.array([1225., 1600.]) 
            self.mode_name = 'G150M' 
            self.R = 30000. 
          
        if 'G180M' in mode_name: 
            print('Setting the spectrograph to mode: ', mode_name) 
            self.bef = self.lumos_table['Med_Res_BEF'] 
            self.delta_lambda = self.wave / 30000. 
            self.lambda_range = np.array([1550., 1900.]) 
            self.mode_name = 'G180M' 
            self.R = 30000. 
          
        if 'G155L' in mode_name: 
            print('Setting the spectrograph to mode: ', mode_name) 
            self.bef = self.lumos_table['Low_Res_BEF'] 
            self.delta_lambda = self.wave / 5000. 
            self.lambda_range = np.array([1000., 2000.]) 
            self.mode_name = 'G155L' 
            self.R = 5000.  

        if 'G145LL' in mode_name: 
            print('Setting the spectrograph to mode: ', mode_name) 
            self.bef = self.lumos_table['LL_mode_BEF'] 
            self.delta_lambda = self.wave / 500. 
            self.lambda_range = np.array([900., 2000.]) 
            self.mode_name = 'G145LL' 
            self.R = 500. 
            """
        # new parameter reader
        if 'G110M' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           g120m = Table.read(cwd+'/data/instrument/'
                              +'UVI_v1temp_filt110M_MCP_performance_071024.txt', 
                              format='ascii',delimiter=',') 
           self.wave = g120m['Wavelength']
           self.aeff = g120m['A_Eff']
           self.bef = g120m['BEF'] 
           self.delta_lambda = self.wave / 40000. 
           self.lambda_range = np.array([950., 2100.]) 
           self.mode_name = 'G110M' 
           self.R = 40000.     
        
        if 'G120M' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           g120m = Table.read(cwd+'/data/instrument/'
                              +'UVI_v1temp_G120M_MCP_performance_071024.txt', 
                              format='ascii',delimiter=',') 
           self.wave = g120m['Wavelength']
           self.aeff = g120m['A_Eff']
           self.bef = g120m['BEF'] 
           self.delta_lambda = self.wave / 40000. 
           self.lambda_range = np.array([950., 2100.]) 
           self.mode_name = 'G120M' 
           self.R = 40000.    

        if 'G140M' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           e145h = Table.read(cwd+'/data/instrument/'+
                              'UVI_v1temp_filt140M_MCP_performance_071024.txt',
                              format='ascii',delimiter=',') 
           self.wave = e145h['Wavelength']
           self.aeff = e145h['A_Eff']
           self.bef = e145h['BEF'] 
           self.delta_lambda = self.wave / 40000. 
           self.lambda_range = np.array([1000., 1600.]) 
           self.mode_name = 'G140M' 
           self.R = 40000. 

        if 'E145H' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           e145h = Table.read(cwd+'/data/instrument/'+
                              'UVI_v1temp_E145H_MCP_performance_071024.txt',
                              format='ascii',delimiter=',') 
           self.wave = e145h['Wavelength']
           self.aeff = e145h['A_Eff']
           self.bef = e145h['BEF'] 
           self.delta_lambda = self.wave / 100000. 
           self.lambda_range = np.array([1000., 1600.]) 
           self.mode_name = 'E145H' 
           self.R = 100000.    
           
        if 'G155L' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           e145h = Table.read(cwd+'/data/instrument/'+
                              'UVI_v1temp_G155L_MCP_performance_071024.txt', 
                              format='ascii',delimiter=',') 
           self.wave = e145h['Wavelength']
           self.aeff = e145h['A_Eff']
           self.bef = e145h['BEF'] 
           self.delta_lambda = self.wave / 8000. 
           self.lambda_range = np.array([1000., 2000.]) 
           self.mode_name = 'G155L' 
           self.R = 8000.           
        
        if 'G160M' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           g120m = Table.read(cwd+'/data/instrument/'+
                              'UVI_v1temp_filt160M_MCP_performance_071024.txt',
                              format='ascii',delimiter=',') 
           self.wave = g120m['Wavelength']
           self.aeff = g120m['A_Eff']
           self.bef = g120m['BEF'] 
           self.delta_lambda = self.wave / 40000. 
           self.lambda_range = np.array([950., 2100.]) 
           self.mode_name = 'G160M' 
           self.R = 40000.     

        if 'G180M' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           e145h = Table.read(cwd+'/data/instrument/'+
                              'UVI_v1temp_filt180M_MCP_performance_071024.txt',
                              format='ascii',delimiter=',') 
           self.wave = e145h['Wavelength']
           self.aeff = e145h['A_Eff']
           self.bef = e145h['BEF'] 
           self.delta_lambda = self.wave / 30000. 
           self.lambda_range = np.array([1550., 1900.]) 
           self.mode_name = 'G180M' 
           self.R = 30000.    

        if 'IFU' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           e145h = Table.read(cwd+'/data/instrument/'+
                              'UVI_v1temp_IFU_MCP_performance_071024.txt', 
                              format='ascii',delimiter=',') 
           self.wave = e145h['Wavelength']
           self.aeff = e145h['A_Eff']
           self.bef = e145h['BEF'] 
           self.delta_lambda = self.wave / 5000. 
           self.lambda_range = np.array([950., 2099.]) 
           self.mode_name = 'IFU' 
           self.R = 5000.    
           
           
class Spectropolarimeter(): 
    
    def __init__(self): 

        self.name = 'POLLUX' 
        pollux = Table.read(cwd+'/data/instrument/POLLUX_fuv.csv',
                            format='ascii',delimiter=',') 
        self.wave = pollux['Wave']*10.0
        self.aeff = pollux['EFF']
        self.bef = pollux['BEF'] 
        self.delta_lambda = self.wave / 120750. #  EXTREMELY ROUGH resel width 
        self.pollux_table = pollux 
        self.lambda_range = np.array([1000., 1200.]) 
        self.mode_name = 'FUV' 
        self.R = 30000. 

    def set_mode(self, mode_name): 
        self.mode_names = mode_name   

        if 'FUVPOL' in mode_name:
           print('Setting the POLLUX to mode: ', mode_name) 
           fuvpol = Table.read(cwd+'/data/instrument/POLLUX_fuv.csv', 
                               format='ascii',delimiter=',') 
           self.wave = fuvpol['Wave']*10.0
           self.aeff = fuvpol['EFF']
           self.bef = fuvpol['BEF'] 
           self.lambda_range = np.array([1000., 1200.]) 
           self.mode_name = 'FUV' 
           self.R = 70000.     
           self.delta_lambda = self.wave / self.R 

         
        if 'MUV' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           muv = Table.read(cwd+'/data/instrument/POLLUX_muvpol.csv', 
                            format='ascii',delimiter=',') 
           self.wave = muv['Wave']*10.0
           self.bef = muv['BEF'] 
           self.lambda_range = np.array([1180., 2360.]) 
           self.mode_name = 'MUV' 
           self.R = 70000.
           self.delta_lambda = self.wave /  self.R
           self.dichroic_eff = 0.8 #muv['dichro_RHO']
           self.coat_R       = muv['coat_R']
           self.echelle_R    = muv['echelle_R']
           self.CD_R         = muv['CD_R']
           self.Det_QE       = muv['dDOPED_QE']
           self.mirror_R     = muv['R_percent']    
           aeff = (self.dichroic_eff*self.echelle_R *self.CD_R*self.Det_QE*
                   (self.coat_R**3)*(self.coat_R**3))
           self.aeff = aeff
           
        if 'MUVPOL' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           muvpol = Table.read(cwd+'/data/instrument/POLLUX_muvpol.csv', 
                               format='ascii',delimiter=',') 
           self.wave = muvpol['Wave']*10.0
           self.aeff = muvpol['EFF']
           self.bef = muvpol['BEF'] 
           self.lambda_range = np.array([1180., 2360.]) 
           self.mode_name = 'MUVPOL' 
           self.R = 70000. 
           self.delta_lambda = self.wave /  self.R

        if 'NUV' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           nuv = Table.read(cwd+'/data/instrument/POLLUX_nuv.csv', 
                            format='ascii',delimiter=',') 
           self.wave = nuv['Wave']*10.0
           self.aeff = nuv['Aeff']
           self.bef = nuv['BEF'] 
           self.lambda_range = np.array([2360., 4720.]) 
           self.mode_name = 'NUV' 
           self.R = 70000.
           self.delta_lambda = self.wave /  self.R
           """

        if 'NUVPOL' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           nuvpol = Table.read(cwd+'/data/instrument/POLLUX_nuvpol.txt', 
                               format='ascii',delimiter=',') 
           self.wave = nuvpol['Wave']*10.0
           self.aeff = nuvpol['A_Eff']
           self.bef = nuvpol['BEF'] 
           self.delta_lambda = self.wave / 70000. 
           self.lambda_range = np.array([2360., 4720.]) 
           self.mode_name = 'NUVPOL' 
           self.R = 100000.    
           """
           
        if 'OPT' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           opt = Table.read(cwd+'/data/instrument/POLLUX_opt.csv', 
                            format='ascii',delimiter=',') 
           self.wave = opt['Wave']*10.0
           self.aeff = opt['Aeff']
           self.bef = opt['BEF'] 
           self.lambda_range = np.array([1000., 1600.]) 
           self.mode_name = 'OPT' 
           self.R = 100000.    
           self.delta_lambda = self.wave /  self.R
           """
        if 'NIROPOL' in mode_name:
           print('Setting the spectrograph to mode: ', mode_name) 
           niropol = Table.read(cwd+'/data/instrument/POLLUX_niropol.txt', 
                                format='ascii',delimiter=',') 
           self.wave = niropol['Wave']*10.0
           self.aeff = niropol['A_Eff']
           self.bef = niropol['BEF'] 
           self.delta_lambda = self.wave / 77042. 
           self.lambda_range = np.array([1000., 1600.]) 
           self.mode_name = 'NIROPOL' 
           self.R = 100000.    
           """
           if 'NIR' in mode_name:
              print('Setting the spectrograph to mode: ', mode_name) 
              nir = Table.read(cwd+'/data/instrument/POLLUX_nir.csv', 
                               format='ascii',delimiter=',') 
              self.wave = nir['Wave']*10.0
              self.aeff = nir['Aeff']
              self.bef = nir['BEF'] 
              self.lambda_range = np.array([1000., 1600.]) 
              self.mode_name = 'NIR' 
              self.R = 100000.    
              self.delta_lambda = self.wave /  self.R            