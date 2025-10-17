"""
Read spectra for HWO SNR Calculator
Author:    A. G. Sreejith
Version:  0.5 18.09.2025     Initial beta release
Version:  0.5 15.10.2025     Replaced pysznphot with synphot
Notes:    Modified from HWO tools for more spectra options    

"""
import os
import numpy as np
import synphot as S 
import stsynphot as St 
import astropy.units as u
from astropy.io import ascii
from synphot.models import Empirical1D

def add_spectrum_to_library(requested=None):
    """
    Parameters
    ----------
    requested : Selected source name, default is None.

    Returns
    -------
    spec_dict : Araay containing the requested spectra.

    """

    cwd = "../"

    # Define the bandpass (GALEX FUV)
    band = St.band('galex,fuv')
    spec_dict = {}
    
    # --- Special cases (fixed spectra) ---
    special_sources = {
      "Classical T Tauri": "CTTS_etc_d140pc_101116.txt",
      "M1 Dwarf": "dM1_etc_d5pc_101116.txt",
      "10 Myr Starburst":"10Myr_Starburst_nodust.dat",
      "Galaxy with f_esc, HI=1, HeI=1":"fesc/fe_lyccontrot1.000hi1.000hei.txt",
      "Galaxy with f_esc, HI=0.001, HeI=1'":"fesc/fe_lyccontrot0.001hi1.000hei.txt"
    }

    for name, fname in special_sources.items():
        if (requested is None) or (requested in name):
            if (name == "Galaxy with f_esc, HI=1, HeI=1" or 
                name == "Galaxy with f_esc, HI=0.001, HeI=1"):
                tab = ascii.read(os.path.join(cwd, "data/source", fname))
                sp = S.SourceSpectrum(Empirical1D, points=tab["wave"], 
                                      lookup_table=tab["flux"])
                #sp  = S.ArraySpectrum(wave=tab['lam'], flux=tab['lh1=17.5'], 
                #                     waveunits='Angstrom', fluxunits='flam')
            else :
                tab = ascii.read(os.path.join(cwd, "data/source", fname),
                                 names=["wave", "flux"])
                sp = S.SourceSpectrum(Empirical1D, points=tab["wave"], 
                                      lookup_table=tab["flux"])
                #sp  = S.ArraySpectrum(wave=tab["wave"], flux=tab["flux"],
                #                 waveunits="Angstrom", fluxunits="flam")
            trgt = sp.normalize(21.0 * u.ABmag, band,force=True)
            #trgt = sp.renorm(21., "abmag", S.ObsBandpass("galex,fuv"))
            spec_dict[name] = trgt
            spec_dict[name].wave = np.array(tab['wave'])
            spec_dict[name].flux = trgt(tab['wave'] * u.AA, 
                                        flux_unit=S.units.FLAM).value
            

    fits_sources = {
        "QSO": os.path.join(cwd, "data/source/models", "grid", "agn", 
                            "qso_template.fits"),
        "O5V Star": os.path.join(cwd, "data/source/models", "grid", "pickles", 
                                 "dat_uvk", "pickles_uk_1.fits"),
        "G2V Star": os.path.join(cwd, "data/source/models", "grid", "pickles",
                                  "dat_uvk", "pickles_uk_26.fits"),
        "Orion Nebula": os.path.join(cwd, "data/source/models", "grid", 
                                     "galactic", "orion_template.fits"),
        "G191B2B (WD)": os.path.join(cwd, "data/source/models","calspec", 
                                     "g191b2b_mod_010.fits"),
        "GD71 (WD)": os.path.join(cwd, "data/source/models","calspec",
                                  "gd71_fos_003.fits"),
        "GD153 (WD)": os.path.join(cwd, "data/source/models","calspec",
                                   "gd153_fos_003.fits"),
        "Starburst, No Dust": os.path.join(cwd, "data/source/models", "grid", 
                                           "kc96", "starb1_template.fits"),
        "Starburst, E(B-V) = 0.6": os.path.join(cwd,"data/source/models","grid", 
                                           "kc96","starb6_template.fits")
                
    }

    for name, path in fits_sources.items():
        if (requested is None) or (requested in name):
            sp = S.SourceSpectrum.from_file(path)
            sp = sp.normalize(21.0 * u.ABmag, band)
            spec_dict[name] = sp
            wave = sp.waveset.to_value(u.AA)
            flux = sp(wave, flux_unit=S.units.FLAM).value
            spec_dict[name].wave = wave
            spec_dict[name].flux = flux


    # --- Regular stars ---
    stars = ["Wasp-33","Procyon","Tau_ceti","HD120411",
             "V-esp-eri","Trappist-1","HD40307","GJ832","GJ729"]

    for star in stars:
        if (requested is None) or (requested in star):
            fname = f"{star}.dat"
            tab = ascii.read(os.path.join(cwd, "data/source", fname),
                             names=["wave","flux"])
            
            #sp = S.ArraySpectrum(wave=tab["wave"], flux=tab["flux"],
            #                     waveunits="Angstrom", fluxunits="flam")
            sp = S.SourceSpectrum(Empirical1D, points=tab["wave"]*u.AA,
                                  lookup_table=tab["flux"]*S.units.FLAM)
            spec_dict[star] = sp
            spec_dict[star].wave = np.array(tab['wave'])
            spec_dict[star].flux = sp(tab['wave'] * u.AA, 
                                      flux_unit=S.units.FLAM).value
            

    # --- Models ---
    models = ['M7','M3','M1','K2','K2V','G8','G2V','G0V','F5','A8']    
    
    for model in models:
        if (requested is None) or (requested in model):
            print(model)
            fname = f"{model}.dat"
            tab = ascii.read(os.path.join(cwd, "data/source", fname),
                             names=["wave","flux"])
            #sp = S.ArraySpectrum(wave=tab["wave"], flux=tab["flux"],
            #                     waveunits="Angstrom", fluxunits="flam")
            sp = S.SourceSpectrum(Empirical1D, points=tab["wave"]*u.AA,
                                  lookup_table=tab["flux"]*S.units.FLAM)
            spec_dict[model] = sp
            spec_dict[model].wave = np.array(tab['wave'])
            spec_dict[model].flux = sp(tab['wave'] * u.AA, 
                                      flux_unit=S.units.FLAM).value


    #flatsp = S.FlatSpectrum(21, fluxunits='flam')
    #flat = flatsp.renorm(21., 'abmag', S.ObsBandpass('galex,fuv'))
    #spec_dict['Flat in F_lambda'] = flat  

    return spec_dict


   