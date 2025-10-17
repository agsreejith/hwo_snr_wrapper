====================================================
README - HWO SNR Calculator Wrapper
====================================================
Author:  A. G. Sreejith
Version: 0.5 (17.09.2025) - Initial beta release

----------------------------------------------------
DESCRIPTION
----------------------------------------------------

This code is a wrapper for the HWO (Habitable Worlds Observatory) 
SNR Calculator modified from HWO tools 
(https://github.com/spacetelescope/hwo-tools)

It calculates:
- Signal-to-noise ratios (SNR) for the Pollux instrument

Outputs:
- SNR text files
- Diagnostic plots

----------------------------------------------------
INSTALLATION
----------------------------------------------------

1. Clone or download the repository:

   git clone https://github.com/agsreejith/hwo_snr_wrapper.git
   cd hwo_snr_wrapper

2. Create a Python environment (recommended):

3. Install required dependencies:
   
    numpy, matplotlib, pandas, astropy, synphot, stsynphot
    
    After following install instructions at stsynphot
    your PYSYN_CDBS should be properly configured.

4. Make sure the following local scripts are available 
   in your working directory (or in PYTHONPATH):

   - get_pollux_spectra.py
   - pollux_snr_util.py
   - Telescope.py
   - planet_sim.py

----------------------------------------------------
DIRECTORY STRUCTURE
----------------------------------------------------

project_root/
|
|-- data/                        (input data)
|
|-- src/
|   |-- get_pollux_spectra.py
|   |-- pollux_snr_util.py
|   |-- Telescope.py
|   |-- planet_sim.py
|   |-- hwo_snr_wrapper.py       (this main script)
|
|-- output/                      (generated SNR text files and plots)

----------------------------------------------------
AVAILABLE MODES (Pollux)
----------------------------------------------------

The Pollux spectropolarimeter modes currently supported are:

- FUVPOL    : Far Ultraviolet with Polarimetry
- MUV       : Mid Ultraviolet spectroscopy
- MUVPOL    : Mid Ultraviolet with Polarimetry
- NUV       : Near Ultraviolet spectroscopy
- NUVPOL    : Near Ultraviolet with Polarimetry
- OPT       : Optical spectroscopy
- OPTPOL    : Optical with Polarimetry
- NIR       : Near Infrared spectroscopy
- NIRPOL    : Near Infrared with Polarimetry

Select the mode by setting:

   mode_selected = "MUV"    # Example

----------------------------------------------------
AVAILABLE TARGETS
----------------------------------------------------

Templates (sources) are provided via get_pollux_spectra. 
Current available sources include:

- Procyon
- Tau Ceti
- HD 120411
- V Eps Eri
- WASP-33
- TRAPPIST-1
- GJ 832
- GJ 729
- HD 40307
- Classical T Tauri
- O5V Star
- O5V Star
- G2V Star
- M1 Dwarf
- G191B2B (WD)
- GD71 (WD)
- GD153 (WD)
- QSO
- 10 Myr Starburst
- Galaxy with f_esc, HI=1, HeI=1
- Galaxy with f_esc, HI=0.001, HeI=1
- Orion Nebula
- Starburst, No Dust
- Starburst, E(B-V) = 0.6

Models
- M7
- M3
- M1
- K2
- K2V
- G8
- G2V
- G0V
- F5
- A8

Set the source by updating: 

   template = "Trappist-1"
   
if model is selected, set:
 
   model = True

and also specify:

   distance = 10 (in parsec)


----------------------------------------------------
USAGE
----------------------------------------------------

Example workflow:

1. Edit the configuration section in hwo_snr_wrapper.py:
   - Select telescope size and temperature
   - Choose instrument mode
   - Define exposure time
   - Pick source template

2. Run the script:

   python hwo_snr_wrapper.py

3. Outputs will be written to the "output/" folder:
   - Text files with calculated SNR values
   - Plots with SNR curves

----------------------------------------------------
SORCE RESOLUTION
----------------------------------------------------
- M7  variable (above 2k, 100K or more in FUV and IR)  
- M3  variable (above 2k, 100K or more in FUV and IR) 
- M1  variable (above 2k, 100K or more in FUV and IR) 
- K2  variable (above 2k, 50k or more in NUV, 100K or more in FUV and IR)
- K2V variable (above 2k, 10k or more above 1200A, above 50K for most of the spectrum)
- G8  variable (above 2k, 10k or more above 1200A, above 50K for most of the spectrum)
- G2V variable (above 1k, above 500K above 1500A) 
- G0V variable (above 1k, 100k or more in IR)
- F5  variable (above 3k, 60k or more for most of the spectrum) 
- A8  variable (above 100k)
----------------------------------------------------
EXTENDING THE CODE
----------------------------------------------------

- Add new source templates in get_pollux_spectra.py
- Add new instrument modes or configurations in Telescope.py

====================================================
END OF README
====================================================


# hwo_snr_wrapper
