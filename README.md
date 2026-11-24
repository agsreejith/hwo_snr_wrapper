----------------------------------------------------
# README - HWO SNR Calculator Wrapper
----------------------------------------------------
Author:  A. G. Sreejith, IWF, ÖAW, Graz, Austria
Version: 0.5 (17.09.2025) - Initial beta release
Version: 0.6 (31.03.2026) - Updated effective area curves & 
                            wavelength range for all modes.
                          - resolving power of all modes 
                            set to 100K
                          - updated source list
                          - renamed OPT & OPTPOL to VIS & VISPOL
                          - aditional check to see if source spectrum
                            is avaliable & warning to check wave coverage  
Version: 0.7 (23.11.2026) - Created makefile and command line runnable
                          - renamed VIS & VISPOL  to OPT & OPTPOL
                          - browser-based front end (NiceGUI + Plotly)  
----------------------------------------------------
DESCRIPTION
----------------------------------------------------

A wrapper, modified from the [HWO tools](https://github.com/spacetelescope/hwo-tools), for the Habitable Worlds Observatory (HWO) signal-to-noise calculator. It estimates the signal-to-noise ratio (S/N) for the **Pollux** spectropolarimeter on an 8 m LUVOIR-class telescope, and includes a browser-based front end (NiceGUI + Plotly).

**Outputs**

- SNR text files (wavelength, flux, counts, S/N)
- Plots: flux, S/N per resolution element and effective area
- In the web app: interactive figures with downloads for the SNR table (`.txt`) and the figure (`.html`)

## Features

- **Spectroscopy or polarimetry:** switch between the two with one button. The mode list changes with it:
  - Spectroscopy: `MUV`, `NUV`, `OPT`, `NIR`
  - Polarimetry: `FUVPOL`, `MUVPOL`, `NUVPOL`, `OPTPOL`, `NIRPOL`

- **Models or templates:** choose the source type.
  - Models (flux scaled to a distance in parsec): `M7`, `M3`, `M1`, `K2`, `K2V`, `G8`, `G2V`, `G0V`, `F5`, `A8`
  - Templates (used as provided by `get_pollux_spectra`): Procyon, Tau Ceti, HD 120411, V Eps Eri, WASP-33, TRAPPIST-1, GJ 832, GJ 729, HD 40307, Classical T Tauri, G191B2B (WD), GD71 (WD), GD153 (WD), QSO, 10 Myr Starburst, Galaxy with f_esc (HI=1, HeI=1), Galaxy with f_esc (HI=0.001, HeI=1), Orion Nebula, Starburst (No Dust), Starburst (E(B-V) = 0.6)

- **Adjustable inputs:** exposure time, distance (models only), and telescope aperture, temperature and diffraction limit (defaults: 8 m, 280 K, 500 nm).
- **Interactive plots:** Flux, S/N ratios and effective areas as a function of wavelength.

## Directory structure

```
project_root/
├── Makefile
├── pyproject.toml            dependencies and the pollux_etc command
├── data/                     input data (created by make install)
├── output/                   generated files (created by make install)
├── petc/                     virtual environment (created by make install)
└── src/
    ├── app3.py               NiceGUI front end
    ├── pollux_etc.py         launcher behind the pollux_etc command
    ├── get_pollux_spectra.py
    ├── pollux_snr_util.py
    ├── Telescope.py
    ├── planet_sim.py
    └── hwo_snr_wrapper.py    command-line script
```

The front end and the local modules must stay together in `src/`, because they import each other by name.

## Installation

1. Clone or download the repository:
   ```
   git clone https://github.com/agsreejith/hwo_snr_wrapper.git
   cd hwo_snr_wrapper
   ```
2. Check that the files listed above are in place, in particular the four local modules in `src/`: `get_pollux_spectra.py`, `pollux_snr_util.py`, `Telescope.py` and `planet_sim.py`.
3. Create the virtual environment and install everything.

   **With make (Linux and macOS):**
   ```
   make install
   source petc/bin/activate
   ```

   **By hand (any system):**
   ```
   python3 -m venv petc
   source petc/bin/activate          # Windows: petc\Scripts\Activate.ps1
   pip install --upgrade pip
   pip install -e .
   ```
   This installs numpy, matplotlib, pandas, astropy, synphot, stsynphot, nicegui and plotly (Python 3.9 or newer).
4. Install the stsynphot reference data by following the [stsynphot installation instructions](https://stsynphot.readthedocs.io/), then point `PYSYN_CDBS` at it:
   ```
   export PYSYN_CDBS=/path/to/data       # Windows PowerShell: $env:PYSYN_CDBS="C:\path\to\data"
   ```
5. Check the setup:
   ```
   make check
   ```

## Usage

### Web app

```
pollux_etc
```

This opens <http://localhost:8000> in your browser. The command is available inside the activated environment. To use it from any terminal, run `make link`, which symlinks it into `~/.local/bin`.

Options:

| Option | Effect |
|---|---|
| `--port 8080` | Use another port (the next free port is used if it is busy) |
| `--no-browser` | Do not open a browser tab |
| `--host 0.0.0.0` | Make the app reachable from other machines on the network (it has no login, so use it only on trusted networks) |

**Without installing the command**, run the app directly:
```
cd src
python app3.py                  # add --reload while developing, --port 8080 to change the port
```

In the app, pick **Spectroscopy** or **Polarimetry**, choose a mode, set the exposure time, pick a **Model** or **Template** source, and press **Calculate SNR**. Then download the SNR table or the interactive plot.

### Command-line script

```
make script
```

This runs `src/hwo_snr_wrapper.py` from `src/`. Edit the settings at the top of that file (mode, exposure time, source). Results are written to `output/` as `hwo_<source>_<mode>_snr.txt` and `.png`.

## Makefile targets

| Command | What it does |
|---|---|
| `make install` | Create `petc/`, install dependencies and the `pollux_etc` command |
| `make link` / `make unlink` | Add or remove `pollux_etc` in `~/.local/bin` |
| `make launch` | Run `pollux_etc` |
| `make check` | Verify the source files, directories and `PYSYN_CDBS` |
| `make app` / `make app-open` | Run `app3.py` in development mode with auto-reload (the second also opens the browser) |
| `make script` | Run the original command-line script |
| `make clean-output` | Delete generated files in `output/` |
| `make clean` | Remove `__pycache__` and build folders |
| `make distclean` | `clean`, plus delete the virtual environment |

Settings can be overridden on the command line, for example `make app PORT=8080` or `make install VENV=myenv`.

## Troubleshooting

- **`ModuleNotFoundError` for `Telescope`, `pollux_snr_util` and similar:** the app was started outside `src/` or the file is missing. Run `make check`.
- **`ModuleNotFoundError` for a package:** the environment is not activated, or the install did not finish. Run `source petc/bin/activate` and `pip install -e .`.
- **"Source ... was not found in the spectral library":** the names in the app's source lists must match the keys used by `get_pollux_spectra.py`. The error message lists the available keys. Update the `TEMPLATES` or `MODELS` list at the top of `src/app3.py`.
- **No spectra or a stsynphot error:** check that `PYSYN_CDBS` is set and points to the reference data.
- **Spectrum coverage warning:** the chosen source does not cover the full wavelength range of the selected mode, so the result is incomplete in that range.
- **`No module named 'h11._connection'`:** reinstall the web server pieces with `pip install --force-reinstall h11 uvicorn nicegui`.
- **Port already in use:** use `pollux_etc --port 8080`.

## Source Resolution
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

##  Extending the code

- Add new source templates in get_pollux_spectra.py
- Add new instrument modes or configurations in Telescope.py

## Credits

Based on the [HWO tools](https://github.com/spacetelescope/hwo-tools) SNR calculator from the Space Telescope Science Institute.
This code is a wrapper for the HWO (Habitable Worlds Observatory) 
SNR Calculator modified from HWO tools 
(https://github.com/spacetelescope/hwo-tools)

----------------------------------------------------
END OF README
----------------------------------------------------


# hwo_snr_wrapper
