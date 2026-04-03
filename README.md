# Pat's ELC Helper Tool
Pat's ELC Helper Tool is a Python-based GUI utility for preparing Kepler and TESS light curve data for use with the ELC (Eclipsing Light Curve) Fortran software \([Orosz & Hauschildt, 2000](https://ui.adsabs.harvard.edu/abs/2000A&A...364..265O/abstract)\). It automates the process of downloading Kepler and TESS data from MAST, identifying eclipses, detrending the light curves, removing gamma rays from the data, and generating some of the input files required by ELC.

The source code can be viewed [here](https://github.com/patspalding/pat-elc-helper-tool/blob/main/pats_elc_helper_tool.py).

## REQUIREMENTS

An internet connection is required at runtime for downloading data from the MAST archive. Sufficient disk space for data and output files (~300 MB for targets with a large amount of Kepler and TESS data) is also required.

## INSTALLATION

Download the appropriate pre-built executable for your operating system:

  Windows  - pats_elc_helper_tool_windows.exe
  
  macOS    - pats_elc_helper_tool_macos

Windows
  Download pats_elc_helper_tool_windows.exe and place it in a working directory of your choice. Double-click to run, or launch from a terminal:
  
    pats_elc_helper_tool_windows.exe
    
  The tool will create subdirectories in the same folder for downloaded MAST data and generated output files.

macOS
  Download pats_elc_helper_tool_mac and place it in a working directory. The first time you run it, macOS may block the executable as it is from an unidentified developer. To allow it, open Terminal, navigate to the directory containing the file, and run:
    
    chmod +x pats_elc_helper_tool_mac
    
    ./pats_elc_helper_tool_mac
    
  Alternatively, right-click the file in Finder, select Open, and confirm when prompted.

If there are any issues with the executables, or if you wish to run on Linux, the Python source file can be downloaded [here](https://github.com/patspalding/pat-elc-helper-tool/blob/main/pats_elc_helper_tool.py), and executed with:

    python pats_elc_helper_tool.py

This method requires the following installations and libraries:

  Python 3.14+
  
  NumPy 2.4.3+ 
  
  SciPy 1.17.1+
  
  pandas 3.0.1+
  
  matplotlib 3.10.8+
  
  astropy 7.2.0+

## USAGE GUIDE / UI WALKTHROUGH

On launch, a GUI window will appear with the following sections:

### DATA SETTINGS

  #### Target Name
  
  Enter the name of the target system (e.g. "KIC 12345678", "TIC 987654321", or a resolvable name like "V404 Cyg"). Kepler targets must be resolvable to a KIC number via SIMBAD; TESS targets must resolve to a TIC number.

  #### Instrument Data to Use
    
  Select one or both of:
    
  ##### Kepler
  Downloads and processes Kepler long-cadence (LC) and short-cadence (SC) light curve FITS files.
    
  ##### TESS
  Downloads and processes TESS light curve FITS files.

### ECLIPSE TIMES
  #### Find Eclipse Times from...
  
  ##### Fitting Light Curve Data
  The tool will search for eclipses in the data using SciPy peak detection and Gaussian fitting. Requires the following bounds to be specified:
  
  Eclipse Depth Min and Max (Fractional Flux) - The minimum and maximum fractional flux depth of eclipses to search for. Values must be between 0 and 1 (e.g. 0.01 and 0.5).

  Eclipse Duration Min and Max (Days) - The minimum and maximum duration (beginning of ingress to end of egress) of eclipses in days (e.g. 0.05 and 0.5).

  ##### Reading ELC Eclipse Times File
  The tool will read eclipse midpoint times from previously generated ELC eclipse times output files. Requires the path to the folder containing those files to be specified in the "Folder Containing Eclipse Data" field.

  #### Eclipse Times to Remove from Data 
  Enter a space-separated list of eclipse midpoint times (in BJD - 2,455,000) to exclude from processing, e.g. for known contaminated eclipses.

### GENERATE FILES

  Select which output files to generate (see OUTPUT FILES section below):

  #### Eclipse Time Files
  Generates a file of eclipse midpoint times.

  #### ELCgap.inp File
  Generates an ELCgap.inp file, which contains a list of start and end times for "gaps" that ELC will ignore when making forward models. 

  #### ELCSC.inp File
  Generates ELCSC.inp file, which contains a list of start and end times where the data are in short-cadence. Only available if Kepler data are being used.

  #### Individual Eclipses to Seperate Files
  Writes the detrended light curve for each eclipse to its own separate data file.

  #### Instrument Data to File
  Writes the full detrended light curve for each instrument to its own file.

### PLOT SETTINGS

  #### Plot Detrended Data
  If checked, displays the light curve plot of the detrended data after processing. 

  #### Save Plots?
  If checked, saves the plot as a PNG to the images/ subdirectory.

### VERBOSITY
  Controls how much information is printed to the console during processing:

  Quiet - Prints only essential status messages and critical errors
  
  Normal - Prints standard progress updates (default)
  
  Verbose - Prints detailed per-eclipse and per-file debug information

### RUN / EXIT

  #### Run
  Validates inputs and begins processing. The button is disabled while the tool is running to prevent duplicate runs.
  
  #### Exit
  Closes the tool.

### CONSOLE OUTPUT
  The scrollable text area at the bottom of the window displays status messages, warnings, and errors generated during processing.

## OUTPUT FILES

All output files are written to a subdirectory named after the target (lowercase, spaces replaced with underscores), e.g. "kic_12345678/".

  ### {target}/{target}_{instrument}_photo.dat
  
  Tab-separated file containing the detrended light curve for a given instrument. Columns: time (BJD - 2,455,000), normalised flux, flux error.

  ### {target}/ELCSC.inp
  
  Tab-separated file listing the start and end times of each Kepler short-cadence data segment. Used by ELC to identify SC time ranges. Only generated when Kepler data is selected.

  ### {target}/ELCgap.inp
    
  Tab-separated file listing time gaps larger than 1 day in the combined light curve. Each row contains the end time of one segment and the start time of the next. Used by ELC to skip over gaps in data for efficiency.

  ### {target}/{target}_{type}time.dat
  
  Tab-separated eclipse time files, one per eclipse type (e.g. prim, sec, primtrans, sectrans). Columns: eclipse number, midpoint time (BJD - 2,455,000), midpoint error. Intended to be used as eclipse time input file to ELC. If an ELC eclipse time file is used to find eclipse times, eclipse type will be based from the ELC eclipse time filename. Otherwise, all eclipses will be designated type "ecl" and output to a single file.
    
  ### {target}/ecl_data/{target}_{type}{N}.dat
  
  Individual eclipse data files, one per eclipse event. {N} is an arbitrary number used to differentiate eclipses and should not be confused with  the eclipse number based on period.  Columns: time, normalised flux, flux error. Intended to be used for detailed midpoint fitting in ELC. 
    
  ### images/{title}.png
  
  Light curve plot images, saved if "Save Plots?" is checked.

MAST data is downloaded to a mast_data/ subdirectory in the working directory and is cached so subsequent runs do not re-download existing files.

## SUMMARY OF TOOL FUNCTIONS

The tool first attempts to pull all available MAST data through astroquery using the input target name and selected instruments. 

If any files or images are selected to be generated, the tool creates the required directories using pathlib.

The tool then loops through each selected instrument. It first identifies all downloaded MAST data in mast_data/ subdirectory for the instrument by using SIMBAD to resolve the given target name to a KIC (if Kepler)/TIC (if TESS) identifier.

If there is no available MAST data for an instrument for the target in the mast_data/ subdirectory, it skips that instrument.

The tool then attempts to find the eclipse midpoint times, either by fitting the light curve, or reading eclipse midpoints from an ELC eclipse times output file.

If the "Fitting Light Curve Data" option is selected, the tool uses scipy.signal.find_peaks on the inverted FITS file data with the eclipse duration and depth bounds given in order to identify local minima. The tool then fits each identified local minimum with an inverted Gaussian multiplied by a polynomial (to account for the non-flat flux baseline) using curve_fit. If then fit fails, or the derivated standard deviation error on the eclipse depth or duration is greater than 5%, the midpoint is rejected.

If the "Reading ELC Eclipse Times File" option is selected, the tool reads the midpoints from ELC eclipse times output files from the input directory. The type for each eclipse is determined by the filename. The tool then derives the duration and midpoint error for each eclipse by fitting an inverted Gaussian multiplied by a polynomial (to account for the non-flat flux baseline) using curve_fit.

If the tool was unable coherent midpoints and durations from the previous step, it skips to the next instrument.

The tool then detrends the data using the eclipse midpoints and durations  from the previous steps. The tool masks out one duration around each eclipse, centered on the midpoint, then fits a polynomial to the remaining flux data within 3 durations of each midpoint. The data are rejected if there is insufficient coverage on one side or if the detrended out-of-eclipse flux data still have a slope magnitude greater than 0.005. On a successful detrend, the in-eclipse flux and error are divided by the best-fit polynomial to the out-of-eclipse data and appended to the flux_data dictionary. The flux_data dictionary is output as a pandas DataFrame with columns "time", "flux", and "flux_err".

If no flux data are returned by the detrender, the tool skips to the next instrument.

The tool then attempts to remove bad data points caused by gamma rays. A data point is determined to have been caused a gamma ray if it is 5 standard deviations above the median out-of-eclipse flux data.

The tool then cuts out any eclipses within 1 duration of any of the times entered by the user in the "Eclipse Times to Remove from Data" field.

The tool then cuts out any non-eclipse data from the flux_data DataFrame, then saves it for plotting later. If the user has selected the "Instrument Data to File?" option, the tool writes the instrument data to its own file.

After the tool finishes these processes for each instrument, it writes individual eclipse data files, the ELCgap.inp file, and the eclipse times files, if the relevant options were selected. If the relevant option was selected, it plots the detrended light curve for each instrument, and saves the plots if the "Save Plots?" option was selected.

## KNOWN ISSUES / LIMITATIONS

### Gaussian eclipse fitting
  The tool uses a Gaussian to approximate the eclipse shape. Because of this, the tool may not accurately estimate eclipse depth or duration, or error on the eclipse midpoint. It is recommended that the user use ELC for more rigorous eclipse time fitting using the individual eclipse data files output by the tool.

### macOS Matplotlib issues
  Currently, the tool uses matplotlib.pyplot.show() when displaying plots, which can cause crashes on macOS.   It is recommended that the user not enable the "Show Final Light Curve Plots?" option when running on macOS, but instead save plots and view them in a compatible image viewer.

### Windows / Python 3.13+ tkinter GC warnings
  On Windows with Python 3.13 and later, benign error messages of the form "Exception ignored in: <function Variable.__del__>" may appear in the terminal window. These are a known CPython bug and do not affect the correctness of results. These should be suppressed within the tool's console output automatically.

### Eclipse type
  If the "Fitting Light Curve Data" option is selected for finding eclipse midpoint times, the tool will not distinguish between eclipses of different types. This mostly causes issues when estimating the eclipse number (see "Eclipse number estimation", below).
  
### Eclipse number estimation
  If only one eclipse is present in the data, or if the "Fitting Light Curve Data" option is selected for finding eclipse midpoint times, eclipse number estimation will be unreliable as it is based on the minimum interval between each eclipse type. Manual review of the output eclipse time files is advised.

### Eclipse detection sensitivity
  When finding eclipse midpoints from fitting data, the tool may miss eclipses, or produce false positives in noisy data. Using the "Reading ELC Eclipse Times File" option if there are known midpoints will typically give more reliable results.

### Detrending coverage requirement
  The detrender requires a minimum of 5 data points on each side of an eclipse midpoint in the out-of-eclipse region. Incomplete eclipses and  eclipses near the start or end of a data segment may be skipped.

### TESS and Kepler data combined runs

  When both instruments are selected, eclipse detection is run independently per instrument. If eclipse depths or durations differ significantly between instruments (due to e.g. contamination), the same bounds may not be optimal for both. Consider running instruments separately if results are inconsistent.
