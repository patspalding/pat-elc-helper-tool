# Imports
import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from glob import glob

import warnings
from scipy.optimize import OptimizeWarning
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

from astropy.io import fits
from astroquery.simbad import Simbad
from astroquery.mast import Observations

import time as t
import pandas as pd

import pathlib

import matplotlib.pyplot as plt

import logging
import threading

import tkinter as tk
from tkinter import ttk, messagebox
from tkinter.scrolledtext import ScrolledText

import sys

import os
os.chdir(pathlib.Path(__file__).parent)

def create_directories(target,gen_indiv_ecls,save_plots):
    """
    Create directories to help organize generated files
    """
    
    target = target.lower().replace(" ", "_")
    directories = [f"{target}"]

    if gen_indiv_ecls:
        directories.append(f"{target}/ecl_data")

    if save_plots:
        directories.append("images")

    for directory in directories:
        if not pathlib.Path(directory).exists():
            pathlib.Path(directory).mkdir(parents=True)
            vprint(f"Created folder '{directory}'",level="normal")

def vprint(message, level="normal"):
    """
    Print message only if level is within current verbosity setting.
    """
    levels = {"quiet": 0, "normal": 1, "verbose": 2}
    if levels.get(level, 1) <= levels.get(VERBOSITY, 1):
        print(message)

def get_mast_files(targetname,instruments,verbosity):
    """
    Pull Kepler and TESS timeseries data from MAST based on target name
    """

    verbose = verbosity == "verbose" 

    if not verbose:
        logging.getLogger("astroquery").setLevel(logging.WARNING)    
    else:
        logging.getLogger("astroquery").setLevel(logging.INFO)    

    # Query observations
    obs = Observations.query_object(
        targetname, 
        radius=0.0
    )


    # Keep only timeseries data products
    obs = obs[obs['dataproduct_type'] == 'timeseries']

    # Keep only selected instrument data
    obs = obs[np.isin(obs['obs_collection'], instruments)]

    # Get all products for those observations
    products = Observations.get_product_list(obs)

    # Keep only light curve FITS files
    products = Observations.filter_products(
        products,
        productType="SCIENCE",
        extension="fits"
    )

    # Make array of filenames
    filenames = np.array(products['productFilename'], dtype=str)

    # Mask out files based on ending
    w_kep_lc = np.char.endswith(filenames, "llc.fits")
    w_kep_sc = np.char.endswith(filenames, "slc.fits")
    w_tess = np.char.endswith(filenames, "_lc.fits")
    
    final_mask = w_kep_lc | w_kep_sc | w_tess
    
    products = products[final_mask]

    # Download
    manifest = Observations.download_products(products, download_dir="mast_data", cache=True)

def get_fits_filenames(targetname,mission,cadence=None):
    """
    Get list of filenames depending on target, mission (Kepler or TESS), and cadence (short or long cadence; for Kepler only)
    """
    targetname_upper = targetname.upper()
    mission_upper = mission.upper()
    if cadence is not None: 
        cadence = cadence.upper()

    
    # If KIC or TIC number provided, use that
    if "KIC" in targetname_upper and mission_upper == "KEPLER":
        name = targetname_upper
    elif "TIC" in targetname_upper and mission_upper == "TESS":
        name = targetname_upper

    # Otherwise, query Simbad with the provided name
    else:
        result_table = Simbad.query_objectids(targetname)
    
        for result in result_table:
            result_str = str(result)
            if "KIC" in result_str and mission_upper == "KEPLER":
                name = result_str
            elif "TIC" in result_str and mission_upper == "TESS":
                name = result_str

    try:
        # Get files from mission, cadence, and KIC/TIC number
        if "KIC" in str(name) and mission_upper == "KEPLER":
            # Short cadence Kepler data
            if cadence == 'SC':
                num = f"{int(str(name).split()[-1]):09d}"
                return list(pathlib.Path('.').rglob(f"kplr{num}*slc.fits"))
            
            # Long cadence Kepler data
            elif cadence == 'LC' or cadence == None:
                num = f"{int(str(name).split()[-1]):09d}"
                return list(pathlib.Path('.').rglob(f"kplr{num}*llc.fits"))

        # TESS data
        if mission_upper == 'TESS' and 'TIC' in str(name):
            num = f'{int(str(name).split(" ")[-1]):016d}'
            return list(pathlib.Path('.').rglob(f"tess*{num}*_lc.fits"))

    except UnboundLocalError:
        return list()


def filter_data_rel(fnames):
    """
    Return filenames determined to be reliable based on the DATA_REL parameter in each file's header
    """
    return fnames[np.array([fits.getheader(fname)['DATA_REL'] == 25 for fname in fnames])]

def func(x,a,b,c,d):
    """
    Third-degree polynomial for fitting
    """
    return a*x**3 + b*x**2 + c*x + d

def gaussian(t,A,b,c):
    """
    Inverted Gaussian with flat base for fitting; A is eclipse depth, b is eclipse time, c is ~eclipse duration
    """
    
    return 1 - A*np.exp(-1/2*((t-b)/c)**2)

def get_instrument_qual_time(img_data,instrument):
    """
    Get time and quality depending on the instrument
    Filter time based on quality and adjust given instrument timestamps to be in BJD - 24550000
    """
    instrument = instrument.upper()
    
    if instrument == 'KEPLER':
        qual = img_data['SAP_QUALITY'] <= 16
        time = img_data['TIME'][qual] - 167
        
    elif instrument == 'TESS':
        qual = img_data['QUALITY'] == 0
        time = img_data['TIME'][qual] + 2000

    return qual, time    

def make_in_eclipse_mask(time,ecl_data):
    """
    Make mask of midpoints within the given time, and in-eclipse time within the given time, based on the data in ecl_data
    """

    time = np.array(time)

    midpoints = ecl_data["midpoint"].values
    durations = ecl_data["duration"].values
        
    # Get eclipses only within time range
    w_midpoints = np.any(
        (midpoints[:, None] >= np.min(time)) &
        (midpoints[:, None] <= np.max(time)),
        axis=1
    )


    # Get in-eclipse mask
    w_in_eclipse = np.any(
        (time[:, None] >= midpoints[w_midpoints] - 0.5*durations[w_midpoints]) &
        (time[:, None] <= midpoints[w_midpoints] + 0.5*durations[w_midpoints]),
        axis=1
    )

    return w_midpoints, w_in_eclipse

def load_flux_data(fname,instrument):
    # Get image data
    img_data = fits.getdata(fname, ext=1)

    # Filter data based on quality; adjust given timestamps to be in BJD - 24550000
    qual, time = get_instrument_qual_time(img_data,instrument)
    flux = img_data['SAP_FLUX'][qual]
    err = img_data['SAP_FLUX_ERR'][qual]

    w_nan = np.isnan(flux) | np.isnan(time) | np.isnan(err)
    time = time[~w_nan]
    flux = flux[~w_nan]
    err = err[~w_nan]

    return time, flux, err

def mdpt_detrend(fnames,ecl_data,instrument):
    """
    Detrend the data based on a list of eclipse midpoints
    ecl_data is the eclipse data DataFrame which is used to mask out in-eclipse data for detrending
    instrument is the mision name (Kepler or TESS)
    """

    midpoints = ecl_data["midpoint"].values
    durations = ecl_data["duration"].values

    # Create dict to hold flux data
    flux_data = {"time":[],"flux":[],"flux_err":[]}

    # Loop through each data file
    for fname in fnames:
        # Get data from FITS file
        time, flux, err = load_flux_data(fname,instrument)

        # Get in-eclipse masks
        w_midpoints, w_in_eclipse = make_in_eclipse_mask(time,ecl_data)
        
        # Skip file if there is no eclipse data
        if time[w_in_eclipse].size == 0:
            vprint(f"No eclipse data. Skipping file {fname}",level="verbose")
            continue

        # Loop through all eclipses within time range
        for idx, row in ecl_data[w_midpoints].iterrows():
            
            # Get data around eclipses within 3 durations of the midpoint
            w_range = (time < (row["midpoint"] + 3*row["duration"])) & (time > (row["midpoint"] - 3*row["duration"]))
            time_range = time[w_range]
            flux_range = flux[w_range]
            err_range = err[w_range]

            # Isolate out-of-eclipse data           
            out_ecl_time = time[~w_in_eclipse & w_range]
            out_ecl_flux = flux[~w_in_eclipse & w_range]

            # Ensure eclipse has coverage on both sides
            left_coverage = np.sum(out_ecl_time < row["midpoint"])
            right_coverage = np.sum(out_ecl_time > row["midpoint"])
            min_pts = 5
            if left_coverage < min_pts or right_coverage < min_pts:
                vprint(f"Skipping eclipse at {row['midpoint']:.4f}: insufficient coverage on one side.",level="verbose")
                continue

            try:
                proc_flux,proc_err = poly_detrender(time_range-row["midpoint"],flux_range,err_range,out_ecl_time-row["midpoint"],out_ecl_flux)

                w_range_midpoints, w_range_in_eclipse = make_in_eclipse_mask(time_range,ecl_data)
                out_ecl_detrended = proc_flux[~w_range_in_eclipse]

                # Slope of detrended out-of-eclipse flux (should be ~0 if flat)
                slope, intercept = np.polyfit(out_ecl_time, out_ecl_detrended, 1)
                if abs(slope) > 0.005:
                    vprint(f"Bad detrend. Data near {row['midpoint']} discarded.",level="verbose")
                    continue

                # Append processed time, flux, and error to dict
                flux_data["time"].extend(time_range)
                flux_data["flux"].extend(proc_flux)
                flux_data["flux_err"].extend(proc_err)

            except Exception as e:
                vprint(f"Detrender failed: {e}",level="quiet")
                continue

    # Output sorted DataFrame with time, flux, and flux error
    return pd.DataFrame.from_dict(flux_data).sort_values(by="time").reset_index(drop=True)

def fit_eclipses(fnames,instrument,ecl_midpoints=[],ecl_bounds=([0.001,0.01],[0.99,2.0])):
    """
    Fit for eclipses using a Gaussian and return parameters from fitting.
    """

    ecl_data = [] # List of eclipse data dicts containing midpoint time, depth, and duration

    naive = len(ecl_midpoints) == 0

    # Get a "data buffer" around eclipses to ensure the entire eclipse is in the fitting data
    if naive:
        data_buffer = 1.5*ecl_bounds[1][1]
    else:
        data_buffer = 0.5
   
    # Loop through each data file
    for fname in fnames:

        # Get image data
        time, flux, err = load_flux_data(fname,instrument)

        if naive:
            # Find midpoints using find_peaks

            # Convert bounds to find_peaks units
            cadence = np.min(np.diff(time))  # days per data point
            flux_lvl = np.median(flux)
        
            min_width_pts = np.max([int(1/2.35*ecl_bounds[0][1] / cadence),1])
            max_width_pts = np.max([int(1/2.35*ecl_bounds[1][1] / cadence),1])
            min_height = abs(ecl_bounds[0][0] * flux_lvl)
            max_height = abs(ecl_bounds[0][1] * flux_lvl)

            # Find local minima using find_peaks
            peaks, props = find_peaks(-flux, prominence=(min_height,max_height), width=(min_width_pts, max_width_pts))
            midpoints = time[peaks]

            if len(midpoints) == 0:
                vprint(f"No eclipses found in {fname}, skipping file.", level="verbose")
                continue

        else:
            midpoints = ecl_midpoints

        # Loop through midpoints
        for midpoint in midpoints:

            # Get data around possible midpoints within the given eclipse bounds
            w_range = (time < (midpoint + data_buffer)) & (time > (midpoint - data_buffer))

            # Skip if there are too few points
            if np.sum(w_range) < 8:
                vprint(f"Too few points around possible midpoint {midpoint:.4f}, skipping.", level="verbose")
                continue

            try:
                popt,pcov = fit_ecl_gaussian(time[w_range], flux[w_range], err[w_range], midpoint,ecl_bounds)
                a, b, c, d, e, f, g = popt

                # Reject if high error on duration or depth
                perr = np.sqrt(np.diag(pcov))

                if perr[0] > 0.05*a:
                    vprint(f"Derived depth error on possible midpoint {midpoint:.4f} too high, skipping.",level="verbose")
                    continue

                if perr[2] > 0.05*c:
                    vprint(f"Derived duration error on possible midpoint {midpoint:.4f} too high, skipping.",level="verbose")
                    continue
                    
                # If eclipse has not been rejected, add to ecl_data
                ecl_dict = {"midpoint":b,"error":pcov[1,1],"duration":6*c,"type":"ecl"}
                ecl_data.append(ecl_dict)

            except Exception as e:
                vprint(f"Fit failed at midpoint {midpoint:.4f}: {e}",level="verbose")
                continue

    # If empty, output empty DataFrame (an error will be thrown outside of the method)
    if len(ecl_data) == 0:
        return pd.DataFrame.from_dict(ecl_data)
    
    # Output dataframe of eclipse times
    return pd.DataFrame.from_dict(ecl_data).drop_duplicates().sort_values(by="midpoint").reset_index(drop=True)

def ecl_gaussian(t,A,b,c,d,e,f,g):
    """
    Inverted Gaussian with flat base for fitting; A is eclipse depth, b is eclipse time, c is the standard deviation
    """

    t_start = b - 3*c
    t_end = b + 3*c
    
    w_in_eclipse = (t > t_start) & (t < t_end)
    
    ecl_gauss = np.ones(len(t))
    ecl_gauss[~w_in_eclipse] = func(t[~w_in_eclipse],d,e,f,g)
    ecl_gauss[w_in_eclipse] = func(t[w_in_eclipse],d,e,f,g)*gaussian(t[w_in_eclipse],A,b,c)

    return ecl_gauss

def fit_ecl_gaussian(time, flux, err, midpoint, ecl_bounds=([0.001,0.01],[0.99,2.0])):
    """
    Attempt to fit an inverted Gaussian + polynomial to a light curve
    Returns (popt, pcov) on success, raises on failure
    """
    
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", OptimizeWarning)

        try:
            popt_poly, _ = curve_fit(func, time, flux, maxfev=100000)
        except (RuntimeError, ValueError) as e:
            raise RuntimeError(f"Baseline poly fit failed: {e}")
    
        d, e, f, g = popt_poly
        depth = 1 - np.min(flux) / np.max(flux)
        width = 0.166 * np.mean([ecl_bounds[0][1], ecl_bounds[1][1]])
    
        bounds = (
            [ecl_bounds[0][0], time[0],  0.166*ecl_bounds[0][1], -np.inf, -np.inf, -np.inf, -np.inf],
            [ecl_bounds[1][0], time[-1], 0.166*ecl_bounds[1][1],  np.inf,  np.inf,  np.inf,  np.inf]
        )
        try:

            popt, pcov = curve_fit(
                ecl_gaussian, time, flux,
                p0=[depth, midpoint, width, d, e, f, g],
                bounds=bounds, maxfev=100000
            )
        except (RuntimeError, ValueError) as e:
            raise RuntimeError(f"Eclipse fit failed: {e}")
    return popt, pcov

def poly_detrender(time, flux, err, usable_time, usable_flux):
    n_pts = len(usable_time)

    if n_pts < 4:
        vprint(f"Too few out-of-eclipse points for fitting: {n_pts}",level="verbose")

    # Set degree based on number of points
    degree = 1 if n_pts < 8 else (2 if n_pts < 20 else 3)

    # Find best-fitting polynomial
    coeffs = np.polyfit(usable_time, usable_flux, degree)
    trend = np.polyval(coeffs, time)
    
    return flux/trend, err/trend

def remove_cosmic_rays(flux_data,ecl_data):
    """
    Remove possible cosmic rays; assumed to be any data 5 standard deviations above the out-of-eclipse median
    """
    
    w_midpoints, w_in_eclipse = make_in_eclipse_mask(flux_data["time"],ecl_data)

    out_eclipse_flux = flux_data["flux"][~w_in_eclipse]
    med = np.median(out_eclipse_flux)
    noise = np.std(out_eclipse_flux)

    w_cosmic = (flux_data["flux"] > (med + 5*noise))

    return flux_data[~w_cosmic]

def get_elc_eclipse_times(ecl_path):
    """
    Get list of eclipse/transit times from an ELC eclipse times output file in the user-selected directory.
    Takes a single argument, ecl_path, a string indicating the directory where the ELC eclipse times files
    are stored.
    Returns a pandas DataFrame with "midpoint" and "type" derived from the file, and "NaN" for the
    eclipse duration and error on a successful file read.
    On an unsuccessful read, returns an empty DataFrame.
    """
    
    ecl_data = []
        
    # Get all ELC "time" files
    ecl_files = list(pathlib.Path(ecl_path).glob("ELC*time.dat"))

    if len(ecl_files) > 0:
        # Loop through eclipse files
        
        for file in ecl_files:

            # Determine eclipse type based on file name
            ecl_file = np.genfromtxt(file)

            ecl_type = "ecl"
            if "prim" in file.name:
                ecl_type = "prim"
            elif "sec" in file.name:
                ecl_type = "sec"
            elif "3tran1" in file.name:
                ecl_type = "primtrans"
            elif "3tran2" in file.name:
                ecl_type = "sectrans"         

            # Read each line of file and get the eclipse midpoint
            for line in ecl_file:
                ecl_dict = {"midpoint":float(line[1]),"error":np.nan,"duration":np.nan,"type":ecl_type}
                ecl_data.append(ecl_dict)
        
        return pd.DataFrame.from_dict(ecl_data).sort_values(by="midpoint").reset_index(drop=True) # Return the DataFrame sorted by midpoint time
    return pd.DataFrame()  # Return an empty DataFrame if no files found


def fill_nans(ecl_data,ecl_data_naive):
    """
    Fill in NaN fields in ecl_data using ecl_data_naive 
    """

    new_ecl_data = {}
        
    for i, row in ecl_data.iterrows():
        # Find closest midpoint
        w_closest_midpoint = np.argmin(np.abs(ecl_data_naive["midpoint"] - row["midpoint"]))

        # Iterate through columns within row
        for column in row.index:
            # Save values
            if column not in new_ecl_data:
                # Create column in new_ecl_data if it doesn't already exist
                new_ecl_data[column] = []

            if not pd.isnull(row[column]):
                # If not NaN, append ecl_data entry to new dataFrame
                new_ecl_data[column].append(row[column])
            elif pd.isnull(row[column]) and not pd.isnull(ecl_data_naive.loc[w_closest_midpoint][column]):
                # If NaN in ecl_data but non-NaN in ecl_data_naive append ecl_data_naive to new dataFrame
                new_ecl_data[column].append(ecl_data_naive.loc[w_closest_midpoint][column])
            else:
                # Else, append NaN
                new_ecl_data[column].append(np.nan)

    return pd.DataFrame(new_ecl_data)

def cut_non_eclipse_data(flux_data,ecl_data):
    """
    Given flux and eclipse data, cut out all non-eclipse data
    """

    w_midpoints, w_in_eclipse = make_in_eclipse_mask(flux_data["time"],ecl_data)
    
    return flux_data[w_in_eclipse].reset_index(drop=True)

def write_ELCSC(fnames,target,instrument):
    """
    Use a list of MAST files to calculate what times short cadence data occur and create an ELCSC.inp file
    """
    
    times = []    

    for fname in fnames:
        # Get data from FITS file
        time, flux, err = load_flux_data(fname,instrument)
    
        # Append initial and final time of file
        times.append({"Initial Time":np.min(time),"Final Time":np.max(time)})

    # Put times into pandas DataFrame and output to TSV
    times_df = pd.DataFrame.from_dict(times).sort_values(by="Initial Time").drop_duplicates()
    times_df.to_csv(f'{target.lower().replace(" ","_")}/ELCSC.inp', sep='\t', header=False, index=False)
    vprint("Created ELCSC.inp",level="normal")

def write_ELCgap(flux_data,target):
    """
    Find "gaps" in the data in order to create an ELCgap.inp file.
    Takes a list of pathnames. Files should be formatted as time, flux, error with some kind of separator.
    """
    
    ELCgap_path = f'{target.lower().replace(" ", "_")}/ELCgap.inp'

    with open(ELCgap_path,"w") as file:
        for i in range(len(flux_data) - 1):
            if (flux_data["time"].iloc[i+1] - flux_data["time"].iloc[i]) > 1.0:
                write_string = f"{(flux_data['time'].iloc[i] + 0.1):14.10f}" + "\t" + f"{(flux_data['time'].iloc[i+1] - 0.1):14.10f}" + "\n"
                file.write(write_string)
    vprint("Created ELCgap.inp",level="normal")

def write_data(flux_data,instrument,target):
    """
    Write processed Kepler or TESS data from MAST to a file for use in ELC 
    """
    target = target.lower().replace(" ", "_")
    path_str = f"{target.lower()}/{target.lower()}_{instrument.lower()}_photo.dat"
    placeholder_df = flux_data.copy()
    placeholder_df.to_csv(path_str, sep='\t',columns=["time","flux","flux_err"], header=False, index=False)
    vprint(f"Saved {instrument} data to {path_str}",level="normal") 

def write_eclipses(flux_data,ecl_data,target):
    """
    Write eclipses to data files for more precise midpoint fitting using ELC
    """

    # Set up eclipse number dict
    ecl_nums = {}

    # Convert the target name to lower case, and replace spaces with underscores
    # for filename safety
    target = target.lower().replace(" ", "_")
    
    # Loop through all eclipses
    for i, row in ecl_data.iterrows():

        # Skip if duration is NaN
        if pd.isna(row["duration"]):
            vprint(f"Skipping eclipse at {row['midpoint']:.4f}: duration is NaN", level="verbose")
            continue

        # Get all data within 3 eclipse durations of the midpoint
        data_range = 3*row["duration"]
        w_range = (flux_data["time"] >= (row["midpoint"] - data_range)) & (flux_data["time"] <= (row["midpoint"] + data_range))
        time_range = flux_data["time"][w_range]

        # Skip if no eclipse data
        if len(time_range) == 0:
            continue
            
        flux_range = flux_data["flux"][w_range]
        err_range = flux_data["flux_err"][w_range]
    
        flag = row["type"]

        # Set eclipse number (arbitary, used to distinguish between eclipse files)
        if flag not in ecl_nums:
            ecl_nums[flag] = 0
        ecl_nums[flag] += 1

        # Path to write to
        path_str = f"{target}/ecl_data/{target}_{flag}{ecl_nums[flag]}.dat"

        # Collect data into DataFrame and write to TSV file 
        ecl_df = pd.DataFrame({
            "time": time_range,
            "flux": flux_range,
            "err": err_range
        })

        ecl_df.to_csv(path_str, sep="\t", header=False, index=False)
    vprint(f"Wrote eclipses to folder {target}/ecl_data",level="normal") 

def gen_eclipsetimes(ecl_data,target):
    """
    Generate eclipse time data files for use in ELC. 
    Takes in a pandas DataFrame ecl_data, which must contain columns 
    "midpoints","error","type", and a string target that indicates the 
    target name.
    Outputs a TSV with the eclipse number, eclipse midpoint, and error
    on the midpoint for each seperate eclipse type in ecl_data.
    Returns nothing.
    """

    # Convert the target name to lower case, and replace spaces with underscores
    # for filename safety
    target = target.lower().replace(" ", "_")

    # Loop through every unique eclipse type
    for type in np.unique(ecl_data["type"]):
        output_path = f'{target}/{target}_{type}time.dat'
        w_type =  ecl_data["type"] == type
        type_only = ecl_data[w_type].reset_index(drop=True)
       
        # Get period (based on minimum difference between eclipse midpoints) and
        # reference time (based on earliest eclipse) in order to calculate eclipse number
        diffs = type_only["midpoint"].sort_values().diff()
        period = diffs[diffs > 0].min()
        reftime = type_only["midpoint"].min()

        # Calculate eclipse numbers
        type_only["num"] = round(1 + (type_only["midpoint"] - reftime)/period)

        # Insert a blank column for making initial tab
        type_only.insert(0, 'blank', '')

        # Drop duplicate entries
        type_only = type_only.drop_duplicates()

        # Output to TSV
        type_only.to_csv(output_path, sep="\t",columns=["blank","num","midpoint","error"], header=False, index=False,na_rep="NaN")
        vprint(f"Wrote eclipse times to {output_path}",level="normal") 
     

def remove_bad_ecls(ecl_data, ecl_subset, flag="bad"):
    """
    Remove bad eclipses using a listlike of eclipses
    """

    subset_indices = [idx for idx, row in ecl_data.iterrows() for ecl in ecl_subset if abs(ecl - row["midpoint"]) < row["duration"]]

    if flag == "bad":
        return ecl_data.loc[~ecl_data.index.isin(subset_indices)].reset_index(drop=True)
    elif flag == "good":
        return ecl_data.loc[ecl_data.index.isin(subset_indices)].reset_index(drop=True)

def plot_lightcurve(flux_data, title, show=None, save=None):
    """
    Plot a lightcurve from flux_data DataFrame
    """

    time = flux_data["time"].values.copy().astype(float)
    flux = flux_data["flux"].values.copy().astype(float)
    
    # Insert NaN where gaps between consecutive points exceed 0.25 days
    gap_threshold = 0.25
    gap_indices = np.where(np.diff(time) > gap_threshold)[0] + 1
    time = np.insert(time, gap_indices, np.nan)
    flux = np.insert(flux, gap_indices, np.nan)

    plt.figure()
    plt.title(title)
    plt.xlabel("Time (BJD - 2,455,000)")
    plt.ylabel("Normalized Flux")
    plt.plot(time,flux, 'k-')
    if save:
        # Convert title in filename
        # Replace Latex formatting in filename
        # Replace spaces with underscores for file safety
        translation_table =  str.maketrans("","",r"$()")
        plot_filename= "images/" + title.translate(translation_table).lower().replace(" ","_") + ".png"
        plt.savefig(plot_filename)
        vprint(f"Saved image to {plot_filename}",level="normal")
    if show:
        plt.show()

VERBOSITY = "normal"  # default before run_pipeline sets it

def run_pipeline(target,instruments,ecl_times,ecl_bounds,ecl_path,bad_ecls,gen_ecl_times,gen_elcgap,gen_elcsc,gen_indiv_ecls,gen_data_file,save_plots,verbosity):

    global VERBOSITY
    VERBOSITY = verbosity

    # Download files from MAST
    vprint(f"Downloading all available {target} data from MAST for selected instruments...",level="normal") 
    get_mast_files(target,instruments,verbosity)

    gen_any = [gen_ecl_times,gen_elcgap,gen_elcsc,gen_indiv_ecls,gen_data_file,save_plots]
    # Create directories
    if any(gen_any):
        vprint("Checking for required directories...",level="normal") 
        create_directories(target,gen_indiv_ecls,save_plots)

    ecl_data = {}
    flux_data = {}
    plot_data = {}

    for instrument in instruments:

        if instrument.upper() == "KEPLER":
            # Get list of Kepler SC and LC filenames
            kepler_sc_names = list(get_fits_filenames(target,'Kepler','SC'))
            kepler_lc_names = list(get_fits_filenames(target,'Kepler','LC'))
            fnames = np.array(kepler_lc_names + kepler_sc_names)
        
            if len(fnames) == 0:
                vprint(f"Could not find {instrument} data for {target}. Skipping to next instrument.",level="quiet")
                continue


            # Get data reliability from Kepler files headers
            fnames = filter_data_rel(fnames)
        
            #write ELCSC
            if gen_elcsc:
                vprint("Writing ELCSC file...",level="normal")
                write_ELCSC(kepler_sc_names,target,instrument)
        
        elif instrument.upper() == "TESS":
            # Get list of TESS filenames
            fnames = get_fits_filenames(target,'TESS')
               
        if len(fnames) == 0:
            vprint(f"Could not find {instrument} data for {target}. Skipping to next instrument.",level="quiet")
            continue

        # If fitting for eclipse midpoint times from data
        if ecl_times == "1":
            vprint("Fitting for eclipse midpoint times...",level="normal") 
            ecl_data[instrument] = fit_eclipses(fnames,instrument,ecl_midpoints=[],ecl_bounds=ecl_bounds)
            if ecl_data[instrument].empty:
                vprint(f"No eclipses found in {instrument} data. Skipping to next instrument.",level="quiet")
                continue

        # If reading ELC eclipse times file
        elif ecl_times == "2":
            vprint("Finding eclipse midpoints from file...",level="normal")
            ecl_data[instrument] = get_elc_eclipse_times(ecl_path)
            if "midpoint" not in ecl_data[instrument]:
                raise ValueError("File read error. Make sure you are using an unedited eclipse times output file from ELC.")
            vprint("Filling in missing data...",level="normal")
            naive_ecl_data = fit_eclipses(fnames,instrument,ecl_midpoints=ecl_data[instrument]["midpoint"].values)
            if naive_ecl_data.empty:
                vprint(f"No eclipses found in {instrument} data. Skipping to next instrument.",level="quiet")
                continue
            ecl_data[instrument] = fill_nans(ecl_data[instrument],naive_ecl_data)
            if ecl_data[instrument].empty:
                vprint(f"No eclipses found in {instrument} data. Skipping to next instrument.",level="quiet")
                continue

        # Detrend
        vprint("Detrending...",level="normal")
        start = t.time()
        flux_data[instrument] = mdpt_detrend(fnames,ecl_data[instrument],instrument)
        end = t.time()

        # If no data returned, skip to next instrument
        if flux_data[instrument].empty:
            vprint(f"Detrending produced no usable data for {instrument} data. Skipping to next instrument.",level="quiet")
            continue
        vprint(f"Detrender finished in {end-start:.2f}s.",level="normal") 

        # Remove cosmic rays
        vprint("Removing cosmic rays...",level="normal") 
        flux_data[instrument] = remove_cosmic_rays(flux_data[instrument],ecl_data[instrument])
        vprint("Cosmic rays removed.",level="normal") 

        # If the user chose to remove any eclipses from the data
        if len(bad_ecls) > 0:
            vprint(f"Cutting out bad data for {instrument} data.",level="normal") 
            ecl_data[instrument] = remove_bad_ecls(ecl_data[instrument],bad_ecls,flag="bad")

        vprint(f"Cutting out non-eclipse data for {instrument} data.",level="normal") 
        # Cut out non-eclipse data
        flux_data[instrument] = cut_non_eclipse_data(flux_data[instrument],ecl_data[instrument])
        plot_data[instrument] = {}
        plot_data[instrument]["data"] = flux_data[instrument]
        plot_data[instrument]["title"] = fr"{target} Light Curve (${instrument}$ Data Only)"

        # Write instrument data to file
        if gen_data_file:
            vprint(f"Writing {instrument} data to file.",level="normal")
            write_data(flux_data[instrument],instrument,target)
        
    # Concat flux and eclipse data; raise error if not data is found or another error occurs
    try:
        ecl_data_total = pd.concat([ecl_data[instrument] for instrument in ecl_data]).reset_index(drop=True)
    except Exception as e:
        raise ValueError("No eclipse data found. Try again with different bounds.")
    if ecl_data_total.empty:
        raise ValueError("No eclipse data found. Try again with different bounds.")

    try:
        flux_data_total = pd.concat([flux_data[instrument] for instrument in flux_data]).reset_index(drop=True)
    except Exception as e:
        raise ValueError("Detrending produced no usable data.")
    if flux_data_total.empty:
        raise ValueError("Detrending produced no usable data.")

    # Write eclipses for more thorough fitting with ELC
    if gen_indiv_ecls:
        vprint("Writing indiviudal eclispes to files.",level="normal")
        write_eclipses(flux_data_total,ecl_data_total,target)

    # Write ELCgap.inp
    if gen_elcgap:
        vprint(f"Writing ELCgap file.",level="normal")
        write_ELCgap(flux_data_total,target)

    # generate eclipse times:
    if gen_ecl_times:
        vprint(f"Generating eclipse times...",level="normal")
        gen_eclipsetimes(ecl_data_total,target)
    vprint("Done!",level="quiet")

    return plot_data

def run_thread(root, run_button, target, instruments, ecl_times, ecl_path, ecl_bounds,
               bad_ecls, gen_ecl_times, gen_elcgap, gen_elcsc,
               gen_indiv_ecls, gen_data_file, save_plots, verbosity,
               plot_data_bool):
    vprint("Running tool...",level="quiet")
    try:
        plot_data = run_pipeline(
                    target=target,
                    instruments=instruments,
                    ecl_times=ecl_times,
                    ecl_path=ecl_path,
                    ecl_bounds=ecl_bounds,
                    bad_ecls=bad_ecls,
                    gen_ecl_times=gen_ecl_times,
                    gen_elcgap=gen_elcgap,
                    gen_elcsc=gen_elcsc,
                    gen_indiv_ecls=gen_indiv_ecls,
                    gen_data_file=gen_data_file,
                    save_plots=save_plots,
                    verbosity=verbosity
                )

        for instrument in plot_data:
            root.after(0, lambda i=instrument: plot_lightcurve(
                plot_data[i]["data"], plot_data[i]["title"], plot_data_bool,
                save_plots
            ))
                
        root.after(0, lambda: messagebox.showinfo("Success", "Done!"))
    except Exception as e:
        root.after(0, lambda err=str(e): messagebox.showerror("Error", err))
    finally:
        root.after(0, lambda: run_button.configure(state="normal"))
        root.after(0, lambda: vprint("",level="quiet"))

def launch_ui():
    """
    UI for tool
    """

    # Set up variables and title
    root = tk.Tk()
    root.title("Pat's ELC Helper Tool")

    target_var = tk.StringVar()
    kepler_var = tk.BooleanVar()
    tess_var = tk.BooleanVar()
    ecl_times_var = tk.StringVar()
    ecl_path_var = tk.StringVar()

    depth_var_lower = tk.StringVar()
    depth_var_upper = tk.StringVar()
    dur_var_lower = tk.StringVar()
    dur_var_upper = tk.StringVar()

    bad_ecls_var = tk.StringVar()
    gen_ecl_times_var = tk.BooleanVar()
    gen_elcgap_var = tk.BooleanVar()
    gen_elcsc_var = tk.BooleanVar()
    indiv_ecls_var = tk.BooleanVar()
    data_file_var = tk.BooleanVar()

    plot_data_bool_var = tk.BooleanVar()
    save_plots_var = tk.BooleanVar()

    elcsc_check = ttk.Checkbutton(root, text="ELCSC.inp File?", variable=gen_elcsc_var)
    elcsc_check.configure(state="disabled")

    depth_entry_lower = ttk.Entry(root, textvariable=depth_var_lower)
    depth_entry_upper = ttk.Entry(root, textvariable=depth_var_upper)
    dur_entry_lower = ttk.Entry(root, textvariable=dur_var_lower)
    dur_entry_upper = ttk.Entry(root, textvariable=dur_var_upper)
    ecl_path_entry = ttk.Entry(root, textvariable=ecl_path_var, width=40)

    depth_entry_lower.configure(state="disabled")
    depth_entry_upper.configure(state="disabled")
    dur_entry_lower.configure(state="disabled")
    dur_entry_upper.configure(state="disabled")
    ecl_path_entry.configure(state="disabled")

    verbosity_var = tk.StringVar()
    verbosity_var.set("normal")

    def ecltimesvar_set():
        """
        Method for handling which fields are enabled or disabled based on which eclipse midpoint method is selected
        """
        
        if ecl_times_var.get() == "2":
            depth_entry_lower.configure(state="disabled")
            depth_entry_upper.configure(state="disabled")
            dur_entry_lower.configure(state="disabled")
            dur_entry_upper.configure(state="disabled")
            ecl_path_entry.configure(state="normal")
            depth_var_lower.set("")
            depth_var_upper.set("")
            dur_var_lower.set("")
            dur_var_upper.set("")


        elif ecl_times_var.get() == "1":
            depth_entry_lower.configure(state="normal")
            depth_entry_upper.configure(state="normal")
            dur_entry_lower.configure(state="normal")
            dur_entry_upper.configure(state="normal")
            ecl_path_entry.configure(state="disabled")
            ecl_path_var.set("")


    def keplersc_set():
        """
        Method for handling if the ELSC file field is enabled or disabled based on if any Kepler data is being used 
        """
        if kepler_var.get():
            elcsc_check.configure(state="normal")
        else:
            gen_elcsc_var.set(False)
            elcsc_check.configure(state="disabled")

    # Organize rows and columns
    # Target name and data settings
    row = 0
    ttk.Label(root, text="Data Settings").grid(row=row, column=0, sticky="w")
    row += 1
    
    ttk.Label(root, text="Target Name:").grid(row=row, column=0, sticky="w")
    ttk.Entry(root, textvariable=target_var).grid(row=row, column=1)
    row += 1

    ttk.Label(root, text="Instrument Data to Use:").grid(row=row, column=0, sticky="w")
    ttk.Checkbutton(root, text="Kepler", variable=kepler_var,command=keplersc_set).grid(row=row, column=1, sticky="w")
    ttk.Checkbutton(root, text="TESS", variable=tess_var).grid(row=row, column=2, sticky="w")
    row += 2

    # Eclipse midpoint times settings
    ttk.Label(root, text="Find Eclipse Times from...").grid(row=row, column=0, sticky="w")
    values = {"Fitting Light Curve Data" : "1", 
            "Reading ELC Eclipse Times File" : "2", } 

    for (text, value) in values.items(): 
        ttk.Radiobutton(root, text = text, variable = ecl_times_var, 
            value = value,command=ecltimesvar_set).grid(row=row,column=int(value), sticky="w")
    row += 1

    ttk.Label(root, text="Eclipse Depth Min and Max (Fractional Flux):").grid(row=row, column=0, sticky="w")
    depth_entry_lower.grid(row=row, column=1)
    depth_entry_upper.grid(row=row, column=2)
    row += 1

    ttk.Label(root, text="Eclipse Duration Min and Max (Days):").grid(row=row, column=0, sticky="w")
    dur_entry_lower.grid(row=row, column=1)
    dur_entry_upper.grid(row=row, column=2)
    row += 1

    ttk.Label(root, text="Folder Containing Eclipse Data:").grid(row=row, column=0, sticky="w")
    ecl_path_entry.grid(row=row, column=1, columnspan=2)
    row += 2

    # Files to be generated
    ttk.Label(root, text="Generate Files:").grid(row=row, column=0, sticky="w")
    row += 1
    ttk.Checkbutton(root, text="Eclipse Time Files?", variable=gen_ecl_times_var).grid(row=row, column=0, sticky="w")
    ttk.Checkbutton(root, text="ELCgap.inp File?", variable=gen_elcgap_var).grid(row=row, column=1, sticky="w")
    elcsc_check.grid(row=row, column=2, sticky="w")
    row += 1

    ttk.Checkbutton(root, text="Individual Eclipses to Separate Files?", variable=indiv_ecls_var).grid(row=row, column=0, sticky="w")
    ttk.Checkbutton(root, text="Instrument Data to File?", variable=data_file_var).grid(row=row, column=1, sticky="w")
    row += 1

    # Remove bad eclipse times
    ttk.Label(root, text="Eclipse Times to Remove from Data:").grid(row=row, column=0, sticky="w")
    ttk.Entry(root, textvariable=bad_ecls_var).grid(row=row, column=1)
    row += 2
    
    # Plot settings
    ttk.Label(root, text="Plot Settings:").grid(row=row, column=0, sticky="w")
    row += 1
    ttk.Checkbutton(root, text="Show Final Light Curve Plots?", variable=plot_data_bool_var).grid(row=row, column=0, sticky="w")
    ttk.Checkbutton(root, text="Save Plots?", variable=save_plots_var).grid(row=row, column=1, sticky="w")
    row += 2

    # Verbosity level
    ttk.Label(root, text="Verbosity Level:").grid(row=row, column=0, sticky="w")
    row += 1
    verb_values = {"Quiet" : "quiet", 
            "Normal" : "normal",
            "Verbose" : "verbose"} 

    n = 0
    for (text, value) in verb_values.items(): 
        ttk.Radiobutton(root, text = text, variable = verbosity_var, 
            value = value).grid(row=row,column=n, sticky="w")
        n += 1
    row += 2


    def on_run():

        # Get instruments
        instruments = []
        if kepler_var.get():
            instruments.append("Kepler")
        if tess_var.get():
            instruments.append("TESS")

        ecl_bounds = ([np.nan, np.nan],[np.nan, np.nan])

        # Set eclipse depth bounds
        if depth_var_lower.get() and depth_var_upper.get():
            ecl_bounds[0][0] = float(depth_var_lower.get())
            ecl_bounds[1][0] = float(depth_var_upper.get())

        # Set eclipse duration bounds
        if dur_var_lower.get() and dur_var_upper.get():
            ecl_bounds[0][1] = float(dur_var_lower.get())
            ecl_bounds[1][1] = float(dur_var_upper.get())

        # Path to ELC eclipse data
        ecl_path = str(ecl_path_var.get()) if ecl_path_var.get() else ""

        # Bad Eclipse Times
        bad_ecls = bad_ecls_var.get().split()
        bad_ecls = [float(x) for x in bad_ecls] if bad_ecls else []

        # Read in other variables
        target=target_var.get().strip()
        ecl_times=ecl_times_var.get()
        gen_ecl_times=gen_ecl_times_var.get()
        gen_elcgap=gen_elcgap_var.get()
        gen_elcsc=gen_elcsc_var.get()
        gen_indiv_ecls=indiv_ecls_var.get()
        gen_data_file=data_file_var.get()
        plot_data_bool=plot_data_bool_var.get()
        save_plots=save_plots_var.get()
        verbosity=verbosity_var.get()

        try:
            # Error handling
            if target == "":
                raise ValueError("Enter a target name.")

            if len(instruments) == 0:
                raise ValueError("Select at least one instrument.")

            if not ecl_times:
                raise ValueError("Select a method of acquiring eclipse times.")

            if ecl_times == "1":
                if np.isnan(ecl_bounds).any():
                    raise ValueError("Empty eclipse depth and/or eclipse duration field(s).")

                if ecl_bounds[1][0] >= 1 or ecl_bounds[0][0] <= 0:
                    raise ValueError("Eclipse depth bounds must be less than 1 and greater than 0.")

                if ecl_bounds[0][1] <= 0:
                    raise ValueError("Minimum eclipse duration must be greater than zero.")
    
                if ecl_bounds[0][1] >= ecl_bounds[1][1]:
                    raise ValueError("Eclipse duration lower bound must be smaller than upper bound.")

                if ecl_bounds[0][0] >= ecl_bounds[1][0]:
                    raise ValueError("Eclipse depth lower bound must be smaller than upper bound.")

            if ecl_times == "2":
                if ecl_path == "":
                    raise ValueError("Enter the path to the folder containg ELC eclipse time data.")

                if not pathlib.Path(ecl_path).exists():
                    raise ValueError("Entered eclipse time data path is not valid. Enter a valid path to the folder containg ELC eclipse time data.")

                if not pathlib.Path(ecl_path).is_dir():
                    raise ValueError("Entered eclipse time data path is not a directory. Enter a valid path to the folder containg ELC eclipse time data.")

            if not verbosity:
                raise ValueError("Select a verbosity level.")

        except ValueError as e:
            messagebox.showerror("Error", str(e))
            return
            
        # Disable [Run] button while code is running
        run_button.configure(state="disabled")

        threading.Thread(target=run_thread,
                         args=(root, run_button, target, instruments, ecl_times, ecl_path, ecl_bounds,
                               bad_ecls, gen_ecl_times, gen_elcgap, gen_elcsc,
                               gen_indiv_ecls, gen_data_file, save_plots, verbosity,
                               plot_data_bool),
                          daemon=True).start()

    def on_close():
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

        # Delete each variable by name
        nonlocal target_var, kepler_var, tess_var, ecl_times_var, ecl_path_var
        nonlocal depth_var_lower, depth_var_upper, dur_var_lower, dur_var_upper
        nonlocal bad_ecls_var, gen_ecl_times_var, gen_elcgap_var, gen_elcsc_var
        nonlocal indiv_ecls_var, data_file_var, plot_data_bool_var, save_plots_var
        nonlocal verbosity_var

        del target_var, kepler_var, tess_var, ecl_times_var, ecl_path_var
        del depth_var_lower, depth_var_upper, dur_var_lower, dur_var_upper
        del bad_ecls_var, gen_ecl_times_var, gen_elcgap_var, gen_elcsc_var
        del indiv_ecls_var, data_file_var, plot_data_bool_var, save_plots_var
        del verbosity_var

        root.update()
        root.destroy()

    # Run and Exit buttons
    run_button = ttk.Button(root, text="Run", command=on_run)
    exit_button = ttk.Button(root, text="Exit", command=on_close)
    root.protocol("WM_DELETE_WINDOW", on_close)  # also handle the X button

    run_button.grid(row=row, column=0)
    exit_button.grid(row=row, column=2)
    row += 1

    console = ScrolledText(root, height=15, state="disabled")
    console.grid(row=row, column=0, columnspan=3, sticky="nsew")
    sys.stdout = ConsoleRedirect(console)
    sys.stderr = ConsoleRedirect(console)
    row += 1

    # Set row and column size in pixels
    col_count, row_count = root.grid_size()

    for col in range(col_count):
        root.grid_columnconfigure(col, minsize=250)

    for row in range(row_count):
        root.grid_rowconfigure(row, minsize=25)

    root.mainloop()

# Console handling to avoid excessive messaging
class ConsoleRedirect:
    def __init__(self, widget):
        self.widget = widget
    def write(self, message):
        if ("main thread is not in main loop" in message or (
            "Variable.__del__" in message or 
            "Image.__del__" in message or
            "self.tk.call" in message or
            "self._tk.call" in message
        )) and sys.platform == "win32":
            return  # suppress known CPython 3.13 Windows tkinter GC noise
        self.widget.after(0, self._append_text, message)

    def _append_text(self, message):
        self.widget.configure(state="normal")
        self.widget.insert("end", message)
        self.widget.see("end")
        self.widget.configure(state="disabled")

    def flush(self):
        pass

if __name__ == "__main__":
    # Filter warnings to avoid excessive messaging
    warnings.filterwarnings("ignore", message=".*RuntimeError:.*")
    warnings.filterwarnings("ignore", message=".*main thread is not in main loop.*")
    launch_ui()
