import bz2
import os
import time
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from astropy.io import fits
import numpy as np
import glob
from astroquery.vizier import Vizier
import sep
import pandas as pd
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy import units as u
from astropy.time import Time
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import logging
import re
import sys
import csv
from matplotlib.colors import Normalize
import argparse

# logging to save all the outputs
logging.basicConfig(filename='processing.log', level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger()

# Function to redirect stdout and stderr to the logger
class StreamToLogger:
    def __init__(self, logger, log_level=logging.INFO):
        self.logger = logger
        self.log_level = log_level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.log_level, line.rstrip())

    def flush(self):
        pass

sys.stdout = StreamToLogger(logger, logging.INFO)
sys.stderr = StreamToLogger(logger, logging.ERROR)

def create_circular_mask(h, w, center=None, radius=None):
    """
    Circular mask for masking out everything other than sky
    """
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return mask

def catalog_query(date_obs):
    """
    Set the longitude, latitude and height of TMO
    Calculate the Local Sidereal Time
    Fetch the catalog query
    Convert the ra and dec to alt-az 
    The width of the query can be changed. 78d is used to match the portion of visible sky
    """
    lon = -1 - 34/60 - 2.2/3600
    lat = 52 + 22/60 + 36.5/3600
    height = 84
    tmo = EarthLocation.from_geodetic(lon=lon*u.deg, lat=lat*u.deg, height=height*u.m)
    time = Time(date_obs, scale='utc', location=tmo)
    lst = time.sidereal_time('apparent')
    
    Vizier.ROW_LIMIT = -1
    MAG_LIMIT = 5
    res = Vizier.query_region(SkyCoord(lst.deg, lat, unit=(u.deg, u.deg), frame='icrs'),
                          width="78d", catalog="I/259/tyc2",
                          column_filters={'VTmag': f'<={MAG_LIMIT}'})
    cat= res[0]
    ra = np.array(cat["RA_ICRS_"])
    dec = np.array(cat["DE_ICRS_"])
    mag = np.array(cat["VTmag"])
    
    cat_radec = SkyCoord(ra, dec, unit=(u.deg, u.deg), frame='icrs')
    cat_altaz = cat_radec.transform_to(AltAz(obstime=date_obs,location=tmo))
    
    df_catalog = pd.DataFrame({
        'Magnitude': mag,
        'Cat_altitude': cat_altaz.alt.deg,
        'Cat_azimuth': cat_altaz.az.deg
    })
    
    return df_catalog

def save_results_to_csv(file_path, date_obs, mean_zero_point, zero_point_error, csv_file_path, zp_quad1, zp_quad2, zp_quad3, zp_quad4, 
                            zp_err_quad1, zp_err_quad2, zp_err_quad3, zp_err_quad4):
    """
    Save all the data to the csv file. 
    """
    header = ['filename', 'date_obs', 'zero_point', 'zero_point_error', 'zero_point_q1', 'zero_point_error_q1', 'zero_point_q2', 'zero_point_error_q2', 'zero_point_q3', 'zero_point_error_q3', 'zero_point_q4', 'zero_point_error_q4' ]
    
    # Collect the data for each file processed
    row = {
        'filename': os.path.basename(file_path),
        'date_obs': date_obs,
        'zero_point': mean_zero_point,
        'zero_point_error': zero_point_error,
        'zero_point_q1': zp_quad1,
        'zero_point_error_q1': zp_err_quad1,
        'zero_point_q2':zp_quad2,
        'zero_point_error_q2': zp_err_quad2,
        'zero_point_q3':zp_quad3,
        'zero_point_error_q3': zp_err_quad3,
        'zero_point_q4':zp_quad4,
        'zero_point_error_q4': zp_err_quad4
    }    
    # Check if file exists to write header or append
    file_exists = os.path.isfile(csv_file_path)
    
    # Write the row to the CSV
    with open(csv_file_path, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=header)
        
        # Write the header only if the file does not exist
        if not file_exists:
            writer.writeheader()
        
        # Write the data row
        writer.writerow(row)

def process_existing_fits_files(src_folder, zp_plot, csv_file_path, save_plot):
    """
    To start with existing files before the new files.
    """
    fits_files = glob.glob(os.path.join(src_folder, '*.fits')) + glob.glob(os.path.join(src_folder, '*.fits.bz2'))
    
    # Extract the datetime part from filenames and sort the files
    def get_datetime_from_filename(filename):
        match = re.search(r'(\d{4}_\d{2}_\d{2}__\d{2}_\d{2}_\d{2})', filename)
        if match:
            return pd.to_datetime(match.group(1), format='%Y_%m_%d__%H_%M_%S')
        return None

    fits_files = sorted(fits_files, key=get_datetime_from_filename)
    for file_path in fits_files:
        zp_plot= process_fits_file(file_path, zp_plot, csv_file_path, save_plot)
    return zp_plot

class FITSFileHandler(FileSystemEventHandler):
    """
    For real-time processing
    """
    def __init__(self, zp_plot, csv_file_path, save_plot):
        self.zp_plot = zp_plot
        self.csv_file_path= csv_file_path
        self.save_plot= save_plot
    
    def on_created(self, event):
        if not event.is_directory and (event.src_path.endswith('.fits') or event.src_path.endswith('.fits.bz2')):
            time.sleep(2)
            process_fits_file(event.src_path, self.zp_plot, self.csv_file_path, self.save_plot)  

def detect_sources(masked_data, h, w, mask_radius, gain, exposure_time):
    """
    Does source detect of objects.
    Can change the aperture size. Currently used: 4,6,8.
    Does gain calibration and normalisation at gain_target: 150
    Model from previous data: y= coeffs[0] * gain**2 + coeffs[1] * gain + coeffs[2] 
    coeffs=[0.59459824, 0.13382553, 0.27424068]
    Also normalises flues with exposure_time
    """
    try:
        # Background subtraction
        bkg = sep.Background(masked_data)
        masked_data_sub = masked_data - bkg

        sources = sep.extract(masked_data_sub, 3.0, err=bkg.globalrms)
        if sources is None or len(sources) == 0:
            logging.error("No sources detected.")
            return None

        cen_x = w / 2
        cen_y = h / 2
        radius = np.sqrt((sources['x'] - cen_x) ** 2 + (sources['y'] - cen_y) ** 2)

        sep_fluxes = []
        for i, source in enumerate(sources):
            try:
                x = np.array([source['x']])
                y = np.array([source['y']])

                flux, fluxerr, flag = sep.sum_circle(masked_data_sub, x, y, 4, bkgann=(6,8), err=bkg.globalrms, gain=gain)
                if flux is None:
                    logging.error(f"Flux computation returned None for source at index {i}.")
                else:
                    sep_fluxes.append(flux[0])

            except Exception as e:
                logging.error(f"Error calculating flux for source at index {i}: {e}")

        # Ensure sep_fluxes is not empty
        if not sep_fluxes:
            logging.error("No flux values were computed.")
            return None

        # Normalize flux based on gain using a tanh function
        gain_target = 150
        coeffs=[0.59459824, 0.13382553, 0.27424068]
        flux_target = coeffs[0] * gain_target**2 + coeffs[1] * gain_target + coeffs[2]
        flux_observed = coeffs[0] * gain**2 + coeffs[1] * gain + coeffs[2]
        flux_normalized = np.array(sep_fluxes) * (flux_target / flux_observed)
        
        npix = sources['npix']
        mask_condition = (radius < mask_radius) & (npix < 20)
        valid_indices = np.where(mask_condition)[0]

        df_sources = pd.DataFrame({
            'x': sources['x'][valid_indices],
            'y': sources['y'][valid_indices],
            'radius': radius[valid_indices],
            'flux_normalized': flux_normalized[valid_indices] / exposure_time
        })
        df_sources = df_sources.sort_values(by='flux_normalized', ascending=False).reset_index(drop=True)
        return df_sources

    except Exception as e:
        logging.error(f"Error detecting sources: {e}")
        return None
    
def divide_into_quadrants(df, center_x, center_y, x_col='x', y_col='y'):
    """
    Divides stars into 4 quadrants using the center x and y pixel positions.
    """
    quad1 = df[(df[x_col] >= center_x) & (df[y_col] >= center_y)]  #Top Right
    quad2 = df[(df[x_col] < center_x) & (df[y_col] >= center_y)]   #Top Left
    quad3 = df[(df[x_col] < center_x) & (df[y_col] < center_y)]    #Bottom Left
    quad4 = df[(df[x_col] >= center_x) & (df[y_col] < center_y)]   #Bottom Right
    return quad1, quad2, quad3, quad4

def calculate_seeing_ratios(matched_stars_quad, catalog_quad):
    """
    Does seeing ratio to detect clouds or unfavourable conditions.
    Ratio: Matched Stars/ Catalog star. 
    Basically: sources that are detected and cataloged matched divided by number of catalog sources present in the quadrant
    Ratio: 1 means good and 0 means bad.
    """
    if catalog_quad.empty:
        return np.nan

    matched_count = len(matched_stars_quad)
    catalog_count = len(catalog_quad)

    if catalog_count > 0:
        return matched_count / catalog_count
    else:
        return np.nan
    
def barghini_model(df_catalog):
    """
    Using the Barghini model to get x and y pixel position of catalog sources.
    5 parameters- R, F, a0, x_c, y_c.
    Model parameters calculated before are below:
    """
    # Model parameters
    R= 733.87070837
    F= 1.36494456
    a0= -21.97627716
    x_c= 777.31085116
    y_c= 511.41254342

    a_cat = df_catalog['Cat_azimuth']
    z_cat = 90 - df_catalog['Cat_altitude']  
    r_cat = R * np.sin(np.radians(z_cat/F))
    x_cat = x_c - r_cat * np.cos(np.radians(a_cat - a0))      
    y_cat = y_c + r_cat * np.sin(np.radians(a_cat - a0))  
    
    df_catalog['x_cat'] = x_cat
    df_catalog['y_cat'] = y_cat
    
    return df_catalog
    
def match_sources(df_sources, df_catalog, tolerance=5.0):
    """
    Match sources detected with catalog sources based on x and y positions.
    """
    matched_sources = []
    matched_catalog = []

    for i, (x, y) in df_sources[['x', 'y']].iterrows():
        distances = np.sqrt((df_catalog['x_cat'] - x) ** 2 + (df_catalog['y_cat'] - y) ** 2)
        min_index = distances.idxmin()
        
        # Check if the distance is within the tolerance
        if distances[min_index] <= tolerance:
            matched_sources.append(i)
            matched_catalog.append(min_index)

    # Create DataFrames for matched sources
    matched_df_sources = df_sources.loc[matched_sources].reset_index(drop=True)
    matched_df_catalog = df_catalog.loc[matched_catalog].reset_index(drop=True)
    matched_stars = pd.concat([matched_df_sources, matched_df_catalog], axis=1)

    return matched_stars

def calculate_zero_point(matched_stars):    
    """
    Calculates zero point
    """
    source_mags = matched_stars['Magnitude']
    safe_flux = np.clip(matched_stars['flux_normalized'], a_min=1e-10, a_max=None)
    inst_mags = -2.5 * np.log10(safe_flux)
    difference = source_mags - inst_mags
    zero_point = np.mean(difference)
    zero_point_std = np.std(difference)
    zero_point_error = zero_point_std / np.sqrt(len(difference))

    return zero_point, zero_point_error

def plot_results(zp_plot, save_plot, masked_data, seeing_ratios):
    """
    Plots all the 4 plots.
    """    
    plt.rcParams.update({
        'font.size': 25,           # Base font size for all text
        'axes.titlesize': 30,      # Title font size
        'axes.labelsize': 25,      # Axis labels font size
        'xtick.labelsize': 25,     # X-axis tick labels font size
        'ytick.labelsize': 25,     # Y-axis tick labels font size
        'legend.fontsize': 25,     # Legend font size
    })
    df = pd.DataFrame(zp_plot)
    df['datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], errors='coerce')
    df['zero_point'] = df['zero_point'].astype(float)
    df['zero_point_error'] = df['zero_point_error'].astype(float)
    
    past_hour_limit = df['datetime'].max() - pd.Timedelta(hours=1)
    df_past_hour = df[df['datetime'] >= past_hour_limit]
    
    # Create the plots
    if masked_data is None or seeing_ratios is None:
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(20, 24), gridspec_kw={'height_ratios': [1, 1, 1]})
    else:
        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(20, 32), gridspec_kw={'height_ratios': [1, 1, 1, 2]})
    
    # Top plot: Data from the past hour
    time_range = df['datetime'].max() - df['datetime'].min()

    # If the data covers more than an hour, apply the past hour limit
    if time_range > pd.Timedelta(hours=1):
        df_past_hour = df[df['datetime'] >= past_hour_limit]
    else:
        df_past_hour = df
    
    ax1.errorbar(df_past_hour['datetime'], df_past_hour['zero_point'], yerr=df_past_hour['zero_point_error'], 
             marker='o', linestyle='-', color='lightsteelblue', markerfacecolor='cornflowerblue', markeredgecolor='cornflowerblue',
             ecolor='lightgray')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Zero Point')
    ax1.set_title('Zero Point Over Time - Past Hour')
    ax1.grid(True)
    ax1.tick_params(axis='x', rotation=45)

    # Set the x-axis limits to start from the earliest data point and extend to the latest plus 5 minutes
    ax1.set_xlim([df_past_hour['datetime'].min() - pd.Timedelta(minutes=5), df['datetime'].max() + pd.Timedelta(minutes=5)])
    ax1.set_ylim(ax1.get_ylim()[::-1])

    if df_past_hour['zero_point'].isna().all():
        ax1.set_ylim(10, 14)  # Example default range
    else:
        ymin, ymax = np.nanmin(df_past_hour['zero_point']), np.nanmax(df_past_hour['zero_point'])
        ax1.set_ylim(ymax + 0.2, ymin - 0.2)

    # Format the x-axis
    locator1 = mdates.AutoDateLocator()
    formatter1 = mdates.DateFormatter('%d %H:%M')
    ax1.xaxis.set_major_locator(locator1)
    ax1.xaxis.set_major_formatter(formatter1)
    
    # Bottom plot: Full dataset with all data points including outliers
    ax2.errorbar(df['datetime'], df['zero_point'], yerr=df['zero_point_error'], 
                 marker='o', linestyle='-', color='lightsteelblue', markerfacecolor='steelblue', markeredgecolor='steelblue',
                 ecolor='lightgray')
    ax2.set_xlabel('Time')
    ax2.set_ylabel('Zero Point')
    ax2.set_title('Zero Point Over Time - Full Dataset')
    ax2.grid(True)
    ax2.tick_params(axis='x', rotation=45)
    ax2.set_ylim(ax2.get_ylim()[::-1])

    if df['zero_point'].isna().all():
        ax2.set_ylim(10, 14)  # Example default range
    else:
        ymin_full, ymax_full = np.nanmin(df['zero_point']), np.nanmax(df['zero_point'])
        ax2.set_ylim(ymax_full + 0.2, ymin_full - 0.2)


    locator2 = mdates.AutoDateLocator()
    formatter2 = mdates.DateFormatter('%d %H:%M')
    ax2.xaxis.set_major_locator(locator2)
    ax2.xaxis.set_major_formatter(formatter2)
    ax2.set_xlim([df['datetime'].min() - pd.Timedelta(minutes=5), df['datetime'].max() + pd.Timedelta(minutes=5)])
    
    # Third plot: Zero point data for each quadrant
    ax3.errorbar(df_past_hour['datetime'], df_past_hour['zero_point_q1'], yerr=df_past_hour['zero_point_error_q1'], 
             marker='o', linestyle='-', label='Quad 1', color='cadetblue', markerfacecolor='cadetblue', markeredgecolor='cadetblue',
             ecolor='lightgray')
    ax3.errorbar(df_past_hour['datetime'], df_past_hour['zero_point_q2'], yerr=df_past_hour['zero_point_error_q2'], 
             marker='o', linestyle='-', label='Quad 2', color='dodgerblue', markerfacecolor='dodgerblue', markeredgecolor='dodgerblue',
             ecolor='lightgray')
    ax3.errorbar(df_past_hour['datetime'], df_past_hour['zero_point_q3'], yerr=df_past_hour['zero_point_error_q3'], 
             marker='o', linestyle='-', label='Quad 3', color='palevioletred', markerfacecolor='palevioletred', markeredgecolor='palevioletred',
             ecolor='lightgray')
    ax3.errorbar(df_past_hour['datetime'], df_past_hour['zero_point_q4'], yerr=df_past_hour['zero_point_error_q4'], 
             marker='o', linestyle='-', label='Quad 4', color='goldenrod', markerfacecolor='goldenrod', markeredgecolor='goldenrod',
             ecolor='lightgray')
    ax3.set_xlabel('Time')
    ax3.set_ylabel('Zero Point')
    ax3.set_title('Zero Point Over Time - Quadrants - Past Hour')
    ax3.grid(True)
    ax3.tick_params(axis='x', rotation=45)
    ax3.set_xlim([df_past_hour['datetime'].min() - pd.Timedelta(minutes=5), df['datetime'].max() + pd.Timedelta(minutes=5)])
    ax3.legend(loc='upper left')
    ax3.set_ylim(ax3.get_ylim()[::-1])

    if df_past_hour[['zero_point_q1', 'zero_point_q2', 'zero_point_q3', 'zero_point_q4']].isna().all().all():
        ax3.set_ylim(19, 21)  # Example default range
    else:
        ymin_quad = np.nanmin(df_past_hour[['zero_point_q1', 'zero_point_q2', 'zero_point_q3', 'zero_point_q4']].min())
        ymax_quad = np.nanmax(df_past_hour[['zero_point_q1', 'zero_point_q2', 'zero_point_q3', 'zero_point_q4']].max())
        ax3.set_ylim(ymax_quad + 0.2, ymin_quad - 0.2)

    locator3 = mdates.AutoDateLocator()
    formatter3 = mdates.DateFormatter('%d %H:%M')
    ax3.xaxis.set_major_locator(locator3)
    ax3.xaxis.set_major_formatter(formatter3)
    
    # Fourth plot: Heatmap overlay on the masked image
    if masked_data is not None and seeing_ratios is not None:
        im = ax4.imshow(masked_data, cmap='gray')
    
        # Overlay heatmap based on seeing ratios
        overlay_colors = {
            'quad1': seeing_ratios['quad1'],
            'quad2': seeing_ratios['quad2'],
            'quad3': seeing_ratios['quad3'],
            'quad4': seeing_ratios['quad4'],
        }

        # Create an overlay image the same size as masked_data with transparency
        overlay = np.zeros((*masked_data.shape, 4))
    
        # Determine center and radius for the circular field
        h, w = masked_data.shape
        center_y, center_x = h // 2, w // 2
        radius = 450
    
        # Create a meshgrid to calculate distances from the center
        y, x = np.ogrid[:h, :w]
        distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)

        # Apply the overlays to each quadrant within the circular region
        mask_quad1 = (x >= center_x) & (y < center_y) & (distance_from_center <= radius)  # Top Right
        mask_quad2 = (x < center_x) & (y < center_y) & (distance_from_center <= radius)   # Top Left
        mask_quad3 = (x < center_x) & (y >= center_y) & (distance_from_center <= radius)  # Bottom Left
        mask_quad4 = (x >= center_x) & (y >= center_y) & (distance_from_center <= radius) # Bottom Right
    
        # Apply the heatmap values to the overlay
        overlay[mask_quad1, :3] = plt.cm.Spectral(overlay_colors['quad1'])[:3]
        overlay[mask_quad2, :3] = plt.cm.Spectral(overlay_colors['quad2'])[:3]
        overlay[mask_quad3, :3] = plt.cm.Spectral(overlay_colors['quad3'])[:3]
        overlay[mask_quad4, :3] = plt.cm.Spectral(overlay_colors['quad4'])[:3]
        overlay[mask_quad1 | mask_quad2 | mask_quad3 | mask_quad4, 3] = 0.4  # Set transparency
    
        # Show the overlay with transparency
        heatmap = ax4.imshow(overlay, extent=(0, w, 0, h), alpha=0.5, cmap='Spectral', norm=Normalize(vmin=0, vmax=1))
        
        # Add quadrant boundaries for visualization
        ax4.axvline(x=center_x, color='white', linestyle='--', linewidth=2)  # Vertical line at center_x
        ax4.axhline(y=center_y, color='white', linestyle='--', linewidth=2)  # Horizontal line at center_y
        ax4.text(center_x + w // 3, center_y - h // 3, f'Quad 4\nRatio: {overlay_colors["quad4"]:.2f}', color='white', fontsize=20, ha='center')  # Bottom Right
        ax4.text(center_x - w // 3, center_y - h // 3, f'Quad 3\nRatio: {overlay_colors["quad3"]:.2f}', color='white', fontsize=20, ha='center')  # Bottom Left
        ax4.text(center_x - w // 3, center_y + h // 3, f'Quad 2\nRatio: {overlay_colors["quad2"]:.2f}', color='white', fontsize=20, ha='center')  # Top Left
        ax4.text(center_x + w // 3, center_y + h // 3, f'Quad 1\nRatio: {overlay_colors["quad1"]:.2f}', color='white', fontsize=20, ha='center')  # Top Right

        
        ax4.set_title('Allsky Images with Seeing Ratio')
        # Add color bar next to the heatmap and ensure it matches the heatmap scale
        norm = Normalize(vmin=0, vmax=1)
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='Spectral'), ax=ax4, orientation='vertical', fraction=0.03, pad=0.08)
        cbar.set_label('Seeing Ratio')

    # Adjust layout to prevent overlap
    plt.tight_layout()
    fig.savefig(save_plot, format='jpeg', dpi=300)
    plt.close(fig)

def plot_full_dataset(zp_plot, save_plot):
    """
    Plots a big full dataset plot at the end of the night
    """
    plt.rcParams.update({
        'font.size': 25,
        'axes.titlesize': 30,
        'axes.labelsize': 25,
        'xtick.labelsize': 25,
        'ytick.labelsize': 25,
        'legend.fontsize': 25,
    })
    
    df = pd.DataFrame(zp_plot)
    df['datetime'] = pd.to_datetime(df['date'] + ' ' + df['time'], errors='coerce')
    df['zero_point'] = df['zero_point'].astype(float)
    df['zero_point_error'] = df['zero_point_error'].astype(float)

    fig, ax = plt.subplots(figsize=(25, 15))
    
    clipped_errors_full = df['zero_point_error'].clip(upper=2.0)  # Clip error bars

    ax.errorbar(df['datetime'], df['zero_point'], yerr=clipped_errors_full, 
                 marker='o', linestyle='-', color='lightsteelblue', markerfacecolor='steelblue', markeredgecolor='steelblue',
                 ecolor='lightgray')
    ax.set_xlabel('Time')
    ax.set_ylabel('Zero Point')
    ax.set_title('Zero Point Over Time - Full Dataset')
    ax.grid(True)
    ax.tick_params(axis='x', rotation=45)
    ax.set_ylim(ax.get_ylim()[::-1])
    
    locator = mdates.AutoDateLocator()
    formatter = mdates.DateFormatter('%d %H:%M')
    ax.xaxis.set_major_locator(locator)
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xlim([df['datetime'].min() - pd.Timedelta(minutes=5), df['datetime'].max() + pd.Timedelta(minutes=5)])

    plt.tight_layout()
    fig.savefig(save_plot, format='jpeg', dpi=300)
    plt.close(fig)

def process_fits_file(file_path, zp_plot, csv_file_path, save_plot):
    """
    Checks if the file is processed before in the csv file.
    If yes, skips the processing and stores the data from csv. The quadrant image plot is skipped.
    If no, does all the calculations and updates the plot.
    """
    try: 
        base_name = os.path.basename(file_path)
        
        # -----Check if the file is already processed before in the csv file----
        if os.path.exists(csv_file_path):
            df_csv = pd.read_csv(csv_file_path)
            if base_name in df_csv['filename'].values:
                logger.info(f"Skipping {base_name}, already processed.")

                row = df_csv[df_csv['filename'] == base_name].iloc[0]
                date_obs = row['date_obs']
                mean_zero_point = row['zero_point']
                zero_point_error = row['zero_point_error']
                zp_quad1= row['zero_point_q1']
                zp_quad2= row['zero_point_q2']
                zp_quad3= row['zero_point_q3']
                zp_quad4= row['zero_point_q4']
                zp_err_quad1= row['zero_point_error_q1']
                zp_err_quad2= row['zero_point_error_q2']
                zp_err_quad3= row['zero_point_error_q3']
                zp_err_quad4= row['zero_point_error_q4']
                
                # Append the zero point to the plot data
                zp_plot.append({
                    'date': date_obs.split('T')[0],
                    'time': date_obs.split('T')[1],
                    'zero_point': mean_zero_point,
                    'zero_point_error': zero_point_error,
                    'zero_point_q1': zp_quad1, 'zero_point_error_q1': zp_err_quad1,
                    'zero_point_q2': zp_quad2, 'zero_point_error_q2': zp_err_quad2,
                    'zero_point_q3': zp_quad3, 'zero_point_error_q3': zp_err_quad3,
                    'zero_point_q4': zp_quad4, 'zero_point_error_q4': zp_err_quad4
                })
                plot_results(zp_plot, save_plot, None, None)
                return zp_plot 
        
        # -----Process the file if not----
        if file_path.endswith('.bz2'):
            with bz2.open(file_path, 'rb') as f:
                with fits.open(f) as hdul:
                    data = hdul[0].data
                    data = data.byteswap().newbyteorder().astype(np.float32)
                    header = hdul[0].header
                    exposure_time = header['EXPOSURE']
                    date_obs = header['DATE-OBS']
                    gain= header['GAIN_ELE']          
        else:
            with fits.open(file_path) as hdul:
                data = hdul[0].data
                data = data.byteswap().newbyteorder().astype(np.float32)
                header = hdul[0].header
                exposure_time = header['EXPOSURE']
                date_obs = header['DATE-OBS']
                gain= header['GAIN_ELE']
                
        if gain < 50:
            logger.info(f"Low gain ({gain}) detected for {base_name}. Setting zero points to NaN.")
            zero_point = np.nan
            zero_point_error = np.nan
            zp_quad1 = np.nan
            zp_quad2 = np.nan
            zp_quad3 = np.nan
            zp_quad4 = np.nan
            zp_err_quad1 = np.nan
            zp_err_quad2 = np.nan
            zp_err_quad3 = np.nan
            zp_err_quad4 = np.nan
            
            h, w = data.shape
            mask_radius = 450
            mask = create_circular_mask(h, w, radius=mask_radius)
            masked_data = data.copy()
            masked_data[~mask] = 0
            seeing_ratios = {
                'quad1': np.nan,
                'quad2': np.nan,
                'quad3': np.nan,
                'quad4': np.nan
            }     
        else:
            # Masking data
            h, w = data.shape
            mask_radius = 450
            mask = create_circular_mask(h, w, radius=mask_radius)
            masked_data = data.copy()
            masked_data[~mask] = 0
                             
            # Barghini Model
            df_sources = detect_sources(masked_data, h, w, mask_radius, gain, exposure_time)
            df_catalog = catalog_query(date_obs)
            df_catalog = barghini_model(df_catalog)
            matched_stars = match_sources(df_sources, df_catalog)
            
            logging.debug(f"df_sources columns: {df_sources.columns}")
            logging.debug(f"df_catalog columns: {df_catalog.columns}")
            logging.debug(f"matched_stars columns: {matched_stars.columns}")
            
            # Quadrant Stuff
            center_x, center_y = w / 2, h / 2
            matched_quad1, matched_quad2, matched_quad3, matched_quad4 = divide_into_quadrants(matched_stars, center_x, center_y)
            catalog_quad1, catalog_quad2, catalog_quad3, catalog_quad4 = divide_into_quadrants(df_catalog, center_x, center_y, x_col='x_cat', y_col='y_cat')
            seeing_ratios = {
                'quad1': calculate_seeing_ratios(matched_quad1, catalog_quad1),
                'quad2': calculate_seeing_ratios(matched_quad2, catalog_quad2),
                'quad3': calculate_seeing_ratios(matched_quad3, catalog_quad3),
                'quad4': calculate_seeing_ratios(matched_quad4, catalog_quad4)
            }
            zp_quad1, zp_err_quad1 = calculate_zero_point(matched_quad1)
            zp_quad2, zp_err_quad2 = calculate_zero_point(matched_quad2)
            zp_quad3, zp_err_quad3 = calculate_zero_point(matched_quad3)
            zp_quad4, zp_err_quad4 = calculate_zero_point(matched_quad4)
            zero_point, zero_point_error = calculate_zero_point(matched_stars)
        
        # Append the zero point to the plot data (will be NaN if gain < 50)
        zp_plot.append({'date': date_obs.split('T')[0], 'time': date_obs.split('T')[1], 
                        'zero_point': zero_point, 'zero_point_error': zero_point_error,
                        'zero_point_q1': zp_quad1, 'zero_point_error_q1': zp_err_quad1,
                        'zero_point_q2': zp_quad2, 'zero_point_error_q2': zp_err_quad2,
                        'zero_point_q3': zp_quad3, 'zero_point_error_q3': zp_err_quad3,
                        'zero_point_q4': zp_quad4, 'zero_point_error_q4': zp_err_quad4})
        plot_results(zp_plot, save_plot, masked_data, seeing_ratios) 
        save_results_to_csv(file_path, date_obs, zero_point, zero_point_error, csv_file_path, 
                            zp_quad1, zp_quad2, zp_quad3, zp_quad4, 
                            zp_err_quad1, zp_err_quad2, zp_err_quad3, zp_err_quad4)
        logger.info(f"File {base_name} processed successfully.")
        return zp_plot
 
    except Exception as e:
        logging.error(f"Error processing file {file_path}: {e}")
    return zp_plot

def main(src_folder, csv_file_path, save_plot):
    zp_plot = []
    # Process existing FITS files
    zp_plot= process_existing_fits_files(src_folder, zp_plot, csv_file_path,save_plot)
    
    # Start real-time monitoring
    event_handler = FITSFileHandler( zp_plot, csv_file_path, save_plot)
    observer = Observer()
    observer.schedule(event_handler, src_folder, recursive=False) 
    observer.start()
    try:
        while True:
            time.sleep(3)
    except KeyboardInterrupt:
        observer.stop()
    observer.join()

    # Plot an end of the night final plot
    plot_full_dataset(zp_plot,save_plot)

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Process and monitor FITS files.')
    parser.add_argument('--src_folder', type=str, required=True, help='Path to the source folder containing FITS files.')
    parser.add_argument('--csv_file', type=str, required=True, help='Path to the CSV file.')
    parser.add_argument('--save_plot', type=str, required=True, help='Path to save the plot image.')

    # Parse arguments
    args = parser.parse_args()

    # Call main with parsed arguments
    main(args.src_folder, args.csv_file, args.save_plot)