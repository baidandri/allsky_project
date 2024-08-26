# Night Sky Monitoring Using Allsky Camera

## Overview
The Marsh Observatory at the University of Warwick houses an all-sky camera that monitors the night sky. This project includes a Python script that processes real-time data from the camera. The script performs the following tasks:

- **Zero-Point Calculation**: Computes the zero-point for each image file.
- **Zero-Point Plot**: Updates the plot of zero-point over time
- **Quadrant Analysis**: Divides the sky into four quadrants to calculate and plot the zero-point.
- **Sky Visualization**: Creates a visual representation of the night sky.
- **Seeing Ratio**: Calculates the ratio of detected stars to catalog stars for each quadrant and overlays a heatmap on the night sky.
- At the end of the night, it produces a final plot of zero-point over time.

## Installation
To set up the project and install all required dependencies, follow these steps:

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/baidandri/allsky_project.git
   cd allsky_project

2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt

### Command-Line Usage
Here is the syntax for running the script. You need to specify the following:

- **`src_folder`**: The path to the folder where all the images are located.
- **`csv_file`**: The path where the script will save the output CSV file.
- **`save_plot`**: The path where the script will save the output plot (image file).

The script will automatically create the specified CSV file and plot image if they do not already exist.

   ```bash
   python zero_point.py --src_folder <path_to_source_folder> --csv_file <path_to_csv_file> --save_plot <path_to_plot_image>
   ```

### Example Usage
For example, for a folder named `allsky_2024-08-23` with source images either in `.fits` or `.fits.bz2` format, and you want to save the results to a CSV file named `zp_results.csv` and the plot to `plot_results.jpeg`, you would run:

   ```bash
   python zero_point.py --src_folder allsky_2024-08-23 --csv_file zp_results.csv --save_plot plot_results.jpeg
   ```


## Important Notes

1. **Projection Model**:
   - The code utilizes the projection model described in the paper: *"Astrometric calibration for all-sky cameras revisited"* by D. Barghini, D. Gardiol, et al., published in *Astronomy & Astrophysics*, Vol. 626, A105 (2019).
   - You can access the paper here: [https://doi.org/10.1051/0004-6361/201935580](https://doi.org/10.1051/0004-6361/201935580).
   - This model with calculated parameters is implemented in the script in the function called `barghini_model` at `260`.

2. **Gain Normalisation**:
   - The camera gain changes constantly throughout the night to accomodate the sky brightness.
   - The measured fluxes are normalised using previous data. The polynomial model used: flux= coeffs[0] * gain**2 + coeffs[1] * gain + coeffs[2] with coeffs=[0.59459824, 0.13382553, 0.27424068]. 
   - This model is implemented in the script in the function called `detect_sources` from `209` to `213`.

3. **Exposure Time Normalisation**:
   - The flux is also normalised by dividing it with exposure time on line `223`.

4. **Catalog Query**:
   - The tycho-2 catalog is used and the width of the sky query is 78d as seen in function `catalog_query`. The magnitude limit is set as 5 at line `77`.
   - The location of the observatory with latitude, longitude and height above sea level is also defined in line from `69` to `71`.

5. **Sleep Timers**:
   - Two sleep timers are used: one at line 162 for 2 seconds to ensure the new file is completely written before processing, and another at line 690 for 3 seconds to reduce CPU usage by delaying each loop iteration.

6. **Processing File**:
   - A log file named `processing.log` is defined on line `25` to track processing steps.

7. **Stopping Real-Time Processing**:
   - To stop real-time processing, you can use a keyboard interrupt (Ctrl+C).
   - If you wish to disable the real-time file check, comment out lines `684` to `693`.

8. **Masking Data**:
   - A circular mask is created to mask the low horizon and trees in the image.
   - A radius of 450 is used in lines `457`,`616`, and `629`.

