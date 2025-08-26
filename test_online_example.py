# Usage: python test_online_example.py <alt_km> <lat> <lon> <aos_deg> <aoa_deg>
# Example: python test_online_example.py 400 0 0 15 5

import os
import sys
import time
import numpy as np
import datetime
from scipy.io import loadmat
from calc.ADBSatImport import ADBSatImport
from calc.pymsis_data import get_atmospheric_data
from calc.calc_coeff import calc_coeff
from postpro.plot_surfq import plot_surfq

# --- Argument Parsing ---
if len(sys.argv) != 6:
    print("Usage: python test_online_example.py <alt_km> <lat_deg> <lon_deg> <aos_deg> <aoa_deg>")
    print("Example: python test_online_example.py 400 0 0 15 5")
    sys.exit(1)

start = time.time()
mod_name = 'Cube'

# --- Configuration ---
base_path = os.path.dirname(os.path.abspath(__file__))
mod_in = os.path.join(base_path, 'inou/obj_files', f"{mod_name}.obj")
mod_out_dir = os.path.join(base_path, 'inou/models')
res_out_dir = os.path.join(base_path, 'inou/results', mod_name)
mesh_path = os.path.join(mod_out_dir, f"{mod_name}.mat")

os.makedirs(mod_out_dir, exist_ok=True)
os.makedirs(res_out_dir, exist_ok=True)

# --- Simulation Parameters from Command Line ---
alt_km = float(sys.argv[1])
lat = float(sys.argv[2])
lon = float(sys.argv[3])
aos_deg = float(sys.argv[4])
aoa_deg = float(sys.argv[5])

# --- Fixed Atmospheric Parameters for pymsis ---
date = datetime.datetime(2010, 3, 15, 12, 0, 0)
f107 = 80.0
f107a = 75.0
ap = 7

print(f"Running online test for {mod_name} with command line parameters:")
print(f"Altitude: {alt_km} km, Latitude: {lat}°, Longitude: {lon}°")
print(f"AoS: {aos_deg}°, AoA: {aoa_deg}°")
print(f"Fixed Date: {date}, F10.7: {f107}, AP: {ap}")

# 1. Load and prepare geometry
# This converts the .obj file to a .mat file used by the calculation engine
mod_mat_path = ADBSatImport(mod_in, mod_out_dir, mod_name, verbose=False)
N_elems = np.shape(loadmat(mod_mat_path)['meshdata']['XData'][0, 0])[1]

# 2. Set initial GSI and SRP model parameters
inparam = {
    "gsi_model": 'Sentman',
    "alpha": np.ones(N_elems) * 0.9,
    "alphaN": np.ones(N_elems),
    "sigmaN": np.ones(N_elems),
    "sigmaT": np.ones(N_elems),
    "Tw": 300,
    "sol_cR": 0.15,
    "sol_cD": 0.25
}

# 3. Get atmospheric conditions using the new online function
inparam = get_atmospheric_data(inparam, date, lon, lat, alt_km, f107, f107a, ap)

print(f"Calculated Tinf: {inparam['Tinf']:.2f} K, Vinf: {inparam['vinf']:.2f} m/s, Speed Ratio s: {inparam['s']:.2f}")

# 4. Run the coefficient calculation
# Note: aos is adjusted by +90 degrees to match the coordinate system expectations if needed
file_out = calc_coeff(
    mod_mat_path, res_out_dir, 
    [np.radians(aoa_deg)], [np.radians(aos_deg + 90)],
    inparam, 
    flag_shad=True, 
    flag_sol=True, 
    dp=False, 
    delete_temp_files=False, 
    verbose=True
)

print(f"\n✅ Py-ADBSat (Online) finished in {time.time() - start:.2f} s")
print(f"Results saved to: {file_out}")

# 5. Visualization
print("\nGenerating plot...")
plot_surfq(file_out, mesh_path, aoa_deg, aos_deg, 'cd', show_normals=True)
