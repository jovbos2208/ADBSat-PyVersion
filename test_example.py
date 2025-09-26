# Usage: python test_example.py <alt[km]> <aos[deg]> <index>
# Example: python test_example.py 400 15 100

import os, sys, time
import numpy as np
import pandas as pd
from scipy.io import loadmat
from calc.ADBSatImport import ADBSatImport
from calc.environment import environment
from calc.calc_coeff import calc_coeff
from postpro.plot_surfq import plot_surfq
from calc.ADBSatConstants import ConstantsData

start = time.time()
mod_name = 'CHAMP'
# Use the script's own path to build absolute paths
base_path = os.path.dirname(os.path.abspath(__file__))

mod_in = os.path.join(base_path, 'inou/obj_files', f"{mod_name}.obj")
mod_out_dir = os.path.join(base_path, 'inou/models')
res_out_dir = os.path.join(base_path, 'inou/results', mod_name)
mesh_path = os.path.join(mod_out_dir, f"{mod_name}.mat")

# Ensure output directories exist
os.makedirs(mod_out_dir, exist_ok=True)
os.makedirs(res_out_dir, exist_ok=True)

# Read command line arguments
if len(sys.argv) != 4:
    print("Usage: script.py <alt[km]> <aos[deg]> <db_index>")
    sys.exit(1)

alt = int(sys.argv[1]) * 1e3
aos_deg = float(sys.argv[2])
idx = int(sys.argv[3])
aoa_deg = 0  # fixed

# Note on wind direction:
# AoA = 0°, AoS = 0° → Freestream (wind vector) points exactly in the -X direction in the body frame.
# (i.e., the satellite is flying in the +X direction; the incoming flow comes from -X)
# Therefore, surfaces with normals pointing in +X experience maximum dynamic pressure.

# Load and prepare geometry
mod_mat_path = ADBSatImport(mod_in, mod_out_dir, mod_name, verbose=False)
N_elems = np.shape(loadmat(mod_mat_path)['meshdata']['XData'][0, 0])[1]

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

# Load atmospheric conditions
database_path = os.path.join(base_path, f"atmos_data/database_{int(alt/1000):03d}km.csv")
database = pd.read_csv(database_path)
inparam = environment(inparam, database, idx, h=alt)

# Run coefficient calculation
print(f'Altitude: {alt*1e-3} km')
print('----------------------')
file_out = calc_coeff(mod_mat_path, res_out_dir, [np.radians(aoa_deg)], [np.radians(aos_deg)],
                      inparam, flag_shad=True, flag_sol=True, dp=False, delete_temp_files=False, verbose=False)
print('----------------------')
print(f"✅ Py-ADBSat finished in {time.time() - start:.2f} s")

# Visualization
plot_surfq(file_out, mesh_path, aoa_deg, aos_deg, 'cd', show_normals=False)
