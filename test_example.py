import os
import numpy as np
import pandas as pd
import sys
import time
from scipy.io import loadmat,savemat
from calc.ADBSatImport import ADBSatImport
from calc.environment import environment
from calc.calc_coeff import calc_coeff
from postpro.plot_surfq import plot_surfq
from calc.ADBSatConstants import ConstantsData

start = time.time()

# Eingabedaten
mod_name = 'Cube'
adbsat_path = '/home/jovan/software/ADBSat-PyVersion/'  # Basispfad anpassen
mod_in = os.path.join(adbsat_path, 'inou', 'obj_files', f"{mod_name}.obj")
mod_out = os.path.join(adbsat_path, 'inou', 'models')
res_out = os.path.join(adbsat_path, 'inou', 'results', mod_name)

verbose = False
delete_temp_files = False

# Importieren des Modells
mod_out = ADBSatImport(mod_in, mod_out,mod_name, verbose)
mesh = loadmat(mod_out)
N_elems = np.shape(mesh['meshdata']['XData'][0,0])[1]

if len(sys.argv) != 4:
    print("Usage: script.py <alt:int> <aos_deg:float> <idx:int>")
    sys.exit(1)

constants = ConstantsData()
# Argumente parsen
alt = int(sys.argv[1]) * 1e3  # Höhe in Metern
aos_deg = float(sys.argv[2])  # Seitenwinkel (AoS)
idx = int(sys.argv[3])        # Datenbank-Index
    
inc = 130  # Inklination in Grad
env = { "h": alt}

aoa_deg = 0  # Angle of attack in degrees
aos_deg = aos_deg+90  # Angle of sideslip in degrees

# Modellparameter
shadow = True
solar = True
inparam = {
    "gsi_model": 'Sentman',
    "alpha": np.ones(N_elems)*0.9,
    "alphaN": np.ones(N_elems),
    "sigmaN" : np.ones(N_elems),
    "sigmaT" : np.ones(N_elems),
    "Tw": 300,
    "sol_cR": 0.15,
    "sol_cD": 0.25
}
# Verbose und Cleanup


database = pd.read_csv(f"atmos_data/database_{int(sys.argv[1]):03d}km.csv")
# Umgebungseigenschaften berechnen
inparam = environment(inparam,database,idx,**env)

# Koeffizienten berechnen
file_out = calc_coeff(mod_out, res_out, [np.radians(aoa_deg)], [np.radians(aos_deg)], inparam, shadow, solar, delete_temp_files, verbose)

mesh_path = os.path.join(adbsat_path, 'inou', 'models',f'{mod_name}.mat')

end = time.time()

print(f'Py-ADBSat finished!! ({end-start} s)')

# Visualisieren der Ergebnisse
if verbose and not delete_temp_files:
    plot_surfq(file_out, mesh_path, aoa_deg, aos_deg, 'cd',show_normals=True)
    

