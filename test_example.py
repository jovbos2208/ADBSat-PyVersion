import os
import numpy as np
import pandas as pd
import sys
from scipy.io import loadmat,savemat
from calc.ADBSatImport import ADBSatImport
from calc.environment import environment
from calc.calc_coeff import calc_coeff
from postpro.plot_surfq import plot_surfq
from calc.ADBSatConstants import ConstantsData

# Eingabedaten
mod_name = 'Cube'
adbsat_path = '/Path/to/ADBSat-PyVersion-main/'  # Basispfad anpassen
mod_in = os.path.join(adbsat_path, 'inou', 'obj_files', f"{mod_name}.obj")
mod_out = os.path.join(adbsat_path, 'inou', 'models')
res_out = os.path.join(adbsat_path, 'inou', 'results', mod_name)

verbose = True
delete_temp_files = False

# Importieren des Modells
mod_out = ADBSatImport(mod_in, mod_out,mod_name, verbose)
mesh = loadmat(mod_out)
N_elems = np.shape(mesh['meshdata']['XData'][0,0])[1]

# Bedingungen
constants = ConstantsData()
alt = sys.argv[1] * 1e3  # Höhe in Metern
if len(sys.argv)==3:
    idx = sys.argv[2]
else:
    idx=np.random.uniform(0,190000)
    
inc = 130  # Inklination in Grad
env = { "h": alt}

aoa_deg = 0  # Angle of attack in degrees
aos_deg = 90  # Angle of sideslip in degrees

# Modellparameter
shadow = True
solar = True
inparam = {
    "gsi_model": 'Maxwell',
    "alpha": np.ones(N_elems),
    "alphaN": np.ones(N_elems),
    "sigmaN" : np.ones(N_elems),
    "sigmaT" : np.ones(N_elems),
    "Tw": 300,
    "sol_cR": 0.15,
    "sol_cD": 0.25
}
# Verbose und Cleanup


database = pd.read_csv(f"atmos_data/msise00_database_{alt/1e3:03d}km.csv")
# Umgebungseigenschaften berechnen
inparam = environment(inparam,database,idx,**env)

print('Dynamic Pressure: \n')
print(0.5*np.sum(inparam['rho'])*inparam['mmean']/1000/constants.NA*np.sqrt(inparam['vinf']**2+inparam['vth']**2)**2)

# Koeffizienten berechnen
file_out = calc_coeff(mod_out, res_out, [np.radians(aoa_deg)], [np.radians(aos_deg)], inparam, shadow, solar, delete_temp_files, verbose)

mesh_path = os.path.join(adbsat_path, 'inou', 'models',f'{mod_name}.mat')
# Visualisieren der Ergebnisse
if verbose and not delete_temp_files:
    plot_surfq(file_out, mesh_path, aoa_deg, aos_deg, 'cd')
