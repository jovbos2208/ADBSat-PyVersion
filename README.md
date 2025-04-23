# ADBSat-PyVersion

## About
ADBSat-PyVersion is a Python translation of the original [ADBSat](https://github.com/nhcrisp/ADBSat) project. ADBSat was originally developed to provide an analysis tool for spacecraft attitude dynamics and control. This repository aims to bring its functionality into Python for better accessibility and further development.

## Repository Structure

```
.
├── atmos_data
│   ├── database100km.zip
│   ├── database150km.zip
│   ├── database200km.zip
│   ├── database250km.zip
│   ├── database300km.zip
│   ├── database350km.zip
│   └── database400km.zip
├── calc
│   ├── ADBSatConstants.py
│   ├── ADBSatFcn.py
│   ├── ADBSatImport.py
│   ├── calc_coeff.py
│   ├── coeff_CLL.py
│   ├── coeff_cook.py
│   ├── coeff_DRIA.py
│   ├── coeff_maxwell.py
│   ├── coeff_newton.py
│   ├── coeff_schaaf.py
│   ├── coeff_sentman.py
│   ├── coeff_solar.py
│   ├── coeff_storchHyp.py
│   ├── environment.py
│   ├── importobjtri.py
│   ├── __init__.py
│   ├── insidetri.py
│   ├── mainCoeff.py
│   ├── obj_fileTri2patch.py
│   ├── shadowAnaly.py
│   └── surfaceNormals.py
├── inou
│   ├── models
│   ├── obj_files
│   │   ├── Cube.obj
│   │   ├── CubeSat.obj
│   │   ├── meshlab_reset_origin.mlx
│   │   └── stl2obj.py
│   ├── results
│   │   ├── Cube
│   │   └── CubeSat
│   └── stl_files
│       └── Cube.STL
├── postpro
│   ├── plotNormals.py
│   └── plot_surfq.py
├── README.md
└── test_example.py

```

## Installation Guide
To download the `ABDSat-PyVersion` repository and navigate into the directory, run the following command in your terminal:

```sh
git clone git@github.com:jovbos2208/ABDSat-PyVersion.git
cd ABDSat-PyVersion/atmos_data
unzip '*.zip'
cd ..
```

Now, a Virtual Environment must be created, and the following dependencies need to be installed:

```sh
pip install numpy scipy matplotlib mplcursors pandas
```

## Testing the Installation
To test the installation, run:

```sh
python test_example.py 250
```

If everything works well, a plot should pop up with a cube with colored sides.

## Customizing the Test Example
You can simply edit the `test_example.py` script. If you want to use your own objects, you have to insert the OBJ file of your object into `inou/obj-files` and change the `mod_name` accordingly. Please add a directory with the object's name in `inou/results` before running it for the first time with the new object.

Please provide one of the heights from 100km to 400km in 10km-steps when executing the Python script. An index for a database entry is optional (e.g., for the comparison of models, etc.).

```sh
python test_example.py 200 4000
```

## References
Sinpetru, L. A., Crisp, N. H., Mostaza-Prieto, D., Livadiotti, S., and Roberts, P. C. E. “ADBSat: Methodology of a Novel Panel Method Tool for Aerodynamic Analysis of Satellites.” Computer Physics Communications, Vol. 275, 2022, p. 108326. <https://doi.org/10.1016/j.cpc.2022.108326>.

Sinpetru, L. A., Crisp, N. H., Roberts, P. C. E., Sulliotti-Linner, V., Hanessian, V., Herdrich, G. H., Romano, F., Garcia-Almiñana, D., Rodríguez-Donaire, S., and Seminari, S. “ADBSat: Verification and Validation of a Novel Panel Method for Quick Aerodynamic Analysis of Satellites.” Computer Physics Communications, Vol. 275, 2022, p. 108327. <https://doi.org/10.1016/j.cpc.2022.108327>.

Mostaza-Prieto, D., Characterisation and Applications of Aerodynamic Torques on Satellites (PhD Thesis), The University of Manchester, 2017.

## Licence
This project is licensed under the GPL-3.0 License - see the LICENSE file for details.




