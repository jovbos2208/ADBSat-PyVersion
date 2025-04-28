# ADBSat-PyVersion

## About
ADBSat-PyVersion is a Python translation of the original [ADBSat](https://github.com/nhcrisp/ADBSat) project. ADBSat was originally developed to provide an analysis tool for spacecraft attitude dynamics and control. This repository aims to bring its functionality into Python for better accessibility and further development.

## Repository Structure

```
.
├── atmos_data
│   ├── database_100km.csv
│   ├── database_110km.csv
│   ├── database_120km.csv
│   ├── database_130km.csv
│   ├── database_140km.csv
│   ├── database_150km.csv
│   ├── database_160km.csv
│   ├── database_170km.csv
│   ├── database_180km.csv
│   ├── database_190km.csv
│   ├── database_200km.csv
│   ├── database_210km.csv
│   ├── database_220km.csv
│   ├── database_230km.csv
│   ├── database_240km.csv
│   ├── database_250km.csv
│   ├── database_260km.csv
│   ├── database_270km.csv
│   ├── database_280km.csv
│   ├── database_290km.csv
│   ├── database_300km.csv
│   ├── database_310km.csv
│   ├── database_320km.csv
│   ├── database_330km.csv
│   ├── database_340km.csv
│   ├── database_350km.csv
│   ├── database_360km.csv
│   ├── database_370km.csv
│   ├── database_380km.csv
│   ├── database_390km.csv
│   ├── database_400km.csv
│   └── data_load.py
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
│   ├── __pycache__
│   │   ├── ADBSatConstants.cpython-312.pyc
│   │   ├── ADBSatConstants.cpython-313.pyc
│   │   ├── ADBSatImport.cpython-312.pyc
│   │   ├── ADBSatImport.cpython-313.pyc
│   │   ├── calc_coeff.cpython-312.pyc
│   │   ├── calc_coeff.cpython-313.pyc
│   │   ├── coeff_CLL.cpython-312.pyc
│   │   ├── coeff_CLL.cpython-313.pyc
│   │   ├── coeff_cook.cpython-312.pyc
│   │   ├── coeff_cook.cpython-313.pyc
│   │   ├── coeff_DRIA.cpython-312.pyc
│   │   ├── coeff_DRIA.cpython-313.pyc
│   │   ├── coeff_maxwell.cpython-312.pyc
│   │   ├── coeff_maxwell.cpython-313.pyc
│   │   ├── coeff_newton.cpython-312.pyc
│   │   ├── coeff_newton.cpython-313.pyc
│   │   ├── coeff_schaaf.cpython-312.pyc
│   │   ├── coeff_schaaf.cpython-313.pyc
│   │   ├── coeff_sentman.cpython-312.pyc
│   │   ├── coeff_sentman.cpython-313.pyc
│   │   ├── coeff_solar.cpython-312.pyc
│   │   ├── coeff_solar.cpython-313.pyc
│   │   ├── coeff_storchHyp.cpython-312.pyc
│   │   ├── coeff_storchHyp.cpython-313.pyc
│   │   ├── environment.cpython-312.pyc
│   │   ├── environment.cpython-313.pyc
│   │   ├── __init__.cpython-312.pyc
│   │   ├── __init__.cpython-313.pyc
│   │   ├── insidetri.cpython-312.pyc
│   │   ├── insidetri.cpython-313.pyc
│   │   ├── mainCoeff.cpython-312.pyc
│   │   ├── mainCoeff.cpython-313.pyc
│   │   ├── obj_fileTri2patch.cpython-312.pyc
│   │   ├── obj_fileTri2patch.cpython-313.pyc
│   │   ├── shadowAnaly.cpython-312.pyc
│   │   ├── shadowAnaly.cpython-313.pyc
│   │   ├── surfaceNormals.cpython-312.pyc
│   │   └── surfaceNormals.cpython-313.pyc
│   ├── shadowAnaly.py
│   └── surfaceNormals.py
├── inou
│   ├── models
│   │   ├── BM_aussen.mat
│   │   ├── BM_VAR1_1.mat
│   │   ├── Box.mat
│   │   ├── cube.mat
│   │   ├── Cube.mat
│   │   ├── CubeSat.mat
│   │   └── Opt_Sat.mat
│   ├── obj_files
│   │   ├── Cube.obj
│   │   ├── meshlab_reset_origin.mlx
│   │   └── stl2obj.py
│   ├── results
│   │   └── Cube
│   │       └── Cube_Cook_aoa0_aos90.0.mat
│   └── stl_files
│       └── Cube.STL
├── postpro
│   ├── plotNormals.py
│   ├── plot_surfq.py
│   └── __pycache__
│       ├── plot_surfq.cpython-312.pyc
│       └── plot_surfq.cpython-313.pyc
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




