# ADBSat-PyVersion

## About
ADBSat-PyVersion is a Python translation of the original [ADBSat](https://github.com/nhcrisp/ADBSat) project. ADBSat was originally developed to provide an analysis tool for spacecraft attitude dynamics and control. This repository aims to bring its functionality into Python for better accessibility and further development.

## Repository Structure

```
.
├── atmos_data
│  ├── database_100km.csv
│  ├── database_110km.csv
│  ├── database_120km.csv
│  ├── database_130km.csv
│  ├── database_140km.csv
│  ├── database_150km.csv
│  ├── database_160km.csv
│  ├── database_170km.csv
│  ├── database_180km.csv
│  ├── database_190km.csv
│  ├── database_200km.csv
│  ├── database_210km.csv
│  ├── database_220km.csv
│  ├── database_230km.csv
│  ├── database_240km.csv
│  ├── database_250km.csv
│  ├── database_260km.csv
│  ├── database_270km.csv
│  ├── database_280km.csv
│  ├── database_290km.csv
│  ├── database_300km.csv
│  ├── database_310km.csv
│  ├── database_320km.csv
│  ├── database_330km.csv
│  ├── database_340km.csv
│  ├── database_350km.csv
│  ├── database_360km.csv
│  ├── database_370km.csv
│  ├── database_380km.csv
│  ├── database_390km.csv
│  ├── database_400km.csv
│  └── data_load.py

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
You can simply edit the `test_example.py` script. If you want to use your own objects, you have to insert the OBJ file of your object into `inou/obj-files` and change the `mod_name` accordingly. 

Please provide a betwenn 100km and 400km in 10km-step (10,110,120,130,...) when executing the Python script. An index for a database entry is optional (e.g., for the comparison of models, etc.).

```sh
python test_example.py 200 4000
```

(Databases for more heights will be added soon.)

