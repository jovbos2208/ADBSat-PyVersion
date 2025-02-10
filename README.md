# ADBSat-PyVersion

## About
ADBSat-PyVersion is a Python translation of the original [ADBSat](https://github.com/nhcrisp/ADBSat) project. ADBSat was originally developed to provide an analysis tool for spacecraft attitude dynamics and control. This repository aims to bring its functionality into Python for better accessibility and further development.

## Repository Structure

```
# Insert tree structure here
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

Please provide one of the heights `[100,150,200,250,300,350,400]` when executing the Python script. An index for a database entry is optional (e.g., for the comparison of models, etc.).

```sh
python test_example.py 200 4000
```

(Databases for more heights will be added soon.)

