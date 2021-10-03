# CDIP Wave-Power Plots  #

This repository contains python code to create CDIP wave power plots. 

### Repository Description ###

CDIP Wave Power Product creates matrix plots of wave energy. These characterize the various combina
tions of Significant Wave Height (Hs) and Energy Period (Te) of a particular CDIP buoy over time. T
he frequency of occurrence is colored and labeled as a percent of the total distribution of wave me
trics Te and Hs. Low occurences (green/yellow) indicate sea states which are less prevalent than hi
gher occurences (orange/red).

### Requirements ###

* Python3
* Libraries: numpy, matplotlib, cdippy
* Plot server vs local
* CDIP NetCDF data

### Usage ###
Default usage (no input params) will create 3 output figures:
* Wave Power Matrix: Frequency of occurence
* Wave Power Matrix: Percent of total power
* Wave Power Matrix: Wave power numbers
`python wave_power.py`


### Command line
* Run code from command line
```bash
python wave_power.py
```

### Django
* Run code from web
```bash
cd powerplot
python manage.py runserver
```
* Check out on web at:
```
http://127.0.0.1:8000/wave_power
http://127.0.0.1:8000/wave_power/?stn=191
http://127.0.0.1:8000/wave_power/?stn=191&date=201702
http://127.0.0.1:8000/wave_power/?stn=191&date=201702&type=percent
```

### Build instructions
From the command line, run the command

```bash
git clone https://randy_bucciarelli@bitbucket.org/randy_bucciarelli/wave-power.git
```

Once you have the repository cloned, you can update it periodically to match master repo

```bash
git pull origin master
```

This repository includes an [environment file](./wave_power_env.yml) which you can use to set up your python environment. To install this environment, type the following

```bash
cd wave-power
conda env create -f wave_power_env.yml
```

This will create a new environment called `wave_power`. To activate this environment, type

```bash
source activate wave_power
```

You will need to install the CDIPPY library to access CDIP data
```bash
# Install cdippy
cd ../cdippy
python setup.py sdist
pip install dist/cdippy-*.tar.gz
```
When you are done working, deactive the environment by typing

```bash
source deactivate
```

