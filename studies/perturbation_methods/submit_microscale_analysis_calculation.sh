#!/bin/bash
#PBS -A NSAP0003    
#PBS -N Micro_anl
#PBS -j oe
#PBS -m abe
#PBS -M hawbecke@ucar.edu
#PBS -q share     
#PBS -l walltime=05:00:00
### Request 10 CPUS for 10 threads
#PBS -l select=1:ncpus=4:mpiprocs=4:mem=200GB

simname='WPert'

# Activate your conda environment
source /glade/u/home/hawbecke/.bashrc
conda activate pyhawbeck

# Set python executable path:
TMP_PYTHON_PATH=$(which python3.7)

./convert_notebook_to_py.sh calc_only

### Run notebook as .py
#/glade/u/home/hawbecke/local/envs/pyhawbeck/bin/python3.7 -u PMIC-Microscale-analysis_template_tslist.py &> log.micro_anl 
$TMP_PYTHON_PATH -u PMIC-Microscale-analysis_template_tslist.py &> log.micro_anl_$simname 
