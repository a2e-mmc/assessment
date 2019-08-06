# assessment/examples

Example notebooks to get started with assessment using the a2e-mmc
repositories. The following examples demonstrate some basic functionality of
the `mmctools` repository and provide a basic analysis framework. 

## data processing

`examples/001-DAP_data_processing.ipynb`

- Download data from the A2e [Data Archive & Portal (DAP)](https://a2e.energy.gov/about/dap)
  within a Python notebook using [`dap-py`](https://github.com/a2edap/dap-py)
- Read Wind Profiler radar data, and optionally, the scan properties using the
  `mmctools.measurements.radar.profiler()`
- Process a directory of downloaded radar data, combining the data into a single
  Pandas dataframe

`examples/002-process_filelist.ipynb`

- Process a list of downloaded radar data files, using the radar profiler
  reader, combining the data into a single Pandas dataframe

`examples/003-process_subdirs_quickplot.ipynb`

- Process subdirectories of downloaded sodar data files, combining the data into
  a single Pandas dataframe
- Quick time-height plot

## advanced data processing 

`datasets/SWiFT/process_TTU_tower.ipynb`

- Define a custom data reader
- Read and standardize data (data-set dependent calculations)
- Calculate standard statistical quantities, including variances and
  covariances using `mmctools.helper_functions.covariance()`
- Visualize data using `mmctools.plotting` module

