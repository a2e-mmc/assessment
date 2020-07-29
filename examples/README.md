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
  Pandas dataframe using `mmctools.dataloaders.read_dir()`

`examples/002-process_filelist.ipynb`

- Process a list of downloaded radar data files, using the radar profiler
  reader, combining the data into a single Pandas dataframe using 
  `mmctools.dataloaders.read_files()`

`examples/003-process_subdirs_quickplot.ipynb`

- Process subdirectories of downloaded sodar data files, combining the data into
  a single Pandas dataframe using `mmctools.dataloaders.read_date_dirs()`
- Quick time-height plot using `mmctools.plotting.plot_timeheight()`

`examples/004-process_met_data.ipynb`

- **Complete example** demonstrating DAP file downloads (using `dap-py`), data
  processing, and calculation of select atmospheric quantities.
- Met mast data are processed using `mmctools.dataloaders.read_files()` in
  conjunction with the `mmctools.measurements.metmast` submodule. The `metmast`
  submodule provides `read_data()` for reading time series (e.g., measured by
  sonic anemometers) at different heights within a single file.
- In this case, the metmast data format is _NOT_ standard, so the fields are
  defined by `metmast.RMYoung_05106`. This may be used as an example for other
  similar met data.
- Atmospheric calculations have been performed with different empirical models
  implemented within `mmctools.helper_functions` wherever applicable. 

## advanced data processing 

`datasets/SWiFT/process_TTU_tower.ipynb`

- Define a custom data reader
- Read and standardize data (data-set dependent calculations)
- Calculate standard statistical quantities, including variances and
  covariances using `mmctools.helper_functions.covariance()`
- Visualize data using `mmctools.plotting` module

