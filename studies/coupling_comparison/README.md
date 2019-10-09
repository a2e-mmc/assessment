# Flat terrain coupling comparison
Lead: Matthew CHurchfield

Set of Jupyter notebooks used for the coupling comparison study based on the SWiFT diurnal cycle case

## Preprocessing

### Internal
Notebooks to preprocess mesoscale data for internal forcing coupling methods, including:

- `wrf_to_sowfa`: Generate SOWFA readable internal forcing files based on WRF data
- `obs_to_sowfa`: Generate SOWFA readable internal forcing files based on observations

### Boundary
Notebooks to preprocess mescoscale data for boundary forcing coupling methods

## Postprocessing
Notebooks for individual postprocessing and data analysis, including:

- `analyse_virtual_tower_data_template`: Read and process data of virtual towers

## Comparison
Notebooks to compare performance of different coupling techniques, including:

- `assess_virtual_tower_data`: Assess virtual tower data of different microscale simulations
