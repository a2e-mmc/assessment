# Flat-terrain coupling comparison

## Preprocessing

### Internal
Notebooks to preprocess mesoscale data for internal forcing coupling methods, including:

- `wrf_to_sowfa`: Generate SOWFA readable internal forcing files based on WRF data
- `obs_to_sowfa`: Generate SOWFA readable internal forcing files based on observations


## Simulations

### WRF/SOWFA
- `internal_bcc_wrf_sowfa`: WRF to SOWFA based on internal forcing with mesoscale budget components (Dries)
- `internal_pat_wrf_sowfa`: WRF to SOWFA based on internal forcing with mesoscale profile assimilation (Dries)
- `internal_pat_obs_sowfa_noT`: Observations to SOWFA based on internal forcing with observed profile assimilation (no temperature forcing) (Dries)
- ...


## Postprocessing
Notebooks for individual postprocessing and data analysis, including:

- `analyse_virtual_tower_data_template`: Read and process data of virtual towers
- `virtual_tower_assessment`: Assess virtual tower data of different microscale simulations
