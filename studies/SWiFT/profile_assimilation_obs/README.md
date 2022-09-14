# Demonstrate the use of the profile assimilation technique with observations
Lead: Eliot Quon

Jupyter notebooks and file with helper functions to produce figures of the paper

"Using observational data to drive large-eddy simulations of a diurnal cycle at the SWiFT site",
by Dries Allaerts, Eliot Quon, Caroline Draxl, Patrick Hawbecker, and Matthew Churchfield


The study is based on the following simulations:

- `internal_pat_obs_sowfa_noT`: Observations to SOWFA based on internal forcing with observed profile assimilation (no temperature forcing)
- `internal_pat_obs_sowfa_wrfT`: Observations to SOWFA based on internal forcing with observed profile assimilation (temperature forcing based on mesoscale profile assimilation)
- `internal_pat_obs_sowfa_wrfTadv`: Observations to SOWFA based on internal forcing with observed profile assimilation (temperature forcing based on mesoscale budget component for temperature)
- `internal_pat_wrf_sowfa`: WRF to SOWFA based on internal forcing with mesoscale profile assimilation


The numerical setup, microscale output and postprocessing is based on the coupling comparison study in `assessment/studies/coupling_comparison`. Hence, preprocessing of WRF and observational data
is based on the notebooks

- `assessment/studies/coupling_comparison/preprocessing/internal/wrf_to_sowfa`
- `assessment/studies/coupling_comparison/preprocessing/internal/obs_to_sowfa`

The output of 9 virtual towers are read and processed by the notebook

- `assessment/studies/coupling_comparison/postprocessing/analyse_virtual_tower_data_template`
