# Demonstrate the use of the profile assimilation technique with observations
Lead: Dries Allaerts

Jupyter notebooks and file with helper functions to produce figures of the paper

"Using observational mean-flow data to drive large-eddy simulations of a diurnal cycle at the SWiFT site",
by Dries Allaerts, Eliot Quon, and Matthew Churchfield (submitted to Wind Energy on September 19, 2022)


The study is based on the following simulations:
| Simulation name | Momentum forcing | Virtual potential temperature forcing |
|--------------|------------------|---------------------------------------|
| internal.ipa.obs.noT | IPA w/ data from tower + radar | no forcing |
| internal.ipa.obs.Tassim | IPA w/ data from tower + radar | forcing at single height |
| internal.ipa.obs.allT | IPA w/ data from tower + radar | IPA w/ data from tower + RASS + sounding |
| internal.dpa.obs.noT | DPA w/ data from tower + radar | no forcing |
| internal.dpa.obs.Tassim | DPA w/ data from tower + radar | forcing at single height |
| internal.dpa.obs.allT | DPA w/ data from tower + radar | DPA w/ data from tower + RASS + sounding |
| internal.ipa.obs.wrfT | IPA w/ data from tower + radar | IPA w/ data from WRF |
| internal.ipa.obs.wrfTadv | IPA w/ data from tower + radar | BCC w/ data from WRF |
| internal.ipa.wrf | IPA w/ data from WRF | IPA w/ data from WRF |

The simulation setup files can be found [here](https://github.com/a2e-mmc/SOWFA-setups/tree/master/SWiFT).

The numerical setup and microscale output is based on the coupling comparison study in `assessment/studies/SWiFT/coupling_comparison`. Hence, preprocessing of WRF and observational data is based on the notebooks

- [assessment/studies/SWiFT/coupling_comparison/preprocessing/internal/wrf_to_sowfa](https://github.com/a2e-mmc/assessment/blob/master/studies/SWiFT/coupling_comparison/preprocessing/internal/wrf_to_sowfa.ipynb)
- [assessment/studies/SWiFT/coupling_comparison/preprocessing/internal/obs_to_sowfa](https://github.com/a2e-mmc/assessment/blob/master/studies/SWiFT/coupling_comparison/preprocessing/internal/obs_to_sowfa.ipynb)

Postprocessing of SOWFA output data is based on the following python scripts:
- [archive_sowfa_probe.py](archive_sowfa_probe.py): Script to convert virtual tower data (so called probes) from SOWFA specific output format to a netcdf file.
- [probe_raw2stats.py](probe_raw2stats.py): Script to calculate 10-min statistics from raw virtual tower data. Output is stored in a netcdf file
- [probe_raw2spectra.py](probe_raw2spectra.py): Script to calculate spectra from raw virtual tower data. Output is stored in a netcdf file

