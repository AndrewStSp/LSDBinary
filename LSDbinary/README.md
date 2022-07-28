# LSDbinary
A module to compute individual LSD profiles and LSD-based model spectra of both binary components, as well as to estimate their individual RVs

## Running the module
The module is ran from the command line with the following command: **./LSDbinary**
* *Input:*
  - 'LSD.conf': configuration file that configures a run
  - '_primary/secondary.corr': per binary component, a file that contains initial guess LSD profiles, local corrections, and estimate of the flux ratio. Output from the *LSDPrepare* module
  - '_primary/secondary.lin': per binary component, a file that contains line list ('line mask'). Output from the *LSDPrepare* module
  - 'vrads.init': a file that contains initial guess RVs (last two columns) for both binary components. The file also includes information about observed spectra and times/phases of observations (columns preceding RVs)
* *Output:* 
  - 'ObsSpectrumName.lsd1': LSD profile of the primary component
  - 'ObsSpectrumName.lsds': LSD profile of the secondary component
  - 'ObsSpectrumName.prof': file that contains LSD-based model spectra. Columns are: (1) RV; (2) composite spectrum (normalized flux); (3) spectrum of the primary component (normalized flux); (3) spectrum of the secondary component (normalized flux)
  - 'Vrads.dat': a file that contains estimated by *LSDbinary* RVs of both binary components and (optionally) their radii ratio 
* *Example (RZ Cas from Tkachenko et al. 2022):*
  - 'LSD.conf': configuration file that configures a run
  - 'input folder':
    - 'RZCas_primary.corr': initial guess LSD profiles and local corrections for the primary component
    - 'RZCas_primary.lin': line mask for the primary component
    - 'RZCas_secondary.corr': initial guess LSD profiles and local corrections for the secondary component
    - 'RZCas_secondary.lin': line mask for the secondary component
    - 'vrads.init': initial guess RVs
  - 'RZ_CAS_test_spectrum.asc': simulated observed spectra of the RZ Cas system