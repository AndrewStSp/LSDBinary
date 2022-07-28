# LSDPrepare
A module to take care of preparations for a run with the *LSDBinary* algorithm. This module covers all aspects of the *LSDInit* module described on the main page

## Running the module
The module is ran from the command line with the following command: **./LSDPrepare**
* *Input:*
  - 'LSDPrepare.conf': configuration file that configures a run
  - atmosphere models for both binary components (e.g. Kurucz models)
* *Output:* 
  - 'StarName_primary/secondary.corr': key output file that will be used by the *LSDBinary* algorithm
  - 'StarName_primary/secondary.lin': line list that has been used to compute the initial guess LSD profiles and will be employed by the *LSDBinary* algorithm
  - 'StarName_primary/secondary.rgs': synthetic spectrum that has been used to compute the initial guess LSD profiles
* *Example (RZ Cas from Tkachenko et al. 2022):*
  - 'LSDPrepare.conf': configuration file that configures a run
  - atmosphere models: 'lp0000_08800_0430_0020_on.mod' (primary) and 'lp0000_04800_0350_0020_on.mod' (secondary)
