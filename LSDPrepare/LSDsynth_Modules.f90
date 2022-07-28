module observed ! Observed spectra data
 character(100) observed_spectrum, line_list
 real(8) Vsini ! Vsini of the star
 real(8) wavelength_step ! Wavelength step in the synthetic spectrum
 real(8), allocatable :: wfirst(:),wlast(:) ! Wavelengths regions for fitting: first and last wavelength for each region of observed spectrum
 real(8), allocatable :: wave_obs(:), R_obs(:) ! given observed spectrum wavelengths and residual flux
 integer number_regions	! number of regions for fitting
 integer nobs ! number of points in observed spectrum
end module 

module mask
 real(8), allocatable :: w_mask(:),r_mask(:) ! wavelengths and intensities
 character(4), allocatable :: el_mask(:) ! chemical element IDs
 integer nlines ! number of lines in the mask
end module

module LSD_data
 real(8) FWHM_tot ! FWHM in velocity space (km/s)
 integer(4) nRV ! Number of pixels across the LSD profile
 real(8), allocatable :: vLSD(:), rLSD(:) ! RVs and intensity
 real(8), allocatable :: model(:) ! model spectrum
end module

