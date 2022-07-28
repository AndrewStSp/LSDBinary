module observed ! Observed spectra data
integer number_spectra ! number of observed spectra
integer if_FITS ! Switch to possible format of observed spectra: 0 - ASCII, 1- FITS (equidistant wavelengths 1D spectrum
real(8), allocatable :: wave_obs(:), R_obs(:) ! given observed spectrum wavelengths and residual flux
integer nobs ! number of points in observed spectrum
character(80), allocatable :: observed_names(:) ! names of files with observed spectra
real(8), allocatable :: wfirst(:),wlast(:) ! Wavelengths regions for fitting: first and last wavelength for each region of observed spectrum
integer number_regions ! number of regions for fitting
! FITS only
integer bitpix ! FITS data format
integer naxis1 ! size of data array
integer nb     ! size of FITS header
real(8) dwave  ! wavelength step
real(8) w0     ! first wavelength
end module

module synthetic ! Data from synthetic spectra
real(8), allocatable :: w_mask_prim(:),r_mask_prim(:),w_mask_sec(:),r_mask_sec(:) ! wavelengths and intensities
character(2), allocatable :: el_mask_prim(:),el_mask_sec(:)                       ! chemical elements ID
integer nlines_prim,nlines_sec                        ! numbers of lines in masks
real(8), allocatable :: wcorr(:), corr_prim(:), corr_sec(:), vLSD_syn1(:),rLSD_syn1(:),vLSD_syn2(:),rLSD_syn2(:) ! Correction to model spectra and synthetic LSD's
integer(8) nLSD_syn1,nLSD_syn2,nsynth ! numbers of points in synthetic LSD's and in correction data 
real(8), allocatable :: xint(:),yint(:) ! arrays for interpolation
real(8) cp_syn(3) ! approximated continuum flux_ratio from synthetic spectra
character(80) file_prim_corr, file_sec_corr, synth_path 
end module

module LSD_data
real(8) Rv1,Rv2  ! Components radial velocities
real(8) BJD ! Julian date
real(8) phase ! Given spectrum phase
real(8), allocatable :: phases(:)
integer nphases
real(8) R2, RR2 ! Square of radia ratio 
integer nprim,nsec  ! Number of points in LSD profiles
real(8), allocatable :: vLSD_prim(:),rLSD_prim(:), vLSD_sec(:),rLSD_sec(:) ! LSD profiles
real(8), allocatable :: model_prim(:),model_sec(:),model(:) ! model spectra of components and total model spectrum
real(8) cp(3) ! light factor polinomial coefficients
real(8) regpar(2) ! Tichonov quadratic regularization parameters
interface
 subroutine model_spectrum(Jacobian,mode)
	character(4), optional :: mode
	real(8), optional :: Jacobian(:,:)
 end subroutine model_spectrum
end interface
end module 

module velocity ! Initial guess for Vrad of components
integer if_calc ! Switch to option: read from table or calculate
real(8) e, K1,K2, Vgamma, omega ! Elements of orbits: eccentricity,semiamplitudes,gamma velocity, longitude of periastron
character(80), allocatable :: f_names(:) ! Names of files from table
real(8), allocatable :: ph(:),JD(:),v1(:),v2(:) ! Data from table
character slash
end module
 


