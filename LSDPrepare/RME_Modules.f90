module synthV_data
real(8), allocatable :: I_c(:,:),I_l(:,:),f1(:),f1c(:),f2(:),f2c(:)  ! Local intensities (total,continuum),fluxes, fluxes of secondary
real(8) wavelength_first,wavelength_last,wavelength_step
integer number_of_wavelengths
parameter (ndeg=3,dmu=0.05d0) ! polinomial degree for approximate mu dependence of intensity
character(80) :: string_name
end module

module coordinates_data ! Coordinates in stellar surface
real(8) :: pi=3.1415925d0
parameter (nbin=128)
real(8) x_gau(nbin),wt_gau(nbin),y_gau(nbin,nbin),wt_y(nbin,nbin) ! coordinates and weights for gauss integration
real(8) mu_s(nbin,nbin) !mu cos(theta) for each point of visible disk
real(8) z(nbin) ! Doppler shift factor for earch x cartesian coordinate
end module

module stellar
real(8) ratio  ! R2/R1 - radii ratio
real(8) V_sini ! Vsini for primary
real(8) inclin ! inclination angle
real(8) phase,V1,V2 ! phase and radial velocities
real(8) ar1 ! ar=a/R1
real(8) Resolution, FWHM ! Resolved power and corresponding Gaussian parameter 
end module