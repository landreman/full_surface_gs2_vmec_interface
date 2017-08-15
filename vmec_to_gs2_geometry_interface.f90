! vmec_to_gs2_geometry_interface.f90
! Written by Matt Landreman, University of Maryland
! Initial code written August 2017.
! Skip down ~25 lines for detailed description of the input and output parameters.

module vmec_to_gs2_geometry_interface_mod

  implicit none

  private

  public :: vmec_to_gs2_geometry_interface

  real :: theta_pest_target, zeta0
  real, dimension(2) :: vmec_radial_weight_full, vmec_radial_weight_half
  integer, dimension(2) :: vmec_radial_index_full, vmec_radial_index_half

contains

  subroutine vmec_to_gs2_geometry_interface(vmec_filename, nalpha, nzgrid, zeta_center, number_of_field_periods_to_include, &
       desired_normalized_toroidal_flux, vmec_surface_option, verbose, &
       normalized_toroidal_flux_used, safety_factor_q, shat, &
       alpha, zeta, bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0)

    use read_wout_mod, nzgrid_vmec => nzgrid  ! VMEC has a variable nzgrid which conflicts with our nzgrid, so rename vmec's version.

    implicit none

    !*********************************************************************
    ! Input parameters
    !*********************************************************************

    ! vmec_filename is the vmec wout_* file that will be read.
!    character(len=2000), intent(in) :: vmec_filename
    character(*), intent(in) :: vmec_filename

    ! nalpha is the number of grid points in the alpha coordinate:
    integer, intent(in) :: nalpha
  
    ! The zeta grid has nzgrid*2+1 points, including the "repeated" point at index -nzgrid and +nzgrid.
    integer, intent(in) :: nzgrid

    ! The zeta domain is centered at zeta_center. Setting zeta_center = 2*pi*N/nfp for any integer N should
    ! yield identical results to setting zeta_center = 0, where nfp is the number of field periods (as in VMEC).
    real, intent(in) :: zeta_center

    ! If number_of_field_periods_to_include is >= 1, then this parameter does what you think:
    ! the extent of the toroidal in zeta will be 2*pi*number_of_field_periods_to_include/nfp.
    ! If number_of_field_periods_to_include is < 1, the entire 2*pi toroidal domain will be included.
    ! If number_of_field_periods_to_include is > nfp, an error will result.
    integer, intent(in) :: number_of_field_periods_to_include

    ! The parameter desired_normalized_toroidal_flux determines which flux surface from the VMEC file will be used
    ! for the computation. This parameter should lie in the interval [0,1].
    real, intent(in) :: desired_normalized_toroidal_flux

    ! If vmec_surface_option = 0, the magnetic surface specified by desired_normalized_toroidal_flux will be used,
    ! by interpolating between the surfaces available in the vmec file.
    ! If vmec_surface_option = 1, the magnetic surface on vmec's HALF radial mesh will be used that is closest to desired_normalized_toroidal_flux.
    ! If vmec_surface_option = 2, the magnetic surface on vmec's FULL radial mesh will be used that is closest to desired_normalized_toroidal_flux.    
    ! Other values of vmec_surface_option will cause the program to abort with an error.
    integer, intent(in) :: vmec_surface_option

    ! If verbose is .true., lots of diagnostic information is printed.
    logical, intent(in) :: verbose

    !*********************************************************************
    ! Output quantities
    !*********************************************************************

    ! On exit, normalized_toroidal_flux_used holds the flux surface that was actually used for the geometry,
    ! as measured by psi_toroidal / psi_{toroidal,edge}
    real, intent(out) :: normalized_toroidal_flux_used

    ! Safety factor q = 1/iota
    real, intent(out) :: safety_factor_q

    ! Magnetic shear shat = (x/q) * (d q / d x) where x = Aminor_p * sqrt(psi_toroidal / psi_{toroidal,edge})
    ! and Aminor_p is the minor radius calculated by VMEC.
    real, intent(out) :: shat

    ! On exit, alpha holds the grid points in alpha = theta_p - iota * zeta, where theta_p is the PEST toroidal angle
    real, dimension(nalpha), intent(out) :: alpha

    ! On exit, zeta holds the grid points in the toroidal angle zeta
    real, dimension(-nzgrid:nzgrid), intent(out) :: zeta

    real, dimension(nalpha, -nzgrid:nzgrid), intent(out) :: bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0
    
    !*********************************************************************
    ! Variables used internally by this subroutine
    !*********************************************************************

!!$    integer :: i, itheta, izeta, imn, tic, toc, countrate, iflag, ierr, iopen, tic1, toc1, iunit
!!$    real(dp) :: angle, sinangle, cosangle, dsinangledtheta, dsinangledzeta, dcosangledtheta, dcosangledzeta
!!$    real(dp) :: angle2, sinangle2, cosangle2, dsinangle2dzeta, dcosangle2dzeta
!!$    real(dp) :: weight1, weight2, theta, r_temp, z_temp, dnorm
!!$    integer :: ntheta_coordTransform, nzeta_coordTransform
!!$    real(dp), dimension(:,:), allocatable :: r_coordTransform, z_coordTransform, major_R_squared
!!$    real(dp), dimension(:), allocatable :: rmnc_vmecLast, zmns_vmecLast
!!$    real(dp) :: rootSolve_abserr, rootSolve_relerr, theta_rootSolve_min, theta_rootSolve_max, theta_rootSolve_soln
!!$    integer :: fzeroFlag, mpol, ntor, jm, jn, index

    integer :: j, index, izeta, ialpha, which_surface, isurf, m, n
    real :: angle, sin_angle, cos_angle, temp
    real, dimension(:,:), allocatable :: theta_vmec
    integer :: ierr, iopen, fzero_flag, number_of_field_periods_to_include_final
    real :: dphi, iota, min_dr2, ds, d_pressure_d_s, scale_factor
    real :: theta_vmec_min, theta_vmec_max
    real, dimension(:), allocatable :: dr2, normalized_toroidal_flux_full_grid, normalized_toroidal_flux_half_grid
    real, dimension(:), allocatable :: d_pressure_d_s_on_half_grid
    real :: rootSolve_abserr, rootSolve_relerr
    logical :: non_Nyquist_mode_available, found_imn

    real, parameter :: pi = 3.1415926535897932d+0
    real, parameter :: zero = 0.0d+0
    real, parameter :: one = 1.0d+0

    !*********************************************************************
    ! VMEC variables of interest:
    ! ns = number of flux surfaces used by VMEC
    ! nfp = number of field periods, e.g. 5 for W7-X, 4 for HSX
    ! iotas = rotational transform (1/q) on the half grid.
    ! presf = pressure on the full grid.
    !
    ! All VMEC quantities (B, pressure, etc) are in SI units.
    ! 
    ! In VMEC, quantities on the half grid have the same number of array elements (ns) as quantities on the full grid,
    ! but the first array element is 0.
    !
    !*********************************************************************

    !*********************************************************************
    ! Beginning of executable statements.
    !*********************************************************************

    if (verbose) print *,"Entering subroutine vmec_to_gs2_geometry_interface."
    
    !*********************************************************************
    ! Do some validation.
    !*********************************************************************

    if (nalpha<1) then
       print *,"Error! nalpha must be >= 1. Instead it is",nalpha
       stop
    end if

    if (nzgrid<1) then
       print *,"Error! nzgrid must be >= 1. Instead it is",nzgrid
       stop
    end if

    if (desired_normalized_toroidal_flux <= 0) then
       print *,"Error! desired_normalized_toroidal_flux must be >0. Instead it is",desired_normalized_toroidal_flux
       stop
    end if

    if (desired_normalized_toroidal_flux > 1) then
       print *,"Error! desired_normalized_toroidal_flux must be <= 1. Instead it is",desired_normalized_toroidal_flux
       stop
    end if

    !*********************************************************************
    ! Read in everything from the vmec wout file using libstell.
    !*********************************************************************

    if (verbose) print *,"  About to read VMEC wout file ",trim(vmec_filename)
    call read_wout_file(vmec_filename, ierr, iopen)
    if (iopen .ne. 0) stop 'error opening wout file'
    if (ierr .ne. 0) stop 'error reading wout file'
    if (verbose) print *,"  Successfully read VMEC data from ",trim(vmec_filename)
       
    if (verbose) print *,"  Number of field periods (nfp):",nfp
    if (verbose) print *,"  Stellarator-asymmetric? (lasym):",lasym
!!$    if (lasym) then
!!$       stop "Error! geometry_option_plasma=4 is not yet implemented for lasym=true"
!!$    end if

    ! There is a bug in libstell read_wout_file for ASCII-format wout files, in which the xm_nyq and xn_nyq arrays are sometimes
    ! not populated. The next few lines here provide a workaround:
    if (maxval(abs(xm_nyq)) < 1 .and. maxval(abs(xn_nyq)) < 1) then
       if (mnmax_nyq == mnmax) then
          if (verbose) print *,"xm_nyq and xn_nyq arrays are not populated in the wout file. Using xm and xn instead."
          xm_nyq = xm
          xn_nyq = xn
       else
          print *,"Error! xm_nyq and xn_nyq arrays are not populated in the wout file, and mnmax_nyq != mnmax."
          stop
       end if
    end if

    ! --------------------------------------------------------------------------------
    ! Do some sanity checking to ensure the VMEC arrays have some expected properties.
    ! --------------------------------------------------------------------------------

    ! 'phi' is vmec's array of the toroidal flux (not divided by 2pi!) on its radial grid.
    if (abs(phi(1)) > 1d-14) then
       print *,"Error! VMEC phi array does not begin with 0."
       print *,"phi:",phi
       stop
    end if

    dphi = phi(2) - phi(1)
    do j=3,ns
       if (abs(phi(j)-phi(j-1)-dphi) > 1d-11) then
          print *,"Error! VMEC phi array is not uniformly spaced."
          print *,"phi:",phi
          stop
       end if
    end do

    ! The variable called 'phips' in the wout file is called just 'phip' in read_wout_mod.F.
    ! phips is on the half-mesh, so skip first point.
    do j=2,ns
       if (abs(phip(j)+phi(ns)/(2*pi)) > 1d-11) then
          print *,"Error! VMEC phips array is not constant and equal to -phi(ns)/(2*pi)."
          print *,"phip(s):",phip
          stop
       end if
    end do

    ! The first mode in the m and n arrays should be m=n=0:
    if (xm(1) .ne. 0) stop "First element of xm in the wout file should be 0."
    if (xn(1) .ne. 0) stop "First element of xn in the wout file should be 0."
    if (xm_nyq(1) .ne. 0) stop "First element of xm_nyq in the wout file should be 0."
    if (xn_nyq(1) .ne. 0) stop "First element of xn_nyq in the wout file should be 0."

    ! --------------------------------------------------------------------------------
    ! End of sanity checks.
    ! --------------------------------------------------------------------------------

    allocate(normalized_toroidal_flux_full_grid(ns))
    normalized_toroidal_flux_full_grid = [( real(j-1)/(ns-1), j=1,ns )]

    ! Build an array of the half grid points:
    allocate(normalized_toroidal_flux_half_grid(ns-1))
    do j = 1,ns-1
       normalized_toroidal_flux_half_grid(j) = (normalized_toroidal_flux_full_grid(j) + normalized_toroidal_flux_full_grid(j+1))*(0.5d+0)
    end do

    !*********************************************************************
    ! Determine which flux surface to use, based on 
    ! desired_normalized_toroidal_flux and vmec_surface_option.
    !*********************************************************************

    ! Possible values of vmec_surface_option:
    ! 0 = Use the exact radius requested.
    ! 1 = Use the nearest value of the VMEC half grid.
    ! 2 = Use the nearest value of the VMEC full grid.
    
    select case (vmec_surface_option)
    case (0)
       ! Use exact radius requested.
       normalized_toroidal_flux_used = desired_normalized_toroidal_flux

    case (1)
       ! Use nearest value of the VMEC half grid

       ! Compute differences
       allocate(dr2(ns-1))
       dr2 = (normalized_toroidal_flux_half_grid - desired_normalized_toroidal_flux) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns-1
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       normalized_toroidal_flux_used = normalized_toroidal_flux_half_grid(index)
       deallocate(dr2)

    case (2)
       ! Use nearest value of the VMEC full grid

       ! Compute differences
       allocate(dr2(ns))
       dr2 = (normalized_toroidal_flux_full_grid - desired_normalized_toroidal_flux) ** 2

       index = 1
       min_dr2 = dr2(1)
       ! Find the index of minimum error:
       do j=2,ns
          if (dr2(j)<min_dr2) then
             index = j
             min_dr2 = dr2(j)
          end if
       end do

       normalized_toroidal_flux_used = normalized_toroidal_flux_full_grid(index)
       deallocate(dr2)

    case default
       print *,"Error! vmec_surface_option must be 0, 1, or 2. It is instead ",vmec_surface_option
       stop
    end select

    ! --------------------------------------------------------------------------------
    ! Done choosing the actual radius to use.
    ! --------------------------------------------------------------------------------

    ! In general, we get quantities for gs2 by linear interpolation, taking a weighted average of the quantity from
    ! 2 surfaces in the VMEC file. Sometimes the weights are 0 and 1, i.e. no interpolation is needed.

    ! For any VMEC quantity Q on the full grid, the value used in GS2 will be
    !  Q_gs2 = Q(vmec_radial_index_full(1))*vmec_radial_weight_full(1) + Q(vmec_radial_index_full(2))*vmec_radial_weight_full(2)

    ! For any VMEC quantity Q on the half grid, the value used in GS2 will be
    !  Q_gs2 = Q(vmec_radial_index_half(1))*vmec_radial_weight_half(1) + Q(vmec_radial_index_half(2))*vmec_radial_weight_half(2)


    ! Handle quantities for the full grid
    if (normalized_toroidal_flux_used>1) then
       stop "Error! normalized_toroidal_flux_used cannot be >1"
    elseif (normalized_toroidal_flux_used<0) then
       stop "Error! normalized_toroidal_flux_used cannot be <0"
    elseif (normalized_toroidal_flux_used==1) then
       vmec_radial_index_full(1) = ns-1
       vmec_radial_index_full(2) = ns
       vmec_radial_weight_full(1) = zero
    else
       ! normalized_toroidal_flux_used is >= 0 and <1
       ! This is the most common case.
       vmec_radial_index_full(1) = floor(normalized_toroidal_flux_used*(ns-1))+1
       vmec_radial_index_full(2) = vmec_radial_index_full(1) + 1
       vmec_radial_weight_full(1) = vmec_radial_index_full(1) - normalized_toroidal_flux_used*(ns-one)
    end if
    vmec_radial_weight_full(2) = one - vmec_radial_weight_full(1)
    
    ! Handle quantities for the half grid
    if (normalized_toroidal_flux_used < normalized_toroidal_flux_half_grid(1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the magnetic axis.) Results are likely to be inaccurate."

       ! We start at element 2 since element 1 is always 0 for quantities on the half grid.
       vmec_radial_index_half(1) = 2
       vmec_radial_index_half(2) = 3
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_used) / (normalized_toroidal_flux_half_grid(2) - normalized_toroidal_flux_half_grid(1))

    elseif (normalized_toroidal_flux_used > normalized_toroidal_flux_half_grid(ns-1)) then
       print *,"Warning: extrapolating beyond the end of VMEC's half grid."
       print *,"(Extrapolating towards the last closed flux surface.) Results may be inaccurate."
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_used) &
            / (normalized_toroidal_flux_half_grid(ns-1) - normalized_toroidal_flux_half_grid(ns-2))

    elseif (normalized_toroidal_flux_used == normalized_toroidal_flux_half_grid(ns-1)) then
       ! We are exactly at the last point of the half grid
       vmec_radial_index_half(1) = ns-1
       vmec_radial_index_half(2) = ns
       vmec_radial_weight_half(1) = zero
    else
       ! normalized_toroidal_flux_used is inside the half grid.
       ! This is the most common case.
       vmec_radial_index_half(1) = floor(normalized_toroidal_flux_used*(ns-1) + 0.5d+0)+1
       if (vmec_radial_index_half(1) < 2) then
          ! This can occur sometimes due to roundoff error.
          vmec_radial_index_half(1) = 2
       end if
       vmec_radial_index_half(2) = vmec_radial_index_half(1) + 1
       vmec_radial_weight_half(1) = vmec_radial_index_half(1) - normalized_toroidal_flux_used*(ns-one) - (0.5d+0)
    end if
    vmec_radial_weight_half(2) = one-vmec_radial_weight_half(1)

    if (verbose) then
       if (abs(vmec_radial_weight_half(1)) < 1e-14) then
          print "(a,i4,a,i4,a)","   Using radial index ",vmec_radial_index_half(2)," of ",ns," from vmec's half mesh."
       elseif (abs(vmec_radial_weight_half(2)) < 1e-14) then
          print "(a,i4,a,i4,a)","   Using radial index ",vmec_radial_index_half(1)," of ",ns," from vmec's half mesh."
       else
          print "(a,i4,a,i4,a,i4,a)", "   Interpolating using radial indices ",vmec_radial_index_half(1)," and ",vmec_radial_index_half(2),&
               " of ",ns," from vmec's half mesh."
          print "(a,f17.14,a,f17.14)", "   Weights for half mesh = ",vmec_radial_weight_half(1)," and ",vmec_radial_weight_half(2)
          print "(a,i4,a,i4,a,i4,a)", "   Interpolating using radial indices ",vmec_radial_index_full(1)," and ",vmec_radial_index_full(2),&
               " of ",ns," from vmec's full mesh."
          print "(a,f17.14,a,f17.14)", "   Weights for full mesh = ",vmec_radial_weight_full(1)," and ",vmec_radial_weight_full(2)
       end if
    end if

    !*********************************************************************
    ! Evaluate several radial-profile functions at the flux surface
    ! we ended up choosing.
    !*********************************************************************

    iota = iotas(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + iotas(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    if (verbose) print *,"  iota =",iota
    safety_factor_q = 1/iota

    allocate(d_pressure_d_s_on_half_grid(ns))
    d_pressure_d_s_on_half_grid = 0
    ds = normalized_toroidal_flux_full_grid(2) - normalized_toroidal_flux_full_grid(1)
    d_pressure_d_s_on_half_grid(2:ns) = (presf(2:ns) - presf(1:ns-1)) / ds
    d_pressure_d_s =  &
         d_pressure_d_s_on_half_grid(vmec_radial_index_half(1)) * vmec_radial_weight_half(1) &
         + d_pressure_d_s_on_half_grid(vmec_radial_index_half(2)) * vmec_radial_weight_half(2)
    if (verbose) print *,"  d pressure / d s =",d_pressure_d_s

    !*********************************************************************
    ! Set up the coordinate grids.
    !*********************************************************************

    alpha = [( ((j-1)*2*pi) / nalpha, j=1, nalpha )]

    if (number_of_field_periods_to_include > nfp) then
       print *,"Error! number_of_field_periods_to_include > nfp"
       print *,"  number_of_field_periods_to_include =",number_of_field_periods_to_include
       print *,"  nfp =",nfp
       stop
    end if
    number_of_field_periods_to_include_final = number_of_field_periods_to_include
    if (number_of_field_periods_to_include <1) then
       number_of_field_periods_to_include_final = nfp
       if (verbose) print *,"  Since number_of_field_periods_to_include was < 1, it is being reset to nfp =",nfp
    end if

    zeta = [( zeta_center + (pi*j*number_of_field_periods_to_include_final)/(nfp*nzgrid), j=-nzgrid, nzgrid )]

    !*********************************************************************
    ! We know theta_pest = alpha + iota * zeta, but we need to determine
    ! theta_vmec = theta_pest - Lambda.
    !*********************************************************************

    allocate(theta_vmec(nalpha, -nzgrid:nzgrid))

    if (verbose) print *,"  Beginning root solves to determine theta_vmec."
    rootSolve_abserr = 1.0d-10
    rootSolve_relerr = 1.0d-10
    !open(unit=5,file="testStraightFieldLines",status='new',form='formatted')
    !write (5,*) ntheta_coordTransform, nzeta_coordTransform
    do izeta = -nzgrid, nzgrid
       zeta0 = zeta(izeta)
       do ialpha = 1,nalpha
          theta_pest_target = alpha(ialpha) + iota * zeta0
          ! Guess that theta_vmec will be within 0.3 radians of theta_pest:
          theta_vmec_min = theta_pest_target - 0.3
          theta_vmec_max = theta_pest_target + 0.3
          
          ! In the 4th argument, we are telling the root-finder (fzero) to use theta_pest as the initial guess for theta_vmec.
          call fzero(fzero_residual, theta_vmec_min, theta_vmec_max, theta_pest_target, &
               rootSolve_relerr, rootSolve_abserr, fzero_flag)
          ! Note: fzero returns its answer in theta_vmec_min
          theta_vmec(ialpha,izeta) = theta_vmec_min
          if (fzero_flag == 4) then
             stop "ERROR: fzero returned error 4: no sign change in residual"
          else if (fzero_flag > 2) then
             print *,"WARNING: fzero returned an error code:",fzero_flag
          end if
       end do
    end do
    if (verbose) then
       print *,"  Done with root solves. Here comes theta_vmec:"
       do j = 1, nalpha
          print *,theta_vmec(j,:)
       end do
    end if
       
    !*********************************************************************
    ! Now that we know the grid points in theta_vmec, we can evaluate
    ! all the geometric quantities on the grid points.
    !*********************************************************************

    do imn_nyq = 1, mnmax_nyq ! All the quantities we need except R and Z use the _nyq mode numbers.
       m = xm_nyq(imn_nyq)
       n = xn_nyq(imn_nyq)/nfp

       if (abs(m) >= mpol .or. abs(n) > ntor) then
          non_Nyquist_mode_available = .false.
       else
          non_Nyquist_mode_available = .true.
          ! Find the imn in the non-Nyquist arrays that corresponds to the same m and n.
          found_imn = .false.
          do imn = 1,mnmax
             if (xm(imn)==m .and. xn(imn)==n*nfp) then
                found_imn = .true.
                exit
             end if
          end do
          if ((xm(imn) .ne. m) .or. (xn(imn) .ne. n*nfp)) stop "Something went wrong!"
          if (.not. found_imn) stop "Error! imn could not be found matching the given imn_nyq."
       end if

       ! -----------------------------------------------------
       ! First, consider just the stellarator-symmetric terms:
       ! -----------------------------------------------------
       
       b = bmnc(imn_nyq,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
            + bmnc(imn_nyq,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)
       
       ! All quantities are multiplied by a variable scale_factor which can in principle depend on m and n.
       ! For now we just set scale_factor = 1. In the future, scale_factor could be used to lower the
       ! symmetry-breaking Fourier components, or filter out certain Fourier components in some way.
       scale_factor = 1
       b = b*scale_factor
	

          ! Evaluate the radial derivatives we will need:
          dpsi = phi(2)/(2*pi)  ! Doesn't need to be in the loops, but here for convenience.

          ! B, B_sub_theta, and B_sub_zeta are on the half mesh, so their radial derivatives are on the full mesh.
          ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.

          vmec_dBHatdpsiHat(2:ns-1) = (bmnc(imn_nyq,3:ns) - bmnc(imn_nyq,2:ns-1)) / dpsi
          ! Simplistic "extrapolation" at the endpoints:
          vmec_dBHatdpsiHat(1) = vmec_dBHatdpsiHat(2)
          vmec_dBHatdpsiHat(ns) = vmec_dBHatdpsiHat(ns-1)

          vmec_dBHat_sub_theta_dpsiHat(2:ns-1) = (bsubumnc(imn_nyq,3:ns) - bsubumnc(imn_nyq,2:ns-1)) / dpsi
          vmec_dBHat_sub_theta_dpsiHat(1) = vmec_dBHat_sub_theta_dpsiHat(2)
          vmec_dBHat_sub_theta_dpsiHat(ns) = vmec_dBHat_sub_theta_dpsiHat(ns-1)

          vmec_dBHat_sub_zeta_dpsiHat(2:ns-1) = (bsubvmnc(imn_nyq,3:ns) - bsubvmnc(imn_nyq,2:ns-1)) / dpsi
          vmec_dBHat_sub_zeta_dpsiHat(1) = vmec_dBHat_sub_zeta_dpsiHat(2)
          vmec_dBHat_sub_zeta_dpsiHat(ns) = vmec_dBHat_sub_zeta_dpsiHat(ns-1)

          if (non_Nyquist_mode_available) then
             vmec_dRdpsiHat(2:ns) = (rmnc(imn,2:ns) - rmnc(imn,1:ns-1)) / dpsi
             vmec_dRdpsiHat(1) = 0

             vmec_dZdpsiHat(2:ns) = (zmns(imn,2:ns) - zmns(imn,1:ns-1)) / dpsi
             vmec_dZdpsiHat(1) = 0
          else
             vmec_dRdpsiHat = 0
             vmec_dZdpsiHat = 0
          end if

          ! End of evaluating radial derivatives.
             
          do ialpha = 1,nalpha
             do izeta = -nzgrid, nzgrid
                angle = m * theta_vmec(ialpha,izeta) - n * nfp * zeta(izeta)
                cos_angle = cos(angle)
                sin_angle = sin(angle)
                
                BHat(ialpha,izeta) = BHat(ialpha,izeta) + b * cos_angle
                
                dbHatdtheta(ialpha,izeta) = dBHatdtheta(ialpha,izeta) - m * b * sin_angle
                
                dbHatdzeta(ialpha,izeta) = dBHatdzeta(ialpha,izeta) + n * nfp * b * sin_angle
                
                do isurf = 1,2
                   ! Handle Jacobian:
                   temp = gmnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scale_factor
                   sqrt_g(ialpha,izeta) = sqrt_g(ialpha,izeta) + temp * cos_angle

                   ! Handle B sup theta:
                   ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                   temp = bsupumnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scale_factor
                   BHat_sup_theta(ialpha,izeta) = BHat_sup_theta(ialpha,izeta) + temp * cos_angle
                   dBHat_sup_theta_dzeta(ialpha,izeta) = dBHat_sup_theta_dzeta(ialpha,izeta) + n * nfp * temp * sin_angle

                   ! Handle B sup zeta:
                   ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or nfp needed.
                   temp = bsupvmnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scale_factor
                   BHat_sup_zeta(ialpha,izeta) = BHat_sup_zeta(ialpha,izeta) + temp * cos_angle
                   dBHat_sup_zeta_dtheta(ialpha,izeta) = dBHat_sup_zeta_dtheta(ialpha,izeta) - m * temp * sin_angle

                   ! Handle B sub theta:
                   ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                   temp = bsubumnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scale_factor
                   BHat_sub_theta(ialpha,izeta) = BHat_sub_theta(ialpha,izeta) + temp * cos_angle
                   dBHat_sub_theta_dzeta(ialpha,izeta) = dBHat_sub_theta_dzeta(ialpha,izeta) + n * nfp * temp * sin_angle

                   ! Handle B sub zeta:
                   ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                   temp = bsubvmnc(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                   temp = temp*scale_factor
                   BHat_sub_zeta(ialpha,izeta) = BHat_sub_zeta(ialpha,izeta) + temp * cos_angle
                   dBHat_sub_zeta_dtheta(ialpha,izeta) = dBHat_sub_zeta_dtheta(ialpha,izeta) - m * temp * sin_angle

                   ! Handle B sub psi.
                   ! Unlike the other components of B, this one is on the full mesh.
                   ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                   temp = bsubsmns(imn_nyq,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                   temp = temp*scale_factor
                   BHat_sub_psi(ialpha,izeta) = BHat_sub_psi(ialpha,izeta) + temp * sin_angle
                   dBHat_sub_psi_dtheta(ialpha,izeta) = dBHat_sub_psi_dtheta(ialpha,izeta) + m * temp * cos_angle
                   dBHat_sub_psi_dzeta(ialpha,izeta)  = dBHat_sub_psi_dzeta(ialpha,izeta) - n * nfp * temp * cos_angle

                   ! Handle dBHatdpsiHat.
                   ! Since bmnc is on the half mesh, its radial derivative is on the full mesh.
                   temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                   temp = temp*scale_factor
                   dBHatdpsiHat(ialpha,izeta) = dBHatdpsiHat(ialpha,izeta) + temp * cos_angle

                   ! Handle dBHat_sub_theta_dpsiHat.
                   ! Since bsubumnc is on the half mesh, its radial derivative is on the full mesh.
                   temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                   temp = temp*scale_factor
                   dBHat_sub_theta_dpsiHat(ialpha,izeta) = dBHat_sub_theta_dpsiHat(ialpha,izeta) + temp * cos_angle

                   ! Handle dBHat_sub_zeta_dpsiHat.
                   ! Since bsubvmnc is on the half mesh, its radial derivative is on the full mesh.
                   temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                   temp = temp*scale_factor
                   dBHat_sub_zeta_dpsiHat(ialpha,izeta) = dBHat_sub_zeta_dpsiHat(ialpha,izeta) + temp * cos_angle

                   ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                   if (non_Nyquist_mode_available) then

                      ! Handle R, which is on the full mesh
                      temp = rmnc(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scale_factor
                      R(ialpha,izeta) = R(ialpha,izeta) + temp * cos_angle
                      dRdtheta(ialpha,izeta) = dRdtheta(ialpha,izeta) - temp * m * sin_angle
                      dRdzeta(ialpha,izeta)  = dRdzeta(ialpha,izeta)  + temp * n * nfp * sin_angle

                      ! Handle Z, which is on the full mesh
                      temp = zmns(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scale_factor
                      !Z(ialpha,izeta) = Z(ialpha,izeta) + temp * sin_angle  ! We don't actually need Z itself, only derivatives of Z.
                      dZdtheta(ialpha,izeta) = dZdtheta(ialpha,izeta) + temp * m * cos_angle
                      dZdzeta(ialpha,izeta)  = dZdzeta(ialpha,izeta)  - temp * n * nfp * cos_angle

                      ! Handle dRdpsiHat.
                      ! Since R is on the full mesh, its radial derivative is on the half mesh.
                      temp = vmec_dRdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor
                      dRdpsiHat(ialpha,izeta) = dRdpsiHat(ialpha,izeta) + temp * cos_angle

                      ! Handle dZdpsiHat.
                      ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                      temp = vmec_dZdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor
                      dZdpsiHat(ialpha,izeta) = dZdpsiHat(ialpha,izeta) + temp * sin_angle

                   end if
                end do
             end do
          end do

       ! -----------------------------------------------------
       ! Now consider the stellarator-asymmetric terms.
       ! NOTE: This functionality has not been tested as thoroughly !!!
       ! -----------------------------------------------------

       if (lasym) then

          b = bmns(imn_nyq,vmecRadialIndex_half(1)) * vmecRadialWeight_half(1) &
               + bmns(imn_nyq,vmecRadialIndex_half(2)) * vmecRadialWeight_half(2)

          ! Set scale_factor to rippleScale for non-axisymmetric or non-quasisymmetric modes
          scale_factor = setScaleFactor(n,m)
          b = b*scale_factor
          
             
             ! Evaluate the radial derivatives we will need:
             dpsi = phi(2)/(2*pi)  ! Doesn't need to be in the loops, but here for convenience.
             
             ! B, B_sub_theta, and B_sub_zeta are on the half mesh, so their radial derivatives are on the full mesh.
             ! R and Z are on the full mesh, so their radial derivatives are on the half mesh.
             
             vmec_dBHatdpsiHat(2:ns-1) = (bmns(imn_nyq,3:ns) - bmns(imn_nyq,2:ns-1)) / dpsi
             ! Simplistic "extrapolation" at the endpoints:
             vmec_dBHatdpsiHat(1) = vmec_dBHatdpsiHat(2)
             vmec_dBHatdpsiHat(ns) = vmec_dBHatdpsiHat(ns-1)
             
             vmec_dBHat_sub_theta_dpsiHat(2:ns-1) = (bsubumns(imn_nyq,3:ns) - bsubumns(imn_nyq,2:ns-1)) / dpsi
             vmec_dBHat_sub_theta_dpsiHat(1) = vmec_dBHat_sub_theta_dpsiHat(2)
             vmec_dBHat_sub_theta_dpsiHat(ns) = vmec_dBHat_sub_theta_dpsiHat(ns-1)
             
             vmec_dBHat_sub_zeta_dpsiHat(2:ns-1) = (bsubvmns(imn_nyq,3:ns) - bsubvmns(imn_nyq,2:ns-1)) / dpsi
             vmec_dBHat_sub_zeta_dpsiHat(1) = vmec_dBHat_sub_zeta_dpsiHat(2)
             vmec_dBHat_sub_zeta_dpsiHat(ns) = vmec_dBHat_sub_zeta_dpsiHat(ns-1)
             
             if (non_Nyquist_mode_available) then
                vmec_dRdpsiHat(2:ns) = (rmns(imn,2:ns) - rmns(imn,1:ns-1)) / dpsi
                vmec_dRdpsiHat(1) = 0
             
                vmec_dZdpsiHat(2:ns) = (zmnc(imn,2:ns) - zmnc(imn,1:ns-1)) / dpsi
                vmec_dZdpsiHat(1) = 0
             else
                vmec_dRdpsiHat = 0
                vmec_dZdpsiHat = 0
             end if
             
             ! End of evaluating radial derivatives.
             
             do ialpha = 1,nalpha
                do izeta = -nzgrid, nzgrid
                   angle = m * theta_vmec(ialpha,izeta) - n * nfp * zeta(izeta)
                   cos_angle = cos(angle)
                   sin_angle = sin(angle)
                      
                   BHat(ialpha,izeta) = BHat(ialpha,izeta) + b * sin_angle
                   dbHatdtheta(ialpha,izeta) = dBHatdtheta(ialpha,izeta) + m * b * cos_angle
                   dbHatdzeta(ialpha,izeta) = dBHatdzeta(ialpha,izeta) - n * nfp * b * cos_angle
                   
                   do isurf = 1,2
                      ! Handle Jacobian:
                      temp = gmns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor
                      sqrt_g(ialpha,izeta) = sqrt_g(ialpha,izeta) + temp * sin_angle
                         
                      ! Handle B sup theta:
                      ! Note that VMEC's bsupumnc and bsupumns are exactly the same as SFINCS's BHat_sup_theta, with no conversion factors of 2pi needed.
                      temp = bsupumns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor	
                      BHat_sup_theta(ialpha,izeta) = BHat_sup_theta(ialpha,izeta) + temp * sin_angle
                      dBHat_sup_theta_dzeta(ialpha,izeta) = dBHat_sup_theta_dzeta(ialpha,izeta) - n * nfp * temp * cos_angle
                      
                      ! Handle B sup zeta:
                      ! Note that VMEC's bsupvmnc and bsupvmns are exactly the same as SFINCS's BHat_sup_zeta, with no conversion factors of 2pi or nfp needed.
                      temp = bsupvmns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor
                      BHat_sup_zeta(ialpha,izeta) = BHat_sup_zeta(ialpha,izeta) + temp * sin_angle
                      dBHat_sup_zeta_dtheta(ialpha,izeta) = dBHat_sup_zeta_dtheta(ialpha,izeta) + m * temp * cos_angle
                      
                      ! Handle B sub theta:
                      ! Note that VMEC's bsubumnc and bsubumns are exactly the same as SFINCS's BHat_sub_theta, with no conversion factors of 2pi needed.
                      temp = bsubumns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor
                      BHat_sub_theta(ialpha,izeta) = BHat_sub_theta(ialpha,izeta) + temp * sin_angle
                      dBHat_sub_theta_dzeta(ialpha,izeta) = dBHat_sub_theta_dzeta(ialpha,izeta) - n * nfp * temp * cos_angle
                      
                      ! Handle B sub zeta:
                      ! Note that VMEC's bsubvmnc and bsubvmns are exactly the same as SFINCS's BHat_sub_zeta, with no conversion factors of 2pi needed.
                      temp = bsubvmns(imn_nyq,vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                      temp = temp*scale_factor
                      BHat_sub_zeta(ialpha,izeta) = BHat_sub_zeta(ialpha,izeta) + temp * sin_angle
                      dBHat_sub_zeta_dtheta(ialpha,izeta) = dBHat_sub_zeta_dtheta(ialpha,izeta) + m * temp * cos_angle
                      
                      ! Handle B sub psi.
                      ! Unlike the other components of B, this one is on the full mesh.
                      ! Notice B_psi = B_s * (d s / d psi), and (d s / d psi) = 1 / psiAHat
                      temp = bsubsmnc(imn_nyq,vmecRadialIndex_full(isurf)) / psiAHat * vmecRadialWeight_full(isurf)
                      temp = temp*scale_factor
                      BHat_sub_psi(ialpha,izeta) = BHat_sub_psi(ialpha,izeta) + temp * cos_angle
                      dBHat_sub_psi_dtheta(ialpha,izeta) = dBHat_sub_psi_dtheta(ialpha,izeta) - m * temp * sin_angle
                      dBHat_sub_psi_dzeta(ialpha,izeta)  = dBHat_sub_psi_dzeta(ialpha,izeta) + n * nfp * temp * sin_angle
                      
                      ! Handle dBHatdpsiHat.
                      ! Since bmns is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHatdpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scale_factor
                      dBHatdpsiHat(ialpha,izeta) = dBHatdpsiHat(ialpha,izeta) + temp * sin_angle
                      
                      ! Handle dBHat_sub_theta_dpsiHat.
                      ! Since bsubumns is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHat_sub_theta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scale_factor
                      dBHat_sub_theta_dpsiHat(ialpha,izeta) = dBHat_sub_theta_dpsiHat(ialpha,izeta) + temp * sin_angle
                      
                      ! Handle dBHat_sub_zeta_dpsiHat.
                      ! Since bsubvmns is on the half mesh, its radial derivative is on the full mesh.
                      temp = vmec_dBHat_sub_zeta_dpsiHat(vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                      temp = temp*scale_factor
                      dBHat_sub_zeta_dpsiHat(ialpha,izeta) = dBHat_sub_zeta_dpsiHat(ialpha,izeta) + temp * sin_angle
                      
                      ! Handle arrays that use xm and xn instead of xm_nyq and xn_nyq.
                      if (non_Nyquist_mode_available) then

                         ! Handle R, which is on the full mesh
                         temp = rmns(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scale_factor
                         R(ialpha,izeta) = R(ialpha,izeta) + temp * sin_angle
                         dRdtheta(ialpha,izeta) = dRdtheta(ialpha,izeta) + temp * m * cos_angle
                         dRdzeta(ialpha,izeta)  = dRdzeta(ialpha,izeta)  - temp * n * nfp * cos_angle

                         ! Handle Z, which is on the full mesh
                         temp = zmnc(imn,vmecRadialIndex_full(isurf)) * vmecRadialWeight_full(isurf)
                         temp = temp*scale_factor
                         ! Z(ialpha,izeta) = Z(ialpha,izeta) + temp * cos_angle   ! We don't actually need Z itself, only derivatives of Z.
                         dZdtheta(ialpha,izeta) = dZdtheta(ialpha,izeta) - temp * m * sin_angle
                         dZdzeta(ialpha,izeta)  = dZdzeta(ialpha,izeta)  + temp * n * nfp * sin_angle

                         ! Handle dRdpsiHat.
                         ! Since R is on the full mesh, its radial derivative is on the half mesh.
                         temp = vmec_dRdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scale_factor
                         dRdpsiHat(ialpha,izeta) = dRdpsiHat(ialpha,izeta) + temp * sin_angle
                         
                         ! Handle dZdpsiHat.
                         ! Since Z is on the full mesh, its radial derivative is on the half mesh.
                         temp = vmec_dZdpsiHat(vmecRadialIndex_half(isurf)) * vmecRadialWeight_half(isurf)
                         temp = temp*scale_factor
                         dZdpsiHat(ialpha,izeta) = dZdpsiHat(ialpha,izeta) + temp * cos_angle
                      end if
                   end do
                end do
             end do
          end if
       end do

    do izeta = 1,Nzeta
       cos_angle = cos(zeta(izeta))
       sin_angle = sin(zeta(izeta))

       ! X = R * cos(zeta)
       dXdtheta(:,izeta) = dRdtheta(:,izeta) * cos_angle
       dXdzeta(:,izeta) = dRdzeta(:,izeta) * cos_angle - R(:,izeta) * sin_angle
       dXdpsiHat(:,izeta) = dRdpsiHat(:,izeta) * cos_angle

       ! Y = R * sin(zeta)
       dYdtheta(:,izeta) = dRdtheta(:,izeta) * sin_angle
       dYdzeta(:,izeta) = dRdzeta(:,izeta) * sin_angle + R(:,izeta) * cos_angle
       dYdpsiHat(:,izeta) = dRdpsiHat(:,izeta) * sin_angle
    end do

    if (verbose) print *,"Leaving vmec_to_gs2_geometry_interface."

  end subroutine vmec_to_gs2_geometry_interface

  ! --------------------------------------------------------------------------

  function fzero_residual(theta_vmec_try)

    use read_wout_mod, only: xm, xn, mnmax, lmns, lmnc, lasym
    ! Note that lmns and lmnc use the non-Nyquist xm, xn, and mnmax.
    ! Also note that lmns and lmnc are on the HALF grid.

    implicit none
    
    real :: theta_vmec_try, fzero_residual
    real :: angle, sinangle, cosangle
    integer :: imn, which_surface

    ! residual = (theta_pest based on theta_vmec_try) - theta_pest_target = theta_vmec_try + Lambda - theta_pest_target
    fzero_residual = theta_vmec_try - theta_pest_target

    do imn = 1, mnmax
       angle = xm(imn)*theta_vmec_try - xn(imn)*zeta0
       sinangle = sin(angle)
       cosangle = cos(angle)
       do which_surface = 1,2
          fzero_residual = fzero_residual + vmec_radial_weight_half(which_surface) * lmns(imn,vmec_radial_index_half(which_surface)) * sinangle
          if (lasym) then
             fzero_residual = fzero_residual + vmec_radial_weight_half(which_surface) * lmnc(imn,vmec_radial_index_half(which_surface)) * cosangle
          end if
       end do
    end do

  end function fzero_residual

end module vmec_to_gs2_geometry_interface_mod
