! vmec_to_gs2_geometry_interface.f90
! Written by Matt Landreman, University of Maryland
! Initial code written August 2017.
! Skip down ~20 lines for detailed description of the input and output parameters.

module vmec_to_gs2_geometry_interface_mod

  implicit none

  private

  public :: vmec_to_gs2_geometry_interface

  real :: theta_rootSolve_target, zeta

contains

  subroutine vmec_to_gs2_geometry_interface(vmec_filename, nalpha, nzgrid, zeta_center, number_of_field_periods_to_include, &
       desired_normalized_toroidal_flux, vmec_surface_option, verbose, &
       normalized_toroidal_flux_used, safety_factor_q, shat, &
       alpha, zeta, bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0)

!    use read_wout_mod, only: read_wout_file, nfp, lasym, xm, xn, xm_nyq, xn_nyq, ns, phi
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
    ! If number_of_field_periods_to_include is <1, the entire 2*pi toroidal domain will be included.
    ! If number_of_field_periods_to_include is > nfp, an error will result.
    integer, intent(in) :: number_of_field_periods_to_include

    ! The parameter desired_normalized_toroidal_flux determines which flux surface from the VMEC file will be used
    ! for the computation. This parameter should lie in the interval [0,1]
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
    real, dimension(nzgrid*2+1), intent(out) :: zeta

    real, dimension(nalpha, nzgrid*2+1), intent(out) :: bmag, gradpar, gds2, gds21, gds22, gbdrift, gbdrift0, cvdrift, cvdrift0
    
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

    integer :: j, index
    real, dimension(:,:), allocatable :: alpha_2D, zeta_2D
    integer :: ierr, iopen
    real :: dphi, iota, min_dr2, ds, d_pressure_d_s
    real, dimension(2) :: vmec_radial_weight_full, vmec_radial_weight_half
    integer, dimension(2) :: vmec_radial_index_full, vmec_radial_index_half
    real, dimension(:), allocatable :: dr2, normalized_toroidal_flux_full_grid, normalized_toroidal_flux_half_grid
    real, dimension(:), allocatable :: d_pressure_d_s_on_half_grid
    
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

       
!!$       ! Average R and Z from the outermost 2 grid points in vmec's full mesh
!!$       ! to get R and Z on the outermost point of vmec's half mesh:
!!$       weight1 = 0.5_dp
!!$       weight2 = 0.5_dp
!!$       allocate(rmnc_vmecLast(mnmax_vmec),stat=iflag)
!!$       if (iflag .ne. 0) stop 'Allocation error!'
!!$       allocate(zmns_vmecLast(mnmax_vmec),stat=iflag)
!!$       if (iflag .ne. 0) stop 'Allocation error!'
!!$       rmnc_vmecLast = rmnc_vmec(:,ns-1) * weight1 + rmnc_vmec(:,ns) * weight2
!!$       zmns_vmecLast = zmns_vmec(:,ns-1) * weight1 + zmns_vmec(:,ns) * weight2
!!$       
!!$       ! Since the "original" vmec poloidal angle is chosen to have a very condensed
!!$       ! Fourier spectrum, we probably need more Fourier modes to represent the surface using the
!!$       ! straight-field-line coordinate.
!!$       mpol = mpol_vmec*mpol_transform_refinement
!!$       ntor = ntor_vmec*ntor_transform_refinement
!!$       
!!$       ! Beginning of coordinate transformation.
!!$       ! Set up high-resolution grid in the "new" theta coordinate:
!!$       ntheta_coordTransform = mpol * 2 
!!$       nzeta_coordTransform = ntor * 2
!!$       allocate(r_coordTransform(ntheta_coordTransform, nzeta_coordTransform), stat=iflag)
!!$       if (iflag .ne. 0) stop 'Allocation error!'
!!$       allocate(z_coordTransform(ntheta_coordTransform, nzeta_coordTransform), stat=iflag)
!!$       if (iflag .ne. 0) stop 'Allocation error!'
!!$       r_coordTransform = 0
!!$       z_coordTransform = 0
!!$       
!!$       call system_clock(tic1)
!!$       rootSolve_abserr = 1.0e-10_dp
!!$       rootSolve_relerr = 1.0e-10_dp
!!$       !open(unit=5,file="testStraightFieldLines",status='new',form='formatted')
!!$       !write (5,*) ntheta_coordTransform, nzeta_coordTransform
!!$       do izeta = 1,nzeta_coordTransform
!!$          zeta = (izeta-1.0_dp)/nzeta_coordTransform
!!$          do itheta = 1,ntheta_coordTransform
!!$             ! For each value of the new coordinates, solve for the old theta:
!!$             theta_rootSolve_target = (itheta-1.0_dp)/ntheta_coordTransform
!!$             theta_rootSolve_min = theta_rootSolve_target - 0.3
!!$             theta_rootSolve_max = theta_rootSolve_target + 0.3
!!$             
!!$             call fzero(fzero_residual, theta_rootSolve_min, theta_rootSolve_max, theta_rootSolve_target, &
!!$                  rootSolve_relerr, rootSolve_abserr, fzeroFlag)
!!$             ! Note: fzero returns its answer in theta_rootSolve_min
!!$             theta_rootSolve_soln = theta_rootSolve_min
!!$             if (fzeroFlag == 4) then
!!$                stop "ERROR: fzero returned error 4: no sign change in residual"
!!$             else if (fzeroFlag > 2) then
!!$                print *,"WARNING: fzero returned an error code:",fzeroFlag
!!$             end if
!!$             ! Now that we have the old theta, evaluate r and z:
!!$             r_temp = 0
!!$             z_temp = 0
!!$             do imn = 1, mnmax_vmec
!!$                !angle = twopi*(xm_vmec(imn)*theta_rootSolve_soln - xn_vmec(imn)*zeta/nfp)
!!$                angle = xm_vmec(imn)*theta_rootSolve_soln - xn_vmec(imn)*zeta
!!$                r_temp = r_temp + rmnc_vmecLast(imn)*cos(angle)
!!$                z_temp = z_temp + zmns_vmecLast(imn)*sin(angle)
!!$             end do
!!$             r_coordTransform(itheta,izeta) = r_temp
!!$             z_coordTransform(itheta,izeta) = z_temp
!!$             !write(5,*) theta_rootSolve_soln
!!$          end do
!!$       end do
!!$       !close(unit=5)
!!$       call system_clock(toc1)
!!$       print *,"  Time for root solving:",real(toc1-tic1)/countrate

    if (verbose) print *,"Leaving vmec_to_gs2_geometry_interface."

  end subroutine vmec_to_gs2_geometry_interface

  ! --------------------------------------------------------------------------

  function fzero_residual(theta_old)

    use read_wout_mod, only: xm_vmec => xm, xn_vmec => xn, mnmax_vmec => mnmax, lmns, ns
    
    implicit none
    
    real :: theta_old, fzero_residual
    integer :: imn

    ! residual = twopi*(u_new - u_new_target) = (twopi*u_old + lambda) - u_new_target*twopi
    fzero_residual = theta_old - theta_rootSolve_target

    do imn = 1, mnmax_vmec
       fzero_residual = fzero_residual + lmns(imn,ns)*sin(xm_vmec(imn)*theta_old - xn_vmec(imn)*zeta)
    end do

  end function fzero_residual

end module vmec_to_gs2_geometry_interface_mod
