program FD2D_ELASTIC
    !####################################################################################################################
    ! FD2D: Finite difference solver for the elastic wave equation in cartesian coordinates with flat free surface.
    ! NOTE: This program has not been validated.
    ! 
    ! Language: Fortran 90, with parralel impelementation using OpenMP API
    ! 
    ! Author: Mus Benziane
    ! Source: Virieux 1986
    !         other sources [2] are mentioned as comments.
    !
    ! Input file: [Example]
    ! testing              ! prefix of model file names, suffixes: _vp _vs _rho
    ! acqui_src            ! sources positions ascii file name
    ! acqui_rcv            ! receiver positions ascii file name
    ! 1996                 ! zmax
    ! 3196                 ! xmax
    ! 500                  ! nz
    ! 800                  ! nx
    ! 4                    ! dz
    ! 4                    ! dx
    ! 0.0004               ! dt
    ! 4000                 ! nt
    ! 6.                   ! Wavelet's peak frequency
    ! 50                   ! Snapshot interval
    ! 1                    ! [0/1] explosive, otherwise vertical
    ! 1                    ! [0/1] Free Surface - if set to 0, rigid BC are used, i.e no displacement at boundary
    ! 1                    ! [0/1] Absorbing
    ! .65                  ! Att constant for sponge layer
    ! 70                   ! Sponge layer width
    !
    ! -> Acquisition file in ASCII Format: 
    !  1) receivers ASCII file [N: Receivers number]
    !
    !  z1 x1
    !  z2 x2
    !  .  .
    !  zN xN
    !
    !  2) Shot points ASCII file [M: Shot points number]
    !
    !  z1 x1
    !  z2 x2
    !  .  .
    !  zM xM
    ! 
    ! -> Model files in C-Style binary floats [doubles]: Vp, Vs, Rho files are needed.
    !                                                  : For simple models, use create2Dmodel_files.f90
    ! 
    ! -> Outputs are created in OUTPUT/ if OUTPUT/ is not created by the user, the program will not handle it.
    ! Output files in OUTPUT directory:
    !
    !####################################################################################################################


    !$ USE OMP_LIB
    use read_acqui

    implicit none
    real(kind=8)                                    :: dz, dx, dt, f0, xmax, zmax, lambdamin, gppw, tol
    real(kind=8)                                    :: CFL, d(2), attConst, taper, time, t_cpu_0, t_cpu_1, t_cpu
    real(kind=8), dimension(:),     allocatable     :: src, gz, gx, zrcv, xrcv, zsrc, xsrc
    real(kind=8), dimension(:,:),   allocatable     :: vp2D, vs2D, rho2D, B, M, L, sigma, tau, xi, u, v, gzx, SGX, SGZ
    real(kind=8), dimension(:,:,:), allocatable     :: Usnap, Vsnap
    integer                                         :: nz,nx, isnap, nt, it, iz, ix, is, k, srctype,i
    integer                                         :: gWidth, FS, ABS, reclsnaps, n_workers, ir, t0, t1, nshot, nrcv
    integer, dimension(:), allocatable              :: isrc, jsrc, ircv, jrcv
    character(len=40)                               :: filename, outname_u, outname_v, filecheck, modnameprefix
    character(len=40)                               :: modname_vp, modname_vs, modname_rho, isstr, acqui_src, acqui_rcv


    write(*,*) "##########################################"
    write(*,*) "############### OpenMP     ###############"

    !$OMP PARALLEL
    !$ n_workers = OMP_GET_NUM_THREADS()
    !$OMP END PARALLEL
    !$ print '(3X,"Number of workers ->  ",i2)',n_workers

    call cpu_time(t_cpu_0)
    call system_clock(count=t0, count_rate=ir)

    filename       = "parameters.in"

    write(*,*) "##########################################"
    write(*,*) "######## Reading parameters file #########"
    write(*,*) "##########################################"

    print*,"Is the parameters input file (parameters.in) [Yes/no]"
    read(*,*) filecheck

    if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
            filecheck=="Y") then
        write(*,*) "Reading simulation parameters..."

    elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
            filecheck=="N") then
        write(*,*) "Enter simulation parameters text file name with extension"
        write(*,*) "40 characters max"
        read(*,*) filename

    else
        write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
        write(*,*) "The program have been terminated, please star over"
        stop
    end if

    open (2, file=filename, status = 'old')
    read(2,*) modnameprefix
    read(2,*) acqui_src
    read(2,*) acqui_rcv
    read(2,*) zmax
    read(2,*) xmax
    read(2,*) nz
    read(2,*) nx
    read(2,*) dz
    read(2,*) dx
    read(2,*) dt
    read(2,*) nt
    read(2,*) f0
    read(2,*) isnap
    read(2,*) srctype
    read(2,*) FS
    read(2,*) ABS
    read(2,*) attConst
    read(2,*) gWidth
    close(2)

    call countpoints(acqui_rcv,nrcv)                    ! Count number of receivers
    call countpoints(acqui_src,nshot)                   ! Count number of shot points

    100 format (A,F6.1,X,F6.1)
    120 format (A,F6.1)
    140 format (A,I4,X,I4)
    160 format (A,I4)

    write(*,100)  "Maximum distance zmax, xmax    -> ",zmax, xmax
    write(*,140)  "Grid size nz, nx               -> ",nz, nx
    write(*,100)  "Grid points spacing dz, dx     -> ",dz, dx
    write(*,120)  "time step                      -> ",dt
    write(*,160)  "Number of time steps           -> ",nt
    write(*,120)  "Ricker's peak frequency        -> ",f0
    write(*,160)  "Snapshot interval              -> ",isnap
    write(*,160)  "Number Shot points             -> ",nshot
    write(*,160)  "Number of receivers            -> ",nrcv
    write(*,160)  "Free Surface BC [1/0]          -> ",FS
    write(*,160)  "Absorbing BC [1/0]             -> ",ABS


    modname_vp  = TRIM(modnameprefix)//'_vp'
    modname_vs  = TRIM(modnameprefix)//'_vs'
    modname_rho = TRIM(modnameprefix)//'_rho'

    allocate(vp2D(nz,nx),vs2D(nz,nx),rho2D(nz,nx))     ! Model parameters matrices
    allocate(B(nz,nx),M(nz,nx),L(nz,nx))               ! Model parameters matrices
    allocate(sigma(nz,nx),tau(nz,nx),xi(nz,nx))        ! Stress tensor parameters matrices
    allocate(u(nz,nx),v(nz,nx))                        ! Solution parameters matrices
    allocate(src(nt))                                  ! Source time function
    allocate(gz(nz),gx(nx))                            ! Abosrbing boundary condition: tmp matrices
    allocate(gzx(nz,nx))                               ! Abosrbing boundary condition: Grid
    allocate(Usnap(NINT(REAL(nt/isnap)),nz,nx), &      ! Wavefield snapshot matrice for binary output ...
             Vsnap(NINT(REAL(nt/isnap)),nz,nx))        ! Fastest dimension: time steps, then nz, lastly, nx
    allocate(xrcv(nrcv),zrcv(nrcv),xsrc(nshot),zsrc(nshot))
    allocate(SGZ(nt,nrcv),SGX(nt,nrcv))                ! Seismograms <<Shot Gathers>> matrices for binary output

    call readmodelfiles(nz,nx,modnameprefix,vp2D,vs2D,rho2D)       ! Subroutine to read model binaries
    call ricker(nt,f0,dt,src)                                      ! Source time function
    call read_acqui_geo(acqui_rcv,nrcv,zrcv,xrcv)                  ! Read acquisition geometry from ASCII File: receivers
    call read_acqui_geo(acqui_src,nshot,zsrc,xsrc)                 ! Same                                     : Shot points

    !###############################################
    !### Acquisition geometry and indices stuff ####
    !### Source indices
    isrc = NINT(zsrc / dz) + 1
    jsrc = NINT(xsrc / dx) + 1
    !### receiver indices
    ircv = NINT(zrcv / dz) + 1
    jrcv = NINT(xrcv / dx) + 1
    !###############################################

    
    !###############################################
    !### Elastic parameters                     ####
    B = 1 / rho2D                       
    M = rho2D * vs2D**2
    L = rho2D * (vp2D**2 - 2*vs2D**2)
    !###############################################


    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"

     CFL = dt  * maxval(vp2D(:,:)) * SQRT(1/(dx**2) + 1/(dz**2))
    if (CFL > .6) then
        print"(a14,f6.3)"," Courant number is ",CFL
        print*,"Decrease time step, the program has been terminated"
        stop
    else
        print"(a14,f6.3)","Courant number is ",CFL
        print*,"Simulation is stable"
    end if
    write(*,*)"##########################################"


    write(*,*)"##########################################"
    write(*,*)"########## Space Sampling check ##########"
    write(*,*)"##########################################"

    lambdamin = MINVAL(vs2D)/(f0*2.5)    ! Shortest wavelength
    d         = (/dz,dx/)
    gppw     = lambdamin / maxval(d)     ! Grid points per wavelength
    tol      = 18                        ! Tolerance set to 25 Grid points / Wavelength

    print"(a32,f6.1)", " Grid points per minimum wavelength ->", gppw

    if ((gppw)<tol) then
        print*,"Grid spacing size is too large"
        print*,"Numerical dispersion might be present"
        write(*,*) "Reduce spacing (dz or dx, whichever the largest) size or frequency"
        write(*,*) "The program has been terminated"
        stop
    else
        print*, "Spatial sampling as OK!"
    end if


    !###############################################
    !### Absorbing boundary condition grid prep ####
    !### See Cerjan 1985
    gz(:)    = 1
    gx(:)    = 1
    gzx(:,:) = 1

    do k=1,gWidth
        taper       = exp(-(attConst*(gWidth-k)/gWidth)**2)
        gz(nz-k+1)  = taper
        gx(k)       = taper
        gx(nx-k+1)  = taper
    end do

    if (FS == 0) then ! If free surface boundary condtions, sponge layers won't be used below the surface
        do k=1,gWidth
            taper       = exp(-(attConst*(gWidth-k)/20)**2)
            gz(k)       = taper
        end do
    end if

    do ix=1,nx
        do iz=1,nz
            gzx(iz,ix) = gz(iz) * gx(ix)
        end do
    end do
    !###############################################

    if (FS==1) then  ! If free surface boundary conditions, source is moved to free surface
        isrc(:) = 3
        ircv(:) = 3
    end if


    write(*,*)"##########################################"
    write(*,*)"########## Begin shot loop      ##########"
    write(*,*)"##########################################"
    do is=1,nshot

        write(*,*)"##########################################"
        print*,   "########### At shot point ->",is, "/",nshot
        write(*,*)"##########################################"
        k = 0
        do it=1,nt
            !###############################################
            !### Begin time loop                       #####
            !###############################################

            if (srctype==1) then    ! Explosive source
                u(isrc(is),jsrc(is)) = u(isrc(is),jsrc(is)) + dt * src(it)
                v(isrc(is),jsrc(is)) = v(isrc(is),jsrc(is)) + dt * src(it)
            else                    ! Vertical source
                v(isrc(is),jsrc(is)) = v(isrc(is),jsrc(is)) + dt * src(it)
            end if

            !$OMP PARALLEL DO PRIVATE(iz,ix) SHARED(u,v,sigma,xi,tau) SCHEDULE(static)
            do iz=3,nz-2
                do ix=3,nx-2
                    ! Staggered grids finite difference scheme << Velocity-Stress formulation >>
                    u(iz,ix)     = u(iz,ix)     +  B(iz,ix) * (dt/dx) * (9./8. * (sigma(iz, ix+1) - sigma(iz,ix))    & 
                                                                       -1./24. * (sigma(iz, ix+2) - sigma(iz,ix-1))) &  
                                                +  B(iz,ix) * (dt/dz) * (9./8. * (xi(iz+1, ix)    - xi(iz,ix))       & 
                                                                       -1./24. * (xi(iz+2, ix)    - xi(iz-1,ix)))

                    v(iz,ix)     = v(iz,ix)     +  B(iz,ix) * (dt/dx) * (9./8. * (xi(iz, ix+1)    - xi(iz,ix))       &
                                                                       -1./24. * (xi(iz, ix+2)    - xi(iz,ix-1)))    &
                                                +  B(iz,ix) * (dt/dz) * (9./8. * (tau(iz+1, ix)   - tau(iz,ix))      &
                                                                       -1./24. * (tau(iz+2, ix)   - tau(iz-1,ix)))                                              

                    sigma(iz,ix) = sigma(iz,ix) + (L(iz,ix) + 2*M(iz,ix)) * (dt/dx) * (9./8. * (u(iz,ix)   - u(iz,ix-1))  &
                                                                                     -1./24. * (u(iz,ix+1) - u(iz,ix-2))) &
                                                +  L(iz,ix)               * (dt/dz) * (9./8. * (v(iz,ix)   - v(iz-1,ix))  &
                                                                                     -1./24. * (v(iz+1,ix) - v(iz-2,ix)))

                    tau(iz,ix)   = tau(iz,ix)   + (L(iz,ix) + 2*M(iz,ix)) * (dt/dz) * (9./8. * (v(iz,ix)   - v(iz-1,ix))  &
                                                                                     -1./24. * (v(iz+1,ix) - v(iz-2,ix))) &
                                                +                L(iz,ix) * (dt/dx) * (9./8. * (u(iz,ix)   - u(iz,ix-1))  &
                                                                                     -1./24. * (u(iz,ix+1) - u(iz,ix-2)))

                    xi(iz,ix)    = xi(iz,ix)    +  M(iz,ix)               * (dt/dz) * (9./8. * (u(iz,ix)   - u(iz-1,ix))  &
                                                                                     -1./24. * (u(iz+1,ix) - u(iz-2,ix))) &        
                                                +                M(iz,ix) * (dt/dx) * (9./8. * (v(iz,ix)   - v(iz,ix-1))  &
                                                                                     -1./24. * (v(iz,ix+1) - v(iz,ix-2)))
                end do
            end do
            !$OMP END PARALLEL DO
            
            !##################################################
            !### Free Surface Boundary conditions      ########
            !### See: Moczo, P., Kristek, J., & Gális, M. 
            ! (2014). The Finite-Difference Modelling of 
            ! Earthquake Motions: Waves and Ruptures. Cambridge: 
            ! Cambridge University Press. 
            ! doi:10.1017/CBO9781139236911     
            ! Chapter 7.5.1: << Stress Imaging >> 
            ! This is a trial for planar FS boundary conditions
            if (FS==1) then
                tau(1,:) = -tau(5,:)
                xi(1,:)  = -xi(5,:)
                u(1,:)   =  u(5,:)
                v(1,:)   =  v(5,:)

                tau(2,:) = -tau(4,:)
                xi(2,:)  = -xi(4,:)
                u(2,:)   =  u(4,:)
                v(2,:)   =  v(4,:)
            end if
            !##################################################


            !##################################################
            !### Abosrbing Boundary conditions         ########
            if (ABS==1) then
                tau   = tau   * gzx
                sigma = sigma * gzx
                xi    = xi    * gzx
                u     = u     * gzx
                v     = v     * gzx
            end if
            !##################################################

            ! Get wavefield snapshot
            if (mod(it,isnap) == 0) then
                k = k + 1
                Usnap(k,:,:) = u
                Vsnap(k,:,:) = v
            end if
            
            ! Display progress
            if (mod(it,NINT(nt/10.)) == 0) then
                print*, "########### At time sample ->",it, "/",nt
            end if

            !##### Extract displacement at receiver point
            do i=1,nrcv
                SGZ(it,i) = v(ircv(i),jrcv(i))
                SGX(it,i) = u(ircv(i),jrcv(i))
            end do
        end do ! time loop ends here
                
        write(isstr,"(I5.5)") is
        outname_u      = "OUTPUT/SnapX_"//TRIM(isstr)
        outname_v      = "OUTPUT/SnapZ_"//TRIM(isstr)

        inquire(iolength=reclsnaps) Vsnap
        
        ! ### Write solution binary: Snapshot
        open(3,file=outname_u,access="direct",recl=reclsnaps)
        write(3,rec=1) Usnap
        close(3)

        open(4,file=outname_v,access="direct",recl=reclsnaps)
        write(4,rec=1) Vsnap
        close(4)


        outname_u      = "OUTPUT/SeismX_"//TRIM(isstr)
        outname_v      = "OUTPUT/SeismZ_"//TRIM(isstr)

        inquire(iolength=reclsnaps) SGZ

        ! ### Write solution binary: Seismograms
        open(5,file=outname_u,access="direct",recl=reclsnaps)
        write(5,rec=1) SGX
        close(5)

        open(7,file=outname_v,access="direct",recl=reclsnaps)
        write(7,rec=1) SGZ
        close(7)

    end do  ! shot loop ends here

    call system_clock(count=t1, count_rate=ir)
    time = real(t1 - t0,kind=8) / real(ir,kind=8)

    call cpu_time(t_cpu_1)
    t_cpu = t_cpu_1 - t_cpu_0

    write(*,*) "##########################################"
    write(*,*) "######### TIME:                 ##########"
    print '(//3X,"Elapsed Time        : ",1PE10.3," [s]",/ &
            &,3X,"CPU Time            : ",1PE10.3," [s]",//)', &
            & time,t_cpu
    write(*,*) "##########################################"

    deallocate(vp2D,vs2D,rho2D)
    deallocate(B,M,L)
    deallocate(sigma,tau,xi)
    deallocate(u,v)
    deallocate(src)
    deallocate(gz,gx)
    deallocate(gzx)
    deallocate(Usnap,Vsnap)
    deallocate(xrcv,zrcv,xsrc,zsrc)
    deallocate(SGZ,SGX)
end program