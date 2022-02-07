program FD2D_VISCOACOUSTIC
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
    real(kind=8)                                    :: dz, dt, f0, lambdamin, gppw, tol, c1, c2, c3, c4, wPML, sigma_m,temp
    real(kind=8)                                    :: CFL,  time, t_cpu_0, t_cpu_1, t_cpu
    real(kind=8), dimension(:),     allocatable     :: src, sigmaz, sigmax, zrcv, xrcv, zsrc, xsrc, tmpS, tmpE
    real(kind=8), dimension(:,:),   allocatable     :: c2D, qf2D, rho2D, K, tau_sigma,tau_epsi, p, pz, px, w, u, &
                                                       ctmp, qftmp, rhotmp, seis_p, seis_w, seis_u, r
    real(kind=8), dimension(:,:,:), allocatable     :: Wsnap, Usnap, Psnap
    integer                                         :: nz,nx, isnap, nt, it, iz, ix, is, i,npml,c
    integer                                         :: is_FS, is_PML, reclsnaps, n_workers, ir, t0, t1, nshot, nrcv
    integer, dimension(:), allocatable              :: isrc, jsrc, ircv, jrcv
    character(len=40)                               :: filename, outname_p, outname_w, outname_u, filecheck, modnameprefix
    character(len=40)                               :: modname_c, modname_qf, modname_rho, isstr, acqui_src, acqui_rcv


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
    read(2,*) nz
    read(2,*) nx
    read(2,*) dz
    read(2,*) dt
    read(2,*) nt
    read(2,*) f0
    read(2,*) isnap
    read(2,*) is_FS
    read(2,*) is_PML
    read(2,*) npml
    close(2)

    call countpoints(acqui_rcv,nrcv)                    ! Count number of receivers
    call countpoints(acqui_src,nshot)                   ! Count number of shot points

    100 format (A,F6.1,X,F6.1)
    120 format (A,F6.1)
    140 format (A,I4,X,I4)
    160 format (A,I4)

    write(*,100)  "Maximum distance zmax, xmax    -> ",(nz-1)*dz, (nx-1)*dz
    write(*,140)  "Grid size nz, nx               -> ",nz, nx
    write(*,100)  "Grid points spacing dz, dx     -> ",dz
    write(*,120)  "time step                      -> ",dt
    write(*,160)  "Number of time steps           -> ",nt
    write(*,120)  "Ricker's peak frequency        -> ",f0
    write(*,160)  "Snapshot interval              -> ",isnap
    write(*,160)  "Number Shot points             -> ",nshot
    write(*,160)  "Number of receivers            -> ",nrcv
    write(*,160)  "Free Surface BC [1/0]          -> ",is_FS
    write(*,160)  "PML  [1/0]                     -> ",is_PML


    modname_c   = TRIM(modnameprefix)//'_c'
    modname_qf  = TRIM(modnameprefix)//'_qf'
    modname_rho = TRIM(modnameprefix)//'_rho'

    allocate(ctmp(nz,nx),qftmp(nz,nx),rhotmp(nz,nx))   ! Model parameters matrices
    allocate(src(nt))                                  ! Source time function
    allocate(xrcv(nrcv),zrcv(nrcv),xsrc(nshot),zsrc(nshot))

    call readmodelfiles(nz,nx,modnameprefix,ctmp,qftmp,rhotmp)     ! Subroutine to read model binaries
    call ricker(nt,f0,dt,src)                                      ! Source time function
    call read_acqui_geo(acqui_rcv,nrcv,zrcv,xrcv)                  ! Read acquisition geometry from ASCII File: receivers
    call read_acqui_geo(acqui_src,nshot,zsrc,xsrc)                 ! Same                                     : Shot points




    !###############################################
    !### Acquisition geometry and indices stuff ####
    !### Source indices
    allocate(isrc(nshot),jsrc(nshot),ircv(nrcv),jrcv(nrcv))
    isrc = NINT(zsrc / dz) + 1
    jsrc = NINT(xsrc / dz) + 1
    !### receiver indices
    ircv = NINT(zrcv / dz) + 1
    jrcv = NINT(xrcv / dz) + 1
    !###############################################

    !###############################################
    !#####       Extend model for PML           ####
    if (is_PML==1) then
        allocate(c2D(nz+npml,nx+2*npml))
        allocate(qf2D(nz+npml,nx+2*npml))
        allocate(rho2D(nz+npml,nx+2*npml))

        c2D(1:nz,npml+1:nx+npml)   = ctmp
        qf2D(1:nz,npml+1:nx+npml)  = qftmp
        rho2D(1:nz,npml+1:nx+npml) = rhotmp

        do i=1,npml
            c2D(:,i)          = ctmp(:,1)
            c2D(:,nx+npml+i)  = ctmp(:,nx) 

            qf2D(:,i)         = qftmp(:,1)
            qf2D(:,nx+npml+i) = qftmp(:,nx) 

            rho2D(:,i)        = rhotmp(:,1)
            rho2D(:,nx+npml+i)= rhotmp(:,nx) 
        end do


        do i=1,npml
            c2D(nz+i,:)      = c2D(nz,:)
            qf2D(nz+i,:)     = qf2D(nz,:)
            rho2D(nz+i,:)    = rho2D(nz,:)
        end do

        ! Adjust nz, nx for PML grid points.
        nz = nz + npml
        nx = nx + 2*npml
    else 
        allocate(c2D(nz,nx))
        allocate(qf2D(nz,nx))
        allocate(rho2D(nz,nx))
        c2D   = ctmp
        qf2D  = qftmp
        rho2D = rhotmp
    end if
    !###############################################

    
    !###############################################
    !#######         Model parameters       ########
    
    allocate(K(nz,nx),tau_sigma(nz,nx),tau_epsi(nz,nx))

    K         = c2D**2*rho2D
    tau_sigma = (sqrt(1. + 1./qf2D**2)-1./qf2D)/f0
    tau_epsi  = 1./(f0**2*tau_sigma)
    !###############################################




    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"

    CFL = dt  * maxval(c2D(:,:)) * SQRT(1/(dz**2) + 1/(dz**2))
    if (CFL > .34) then
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

    lambdamin = MINVAL(c2D)/(f0*2.5)     ! Shortest wavelength
    gppw     = lambdamin / dz            ! Grid points per wavelength
    tol      = 12                        ! Tolerance set to 12 Grid points / Wavelength

    print"(a32,f6.1)", " Grid points per minimum wavelength ->", gppw

    if ((gppw)<tol) then
        print*,"Grid spacing size is too large"
        print*,"Numerical dispersion might be present"
        write(*,*) "Reduce spacing (dz or dz, whichever the largest) size or frequency"
        write(*,*) "The program has been terminated"
        stop
    else
        print*, "Spatial sampling as OK!"
    end if

    !###############################################
    !#######         Sentcils coefficient   ########
    c1 = 1225.0/1024.0 
    c2 =-245.0 /3072.0
    c3 = 49.0  /5120.0   
    c4 = -5.0  /7168.0
    !###############################################
    
    !###############################################
    !#######            PML profile         ########
    allocate(tmpS(npml),tmpE(npml),sigmaz(nz),sigmax(nx))
    sigma_m   = 450
    wPML      = dz * (npml-1)
    sigmaz(:) = 0
    sigmax(:) = 0
    do i=1,npml
        tmpS(i) = sigma_m * ((npml-i)*dz/wPML)**2
        tmpE(i) = sigma_m * ((i-1)*dz/wPML)**2
    end do 
    
    sigmaz(nz-npml+1:nz) = tmpE
    sigmax(nx-npml+1:nx) = tmpE
    sigmax(1:npml)       = tmpS

    !###############################################

    allocate(w(nz,nx),u(nz,nx),r(nz,nx))               ! Solution  matrices
    allocate(p(nz,nx),pz(nz,nx),px(nz,nx))
    allocate(Wsnap(NINT(REAL(nt/isnap)),nz,nx),    &   ! Wavefield snapshot matrice for binary output ...
             Usnap(NINT(REAL(nt/isnap)),nz,nx),    &
             Psnap(NINT(REAL(nt/isnap)),nz,nx))        ! Fastest dimension: time steps, then nz, lastly, nx
    allocate(seis_p(nt,nrcv),seis_w(nt,nrcv),seis_u(nt,nrcv))


    write(*,*)"##########################################"
    write(*,*)"########## Begin shot loop      ##########"
    write(*,*)"##########################################"
    
    do is=1,nshot
        write(*,*)"##########################################"
        print*,   "########### At shot point ->",is, "/",nshot
        write(*,*)"##########################################"
        c = 0 
        !###############################################
        !### Begin time loop                       #####
        !###############################################
        
        do it=1,nt
            do iz=5,nz-4
                do ix=5,nx-4
                    temp      = K(iz,ix)*dt*(tau_epsi(iz,ix)/tau_sigma(iz,ix)-1.0)/dz/tau_sigma(iz,ix)
                    r(iz,ix)  = (1.0-dt/tau_sigma(iz,ix))*r(iz,ix)-temp                 * &
                                (c1*(u(iz,ix+1)-u(iz,ix)+w(iz,ix)-w(iz-1,ix))           + &
                                 c2*(u(iz,ix+2)-u(iz,ix-1)+w(iz+1,ix)-w(iz-2,ix))       + &
                                 c3*(u(iz,ix+3)-u(iz,ix-2)+w(iz+2,ix)-w(iz-3,ix))       + &
                                 c4*(u(iz,ix+4)-u(iz,ix-3)+w(iz+3,ix)-w(iz-4,ix)))


                    temp      = K(iz,ix)*dt*tau_epsi(iz,ix)/dz/tau_sigma(iz,ix)
                    pz(iz,ix) = (1-dt*sigmaz(iz)) * pz(iz,ix)-r(iz,ix)*dt-temp          * &    
                                (c1*(w(iz,ix)-w(iz-1,ix))                               + &
                                 c2*(w(iz+1,ix)-w(iz-2,ix))                             + &
                                 c3*(w(iz+2,ix)-w(iz-3,ix))                             + &
                                 c4*(w(iz+3,ix)-w(iz-4,ix)))

                    px(iz,ix) = (1-dt*sigmax(ix)) * px(iz,ix)-r(iz,ix)*dt-temp          * &    
                                (c1*(u(iz,ix+1)-u(iz,ix))                               + &
                                 c2*(u(iz,ix+2)-u(iz,ix-1))                             + &
                                 c3*(u(iz,ix+3)-u(iz,ix-2))                             + &
                                 c4*(u(iz,ix+4)-u(iz,ix-3)))
     
                    p(iz,ix)  = pz(iz,ix) + px(iz,ix)     
                    

                    temp      = dt/(rho2D(iz,ix)*dz) 
                    w(iz,ix)  = (1-dt*sigmaz(iz)) * w(iz,ix)-temp                       * &
                                  (c1*(p(iz+1,ix)-p(iz,ix))                             + &
                                   c2*(p(iz+2,ix)-p(iz-1,ix))                           + &  
                                   c3*(p(iz+3,ix)-p(iz-2,ix))                           + &
                                   c4*(p(iz+4,ix)-p(iz-3,ix)))

                    u(iz,ix) = (1-dt*sigmaz(ix)) * u(iz,ix) - temp                      * &
                                    (c1*(p(iz,ix)-p(iz,ix-1))                           + &
                                     c2*(p(iz,ix+1)-p(iz,ix-2))                         + &  
                                     c3*(p(iz,ix+2)-p(iz,ix-3))                         + &
                                     c4*(p(iz,ix+3)-p(iz,ix-4)))
        
                end do
            end do

            p(isrc(is),jsrc(is))  = p(isrc(is),jsrc(is)) + dt * src(it) 

            if (is_FS==1) then
                p(5,:)  = 0;
                pz(5,:) = 0;
                px(5,:) = 0;
            end if

            do ir=1,nrcv
                seis_p(it,ir) = p(ircv(ir),jrcv(ir))
                seis_w(it,ir) = w(ircv(ir),jrcv(ir))
                seis_u(it,ir) = u(ircv(ir),jrcv(ir))
            end do

            if (mod(it,isnap) == 0) then
                Psnap(c,:,:) = p
                Wsnap(c,:,:) = w
                Usnap(c,:,:) = u
                c = c + 1
            end if

                    ! Display progress
            if (mod(it,NINT(nt/10.)) == 0) then
                print*, "########### At time sample ->",it, "/",nt
            end if
        end do 


        write(isstr,"(I5.5)") is
        outname_p      = "OUTPUT/field_P_"//TRIM(isstr)
        outname_w      = "OUTPUT/field_W_"//TRIM(isstr)
        outname_u      = "OUTPUT/field_U_"//TRIM(isstr)

        inquire(iolength=reclsnaps) Psnap

        ! ### Write solution binary: Snapshot
        open(3,file=outname_p,access="direct",recl=reclsnaps)
        write(3,rec=1) Psnap
        close(3)

        open(4,file=outname_w,access="direct",recl=reclsnaps)
        write(4,rec=1) Wsnap
        close(4)


        open(5,file=outname_u,access="direct",recl=reclsnaps)
        write(5,rec=1) Usnap
        close(5)


        outname_p      = "OUTPUT/seism_P_"//TRIM(isstr)
        outname_w      = "OUTPUT/Seism_W_"//TRIM(isstr)
        outname_u      = "OUTPUT/Seism_U_"//TRIM(isstr)
        

        inquire(iolength=reclsnaps) seis_p

        ! ### Write solution binary: Seismograms
        open(6,file=outname_p,access="direct",recl=reclsnaps)
        write(6,rec=1) seis_p
        close(6)

        open(7,file=outname_w,access="direct",recl=reclsnaps)
        write(7,rec=1) seis_w
        close(7)

        open(8,file=outname_w,access="direct",recl=reclsnaps)
        write(8,rec=1) seis_w
        close(8)

    end do



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

end program