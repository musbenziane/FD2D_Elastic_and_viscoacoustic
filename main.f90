program FD2D_ELASTIC
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

    call countpoints(acqui_rcv,nrcv)
    call countpoints(acqui_src,nshot)

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

    modname_vp  = TRIM(modnameprefix)//'_vp'
    modname_vs  = TRIM(modnameprefix)//'_vs'
    modname_rho = TRIM(modnameprefix)//'_rho'

    allocate(vp2D(nz,nx),vs2D(nz,nx),rho2D(nz,nx))
    allocate(B(nz,nx),M(nz,nx),L(nz,nx))
    allocate(sigma(nz,nx),tau(nz,nx),xi(nz,nx))
    allocate(u(nz,nx),v(nz,nx))
    allocate(src(nt))
    allocate(gz(nz),gx(nx))
    allocate(gzx(nz,nx))
    allocate(Usnap(NINT(REAL(nt/isnap)),nz,nx),Vsnap(NINT(REAL(nt/isnap)),nz,nx))
    allocate(xrcv(nrcv),zrcv(nrcv),xsrc(nshot),zsrc(nshot))
    allocate(SGZ(nt,nrcv),SGX(nt,nrcv))

    call readmodelfiles(nz,nx,modnameprefix,vp2D,vs2D,rho2D)
    call ricker(nt,f0,dt,src)                                      ! Source time function
    call read_acqui_geo(acqui_rcv,nrcv,zrcv,xrcv)
    call read_acqui_geo(acqui_src,nshot,zsrc,xsrc)

    !###############################################
    !### Acquisition geometry and indices stuff ####
    !### Source indices
    isrc = NINT(zsrc / dz )+1
    jsrc = NINT(xsrc / dx )+1
    !### receiver indices
    ircv = NINT(zrcv / dz ) + 1
    jrcv = NINT(xrcv / dx ) + 1


    !###############################################

    B = 1 / rho2D
    M = rho2D * vs2D**2
    L = rho2D * (vp2D**2 - 2*vs2D**2)

    write(*,*)"##########################################"
    write(*,*)"############### CFL Check ################"

     CFL = dt  * maxval(vp2D(:,:)) * SQRT(1/(dx**2) + 1/(dz**2))
    if (CFL > .9) then
        print"( a14,f6.3)"," Courant number is ",CFL
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


    lambdamin = MINVAL(vs2D)/(f0*2.5)
    d         = (/dz,dx/)
    gppw     = lambdamin / maxval(d)
    tol      = 13.9

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

    gz(:)    = 1
    gx(:)    = 1
    gzx(:,:) = 1

    do k=1,gWidth
        taper       = exp(-(attConst*(gWidth-k)/gWidth)**2)
        gz(nz-k+1)  = taper
        gx(k)       = taper
        gx(nx-k+1)  = taper
    end do


    if (FS == 0) then
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

    write(*,*)"##########################################"
    write(*,*)"########## Begin shot loop      ##########"
    write(*,*)"##########################################"
    print*,nshot
    do is=1,nshot

        write(*,*)"##########################################"
        print*,   "########### At shot point ->",is, "/",nshot
        write(*,*)"##########################################"
        k = 0
        do it=1,nt
            if (srctype==1) then
                u(isrc(is),jsrc(is)) = u(isrc(is),jsrc(is)) + dt * src(it)
                v(isrc(is),jsrc(is)) = v(isrc(is),jsrc(is)) + dt * src(it)
            else
                v(isrc(is),jsrc(is)) = v(isrc(is),jsrc(is)) + dt * src(it)
            end if

            !$OMP PARALLEL DO PRIVATE(iz,ix) SHARED(u,v,sigma,xi,tau) SCHEDULE(static)
            do iz=2,nz-1
                do ix=2,nx-1

                    u(iz,ix)     = u(iz,ix) +      B(iz,ix) * (dt/dx) * (sigma(iz, ix) - sigma(iz,ix-1)) + &
                                                   B(iz,ix) * (dt/dz) * (xi(iz, ix) - xi(iz-1,ix))

                    v(iz,ix)     = v(iz,ix) +      B(iz,ix) * (dt/dx) * (xi(iz, ix) - xi(iz,ix-1)) +  &
                                                   B(iz,ix) * (dt/dz) * (tau(iz, ix) - tau(iz-1,ix))

                    sigma(iz,ix) = sigma(iz,ix) + (L(iz,ix) + 2*M(iz,ix)) * (dt/dx) * (u(iz,ix+1) - u(iz,ix)) + &
                                                   L(iz,ix) * (dt/dz) * (v(iz+1,ix) - v(iz,ix))

                    tau(iz,ix)   = tau(iz,ix)   + (L(iz,ix) + 2*M(iz,ix)) * (dt/dz) * (v(iz+1,ix) - v(iz,ix)) + &
                                                   L(iz,ix) * (dt/dx) * (u(iz,ix+1) - u(iz,ix))

                    xi(iz,ix)    = xi(iz,ix)    +  M(iz,ix) * (dt/dz) * (u(iz+1,ix) - u(iz,ix))  + &
                                                   M(iz,ix) * (dt/dx) * (v(iz,ix+1) - v(iz,ix))
                end do
            end do
            !$OMP END PARALLEL DO

            if (FS==1) then
                tau(1,:) = -tau(3,:);
                xi(1,:)  = -xi(3,:);
                u(1,:)   = u(3,:);
                v(1,:)   = v(3,:);
            end if

            if (ABS==1) then
                tau   = tau   * gzx
                sigma = sigma * gzx
                xi    = xi    * gzx
                u     = u     * gzx
                v     = v     * gzx
            end if


            if (mod(it,isnap) == 0) then
                k = k + 1
                Usnap(k,:,:) = u
                Vsnap(k,:,:) = v
            end if
            if (mod(it,NINT(nt/20.)) == 0) then
                print*, "########### At time sample ->",it, "/",nt
            end if

            !##### Extract displacement at receiver point

            do i=1,nrcv
                SGZ(it,i) = v(ircv(i),jrcv(i))
                SGX(it,i) = u(ircv(i),jrcv(i))
            end do

        end do

        write(isstr,"(I5.5)") is
        outname_u      = "OUTPUT/SnapX_"//TRIM(isstr)
        outname_v      = "OUTPUT/SnapZ_"//TRIM(isstr)

        inquire(iolength=reclsnaps) Vsnap

        open(3,file=outname_u,access="direct",recl=reclsnaps)
        write(3,rec=1) Usnap
        close(3)

        open(4,file=outname_v,access="direct",recl=reclsnaps)
        write(4,rec=1) Vsnap
        close(4)


        outname_u      = "OUTPUT/SeismX_"//TRIM(isstr)
        outname_v      = "OUTPUT/SeismZ_"//TRIM(isstr)

        inquire(iolength=reclsnaps) SGZ

        open(5,file=outname_u,access="direct",recl=reclsnaps)
        write(5,rec=1) SGX
        close(5)
        open(6,file=outname_v,access="direct",recl=reclsnaps)
        write(6,rec=1) SGZ
        close(6)

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

