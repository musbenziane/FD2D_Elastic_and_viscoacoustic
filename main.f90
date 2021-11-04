program FD2D_ELASTIC
    implicit none
    real(kind=8)                                    :: dz, dx, dt, f0, xmax, zmax, xsrc, lambdamin, gppw,  &
                                                        CFL, d(2), attConst, taper
    real(kind=8), dimension(:), allocatable         :: src, gz, gx
    real(kind=8), dimension(:,:), allocatable       :: vp2D, vs2D, rho2D, B, M, L, sigma, tau, xi, u, v, gzx
    real(kind=8), dimension(:,:,:), allocatable     :: Usnap, Vsnap
    integer                                         :: nz,nx, isnap, nt, tol, it, iz, ix, k, srctype, isrc, jsrc, &
                                                       gWidth, FS, ABS, reclsnaps
    character(len=40)                               :: filename, outname_u, outname_v, filecheck, modnameprefix
    character(len=40)                               :: modname_vp, modname_vs, modname_rho

    filename       = "parameters.in"
    outname_u      = "OUTPUT/snapshots_u.bin"
    outname_v      = "OUTPUT/snapshots_v.bin"

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
    read(2,*) zmax
    read(2,*) xmax
    read(2,*) nz
    read(2,*) nx
    read(2,*) dz
    read(2,*) dx
    read(2,*) dt
    read(2,*) nt
    read(2,*) f0
    read(2,*) jsrc
    read(2,*) isrc
    read(2,*) isnap
    read(2,*) srctype
    read(2,*) FS
    read(2,*) ABS
    read(2,*) attConst
    read(2,*) gWidth
    close(2)

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
    write(*,120)  "Source position                -> ",xsrc
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

    call readmodelfiles(nz,nx,modnameprefix,vp2D,vs2D,rho2D)
    call ricker(nt,f0,dt,src)                                      ! Source time function

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
    tol      = 15;

    print"(a32,f6.1)", " Grid points per minimum wavelength ->", gppw

    if ((gppw)<tol) then
        print*,"Grid spacing size is too large"
        print*,"Numerical dispersion might be present"
        print*,"Do you wish to continue anyways? Yes/no"
        read*, filecheck

        if (filecheck=="Yes" .or. filecheck=="yes" .or. filecheck=="y" .or. &
                filecheck=="Y") then
            print*,"Proceeding..."

        elseif  (filecheck=="No" .or. filecheck=="no" .or. filecheck=="n" .or. &
                filecheck=="N") then
            write(*,*) "Reduce spacing (dz or dx, whichever the largest) size or frequency"
            write(*,*) "The program has been terminated"
            stop
        else
            write(*,*) "Only: Yes/yes/Y/y & No/no/N/n are handled"
            write(*,*) "The program have been terminated, please star over"
            stop
        end if
    else
        print*, "Spatial sampling as OK!"
    end if


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

    k = 0
    do it=1,nt
        if (srctype==1) then
            u(isrc,jsrc) = u(isrc,jsrc) + dt * src(it)
            v(isrc,jsrc) = v(isrc,jsrc) + dt * src(it)
        else
            v(isrc,jsrc) = v(isrc,jsrc) + dt * src(it)
        end if

        do iz=2,nz-1
            do ix=2,nx-1
                u(iz,ix)     = u(iz,ix) + B(iz,ix) * (dt/dx) * (sigma(iz, ix) - sigma(iz,ix-1)) + &
                                          B(iz,ix) * (dt/dz) * (xi(iz, ix) - xi(iz-1,ix))

                v(iz,ix)     = v(iz,ix) + B(iz,ix) * (dt/dx) * (xi(iz, ix) - xi(iz,ix-1)) +  &
                                          B(iz,ix) * (dt/dz) * (tau(iz, ix) - tau(iz-1,ix))

                sigma(iz,ix) = sigma(iz,ix) + (L(iz,ix) + 2*M(iz,ix)) * (dt/dx) * (u(iz,ix+1) - u(iz,ix)) + &
                                               L(iz,ix) * (dt/dz) * (v(iz+1,ix) - v(iz,ix))

                tau(iz,ix)   = tau(iz,ix)   + (L(iz,ix) + 2*M(iz,ix)) * (dt/dz) * (v(iz+1,ix) - v(iz,ix)) + &
                                               L(iz,ix) * (dt/dx) * (u(iz,ix+1) - u(iz,ix))

                xi(iz,ix)    = xi(iz,ix)    + M(iz,ix) * (dt/dz) * (u(iz+1,ix) - u(iz,ix))  + &
                                              M(iz,ix) * (dt/dx) * (v(iz,ix+1) - v(iz,ix))

            end do
        end do

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
        if (mod(it,NINT(nt/100.)) == 0) then
            print*, "########### At time sample ->",it, "/",nt
        end if

    end do

    inquire(iolength=reclsnaps) Vsnap

    open(3,file=outname_u,access="direct",recl=reclsnaps)
    write(3,rec=1) Usnap
    close(3)

    open(4,file=outname_v,access="direct",recl=reclsnaps)
    write(4,rec=1) Vsnap
    close(4)



    deallocate(vp2D,vs2D,rho2D)
    deallocate(B,M,L)
    deallocate(sigma,tau,xi)
    deallocate(u,v)
    deallocate(src)
    deallocate(gz,gx)
    deallocate(gzx)
    deallocate(Usnap,Vsnap)
end program
