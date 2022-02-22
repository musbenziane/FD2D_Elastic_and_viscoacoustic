! Created by mus on 29/07/2021.

program create2Dmodel_files
    implicit none
    integer :: modtype, ns, nend, nz, nx
    real (kind=8) :: c, qf, rho, rat, ratqf
    real (kind=8), dimension(:,:), allocatable :: c2D, qf2D, rho2D

    write(*,*) "Simple model -> 0 : ,c = 2000 m/s, qf = 80, rho = 1800 kg/m³"
    write(*,*) "Homogeneous -> 1 | Heterogenous  [1 diffrent layer ] -> 2"
    read (*,*) modtype
    write(*,*) "Give the number of elements: nz  nx"
    read (*,*) nz, nx

    allocate(c2D(nz,nx))
    allocate(qf2D(nz,nx))
    allocate(rho2D(nz,nx))


    select case (modtype)

        case (0)
        c2D(:,:)   = 2000
        qf2D(:,:)  = 2000
        rho2D(:,:) = 1800
        
        open(8,file="testing_c",access='direct',recl=nz*nx*8)
        write(8,rec=1) c2D
        close(8)

        open(9,file="testing_qf",access='direct',recl=nz*nx*8)
        write(9,rec=1) qf2D
        close(9)

        open(10,file="testing_rho",access='direct',recl=nz*nx*8)
        write(10,rec=1) rho2D
        close(10)

        deallocate(c2D)
        deallocate(qf2D)
        deallocate(rho2D)

        case (1)
        write(*,*) "Give velocity, Q  & density values"
        read (*,*) c, qf, rho
        c2D(:,:)   = c
        qf2D(:,:)  = qf
        rho2D(:,:) = rho

        open(11,file="testing_c",access='direct',recl=nz*nx*8)
        write(11,rec=1) c2D
        close(11)

        open(12,file="testing_qf",access='direct',recl=nz*nx*8)
        write(12,rec=1) qf2D
        close(12)

        open(13,file="testing_rho",access='direct',recl=nz*nx*8)
        write(13,rec=1) rho2D
        close(13)

        deallocate(c2D)
        deallocate(qf2D)
        deallocate(rho2D)

        case (2)
        write(*,*) "Give velocity, Q & density values of the background"
        read (*,*) c, qf, rho
        write(*,*) "Give ratio (for c & rho) of the low/high velocity zone not in %"
        read (*,*) rat
        write(*,*) "Give ratio (for Q) of the low/high velocity zone not in %"
        read (*,*) ratqf
        write(*,*) "give position in element number of the low/high velocity zone [nstart, nend]"
        read (*,*) ns, nend

        c2D(:,:)   = c
        qf2D(:,:)  = qf
        rho2D(:,:) = rho

        c2D(ns:nend,:)  = c * rat
        qf2D(ns:nend,:) = qf * ratqf
        rho2D(ns:nend,:)= rho * rat

        open(14,file="testing_c",access='direct',recl=nz*nx*8)
        write(14,rec=1) c2D
        close(14)

        open(15,file="testing_qf",access='direct',recl=nz*nx*8)
        write(15,rec=1) qf2D
        close(15)

        open(16,file="testing_rho",access='direct',recl=nz*nx*8)
        write(16,rec=1) rho2D
        close(16)

        deallocate(c2D)
        deallocate(qf2D)
        deallocate(rho2D)


    case default
        write (*,*) "There has been an issue, start over"

    end select

    write(*,*) "Model files names are: testing_c, testing_qf & testing_rho"

end program create2Dmodel_files