! Created by mus on 29/07/2021.

program create2Dmodel_files
    implicit none
    integer :: modtype, ns, nend, nz, nx
    real (kind=8) :: vs, vp, rho, rat
    real (kind=8), dimension(:,:), allocatable :: vp2D, vs2D, rho2D

    write(*,*) "Simple model -> 0 : ,vp = 3500 m/s, vs = 2000 m/s, rho = 1800 kg/m³"
    write(*,*) "Homogeneous -> 1 | Heterogenous  [1 diffrent layer ] -> 2"
    read (*,*) modtype
    write(*,*) "Give the number of elements: nz  nx"
    read (*,*) nz, nx

    allocate(vp2D(nz,nx))
    allocate(vs2D(nz,nx))
    allocate(rho2D(nz,nx))


    select case (modtype)

        case (0)
        vp2D(:,:)  = 3500
        vs2D(:,:)  = 2000
        rho2D(:,:) = 1800
        open(8,file="testing_vp",access='direct',recl=nz*nx*8)
        write(8,rec=1) vp2D
        close(8)

        open(9,file="testing_vs",access='direct',recl=nz*nx*8)
        write(9,rec=1) vs2D
        close(9)

        open(10,file="testing_rho",access='direct',recl=nz*nx*8)
        write(10,rec=1) rho2D
        close(10)

        deallocate(vp2D)
        deallocate(vs2D)
        deallocate(rho2D)

        case (1)
        write(*,*) "Give velocity (P a S) value & density value"
        read (*,*) vp, vs, rho
        vp2D(:,:)  = vp
        vs2D(:,:)  = vs
        rho2D(:,:) = rho

        open(11,file="testing_vp",access='direct',recl=nz*nx*8)
        write(11,rec=1) vp2D
        close(11)

        open(12,file="testing_vs",access='direct',recl=nz*nx*8)
        write(12,rec=1) vs2D
        close(12)

        open(13,file="testing_rho",access='direct',recl=nz*nx*8)
        write(13,rec=1) rho2D
        close(13)

        deallocate(vp2D)
        deallocate(vs2D)
        deallocate(rho2D)

        case (2)
        write(*,*) "Give velocity (P & S) value & density value of the background"
        read (*,*) vp, vs, rho
        write(*,*) "Give ratio of the low/high velocity zone not in %"
        read (*,*) rat
        write(*,*) "give position in element number of the low/high velocity zone [nstart, nend]"
        read (*,*) ns, nend

        vp2D(:,:)  = vp
        vs2D(:,:)  = vs
        rho2D(:,:) = rho

        vp2D(ns:nend,:) = vp  * rat
        vs2D(ns:nend,:) = vs  * rat
        rho2D(ns:nend,:)= rho * rat

        open(14,file="testing_vp",access='direct',recl=nz*nx*8)
        write(14,rec=1) vp2D
        close(14)

        open(15,file="testing_vs",access='direct',recl=nz*nx*8)
        write(15,rec=1) vs2D
        close(15)

        open(16,file="testing_rho",access='direct',recl=nz*nx*8)
        write(16,rec=1) rho2D
        close(16)

        deallocate(vp2D)
        deallocate(vs2D)
        deallocate(rho2D)


    case default
        write (*,*) "There has been an issue, start over"

    end select

    write(*,*) "Model files names are: testing_vp, testing_vs & testing_rho"

end program create2Dmodel_files