! Created by mus on 30/07/2021.

subroutine readmodelfiles(nz, nx, prefix, vp2D,vs2D, rho2D)
    implicit none
    integer, intent(in)                            :: nz,nx
    real (kind=8), dimension(nz,nx), intent(out)   :: vp2D,vs2D, rho2D
    character(len=40)                              :: prefix, modname_vp, modname_vs,modname_rho

    modname_vp  = TRIM(prefix)//'_vp'
    modname_vs  = TRIM(prefix)//'_vs'
    modname_rho = TRIM(prefix)//'_rho'

    open(17,file=modname_vp,access='direct',recl=nz*nx*8)
    read(17,rec=1) vp2D
    close(17)

    open(18,file=modname_vs,access='direct',recl=nz*nx*8)
    read(18,rec=1) vs2D
    close(18)

    open(19,file=modname_rho,access='direct',recl=nz*nx*8)
    read(19,rec=1) rho2D
    close(19)

end subroutine readmodelfiles