! Ricker wavelet
subroutine ricker(nt,f0,dt,source)
    implicit none
    integer                      :: nt, it
    real (kind=8)                :: f0, pt, a_ricker, dt, t0,t
    real (kind=8),dimension(nt)  :: source, temp

    source = 0
    temp   = 0
    pt = 1 / f0
    t0=pt/dt
    a_ricker=4./pt

    do it=1,nt+1
        t=(it-t0)*dt
        temp(it) = -2.*a_ricker*t*exp(-(a_ricker*t)**2)
    enddo

    source = temp
    do it=1,nt
        source(it) = temp(it+1) - temp(it)
        source(it) = -source(it)
    end do

END subroutine ricker
