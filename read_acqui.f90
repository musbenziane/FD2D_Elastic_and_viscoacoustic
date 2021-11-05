module read_acqui
    implicit none
    contains
    subroutine countpoints(filename,c)
        implicit none
        character(len=40), intent(in)      :: filename
        real(kind=8)                       :: tmpz, tmpx
        integer, intent(out)               :: c

        c = 0
        open(5,file=filename)
        do while(.true.)
            read(5,*,end=999) tmpz,tmpx
            c = c + 1
        end do
        999 continue
        close(5)
    end subroutine countpoints

    subroutine read_acqui_geo(filename,n,zz,xx)
        implicit none
        character(len=40), intent(in)            :: filename
        integer, intent(in)                      :: n
        real(kind=8), dimension(n), intent(out)  :: xx, zz
        integer                                  :: c

        open(65,file=filename)
        c = 1
        do while(.true.)
            read(65,*,end=994) zz(c),xx(c)
            c = c + 1
            if (c==n) then
                exit
            end if
        end do
        994 continue
        close(65)

    end subroutine read_acqui_geo

end module read_acqui
