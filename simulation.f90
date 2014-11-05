module simulation
    use constants
    use tools
    use math
    use material
    use distributions
    use domain
    use particle
    implicit none
    
    public
    
contains

subroutine initisot(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
    character(len=*), intent(in) :: disp, relax
    logical, intent(in) :: one, mt, vol, gf
    integer, intent(in) :: nemit, ncell, ntime
    real(8), intent(in) :: length, side, tend, T, Thot, Tcold
    integer, parameter :: nbdry = 6
!     type(axis) :: grid, vgen
    type(boundary) :: bdry_arr(nbdry)
    real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
    integer :: emit_arr(nbdry)
    
    if (one) then
        call initomp(1)
    else
        call initomp()
    end if
    call initrand(mt)

    call initmat(disp, relax, T)
    call initdist()

!     grid = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, ncell)
!     call setgrid(tend, ntime, grid)
    call settime(tend, ntime)
    call setgrid((/0d0, 0d0, 1d0/), 0d0, length, ncell, product(corner))
    call initrecord((/0d0, 0d0, 1d0/))

    origin = (/0d0, 0d0, 0d0/)
    xvec = (/side, 0d0, 0d0/)
    yvec = (/0d0, side, 0d0/)
    zvec = (/0d0, 0d0, length/)
    corner = xvec + yvec + zvec
    bdry_arr(1) = makebdry(ISOT_BC, origin, xvec, yvec, Thot)
    bdry_arr(2) = makebdry(SPEC_BC, origin, yvec, zvec)
    bdry_arr(3) = makebdry(SPEC_BC, origin, zvec, xvec)
    bdry_arr(4) = makebdry(ISOT_BC, corner, -yvec, -xvec, Tcold)
    bdry_arr(5) = makebdry(SPEC_BC, corner, -zvec, -yvec)
    bdry_arr(6) = makebdry(SPEC_BC, corner, -xvec, -zvec)
!     call setbdry(bdry_arr, product(corner))
    call setbdry(bdry_arr)
    call calculateemit(nemit)
    
!     vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
!     call setvolumetric(vol, vgen)
    if (vol) then
        call setvolumetric((/0d0, 0d0, 1d0/), (/0d0, length/))
    end if
    
    emit_arr = getemit()
    print ('(/,A)'), 'Isothermal'
    print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
    print ('(A12,ES15.8)'), 'ktheory = ', ktheory
    print ('(A12,16I10)'),  'emit_arr = ', emit_arr
end subroutine initisot

subroutine initbulk(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
    character(len=*), intent(in) :: disp, relax
    logical, intent(in) :: one, mt, vol, gf
    integer, intent(in) :: nemit, ncell, ntime
    real(8), intent(in) :: length, side, tend, T, Thot, Tcold
    real(8) :: zvec(3)

    call initisot(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)

    zvec = (/0d0, 0d0, length/)
    call setbdrypair(1, 4, zvec)
    call calculateemit(nemit)

    print ('(/,A)'), 'Bulk'
end subroutine initbulk

subroutine initfilm(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
    character(len=*), intent(in) :: disp, relax
    logical, intent(in) :: one, mt, vol, gf
    integer, intent(in) :: nemit, ncell, ntime
    real(8), intent(in) :: length, side, tend, T, Thot, Tcold
    integer, parameter :: nbdry = 6
!     type(axis) :: grid, vgen
    type(boundary) :: bdry_arr(nbdry)
    real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
    integer :: emit_arr(nbdry)

    if (one) then
        call initomp(1)
    else
        call initomp()
    end if
    call initrand(mt)

    call initmat(disp, relax, T)
    call initdist()

!     grid = makeaxis((/0d0, 1d0, 0d0/), 0d0, side, ncell)
!     call setgrid(tend, ntime, grid)
    call settime(tend, ntime)
    call setgrid((/0d0, 1d0, 0d0/), 0d0, side, ncell, product(corner))
    call initrecord((/0d0, 0d0, 1d0/))

    origin = (/0d0, 0d0, 0d0/)
    xvec = (/length, 0d0, 0d0/)
    yvec = (/0d0, side, 0d0/)
    zvec = (/0d0, 0d0, length/)
    corner = xvec + yvec + zvec

    bdry_arr(1) = makebdry(PERI_BC, origin, xvec, yvec, Thot)
    bdry_arr(2) = makebdry(SPEC_BC, origin, yvec, zvec)
    bdry_arr(3) = makebdry(DIFF_BC, origin, zvec, xvec)
    bdry_arr(4) = makebdry(PERI_BC, corner, -yvec, -xvec, Tcold)
    bdry_arr(5) = makebdry(SPEC_BC, corner, -zvec, -yvec)
    bdry_arr(6) = makebdry(DIFF_BC, corner, -xvec, -zvec)
!     call setbdry(bdry_arr, product(corner))
    call setbdry(bdry_arr)
    call setbdrypair(1, 4, zvec)
    call calculateemit(nemit)

!     vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
!     call setvolumetric(vol, vgen)
    if (vol) then
        call setvolumetric((/0d0, 0d0, 1d0/), (/0d0, length/))
    end if

    emit_arr = getemit()
    print ('(/,A)'), 'Film'
    print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
    print ('(A12,ES15.8)'), 'ktheory = ', ktheory
    print ('(A12,16I10)'),   'emit_arr = ', emit_arr
end subroutine initfilm

subroutine inithollow(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, wall, tend, T, Thot, Tcold)
    character(len=*), intent(in) :: disp, relax
    logical, intent(in) :: one, mt, vol, gf
    integer, intent(in) :: nemit, ncell, ntime
    real(8), intent(in) :: length, side, wall, tend, T, Thot, Tcold
    integer, parameter :: nbdry = 16
!     type(axis) :: grid, vgen
    type(boundary) :: bdry_arr(nbdry)
    real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3), xwall(3), ywall(3)
    integer :: emit_arr(nbdry)

    if (one) then
        call initomp(1)
    else
        call initomp()
    end if
    call initrand(mt)

    call initmat(disp, relax, T)
    call initdist()

!     grid = makeaxis((/0d0, 1d0, 0d0/), 0d0, side, ncell)
!     call setgrid(tend, ntime, grid)
    call settime(tend, ntime)
    call setgrid((/0d0, 1d0, 0d0/), 0d0, side, ncell, &
        product(corner) - product(corner - 2*xwall - 2*ywall))
    call initrecord((/0d0, 0d0, 1d0/))

    origin = (/0d0, 0d0, 0d0/)
    xvec = (/side, 0d0, 0d0/)
    yvec = (/0d0, side, 0d0/)
    zvec = (/0d0, 0d0, length/)
    corner = xvec + yvec + zvec
    xwall = (/wall, 0d0, 0d0/)
    ywall = (/0d0, wall, 0d0/)

    bdry_arr(1) = makebdry(PERI_BC, origin, xvec-xwall, ywall, Thot)
    bdry_arr(2) = makebdry(PERI_BC, origin + xvec-xwall, xwall, yvec-ywall, Thot)
    bdry_arr(3) = makebdry(PERI_BC, origin + xwall + yvec-ywall, xvec-xwall, ywall, Thot)
    bdry_arr(4) = makebdry(PERI_BC, origin + ywall, xwall, yvec-ywall, Thot)

    bdry_arr(5) = makebdry(PERI_BC, corner, -ywall, -xvec+xwall, Tcold)
    bdry_arr(6) = makebdry(PERI_BC, corner - xvec+xwall, -yvec+ywall, -xwall, Tcold)
    bdry_arr(7) = makebdry(PERI_BC, corner - xwall - yvec+ywall, -ywall, -xvec+xwall, Tcold)
    bdry_arr(8) = makebdry(PERI_BC, corner - ywall, -yvec+ywall, -xwall, Tcold)

    bdry_arr(9) = makebdry(DIFF_BC, origin, yvec, zvec)
    bdry_arr(10) = makebdry(DIFF_BC, origin + xwall+ywall, zvec, yvec-2*ywall)
    bdry_arr(11) = makebdry(DIFF_BC, origin, zvec, xvec)
    bdry_arr(12) = makebdry(DIFF_BC, origin + xwall+ywall, xvec-2*xwall, zvec)

    bdry_arr(13) = makebdry(DIFF_BC, corner, -zvec, -yvec)
    bdry_arr(14) = makebdry(DIFF_BC, corner - xwall-ywall, -yvec+2*ywall, -zvec)
    bdry_arr(15) = makebdry(DIFF_BC, corner, -xvec, -zvec)
    bdry_arr(16) = makebdry(DIFF_BC, corner - xwall-ywall, -zvec, -xvec+2*xwall)

!     call setbdry(bdry_arr, product(corner) - product(corner - 2*xwall - 2*ywall))
    call setbdry(bdry_arr)
    call setbdrypair(1, 7, zvec)
    call setbdrypair(2, 8, zvec)
    call setbdrypair(3, 5, zvec)
    call setbdrypair(4, 6, zvec)
    call calculateemit(nemit)

!     vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
!     call setvolumetric(vol, vgen)
    if (vol) then
        call setvolumetric((/0d0, 0d0, 1d0/), (/0d0, length/))
    end if

    emit_arr = getemit()
    print ('(/,A)'), 'Hollow'
    print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
    print ('(A12,ES15.8)'), 'ktheory = ', ktheory
    print ('(A12,16I10)'),   'emit_arr = ', emit_arr
end subroutine inithollow

subroutine initunit(disp, relax, one, mt, nemit, ntime, a, b, c, d, tend, T, Thot, Tcold)
    character(len=*), intent(in) :: disp, relax
    logical, intent(in) :: one, mt
    integer, intent(in) :: nemit, ntime
    real(8), intent(in) :: a, b, c, d, tend, T, Thot, Tcold
    integer, parameter :: nvert = 50
    integer, parameter :: nbdry = 34
!     type(axis) :: grid, vgen
    type(boundary) :: bdry_arr(nbdry)
    real(8) :: l, vertshi(3, nvert/2), verts(3, nvert), temps(nbdry), mv(3)
    integer :: i, bdrys(4, nbdry), emit_arr(nbdry)
    logical :: tris(nbdry)

    if (one) then
        call initomp(1)
    else
        call initomp()
    end if
    call initrand(mt)

    call initmat(disp, relax, T)
    call initdist()

    l = c/4d0*sqrt(2d0)

!     grid = makeaxis((/0d0, 0d0, 1d0/), -l, l, 1) ! ncell = 1
!     call setgrid(tend, ntime, grid)
    call settime(tend, ntime)
    call setgrid((/0d0, 0d0, 1d0/), (/-l, -(l+a)/2, (l+a)/2, l/), (/1d0,1d0,1d0/))
    call initrecord((/0d0, 0d0, 1d0/))
    
    vertshi = reshape((/ &
          0d0, 0d0,   l, & 
          0d0,   b,   l, &
          2*a,   b,   l, &
          2*a, 0d0,   l, &
        2*a+d, 0d0,   l, &
        2*a+d,   b,   l, &
        2*a+d, b+d,   l, &
           -d, b+d,   l, &
           -d,   b,   l, &
           -d, 0d0,   l, &
          2*a,   b,   a, &
            a,   b,   a, &
          0d0,   b,   a, &
          0d0, 0d0,   a, &
           -d, 0d0,   a, &
           -d,   b,   a, &
           -d, b+d,   a, &
           -d, b+d, a+d, &
          a+d, b+d, a+d, &
        2*a+d, b+d, a+d, &
            a,   l,   a, &
           -d,   l,   a, &
           -d,   l, a+d, &
            a,   l, a+d, &
          a+d,   l, a+d &
        /), (/3, nvert/2/))
    
    do i = 1,nvert/2
        verts(:,i) = vertshi(:,i)
        verts(1:2,nvert-i+1) = vertshi(1:2,i)
        verts(3,nvert-i+1) = -vertshi(3,i)
    end do
    
    bdrys = reshape((/ &
         1,  2,  9, PERI_BC, & !01
         3,  4,  5, PERI_BC, & !02
         6,  7,  8, PERI_BC, & !03
         3,  2, 13, DIFF_BC, & !04
         2,  1, 14, DIFF_BC, & !05
         1, 10, 15, SPEC_BC, & !06
        10,  8, 17, DIFF_BC, & !07
         8,  7, 20, DIFF_BC, & !08
        13, 14, 15, SPEC_BC, & !09
        12, 16, 22, DIFF_BC, & !10
        17, 18, 23, SPEC_BC, & !11
        21, 22, 23, SPEC_BC, & !12
        18, 19, 25, DIFF_BC, & !13
        11, 12, 39, DIFF_BC, & !14
        12, 21, 30, DIFF_BC, & !15
        24, 25, 26, SPEC_BC, & !16
        25, 19, 32, DIFF_BC, & !17
        19, 20, 31, DIFF_BC, & !18
         7,  5, 46, DIFF_BC, & !19
         5,  4, 47, SPEC_BC, & !20
         4,  3, 48, DIFF_BC, & !21
        32, 33, 28, DIFF_BC, & !22
        27, 28, 29, SPEC_BC, & !23
        33, 34, 29, SPEC_BC, & !24
        35, 39, 30, DIFF_BC, & !25
        35, 36, 37, SPEC_BC, & !26
        33, 31, 44, DIFF_BC, & !27
        36, 34, 43, DIFF_BC, & !28
        37, 36, 41, SPEC_BC, & !29
        38, 37, 50, DIFF_BC, & !30
        40, 38, 49, DIFF_BC, & !31
        42, 43, 44, PERI_BC, & !32
        45, 46, 47, PERI_BC, & !33
        41, 42, 49, PERI_BC &  !34
        /), (/4, nbdry/))
    
    temps = 0d0
    temps(1:3) = Tcold
    temps(nbdry-2:nbdry) = Thot
    tris = .false.
    
    bdry_arr = makebdry_arr(verts, bdrys(4,:), bdrys(1:3,:), temps, tris)
    
!     call setbdry(bdry_arr, 2*d*(a+b+d)*(2*l-a) + 4*d*(a+d)*(l-b-d))
    call setbdry(bdry_arr)
    
    mv = (/0d0, 0d0, -2*l/)
    call setbdrypair(1, nbdry, mv)
    call setbdrypair(2, nbdry-1, mv)
    call setbdrypair(3, nbdry-2, mv)
    
    call calculateemit(nemit)

    ! No volumetric generation
!     vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, 1d0, 0)
!     call setvolumetric(.false., vgen)
    
    emit_arr = getemit()
    print ('(/,A)'), 'Unit'
    print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
    print ('(A12,ES15.8)'), 'ktheory = ', ktheory
    print ('(A12,34I10)'),   'emit_arr = ', emit_arr
end subroutine initunit

subroutine simulate(maxscat)
    integer, intent(in) :: maxscat
    integer :: i, nemit, nscat, progress
    real(8) :: t
    type(phonon) :: phn
    
    progress = 0
    nemit = getnemit()
    !$omp parallel do private(phn, t, nscat) shared(progress)
    do i = 1, nemit
        phn = emit()

        call drawemittime(t)
        nscat = 0
        do while (isalive(phn))
            call advect(phn, t)
            call scatter(phn, nscat)
            if (maxscat > 0 .and. nscat >= maxscat) then
                call remove(phn)
            end if
        end do
        !$omp critical
        progress = progress + 1
        call showprogress(progress, nemit, min(nemit, 20))
        !$omp end critical
    end do
    !$omp end parallel do
end subroutine simulate

subroutine simulateone(maxscat, maxcoll)
    integer, intent(in) :: maxscat, maxcoll
    integer :: ncoll, nscat, ntraj
    real(8) :: t, traj(3, 0:2*maxcoll)
    type(phonon) :: phn
    
    phn = emit()
    
    ntraj = 2*maxcoll
    call inittraj( ntraj, getpos(phn) )

    t = 0
    ncoll = 0
    nscat = 0
    do while (isalive(phn))
        call advect(phn, t)
        call scatter(phn, nscat)
        ncoll = ncoll + 1
        if (ncoll >= maxcoll .or. (maxscat > 0 .and. nscat >= maxscat)) then
            call remove(phn)
        end if
    end do
    
    call gettraj( ntraj, traj )
    call writematlab( transpose(traj(:, 0:ntraj)), '(ES16.8)', 2, 'traj', 'x' )
end subroutine simulateone

subroutine writetemp(ncell, ntime)
    integer, intent(in) :: ncell, ntime
    real(8) :: T(ncell, 0:ntime)
    
    if (ntime == 0) then
        T(:,0) = getsteadytemp()
        call writematlab(T(:,0), '(ES16.8)', 2, 'temp', 'T')
    else
        T = gettranstemp()
        call writematlab(T, '(ES16.8)', 2, 'temp', 'T')
    end if
end subroutine writetemp

subroutine writeflux(ncell)
    integer, intent(in) :: ncell
    real(8) :: flux(ncell)
    
    flux = getgridflux()
    call writematlab(flux, '(ES16.8)', 3, 'flux', 'j')
end subroutine writeflux

real(8) function getcond(deltaT, length)
    real(8), intent(in) :: deltaT, length
    real(8) :: flux
    
    flux = getflux()
    getcond = -flux *length/deltaT
end function getcond

end module simulation