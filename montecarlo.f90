program montecarlo
    use omp_lib
    use tools
    use simulation
    implicit none
    
    character(128) :: disp, relax
    logical :: one, mt, vol, gf
    integer :: nemit, ncell, ntime, maxscat, maxcoll
    real(8) :: length, side, wall, a, b, c, d, tend, T, Thot, Tcold
    
    character(128) :: whichsim
    
    print ('(A)'), timestamp()
    call starttimer()
    
    call parsecmdargs()
    
    call get_command_argument(1, whichsim)
    select case (whichsim)
        case ('isot')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1/) == 1)
            call initisot(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
        case ('bulk')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1/) == 1)
            call initbulk(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
        case ('film')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1/) == 1)
            call initfilm(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
        case ('hollow')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1/) == 1)
            call inithollow(disp, relax, one, mt, vol, gf, nemit, ncell, ntime, length, side, wall, tend, T, Thot, Tcold)
        case ('unit')
            vol = .false.
            gf  = .false.
            ncell = 1
            length = c/2d0*sqrt(2d0)
            call printargs((/1,1,1,1,0,0,1,0,1,1,1,1,0,0,0,1,1,1,1,1,1,1/) == 1)
            call initunit(disp, relax, one, mt, nemit, ntime, a, b, c, d, tend, T, Thot, Tcold)
        case default
            print *, 'Invalid sim type: "', trim(whichsim), '"'
            call exit
    end select
    
    if (one) then
        call simulateone(maxscat, maxcoll)
    else
        call simulate(maxscat)
        
        if (ncell > 1) call writetemp(ncell, ntime)
        
        if (gf) then
            call writeflux(ncell)
        else
            print ('(A,ES16.8,A)'), 'j = ', getflux(), ' W/m^2'
            print ('(A,ES16.8,A)'), 'k = ', getcond(Tcold - Thot, length), ' W/m-K'
        end if
    end if
    
contains

subroutine parsecmdargs()
    integer :: i, eqsign, inplen
    character(128) :: input, arg, val
    character(130) :: charval
    
    disp  = './input/Si_disp.txt'
    relax = './input/Si_relax.txt'
    
    one = .false.
    mt  = .true.
    vol = .false.
    gf  = .false.

    nemit = 10000
    ncell = 1
    ntime = 0               !zero for steady
    maxscat = 1             !zero for time limit only
    maxcoll = 10*maxscat    !for single phonon simulation
    
    tend = 100d-9
    
    length = 1d-6
    side = 100d-9
    wall = 10d-9
    
    a = 0.8d-6
    b = 0.2d-6
    c = 8.0d-6
    d = 40.d-9
    
    T = 300d0
    Thot = 3d0
    Tcold = -3d0
    
    do i = 2, command_argument_count()
        call get_command_argument(i, input)
        eqsign = index(input, '=')
        inplen = len_trim(input)
        
        if (eqsign > 1 .and. eqsign < inplen) then
            arg = input(1:eqsign-1)
            val = input(eqsign+1:inplen)
        
            select case (arg)
                case ('disp')
                    charval = '"'// val //'"'
                    read(charval,*) disp
                case ('relax')
                    charval = '"'// val //'"'
                    read(charval,*) relax
                case ('one')
                    read(val,*) one
                case ('mt')
                    read(val,*) mt
                case ('vol')
                    read(val,*) vol
                case ('gf')
                    read(val,*) gf
                case ('nemit')
                    read(val,*) nemit
                case ('ncell')
                    read(val,*) ncell
                case ('ntime')
                    read(val,*) ntime
                case ('maxscat')
                    read(val,*) maxscat
                case ('maxcoll')
                    read(val,*) maxcoll
                case ('tend')
                    read(val,*) tend
                case ('length')
                    read(val,*) length
                case ('side')
                    read(val,*) side
                case ('wall')
                    read(val,*) wall
                case ('a')
                    read(val,*) a
                case ('b')
                    read(val,*) b
                case ('c')
                    read(val,*) c
                case ('d')
                    read(val,*) d
                case ('T')
                    read(val,*) T
                case ('Thot')
                    read(val,*) Thot
                case ('Tcold')
                    read(val,*) Tcold
                case default
                    print *, 'Invalid input: "', trim(input), '"'
            end select
        else
            print *, 'Invalid input: "', trim(input), '"'
        end if
    end do
end subroutine parsecmdargs

subroutine printargs(isused)
    logical, intent(in) :: isused(22)
    
    if (isused(1))  write(*,'(A12,A)') 'disp = ', trim(disp)
    if (isused(2))  write(*,'(A12,A)') 'relax = ', trim(relax)
    if (isused(3))  write(*,'(A12,L2)') 'one = ', one
    if (isused(4))  write(*,'(A12,L2)') 'mt = ', mt
    if (isused(5))  write(*,'(A12,L2)') 'vol = ', vol
    if (isused(6))  write(*,'(A12,L2)') 'gf = ', gf
    if (isused(7))  write(*,'(A12,I10)') 'nemit = ', nemit
    if (isused(8))  write(*,'(A12,I10)') 'ncell = ', ncell
    if (isused(9))  write(*,'(A12,I10)') 'ntime = ', ntime
    if (isused(10)) write(*,'(A12,I10)') 'maxscat = ', maxscat
    if (isused(11)) write(*,'(A12,I10)') 'maxcoll = ', maxcoll
    if (isused(12)) write(*,'(A12,ES10.3)') 'tend = ', tend
    if (isused(13)) write(*,'(A12,ES10.3)') 'length = ', length
    if (isused(14)) write(*,'(A12,ES10.3)') 'side = ', side
    if (isused(15)) write(*,'(A12,ES10.3)') 'wall = ', wall
    if (isused(16)) write(*,'(A12,ES10.3)') 'a = ', a
    if (isused(17)) write(*,'(A12,ES10.3)') 'b = ', b
    if (isused(18)) write(*,'(A12,ES10.3)') 'c = ', c
    if (isused(19)) write(*,'(A12,ES10.3)') 'd = ', d
    if (isused(20)) write(*,'(A12,F10.3)') 'T = ', T
    if (isused(21)) write(*,'(A12,F10.3)') 'Thot = ', Thot
    if (isused(22)) write(*,'(A12,F10.3)') 'Tcold = ', Tcold
    
end subroutine printargs

end program montecarlo