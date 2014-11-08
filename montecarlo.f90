program montecarlo
    use omp_lib
    use tools
    use simulation
    implicit none
    
    character(128) :: disp, relax, output
    logical :: one, mt, vol
    integer :: maxthrd, nemit, ngrid, ntime, maxscat, maxcoll
    real(8) :: tend, length, side, wall, a, b, c, d, T, Thot, Tcold
    
    character(128) :: stamp, whichsim
    
    stamp = timestamp()
    print ('(A,/,A)'), '----------------------------------------', trim(stamp)
    
    call parsecmdargs()
    
    call preinit(disp, relax, one, mt, maxthrd, ntime, tend, T)
    
    call get_command_argument(1, whichsim)
    print ('(/,A)'), trim(whichsim)
    select case (whichsim)
        case ('isot')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1/)==1)
            call initisot(vol, ngrid, length, side, Thot, Tcold)
        case ('bulk')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1/)==1)
            call initbulk(vol, ngrid, length, side, Thot, Tcold)
        case ('film')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,1,1,1/)==1)
            call initfilm(vol, ngrid, length, side, Thot, Tcold)
        case ('hollow')
            call printargs((/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,1,1/)==1)
            call inithollow(vol, ngrid, length, side, wall, Thot, Tcold)
        case ('unit')
            vol = .false.
            ngrid = 3
            call printargs((/1,1,1,1,1,0,1,1,0,1,1,1,1,0,0,0,1,1,1,1,1,1,1/)==1)
            call initunit(a, b, c, d, Thot, Tcold)
        case default
            print *, 'Invalid sim type: "', trim(whichsim), '"'
            call exit
    end select
    
    call postinit(nemit)
    
    if (one) then
        call simulateone(maxscat, maxcoll)
    else
        call simulate(maxscat)
        
!         if (ngrid > 1) then
!             call writetemp(ngrid, ntime)
!             call writeflux(ngrid)
!         else
!             print ('(A,ES16.8,A)'), 'j = ', getflux(), ' W/m^2'
!             print ('(A,ES16.8,A)'), 'k = ', getcond(Tcold - Thot, length), ' W/m-K'
!         end if
    end if
    
    call postsimulate()
    
    if (output == '') then
        call writeresults(one, maxcoll, ngrid, ntime)
    else
        call writeresults(one, maxcoll, ngrid, ntime, 2, output, stamp)
    end if
    
contains

subroutine parsecmdargs()
    integer :: i, eqsign, inplen
    integer, parameter :: collfactor = 10
    character(128) :: input, arg, val
    character(130) :: charval
    
    disp  = './input/Si_disp.txt'
    relax = './input/Si_relax.txt'
    output = ''
    
    one = .false.
    mt  = .true.
    vol = .false.
    
    maxthrd = 0             !zero for no limit
    nemit = 10000
    ngrid = 0
    ntime = 0               !zero for steady
    maxscat = 1             !zero for time limit only
    maxcoll = collfactor    !for single phonon simulation
    
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
                case ('output')
                    charval = '"'// val //'"'
                    read(charval,*) output
                case ('one')
                    read(val,*) one
                case ('mt')
                    read(val,*) mt
                case ('vol')
                    read(val,*) vol
                case ('maxthrd')
                    read(val,*) maxthrd
                case ('nemit')
                    read(val,*) nemit
                case ('ngrid')
                    read(val,*) ngrid
                case ('ntime')
                    read(val,*) ntime
                case ('maxscat')
                    read(val,*) maxscat
                    if (maxcoll == collfactor) maxcoll = collfactor*maxscat
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
    logical, intent(in) :: isused(23)
    
    if (isused(1))  write(*,'(A12,A)')      '   disp = ', trim(disp)
    if (isused(2))  write(*,'(A12,A)')      '  relax = ', trim(relax)
    if (isused(3))  write(*,'(A12,A)')      ' output = ', trim(output)
    if (isused(4))  write(*,'(A12,L2)')     '    one = ', one
    if (isused(5))  write(*,'(A12,L2)')     '     mt = ', mt
    if (isused(6))  write(*,'(A12,L2)')     '    vol = ', vol
    if (isused(7))  write(*,'(A12,I11)')    'maxthrd = ', maxthrd
    if (isused(8))  write(*,'(A12,I11)')    '  nemit = ', nemit
    if (isused(9))  write(*,'(A12,I11)')    '  ngrid = ', ngrid
    if (isused(10)) write(*,'(A12,I11)')    '  ntime = ', ntime
    if (isused(11)) write(*,'(A12,I11)')    'maxscat = ', maxscat
    if (isused(12)) write(*,'(A12,I11)')    'maxcoll = ', maxcoll
    if (isused(13)) write(*,'(A12,ES10.3)') '   tend = ', tend
    if (isused(14)) write(*,'(A12,ES10.3)') ' length = ', length
    if (isused(15)) write(*,'(A12,ES10.3)') '   side = ', side
    if (isused(16)) write(*,'(A12,ES10.3)') '   wall = ', wall
    if (isused(17)) write(*,'(A12,ES10.3)') '      a = ', a
    if (isused(18)) write(*,'(A12,ES10.3)') '      b = ', b
    if (isused(19)) write(*,'(A12,ES10.3)') '      c = ', c
    if (isused(20)) write(*,'(A12,ES10.3)') '      d = ', d
    if (isused(21)) write(*,'(A12,F10.3)')  '      T = ', T
    if (isused(22)) write(*,'(A12,F10.3)')  '   Thot = ', Thot
    if (isused(23)) write(*,'(A12,F10.3)')  '  Tcold = ', Tcold
    write(*,*)
    
end subroutine printargs

end program montecarlo