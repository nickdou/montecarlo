program montecarlo
	use tools
	use simulation
	implicit none
	
	logical :: one, vol, gf
	
	integer :: nemit, ncell, ntime, maxscat, maxcoll
	integer :: i, j
	integer, parameter :: nside = 1, nwall = 1
	
	real(8) :: a, b, c, d, length, side, wall, tend, T, Thot, Tcold
	real(8) :: k, cond(nside, nwall)
	
	character(len=*), parameter :: disp  = 'input/Si_disp.txt'
	character(len=*), parameter :: relax = 'input/Si_relax.txt'
	
	one = .false.
	vol = .true.
	gf  = .true.
	
	nemit = 10000000
	ncell = 100
	ntime = 0              !zero for steady
	length = 1e-6
	side = 100e-9
	wall = 10d-9
	tend = 100d-9
	T = 300d0
	Thot = 3d0
	Tcold = -3d0
	
	maxscat = 1             !zero for time limit only
	maxcoll = 2*maxscat     !for single phonon simulation
	
	! isot: maxscat high
	! bulk: maxscat one
	! film: maxscat one, gf true
	
!	do j = 1,nwall
!		do i = 1,nside
!			call run(disp, relax, one, vol, gf, nemit, ncell, ntime, length, side(i), wall(j), tend, T, Thot, Tcold, maxscat, maxcoll, k)
!			cond(i,j) = k
!		end do
!	end do
!	print ('(A12)'), 'k = '
!	call printarray(transpose(cond), '(ES16.8)')
	
! 	call run(disp, relax, one, vol, gf, nemit, ncell, ntime, length, side, wall, tend, T, Thot, Tcold, maxscat, maxcoll)
	
	one = .false.
	nemit = 10000000
	ntime = 0
	a = 0.8d-6
	b = 0.2d-6
	c = 8.0d-6
	d = 40.d-9
	tend = 100d-9
	T = 300d0
	Thot = 3d0
	Tcold = -3d0
	
	maxscat = 100
	maxcoll = 2*maxscat
	
	call rununit(disp, relax, one, nemit, ntime, a, b, c, d, tend, T, Thot, Tcold, maxscat, maxcoll, k)
	
contains

subroutine run(disp, relax, one, vol, gf, nemit, ncell, ntime, length, side, wall, tend, T, Thot, Tcold, maxscat, maxcoll, k)
	character(len=*), intent(in) :: disp, relax
	logical, intent(in) :: one, vol, gf
	integer, intent(in) :: nemit, ncell, ntime, maxscat, maxcoll
	real(8), intent(in) :: length, side, wall, tend, T, Thot, Tcold
	real(8), intent(out), optional :: k
	real(8) :: cond
	
	print ('(A12,A)'), 'disp = ', disp
	print ('(A12,A)'), 'relax = ', relax
	print ('(A12,L)'), 'one = ', one
	print ('(A12,L)'), 'vol = ', vol
	print ('(A12,L)'), 'gf = ', gf
	print ('(A12,I10)'), 'nemit = ', nemit
	print ('(A12,I10)'), 'ncell = ', ncell
	print ('(A12,I10)'), 'ntime = ', ntime
	print ('(A12,I10)'), 'maxscat = ', maxscat
	print ('(A12,I10)'), 'maxcoll = ', maxcoll
	print ('(A12,ES10.3)'), 'tend = ', tend
	print ('(A12,ES10.3)'), 'length = ', length
	print ('(A12,ES10.3)'), 'side = ', side
	print ('(A12,ES10.3)'), 'wall = ', wall
	print ('(A12,F10.3)'), 'T = ', T
	print ('(A12,F10.3)'), 'Thot = ', Thot
	print ('(A12,F10.3)'), 'Tcold = ', Tcold
	
! 	call initisot1d(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
! 	call initbulk1d(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
! 	call initfilm(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
	call inithollow(disp, relax, vol, gf, nemit, ncell, ntime, length, side, wall, tend, T, Thot, Tcold)
	
	if (one) then
		call simulateone(maxscat, maxcoll)
		
	else
		call simulate(maxscat)
		
		call writetemp(ncell, ntime)
		if (gf) then
			call writeflux(ncell)
		else
			cond = getcond(Tcold - Thot, length)
			if (present(k)) then
				k = cond
			end if
		end if
	end if
end subroutine run

subroutine rununit(disp, relax, one, nemit, ntime, a, b, c, d, tend, T, Thot, Tcold, maxscat, maxcoll, k)
	character(len=*), intent(in) :: disp, relax
	logical, intent(in) :: one
	integer, intent(in) :: nemit, ntime, maxscat, maxcoll
	real(8), intent(in) :: a, b, c, d, tend, T, Thot, Tcold
	real(8), intent(out), optional :: k
	real(8) :: cond
	
	print ('(A12,A)'), 'disp = ', disp
	print ('(A12,A)'), 'relax = ', relax
	print ('(A12,L)'), 'one = ', one
	print ('(A12,I10)'), 'nemit = ', nemit
	print ('(A12,I10)'), 'ntime = ', ntime
	print ('(A12,I10)'), 'maxscat = ', maxscat
	print ('(A12,I10)'), 'maxcoll = ', maxcoll
	print ('(A12,ES10.3)'), 'tend = ', tend 
	print ('(A12,ES10.3)'), 'a = ', a
	print ('(A12,ES10.3)'), 'b = ', b
	print ('(A12,ES10.3)'), 'c = ', c
	print ('(A12,ES10.3)'), 'd = ', d
	print ('(A12,F10.3)'), 'T = ', T
	print ('(A12,F10.3)'), 'Thot = ', Thot
	print ('(A12,F10.3)'), 'Tcold = ', Tcold
	
	call initunit(disp, relax, nemit, ntime, a, b, c, d, tend, T, Thot, Tcold)
	
	if (one) then
		call simulateone(maxscat, maxcoll)
	
	else
		call simulate(maxscat)
		
		cond = getcond(Tcold - Thot, length)
		if (present(k)) then
			k = cond
		end if
	end if
	
end subroutine rununit

end program montecarlo