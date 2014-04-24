module simulation
	use constants
	use tools
	use distributions
	use domain
	use particle
	implicit none
	
	public
	
contains

subroutine initisot1d(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, &
	Teq, Thot, Tcold)
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nomega, nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, Teq, Thot, Tcold
	type(axis) :: grid, vgen
	type(rectbdry) :: bdry_arr(6)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
	
	call initrand()
	
	call initomega(nomega)
	call initpropcdf(Teq)
	
	grid = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, ncell)
	call setgrid(tend, ntime, grid)
	call initrecord(gf, (/0d0, 0d0, 1d0/))
	
	origin = (/0d0, 0d0, 0d0/)
	corner = (/side, side, length/)
	xvec = (/side, 0d0, 0d0/)
	yvec = (/0d0, side, 0d0/)
	zvec = (/0d0, 0d0, length/)
	bdry_arr(1) = makebdry(ISOT_BC, origin, xvec, yvec, Thot)
	bdry_arr(2) = makebdry(SPEC_BC, origin, yvec, zvec)
	bdry_arr(3) = makebdry(SPEC_BC, origin, zvec, xvec)
	bdry_arr(4) = makebdry(ISOT_BC, corner, -yvec, -xvec, Tcold)
	bdry_arr(5) = makebdry(SPEC_BC, corner, -zvec, -yvec)
	bdry_arr(6) = makebdry(SPEC_BC, corner, -xvec, -zvec)
	call setbdry(bdry_arr, side**2*length)
	call calculateemit(nemit)
	
	vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
	call setvolumetric(vol, vgen)
	
	print ('(A12,ES15.8)'), 'Kn = ', tau(0d0, 1, 300d0)*vel(0d0, 1)/length
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	print ('(A12,ES15.8)'), 'ktheory = ', calculatek(Teq)
	print ('(A12,6I10)'),   'emit_arr = ', getemit()
	print ('(A12,ES15.8)'), 'xmax = ', hbar*omegamax/(kb*Teq)
end subroutine initisot1d

subroutine initbulk1d(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, &
	Teq, Thot, Tcold)
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nomega, nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, Teq, Thot, Tcold
	real(8) :: zvec(3)
	
	call initisot1d(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	
	zvec = (/0d0, 0d0, length/)
	call setbdrypair(1, 4, zvec)
	call calculateemit(nemit)
end subroutine initbulk1d

subroutine initfilm(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, &
	Teq, Thot, Tcold)
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nomega, nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, Teq, Thot, Tcold
	type(axis) :: grid, vgen
	type(rectbdry) :: bdry_arr(6)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
	
	call initrand()
	
	call initomega(nomega)
	call initpropcdf(Teq)
	
	grid = makeaxis((/0d0, 1d0, 0d0/), 0d0, side, ncell)
	call setgrid(tend, ntime, grid)
	call initrecord(gf, (/0d0, 0d0, 1d0/))
	
	origin = (/0d0, 0d0, 0d0/)
	corner = (/side, side, length/)
	xvec = (/side, 0d0, 0d0/)
	yvec = (/0d0, side, 0d0/)
	zvec = (/0d0, 0d0, length/)
	bdry_arr(1) = makebdry(PERI_BC, origin, xvec, yvec, Thot)
	bdry_arr(2) = makebdry(SPEC_BC, origin, yvec, zvec)
	bdry_arr(3) = makebdry(DIFF_BC, origin, zvec, xvec)
	bdry_arr(4) = makebdry(PERI_BC, corner, -yvec, -xvec, Tcold)
	bdry_arr(5) = makebdry(SPEC_BC, corner, -zvec, -yvec)
	bdry_arr(6) = makebdry(DIFF_BC, corner, -xvec, -zvec)
	call setbdry(bdry_arr, side**2*length)
	call setbdrypair(1, 4, zvec)
	call calculateemit(nemit)
	
	vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
	call setvolumetric(vol, vgen)
	
	print ('(A12,ES15.8)'), 'Kn = ', tau(0d0, 1, 300d0)*vel(0d0, 1)/length
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	print ('(A12,ES15.8)'), 'ktheory = ', calculatek(Teq)
	print ('(A12,6I10)'),   'emit_arr = ', getemit()
	print ('(A12,ES15.8)'), 'xmax = ', hbar*omegamax/(kb*Teq)
end subroutine initfilm

subroutine simulate(maxscat)
	integer, intent(in) :: maxscat
	integer :: i, nemit, nscat
	real(8) :: t
	type(phonon) :: phn
	
	nemit = sum( getemit() )
	do i = 1, nemit
!		print ('(/,A,I2)'), 'Particle ', i 
		phn = emit()
		
		call drawemittime(t)
		nscat = 0
		do while (isalive(phn))
			call advect(phn, t)
			call scatter(phn, nscat)
			if (maxscat > 0 .and. nscat >= maxscat) then
				call kill(phn)
			end if
		end do
		call showprogress(i, nemit, min(nemit, 20))
	end do
end subroutine simulate

subroutine simulateone(maxscat, maxcoll)
	integer, intent(in) :: maxscat, maxcoll
	integer :: ncoll, nscat
	real(8) :: t
	type(phonon) :: phn
	character(len=8), parameter :: fmt = '(ES16.8)'
	integer, parameter :: unit = 2
	
	phn = emit()
	call inittraj( 2*maxcoll, getpos(phn) )
	
	t = 0
	ncoll = 0
	nscat = 0
	do while (isalive(phn))
		call advect(phn, t)
		call scatter(phn, nscat)
		ncoll = ncoll + 1
		if (ncoll >= maxcoll .or. (maxscat > 0 .and. nscat >= maxscat)) then
			call kill(phn)
		end if
	end do
	
	call writematlab( transpose(gettraj()), '(ES16.8)', 2, 'traj', 'x' )
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
	call writematlab(flux, '(ES16.8)', 3, 'flux', 'q')
end subroutine writeflux

subroutine printcond(deltaT, length)
	real(8), intent(in) :: deltaT, length
	real(8) :: flux, conductivity
	
	flux = getflux()
	conductivity = -flux *length/deltaT
	
	print ('(A,ES16.8,A)'), 'j = ', flux, ' W/m^2'
	print ('(A,ES16.8,A)'), 'k = ', conductivity, ' W/m-K'
end subroutine printcond

end module simulation