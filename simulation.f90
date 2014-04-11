module simulation
	use tools
	use distributions
	use domain
	use particle
	implicit none
	
	public
	
contains

subroutine initboxisot(nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	integer, intent(in) :: nomega, nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, Teq, Thot, Tcold
	type(rectbdry) :: bdry_arr(6)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
	
	call initrand()
	
	call initomega(nomega)
	call initpropcdf(Teq)
	
	call setspatgrid((/0d0, 0d0, 1d0/), 0d0, length, ncell)
	call settimegrid(tend, ntime)
	
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
	call setbdry(bdry_arr, nemit, side**2*length)
	
	print ('(A12,6I10)'), 'emit_arr = ', getemit()
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	
end subroutine initboxisot

subroutine initboxperi(nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	integer, intent(in) :: nomega, nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, Teq, Thot, Tcold
	real(8) :: zvec(3), eye(3,3)
	
	call initboxisot(nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	
	zvec = (/0d0, 0d0, length/)
	eye = reshape((/1d0, 0d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0, 1d0/), (/3, 3/))
	call setpair(1, 4, eye, eye, zvec)
	
end subroutine initboxperi

subroutine simulate(maxscat)
	integer, intent(in) :: maxscat
	integer :: i, nemit, nscat
	real(8) :: t
	type(phonon) :: phn
	
	nemit = sum( getemit() )
	do i = 1, nemit
		call showprogress(i, nemit, min(nemit, 20))
!		print ('(/,A,I2)'), 'Particle ', i 
		phn = emitbdry()
		
		call getemittime(t)
		nscat = 0
		do while (isalive(phn))
			call advect(phn, t)
			call scatter(phn, nscat)
			if (maxscat > 0 .and. nscat >= maxscat) then
				call kill(phn)
			end if
		end do
	end do
end subroutine simulate

subroutine simulateone(maxscat, maxcoll)
	integer, intent(in) :: maxscat, maxcoll
	integer :: ncoll, nscat
	real(8) :: t
	type(phonon) :: phn
	character(len=8), parameter :: fmt = '(ES16.8)'
	integer, parameter :: unit = 2
	
	phn = emitbdry()
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

end module simulation