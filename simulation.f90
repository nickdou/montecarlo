module simulation
	use distributions
	use domain
	use particle
	implicit none
	
	public
	
contains

subroutine initboxisot(nomega, nemit, ncell, length, side, tend, Teq, Thot, Tcold)
	integer, intent(in) :: nomega, nemit, ncell
	real(8), intent(in) :: length, side, tend, Teq, Thot, Tcold
	type(rectbdry) :: bdry_arr(6)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
	
	call initrand()
	
	call initomega(nomega)
	call initpropcdf(Teq)
	
	call setend(tend)
	call setgrid((/0d0, 0d0, 1d0/), 0d0, length, ncell)
	
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

subroutine simulate(len)
	integer, intent(in) :: len
	integer :: i, nemit
	real(8) :: t
	type(phonon) :: phn
	
	nemit = sum( getemit() )
	do i = 1, nemit
		call showprogress(i, nemit, len)
!		print ('(/,A,I2)'), 'Particle ', i 
		phn = emitbdry()
		t = 0
		do
			if (.not. isalive(phn)) exit
			call scatter(phn)
			call advect(phn, t)
		end do
	end do
end subroutine simulate

subroutine showprogress(i, iend, len)
	integer, intent(in) :: i, iend, len
	real(8) :: x, unit
	integer :: bars
	
	unit = dble(len)/iend
	x = i*unit
	bars = floor(x)
	if ((x - bars)/unit < 0.999999) then
		print *, nint(100*dble(i)/iend), '%  ', &
			'[', repeat('|',bars), repeat('-',len-bars), ']'
			
	end if
end subroutine showprogress

subroutine printtemp(ncell)
	integer, intent(in) :: ncell
	real(8) :: T(ncell)
	
	T = getsteadytemp()
	open(unit=2, file='temp.m', action='write', status='replace')
	write(2,'(A)') 'T = ['
	write(2,'(ES15.8)') T
	write(2,'(A)') '];'
	write(2,'(A)') 'plot(T)'
!	print ('(ES15.8)'), T
end subroutine printtemp

end module simulation