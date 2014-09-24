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

subroutine initisot1d(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
	character(len=*), intent(in) :: disp, relax
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, T, Thot, Tcold
	integer, parameter :: nbdry = 6
	type(axis) :: grid, vgen
	type(boundary) :: bdry_arr(nbdry)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
	integer :: emit_arr(nbdry)

	call initrand()

	call initmat(disp, relax, T)
	call initdist()

	grid = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, ncell)
	call setgrid(tend, ntime, grid)
	call initrecord(gf, (/0d0, 0d0, 1d0/))

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
	call setbdry(bdry_arr, product(corner))
	call calculateemit(nemit)

	vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
	call setvolumetric(vol, vgen)
	
	emit_arr = getemit()
	print ('(/,A)'), 'Isothermal'
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	print ('(A12,ES15.8)'), 'ktheory = ', ktheory
	print ('(A12,16I10)'),  'emit_arr = ', emit_arr
end subroutine initisot1d

subroutine initbulk1d(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
	character(len=*), intent(in) :: disp, relax
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, T, Thot, Tcold
	real(8) :: zvec(3)

	call initisot1d(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)

	zvec = (/0d0, 0d0, length/)
	call setbdrypair(1, 4, zvec)
	call calculateemit(nemit)

	print ('(/,A)'), 'Bulk'
end subroutine initbulk1d

subroutine initfilm(disp, relax, vol, gf, nemit, ncell, ntime, length, side, tend, T, Thot, Tcold)
	character(len=*), intent(in) :: disp, relax
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nemit, ncell, ntime
	real(8), intent(in) :: length, side, tend, T, Thot, Tcold
	integer, parameter :: nbdry = 6
	type(axis) :: grid, vgen
	type(boundary) :: bdry_arr(nbdry)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3)
	integer :: emit_arr(nbdry)

	call initrand()

	call initmat(disp, relax, T)
	call initdist()

	grid = makeaxis((/0d0, 1d0, 0d0/), 0d0, side, ncell)
	call setgrid(tend, ntime, grid)
	call initrecord(gf, (/0d0, 0d0, 1d0/))

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
	call setbdry(bdry_arr, product(corner))
	call setbdrypair(1, 4, zvec)
	call calculateemit(nemit)

	vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
	call setvolumetric(vol, vgen)

	emit_arr = getemit()
	print ('(/,A)'), 'Film'
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	print ('(A12,ES15.8)'), 'ktheory = ', ktheory
	print ('(A12,16I10)'),   'emit_arr = ', emit_arr
end subroutine initfilm

subroutine inithollow(disp, relax, vol, gf, nemit, ncell, ntime, length, side, wall, tend, T, Thot, Tcold)
	character(len=*), intent(in) :: disp, relax
	logical, intent(in) :: vol, gf
	integer, intent(in) :: nemit, ncell, ntime
	real(8), intent(in) :: length, side, wall, tend, T, Thot, Tcold
	integer, parameter :: nbdry = 16
	type(axis) :: grid, vgen
	type(boundary) :: bdry_arr(nbdry)
	real(8) :: origin(3), corner(3), xvec(3), yvec(3), zvec(3), xwall(3), ywall(3)
	integer :: emit_arr(nbdry)

	call initrand()

	call initmat(disp, relax, T)
	call initdist()

	grid = makeaxis((/0d0, 1d0, 0d0/), 0d0, side, ncell)
	call setgrid(tend, ntime, grid)
	call initrecord(gf, (/0d0, 0d0, 1d0/))

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

	call setbdry(bdry_arr, product(corner) - product(corner - 2*xwall - 2*ywall))
	call setbdrypair(1, 7, zvec)
	call setbdrypair(2, 8, zvec)
	call setbdrypair(3, 5, zvec)
	call setbdrypair(4, 6, zvec)
	call calculateemit(nemit)

	vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, length, 0)
	call setvolumetric(vol, vgen)

	emit_arr = getemit()
	print ('(/,A)'), 'Hollow'
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	print ('(A12,ES15.8)'), 'ktheory = ', ktheory
	print ('(A12,16I10)'),   'emit_arr = ', emit_arr
end subroutine inithollow

subroutine initunit(disp, relax, nemit, ntime, a, b, c, d, tend, T, Thot, Tcold)
	character(len=*), intent(in) :: disp, relax
	integer, intent(in) :: nemit, ntime
	real(8), intent(in) :: a, b, c, d, tend, T, Thot, Tcold
	integer, parameter :: nbdry = 24
	type(axis) :: grid, vgen
	type(boundary) :: bdry_arr(nbdry)
	real(8) :: l, bdrydata_arr(nbdry, 10)
	integer :: i, bc_arr(nbdry), emit_arr(nbdry)

	call initrand()

	call initmat(disp, relax, T)
	call initdist()

	l = c/4d0*sqrt(2d0)

	grid = makeaxis((/0d0, 0d0, 1d0/), -l, 0d0, 1) ! ncell = 1
	call setgrid(tend, ntime, grid)
	call initrecord(.false., (/0d0, 0d0, 1d0/)) ! gf = .false.

	bdrydata_arr = transpose(reshape((/ &
		   -d, 0d0,   -l,       d,    0d0,  0d0,   0d0,      b,   0d0,  Thot, &
		   -d,   b,   -l, 2*a+2*d,    0d0,  0d0,   0d0,      d,   0d0,  Thot, &
		  2*a, 0d0,   -l,       d,    0d0,  0d0,   0d0,      b,   0d0,  Thot, &
		   -d, 0d0,   -l,     0d0,    b+d,  0d0,   0d0,    0d0,   l-a,   0d0, &
		   -d, b+d,   -l,   a+2*d,    0d0,  0d0,   0d0,    0d0, l-a-d,   0d0, &
		  a+d, b+d,   -l,       a,    0d0,  0d0,   0d0,    0d0,     l,   0d0, &
		2*a+d, b+d,   -l,     0d0,   -b-d,  0d0,   0d0,    0d0,     l,   0d0, &
		2*a+d, 0d0,   -l,      -d,    0d0,  0d0,   0d0,    0d0,     l,   0d0, &
		  2*a, 0d0,   -l,     0d0,      b,  0d0,   0d0,    0d0,     l,   0d0, &
		  2*a,   b,   -l,      -a,    0d0,  0d0,   0d0,    0d0,     l,   0d0, &
		    a,   b,   -l,      -a,    0d0,  0d0,   0d0,    0d0,   l-a,   0d0, &
		  0d0,   b,   -l,     0d0,     -b,  0d0,   0d0,    0d0,   l-a,   0d0, &
		  0d0, 0d0,   -l,      -d,    0d0,  0d0,   0d0,    0d0,   l-a,   0d0, &
		  0d0, 0d0,   -a,      -d,    0d0,  0d0,   0d0,      b,   0d0,   0d0, &
		   -d,   l,   -a,     0d0, -l+b+d,  0d0,   0d0,    0d0,    -d,   0d0, &
		   -d,   l, -a-d,     0d0, -l+b+d,  0d0, a+2*d,    0d0,   0d0,   0d0, &
		   -d,   l, -a-d,     a+d,    0d0,  0d0,   0d0,    0d0,     d,   0d0, &
		   -d,   l,   -a,     a+d,    0d0,  0d0,   0d0,   -l+b,   0d0,   0d0, &
		    a,   l,   -a,     0d0,    0d0,    a,   0d0,   -l+b,   0d0,   0d0, &
		    a,   l,  0d0,     0d0,    0d0, -a-d,     d,    0d0,   0d0,   0d0, &
		  a+d,   l,  0d0,     0d0,    0d0, -a-d,   0d0, -l+b+d,   0d0,   0d0, &
		    a,   l,  0d0,       d,    0d0,  0d0,   0d0,   -l+b,   0d0, Tcold, &
		  a+d,   b,  0d0,     0d0,      d,  0d0,     a,    0d0,   0d0, Tcold, &
		  2*a,   b,  0d0,       d,    0d0,  0d0,   0d0,     -b,   0d0, Tcold  &
		/), (/10, nbdry/)))

	bc_arr = (/ ISOT_BC, ISOT_BC, ISOT_BC, &
		DIFF_BC, DIFF_BC, DIFF_BC, DIFF_BC, SPEC_BC, &
		DIFF_BC, DIFF_BC, DIFF_BC, DIFF_BC, SPEC_BC, SPEC_BC, &
		SPEC_BC, DIFF_BC, SPEC_BC, DIFF_BC, DIFF_BC, SPEC_BC, DIFF_BC, &
		ISOT_BC, ISOT_BC, ISOT_BC /)

	do i = 1,nbdry
		bdry_arr(i) = makebdry(bc_arr(i), bdrydata_arr(i,1:3), bdrydata_arr(i,4:6), bdrydata_arr(i,7:9), bdrydata_arr(i,10))
	end do

	call setbdry(bdry_arr, d*(a+b+d)*(2*l-a) + 2*d*(a+d)*(l-b-d))

	call calculateemit(nemit)

	! No volumetric generation
	vgen = makeaxis((/0d0, 0d0, 1d0/), 0d0, 1d0, 0)
	call setvolumetric(.false., vgen)
	
	emit_arr = getemit()
	print ('(/,A)'), 'Unit'
	print ('(A12,ES15.8)'), 'Eeff = ', geteeff()
	print ('(A12,ES15.8)'), 'ktheory = ', ktheory
	print ('(A12,16I10)'),   'emit_arr = ', emit_arr
end subroutine initunit

subroutine simulate(nemit, maxscat)
	integer, intent(in) :: nemit, maxscat
	integer :: i, nscat
	real(8) :: t
	type(phonon) :: phn
	
	do i = 1, nemit
! 		print ('(/,A,I2)'), 'Particle ', i
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

! 	call writematlab( transpose(gettraj()), '(ES16.8)', 2, 'traj', 'x' )
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

real(8) function getcond(deltaT, length) result(k)
	real(8), intent(in) :: deltaT, length
	real(8) :: flux, cond
	
	flux = getflux()
	cond = -flux *length/deltaT
	k = cond
	
	print ('(A,ES16.8,A)'), 'j = ', flux, ' W/m^2'
	print ('(A,ES16.8,A)'), 'k = ', cond, ' W/m-K'
end function getcond

end module simulation