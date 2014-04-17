module domain
	use math
	use distributions
	implicit none
	
	private :: tend, tstep, ntime, gridnorm, gridloc, cellsize, ncell, nbdry, nemit, &
		emitind, Eeff, volume, gridtime_arr, incgridtime, gridnum_arr, cumflux, &
		bdry_arr, emit_arr, getcoll, ntraj, maxtraj, trajectory, &
		addcumflux_real, addcumflux_int
	public  :: rectbdry, setgrid, initrecord, getcoord, getxind, coordtoind, gettind, &
		makebdry, setbdry, setbdrypair, getemit, geteeff, getemittime, getemitstate, &
		updatestate, recordtime, recordnum, getsteadytemp, inittraj, appendtraj, gettraj, &
		addcumflux, getflux
	
	type rectbdry
		private
		integer :: bc, pair
		real(8) :: o(3), n(3), loc, dx, dy, area, rot(3,3), &
			temp, pairrot(3,3), pairmv(3)
	end type rectbdry
	
	interface addcumflux
		module procedure addcumflux_real, addcumflux_int
	end interface addcumflux
	
	integer, parameter :: SPEC_BC = 1, &
	                      DIFF_BC = 2, &
	                      ISOT_BC = 3, &
	                      PERI_BC = 4
	
	real(8) :: tend, tstep, gridnorm(3), gridloc, cellsize, Eeff, volume, cumflux
	integer :: ntime, ncell, nbdry, nemit, emitind, ntraj, maxtraj
	real(8), allocatable :: trajectory(:,:)
	real(8), allocatable :: gridtime_arr(:)
	integer, allocatable :: gridnum_arr(:,:)
	type(rectbdry), allocatable :: bdry_arr(:)
	integer, allocatable :: emit_arr(:)

contains

!subroutine setend(t)
!	real(8), intent(in) :: t
!	
!	tend = t
!end subroutine setend
!
!subroutine setspatgrid(norm, lower, upper, num)
!	real(8), intent(in) :: norm(3), lower, upper
!	integer, intent(in) :: num
!	
!	gridnorm = norm
!	gridloc = lower
!	ncell = num
!	cellsize = (upper - lower)/num
!	
!	cumflux = 0
!	
!	allocate( gridtime_arr(num) )
!	gridtime_arr = 0
!end subroutine setspatgrid
!
!subroutine settimegrid(t, num)
!	real(8), intent(in) :: t
!	integer, intent(in) :: num
!	
!	tend = t
!	ntime = num
!	
!	if (num /= 0) then
!		tstep = t/num
!	
!		allocate( gridnum_arr(ncell, 0:num) )
!		gridnum_arr = 0
!	end if
!end subroutine settimegrid

subroutine setgrid(t, tnum, norm, lower, upper, xnum)
	real(8), intent(in) :: t, norm(3), lower, upper
	integer, intent(in) :: tnum, xnum
	
	tend = t
	ntime = tnum
	
	gridnorm = norm
	gridloc = lower
	ncell = xnum
	cellsize = (upper - lower)/xnum
end subroutine setgrid

subroutine initrecord()
	if (ntime /= 0) then
		tstep = tend/ntime
		
		allocate( gridnum_arr(ncell, 0:ntime) )
		gridnum_arr = 0
	end if
	
	cumflux = 0
	
	allocate( gridtime_arr(ncell) )
	gridtime_arr = 0
end subroutine initrecord

subroutine inittraj(num, pos)
	integer, intent(in) :: num
	real(8), intent(in) :: pos(3)
	
	maxtraj = num
	ntraj = 0
	
	allocate( trajectory(3, 0:num) )
	trajectory(:,0) = pos
end subroutine inittraj

subroutine appendtraj(pos)
	real(8), intent(in) :: pos(3)
	
	if (allocated(trajectory)) then
		if (ntraj < maxtraj) then
			ntraj = ntraj + 1
			trajectory(:,ntraj) = pos
		end if
	end if
end subroutine appendtraj

function gettraj() result(traj)
	real(8), allocatable :: traj(:,:)
	
	allocate( traj(3, 0:ntraj) )
	traj = trajectory(:, 0:ntraj)
end function gettraj

real(8) pure function getcoord(pos) result(coord)
	real(8), intent(in) :: pos(3)
	
	coord = (dot_product(pos, gridnorm) - gridloc)/cellsize
end function getcoord

integer pure function coordtoind(coord) result(ind)
	real(8), intent(in) :: coord
	
	ind = min(ncell, max(1, ceiling(coord)))
end function coordtoind

integer pure function getxind(pos) result(ind)
	real(8), intent(in) :: pos(3)
	
	ind = coordtoind( getcoord(pos) )
end function getxind

integer pure function gettind(time) result(tind)
	real(8), intent(in) :: time
	
	tind = min(ntime, max(0, floor(time/tstep)))
end function gettind

type(rectbdry) pure function makebdry(bc, o, xvec, yvec, temp) result(bdry)
	integer, intent(in) :: bc
	real(8), intent(in) :: o(3), xvec(3), yvec(3)
	real(8), intent(in), optional :: temp
	real(8) :: cross(3), xaxis(3), yaxis(3), zaxis(3)
	
	if (dot_product(xvec, yvec) == 0) then
		bdry%bc = bc
		bdry%o = o
		bdry%dx = normtwo(xvec)
		bdry%dy = normtwo(yvec)
		
		cross = cross_product(xvec, yvec)
		xaxis = unitvec(xvec)
		yaxis = unitvec(yvec)
		zaxis = unitvec(cross)
	
		bdry%area = normtwo(cross)
		bdry%n = zaxis
		bdry%loc = dot_product(o, zaxis)
		bdry%rot = rotmatrix(xaxis, yaxis, zaxis)
	
		if (bc == ISOT_BC .or. bc == PERI_BC) then
			bdry%temp = temp
		else
			bdry%temp = 0
		end if
	end if
end function makebdry

subroutine setbdry(arr, vol)
	type(rectbdry), intent(in) :: arr(:)
	real(8), intent(in) :: vol
	
	nbdry = size(arr)
	allocate( bdry_arr(nbdry) )
	bdry_arr = arr
	
	volume = vol
end subroutine setbdry

subroutine setbdrypair(ind1, ind2, rot1, rot2, mv)
! TODO: check that rot2 = inv(rot1) or calculate one from the other
	integer, intent(in) :: ind1, ind2
	real(8), intent(in) :: rot1(3,3), rot2(3,3), mv(3)
	real(8) :: temp1, temp2
	
	temp1 = bdry_arr(ind1)%temp
	temp2 = bdry_arr(ind2)%temp
	
	bdry_arr(ind1)%bc = PERI_BC
	bdry_arr(ind1)%pair = ind2
	bdry_arr(ind1)%pairrot = rot1
	bdry_arr(ind1)%pairmv = mv
	bdry_arr(ind1)%temp = temp1 - temp2
	
	bdry_arr(ind2)%bc = PERI_BC
	bdry_arr(ind2)%pair = ind1
	bdry_arr(ind2)%pairrot = rot2
	bdry_arr(ind2)%pairmv = -mv
	bdry_arr(ind2)%temp = temp2 - temp1
end subroutine setbdrypair

subroutine calculateemit(num)
	integer, intent(in) :: num
	real(8) :: areatemp_arr(nbdry)
	real(8) :: sumareatemp
	integer :: i
	
	areatemp_arr = (/(bdry_arr(i)%area*abs(bdry_arr(i)%temp), i=1,nbdry)/)
	sumareatemp = sum(areatemp_arr)
	
	if (.not. allocated(emit_arr)) then
		allocate( emit_arr(nbdry) )
	else if (size(emit_arr, 1) /= nbdry) then
		deallocate( emit_arr )
		allocate( emit_arr(nbdry) )
	end if
	emit_arr = nint( num*areatemp_arr/sumareatemp )
	nemit = sum(emit_arr)
	emitind = 1
	
	Eeff = (sumareatemp*tend/nemit) * getpseudoflux()
end subroutine calculateemit

pure function getemit() result(arr)
	integer :: arr(nbdry)
	
	arr = emit_arr
end function getemit

real(8) pure function geteeff() result(E)
	E = Eeff
end function geteeff

subroutine getemittime(time)
	real(8), intent(out) :: time
	real(8) :: r
	
	if (ntime == 0) then
		time = 0
	else
		call random_number(r)
		time = r*tend
	end if
end subroutine getemittime

subroutine getemitstate(pos, dir, sign, volumetric)
	real(8), intent(out) :: pos(3), dir(3)
	logical, intent(out) :: sign
	logical, intent(in) :: volumetric
	integer :: ind
	real(8) :: lin(3)
	
	if (nemit > 0) then
		do
			if (emit_arr(emitind) > 0) then
				emit_arr(emitind) = emit_arr(emitind) - 1
				nemit = nemit - 1
				exit
			else
				emitind = emitind + 1
			end if
		end do
	
		call drawposrect(pos)
		pos = (/bdry_arr(emitind)%dx, bdry_arr(emitind)%dy, 1d0/)*pos
		pos = bdry_arr(emitind)%o + matmul(bdry_arr(emitind)%rot, pos)
		
		if (volumetric) then
			call drawposlin(lin)
			pos = pos - project(pos, gridnorm)
			pos = pos + (gridloc + ncell*cellsize*lin(3))*gridnorm
		end if
		
		call drawanghalf(dir)
		dir = matmul(bdry_arr(emitind)%rot, dir)
		
		sign = (bdry_arr(emitind)%temp >= 0)
	end if
end subroutine getemitstate

pure subroutine getcoll(pierced, coll, pos1, pos2, dx, dy)
	real(8), intent(in) :: pos1(3), pos2(3), dx, dy
	logical, intent(out) :: pierced
	real(8), intent(out) :: coll(3)
	real(8) :: z1, z2
	
	z1 = pos1(3)
	z2 = pos2(3)
	coll = (z2*pos1 - z1*pos2)/(z2 - z1)
	
	pierced = (z1*z2 < 0) .and. &
		(coll(1) >= 0 .and. coll(1) <= dx) .and. &
		(coll(2) >= 0 .and. coll(2) <= dy)
end subroutine getcoll

subroutine updatestate(ind, bc, x, dir, deltax, t, deltat)
	real(8), intent(inout) :: x(3), dir(3), deltax, t, deltat
	integer, intent(out) :: ind, bc
	integer :: i
	real(8) :: v, aim(3), o(3), dx, dy, rot(3,3), rotT(3,3), pos1(3), pos2(3)
	logical :: piercedi, pierced(nbdry)
	real(8) :: colli(3), coll(3,nbdry), dist(nbdry)
	
	v = deltax/deltat
	if (t + deltat > tend) then
		deltat = tend - t
		deltax = v*deltat
	end if
	aim = x + deltax*dir
	
	do i = 1,nbdry
		o = bdry_arr(i)%o
		dx = bdry_arr(i)%dx
		dy = bdry_arr(i)%dy
		rot = bdry_arr(i)%rot
		rotT = transpose(rot)
		
		pos1 = matmul(rotT, x - o)
		pos2 = matmul(rotT, aim - o)
		call getcoll(piercedi, colli, pos1, pos2, dx, dy)
		
		pierced(i) = piercedi
		coll(:,i) = o + matmul(rot, colli)
		dist(i) = normtwo(colli - pos1)
	end do
	
	ind = minloc(dist, 1, pierced)
	bc = 0
	if (ind == 0) then
		x = aim
	else
		bc = bdry_arr(ind)%bc
		deltax = dist(ind)
		deltat = deltax / v
		x = coll(:,ind)
	end if
	
	t = t + deltat
	if (t >= tend) then
		bc = 0
		dir = 0
	end if
end subroutine updatestate

subroutine applybc(ind, bc, x, dir)
	integer, intent(in) :: ind, bc
	real(8), intent(inout) :: x(3), dir(3)
	
	if (bc /= 0 .and. any(dir /= 0, 1)) then
		select case (bc)
		case (SPEC_BC)
			dir = dir - 2*project(dir, bdry_arr(ind)%n)
		case (DIFF_BC)
			call drawanghalf(dir)
			dir = matmul(bdry_arr(ind)%rot, dir)
		case (ISOT_BC)
			dir = 0
		case (PERI_BC)
			x = x + bdry_arr(ind)%pairmv
			dir = matmul(bdry_arr(ind)%pairrot, dir)
		end select
	end if
end subroutine applybc

subroutine incgridtime(time, ind1, ind2)
	real(8), intent(in) :: time
	integer, intent(in) :: ind1
	integer, intent(in), optional :: ind2
	
	if (present(ind2)) then
		if (ind1 <= ind2) then
			gridtime_arr(ind1:ind2) = gridtime_arr(ind1:ind2) + time
		end if
	else
		gridtime_arr(ind1) = gridtime_arr(ind1) + time
	end if
end subroutine incgridtime

subroutine recordtime(sign, deltat, xold, xnew)
	logical, intent(in) :: sign
	real(8), intent(in) :: deltat, xold(3), xnew(3)
	real(8) :: coordold, coordnew, dcoord
	integer :: pm, indold, indnew
	
	if (ntime == 0) then
		pm = signtoint(sign)
	
		coordold = getcoord(xold)
		coordnew = getcoord(xnew)
		dcoord = abs(coordnew - coordold)
	
		indold = coordtoind(coordold)
		indnew = coordtoind(coordnew)
	
		if (indold < indnew) then
			call incgridtime(pm*deltat*(indold - coordold)/dcoord, indold)
			call incgridtime(pm*deltat/dcoord, indold+1, indnew-1)
			call incgridtime(pm*deltat*(1 - (indnew - coordnew))/dcoord, indnew)
		else if (indnew < indold) then
			call incgridtime(pm*deltat*(indnew - coordnew)/dcoord, indnew)
			call incgridtime(pm*deltat/dcoord, indnew+1, indold-1)
			call incgridtime(pm*deltat*(1 - (indold - coordold))/dcoord, indold)
		else
			call incgridtime(pm*deltat, indold)
		end if
	
!		print ('(2I10)'), indold, indnew
!		print ('(I4,F8.3,I4)'), floor(coordold), coordold, ceiling(coordold)
!		print ('(2F10.3,A,2ES11.3,A,2ES11.3)'), coordold, coordnew, ' :', &
!			gridtime_arr(1:2), ' ... ', gridtime_arr(ncell-1:ncell)
	end if
end subroutine recordtime

subroutine recordnum(sign, t, deltat, xold, xnew)
	logical, intent(in) :: sign
	real(8), intent(in) :: t, deltat, xold(3), xnew(3)
	real(8) :: xi, coord
	integer :: pm, told, tnew, tind, xind
	
	if (ntime /= 0) then
		pm = signtoint(sign)
		
		told = gettind(t - deltat)
		tnew = gettind(t)
		if (told /= 0) then
			told = told + 1
		end if
		
		if (told <= tnew) then
			do tind = told, tnew
				xi = (t - tind*tstep)/deltat
				xind = getxind(xold*xi + xnew*(1 - xi))
				gridnum_arr(xind, tind) = gridnum_arr(xind, tind) + pm
			end do
		end if
	end if
end subroutine recordnum

subroutine addcumflux_real(sign, deltax)
	logical, intent(in) :: sign
	real(8), intent(in) :: deltax(3)
	integer :: pm
	
	pm = signtoint(sign)
	cumflux = cumflux + pm*dot_product(gridnorm, deltax)
end subroutine addcumflux_real

subroutine addcumflux_int(sign, ind)
	logical, intent(in) :: sign
	integer, intent(in) :: ind
	integer :: pm
	
	pm = signtoint(sign)
	cumflux = cumflux - pm*dot_product(gridnorm, bdry_arr(ind)%pairmv)
end subroutine addcumflux_int

pure function getsteadytemp() result(T)
	real(8) :: T(ncell)
	
	T = Eeff*ncell/(volume*tend*getpseudoenergy()) * gridtime_arr
end function getsteadytemp

pure function gettranstemp() result(T)
	real(8) :: T(ncell, 0:ntime)
	
	T = Eeff*ncell/(volume*getpseudoenergy()) * gridnum_arr
end function gettranstemp

real(8) pure function getflux() result(flux)
	flux = Eeff/(volume*tend) * cumflux
end function getflux

end module domain