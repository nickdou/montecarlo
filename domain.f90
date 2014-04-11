module domain
	use math
	use distributions
	implicit none
	
	private :: tend, tstep, ntime, gridnorm, gridloc, cellsize, ncell, nbdry, nemit, &
		emitind, Eeff, volume, gridtime_arr, gridnum_arr, bdry_arr, emit_arr, getcoll, &
		incgridtime, ntraj, maxtraj, trajectory, appendtraj
	public  :: rectbdry, setend, setspatgrid, getcoord, getind, coordtoind, gettind, &
		makebdry, setbdry, setpair, getemit, geteeff, getemittime, getemitstate, &
		updatestate, recordtime, recordnum, getsteadytemp, inittraj, gettraj
	
	type rectbdry
		private
		integer :: bc, pair
		real(8) :: o(3), n(3), loc, dx, dy, area, rot(3,3), &
			temp, pairrot(3,3), pairmv(3)
	end type rectbdry
	
	integer, parameter :: SPEC_BC = 1, &
	                      DIFF_BC = 2, &
	                      ISOT_BC = 3, &
	                      PERI_BC = 4
	
	real(8) :: tend, tstep, gridnorm(3), gridloc, cellsize, Eeff, volume
	integer :: ntime, ncell, nbdry, nemit, emitind, ntraj, maxtraj
	real(8), allocatable :: trajectory(:,:)
	real(8), allocatable :: gridtime_arr(:)
	integer, allocatable :: gridnum_arr(:,:)
	type(rectbdry), allocatable :: bdry_arr(:)
	integer, allocatable :: emit_arr(:)

contains

subroutine setend(t)
	real(8), intent(in) :: t
	
	tend = t
end subroutine setend

subroutine setspatgrid(norm, lower, upper, num)
	real(8), intent(in) :: norm(3), lower, upper
	integer, intent(in) :: num 
	
	gridnorm = norm
	gridloc = lower
	ncell = num
	cellsize = (upper - lower)/num
	
	allocate( gridtime_arr(num) )
	gridtime_arr = 0
end subroutine setspatgrid

subroutine settimegrid(t, num)
	real(8), intent(in) :: t
	integer, intent(in) :: num
	
	tend = t
	ntime = num
	
	if (num /= 0) then
		tstep = t/num
	
		allocate( gridnum_arr(ncell, 0:num) )
		gridnum_arr = 0
	end if
end subroutine settimegrid

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

integer pure function getind(pos) result(ind)
	real(8), intent(in) :: pos(3)
	
	ind = coordtoind( getcoord(pos) )
end function getind

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

subroutine setbdry(arr, num, vol)
	type(rectbdry), intent(in) :: arr(:)
	integer, intent(in) :: num
	real(8), intent(in) :: vol
	real(8), allocatable :: areatemp_arr(:)
	real(8) :: sumareatemp
	integer :: i
	
	nbdry = size(arr)
	allocate( bdry_arr(nbdry) )
	bdry_arr = arr
	
	allocate ( areatemp_arr(nbdry) )
	areatemp_arr = (/(arr(i)%area*abs(arr(i)%temp), i=1,nbdry)/)
	sumareatemp = sum(areatemp_arr)
	
	allocate( emit_arr(nbdry) )
	emit_arr = nint( num*areatemp_arr/sumareatemp )
!	print *, emit_arr
	nemit = sum(emit_arr)
	emitind = 1
	
	Eeff = (sumareatemp*tend/nemit) * getpseudoflux()
	deallocate(areatemp_arr)
	
	volume = vol
end subroutine setbdry

subroutine setpair(ind1, ind2, rot1, rot2, mv)
! TODO: check that rot2 = inv(rot1) or calculate one from the other
	integer, intent(in) :: ind1, ind2
	real(8), intent(in) :: rot1(3,3), rot2(3,3), mv(3)
	
	bdry_arr(ind1)%bc = PERI_BC
	bdry_arr(ind1)%pair = ind2
	bdry_arr(ind1)%pairrot = rot1
	bdry_arr(ind1)%pairmv = mv
	
	bdry_arr(ind2)%bc = PERI_BC
	bdry_arr(ind2)%pair = ind1
	bdry_arr(ind2)%pairrot = rot2
	bdry_arr(ind2)%pairmv = -mv
end subroutine setpair

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

subroutine getemitstate(pos, dir, sign)
	real(8), intent(out) :: pos(3), dir(3)
	logical, intent(out) :: sign
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
		
		call drawposlin(lin)
		pos = pos - project(pos, gridnorm)
		pos = pos + (gridloc + ncell*cellsize*lin(3))*gridnorm
		
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

subroutine updatestate(x, dir, deltax, t, deltat)
	real(8), intent(inout) :: x(3), dir(3), deltax, t, deltat
	integer :: i, ind, bc
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
	if (ind == 0) then
		x = aim
		call appendtraj(x)
!		print ('(A,F10.5,A,F10.3,A,3F10.5)'), &
!			'dx = ', deltax*1e6, '  dt = ', deltat*1e12, '  x = ', x*1e6
	else
		deltax = dist(ind)
		deltat = deltax / v
		x = coll(:,ind)
		call appendtraj(x)
!		print ('(A,F10.5,A,F10.3,A,3F10.5)'), &
!			'dx = ', deltax*1e6, '  dt = ', deltat*1e12, '  x = ', x*1e6
		
		bc = bdry_arr(ind)%bc
		select case (bc)
		case (SPEC_BC)
			dir = dir - 2*project(dir, bdry_arr(ind)%n)
!			print ('(A,I1)'), 'spec ', ind
		case (DIFF_BC)
			call drawanghalf(dir)
			dir = matmul(bdry_arr(ind)%rot, dir)
		case (ISOT_BC)
!			print ('(A,I1)'), 'isot ', ind
			dir = 0
		case (PERI_BC)
			x = x + bdry_arr(ind)%pairmv
			dir = matmul(bdry_arr(ind)%pairrot, dir)
			call appendtraj(x)
		end select
	end if
	
	t = t + deltat
	if (t >= tend) then
!		print *, 'tend'
		dir = 0
	end if
end subroutine updatestate

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
	integer :: pm, tind, xind, told, tnew
	
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
				xind = getind(xold*xi + xnew*(1 - xi))
				gridnum_arr(xind, tind) = gridnum_arr(xind, tind) + pm
			end do
		end if
		
!		if (t >= tend) then
!			tind = ntime
!		else
!			tind = ceiling(t/tstep) - 1
!		end if
!		xi = (t - tind*tstep)/deltat
!	
!		do while (xi <= 1.000001)
!			xind = getind(xold*xi + xnew*(1 - xi))
!			gridnum_arr(xind, tind) = gridnum_arr(xind, tind) + pm
!		
!			tind = tind - 1
!			xi = (t - tind*tstep)/deltat
!		end do
	end if
end subroutine recordnum

function getsteadytemp() result(T)
	real(8) :: T(ncell)
	
	T = Eeff*ncell/(volume*tend*getpseudoenergy()) * gridtime_arr
end function getsteadytemp

function gettranstemp() result(T)
	real(8) :: T(ncell, 0:ntime)
	
	T = Eeff*ncell/(volume*getpseudoenergy()) * gridnum_arr
end function gettranstemp

end module domain