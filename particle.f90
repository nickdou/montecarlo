module particle
	use material
	use distributions
	use domain
	implicit none
	
	public
	
	type phonon
		private
		logical :: sign, alive
		real(8) :: omega, xinitial(3), x(3), dir(3), tscat
		integer :: p
	end type phonon
contains

logical pure function isalive(phn) result(alive)
	type(phonon), intent(in) :: phn
	
	alive = phn%alive
end function isalive

subroutine kill(phn)
	type(phonon), intent(inout) :: phn
	
	call addcumdisp(phn%sign, phn%x - phn%xinitial)
	phn%alive = .false.
end subroutine kill

pure function getpos(phn) result(x)
	type(phonon), intent(in) :: phn
	real(8) :: x(3)
	
	x = phn%x
end function getpos

type(phonon) function emit() result(phn)
	call drawfluxprop(phn%omega, phn%p)
	call drawemitstate(phn%x, phn%dir, phn%sign)
	call drawscattime(phn%tscat, phn%omega, phn%p)
	phn%xinitial = phn%x
	phn%alive = .true.
end function emit

subroutine scatter(phn, nscat)
	type(phonon), intent(inout) :: phn
	integer, intent(inout) :: nscat
	
	if (phn%tscat == 0) then
		call drawscatterprop(phn%omega, phn%p)
		call drawangiso(phn%dir)
		call drawscattime(phn%tscat, phn%omega, phn%p)
		
		nscat = nscat + 1
	end if
end subroutine scatter

subroutine advect(phn, t)
	type(phonon), intent(inout) :: phn
	real(8), intent(inout) :: t
	real(8) :: x(3), dir(3), deltax, deltat
	integer :: ind, bc
	
	x = phn%x
	dir = phn%dir
	deltax = vel(phn%omega, phn%p)*phn%tscat
	deltat = phn%tscat
	
	call updatestate(ind, bc, x, dir, deltax, t, deltat)
	call appendtraj(x)
	
	call recordtime(phn%sign, deltat, phn%x, x)
	call recorddisp(phn%sign, phn%x, x)
	call recordloc(phn%sign, t, deltat, phn%x, x)
	
	call applybc(ind, bc, x, dir)
	if (bc == PERI_BC) then
		call appendtraj(x)
		call addcumdisp(phn%sign, ind)
	end if
	
	phn%x = x
	phn%dir = dir
	phn%tscat = phn%tscat - deltat
	
	if (all(dir == 0, 1)) then
		call kill(phn)
	end if
end subroutine advect

end module particle