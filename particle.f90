module particle
	use material
	use distributions
	use domain
	implicit none
	
	public
	
	type phonon
		private
		logical :: sign, alive
		real(8) :: omega, x(3), dir(3), tscat
		integer :: p
	end type phonon
contains

logical pure function isalive(phn) result(alive)
	type(phonon), intent(in) :: phn
	
	alive = phn%alive
end function isalive

pure subroutine kill(phn)
	type(phonon), intent(inout) :: phn
	
	phn%alive = .false.
end subroutine kill

real(8) pure function getcoord_phn(phn) result(coord)
	type(phonon), intent(in) :: phn
	
	coord = getcoord(phn%x)
end function getcoord_phn

type(phonon) function emitbdry() result(phn)
	
	call drawfluxprop(phn%omega, phn%p)
	call getemitstate(phn%x, phn%dir, phn%sign)
	call drawscattime(phn%tscat, phn%omega, phn%p)
	phn%alive = .true.
end function emitbdry

subroutine scatter(phn)
	type(phonon), intent(inout) :: phn
	
	if (phn%tscat == 0) then
		call drawscatterprop(phn%omega, phn%p)
		call drawangiso(phn%dir)
		call drawscattime(phn%tscat, phn%omega, phn%p)
	end if
end subroutine scatter

subroutine advect(phn, t)
	type(phonon), intent(inout) :: phn
	real(8), intent(inout) :: t
	real(8) :: x(3), dir(3), deltax, deltat
	
	x = phn%x
	dir = phn%dir
	deltax = vel(phn%omega, phn%p)*phn%tscat
	deltat = phn%tscat
	call updatestate(x, dir, deltax, t, deltat)
	
	if (phn%sign) then
		call recordtime(deltat, phn%x, x)
	else
		call recordtime(-deltat, phn%x, x)
	end if
	
	phn%x = x
	phn%dir = dir
	phn%tscat = phn%tscat - deltat
	
	if (all(phn%dir == 0, 1)) then
		call kill(phn)
	end if
end subroutine advect

end module particle