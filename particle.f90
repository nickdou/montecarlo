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

pure function getpos(phn) result(x)
	type(phonon), intent(in) :: phn
	real(8) :: x(3)
	
	x = phn%x
end function getpos

!real(8) pure function getcoord_phn(phn) result(coord)
!	type(phonon), intent(in) :: phn
!	
!	coord = getcoord(phn%x)
!end function getcoord_phn

type(phonon) function emitbdry() result(phn)
	
	call drawfluxprop(phn%omega, phn%p)
	call getemitstate(phn%x, phn%dir, phn%sign)
	call drawscattime(phn%tscat, phn%omega, phn%p)
	phn%alive = .true.
end function emitbdry

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
	
	x = phn%x
	dir = phn%dir
	deltax = vel(phn%omega, phn%p)*phn%tscat
	deltat = phn%tscat
	call updatestate(x, dir, deltax, t, deltat)
	
	call recordtime(phn%sign, deltat, phn%x, x)
	call recordnum(phn%sign, t, deltat, phn%x, x)
	
	phn%x = x
	phn%dir = dir
	phn%tscat = phn%tscat - deltat
	
	if (all(phn%dir == 0, 1)) then
		call kill(phn)
	end if
end subroutine advect

end module particle