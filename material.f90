module material
	use constants
	implicit none
	
	private :: c, tauconst
	public  :: npol, omegamax, vel, dos, tau
	
	integer, parameter :: npol = 1
	real(8), parameter :: c = 6000 ! m/s for Si
	real(8), parameter :: omegamax = 2*pi*10d12 ! 10 THz for Si
	real(8), parameter :: tauconst = 10d-12 ! 10 ps
	
contains

real(8) pure function vel(omega, p)
	real(8), intent(in) :: omega
	integer, intent(in) :: p
	
	vel = c
end function vel

real(8) pure function dos(omega, p)
	real(8), intent(in) :: omega
	integer, intent(in) :: p
	
	if (omega <= omegamax) then
		dos = 3*omega**2 / (2*pi**2*c**3)
	else
		dos = 0
	end if
end function dos

real(8) pure function tau(omega, p, T)
	real(8), intent(in) :: omega, T
	integer, intent(in) :: p
	
	tau = tauconst
end function tau

end module material