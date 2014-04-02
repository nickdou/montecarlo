program montecarlo
	use constants
	use math
	use material
	use distributions
	use domain
	use particle
	use simulation
	implicit none
	
	integer :: nomega, nemit, ncell, len
	real(8) :: length, side, tend, Teq, Thot, Tcold
	
	nomega = 1000
	nemit = 10000000
	ncell = 1000
	length = 1d-6
	side = 1d-6
	tend = 10d-9
	Teq = 300d0
	Thot = 1d0
	Tcold = -1d0
	
	len = 20
	
	call initboxisot(nomega, nemit, ncell, length, side, tend, Teq, Thot, Tcold)
	call simulate(len)
	call printtemp(ncell)
end program montecarlo