program montecarlo
	use constants
	use math
	use material
	use distributions
	use domain
	use particle
	use simulation
	implicit none
	
	integer :: nomega, nemit, ncell, ntime, maxscat, maxcoll
	real(8) :: tauconst, length, side, tend, Teq, Thot, Tcold
	
	tauconst = 10d-12
	
	nomega = 1000
	nemit = 100000
	ncell = 100
	ntime = 0              !zero for steady
	length = 1d-8
	side = 1d-8
	tend = 100d-9
	Teq = 300d0
	Thot = 1d0
	Tcold = -1d0
	
	maxscat = 1             !zero for time limit only
	maxcoll = 100
	
	call settau(tauconst)
	call initboxperi(nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	
	call simulate(maxscat)
	call writetemp(ncell, ntime)

!	call simulateone(maxscat, maxcoll)
	
end program montecarlo