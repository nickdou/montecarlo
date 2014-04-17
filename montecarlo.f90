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
	logical :: volumetric
	
	tauconst = 10d-12
	
	nomega = 10000
	nemit = 10000000
	ncell = 100
	ntime = 0              !zero for steady
	length = 1d-6
	side = 1d-6
	tend = 100d-9
	Teq = 300d0
	Thot = 3d0
	Tcold = -3d0
	
	volumetric = .true.
	
	maxscat = 1             !zero for time limit only
	maxcoll = 100           !for single phonon simulation
	
	call settau(tauconst)
!	call initboxisot(nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	call initboxperi(nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	
	call simulate(volumetric, maxscat)
	call writetemp(ncell, ntime)
	
!	call simulateone(volumetric, maxscat, maxcoll)
	
	call printresults(Tcold - Thot, length)
		
end program montecarlo