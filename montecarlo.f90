program montecarlo
	use material
	use simulation
	implicit none
	
	logical :: one, vol, gf
	integer :: nomega, nemit, ncell, ntime, maxscat, maxcoll
	real(8) :: tauconst, length, side, tend, Teq, Thot, Tcold
	
	one = .false.
	vol = .true.
	gf = .true.
	
	tauconst = 10d-12
	
	nomega = 1000
	nemit = 1000000
	ncell = 100
	ntime = 0              !zero for steady
	length = 10000d-9
	side = 10000d-9
	tend = 100d-9
	Teq = 300d0
	Thot = 3d0
	Tcold = -3d0
	
	maxscat = 1             !zero for time limit only
	maxcoll = 100           !for single phonon simulation
	
	call settau(tauconst)
!	call initisot1d(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
!	call initbulk1d(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
!	call initfilm(vol, gf, nomega, nemit, ncell, ntime, length, side, tend, Teq, Thot, Tcold)
	
	if (one) then
		call simulateone(maxscat, maxcoll)
		
	else
		call simulate(maxscat)
		
		call writetemp(ncell, ntime)
		if (gf) then
			call writeflux(ncell)
		else
			call printcond(Tcold - Thot, length)
		end if
	end if		
	
end program montecarlo