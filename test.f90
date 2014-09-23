program test
	implicit none
	real(8) :: array(2) = (/1d0, 1d0/)
	
	print *, 0 <= array .and. array <= 1
	
contains

end program test