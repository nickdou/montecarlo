program test
	implicit none
	integer, parameter :: N = 5
	real(8), allocatable :: array(:)
	integer :: i
	
	allocate( array(N) )
	array = (/(dble(i), i=1,N)/)
	print *, array
	
contains

end program test