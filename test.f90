program test
	implicit none
	
	integer, allocatable :: array(:)
	integer :: n
	
	call initarray((/ 0, 1, 2, 3, 4 /))
	do n = 1,10
		call decrement()
		print *, array
	end do
	
contains

subroutine initarray(arr)
	integer, intent(in) :: arr(:)
	
	allocate( array(size(arr, 1)) )
	array = arr
end subroutine initarray

subroutine decrement()
	integer :: ind
	
	if (sum(array) > 0) then
		ind = minloc(array > 0, 1)
		array(ind) = array(ind) - 1
	end if
end subroutine decrement

end program test