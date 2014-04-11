program test
	implicit none
	
	call testdo(3)
	
contains

subroutine testdo(n)
	integer, intent(in) :: n
	integer :: i
	
	do i = 1,100
		if (i == n) exit
	end do
	print *, i
end subroutine testdo

end program test