module tools
	implicit none
	
	public
	
	interface makearray
		procedure makearray1, makearray2
	end interface makearray
	
	interface printarray
		module procedure printarray1, printarray2
	end interface printarray
	
	interface writematlab
		module procedure writematlab1, writematlab2
	end interface writematlab

contains

function makearray1(m) result(array)
	integer, intent(in) :: m
	real(8) :: array(m)
	integer :: i
	
	array = (/(i, i=1,m)/)
end function makearray1

function makearray2(m, n) result(array)
	integer, intent(in) :: m, n
	real(8) :: array(m, n)
	
	array = reshape(makearray1(m*n), (/m, n/))
end function makearray2

subroutine printarray2(array, fmt, un)
	real(8), intent(in) :: array(:,:)
	character(len=*), intent(in) :: fmt
	integer, intent(in), optional :: un
	integer :: unit, i, j
	
	if (present(un)) then
		unit = un
	else
		unit = 6
	end if
	
	do i = 1,size(array, 1)
		do j = 1,size(array, 2)
			write(unit, fmt, advance='no') array(i,j)
		end do
		write(unit, *)
	end do
end subroutine printarray2

subroutine printarray1(vector, fmt, un, row)
	real(8), intent(in) :: vector(:)
	character(len=*), intent(in) :: fmt
	integer, intent(in), optional :: un
	logical, intent(in), optional :: row
	integer :: unit
	logical :: transpose = .false.
	
	if (present(un)) then
		unit = un
	else
		unit = 6
	end if
	
	if (present(row)) then
		transpose = row
	end if
	
	if (transpose) then
		call printarray2( reshape(vector, (/1, size(vector, 1)/)), fmt, un )
	else
		write(unit, fmt) vector
	end if
end subroutine printarray1

subroutine writematlab2(array, fmt, unit, filename, var)
	real(8), intent(in) :: array(:,:)
	character(len=*), intent(in) :: fmt, filename
	integer, intent(in) :: unit
	character(len=*), intent(in), optional :: var
	
	open(unit, file=filename//'.m', action='write', status='replace')
	
	if (present(var)) then
		write(unit,'(A, A)') var, ' = ['
	else
		write(unit,'(A)') 'A = ['
	end if
	
	call printarray2(array, fmt, unit)
	write(unit,'(A)') '];'
	
end subroutine writematlab2

subroutine writematlab1(vector, fmt, unit, filename, var)
	real(8), intent(in) :: vector(:)
	character(len=*), intent(in) :: fmt, filename
	integer, intent(in) :: unit
	character(len=*), intent(in), optional :: var
	
	open(unit, file=filename//'.m', action='write', status='replace')
	
	if (present(var)) then
		write(unit,'(A, A)') var, ' = ['
	else
		write(unit,'(A)') 'x = ['
	end if
	
	call printarray1(vector, fmt, unit)
	write(unit,'(A)') '];'
	
end subroutine writematlab1

subroutine showprogress(i, iend, num)
	integer, intent(in) :: i, iend
	integer, intent(in), optional :: num
	real(8) :: x, unit
	integer :: len, bars
	
	if (present(num)) then
		len = num
	else
		len = 20
	end if
		
	unit = dble(len)/iend
	x = i*unit
	bars = floor(x)
	if ((x - bars)/unit < 0.999999) then
		print *, nint(100*dble(i)/iend), '%  ', &
			'[', repeat('|',bars), repeat('-',len-bars), ']'
	end if
	
end subroutine showprogress

end module tools