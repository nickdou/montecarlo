module tools
	implicit none
	
	public
	
	interface makearray
		procedure makearray1, makearray2
	end interface makearray
	
	interface printarray
		module procedure printarray1, printarray2
	end interface printarray
	
	interface alloc
		module procedure alloc1_real, alloc1_int, alloc2_real, alloc2_int
	end interface alloc
	
	interface writematlab
		module procedure writematlab1, writematlab2
	end interface writematlab

contains

subroutine initrand()
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
	
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid, t(2), s
	integer(8) :: count, tms

	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(newunit=un, file="/dev/urandom", access="stream", &
		 form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
	   read(un) seed
	   close(un)
	else
	   ! Fallback to XOR:ing the current time and pid. The PID is
	   ! useful in case one launches multiple instances of the same
	   ! program in parallel.
	   call system_clock(count)
	   if (count /= 0) then
		  t = transfer(count, t)
	   else
		  call date_and_time(values=dt)
		  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
			   + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
			   + dt(3) * 24 * 60 * 60 * 60 * 1000 &
			   + dt(5) * 60 * 60 * 1000 &
			   + dt(6) * 60 * 1000 + dt(7) * 1000 &
			   + dt(8)
		  t = transfer(tms, t)
	   end if
	   s = ieor(t(1), t(2))
	   pid = getpid() + 1099279 ! Add a prime
	   s = ieor(s, pid)
	   if (n >= 3) then
		  seed(1) = t(1) + 36269
		  seed(2) = t(2) + 72551
		  seed(3) = pid
		  if (n > 3) then
			 seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
		  end if
	   else
		  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
	   end if
	end if
	call random_seed(put=seed)
end subroutine initrand

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

pure subroutine alloc1_real(arr, n, zero)
	real(8), intent(inout), allocatable :: arr(:)
	integer, intent(in) :: n
	logical, intent(in), optional :: zero
	
	if (allocated(arr)) then
		deallocate( arr )
	end if
	
	if (present(zero)) then
		if (zero) then
			allocate( arr(0:n) )
		else
			allocate( arr(n) )
		end if
	else
		allocate( arr(n) )
	end if
	
	arr = 0
end subroutine alloc1_real

pure subroutine alloc1_int(arr, n, zero)
	integer, intent(inout), allocatable :: arr(:)
	integer, intent(in) :: n
	logical, intent(in), optional :: zero
	
	if (allocated(arr)) then
		deallocate( arr )
	end if
	
	if (present(zero)) then
		if (zero) then
			allocate( arr(0:n) )
		else
			allocate( arr(n) )
		end if
	else
		allocate( arr(n) )
	end if
	
	arr = 0
end subroutine alloc1_int

pure subroutine alloc2_real(arr, m, n, zero)
	real(8), intent(inout), allocatable :: arr(:,:)
	integer, intent(in) :: m, n
	logical, intent(in), optional :: zero(2)
	
	if (allocated(arr)) then
		deallocate( arr )
	end if
	
	if (present(zero)) then
		if (all(zero)) then
			allocate( arr(0:m, 0:n) )
		else if (zero(1)) then
			allocate( arr(0:m, n) )
		else if (zero(2)) then
			allocate( arr(m, 0:n) )
		else
			allocate( arr(m,n) )
		end if
	else
		allocate( arr(m,n) )
	end if
	
	arr = 0
end subroutine alloc2_real

pure subroutine alloc2_int(arr, m, n, zero)
	integer, intent(inout), allocatable :: arr(:,:)
	integer, intent(in) :: m, n
	logical, intent(in), optional :: zero(2)
	
	if (allocated(arr)) then
		deallocate( arr )
	end if
	
	if (present(zero)) then
		if (all(zero)) then
			allocate( arr(0:m, 0:n) )
		else if (zero(1)) then
			allocate( arr(0:m, n) )
		else if (zero(2)) then
			allocate( arr(m, 0:n) )
		else
			allocate( arr(m,n) )
		end if
	else
		allocate( arr(m,n) )
	end if
	
	arr = 0
end subroutine alloc2_int

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