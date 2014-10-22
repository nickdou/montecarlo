module tools
	use omp_lib
	use mt_stream
	implicit none
	
	public
	
! 	interface zeros
! 		module procedure zeros1, zeros2
! 	end interface zeros
	
! 	interface alloc
! 		module procedure alloc1_real, alloc1_int, alloc2_real, alloc2_int
! 	end interface alloc

	interface randnum
		module procedure randnum_dbl, randnum_arr
	end interface randnum
	
	interface printarray
		module procedure printarray1, printarray2
	end interface printarray
	
	interface writematlab
		module procedure writematlab1, writematlab2
	end interface writematlab
	
	integer :: nthreads = 1
	logical :: usemt
	type(mt_state), allocatable :: randstate(:)
	
contains

subroutine initomp(n)
	integer, intent(in), optional :: n
	
!$ 	nthreads = omp_get_max_threads()
!$ 	if (present(n)) then
!$ 		nthreads = min(nthreads, n)
!$ 	end if
!$ 	print ('(A,I2,A,/)'), 'OMP enabled, ', nthreads, ' threads available'
end subroutine

subroutine initrand(mt)
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html

	logical, intent(in) :: mt
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid, t(2), s
	integer(8) :: count, tms
	type(mt_state) :: parentstate
	
	usemt = mt
	if (mt) then
		n = 1
	else
		call random_seed(size = n)
	end if
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(unit=un, file="/dev/urandom", access="stream", &
		form="unformatted", action="read", status="old", iostat=istat)
	
	if (istat == 0) then
		print ('(A)'), 'Random seed from urandom'
		read(un) seed
		close(un)
	   
	else
		print ('(A)'), 'Random seed from time and pid'
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
	print *, seed
	if (usemt) then
		call set_mt19937
		call new(parentstate)
		call init(parentstate, seed)
		
		allocate( randstate(0:nthreads-1) )
		do i = 0, nthreads-1
			call create_stream(parentstate, randstate(i), i)
		end do
	else
		call random_seed(put=seed)
	end if
end subroutine initrand

subroutine randnum_dbl(x)
	real(8), intent(out) :: x
	integer :: threadi = 0
	
!$  threadi = omp_get_thread_num()
	if (usemt) then
		x = genrand_double1(randstate(threadi))
	else
		call random_number(x)
	end if
end subroutine randnum_dbl

subroutine randnum_arr(x)
	real(8), intent(out) :: x(:)
	integer :: i, threadi = 0
	
!$  threadi = omp_get_thread_num()
	if (usemt) then
		do i = 1,size(x)
			x(i) = genrand_double1(randstate(threadi))
		end do
	else
		call random_number(x)
	end if
end subroutine randnum_arr

! function zeros1(m) result(array)
! 	integer, intent(in) :: m
! 	real(8) :: array(m)
! 	integer :: i
!
! 	array = (/(0d0, i=1,m)/)
! end function zeros1
!
! function zeros2(m, n) result(array)
! 	integer, intent(in) :: m, n
! 	real(8) :: array(m, n)
!
! 	array = reshape(zeros1(m*n), (/m, n/))
! end function zeros2

subroutine printarray2(array, fmt, un)
	real(8), intent(in) :: array(:,:)
	character(len=*), intent(in) :: fmt
	integer, intent(in), optional :: un
	integer :: i, j, unit = 6
	
	if (present(un)) then
		unit = un
	end if
	
	do i = 1, size(array,1)
		do j = 1, size(array,2)
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
	integer :: unit = 6
	logical :: transpose = .false.
	
	if (present(un)) then
		unit = un
	end if
	
	if (present(row)) then
		transpose = row
	end if
	
	if (transpose) then
		call printarray2( reshape(vector, (/1, size(vector,1)/)), fmt, un )
	else
		write(unit, fmt) vector
	end if
end subroutine printarray1

! pure subroutine alloc1_real(arr, n, zero)
! 	real(8), intent(inout), allocatable :: arr(:)
! 	integer, intent(in) :: n
! 	logical, intent(in), optional :: zero
!
! 	if (allocated(arr)) then
! 		deallocate( arr )
! 	end if
!
! 	if (present(zero)) then
! 		if (zero) then
! 			allocate( arr(0:n) )
! 		else
! 			allocate( arr(n) )
! 		end if
! 	else
! 		allocate( arr(n) )
! 	end if
!
! 	arr = 0
! end subroutine alloc1_real
!
! pure subroutine alloc1_int(arr, n, zero)
! 	integer, intent(inout), allocatable :: arr(:)
! 	integer, intent(in) :: n
! 	logical, intent(in), optional :: zero
!
! 	if (allocated(arr)) then
! 		deallocate( arr )
! 	end if
!
! 	if (present(zero)) then
! 		if (zero) then
! 			allocate( arr(0:n) )
! 		else
! 			allocate( arr(n) )
! 		end if
! 	else
! 		allocate( arr(n) )
! 	end if
!
! 	arr = 0
! end subroutine alloc1_int
!
! pure subroutine alloc2_real(arr, m, n, zero)
! 	real(8), intent(inout), allocatable :: arr(:,:)
! 	integer, intent(in) :: m, n
! 	logical, intent(in), optional :: zero(2)
!
! 	if (allocated(arr)) then
! 		deallocate( arr )
! 	end if
!
! 	if (present(zero)) then
! 		if (all(zero)) then
! 			allocate( arr(0:m, 0:n) )
! 		else if (zero(1)) then
! 			allocate( arr(0:m, n) )
! 		else if (zero(2)) then
! 			allocate( arr(m, 0:n) )
! 		else
! 			allocate( arr(m,n) )
! 		end if
! 	else
! 		allocate( arr(m,n) )
! 	end if
!
! 	arr = 0
! end subroutine alloc2_real
!
! pure subroutine alloc2_int(arr, m, n, zero)
! 	integer, intent(inout), allocatable :: arr(:,:)
! 	integer, intent(in) :: m, n
! 	logical, intent(in), optional :: zero(2)
!
! 	if (allocated(arr)) then
! 		deallocate( arr )
! 	end if
!
! 	if (present(zero)) then
! 		if (all(zero)) then
! 			allocate( arr(0:m, 0:n) )
! 		else if (zero(1)) then
! 			allocate( arr(0:m, n) )
! 		else if (zero(2)) then
! 			allocate( arr(m, 0:n) )
! 		else
! 			allocate( arr(m,n) )
! 		end if
! 	else
! 		allocate( arr(m,n) )
! 	end if
!
! 	arr = 0
! end subroutine alloc2_int

subroutine writematlab2(array, fmt, unit, filename, var)
	real(8), intent(in) :: array(:,:)
	character(len=*), intent(in) :: fmt, filename, var
	integer, intent(in) :: unit
	
	open(unit, file=filename//'.m', action='write', status='replace')
	
	write(unit,'(A, A)') var, ' = ['
	call printarray2(array, fmt, unit)
	write(unit,'(A)') '];'
	
end subroutine writematlab2

subroutine writematlab1(vector, fmt, unit, filename, var)
	real(8), intent(in) :: vector(:)
	character(len=*), intent(in) :: fmt, filename, var
	integer, intent(in) :: unit
	
	open(unit, file=filename//'.m', action='write', status='replace')
	
	write(unit,'(A, A)') var, ' = ['
	call printarray1(vector, fmt, unit)
	write(unit,'(A)') '];'
	
end subroutine writematlab1

subroutine showprogress(i, itot, num)
	integer, intent(in) :: i, itot
	integer, intent(in), optional :: num
	real(8) :: xstep, x
	integer :: nbars = 20, bars
	
	if (present(num)) then
		nbars = num
	end if
	
	xstep = dble(nbars)/itot
	x = i*xstep
	bars = floor(x)
	
	if ((x - floor(x))/xstep < 1d0 - 1d-6) then
		print *, nint(100*dble(i)/itot), '%  ', &
			'[', repeat('|',bars), repeat('-',nbars-bars), ']'
	end if
	
end subroutine showprogress

end module tools