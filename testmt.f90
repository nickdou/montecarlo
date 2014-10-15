program testmt
	use omp_lib
	use mt_stream
	implicit none

	logical, parameter :: USE_MT = .true.
	integer, parameter :: seed = 2543087
	integer :: nthreads = 1
	type(mt_state) :: mtsparent

	integer, parameter :: niters = 100
	integer, parameter :: nx = 100000

!$ 	nthreads = omp_get_max_threads()
!$ 	print ('(A,I2,A,/)'), 'OMP enabled, ', nthreads, ' threads available'

	if (USE_MT) then
		print *, 'Using Mersenne Twister PRNG'
		call initmt()
	else
		print *, 'Using intrinsic PRNG'
		call random_seed()
	end if
	
	print *, 'niters = ', niters
	print *, '    nx = ', nx
	
	call loop(niters, nx)

contains

subroutine initmt()
	call set_mt19937
	call new(mtsparent)
	call init(mtsparent, seed)
end subroutine initmt

subroutine loop(niters, nx)
	integer, intent(in) :: niters, nx
	type(mt_state) :: mts
	integer :: i, mtsid
	real(8) :: mean(niters), mu, sigma 

	!$omp parallel private(mtsid, mts)
	if (USE_MT) then
		mtsid = 0
!$ 		mtsid = omp_get_thread_num()
		!$omp critical
		call create_stream(mtsparent, mts, mtsid)
		!$omp end critical
	end if
	!$omp do
	do i = 1,niters
		mean(i) = mc(mts, nx)
	end do
	!$omp end do
	!$omp end parallel

	mu = sum(mean)/niters
	sigma = sqrt(sum((mean - mu)**2))/niters

	print *, '    mu = ', mu
	print *, ' sigma = ', sigma
end subroutine loop

real(8) function mc(mts, nx) result(mean)
	type(mt_state), intent(inout) :: mts
	integer, intent(in) :: nx
	integer :: i
	real(8) :: x, cum

	cum = 0d0
	if (USE_MT) then
		do i = 1,nx
			x = genrand_double1(mts)
			cum = cum + x
		end do
	else
		do i = 1,nx
			call random_number(x)
			cum = cum + x
		end do
	end if
	mean = cum/nx
end function mc

end program testmt
