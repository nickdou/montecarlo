program testomp
	use omp_lib
	implicit none
	
! 	call hello()
! 	call threads()
! 	call loop()
	call array()
contains

subroutine hello()
	!$omp parallel
	print *, 'Hello World!'
	!$omp end parallel
end subroutine hello

subroutine threads()
	!$omp parallel
	print *, 'Thread', omp_get_thread_num()
	!$omp barrier
	!$omp single
	print *, 'There are', omp_get_num_threads(), 'threads'
	!$omp end single
	!$omp end parallel
end subroutine threads

subroutine loop()
	integer :: i
	!$omp parallel do
	do i = 1, omp_get_max_threads() !OMP_GET_NUM_THREADS doesn't work here
		print *, 'Thread', omp_get_thread_num(), 'Iteration', i
	end do
	!$omp end parallel do
end subroutine loop

subroutine array()
	integer :: i, N
	integer, allocatable :: arr(:)
	
	N = 100000000 * 12
!$ 	N = 100000000 * omp_get_max_threads() !conditionally compiled
!$  print *, 'OMP enabled'
	allocate( arr(N) )
	!$omp parallel do
	do i = 1, N
		arr(i) = i
	end do
	!$omp end parallel do
end subroutine array

end program testomp