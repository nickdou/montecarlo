program testomp
	use omp_lib
	implicit none
	
	integer, parameter :: iters = 10000000
	integer, parameter :: maxthreads = 8
	integer :: nthreads = 1
	
!$  print *, 'OMP enabled'
	
! 	call hello()
! 	call threadnum()
! 	call loop()
! 	call array()
	call race_reduction()
! 	call race_atomic()

contains

subroutine hello()
	!$omp parallel
	print *, 'Hello World!'
	!$omp end parallel
end subroutine hello

subroutine threadnum()
	integer :: threadi = 1
	
!$ 	nthreads = omp_get_max_threads()
	!$omp parallel private(threadi)
!$ 	threadi = omp_get_thread_num()
 	print *, 'Thread', threadi
	!$omp barrier
	!$omp single
 	print *, 'There are', nthreads, 'threads'
	!$omp end single
	!$omp end parallel
end subroutine threadnum

subroutine loop()
	integer :: i, N, threadi = 1
	
!$ 	nthreads = omp_get_max_threads()
	N = 2*nthreads
	!$omp parallel private(threadi)
!$	threadi = omp_get_thread_num()
	!$omp do
	do i = 1, N
		print *, 'Thread', threadi, 'Iteration', i
	end do
	!$omp end do
	!$omp end parallel
	print *, 'There are', nthreads, 'threads'
end subroutine loop

subroutine array()
	integer :: i, N
	integer, allocatable :: arr(:)
	
	N = iters * maxthreads
!$ 	N = iters * omp_get_max_threads()
	allocate( arr(N) )
	!$omp parallel do
	do i = 1, N
		arr(i) = i
	end do
	!$omp end parallel do
end subroutine array

subroutine race_reduction()
	integer, parameter :: bins = 3
	integer :: i, N, ind
	integer :: arr(bins) = 0
	
	N = iters * maxthreads
!$ 	N = iters * omp_get_max_threads()
	print *, '  N =', N
	!$omp parallel do private(ind) reduction(+:arr)
	do i = 1, N
		ind = mod(i, bins) + 1
		arr(ind) = arr(ind) + 1
	end do
	!$omp end parallel do
	print *, 'arr =', arr
	print *, 'sum =', sum(arr)
end subroutine race_reduction

subroutine race_atomic()
	integer, parameter :: bins = 3
	integer :: i, N, ind
	integer :: arr(bins) = 0
	
	N = iters * maxthreads
!$ 	N = iters * omp_get_max_threads()
	print *, '  N =', N
	!$omp parallel do private(ind)
	do i = 1, N
		ind = mod(i, bins) + 1
		!$omp atomic
		arr(ind) = arr(ind) + 1
	end do
	!$omp end parallel do
	print *, 'arr =', arr
	print *, 'sum =', sum(arr)
end subroutine race_atomic

end program testomp