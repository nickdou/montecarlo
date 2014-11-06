program testmath
	use constants
	use tools
	use math
	implicit none
	
!     call initomp()
	call initrand(.true.)
!     call printconst()
!     call testrand(100, 100000)
!     call testprintarray(5, 5)
!     call testwritematlab(5, 5)
!     call testcumsum(5, 5)
!     call testpdftocdf()
    call testsearchbin()
!     call testnormtwo()
!     call testunitvec()
!     call testproject()
!     call testcrossproduct()
!     call testinverse2()
!     call testangdir()
contains

subroutine printconst()
	print ('(A10, ES25.17)'), 'pi = ', pi
	print ('(A10, ES25.17)'), 'hbar = ', hbar
	print ('(A10, ES25.17)'), 'kb = ', kb
end subroutine printconst

subroutine testrand(niters, nx)
	integer, intent(in) :: niters, nx
	integer :: i
	real(8) :: mean(niters), mu, sigma
	
	!$omp parallel do
	do i = 1, niters
		mean(i) = mc(nx)
	end do
	!$omp end parallel do
	
	mu = sum(mean)/niters
	sigma = sqrt(sum((mean - mu)**2))/niters
	
	print *, '    mu = ', mu
	print *, ' sigma = ', sigma
end subroutine testrand

real(8) function mc(nx) result(mean)
	integer, intent(in) :: nx
	integer :: i
	real(8) :: x_arr(nx), cum

! 	cum = 0d0
! 	do i = 1,nx
! 		call randnum(x)
! 		cum = cum + x
! 	end do
	
	call randnum(x_arr)
	cum = sum(x_arr)
	mean = cum/nx
end function mc

subroutine testprintarray(m, n)
	integer, intent(in) :: m, n
	real(8) :: array1(m), array2(m,n)
	character(len=6), parameter :: fmt = '(F8.0)'
	integer :: i
	
	array1 = (/(dble(i), i=1,m)/)
	array2 = reshape((/(dble(i), i=1,m*n)/), (/m,n/))

	call printarray(array1, fmt)
	print *
	call printarray(array2, fmt)
end subroutine testprintarray

subroutine testwritematlab(m, n)
	integer, intent(in) :: m, n
	real(8) :: array1(m), array2(m,n)
	character(len=6), parameter :: fmt = '(F8.0)'
	integer :: i
	
	array1 = (/(dble(i), i=1,m)/)
	array2 = reshape((/(dble(i), i=1,m*n)/), (/m,n/))
	
	call writematlab(array1, fmt, 2, 'vector', 'vec')
	call writematlab(array2, fmt, 2, 'array', 'arr')
end subroutine testwritematlab

subroutine testcumsum(m, n)
	integer, intent(in) :: m, n
	real(8) :: array1(m)
	character(len=6), parameter :: fmt = '(F8.0)'
	integer :: i
	
	array1 = (/(dble(i), i=1,m)/)
	
	call printarray(array1, fmt, row=.true.)
	call printarray( cumsum(array1), fmt, row=.true.)
end subroutine testcumsum

subroutine testpdftocdf()
	integer, parameter :: m = 5
	real(8) :: pdf(m), cdf(m)
	integer :: i
	
	pdf = (/(0d0, i=1,m)/)
	cdf = pdftocdf(pdf)
	
	call printarray(pdf, '(F8.3)')
	print *
	call printarray(cdf, '(F8.3)')
end subroutine testpdftocdf

subroutine testsearchbin()
	integer, parameter :: m = 20
	integer, parameter :: N = 10
	real(8) :: cdf(m), r
	integer :: i
	
	cdf = (/(dble(i)/m, i=1,m)/)
    call printarray(cdf, '(F8.3)', row=.true.)
    
!     print *, searchbin(cdf, dble(m))
!     print *, searchbin(cdf, dble(m+1))
!     print *, searchbin(cdf, -1d0)
    
    do i = 1, N
        call randnum(r)
        print ('(F8.3,2I6)'), r, ceiling(r*m), searchbin(cdf, r)
    end do

    do i = 0, m-1
        r = dble(i)/m
        if (searchbin(cdf, r) /= nint(r*m+1)) then
            print ('(I6)'), i
        end if
    end do
end subroutine testsearchbin

subroutine testnormtwo()
	real(8) :: vec(2)
	
	vec =  (/1, 1/)
	call printarray(vec, '(F8.3)', row=.true.)
	print ('(F8.3)'), normtwo(vec)
end subroutine testnormtwo

subroutine testunitvec()
	real(8) :: vec(2)
	
	vec = (/3, 4/)
! 	vec = (/0, 0/)
	call printarray(vec, '(F8.3)', row=.true.)
	call printarray(unitvec(vec), '(F8.3)', row=.true.)
end subroutine testunitvec

subroutine testproject()
	real(8) :: a(3), b(3)
	
!	call randnum(a)
	a = (/1, 0, 0/)
	b = (/0, 0, 0/)
	call printarray(a, '(F8.3)', row=.true.)
	call printarray(b, '(F8.3)', row=.true.)
	call printarray(project(a, b), '(F8.3)', row=.true.)
end subroutine testproject

subroutine testcrossproduct()
	real(8) :: a(3), b(3)
	
!	call randnum(a)
	a = (/-1, 0, 0/)
	b = (/0, 1, 0/)
	call printarray(a, '(F8.3)', row=.true.)
	call printarray(b, '(F8.3)', row=.true.)
	call printarray(cross_product(a, b), '(F8.3)', row=.true.)
end subroutine testcrossproduct

subroutine testangdir()
	real(8) :: ang(2), ang2(2), dir(3)
	
	call randnum(ang)
	ang = 1 - 2*ang
!	ang = (/0, 0/)
	dir = angtodir(ang)
	ang2 = dirtoang(dir)
	
	call printarray(ang,  '(F8.3)', row=.true.)
	call printarray(dir,  '(F8.3)', row=.true.)
	print ('(F8.3)'), normtwo(dir)
	call printarray(ang2, '(F8.3)', row=.true.)
end subroutine testangdir

end program testmath