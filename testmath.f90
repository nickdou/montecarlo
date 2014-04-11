program testmath
	use constants
	use tools
	use math
	implicit none
	
	call initrand()
!	call printconst()
!	call testprintarray(5, 5)
	call testwritematlab(5, 5)
!	call testcumsum(5, 5)
!	call testpdftocdf()
!	call testsearchintvl()
!	call testnormtwo()
!	call testunitvec()
!	call testproject()
!	call testcrossproduct()
!	call testangdir()
contains

subroutine printconst()
	print ('(A10, ES25.17)'), 'pi = ', pi
	print ('(A10, ES25.17)'), 'hbar = ', hbar
	print ('(A10, ES25.17)'), 'kb = ', kb
end subroutine printconst

subroutine testprintarray(m, n)
	integer, intent(in) :: m, n
	character(len=6), parameter :: fmt = '(F8.0)'
	
	call printarray(makearray(m), fmt)
	print *
	call printarray(makearray(m, n), fmt)
end subroutine testprintarray

subroutine testwritematlab(m, n)
	integer, intent(in) :: m, n
	character(len=6), parameter :: fmt = '(F8.0)'
	
	call writematlab(makearray(m), fmt, 2, 'vector', 'vec')
	call writematlab(makearray(m,n), fmt, 2, 'array', 'arr')
	
end subroutine testwritematlab

subroutine testcumsum(m, n)
	integer, intent(in) :: m, n
	character(len=6), parameter :: fmt = '(F8.0)'
	real(8) :: array1(m), array2(m, n)
	integer :: i
	
	array1 = makearray(m)
	call printarray(array1, fmt, row=.true.)
	call printarray( cumsum(array1), fmt, row=.true.)
	print *
	
	array2 = makearray(m, n)
	call printarray(array2, fmt)
	call printarray( cumsum(array2, 1), fmt )
	call printarray( cumsum(array2, 2), fmt )
end subroutine testcumsum

subroutine testpdftocdf()
	integer, parameter :: m = 5
	real(8) :: cdf(m), pdf(m)
	
	pdf = makearray(m)
!	pdf = (/0.2, 0.2, 0.2, 0.2, 0.2/)
	cdf = pdftocdf(pdf)
	
	call printarray(pdf, '(F8.0)')
	print *
	call printarray(cdf, '(F8.3)')
end subroutine testpdftocdf

subroutine testsearchintvl()
	integer, parameter :: m = 100
	integer, parameter :: N = 10
	real(8) :: cdf(m), r
	integer :: i
	
	cdf = (/(dble(i)/m, i=1,m)/)
	call printarray(cdf, '(F8.3)', row=.true.)
	
	do i = 1, N
		call random_number(r)
		print ('(F8.3,2I6)'), r, ceiling(r*m), searchintvl(cdf, r)
	end do
	
	do i = 0, m-1
		r = dble(i)/m
		if (searchintvl(cdf, r) /= nint(r*m+1)) then
			print ('(I6)'), i
		end if
	end do
end subroutine testsearchintvl

subroutine testnormtwo()
	real(8) :: vec(2)
	
	vec =  (/1, 1/)
	call printarray(vec, '(F8.3)', row=.true.)
	print ('(F8.3)'), normtwo(vec)
end subroutine testnormtwo

subroutine testunitvec()
	real(8) :: vec(2)
	
	vec = (/3, 4/)
!	vec = (/0, 0/)
	call printarray(vec, '(F8.3)', row=.true.)
	call printarray(unitvec(vec), '(F8.3)', row=.true.)
end subroutine testunitvec

subroutine testproject()
	real(8) :: a(3), b(3)
	
!	call random_number(a)
	a = (/1, 0, 0/)
	b = (/0, 0, 0/)
	call printarray(a, '(F8.3)', row=.true.)
	call printarray(b, '(F8.3)', row=.true.)
	call printarray(project(a, b), '(F8.3)', row=.true.)
end subroutine testproject

subroutine testcrossproduct()
	real(8) :: a(3), b(3)
	
!	call random_number(a)
	a = (/-1, 0, 0/)
	b = (/0, 1, 0/)
	call printarray(a, '(F8.3)', row=.true.)
	call printarray(b, '(F8.3)', row=.true.)
	call printarray(cross_product(a, b), '(F8.3)', row=.true.)
end subroutine testcrossproduct

subroutine testangdir()
	real(8) :: ang(2), ang2(2), dir(3)
	
	call random_number(ang)
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