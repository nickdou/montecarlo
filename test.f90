program test
	use tools
	implicit none
	
	call makeunit(0.8d0, 0.2d0, 8d0, 0.04d0)

contains

subroutine makeunit(a, b, c, d)
	real(8), intent(in) :: a, b, c, d
	real(8) :: l, data_arr(9,24)
	
	l = c/4*sqrt(2d0)
	
	data_arr = reshape((/ &
	0d0,   -l,    -d,    0d0,  0d0,       d,      b,   0d0,   0d0, &
	  b,   -l,    -d,    0d0,  0d0, 2*a+2*d,      d,   0d0,   0d0, &
	0d0,   -l,   2*a,    0d0,  0d0,       d,      b,   0d0,   0d0, &
	0d0,   -l,    -d,    b+d,  0d0,     0d0,    0d0,   l-a,   0d0, &
	b+d,   -l,    -d,    0d0,  0d0,   a+2*d,    0d0, l-a-d,   0d0, &
	b+d,   -l,   a+d,    0d0,  0d0,       a,    0d0,     l,   0d0, &
	b+d,   -l, 2*a+d,   -b-d,  0d0,     0d0,    0d0,     l,   0d0, &
	0d0,   -l, 2*a+d,    0d0,  0d0,      -d,    0d0,     l,   0d0, &
	0d0,   -l,   2*a,      b,  0d0,     0d0,    0d0,     l,   0d0, &
	  b,   -l,   2*a,    0d0,  0d0,      -a,    0d0,     l,   0d0, &
	  b,   -l,     a,    0d0,  0d0,      -a,    0d0,   l-a,   0d0, &
	  b,   -l,   0d0,     -b,  0d0,     0d0,    0d0,   l-a,   0d0, &
	0d0,   -l,   0d0,    0d0,  0d0,      -d,    0d0,   l-a,   0d0, &
	0d0,   -a,   0d0,    0d0,  0d0,      -d,      b,   0d0,   0d0, &
	  l,   -a,    -d, -l+b+d,  0d0,     0d0,    0d0,    -d,   0d0, &
	  l, -a-d,    -d, -l+b+d,  0d0,     0d0,    0d0,   0d0, a+2*d, &
	  l, -a-d,    -d,    0d0,  0d0,     a+d,    0d0,     d,   0d0, &
	  l,   -a,    -d,    0d0,  0d0,     a+d,   -l+b,   0d0,   0d0, &
	  l,   -a,     a,    0d0,    a,     0d0,   -l+b,   0d0,   0d0, &
	  l,  0d0,     a,    0d0, -a-d,     0d0,    0d0,   0d0,     d, &
	  l,  0d0,   a+d,    0d0, -a-d,     0d0, -l+b+d,   0d0,   0d0, &
	  l,  0d0,     a,    0d0,  0d0,       d,   -l+b,   0d0,   0d0, &
	  b,  0d0,   a+d,      d,  0d0,     0d0,    0d0,   0d0,     a, &
	  b,  0d0,   2*a,    0d0,  0d0,       d,     -b,   0d0,   0d0 /), &
	(/9, 24/))
	
	call printarray(transpose(data_arr), '(F8.3)')
	
end subroutine makeunit

end program test