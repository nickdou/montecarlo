program testprecision
    implicit none
    
    real(8), parameter :: zero = 0d0
    real(8), parameter :: small = 2d0**(-1022-52)
    real(8) :: large = 2d0**(1023)
    integer :: i
    
    print ('(A,ES10.3)'), 'epsilon = ', epsilon(zero)
    print ('(A,ES10.3)'), 'tiny = ', tiny(zero)
    print ('(A,I5)'), 'precision = ', precision(zero)
    print ('(A,I5)'), 'range = ', range(zero)
    print ('(A,I5)'), 'radix = ', radix(zero)
    print ('(A,I5)'), 'minexp = ', minexponent(zero)
    print ('(A,I5)'), 'maxexp = ', maxexponent(zero)
    
    do i = 1,52
        large = large + 2d0**(1023-i)
    end do
    
    print ('(Z16)'), small
    print ('(Z16)'), large
    print ('(L2)'), small == 0d0
    print ('(L2)'), 1d0 + small > 1d0
    print ('(L2)'), 1d0 + epsilon(1d0) > 1d0

end program testprecision