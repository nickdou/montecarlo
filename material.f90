module material
    use constants
    use tools
    implicit none
    
    public
    
    integer :: npol, nomega
    real(8) :: Teq, ktheory
    real(8), allocatable :: omega_arr(:), domega_arr(:,:), &
        vel_arr(:,:), dos_arr(:,:), tau_arr(:,:), dedT_arr(:,:)
    
contains

real(8) pure function tauinv(omega, T, coeffs)
    real(8), intent(in) :: omega, T, coeffs(4)
    
    tauinv = coeffs(1)*omega**coeffs(2)*T**coeffs(3)*exp(-coeffs(4)/T)
end function tauinv

subroutine initmat(disp, relax, T)
    character(len=*), intent(in) :: disp, relax
    real(8), intent(in) :: T
    integer, parameter :: un1 = 2, un2 = 3
    integer :: q, p
    real(8) :: coeffspp(4), coeffsimp(4), x
    real(8), allocatable :: row(:), x_arr(:)
    
    Teq = T
    
    open(unit=un1, file=disp, action='read', status='old')
    read(un1, *) nomega, npol
    
    allocate( omega_arr(nomega) )
    allocate( domega_arr(npol, nomega) )
    allocate( vel_arr(npol, nomega) )
    allocate( dos_arr(npol, nomega) )
    allocate( tau_arr(npol, nomega) )
    allocate( dedT_arr(npol, nomega) )
    
!   call alloc(omega_arr, nomega)
!   call alloc(domega_arr, npol, nomega)
!   call alloc(vel_arr, npol, nomega)
!   call alloc(dos_arr, npol, nomega)
!   call alloc(tau_arr, npol, nomega)
!   call alloc(dedT_arr, npol, nomega)
    
    allocate( row(2*npol) )
    do q = 1,nomega
        read(un1, *) omega_arr(q), domega_arr(1,q), row
        domega_arr(2:,q) = domega_arr(1,q)
        do p = 1,npol
            vel_arr(p,q) = row(2*p-1)
            dos_arr(p,q) = row(2*p)
        end do
    end do
    deallocate( row )
    close(un1)
    
    open(unit=un2, file=relax, action='read', status='old')
    do p = 1,npol
        read(un2, *) coeffspp
        read(un2, *) coeffsimp
        do q = 1,nomega
            tau_arr(p,q) = 1/(tauinv(omega_arr(q), T, coeffspp) + &
                tauinv(omega_arr(q), T, coeffsimp))
        end do
    end do
    close(un2)
    
    allocate( x_arr(nomega) )
    x_arr = hbar/(kb*T) * omega_arr
    do q = 1,nomega
        x = x_arr(q)
        if (abs(x) < eps) then
            dedT_arr(:,q) = 0
        else
            dedT_arr(:,q) = kb*(x/(2*sinh(x/2)))**2
        end if
    end do
    deallocate( x_arr )
    
    ktheory = sum(tau_arr*vel_arr**2*dedT_arr*dos_arr*domega_arr)/3
    
end subroutine initmat

end module material