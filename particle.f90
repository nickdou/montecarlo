module particle
    use material
    use distributions
    use domain
    implicit none
    
    public
    
    type phonon
        private
        logical :: sign, alive
        real(8) :: x0(3), x(3), dir(3), tscat
        integer :: p, q, ind
    end type phonon
contains

logical pure function isalive(phn) result(alive)
    type(phonon), intent(in) :: phn
    
    alive = phn%alive
end function isalive

subroutine remove(phn)
    type(phonon), intent(inout) :: phn
    
    call addcumdisp(phn%sign, phn%x - phn%x0)
    phn%alive = .false.
end subroutine remove

pure function getpos(phn) result(x)
    type(phonon), intent(in) :: phn
    real(8) :: x(3)
    
    x = phn%x
end function getpos

type(phonon) function emit() result(phn)
    call drawfluxprop(phn%p, phn%q)
    call drawemitstate(phn%ind, phn%x, phn%dir, phn%sign)
    call drawscattime(phn%tscat, phn%p, phn%q)
    phn%x0 = phn%x
    phn%alive = .true.
    
!     print ('(A1,3X,I3,7ES12.4)'), 'e', phn%ind, phn%x, phn%dir, phn%tscat
end function emit

subroutine scatter(phn, nscat)
    type(phonon), intent(inout) :: phn
    integer, intent(inout) :: nscat
    
!     if (phn%tscat < eps) then
    if (phn%ind == 0) then
        call drawscatterprop(phn%p, phn%q)
        call drawangiso(phn%dir)
        call drawscattime(phn%tscat, phn%p, phn%q)
        
!         phn%ind = 0
        nscat = nscat + 1
        
!         print ('(A1,3X,I3,7ES12.4)'), 's', phn%ind, phn%x, phn%dir, phn%tscat
    end if
end subroutine scatter

subroutine advect(phn, t)
    type(phonon), intent(inout) :: phn
    real(8), intent(inout) :: t
    real(8) :: x(3), dir(3), v, deltat
    integer :: ind, bc
    
    x = phn%x
    dir = phn%dir
    ind = phn%ind
    v = vel_arr(phn%p, phn%q)
    deltat = phn%tscat
    
    call updatestate(bc, ind, x, dir, v, deltat, t)
    call appendtraj(x)
    
    call recordtime(phn%sign, deltat, phn%x, x)
    call recorddisp(phn%sign, phn%x, x)
    call recordloc(phn%sign, t, deltat, phn%x, x)
    
!     print ('(A1,3X,I3,6ES12.4)'), 'a', ind, x, dir
    
    if (bc == PERI_BC) then
        call addcumdisp(phn%sign, ind)
        call applybc(bc, ind, x, dir)
        call appendtraj(x)
    else
        call applybc(bc, ind, x, dir)
    end if
    
    phn%x = x
    phn%dir = dir
    phn%ind = ind
    phn%tscat = max(0d0, phn%tscat - deltat)
    
!     print ('(A1,2I3,7ES12.4)'), 'b', bc, ind, x, dir, phn%tscat
    
    if (all(abs(dir) < eps, 1)) then
        call remove(phn)
    end if
end subroutine advect

end module particle
