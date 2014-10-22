module distributions
    use constants
    use tools
    use math
    use material
    implicit none
    
    public
    
!   private 
!   public  :: getnomega, getteq, dedT, calculatek, energypdf, fluxpdf, scatterpdf, &
!       initomega, initpropcdf, getpseudoenergy, getpseudoflux, &
!       drawenergyprop, drawfluxprop, drawscatterprop, &
!       drawposrect, drawangiso, drawanghalf, drawscattime
!   
!   interface initomega
!       module procedure initomega_int, initomega_real
!   end interface initomega
    
!   integer :: nomega
!   real(8) :: domega, Teq
!   real(8), allocatable :: omega_arr(:)
    real(8), allocatable :: energypdf_arr(:,:), fluxpdf_arr(:,:), scatterpdf_arr(:,:)
    real(8), allocatable :: energycdf_arr(:,:), fluxcdf_arr(:,:), scattercdf_arr(:,:)
    
contains

!real(8) pure function getnomega() result(n)
!   n = nomega
!end function getnomega
!
!real(8) pure function getteq() result(T)
!   T = Teq
!end function getteq

!real(8) pure function dedT(omega, T)
!   real(8), intent(in) :: omega, T
!   real(8) :: x
!   
!   x = hbar*omega/(kb*T)
!   if (x == 0) then
!       dedT = 0
!   else
!       dedT = kb*(x/(2*sinh(x/2)))**2
!   end if
!end function dedT
!
!real(8) pure function calculatek(T) result(k)
!   real(8), intent(in) :: T
!   real(8) :: integrand(npol*nomega)
!   integer :: p, i
!   
!   integrand = (/ ((tau(omega_arr(i), p, T)*vel(omega_arr(i), p)**2*&
!       dedT(omega_arr(i), T)*dos(omega_arr(i), p), p=1,npol), i=1,nomega) /)
!   k = domega/3 * sum(integrand)
!end function calculatek

!real(8) pure function energypdf(omega, p, T) result(pdf)
!   real(8), intent(in) :: omega, T
!   integer, intent(in) :: p
!   
!   pdf = dos(omega, p)*dedT(omega, T)
!end function energypdf
!
!real(8) pure function fluxpdf(omega, p, T) result(pdf)
!   real(8), intent(in) :: omega, T
!   integer, intent(in) :: p
!   
!   pdf = dos(omega, p)*vel(omega, p)*dedT(omega, T)
!end function fluxpdf
!
!real(8) pure function scatterpdf(omega, p, T) result(pdf)
!   real(8), intent(in) :: omega, T
!   integer, intent(in) :: p
!   
!   pdf = dos(omega, p)/tau(omega, p, T)*dedT(omega, T)
!end function scatterpdf

subroutine initdist()
    allocate( energypdf_arr(npol, nomega) )
    allocate( fluxpdf_arr(npol, nomega) )
    allocate( scatterpdf_arr(npol, nomega) )
!   call alloc(energypdf_arr, npol, nomega)
!   call alloc(fluxpdf_arr, npol, nomega)
!   call alloc(scatterpdf_arr, npol, nomega)
    
    energypdf_arr = dedT_arr*dos_arr*domega_arr
    fluxpdf_arr = vel_arr*dedT_arr*dos_arr*domega_arr
    scatterpdf_arr = dedT_arr/tau_arr*dos_arr*domega_arr
    
    allocate( energycdf_arr(0:npol, nomega) )
    allocate( fluxcdf_arr(0:npol, nomega) )
    allocate( scattercdf_arr(0:npol, nomega) )
    
    energycdf_arr = calculatecdf(energypdf_arr)
    fluxcdf_arr = calculatecdf(fluxpdf_arr)
    scattercdf_arr = calculatecdf(scatterpdf_arr)
end subroutine initdist

!subroutine initomega_int(n)
!   integer, intent(in) :: n
!   integer :: i
!
!   nomega = n
!   domega = omegamax/n
!   allocate( omega_arr(n) )
!   omega_arr = (/ ((i - 0.5)*domega, i=1,n) /)
!end subroutine initomega_int
!
!subroutine initomega_real(delta)
!   real(8), intent(in) :: delta
!   
!   call initomega_int( nint(omegamax/delta) )
!end subroutine initomega_real

!subroutine calculatecdf(pdf_arr, cdf_arr, pdf, T)
!   real(8), external :: pdf
!   real(8), intent(in) :: T
!   real(8), intent(out), allocatable :: pdf_arr(:,:), cdf_arr(:,:)
!   integer :: i, p
!   real(8) :: flat_arr( npol*nomega )
!   
!   allocate( pdf_arr(npol, nomega) )
!   flat_arr = (/((pdf(omega_arr(i), p, T), p=1,npol), i=1,nomega)/)
!   pdf_arr = reshape(flat_arr, (/npol, nomega/))
!   
!   allocate( cdf_arr(0:npol, nomega) )
!   cdf_arr(0,:) = pdftocdf( sum(pdf_arr, 1) )
!   do i = 1,nomega
!       cdf_arr(1:npol,i) = pdftocdf( pdf_arr(:,i) )
!   end do
!end subroutine calculatecdf

function calculatecdf(pdf_arr) result(cdf_arr)
    real(8), intent(in) :: pdf_arr(npol, nomega)
    real(8) :: cdf_arr(0:npol, nomega)
    integer :: q
    
    cdf_arr(0,:) = pdftocdf( sum(pdf_arr, 1) )
    do q = 1,nomega
        cdf_arr(1:npol,q) = pdftocdf( pdf_arr(:,q) )
    end do
end function calculatecdf

pure function getomegacdf(cdf_arr) result(omegacdf)
    real(8), intent(in) :: cdf_arr(0:npol, nomega)
    real(8) :: omegacdf(nomega)
    
    omegacdf = cdf_arr(0,:)
end function getomegacdf

pure function getpolcdf(cdf_arr, ind) result(polcdf)
    real(8), intent(in) :: cdf_arr(0:npol, nomega)
    integer, intent(in) :: ind
    real(8) :: polcdf(npol)
    
    polcdf = cdf_arr(1:npol, ind)
end function getpolcdf

!subroutine initpropcdf(T)
!   real(8), intent(in) :: T
!   
!   Teq = T
!   call calculatecdf(energypdf_arr, energycdf_arr, energypdf, T)
!   call calculatecdf(fluxpdf_arr, fluxcdf_arr, fluxpdf, T)
!   call calculatecdf(scatterpdf_arr, scattercdf_arr, scatterpdf, T)
!end subroutine initpropcdf

real(8) pure function getpseudoenergy() result(energy)
    energy = sum(energypdf_arr)
end function getpseudoenergy

real(8) pure function getpseudoflux() result(flux)
    flux = sum(fluxpdf_arr)/4
end function getpseudoflux

subroutine drawprop(p, q, cdf_arr)
    real(8), intent(in) :: cdf_arr(:,:)
    integer, intent(out) :: p, q
    real(8) :: r(2)
    
    call randnum(r)
    q = searchintvl( getomegacdf(cdf_arr), r(1) )
    p = searchintvl( getpolcdf(cdf_arr, q), r(2) )
end subroutine drawprop

subroutine drawenergyprop(p, q)
    integer, intent(out) :: p, q
    
    call drawprop(p, q, energycdf_arr)
end subroutine drawenergyprop

subroutine drawfluxprop(p, q)
    integer, intent(out) :: p, q
    
    call drawprop(p, q, fluxcdf_arr)
end subroutine drawfluxprop

subroutine drawscatterprop(p, q)
    integer, intent(out) :: p, q
    
    call drawprop(p, q, scattercdf_arr)
end subroutine drawscatterprop

subroutine drawposrect(pos)
    real(8), intent(out) :: pos(3)
    real(8) :: r(2)
    
    pos = 0
    call randnum(r)
    pos(1:2) = r
end subroutine drawposrect

subroutine drawangiso(dir)
    real(8), intent(out) :: dir(3)
    real(8) :: r(2), mu, phi
    
    call randnum(r)
    mu = 1 - 2*r(1)
    phi = pi*(1 - 2*r(2))
    dir = angtodir((/ mu, phi /))
end subroutine drawangiso

subroutine drawanghalf(dir)
    real(8), intent(out) :: dir(3)
    real(8) :: r(2), mu, phi
    
    call randnum(r)
    mu = sqrt(1 - r(1))
    phi = pi*(1 - 2*r(2))
    dir = angtodir((/ mu, phi /))
end subroutine drawanghalf

subroutine drawscattime(tscat, p, q)
    real(8), intent(out) :: tscat
    integer, intent(in) :: p, q
    real(8) :: r
    
    call randnum(r)
    tscat = -tau_arr(p,q)*log(r)
end subroutine drawscattime

end module distributions