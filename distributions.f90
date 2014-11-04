module distributions
    use constants
    use tools
    use math
    use material
    implicit none
    
    public
    
    real(8), allocatable :: energypdf_arr(:,:), fluxpdf_arr(:,:), scatterpdf_arr(:,:)
    real(8), allocatable :: energycdf_arr(:,:), fluxcdf_arr(:,:), scattercdf_arr(:,:)
    
contains

subroutine initdist()
    allocate( energypdf_arr(npol, nomega) )
    allocate( fluxpdf_arr(npol, nomega) )
    allocate( scatterpdf_arr(npol, nomega) )
    
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