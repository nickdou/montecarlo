module math
	implicit none
	
	public
	
! 	interface cumsum
! 		module procedure cumsum_real, cumsum_arr
! 	end interface cumsum

contains

pure function cumsum(arr) result(cum_arr)
	real(8), intent(in) :: arr(:)
	real(8) :: cum_arr(size(arr))
	real(8) :: cum
	integer :: i
	
	cum = 0d0
	do i = 1,size(arr)
		cum = cum + arr(i)
		cum_arr(i) = cum
	end do
end function cumsum

! pure function cumsum_arr(arr, dim) result(cum_arr)
! 	real(8), intent(in) :: arr(:,:)
! 	integer, intent(in) :: dim
! 	real(8), allocatable :: cum_arr(:,:)
! 	integer :: m, n, l, u, k
!
! 	m = size(arr, 1);
! 	n = size(arr, 2);
! !	allocate( cum_arr(m,n) )
! 	call alloc(cum_arr, m, n)
!
! 	if (dim == 1) then
! 		l = lbound(arr, 2)
! 		u = ubound(arr, 2)
! 		do k = l,u
! 			cum_arr(:,k) = cumsum_real( arr(:,k) )
! 		end do
!
! 	else if (dim == 2) then
! 		l = lbound(arr, 1)
! 		u = ubound(arr, 1)
! 		do k = l,u
! 			cum_arr(k,:) = cumsum_real( arr(k,:) )
! 		end do
!
! 	else
! 		cum_arr = arr
! 	end if
! end function cumsum_arr

pure function pdftocdf(pdf) result(cdf)
	real(8), intent(in) :: pdf(:)
	real(8) :: cdf(size(pdf))
	integer :: n
	
	cdf = cumsum(pdf)
	n = size(cdf)
	if (cdf(n) /= 0) then
		cdf = cdf/cdf(n)
	end if
end function pdftocdf

integer pure function searchintvl(arr, val) result(ind)
	real(8), intent(in) :: arr(:)
	real(8), intent(in) :: val
	integer :: low, mid, high
	
	low = 1
	high = size(arr)
	do while (low < high)
		mid = (low + high)/2
		if (arr(mid) > val) then
			high = mid
		else
			low = mid+1
		end if
	end do
	ind = low
end function searchintvl

real(8) pure function normtwo(vec)
	real(8), intent(in) :: vec(:)
	
	normtwo = sqrt( sum(vec**2) )
end function normtwo

pure function unitvec(vec)
	real(8), intent(in) :: vec(:)
	real(8) :: unitvec(size(vec))
	real(8) :: norm
	
	norm = normtwo(vec)
	
	unitvec = vec
	if (norm /= 0) then
		unitvec = vec / norm
	end if
end function unitvec

pure function project(vec, dir) result(proj)
	real(8), intent(in) :: vec(3), dir(3)
	real(8) :: unit(3), proj(3)
	
	unit = unitvec(dir)
	proj = dot_product(vec, unit)*unit
end function project

pure function cross_product(a, b) result(c)
	real(8), intent(in) :: a(3), b(3)
	real(8) :: c(3)
	
	c = (/ a(2)*b(3) - a(3)*b(2), &
	       a(3)*b(1) - a(1)*b(3), &
	       a(1)*b(2) - a(2)*b(1) /)
end function cross_product

pure function angtodir(ang) result(dir)
	real(8), intent(in) :: ang(2)
	real(8) :: costheta, sintheta, phi
	real(8) :: dir(3)
	
	costheta = ang(1)
	sintheta = sqrt(1 - costheta**2)
	phi = ang(2)
	
	dir = (/ sintheta*cos(phi), sintheta*sin(phi), costheta /)
end function angtodir

pure function dirtoang(dir) result(ang)
	real(8), intent(in) :: dir(3)
	real(8) :: x, y, z, r, ang(2)
	
	x = dir(1)
	y = dir(2)
	z = dir(3)
	r = normtwo(dir)
	
	if (r == 0 .or. r == z) then
		ang = (/ 1, 0 /)
	else if (r == -z) then
		ang = (/ -1, 0 /)
	else
		ang = (/ z/r, atan2(y, x) /)
	end if
end function dirtoang

! pure function rotmatrix(xaxis, yaxis, zaxis) result(mat)
! 	real(8), intent(in) :: xaxis(3), yaxis(3), zaxis(3)
! 	real(8) :: mat(3,3)
!
! 	mat(:,1) = xaxis
! 	mat(:,2) = yaxis
! 	mat(:,3) = zaxis
! end function rotmatrix

pure function inverse2(mat) result(inv)
	real(8), intent(in) :: mat(2,2)
	real(8) :: det, inv(2,2)
	
	det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
	inv = 1/det * reshape((/mat(2,2), -mat(2,1), -mat(1,2), mat(1,1)/), (/2,2/))
end function inverse2

integer pure function signtoint(sign) result(pm)
	logical, intent(in) :: sign
	
	if (sign) then
		pm = 1
	else
		pm = -1
	end if
end function signtoint

end module math