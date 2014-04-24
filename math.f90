module math
	implicit none
	
	public
	
	interface cumsum
		module procedure cumsum_real, cumsum_arr
	end interface cumsum

contains

pure function cumsum_real(arr) result(cum_arr)
	real(8), intent(in) :: arr(:)
	real(8), allocatable :: cum_arr(:)
	real(8) :: cum
	integer :: n, l, u, i
	
	n = size(arr)
	allocate( cum_arr(n) )
	
	l = lbound(arr, 1)
	u = ubound(arr, 1)
	cum = 0d0
	do i = l,u
		cum = cum + arr(i)
		cum_arr(i) = cum
	end do
end function cumsum_real

pure function cumsum_arr(arr, dim) result(cum_arr)
	real(8), intent(in) :: arr(:,:)
	integer, intent(in) :: dim
	real(8), allocatable :: cum_arr(:,:)
	integer :: m, n, l, u, k
	
	m = size(arr, 1);
	n = size(arr, 2);
	allocate( cum_arr(m,n) )
	
	if (dim == 1) then
		l = lbound(arr, 2)
		u = ubound(arr, 2)
		do k = l,u
			cum_arr(:,k) = cumsum_real( arr(:,k) )
		end do
		
	else if (dim == 2) then
		l = lbound(arr, 1)
		u = ubound(arr, 1)
		do k = l,u
			cum_arr(k,:) = cumsum_real( arr(k,:) )
		end do
		
	else
		cum_arr = arr
	end if
end function cumsum_arr

pure function pdftocdf(pdf) result(cdf)
	real(8), intent(in) :: pdf(:)
	real(8), allocatable :: cdf(:)
	integer :: n
	
	n = size(pdf)
	allocate( cdf(n) )
	cdf = cumsum(pdf)
	if (cdf(n) /= 0) then
		cdf = cdf/cdf(n)
	end if
end function pdftocdf

integer pure function searchintvl(arr, val) result(ind)
	real(8), intent(in) :: arr(:)
	real(8), intent(in) :: val
	integer :: low, high, mid
	
	low = lbound(arr, 1)
	high = ubound(arr, 1)
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
	if (norm /= 0) then
		unitvec = vec / norm
	else
		unitvec = vec
	end if
end function unitvec

pure function project(a, b) result(proj)
	real(8), intent(in) :: a(3), b(3)
	real(8) :: dir(3), proj(3)
	
	dir = unitvec(b)
	proj = dot_product(a, dir)*dir
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
	real(8) :: mu, lambda, phi, dir(3)
	
	mu = ang(1)
	lambda = sqrt(1 - mu**2)
	phi = ang(2)
	dir = (/ lambda*cos(phi), lambda*sin(phi), mu /)
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

pure function rotmatrix(xaxis, yaxis, zaxis) result(mat)
	real(8), intent(in) :: xaxis(3), yaxis(3), zaxis(3)
	real(8) :: mat(3,3)
	
	mat(:,1) = xaxis
	mat(:,2) = yaxis
	mat(:,3) = zaxis
end function rotmatrix

integer pure function signtoint(sign) result(pm)
	logical, intent(in) :: sign
	
	if (sign) then
		pm = 1
	else
		pm = -1
	end if
end function signtoint

end module math