module math
	implicit none
	
	public
	
	interface cumsum
		module procedure cumsum_real, cumsum_arr
	end interface cumsum

contains

subroutine initrand()
! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
	
	integer, allocatable :: seed(:)
	integer :: i, n, un, istat, dt(8), pid, t(2), s
	integer(8) :: count, tms

	call random_seed(size = n)
	allocate(seed(n))
	! First try if the OS provides a random number generator
	open(newunit=un, file="/dev/urandom", access="stream", &
		 form="unformatted", action="read", status="old", iostat=istat)
	if (istat == 0) then
	   read(un) seed
	   close(un)
	else
	   ! Fallback to XOR:ing the current time and pid. The PID is
	   ! useful in case one launches multiple instances of the same
	   ! program in parallel.
	   call system_clock(count)
	   if (count /= 0) then
		  t = transfer(count, t)
	   else
		  call date_and_time(values=dt)
		  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
			   + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
			   + dt(3) * 24 * 60 * 60 * 60 * 1000 &
			   + dt(5) * 60 * 60 * 1000 &
			   + dt(6) * 60 * 1000 + dt(7) * 1000 &
			   + dt(8)
		  t = transfer(tms, t)
	   end if
	   s = ieor(t(1), t(2))
	   pid = getpid() + 1099279 ! Add a prime
	   s = ieor(s, pid)
	   if (n >= 3) then
		  seed(1) = t(1) + 36269
		  seed(2) = t(2) + 72551
		  seed(3) = pid
		  if (n > 3) then
			 seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
		  end if
	   else
		  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
	   end if
	end if
	call random_seed(put=seed)
end subroutine initrand

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