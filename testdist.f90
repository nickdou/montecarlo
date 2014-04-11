program testdist
	use material
	use distributions
	implicit none
	
	integer :: i
	real(8) :: corner(3), dir1(3), dir2(3), pos(3)
	real(8) :: omega, dir(3), xaxis(3), yaxis(3), zaxis(3), rot(3,3), eye(3,3)
	integer :: nomega, p
	
	print *, 'omegamax = ', omegamax
	
	print *, 'velocity = ', vel(omegamax, npol)
	
	print *, 'dos(omega = 0) = ', dos(0d0, 1)
	
	print *, 'dos(omega = omegamax) = ', dos(omegamax, 1)
	
	print *, 'tau = ', tau(omegamax, 1, 300d0)
	
	call initomega(1d12)
	nomega = getnomega()
	print *, 'nomega = ', nomega
	
	call initpropcdf(300d0)
	print *, 'Teq = ', getteq()
	
!	do i = 1,nomega
!		print ('(E12.5, E12.5, E12.5, E12.5)'), omega_arr(i), polcdf_arr(i,:)
!	end do
	
	call initrand()
!	call testpos()
!	call testprop()
	call testang()
	
contains

subroutine testpos()
	corner = (/ 0, 0, 0 /)
	dir1 = (/ 1, 0, 0 /)
	dir2 = (/ 0, 1, 1 /)
	print *, 'Draw 1000 random positions:'
	open(unit=2, file='drawpos.m', action='write', status='replace')
	write(2,'(A)') 'data = ['
	do i=1,1000
		call drawposquad(pos, corner, dir1, dir2)
		write(2,'(3F15.8)') pos
	end do
	write(2,'(A)') '];'
	write(2,'(A)') 'figure(1)'
	write(2,'(A)') 'scatter3(data(:,1), data(:,2), data(:,3))'
	write(2,'(A,/,A,/,A)') 'xlabel(''x'')', 'ylabel(''y'')', 'zlabel(''z'')'
end subroutine testpos

subroutine testprop()
	print *, 'Draw 1000 random (omega, p):'
	open(unit=2, file='drawomega.m', action='write', status='replace')
	write(2,'(A)') 'function drawomega(nbins)'
	write(2,'(A,/,4X,A,/,A)') 'if nargin < 1', 'nbins = 3;', 'end'
	write(2,'(A)') 'data = ['
	do i=1,1000
		call drawprop(omega, p)
		write(2,'(I5)') p
!		write(2,'(E15.8)') omega
	end do
	write(2,'(A)') '];'
	write(2,'(A)') 'hist(data, nbins)'
end subroutine testprop

subroutine testang()
	xaxis = unitvec((/ 0d0, 1d0, 0d0 /))
	yaxis = (/ 0d0, 0d0, 1d0 /)
	yaxis = unitvec(yaxis - project(yaxis, xaxis))
	zaxis = unitvec( cross_product(xaxis, yaxis) )
	rot = rotmatrix(xaxis, yaxis, zaxis)
	eye = matmul(transpose(rot), rot)
	print *, 'rot = '
	print ('(3F15.8)'), (rot(i,:), i=1,3)
	print *, 'Verify rotation matrix'
	print ('(3F15.8)'), (eye(i,:), i=1,3)
	print *, 'Draw 1000 random directions:'
	open(unit=2, file='drawang.m', action='write', status='replace')
	write(2,'(A)') 'data = ['
	do i=1,1000
!		call drawangiso(dir)
		call drawanghalf(dir)
		dir = matmul(rot, dir)
		write(2,'(3F15.8)') dir
!		print ('(3F15.8)'), dir
	end do
	write(2,'(A)') '];'
	write(2,'(A)') 'figure(1)'
	write(2,'(A)') 'scatter3(data(:,1), data(:,2), data(:,3))'
	write(2,'(A,/,A,/,A)') 'xlabel(''x'')', 'ylabel(''y'')', 'zlabel(''z'')'
end subroutine testang

end program testdist