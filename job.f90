program job
    implicit none
    
    character(7),  parameter :: filename = 'job.pbs'
    character(11), parameter :: hashbang = '#!/bin/bash'
    character(5),  parameter :: pbs = '#PBS '
    
    integer :: unit = 1
    integer :: i, j
    integer :: iters
    character(128) :: logdir, jobname, exec
    character(32) :: c(3), d(4)
    
    c = (/'4 ', '8 ', '16'/)
    d = (/'0.01', '0.02', '0.04', '0.08'/)
    iters = 10 
    
    logdir = ''
    do i = 1,3
        do j = 1,4
            jobname = 'octet_unit_08_02_'// trim(c(i)) //'_0'// d(j)(3:4) //'_1e6_1e5_5d-6'
            exec = 'time ./run unit a=0.8d-6 b=0.2d-6 c='// &
                trim(c(i)) //'d-6 d='// trim(d(j)) //&
                'd-6 nemit=1000000 maxscat=100000 tend=5d-6'
            call writejob()
            print ('(A,/,A)'), trim(jobname), trim(exec)
            call system('qsub ' // filename)
        end do
    end do
    
contains
    
subroutine writejob()
    integer :: k
    
    open(unit, file=filename, action='write', status='replace')
    write(unit, '(A)') hashbang
    write(unit, '(A)') pbs // '-N ' // trim(jobname)
    write(unit, '(A)') pbs // '-o ' // trim(logdir) // trim(jobname) // '.log'
    write(unit, '(A)') pbs // '-j ' // 'oe'
    write(unit, '(A)') pbs // '-l nodes=1:ppn=12'
    write(unit, '(A)') pbs // '-l walltime=24:00:00'
    write(unit, '(A)') pbs // '-q default'
    write(unit, '(/,A)') 'cd $PBS_O_WORKDIR'
    do k = 1, iters
        write(unit, '(A)') trim(exec)
    end do
    close(unit)
end subroutine writejob

end program job
