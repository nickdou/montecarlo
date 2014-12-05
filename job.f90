program job
    implicit none
    
    character(*), parameter :: pbsfile = 'job.pbs'
    character(32) :: input
    logical :: dosubmit
    
    dosubmit = .false.
    if (command_argument_count() > 0) then
        call get_command_argument(1, input)
        read(input,'(L1)') dosubmit
    else
        write(*,'(A)',advance='no') 'Enable job submission? (T/F) '
        read(*,'(L1)') dosubmit
    end if
    
    if (dosubmit) then
        print *, 'Job submission enabled'
    else
        print *, 'Job submission disabled'
    end if
    
    call sweepscat()
!     call sweepcd()
    
contains

subroutine sweepscat()
    character(128) :: jobname, logdir, exec(1)
    character(32) :: iter, scat
    integer :: i, j, maxscat(9)
    
    maxscat = (/(i, i=1,size(maxscat))/)
    
    logdir = ''
    do j = 1, size(maxscat)
        write(scat,*) maxscat(j)
        scat = adjustl(scat)
        
        do i = 11,20
            write(iter,'(I2.2)') i
            iter = adjustl(iter)
            
            jobname = 'octet_unit_08_02_4_004_1e5_'// trim(scat) // &
                'e4_5d-6_'// trim(iter)
            exec(1) = 'time ./run unit a=0.8d-6 b=0.2d-6 c=4d-6 d=0.04d-6 '// &
                'nemit=100000 maxscat='// trim(scat) //'0000 tend=5d-6'
        
            call writejob(jobname, logdir, exec)
            call submitjob()
        end do
    end do
    
end subroutine sweepscat

subroutine sweepcd()
    character(128) :: jobname, logdir, exec(1)
    character(32) :: c(3), d(4)
    integer :: i, j
    
    c = (/'4 ', '8 ', '16'/)
    d = (/'0.01', '0.02', '0.04', '0.08'/)
    
    logdir = ''
    do i = 1, size(c)
        do j = 1, size(d)
            jobname = 'octet_unit_08_02_'// trim(c(i)) //'_0'// d(j)(3:4) //'_1e6_1e5_5d-6'
            exec(1) = 'time ./run unit a=0.8d-6 b=0.2d-6 c='// &
                trim(c(i)) //'d-6 d='// trim(d(j)) // &
                'd-6 nemit=1000000 maxscat=100000 tend=5d-6'
            
            call writejob(jobname, logdir, exec)
            call submitjob()
        end do
    end do
end subroutine sweepcd

subroutine writejob(jobname, logdir, exec)
    character(128), intent(in) :: jobname, logdir, exec(:)
    integer, parameter :: unit = 1
    integer :: i
    
    write(*, '(A)') trim(jobname)
    
    open(unit, file=pbsfile, action='write', status='replace')
    write(unit, '(A)') '#!/bin/bash'
    write(unit, '(A)') '#PBS -N ' // trim(jobname)
    write(unit, '(A)') '#PBS -o ' // trim(logdir) // trim(jobname) // '.log'
    write(unit, '(A)') '#PBS -j oe'
    write(unit, '(A)') '#PBS -l nodes=1:ppn=12'
    write(unit, '(A)') '#PBS -l walltime=24:00:00'
    write(unit, '(A)') '#PBS -q default'
    write(unit, '(/,A)') 'cd $PBS_O_WORKDIR'
    
    do i = 1, size(exec)
        write(*, '(A)') trim(exec(i))
        write(unit, '(A)') trim(exec(i))
    end do
    close(unit)
end subroutine writejob

subroutine submitjob()
    if (dosubmit) then
        call system('qsub ' // pbsfile)
    end if
end subroutine submitjob

end program job
