program test
    use tools
    implicit none
    
    real(8) :: a, b
    character(64) :: cmd

!     call testsystem()
    
    call testparse()

!     print ('(A)'), timestamp()
!     call starttimer()
!     call testprogress()
contains


subroutine testsystem()
    integer :: status
    
    call parsecmdargs()
    print ('(A)'), cmd
    call system(cmd, status)
    if (status /= 0) call exit(status)
end subroutine testsystem

subroutine testparse()
    a = 0d0
    b = 0d0
    call parsecmdargs()
    print *, '   a = ', a
    print *, '   b = ', b
end subroutine testparse

subroutine parsecmdargs()
    integer :: i, eqsign, length
    character(64) :: input, arg, val
    character(66) :: charval
    
    do i = 1, command_argument_count()
        call get_command_argument(i, input)
        input = adjustl(input)
        eqsign = index(input, '=')
        length = len_trim(input)
        if (eqsign > 1 .and. eqsign < length) then
            arg = input(1:eqsign-1)
            val = input(eqsign+1:length)
            
            select case (arg)
                case ('a')
                    read(val,*) a
                case ('b')
                    read(val,*) b
                case ('cmd')
                    charval = '"'// val //'"'
                    read(charval,*) cmd
                case default
                    print *, 'Invalid input: "', trim(input), '"'
            end select
        else
            print *, 'Invalid input: "', trim(input), '"'
        end if
    end do
end subroutine parsecmdargs

subroutine testprogress()
    integer :: i
    integer, parameter :: N = 20
    
    do i = 1, N
        call showprogress(i, N)
        call sleep(1)
    end do
end subroutine testprogress

end program test