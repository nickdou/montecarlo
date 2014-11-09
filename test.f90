program test
    use tools
    implicit none
    
    real(8) :: a, b
    character(64) :: cmd
    integer :: i
    integer, parameter :: l=0, u=4
    real(8) :: array(l:u), array2(size(array))
    logical :: mask(l:u)
    
    array = (/(dble(i), i=l,u)/)
    array2 = array
    mask = .true.
    mask(l) = .true.
    print *, array
    print *, array2
    print *, mask
    print *, minloc(array, 1, mask), minloc(array2, 1, mask)
    print *, minloc(array, 1), minloc(array2, 1)
    print *, minval(array, 1, mask)
    
!     open(2, file='test.txt', action='write', position='append')
!     call writearr(reshape((/(dble(i), i=1,20)/),(/5,4/)), 2)
!     close(2)
!     call testsystem()
    
!     call testparse()

!     print ('(A)'), timestamp()
!     call starttimer()
!     call testprogress()
contains

subroutine writearr(arr, un)
    real(8), intent(in) :: arr(:,0:)
    integer, intent(in), optional :: un
    character(16) :: col
    character(32) :: fmt, row
    integer :: unit
    
    if (present(un)) then
        unit = un
    else
        unit = 6
    end if
    
    write (col, *) size(arr,1)
    fmt = '(A4,'// trim(adjustl(col)) //'F8.1)'
    row = '(4X,'// trim(adjustl(col)) //'F8.1)'
    write (unit, '(A,/,A)') fmt, row
    write (unit, fmt) 'A = ', arr(:,0)
    write (unit, row) arr(:,1:ubound(arr,2))
end subroutine writearr

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
