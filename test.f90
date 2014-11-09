program test
    use tools
    implicit none
    
    real(8) :: a, b
    character(64) :: cmd
    real(8) :: xold(3), xnew(3), bounds(2,3)
    
    xold = (/1d0, 1d0, 1d0/)
    xnew = (/1d0, 1d0, 1d0/)
    bounds = reshape((/0d0, 2d0, 0d0, 2d0, 0d0, 2d0/), (/2, 3/))
    
    print *, xold
    print *, bounds(1,:)
    print *, bounds(2,:)
    print *, all(xold >= bounds(1,:) .and. xold <= bounds(2,:))
    
    
contains

subroutine testminloc(l, u)
    integer, intent(in) :: l, u
    real(8) :: arr(l:u), arr1(size(arr))
    logical :: mask(l:u)
    
    call initrand(.false.)
    call randnum(arr)
    arr(l) = -1d0
    arr1 = arr
    mask = .true.
    mask(l) = .false.
    print *, arr
    print *, arr1
    print *, mask
    print *, minloc(arr, 1, mask), minloc(arr1, 1, mask)
    print *, minloc(arr, 1), minloc(arr1, 1)
    print *, minval(arr, 1, mask)
end subroutine testminloc

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
