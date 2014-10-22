module domain
    use omp_lib
    use tools
    use math
    use distributions
    implicit none
    
    private
    public  :: SPEC_BC, DIFF_BC, ISOT_BC, PERI_BC, axis, boundary, &
        makeaxis, setgrid, initrecord, inittraj, appendtraj, gettraj, &
        makebdry, makebdry_arr, setbdry, setbdrypair, &
        calculateemit, setvolumetric, getemit, geteeff, &
        drawemittime, drawemitstate, updatestate, applybc, &
        recordtime, recorddisp, recordloc, addcumdisp, &
        getsteadytemp, gettranstemp, getflux, getgridflux
    
    type axis
        private
        integer :: n
        real(8) :: dir(3), lo, hi, h
    end type axis
    
    type boundary
        private
        integer :: bc, pair
        real(8) :: origin(3), normal(3), ivec(3), jvec(3), area, rot(3,3), inv(3,3), temp, pairmv(3)
        logical :: tri
    end type boundary
    
    interface addcumdisp
        module procedure addcumdisp_real, addcumdisp_int
    end interface addcumdisp
    
    integer, parameter :: SPEC_BC = 1, &
                          DIFF_BC = 2, &
                          ISOT_BC = 3, &
                          PERI_BC = 4
    
    real(8) :: tend, tstep, flow(3), Eeff, volume
    integer :: ntime, ncell, nbdry, nemit, emitind, ntraj, maxtraj
    logical :: gridflux, volumetric
    type(axis) :: grid, vgen
    real(8), allocatable :: trajectory(:,:)
    real(8), allocatable :: cumdisp(:)
    real(8), allocatable :: gridtime_arr(:,:), griddisp_arr(:,:)
    integer, allocatable :: gridloc_arr(:,:,:)
    type(boundary), allocatable :: bdry_arr(:)
    integer, allocatable :: emit_arr(:)

contains

type(axis) pure function makeaxis(dir, lo, hi, n) result(ax)
    real(8), intent(in) :: dir(3), lo, hi
    integer, intent(in) :: n
    
    ax%dir = dir
    ax%lo = lo
    ax%hi = hi
    ax%n = n
    if (n /= 0) then
        ax%h = (hi - lo)/n
    else
        ax%h = 0
    end if
end function makeaxis

subroutine setgrid(t, nt, ax)
    real(8), intent(in) :: t
    integer, intent(in) :: nt
    type(axis), intent(in) :: ax
    
    tend = t
    ntime = nt
    
    grid = ax
    ncell = ax%n
end subroutine setgrid

subroutine initrecord(gf, dir)
    logical, intent(in) :: gf
    real(8), intent(in) :: dir(3)
    
    if (ntime == 0) then
        allocate( gridtime_arr(0:nthreads-1, ncell) )
!       call alloc(gridtime_arr, ncell)
        gridtime_arr = 0
    else
        tstep = tend/ntime
        
        allocate( gridloc_arr(0:nthreads-1, ncell, 0:ntime) )
!       call alloc(gridloc_arr, ncell, ntime, (/.false.,.true./))
        
        gridloc_arr = 0
    end if
    
    gridflux = gf
    flow = dir
    if (gf) then
        allocate( griddisp_arr(0:nthreads-1, ncell) )
!       call alloc(griddisp_arr, ncell)
        griddisp_arr = 0
    else
        allocate( cumdisp(0:nthreads-1) )
        cumdisp = 0
    end if
end subroutine initrecord

subroutine inittraj(num, pos)
    integer, intent(in) :: num
    real(8), intent(in) :: pos(3)
    
    maxtraj = num
    ntraj = 0
    
    allocate( trajectory(3, 0:num) )
!   call alloc(trajectory, 3, num, (/.false.,.true./))
    trajectory(:,0) = pos
end subroutine inittraj

subroutine appendtraj(pos)
    real(8), intent(in) :: pos(3)
    
    if (allocated(trajectory)) then
        if (ntraj < maxtraj) then
            ntraj = ntraj + 1
            trajectory(:,ntraj) = pos
        end if
    end if
end subroutine appendtraj

function gettraj() result(traj)
    real(8) :: traj(3, 0:ntraj)
    
!   allocate( traj(3, 0:ntraj) )
    traj = trajectory(:, 0:ntraj)
end function gettraj

real(8) pure function getcoord(pos) result(coord)
    real(8), intent(in) :: pos(3)
    
    coord = (dot_product(pos, grid%dir) - grid%lo)/grid%h
end function getcoord

integer pure function coordtoind(coord) result(ind)
    real(8), intent(in) :: coord
    
    ind = min(ncell, max(1, ceiling(coord)))
end function coordtoind

integer pure function getxind(pos) result(ind)
    real(8), intent(in) :: pos(3)
    
    ind = coordtoind( getcoord(pos) )
end function getxind

integer pure function gettind(time) result(tind)
    real(8), intent(in) :: time
    
    tind = min(ntime, max(0, floor(time/tstep)))
end function gettind

type(boundary) pure function makebdry(bc, origin, ivec, jvec, temp, tri) result(bdry)
    integer, intent(in) :: bc
    real(8), intent(in) :: origin(3), ivec(3), jvec(3)
    real(8), intent(in), optional :: temp
    logical, intent(in), optional :: tri
    real(8) :: kvec(3), dotij, normi, normk, rot(3,3), inv(3,3)
    
    bdry%bc = bc
    bdry%origin = origin
    bdry%ivec = ivec
    bdry%jvec = jvec
    
    dotij = dot_product(ivec, jvec)
    kvec = cross_product(ivec, jvec)
    normi = normtwo(ivec)
    normk = normtwo(kvec)
    
    rot(:,1) = ivec/normi
    rot(:,2) = (jvec - ivec*(dotij/normi**2))*(normi/normk)
    rot(:,3) = kvec/normk
    
    inv(1,:) = rot(:,1)/normi - rot(:,2)*(dotij/(normi*normk))
    inv(2,:) = rot(:,2)*(normi/normk)
    inv(3,:) = rot(:,3)/normk
    
    bdry%area = normk
    bdry%normal = rot(:,3)
    bdry%rot = rot
    bdry%inv = inv
    
    if (present(tri)) then
        bdry%tri = tri
    else
        bdry%tri = .false.
    end if
    
    if (bc == ISOT_BC .or. bc == PERI_BC) then
        bdry%temp = temp
    else
        bdry%temp = 0d0
    end if
    
    bdry%pairmv = 0d0
end function makebdry

pure function makebdry_arr(verts, bcs, faces, temps, tris) result(bdry_arr)
    real(8), intent(in) :: verts(:,:), temps(:)
    integer, intent(in) :: bcs(:), faces(:,:)
    logical, intent(in) :: tris(:)
    type(boundary), allocatable :: bdry_arr(:)
    integer :: i, n
    real(8) :: origin(3), ivec(3), jvec(3)
    
    n = size(bcs,1)
    allocate( bdry_arr(n) )
    do i = 1,n
        origin = verts(:,faces(2,i))
        ivec = verts(:,faces(1,i)) - origin
        jvec = verts(:,faces(3,i)) - origin
        bdry_arr(i) = makebdry(bcs(i), origin, ivec, jvec, temps(i), tris(i))
    end do
    
end function makebdry_arr

subroutine setbdry(arr, vol)
    type(boundary), intent(in) :: arr(:)
    real(8), intent(in) :: vol
    
    nbdry = size(arr)
    if (allocated(bdry_arr)) then
        deallocate( bdry_arr )
    end if
    allocate( bdry_arr(nbdry) )
    bdry_arr = arr
    
    volume = vol
end subroutine setbdry

subroutine setbdrypair(ind1, ind2, mv)
! TODO: check that rot2 = inv(rot1) or calculate one from the other
    integer, intent(in) :: ind1, ind2
    real(8), intent(in) :: mv(3)
    real(8) :: temp1, temp2
    
    temp1 = bdry_arr(ind1)%temp
    temp2 = bdry_arr(ind2)%temp
    
    bdry_arr(ind1)%bc = PERI_BC
    bdry_arr(ind1)%pair = ind2
    bdry_arr(ind1)%pairmv = mv
    bdry_arr(ind1)%temp = temp1 - temp2
    
    bdry_arr(ind2)%bc = PERI_BC
    bdry_arr(ind2)%pair = ind1
    bdry_arr(ind2)%pairmv = -mv
    bdry_arr(ind2)%temp = temp2 - temp1
end subroutine setbdrypair

subroutine calculateemit(num)
    integer, intent(in) :: num
    real(8) :: areatemp_arr(nbdry)
    real(8) :: sumareatemp
    integer :: i
    
    areatemp_arr = (/(bdry_arr(i)%area*abs(bdry_arr(i)%temp), i=1,nbdry)/)
    sumareatemp = sum(areatemp_arr)
    
    if (.not. allocated(emit_arr)) then
        allocate( emit_arr(nbdry) )
    else if (size(emit_arr, 1) /= nbdry) then
        deallocate( emit_arr )
        allocate( emit_arr(nbdry) )
    end if
!   call alloc(emit_arr, nbdry)
    
    emit_arr = nint( num*areatemp_arr/sumareatemp )
    nemit = sum(emit_arr)
    emitind = 1
    
    Eeff = (sumareatemp*tend/nemit) * getpseudoflux()
end subroutine calculateemit

pure function getemit() result(arr)
    integer :: i, arr(nbdry)
    
    arr = 0
    do i = 1,nbdry
        arr(i) = emit_arr(i)
    end do
end function getemit

real(8) pure function geteeff() result(E)
    E = Eeff
end function geteeff

subroutine setvolumetric(vol, ax)
    logical, intent(in) :: vol
    type(axis), intent(in), optional :: ax
    
    volumetric = vol
    if (vol) then
        vgen = ax
    end if
end subroutine setvolumetric

subroutine drawemittime(time)
    real(8), intent(out) :: time
    real(8) :: r
    
    if (ntime == 0) then
        time = 0
    else
        call randnum(r)
        time = r*tend
    end if
end subroutine drawemittime

subroutine drawemitstate(ind, pos, dir, sign)
    integer, intent(out) :: ind
    real(8), intent(out) :: pos(3), dir(3)
    logical, intent(out) :: sign
    real(8) :: r
    
    if (nemit > 0) then
        !$omp critical
        do
            if (emit_arr(emitind) > 0) then
                emit_arr(emitind) = emit_arr(emitind) - 1
                nemit = nemit - 1
                exit
            else
                emitind = emitind + 1
            end if
        end do
        !$omp end critical
        ind = emitind
        
        call drawposrect(pos)
!       pos = (/bdry_arr(emitind)%dx, bdry_arr(emitind)%dy, 1d0/)*pos
!       pos = bdry_arr(emitind)%origin + matmul(bdry_arr(emitind)%rot, pos)
        pos = bdry_arr(emitind)%origin + pos(1)*bdry_arr(emitind)%ivec + &
            pos(2)*bdry_arr(emitind)%jvec
        
        if (volumetric) then
            call randnum(r)
!           r = 0.99d0 !DEBUG
            pos = pos - project(pos, vgen%dir)
            pos = pos + ((1-r)*vgen%lo + r*vgen%hi) * vgen%dir
        end if
        
        call drawanghalf(dir)
        dir = matmul(bdry_arr(emitind)%rot, dir)
!       dir = (/0d0, 0d0, 1d0/) !DEBUG
        
        sign = (bdry_arr(emitind)%temp >= 0)
    end if
end subroutine drawemitstate

! pure subroutine getcoll(pierced, coll, pos1, pos2, tri)
!   real(8), intent(in) :: pos1(3), pos2(3)
!   logical, intent(in) :: tri
!   logical, intent(out) :: pierced
!   real(8), intent(out) :: coll(3)
!   real(8) :: k1, k2
!
!   k1 = pos1(3)
!   k2 = pos2(3)
!   coll = (k2*pos1 - k1*pos2)/(k2 - k1)
!
!   if (k1*k2 > 0) then
!       pierced = .false.
!   else if (tri) then
!       pierced = all(coll(1:2) >= 0) .and. sum(coll(1:2)) <= 1
!   else
!       pierced = all(coll(1:2) >= 0 .and. coll(1:2) <= 1)
!   end if
! end subroutine getcoll

pure subroutine getcoll(pierced, coll, x, dir, tri)
    real(8), intent(in) :: x(3), dir(3)
    logical, intent(in) :: tri
    logical, intent(out) :: pierced
    real(8), intent(out) :: coll(2)
    real(8) :: dist
    
    dist = -x(3)/dir(3)
    coll = x(1:2) + dist*dir(1:2)
    
    if (dist < 0) then
        pierced = .false.
    else if (tri) then
        pierced = all(coll >= 0) .and. sum(coll) <= 1
    else
        pierced = all(coll >= 0 .and. coll <= 1)
    end if
end subroutine getcoll

subroutine updatestate(bc, ind, x, dir, v, deltat, t)
    integer, intent(inout) :: ind
    real(8), intent(inout) :: x(3), dir(3), v, deltat, t
    integer, intent(out) :: bc
    integer :: i
!   real(8) :: aim(3), origin(3), ivec(3), jvec(3), inv(3,3), pos1(3), pos2(3)
    real(8) :: xi(3), diri(3)
!   logical :: tri, piercedi, pierced(nbdry)
    logical :: pierced(0:nbdry)
!   real(8) :: colli(3), coll(3,nbdry), dist(nbdry)
    real(8) :: colli(2), coll(3, 0:nbdry), dt(0:nbdry)
    
!   if (t + deltat > tend) then
!       deltat = tend - t
!   end if
!   aim = x + v*deltat*dir
    dt(0) = min(deltat, tend - t)
    coll(:,0) = x + v*dt(0)*dir
    
    do i = 1,nbdry
!       origin = bdry_arr(i)%origin
!       ivec = bdry_arr(i)%ivec
!       jvec = bdry_arr(i)%jvec
!       inv = bdry_arr(i)%inv
!       tri = bdry_arr(i)%tri
!
!       pos1 = matmul(inv, x - origin)
!       pos2 = matmul(inv, aim - origin)
!       call getcoll(piercedi, colli, pos1, pos2, tri)
!
!       pierced(i) = piercedi
!       coll(:,i) = origin + colli(1)*ivec + colli(2)*jvec
!       dist(i) = normtwo(coll(:,i) - x)
        
        xi = matmul(bdry_arr(i)%inv, x - bdry_arr(i)%origin)
        diri = matmul(bdry_arr(i)%inv, dir)
        call getcoll(pierced(i), colli, xi, diri, bdry_arr(i)%tri)
        
        coll(:,i) = bdry_arr(i)%origin + &
            colli(1)*bdry_arr(i)%ivec + colli(2)*bdry_arr(i)%jvec
        dt(i) = normtwo(coll(:,i) - x)/v
    end do
    
    pierced(ind) = .false.
    pierced(0) = .true.
    ind = minloc(dt, 1, pierced) - 1 !subtract one because of zero index
!   print *, dt
!   print *, pierced
    
    if (ind == 0) then
        bc = 0
    else
        bc = bdry_arr(ind)%bc
    end if
    x = coll(:,ind)
    deltat = dt(ind)
    
    t = t + deltat
    if (t >= tend) then
        bc = 0
        dir = 0
    end if
end subroutine updatestate

subroutine applybc(bc, ind, x, dir)
    integer, intent(in) :: bc
    integer, intent(inout) :: ind
    real(8), intent(inout) :: x(3), dir(3)
    
    if (bc /= 0 .and. any(abs(dir) >= eps, 1)) then
        select case (bc)
        case (SPEC_BC)
            dir = dir - 2*project(dir, bdry_arr(ind)%normal)
        case (DIFF_BC)
            call drawanghalf(dir)
            dir = matmul(bdry_arr(ind)%rot, dir)
        case (ISOT_BC)
            dir = 0
        case (PERI_BC)
            x = x + bdry_arr(ind)%pairmv
            ind = bdry_arr(ind)%pair
        end select
    end if
end subroutine applybc

!subroutine incgridtime(time, ind1, ind2)
!   real(8), intent(in) :: time
!   integer, intent(in) :: ind1
!   integer, intent(in), optional :: ind2
!   
!   if (present(ind2)) then
!       if (ind1 <= ind2) then
!           gridtime_arr(ind1:ind2) = gridtime_arr(ind1:ind2) + time
!       end if
!   else
!       gridtime_arr(ind1) = gridtime_arr(ind1) + time
!   end if
!end subroutine incgridtime

!subroutine recordtime(sign, deltat, xold, xnew)
!   logical, intent(in) :: sign
!   real(8), intent(in) :: deltat, xold(3), xnew(3)
!   real(8) :: coordold, coordnew, dcoord
!   integer :: pm, indold, indnew
!   
!   if (ntime == 0) then
!       pm = signtoint(sign)
!   
!       coordold = getcoord(xold)
!       coordnew = getcoord(xnew)
!       dcoord = abs(coordnew - coordold)
!   
!       indold = coordtoind(coordold)
!       indnew = coordtoind(coordnew)
!   
!       if (indold < indnew) then
!           call incgridtime(pm*deltat*(indold - coordold)/dcoord, indold)
!           call incgridtime(pm*deltat/dcoord, indold+1, indnew-1)
!           call incgridtime(pm*deltat*(1 - (indnew - coordnew))/dcoord, indnew)
!       else if (indnew < indold) then
!           call incgridtime(pm*deltat*(indnew - coordnew)/dcoord, indnew)
!           call incgridtime(pm*deltat/dcoord, indnew+1, indold-1)
!           call incgridtime(pm*deltat*(1 - (indold - coordold))/dcoord, indold)
!       else
!           call incgridtime(pm*deltat, indold)
!       end if
!   
!!      print ('(2I10)'), indold, indnew
!!      print ('(I4,F8.3,I4)'), floor(coordold), coordold, ceiling(coordold)
!!      print ('(2F10.3,A,2ES11.3,A,2ES11.3)'), coordold, coordnew, ' :', &
!!          gridtime_arr(1:2), ' ... ', gridtime_arr(ncell-1:ncell)
!   end if
!end subroutine recordtime

subroutine incgrid(cum_arr, val, ind1, ind2)
    real(8), intent(inout) :: cum_arr(0:,:)
    real(8), intent(in) :: val
    integer, intent(in) :: ind1
    integer, intent(in), optional :: ind2
    integer :: threadi = 0
    
!$  threadi = omp_get_thread_num()
    if (present(ind2)) then
        if (ind1 <= ind2) then
            cum_arr(threadi, ind1:ind2) = cum_arr(threadi, ind1:ind2) + val
        end if
    else
        cum_arr(threadi, ind1) = cum_arr(threadi, ind1) + val
    end if
end subroutine incgrid

subroutine recordgrid(cum_arr, sign, delta, coordold, coordnew)
    real(8), intent(inout) :: cum_arr(0:,:)
    logical, intent(in) :: sign
    real(8), intent(in) :: delta, coordold, coordnew
    real(8) :: dcoord
    integer :: pm, indold, indnew
    
    pm = signtoint(sign)
    dcoord = abs(coordnew - coordold)

    indold = coordtoind(coordold)
    indnew = coordtoind(coordnew)

    if (indold < indnew) then
        call incgrid(cum_arr, pm*delta*(indold - coordold)/dcoord, indold)
        call incgrid(cum_arr, pm*delta/dcoord, indold+1, indnew-1)
        call incgrid(cum_arr, pm*delta*(1 - (indnew - coordnew))/dcoord, indnew)
    else if (indnew < indold) then
        call incgrid(cum_arr, pm*delta*(indnew - coordnew)/dcoord, indnew)
        call incgrid(cum_arr, pm*delta/dcoord, indnew+1, indold-1)
        call incgrid(cum_arr, pm*delta*(1 - (indold - coordold))/dcoord, indold)
    else
        call incgrid(cum_arr, pm*delta, indold)
    end if
end subroutine recordgrid

subroutine recordtime(sign, deltat, xold, xnew)
    logical, intent(in) :: sign
    real(8), intent(in) :: deltat, xold(3), xnew(3)
    
    if (ntime == 0) then
        call recordgrid(gridtime_arr, sign, deltat, getcoord(xold), getcoord(xnew))
    end if
end subroutine recordtime

subroutine recorddisp(sign, xold, xnew)
    logical, intent(in) :: sign
    real(8), intent(in) :: xold(3), xnew(3)
    real(8) :: disp
    
    if (gridflux) then
        disp = dot_product(xnew - xold, flow) 
        call recordgrid(griddisp_arr, sign, disp, getcoord(xold), getcoord(xnew))
    end if
end subroutine recorddisp

subroutine recordloc(sign, t, deltat, xold, xnew)
    logical, intent(in) :: sign
    real(8), intent(in) :: t, deltat, xold(3), xnew(3)
    real(8) :: xi
    integer :: pm, told, tnew, tind, xind
    integer :: threadi = 0
    
!$  threadi = omp_get_thread_num()
    if (ntime /= 0) then
        pm = signtoint(sign)
        
        told = gettind(t - deltat)
        tnew = gettind(t)
        if (told /= 0) then
            told = told + 1
        end if
        
        if (told <= tnew) then
            do tind = told, tnew
                xi = (t - tind*tstep)/deltat
                xind = getxind(xold*xi + xnew*(1 - xi))
                gridloc_arr(threadi, xind, tind) = gridloc_arr(threadi, xind, tind) + pm
            end do
        end if
    end if
end subroutine recordloc

subroutine addcumdisp_real(sign, deltax)
    logical, intent(in) :: sign
    real(8), intent(in) :: deltax(3)
    integer :: pm, threadi = 0
    
!$  threadi = omp_get_thread_num()
    if (.not. gridflux) then
        pm = signtoint(sign)
        cumdisp(threadi) = cumdisp(threadi) + pm*dot_product(flow, deltax)
!       print *, 'addcumdisp', deltax
    end if
end subroutine addcumdisp_real

subroutine addcumdisp_int(sign, ind)
    logical, intent(in) :: sign
    integer, intent(in) :: ind
    integer :: pm, threadi = 0
    
!$  threadi = omp_get_thread_num()
    if (.not. gridflux) then
        pm = signtoint(sign)
        cumdisp(threadi) = cumdisp(threadi) - pm*dot_product(flow, bdry_arr(ind)%pairmv)
    end if
end subroutine addcumdisp_int

pure function getsteadytemp() result(T)
    real(8) :: T(ncell)
    
    T = Eeff*ncell/(volume*tend*getpseudoenergy()) * sum(gridtime_arr, 1)
end function getsteadytemp

pure function gettranstemp() result(T)
    real(8) :: T(ncell, 0:ntime)
    
    T = Eeff*ncell/(volume*getpseudoenergy()) * sum(gridloc_arr, 1)
end function gettranstemp

real(8) pure function getflux() result(flux)
    
    flux = Eeff/(volume*tend) * sum(cumdisp)
end function getflux

pure function getgridflux() result(flux)
    real(8) :: flux(ncell)
    
    flux = Eeff*ncell/(volume*tend) * sum(griddisp_arr, 1)
end function getgridflux

end module domain