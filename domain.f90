module domain
    use omp_lib
    use tools
    use math
    use distributions
    implicit none
    
    private
    public  :: SPEC_BC, DIFF_BC, ISOT_BC, PERI_BC, boundary, &
        settime, setgrid, initrecord, inittraj, appendtraj, gettraj, &
        makebdry, makebdry_arr, setbdry, setbdrypair, &
        calculateemit, setvolumetric, getemit, getnemit, geteeff, &
        drawemittime, drawemitstate, updatestate, applybc, &
        recordtime, recorddisp, recordloc, addcumdisp, &
        getsteadytemp, gettranstemp, getflux, getgridflux, getcond, getnstop
    
!     type axis
!         private
!         integer :: n
!         real(8) :: dir(3), lo, hi, h
!     end type axis
    
    type boundary
        private
        integer :: bc, pair
        real(8) :: origin(3), normal(3), ivec(3), jvec(3), area, rot(3,3), inv(3,3), temp, pairmv(3)
        logical :: tri
    end type boundary
    
    interface setgrid
        module procedure setgrid_none, setgrid_arr, setgrid_int
    end interface setgrid
    
    interface addcumdisp
        module procedure addcumdisp_real, addcumdisp_int
    end interface addcumdisp
    
    integer, parameter :: SPEC_BC = 1, &
                          DIFF_BC = 2, &
                          ISOT_BC = 3, &
                          PERI_BC = 4
    
!     logical :: gridflux = .false.
    logical :: volumetric = .false.
    real(8) :: tend, tstep, Eeff, gradT
    real(8) :: flow(3), voldir(3), vollim(2), griddir(3) 
    integer :: ntime, ngrid, nbdry, ntraj, maxtraj
!     type(axis) :: grid, vgen

    real(8), allocatable :: grid(:), dgrid(:), volume(:)
    logical, allocatable :: fluxmask(:)
    type(boundary), allocatable :: bdry_arr(:)
    integer, allocatable :: emit_arr(:)
    
    real(8), allocatable :: trajectory(:,:)
    real(8), allocatable :: cumdisp(:)
    real(8), allocatable :: gridtime_arr(:,:), griddisp_arr(:,:)
    integer, allocatable :: gridloc_arr(:,:,:)
    integer, allocatable :: nstop(:)

contains

! type(axis) pure function makeaxis(dir, lo, hi, n) result(ax)
!     real(8), intent(in) :: dir(3), lo, hi
!     integer, intent(in) :: n
!
!     ax%dir = dir
!     ax%lo = lo
!     ax%hi = hi
!     ax%n = n
!     if (n /= 0) then
!         ax%h = (hi - lo)/n
!     else
!         ax%h = 0
!     end if
! end function makeaxis

! subroutine setgrid(t, nt, ax)
!     real(8), intent(in) :: t
!     integer, intent(in) :: nt
!     type(axis), intent(in) :: ax
!
!     tend = t
!     ntime = nt
!
!     grid = ax
!     ngrid = ax%n
! end subroutine setgrid

subroutine settime(t, nt)
    real(8), intent(in) :: t
    integer, intent(in) :: nt
    
    tend = t
    ntime = nt
end subroutine settime

subroutine setgrid_none(vol, gT)
    real(8), intent(in) :: vol, gT
    
!     gridflux = .false.
    ngrid = 0
    allocate( volume(0:0) )
    volume(0) = vol
    gradT = gT
end subroutine setgrid_none

subroutine setgrid_arr(dir, arr, vol, mask, gT)
    real(8), intent(in) :: dir(3), arr(0:), vol(:), gT
    logical, intent(in) :: mask(:)
    
!     gridflux = .true.
    ngrid = size(arr,1) - 1 !number of grid cells
    griddir = dir
    
    allocate( grid(0:ngrid) )
    allocate( dgrid(ngrid) )
    allocate( volume(ngrid) )
    allocate( fluxmask(ngrid) )
    grid = arr
    dgrid = arr(1:ngrid) - arr(0:ngrid-1)
    volume = vol
    fluxmask = mask
    gradT = gT
end subroutine setgrid_arr

subroutine setgrid_int(dir, lo, hi, num, vol, gT)
    real(8), intent(in) :: dir(3), lo, hi, vol, gT
    integer, intent(in) :: num
    real(8) :: arr(0:num), vol_arr(num)
    logical :: mask(num)
    integer :: i
    
    if (num == 0) then
        call setgrid_none(vol, gT)
    else if (lo >= hi) then
        print *, 'Error: setgrid: lo >= hi'
    else
        arr = (/(lo*(num-i)/num + hi*i/num, i=0,num)/)
        vol_arr = vol/num
        mask = .true.
        call setgrid_arr(dir, arr, vol_arr, mask, gT)
    end if
end subroutine setgrid_int

subroutine initrecord(dir)
!     logical, intent(in) :: gf
    real(8), intent(in) :: dir(3)
    
    if (ntime == 0) then
        allocate( gridtime_arr(0:nthreads-1, ngrid) )
        gridtime_arr = 0
    else
        tstep = tend/ntime
        
        allocate( gridloc_arr(0:nthreads-1, ngrid, 0:ntime) )
        gridloc_arr = 0
    end if
    
!     gridflux = gf
    flow = dir
    if (ngrid /= 0) then
        allocate( griddisp_arr(0:nthreads-1, ngrid) )
        griddisp_arr = 0
    else
        allocate( cumdisp(0:nthreads-1) )
        cumdisp = 0
    end if
    
    allocate( nstop(0:nthreads-1) )
    nstop = 0
end subroutine initrecord

subroutine inittraj(num, pos)
    integer, intent(in) :: num
    real(8), intent(in) :: pos(3)
    
    maxtraj = num
    ntraj = 0
    
    allocate( trajectory(3, 0:num) )
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

pure subroutine gettraj(num, traj)
    integer, intent(inout) :: num
    real(8), intent(inout) :: traj(3, 0:num)
    
    num = min(num, ntraj)
    traj(:, 0:num) = trajectory(:, 0:num)
end subroutine gettraj

real(8) pure function getcoord(pos) result(coord)
    real(8), intent(in) :: pos(3)
    
!     coord = (dot_product(pos, grid%dir) - grid%lo)/grid%h
    coord = dot_product(pos, griddir)
end function getcoord

integer pure function coordtoind(coord) result(ind)
    real(8), intent(in) :: coord
    
!     ind = min(ngrid, max(1, ceiling(coord)))
    ind = searchbin(grid(1:ngrid), coord)
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
    type(boundary) :: bdry_arr(size(bcs,1))
    integer :: i, n
    real(8) :: origin(3), ivec(3), jvec(3)
    
    n = size(bcs,1)
    do i = 1,n
        origin = verts(:,faces(2,i))
        ivec = verts(:,faces(1,i)) - origin
        jvec = verts(:,faces(3,i)) - origin
        bdry_arr(i) = makebdry(bcs(i), origin, ivec, jvec, temps(i), tris(i))
    end do
    
end function makebdry_arr

subroutine setbdry(arr)
    type(boundary), intent(in) :: arr(:)
!     real(8), intent(in) :: vol
    
    nbdry = size(arr)
    if (allocated(bdry_arr)) then
        deallocate( bdry_arr )
    end if
    allocate( bdry_arr(nbdry) )
    bdry_arr = arr
    
!     volume = vol
end subroutine setbdry

subroutine setbdrypair(ind1, ind2, mv)
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
    real(8) :: areatemp_arr(nbdry), sumareatemp
    integer :: i
    
    areatemp_arr = (/(bdry_arr(i)%area*abs(bdry_arr(i)%temp), i=1,nbdry)/)
    sumareatemp = sum(areatemp_arr)
    
    if (.not. allocated(emit_arr)) then
        allocate( emit_arr(nbdry) )
    else if (size(emit_arr, 1) /= nbdry) then
        deallocate( emit_arr )
        allocate( emit_arr(nbdry) )
    end if
    
    emit_arr = ceiling( num*areatemp_arr/sumareatemp )
    
    Eeff = (sumareatemp*tend/sum(emit_arr)) * getpseudoflux()
end subroutine calculateemit

pure function getemit() result(arr)
    integer :: arr(nbdry)
    
    arr = emit_arr
end function getemit

integer pure function getnemit() result(num)
    num = sum(emit_arr)
end function getnemit

real(8) pure function geteeff() result(E)
    E = Eeff
end function geteeff

! subroutine setvolumetric(vol, ax)
!     logical, intent(in) :: vol
!     type(axis), intent(in), optional :: ax
!
!     volumetric = vol
!     if (vol) then
!         vgen = ax
!     end if
! end subroutine setvolumetric

subroutine setvolumetric(dir, lim)
    real(8), intent(in) :: dir(3), lim(2)
    
    volumetric = .true.
    voldir = dir
    vollim = lim
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
    logical :: isdone
    real(8) :: r
    
    !$omp critical
    isdone = ( sum(emit_arr) == 0 )
    if (.not. isdone) then
        ind = maxloc(emit_arr, 1)
        emit_arr(ind) = emit_arr(ind) - 1
    end if
    !$omp end critical
    
    if (.not. isdone) then
        call drawposrect(pos)
        pos = bdry_arr(ind)%origin + pos(1)*bdry_arr(ind)%ivec + &
            pos(2)*bdry_arr(ind)%jvec
    
        if (volumetric) then
            call randnum(r)
            pos = pos - project(pos, voldir)
            pos = pos + ((1-r)*vollim(1) + r*vollim(2))*voldir
        end if
    
        call drawanghalf(dir)
        dir = matmul(bdry_arr(ind)%rot, dir)
    
        sign = (bdry_arr(ind)%temp >= 0)
    end if
end subroutine drawemitstate

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
    real(8) :: xi(3), diri(3)
    logical :: pierced(0:nbdry)
    real(8) :: colli(2), coll(3, 0:nbdry), dt(0:nbdry)
    integer :: i, indold
    integer :: threadi = 0
    
    dt(0) = min(deltat, tend - t)
    coll(:,0) = x + v*dt(0)*dir
    
    do i = 1,nbdry
        xi = matmul(bdry_arr(i)%inv, x - bdry_arr(i)%origin)
        diri = matmul(bdry_arr(i)%inv, dir)
        call getcoll(pierced(i), colli, xi, diri, bdry_arr(i)%tri)
        
        coll(:,i) = bdry_arr(i)%origin + &
            colli(1)*bdry_arr(i)%ivec + colli(2)*bdry_arr(i)%jvec
        dt(i) = normtwo(coll(:,i) - x)/v
    end do
    
    pierced(ind) = .false.
    pierced(0) = .true.
    indold = ind
    ind = minloc(dt(1:nbdry), 1, pierced(1:nbdry))
    if (dt(ind) >= dt(0)) then
        ind = 0
    end if
!     print('(1X,I3,7L2,7ES10.2)'), ind, pierced, dt1
    
    if (ind == 0) then
        bc = 0
    else
        bc = bdry_arr(ind)%bc
    end if
    
    xi = coll(:,ind)
!     if (abs(getcoord(xi) - getcoord(x)) < 0d0) then
!     if (x(3) < grid(0) .or. x(3) > grid(ngrid) .or. &
!         xi(3) < grid(0) .or. x(3) > grid(ngrid)) then
!         print *, 'Warning: updatestate: phonon did not move'
!         if (ntraj < 1000) then
!             print *, 'traj = '
!             print ('(3ES16.8)'), trajectory(:,0:ntraj), x
!         else
!             print *, 'xold = ', x
!             print *, 'xnew = ', xi
!         end if
!         print *, 'coord = ', getcoord(xi), getcoord(x)
!         print *, 'indold, ind  = ', indold, ind
!         call exit
!     end if
    x = xi
    x = coll(:,ind)
    
    deltat = dt(ind)
    
    t = t + deltat
    if (t >= tend) then
        bc = 0
        dir = 0
        
!$      threadi = omp_get_thread_num()
        nstop(threadi) = nstop(threadi) + 1
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

! subroutine incgrid(cum_arr, val, ind1, ind2)
!     real(8), intent(inout) :: cum_arr(0:,:)
!     real(8), intent(in) :: val
!     integer, intent(in) :: ind1
!     integer, intent(in), optional :: ind2
!     integer :: threadi = 0
!
! !$  threadi = omp_get_thread_num()
!     if (present(ind2)) then
!         if (ind1 <= ind2) then
!             cum_arr(threadi, ind1:ind2) = cum_arr(threadi, ind1:ind2) + val
!         end if
!     else
!         cum_arr(threadi, ind1) = cum_arr(threadi, ind1) + val
!     end if
! end subroutine incgrid

subroutine recordgrid(cum_arr, sign, delta, coordold, coordnew)
    real(8), intent(inout) :: cum_arr(0:nthreads-1, ngrid)
    logical, intent(in) :: sign
    real(8), intent(in) :: delta, coordold, coordnew
    real(8) :: coordlo, coordhi, dcoord, total, add(ngrid)
    integer :: pm, indlo, indhi
    integer :: threadi = 0
    
!$  threadi = omp_get_thread_num()
    
    pm = signtoint(sign)
    coordlo = max(min(coordold, coordnew), grid(0))
    coordhi = min(max(coordold, coordnew), grid(ngrid))
    dcoord = coordhi - coordlo
    if (dcoord < 0d0) then
        print *, 'Warning: recordgrid: dcoord = ', dcoord, ' < eps'
        print *, 'coordold, coordnew = ', coordold, coordnew
        print *, 'gridlo,   gridhi   = ', grid(0), grid(ngrid)
        print *, 'coordlo,  coordlo  = ', coordlo, coordhi
        if (ntraj < 1000) then
            print *, 'traj = '
            print ('(3ES16.8)'), trajectory(:,0:ntraj)
        end if
        call exit
    end if
    total = dcoord/abs(coordnew - coordold)*delta
    indlo = coordtoind(coordlo)
    indhi = coordtoind(coordhi)
!     print *, indlo, coordlo
!     print *, indhi, coordhi
    
    add = 0d0
    if (indlo == indhi) then
        add(indlo) = coordhi - coordlo
    else
        add(indlo:indhi) = dgrid(indlo:indhi)
        add(indlo) = grid(indlo) - coordlo
        add(indhi) = coordhi - grid(indhi-1)
    end if
!     if (abs(sum(add) - dcoord) > 1d-16) then
!         print *, 'Warning: recordgrid: sum(add) - dcoord = ', sum(add) - dcoord
!     end if
    add = pm*total*add/dcoord
    
    cum_arr(threadi,:) = cum_arr(threadi,:) + add
    
!     if (indold < indnew) then
!         call incgrid(cum_arr, pm*delta*(grid(indold) - coordold)/dcoord, indold)
!         call incgrid(cum_arr, pm*delta/dcoord, indold+1, indnew-1)
!         call incgrid(cum_arr, pm*delta*(1 - (indnew - coordnew))/dcoord, indnew)
!     else if (indnew < indold) then
!         call incgrid(cum_arr, pm*delta*(indnew - coordnew)/dcoord, indnew)
!         call incgrid(cum_arr, pm*delta/dcoord, indnew+1, indold-1)
!         call incgrid(cum_arr, pm*delta*(1 - (indold - coordold))/dcoord, indold)
!     else
!         call incgrid(cum_arr, pm*delta, indold)
!     end if
end subroutine recordgrid

subroutine recordtime(sign, deltat, xold, xnew)
    logical, intent(in) :: sign
    real(8), intent(in) :: deltat, xold(3), xnew(3)
    
    if (ntime == 0 .and. ngrid /= 0) then
        call recordgrid(gridtime_arr, sign, deltat, getcoord(xold), getcoord(xnew))
    end if
end subroutine recordtime

subroutine recorddisp(sign, xold, xnew)
    logical, intent(in) :: sign
    real(8), intent(in) :: xold(3), xnew(3)
    real(8) :: disp
    
    if (ngrid /= 0) then
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
    if (ntime /= 0 .and. ngrid /= 0) then
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
    if (ngrid == 0) then
        pm = signtoint(sign)
        cumdisp(threadi) = cumdisp(threadi) + pm*dot_product(flow, deltax)
    end if
end subroutine addcumdisp_real

subroutine addcumdisp_int(sign, ind)
    logical, intent(in) :: sign
    integer, intent(in) :: ind
    integer :: pm, threadi = 0
    
!$  threadi = omp_get_thread_num()
    if (ngrid == 0) then
        pm = signtoint(sign)
        cumdisp(threadi) = cumdisp(threadi) - pm*dot_product(flow, bdry_arr(ind)%pairmv)
    end if
end subroutine addcumdisp_int

function getsteadytemp() result(T)
    real(8) :: T(ngrid)
    
!     T = Eeff*ngrid/(volume*tend*getpseudoenergy()) * sum(gridtime_arr, 1)
    if (ntime == 0 .and. ngrid /= 0) then
        T = Eeff/(tend*getpseudoenergy()) * sum(gridtime_arr, 1)/volume
    else
        T = 0d0
    end if
end function getsteadytemp

pure function gettranstemp() result(T)
    real(8) :: T(ngrid, 0:ntime)
    real(8) :: sumgridloc(ngrid, 0:ntime)
    integer :: i

!     T = Eeff*ngrid/(volume*getpseudoenergy()) * sum(gridloc_arr, 1)
    if (ntime /= 0 .and. ngrid /= 0) then
        sumgridloc = sum(gridloc_arr, 1)
        do i = 0,ntime
            T(:,i) = Eeff/getpseudoenergy() * sumgridloc(:,i)/volume
        end do
    else
        T = 0d0
    end if
end function gettranstemp

real(8) pure function getflux() result(flux)
    
    if (ngrid == 0) then
        flux = Eeff/tend * sum(cumdisp)/volume(0)
    else 
        flux = Eeff/tend * sum(sum(griddisp_arr, 1), 1, fluxmask)/sum(volume, 1, fluxmask)
    end if
end function getflux

pure function getgridflux() result(flux)
    real(8) :: flux(ngrid)
    
    if (ngrid == 0) then
        flux = getflux()
    else
        flux = Eeff/tend * sum(griddisp_arr, 1)/volume
    end if
end function getgridflux

real(8) pure function getcond()
    getcond = -getflux()/gradT
end function getcond

integer pure function getnstop()
    getnstop = sum(nstop)
end function getnstop

end module domain
