!  1D particle in cell code by Y.Lv
!  2/18/2016

program pic

!-------------------------------------------------------------------------------

use poisson1d

implicit none

! Read-in parameters: particle loading, k

integer		::	npg,k
real, parameter	::	pi = 3.1415926535898
real, parameter ::	L = 2.*pi ! box size
integer, parameter	::	ng = 64 ! number of mesh zones
real, parameter	::	q1 = -1., q2 = -1.
real, parameter	::	m1 = 1., m2 = 1.
real, parameter	::	v1 = 4., v2 = -4. 	

! Others

integer		::	np
real		::	omegap, tmax
real		::	dt, dxg, dxp 
real, allocatable :: pos1(:), pos2(:), vel1(:), vel2(:), dv1(:), dv2(:)
real		::  rho(0:ng-1), phi(0:ng-1), mesh(0:ng-1)
character(len=30)  :: deb, deb2, grow, p4t, p2t, pt, f4t, f2t, ft
character(len=30)  ::  arg1, arg2

integer :: i, icell, ierr, iter=0 
integer :: ileft, iright, inext
real :: t=0, dist, accl, vmax, dvmax
real :: smallt=0.1, rho_b

!-------------------------------------------------------------------------------

! Assigning values

call getarg(1,arg1)
read(arg1,*) npg 
call getarg(2,arg2)
read(arg2,*) k

np = npg * ng	! number of particles per stream
omegap = sqrt(2.*np/L)	! plasma frequency
tmax = 16.*pi/omegap	! maximum simulation time
dxg = L/ng	! mesh spacing
dt = dxg/4.	! initial time step interval
dxp = L/np	! initial zone spacing
rho_b = -np*(q1+q2)/L

f4t = 'out4t.dat'
f2t = 'out2t.dat'
ft = 'outt.dat'
p4t = 'part4t.dat'
p2t = 'part2t.dat'
pt = 'partt.dat'
deb = 'debug.dat'
deb2 = 'debugf.dat'
grow = 'growth.dat'

! Initialize particles

allocate(pos1(0:np-1), stat=ierr)
allocate(pos2(0:np-1), stat=ierr)
allocate(vel1(0:np-1), stat=ierr)
allocate(vel2(0:np-1), stat=ierr)
allocate(dv1(0:np-1), stat=ierr)
allocate(dv2(0:np-1), stat=ierr)

vel1(:) = v1
vel2(:) = v2

do i = 0, np-1
  pos1(i) = (i+1./2.)*dxp
  pos1(i) = pos1(i)+ 0.01*L*sin(2.*pi*k*pos1(i)/L)
  pos2(i) = (i+1./2.)*dxp
  pos2(i) = pos2(i)- 0.01*L*sin(2.*pi*k*pos2(i)/L)
enddo

! Initialize mesh
do i = 0, ng-1
   mesh(i) = (i+1./2.)*dxg
end do
print *, "Done initialization. Move now..." 

! Main loop
do while ( t .le. tmax)

    ! Debugging/movie
    open(9, file=deb, form='formatted', position='append')
    do i = 0, np-1
        write(9,'(4F10.3)') pos1(i), vel1(i), pos2(i), vel2(i)
    enddo
    write(9,*)
    close(9)

    rho(:) = rho_b
    phi(:) = 0
    do i = 0, np-1 ! Loop over particles

	! First advancement of leapfrog
 	pos1(i) = pos1(i) + 1./2.*vel1(i)*dt
        pos2(i) = pos2(i) + 1./2.*vel2(i)*dt

        ! Check if any particle leaves the box
        if(pos1(i)>L) pos1(i) = pos1(i)-L
        if(pos1(i)<0) pos1(i) = L+pos1(i)
	if(pos2(i)>L) pos2(i) = pos2(i)-L
        if(pos2(i)<0) pos2(i) = L+pos2(i)
    
	! Calculating charge density using CIC
  	icell = floor(pos1(i)/dxg)
  	dist = mesh(icell)-pos1(i)
  	if (dist > 0) then
	    if (icell .ne. 0) then
		inext = icell -1
	    else
		inext = ng -1
	    end if 
	    rho(inext) = rho(inext)+q1*dist/dxg**2.
	    rho(icell) = rho(icell)+q1*(dxg-dist)/dxg**2.
  	else 	!dist<0
	    if (icell .ne. ng-1) then
		inext = icell+1
	    else
		inext = 0
	    end if 
	    rho(inext) = rho(inext)-q1*dist/dxg**2.
	    rho(icell) = rho(icell)+q1*(dxg+dist)/dxg**2.
  	end if

	icell = floor(pos2(i)/dxg)
        dist = mesh(icell)-pos2(i)
        if (dist > 0) then
            if (icell .ne. 0) then
                inext =	icell -1
            else
                inext =	ng -1
            end	if
            rho(inext) = rho(inext)+q2*dist/dxg**2.
       	    rho(icell) = rho(icell)+q2*(dxg-dist)/dxg**2.
        else    !dist<0
            if (icell .ne. ng-1) then
                inext =	icell+1
            else
                inext =	0 
            end	if
            rho(inext) = rho(inext)-q2*dist/dxg**2.
            rho(icell) = rho(icell)+q2*(dxg+dist)/dxg**2.
        end if

    end do
!    rho(0) = rho(1)
!    rho(ng-1)= rho(ng-2)

    ! Solve Poisson equation
    if(iter == 0) then
        call solve_poisson(rho,phi,ng,dxg,1)
    else
        call solve_poisson(rho,phi,ng,dxg,0)
    end if

    ! Debugging again
    open(30, file=deb2, form='formatted', position='append')
    do i = 0, ng-1
        write(30,'(3F10.3)') mesh(i),rho(i), phi(i)
    enddo
    write(30,*)
    close(30)

    ! Second advancement of leapfrog: update velocity and position

    do i = 0, np-1
        icell = floor(pos1(i)/dxg)
	if(icell == 0 ) then
	    ileft = ng-1
	else
	    ileft = icell-1
	end if
	if(icell == ng-1 ) then
           iright = 0
        else
            iright = icell+1
	end if
!	if(icell == 0) then
!	    ileft = icell
!	    iright = icell+2
!	else if(icell == ng-1) then
!	    ileft = icell-2
!	    iright = icell
!	else
!	    ileft = icell-1
!	    iright = icell+1
!	end if 
        accl = -q1/m1*(phi(iright)-phi(ileft))/2./dxg

	vel1(i) = vel1(i)+accl*dt
	pos1(i) = pos1(i)+1./2.*vel1(i)*dt

	! Do the same thing for the second stream
        icell = floor(pos2(i)/dxg)
	if(icell == 0) then
            ileft = ng-1
        else
            ileft = icell-1
        end if
	if(icell == ng-1) then
            iright = 0
        else
            iright = icell+1
        end if
!	if(icell == 0) then
!            ileft = icell
!            iright = icell+2
!        else if(icell == ng-1) then
!            ileft = icell-2
!            iright = icell
!        else
!            ileft = icell-1
!            iright = icell+1
!        end if

        accl = -q2/m2*(phi(iright)-phi(ileft))/2./dxg

        vel2(i)	= vel2(i) + accl*dt
        pos2(i) = pos2(i) + 1./2.*vel2(i)*dt

	dv1(i) = abs(vel1(i)-v1)
	dv2(i) = abs(vel2(i)-v2)
    end do

    ! Output

    if (abs(t-tmax/4) < smallt) then 
	call outputpart(p4t,vel1, pos1, vel2, pos2, np)
	call outputfield(f4t, mesh, rho, phi, ng)
    else if (abs(t-tmax/2) < smallt) then
        call outputpart(p2t,vel1, pos1, vel2, pos2, np)
        call outputfield(f2t, mesh, rho, phi, ng)
    else if (abs(t-tmax) < smallt) then
        call outputpart(pt,vel1, pos1, vel2, pos2, np)
        call outputfield(ft, mesh, rho, phi, ng)
    end if

    iter = iter + 1
    vmax = max(maxval(abs(vel1)),maxval(abs(vel2)))
    dvmax = max(maxval(dv1),maxval(dv2))
    dt = dxg/vmax
    t = t + dt
  
    open(40, file=grow, form='formatted', position='append')
    write(40,'(2F10.3)') t, dvmax

end do

close(40)

print *, "reached maximum simulation time."

end program pic

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine outputpart(fname, vel1, pos1, vel2, pos2, np)

character(len=30)  :: fname
integer :: np
real	:: vel1(0:np-1), pos1(0:np-1), vel2(0:np-1), pos2(0:np-1) 
integer :: i

open(10, file=fname, form='formatted', position='append')

do i = 0, np-1
  write(10,'(4F10.3)') pos1(i), vel1(i), pos2(i), vel2(i)
enddo
write(10,*)

close(10)

end subroutine

!-------------------------------------------------------------------------------
subroutine outputfield(fname, mesh, rho, phi, ng)

character(len=30)  :: fname
integer :: ng
real    :: mesh(0:ng-1), rho(0:ng-1), phi(0:ng-1)
integer	:: i

open(20, file=fname, form='formatted', position='append')

do i = 0, ng-1
  write(20,'(3F10.3)') mesh(i), rho(i), phi(i)
enddo
write(20,*)

close(20)

end subroutine

