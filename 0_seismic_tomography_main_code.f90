!----------------------------------------------
! Metropolis Algorithm
! Geophysical Inversion | Seismic Tomography
! Author: Yan Viegas
!----------------------------------------------

allocatable tobs(:), xobs(:), coe0(:), coe(:)
m = 220; n = 28
allocate(tobs(m), xobs(m), coe0(n), coe(n))

open(5,file='model.dat')
open(10,file='traveltime1.dat')
open(20,file='traveltime2.dat')

do j = 1,m/2
  read(20,*) xobs(j), tobs(j)
enddo
do k = 1,m/2
  xobs((m/2)+k) = xobs(k)
  tobs((m/2)+k) = tobs(k)
enddo
do l = 1,m/2
  read(10,*) xobs(l), tobs(l)
enddo

do i = 1,n
  call random_number(s)
  z = 6.0*s - 3.0
  coe0(i) = z
enddo

call inv_metro(m,n,coe0,xobs,tobs,coe)

write(*,*) 'Parameters of the model:'
write(5,*) 'Parameters of the model:'
do i = 1,n
  write(*,*) 'coe(', i,') =', coe(i)
  write(5,*) 'p(', i,') =', coe(i)
enddo

close(5); close(10); close(20)

deallocate(tobs,xobs,coe0,coe)
stop
end

!----------------------------------------------------------------------
