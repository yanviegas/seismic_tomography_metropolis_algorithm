subroutine int_shep(l1,xcalc1,xiobs,tcalc1,ticalc)
!-----------------------------------------------------------------------
! Interpolacao de valores discretos usando o algoritmo de Shepard.
!-----------------------------------------------------------------------
dimension xcalc1(l1), tcalc1(l1)
allocatable w(:)
allocate(w(l1))

s = 0.0
do i = 1,l1
  w1 = (xiobs - xcalc1(i))**2
  w(i) = 1.0/w1
  s = s + w(i)
enddo

f = 0.0
do j = 1,l1
  f = f + tcalc1(j)*w(j)
enddo

ticalc = f/s

deallocate(w)
return
end
