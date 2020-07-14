subroutine inversion_metropolis(m,n,coe0,xobs,tobs,coe)
!-----------------------------------------------------------------------
! Inversao de parametros polinomiais com o Metodo Metropolis.
! Autor: Yan Carlos Viegas De Jesus
!
! Variaveis: m (= numero de medidas realizadas), n (= numero de
!            coeficientes do polinomio 2D), coe0 (= coeficientes do
!            modelo inicial), xobs (= valores dos espa‡amentos entre
!            fonte sismica e cada geofone), tobs (= valores observados
!            dos tempos de transito), coe (= coeficientes do modelo
!            invertido).
!-----------------------------------------------------------------------
dimension coe0(n), xobs(m), tobs(m), coe(n)    ! Variaveis principais.
allocatable tcalc0(:), tcalc(:)     ! Variaveis auxiliares.
allocate(tcalc0(m), tcalc(m))

call random_number(s)     ! Perturbacao aleatoria do modelo inicial.
z = 2.0*s - 1.0
do i = 1,n
  coe(i) = coe0(i) + (3.0/10.0)*z*coe0(i)
enddo
niter = 1000000      ! Numero de iteracoes.
k = 0; key = 1
do while (key .eq. 1)
!-----------------------------------------------------------------------
  call raytrace_time_inv(m,n,coe,xobs,tcalc)     ! Tracamento de raios.
!-----------------------------------------------------------------------
  error1 = dot_product(tobs - tcalc,tobs - tcalc)
  if (error1 .gt. 0.001) then
!-----------------------------------------------------------------------
    call raytrace_time_inv(m,n,coe0,xobs,tcalc0)
!-----------------------------------------------------------------------
    error2 = dot_product(tobs - tcalc0,tobs - tcalc0)
    if (error1 .lt. error2) then
      coe0(:) = coe(:)
      do j = 1,n
        errorm = (tobs(j) - tcalc(j))/sqrt((tobs(j))**2)
        coe(j) = coe0(j) + (3.0/10.0)*errorm*coe0(j)
      enddo
    else
      call random_number(a)
      w = (exp(-error1))/(exp(-error2))    ! Calculo das energias.
      if (a .lt. w) then
        coe0(:) = coe(:)
        do j = 1,n
          errorm = (tobs(j) - tcalc(j))/sqrt((tobs(j))**2)
          coe(j) = coe0(j) + (3.0/10.0)*errorm*coe0(j)
        enddo
      else
        do j = 1,n
          errorm = (tobs(j) - tcalc0(j))/sqrt((tobs(j))**2)
          coe(j) = coe0(j) + (3.0/10.0)*errorm*coe0(j)
        enddo
      endif
    endif
    key = 1
  else
    key = 0
  endif
  k = k + 1
  write(*,*) 'Number of iteractions =', k, '/', niter - k, 'left.'
  if (k .eq. niter) key = 0
enddo
deallocate(tcalc0,tcalc)
return
end

!-----------------------------------------------------------------------
