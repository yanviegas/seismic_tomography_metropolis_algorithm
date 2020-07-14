subroutine raytrace_time_inv(m,n,coe,xobs,tcalc)
!-----------------------------------------------------------------------
! Tracamento de raios para determinar os tempos de transito na Inversao.
!
! Autores: Yan Carlos Viegas De Jesus & Dr. Wilson Mouzer Figueir¢
!
! Variaveis: Variaveis: m (= numero de medidas realizadas), n (= numero de
!            coeficientes do polinomio 2D), coe0 (= coeficientes do
!            modelo atual), xobs (= valores dos espa‡amentos entre
!            fonte sismica e cada geofone), tobs (= valores calculados
!            dos tempos de transito).
!-----------------------------------------------------------------------
dimension coe(n), xobs(m), tcalc(m)
allocatable xcalc1(:), tcalc1(:)
l1 = 360
allocate(xcalc1(l1), tcalc1(l1))
!-----------------------------------------------------------------------
do l = 1,l1
  THETA = (3.1415926536/360.0)*l
  DTAU = 0.012500
  S = 16.0
  V = coe(1) + coe(2)*S + coe(4)*(S**2) + coe(7)*(S**3) + coe(11)*(S**4) &
  & + coe(16)*(S**5) + coe(22)*(S**6)
  P0 = 1/V
  PX0 = P0*COS(THETA)
  PZ0 = P0*SIN(THETA)
  X = S
  Z = 0.0
  T = 0.0
  key = 1
  do while (key .eq. 1)
      UMGRADX = coe(2) + 2*coe(4)*X + coe(5)*Z + 3*coe(7)*(X**2) &
    &           + 2*coe(8)*X*Z + coe(9)*(Z**2) + 4*coe(11)*(X**3) &
    &           + 3*coe(12)*(X**2)*Z + 2*coe(13)*X*(Z**2) &
    &           + coe(14)*(Z**3) + 5*coe(16)*(X**4) + 4*coe(17)*(X**3)*Z &
    &           + 3*coe(18)*(X**2)*(Z**2) + 2*coe(19)*X*(Z**3) &
    &           + coe(20)*(Z**4) + 6*coe(22)*(X**5) + 5*coe(23)*(X**4)*Z &
    &           + 4*coe(24)*(X**3)*(Z**2) + 3*coe(25)*(X**2)*(Z**3) &
    &           + 2*coe(26)*X*(Z**4) + coe(27)*(Z**5)
      UMGRADZ = coe(3) + coe(5)*X + 2*coe(6)*Z + coe(8)*(X**2) &
    &           + 2*coe(9)*X*Z + 3*coe(10)*(Z**2) + coe(12)*(X**3) &
    &           + 2*coe(13)*(X**2)*Z + 3*coe(14)*X*(Z**2) &
    &           + 4*coe(15)*(Z**3) + coe(17)*(X**4) + 2*coe(18)*(X**3)*Z &
    &           + 3*coe(19)*(X**2)*(Z**2) + 4*coe(20)*X*(Z**3) &
    &           + 5*coe(21)*(Z**4) + coe(23)*(X**5) + 2*coe(24)*(X**4)*Z &
    &           + 3*coe(25)*(X**3)*(Z**2) + 4*coe(26)*(X**2)*(Z**3) &
    &           + 5*coe(27)*X*(Z**4) + 6*coe(28)*(Z**5)
      GRADX = -UMGRADX/(V**3)
      GRADZ = -UMGRADZ/(V**3)
      XANT = X
      ZANT = Z
      X = X + PX0*DTAU
      Z = Z + PZ0*DTAU
      IF ((X.GT.32.0).OR.(X.LT.0.0)) THEN
        IF (Z.LE.0.0) THEN
          tcalc1(l) = T
          xcalc1(l) = XANT
          key = 0
        ENDIF
      ENDIF
      IF ((Z.GT.4.0).OR.(Z.LT.0.0)) THEN
        IF (Z.LE.0.0) THEN
          tcalc1(l) = T
          xcalc1(l) = XANT
          key = 0
        ENDIF
        IF (Z.GE.4.0) THEN
          key = 0
        ENDIF
      ENDIF
      D = SQRT((X-XANT)**2 + (Z-ZANT)**2)
      T = T + D/V
      PX0 = PX0 + GRADX*DTAU
      PZ0 = PZ0 + GRADZ*DTAU
      V = coe(1) + coe(2)*X + coe(3)*Z + coe(4)*(X**2) + coe(5)*X*Z &
    &   + coe(6)*(Z**2) + coe(7)*(X**3) + coe(8)*(X**2)*Z &
    &   + coe(9)*X*(Z**2) + coe(10)*(Z**3) + coe(11)*(X**4) &
    &   + coe(12)*(X**3)*Z + coe(13)*(X**2)*(Z**2) + coe(14)*X*(Z**3) &
    &   + coe(15)*(Z**4) + coe(16)*(X**5) + coe(17)*(X**4)*Z &
    &   + coe(18)*(X**3)*(Z**2) + coe(19)*(X**2)*(Z**3) &
    &   + coe(20)*X*(Z**4) + coe(21)*(Z**5) + coe(22)*(X**6) &
    &   + coe(23)*(X**5)*Z + coe(24)*(X**4)*(Z**2) &
    &   + coe(25)*(X**3)*(Z**3) + coe(26)*(X**2)*(Z**4) &
    &   + coe(27)*X*(Z**5) + coe(28)*(Z**6)
      PXN = V*SQRT(PX0**2 + PZ0**2)
      PX0 = PX0/PXN
      PZ0 = PZ0/PXN
  enddo
enddo

do j = 1,m
  xiobs = xobs(j)
  call int_shep(l1,xcalc1,xiobs,tcalc1,ticalc)
  tcalc(j) = ticalc
enddo

deallocate(xcalc1, tcalc1)
return
end

!-----------------------------------------------------------------------
