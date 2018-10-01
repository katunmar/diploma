PROGRAM SURFFLUX

implicit none


real, parameter :: cp = 1013., ro0 = 0.923, R = 287., z1 = 0.01, niu = 10.4e-6, g=9.814
real :: AR1(6)
real :: AR2(11)
integer :: it, n, m, i, j
real, dimension (10000,10) :: a

real :: TETSURF(10000) , TETAIR(10000) , HUMSURF, HUMAIR, P0(10000), P(10000)
real :: tau, u_star, zold, znew, cnt

OPEN (10, FILE='12', STATUS='old')
Read (10,*) n, m
do i=1,n
	read(10, *) (a(i,j),j=1,m)
	!print *, i
end do 

print *, "hello world!"

call PBLDAT

!*====================================================================
!*     .....DEFENITION OF DRAG AND HEAT EXCHANGE COEFFICIENTS......  =
!*     DETAILS OF ALGORITM ARE GIVEN IN:                             =
!*     A.L.KAZAKOV,V.N.LYKOSSOV,"TRUDY ZAP.SIB.NII",1982,N.55,3-20   =
!*                    INPUT DATA:
!*     AR1(1) - ABS(WIND VELOCITY) AT CONSTANT FLUX LAYER            =
!*                                 (CFL) HIGHT (M/S)                 =
!*     AR1(2) - DIFFERENCE BETWEEN POTENTIAL TEMPERATURE AT CFL HIGHT=
!*                             AND AT SURFACE  ( DEG. K)             =
!*     AR1(3) - SEMI-SUM OF POTENTIAL TEMPERATURE AT CFL HIGHT AND   =
!*                             AND AT SURFACE  ( DEG. K)             =
!*     AR1(4) - DIFFERENCE BETWEEN HUMIDITY AT CFL HIGHT             =
!*                             AND A SURFACE   ( GR/GR )             =
!*     AR1(5) - CFL HIGHT ( M )                                      =
!*     AR1(6) - ROUGHNESS OF SURFACE ( M ); FOR SEA SURFACE PUT -1   =
!*      IT    - NUMBER OF ITERATIONS                                 =
!*                     OUTPUT DATA:
!*     AR2(1) - NON-DIMENSIONAL CFL HIGHT                            =
!*     AR2(2) - RICHARDSON NUMBER                                    =
!*     AR2(3) - REYNODS NUMBER                                       =
!*     AR2(4) - LN(ZU/ZT)                                            =
!*     AR2(5) - DYNAMICAL ROUGHNESS ZU (M)                           =
!*     AR2(6) - THERMAL   ROUGHNESS ZT (M)                           =
!*     AR2(7) - CRITICAL RICHARDSON NUMBER                           =
!*     AR2(8) - TRANSFER COEFFICIENT FOR MOMENTUM                    =
!*     AR2(9) - TRANSFER COEFFICIENT FR HEAT                         =
!*     AR2(10)- COEFFICIENT OF TURBULENCE (KM) AT CFL HIGHT (M**2/S) =
!*     AR2(11)- ALFT=KT/KM ( KT-COEFFICIENT OF TURBULENCE FOR HEAT)  =
!*     COMMENT: DRAG COEFFICIENT =          AR2(8)*AR2(8)            =
!*              HEAT EXCHANGE COEFFICIENT = AR2(8)*AR2(9)            =
!*====================================================================


open(9, file='z0T_new_out_max.txt')

!zold=0.3
IT = 100
do i=1,n
    P0(i) = a(i,10)
    P(i) = a(i,9)
    
	TETSURF(i) = (a(i,3)+273)*(P0(i)/P(i))**(R/cp)
	TETAIR(i)= (a(i,2)+273)*(P0(i)/P(i))**(R/cp)

	HUMSURF = 1.e-4
	HUMAIR = 1.e-5

	AR1(1) = a(i,1)                          !AR1(1)
	AR1(2) = TETAIR(i) - TETSURF(i)          !AR1(2)
	AR1(3) = 0.5*(TETAIR(i) + TETSURF(i))    !AR1(3)
	AR1(4) = HUMAIR - HUMSURF     
	AR1(5) = 2.
	!AR1(6) = 0.01
	!call DRAG3(AR1,IT,AR2)
    !print*, sign(1.,AR1(2))
	
	zold = 0.
	znew = 0.0004
	cnt = 0
	do while(abs(zold-znew).GT.0.00001)
		!print *, zold - znew
		zold = znew
		AR1(6) = zold
		call DRAG3(AR1,IT,AR2)
		tau=ro0*AR2(8)*AR2(8)*AR1(1)*AR1(1)
		u_star=sqrt(tau/ro0)
		znew = ((0.135*niu)/u_star)+2.0*0.0001*exp(-((u_star-0.25)/0.15)**2.)+((0.03*u_star*u_star)/g)
		
		cnt = cnt + 1
		
	end do
		
	write(9,*) 'ZU = ', AR2(5), 'ZT = ', AR2(6), 'H = ',- cp*ro0*AR2(8)*AR2(9)*AR1(1)*AR1(2) 
	!'H = ',- cp*ro0*AR2(8)*AR2(9)*AR1(1)*AR1(2), 'tau =', ro0*AR2(8)*AR2(8)*AR1(1)*AR1(1), 'u* = ', u_star
	!'H = ',- cp*ro0*AR2(8)*AR2(9)*AR1(1)*AR1(2),'H(изм.) = ', a(i,7), AR2(8), AR2(9), 'ZU = ', AR2(5), 'ZT = ', AR2(6)
	!'z0 = ', znew,  'u* = ', u_star , cnt
	
	!'H = ',- cp*ro0*AR2(8)*AR2(9)*AR1(1)*AR1(2),'H(изм.) = ', a(i,7), 'u* = ', u, 'u*(изм.) = ', a(i,8), 'z0 = ', znew
enddo
!print*, 'LE = ', - L*ro0*AR2(8)*AR2(9)*AR1(1)*AR1(4)
!print*, 'tau = ', ro0*AR2(8)*AR2(8)*AR1(1)*AR1(1)


END PROGRAM SURFFLUX
