Module Distributions
  Use iso_fortran_env, only: wp => real64
  Use Parameters

Contains

Subroutine Proportionality(Ee_max,Ep_max,A_e,A_p)
     USE iso_fortran_env, ONLY: wp => REAL64
     USE Parameters
!!$     ---------------------------------------------------------------
  IMPLICIT NONE
  
!!$     ---------------------------------------------------------------
!!$     Programa que calcula las constantes de proporcionalidad Ae y Ap
!!$     para las distribuciones ne(E) y np(E)
!!$     ---------------------------------------------------------------

!!$-----------------------------------------------------------------------  
!!$    Declaracion de variables
!!$-----------------------------------------------------------------------  

  INTEGER :: k,i                  ! Contador para los do's
  INTEGER, parameter :: nmax=5000 ! Cantidad de puntos de integracion

  REAL(WP), allocatable :: Ee(:)  ! Energias de los electrones
  REAL(WP) :: deltaEe             ! Paso  de los electrones

  REAL(WP), allocatable :: Ep(:)  ! Energias de los protones
  REAL(WP) :: deltaEp             ! Paso  de los protones

  REAL(WP) :: integralp=0.d0      ! Integral para los protones
  REAL(WP) :: integrale=0.d0      ! Integral para los electrones
  Real(WP) :: w=0_WP

  Real(wp), intent(IN) :: Ee_max, Ep_max

  Real(WP), intent(out) :: A_p                     ! Constante de proporcionalidad para los protones
  Real(WP), intent(out) :: A_e                     ! Constante de proporcionalidad para los electrones

  REAL(WP), external :: f

  
!!$-----------------------------------------------------------------------  
!!$    Calculamos el paso para las energias de electrones y protones
!!$-----------------------------------------------------------------------


  deltaEe=(Ee_max/Eemin)**(1.0/(real(nmax)-1.0))
  Allocate(Ee(nmax))  

  Ee(1)=Eemin
  
  Do k=2,nmax
     Ee(k)=Ee(k-1)*deltaEe
  EndDo
 
  deltaEp=(Ep_max/Epmin)**(1.0/(real(nmax)-1.0))
  Allocate(Ep(nmax))  

  Ep(1)=Epmin
  
  Do k=2,nmax
     Ep(k)=Ep(k-1)*deltaEp
  EndDo
  
!!$-----------------------------------------------------------------------  
!!$    Calculamos las integrales. Usamos la regla del rectangulo
!!$-----------------------------------------------------------------------
  
     Do i=1,nmax-1                            ! Barremos las energias de los elect.
        w=Ee(i+1)-Ee(i)                       ! Calculamos la base del rectangulo
        integrale=integrale+f(Ee(i),Ee_max)*w  ! Calculamos la altura como f(Ek,Ei)
     EndDo

     w=0.d0

     Do i=1,nmax-1                            ! Barremos las energias de los elect.
        w=Ep(i+1)-Ep(i)                       ! Calculamos la base del rectangulo
        integralp=integralp+f(Ep(i),Ep_max)*w  ! Calculamos la altura como f(Ek,Ei)
     EndDo

!     integralp=20.4938
!     integrale=123.225
     
!!$-----------------------------------------------------------------------  
!!$    Calculamos las constantes de proporcionalidad
!!$-----------------------------------------------------------------------

     A_p=(Xrel*Esnr)/((Kep+1.0)*integralp)

     A_e=(Xrel*Esnr)/((Kep+1.0/Kep)*integrale)

!!$-----------------------------------------------------------------------  
!!$    Escribimos los resultados
!!$-----------------------------------------------------------------------
!     Write(20,*) "p= ", p
!     Write(20,*) "Epmax= ", Epmax
!     Write(20,*) "Epmin= ", Epmin
!     Write(20,*) "Eemax= ", Eemax
!     Write(20,*) "Eemin= ", Eemin
!     Write(20,*)   
!     Write(20,*) "integralp= ", integralp
!    Write(20,*) "integrale= ", integrale
!     Write(20,*)
!     Write(20,*) "Ap= ", Ap
!     Write(20,*) "Ae= ", Ae

   End Subroutine Proportionality

 End Module Distributions
 
   
!!$-----------------------------------------------------------------------
!!$ Funcion a integrar
!!$-----------------------------------------------------------------------
REAL(wp) FUNCTION  f(E,Emax)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  IMPLICIT NONE
  
  Real(WP), INTENT(IN) :: E    ! Energia de las particulas
  Real(WP), INTENT(IN) :: Emax ! Energia maxima de las particulas
  
  f=E**(-p+1.0)*exp(-1.0d0*(E/Emax))
  
END function f
!!$-----------------------------------------------------------------------
  
