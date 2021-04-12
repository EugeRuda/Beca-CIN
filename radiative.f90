Module Radiative
  Use iso_fortran_env, only: wp => real64
  Use Parameters
  Use Functions

Contains       !Secciones eficaces + potencia emitida por distribuciones para cada proceso físico

  !-------------------------------------------------------------------------------
  !Bremsstrahlung Relativista
  !-------------------------------------------------------------------------------

  REAL(wp) FUNCTION  bremss(Eph,Ee)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE parameters
  Use Functions
  IMPLICIT NONE

  Real(WP), INTENT(IN) :: Ee  ! Energia de los electrones
  Real(WP), INTENT(IN) :: Eph ! Energia de los fotones

  Real(WP) :: x,Phi, SecEdif    !Seccion eficaz diferencial

  !Caso apantallamiento

  x = 1. - (Eph/Ee)
 
  phi=1. + x**2 - (2./3.)*x
    
  phi = (phi*log(191./(Z**(1./3.)))) + (1./9.)*x
   
  SecEdif = (4.*alfa*(re**2)*(Z**2))*Phi/Eph 

  bremss = SecEdif*c*Eph*na 

END function bremss

REAL(wp) FUNCTION  fbremss(Eph,Ee,A_e,Ee_max)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE parameters
  Use Functions
  IMPLICIT NONE
  
  Real(WP), INTENT(IN) :: Ee,Ee_max,A_e  ! Energia de los electrones
  Real(WP), INTENT(IN) :: Eph ! Energia de los fotones
  
  fbremss = ne(Ee,A_e,Ee_max)*bremss(Eph,Ee)

END function fbremss


!-------------------------------------------------------------------------------
!Inverse Compton
!-------------------------------------------------------------------------------

REAL(wp) FUNCTION  SigmaIC(Ee,Eph,Egamma)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE

  Real(WP), INTENT(IN) :: Ee     ! Energia de las particulas
  Real(WP), INTENT(IN) :: Eph    ! Energia de los fotones interactuantes (baja energia)
  Real(WP), INTENT(IN) :: Egamma ! Energia de los fotones gamma resultantes

  Real(WP) :: e_ph, e_gamma, gamma, x, f ! e minuscula = epsilon

  e_ph = Eph/Eer
  e_gamma = Egamma/Eer
  gamma = Ee/Eer
  ! Eer es la energia en reposo del electron, mec^2, y esta definida en Parameters
 
  x = e_gamma/(4._WP*e_ph*(gamma**2._wp)*(1._wp-(e_gamma/gamma)))
    
  If ((1._WP >= x).AND.(x >= (1._WP/(4._WP*gamma**2._WP)))) then

       f = (4._wp*e_ph*gamma*x)**2._wp
       f = f*(1._wp - x)
       f = f/(2._wp*(1._wp + 4._wp*e_ph*gamma*x))
       f = f + 2._wp*x*log(x) + x + 1._wp - 2._wp*x**2._wp

  Else

       f = 0._WP

  End If
         
!sigmath (en cm2)

  SigmaIC = 3._wp*sigmath*f/(4._wp*e_ph*(gamma**2._wp))

END function SigmaIC

REAL(wp) FUNCTION  fic(Ee,Eph,Egamma,A_e,Ee_max)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE
  
  Real(WP), INTENT(IN) :: Ee,Ee_max,A_e  ! Energia de los electrones
  Real(WP), INTENT(IN) :: Eph ! Energia de los fotones interactuantes
  Real(WP), INTENT(IN) :: Egamma ! Energia de los fotones resultantes
  
  fic = (Egamma/Eer)*ne(Ee,A_e,Ee_max)*SigmaIC(Ee,Eph,Egamma)*c*(Ucmb/Eph) ! n_cmb = Ucmb/Ecmb , n = densidad de particulas [cm^-3] , U = densidad de energia [erg/cm^3] , E = energia de las particulas [erg] 

END function fic


!-------------------------------------------------------------------------------
!Proton-Proton
!-------------------------------------------------------------------------------

REAL(wp) FUNCTION  sigma_PP(E)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE

  Real(WP), INTENT(IN) :: E  ! Energia
  Real(WP) :: L1
  Real(WP) :: E_th=280*Mev   ! Umbral de energia cinetica de los protones
                             ! para que se produzcan las reacciones

  L1 = log(E/TeV)                                         ! L = ln(Ep/1TeV)
  sigma_PP = (34.3d0 + 1.88d0*L1 + 0.25d0*L1**2)*1.d-27   !cm²
  
  ! Agregamos un factor para que ajuste mejor a bajas energías
  sigma_PP = sigma_PP*(1.d0-(E_th/E)**4.)**2.
  
END function SIGMA_PP

REAL(wp) FUNCTION  emisividadpp(Epi,A_p,Ep_max)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE

  Real(WP), INTENT(IN) :: Epi,Ep_max,A_p ! Energia de los piones
  
  ! Calculamos la emisividad
  emisividadpp = c*na*sigma_PP(Epr+(Epi/0.17_wp))*Np(Epr+(Epi/0.17_wp),A_p,Ep_max)/(0.17_wp)

END function emisividadpp

REAL(wp) FUNCTION  fpp(Epi,A_p,Ep_max)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE
  
  Real(WP), INTENT(IN) :: Epi,Ep_max,A_p  ! Energia de los piones
  
  fpp=(2.0d0*emisividadpp(Epi,A_p,Ep_max))/(sqrt((Epi**2)-(Epir**2)))
  
END function fpp

!-------------------------------------------------------------------------------
!Sincrotron
!-------------------------------------------------------------------------------

REAL(wp) FUNCTION  Potsyn(Ee,Eph)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE

  Real(WP), INTENT(IN) :: Ee  ! Energia de los electrones
  Real(WP), INTENT(IN) :: Eph ! Energia de los fotones

  Real(WP) :: gamma            ! Factor de Lorentz de los electrones
  Real(WP) :: E_c              ! Energia critica de los electrones
  Real(WP) :: x
  
  gamma=Ee/Eer

  E_c=(0.75/pi)*((qe*h*B)/(me*c))*(gamma**2)

  x=Eph/E_c
  
!  Write(*,*) x
  
  Potsyn=(sqrt(3.0)*(qe**3))/(h*Eer)*B*1.85*(x**(1.0/3.0))*exp(-1.0*x)

END function POTsyn

REAL(wp) FUNCTION  fsyn(Ee,Eph,A_e,Ee_max)
  USE iso_fortran_env, ONLY: wp => REAL64
  USE Parameters
  Use Functions
  IMPLICIT NONE
  
  Real(WP), INTENT(IN) :: Ee,Ee_max,A_e  ! Energia de los electrones
  Real(WP), INTENT(IN) :: Eph ! Energia de los fotones
  
  fsyn=Potsyn(Ee,Eph)*Ne(Ee,A_e,Ee_max)
  
END function fsyn

End Module Radiative
