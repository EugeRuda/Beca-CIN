Module Functions
  Use iso_fortran_env, only:wp => real64
  Use Parameters

Contains

  !==============================================================================
  !Funciones para el ejercicio 1 y 2
  !==============================================================================
  !Tiempo de enfriamiento
  !----------------------
  
  Real(wp) Function Tcool(E,dEdt)
    Use iso_fortran_env, only: wp => real64
    Use Parameters
    Implicit None
    Real(wp), intent(in) :: E
    Real(wp), intent(in) :: dEdt
    
    Tcool = -(1.d0/dEdt) * E    !Donde b = dE/dt.Las próximas funciones dan el b para cada proceso.
    
    Return
  End Function Tcool
  
  !Radiacion sincrotron
  !--------------------
  
  Real(wp) Function b_sinc(E)
    Use Parameters
    Use iso_fortran_env, only: wp => real64
    Implicit None
    Real(wp), intent(in) :: E
    
    b_sinc = -b0sinc * E**2
    
    Return
  End Function b_sinc
  
  !Radiacion Inverse Compton
  !-------------------------
  
  Real(wp) Function b_ic(E)
    Use iso_fortran_env, only: wp => real64
    Use Parameters
    Implicit None
    Real(wp), intent(in) :: E
    Real(wp) :: gamma
    
    gamma = E/Eer
    
    b_ic =-Eer * (5.5D17) * (Ttilde**3) * gamma * &
         &(log(1.D0 + 0.55D0 * gamma * Ttilde) / (1.D0 + 25.D0 * Ttilde * gamma))&
         &* (1.D0 + (1.4D0 * gamma * Ttilde) / (1.D0 + 12.D0 * gamma**2 * Ttilde**2)) &
         &* (Ucmb/Ubb)
    
    Return
  End Function b_ic
  
  !Bremsstrahlung Relativista (con apantallamiento)
  !------------------------------------------------
  
  Real(wp) Function b_br(E)
    Use iso_fortran_env, only: wp => real64
    Use Parameters
    Implicit None
    Real(wp), intent(in) :: E
    
    b_br = b0br * E
    
    Return
  End Function b_br
  
  !Proton - proton
  !---------------
  
  Real(wp) Function b_pp(E)
    Use iso_fortran_env, only: wp => real64
    Use Parameters
    Implicit None
    Real(wp),intent(in) :: E
    Real(wp) :: L, sigma
    
    If (E>=Eth) then
       L = log(E / (1.D12*eV))
       sigma = (34.3D0 + 1.88D0*L + 0.25D0*L**2) * (1.D0 - (Eth/E)**4.D0)**2 * mb
       b_pp = -c * na * Kpp * sigma * E
    Else
       b_pp = 0.D0
    End If
    
    Return
  End Function b_pp
  
  !Tiempo de aceleracion
  !---------------------
  
  Real(wp) Function Tacel(E)
    Use iso_fortran_env, only: wp => real64
    Use Parameters
    Implicit None
    Real(wp), intent(in) :: E
    
    Tacel = E / (eta * qe * c * B)
    
    Return
  End Function Tacel
  
  !===============================================================================

  real(wp) function ne(E,A_e,Ee_max)
    
    use iso_fortran_env , ONLY: wp=>REAL64
    use Parameters
    implicit none
    
    real(wp),intent(in) :: E         !E va tomando los valores de Ee(i)
    real(WP),intent(in) :: Ee_max,A_e 
    
    
    ne=A_e*(E**(-p))*exp(-1.0*(E/Ee_max))   
    
  end function ne
  
  REAL(wp) FUNCTION  Np(E,A_p,Ep_max)
    USE iso_fortran_env, ONLY: wp => REAL64
    USE Parameters
    IMPLICIT NONE
    
    Real(WP), INTENT(IN) :: E,Ep_max,A_p     ! Energia de las particulas
    
    Np=A_p*(E**(-p))*exp(-1.0*(E/Ep_max)) ! Tomamos Emax=Epr+Epimax/0.17
    
  END function NP
  
End Module Functions
