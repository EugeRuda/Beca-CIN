Program TrabajoFinal

  !===============================================================================
  !TRABAJO FINAL DE INTRODUCCION A LA ASTROFISICA RELATIVISTA - 2018
  !===============================================================================
  !En este trabajo estudiaremos el Remanente de Supernova RX J1713.7-3946
  !-------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  !Declaracion de modulos a utilizar
  !-------------------------------------------------------------------------------
  
  Use iso_fortran_env, only: wp => real64
  Use Parameters
  Use Radiative
  Use Functions
  Use Distributions

  !-------------------------------------------------------------------------------
  !Declaracion de variables
  !-------------------------------------------------------------------------------
  
  Implicit None

  Real(wp) :: Ee, Eemax, Ae                          !Energia de los electrones
  Real(wp) :: Ep, Epmax, Ap                          !Energia de los protones
  Real(wp) :: Eph, Ephmax_e, Ephmax_p                !Energia de los fotones
  Real(wp) :: Epi, Epimin, Epimax                    !Energia de los piones
  Real(wp) :: ke, kp, kph_e, kph_p , kpi             !Pasos de integracion
  Real(wp) :: base                                   !Base de los rectangulos para integrar

  Integer :: i, k                          !contadores en los ciclos do
  Integer, parameter :: nmax = 1000        !cantidad de puntos de integracion

  Real(wp) :: Lsy, Lic, Lbremss, Lpp       !Luminosidades 

  !-------------------------------------------------------------------------------
  !Preparamos los archivos
  !-------------------------------------------------------------------------------

  Open(10, File='Tiempos-Cool-Acel.dat')
  Write(10,*)  '# Ee                        Ep                        Tcool_sinc  &
         &              Tcool_bremss              Tcool_ic                  Tcool_pp &
         &                 Tacel_e                   Tacel_p                   Tcool_TOT'

  Open(20, File='Bremss.dat')
  Write(20,*) '# log(E[eV])              log(E*L(E)[erg/s])   '
  Open(30, File='IC.dat')
  Write(30, *) '#    log(E[eV])          log(E*L(E)[erg/s])'
  Open(40, File="PP.dat")
  Write(40,*) "#  log(E[eV])                 log(E*L(E)[erg/s])"
  Open(50, File="Sincrotron.dat")
  Write(50,*) "#  log(E[eV])              log(E*L(E)[erg/s])"

  !-------------------------------------------------------------------------------
  !Calculamos los pasos de integracion para los distintos procesos
  !-------------------------------------------------------------------------------
  
  ke = (Eemax1/Eemin)**(0.001d0)
  
  kp = (Epmax1/Epmin)**(0.001d0)

  !-------------------------------------------------------------------------------
  !Calculamos los tiempos de enfriamiento y el tiempo de aceleracion (Ejercicios 1 y 2)
  !-------------------------------------------------------------------------------

  Ee = Eemin
  Ep = Epmin
  
  Do i=1,1000
     
     Write(10,*) log10(Ee/eV) , log10(Ep/eV), log10(Tcool(Ee,b_sinc(Ee))),&
          & log10(Tcool(Ee,b_br(Ee))), log10(Tcool(Ee,b_ic(Ee))), &
          & log10(Tcool(Ep,b_pp(Ep))),log10(Tacel(Ee)), log10(Tacel(Ep)),&
          & log10(Tcool(Ee,b_sinc(Ee)+b_br(Ee)+b_ic(Ee)))

     If(abs(log10(Tacel(Ee))-log10(Tsnr))<=1.d-2)then
        Eemax = Ee
     End If
     If(abs(log10(Tacel(Ep))-log10(Tsnr))<=1.d-2)then
        Epmax = Ep
     End If
     
     Ep = Ep * kp

     Ee = Ee * ke

  End Do

  Close(10)

  Ee = Eemin
  Ep = Epmin
  Eph = Ephmin_e
  Ephmax_e = Eemax - Eer
  Ephmax_p = Epmax

  kph_e = (Ephmax_e/Ephmin_e)**(0.001d0)

  kph_p = (Ephmax_p/Ephmin_p)**(0.001d0)

  ke = (Eemax/Eemin)**(0.001d0)
  
  kp = (Epmax/Epmin)**(0.001d0)

  Call Proportionality(Eemax,Epmax,Ae,Ap)
     
     !-------------------------------------------------------------------------------
     !Calculamos la SED por Bremsstrahlung
     !-------------------------------------------------------------------------------
  Do k=1,nmax !barro las energias de los fotones
     
     Lbremss = 0.d0
     Ee = Eemin
     
     Do i = 1, nmax !barro las energias de los electrones
        
        base = Ee*ke-Ee
        
        If(Ee >= Eph)then
           Lbremss = Lbremss + fbremss(Eph,Ee,Ae,Eemax)*base
        End If
        
        Ee = Ee*ke
        
     End Do
     
     Write(20,*) log10(Eph*erg), log10(Eph*Lbremss)

     Eph = Eph*kph_e
     
  End Do
  
     !-------------------------------------------------------------------------------
     !Calculamos la SED para Inverse Compton
     !-------------------------------------------------------------------------------
  Eph = Ephmin_e

  Do k=1,nmax

     Lic = 0.d0
     Ee = Eemin
     
     Do i = 1, nmax-1
        
        base = Ee*ke-Ee
        
        Lic = Lic + fic(Ee,Ephcmb,Eph,Ae,Eemax)*base
        
        Ee = Ee*ke
        
     End Do
     
     Write(30,*) log10(Eph*erg), log10(Eph*Lic)
     
     Eph = Eph*kph_e
     
  End Do
  
  !-------------------------------------------------------------------------------
  !Calculamos la SED por Proton-Proton
  !-------------------------------------------------------------------------------
  Eph = Ephmin_p
  Ep = Epmin
  
  Do k=1,nmax
     
     Lpp = 0.d0
     Epimin = Eph + (Epir**2)/(4.d0*Eph)
     Epimax = 0.1d0*Epmax
     Epi = Epimin
     
     kpi = (Epimax/Epimin)**(1.d0/(real(nmax)-1.d0))
     
     Do i=1,nmax-1
        If(Epi<Ep)then
           base = Epi*kpi-Epi
           Lpp = Lpp + fpp(Epi,Ap,Epmax)*base
        End If
        Epi=Epi*kpi
     End Do
     
     Lpp = Eph * Lpp
     
     Write(40,*) log10(Eph*erg), log10(Eph*Lpp)
     
     Eph = Eph*kph_p
     Ep=Ep*kp
     
  End Do
    
    !-------------------------------------------------------------------------------
  !Calculamos la SED por Sincrotron
  !-------------------------------------------------------------------------------
  Eph = Ephmin_e

  
  Do k=1,nmax
     
     Ee = Eemin
     Lsy = 0.d0
     
     Do i=1,nmax-1
        base = Ee*ke-Ee
        Lsy = Lsy + fsyn(Ee,Eph,Ae,Eemax)*base
        Ee=Ee*ke
     End Do
     
     Write(50,*) log10(Eph*erg),log10(Eph*Lsy)
     
     Eph = Eph*kph_e
     
  End Do
  
  Close(20)
  Close(30)
  Close(40)
  Close(50)
  
End Program TrabajoFinal
