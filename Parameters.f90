Module Parameters
  Use iso_fortran_env, only:wp => real64
  Implicit None
  
  !------------------------------------------------------------------------------------------------------------
  !  Factores de conversion (cgs)
  !------------------------------------------------------------------------------------------------------------

  Real(wp), parameter :: eV = 1.602176d-12      ! 1 electron-volt

  Real(wp), parameter :: erg = 6.2415091d11     ! 1 ergio (en eV)

  Real(wp), parameter :: mb = 1.d-27            ! 1 milibarn

  Real(wp), parameter :: uG = 1.0d-6            ! 1 microGauss

  Real(wp), parameter :: pc = 3.086d18          ! 1 parsec

  Real(wp), parameter :: kpc = 3.086d21         ! 1 kiloparsec

  Real(wp), parameter :: yr = 3.154d7           ! 1 año

  Real(wp), parameter :: UA = 1.496d13          ! 1 Unidad astronomica

  Real(wp), parameter :: MeV = 1.6021766d-6     ! 1 MeV
  
  Real(wp), parameter :: Tev = 1.6021766        ! 1 TeV
  
  !------------------------------------------------------------------------------------------------------------
  ! Constantes universales (cgs)
  !------------------------------------------------------------------------------------------------------------

  Real(wp), parameter :: pi = 4.d0*atan(1.d0)   ! pi
  
  Real(wp), parameter :: c = 2.99792458d10      ! Velocidad de la luz en el vacío

  Real(wp), parameter :: me = 9.1093826d-28     ! Masa del electron

  Real(wp), parameter :: mp = 1.6726231d-24     ! Masa del proton
  
  Real(wp), parameter :: Eer = 8.18714e-7       ! Energia del electron en reposo

  Real(wp), parameter :: Epr = 1.503284d-3      ! Energía del protón en reposo

  Real(wp), parameter :: Epir = 2.162568d-4     ! Energía del pión neutro en rep

  Real(wp), parameter :: alfa = 1.d0/137.d0     ! Constante de estructura fina
  
  Real(wp), parameter :: qe = 4.8032047d-10     ! Carga del electron
  
  Real(wp), parameter :: re = qe**2/Eer        ! Radio clasico del elctron

  Real(wp), parameter :: h = 6.626069d-27       ! Constante de Planck

  Real(wp), parameter :: kb = 1.38065d-16       ! Constante de Boltzmann

  Real(wp), parameter :: sigmaSB = 5.6704d-5    ! Constante de Stefan-Boltzmann

  Real(wp), parameter :: sigmath = 6.65d-25     ! Seccion eficaz de Thomson

  Real(wp), parameter :: Ucmb= 4.01d-13         ! Densidad de energia del campo de fotones del CMB

  Real(wp), parameter :: Tcmb=2.725d0           ! Temperatura del CMB
  
  Real(wp), parameter :: Ubb = (4.d0 * sigmaSB * Tcmb**4d0)/c  !Densidad de energia de un cuerpo negro emitiendo a Tcmb

  Real(wp), parameter :: Ephcmb = 2.7*kb*Tcmb   ! Energia promedio de los fotones del CMB

  !------------------------------------------------------------------------------------------------------------
  ! Datos del problema
  !------------------------------------------------------------------------------------------------------------
  
  Real(wp), parameter :: Dsnr = 1.d0*kpc        ! Distancia al remanente

  Real(wp), parameter :: Esnr = 1.d51           ! Energia liberada en la explosion

  Real(wp), parameter :: Tsnr = 1600.d0*yr      ! Edad del remanente

  Real(wp), parameter :: Rout = 10.d0*pc        ! Radio externo del shell

  Real(wp), parameter :: Rin = 7.5d0*pc         ! Radio interno del shell

  Real(wp), parameter :: Vol = (4.d0/3.d0)*pi*(Rout**3.d0-Rin**3.d0) ! Volumen de la cascara
  
  Real(wp), parameter :: fvol = 0.16d0          ! Fraccion emisora del volumen total

  Real(wp), parameter :: na = 0.2d0           ! Densidad del medio interestelar

  Real(wp), parameter :: B = 7.5*uG              ! Campo magnetico en la region de aceleracion

  Real(wp), parameter :: Xrel = 0.23d0           ! Fraccion de energia en particulas relativistas

  Real(wp), parameter :: eta = 3.d-5            ! Eficiencia del mecanismo de aceleracion

  Real(wp), parameter :: Kep = 1.4d-3            ! Cociente de energia inyectada en electrones y protones (E/Ep)

  Real(wp), parameter :: Z = 1.5d0              ! Numero atomico de los atomos del medio interestelar

  !------------------------------------------------------------------------------------------------------------
  !  Parametros generales del problema
  !------------------------------------------------------------------------------------------------------------

  Real(wp), parameter :: Eemax1 = 1.D16*eV       ! Energia maxima de los electrones para el ejercicio 1

  Real(wp), parameter :: Epmax1 = 1.D20*eV       ! Energia maxima de los protones para el ejercicio 1

  
  Real(wp), parameter :: Eemin = 1.d0*MeV       ! Energia minima de los electrones

  Real(wp), parameter :: Epmin = 1.d9*eV        ! Energia minima de los protones
  
!  Real(wp), parameter :: Eemax = 3.1754d13*eV   ! Energia maxima de los electrones

!  Real(wp), parameter :: Epmax = 3.1747d13*eV   ! Energia maxima de los protones

  Real(wp), parameter :: Ephmin_e = 1.d-7*eV    !Energia minima de los fotones

!  Real(wp), parameter :: Ephmax_e = Eemax - Eer !Energia maxima de los fotones

  Real(wp), parameter :: Ephmin_p = 1.d7*ev

!  Real(wp), parameter :: Ephmax_p = Epmax

  Real(wp), parameter :: Alpha = 0.579469       !Indice espectral de la distribucion de fotones emitida por Sincrotron

!  Real(wp), parameter :: p = 2.d0*alpha+1.d0    !Indice espectral de la distribucion de electrones
  Real(wp), parameter :: p=2.d0
  
!  Real(wp), parameter :: Ae = 6.163218d45       !Constante de normalizacion de la distribucion de electrones

!  Real(wp), parameter :: Ap = 7.314683d48       !Constante de normalizacion de la distribucion de protones
  
  !------------------------------------------------------------------------------------------------------------
  !  Sincrotron
  !------------------------------------------------------------------------------------------------------------
  
  Real(wp), parameter :: b0sinc = 6.6d-4 *( (B**2) / (Eer**2) ) * eV ! bsinc(E) = b0sinc * E**2

  !------------------------------------------------------------------------------------------------------------
  !  Inverse Compton
  !------------------------------------------------------------------------------------------------------------

  Real(wp), parameter :: Ttilde = kb * Tcmb / Eer

  !------------------------------------------------------------------------------------------------------------
  !  Bremsstrahlung Relativista
  !------------------------------------------------------------------------------------------------------------

  Real(wp), parameter :: b0br = -4.d0 * na * Z**2 * re**2 * alpha * c *(log(183.d0*Z**(-1.d0/3.d0))-1.d0/18.d0)

                         !bbr(E) = b0br * E
  
  !------------------------------------------------------------------------------------------------------------
  !  Proton-proton
  !------------------------------------------------------------------------------------------------------------

  Real(wp), parameter :: Eth = 1.22d9 * eV      !Energia umbral para la interaccion p-p

  Real(wp), parameter :: Kpp = 0.5d0            !

  End Module Parameters
