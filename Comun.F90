!##############################################################################
!# Copyright 2011,2012 Ignacio Fdez. Galván, M. Luz Sánchez,                  #
!#                     Aurora Muñoz Losa, M. Elena Martín, Manuel A. Aguilar  #
!#                                                                            #
!# This file is part of ASEP-MD.                                              #
!#                                                                            #
!# ASEP-MD is free software: you can redistribute it and/or modify it under   #
!# the terms of the GNU General Public License as published by the Free       #
!# Software Foundation, either version 3 of the License, or (at your option)  #
!# any later version.                                                         #
!#                                                                            #
!# ASEP-MD is distributed in the hope that it will be useful, but WITHOUT ANY #
!# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  #
!# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      #
!# details.                                                                   #
!#                                                                            #
!# You should have received a copy of the GNU General Public License along    #
!# with ASEP-MD. If not, see <http://www.gnu.org/licenses/>.                  #
!##############################################################################

!-------------------------------------------------------------------------------
! Se define el tipo de datos "Atomo"
!-------------------------------------------------------------------------------
MODULE TipoAtomo
  TYPE Atomo
    DOUBLE PRECISION, DIMENSION(3) :: pos
    DOUBLE PRECISION :: q,m
    INTEGER :: id,z
    CHARACTER(LEN=16) :: nom
  END TYPE Atomo
END MODULE TipoAtomo

!-------------------------------------------------------------------------------
! En este módulo se define el sistema que se calculará
!-------------------------------------------------------------------------------
MODULE Sistema

  USE TipoAtomo
  TYPE(Atomo), DIMENSION(:), ALLOCATABLE :: Soluto,Disolvente,Disolvente2
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MallaSoluto,CargasDisolvente
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: InterAtom
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: QInter
  INTEGER :: MoleculasSoluto,MoleculasDisolvente,MoleculasDisolvente2, &
             TipoPotencial
  CHARACTER(LEN=32) :: NombreSoluto,NombreDisolvente,NombreDisolvente2

END MODULE Sistema

!-------------------------------------------------------------------------------
! En este modulo se encuentran algunos resultados del cálculo cuántico
!-------------------------------------------------------------------------------
MODULE DatosQM

  USE TipoAtomo
  TYPE(Atomo), DIMENSION(:), ALLOCATABLE :: MolQM
  INTEGER :: CargaQM,MultipQM
  DOUBLE PRECISION :: EnergiaQM,Energia2QM
  DOUBLE PRECISION, DIMENSION(3) :: DipoloQM
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GradQM,Grad2QM
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: HessQM,DisolvQM,PotQM

END MODULE DatosQM

!-------------------------------------------------------------------------------
! En este módulo se definen las interacciones QM/MM
!-------------------------------------------------------------------------------
MODULE DatosQMMM

  DOUBLE PRECISION :: EnergiaEQM,EnergiaEMM,EnergiaVdW
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GradVdW
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: HessVdW

END MODULE DatosQMMM

!-------------------------------------------------------------------------------
! En este módulo se definen algunas constantes y factores de conversión
!-------------------------------------------------------------------------------
MODULE Unidades

  !Factores de conversion
  DOUBLE PRECISION, PARAMETER :: &
    SIMasa         =1.09776929D30, &      !1 kg [me]
    SILongitud     =1.88972613289D10, &   !1 m [a0]
    SICarga        =6.241509647D18, &     !1 C [e]
    SITiempo       =4.134137333656D16, &  !1 s [hbar/Eh]
    SIEnergia      =2.29371269D17, &      !1 J [Eh]
    DebyeAtomica   =3.934303073D-1, &     !1 D [e*a0]
    AngstromAtomica=1.88972613289D0, &    !1 A [a0]
    KcalmolAtomica =1.5936014507D-3, &    !1 kcal/mol [Eh]
    AmuAtomica     =1.82288848D3, &       !1 u [me]
    KBoltzmann     =3.1668153D-6, &       !kB [Eh/K]
    Pi             =3.14159265358979D0, & !Pi
    Grado          =1.74532925199433D-2   !Pi/180
  
  !Símbolos atómicos
  CHARACTER(LEN=2), DIMENSION(0:92), PARAMETER :: Simbolo=(/'Bq', &
    'H ',                              'He', &
    'Li','Be','B ','C ','N ','O ','F ','Ne', &
    'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
    'K ','Ca', &
    'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
              'Ga','Ge','As','Se','Br','Kr', &
    'Rb','Sr', &
    'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', &
              'In','Sn','Sb','Te','I ','Xe', &
    'Cs','Ba', &
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
    'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
              'Tl','Pb','Bi','Po','At','Rn', &
    'Fr','Ra', &
    'Ac','Th','Pa','U '/)

  !Ver también:
  ! "Bonding and Structure", N. W. Alcock, 1990

  !Radios covalentes
  ! http://www.webofelements.com
  ! Sin datos -> NDRC
  REAL, PARAMETER :: NDRC=1.50
  DOUBLE PRECISION, DIMENSION(0:92), PARAMETER :: RadiosCov=(/0.00, &
    0.37,                              0.32, &
    1.34,0.90,0.82,0.77,0.75,0.73,0.71,0.69, &
    1.54,1.30,1.18,1.11,1.06,1.02,0.99,0.97, &
    1.96,1.74, &
    1.44,1.36,1.25,1.27,1.39,1.25,1.26,1.21,1.38,1.31, &
              1.26,1.22,1.19,1.16,1.14,1.10, &
    2.11,1.92, &
    1.62,1.48,1.37,1.45,1.56,1.26,1.35,1.31,1.53,1.48, &
              1.44,1.41,1.38,1.35,1.33,1.30, &
    2.25,1.98, &
    1.69,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC,NDRC, &
    1.60,1.50,1.38,1.46,1.59,1.28,1.37,1.28,1.44,1.49, &
              1.48,1.47,1.46,NDRC,NDRC,1.45, &
    NDRC,NDRC, &
    NDRC,NDRC,NDRC,NDRC/) * AngstromAtomica

  !Radios de van der Waals
  ! Bondi, J. Phys. Chem. 68 (1964), 441
  ! Sin datos -> NDVW
  REAL, PARAMETER :: NDVW=2.00
  DOUBLE PRECISION, DIMENSION(0:92), PARAMETER :: RadiosVdW=(/0.00, &
    1.20,                              1.40, &
    1.82,NDVW,NDVW,1.70,1.55,1.52,1.47,1.54, &
    2.27,1.73,NDVW,2.10,1.80,1.80,1.75,1.88, &
    2.75,NDVW, &
    NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,1.63,1.40,1.39, &
              1.87,NDVW,1.85,1.90,1.85,2.02, &
    NDVW,NDVW, &
    NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,1.63,1.72,1.58, &
              1.93,2.17,NDVW,2.06,1.98,2.16, &
    NDVW,NDVW, &
    NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW, &
    NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,NDVW,1.75,1.66,1.55, &
              1.96,2.02,NDVW,NDVW,NDVW,NDVW, &
    NDVW,NDVW, &
    NDVW,NDVW,NDVW,NDVW/) * AngstromAtomica

END MODULE Unidades

!-------------------------------------------------------------------------------
! Mensajes y variables en español
!-------------------------------------------------------------------------------
MODULE Espanol

#ifndef LLL
#define LLL 256
#endif

  CHARACTER(LEN=LLL), DIMENSION(6) :: Titulo
    DATA Titulo(1) /'========================================================='/
    DATA Titulo(2) /'|                    ASEP-MD v. 2.0                     |'/
    DATA Titulo(3) /'|                  ------------------                   |'/
    DATA Titulo(4) /'| Copyright 2011,2012 I.F.G., M.L.S., A.M.L., M.E.M.,   |'/
    DATA Titulo(5) /'|                     M.A.A.                            |'/
    DATA Titulo(6) /'========================================================='/

  CHARACTER(LEN=32), DIMENSION(44) :: Variables
    DATA Variables( 1) /'TipoCargas'/       !TipoCargas
    DATA Variables( 2) /'EntradaMM'/        !EntradaMM
    DATA Variables( 3) /'EntradaQM'/        !EntradaQM
    DATA Variables( 4) /'ProgramaQM'/       !ProgramaQM
    DATA Variables( 5) /'EjecutableQM'/     !EjecutableQM
    DATA Variables( 6) /'SalidaQM'/         !SalidaQM
    DATA Variables( 7) /'CalcHessiana'/     !CalcHessiana
    DATA Variables( 8) /'Actualizacion'/    !Actualizacion
    DATA Variables( 9) /'EstadoTransicion'/ !EstadoTransicion
    DATA Variables(10) /'MetodoOptim'/      !MetodoOptim
    DATA Variables(11) /'TipoCoordenadas'/  !TipoCoordenadas
    DATA Variables(12) /'AltCoordenadas'/   !AltCoordenadas
    DATA Variables(13) /'MaxIterOpt'/       !MaxIterOpt
    DATA Variables(14) /'ConvGradOpt'/      !ConvGradOpt
    DATA Variables(15) /'ConvEnerOpt'/      !ConvEnerOpt
    DATA Variables(16) /'ConvPasoOpt'/      !ConvPasoOpt
    DATA Variables(17) /'MaxPasoOpt'/       !MaxPasoOpt
    DATA Variables(18) /'HessInicial'/      !HessInicial
    DATA Variables(19) /'BusquedaLineal'/   !BusquedaLineal
    DATA Variables(20) /'FicheroCargas'/    !FicheroCargas
    DATA Variables(21) /'RadioCavidad'/     !RadioCavidad
    DATA Variables(22) /'Subdivisiones'/    !Subdivisiones
    DATA Variables(23) /'RadioDisolvente'/  !RadioDisolvente
    DATA Variables(24) /'OptContinuacion'/  !OptContinuacion
    DATA Variables(25) /'FChkGaussian'/     !FChkGaussian
    DATA Variables(26) /'CargasExternas'/   !CargasExternas
    DATA Variables(27) /'ProgramaMM'/       !ProgramaMM
    DATA Variables(28) /'EjecutableMM'/     !EjecutableMM
    DATA Variables(29) /'DumpextMoldy'/     !DumpextMoldy
    DATA Variables(30) /'NumConfig'/        !NumConfig
    DATA Variables(31) /'Dielectrica'/      !Dielectrica
    DATA Variables(32) /'DistMalla'/        !DistMalla
    DATA Variables(33) /'FactorMalla'/      !FactorMalla
    DATA Variables(34) /'TipoCavidad'/      !TipoCavidad
    DATA Variables(35) /'TipoMalla'/        !TipoMalla
    DATA Variables(36) /'CorteRF'/          !CorteRF
    DATA Variables(37) /'DistCargas'/       !DistCargas
    DATA Variables(38) /'TipoReduccion'/    !TipoReduccion
    DATA Variables(39) /'SalidaMM'/         !SalidaMM
    DATA Variables(40) /'SalidaOpt'/        !SalidaOpt
    DATA Variables(41) /'MaxIter'/          !MaxIter
    DATA Variables(42) /'InicioVacio'/      !InicioVacio
    DATA Variables(43) /'Inicio'/           !Inicio
    DATA Variables(44) /'TrayectoriaMM'/    !TrayectoriaMM
  CHARACTER(LEN=LLL), DIMENSION(0:40) :: Errores
    DATA Errores( 0) /'Error no definido.'/
    DATA Errores( 1) /'La matriz no es simetrica.'/
    DATA Errores( 2) /'Superado el numero maximo de iteraciones.'/
    DATA Errores( 3) /'La matriz es singular.'/
    DATA Errores( 4) /'Error en la plantilla de Gaussian.'/
    DATA Errores( 5) /'Error en la plantilla de Moldy.'/
    DATA Errores( 6) /'Dos moleculas con el mismo nombre.'/
    DATA Errores( 7) /'Tipo de potencial interatomico no valido.'/
    DATA Errores( 8) /'Variable inexistente.'/
    DATA Errores( 9) /'Error en la salida de Gaussian.'/
    DATA Errores(10) /'Tipo de cargas atomicas no reconocido.'/
    DATA Errores(11) /'Programa de calculo cuantico no reconocido.'/
    DATA Errores(12) /'Metodo de optimizacion no reconocido.'/
    DATA Errores(13) /'Fallo en el calculo del factor lambda.'/
    DATA Errores(14) /'Paso Newton-Raphson.'/
    DATA Errores(15) /'Formula de actualizacion de la hessiana no reconocida.'/
    DATA Errores(16) /'No hay potencial al que ajustar las cargas.'/
    DATA Errores(17) /'Formato incorrecto.'/
    DATA Errores(18) /'Atomo inexistente.'/
    DATA Errores(19) /'La conversion de coordenadas diverge.'/
    DATA Errores(20) /'Tipo de coordenadas no reconocido.'/
    DATA Errores(21) /'Tipo de hessiana inicial no reconocido.'/
    DATA Errores(22) /'Tipo de busqueda lineal no reconocido.'/
    DATA Errores(23) /'No se puede buscar un ET con metodos de 1er orden.'/
    DATA Errores(24) /'No se puede usar busqueda lineal para un ET.'/
    DATA Errores(25) /'FicheroCargas es necesario para cargas de tipo externo.'/
    DATA Errores(26) /'Error en el fichero de continuacion de la optimizacion.'/
    DATA Errores(27) /'Error en la salida de Molcas.'/
    DATA Errores(28) /'Molcas no termino correctamente.'/
    DATA Errores(29) /'Molcas no calcula cargas ESP.'/
    DATA Errores(30) /'No se puede repetir indice en moleculas diferentes.'/
    DATA Errores(31) /'No se pudo encontrar una unidad libre.'/
    DATA Errores(32) /'Programa de dinamica molecular no reconocido.'/
    DATA Errores(33) /'Numero de configuraciones mayor que el disponible.'/
    DATA Errores(34) /'Tipo de cavidad no reconocido.'/
    DATA Errores(35) /'Tipo de malla no reconocido.'/
    DATA Errores(36) /'El numero de cargas es mayor que el numero de puntos con potencial conocido.'/
    DATA Errores(37) /'Error en la definicion del sistema.'/
    DATA Errores(38) /'El MM generico requiere una sola molecula de soluto.'/
    DATA Errores(39) /'Numero incorrecto de atomos en las configuraciones.'/
    DATA Errores(40) /'La primera iteracion no puede ser mayor que el maximo.'/
  CHARACTER(LEN=LLL), DIMENSION(77) :: Textos
    DATA Textos( 1) /'Maxima componente del gradiente:'/
    DATA Textos( 2) /'Diferencia de energia:'/
    DATA Textos( 3) /'Maxima componente del incremento:'/
    DATA Textos( 4) /'Se ha alcanzado la convergencia'/
    DATA Textos( 5) /'Se ha superado el numero maximo de iteraciones'/
    DATA Textos( 6) /'Numero de coordenadas:'/
    DATA Textos( 7) /'COORDENADAS Y GRADIENTE INICIALES'/
    DATA Textos( 8) /'ITERACION'/
    DATA Textos( 9) /'Cartesianas (angstrom, Eh)'/
    DATA Textos(10) /'Cartesianas ponderadas (angstrom, 1 uma, Eh)'/
    DATA Textos(11) /'Internas (angstrom, grados, Eh)    coordenada        gradiente'/
    DATA Textos(12) /'Diedro'/
    DATA Textos(13) /'Angulo'/
    DATA Textos(14) /'Distancia'/
    DATA Textos(15) /'Energia total:'/
    DATA Textos(16) /'Modulo del gradiente:'/
    DATA Textos(17) /'**Calculo de la hessiana exacta**'/
    DATA Textos(18) /'**Busqueda lineal del minimo**'/
    DATA Textos(19) /'Hessiana exacta (Eh/(a0^2))'/
    DATA Textos(20) /'Hessiana aproximada y proyectada (Eh/(a0^2))'/
    DATA Textos(21) /'Combinacion'/
    DATA Textos(22) /'Continuacion de un calculo interrumpido'/
    DATA Textos(23) /'Energia QM:'/
    DATA Textos(24) /'Despues de la simulacion MM'/
    DATA Textos(25) /'Despues del calculo QM'/
    DATA Textos(26) /'Interaccion <q-q> soluto-disolvente:'/
    DATA Textos(27) /'Interaccion <VdW> soluto-disolvente:'/
    DATA Textos(28) /'Momento dipolar:'/
    DATA Textos(29) /'Energia total:'/
    DATA Textos(30) /'Diferencia de energia:'/
    DATA Textos(31) /'Interaccion QM-<q> soluto-disolvente:'/
    DATA Textos(32) /'Energia interna:'/
    DATA Textos(33) /'Ajuste del disolvente'/
    DATA Textos(34) /'Cargas iniciales:'/
    DATA Textos(35) /'Despues de reducir:'/
    DATA Textos(36) /'RErr:'/
    DATA Textos(37) /'Rms:'/
    DATA Textos(38) /'RRms:'/
    DATA Textos(39) /'Error maximo en el potencial:'/
    DATA Textos(40) /'Puntos para ajustar el potencial:'/
    DATA Textos(41) /'SI'/
    DATA Textos(42) /'NO'/
    DATA Textos(43) /'Mulliken'/
    DATA Textos(44) /'ESP'/
    DATA Textos(45) /'Potencial'/
    DATA Textos(46) /'Externo'/
    DATA Textos(47) /'Generico'/
    DATA Textos(48) /'Gaussian'/
    DATA Textos(49) /'Molcas'/
    DATA Textos(50) /'Ninguna'/
    DATA Textos(51) /'BFGS'/
    DATA Textos(52) /'Murtagh-Sargent'/
    DATA Textos(53) /'Powell'/
    DATA Textos(54) /'Bofill'/
    DATA Textos(55) /'Gradiente'/
    DATA Textos(56) /'Gradiente-conjugado'/
    DATA Textos(57) /'Newton-Raphson'/
    DATA Textos(58) /'RFO'/
    DATA Textos(59) /'Cartesianas'/
    DATA Textos(60) /'Ponderadas'/
    DATA Textos(61) /'Internas'/
    DATA Textos(62) /'Exacta'/
    DATA Textos(63) /'Unidad'/
    DATA Textos(64) /'Simple'/
    DATA Textos(65) /'Modelo'/
    DATA Textos(66) /'Completa'/
    DATA Textos(67) /'Moldy'/
    DATA Textos(68) /'Esferica'/
    DATA Textos(69) /'Proporcional'/
    DATA Textos(70) /'Constante'/
    DATA Textos(71) /'Cubica'/
    DATA Textos(72) /'FCC'/
    DATA Textos(73) /'BCC'/
    DATA Textos(74) /'Normal'/
    DATA Textos(75) /'Traslacion'/
    DATA Textos(76) /'Rotacion'/
    DATA Textos(77) /'Cartesiana'/

END MODULE Espanol

!-------------------------------------------------------------------------------
! Mensajes y variables en ingles
!-------------------------------------------------------------------------------
MODULE English

#ifndef LLL
#define LLL 256
#endif

  CHARACTER(LEN=LLL), DIMENSION(6) :: Titulo
    DATA Titulo(1) /'========================================================='/
    DATA Titulo(2) /'|                    ASEP-MD v. 2.0                     |'/
    DATA Titulo(3) /'|                  ------------------                   |'/
    DATA Titulo(4) /'| Copyright 2011,2012 I.F.G., M.L.S., A.M.L., M.E.M.,   |'/
    DATA Titulo(5) /'|                     M.A.A.                            |'/
    DATA Titulo(6) /'========================================================='/

  CHARACTER(LEN=32), DIMENSION(44) :: Variables
    DATA Variables( 1) /'ChargesType'/     !TipoCargas
    DATA Variables( 2) /'MMInput'/         !EntradaMM
    DATA Variables( 3) /'QMInput'/         !EntradaQM
    DATA Variables( 4) /'QMProgram'/       !ProgramaQM
    DATA Variables( 5) /'QMExec'/          !EjecutableQM
    DATA Variables( 6) /'QMOutput'/        !SalidaQM
    DATA Variables( 7) /'HessianCalc'/     !CalcHessiana
    DATA Variables( 8) /'Update'/          !Actualizacion
    DATA Variables( 9) /'TransitionState'/ !EstadoTransicion
    DATA Variables(10) /'OptimMethod'/     !MetodoOptim
    DATA Variables(11) /'CoordinatesType'/ !TipoCoordenadas
    DATA Variables(12) /'AltCoordinates'/  !AltCoordenadas
    DATA Variables(13) /'MaxIterOpt'/      !MaxIterOpt
    DATA Variables(14) /'GradConvOpt'/     !ConvGradOpt
    DATA Variables(15) /'EnerConvOpt'/     !ConvEnerOpt
    DATA Variables(16) /'StepConvOpt'/     !ConvPasoOpt
    DATA Variables(17) /'MaxStepOpt'/      !MaxPasoOpt
    DATA Variables(18) /'InitialHess'/     !HessInicial
    DATA Variables(19) /'LinearSearch'/    !BusquedaLineal
    DATA Variables(20) /'ChargesFile'/     !FicheroCargas
    DATA Variables(21) /'CavityRadius'/    !RadioCavidad
    DATA Variables(22) /'Subdivisions'/    !Subdivisiones
    DATA Variables(23) /'SolventRadius'/   !RadioDisolvente
    DATA Variables(24) /'OptResume'/       !OptContinuacion
    DATA Variables(25) /'GaussianFChk'/    !FChkGaussian
    DATA Variables(26) /'ExternalCharges'/ !CargasExternas
    DATA Variables(27) /'MMProgram'/       !ProgramaMM
    DATA Variables(28) /'MMExec'/          !EjecutableMM
    DATA Variables(29) /'DumpextMoldy'/    !DumpextMoldy
    DATA Variables(30) /'ConfigNum'/       !NumConfig
    DATA Variables(31) /'Dielectric'/      !Dielectrica
    DATA Variables(32) /'GridDist'/        !DistMalla
    DATA Variables(33) /'GridFactor'/      !FactorMalla
    DATA Variables(34) /'CavityType'/      !TipoCavidad
    DATA Variables(35) /'GridType'/        !TipoMalla
    DATA Variables(36) /'RFCutoff'/        !CorteRF
    DATA Variables(37) /'ChargesDist'/     !DistCargas
    DATA Variables(38) /'ReductionType'/   !TipoReduccion
    DATA Variables(39) /'MMOutput'/        !SalidaMM
    DATA Variables(40) /'OptOutput'/       !SalidaOpt
    DATA Variables(41) /'MaxIter'/         !MaxIter
    DATA Variables(42) /'VacuumStart'/     !InicioVacio
    DATA Variables(43) /'Start'/           !Inicio
    DATA Variables(44) /'MMTrajectory'/    !TrayectoriaMM
  CHARACTER(LEN=LLL), DIMENSION(0:40) :: Errores
    DATA Errores( 0) /'Undefined error.'/
    DATA Errores( 1) /'The matrix is not symmetric.'/
    DATA Errores( 2) /'Maximum number of iterations exceeded.'/
    DATA Errores( 3) /'The matrix is singular.'/
    DATA Errores( 4) /'Error in the Gaussian template.'/
    DATA Errores( 5) /'Error in the Moldy template.'/
    DATA Errores( 6) /'Two molecules have the same name.'/
    DATA Errores( 7) /'Invalid interatomic potential.'/
    DATA Errores( 8) /'Non-existent variable.'/
    DATA Errores( 9) /'Error in the Gaussian output.'/
    DATA Errores(10) /'Unrecognized atomic charges type.'/
    DATA Errores(11) /'Unrecognized quantum calculation program.'/
    DATA Errores(12) /'Unrecognized optimization method.'/
    DATA Errores(13) /'Calculation of lambda factor failed.'/
    DATA Errores(14) /'Newton-Raphson step.'/
    DATA Errores(15) /'Unrecognized Hessian update formula.'/
    DATA Errores(16) /'No potential to fit charges to.'/
    DATA Errores(17) /'Wrong format.'/
    DATA Errores(18) /'Atom does not exist.'/
    DATA Errores(19) /'Divergent coordinate conversion.'/
    DATA Errores(20) /'Unrecognized coordinates type.'/
    DATA Errores(21) /'Unrecognized initial Hessian type.'/
    DATA Errores(22) /'Unrecognized linear search type.'/
    DATA Errores(23) /'Cannot search TS with 1st order method.'/
    DATA Errores(24) /'Cannot use linear search for a TS.'/
    DATA Errores(25) /'ChargesFile needed for charges of external type.'/
    DATA Errores(26) /'Error in resume file for optimization.'/
    DATA Errores(27) /'Error in the Molcas output.'/
    DATA Errores(28) /'Molcas did not finish correctly.'/
    DATA Errores(29) /'Molcas does not calculate ESP charges.'/
    DATA Errores(30) /'Must not repeat id in different molecules.'/
    DATA Errores(31) /'Could not find a free unit.'/
    DATA Errores(32) /'Unrecognized molecular dynamics program.'/
    DATA Errores(33) /'Number of configurations higher than available.'/
    DATA Errores(34) /'Unrecognized cavity type.'/
    DATA Errores(35) /'Unrecognized grid type.'/
    DATA Errores(36) /'Number of charges larger than number of points with known potential.'/
    DATA Errores(37) /'Error in the system definition.'/
    DATA Errores(38) /'Generic MM requires a single solute molecule.'/
    DATA Errores(39) /'Wrong number of atoms in configurations.'/
    DATA Errores(40) /'The first iteration cannot be larger than the maximum.'/
  CHARACTER(LEN=LLL), DIMENSION(77) :: Textos
    DATA Textos( 1) /'Maximum component of the gradient:'/
    DATA Textos( 2) /'Energy difference:'/
    DATA Textos( 3) /'Maximum component of the increment:'/
    DATA Textos( 4) /'Convergence reached'/
    DATA Textos( 5) /'Maximum number of iterations reached'/
    DATA Textos( 6) /'Number of coordinates:'/
    DATA Textos( 7) /'INITIAL COORDINATES AND GRADIENT'/
    DATA Textos( 8) /'ITERATION'/
    DATA Textos( 9) /'Cartesians (angstrom, Eh)'/
    DATA Textos(10) /'Mass-weighted cartesians (angstrom, 1 uma, Eh)'/
    DATA Textos(11) /'Internals (angstrom, degrees, Eh)  coordinate         gradient'/
    DATA Textos(12) /'Dihedral'/
    DATA Textos(13) /'Angle'/
    DATA Textos(14) /'Distance'/
    DATA Textos(15) /'Total energy:'/
    DATA Textos(16) /'Gradient norm:'/
    DATA Textos(17) /'**Exact Hessian calculation**'/
    DATA Textos(18) /'**Linear minimum search**'/
    DATA Textos(19) /'Exact Hessian (Eh/(a0^2))'/
    DATA Textos(20) /'Approximate and projected Hessian (Eh/(a0^2))'/
    DATA Textos(21) /'Combination'/
    DATA Textos(22) /'Resuming an interrupted calculation'/
    DATA Textos(23) /'QM Energy:'/
    DATA Textos(24) /'After the MM simulation'/
    DATA Textos(25) /'After the QM calculation'/
    DATA Textos(26) /'Solute-solvent <q-q> interaction:'/
    DATA Textos(27) /'Solute-solvent <VdW> interaction:'/
    DATA Textos(28) /'Dipole moment:'/
    DATA Textos(29) /'Total energy:'/
    DATA Textos(30) /'Energy difference:'/
    DATA Textos(31) /'Solute-solvent QM-<q> interaction:'/
    DATA Textos(32) /'Internal energy:'/
    DATA Textos(33) /'Solvent fit'/
    DATA Textos(34) /'Initial charges:'/
    DATA Textos(35) /'After reduction:'/
    DATA Textos(36) /'RErr:'/
    DATA Textos(37) /'Rms:'/
    DATA Textos(38) /'RRms:'/
    DATA Textos(39) /'Maximum error in the potential:'/
    DATA Textos(40) /'Points for potential fitting:'/
    DATA Textos(41) /'YES'/
    DATA Textos(42) /'NO'/
    DATA Textos(43) /'Mulliken'/
    DATA Textos(44) /'ESP'/
    DATA Textos(45) /'Potential'/
    DATA Textos(46) /'External'/
    DATA Textos(47) /'Generic'/
    DATA Textos(48) /'Gaussian'/
    DATA Textos(49) /'Molcas'/
    DATA Textos(50) /'None'/
    DATA Textos(51) /'BFGS'/
    DATA Textos(52) /'Murtagh-Sargent'/
    DATA Textos(53) /'Powell'/
    DATA Textos(54) /'Bofill'/
    DATA Textos(55) /'Steepest-descent'/
    DATA Textos(56) /'Conjugate-gradient'/
    DATA Textos(57) /'Newton-Raphson'/
    DATA Textos(58) /'RFO'/
    DATA Textos(59) /'Cartesian'/
    DATA Textos(60) /'Mass-weighted'/
    DATA Textos(61) /'Internal'/
    DATA Textos(62) /'Exact'/
    DATA Textos(63) /'Unit'/
    DATA Textos(64) /'Simple'/
    DATA Textos(65) /'Model'/
    DATA Textos(66) /'Full'/
    DATA Textos(67) /'Moldy'/
    DATA Textos(68) /'Spherical'/
    DATA Textos(69) /'Proportional'/
    DATA Textos(70) /'Constant'/
    DATA Textos(71) /'Cubic'/
    DATA Textos(72) /'FCC'/
    DATA Textos(73) /'BCC'/
    DATA Textos(74) /'Normal'/
    DATA Textos(75) /'Translation'/
    DATA Textos(76) /'Rotation'/
    DATA Textos(77) /'Cartesian'/
  
END MODULE English
