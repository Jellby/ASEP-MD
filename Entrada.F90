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
! En este módulo se definen las variables globales del cálculo
!-------------------------------------------------------------------------------
MODULE Parametros

#ifndef LLL
#define LLL 256
#endif

CHARACTER(LEN=LLL) :: EjecutableMM,EjecutableQM,EntradaMM,EntradaQM,SalidaQM, &
                      FChkGaussian,AltCoordenadas,FicheroCargas,Exten, &
                      Extension,OptContinuacion,CargasExternas,DumpextMoldy, &
                      SalidaMM,SalidaOpt,TrayectoriaMM
INTEGER :: ProgramaMM,ProgramaQM,TipoCargas,TipoCoordenadas,CalcHessiana, &
           Actualizacion,MetodoOptim,MaxIterOpt,HessInicial,BusquedaLineal, &
           Subdivisiones,TipoCavidad,NumConfig,TipoMalla,TipoReduccion, &
           MaxIter,Inicio
DOUBLE PRECISION :: ConvGradOpt,ConvEnerOpt,ConvPasoOpt,MaxPasoOpt, &
                    RadioCavidad,RadioDisolvente,Dielectrica,DistMalla, &
                    FactorMalla,CorteRF,DistCargas
LOGICAL :: EstadoTransicion,InicioVacio

END MODULE Parametros

!-------------------------------------------------------------------------------
! Este módulo tiene subrutinas para leer la entrada general
!-------------------------------------------------------------------------------
MODULE Entrada

#ifndef LLL
#define LLL 256
#endif

#ifndef LANG
#define LANG Espanol
#endif

USE LANG

!-------------------------------------------------------------------------------
! VarComp: Nombres de las variables en minúscula
! VarEnt:  Vector para marcar las variables especificadas en la entrada
!-------------------------------------------------------------------------------
#ifdef __PORTLAND__
CHARACTER(LEN=32), DIMENSION(SIZE(Variables,1)) :: VarComp
#else
CHARACTER(LEN=LEN(Variables(1))), DIMENSION(SIZE(Variables,1)) :: VarComp
#endif
LOGICAL, DIMENSION(SIZE(Variables,1)) :: VarEnt

CONTAINS
!LeerEntrada(Fich)
!LeerVariable(Lin)
!ValoresDefecto

!-------------------------------------------------------------------------------
! Lee el fichero de entrada general
!-------------------------------------------------------------------------------
! Fich:    Fichero de entrada
! Fin:     Variable para controlar el fin del fichero
! Linea:   La línea que se va leyendo
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE LeerEntrada(Fich)
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich

  CHARACTER(LEN=LLL) :: Linea
  INTEGER :: i,Fin

  !Para empezar, ninguna variable está en la entrada
  VarEnt(:)=.FALSE.

  !Se construyen los nombres de las variables en minúscula
  VarComp(:)=Variables(:)
  DO i=1,SIZE(VarComp,1)
    CALL PasarMinusculas(VarComp(i))
  END DO

  DO
    CALL LeerSiguienteLinea(Fich,Linea,Fin)
    IF (Fin /= 0) EXIT
    CALL LeerVariable(Linea)
  END DO

  CALL ValoresDefecto

END SUBROUTINE LeerEntrada

!-------------------------------------------------------------------------------
! Esta subrutina asigna el valor correspondiente a la variable indicada en cada
! línea
!-------------------------------------------------------------------------------
! Lin:     La línea que se está leyendo
! Var:     El nombre de la variable (antes de '=')
! Val:     El valor de la variable (después de '=')
! Pos:     Posición del signo '='
!-------------------------------------------------------------------------------
SUBROUTINE LeerVariable(Lin)
  USE Parametros
  USE Utilidades
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: Lin

  CHARACTER(LEN=LLL) :: Var,Val
  INTEGER :: Pos

  !Se separan las partes antes y después del signo '='
  Pos=INDEX(Lin,'=')
  Var=Lin(:Pos-1)
  Var=TRIM(ADJUSTL(Var))
  CALL PasarMinusculas(Var)
  Val=Lin(Pos+1:)
  Val=TRIM(ADJUSTL(Val))
  IF (INDEX(Val,' ') > 0) Val=Val(:INDEX(Val,' '))

  !TipoCargas
  IF (TRIM(Var) == TRIM(VarComp(1))) THEN
    VarEnt(1)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('mulliken','0')
      TipoCargas=0
     CASE ('esp','qm','1')
      TipoCargas=1
     CASE ('potencial','potential','2')
      TipoCargas=2
     CASE ('externo','external','3')
      TipoCargas=3
     CASE DEFAULT
      CALL Mensaje('LeerVariable',10,.TRUE.)
    END SELECT

  !EntradaMM
  ELSE IF (TRIM(Var) == TRIM(VarComp(2))) THEN
    VarEnt(2)=.TRUE.
    EntradaMM=Val

  !EntradaQM
  ELSE IF (TRIM(Var) == TRIM(VarComp(3))) THEN
    VarEnt(3)=.TRUE.
    EntradaQM=Val

  !ProgramaQM
  ELSE IF (TRIM(Var) == TRIM(VarComp(4))) THEN
    VarEnt(4)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('generico','generic','0')
      ProgramaQM=0
     CASE ('gaussian','1')
      ProgramaQM=1
     CASE ('molcas','2')
      ProgramaQM=2
     CASE DEFAULT
      CALL Mensaje('LeerVariable',11,.TRUE.)
    END SELECT

  !EjecutableQM
  ELSE IF (TRIM(Var) == TRIM(VarComp(5))) THEN
    VarEnt(5)=.TRUE.
    EjecutableQM=Val

  !SalidaQM
  ELSE IF (TRIM(Var) == TRIM(VarComp(6))) THEN
    VarEnt(6)=.TRUE.
    SalidaQM=Val

  !CalcHessiana
  ELSE IF (TRIM(Var) == TRIM(VarComp(7))) THEN
    VarEnt(7)=.TRUE.
    READ(Val,*) CalcHessiana
    IF (CalcHessiana < 0) CalcHessiana=0

  !Actualizacion
  ELSE IF (TRIM(Var) == TRIM(VarComp(8))) THEN
    VarEnt(8)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('no','ninguna','none','0')
      Actualizacion=0
     CASE ('bfgs','1')
      Actualizacion=1
     CASE ('ms','murtagh-sargent','2')
      Actualizacion=2
     CASE ('psb','powell','3')
      Actualizacion=3
     CASE ('bofill','4')
      Actualizacion=4
     CASE DEFAULT
      CALL Mensaje('LeerVariable',15,.TRUE.)
    END SELECT

  !EstadoTransicion
  ELSE IF (TRIM(Var) == TRIM(VarComp(9))) THEN
    VarEnt(9)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('si','yes','1')
      EstadoTransicion=.TRUE.
     CASE DEFAULT
      EstadoTransicion=.FALSE.
    END SELECT

  !MetodoOptim
  ELSE IF (TRIM(Var) == TRIM(VarComp(10))) THEN
    VarEnt(10)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('gradiente','gradient','sd','steepest-descent','0')
      MetodoOptim=0
     CASE ('gradiente-conjugado','cg','conjugated-gradient','1')
      MetodoOptim=1
     CASE ('newton','newton-raphson','quasi-newton','cuasi-newton','nr','qn', &
           '2')
      MetodoOptim=2
     CASE ('rfo','prfo','p-rfo','3')
      MetodoOptim=3
     CASE DEFAULT
      CALL Mensaje('LeerVariable',12,.TRUE.)
    END SELECT

  !TipoCoordenadas
  ELSE IF (TRIM(Var) == TRIM(VarComp(11))) THEN
    VarEnt(11)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('cartesianas','cartesian','0')
      TipoCoordenadas=0
     CASE ('ponderadas','mass-weighted','1')
      TipoCoordenadas=1
     CASE ('internas','internal','2')
      TipoCoordenadas=2
     CASE DEFAULT
      CALL Mensaje('LeerVariable',20,.TRUE.)
    END SELECT

  !AltCoordenadas
  ELSE IF (TRIM(Var) == TRIM(VarComp(12))) THEN
    VarEnt(12)=.TRUE.
    AltCoordenadas=Val

  !MaxIterOpt
  ELSE IF (TRIM(Var) == TRIM(VarComp(13))) THEN
    VarEnt(13)=.TRUE.
    READ(Val,*) MaxIterOpt
    IF (MaxIterOpt < 0) MaxIterOpt=0

  !ConvGradOpt
  ELSE IF (TRIM(Var) == TRIM(VarComp(14))) THEN
    VarEnt(14)=.TRUE.
    READ(Val,*) ConvGradOpt
    IF (ConvGradOpt < 0.0D0) ConvGradOpt=0.0D0

  !ConvEnerOpt
  ELSE IF (TRIM(Var) == TRIM(VarComp(15))) THEN
    VarEnt(15)=.TRUE.
    READ(Val,*) ConvEnerOpt
    IF (ConvEnerOpt < 0.0D0) ConvEnerOpt=0.0D0

  !ConvPasoOpt
  ELSE IF (TRIM(Var) == TRIM(VarComp(16))) THEN
    VarEnt(16)=.TRUE.
    READ(Val,*) ConvPasoOpt
    IF (ConvPasoOpt < 0.0D0) ConvPasoOpt=0.0D0

  !MaxPasoOpt
  ELSE IF (TRIM(Var) == TRIM(VarComp(17))) THEN
    VarEnt(17)=.TRUE.
    READ(Val,*) MaxPasoOpt
    MaxPasoOpt=ABS(MaxPasoOpt)

  !HessInicial
  ELSE IF (TRIM(Var) == TRIM(VarComp(18))) THEN
    VarEnt(18)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('exacta','exact','0')
      HessInicial=0
     CASE ('unidad','unit','1')
      HessInicial=1
     CASE ('simple','2')
      HessInicial=2
     CASE ('modelo','model','3')
      HessInicial=3
     CASE DEFAULT
      CALL Mensaje('LeerVariable',21,.TRUE.)
    END SELECT

  !BusquedaLineal
  ELSE IF (TRIM(Var) == TRIM(VarComp(19))) THEN
    VarEnt(19)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('no','ninguna','none','0')
      BusquedaLineal=0
     CASE ('simple','1')
      BusquedaLineal=1
     CASE ('completa','full','2')
      BusquedaLineal=2
     CASE DEFAULT
      CALL Mensaje('LeerVariable',22,.TRUE.)
    END SELECT

  !FicheroCargas
  ELSE IF (TRIM(Var) == TRIM(VarComp(20))) THEN
    VarEnt(20)=.TRUE.
    FicheroCargas=Val

  !RadioCavidad
  ELSE IF (TRIM(Var) == TRIM(VarComp(21))) THEN
    VarEnt(21)=.TRUE.
    READ(Val,*) RadioCavidad
    IF (RadioCavidad < 0.0D0) RadioCavidad=-RadioCavidad

  !Subdivisiones
  ELSE IF (TRIM(Var) == TRIM(VarComp(22))) THEN
    VarEnt(22)=.TRUE.
    READ(Val,*) Subdivisiones
    IF (Subdivisiones < 0) Subdivisiones=0
    IF (Subdivisiones > 5) Subdivisiones=5

  !RadioDisolvente
  ELSE IF (TRIM(Var) == TRIM(VarComp(23))) THEN
    VarEnt(23)=.TRUE.
    READ(Val,*) RadioDisolvente
    IF (RadioDisolvente < 0.0D0) RadioDisolvente=0.0D0

  !OptContinuacion
  ELSE IF (TRIM(Var) == TRIM(VarComp(24))) THEN
    VarEnt(24)=.TRUE.
    OptContinuacion=Val

  !FChkGaussian
  ELSE IF (TRIM(Var) == TRIM(VarComp(25))) THEN
    VarEnt(25)=.TRUE.
    FChkGaussian=Val

  !CargasExternas
  ELSE IF (TRIM(Var) == TRIM(VarComp(26))) THEN
    VarEnt(26)=.TRUE.
    CargasExternas=Val

  !ProgramaMM
  ELSE IF (TRIM(Var) == TRIM(VarComp(27))) THEN
    VarEnt(27)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('generico','generic','0')
      ProgramaMM=0
     CASE ('moldy','1')
      ProgramaMM=1
     CASE DEFAULT
      CALL Mensaje('LeerVariable',32,.TRUE.)
    END SELECT

  !EjecutableMM
  ELSE IF (TRIM(Var) == TRIM(VarComp(28))) THEN
    VarEnt(28)=.TRUE.
    EjecutableMM=Val

  !DumpextMoldy
  ELSE IF (TRIM(Var) == TRIM(VarComp(29))) THEN
    VarEnt(29)=.TRUE.
    DumpextMoldy=Val

  !NumConfig
  ELSE IF (TRIM(Var) == TRIM(VarComp(30))) THEN
    VarEnt(30)=.TRUE.
    READ(Val,*) NumConfig
    IF (NumConfig < 0) NumConfig=0

  !Dielectrica
  ELSE IF (TRIM(Var) == TRIM(VarComp(31))) THEN
    VarEnt(31)=.TRUE.
    READ(Val,*) Dielectrica
    IF (Dielectrica < 1.0D0) Dielectrica=1.0D0

  !DistMalla
  ELSE IF (TRIM(Var) == TRIM(VarComp(32))) THEN
    VarEnt(32)=.TRUE.
    READ(Val,*) DistMalla
    IF (DistMalla < 0.0D0) DistMalla=ABS(DistMalla)

  !FactorMalla
  ELSE IF (TRIM(Var) == TRIM(VarComp(33))) THEN
    VarEnt(33)=.TRUE.
    READ(Val,*) FactorMalla

  !TipoCavidad
  ELSE IF (TRIM(Var) == TRIM(VarComp(34))) THEN
    VarEnt(34)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('esferica','spherical','0')
      TipoCavidad=0
     CASE ('proporcional','proportional','1')
      TipoCavidad=1
     CASE ('constante','constant','2')
      TipoCavidad=2
     CASE DEFAULT
      CALL Mensaje('LeerVariable',34,.TRUE.)
    END SELECT

  !TipoMalla
  ELSE IF (TRIM(Var) == TRIM(VarComp(35))) THEN
    VarEnt(35)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('cubica','cubic','sc','0')
      TipoMalla=0
     CASE ('compacta','face-centered','fcc','1')
      TipoMalla=1
     CASE ('body-centered','bcc','2')
      TipoMalla=2
     CASE DEFAULT
      CALL Mensaje('LeerVariable',35,.TRUE.)
    END SELECT

  !CorteRF
  ELSE IF (TRIM(Var) == TRIM(VarComp(36))) THEN
    VarEnt(36)=.TRUE.
    READ(Val,*) CorteRF
    IF (CorteRF <= 0.0D0) CorteRF=HUGE(CorteRF)

  !DistCargas
  ELSE IF (TRIM(Var) == TRIM(VarComp(37))) THEN
    VarEnt(37)=.TRUE.
    READ(Val,*) DistCargas
    IF (DistCargas < 0.0D0) DistCargas=ABS(DistCargas)

  !TipoReduccion
  ELSE IF (TRIM(Var) == TRIM(VarComp(38))) THEN
    VarEnt(38)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('normal','-1')
      TipoReduccion=-1
     CASE ('cubica','cubic','sc','0')
      TipoReduccion=0
     CASE ('compacta','face-centered','fcc','1')
      TipoReduccion=1
     CASE ('body-centered','bcc','2')
      TipoReduccion=2
     CASE DEFAULT
      CALL Mensaje('LeerVariable',35,.TRUE.)
    END SELECT

  !SalidaMM
  ELSE IF (TRIM(Var) == TRIM(VarComp(39))) THEN
    VarEnt(39)=.TRUE.
    SalidaMM=Val

  !SalidaOpt
  ELSE IF (TRIM(Var) == TRIM(VarComp(40))) THEN
    VarEnt(40)=.TRUE.
    SalidaOpt=Val

  !MaxIter
  ELSE IF (TRIM(Var) == TRIM(VarComp(41))) THEN
    VarEnt(41)=.TRUE.
    READ(Val,*) MaxIter
    IF (MaxIter < 0) MaxIter=0

  !InicioVacio
  ELSE IF (TRIM(Var) == TRIM(VarComp(42))) THEN
    VarEnt(42)=.TRUE.
    CALL PasarMinusculas(Val)
    SELECT CASE (TRIM(Val))
     CASE ('si','yes','1')
      InicioVacio=.TRUE.
     CASE DEFAULT
      InicioVacio=.FALSE.
    END SELECT

  !Inicio
  ELSE IF (TRIM(Var) == TRIM(VarComp(43))) THEN
    VarEnt(43)=.TRUE.
    READ(Val,*) Inicio
    IF (Inicio < 1) Inicio=1

  !TrayectoriaMM
  ELSE IF (TRIM(Var) == TRIM(VarComp(44))) THEN
    VarEnt(44)=.TRUE.
    TrayectoriaMM=Val

  !Variable inexistente
  ELSE
    CALL Mensaje('LeerVariable',8,.TRUE.)

  END IF

END SUBROUTINE LeerVariable

!-------------------------------------------------------------------------------
! Asigna valores por defecto a las variables
!-------------------------------------------------------------------------------
SUBROUTINE ValoresDefecto
  USE Parametros
  USE Utilidades
  IMPLICIT NONE

  !Se asignan las valores a las variables que no están en la entrada
  IF (.NOT. VarEnt(1)) TipoCargas=1
  IF (.NOT. VarEnt(2)) EntradaMM='input.mm'
  IF (.NOT. VarEnt(3)) EntradaQM='input.qm'
  IF (.NOT. VarEnt(4)) ProgramaQM=1
  IF (.NOT. VarEnt(5)) EjecutableQM='true'
  IF (.NOT. VarEnt(6)) SalidaQM='output.qm'
  IF (.NOT. VarEnt(7)) CalcHessiana=0
  IF (.NOT. VarEnt(8)) Actualizacion=1
  IF (.NOT. VarEnt(9)) EstadoTransicion=.FALSE.
  IF (.NOT. VarEnt(10)) MetodoOptim=3
  IF (.NOT. VarEnt(11)) TipoCoordenadas=2
  IF (.NOT. VarEnt(12)) AltCoordenadas=''
  IF (.NOT. VarEnt(13)) MaxIterOpt=0
  IF (.NOT. VarEnt(14)) ConvGradOpt=3.0D-4
  IF (.NOT. VarEnt(15)) ConvEnerOpt=1.0D-6
  IF (.NOT. VarEnt(16)) ConvPasoOpt=3.0D-4
  IF (.NOT. VarEnt(17)) MaxPasoOpt=3.0D-1
  IF (.NOT. VarEnt(18)) HessInicial=2
  IF (.NOT. VarEnt(19)) BusquedaLineal=0
  IF (.NOT. VarEnt(20)) FicheroCargas=''
  IF (.NOT. VarEnt(21)) RadioCavidad=4.5D0
  IF (.NOT. VarEnt(22)) Subdivisiones=0
  IF (.NOT. VarEnt(23)) RadioDisolvente=0.0D0
  IF (.NOT. VarEnt(24)) OptContinuacion=''
  IF (.NOT. VarEnt(25)) FChkGaussian='fchk'
  IF (.NOT. VarEnt(26)) CargasExternas=''
  IF (.NOT. VarEnt(27)) ProgramaMM=1
  IF (.NOT. VarEnt(28)) EjecutableMM='true' 
  IF (.NOT. VarEnt(29)) DumpextMoldy='dumpext' 
  IF (.NOT. VarEnt(30)) NumConfig=0
  IF (.NOT. VarEnt(31)) Dielectrica=1.0D0
  IF (.NOT. VarEnt(32)) DistMalla=0.5D0
  IF (.NOT. VarEnt(33)) FactorMalla=0.75D0
  IF (.NOT. VarEnt(34)) TipoCavidad=0
  IF (.NOT. VarEnt(35)) TipoMalla=0
  IF (.NOT. VarEnt(36)) CorteRF=HUGE(CorteRF)
  IF (.NOT. VarEnt(37)) DistCargas=1.0D0
  IF (.NOT. VarEnt(38)) TipoReduccion=-1
  IF (.NOT. VarEnt(39)) SalidaMM='output.mm'
  IF (.NOT. VarEnt(40)) SalidaOpt='optim'
  IF (.NOT. VarEnt(41)) MaxIter=0
  IF (.NOT. VarEnt(42)) InicioVacio=.TRUE.
  IF (.NOT. VarEnt(43)) Inicio=1
  IF (.NOT. VarEnt(44)) TrayectoriaMM='traj.dcd'

  !Algunos valores dependen de otros, pero sólo si no están en la entrada
  IF (EstadoTransicion) THEN
    IF (.NOT. VarEnt(8)) Actualizacion=3
    IF (.NOT. VarEnt(18)) HessInicial=0
  END IF

  IF (MetodoOptim < 2) THEN
    IF (.NOT. VarEnt(19)) BusquedaLineal=1
  END IF

  IF (TipoCoordenadas == 1) THEN
    IF (.NOT. VarEnt(17)) MaxPasoOpt=3.0D0
  END IF

  IF (ProgramaQM == 2) THEN
    IF (.NOT. VarEnt(1)) TipoCargas=3
  END IF

  IF ((TipoCavidad == 0) .OR. (TipoCavidad == 2)) THEN
    IF (.NOT. VarEnt(21)) RadioCavidad=9.0D0
  END IF

  !Algunas opciones son incompatibles
  IF (EstadoTransicion) THEN
    IF (MetodoOptim < 2) CALL Mensaje('ValoresDefecto',23,.TRUE.)
    IF (BusquedaLineal > 0) CALL Mensaje('ValoresDefecto',24,.TRUE.)
  END IF

  IF (TipoCargas == 3) THEN
    IF (TRIM(FicheroCargas) == '') CALL Mensaje('ValoresDefecto',25,.TRUE.)
  END IF

  IF (ProgramaQM == 2) THEN
    IF (TipoCargas == 1) CALL Mensaje('ValoresDefecto',29,.TRUE.)
  END IF

  IF ((MaxIter < Inicio) .AND. VarEnt(43)) THEN
    CALL Mensaje('ValoresDefecto',40,.TRUE.)
  END IF

END SUBROUTINE ValoresDefecto

END MODULE Entrada
