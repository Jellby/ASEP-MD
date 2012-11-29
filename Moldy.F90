!##############################################################################
!# Copyright 2011 Ignacio Fdez. Galván, M. Luz Sánchez, Aurora Muñoz Losa,    #
!#                M. Elena Martín, Manuel A. Aguilar                          #
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

MODULE Moldy

#ifndef LLL
#define LLL 256
#endif

USE Unidades

!-------------------------------------------------------------------------------
! MoldyInput:    Nombre del fichero de entrada de Moldy
! MoldyDump:     Nombre de los ficheros "dump" de Moldy
! MoldySave:     Nombre de los ficheros "save" de Moldy
! MoldyBackup:   Nombre de los ficheros "MDBCK" de Moldy
! MoldyTemp:     Temperatura de la simulación
! MoldyMasa:     Unidad de masa en el fichero de entrada de Moldy, en m_e
! MoldyLongitud: Unidad de longitud en el fichero de entrada de Moldy, en a_0
! MoldyTiempo:   Unidad de tiempo en el fichero de entrada de Moldy, en hbar/Eh
! MoldyCarga:    Unidad de carga en el fichero de entrada de Moldy, en e
! MoldyEnergia:  Unidad de energía en el fichero de entrada de Moldy, en Eh
!-------------------------------------------------------------------------------
CHARACTER(LEN=LLL) :: MoldyInput,MoldyDump,MoldySave,MoldyBackup
INTEGER :: MoldyConfigs
DOUBLE PRECISION :: MoldyTemp, &
  MoldyMasa=AmuAtomica, &
  MoldyLongitud=AngstromAtomica, &
  MoldyTiempo=1D-13*SITiempo, &
  MoldyCarga=1.0D0, &
  MoldyEnergia=AmuAtomica*(AngstromAtomica/(1D-13*SITiempo))**2

CONTAINS
!LeerSistemaMoldy(Fich)
!LeerMoleculaMoldy(Fich,Nom,Mols,Dat)
!LeerControlMoldy(Fich)
!LeerConfigsMoldy(Centrar)
!CentrarConfigMoldy(Celda,Cent,Cuat,Centrar,U)
!ModificarMoldy(Ent,Sal)

!-------------------------------------------------------------------------------
! Esta subrutina lee el sistema que se va a estudiar a partir del fichero de
! entrada de moldy
!-------------------------------------------------------------------------------
! Fich:     Numero del fichero de entrada
! Error:    Variable para controlar los errores
! Linea:    La línea que se va leyendo
! Especie*: Las moléculas que se encuentran
! Nombre*:  El nombre de cada molécula
! Molec*:   El numero de moléculas de cada clase
! Renum:    Matriz para renumerar los átomos y que sean seguidos
! Aux:      Vector auxiliar para leer los parámetros
! i,j,k:    Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerSistemaMoldy(Fich)
  USE Parametros
  USE Sistema
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich

  TYPE (Atomo), DIMENSION(:), POINTER :: Especie1,Especie2,Especie3
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Aux
  INTEGER, DIMENSION(:), ALLOCATABLE :: Renum
  CHARACTER(LEN=LLL) :: Linea,Nombre1,Nombre2,Nombre3
  INTEGER :: Error,Molec1,Molec2,Molec3,i,j,k

  !Se lee la primera molécula
  CALL LeerMoleculaMoldy(Fich,Nombre1,Molec1,Especie1)
  !Se lee la segunda molécula
  CALL LeerMoleculaMoldy(Fich,Nombre2,Molec2,Especie2)
  !Se lee la tercera molécula
  CALL LeerMoleculaMoldy(Fich,Nombre3,Molec3,Especie3)

  !Se ponen los datos en su sitio ("Soluto" es la primera especie)
  IF (SIZE(Especie2,1) > 0) THEN
    IF (TRIM(Nombre1) == TRIM(Nombre2)) &
      CALL Mensaje('LeerSistemaMoldy',6,.TRUE.)
    IF (SIZE(Especie3,1) > 0) THEN
      IF (TRIM(Nombre1) == TRIM(Nombre3)) &
        CALL Mensaje('LeerSistemaMoldy',6,.TRUE.)
      IF (TRIM(Nombre2) == TRIM(Nombre3)) &
        CALL Mensaje('LeerSistemaMoldy',6,.TRUE.)
    END IF
  END IF
  IF (ALLOCATED(Soluto)) DEALLOCATE(Soluto)
  IF (ALLOCATED(Disolvente)) DEALLOCATE(Disolvente)
  IF (ALLOCATED(Disolvente2)) DEALLOCATE(Disolvente2)
  ALLOCATE(Soluto(SIZE(Especie1,1)),Disolvente(SIZE(Especie2,1)), &
    Disolvente2(SIZE(Especie3,1)))
  Soluto(:)=Especie1(:)
  Disolvente(:)=Especie2(:)
  Disolvente2(:)=Especie3(:)
  NombreSoluto=TRIM(Nombre1)
  NombreDisolvente=TRIM(Nombre2)
  NombreDisolvente2=TRIM(Nombre3)
  MoleculasSoluto=Molec1
  MoleculasDisolvente=Molec2
  MoleculasDisolvente2=Molec3
  DEALLOCATE(Especie1,Especie2,Especie3)

  !Pasa el resto hasta encontrar 'end'
  DO
    CALL LeerSiguienteLinea(Fich,Linea,Error)
    IF (Error /= 0) CALL Mensaje('LeerSistemaMoldy',5,.TRUE.)
    IF (Linea(1:3) == 'end') EXIT
  END DO

  !Se renumeran los átomos para que estén seguidos
  k=MAX(MAXVAL(Soluto(:)%id),MAXVAL(Disolvente(:)%id),MAXVAL(Disolvente2(:)%id))
  ALLOCATE(Renum(k))
  Renum(:)=0
  j=0
  k=0
  DO i=1,SIZE(Soluto,1)
    IF (Renum(Soluto(i)%id) == 0) THEN
      k=k+1
      Renum(Soluto(i)%id)=k
     ELSE
      IF (Renum(Soluto(i)%id) <= j) CALL Mensaje('LeerSistemaMoldy',30,.TRUE.)
    END IF
    Soluto(i)%id=Renum(Soluto(i)%id)
  END DO
  j=MAXVAL(Renum)
  DO i=1,SIZE(Disolvente,1)
    IF (Renum(Disolvente(i)%id) == 0) THEN
      k=k+1
      Renum(Disolvente(i)%id)=k
     ELSE
      IF (Renum(Disolvente(i)%id) <= j) &
        CALL Mensaje('LeerSistemaMoldy',30,.TRUE.)
    END IF
    Disolvente(i)%id=Renum(Disolvente(i)%id)
  END DO
  j=MAXVAL(Renum)
  DO i=1,SIZE(Disolvente2,1)
    IF (Renum(Disolvente2(i)%id) == 0) THEN
      k=k+1
      Renum(Disolvente2(i)%id)=k
     ELSE
      IF (Renum(Disolvente2(i)%id) <= j) &
        CALL Mensaje('LeerSistemaMoldy',30,.TRUE.)
    END IF
    Disolvente2(i)%id=Renum(Disolvente2(i)%id)
  END DO
  j=MAXVAL(Renum)

  !Lee el potencial interatómico
  CALL LeerSiguienteLinea(Fich,Linea,Error)
  IF (Error /= 0) CALL Mensaje('LeerSistemaMoldy',5,.TRUE.)
  SELECT CASE (TRIM(Linea))
   CASE ('lennard-jones')
    TipoPotencial=1
    i=2
   CASE ('generic')
    TipoPotencial=2
    i=6
   CASE DEFAULT
    CALL Mensaje('LeerSistemaMoldy',7,.TRUE.)
  END SELECT
  IF (ALLOCATED(InterAtom)) DEALLOCATE (InterAtom,QInter)
  ALLOCATE(QInter(j,j),InterAtom(j,j,i),Aux(i))
  QInter(:,:)=.FALSE.
  InterAtom(:,:,:)=0.0D0
  DO
    CALL LeerSiguienteLinea(Fich,Linea,Error)
    IF (Error /= 0) CALL Mensaje('LeerSistemaMoldy',5,.TRUE.)
    IF (Linea(1:3) == 'end') EXIT
    READ(Linea,*) i,j,(Aux(k),k=1,SIZE(Aux,1))
    IF ((Renum(i) == 0) .OR. (Renum(j) == 0) .OR. (MAX(i,j) > SIZE(Renum,1))) &
      CALL Mensaje('LeerSistemaMoldy',18,.TRUE.)
    i=Renum(i)
    j=Renum(j)
    QInter(i,j)=.TRUE.
    QInter(j,i)=.TRUE.
    InterAtom(i,j,:)=Aux(:)
    InterAtom(j,i,:)=InterAtom(i,j,:)
  END DO
  DEALLOCATE(Renum,Aux)

  !Se cambian las unidades
  SELECT CASE (TipoPotencial)
   CASE (1) !lennard-jones
    InterAtom(:,:,1)=InterAtom(:,:,1)*MoldyEnergia
    InterAtom(:,:,2)=InterAtom(:,:,2)*MoldyLongitud
   CASE (2) !generic
    InterAtom(:,:,1)=InterAtom(:,:,1)*MoldyEnergia
    InterAtom(:,:,2)=InterAtom(:,:,2)/MoldyLongitud
    InterAtom(:,:,3)=InterAtom(:,:,3)*MoldyEnergia*MoldyLongitud**12
    InterAtom(:,:,4)=InterAtom(:,:,4)*MoldyEnergia*MoldyLongitud**4
    InterAtom(:,:,5)=InterAtom(:,:,5)*MoldyEnergia*MoldyLongitud**6
    InterAtom(:,:,6)=InterAtom(:,:,6)*MoldyEnergia*MoldyLongitud**8
  END SELECT

END SUBROUTINE LeerSistemaMoldy

!-------------------------------------------------------------------------------
! Lee los datos de una molécula del fichero de entrada de moldy
!-------------------------------------------------------------------------------
! Fich:    Numero del fichero que se lee
! Nom:     Nombre de la molécula
! Mols:    Numero de moléculas de esta clase
! Dat:     Matriz donde se guardaran los datos
! Error:   Variable para controlar los errores
! Linea:   Línea que se va leyendo
! Simb:    Símbolo atómico
! Atom:    Número de átomos en la molécula
! UTmp:    Unidad temporal
! i,j,k:   Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerMoleculaMoldy(Fich,Nom,Mols,Dat)
  USE TipoAtomo
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich
  CHARACTER(LEN=*), INTENT(OUT) :: Nom
  INTEGER, INTENT(OUT) :: Mols
  TYPE(Atomo), DIMENSION(:), POINTER :: Dat

  CHARACTER(LEN=LLL) :: Linea
  CHARACTER(LEN=2) :: Simb
  INTEGER :: Error,Atom,UTmp,i,j,k

  !Lee el nombre y número de moléculas
  CALL LeerSiguienteLinea(Fich,Linea,Error)
  IF (Error /= 0) CALL Mensaje('LeerMoleculaMoldy',5,.TRUE.)
  !Si no hay más moléculas, se devuelve una molécula vacía
  IF (TRIM(Linea) == 'end') THEN
    Nom=''
    Mols=0
    ALLOCATE(Dat(0))
    BACKSPACE(Fich)
    RETURN
  END IF
  READ(Linea,*,IOSTAT=Error) Nom,Mols
  IF (Error /= 0) CALL Mensaje('LeerMoleculaMoldy',5,.TRUE.)

  !Cuenta el número de átomos
  Atom=0
  UTmp=NuevaUnidad()
  OPEN(UTmp,STATUS='SCRATCH',ACTION='READWRITE')
  DO
    CALL LeerSiguienteLinea(Fich,Linea,Error)
    IF (Error /= 0) CALL Mensaje('LeerMoleculaMoldy',5,.TRUE.)
    SELECT CASE (Linea(1:1))
     CASE ('0','1','2','3','4','5','6','7','8','9')
      Atom=Atom+1
      WRITE(UTmp,'(A)') TRIM(Linea)
     CASE DEFAULT
      EXIT
    END SELECT
  END DO
  BACKSPACE(Fich)

  !Lee los datos de los átomos
  REWIND(UTmp)
  ALLOCATE(Dat(Atom))
  DO i=1,Atom
    READ(UTmp,'(A)') Linea
    READ(Linea,*) Dat(i)%id
    k=0
    DO j=1,i-1
      IF (Dat(j)%id == Dat(i)%id) k=j
    END DO
    IF (k == 0) THEN
      READ(Linea,*) Dat(i)%id,Dat(i)%pos(:),Dat(i)%m,Dat(i)%q,Dat(i)%nom
      Dat(i)%m=Dat(i)%m*MoldyMasa
      Dat(i)%q=Dat(i)%q*MoldyCarga
     ELSE
      READ(Linea,*) Dat(i)%id,Dat(i)%pos(:)
      Dat(i)%m=Dat(k)%m
      Dat(i)%q=Dat(k)%q
      Dat(i)%nom=Dat(k)%nom
    END IF
    Dat(i)%pos(:)=Dat(i)%pos(:)*AngstromAtomica
    Simb=Dat(i)%nom(1:2)
    IF ((ICHAR(Simb(2:2)) < 97) .OR. (ICHAR(Simb(2:2)) > 122)) Simb(2:2)=' '
    Dat(i)%z=0
    DO j=0,SIZE(Simbolo,1)
      IF (Simb == Simbolo(j)) Dat(i)%z=j
    END DO
  END DO
  CLOSE(UTmp)

END SUBROUTINE LeerMoleculaMoldy

!-------------------------------------------------------------------------------
! Lee algunos datos del fichero de control de Moldy
!-------------------------------------------------------------------------------
! Fich:    Numero del fichero de control
! Linea:   Línea que se va leyendo
! Ind:     Posición del signo =
! Var:     Nombre de la variable que se lee
! Val:     Valor de la variable
! Error:   Variable para controlar los errores
!-------------------------------------------------------------------------------
SUBROUTINE LeerControlMoldy(Fich)
  USE Parametros
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich

  CHARACTER(LEN=LLL) :: Linea,Var,Val
  INTEGER :: Ind,Error,NSteps,BeginDump,DumpInterval

  !Inicia las variables
  DumpInterval=0

  !Lee las variables del fichero
  Error=0
  DO
    CALL LeerSiguienteLinea(Fich,Linea,Error)
    IF (Error /= 0) EXIT
    Ind=INDEX(Linea,'=')
    Var=Linea(:MAX(Ind-1,0))
    Val=Linea(MIN(Ind+1,LEN(Linea)):)
    SELECT CASE (TRIM(Var))
     CASE ('sys-spec-file')
      MoldyInput=TRIM(ADJUSTL(Val))
     CASE ('save-file')
      MoldySave=TRIM(ADJUSTL(Val))
     CASE ('temperature')
      READ(Val,*) MoldyTemp
     CASE ('mass-unit')
      READ(Val,*) MoldyMasa
      MoldyMasa=MoldyMasa*SIMasa
     CASE ('length-unit')
      READ(Val,*) MoldyLongitud
      MoldyLongitud=MoldyLongitud*SILongitud
     CASE ('time-unit')
      READ(Val,*) MoldyTiempo
      MoldyTiempo=MoldyTiempo*SITiempo
     CASE ('charge-unit')
      READ(Val,*) MoldyCarga
      MoldyCarga=MoldyCarga*SICarga
     CASE ('dump-file')
      MoldyDump=TRIM(ADJUSTL(Val))
     CASE ('backup-file')
      MoldyBackup=TRIM(ADJUSTL(Val))
     CASE ('nsteps')
      READ(Val,*) NSteps
     CASE ('begin-dump')
      READ(Val,*) BeginDump
     CASE ('dump-interval')
      READ(Val,*) DumpInterval
    END SELECT
  END DO

  !Calcula la unidad de energía en función de las otras unidades
  MoldyEnergia=MoldyMasa*(MoldyLongitud/MoldyTiempo)**2

  !Calcula el número de configuraciones en los dump
  IF (DumpInterval > 0) THEN
    MoldyConfigs=(NSteps-BeginDump)/DumpInterval+1
   ELSE
    MoldyConfigs=0
  END IF
  IF (NumConfig > MoldyConfigs) THEN
    NumConfig=MoldyConfigs
    CALL Mensaje('LeerControlMoldy',33,.FALSE.)
  END IF

END SUBROUTINE LeerControlMoldy

!-------------------------------------------------------------------------------
! Genera un fichero binario con las configuraciones del moldy
!-------------------------------------------------------------------------------
! Centrar:  Numero de la especie sobre la que se centran las configuraciones
!           (0: Soluto, 1: Disolvente, 2: Disolvente2...)
! Cent:     Centros de masas de las moléculas
! Cuat:     Cuaterniones de las moléculas poliatómicas
! Celda:    Vectores de la celda de simulación
! Paso:     Intervalo entre las configuraciones
! UL,UC,UQ: Ficheros de donde se leen las configuraciones
! Mols:     Número de moléculas
! Cuats:    Número de moléculas poliatómicas
! Conf:     Número de la configuración
! i,j:      Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerConfigsMoldy(Centrar)
  USE Configuraciones
  USE Parametros
  USE Utilidades
  USE Sistema
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Centrar

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Cent,Cuat
  DOUBLE PRECISION, DIMENSION(3,3) :: Celda
  DOUBLE PRECISION :: Paso
  INTEGER :: UL,UC,UQ,Mols,Cuats,Conf,i,j

  !Se calcula el número de moléculas y de cuaterniones
  Mols=0
  Cuats=0
  Mols=Mols+MoleculasSoluto
  IF (SIZE(Soluto,1) > 1) Cuats=Cuats+MoleculasSoluto
  Mols=Mols+MoleculasDisolvente
  IF (SIZE(Disolvente,1) > 1) Cuats=Cuats+MoleculasDisolvente
  Mols=Mols+MoleculasDisolvente2
  IF (SIZE(Disolvente2,1) > 1) Cuats=Cuats+MoleculasDisolvente2

  ALLOCATE(Cent(Mols,3),Cuat(Cuats,4))

  !Se abren los ficheros de donde se leen las configuraciones
  UL=NuevaUnidad()
  OPEN(UL,FILE='l.tmp',STATUS='OLD')
  UC=NuevaUnidad()
  OPEN(UC,FILE='c.tmp',STATUS='OLD')
  UQ=NuevaUnidad()
  OPEN(UQ,FILE='q.tmp',STATUS='OLD')

  !Se abre el fichero binario donde se escribe
  CALL AbrirUConf()

  !Se escriben todas las configuraciones que se piden
  Paso=DBLE(MoldyConfigs-1)/DBLE(MAX(NumConfig-1,1))
  Conf=0
  DO i=1,NumConfig
    !Se pasan configuraciones hasta llegar a la deseada
    DO WHILE (Conf < NINT((i-1)*Paso+1)-1)
      Conf=Conf+1
      READ(UL,*)
      READ(UC,*)
      READ(UQ,*)
    END DO
    Conf=Conf+1

    !Se leen los datos
    READ(UL,*) Celda(:,:)
    Celda(:,:)=Celda(:,:)*AngstromAtomica
    READ(UC,*) (Cent(j,:),j=1,Mols)
    Cent(:,:)=Cent(:,:)*AngstromAtomica
    READ(UQ,*) (Cuat(j,:),j=1,Cuats)

    !Se escribe la configuracion centrada sucesivamente en cada molécula de
    !la especie que se indique
    SELECT CASE (Centrar)
     CASE (1)
      DO j=1,MoleculasSoluto
        WRITE(UConf) Conf
        CALL CentrarConfigMoldy(Celda(:,:),Cent(:,:),Cuat(:,:),j,UConf)
      END DO
     CASE (2)
      DO j=MoleculasSoluto+1,MoleculasSoluto+MoleculasDisolvente
        WRITE(UConf) Conf
        CALL CentrarConfigMoldy(Celda(:,:),Cent(:,:),Cuat(:,:),j,UConf)
      END DO
     CASE (3)
      DO j=1,MoleculasSoluto+MoleculasDisolvente+1, &
             MoleculasSoluto+MoleculasDisolvente+MoleculasDisolvente2
        WRITE(UConf) Conf
        CALL CentrarConfigMoldy(Celda(:,:),Cent(:,:),Cuat(:,:),j,UConf)
      END DO
    END SELECT
    Conf=Conf+1
  END DO

  !Se borran los ficheros
  CLOSE(UL,STATUS='DELETE')
  CLOSE(UC,STATUS='DELETE')
  CLOSE(UQ,STATUS='DELETE')

  DEALLOCATE(Cent,Cuat)

END SUBROUTINE LeerConfigsMoldy

!-------------------------------------------------------------------------------
! Centra una configuración en una molécula y la escribe en un fichero binario
!-------------------------------------------------------------------------------
! Celda:       Vectores de la celda de simulación
! Cent:        Centros de masa de las moléculas
! Cuat:        Cuaterniones de las moléculas poliatómicas
! Centrar:     Molécula sobre la que se centrara la configuración
! U:           Fichero donde se escribe
! CentroMasa:  Centros de masa modificados
! Cuaternion:  Cuaterniones modificados
! CoordSol:    Coordenadas de la molécula prototipo de soluto
! CoordDis:    Coordenadas de la molécula prototipo de disolvente
! CoordDis2:   Coordenadas de la molécula prototipo de disolvente 2
! CentroCelda: Coordenadas originales del centro de la celda
! CuatCelda:   Cuaternión original del centro de la celda
! Mols:        Número de moléculas
! i,j:         Contadores
!-------------------------------------------------------------------------------
SUBROUTINE CentrarConfigMoldy(Celda,Cent,Cuat,Centrar,U)
  USE Utilidades
  USE Sistema
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Cent,Cuat
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: Celda
  INTEGER, INTENT(IN) :: Centrar,U

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CentroMasa,Cuaternion, &
                                                   CoordSol,CoordDis,CoordDis2
  DOUBLE PRECISION, DIMENSION(3) :: CentroCelda
  DOUBLE PRECISION, DIMENSION(4) :: CuatCelda
  INTEGER :: Mols,i,j

  ALLOCATE(CentroMasa(SIZE(Cent,1),3),Cuaternion(SIZE(Cent,1),4))

  !Se copian los datos iniciales en las matrices de trabajo
  !Se añaden cuaterniones unidad a las moléculas monoatómicas
  Mols=0
  i=0
  CentroMasa(:,:)=Cent(:,:)
  IF (SIZE(Soluto,1) > 1) THEN
    Cuaternion(Mols+1:Mols+MoleculasSoluto,:)=Cuat(i+1:i+MoleculasSoluto,:)
    i=i+MoleculasSoluto
   ELSE
    Cuaternion(Mols+1:Mols+MoleculasSoluto,:)= &
      SPREAD((/1.0D0,0.0D0,0.0D0,0.0D0/),DIM=1,NCOPIES=MoleculasSoluto)
  END IF
  Mols=Mols+MoleculasSoluto
  IF (SIZE(Disolvente,1) > 1) THEN
    Cuaternion(Mols+1:Mols+MoleculasDisolvente,:)= &
      Cuat(i+1:i+MoleculasDisolvente,:)
    i=i+MoleculasDisolvente
   ELSE
    Cuaternion(Mols+1:Mols+MoleculasDisolvente,:)= &
      SPREAD((/1.0D0,0.0D0,0.0D0,0.0D0/),DIM=1,NCOPIES=MoleculasDisolvente)
  END IF
  Mols=Mols+MoleculasDisolvente
  IF (SIZE(Disolvente2,1) > 1) THEN
    Cuaternion(Mols+1:Mols+MoleculasDisolvente2,:)= &
      Cuat(i+1:i+MoleculasDisolvente2,:)
    i=i+MoleculasDisolvente2
   ELSE
    Cuaternion(Mols+1:Mols+MoleculasDisolvente2,:)= &
      SPREAD((/1.0D0,0.0D0,0.0D0,0.0D0/),DIM=1,NCOPIES=MoleculasDisolvente2)
  END IF
  Mols=Mols+MoleculasDisolvente2

  !Se crean las matrices de coordenadas de soluto y disolvente
  ALLOCATE(CoordSol(SIZE(Soluto,1),3),CoordDis(SIZE(Disolvente,1),3), &
           CoordDis2(SIZE(Disolvente2,1),3))
  DO i=1,SIZE(CoordSol,1)
    CoordSol(i,:)=Soluto(i)%pos(:)
  END DO
  DO i=1,SIZE(CoordDis,1)
    CoordDis(i,:)=Disolvente(i)%pos(:)
  END DO
  DO i=1,SIZE(CoordDis2,1)
    CoordDis2(i,:)=Disolvente2(i)%pos(:)
  END DO

  !Si Centrar=0, no se modifica la configuracion
  IF ((Centrar < 1) .OR. (Centrar > Mols)) THEN
    CentroCelda(:)=0.0D0
    CuatCelda(:)=(/-1.0D0,0.0D0,0.0D0,0.0D0/)
   ELSE
    CentroCelda(:)=CentroMasa(Centrar,:)
    CuatCelda(:)=Cuaternion(Centrar,:)
    !Se invierte el cuaternión para obtener la rotación inversa
    CuatCelda(1)=-CuatCelda(1)
  END IF

!$OMP PARALLEL DO PRIVATE(i,j)
  DO i=1,Mols
    !Se desplaza cada molécula para que entre en la celda unidad
    DO j=1,3
      CentroMasa(i,j)=CentroMasa(i,j)-CentroCelda(j)- &
        Celda(j,j)*NINT((CentroMasa(i,j)-CentroCelda(j))/Celda(j,j))
    END DO
    !Se modifica el cuaternión de cada molécula para considerar la rotación de
    !la celda
    Cuaternion(i,:)=ProductoCuat(CuatCelda(:),Cuaternion(i,:))
  END DO
!$OMP END PARALLEL DO
  !Se giran los centros de masa según la rotación de la celda
  CentroMasa(:,:)=RotarCuaternion(CentroMasa(:,:),CuatCelda(:))

  !Se escriben los datos del centro de la celda
  WRITE(U) Centrar
  WRITE(U) Celda(:,:)
  WRITE(U) CentroCelda(:)
  WRITE(U) -CuatCelda(1),CuatCelda(2:)
  !Y los de todas las moléculas
  Mols=0
  DO i=1,MoleculasSoluto
    WRITE(U) Norma(CentroMasa(Mols+i,:))
    WRITE(U) CentroMasa(Mols+i,:)
    WRITE(U) Cuaternion(Mols+i,:)
    WRITE(U) SPREAD(CentroMasa(Mols+i,:),DIM=1,NCOPIES=SIZE(CoordSol,1))+ &
               RotarCuaternion(CoordSol(:,:),Cuaternion(Mols+i,:))
  END DO
  Mols=Mols+MoleculasSoluto
  DO i=1,MoleculasDisolvente
    WRITE(U) Norma(CentroMasa(Mols+i,:))
    WRITE(U) CentroMasa(Mols+i,:)
    WRITE(U) Cuaternion(Mols+i,:)
    WRITE(U) SPREAD(CentroMasa(Mols+i,:),DIM=1,NCOPIES=SIZE(CoordDis,1))+ &
               RotarCuaternion(CoordDis(:,:),Cuaternion(Mols+i,:))
  END DO
  Mols=Mols+MoleculasDisolvente
  DO i=1,MoleculasDisolvente2
    WRITE(U) Norma(CentroMasa(Mols+i,:))
    WRITE(U) CentroMasa(Mols+i,:)
    WRITE(U) Cuaternion(Mols+i,:)
    WRITE(U) SPREAD(CentroMasa(Mols+i,:),DIM=1,NCOPIES=SIZE(CoordDis2,1))+ &
               RotarCuaternion(CoordDis2(:,:),Cuaternion(Mols+i,:))
  END DO
  Mols=Mols+MoleculasDisolvente2

  DEALLOCATE(CentroMasa,Cuaternion,CoordSol,CoordDis,CoordDis2)

END SUBROUTINE CentrarConfigMoldy

!-------------------------------------------------------------------------------
! Escribe una nueva entrada de moldy, incluyendo fichero de control y de sistema
!-------------------------------------------------------------------------------
! Ent:      Unidad del fichero de control de entrada
! Sal:      Unidad del fichero de control de salida
! Linea:    Línea que se va leyendo
! USis:     Unidad del nuevo fichero de sistema de moldy
! Repetido: Variable que determina si el tipo de átomo ya ha sido usado
! Error:    Variable para controlar los errores
! i,j:      Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ModificarMoldy(Ent,Sal)
  USE Parametros
  USE Sistema
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Ent,Sal

  CHARACTER(LEN=LLL) :: Linea
  INTEGER :: i,j,Error,USis,Repetido

  DO
    READ(Ent,'(A)',IOSTAT=Error) Linea
    IF (Error /= 0) EXIT
    IF (INDEX(Linea,'sys-spec-file') /= 0) THEN
      WRITE(Sal,'(A)') TRIM(Linea)//TRIM(Extension)
     ELSE IF (INDEX(Linea,'save-file') /= 0) THEN
      WRITE(Sal,'(A)') TRIM(Linea)//TRIM(Extension)
     ELSE IF (INDEX(Linea,'dump-file') /= 0) THEN
      WRITE(Sal,'(A)') TRIM(Linea)//TRIM(Extension)//'-'
     ELSE IF (INDEX(Linea,'backup-file') /= 0) THEN
      WRITE(Sal,'(A)') TRIM(Linea)//TRIM(Extension)
     ELSE
      WRITE(Sal,'(A)') TRIM(Linea)
    END IF
  END DO

  USis=NuevaUnidad()
  OPEN(USis,FILE=TRIM(MoldyInput)//TRIM(Extension),STATUS='REPLACE', &
       ACTION='WRITE')

  !Escribe la molécula de soluto
  IF (MoleculasSoluto > 0) THEN
    WRITE(USis,100) TRIM(NombreSoluto),MoleculasSoluto
    DO i=1,SIZE(Soluto,1)
      Repetido=0
      DO j=1,i-1
        IF (Soluto(j)%id == Soluto(i)%id) THEN
          Repetido=j
          EXIT
        END IF
      END DO
      IF (Repetido == 0) THEN
        WRITE(USis,101) Soluto(i)%id,Soluto(i)%pos(:)/MoldyLongitud, &
                     Soluto(i)%m/MoldyMasa,Soluto(i)%q/MoldyCarga, &
                     TRIM(Soluto(i)%nom)
       ELSE
        WRITE(USis,101) Soluto(i)%id,Soluto(i)%pos(:)/MoldyLongitud
      END IF
    END DO
  END IF

  !Escribe la molécula de disolvente
  IF (MoleculasDisolvente > 0) THEN
    WRITE(USis,100) TRIM(NombreDisolvente),MoleculasDisolvente
    DO i=1,SIZE(Disolvente,1)
      Repetido=0
      DO j=1,i-1
        IF (Disolvente(j)%id == Disolvente(i)%id) THEN
          Repetido=j
          EXIT
        END IF
      END DO
      IF (Repetido == 0) THEN
        WRITE(USis,101) Disolvente(i)%id,Disolvente(i)%pos(:)/MoldyLongitud, &
                     Disolvente(i)%m/MoldyMasa,Disolvente(i)%q/MoldyCarga, &
                     TRIM(Disolvente(i)%nom)
       ELSE
        WRITE(USis,101) Disolvente(i)%id,Disolvente(i)%pos(:)/MoldyLongitud
      END IF
    END DO
  END IF

  !Escribe la molécula del segundo disolvente
  IF (MoleculasDisolvente2 > 0) THEN
    WRITE(USis,100) TRIM(NombreDisolvente2),MoleculasDisolvente2
    DO i=1,SIZE(Disolvente2,1)
      Repetido=0
      DO j=1,i-1
        IF (Disolvente2(j)%id == Disolvente2(i)%id) THEN
          Repetido=j
          EXIT
        END IF
      END DO
      IF (Repetido == 0) THEN
        WRITE(USis,101) Disolvente2(i)%id,Disolvente2(i)%pos(:)/MoldyLongitud, &
                     Disolvente2(i)%m/MoldyMasa,Disolvente2(i)%q/MoldyCarga, &
                     TRIM(Disolvente2(i)%nom)
       ELSE
        WRITE(USis,101) Disolvente2(i)%id,Disolvente2(i)%pos(:)/MoldyLongitud
      END IF
    END DO
  END IF

  WRITE(USis,'(A)') 'end'

  !Escribe el potencial de vdW
  SELECT CASE (TipoPotencial)
   CASE (1)
    WRITE(USis,'(A)') 'lennard-jones'
    DO i=1,SIZE(InterAtom,1)
      DO j=1,i
        IF (QInter(i,j)) WRITE(USis,102) i,j,InterAtom(i,j,1)/MoldyEnergia, &
                         InterAtom(i,j,2)/MoldyLongitud
      END DO
    END DO
   CASE (2)
    WRITE(USis,'(A)') 'generic'
    DO i=1,SIZE(InterAtom,1)
      DO j=1,i
        IF (QInter(i,j)) WRITE(USis,102) i,j,InterAtom(i,j,1)/MoldyEnergia, &
                         InterAtom(i,j,2)*MoldyLongitud, &
                         InterAtom(i,j,3)/MoldyEnergia/MoldyLongitud**12, &
                         InterAtom(i,j,4)/MoldyEnergia/MoldyLongitud**4, &
                         InterAtom(i,j,5)/MoldyEnergia/MoldyLongitud**6, &
                         InterAtom(i,j,6)/MoldyEnergia/MoldyLongitud**8
      END DO
    END DO
  END SELECT
  WRITE(USis,'(A)') 'end'

  CLOSE(USis)

100 FORMAT (A,1X,I3)
101 FORMAT (1X,I3,1X,3(F12.8,1X),F7.4,1X,F12.8,1X,A)
102 FORMAT (2(1X,I3),6(F14.6,1X))

END SUBROUTINE ModificarMoldy

END MODULE Moldy
