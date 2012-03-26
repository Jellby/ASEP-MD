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

MODULE Ejecutar

#ifndef LLL
#define LLL 256
#endif

#ifdef __NAG__
USE F90_UNIX_PROC
#endif

USE GenericoQM
USE Gaussian
USE Molcas
USE GenericoMM
USE Moldy

CONTAINS
!EjecutarQM(Der)
!EjecutarMM(Ejec)

!-------------------------------------------------------------------------------
! Ejecuta el programa de cálculo cuántico y lee los resultados
!-------------------------------------------------------------------------------
! Der:                  Orden de las derivadas a calcular
! Linea:                Línea que se ejecuta para el cálculo
! Ent:                  Fichero de entrada
! Sal:                  Fichero de salida
! Car:                  Fichero de cargas externas
! UEnt,USal,UFchk,UCar: Unidades de los ficheros
!-------------------------------------------------------------------------------
SUBROUTINE EjecutarQM(Der)
  USE Parametros
  USE Utilidades
  USE Configuraciones
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Der

  CHARACTER(LEN=LLL) :: Linea,Ent,Sal,Car
  INTEGER :: UEnt,USal,UFchk,UCar

  IF (.NOT. ALLOCATED(DisolvQM)) ALLOCATE(DisolvQM(0,5))
  IF (.NOT. ALLOCATED(PotQM)) ALLOCATE(PotQM(0,4))

  !Los nombres de los ficheros
  Ent=TRIM(EntradaQM)//TRIM(Extension)
  Sal=TRIM(SalidaQM)//TRIM(Extension)
  IF (TRIM(CargasExternas) /= '') THEN
    Car=TRIM(CargasExternas)//TRIM(Extension)
   ELSE
    Car=''
  END IF

  !Según el programa utilizado, se construye la entrada, se ejecuta el programa
  !y se lee la salida correspondiente
  SELECT CASE (ProgramaQM)
   CASE (0) !Genérico
    UEnt=NuevaUnidad()
    OPEN(UEnt,FILE=TRIM(Ent),STATUS='REPLACE',ACTION='WRITE')
    IF (TRIM(Car) /= '') THEN
      UCar=NuevaUnidad()
      OPEN(UCar,FILE=TRIM(Car),STATUS='REPLACE',ACTION='WRITE')
      CALL EntradaGenericoQM(UEnt,Der,UCar)
      CLOSE(UCar)
     ELSE
      CALL EntradaGenericoQM(UEnt,Der)
    END IF
    CLOSE(UEnt)

    Linea=TRIM(EjecutableQM)//' '//TRIM(Ent)//' '//TRIM(Sal)// &
          ' '//TRIM(Extension)
    CALL SYSTEM(TRIM(Linea))

    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Sal),STATUS='OLD',ACTION='READ')
    CALL LeerSalidaGenericoQM(USal)
    CLOSE(USal)

   CASE (1) !Gaussian
    UEnt=NuevaUnidad()
    OPEN(UEnt,FILE=TRIM(EntradaQM),STATUS='OLD',ACTION='READ')
    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Ent),STATUS='REPLACE',ACTION='WRITE')
    IF (TRIM(Car) /= '') THEN
      UCar=NuevaUnidad()
      OPEN(UCar,FILE=TRIM(Car),STATUS='REPLACE',ACTION='WRITE')
      CALL ModificarGaussian(UEnt,USal,Der,UCar)
      CLOSE(UCar)
     ELSE
      CALL ModificarGaussian(UEnt,USal,Der)
    END IF
    CLOSE(UEnt)
    CLOSE(USal)

    Linea=TRIM(EjecutableQM)//' '//TRIM(Ent)//' '//TRIM(Sal)// &
          ' '//TRIM(Extension)
    CALL SYSTEM(TRIM(Linea))

    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Sal),STATUS='OLD',ACTION='READ')
    UFchk=NuevaUnidad()
    OPEN(UFchk,FILE=TRIM(FChkGaussian)//TRIM(Extension),STATUS='OLD',ACTION='READ')
    CALL LeerSalidaGaussian(USal,UFchk)
    CLOSE(USal)
    CLOSE(UFchk)

   CASE (2) !Molcas
    UEnt=NuevaUnidad()
    OPEN(UEnt,FILE=TRIM(EntradaQM),STATUS='OLD',ACTION='READ')
    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Ent),STATUS='REPLACE',ACTION='WRITE')
    IF (TRIM(Car) /= '') THEN
      UCar=NuevaUnidad()
      OPEN(UCar,FILE=TRIM(Car),STATUS='REPLACE',ACTION='WRITE')
      CALL ModificarMolcas(UEnt,USal,Der,UCar)
      CLOSE(UCar)
     ELSE
      CALL ModificarMolcas(UEnt,USal,Der)
    END IF
    CLOSE(UEnt)
    CLOSE(USal)

    Linea=TRIM(EjecutableQM)//' '//TRIM(Ent)//' '//TRIM(Sal)// &
          ' '//TRIM(Extension)
    CALL SYSTEM(TRIM(Linea))

    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Sal),STATUS='OLD',ACTION='READ')
    CALL LeerSalidaMolcas(USal)
    CLOSE(USal)
  END SELECT

  CALL Promedios(Der,MolQM)

END SUBROUTINE EjecutarQM

!-------------------------------------------------------------------------------
! Ejecuta el programa de dinámica molecular y lee las configuraciones
!-------------------------------------------------------------------------------
! Linea:        Línea que se ejecuta para el cálculo
! Ent,Sal,Dump: Ficheros de entrada, salida y dumps
! NumMol:       Número de moléculas
! NumCuat:      Número de moléculas poliatómicas
! UEnt,USal:    Unidades de los ficheros
!-------------------------------------------------------------------------------
SUBROUTINE EjecutarMM(Ejec)
  USE Parametros
  USE Sistema
  USE Utilidades
  IMPLICIT NONE
  LOGICAL, INTENT(IN) :: Ejec

  CHARACTER(LEN=LLL) :: Linea,Ent,Sal,Dump
  INTEGER :: UEnt,USal,NumMol,NumCuat

  Ent=TRIM(EntradaMM)//TRIM(Extension)
  Sal=TRIM(SalidaMM)//TRIM(Extension)
  Dump=TRIM(MoldyDump)//TRIM(Extension)//'-'

  !Se construye la entrada para la dinámica molecular
  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Ent),STATUS='REPLACE',ACTION='WRITE')
    CALL EntradaGenericoMM(USal)
    CLOSE(USal)

   CASE (1) !Moldy
    UEnt=NuevaUnidad()
    OPEN(UEnt,FILE=TRIM(EntradaMM),STATUS='OLD',ACTION='READ')
    USal=NuevaUnidad()
    OPEN(USal,FILE=TRIM(Ent),STATUS='REPLACE',ACTION='WRITE')
    CALL ModificarMoldy(UEnt,USal)
    CLOSE(UEnt)
    CLOSE(USal)
  END SELECT

  IF (.NOT. Ejec) RETURN

  Linea=TRIM(EjecutableMM)//' '//TRIM(Ent)//' '//TRIM(Sal)// &
        ' '//TRIM(Extension)//' '//TRIM(TrayectoriaMM)
  CALL SYSTEM(TRIM(Linea))

  !Se calcula el número de moléculas y de cuaterniones
  NumMol=0
  NumCuat=0
  NumMol=NumMol+MoleculasSoluto
  IF (SIZE(Soluto,1) > 1) NumCuat=NumCuat+MoleculasSoluto
  NumMol=NumMol+MoleculasDisolvente
  IF (SIZE(Disolvente,1) > 1) NumCuat=NumCuat+MoleculasDisolvente
  NumMol=NumMol+MoleculasDisolvente2
  IF (SIZE(Disolvente2,1) > 1) NumCuat=NumCuat+MoleculasDisolvente2

  !Se leen las configuraciones de la dinámica
  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    CALL LeerConfigsGenerico

   CASE (1) !Moldy
    WRITE(Linea,100) &
      TRIM(DumpextMoldy),NumMol,NumCuat,1,TRIM(Dump),'c.tmp'
    CALL SYSTEM(TRIM(Linea))
    WRITE(Linea,100) &
      TRIM(DumpextMoldy),NumMol,NumCuat,2,TRIM(Dump),'q.tmp'
    CALL SYSTEM(TRIM(Linea))
    WRITE(Linea,100) &
      TRIM(DumpextMoldy),NumMol,NumCuat,3,TRIM(Dump),'l.tmp'
    CALL SYSTEM(TRIM(Linea))

    CALL LeerConfigsMoldy(1)
  END SELECT

100 FORMAT(A,' -R ',I5,' -Q ',I5,' -c ',I1,' ',A,'* > ',A)

END SUBROUTINE EjecutarMM

END MODULE Ejecutar
