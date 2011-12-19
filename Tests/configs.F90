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

PROGRAM prueba
  USE Parametros
  USE Entrada
  USE Sistema
  USE Moldy
  USE GenericoMM
  USE DatosQMMM
  USE Configuraciones
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  INTEGER :: U,Error

  CALL LeerEntrada(5)

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),STATUS='OLD',ACTION='READ')
    CALL LeerSistemaGenerico(U)
    CLOSE(U)

   CASE (1) !Moldy
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),STATUS='OLD',ACTION='READ')
    CALL LeerControlMoldy(U)
    CLOSE(U)

    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(MoldyInput),STATUS='OLD',ACTION='READ')
    CALL LeerSistemaMoldy(U)
    CLOSE(U)
  END SELECT

  CALL LeerConfs()
  CALL AbrirUConf()
  DO
    CALL LeerConfig(Error)
    IF (Error /= 0) EXIT
    CALL EscribirXYZ(6,Centros=.TRUE.)
  END DO

CONTAINS

SUBROUTINE LeerConfs()
  IMPLICIT NONE
  CHARACTER (LEN=LLL) :: Linea,Dump
  INTEGER :: NumMol,NumCuat

  Dump=TRIM(MoldyDump)

  NumMol=0
  NumCuat=0
  NumMol=NumMol+MoleculasSoluto
  IF (SIZE(Soluto,1) > 1) NumCuat=NumCuat+MoleculasSoluto
  NumMol=NumMol+MoleculasDisolvente
  IF (SIZE(Disolvente,1) > 1) NumCuat=NumCuat+MoleculasDisolvente
  NumMol=NumMol+MoleculasDisolvente2
  IF (SIZE(Disolvente2,1) > 1) NumCuat=NumCuat+MoleculasDisolvente2

  SELECT CASE (ProgramaMM)
   CASE (0) !Generico
    CALL LeerConfigsGenerico()

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

END SUBROUTINE LeerConfs

END PROGRAM prueba
