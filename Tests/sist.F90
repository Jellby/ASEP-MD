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
  USE GenericoMM
  USE Moldy
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  INTEGER :: U,i,j

  !READ(5,*) Fich
  CALL LeerEntrada(5)

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM))
    CALL LeerSistemaGenerico(U)
    CLOSE(U)

   CASE (1) !Moldy
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM))
    CALL LeerControlMoldy(U)
    CLOSE(U)

    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(MoldyInput))
    CALL LeerSistemaMoldy(U)
    CLOSE(U)
  END SELECT

  CALL OrientarMolecula(Soluto,1)
  CALL OrientarMolecula(Disolvente,1)
  CALL OrientarMolecula(Disolvente2,1)

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    WRITE(6,100) NombreSoluto,MoleculasSoluto
    DO i=1,SIZE(Soluto)
      WRITE(6,101) Soluto(i)%id,Soluto(i)%pos(:), &
                   Soluto(i)%m/AmuAtomica,Soluto(i)%q, &
                   Soluto(i)%nom,Soluto(i)%z
    END DO
   CASE (1) !Moldy
    WRITE(6,100) NombreSoluto,MoleculasSoluto
    DO i=1,SIZE(Soluto)
      WRITE(6,101) Soluto(i)%id,Soluto(i)%pos(:), &
                   Soluto(i)%m/MoldyMasa,Soluto(i)%q/MoldyCarga, &
                   Soluto(i)%nom,Soluto(i)%z
    END DO
  END SELECT

  IF (MoleculasDisolvente > 0) THEN
    SELECT CASE (ProgramaMM)
     CASE (0) !Genérico
      WRITE(6,100) NombreDisolvente,MoleculasDisolvente
      DO i=1,SIZE(Disolvente)
        WRITE(6,104) Disolvente(i)%id,Disolvente(i)%q,Disolvente(i)%nom
      END DO
     CASE (1) !Moldy
      WRITE(6,100) NombreDisolvente,MoleculasDisolvente
      DO i=1,SIZE(Disolvente)
        WRITE(6,101) Disolvente(i)%id,Disolvente(i)%pos(:), &
                     Disolvente(i)%m/MoldyMasa,Disolvente(i)%q/MoldyCarga, &
                     Disolvente(i)%nom,Disolvente(i)%z
      END DO
    END SELECT
  END IF

  IF (MoleculasDisolvente2 > 0) THEN
    WRITE(6,100) NombreDisolvente2,MoleculasDisolvente2
    DO i=1,SIZE(Disolvente2)
      WRITE(6,101) Disolvente2(i)%id,Disolvente2(i)%pos(:), &
                   Disolvente2(i)%m/MoldyMasa,Disolvente2(i)%q/MoldyCarga, &
                   Disolvente2(i)%nom,Disolvente2(i)%z
    END DO
  END IF

  IF (ProgramaMM==1) THEN
    WRITE(6,'(A)') 'end'
  END IF

  SELECT CASE (TipoPotencial)
   CASE (1)
    WRITE(6,'(A)') 'lennard-jones'
   CASE (2)
    WRITE(6,'(A)') 'generic'
  END SELECT
  DO i=1,SIZE(InterAtom,1)
    DO j=i,SIZE(InterAtom,2)
      IF (.NOT. QInter(i,j)) CYCLE
      SELECT CASE (ProgramaMM)
       CASE (0) !Genérico
        WRITE(6,102) i,j,InterAtom(i,j,:)
       CASE (1) !Moldy
        SELECT CASE (TipoPotencial)
         CASE (1)
          WRITE(6,102) i,j,InterAtom(i,j,1)/MoldyEnergia, &
                           InterAtom(i,j,2)/MoldyLongitud
         CASE (2)
          WRITE(6,102) i,j,InterAtom(i,j,1)/MoldyEnergia, &
                           InterAtom(i,j,2)*MoldyLongitud, &
                           InterAtom(i,j,3)/(MoldyEnergia*MoldyLongitud**12), &
                           InterAtom(i,j,4)/(MoldyEnergia*MoldyLongitud**4), &
                           InterAtom(i,j,5)/(MoldyEnergia*MoldyLongitud**6), &
                           InterAtom(i,j,6)/(MoldyEnergia*MoldyLongitud**8)
        END SELECT
      END SELECT
    END DO
  END DO

  IF (ProgramaMM==1) THEN
    WRITE(6,'(A)') 'end'
  END IF

  IF (ProgramaMM==1) THEN
    WRITE(6,*)
    WRITE(6,103) 'Unidad de masa:',    MoldyMasa,    'me',MoldyMasa/AmuAtomica,'amu'
    WRITE(6,103) 'Unidad de longitud:',MoldyLongitud,'a0',MoldyLongitud/AngstromAtomica,'A'
    WRITE(6,103) 'Unidad de tiempo:',  MoldyTiempo,  'hbar/Eh',1.0D15*MoldyTiempo/SITiempo,'fs'
    WRITE(6,103) 'Unidad de carga:',   MoldyCarga,   'e',MoldyCarga,'e'
    WRITE(6,103) 'Unidad de energia:', MoldyEnergia, 'Eh',MoldyEnergia/KcalmolAtomica,'kcal/mol'
  END IF

100 FORMAT(A16,1X,I5)
101 FORMAT(I3,1X,3(F10.6,1X),F7.4,1X,F9.6,1X,A16,1X,I3)
102 FORMAT(2(I3,1X),8(F8.4,:,1X))
103 FORMAT(A,T23,F12.6,1X,A,T44,'(',ES12.6,1X,A,')')
104 FORMAT(I3,1X,F9.6,1X,A16)

END PROGRAM prueba
