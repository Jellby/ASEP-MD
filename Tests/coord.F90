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

PROGRAM prueba
  USE Parametros
  USE Sistema
  USE Entrada
  USE GenericoMM
  USE Moldy
  USE Coordenadas
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  LOGICAL :: ModCoord
  INTEGER :: U,i

  CALL LeerEntrada(5)

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaGenerico(U)
    CLOSE(U)

   CASE (1) !Moldy
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),ACTION='READ',STATUS='OLD')
    CALL LeerControlMoldy(U)
    CLOSE(U)

    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(MoldyInput),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaMoldy(U)
    CLOSE(U)
  END SELECT

  TipoCoord=TipoCoordenadas
  CALL DefinirCoordenadas(Soluto)
  INQUIRE(FILE=TRIM(AltCoordenadas),EXIST=ModCoord)
  IF (ModCoord) THEN
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(AltCoordenadas))
    CALL ModificarCoordenadas(U)
    CLOSE(U)
  END IF
  CALL ConvertirCoordenadas(Cartesianas(Soluto),1)

  SELECT CASE (TipoCoord)
   CASE (0)
    WRITE(6,*) 'Coordenadas cartesianas'
   CASE (1)
    WRITE(6,*) 'Coordenadas cartesianas ponderadas'
   CASE (2)
    WRITE(6,*) 'Coordenadas internas'
  END SELECT
  WRITE(6,*) '( Angstrom, grados, amu )'
  WRITE(6,100) SIZE(Geometria,1)

  DO i=1,SIZE(DefCoord,1)
    SELECT CASE (TipoCoord)
     CASE (0)
      WRITE(6,102) (i-1)/3+1,MOD(i-1,3)+1,Geometria(i)/AngstromAtomica
     CASE (1)
      WRITE(6,102) (i-1)/3+1,MOD(i-1,3)+1,Geometria(i)/AngstromAtomica/SQRT(AmuAtomica)
     CASE (2)
      IF (DefCoord(i,3) == 0) THEN
        WRITE(6,101) DefCoord(i,:),Geometria(i)/AngstromAtomica
       ELSE
        WRITE(6,101) DefCoord(i,:),Geometria(i)/Grado
      END IF
    END SELECT
  END DO

100 FORMAT('N. coordenadas:',1X,I4)
101 FORMAT(4(1X,I4),F12.6)
102 FORMAT(2(1X,I4),F12.6)

END PROGRAM prueba
