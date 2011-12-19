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
  USE Malla
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE

  INTEGER :: U,i

  CALL LeerEntrada(5)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaMM))
  CALL LeerControlMoldy(U)
  CLOSE(U)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(MoldyInput))
  CALL LeerSistemaMoldy(U)
  CLOSE(U)

  CALL OrientarMolecula(Soluto,1)

  CALL MallaMolecula(Soluto,MallaSoluto)

  WRITE(6,100) SIZE(MallaSoluto,1)
  WRITE(6,*)
  DO i=1,SIZE(MallaSoluto,1)
    WRITE(6,101) 'X',MallaSoluto(i,:)
  END DO

100 FORMAT(I6)
101 FORMAT(A2,1X,3(F20.12,1X))

END PROGRAM prueba
