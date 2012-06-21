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
  USE Optimizacion
  USE Utilidades
  IMPLICIT NONE
  INTEGER :: U

  CALL LeerEntrada(5)

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM))
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

  Extension=''
  CALL OptimizarGeometria(Soluto,6)

END PROGRAM prueba
