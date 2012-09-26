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

MODULE Malla

CONTAINS

!-------------------------------------------------------------------------------
! Calcula una malla o red de puntos en la zona ocupada por una molécula
! (definida por esferas)
!-------------------------------------------------------------------------------
! Mol:       Molécula para la que se calcula la malla
! MallaMol:  Malla (matriz de puntos) calculada
! Pmin,Pmax: Coordenadas mínima y máxima de la zona ocupada por la molécula
! Pun:       Coordenadas del punto que se va considerando
! Dist:      Distancia de la red cúbica base
! Lim:       Índices mínimo y máximo de los puntos de la malla
! U:         Fichero temporal para guardar los puntos
! Num:       Número de puntos en la malla
! i,j,k,l:   Contadores
!-------------------------------------------------------------------------------
SUBROUTINE MallaMolecula(Mol,MallaMol)
  USE TipoAtomo
  USE Parametros
  USE Unidades
  USE Sistema
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(IN) :: Mol
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: MallaMol

  DOUBLE PRECISION, DIMENSION(3) :: PMin,Pmax,Pun
  DOUBLE PRECISION :: Dist
  INTEGER, DIMENSION(3,2) :: Lim
  INTEGER :: U,Num,i,j,k,l

  !Se calculan los límites máximos de la región ocupada por la molécula
  PMin(:)=HUGE(PMin)
  PMax(:)=-HUGE(PMax)
  DO i=1,SIZE(Mol,1)
    PMin(:)=MIN(PMin(:),Mol(i)%pos(:)-FactorMalla*RadiosVdW(Mol(i)%z))
    PMax(:)=MAX(PMax(:),Mol(i)%pos(:)+FactorMalla*RadiosVdW(Mol(i)%z))
  END DO

  !Se escoge la distancia de la red cúbica base, según el tipo de malla
  !Todas las mallas están basadas en una red cúbica subyacente
  SELECT CASE (TipoMalla)
   !Malla cúbica simple
   CASE (1)
    Dist=DistMalla
   !Malla cúbica compacta (centrada en las caras)
   CASE (2)
    Dist=DistMalla/SQRT(2.0D0)
   !Malla cúbica centrada en el cuerpo
   CASE (3)
    Dist=DistMalla/SQRT(3.0D0)
  END SELECT

  !Se calculan los índices mínimos y máximos
  DO i=1,3
    Lim(i,1)=FLOOR(PMin(i)/Dist)
    Lim(i,2)=CEILING(PMax(i)/Dist)
  END DO

  U=NuevaUnidad()
  OPEN(U,STATUS='SCRATCH',ACTION='READWRITE',FORM='UNFORMATTED')

  !Para cada índice, se calcula el punto correspondiente y se guarda si está
  !dentro de la región ocupada por la molécula
  Num=0
  DO i=Lim(1,1),Lim(1,2)
    Pun(1)=i*Dist
    DO j=Lim(2,1),Lim(2,2)
      !Se salta el punto si no corresponde al tipo de malla
      SELECT CASE (TipoMalla)
       CASE (3)
        IF (MODULO(j,2) /= MODULO(i,2)) CYCLE
       CASE DEFAULT
      END SELECT
      Pun(2)=j*Dist
      DO k=Lim(3,1),Lim(3,2)
        !Se salta el punto si no corresponde al tipo de malla
        SELECT CASE (TipoMalla)
         !En una red cúbica compacta sólo son puntos de la red los nodos de la
         !red cúbica base que tienen índices tales que i+j+k=0 (mod 2)
         CASE (2)
          IF (MODULO(i+j+k,2) /= 0) CYCLE
         !En una red cúbica centrada en el cuerpo sólo son puntos de la red los
         !nodos de la red cúbica base tales que i=j=k (mod 2)
         CASE (3)
          IF (MODULO(k,2) /= MODULO(i,2)) CYCLE
         CASE DEFAULT
        END SELECT
        Pun(3)=k*Dist
        DO l=1,SIZE(Mol,1)
          IF (Distancia(Pun(:),Mol(l)%pos(:)) <= &
              FactorMalla*RadiosVdW(Mol(l)%z)) THEN
            WRITE(U) Pun(:)
            Num=Num+1
            EXIT
          END IF
        END DO
      END DO
    END DO
  END DO

  REWIND(U)

  !Finalmente se dimensiona la matriz de la malla y se leen los puntos
  IF (ALLOCATED(MallaMol)) DEALLOCATE(MallaMol)
  ALLOCATE(MallaMol(Num,3))
  DO i=1,Num
    READ(U) MallaMol(i,:)
  END DO

  CLOSE(U)

END SUBROUTINE

END MODULE Malla
