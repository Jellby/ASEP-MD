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
  USE Utilidades
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Matriz,Matriz2,Vectores
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Valores
  INTEGER :: N,Rand,i,j

  READ(5,*) N
  Rand=0
  IF (N < 0) Rand=1
  N=ABS(N)
  ALLOCATE(Matriz(N,N),Matriz2(N,N),Vectores(N,N),Valores(N))
  IF (Rand == 0) THEN
    READ(5,*) (Matriz(i,:),i=1,N)
   ELSE
    CALL IniciarAleatorio()
    DO i=1,N
      DO j=i,N
        CALL RANDOM_NUMBER(Matriz(i,j))
        Matriz(j,i)=Matriz(i,j)
      END DO
    END DO
  END IF

  CALL Diagonalizar(Matriz,Vectores,Valores,1)

  Matriz2=0.0D0
  DO i=1,N
    Matriz2(i,i)=Valores(i)
  END DO
  Matriz2=MATMUL(MATMUL(Vectores,Matriz2),TRANSPOSE(Vectores))

  WRITE(6,100) Valores(:)
  WRITE(6,*)
  DO i=1,N
    WRITE(6,100) Vectores(i,:)
  END DO
  WRITE(6,*)
  DO i=1,N
    WRITE(6,101) ABS(Matriz2(i,:)-Matriz(i,:))
  END DO
  WRITE(6,*)
  WRITE(6,*) MAXVAL(ABS(Matriz2(i,:)-Matriz(i,:)))

  DEALLOCATE(Matriz,Matriz2,Vectores,Valores)

100 FORMAT(10F8.3)
101 FORMAT(10ES8.1)

END PROGRAM prueba
