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

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Matriz,Matriz2,Inversa, &
                                                   Vectores,Inv2,Aux
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Valores
  INTEGER :: M,N,Rand,i,j
      
  READ(5,*) M,N
  Rand=0
  IF ((M < 0) .OR. (N < 0)) Rand=1
  M=ABS(M)
  N=ABS(N)
  ALLOCATE(Matriz(M,N),Matriz2(M,N),Inversa(N,M),Vectores(N,N),Valores(N), &
           Inv2(N,M),Aux(N,N))
  IF (Rand == 0) THEN
    READ(5,*) (Matriz(i,:),i=1,M)
   ELSE
    CALL IniciarAleatorio()
    DO i=1,M
      DO j=1,N
        CALL RANDOM_NUMBER(Matriz(i,j))
      END DO
    END DO
  END IF

  Matriz2=Matriz
  CALL Pseudoinversa(Matriz,Valores,Vectores,Inversa)
  Matriz=Matriz2

  WRITE(6,100) Valores(:)
  WRITE(6,*)
  DO i=1,N
    WRITE(6,100) Inversa(i,:)
  END DO
  WRITE(6,*)
  Matriz2=MATMUL(MATMUL(Matriz,Inversa),Matriz)
  WRITE(6,*) 'A - (A * A^+ * A)'
  DO i=1,M
    WRITE(6,101) ABS(Matriz2(i,:)-Matriz(i,:))
  END DO
  WRITE(6,*)
  WRITE(6,*) MAXVAL(ABS(Matriz2(i,:)-Matriz(i,:)))
  WRITE(6,*)
  WRITE(6,*) 'A^+ - (A^+ * A * A^+)'
  Inv2=MATMUL(MATMUL(Inversa,Matriz),Inversa)
  DO i=1,N
    WRITE(6,101) ABS(Inv2(i,:)-Inversa(i,:))
  END DO
  WRITE(6,*)
  WRITE(6,*) MAXVAL(ABS(Inv2(i,:)-Inversa(i,:)))

  DEALLOCATE(Matriz,Matriz2,Vectores,Valores,Inversa,Inv2,Aux)

100 FORMAT(10F8.3)
101 FORMAT(10ES8.1)

END PROGRAM prueba
