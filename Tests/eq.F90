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

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ecuaciones,Eq
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Soluciones,Sol
  INTEGER :: N,Rand,i,j
      
  READ(5,*) N
  Rand=0
  IF (N < 0) THEN
    N=-N
    Rand=1
  END IF
  ALLOCATE(Ecuaciones(N,N),Soluciones(N),Eq(N,N),Sol(N))
  IF (Rand == 0) THEN
    READ(5,*) (Ecuaciones(i,:),i=1,N)
    READ(5,*) Soluciones(:)
   ELSE
    CALL IniciarAleatorio()
    DO i=1,N
      CALL RANDOM_NUMBER(Soluciones(i))
      DO j=i,N
        CALL RANDOM_NUMBER(Ecuaciones(i,j))
      END DO
    END DO
  END IF

  Eq(:,:)=Ecuaciones(:,:)
  Sol(:)=Soluciones(:)
  CALL ResolverLU(Ecuaciones,Soluciones)

  WRITE(6,100) Soluciones(:)
  WRITE(6,*)
  WRITE(6,101) ABS(MATMUL(Eq(:,:),Soluciones(:))-Sol(:))
  WRITE(6,*)
  WRITE(6,*) MAXVAL(ABS(MATMUL(Eq(:,:),Soluciones(:))-Sol(:)))

  DEALLOCATE(Ecuaciones,Soluciones,Eq,Sol)

100 FORMAT(10F8.3)
101 FORMAT(10ES8.1)

END PROGRAM prueba
