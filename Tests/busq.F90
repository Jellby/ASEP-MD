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

MODULE Functions

DOUBLE PRECISION, DIMENSION(2) :: Punto,Vec

CONTAINS

DOUBLE PRECISION FUNCTION MuellerBrown(Factor)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Factor
  DOUBLE PRECISION :: X,Y,F
  INTEGER :: i
  DOUBLE PRECISION, DIMENSION(4), PARAMETER :: &
    A = (/ -200.0D0, -100.0D0, -170.0D0, 15.0D0 /), &
    b = (/ -1.0D0, -1.0D0, -6.5D0, 0.7D0 /), &
    c = (/ 0.0D0, 0.0D0, 11.0D0, 0.6D0 /), &
    d = (/ -10.0D0, -10.0D0, -6.5D0, 0.7D0 /), &
    x0 = (/ 1.0D0, 0.0D0, -0.5D0, -1.0D0 /), &
    y0 = (/ 0.0D0, 0.5D0, 1.5D0, 1.0D0 /)

  F=0.0D0
  DO i=1,4
    X=Punto(1)+Factor*Vec(1)-x0(i)
    Y=Punto(2)+Factor*Vec(2)-y0(i)
    F=F + A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y )
  END DO

  MuellerBrown=F

END FUNCTION MuellerBrown

END MODULE Functions

PROGRAM prueba
  USE Optimizacion
  USE Functions
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(2) :: Punto1,Punto2
  DOUBLE PRECISION :: Minimo
      
  READ(5,*) Punto1(:),Punto2(:)
  Punto(:)=Punto1(:)
  Vec(:)=Punto2(:)-Punto1(:)

  CALL BuscarMinimo(MuellerBrown,0.0D0,1.0D0,Minimo,1.0D-6)

  WRITE(6,100) Punto1(:),MuellerBrown(0.0D0)
  WRITE(6,100) Punto2(:),MuellerBrown(1.0D0)
  WRITE(6,100) Punto1(1)+Minimo*Vec(1),Punto1(2)+Minimo*Vec(2),MuellerBrown(Minimo)

100 FORMAT(5F16.8)

END PROGRAM prueba
