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
  USE Optimizacion
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Hessiana
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Gradiente,Paso,Punto
  DOUBLE PRECISION :: MaxPaso
  INTEGER :: N,Tipo,Obj,Newt,Iter,i,Error
      
  READ(5,*) N,Tipo,Obj,Newt,Iter
  ALLOCATE(Hessiana(N,N),Gradiente(N),Paso(N),Punto(N))
  READ(5,*) MaxPaso,Punto(:)

  WRITE(6,100) Punto(:)
  DO i=1,Iter
    CALL CalcularGradHess(Tipo,Punto,Gradiente,Hessiana)
    CALL IncrementoRFO(Hessiana,Gradiente,Obj,Newt,Paso,Error)
    Paso(:)=Paso(:)*MIN(1.0D0,MaxPaso/SQRT(DOT_PRODUCT(Paso,Paso)))
    Punto(:)=Punto(:)+Paso(:)
    WRITE(6,100) Punto(:)
  END DO

  WRITE(6,*)
  WRITE(6,100) Gradiente(:)

  DEALLOCATE(Hessiana,Gradiente,Paso)

100 FORMAT(5F16.8)

CONTAINS

SUBROUTINE CalcularGradHess(Tipo,Pun,Grad,Hess)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Tipo
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Pun
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: Hess
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Grad
  DOUBLE PRECISION :: X,Y
  DOUBLE PRECISION, DIMENSION(4), PARAMETER :: &
    A = (/ -200.0D0, -100.0D0, -170.0D0, 15.0D0 /), &
    b = (/ -1.0D0, -1.0D0, -6.5D0, 0.7D0 /), &
    c = (/ 0.0D0, 0.0D0, 11.0D0, 0.6D0 /), &
    d = (/ -10.0D0, -10.0D0, -6.5D0, 0.7D0 /), &
    x0 = (/ 1.0D0, 0.0D0, -0.5D0, -1.0D0 /), &
    y0 = (/ 0.0D0, 0.5D0, 1.5D0, 1.0D0 /)

  INTEGER :: N,i

  N=SIZE(Pun)

  SELECT CASE (Tipo)
   CASE (1) ! x^2+y^2
    Grad(:)=2*Pun(:)
    Hess(:,:)=0.0D0
    DO i=1,N
      Hess(i,i)=2.0D0
    END DO
   CASE (2) ! x^3+y^3
    Grad(:)=3*Pun(:)*Pun(:)
    Hess(:,:)=0.0D0
    DO i=1,N
      Hess(i,i)=6*Pun(i)
    END DO
   CASE (3) ! x^2-y^2
    Grad(1)=2*Pun(1)
    Grad(2)=-2*Pun(2)
    Hess(:,:)=0.0D0
    Hess(1,1)=2.0D0
    Hess(2,2)=-2.0D0
   CASE (4) ! x^3+y^3-6xy
    Grad(:)=3*Pun(:)*Pun(:)
    Grad(1)=Grad(1)-6*Pun(2)
    Grad(2)=Grad(2)-6*Pun(1)
    Hess(1,1)=6*Pun(1)
    Hess(2,2)=6*Pun(2)
    Hess(1,2)=-6
    Hess(2,1)=-6
   CASE (5) ! Mueller-Brown
    Grad(:)=0.0D0
    Hess(:,:)=0.0D0; Hess(1,1)=1.0D0; Hess(2,2)=1.0D0
    DO i=1,4
      X=Pun(1)-x0(i)
      Y=Pun(2)-y0(i)
      Grad(1)=Grad(1) + A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                               (2.0D0*b(i)*X + c(i)*Y)
      Grad(2)=Grad(2) + A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                               (2.0D0*d(i)*Y + c(i)*X)
      Hess(1,1)=Hess(1,1) + A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                                   2.0D0*b(i) + &
                            A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                                   (2.0D0*b(i)*X + c(i)*Y)**2
      Hess(2,2)=Hess(2,2) + A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                                   2.0D0*d(i) + &
                            A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                                   (2.0D0*d(i)*Y + c(i)*X)**2
      Hess(1,2)=Hess(1,2) + A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                                   c(i) + &
                            A(i) * EXP( b(i)*X*X + c(i)*X*Y + d(i)*Y*Y ) * &
                                   (2.0D0*b(i)*X + c(i)*Y) * &
                                   (2.0D0*d(i)*Y + c(i)*X)
    END DO
    Hess(2,1)=Hess(1,2)
  END SELECT

END SUBROUTINE CalcularGradHess

END PROGRAM prueba
