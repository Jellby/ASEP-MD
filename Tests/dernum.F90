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
  USE Ejecutar
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION :: Delta
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Energias
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Energias2,Gradientes
  INTEGER :: Error,Der,Num,U,i,ii,j,jj,k,kk
  CHARACTER(LEN=LLL) :: Linea

  Der=0
  Delta=1.0D-3

  CALL LeerEntrada(5)

  REWIND(5)
  READ(5,'(A)',IOSTAT=Error) Linea
  DO WHILE (Error == 0)
    IF (INDEX(Linea,'#Delta') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Delta
    IF (INDEX(Linea,'#Derivada') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Der
    READ(5,'(A)',IOSTAT=Error) Linea
  END DO

  Der=MIN(MAX(Der,0),1)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaMM))
  CALL LeerControlMoldy(U)
  CLOSE(U)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(MoldyInput))
  CALL LeerSistemaMoldy(U)
  CLOSE(U)

  Extension='.mod'

  Num=SIZE(Soluto,1)*3
  ALLOCATE(MolQM(SIZE(Soluto,1)),Energias(Num,2),Energias2(Num,Num,4),&
           Gradientes(Num,Num,2))
  MolQM=Soluto

  DO i=1,Num
    j=(i-1)/3+1
    k=MOD(i-1,3)+1
    MolQM(j)%pos(k)=Soluto(j)%pos(k)-Delta
    CALL EjecutarQM(Der)
    Energias(i,1)=EnergiaQM
    Gradientes(i,:,1)=GradQM
    IF ((CalcHessiana > 0) .AND. (Der == 0)) THEN
      DO ii=1,i-1
        jj=(ii-1)/3+1
        kk=MOD(ii-1,3)+1
        MolQM(jj)%pos(kk)=Soluto(jj)%pos(kk)-Delta
        CALL EjecutarQM(Der)
        Energias2(i,ii,1)=EnergiaQM
        MolQM(jj)%pos(kk)=Soluto(jj)%pos(kk)+Delta
        CALL EjecutarQM(Der)
        Energias2(i,ii,2)=EnergiaQM
        MolQM(jj)%pos(kk)=Soluto(jj)%pos(kk)
      END DO
    END IF
    MolQM(j)%pos(k)=Soluto(j)%pos(k)+Delta
    CALL EjecutarQM(Der)
    Energias(i,2)=EnergiaQM
    Gradientes(i,:,2)=GradQM
    IF ((CalcHessiana > 0) .AND. (Der == 0)) THEN
      DO ii=1,i-1
        jj=(ii-1)/3+1
        kk=MOD(ii-1,3)+1
        MolQM(jj)%pos(kk)=Soluto(jj)%pos(kk)-Delta
        CALL EjecutarQM(Der)
        Energias2(i,ii,3)=EnergiaQM
        MolQM(jj)%pos(kk)=Soluto(jj)%pos(kk)+Delta
        CALL EjecutarQM(Der)
        Energias2(i,ii,4)=EnergiaQM
        MolQM(jj)%pos(kk)=Soluto(jj)%pos(kk)
      END DO
    END IF
    MolQM(j)%pos(k)=Soluto(j)%pos(k)
  END DO

  CALL EjecutarQM(Der)

  IF (Der == 0) THEN
    DO i=1,Num
      GradQM(i)=(Energias(i,2)-Energias(i,1))/(2.0D0*Delta)
      IF (CalcHessiana == 0) CYCLE
      HessQM(i,i)=(Energias(i,1)+Energias(i,2)-2.0D0*EnergiaQM)/(Delta*Delta)
      DO j=1,i-1
        HessQM(i,j)= &
          (Energias2(i,j,1)-Energias2(i,j,2)-Energias2(i,j,3)+Energias2(i,j,4))/ &
          (4.0D0*Delta*Delta)
        HessQM(j,i)=HessQM(i,j)
      END DO
    END DO
   ELSE
    DO i=1,Num
      HessQM(i,:)=(Gradientes(i,:,2)-Gradientes(i,:,1))/(2.0D0*Delta)
    END DO
  END IF

  WRITE(6,*) 'Energia:'
  WRITE(6,100) EnergiaQM
  WRITE(6,*)
  WRITE(6,*) 'Gradiente:'
  DO i=1,SIZE(Soluto,1)
    WRITE(6,100) GradQM((i-1)*3+1:i*3)
  END DO
  IF (CalcHessiana > 0) THEN
    WRITE(6,*)
    WRITE(6,*) 'Hessiana:'
    DO i=1,Num
      WRITE(6,100) HessQM(i,1:i)
    END DO
  END IF

100 FORMAT(5ES16.8)

END PROGRAM prueba
