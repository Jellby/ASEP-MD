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
  DOUBLE PRECISION, DIMENSION(3) :: MinP,MaxP
  INTEGER :: Num,U,i,j,k,l

  Num=2
  MaxP=(/15.0D0,15.0D0,15.0D0/)
  MinP=-MaxP
  !ALLOCATE(PotQM(Num**3,4))
  ALLOCATE(DisolvQM(Num**3,5))
  l=1
  DO i=1,Num
    DO j=1,Num
      DO k=1,Num
        !PotQM(l,1)=MinP(1)+(i-1)*(MaxP(1)-MinP(1))/DBLE(Num-1)
        !PotQM(l,2)=MinP(2)+(j-1)*(MaxP(2)-MinP(2))/DBLE(Num-1)
        !PotQM(l,3)=MinP(3)+(k-1)*(MaxP(3)-MinP(3))/DBLE(Num-1)
        DisolvQM(l,1)=MinP(1)+(i-1)*(MaxP(1)-MinP(1))/DBLE(Num-1)
        DisolvQM(l,2)=MinP(2)+(j-1)*(MaxP(2)-MinP(2))/DBLE(Num-1)
        DisolvQM(l,3)=MinP(3)+(k-1)*(MaxP(3)-MinP(3))/DBLE(Num-1)
        DisolvQM(l,4)=0.0D0
        l=l+1
      END DO
    END DO
  END DO

  CALL LeerEntrada(5)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaMM))
  CALL LeerControlMoldy(U)
  CLOSE(U)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(MoldyInput))
  CALL LeerSistemaMoldy(U)
  CLOSE(U)

  Extension='.mod'

  ALLOCATE(MolQM(SIZE(Soluto,1)))
  MolQM=Soluto

  IF (CalcHessiana > 0) THEN
    CALL EjecutarQM(2)
   ELSE
    CALL EjecutarQM(1)
  END IF

  Soluto=MolQM

  DEALLOCATE(PotQM)

  WRITE(6,103) 'Carga','Multiplicidad'
  WRITE(6,101) CargaQM,MultipQM
  WRITE(6,100) 'Energia:',  EnergiaQM,'Eh'
  WRITE(6,100) 'Dipolo:',   Norma(DipoloQM)/DebyeAtomica,'D'
  WRITE(6,100) 'Gradiente:',Norma(GradQM),'Eh/a0'
  WRITE(6,100) '    (rms):',Norma(GradQM)/SQRT(DBLE(SIZE(GradQM,1))),'Eh/a0'
  WRITE(6,100) 'Cargas:'
  DO i=1,SIZE(Soluto,1)
    WRITE(6,100) '',        Soluto(i)%q,'e'
  END DO

  IF (CalcHessiana > 0 ) THEN
    WRITE(6,103) 'Hessiana:'
    WRITE(6,102) ((HessQM(i,j),j=1,i),i=1,SIZE(HessQM,1))
  END IF

100 FORMAT(A,T12,F12.6,1X,A)
101 FORMAT(I4,T8,I4)
102 FORMAT(8F10.6)
103 FORMAT(A,:,2X,A)

END PROGRAM prueba
