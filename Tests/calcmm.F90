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
  USE Moldy
  USE Cavidad
  USE Configuraciones
  USE Sistema
  USE Ejecutar
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  DOUBLE PRECISION :: Elec,VdW,Eelec,EvdW
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Esf
  INTEGER :: U,i,Num,Error

  CALL LeerEntrada(5)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaMM))
  CALL LeerControlMoldy(U)
  CLOSE(U)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(MoldyInput))
  CALL LeerSistemaMoldy(U)
  CLOSE(U)

  CALL OrientarMolecula(Soluto(:),0)
  CALL OrientarMolecula(Disolvente(:),0)
  CALL OrientarMolecula(Disolvente2(:),0)

  IF (TipoCavidad == 0) THEN
    ALLOCATE(Esf(1,4))
    Esf(1,:)=(/0.0D0,0.0D0,0.0D0,RadioCavidad/)
   ELSE
    ALLOCATE(Esf(SIZE(Soluto,1),4))
    DO i=1,SIZE(Esf,1)
      Esf(i,1:3)=Soluto(i)%pos(:)
      IF (TipoCavidad == 1) THEN 
        Esf(i,4)=RadiosVdW(Soluto(i)%z)
       ELSE
        Esf(:,4)=1.0D0
      END IF
    END DO
    Esf(:,4)=RadioCavidad*Esf(:,4)
  END IF

  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4),RadioDisolvente)

  Extension='.mod'

  CALL EjecutarMM

  U=NuevaUnidad()
  OPEN(U,FILE='configs.xyz')
  Eelec=0.0D0
  EvdW=0.0D0
  Num=0
  CALL AbrirUConf()
  DO
    CALL LeerConfig(Error)
    IF (Error /= 0) EXIT
    CALL EscribirXYZ(U)
    CALL Interacciones(Elec,VdW)
    Eelec=Eelec+Elec
    EvdW=EvdW+VdW
    Num=Num+1
    WRITE(6,101) Elec/KcalmolAtomica,VdW/KcalmolAtomica
  END DO
  CLOSE(U)

  CALL CerrarUConf()

  WRITE(6,100) 'PROMEDIO:'
  WRITE(6,101) Eelec/(Num*KcalmolAtomica),EvdW/(Num*KcalmolAtomica)

100 FORMAT (A)
101 FORMAT ('Eelec: ',F12.6,'; Evdw: ',F12.6)

END PROGRAM prueba
