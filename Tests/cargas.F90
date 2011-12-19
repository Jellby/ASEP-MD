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
  USE Unidades
  USE UtilidadesFis
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Cargas,Puntos
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Potencial
  DOUBLE PRECISION, DIMENSION(3) :: ErrAjuste,Dipolo
  DOUBLE PRECISION :: CargaTotal,Rest
  INTEGER :: NumCargas,NumPuntos,i

  READ(5,*) NumCargas,CargaTotal,Rest

  ALLOCATE(Cargas(NumCargas,4))

  DO i=1,NumCargas
    READ(5,*) Cargas(i,1:3)
  END DO
  Cargas(:,1:3)=Cargas(:,1:3)*AngstromAtomica
  Cargas(:,4)=0.0D0

  READ(5,*) Dipolo(:)

  READ(5,*) NumPuntos

  ALLOCATE(Puntos(NumPuntos,4),Potencial(NumPuntos))

  DO i=1,NumPuntos
    READ(5,*) Puntos(i,1:4)
  END DO
  Puntos(:,1:3)=Puntos(:,1:3)*AngstromAtomica

  IF (ANY(Dipolo(:) /= 0.0D0)) THEN
    CALL AjusteCargas(Cargas(:,1:3),Puntos(:,:),CargaTotal,Cargas(:,4),ErrAjuste,Rest,Dipolo(:))
   ELSE
    CALL AjusteCargas(Cargas(:,1:3),Puntos(:,:),CargaTotal,Cargas(:,4),ErrAjuste,Rest)
  END IF

  CALL PotencialCargas(Cargas(:,:),Puntos(:,1:3),Potencial(:))

  WRITE(6,100) ErrAjuste
  WRITE(6,102) MAXVAL(ABS(Potencial(:)-Puntos(:,4)))
  WRITE(6,*)
  DO i=1,NumCargas
    WRITE(6,101) Cargas(i,1:3)/AngstromAtomica,Cargas(i,4)
  END DO

  DO i=1,3
    Dipolo(i)=SUM(Cargas(:,i)*Cargas(:,4))
  END DO
  WRITE(6,*)
  WRITE(6,103) Dipolo(:)/DebyeAtomica

  DEALLOCATE(Cargas,Puntos,Potencial)

100 FORMAT ('RErr: ',F10.6,'; Rms: ',F10.6,'; RRms: ',F10.6)
101 FORMAT (4(F12.8,1X))
102 FORMAT ('Error maximo en el potencial: ',F10.6)
103 FORMAT ('Momento dipolar: ',F10.6,', ',F10.6,', ',F10.6)

END PROGRAM prueba
