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
  USE Parametros
  USE Entrada
  USE Sistema
  USE DatosQM
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mulliken,ESP,Ajust1,Ajust2,Ajust3,Ajust4
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Puntos
  DOUBLE PRECISION :: InteraccionQM
  CHARACTER(LEN=LLL) :: Linea,Aux
  INTEGER :: U,Num,i,Error

  CALL LeerEntrada(5)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaQM))
  DO
    READ (U,100,IOSTAT=Error) Linea
    IF (Error /= 0) EXIT
    SELECT CASE (TRIM(Linea))
     CASE ('Geometry')
      READ(U,*) Num
      ALLOCATE(Soluto(Num))
      DO i=1,Num
        READ (U,*) Soluto(i)%nom,Aux,Soluto(i)%z,Soluto(i)%pos(:)
        Soluto(i)%pos(:)=Soluto(i)%pos(:)*AngstromAtomica
      END DO
     CASE ('External charges')
      READ(U,*) Num
      ALLOCATE(DisolvQM(Num,5))
      DO i=1,Num
        READ(U,*) DisolvQM(i,1:4)
      END DO
      DisolvQM(:,1:3)=DisolvQM(:,1:3)*AngstromAtomica
    END SELECT
  END DO
  CLOSE(U)

  Num=SIZE(Soluto)
  ALLOCATE(Puntos(Num,3),Mulliken(Num),ESP(Num),Ajust1(Num),Ajust2(Num),Ajust3(Num),Ajust4(Num))

  DO i=1,Num
    Puntos(i,:)=Soluto(i)%pos(:)
  END DO
  CALL PotencialCargas(DisolvQM(:,1:4),Puntos(:,:),Soluto(:)%q)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(SalidaQM))
  DO 
    READ (U,100,IOSTAT=Error) Linea
    IF (Error /= 0) EXIT
    SELECT CASE (TRIM(Linea))
     CASE ('Charge')
      READ(U,*) CargaQM
     CASE ('Dipole moment')
      READ(U,*) DipoloQM(:)
     CASE ('Electrostatic potential')
      READ(U,*) DisolvQM(:,5)
     CASE ('Mulliken charges')
      READ(U,*) Mulliken(:)
     CASE ('ESP charges')
      READ(U,*) ESP(:)
    END SELECT
  END DO

  InteraccionQM=SUM(DisolvQM(:,4)*DisolvQM(:,5))

  CALL Intercambiar(DisolvQM(:,4),DisolvQM(:,5))

  CALL AjusteCargas(Puntos(:,:),DisolvQM(:,1:4),DBLE(CargaQM),Ajust1(:))

  CALL AjusteCargas(Puntos(:,:),DisolvQM(:,1:4),DBLE(CargaQM),Ajust2(:), &
                    Dip=DipoloQM(:))

  CALL AjusteCargas(Puntos(:,:),DisolvQM(:,1:4),DBLE(CargaQM),Ajust3(:), &
                    PotQ=Soluto(:)%q,Energia=InteraccionQM)

  CALL AjusteCargas(Puntos(:,:),DisolvQM(:,1:4),DBLE(CargaQM),Ajust4(:), &
                    Dip=DipoloQM(:),PotQ=Soluto(:)%q,Energia=InteraccionQM)

  CALL Intercambiar(DisolvQM(:,4),DisolvQM(:,5))

  WRITE(6,101) 'Carga QM: ',DBLE(CargaQM)
  WRITE(6,101) 'Dipolo QM: ',DipoloQM(:)/DebyeAtomica
  WRITE(6,101) 'Energia QM: ',InteraccionQM/KcalmolAtomica
  WRITE(6,*)
  WRITE(6,101) 'Carga Mulliken: ',CargaTotal(Mulliken)
  WRITE(6,101) 'Dipolo Mulliken: ',Dipolo(Mulliken)/DebyeAtomica
  WRITE(6,101) 'Energia Mulliken: ',Interaccion(Mulliken)/KcalmolAtomica
  WRITE(6,101) 'Errores: ',ErrAjuste(Mulliken)
  WRITE(6,*)
  WRITE(6,101) 'Carga ESP: ',CargaTotal(ESP)
  WRITE(6,101) 'Dipolo ESP: ',Dipolo(ESP)/DebyeAtomica
  WRITE(6,101) 'Energia ESP: ',Interaccion(ESP)/KcalmolAtomica
  WRITE(6,101) 'Errores: ',ErrAjuste(ESP)
  WRITE(6,*)
  WRITE(6,101) 'Carga Fit: ',CargaTotal(Ajust1)
  WRITE(6,101) 'Dipolo Fit: ',Dipolo(Ajust1)/DebyeAtomica
  WRITE(6,101) 'Energia Fit: ',Interaccion(Ajust1)/KcalmolAtomica
  WRITE(6,101) 'Errores: ',ErrAjuste(Ajust1)
  WRITE(6,*)
  WRITE(6,101) 'Carga Fit(D): ',CargaTotal(Ajust2)
  WRITE(6,101) 'Dipolo Fit(D): ',Dipolo(Ajust2)/DebyeAtomica
  WRITE(6,101) 'Energia Fit(D): ',Interaccion(Ajust2)/KcalmolAtomica
  WRITE(6,101) 'Errores: ',ErrAjuste(Ajust2)
  WRITE(6,*)
  WRITE(6,101) 'Carga Fit(E): ',CargaTotal(Ajust3)
  WRITE(6,101) 'Dipolo Fit(E): ',Dipolo(Ajust3)/DebyeAtomica
  WRITE(6,101) 'Energia Fit(E): ',Interaccion(Ajust3)/KcalmolAtomica
  WRITE(6,101) 'Errores: ',ErrAjuste(Ajust3)
  WRITE(6,*)
  WRITE(6,101) 'Carga Fit(D,E): ',CargaTotal(Ajust4)
  WRITE(6,101) 'Dipolo Fit(D,E): ',Dipolo(Ajust4)/DebyeAtomica
  WRITE(6,101) 'Energia Fit(D,E): ',Interaccion(Ajust4)/KcalmolAtomica
  WRITE(6,101) 'Errores: ',ErrAjuste(Ajust4)

  WRITE(6,*)
  WRITE(6,103) 'Pot','Mulliken','ESP','Fit','Fit(D)','Fit(E)','Fit(D,E)'
  DO i=1,SIZE(Soluto)
    WRITE(6,102) Simbolo(Soluto(i)%z),Soluto(i)%q,Mulliken(i),ESP(i),Ajust1(i),Ajust2(i),Ajust3(i),Ajust4(i)
  END DO

  DEALLOCATE(Soluto,DisolvQM,Puntos,Mulliken,ESP,Ajust1,Ajust2,Ajust3,Ajust4)

100 FORMAT (A)
101 FORMAT (A20,3F10.6)
102 FORMAT (A5,7F10.6)
103 FORMAT (5X,7A10)

CONTAINS

FUNCTION CargaTotal(Cargas)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(Soluto)) :: Cargas
  DOUBLE PRECISION :: CargaTotal

  CargaTotal=SUM(Cargas(:))

END FUNCTION CargaTotal

FUNCTION Interaccion(Cargas)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(Soluto)) :: Cargas
  DOUBLE PRECISION :: Interaccion

  Interaccion=SUM(Cargas(:)*Soluto(:)%q)

END FUNCTION Interaccion

FUNCTION Dipolo(Cargas)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(Soluto)) :: Cargas
  DOUBLE PRECISION, DIMENSION(3) :: Dipolo

  Dipolo(1)=SUM(Cargas(:)*Soluto(:)%pos(1))
  Dipolo(2)=SUM(Cargas(:)*Soluto(:)%pos(2))
  Dipolo(3)=SUM(Cargas(:)*Soluto(:)%pos(3))

END FUNCTION Dipolo

FUNCTION ErrAjuste(Cargas)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(SIZE(Soluto)) :: Cargas
  DOUBLE PRECISION, DIMENSION(3) :: ErrAjuste
  DOUBLE PRECISION, DIMENSION(SIZE(Soluto),4) :: Puntos
  DOUBLE PRECISION, DIMENSION(SIZE(DisolvQM,1)) :: Potencial
  INTEGER :: i

  DO i=1,SIZE(Puntos,1)
    Puntos(i,1:3)=Soluto(i)%pos(:)
    Puntos(i,4)=Cargas(i)
  END DO

  CALL PotencialCargas(Puntos(:,:),DisolvQM(:,1:3),Potencial(:))
  Potencial(:)=Potencial(:)-DisolvQM(:,5)
  i=SIZE(Potencial,1)

  !1. RErr=SUMA[|(V-V')/V|)]/n           n -> numero de puntos en el potencial
  !2. Rms=SQRT[SUMA[(V-V')**2]/n]        V -> potencial real
  !3. RRms=SQRT[SUMA[((V-V')/V)**2]/n]   V'-> potencial calculado
  ErrAjuste(1)=SUM(ABS(Potencial(:)/DisolvQM(:,5)))/i
  ErrAjuste(2)=SQRT(SUM(Potencial(:)**2)/i)
  ErrAjuste(3)=SQRT(SUM(Potencial(:)**2/(DisolvQM(:,5)**2))/i)

END FUNCTION ErrAjuste

END PROGRAM prueba
