!##############################################################################
!# Copyright 2011,2012 Ignacio Fdez. Galván, M. Luz Sánchez,                  #
!#                     Aurora Muñoz Losa, M. Elena Martín, Manuel A. Aguilar  #
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

MODULE GenericoQM

#ifndef LLL
#define LLL 256
#endif

USE DatosQM

CONTAINS
!LeerSalidaGenericoQM(Sal)
!EntradaGenericoQM(Sal,Der,Car)

!-------------------------------------------------------------------------------
! Subrutina que lee una salida genérica de cálculo cuántico
!-------------------------------------------------------------------------------
! Sal:       Fichero de salida
! Coord:     Matriz con las coordenadas para el ajuste de cargas
! Potencial: Matriz temporal para el potencial generado por el soluto
! PotQ:      Vector temporal para el potencial generado por el disolvente
! Linea:     La línea que se va leyendo
! Num:       Número de átomos
! Fin:       Variable para controlar el final del fichero
! Aux:       Variable auxiliar para lectura
! UCar:      Unidad para leer las cargas atómicas
! i,j:       Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerSalidaGenericoQM(Sal)
  USE Parametros
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Sal

  DOUBLE PRECISION, DIMENSION(SIZE(MolQM,1),3) :: Coord
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Potencial
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PotQ
  CHARACTER(LEN=LLL) :: Linea
  INTEGER :: Num,Fin,UCar,i,j

  !Lee los datos del fichero de salida
  DO
    CALL LeerSiguienteLinea(Sal,Linea,Fin)
    IF (Fin /= 0) EXIT
    Linea=ADJUSTL(Linea)
    CALL PasarMinusculas(Linea)

    SELECT CASE (TRIM(Linea))

     !Número de átomos (y dimensiona las matrices)
     CASE ('number of atoms')
      READ(Sal,*) Num
      IF (ALLOCATED(GradQM)) DEALLOCATE(GradQM)
      IF (ALLOCATED(Grad2QM)) DEALLOCATE(Grad2QM)
      IF (ALLOCATED(HessQM)) DEALLOCATE(HessQM)
      ALLOCATE(GradQM(3*Num),Grad2QM(3*Num),HessQM(3*Num,3*Num))
      GradQM(:)=0.0D0
      HessQM(:,:)=0.0D0

     !Carga de la molécula
     CASE ('charge')
      READ(Sal,*) CargaQM

     !Multiplicidad de espín
     CASE ('multiplicity')
      READ(Sal,*) MultipQM

     !Energía
     CASE ('energy')
      READ(Sal,*) EnergiaQM

     !Energía del estado inferior
     CASE ('energy (lower state)')
      READ(Sal,*) Energia2QM

     !Momento dipolar
     CASE ('dipole moment')
      READ(Sal,*) DipoloQM(:)

     !Cargas de Mulliken
     CASE ('mulliken charges')
      IF (TipoCargas == 0) READ(Sal,*) MolQM(:)%q

     !Cargas ESP
     CASE ('esp charges')
      IF (TipoCargas == 1) READ(Sal,*) MolQM(:)%q

     !Gradiente
     CASE ('cartesian gradient')
      READ(Sal,*) GradQM(:)

     !Gradiente del estado inferior
     CASE ('cartesian gradient (lower state)')
      READ(Sal,*) Grad2QM(:)

     !Hessiana
     CASE ('cartesian hessian')
      READ(Sal,*) ((HessQM(i,j),j=1,i),i=1,SIZE(HessQM,1))
      DO i=2,SIZE(HessQM,1)
        HessQM(1:i-1,i)=HessQM(i,1:i-1)
      END DO

     !Potencial electrostático
     CASE ('electrostatic potential')
      IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
        READ(Sal,*) DisolvQM(:,5),PotQM(:,4)
      END IF
 
     CASE DEFAULT
    END SELECT
  END DO

  !Si las cargas son ajustadas al potencial, se calculan
  IF (TipoCargas == 2) THEN
    IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
      !Se crea una matriz temporal con el potencial conjunto de DisolvQM y PotQM
      ALLOCATE(Potencial(SIZE(DisolvQM,1)+SIZE(PotQM,1),4))
      Potencial(1:SIZE(DisolvQM,1),1:3)=DisolvQM(:,1:3)
      Potencial(1:SIZE(DisolvQM,1),4)=DisolvQM(:,5)
      Potencial(SIZE(DisolvQM,1)+1:,:)=PotQM(:,:)
      Coord(:,1)=MolQM(:)%pos(1)
      Coord(:,2)=MolQM(:)%pos(2)
      Coord(:,3)=MolQM(:)%pos(3)
      IF (SIZE(DisolvQM,1) > 0) THEN
        ALLOCATE(PotQ(SIZE(MolQM)))
        CALL PotencialCargas(DisolvQM(:,1:4),Coord,PotQ)
        CALL AjusteCargas(Coord,Potencial,DBLE(CargaQM),MolQM(:)%q, &
                          PotQ=PotQ,Energia=SUM(DisolvQM(:,4)*DisolvQM(:,5)))
        DEALLOCATE(PotQ)
       ELSE
        CALL AjusteCargas(Coord,Potencial,DBLE(CargaQM),MolQM(:)%q)
      END IF
      DEALLOCATE(Potencial)
     ELSE
      !Si no hay potencial externo para ajustar, se dejan las cargas originales
      CALL Mensaje('LeerSalidaGenericoQM',16,.FALSE.)
    END IF
  END IF

  !Si las cargas son externas, se leen del fichero
  IF (TipoCargas == 3) THEN
    UCar=NuevaUnidad()
    OPEN(UCar,FILE=TRIM(FicheroCargas),STATUS='OLD',ACTION='READ')
    READ(UCar,*) MolQM(:)%q
    CLOSE(UCar)
  END IF

  !Se igualan las cargas equivalentes
  CALL IgualarCargas(MolQM(:)%id,MolQM(:)%q)

END SUBROUTINE LeerSalidaGenericoQM

!-------------------------------------------------------------------------------
! Subrutina que genera una entrada genérica de cálculo cuántico
!-------------------------------------------------------------------------------
! Sal:      Fichero de salida (entrada genérica)
! Der:      Orden de la derivada de la energía que se quiere calcular
! Car:      Fichero de cargas externas (opcional)
! FichTemp: Nombre del fichero donde se escriben las cargas externas
! i:        Contador
!-------------------------------------------------------------------------------
SUBROUTINE EntradaGenericoQM(Sal,Der,Car)
  USE Unidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Sal,Der
  INTEGER, INTENT(IN), OPTIONAL :: Car

  CHARACTER(LEN=LLL) :: FichTemp
  INTEGER :: i

  !Escribe la geometría
  WRITE(Sal,'(A)') 'Geometry'
  WRITE(Sal,101) SIZE(MolQM,1)
  DO i=1,SIZE(MolQM,1)
    WRITE(Sal,100) MolQM(i)%nom,Simbolo(MolQM(i)%z),MolQM(i)%z, &
                   MolQM(i)%pos(:)/AngstromAtomica
  END DO
  WRITE(Sal,*)

  !Escribe el orden de la derivada que hay que calcular
  WRITE(Sal,'(A)') 'Derivative'
  WRITE(Sal,101) Der
  WRITE(Sal,*)

  !Escribe las cargas externas
  IF (PRESENT(Car)) THEN
    IF (SIZE(DisolvQM,1) > 0) THEN
      DO i=1,SIZE(DisolvQM,1)
        WRITE(Car,102) DisolvQM(i,1:3)/AngstromAtomica,DisolvQM(i,4)
      END DO
    END IF
    INQUIRE(Car,NAME=FichTemp)
    WRITE(Sal,'(A)') 'External charge file'
    WRITE(Sal,103) TRIM(FichTemp)
    WRITE(Sal,*)
   ELSE
    WRITE(Sal,'(A)') 'External charges'
    WRITE(Sal,101) SIZE(DisolvQM,1)
    IF (SIZE(DisolvQM,1) > 0) THEN
      DO i=1,SIZE(DisolvQM,1)
        WRITE(Sal,102) DisolvQM(i,1:3)/AngstromAtomica,DisolvQM(i,4)
      END DO
    END IF
    WRITE(Sal,*)
  END IF

  !Escribe los puntos donde se calcula el potencial
  WRITE(Sal,'(A)') 'Potential points'
  WRITE(Sal,101) SIZE(DisolvQM,1)+SIZE(PotQM,1)
  IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
    DO i=1,SIZE(DisolvQM,1)
      WRITE(Sal,102) DisolvQM(i,1:3)/AngstromAtomica
    END DO
    DO i=1,SIZE(PotQM,1)
      WRITE(Sal,102) PotQM(i,1:3)/AngstromAtomica
    END DO
  END IF
  WRITE(Sal,*)

100 FORMAT (1X,A16,1X,A2,1X,I3,3(1X,F19.12))
101 FORMAT (4(1X,I5))
102 FORMAT (4(1X,F19.12))
103 FORMAT (1X,A)

END SUBROUTINE EntradaGenericoQM

END MODULE GenericoQM
