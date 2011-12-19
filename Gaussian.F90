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

MODULE Gaussian

#ifndef LLL
#define LLL 256
#endif

USE DatosQM

!-------------------------------------------------------------------------------
! GausFichPot*: Ficheros para el cálculo del potencial electrostático
!-------------------------------------------------------------------------------
INTEGER, PARAMETER :: GausFichPot1=63,GausFichPot2=64

CONTAINS
!LeerSalidaGaussian(Sal,Fchk)
!ModificarGaussian(Ent,Sal,Der)

!-------------------------------------------------------------------------------
! Esta subrutina lee la salida del programa Gaussian
!-------------------------------------------------------------------------------
! Sal:           Fichero de salida de Gaussian
! Fchk:          Fichero fchk de Gaussian
! Coord:         Matriz con las coordenadas para el ajuste de cargas
! Potencial:     Matriz temporal para el potencial generado por el soluto
! PotQ:          Vector temporal para el potencial generado por el disolvente
! EnergiaCargas: Energía de interacción de las cargas externas
! Linea:         La línea que se va leyendo
! FichTemp:      Nombre del fichero de donde se lee el potencial
! Num:           Número de átomos
! UPot,UCar:     Unidades para leer el potencial y las cargas atómicas
! Error:         Variable para controlar los errores
! Aux:           Variable auxiliar para lectura
! i,j:           Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerSalidaGaussian(Sal,Fchk)
  USE Parametros
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Sal,Fchk

  DOUBLE PRECISION, DIMENSION(SIZE(MolQM,1),3) :: Coord
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Potencial
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PotQ
  DOUBLE PRECISION :: EnergiaCargas,Aux
  CHARACTER(LEN=LLL) :: Linea,FichTemp
  INTEGER :: Num,UPot,UCar,Error,i,j

  !Lee los datos presentes en el fichero fchk
  DO
    READ(Fchk,'(A)',IOSTAT=Error) Linea
    IF (Error /= 0) EXIT
    SELECT CASE (TRIM(ADJUSTL(Linea(1:40))))

     CASE ('Number of atoms')
      READ(Linea(45:),*) Num
      IF (ALLOCATED(GradQM)) DEALLOCATE(GradQM)
      IF (ALLOCATED(HessQM)) DEALLOCATE(HessQM)
      ALLOCATE(GradQM(3*Num),HessQM(3*Num,3*Num))
      GradQM(:)=0.0D0
      HessQM(:,:)=0.0D0

     CASE ('Charge')
      READ(Linea(45:),*) CargaQM

     CASE ('Multiplicity')
      READ(Linea(45:),*) MultipQM

     CASE ('Total Energy')
      READ(Linea(45:),*) EnergiaQM

     CASE ('Dipole Moment')
      READ(Fchk,*) DipoloQM(:)

     CASE ('Cartesian Gradient')
      READ(Fchk,*) GradQM(:)

     CASE ('Cartesian Force Constants')
      READ(Fchk,*) ((HessQM(i,j),j=1,i),i=1,SIZE(HessQM,1))
      DO i=2,SIZE(HessQM,1)
        HessQM(1:i-1,i)=HessQM(i,1:i-1)
      END DO

     CASE DEFAULT
    END SELECT
  END DO

  !Se lee el potencial generado por el soluto
  IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
    WRITE(FichTemp,*) GausFichPot2
    FichTemp='fort.'//ADJUSTL(TRIM(FichTemp))
    UPot=NuevaUnidad()
    OPEN(UPot,FILE=TRIM(FichTemp),STATUS='OLD',ACTION='READ')
    DO i=1,SIZE(DisolvQM,1)
      READ(UPot,*) Aux,Aux,Aux,DisolvQM(i,5)
    END DO
    DO i=1,SIZE(PotQM,1)
      READ(UPot,*) Aux,Aux,Aux,PotQM(i,4)
    END DO
    CLOSE(UPot)
  END IF

  !Se lee la energía de las cargas y se resta
  IF (SIZE(DisolvQM,1) > 0) THEN
    DO
      READ(Sal,'(A)',IOSTAT=Error) Linea
      IF (Error /= 0) CALL Mensaje('LeerSalidaGaussian',9,.TRUE.)
      IF (INDEX(Linea,'Self energy of the charges =') /= 0) THEN
        Linea=Linea(INDEX(Linea,'=')+1:)
        READ(Linea,*) EnergiaCargas
        EXIT
      END IF
    END DO
    EnergiaQM=EnergiaQM-EnergiaCargas
  END IF

  !Se leen las cargas atómicas del fichero de salida
  SELECT CASE (TipoCargas)

   CASE (0) !Mulliken
    DO
      READ(Sal,'(A)',IOSTAT=Error) Linea
      IF (Error /= 0) CALL Mensaje('LeerSalidaGaussian',9,.TRUE.)
      IF ((ADJUSTL(Linea) == 'Total atomic charges:') .OR. &
          (ADJUSTL(Linea) == 'Mulliken atomic charges:')) THEN
        READ(Sal,'(A)') Linea
        EXIT
      END IF
    END DO
    DO i=1,SIZE(MolQM,1)
      READ(Sal,*) j,Linea,MolQM(i)%q
    END DO

   CASE (1) !CHELP y similares
    DO
      READ(Sal,'(A)',IOSTAT=Error) Linea
      IF (Error /= 0) CALL Mensaje('LeerSalidaGaussian',9,.TRUE.)
      IF (INDEX(Linea,'Charges from ESP fit') /= 0) THEN
        READ(Sal,'(A)') Linea
        READ(Sal,'(A)') Linea
        EXIT
      END IF
    END DO
    DO i=1,SIZE(MolQM,1)
      READ(Sal,*) j,Linea,MolQM(i)%q
    END DO

   !O se calculan por ajuste al potencial
   CASE(2) !Potencial
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
      CALL Mensaje('LeerSalidaGaussian',16,.FALSE.)
    END IF

   CASE(3) !Externo
    UCar=NuevaUnidad()
    OPEN(UCar,FILE=TRIM(FicheroCargas),STATUS='OLD',ACTION='READ')
    READ(UCar,*) MolQM(:)%q
    CLOSE(UCar)
  END SELECT

  !Se igualan las cargas equivalentes
  CALL IgualarCargas(MolQM(:)%id,MolQM(:)%q)

END SUBROUTINE LeerSalidaGaussian

!-------------------------------------------------------------------------------
! Esta subrutina modifica la entrada de Gaussian para incluir la geometría del
! soluto, las cargas que representan al disolvente, etc.
!-------------------------------------------------------------------------------
! Ent:      Fichero de entrada (plantilla)
! Sal:      Fichero de salida (entrada de Gaussian)
! Der:      Orden de la derivada de la energía que se quiere calcular
! Car:      Fichero de cargas externas (opcional)
! Linea:    La línea que se va leyendo
! FichTemp: Nombre del fichero donde se escriben los puntos para el potencial
! UPot:     Unidad donde se escriben los puntos para el potencial
! Error:    Variable para controlar los errores
! i:        Contador
!-------------------------------------------------------------------------------
SUBROUTINE ModificarGaussian(Ent,Sal,Der,Car)
  USE Unidades
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Ent,Sal,Der
  INTEGER, INTENT(IN), OPTIONAL :: Car

  CHARACTER(LEN=LLL) :: Linea,FichTemp
  INTEGER :: UPot,Error,i

  !Copia el preámbulo (hasta la primera línea en blanco)
  DO
    READ(Ent,'(A)',IOSTAT=Error) Linea
    IF (Error /= 0) CALL Mensaje('ModificarGaussian',4,.TRUE.)
    IF (TRIM(Linea) == '') EXIT
    WRITE(Sal,'(A)') TRIM(Linea)
  END DO

  !Añade los datos necesarios
  WRITE(Sal,'(A)') 'FCheck=All NoSymm'
  SELECT CASE (Der)
   CASE (1)
    WRITE(Sal,'(A)') 'Force'
   CASE (2)
    WRITE(Sal,'(A)') 'Freq'
   CASE DEFAULT
  END SELECT
  IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) &
     WRITE(Sal,'(A)') 'Prop=(Grid,Potential)'
  IF (SIZE(DisolvQM,1) > 0) WRITE(Sal,'(A)') 'Charge=Angstrom'
  WRITE(Sal,'(A)')

  !Copia línea a linea, sustituyendo:
  ! ###(G)### -> geometría de la molécula
  ! ###(C)### -> cargas externas
  ! ###(P)### -> puntos donde se calcula el potencial
  DO
    READ(Ent,'(A)',IOSTAT=Error) Linea
    IF (Error /= 0) EXIT
    SELECT CASE (TRIM(ADJUSTL(Linea)))

     CASE ('###(G)###')
      DO i=1,SIZE(MolQM,1)
        WRITE(Sal,100) Simbolo(MolQM(i)%z),MolQM(i)%pos(:)/AngstromAtomica
      END DO

     CASE ('###(C)###')
      IF (SIZE(DisolvQM,1) > 0) THEN
        IF (PRESENT(Car)) THEN
          DO i=1,SIZE(DisolvQM,1)
            WRITE(Car,101) DisolvQM(i,1:3)/AngstromAtomica,DisolvQM(i,4)
          END DO
          INQUIRE(Car,NAME=FichTemp)
          WRITE(Sal,'(A)') '@'//TRIM(FichTemp)//'/N'
         ELSE
          DO i=1,SIZE(DisolvQM,1)
            WRITE(Sal,101) DisolvQM(i,1:3)/AngstromAtomica,DisolvQM(i,4)
          END DO
        END IF
       ELSE
        READ(Ent,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0) CALL Mensaje('ModificarGaussian',4,.TRUE.)
        IF (TRIM(Linea) /= '') BACKSPACE(Ent)
      END IF

     CASE ('###(P)###')
      IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
        WRITE(Sal,102) SIZE(DisolvQM,1)+SIZE(PotQM,1),3,GausFichPot1, &
                       GausFichPot2
       ELSE
        READ(Ent,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0) CALL Mensaje('ModificarGaussian',4,.TRUE.)
        IF (TRIM(Linea) /= '') BACKSPACE(Ent)
      END IF

     CASE DEFAULT
      WRITE(Sal,'(A)') TRIM(Linea)
    END SELECT
  END DO

  !Escribe los puntos donde se calcula el potencial en un fichero externo
  IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
    WRITE(FichTemp,*) GausFichPot1
    FichTemp='fort.'//ADJUSTL(TRIM(FichTemp))
    UPot=NuevaUnidad()
    OPEN(UPot,FILE=TRIM(FichTemp),STATUS='REPLACE',ACTION='WRITE')
    DO i=1,SIZE(DisolvQM,1)
      WRITE(UPot,101) DisolvQM(i,1:3)/AngstromAtomica
    END DO
    DO i=1,SIZE(PotQM,1)
      WRITE(UPot,101) PotQM(i,1:3)/AngstromAtomica
    END DO
    CLOSE(UPot)
  END IF

100 FORMAT (A2,3(1X,F19.12))
101 FORMAT (4(F19.12,:,1X))
102 FORMAT (I7,3(1X,I4))

END SUBROUTINE ModificarGaussian

END MODULE Gaussian
