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

MODULE Molcas

#ifndef LLL
#define LLL 256
#endif

USE DatosQM

CONTAINS
!LeerSalidaMolcas(Sal)
!ModificarMolcas(Ent,Sal,Der)

!-------------------------------------------------------------------------------
! Esta subrutina lee la salida del programa Molcas
!-------------------------------------------------------------------------------
! Sal:           Fichero de salida de Molcas
! Coord:         Matriz con las coordenadas para el ajuste de cargas
! Potencial:     Matriz temporal para el potencial generado por el soluto
! PotQ:          Vector temporal para el potencial generado por el disolvente
! Linea:         La línea que se va leyendo
! Num:           Número de átomos
! CI:            Número de la raíz que se calcula
! UCar:          Unidad para leer las cargas atómicas
! Error:         Variable para controlar los errores
! Aux:           Variable auxiliar para lectura
! i,j:           Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerSalidaMolcas(Sal)
  USE Parametros
  USE Unidades
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Sal
 
  DOUBLE PRECISION, DIMENSION(SIZE(MolQM,1),3) :: Coord
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Potencial
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: PotQ
  DOUBLE PRECISION :: Aux
  CHARACTER(LEN=LLL) :: Linea
  INTEGER :: Num,CI,UCar,Error,i,j

  MultipQM=0
  Num=SIZE(MolQM,1)

  !Se inician las matrices del gradiente y hessiana
  IF (ALLOCATED(GradQM)) DEALLOCATE(GradQM)
  IF (ALLOCATED(HessQM)) DEALLOCATE(HessQM)
  ALLOCATE(GradQM(3*Num),HessQM(3*Num,3*Num))
  GradQM(:)=0.0D0
  HessQM(:,:)=0.0D0

  !Lee los datos del fichero de salida
  Lectura: DO
    READ(Sal,'(A)',IOSTAT=Error) Linea
    IF (Error /= 0 ) EXIT
    Linea=ADJUSTL(Linea)

    !Carga del cálculo SCF
    IF (Linea(1:16) == 'Molecular charge') THEN
      READ(Linea(17:),*) Aux
      CargaQM=INT(Aux)
    !Carga total a partir de la población de Mulliken
    ELSE IF (Linea(1:23) == 'Total            charge') THEN
      READ(Linea(25:),*) Aux
      CargaQM=INT(Aux)

    !Multiplicidad del cálculo SCF
    ELSE IF (Linea(1:33) == 'Fermi aufbau procedure completed!') THEN
      READ(Sal,'(A)') Linea
      IF (ADJUSTL(Linea(1:13)) == 'nOcc(alpha)=') THEN
        READ(Linea(14:),*) i
        READ(Sal,'(A)') Linea
        READ(Linea(14:),*) j
       ELSE
        READ(Linea(6:),*) i
        j=i
      END IF
      MultipQM=ABS(i-j)+1
    !Multiplicidad del cálculo RASSCF
    ELSE IF (Linea(1:19) == 'Spin quantum number') THEN
      READ(Linea(26:),*) Aux
      MultipQM=INT(2*Aux)+1

    !Raíz que se usa en un cálculo CI/RASSCF
    ELSE IF ((Linea(1:12) == 'CI root used') .OR. &
             (Linea(1:13) == 'CI roots used') .OR. &
             (Linea(1:11) == 'Root passed')) THEN
      READ(Linea(40:),*) CI

    !Energía del cálculo SCF
    ELSE IF (Linea(1:16) == 'Total SCF energy') THEN
      READ(Linea(17:),*) EnergiaQM
    !Energía del cálculo DFT
    ELSE IF (Linea(1:19) == 'Total KS-DFT energy') THEN
      READ(Linea(20:),*) EnergiaQM
    !Energías del cálculo RASSCF
    ELSE IF (Linea(1:23) == 'Final state energy(ies)') THEN
      READ(Sal,*)
      READ(Sal,*)
      DO
        READ(Sal,'(A)') Linea
        IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
        IF (ADJUSTL(Linea(1:17)) /= 'root number') EXIT
        READ(Linea(18:21),*) i
        IF (i == CI) READ(Linea(25:44),*) EnergiaQM
      END DO
    !Energía del cálculo MBPT2 o CASPT2
    ELSE IF (Linea(1:13) == 'Total energy') THEN
      IF (INDEX(Linea,'=') == 0) THEN
        READ(Linea(14:),*) EnergiaQM
       ELSE
        READ(Linea(39:59),*) EnergiaQM
      END IF
    !Energía del cálculo CCSD
    ELSE IF (Linea(1:21) == 'Correlation energy  :') THEN
      READ(Linea(22:),*) Aux
      READ(Sal,'(A)') Linea
      Linea=ADJUSTL(Linea)
      READ(Linea(22:),*) EnergiaQM
      EnergiaQM=EnergiaQM+Aux
    !Energía del cálculo CCSD(T)
    ELSE IF (Linea(1:9) == 'CCSD + T3') THEN
      READ(Linea(18:),*) EnergiaQM
    !Energía del cálculo MRCI
    ELSE IF (Linea(1:17) == 'CORRECTED ENERGY:') THEN
      READ(Linea(18:),*) EnergiaQM

    !Propiedades (cargas, dipolo, potencial) de otras raíces RASSCF
    ELSE IF (Linea(30:45) == 'for root number:') THEN
      READ(Linea(46:),*) i
      IF (i == CI) CYCLE
      DO
        READ(Sal,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
        IF (INDEX(Linea,'Stop Module:') /= 0) EXIT
        IF (Linea(36:51) == 'for root number:') THEN
          BACKSPACE(Sal)
          EXIT
        END IF
      END DO

    !Saltar la salida de MCLR
    ELSE IF (INDEX(Linea,'MOLCAS executing module MCLR') /= 0) THEN
      DO
        READ(Sal,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
        IF (INDEX(Linea,'Stop Module:') /= 0) EXIT
      END DO

    !Momento dipolar
    ELSE IF (Linea(1:13) == 'Dipole Moment') THEN
      READ(Sal,*)
      READ(Sal,100) DipoloQM(:)
      DipoloQM(:)=DipoloQM(:)*DebyeAtomica
    !Momento dipolar del cálculo MRCI
    ELSE IF (Linea(1:18) == 'PROPERTY :MLTPL  1') THEN
      READ(Linea(32:),*) i
      DO
        READ(Sal,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
        IF (ADJUSTL(Linea(1:18)) == 'TOTAL:') EXIT
      END DO
      READ(Linea(19:),*) DipoloQM(i)
      !??? DipoloQM(i)=DipoloQM(i)*DebyeAtomica

    !Cargas de Mulliken
    ELSE IF (Linea(1:16) == 'Mulliken charges') THEN
      IF (TipoCargas /= 0) CYCLE !Mulliken
      DO i=1,Num,12
        DO
          READ(Sal,'(A)',IOSTAT=Error) Linea
          IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
          IF (ADJUSTL(Linea(1:11)) == 'N-E') EXIT
          IF (ADJUSTL(Linea(1:11)) == 'Charge') EXIT
          IF (INDEX(Linea,'spin') /= 0) CYCLE Lectura
        END DO
        READ(Linea(12:),*) MolQM(i:MIN(i+11,Num))%q
      END DO

    !Gradiente analítico
    ELSE IF (Linea(16:34) == 'Molecular gradients') THEN
      DO
        READ(Sal,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
        IF (ADJUSTL(Linea) == 'Irreducible representation: a') EXIT
      END DO
      READ(Sal,*)
      DO i=1,SIZE(GradQM,1)
        READ(Sal,'(A)') Linea
        READ(Linea(27:),*) GradQM(i)
      END DO
    !Gradiente numérico
    ELSE IF (Linea(1:42) == 'MOLCAS executing module NUMERICAL_GRADIENT') THEN
      DO
        READ(Sal,'(A)',IOSTAT=Error) Linea
        IF (Error /= 0 ) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
        Linea=ADJUSTL(Linea)
        IF (Linea(1:19) == 'Numerical gradients') EXIT
      END DO
      IF (INDEX(Linea,'not implemented') /= 0) CYCLE
      READ(Sal,*)
      DO i=1,SIZE(GradQM,1)
        READ(Sal,*)
        READ(Sal,*) GradQM(i)
      END DO

    !Hessiana analítica
    ELSE IF (Linea(1:14) == '*BEGIN HESSIAN') THEN
      READ(Sal,*)
      DO i=1,SIZE(HessQM,1)
        READ(Sal,*)
        READ(Sal,*) HessQM(:,i)
      END DO

    !Potencial electrostático
    ELSE IF (Linea(1:19) == 'Electric Potential:') THEN
      DO i=1,SIZE(DisolvQM,1)
        DO
          READ(Sal,'(A)',IOSTAT=Error) Linea
          IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
          IF (ADJUSTL(Linea(1:20)) == 'Total') EXIT
        END DO
        READ(Linea(21:),*) DisolvQM(i,5)
      END DO
      DO i=1,SIZE(PotQM,1)
        DO
          READ(Sal,'(A)',IOSTAT=Error) Linea
          IF (Error /= 0) CALL Mensaje('LeerSalidaMolcas',27,.TRUE.)
          IF (ADJUSTL(Linea(1:20)) == 'Total') EXIT
        END DO
        READ(Linea(21:),*) PotQM(i,4)
      END DO

    ELSE IF (Linea(1:20) == 'Non-zero return code') THEN
      CALL Mensaje('LeerSalidaMolcas',28,.TRUE.)
    END IF
  END DO Lectura

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
      CALL Mensaje('LeerSalidaMolcas',16,.FALSE.)
    END IF
  END IF

  !O se leen de un fichero externo
  IF (TipoCargas == 3) THEN
    UCar=NuevaUnidad()
    OPEN(UCar,FILE=TRIM(FicheroCargas),STATUS='OLD',ACTION='READ')
    READ(UCar,*) MolQM(:)%q
    CLOSE(UCar)
  END IF

  !Se igualan las cargas equivalentes para el Moldy
  IF (ProgramaMM == 1) CALL IgualarCargas(MolQM(:)%id,MolQM(:)%q)

100 FORMAT(3(17X,F10.4))

END SUBROUTINE LeerSalidaMolcas

!-------------------------------------------------------------------------------
! Esta subrutina modifica la entrada de Molcas para incluir la geometría del
! soluto, las cargas que representan al disolvente, etc.
!-------------------------------------------------------------------------------
! Ent:      Fichero de entrada (plantilla)
! Sal:      Fichero de salida (entrada de Molcas)
! Der:      Orden de la derivada de la energía que se quiere calcular
! Car:      Fichero de cargas externas (opcional)
! Linea:    La línea que se va leyendo
! FichTemp: Nombre del fichero donde se escriben las cargas externas
! Error:    Variable para controlar los errores
! i:        Contador
!-------------------------------------------------------------------------------
SUBROUTINE ModificarMolcas(Ent,Sal,Der,Car)
  USE Unidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Ent,Sal,Der
  INTEGER, INTENT(IN), OPTIONAL :: Car

  CHARACTER(LEN=LLL) :: Linea,FichTemp
  INTEGER :: Error,i

  !Copia línea a línea, sustituyendo:
  ! ###(G)### -> geometría de la molécula
  ! ###(F)### -> cálculo de derivadas
  ! ###(C)### -> cargas externas
  ! ###(P)### -> puntos donde se calcula el potencial
  DO
    READ(Ent,'(A)',IOSTAT=Error) Linea
    IF (Error /= 0) EXIT
    SELECT CASE (TRIM(ADJUSTL(Linea)))

     CASE ('###(G)###')
      WRITE(Sal,'(A)') '  NoSymm'
      WRITE(Sal,*) SIZE(MolQM,1)
      WRITE(Sal,*)
      DO i=1,SIZE(MolQM,1)
        WRITE(Sal,100) Simbolo(MolQM(i)%z),MolQM(i)%pos(:)/AngstromAtomica
      END DO

     CASE ('###(F)###')
      SELECT CASE (Der)
       CASE (1)
        WRITE(Sal,'(A)') '&ALASKA &END'
        WRITE(Sal,'(A)') 'End of input'
       CASE (2)
        WRITE(Sal,'(A)') '&ALASKA &END'
        WRITE(Sal,'(A)') 'End of input'
        WRITE(Sal,'(A)')
        WRITE(Sal,'(A)') '&MCKINLEY &END'
        WRITE(Sal,'(A)') 'End of input'
       CASE DEFAULT
      END SELECT

     CASE ('###(C)###')
      IF (SIZE(DisolvQM,1) > 0) THEN
        WRITE(Sal,'(A)') '  XField'
        IF (PRESENT(Car)) THEN
          WRITE(Car,101) SIZE(DisolvQM,1), 1
          DO i=1,SIZE(DisolvQM,1)
            WRITE(Car,102) DisolvQM(i,1:4)
          END DO
          INQUIRE(Car,NAME=FichTemp)
          WRITE(Sal,'(A)') '    @$MOLCAS_SUBMIT_PWD/'//TRIM(FichTemp)
         ELSE
          WRITE(Sal,101) SIZE(DisolvQM,1), 1
          DO i=1,SIZE(DisolvQM,1)
            WRITE(Sal,102) DisolvQM(i,1:4)
          END DO
        END IF
      END IF

     CASE ('###(P)###')
      IF (SIZE(DisolvQM,1)+SIZE(PotQM,1) > 0) THEN
        WRITE(Sal,'(A)') '  EPot'
        WRITE(Sal,101) SIZE(DisolvQM,1)+SIZE(PotQM,1)
       DO i=1,SIZE(DisolvQM,1)
         WRITE(Sal,102) DisolvQM(i,1:3)
       END DO
       DO i=1,SIZE(PotQM,1)
         WRITE(Sal,102) PotQM(i,1:3)
       END DO
      END IF

     CASE DEFAULT
      WRITE(Sal,'(A)') TRIM(Linea)
    END SELECT
  END DO

100 FORMAT (2X,A2,3(1X,F19.12))
101 FORMAT (3X,5(1X,I6))
102 FORMAT (3X,4(1X,F19.12),'  0.0 0.0 0.0')

END SUBROUTINE ModificarMolcas

END MODULE Molcas
