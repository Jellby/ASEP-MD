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
  USE Coordenadas
  USE GenericoMM
  USE Moldy
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), ALLOCATABLE :: Inicial,Final,Previo
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CoordFin,CoordDif,Coord,Paso, &
                                                 Pesos
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: InterAtomIni,InterAtomFin
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: QInterIni,QInterFin
  LOGICAL :: ModCoord,Sup
  CHARACTER(LEN=LLL) :: Linea,Fich1,Fich2
  DOUBLE PRECISION :: Lambda
  INTEGER :: U,i,j,Error,UEnt,USal

  CALL LeerEntrada(5)
  Sup=.FALSE.

  REWIND(5)
  READ(5,'(A)',IOSTAT=Error) Linea
  DO WHILE (Error == 0)
    IF (INDEX(Linea,'#Inicial') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Fich1
    IF (INDEX(Linea,'#Final') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Fich2
    IF (INDEX(Linea,'#Superponer') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Sup
    READ(5,'(A)',IOSTAT=Error) Linea
  END DO

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(Fich2),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaGenerico(U)
    CLOSE(U)

   CASE (1) !Moldy
    UEnt=NuevaUnidad()
    OPEN(UEnt,FILE=TRIM(EntradaMM),STATUS='OLD',ACTION='READ')
    CALL LeerControlMoldy(UEnt)

    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(Fich2),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaMoldy(U)
    CLOSE(U)
  END SELECT

  ALLOCATE(Final(SIZE(Soluto,1)))
  Final(:)=Soluto(:)
  ALLOCATE(InterAtomFin(SIZE(InterAtom,1),SIZE(InterAtom,2),SIZE(InterAtom,3)))
  ALLOCATE(QInterFin(SIZE(QInter,1),SIZE(QInter,2)))
  InterAtomFin(:,:,:)=InterAtom(:,:,:)
  QInterFin(:,:)=QInter(:,:)

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(Fich1),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaGenerico(U)
    CLOSE(U)

   CASE (1) !Moldy
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(Fich1),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaMoldy(U)
    CLOSE(U)
  END SELECT

  ALLOCATE(Inicial(SIZE(Soluto,1)))
  Inicial(:)=Soluto(:)
  ALLOCATE(InterAtomIni(SIZE(InterAtom,1),SIZE(InterAtom,2),SIZE(InterAtom,3)))
  ALLOCATE(QInterIni(SIZE(QInter,1),SIZE(QInter,2)))
  InterAtomIni(:,:,:)=InterAtom(:,:,:)
  QInterIni(:,:)=QInter(:,:)

  ALLOCATE(Pesos(SIZE(Soluto,1)),Previo(SIZE(Soluto,1)))
  Pesos(:)=1.0D0
  REWIND(5)
  READ(5,'(A)',IOSTAT=Error) Linea
  DO WHILE (Error == 0)
    IF (INDEX(Linea,'#Pesos') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Pesos(:)
    READ(5,'(A)',IOSTAT=Error) Linea
  END DO

  IF (Sup) CALL SuperponerMoleculas(Inicial,Final,Pesos)

  TipoCoord=TipoCoordenadas
  CALL DefinirCoordenadas(Inicial)
  INQUIRE(FILE=TRIM(AltCoordenadas),EXIST=ModCoord)
  IF (ModCoord) THEN
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(AltCoordenadas))
    CALL ModificarCoordenadas(U)
    CLOSE(U)
  END IF

  i=SIZE(DefCoord,1)
  ALLOCATE(CoordFin(i),CoordDif(i),Coord(i),Paso(3*SIZE(Soluto,1)))
  CALL ConvertirCoordenadas(Cartesianas(Final),1)
  CoordFin(:)=Geometria(:)

  QInter(:,:)=(QInterIni(:,:) .OR. QInterFin(:,:))

  Soluto(:)=Inicial(:)
  DO i=0,MaxIter
    Lambda=i/DBLE(MaxIter)
    Previo(:)=Soluto(:)
    CALL ConvertirCoordenadas(Cartesianas(Soluto),1)
    IF (i > 0) THEN
      CoordDif(:)=CoordFin(:)-Geometria(:)
      CALL CorregirDiedros(CoordDif(:))
      CoordDif(:)=CoordDif(:)/DBLE(MaxIter-i+1)
      CALL ConvertirIncremento(CoordDif,Cartesianas(Soluto),Paso)
      CALL CambiarCartesianas(Soluto,Cartesianas(Soluto)+Paso)
    END IF
    Soluto(:)%q=Inicial(:)%q*(1.0D0-Lambda)+Final(:)%q*Lambda
    Soluto(:)%m=Inicial(:)%m*(1.0D0-Lambda)+Final(:)%m*Lambda
    InterAtom(:,:,:)=InterAtomIni(:,:,:)*(1.0D0-Lambda)+ &
                     InterAtomFin(:,:,:)*Lambda

    CALL OrientarMolecula(Soluto,0)

    WRITE(Extension,*) i
    Extension='.iter.'//ADJUSTL(Extension)

    SELECT CASE (ProgramaMM)
     CASE (0) !Genérico
      USal=NuevaUnidad()
      OPEN(USal,FILE=TRIM(EntradaMM)//TRIM(Extension),STATUS='REPLACE',ACTION='WRITE')
      CALL EntradaGenericoMM(USal)
      CLOSE(USal)
 
     CASE (1) !Moldy
      REWIND(UEnt)
      USal=NuevaUnidad()
      OPEN(USal,FILE=TRIM(EntradaMM)//TRIM(Extension),STATUS='REPLACE',ACTION='WRITE')
      CALL ModificarMoldy(UEnt,USal)
      CLOSE(USal)
    END SELECT

    WRITE(6,*) SIZE(Soluto,1)
    WRITE(6,*)
    DO j=1,SIZE(Soluto,1)
      WRITE(6,*) Simbolo(Soluto(j)%z),Soluto(j)%pos/AngstromAtomica
    END DO
  END DO

  IF (ProgramaMM == 1) CLOSE(UEnt)

  DEALLOCATE(Inicial,Final,Previo,CoordFin,CoordDif,Coord,Paso,Pesos, &
             InterAtomIni,InterAtomFin,QInterIni,QInterFin)

CONTAINS

SUBROUTINE CorregirDiedros(Vec)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Vec
  INTEGER :: i

  DO i=1,SIZE(DefCoord,1)
    IF (DefCoord(i,4) == 0) CYCLE
    IF (ABS(Vec(i)) > Pi) THEN
      Vec(i)=MODULO(Vec(i),-SIGN(2.0D0*Pi,Vec(i)))
    END IF
  END DO

END SUBROUTINE CorregirDiedros

END PROGRAM prueba
