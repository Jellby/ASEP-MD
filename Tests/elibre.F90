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

PROGRAM prueba
  USE Parametros
  USE Entrada
  USE Sistema
  USE GenericoMM
  USE Moldy
  USE DatosQMMM
  USE Configuraciones
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:,:), ALLOCATABLE :: Sol
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ELibre
  DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: IntAt
  LOGICAL, DIMENSION(:,:), ALLOCATABLE :: QInterTot
  CHARACTER(LEN=LLL) :: Linea
  DOUBLE PRECISION :: KT,Beta,Temp=0.0D0
  INTEGER :: U,i,Error

  CALL LeerEntrada(5)

  Extension='.iter.0'
  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM)//TRIM(Extension),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaGenerico(U)
    CLOSE(U)
   CASE (1) !Moldy
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),ACTION='READ',STATUS='OLD')
    CALL LeerControlMoldy(U)
    CLOSE(U)
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(MoldyInput)//TRIM(Extension),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaMoldy(U)
    CLOSE(U)
    Temp=MoldyTemp
  END SELECT

  REWIND(5)
  READ(5,'(A)',IOSTAT=Error) Linea
  DO WHILE (Error == 0)
    IF (INDEX(Linea,'#Temperatura') /= 0) READ(Linea(INDEX(Linea,'=')+1:),*) Temp
    READ(5,'(A)',IOSTAT=Error) Linea
  END DO

  IF (Temp == 0.0D0) STOP "Temperatura = 0"
  KT=KBoltzmann*Temp
  Beta=1.0D0/KT

  ALLOCATE(Sol(0:MaxIter,SIZE(Soluto)))
  ALLOCATE(IntAt(0:MaxIter,SIZE(InterAtom,1),SIZE(InterAtom,2),SIZE(InterAtom,3)))
  ALLOCATE(QInterTot(SIZE(QInter,1),SIZE(QInter,2)))

  ALLOCATE(ELibre(0:MaxIter+1,2))

  QInterTot(:,:)=.FALSE.
  DO i=0,MaxIter
    WRITE(Extension,*) i
    Extension='.iter.'//ADJUSTL(Extension)
    SELECT CASE (ProgramaMM)
     CASE (0) !Generico
      U=NuevaUnidad()
      OPEN(U,FILE=TRIM(EntradaMM)//TRIM(Extension),ACTION='READ',STATUS='OLD')
      CALL LeerSistemaGenerico(U)
      CLOSE(U)
     CASE (1) !Moldy
      U=NuevaUnidad()
      OPEN(U,FILE=TRIM(MoldyInput)//TRIM(Extension),ACTION='READ',STATUS='OLD')
      CALL LeerSistemaMoldy(U)
      CLOSE(U)
    END SELECT
    Sol(i,:)=Soluto(:)
    IntAt(i,:,:,:)=InterAtom(:,:,:)
    QInterTot(:,:)=QInterTot(:,:) .OR. QInter(:,:)
  END DO
  QInter(:,:)=QInterTot(:,:)

  DO i=0,MaxIter
    WRITE(Extension,*) i
    Extension='.iter.'//ADJUSTL(Extension)
    WRITE(6,*)
    WRITE(6,'(A)') TRIM(EntradaMM)//TRIM(Extension)

    Soluto(:)=Sol(i,:)
    CALL LeerConfs()
    CALL CalcularInteracciones(i)
  END DO

  WRITE(6,*)
  WRITE(6,'(A)') 'FEP'
  DO i=1,MaxIter
    WRITE(6,100) ELibre(i,:)/KcalmolAtomica
  END DO
  WRITE(6,*)
  WRITE(6,'(A)') 'Total'
  WRITE(6,100) SUM(ELibre(1:MaxIter,:),DIM=1)/KcalmolAtomica

  DEALLOCATE(Sol,IntAt,QInterTot,ELibre)

100 FORMAT(2(1X,F15.8))

CONTAINS

SUBROUTINE LeerConfs()
  IMPLICIT NONE
  CHARACTER (LEN=LLL) :: Linea,Dump
  INTEGER :: NumMol,NumCuat

  SELECT CASE (ProgramaMM)
   CASE (0) !Generico
    CALL LeerConfigsGenerico(TRIM(TrayectoriaMM)//TRIM(Extension))

   CASE (1) !Moldy
    Dump=TRIM(MoldyDump)//TRIM(Extension)//'-'

    NumMol=0
    NumCuat=0
    NumMol=NumMol+MoleculasSoluto
    IF (SIZE(Soluto,1) > 1) NumCuat=NumCuat+MoleculasSoluto
    NumMol=NumMol+MoleculasDisolvente
    IF (SIZE(Disolvente,1) > 1) NumCuat=NumCuat+MoleculasDisolvente
    NumMol=NumMol+MoleculasDisolvente2
    IF (SIZE(Disolvente2,1) > 1) NumCuat=NumCuat+MoleculasDisolvente2

    WRITE(Linea,100) &
      TRIM(DumpextMoldy),NumMol,NumCuat,1,TRIM(Dump),'c.tmp'
    CALL SYSTEM(TRIM(Linea))
    WRITE(Linea,100) &
      TRIM(DumpextMoldy),NumMol,NumCuat,2,TRIM(Dump),'q.tmp'
    CALL SYSTEM(TRIM(Linea))
    WRITE(Linea,100) &
      TRIM(DumpextMoldy),NumMol,NumCuat,3,TRIM(Dump),'l.tmp'
    CALL SYSTEM(TRIM(Linea))

    CALL LeerConfigsMoldy(1)
  END SELECT

100 FORMAT(A,' -R ',I5,' -Q ',I5,' -c ',I1,' ',A,'* > ',A)

END SUBROUTINE LeerConfs

SUBROUTINE CalcularInteracciones(Iter)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Iter
  TYPE(Atomo), DIMENSION(3,SIZE(Soluto)) :: Mol
  DOUBLE PRECISION, DIMENSION(3,2) :: Prom
  DOUBLE PRECISION, DIMENSION(2) :: ELib, EDip
  DOUBLE PRECISION :: Elec,VdW,ERef
  INTEGER :: Num,Error

  Mol(2,:)=Sol(Iter,:)
  IF (Iter > 0) THEN
    Mol(1,:)=Sol(Iter-1,:)
    CALL SuperponerMoleculas(Mol(2,:),Mol(1,:))
  END IF
  IF (Iter < MaxIter) THEN
    Mol(3,:)=Sol(Iter+1,:)
    CALL SuperponerMoleculas(Mol(2,:),Mol(3,:))
  END IF

  Num=0
  Prom(:,:)=0.0D0
  ELib(:)=0.0D0
  CALL AbrirUConf()
  DO
    CALL LeerConfig(Error)
    IF (Error /= 0) EXIT
    InterAtom(:,:,:)=IntAt(Iter,:,:,:)
    CALL Interacciones(Elec,VdW,EDip,Mol(2,:))
    ERef=Elec+VdW
    Prom(2,1)=Prom(2,1)+Elec
    Prom(2,2)=Prom(2,2)+VdW
    IF (Iter > 0) THEN
      InterAtom(:,:,:)=IntAt(Iter-1,:,:,:)
      CALL Interacciones(Elec,VdW,EDip,Mol(1,:))
      Prom(1,1)=Prom(1,1)+Elec
      Prom(1,2)=Prom(1,2)+VdW
      ELib(1)=ELib(1)+EXP(Beta*(ERef-Elec-VdW))
    END IF
    IF (Iter < MaxIter) THEN
      InterAtom(:,:,:)=IntAt(Iter+1,:,:,:)
      CALL Interacciones(Elec,VdW,EDip,Mol(3,:))
      Prom(3,1)=Prom(3,1)+Elec
      Prom(3,2)=Prom(3,2)+VdW
      ELib(2)=ELib(2)+EXP(Beta*(ERef-Elec-VdW))
    END IF
    Num=Num+1
  END DO
  CALL CerrarUConf()

  IF (Num > 0) THEN
    Prom(:,:)=Prom(:,:)/DBLE(Num)
    ELibre(Iter,1)=KT*LOG(ELib(1)/DBLE(Num))
    ELibre(Iter+1,2)=-KT*LOG(ELib(2)/DBLE(Num))
  END IF
  WRITE(6,100) 'EElec:',Prom(:,1)/KcalmolAtomica
  WRITE(6,100) 'EVdW:',Prom(:,2)/KcalmolAtomica

100 FORMAT(A,T7,3(1X,F15.8))

END SUBROUTINE CalcularInteracciones

END PROGRAM prueba
