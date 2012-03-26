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

MODULE Configuraciones

USE Sistema

!-------------------------------------------------------------------------------
! Coords:     Coordenadas de todos los átomos
! CentroMasa: Coordenadas de todos los centros de masas
! MolSol:     Índices de las moléculas de soluto
! MolDis:     Índices de las moléculas de disolvente
! MolDis2:    Índices de las moléculas de disolvente 2
! DefCelda:   Definición (vectores) de la celda de simulación
! Krf,Crf:    Constantes para el cálculo con "reaction field"
! Config:     Número de configuración
! Centro:     Molécula donde se centra la configuración
! UConf:      Fichero binario donde se escriben las configuraciones
! Fuera:      Vector que define para cada átomo si se considera en la config.
! Fichero:    Variable que indica si el fichero binario está abierto o no
!-------------------------------------------------------------------------------
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Coords,CentroMasa
DOUBLE PRECISION, DIMENSION(3,3) :: DefCelda
DOUBLE PRECISION :: Krf,Crf
INTEGER, DIMENSION(:), ALLOCATABLE :: MolSol,MolDis,MolDis2
INTEGER :: Config,Centro,UConf
LOGICAL, DIMENSION(:), ALLOCATABLE :: Fuera
LOGICAL :: FicheroConf=.FALSE.

CONTAINS
!AbrirUConf
!CerrarUConf
!LeerConfig(Error)
!EscribirPDB(Updb)
!EscribirXYZ(Uxyz,Centros)
!Interacciones(Elec,VdW,MolCen,GrVdW,HsVdW)
!InteraccionPar(At1,At2,Elec,VdW,QGrad,Grad,QHess,Hess)
!Promedios(Der,Mol)
!PotencialMM(Puntos,Potencial)
!SeleccionarMols(U,NumCargas)
!ReducirCargas(U,Num,Cargas)

!-------------------------------------------------------------------------------
! Abrir el fichero binario de configuraciones
!-------------------------------------------------------------------------------
SUBROUTINE AbrirUConf
  USE Utilidades
  IMPLICIT NONE

  IF (.NOT. FicheroConf) THEN
    UConf=NuevaUnidad()
    OPEN(UConf,STATUS='SCRATCH',FORM='UNFORMATTED')
    FicheroConf=.TRUE.
   ELSE
    REWIND(UConf)
  END IF

END SUBROUTINE

!-------------------------------------------------------------------------------
! Cerrar el fichero binario de configuraciones
!-------------------------------------------------------------------------------
SUBROUTINE CerrarUConf
  USE Utilidades
  IMPLICIT NONE

  IF (FicheroConf) CLOSE(UConf)
  FicheroConf=.FALSE.

END SUBROUTINE

!-------------------------------------------------------------------------------
! Lee una configuración del fichero binario
!-------------------------------------------------------------------------------
! Error:   Variable para controlar los errores
! Dist:    Distancia de la molécula al origen
! Corte:   Radio de corte
! Num:     Número de átomos
! NumMol:  Número de moléculas
! i,j,k:   Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerConfig(Error)
  USE Parametros
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: Error

  DOUBLE PRECISION :: Dist,Corte
  INTEGER :: Num,NumMol,i,j,k

  !Se calcula el número de átomos
  Num=0
  Num=Num+MoleculasSoluto*SIZE(Soluto,1)
  Num=Num+MoleculasDisolvente*SIZE(Disolvente,1)
  Num=Num+MoleculasDisolvente2*SIZE(Disolvente2,1)
  NumMol=0
  NumMol=NumMol+MoleculasSoluto
  NumMol=NumMol+MoleculasDisolvente
  NumMol=NumMol+MoleculasDisolvente2
  IF (.NOT. ALLOCATED(MolSol)) ALLOCATE(MolSol(MoleculasSoluto))
  IF (.NOT. ALLOCATED(MolDis)) ALLOCATE(MolDis(MoleculasDisolvente))
  IF (.NOT. ALLOCATED(MolDis2)) ALLOCATE(MolDis2(MoleculasDisolvente2))
  IF (.NOT. ALLOCATED(Coords)) ALLOCATE(Coords(Num,3))
  IF (.NOT. ALLOCATED(CentroMasa)) ALLOCATE(CentroMasa(NumMol,3))
  IF (.NOT. ALLOCATED(Fuera)) ALLOCATE(Fuera(Num))

  !Se leen los datos globales de la configuración
  READ(UConf,IOSTAT=Error) Config
  IF (Error /= 0) RETURN
  READ(UConf) Centro
  READ(UConf) DefCelda(:,:)
  READ(UConf)
  READ(UConf)

  !Se calcula el radio de corte (si se usa Moldy)
  IF ((Centro > 0) .AND. (ProgramaMM == 1)) THEN
    Corte=0.5D0*MIN(DefCelda(1,1),DefCelda(2,2),DefCelda(3,3))
   ELSE
    Corte=HUGE(Corte)
  END IF
  Fuera(:)=.FALSE.

  !Se leen las coordenadas de cada molécula, eliminando las que estén fuera
  !del radio de corte
  Num=0
  NumMol=0
  !Un número negativo significa que la molécula se ha eliminado
  MolSol(:)=-1
  j=SIZE(Soluto,1)
  DO i=1,MoleculasSoluto
    READ(UConf) Dist
    READ(UConf) CentroMasa(NumMol+i,:)
    READ(UConf)
    READ(UConf) Coords(Num+1:Num+j,:)
    IF (Dist > Corte) CYCLE
    MolSol(i)=Num
    !Si el MM no es Moldy, se establece el corte por átomos
    IF (ProgramaMM /= 1) THEN
      DO k=Num+1,Num+j
        IF (DOT_PRODUCT(Coords(k,:),Coords(k,:)) > Corte*Corte) Fuera(k)=.TRUE.
      END DO
    END IF
    Num=Num+j
  END DO
  NumMol=NumMol+MoleculasSoluto
  MolDis(:)=-1
  j=SIZE(Disolvente,1)
  DO i=1,MoleculasDisolvente
    READ(UConf) Dist
    READ(UConf) CentroMasa(NumMol+i,:)
    READ(UConf)
    READ(UConf) Coords(Num+1:Num+j,:)
    IF (Dist > Corte) CYCLE
    MolDis(i)=Num
    IF (ProgramaMM /= 1) THEN
      DO k=Num+1,Num+j
        IF (DOT_PRODUCT(Coords(k,:),Coords(k,:)) > Corte*Corte) Fuera(k)=.TRUE.
      END DO
    END IF
    Num=Num+j
  END DO
  NumMol=NumMol+MoleculasDisolvente
  MolDis2(:)=-1
  j=SIZE(Disolvente2,1)
  DO i=1,MoleculasDisolvente2
    READ(UConf) Dist
    READ(UConf) CentroMasa(NumMol+i,:)
    READ(UConf)
    READ(UConf) Coords(Num+1:Num+j,:)
    IF (Dist > Corte) CYCLE
    MolDis2(i)=Num
    IF (ProgramaMM /= 1) THEN
      DO k=Num+1,Num+j
        IF (DOT_PRODUCT(Coords(k,:),Coords(k,:)) > Corte*Corte) Fuera(k)=.TRUE.
      END DO
    END IF
    Num=Num+j
  END DO
  NumMol=NumMol+MoleculasDisolvente2

END SUBROUTINE LeerConfig

!-------------------------------------------------------------------------------
! Escribe una configuración en formato PDB
!-------------------------------------------------------------------------------
! Updb:    Fichero donde se escribe la configuración
! Centros: Variable que indica si se incluyen los centros de masa en la salida
! QCen:    Igual que Centros
! Linea:   Variable de texto para el título
! Num:     Número de átomo
! NumMol:  Número de molécula (residuo)
! NumMol2: Número de molécula total
! i,j:     Contadores
!-------------------------------------------------------------------------------
SUBROUTINE EscribirPDB(Updb,Centros)
  USE Parametros
  USE Unidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Updb
  LOGICAL, INTENT(IN), OPTIONAL :: Centros

  CHARACTER(LEN=LLL) :: Linea
  LOGICAL :: QCen
  INTEGER :: Num,NumMol,NumMol2,i,j

  QCen=.FALSE.
  IF ((ProgramaMM == 1) .AND. PRESENT(Centros)) QCen=Centros

  !Escribe el título (la configuración y la molécula central)
  WRITE(Linea,100) Config,Centro
  WRITE(Updb,101) '',Linea

  !Escribe los átomos del soluto
  Num=0
  NumMol=0
  NumMol2=0
  DO i=1,MoleculasSoluto
    NumMol=NumMol+1
    IF (MolSol(i) < 0) CYCLE
    IF (Qcen) THEN
      Num=Num+1
      WRITE(Updb,103) Num,'X','','SOL','X',NumMol,'', &
                      CentroMasa(NumMol2+i,:)/AngstromAtomica,1.0,0.0,'X',''
    END IF
    DO j=1,SIZE(Soluto,1)
      IF (Fuera(MolSol(i)+j)) CYCLE
      Num=Num+1
      WRITE(Updb,103) Num,Soluto(j)%nom,'','SOL','X',NumMol,'', &
                      Coords(MolSol(i)+j,:)/AngstromAtomica,1.0,0.0, &
                      ADJUSTR(Simbolo(Soluto(j)%z)),''
    END DO
  END DO
  NumMol2=NumMol2+MoleculasSoluto
  !Escribe los átomos del disolvente
  DO i=1,MoleculasDisolvente
    NumMol=NumMol+1
    IF (MolDis(i) < 0) CYCLE
    IF (Qcen) THEN
      Num=Num+1
      WRITE(Updb,103) Num,'X','','DIS','X',NumMol,'', &
                      CentroMasa(NumMol2+i,:)/AngstromAtomica,1.0,0.0,'X',''
    END IF
    DO j=1,SIZE(Disolvente,1)
      IF (Fuera(MolDis(i)+j)) CYCLE
      Num=Num+1
      WRITE(Updb,103) Num,Disolvente(j)%nom,'','DIS','X',NumMol,'', &
                      Coords(MolDis(i)+j,:)/AngstromAtomica,1.0,0.0, &
                      ADJUSTR(Simbolo(Disolvente(j)%z)),''
    END DO
  END DO
  NumMol2=NumMol2+MoleculasDisolvente
  !Escribe los átomos del disolvente 2 ...
  DO i=1,MoleculasDisolvente2
    NumMol=NumMol+1
    IF (MolDis2(i) < 0) CYCLE
    IF (Qcen) THEN
      Num=Num+1
      WRITE(Updb,103) Num,'X','','DI2','X',NumMol,'', &
                      CentroMasa(NumMol2+i,:)/AngstromAtomica,1.0,0.0,'X',''
    END IF
    DO j=1,SIZE(Disolvente2,1)
      IF (Fuera(MolDis2(i)+j)) CYCLE
      Num=Num+1
      WRITE(Updb,103) Num,Disolvente2(j)%nom,'','DI2','X',NumMol,'', &
                      Coords(MolDis2(i)+j,:)/AngstromAtomica,1.0,0.0, &
                      ADJUSTR(Simbolo(Disolvente2(j)%z)),''
    END DO
  END DO
  NumMol2=NumMol2+MoleculasDisolvente2

  !Fin
  WRITE(Updb,102)

100 FORMAT('Conf = ',I5,' Cent = ',I5)
101 FORMAT('TITLE ',T9,A2,T11,A60)
102 FORMAT('END   ')
103 FORMAT('HETATM',I5,T13,A4,A1,A3,T22,A1,I3,A1,T31,3(F8.3),2(F6.2),T77,A2,A2)

END SUBROUTINE EscribirPDB

!-------------------------------------------------------------------------------
! Escribe una configuración en formato XYZ
!-------------------------------------------------------------------------------
! Uxyz:    Fichero donde se escribe la configuración
! Centros: Variable que indica si se incluyen los centros de masa en la salida
! QCen:    Igual que Centros
! Num:     Número de átomos
! NumMol:  Número de moléculas
! i,j:     Contadores
!-------------------------------------------------------------------------------
SUBROUTINE EscribirXYZ(Uxyz,Centros)
  USE Parametros
  USE Unidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Uxyz
  LOGICAL, INTENT(IN), OPTIONAL :: Centros

  LOGICAL :: QCen
  INTEGER :: Num,NumMol,i,j

  QCen=.FALSE.

  !Se calcula el número real de átomos y se escribe la cabecera
  Num=0
  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    Num=COUNT(.NOT. Fuera(:))

   CASE (1) !Moldy
    IF (PRESENT(Centros)) QCen=Centros
    Num=MAX(Num,MAXVAL(MolSol(:))+SIZE(Soluto,1))
    Num=MAX(Num,MAXVAL(MolDis(:))+SIZE(Disolvente,1))
    Num=MAX(Num,MAXVAL(MolDis2(:))+SIZE(Disolvente2,1))
    IF (QCen) THEN
      Num=Num+COUNT(MolSol(:)>=0)
      Num=Num+COUNT(MolDis(:)>=0)
      Num=Num+COUNT(MolDis2(:)>=0)
    END IF
  END SELECT
  WRITE(Uxyz,100) Num
  WRITE(Uxyz,102) Config,Centro

  !Se escriben las coordenadas de los átomos
  NumMol=0
  DO i=1,MoleculasSoluto
    IF (MolSol(i) < 0) CYCLE
    IF (Qcen) WRITE(Uxyz,101) 'X',CentroMasa(NumMol+i,:)/AngstromAtomica
    DO j=1,SIZE(Soluto,1)
      IF (Fuera(MolSol(i)+j)) CYCLE
      WRITE(Uxyz,101) Soluto(j)%nom,Coords(MolSol(i)+j,:)/AngstromAtomica, &
                      Soluto(j)%q
    END DO
  END DO
  NumMol=NumMol+MoleculasSoluto
  DO i=1,MoleculasDisolvente
    IF (MolDis(i) < 0) CYCLE
    IF (Qcen) WRITE(Uxyz,101) 'X',CentroMasa(NumMol+i,:)/AngstromAtomica
    DO j=1,SIZE(Disolvente,1)
      IF (Fuera(MolDis(i)+j)) CYCLE
      WRITE(Uxyz,101) Disolvente(j)%nom,Coords(MolDis(i)+j,:)/AngstromAtomica, &
                      Disolvente(j)%q
    END DO
  END DO
  NumMol=NumMol+MoleculasDisolvente
  DO i=1,MoleculasDisolvente2
    IF (MolDis2(i) < 0) CYCLE
    IF (Qcen) WRITE(Uxyz,101) 'X',CentroMasa(NumMol+i,:)/AngstromAtomica
    DO j=1,SIZE(Disolvente2,1)
      IF (Fuera(MolDis2(i)+j)) CYCLE
      WRITE(Uxyz,101) Disolvente2(j)%nom,Coords(MolDis2(i)+j,:)/AngstromAtomica, &
                      Disolvente2(j)%q
    END DO
  END DO
  NumMol=NumMol+MoleculasDisolvente2

100 FORMAT(I5)
101 FORMAT(A8,4(1X,F10.6))
102 FORMAT('Conf = ',I5,' Cent = ',I5)

END SUBROUTINE EscribirXYZ

!-------------------------------------------------------------------------------
! Calcula las interacciones entre la molécula central y el resto
!-------------------------------------------------------------------------------
! Elec:    Interacción electrostática
! VdW:     Interacción de van der Waals
! MolCen:  Molécula central (sustituye a la real de cada configuración)
! GrVdW:   Componente de vdW del gradiente
! HsVdW:   Componente de vdW de la hessiana
! QGrad:   Indica si se ha de calcular el gradiente
! QHess:   Indica si se ha de calcular la hessiana
! Mol:     La molécula central
! At1,At2: Átomos entre los que se calcula la interacción
! ElecPar: Interacción electrostática entre dos átomos
! VdWPar:  Interacción de van der Waals entre dos átomos
! GradPar: Componente de vdW del gradiente para dos átomos
! HessPar: Componente de vdW de la hessiana para dos átomos
! Num:     Variable auxiliar para llevar las cuentas de las moléculas
! i,j,k:   Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Interacciones(Elec,VdW,MolCen,GrVdW,HsVdW)
  USE Parametros
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT) :: Elec,VdW
  TYPE(Atomo), DIMENSION(:), INTENT(IN), OPTIONAL :: MolCen
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT), OPTIONAL :: GrVdW
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT), OPTIONAL :: HsVdW

  TYPE(Atomo), DIMENSION(:), ALLOCATABLE :: Mol
  TYPE(Atomo) :: At1,At2
  DOUBLE PRECISION :: ElecPar,VdWPar
  DOUBLE PRECISION, DIMENSION(3) :: GradPar
  DOUBLE PRECISION, DIMENSION(3,3) :: HessPar
  LOGICAL :: QGrad,QHess
  INTEGER :: Num,i,j,k

  Elec=0.0D0
  VdW=0.0D0
  IF (PRESENT(GrVdW)) THEN
    GrVdW(:)=0.0D0
    QGrad=.TRUE.
   ELSE
    QGrad=.FALSE.
  END IF
  IF (PRESENT(HsVdW)) THEN
    HsVdW(:,:)=0.0D0
    QHess=.TRUE.
   ELSE
    QHess=.FALSE.
  END IF

  !Se calculan las constantes del "reaction field"
  IF (CorteRF < 1.0D3) THEN
    Krf=(Dielectrica-1.0D0)/(CorteRF**3*(2.0D0*Dielectrica+1.0D0))
    Crf=(3.0D0*Dielectrica)/(CorteRF*(2.0D0*Dielectrica+1.0D0))
   ELSE
    Krf=0.0D0
    Crf=0.0D0
  END IF

  !Si la configuración no esta centrada, no se calcula nada
  IF (Centro < 1) THEN
    RETURN
   !Si está centrada, se copia la molécula central
   ELSE IF (PRESENT(MolCen)) THEN
    ALLOCATE(Mol(SIZE(MolCen,1)))
    Mol(:)=MolCen(:)
   ELSE IF (Centro <= MoleculasSoluto) THEN
    ALLOCATE(Mol(SIZE(Soluto,1)))
    Mol(:)=Soluto(:)
   ELSE IF (Centro <= MoleculasSoluto+MoleculasDisolvente) THEN
    ALLOCATE(Mol(SIZE(Disolvente,1)))
    Mol(:)=Disolvente(:)
   ELSE IF (Centro <= MoleculasSoluto+MoleculasDisolvente+ &
                      MoleculasDisolvente2) THEN
    ALLOCATE(Mol(SIZE(Disolvente2,1)))
    Mol(:)=Disolvente2(:)
  END IF

  !Para cada átomo de la molécula centrada
  DO k=1,SIZE(Mol,1)
    At1=Mol(k)
    Num=0
    !Se calcula la interacción con los átomos de las otras moléculas
!$OMP PARALLEL PRIVATE(i,j,At2,ElecPar,VdWPar,GradPar,HessPar) &
!$OMP          REDUCTION(+:Elec,VdW)
!$OMP DO
    DO i=1,MoleculasSoluto
      !La molécula centrada no se considera
      !Tampoco las que están fuera del corte
      IF ((Num+i == Centro) .OR. (MolSol(i) < 0)) CYCLE
      DO j=1,SIZE(Soluto,1)
        !Si el átomo está omitido (fuera de la esfera) se salta
        IF (Fuera(MolSol(i)+j)) CYCLE
        At2=Soluto(j)
        At2%pos(:)=Coords(MolSol(i)+j,:)
        CALL InteraccionPar(At1,At2,ElecPar,VdWPar,QGrad,GradPar,QHess,HessPar)
        Elec=Elec+ElecPar
        VdW=VdW+VdWPar
!$OMP CRITICAL
        IF (QGrad) GrVdW(k*3-2:k*3)=GrVdW(k*3-2:k*3)+GradPar
        IF (QHess) HsVdW(k*3-2:k*3,k*3-2:k*3)=HsVdW(k*3-2:k*3,k*3-2:k*3)+HessPar
!$OMP END CRITICAL
      END DO
    END DO
!$OMP END DO
!$OMP SINGLE
    Num=Num+MoleculasSoluto
!$OMP END SINGLE
!$OMP DO
    DO i=1,MoleculasDisolvente
      IF ((Num+i == Centro) .OR. (MolDis(i) < 0)) CYCLE
      DO j=1,SIZE(Disolvente,1)
        IF (Fuera(MolDis(i)+j)) CYCLE
        At2=Disolvente(j)
        At2%pos(:)=Coords(MolDis(i)+j,:)
        CALL InteraccionPar(At1,At2,ElecPar,VdWPar,QGrad,GradPar,QHess,HessPar)
        Elec=Elec+ElecPar
        VdW=VdW+VdWPar
!$OMP CRITICAL
        IF (QGrad) GrVdW(k*3-2:k*3)=GrVdW(k*3-2:k*3)+GradPar
        IF (QHess) HsVdW(k*3-2:k*3,k*3-2:k*3)=HsVdW(k*3-2:k*3,k*3-2:k*3)+HessPar
!$OMP END CRITICAL
      END DO
    END DO
!$OMP END DO
!$OMP SINGLE
    Num=Num+MoleculasDisolvente
!$OMP END SINGLE
!$OMP DO
    DO i=1,MoleculasDisolvente2
      IF ((Num+i == Centro) .OR. (MolDis2(i) < 0)) CYCLE
      DO j=1,SIZE(Disolvente2,1)
        IF (Fuera(MolDis2(i)+j)) CYCLE
        At2=Disolvente2(j)
        At2%pos(:)=Coords(MolDis2(i)+j,:)
        CALL InteraccionPar(At1,At2,ElecPar,VdWPar,QGrad,GradPar,QHess,HessPar)
        Elec=Elec+ElecPar
        VdW=VdW+VdWPar
!$OMP CRITICAL
        IF (QGrad) GrVdW(k*3-2:k*3)=GrVdW(k*3-2:k*3)+GradPar
        IF (QHess) HsVdW(k*3-2:k*3,k*3-2:k*3)=HsVdW(k*3-2:k*3,k*3-2:k*3)+HessPar
!$OMP END CRITICAL
      END DO
    END DO
!$OMP END DO
!$OMP SINGLE
    Num=Num+MoleculasDisolvente2
!$OMP END SINGLE
!$OMP END PARALLEL
  END DO

END SUBROUTINE Interacciones

!-------------------------------------------------------------------------------
! Calcula la interacción entre dos átomos
!-------------------------------------------------------------------------------
! At1,At2: Los dos átomos implicados
! Elec:    Energía electrostática entre los dos átomos
! VdW:     Energía de van der Waals entre los dos átomos
! QGrad:   Indica si se ha de calcular el gradiente
! Grad:    Componente de vdW del gradiente sobre At1
! QHess    Indica si se ha de calcular la hessiana
! Hess:    Componente de vdW de la hessiana sobre At1
! Vect:    Vector que une los dos átomos
! Dist:    Distancia entre los dos átomos
! Dist2:   Distancia al cuadrado
! Eps,Sig6,Sig12,Par1..6,Aux1,Aux2: Variables auxiliares para el cálculo vdW
! Unidad:  Matriz unidad 3x3
!-------------------------------------------------------------------------------
SUBROUTINE InteraccionPar(At1,At2,Elec,VdW,QGrad,Grad,QHess,Hess)
  USE Parametros
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), INTENT(IN) :: At1,At2
  DOUBLE PRECISION, INTENT(OUT) :: Elec,VdW
  LOGICAL, INTENT(IN) :: QGrad,QHess
  DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT) :: Grad
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(INOUT) :: Hess

  DOUBLE PRECISION, DIMENSION(3) :: Vect
  DOUBLE PRECISION :: Dist,Dist2,Eps,Sig6,Sig12,Par1,Par2,Par3,Par4,Par5, &
                      Par6,Aux1,Aux2
  DOUBLE PRECISION, DIMENSION(3,3), PARAMETER :: &
    Unidad=RESHAPE((/1.0D0,0.0D0,0.0D0, &
                     0.0D0,1.0D0,0.0D0, &
                     0.0D0,0.0D0,1.0D0 /), (/3,3/))

  !Calcula el vector y la distancia de separación
  Vect(:)=At1%pos(:)-At2%pos(:)
  Dist2=DOT_PRODUCT(Vect(:),Vect(:))
  Dist=SQRT(Dist2)

  !La interacción electrostática se calcula con "reaction field"
  Elec=0.0D0
  IF (Dist < CorteRF) Elec=At1%q*At2%q*(1.0D0/Dist+Krf*Dist2-Crf)

  !La interacción de van der Waals se calcula según el tipo de potencial
  !y sólo si está definido
  VdW=0.0D0
  IF (QGrad) Grad(:)=0.0D0
  IF (QHess) Hess(:,:)=0.0D0
  IF (QInter(At1%id,At2%id)) THEN
    SELECT CASE (TipoPotencial)
     CASE (1) !lennard-jones
      Eps=InterAtom(At1%id,At2%id,1)
      Sig6=(InterAtom(At1%id,At2%id,2)**2/Dist2)**3
      Sig12=Sig6*Sig6
      VdW=Eps*(Sig12-Sig6)
      IF (QGrad .OR. QHess) THEN
        Aux1=Eps*(6.0D0*Sig6-12.0D0*Sig12)/Dist2
        IF (QHess) Aux2=Eps*(168.0D0*Sig12-48.0D0*Sig6)/(Dist2*Dist2)
      END IF
     CASE (2) !generic
      Par1=InterAtom(At1%id,At2%id,1)
      Par2=EXP(-InterAtom(At1%id,At2%id,2)*Dist)
      Par3=InterAtom(At1%id,At2%id,3)/(Dist2**6)
      Par4=InterAtom(At1%id,At2%id,4)/(Dist2**2)
      Par5=InterAtom(At1%id,At2%id,5)/(Dist2**3)
      Par6=InterAtom(At1%id,At2%id,6)/(Dist2**4)
      VdW=Par1*Par2+Par3-Par4-Par5-Par6
      IF (QGrad .OR. QHess) THEN
        Aux1=-Par1*Par2*InterAtom(At1%id,At2%id,2)/Dist- &
             (12.0D0*Par3-4.0D0*Par4-6.0D0*Par5-8.0D0*Par6)/Dist2
        IF (QHess) Aux2=-Par1*Par2*InterAtom(At1%id,At2%id,2)**2/Dist**3+ &
             (168.0D0*Par3-24.0D0*Par4-48.0D0*Par5-80D0*Par6)/(Dist2*Dist2)
      END IF
    END SELECT
    !La fórmula básica del gradiente y la hessiana es igual
    IF (QGrad) THEN
      Grad(:)=Aux1*Vect(:)
    END IF
    IF (QHess) THEN
      Hess(:,:)=Aux2*ProductoTens(Vect(:),Vect(:))+Aux1*Unidad(:,:)
    END IF
  END IF

END SUBROUTINE InteraccionPar

!-------------------------------------------------------------------------------
! Calcula las energías de interacción promedio
!-------------------------------------------------------------------------------
! Der:   Orden de la derivada de la energía de vdW que se quiere calcular
! Mol:   Molécula con respecto a la cual se calculan las interacciones
! Elec:  Energía electrostática temporal
! VdW:   Energía de vdW temporal
! GVdW:  Gradiente de vdW temporal
! HVdW:  Hessiana de vdW temporal
! Num:   Número de configuraciones
! Error: Variable para controlar los errores
!-------------------------------------------------------------------------------
SUBROUTINE Promedios(Der,Mol)
  USE DatosQM
  USE DatosQMMM
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: Der
  TYPE(Atomo), DIMENSION(:), INTENT(IN), OPTIONAL :: Mol
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: GVdW
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: HVdW
  DOUBLE PRECISION :: Elec,VdW
  INTEGER :: Num,Error

  ! Se inician las variables
  IF (ALLOCATED(GradVdW)) DEALLOCATE(GradVdW)
  IF (ALLOCATED(HessVdW)) DEALLOCATE(HessVdW)
  Num=3*SIZE(Mol,1)
  ALLOCATE(GradVdW(Num),HessVdW(Num,Num))
  GradVdW(:)=0.0D0
  HessVdW(:,:)=0.0D0
  ALLOCATE(GVdW(Num),HVdW(Num,Num))

  ! La energía de interacción QM con las cargas es fácil de calcular
  IF (ALLOCATED(DisolvQM)) EnergiaEQM=SUM(DisolvQM(:,4)*DisolvQM(:,5))
  EnergiaEMM=0.0D0
  EnergiaVdW=0.0D0

  ! Se van leyendo configuraciones y se calculan las interacciones
  Num=0
  CALL AbrirUConf()
  DO
    CALL LeerConfig(Error)
    IF (Error /= 0) EXIT
    IF (PRESENT(Mol)) THEN
      SELECT CASE(Der)
       CASE (0)
        CALL Interacciones(Elec,VdW,Mol)
       CASE (1)
        CALL Interacciones(Elec,VdW,Mol,GVdW(:))
       CASE (2)
        CALL Interacciones(Elec,VdW,Mol,GVdW(:),HVdW(:,:))
      END SELECT
     ELSE
      SELECT CASE(Der)
       CASE (0)
        CALL Interacciones(Elec,VdW)
       CASE (1)
        CALL Interacciones(Elec,VdW,GrVdW=GVdW(:))
       CASE (2)
        CALL Interacciones(Elec,VdW,GrVdW=GVdW(:),HsVdW=HVdW(:,:))
      END SELECT
    END IF
    ! Se acumulan los resultados
    EnergiaEMM=EnergiaEMM+Elec
    EnergiaVdW=EnergiaVdW+VdW
    GradVdW(:)=GradVdW(:)+GVdW(:)
    HessVdW(:,:)=HessVdW(:,:)+HVdW(:,:)
    Num=Num+1
  END DO

  ! Finalmente se calculan los promedios y se almacenan en sus variables
  IF (Num > 0) THEN
    EnergiaEMM=EnergiaEMM/DBLE(Num)
    EnergiaVdW=EnergiaVdW/DBLE(Num)
    GradVdW(:)=GradVdW(:)/DBLE(Num)
    HessVdW(:,:)=HessVdW(:,:)/DBLE(Num)
  END IF

  DEALLOCATE(GVdW,HVdW)

END SUBROUTINE Promedios

!-------------------------------------------------------------------------------
! Calcula el potencial electrostático de una configuración en una serie de
! puntos
!-------------------------------------------------------------------------------
! Puntos:    Matriz de puntos donde se ha de calcular el potencial
! Potencial: Vector donde se guarda el potencial calculado
! Dist:      Distancia entre cada par átomo-punto
! Num:       Número de átomos
! i,j,k:     Contadores
!-------------------------------------------------------------------------------
SUBROUTINE PotencialMM(Puntos,Potencial)
  USE Parametros
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Puntos
  DOUBLE PRECISION, DIMENSION(SIZE(Puntos,1)), INTENT(OUT) :: Potencial

  DOUBLE PRECISION :: Dist
  INTEGER :: Num,i,j,k

  !Se calculan las constantes del "reaction field"
  IF (CorteRF < 1.0D3) THEN
    Krf=(Dielectrica-1.0D0)/(CorteRF**3*(2.0D0*Dielectrica+1.0D0))
    Crf=(3.0D0*Dielectrica)/(CorteRF*(2.0D0*Dielectrica+1.0D0))
   ELSE
    Krf=0.0D0
    Crf=0.0D0
  END IF

  Potencial(:)=0.0D0
!$OMP PARALLEL DO PRIVATE(Num,Dist,i,j,k)
  DO k=1,SIZE(Potencial)
    Num=0
    DO i=1,MoleculasSoluto
      !La molécula centrada no se considera
      !Tampoco las que están fuera del corte
      IF ((Num+i == Centro) .OR. (MolSol(i) < 0)) CYCLE
      DO j=1,SIZE(Soluto,1)
        !Si el átomo está omitido (fuera de la esfera) se salta
        IF (Fuera(MolSol(i)+j)) CYCLE
        Dist=Distancia(Coords(MolSol(i)+j,:),Puntos(k,:))
        IF (Dist < CorteRF) Potencial(k)=Potencial(k)+ &
                            Soluto(j)%q*(1.0D0/Dist+Krf*Dist*Dist-Crf)
      END DO
    END DO
    Num=Num+MoleculasSoluto
    DO i=1,MoleculasDisolvente
      IF ((Num+i == Centro) .OR. (MolDis(i) < 0)) CYCLE
      DO j=1,SIZE(Disolvente,1)
        IF (Fuera(MolDis(i)+j)) CYCLE
        Dist=Distancia(Coords(MolDis(i)+j,:),Puntos(k,:))
        IF (Dist < CorteRF) Potencial(k)=Potencial(k)+ &
                            Disolvente(j)%q*(1.0D0/Dist+Krf*Dist*Dist-Crf)
      END DO
    END DO
    Num=Num+MoleculasDisolvente
    DO i=1,MoleculasDisolvente2
      IF ((Num+i == Centro) .OR. (MolDis2(i) < 0)) CYCLE
      DO j=1,SIZE(Disolvente2,1)
        IF (Fuera(MolDis2(i)+j)) CYCLE
        Dist=Distancia(Coords(MolDis2(i)+j,:),Puntos(k,:))
        IF (Dist < CorteRF) Potencial(k)=Potencial(k)+ &
                            Disolvente2(j)%q*(1.0D0/Dist+Krf*Dist*Dist-Crf)
      END DO
    END DO
    Num=Num+MoleculasDisolvente2
  END DO
!$OMP END PARALLEL DO

END SUBROUTINE PotencialMM

!-------------------------------------------------------------------------------
! Selecciona (escribe en un fichero) las cargas correspondientes a las moléculas
! que quedan dentro de la cavidad para una configuración dada
!-------------------------------------------------------------------------------
! U:         Unidad donde se escriben las cargas
! NumCargas: Número acumulativo de cargas escritas en U
! Num:       Número de átomos
! Aux:       Vector con las coordenadas de un átomo
! i,j:       Contadores
!-------------------------------------------------------------------------------
SUBROUTINE SeleccionarMols(U,NumCargas)
  USE Cavidad
  USE Sistema
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: U
  INTEGER, INTENT(INOUT) :: NumCargas

  DOUBLE PRECISION, DIMENSION(3) :: Aux
  INTEGER :: Num,i,j

  Num=0
  DO i=1,MoleculasSoluto
    IF ((Num+i == Centro) .OR. (MolSol(i) < 0)) CYCLE
    IF (.NOT. Interior(CentroMasa(Num+i,:))) CYCLE
    DO j=1,SIZE(Soluto,1)
      !IF (Fuera(MolSol(i)+j)) CYCLE
      Aux(:)=Coords(MolSol(i)+j,:)
      IF (.NOT. Interior(Aux(:))) CYCLE
      WRITE(U) Aux(:),Soluto(j)%q
      NumCargas=NumCargas+1
    END DO
  END DO
  Num=Num+MoleculasSoluto
  DO i=1,MoleculasDisolvente
    IF ((Num+i == Centro) .OR. (MolDis(i) < 0)) CYCLE
    IF (.NOT. Interior(CentroMasa(Num+i,:))) CYCLE
    DO j=1,SIZE(Disolvente,1)
      !IF (Fuera(MolDis(i)+j)) CYCLE
      Aux(:)=Coords(MolDis(i)+j,:)
      IF (.NOT. Interior(Aux(:))) CYCLE
      WRITE(U) Aux(:),Disolvente(j)%q
      NumCargas=NumCargas+1
    END DO
  END DO
  Num=Num+MoleculasDisolvente
  DO i=1,MoleculasDisolvente2
    IF ((Num+i == Centro) .OR. (MolDis2(i) < 0)) CYCLE
    IF (.NOT. Interior(CentroMasa(Num+i,:))) CYCLE
    DO j=1,SIZE(Disolvente2,1)
      !IF (Fuera(MolDis2(i)+j)) CYCLE
      Aux(:)=Coords(MolDis2(i)+j,:)
      IF (.NOT. Interior(Aux(:))) CYCLE
      WRITE(U) Aux(:),Disolvente2(j)%q
      NumCargas=NumCargas+1
    END DO
  END DO
  Num=Num+MoleculasDisolvente2

END SUBROUTINE SeleccionarMols

!-------------------------------------------------------------------------------
! Reduce el número de cargas juntándolas o ajustándolas a una red tridimensional
!-------------------------------------------------------------------------------
! U:         Unidad donde están las cargas de entrada
! Num:       Número de cargas de entrada
! Cargas:    Matriz con las cargas de salida
! CargasIni: Matriz con las cargas de entrada
! Puntos:    Matriz que asigna una carga a cada punto de la red
! Quitar:    Matriz que establece las cargas que se conservan o se eliminan
! NumRed:    Número de cargas una vez reducidas
! Preci:     Valor mínimo de las cargas que se conservan
! Punto:     Posición de cada carga
! Dif:       Vector diferencia entre la carga y el punto de la red
! Dist:      Distancia entre dos puntos
! R:         Distancia entre los puntos de la red cúbica base
! Vecinos:   Matriz que establece los puntos vecinos a uno dado
! Lim:       Índices mínimo y máximo de la red
! Indice:    Índices del punto de la red más cercano a una carga
! i,j,k,l,ii,jj,kk,n,m: Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ReducirCargas(U,Num,Cargas)
  USE Parametros
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: U,Num
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Cargas

  DOUBLE PRECISION, DIMENSION(Num,5) :: CargasIni
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Puntos
  LOGICAL, DIMENSION(Num) :: Quitar
  DOUBLE PRECISION, DIMENSION(3) :: Dif,Punto
  DOUBLE PRECISION :: Dist,R
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-8
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Vecinos
  INTEGER, DIMENSION(3,2) :: Lim
  INTEGER, DIMENSION(3) :: Indice
  INTEGER :: i,j,k,l,ii,jj,kk,n,m,NumRed

  Quitar(:)=.FALSE.

  !Se leen las cargas iniciales
  REWIND(U)
  DO i=1,Num
    READ(U) CargasIni(i,1:4)
  END DO

  !La reducción es distinta si se hace con ajuste a una red o no
  IF (TipoReduccion < 0) THEN

    !Para cada carga
    DO i=1,Num
      IF (Quitar(i)) CYCLE
      Punto(:)=CargasIni(i,1:3)
      !Une todas las siguientes cargas que disten menos de un valor dado
!$OMP PARALLEL DO PRIVATE(j)
      DO j=i+1,Num
        IF (Quitar(j)) CYCLE
        IF (Distancia(CargasIni(j,1:3),Punto(:)) < DistCargas) THEN
          Quitar(j)=.TRUE.
!$OMP CRITICAL
          CargasIni(i,4)=CargasIni(i,4)+CargasIni(j,4)
!$OMP END CRITICAL
        END IF
      END DO
!$OMP END PARALLEL DO
      !Si el valor de la carga es muy pequeño, se elimina
      IF (ABS(CargasIni(i,4)) < Preci) Quitar(i)=.TRUE.
    END DO

   ELSE

    !Se calculan la distancia de la red cubica base y los vecinos más próximos
    !según el tipo de malla
    SELECT CASE (TipoReduccion)
     CASE (0) !Red cúbica simple
      R=DistCargas
      ALLOCATE(Vecinos(13,3))
      !Caras
      Vecinos(1,:)=(/ 1, 0, 0/)
      Vecinos(2,:)=(/ 0, 1, 0/)
      Vecinos(3,:)=(/ 0, 0, 1/)
      !Aristas
      Vecinos(4,:)=(/ 1, 1, 0/)
      Vecinos(5,:)=(/ 1,-1, 0/)
      Vecinos(6,:)=(/ 1, 0, 1/)
      Vecinos(7,:)=(/ 1, 0,-1/)
      Vecinos(8,:)=(/ 0, 1, 1/)
      Vecinos(9,:)=(/ 0, 1,-1/)
      !Vértices
      Vecinos(10,:)=(/ 1, 1, 1/)
      Vecinos(11,:)=(/-1, 1, 1/)
      Vecinos(12,:)=(/ 1,-1, 1/)
      Vecinos(13,:)=(/ 1, 1,-1/)
     CASE (1) !Red cúbica compacta
      R=DistCargas/SQRT(2.0D0)
      ALLOCATE(Vecinos(9,3))
      !Caras
      Vecinos(1,:)=(/ 1, 1, 0/)
      Vecinos(2,:)=(/ 1,-1, 0/)
      Vecinos(3,:)=(/ 1, 0, 1/)
      Vecinos(4,:)=(/ 1, 0,-1/)
      Vecinos(5,:)=(/ 0, 1, 1/)
      Vecinos(6,:)=(/ 0, 1,-1/)
      !Vértices
      Vecinos(7,:)=(/ 2, 0, 0/)
      Vecinos(8,:)=(/ 0, 2, 0/)
      Vecinos(9,:)=(/ 0, 0, 2/)
     CASE (2) !Red cúbica centrada en el cuerpo
      R=DistCargas/SQRT(3.0D0)
      ALLOCATE(Vecinos(7,3))
      !Caras
      Vecinos(1,:)=(/ 2, 0, 0/)
      Vecinos(2,:)=(/ 0, 2, 0/)
      Vecinos(3,:)=(/ 0, 0, 2/)
      Vecinos(4,:)=(/ 1, 1, 1/)
      Vecinos(5,:)=(/-1, 1, 1/)
      Vecinos(6,:)=(/ 1,-1, 1/)
      Vecinos(7,:)=(/ 1, 1,-1/)
     CASE DEFAULT !Sin red
    END SELECT

    !Se aumenta la distancia para dar cierta holgura
    R=1.75D0*R

    !Se calculan los límites de la red
    DO i=1,3
      Lim(i,1)=FLOOR(MINVAL(CargasIni(:,i))/R)
      Lim(i,2)=CEILING(MAXVAL(CargasIni(:,i))/R)
    END DO

    ALLOCATE(Puntos(Lim(1,1):Lim(1,2),Lim(2,1):Lim(2,2),Lim(3,1):Lim(3,2)))
    Puntos(:,:,:)=0

    !Se asigna cada carga a una celda, se suman todas las cargas de cada celda y
    !se mantiene la posición de la más cercana al centro
!$OMP PARALLEL DO PRIVATE(i,j,Indice,Dif,Dist)
    DO i=1,Num
      !Se halla el punto de la red más cercano a cada carga
      SELECT CASE (TipoReduccion)
       CASE (1) !Red cúbica compacta
        Indice(:)=NINT(CargasIni(i,1:3)/R)
        IF (MOD(SUM(Indice(:)),2) /= 0) THEN
          Dif(:)=CargasIni(i,1:3)/R-Indice(:)
          j=SUM(MAXLOC(ABS(Dif(:))))
          Indice(j)=Indice(j)+NINT(SIGN(1.0D0,Dif(j)))
        END IF
       CASE (2) !Red cúbica centrada en el cuerpo
        Indice(:)=2*NINT(0.5D0*CargasIni(i,1:3)/R)
        Dif(:)=CargasIni(i,1:3)/R-Indice(:)
        IF (SUM(ABS(Dif(:))) > 1.5D0) &
          Indice(:)=2*NINT(0.5D0*(CargasIni(i,1:3)/R-(/1.0D0,1.0D0,1.0D0/))) + &
                    (/1,1,1/)
       CASE DEFAULT !Red cúbica simple
        Indice(:)=NINT(CargasIni(i,1:3)/R)
      END SELECT
      Dist=Distancia(CargasIni(i,1:3),Indice(:)*R)
      j=Puntos(Indice(1),Indice(2),Indice(3))
!$OMP CRITICAL
      IF (j == 0) THEN
        !Si no hay carga en esa celda, se pone esta
        Puntos(Indice(1),Indice(2),Indice(3))=i
        CargasIni(i,5)=Dist
       ELSE IF (Dist < CargasIni(j,5)) THEN
        !Si hay una carga, pero esta está más cerca, se sustituye y se suman los
        !valores
        Puntos(Indice(1),Indice(2),Indice(3))=i
        CargasIni(i,5)=Dist
        CargasIni(i,4)=CargasIni(i,4)+CargasIni(j,4)
        Quitar(j)=.TRUE.
       ELSE
        !Si hay una carga y es la más cercana, se suman los valores
        CargasIni(j,4)=CargasIni(j,4)+CargasIni(i,4)
        Quitar(i)=.TRUE.
      END IF
!$OMP END CRITICAL
    END DO
!$OMP END PARALLEL DO

    !Se eliminan las cargas que puedan quedar demasiado cercanas a alguna otra
    DO i=Lim(1,1),Lim(1,2)
      DO j=Lim(2,1),Lim(2,2)
        DO k=Lim(3,1),Lim(3,2)
          n=Puntos(i,j,k)
          IF (n == 0) CYCLE
          !Para cada punto de la red, se cumprueban las distancias con las
          !cargas de los puntos vecinos
          DO l=1,SIZE(Vecinos,1)
            ii=i+Vecinos(l,1)
            jj=j+Vecinos(l,2)
            kk=k+Vecinos(l,3)
            IF ((ii < Lim(1,1)) .OR. (ii > Lim(1,2)) .OR. &
                (jj < Lim(2,1)) .OR. (jj > Lim(2,2)) .OR. &
                (kk < Lim(3,1)) .OR. (kk > Lim(3,2))) CYCLE
            m=Puntos(ii,jj,kk)
            IF (m == 0) CYCLE
            Dist=Distancia(CargasIni(n,1:3),CargasIni(m,1:3))
            !Si la distancia es menor de lo deseado, se elimina la carga más
            !alejada de su centro
           IF (Dist < DistCargas) THEN
              IF (CargasIni(n,5) < CargasIni(m,5)) THEN
                CargasIni(n,4)=CargasIni(n,4)+CargasIni(m,4)
                Quitar(m)=.TRUE.
               ELSE
                CargasIni(m,4)=CargasIni(m,4)+CargasIni(n,4)
                Quitar(n)=.TRUE.
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO

    !Por ultimo se eliminan las cargas que puedan quedar con valor mínimo
    DO i=1,Num
      IF (Quitar(i)) CYCLE
      IF (ABS(CargasIni(i,4)) < Preci) Quitar(i)=.TRUE.
    END DO

    DEALLOCATE(Puntos,Vecinos)

  END IF

  !La salida contiene sólo las cargas que no se han eliminado
  NumRed=COUNT(.NOT. Quitar)
  IF (ALLOCATED(Cargas)) DEALLOCATE(Cargas)
  ALLOCATE(Cargas(NumRed,4))
  j=0
  DO i=1,Num
    IF (Quitar(i)) CYCLE
    j=j+1
    Cargas(j,:)=CargasIni(i,1:4)
  END DO

END SUBROUTINE ReducirCargas

END MODULE Configuraciones
