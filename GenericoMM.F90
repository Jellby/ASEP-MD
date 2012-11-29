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

MODULE GenericoMM
USE DCD

CONTAINS
!LeerSistemaGenerico(Fich)
!LeerConfigsGenerico()
!EntradaGenericoMM(Sal)

!-------------------------------------------------------------------------------
! Esta subrutina lee los datos del sistema en un formato específico
!-------------------------------------------------------------------------------
! Fich:     Numero del fichero de entrada
! Aux:      Vector auxiliar para leer los parámetros
! Renum:    Matriz para renumerar los átomos y que sean seguidos
! Linea:    La línea que se va leyendo
! Num:      Número de elementos en cada bloque
! Fin:      Variable para controlar el fin del fichero
! i,j,k,l:  Contadores
!-------------------------------------------------------------------------------
SUBROUTINE LeerSistemaGenerico(Fich)
  USE Parametros
  USE Sistema
  USE Unidades
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Aux
  INTEGER, DIMENSION(:), ALLOCATABLE :: Renum
  CHARACTER (LEN=LLL) :: Linea
  INTEGER :: Num,Fin,i,j,k,l

  IF (ALLOCATED(Soluto)) DEALLOCATE(Soluto)
  IF (ALLOCATED(Disolvente)) DEALLOCATE(Disolvente)
  IF (ALLOCATED(Disolvente2)) DEALLOCATE(Disolvente2)
  IF (ALLOCATED(InterAtom)) DEALLOCATE(InterAtom,QInter)
  ALLOCATE(Disolvente2(0))

  DO
    CALL LeerSiguienteLinea(Fich,Linea,Fin)
    IF (Fin /= 0) EXIT
    Linea=ADJUSTL(Linea)
    CALL PasarMinusculas(Linea)

    SELECT CASE (TRIM(Linea))

     !Primero debe aparecere la palabra "Solute", seguida de un nombre, el
     !número de moléculas y el número de átomos
     CASE ('solute')
      CALL LeerSiguienteLinea(Fich,Linea,Fin)
      IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
      READ(Linea,*) NombreSoluto
      CALL LeerSiguienteLinea(Fich,Linea,Fin)
      IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
      READ(Linea,*) MoleculasSoluto
      IF (MoleculasSoluto > 1) CALL Mensaje('LeerSistemaGenerico',38,.TRUE.)
      CALL LeerSiguienteLinea(Fich,Linea,Fin)
      IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
      READ(Linea,*) Num
      ALLOCATE(Soluto(Num))
      !Después los átomos del soluto, para cada uno: id, nombre
      !número atómico, masa (en uma), coordenadas (en angstrom) y
      !carga (en electrones)
      DO i=1,Num
        CALL LeerSiguienteLinea(Fich,Linea,Fin)
        IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
        READ(Linea,*) Soluto(i)%id,Soluto(i)%nom,Soluto(i)%z,Soluto(i)%m, &
                      Soluto(i)%pos(:),Soluto(i)%q
        Soluto(i)%m=Soluto(i)%m*AmuAtomica
        Soluto(i)%pos(:)=Soluto(i)%pos(:)*AngstromAtomica
      END DO

     !Después debe aparecer la palabra "Solvent", seguida del número de
     !átomos total del resto del sistema
     CASE ('solvent')
      NombreDisolvente='Solvent'
      MoleculasDisolvente=1
      CALL LeerSiguienteLinea(Fich,Linea,Fin)
      IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
      READ(Linea,*) Num
      ALLOCATE(Disolvente(Num))
      !Para cada átomo: id, nombre y carga (en electrones)
      DO i=1,Num
        READ(Fich,*) Disolvente(i)%id,Disolvente(i)%nom,Disolvente(i)%q
        Disolvente(i)%pos(:)=0.0D0
      END DO
      Disolvente(:)%m=1.0D0
      Disolvente(:)%z=0

     !Finalmente, debe aparecer la palabra "Non-Bonded", seguida del tipo
     !de potencial de Van der Waals y el número de interacciones
     CASE ('non-bonded')
      !Se renumeran los átomos para que estén seguidos
      k=MAX(MAXVAL(Soluto(:)%id),MAXVAL(Disolvente(:)%id))
      ALLOCATE(Renum(k))
      Renum(:)=0
      j=0
      k=0
      DO i=1,SIZE(Soluto,1)
        IF (Renum(Soluto(i)%id) == 0) THEN
          k=k+1
          Renum(Soluto(i)%id)=k
        END IF
        Soluto(i)%id=Renum(Soluto(i)%id)
      END DO
      DO i=1,SIZE(Disolvente,1)
        IF (Renum(Disolvente(i)%id) == 0) THEN
          k=k+1
          Renum(Disolvente(i)%id)=k
        END IF
        Disolvente(i)%id=Renum(Disolvente(i)%id)
      END DO
      j=MAXVAL(Renum)

      CALL LeerSiguienteLinea(Fich,Linea,Fin)
      IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
      CALL PasarMinusculas(Linea)
      SELECT CASE (TRIM(Linea))
       CASE ('lennard-jones')
        TipoPotencial=1
        i=2
       CASE ('generic')
        TipoPotencial=2
        i=6
       CASE DEFAULT
        CALL Mensaje('LeerSistemaGenerico',7,.TRUE.)
      END SELECT
      ALLOCATE(QInter(j,j),InterAtom(j,j,i),Aux(i))
      QInter(:,:)=.FALSE.
      InterAtom(:,:,:)=0.0D0
      CALL LeerSiguienteLinea(Fich,Linea,Fin)
      IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
      READ(Linea,*) Num
      !Para cada interacción: los dos id de los átomos y los parámetros
      !correspondientes (en unidades atómicas)
      DO i=1,Num
        CALL LeerSiguienteLinea(Fich,Linea,Fin)
        IF (Fin /= 0) CALL Mensaje('LeerSistemaGenerico',37,.TRUE.)
        READ(Linea,*) j,k,(Aux(l),l=1,SIZE(Aux,1))
        IF ((Renum(j)*Renum(k) == 0) .OR. (MAX(j,k) > SIZE(Renum,1))) &
          CALL Mensaje('LeerSistemaGenerico',18,.TRUE.)
        j=Renum(j)
        k=Renum(k)
        QInter(j,k)=.TRUE.
        QInter(k,j)=.TRUE.
        InterAtom(j,k,:)=Aux(:)
        InterAtom(k,j,:)=InterAtom(j,k,:)
      END DO
      DEALLOCATE(Renum,Aux)
    END SELECT
  END DO

END SUBROUTINE LeerSistemaGenerico

!-------------------------------------------------------------------------------
! Lee las configuraciones de un fichero DCD y las escribe en el fichero UConf
!-------------------------------------------------------------------------------
! Fich:    Nombre del fichero DCD que se lee (opcional)
!-------------------------------------------------------------------------------
SUBROUTINE LeerConfigsGenerico(Fich)
  USE Configuraciones
  USE Parametros
  USE Sistema
  USE Unidades
  USE Utilidades
  USE UtilidadesFis
  USE TipoAtomo
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: Fich

  TYPE(TipoDCD) :: Tray
  TYPE(Atomo), DIMENSION(:), ALLOCATABLE :: SolConf
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Pos
  DOUBLE PRECISION, DIMENSION(7) :: Mover
  DOUBLE PRECISION :: Paso
  INTEGER :: Num,Conf,i,j,k

  Num=0
  Num=Num+MoleculasSoluto*SIZE(Soluto,1)
  Num=Num+MoleculasDisolvente*SIZE(Disolvente,1)
  Num=Num+MoleculasDisolvente2*SIZE(Disolvente2,1)

  !Si no se especifica ningún fichero se usa TrayectoriaMM
  IF (PRESENT(Fich)) THEN
    Tray%Nombre=TRIM(Fich)
   ELSE
    Tray%Nombre=TRIM(TrayectoriaMM)
  END IF
  CALL AbrirDCD(Tray)
  IF (Tray%NAtomos /= Num) CALL Mensaje('LeerConfigsGenerico',39,.TRUE.)
  IF (Tray%Configs < NumConfig) CALL Mensaje('LeerConfigsGenerico',33,.TRUE.)

  ALLOCATE(SolConf(SIZE(Soluto,1)),Pos(Num,3))
  SolConf(:)=Soluto(:)

  CALL AbrirUConf()

  Paso=DBLE(Tray%Configs-1)/DBLE(MAX(NumConfig-1,1))
  DO i=1,NumConfig
    Conf=NINT((i-1)*Paso+1)
    CALL LeerDCD(Tray,Conf)
    Pos(:,:)=Tray%Coords(:,:)*AngstromAtomica

    !Coordenadas del soluto en la configuración
    DO j=1,SIZE(Soluto,1)
      SolConf(j)%pos(:)=Pos(j,:)
    END DO
    !Transforma la configuración para superponer el soluto
    CALL SuperponerMoleculas(Soluto,SolConf,Trans=Mover)
    Pos(:,:)=RotarCuaternion(Pos(:,:),Mover(1:4))
    Pos(:,:)=Pos(:,:)+SPREAD(Mover(5:7),DIM=1,NCOPIES=Num)

    !Se escribe la configuracion en el fichero binario
    WRITE(UConf) Conf
    WRITE(UConf) 1
    WRITE(UConf) (/0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/)
    WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
    WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
    k=1
    DO j=1,MoleculasSoluto
      WRITE(UConf) 0.0D0
      WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) Pos(k:k+SIZE(Soluto,1)-1,:)
      k=k+SIZE(Soluto,1)
    END DO
    DO j=1,MoleculasDisolvente
      WRITE(UConf) 0.0D0
      WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) Pos(k:k+SIZE(Disolvente,1)-1,:)
      k=k+SIZE(Disolvente,1)
    END DO
    DO j=1,MoleculasDisolvente2
      WRITE(UConf) 0.0D0
      WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) Pos(k:k+SIZE(Disolvente2,1)-1,:)
      k=k+SIZE(Disolvente2,1)
    END DO
  END DO

  CALL CerrarDCD(Tray)

  DEALLOCATE(SolConf,Pos)

END SUBROUTINE LeerConfigsGenerico

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE EntradaGenericoMM(Sal)
  USE Sistema
  USE Unidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Sal

  INTEGER :: i,j,Num

  Num=0
  Num=Num+MoleculasDisolvente*SIZE(Disolvente,1)
  Num=Num+MoleculasDisolvente2*SIZE(Disolvente2,1)

  WRITE(Sal,'(A)') 'Solute'
  WRITE(Sal,100) TRIM(NombreSoluto)
  WRITE(Sal,101) MoleculasSoluto
  WRITE(Sal,101) SIZE(Soluto,1)
  DO i=1,SIZE(Soluto,1)
    WRITE(Sal,102) Soluto(i)%id,Soluto(i)%nom,Soluto(i)%z, &
                   Soluto(i)%m/AmuAtomica,Soluto(i)%pos(:)/AngstromAtomica, &
                   Soluto(i)%q
  END DO

  WRITE(Sal,*)
  WRITE(Sal,'(A)') 'Solvent'
  WRITE(Sal,101) Num
  DO i=1,MoleculasDisolvente
    DO j=1,SIZE(Disolvente,1)
      WRITE(Sal,103) Disolvente(j)%id,Disolvente(j)%nom,Disolvente(j)%q
    END DO
  END DO
  DO i=1,MoleculasDisolvente2
    DO j=1,SIZE(Disolvente2,1)
      WRITE(Sal,103) Disolvente2(j)%id,Disolvente2(j)%nom,Disolvente2(j)%q
    END DO
  END DO

  WRITE(Sal,*)
  WRITE(Sal,'(A)') 'Non-Bonded'
  SELECT CASE (TipoPotencial)
   CASE (1)
    WRITE(Sal,100) 'lennard-jones'
   CASE (2)
    WRITE(Sal,100) 'generic'
  END SELECT
  Num=0
  DO i=1,SIZE(QInter,1)
    DO j=i,SIZE(QInter,2)
      IF (QInter(i,j)) Num=Num+1
    END DO
  END DO
  WRITE(Sal,101) Num
  DO i=1,SIZE(QInter,1)
    DO j=i,SIZE(QInter,2)
      IF (QInter(i,j)) WRITE(Sal,104) i,j,InterAtom(i,j,:)
    END DO
  END DO

100 FORMAT(2X,A)
101 FORMAT(I10)
102 FORMAT(2X,I4,1X,A16,1X,I3,1X,F7.4,4(1X,F12.6))
103 FORMAT(2X,I4,1X,A16,1X,F12.6)
104 FORMAT(2X,I4,1X,I4,7(1X,F14.10))

END SUBROUTINE EntradaGenericoMM

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
SUBROUTINE LeerConfigsGenerico_Borrar
  USE Configuraciones
  USE Parametros
  USE Unidades
  USE Utilidades
  USE UtilidadesFis
  USE Sistema
  USE TipoAtomo
  IMPLICIT NONE

  TYPE(Atomo), DIMENSION(:), ALLOCATABLE :: SolConf
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Pos
  DOUBLE PRECISION, DIMENSION(7) :: Mover
  DOUBLE PRECISION :: Paso,MaxDist
  CHARACTER(LEN=16) :: Nombre
  INTEGER :: UC,Conf,ConfTotal,NumAt,i,j,k

  NumAt=SIZE(Soluto,1)*MoleculasSoluto
  NumAt=NumAt+SIZE(Disolvente,1)*MoleculasDisolvente
  NumAt=NumAt+SIZE(Disolvente2,1)*MoleculasDisolvente2
  ALLOCATE(Pos(NumAt,3))

  ALLOCATE(SolConf(SIZE(Soluto,1)))
  SolConf(:)=Soluto(:)

  !Se abre el fichero de donde se leen las configuraciones
  UC=NuevaUnidad()
  OPEN(UC,FILE='c.tmp',STATUS='OLD')

  !Se abre el fichero binario donde se escribe
  CALL AbrirUConf()

  !Se escriben todas las configuraciones que se piden
  READ(UC,*) ConfTotal
  Paso=DBLE(ConfTotal)/DBLE(NumConfig)
  Conf=1
  DO i=1,NumConfig
    !Se pasan configuraciones hasta llegar a la deseada
    DO WHILE (Conf < NINT(i*Paso))
      READ(UC,*) NumAt
      READ(UC,*)
      DO j=1,NumAt
        READ(UC,*)
      END DO
      Conf=Conf+1
    END DO

    !Se leen los datos
    READ(UC,*) NumAt
    IF (NumAt /= SIZE(Pos,1)) CALL Mensaje('LeerConfigsGenerico',39,.TRUE.)
    READ(UC,*) MaxDist
    MaxDist=MaxDist*AngstromAtomica
    DO j=1,NumAt
      READ(UC,*) Nombre,Pos(j,1:3)
    END DO
    Pos(:,:)=Pos(:,:)*AngstromAtomica

    !Se gira y centra el soluto
    DO j=1,SIZE(Soluto,1)
      SolConf(j)%pos(:)=Pos(j,:)
    END DO
    CALL SuperponerMoleculas(Soluto,SolConf,Trans=Mover)
    Pos(:,:)=RotarCuaternion(Pos(:,:),Mover(1:4))
    Pos(:,:)=Pos(:,:)+SPREAD(Mover(5:7),DIM=1,NCOPIES=SIZE(Pos,1))

    !Se escribe la configuracion en el fichero binario
    WRITE(UConf) Conf
    WRITE(UConf) 1
    WRITE(UConf) (/MaxDist,0.0D0,0.0D0,0.0D0,MaxDist,0.0D0,0.0D0,0.0D0,MaxDist/)
    WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
    WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
    k=1
    DO j=1,MoleculasSoluto
      WRITE(UConf) 0.0D0
      WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) Pos(k:k+SIZE(Soluto,1)-1,:)
      k=k+SIZE(Soluto,1)
    END DO
    DO j=1,MoleculasDisolvente
      WRITE(UConf) 0.0D0
      WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) Pos(k:k+SIZE(Disolvente,1)-1,:)
      k=k+SIZE(Disolvente,1)
    END DO
    DO j=1,MoleculasDisolvente2
      WRITE(UConf) 0.0D0
      WRITE(UConf) (/0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) (/1.0D0,0.0D0,0.0D0,0.0D0/)
      WRITE(UConf) Pos(k:k+SIZE(Disolvente2,1)-1,:)
      k=k+SIZE(Disolvente2,1)
    END DO
    Conf=Conf+1

  END DO

  !Se borran los ficheros
  CLOSE(UC,STATUS='DELETE')

  DEALLOCATE(Pos,SolConf)

END SUBROUTINE LeerConfigsGenerico_Borrar

END MODULE GenericoMM
