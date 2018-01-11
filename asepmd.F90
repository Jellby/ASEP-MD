!##############################################################################
!# Copyright 2011,2012,2018 Ignacio Fdez. Galván, M. Luz Sánchez,             #
!#                     Aurora Muñoz Losa, M. Elena Martín, Manuel A. Aguilar  #
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

PROGRAM ASEPMD
  USE Parametros
  USE Entrada
  USE Moldy
  USE Malla
  USE Cavidad
  USE Configuraciones
  USE Sistema
  USE Ejecutar
  USE Optimizacion
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Potencial,ASEP
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CargasAjustadas
  DOUBLE PRECISION, DIMENSION(4) :: ErrAjuste
  DOUBLE PRECISION, DIMENSION(4) :: IntMD
  DOUBLE PRECISION :: EnergiaAnt,CargaTotal,DifCargas
  INTEGER :: U,i,Iter,IterIni,Num,NumCargas,Error

!>>> Iniciar datos

  EnergiaAnt=0.0D0
  IntMD(:)=0.0D0

!>>> Leer la entrada

  CALL EscribirTitulo()

  CALL LeerEntrada(5)

  CALL EscribirValores()

  SELECT CASE (ProgramaMM)
   CASE (0) !Genérico
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaGenerico(U)
    CLOSE(U)

    !Si no se usa Moldy, sólo se centra la molécula
    CALL OrientarMolecula(Soluto(:),-1)
    CALL OrientarMolecula(Disolvente(:),-1)
    CALL OrientarMolecula(Disolvente2(:),-1)

   CASE (1) !Moldy
    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(EntradaMM),ACTION='READ',STATUS='OLD')
    CALL LeerControlMoldy(U)
    CLOSE(U)

    U=NuevaUnidad()
    OPEN(U,FILE=TRIM(MoldyInput),ACTION='READ',STATUS='OLD')
    CALL LeerSistemaMoldy(U)
    CLOSE(U)

    CALL OrientarMolecula(Soluto(:),0)
    CALL OrientarMolecula(Disolvente(:),0)
    CALL OrientarMolecula(Disolvente2(:),0)
  END SELECT

!>>> Cálculo cuántico inicial

  Iter=0
  WRITE(Extension,*) Iter
  Extension='.cic.'//ADJUSTL(Extension)

  IF (ALLOCATED(MolQM)) DEALLOCATE(MolQM)
  ALLOCATE(MolQM(SIZE(Soluto,1)))
  MolQM=Soluto

  IF (InicioVacio) THEN
    CALL EjecutarQM(0)

    Soluto=MolQM

    CALL EscribirASEPMD()
   ELSE
    ALLOCATE(DisolvQM(0,5),PotQM(0,4))
  END IF

!>>> Empieza el ciclo ASEP/MD

  IF (Polarizable) THEN
    IterIni = 1
   ELSE
    IterIni = Inicio
  END IF

  DO Iter=IterIni,MaxIter

    WRITE(Extension,*) Iter
    Extension='.cic.'//ADJUSTL(Extension)

    EnergiaAnt=EnergiaQM+EnergiaVdW

!>>> Generar la "cavidad" para el soluto y las cargas explícitas del disolvente

    CALL GenerarCavidad

!>>> Lanzar la dinámica

    ! Sólo se hace la dinámica si no es un cálculo polarizable
    ! o en la primera iteración
    IF ((Iter == IterIni) .OR. (.NOT. Polarizable)) THEN
      IF (Polarizable) THEN
        WRITE(Extension,*) Inicio
       ELSE
        WRITE(Extension,*) Iter
      END IF
      Extension='.cic.'//ADJUSTL(Extension)
      CALL EjecutarMM(.TRUE.)

      CALL Promedios(0,Soluto) !Abre el fichero UConf
      IntMD(1)=EnergiaEMM
      IntMD(2)=EnergiaVdW
      IntMD(3)=EnergiaQDip(1)
      IntMD(4)=EnergiaQDip(2)
    END IF

    IF (Polarizable) THEN
      WRITE(Extension,*) Iter
      Extension='.pol.'//ADJUSTL(Extension)
    END IF

!>>> Calcular el potencial promedio y cargas explícitas

    CALL MallaMolecula(Soluto,MallaSoluto)

    CALL PotencialYCargas

!>>> Ajustar las cargas externas

    CALL AjusteCargasExternas

!>>> Cálculo cuántico con disolvente

    DifCargas=0.0D0

    IF (MaxIterOpt > 0) THEN

      U=NuevaUnidad()
      OPEN(U,FILE=TRIM(SalidaOpt)//TRIM(Extension),ACTION='WRITE',STATUS='REPLACE')
      CALL OptimizarGeometria(Soluto,U)
      CLOSE(U)

      IF (ProgramaMM == 1) THEN
        CALL OrientarMolecula(Soluto(:),0)
       ELSE
        !Si no se usa Moldy, sólo se centra la molécula
        CALL OrientarMolecula(Soluto(:),-1)
      END IF

     ELSE

      IF (ALLOCATED(MolQM)) DEALLOCATE(MolQM)
      ALLOCATE(MolQM(SIZE(Soluto,1)))
      MolQM=Soluto
      CALL EjecutarQM(0)
      DO i=1,SIZE(Soluto)
        DifCargas=MAX(DifCargas,ABS(Soluto(i)%q-MolQM(i)%q))
      END DO
      Soluto=MolQM

    END IF

    IF (.NOT. Polarizable) CALL CerrarUConf()

    CALL EscribirASEPMD()

    IF (Polarizable .AND. (DifCargas < ConvCargasPol)) EXIT

    IF (Iter == MaxIter) CALL Mensaje('ASEPMD',44,.FALSE.)

  END DO

!>>> Escribe el sistema final (útil para las cargas del soluto)
  IF (Polarizable) THEN
    Extension='.pol.final'
  ELSE
    Extension='.final'
  END IF
  CALL EjecutarMM(.FALSE.)

CONTAINS

! Escribe la cabecera en la salida
SUBROUTINE EscribirTitulo

#ifndef LANG
#define LANG Espanol
#endif

  USE LANG
  INTEGER :: i

  DO i=1,SIZE(Titulo)
    WRITE(6,101) TRIM(Titulo(i))
  END DO
  WRITE(6,*)

101 FORMAT(A)

END SUBROUTINE EscribirTitulo

! Escribe los datos correspondientes de un ciclo ASEP/MD
SUBROUTINE EscribirASEPMD

#ifndef LANG
#define LANG Espanol
#endif

  USE LANG

  WRITE(6,*)
  WRITE(6,101) TRIM(Textos(8)),Iter
  WRITE(6,20)
  IF (Iter > 0) THEN
    IF ((Iter == IterIni) .OR. (.NOT. Polarizable)) THEN
      WRITE(6,100) TRIM(Textos(24))
      WRITE(6,10)
      WRITE(6,103) TRIM(Textos(26)),IntMD(1),'Eh'
      WRITE(6,104) '',IntMD(1)/KcalmolAtomica,'kcal/mol'
      WRITE(6,103) TRIM(Textos(27)),IntMD(2),'Eh'
      WRITE(6,104) '',IntMD(2)/KcalmolAtomica,'kcal/mol'
      WRITE(6,103) TRIM(Textos(78)),IntMD(3),'Eh'
      WRITE(6,104) '',IntMD(3)/KcalmolAtomica,'kcal/mol'
      WRITE(6,103) TRIM(Textos(79)),IntMD(4),'Eh'
      WRITE(6,104) '',IntMD(4)/KcalmolAtomica,'kcal/mol'
      WRITE(6,*)
    END IF
    WRITE(6,100) TRIM(Textos(33))
    WRITE(6,10)
    WRITE(6,105) TRIM(Textos(34)),NumCargas
    WRITE(6,105) TRIM(Textos(35)),SIZE(CargasDisolvente,1)
    WRITE(6,105) TRIM(Textos(40)),SIZE(ASEP,1)
    WRITE(6,106) TRIM(Textos(36)),ErrAjuste(1)
    WRITE(6,106) TRIM(Textos(37)),ErrAjuste(2),'Eh/e'
    WRITE(6,106) TRIM(Textos(38)),ErrAjuste(3)
    WRITE(6,106) TRIM(Textos(39)),ErrAjuste(4),'Eh/e'
  END IF
  WRITE(6,*)
  WRITE(6,100) TRIM(Textos(25))
  WRITE(6,10)
  WRITE(6,102) TRIM(Textos(23)),EnergiaQM,'Eh'
  IF (Iter > 0) THEN
    WRITE(6,102) TRIM(Textos(29)),EnergiaQM+EnergiaVdW,'Eh'
    WRITE(6,102) TRIM(Textos(30)),EnergiaQM+EnergiaVdW-EnergiaAnt,'Eh'
  END IF
  WRITE(6,103) TRIM(Textos(28)),SQRT(SUM(DipoloQM(:)*DipoloQM(:))),'e*a0'
  WRITE(6,104) '',SQRT(SUM(DipoloQM(:)*DipoloQM(:)))/DebyeAtomica,'D'
  IF (Iter > 0) THEN
    WRITE(6,103) TRIM(Textos(31)),EnergiaEQM,'Eh'
    WRITE(6,104) '',EnergiaEQM/KcalmolAtomica,'kcal/mol'
    WRITE(6,103) TRIM(Textos(26)),EnergiaEMM,'Eh'
    WRITE(6,104) '',EnergiaEMM/KcalmolAtomica,'kcal/mol'
    WRITE(6,103) TRIM(Textos(27)),EnergiaVdW,'Eh'
    WRITE(6,104) '',EnergiaVdW/KcalmolAtomica,'kcal/mol'
    WRITE(6,103) TRIM(Textos(78)),EnergiaQDip(1),'Eh'
    WRITE(6,104) '',EnergiaQDip(1)/KcalmolAtomica,'kcal/mol'
    WRITE(6,103) TRIM(Textos(79)),EnergiaQDip(2),'Eh'
    WRITE(6,104) '',EnergiaQDip(2)/KcalmolAtomica,'kcal/mol'
    WRITE(6,102) TRIM(Textos(32)),EnergiaQM-EnergiaEQM,'Eh'
  END IF
  WRITE(6,20)

 10 FORMAT(80('-'))
 20 FORMAT(80('='))
100 FORMAT(T10,A)
101 FORMAT(A,1X,I3)
102 FORMAT(A,T54,F18.10,1X,A)
103 FORMAT(A,T50,'/',T54,F18.10,1X,A)
104 FORMAT(A,T50,'\',T54,F18.10,1X,A)
105 FORMAT(T15,A,T64,I8,1X,A)
106 FORMAT(T15,A,T54,F18.10,1X,A)

END SUBROUTINE

! Construye la "cavidad" y calcula las posiciones de las cargas ajustadas en la
! superficie
SUBROUTINE GenerarCavidad
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Esf

  ! Según el tipo de cavidad, se define el número y tamaño adecuado de
  ! las esferas que forman la cavidad
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

  ! Una superficie un poco más grande donde se sitúan las cargas ajustadas
  ! (en dos capas)
  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4)+RadioDisolvente,RadioDisolvente, &
                        Subdivisiones)

  IF (ALLOCATED(CargasAjustadas)) DEALLOCATE(CargasAjustadas)
  ALLOCATE(CargasAjustadas(SIZE(CavVert,1)+SIZE(CavCent,1),4))

  ! Las caras están más cerca
  CargasAjustadas(1:SIZE(CavCent,1),1:3)=CavCent(:,1:3)

  ! Los vértices, que son menos, están más lejos
  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4)+2.0D0*RadioDisolvente,RadioDisolvente, &
                        Subdivisiones)

  CargasAjustadas(SIZE(CavCent,1)+1:,1:3)=CavVert(:,1:3)

  ! Otra final, para separarar las moléculas que están dentro, y es la
  ! que permanece en el resto del cálculo
  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4),RadioDisolvente)

  DEALLOCATE(Esf)

END SUBROUTINE GenerarCavidad

! Calcula el potencial sobre el soluto y las cargas del disolvente
SUBROUTINE PotencialYCargas
  DOUBLE PRECISION, DIMENSION(3) :: Pos
  INTEGER :: Utmp,i,j

  IF (ALLOCATED(Potencial)) DEALLOCATE(Potencial)
  IF (ALLOCATED(ASEP)) DEALLOCATE(ASEP)
  ALLOCATE(Potencial(SIZE(MallaSoluto,1)),ASEP(SIZE(MallaSoluto,1)))

  Utmp=NuevaUnidad()
  OPEN(Utmp,STATUS='SCRATCH',ACTION='READWRITE',FORM='UNFORMATTED')

  U=NuevaUnidad()
  ASEP(:)=0.0D0
  Num=0
  NumCargas=0
  CALL AbrirUConf()
  CALL AbrirUDip()
  DO
    CALL LeerConfig(Error)
    IF (Error /= 0) EXIT
    CALL PotencialMM(MallaSoluto,Potencial)
    ASEP(:)=ASEP(:)+Potencial(:)
    CALL SeleccionarMols(Utmp,NumCargas)
    Num=Num+1
  END DO
  CALL CerrarUDip()

  CALL ReducirCargas(Utmp,NumCargas,CargasDisolvente)
  CargasDisolvente(:,4)=CargasDisolvente(:,4)/DBLE(Num)

  ! Se hacen cero las cargas explícitas que están muy cerca de las ajustadas
  DO i=1,SIZE(CargasDisolvente,1)
    Pos=CargasDisolvente(i,1:3)
    DO j=1,SIZE(CargasAjustadas,1)
      IF (Distancia(Pos(:),CargasAjustadas(j,1:3)) < DistCargas) THEN
        CargasDisolvente(i,4)=0.0D0
        EXIT
      END IF
    END DO
  END DO

  CALL PotencialCargas(CargasDisolvente(:,:),MallaSoluto(:,:),Potencial(:))

  ASEP(:)=ASEP(:)/DBLE(Num)

  CLOSE(Utmp)

END SUBROUTINE PotencialYCargas

! Calcula el valor de las cargas del disolvente por ajuste al potencial sobre
! el soluto
SUBROUTINE AjusteCargasExternas
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Pot

  ALLOCATE(Pot(SIZE(ASEP),4))

  Pot(:,1:3)=MallaSoluto(:,1:3)
  Pot(:,4)=ASEP(:)-Potencial(:)
  CargaTotal=SUM(CargasDisolvente(:,4))+SUM(MolQM(:)%q)
  CALL AjusteCargas(CargasAjustadas(:,1:3),Pot(:,:),-CargaTotal, &
                    CargasAjustadas(:,4),ErrAjuste(1:3),1.0D-6)

  Pot(:,4)=Potencial(:)
  CALL PotencialCargas(CargasAjustadas(:,:),MallaSoluto(:,:),Potencial(:))
  Potencial(:)=Potencial(:)+Pot(:,4)

  ErrAjuste(4)=MAXVAL(ABS(Potencial(:)-ASEP(:)))

  DEALLOCATE(Pot)

  IF (ALLOCATED(DisolvQM)) DEALLOCATE(DisolvQM)
  Num=SIZE(CargasAjustadas,1)+SIZE(CargasDisolvente,1)
  ALLOCATE(DisolvQM(Num,5))

  DisolvQM(1:SIZE(CargasAjustadas,1),1:4)=CargasAjustadas(:,:)
  DisolvQM(SIZE(CargasAjustadas,1)+1:,1:4)=CargasDisolvente(:,:)
  DisolvQM(:,5)=0.0D0

END SUBROUTINE AjusteCargasExternas

END PROGRAM ASEPMD
