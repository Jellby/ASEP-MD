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

MODULE Optimizacion

#ifndef LANG
#define LANG Espanol
#endif

#ifdef __NAG__
USE F90_UNIX_IO
#endif

USE Ejecutar
USE DatosQMMM
USE LANG

!-------------------------------------------------------------------------------
! Inc:       Incremento de coordenadas (dirección) en coordenadas de trabajo
! CoordPrev: Coordenadas cartesianas del paso anterior (origen)
!-------------------------------------------------------------------------------
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Inc,CoordPrev

CONTAINS
!OptimizarGeometria(Mol,Fich)
!EnergPunto(Factor)
!MinimoInterpolado(Func0,Func1,Der0,Der1,Minim,Error)
!BuscarMinimo(EnergPunto,A,B,C,Preci)
!IncrementoRFO(Hess,Grad,Modo,Newt,Paso,Error)
!ActualizarHessiana(Inc,IncGrad,Hess,Tipo)
!HessianaInicial(Tipo)

!-------------------------------------------------------------------------------
! Subrutina para optimizar la geometría de una molécula
!-------------------------------------------------------------------------------
! Mol:               Molécula que se optimiza
! Fich:              Número del fichero donde se escribe la salida
! CoordCart:         Coordenadas cartesianas
! GradCart:          Gradiente en coordenadas cartesianas
! IncCart:           Variación de la geometría en coordenadas cartesianas
! IncPrev:           Variación de la iteración anterior
! GradPrev:          Gradiente de la iteración anterior
! CargasPrev:        Cargas atómicas de la iteración anterior
! GradTot,HessTot:   Gradiente y hessiana totales (con contribucion vdW sumada)
! DipoloPrev:        Momento dipolar de la iteración anterior
! Energ:             Energía
! EnergPrev:         Energía de la iteración anterior
! Paso:              Factor por el que se multiplica la variación de geometría
! ModCoord:          Variable para ver si se modifican las coordenadas internas
! Ext:               Extensión de los ficheros que se crean
! Num:               Número de átomos en la molécula
! Iter:              Número de iteración de optimización
! Conv1,Conv2,Conv3: Variables para comprobar la convergencia
! Conv:              Variable que indica si ha convergido la optimización
! Modo:              Modo vibracional que se maximiza (estado de transición)
! Calc:              Variable que indica si es necesario hacer un nuevo cálculo
! Hess:              Variable que indica si se calcula la hessiana exacta
! UMod:              Unidad con las modificaciones de coordenadas
! i,j:               Contadores
!-------------------------------------------------------------------------------
SUBROUTINE OptimizarGeometria(Mol,Fich)
  USE TipoAtomo
  USE Parametros
  USE Coordenadas
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(INOUT) :: Mol
  INTEGER, INTENT(IN) :: Fich

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CoordCart,GradCart,IncCart, &
                                                 IncPrev,GradPrev,CargasPrev, &
                                                 GradTot
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: HessTot
  DOUBLE PRECISION, DIMENSION(3) :: DipoloPrev
  DOUBLE PRECISION :: Energ,EnergPrev,Paso,Conv1,Conv2,Conv3
  LOGICAL :: ModCoord
  CHARACTER(LEN=LEN(Extension)) :: Ext
  INTEGER :: Num,Iter,Conv,Modo,Calc,Hess,UMod,i,j

  !Se inician los datos para los cálculos QM
  Num=SIZE(Mol,1)
  IF (ALLOCATED(MolQM)) DEALLOCATE(MolQM)
  ALLOCATE(MolQM(Num),CoordCart(3*Num),GradCart(3*Num),IncCart(3*Num), &
           CoordPrev(3*Num),CargasPrev(Num),GradTot(3*Num),HessTot(3*Num,3*Num))
  MolQM=Mol
  Ext=Extension

  !Se inician las coordenadas cartesianas y las de trabajo
  CoordCart(:)=Cartesianas(Mol)
  TipoCoord=TipoCoordenadas
  CALL DefinirCoordenadas(Mol)
  INQUIRE(FILE=TRIM(AltCoordenadas),EXIST=ModCoord)
  IF (ModCoord) THEN
    UMod=NuevaUnidad()
    OPEN(UMod,FILE=TRIM(AltCoordenadas),STATUS='OLD',ACTION='READ')
    CALL ModificarCoordenadas(UMod)
    CLOSE(UMod)
  END IF
  ALLOCATE(Inc(SIZE(DefCoord,1)),IncPrev(SIZE(DefCoord,1)), &
           GradPrev(SIZE(DefCoord,1)))

  !Se hace el cálculo inicial
  Iter=0
  WRITE(Extension,*) Iter
  Extension=TRIM(Ext)//'.opt.'//ADJUSTL(Extension)

  CALL LeerContinuacion

  IF (MetodoOptim < 2) THEN
    !Si se emplea un método de primer orden, sólo es necesario el gradiente
    CALL ConvertirCoordenadas(CoordCart(:),1)
   ELSE
    !Si se emplea un método de segundo orden, se usará la hessiana
    CALL ConvertirCoordenadas(CoordCart(:),2)
  END IF

  IF (Iter == 0) THEN
    Hess=0
    IF (MetodoOptim < 2) THEN
      CALL EjecutarQM(1)
      GradTot=GradQM(:)+GradVdW(:)
      CALL ConvertirGradHess(1,GradTot(:))
      CALL Proyectar(Grad=Gradiente(:))
     ELSE
      !Se calcula o no la hessiana inicial según el valor de CalcHessiana
      IF ((CalcHessiana > 0) .OR. (HessInicial == 0)) THEN
        Hess=1
        CALL EjecutarQM(2)
        GradTot(:)=GradQM(:)+GradVdW(:)
        HessTot(:,:)=HessQM(:,:)+HessVdW(:,:)
        CALL ConvertirGradHess(1,GradTot(:),HessTot(:,:))
       ELSE
        CALL EjecutarQM(1)
        GradTot=GradQM(:)+GradVdW(:)
        CALL ConvertirGradHess(1,GradTot(:))
        CALL HessianaInicial(HessInicial)
      END IF
      CALL Proyectar(Grad=Gradiente(:),Hess=Hessiana(:,:))
    END IF

    !Se convierte el gradiente a cartesianas
    CALL ConvertirGradHess(-1,GradCart(:))
    Energ=0.0D0
    Inc(:)=0.0D0
  END IF

  CALL EscribirGeometria

  Conv=0
  DO
    Iter=Iter+1

    !Se actualizan energías y coordenadas
    EnergPrev=Energ
    Energ=EnergiaQM+EnergiaVdW
    CoordPrev(:)=CoordCart(:)
    IncPrev(:)=Inc(:)

    !Se comprueba la convergencia (gradiente y energía)
    Conv1=MAXVAL(ABS(GradCart(:)))
    Conv3=Energ-EnergPrev
    WRITE(Fich,101) TRIM(Textos(1)),Conv1,'Eh/a0'
    WRITE(Fich,101) TRIM(Textos(2)),Conv3,'Eh'
    CALL FLUSH(Fich)
    IF ((Conv1 <= ConvGradOpt) .AND. &
        ((ABS(Conv3) <= ConvEnerOpt) .OR. (Iter == 1))) THEN
      Conv=1
      EXIT
    END IF

    !Se calcula el paso de optimización con el método escogido
    Modo=0
    IF (EstadoTransicion) Modo=1
    SELECT CASE (MetodoOptim)
     CASE (0) !Pendiente máxima
      Inc(:)=-Gradiente(:)
     CASE (1) !Gradiente conjugado (Fletcher-Reeves)
      IF (Iter == 1) GradPrev(:)=Gradiente(:)
      Inc(:)=-Gradiente(:) + IncPrev(:) * &
             DOT_PRODUCT(Gradiente(:),Gradiente(:))/ &
             DOT_PRODUCT(GradPrev(:),GradPrev(:))
     CASE (2) !Newton-Raphson
      CALL IncrementoRFO(Hessiana(:,:),Gradiente(:),Modo,1,Inc(:),Num)
     CASE (3) !RFO
      CALL IncrementoRFO(Hessiana(:,:),Gradiente(:),Modo,0,Inc(:),Num)
    END SELECT

    !Se actualiza el gradiente del paso anterior
    GradPrev(:)=Gradiente(:)

    !Se va a hacer un nuevo cálculo
    WRITE(Extension,*) Iter
    Extension=TRIM(Ext)//'.opt.'//ADJUSTL(Extension)

    Calc=1
    SELECT CASE (BusquedaLineal)
     CASE (0)
      !No hay búsqueda lineal, sólo se recorta el paso
      Paso=MIN(1.0D0,MaxPasoOpt/Norma(Inc))
      Inc(:)=Paso*Inc(:)
     CASE (1)
      !Búsqueda interpolada
      !Se guardan los datos previos
      DipoloPrev(:)=DipoloQM(:)
      CargasPrev(:)=MolQM(:)%q
      GradCart(:)=GradQM(:)+GradVdW(:)
      !Se hace un cálculo con el paso completo
      Paso=MIN(1.0D0,MaxPasoOpt/Norma(Inc))
      Inc(:)=Paso*Inc(:)
      CALL ConvertirIncremento(Inc(:),CoordPrev(:),IncCart(:))
      CoordCart(:)=CoordPrev(:)+IncCart(:)
      CALL CambiarCartesianas(MolQM,CoordCart(:))
      CALL ConvertirCoordenadas(CoordCart(:),1)
      CALL EjecutarQM(1)
      GradTot(:)=GradQM(:)+GradVdW(:)
      CALL ConvertirGradHess(1,GradTot(:))
      CALL Proyectar(Grad=Gradiente(:))
      !Se halla el mínimo del polinomio de interpolación
      CALL MinimoInterpolado(Energ,EnergiaQM,DOT_PRODUCT(GradPrev(:),Inc(:)), &
                             DOT_PRODUCT(Gradiente(:),Inc(:)),Paso,Num)
      !Calc=0 significa que no es necesario hacer un nuevo cálculo
      !Si el paso es muy cercano a 1.0, se toma 1.0
      IF ((Paso > 0.9D0) .AND. (Paso < 1.1D0)) THEN
        Calc=0
        Paso=1.0D0
      END IF
      !Si el paso es muy pequeño o muy grande no es fiable y se toma 1.0
      IF (MetodoOptim > 1) THEN
        IF ((Paso < 0.1D0) .OR. (Paso > 1.1D0)) THEN
          Calc=0
          Paso=1.0D0
        END IF
      END IF
      Inc(:)=Paso*Inc(:)
      CALL ConvertirCoordenadas(CoordPrev(:),1)
     CASE (2)
      !Búsqueda completa
      Paso=MIN(1.0D0,MaxPasoOpt/Norma(Inc))
      CALL BuscarMinimo(EnergPunto,0.0D0,0.4D0,Paso,ConvEnerOpt)
      Inc=Paso*Inc
      CALL ConvertirCoordenadas(CoordPrev(:),1)
    END SELECT

    !Se convierte el nuevo incremento a cartesianas si es necesario
    IF (Calc == 1) THEN
      CALL ConvertirIncremento(Inc(:),CoordPrev(:),IncCart(:))
    END IF

    !Se comprueba la convergencia (gradiente e incremento)
    Conv2=MAXVAL(ABS(IncCart))
    WRITE(Fich,101) TRIM(Textos(3)),Conv2,'a0'
    CALL FLUSH(Fich)
    IF ((Conv1 <= ConvGradOpt) .AND. (Conv2 <= ConvPasoOpt)) THEN
      !Si ha habido busqueda lineal, se recuperan los datos iniciales
      IF (BusquedaLineal > 0) THEN
        CALL CambiarCartesianas(MolQM,CoordPrev(:))
        CALL ConvertirCoordenadas(CoordPrev(:),0)
        MolQM(:)%q=CargasPrev(:)
        EnergiaQM=Energ-EnergiaVdW
        GradQM(:)=GradCart(:)+GradVdW(:)
        Gradiente(:)=GradPrev(:)
        DipoloQM(:)=DipoloPrev(:)
      END IF
      Conv=1
      EXIT
    END IF
    IF (Iter > MaxIterOpt) EXIT

    !Se modifican las coordenadas
    CoordCart(:)=CoordPrev(:)+IncCart(:)
    CALL CambiarCartesianas(MolQM,CoordCart(:))

    !Si Hess=1, hay que calcular la hessiana
    !Si Calc=0, no hay que hacer un nuevo cálculo
    Hess=0
    IF (CalcHessiana > 0) THEN
      IF (MOD(Iter,CalcHessiana) == 0) THEN
        Hess=1
      END IF
    END IF
    IF (MetodoOptim < 2) THEN
      !Métodos de primer orden
      IF (Calc == 1) THEN
        CALL ConvertirCoordenadas(CoordCart(:),1)
        CALL EjecutarQM(1)
        GradTot(:)=GradQM(:)+GradVdW(:)
        CALL ConvertirGradHess(1,GradTot(:))
        CALL Proyectar(Grad=Gradiente(:))
      END IF
     ELSE
      !Métodos de segundo orden
      CALL ConvertirCoordenadas(CoordCart(:),2)
      IF (Hess == 1) THEN
        !Se calcula la hessiana
        CALL EjecutarQM(2)
        GradTot(:)=GradQM(:)+GradVdW(:)
        HessTot(:,:)=HessQM(:,:)+HessVdW(:,:)
        CALL ConvertirGradHess(1,GradTot(:),HessTot(:,:))
        CALL Proyectar(Grad=Gradiente(:),Hess=Hessiana(:,:))
       ELSE
        !Se actualiza la hessiana
        IF (Calc == 1) THEN
          CALL EjecutarQM(1)
          GradTot(:)=GradQM(:)+GradVdW(:)
          CALL ConvertirGradHess(1,GradTot(:))
          CALL Proyectar(Grad=Gradiente(:))
        END IF
        CALL ActualizarHessiana(Inc(:),Gradiente(:)-GradPrev(:),Hessiana(:,:), &
                                Actualizacion)
        CALL Proyectar(Hess=Hessiana(:,:))
      END IF
    END IF

    !Se convierte el gradiente a cartesianas
    CALL ConvertirGradHess(-1,GradCart(:))
    CALL EscribirGeometria
  END DO

  !Para terminar, se escribe si ha convergido la optimización o no
  WRITE(Fich,10)
  WRITE(Fich,*)
  IF (Conv == 1) THEN
    WRITE(Fich,100) TRIM(Textos(4))
   ELSE
    WRITE(Fich,100) TRIM(Textos(5))
  END IF

  !Escribe la hessiana final si se ha usado un método de segundo orden
  IF (MetodoOptim > 1) THEN
    WRITE(Fich,*)
    IF (Hess == 1) THEN
      WRITE(Fich,100) TRIM(Textos(19))
      WRITE(Fich,102) ((HessTot(i,j),j=1,i),i=1,SIZE(HessTot,1))
     ELSE
      WRITE(Fich,100) TRIM(Textos(20))
      CALL ConvertirGradHess(-1,Hess=HessTot(:,:))
      WRITE(Fich,102) ((HessTot(i,j),j=1,i),i=1,SIZE(HessTot,1))
    END IF
  END IF

  Mol=MolQM
  Extension=Ext

  DEALLOCATE(CoordCart,GradCart,IncCart,Inc,CoordPrev,IncPrev,GradPrev, &
             CargasPrev,GradTot,HessTot)

 10 FORMAT(80('='))
100 FORMAT(A)
101 FORMAT(A,T37,F18.10,1X,A)
102 FORMAT(5(1X,ES15.8))

  CONTAINS

  !Subrutina que escribe la geometría y algunos datos durante la optimización
  !UCon: Unidad del fichero de continuación
  !i:    Contador
  SUBROUTINE EscribirGeometria
    USE Unidades
    IMPLICIT NONE

    INTEGER :: UCon,i

    !Escribe los datos necesarios para la continuación del cálculo
    UCon=NuevaUnidad()
    OPEN(UCon,FILE='optim.scratch',STATUS='REPLACE',ACTION='WRITE', &
         FORM='UNFORMATTED')
    WRITE(UCon) Iter
    WRITE(UCon) SIZE(CoordCart)
    WRITE(UCon) CoordCart(:)
    WRITE(UCon) GradCart(:)
    WRITE(UCon) SIZE(DefCoord,1)
    WRITE(UCon) DefCoord(:,:)
    WRITE(UCon) Geometria(:)
    WRITE(UCon) Inc(:)
    WRITE(UCon) Gradiente(:)
    WRITE(UCon) Hessiana(:,:)
    WRITE(UCon) Energ,EnergiaQM,EnergiaVdW,DipoloQM,Hess
    WRITE(UCon) SIZE(MolQM)
    WRITE(UCon) MolQM(:)%q
    WRITE(UCon) (MolQM(i)%pos(:),i=1,SIZE(MolQM))
    WRITE(UCon) GradQM(:)
    WRITE(UCon) GradVdW(:)
    WRITE(UCon) HessVdW(:,:)
    CLOSE(UCon)

    !Escribe la cabecera de cada iteración
    WRITE(Fich,20)
    IF (Iter == 0) THEN
      WRITE(Fich,106) TRIM(Textos(6)),SIZE(DefCoord,1)
      WRITE(Fich,*)
      WRITE(Fich,100) TRIM(Textos(7))
     ELSE
      WRITE(Fich,*)
      WRITE(Fich,106) TRIM(Textos(8)),Iter
    END IF

    !Escribe la geometría y el gradiente en coordenadas cartesianas
    WRITE(Fich,20)
    WRITE(Fich,100) TRIM(Textos(9))
    DO i=1,SIZE(MolQM,1)
      WRITE(Fich,101) MolQM(i)%nom,MolQM(i)%pos(:)/AngstromAtomica
    END DO
    WRITE(Fich,10)
    DO i=1,SIZE(MolQM,1)
      WRITE(Fich,101) MolQM(i)%nom,GradCart(3*i-2:3*i)*AngstromAtomica
    END DO
    WRITE(Fich,20)

    !Escribe la geometría y el gradiente en coordenadas de trabajo
    SELECT CASE (TipoCoord)
     CASE (1)
      WRITE(Fich,100) TRIM(Textos(10))
      DO i=1,SIZE(MolQM,1)
        WRITE(Fich,101) MolQM(i)%nom,Geometria(3*i-2:3*i)/AngstromAtomica
      END DO
      WRITE(Fich,10)
      DO i=1,SIZE(MolQM,1)
        WRITE(Fich,101) MolQM(i)%nom,Gradiente(3*i-2:3*i)*AngstromAtomica
      END DO
      WRITE(Fich,20)
     CASE (2)
      WRITE(Fich,100) TRIM(Textos(11))
      DO i=1,SIZE(DefCoord,1)
        IF (DefCoord(i,4) > 0) THEN
          WRITE(Fich,102) TRIM(Textos(12)),DefCoord(i,1:4),Geometria(i)/Grado, &
                          Gradiente(i)*Grado
        ELSE IF (DefCoord(i,3) > 0) THEN
          WRITE(Fich,103) TRIM(Textos(13)),DefCoord(i,1:3),Geometria(i)/Grado, &
                          Gradiente(i)*Grado
        ELSE IF (DefCoord(i,1) > 0) THEN
          WRITE(Fich,104) TRIM(Textos(14)),DefCoord(i,1:2), &
                          Geometria(i)/AngstromAtomica, &
                          Gradiente(i)*AngstromAtomica
        ELSE IF (DefCoord(i,3) /= 0) THEN
          WRITE(Fich,105) TRIM(Textos(21)),-DefCoord(i,1),Geometria(i)/Grado, &
                          Gradiente(i)*Grado
         ELSE
          WRITE(Fich,105) TRIM(Textos(21)),-DefCoord(i,1), &
                          Geometria(i)/AngstromAtomica, &
                          Gradiente(i)*AngstromAtomica
        END IF
      END DO
      WRITE(Fich,20)
     CASE DEFAULT
    END SELECT

    !Escribe energía y gradiente
    WRITE(Fich,107) TRIM(Textos(23)),EnergiaQM,'Eh'
    WRITE(Fich,107) TRIM(Textos(15)),EnergiaQM+EnergiaVdW,'Eh'
    WRITE(Fich,107) TRIM(Textos(16)),Norma(GradCart),'Eh/a0'
    IF ((MetodoOptim > 1) .AND. (Hess == 1)) WRITE(Fich,100) TRIM(Textos(17))
    IF (BusquedaLineal == 2) WRITE(Fich,100) TRIM(Textos(18))
    WRITE(Fich,10)

   10 FORMAT(80('-'))
   20 FORMAT(80('='))
  100 FORMAT(A)
  101 FORMAT(1X,A,3(1X,F16.10))
  102 FORMAT(1X,A,T13,4(1X,I3),T30,2(1X,F16.10))
  103 FORMAT(1X,A,T13,3(1X,I3),T30,2(1X,F16.10))
  104 FORMAT(1X,A,T13,2(1X,I3),T30,2(1X,F16.10))
  105 FORMAT(1X,A,T13,1(1X,I3),T30,2(1X,F16.10))
  106 FORMAT(A,1X,I3)
  107 FORMAT(A,T23,F18.10,1X,A)

  END SUBROUTINE EscribirGeometria

  !Subrutina que lee los datos intermedios para continuar una optimización
  !Cont: Variable que indica si se lee la continuación o no
  !N:    Tamaño de las matrices
  !UCon: Unidad del fichero de continuación
  SUBROUTINE LeerContinuacion
    IMPLICIT NONE

    LOGICAL :: Cont
    INTEGER :: N,UCon

    INQUIRE(FILE=TRIM(OptContinuacion),EXIST=Cont)

    IF (.NOT. Cont) RETURN

    UCon=NuevaUnidad()
    OPEN(UCon,FILE=TRIM(OptContinuacion),STATUS='OLD',ACTION='READ', &
         FORM='UNFORMATTED')
    READ(UCon) Iter

    READ(UCon) N
    IF ((N /= SIZE(CoordCart)) .OR. (N /= SIZE(GradCart))) &
       CALL Mensaje('OptimizarGeometria',26,.TRUE.)
    READ(UCon) CoordCart(:)
    READ(UCon) GradCart(:)
    IF (ALLOCATED(GradQM)) DEALLOCATE(GradQM)
    IF (ALLOCATED(GradVdW)) DEALLOCATE(GradVdW)
    IF (ALLOCATED(HessVdW)) DEALLOCATE(HessVdW)
    ALLOCATE(GradQM(N),GradVdW(N),HessVdW(N,N))

    READ(UCon) N
    IF ((N /= SIZE(DefCoord,1)) .OR. (N /= SIZE(Geometria)) .OR. &
        (N /=SIZE(Inc))) &
       CALL Mensaje('OptimizarGeometria',26,.TRUE.)
    READ(UCon) DefCoord(:,:)
    READ(UCon) Geometria(:)
    READ(UCon) Inc(:)
    IF (ALLOCATED(Gradiente)) DEALLOCATE(Gradiente)
    IF (ALLOCATED(Hessiana)) DEALLOCATE(Hessiana)
    ALLOCATE(Gradiente(N),Hessiana(N,N))

    READ(UCon) Gradiente(:)
    READ(UCon) Hessiana(:,:)

    READ(UCon) Energ,EnergiaQM,EnergiaVdW,DipoloQM,Hess

    READ(UCon) N
    IF (N /= SIZE(MolQM)) CALL Mensaje('OptimizarGeometria',26,.TRUE.)
    READ(UCon) MolQM(:)%q
    READ(UCon) (MolQM(i)%pos(:),i=1,SIZE(MolQM))

    READ(UCon) GradQM(:)
    READ(UCon) GradVdW(:)
    READ(UCon) HessVdW(:,:)
    CLOSE(UCon)

    WRITE(Fich,100) TRIM(Textos(22))

  100 FORMAT(A)

  END SUBROUTINE LeerContinuacion

END SUBROUTINE OptimizarGeometria

!-------------------------------------------------------------------------------
! Función que obtiene la energía para un cierto desplazamiento en la dirección
! dada por Inc
!-------------------------------------------------------------------------------
! Factor:    Factor que multiplica a Inc
! Incr:      Desplazamiento total en coordenadas de trabajo
! IncCart:   Desplazamiento total en coordenadas cartesianas
! CoordCart: Nuevas coordenadas cartesianas
!-------------------------------------------------------------------------------
FUNCTION EnergPunto(Factor)
  USE Coordenadas
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Factor
  DOUBLE PRECISION :: EnergPunto

  DOUBLE PRECISION, DIMENSION(SIZE(Inc,1)) :: Incr
  DOUBLE PRECISION, DIMENSION(3*SIZE(MolQM,1)) :: IncCart,CoordCart

  !Convierte las coordenadas del origen
  CALL ConvertirCoordenadas(CoordPrev(:),1)
  !Calcula y convierte el desplazamiento total
  Incr(:)=Factor*Inc(:)
  CALL ConvertirIncremento(Incr(:),CoordPrev(:),IncCart(:))
  !Modifica las coordenadas cartesianas
  CoordCart(:)=CoordPrev(:)+IncCart(:)
  CALL CambiarCartesianas(MolQM,CoordCart(:))
  !Hace el cálculo de energía
  CALL EjecutarQM(0)
  EnergPunto=EnergiaQM+EnergiaVdW

END FUNCTION EnergPunto

!-------------------------------------------------------------------------------
! Se calcula en una interpolación polinómica de las energías y derivadas en dos
! puntos (x=0 y x=1)
!-------------------------------------------------------------------------------
! Func0,Func1:    Energías en los puntos x=0 y x=1
! Der0,Der1:      Derivadas en los puntos x=0 y x=1
! Minim:          Mínimo de la interpolación
! Error:          Variable para controlar los errores
! Preci:          Precisión para comparar con cero
! Ca,Cb,Cc,Cd,Ce: Coeficientes del polinomio interpolador
! Aux1,Aux2:      Variables auxiliares
! DerM,DDerM:     Primera y segunda derivada en el mínimo
!-------------------------------------------------------------------------------
! Error=-1: Se realizó una interpolación cubica
! Error=0:  Se realizó una interpolación de cuarto grado
! Error=1:  Ninguna interpolación fue válida
!-------------------------------------------------------------------------------
SUBROUTINE MinimoInterpolado(Func0,Func1,Der0,Der1,Minim,Error)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: Func0,Func1,Der0,Der1
  DOUBLE PRECISION, INTENT(OUT) :: Minim
  INTEGER, INTENT(OUT) :: Error

  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-12
  DOUBLE PRECISION Ca,Cb,Cc,Cd,Ce,Aux1,Aux2,DerM,DDerM

  !Valores iniciales
  Error=1
  Minim=1.0D0

  !Se ajustan los valores a un polinomio de cuarto grado con la condición de
  !que la segunda derivada se anule en un solo punto
  Aux1=Func1-Func0-Der0
  Aux2=Der1-Der0
  Ce=Func0
  Cd=Der0
  !Para que se anule en un solo punto, tiene que cumplirse que:
  ! 9B^2-24AC = 0  (A=Y-3X+C; B=4X-Y-2C)
  !y por lo tanto:
  ! 6XY-6X^2-Y^2 > 0
  Cb=6.0D0*Aux1*(Aux2-Aux1)-Aux2*Aux2
  IF (Cb >= 0.0D0) THEN
    Cc=0.5D0*(6.0D0*Aux1-Aux2+SQRT(2.0D0*Cb))
    Ca=Aux2-3.0D0*Aux1+Cc
    !Hay dos soluciones, se escoge la que da menor A, con A y C positivos
    Cb=0.5D0*(6.0D0*Aux1-Aux2-SQRT(2.0D0*Cb))
    IF ((Aux2-3.0D0*Aux1+Cb > 0.0D0) .AND. (Cb >= 0.0D0) .AND. &
        ((Cc < 0.0D0) .OR. (ABS(Aux2-3.0D0*Aux1+Cb) < ABS(Cc)))) THEN
      Cc=Cb
      Ca=Aux2-3.0D0*Aux1+Cc
    END IF
    Cb=4.0D0*Aux1-Aux2-2.0D0*Cc
    !Se asegura de que la segunda derivada es efectivamente positiva
    IF ((Ca > Preci) .AND. (Cc >= 0.0D0)) THEN
      !Se resuelve una ecuación cúbica con una sola raíz real
      Aux1=(Cc/(2.0D0*Ca)-3.0D0*Cb*Cb/(16.0D0*Ca*Ca))/3.0D0
      Aux2=(Cd/(4.0D0*Ca)+Cb*Cb*Cb/(32.0D0*Ca*Ca*Ca)-Cb*Cc/(8.0D0*Ca*Ca))/2.0D0
      Aux1=SQRT(Aux2*Aux2+Aux1*Aux1*Aux1)
      Minim=(Aux1-Aux2)**(1.0D0/3.0D0)-(Aux1+Aux2)**(1.0D0/3.0D0)-Cb/(4.0D0*Ca)
      !Se comprueban las derivadas en el mínimo
      DerM=4.0D0*Ca*Minim**3+3.0D0*Cb*Minim**2+2.0D0*Cc*Minim+Cd
      DDerM=12.0D0*Ca*Minim**2+6.0D0*Cb*Minim+2.0D0*Cc
      IF ((DerM <= Preci) .AND. (DDerM > 0.0D0)) Error=0
    END IF
  END IF

  !Si el polinomio de cuarto grado no es válido, se hace un ajuste cúbico
  IF (Error == 1) THEN
    Cd=Func0
    Cc=Der0
    Cb=3.0D0*(Func1-Func0)-2.0D0*Der0-Der1
    Ca=2.0D0*(Func0-Func1)+Der0+Der1
    !Se comprueba si existen puntos críticos
    Ce=Cb*Cb-3.0D0*Ca*Cc
    IF (Ce >= 0.0D0) THEN
      IF (ABS(Ca) > Preci) THEN
        !Polinomio de tercer grado, el mínimo es siempre +
        Minim=(-Cb+SQRT(Ce))/(3.0D0*Ca)
      ELSE IF (ABS(Cb) > Preci) THEN
        !Polinomio de segundo grado
        Minim=-0.5D0*Cc/Cb
      END IF
      !Se comprueban las derivadas en el mínimo
      DerM=3.0D0*Ca*Minim**2+2.0D0*Cb*Minim+Cc
      DDerM=6.0D0*Ca*Minim+2.0D0*Cb
      IF ((DerM <= Preci) .AND. (DDerM > 0.0D0)) Error=-1
    END IF
  END IF

  !Si ninguno de los ajustes ha tenido éxito
  IF (Error == 1) Minim=1.0D0

END SUBROUTINE MinimoInterpolado

!-------------------------------------------------------------------------------
! Busca el mínimo de energía con una sola variable por el método de Brent.
!-------------------------------------------------------------------------------
! EnergPunto:  Función que devuelve el valor de la energía con una variable
! A,B:         El intervalo de la búsqueda.
! C:           La posición del mínimo
! Preci:       Precisión relativa a alcanzar en el mínimo
! Error:       Variable para controlar los errores
! Aureo:       Número áureo, (1+SQRT(5))/2
! IterMax:     Número máximo de iteraciones para encontrar el mínimo
! Valor:       Seis puntos con sus valores de energía correspondientes
! Aux:         Variable auxiliar
! Med:         Posición media entre dos valores
! D,E,F,P,Q,R: Variables utilizadas por los algoritmos
! Iter:        Número de iteración
! Fin:         Variable para controlar cuándo ha de finalizar la búsqueda
!-------------------------------------------------------------------------------
SUBROUTINE BuscarMinimo(EnergPunto,A,B,C,Preci)
  USE Utilidades
  IMPLICIT NONE
  INTERFACE
    DOUBLE PRECISION FUNCTION EnergPunto(Factor)
      DOUBLE PRECISION, INTENT(IN) :: Factor
    END FUNCTION EnergPunto
  END INTERFACE
  DOUBLE PRECISION, INTENT(IN) :: A,B
  DOUBLE PRECISION, INTENT(OUT) :: C
  DOUBLE PRECISION, INTENT(IN) :: Preci

  DOUBLE PRECISION, PARAMETER :: Aureo=1.618034D0
  INTEGER, PARAMETER :: IterMax=100
  DOUBLE PRECISION, DIMENSION(6,2) :: Valor
  DOUBLE PRECISION, DIMENSION(2) :: Aux
  DOUBLE PRECISION :: Med,D,E,F,P,Q,R
  INTEGER :: Iter,Fin

  !Se toman dos valores iniciales para comenzar la búsqueda
  Valor(1,1)=A
  Valor(2,1)=B
  Valor(1,2)=EnergPunto(Valor(1,1))
  Valor(2,2)=EnergPunto(Valor(2,1))

  !Se buscan dos puntos (1 y 3) entre los cuales esta el mínimo.
  IF (Valor(2,2) > Valor(1,2)) THEN
    Aux(:)=Valor(1,:)
    Valor(1,:)=Valor(2,:)
    Valor(2,:)=Aux(:)
  END IF
  Valor(3,1)=Valor(2,1)+Aureo*(Valor(2,1)-Valor(1,1))
  Valor(3,2)=EnergPunto(Valor(3,1))
  DO WHILE (Valor(2,2) >= Valor(3,2))
    !El punto 4 es el mínimo de la parábola que pasa por 1, 2 y 3
    R=(Valor(2,1)-Valor(1,1))*(Valor(2,2)-Valor(3,2))
    Q=(Valor(2,1)-Valor(3,1))*(Valor(2,2)-Valor(1,2))
    Valor(4,1)=(Valor(2,1)-Valor(1,1))*R-(Valor(2,1)-Valor(3,1))*Q
    Valor(4,1)=Valor(2,1)-Valor(4,1)/(2*(R-Q))
    Valor(4,2)=EnergPunto(Valor(4,1))
    Fin=0
    !Si 4 está entre 2 y 3
    IF ((Valor(4,1)-Valor(2,1))*(Valor(4,1)-Valor(3,1)) < 0.0D0) THEN
      IF (Valor(4,2) < Valor(3,2)) THEN
        Valor(1,:)=Valor(2,:)
        Valor(2,:)=Valor(4,:)
        Fin=1
      ELSE IF (Valor(4,2) > Valor(2,2)) THEN
        Valor(3,:)=Valor(4,:)
        Fin=1
       ELSE
        Valor(4,1)=Valor(3,1)+Aureo*(Valor(3,1)-Valor(2,1))
        Valor(4,2)=EnergPunto(Valor(4,1))
      END IF
    !Si 4 no está más allá de 3
    ELSE IF ((Valor(4,1)-Valor(3,1))*(Valor(3,1)-Valor(1,1)) < 0.0D0) THEN
      Valor(4,1)=Valor(3,1)+Aureo*(Valor(3,1)-Valor(2,1))
      Valor(4,2)=EnergPunto(Valor(4,1))
    END IF
    !Se toman los puntos 2, 3 y 4
    IF (Fin == 0) THEN
      Valor(1,:)=Valor(2,:)
      Valor(2,:)=Valor(3,:)
      Valor(3,:)=Valor(4,:)
    END IF
  END DO

  !Se localiza con mayor precisión la posición del mínimo
  IF (Valor(1,1) > Valor(3,1)) THEN
    Aux(:)=Valor(1,:)
    Valor(1,:)=Valor(3,:)
    Valor(3,:)=Aux(:)
  END IF
  IF (Valor(1,2) < Valor(3,2)) THEN
    Valor(5,:)=Valor(3,:)
    Valor(6,:)=Valor(1,:)
   ELSE
    Valor(5,:)=Valor(1,:)
    Valor(6,:)=Valor(3,:)
  END IF
  Fin=0
  E=Valor(3,1)-Valor(1,1)
  D=E
  DO Iter=1,IterMax
    F=E
    E=D
    Med=0.5D0*(Valor(1,1)+Valor(3,1))
    !IF (ABS(Valor(3,1)-Valor(1,1)) < 2.0D0*Preci) THEN
    IF (ABS(MAX(Valor(3,2),Valor(1,2))-Valor(2,2)) < Preci) THEN
      Valor(2,1)=Med
      Fin=1
      EXIT
    END IF
    !Se calcula la aproximación parabólica
    R=(Valor(2,1)-Valor(5,1))*(Valor(2,2)-Valor(6,2))
    Q=(Valor(2,1)-Valor(6,1))*(Valor(2,2)-Valor(5,2))
    P=(Valor(2,1)-Valor(5,1))*R-(Valor(2,1)-Valor(6,1))*Q
    Q=2*(Q-R)
    IF (Q == 0.0D0) THEN
      D=-SIGN(Preci,E)
     ELSE
      D=P/Q
    END IF
    !Si no es adecuada, se da un paso "áureo"
    IF ((ABS(D) > ABS(0.5D0*F)) .OR. &
        (((Valor(2,1)+D)-Valor(1,1))*((Valor(2,1)+D)-Valor(3,1)) > 0.0D0)) THEN
      IF (Valor(2,1) < Med) THEN
        D=(2.0D0-Aureo)*(Valor(3,1)-Valor(2,1))
       ELSE
        D=(2.0D0-Aureo)*(Valor(1,1)-Valor(2,1))
      END IF
    END IF
    !Se comprueba la energía y se actúa en consecuencia
    Valor(4,1)=Valor(2,1)+D
    Valor(4,2)=EnergPunto(Valor(4,1))
    IF (Valor(4,2) <= Valor(2,2)) THEN
      IF (D < 0.0D0) THEN
        Valor(3,:)=Valor(2,:)
       ELSE
        Valor(1,:)=Valor(2,:)
      END IF
      Valor(5,:)=Valor(6,:)
      Valor(6,:)=Valor(2,:)
      Valor(2,:)=Valor(4,:)
     ELSE
      IF (D < 0.0D0) THEN
        Valor(1,:)=Valor(4,:)
       ELSE
        Valor(3,:)=Valor(4,:)
      END IF
      IF (Valor(4,2) <= Valor(6,2)) THEN
        Valor(5,:)=Valor(6,:)
        Valor(6,:)=Valor(4,:)
      ELSE IF (Valor(4,2) <= Valor(5,2)) THEN
        Valor(5,:)=Valor(4,:)
      END IF
    END IF
  END DO
  IF (Fin == 0) CALL Mensaje('BuscarMinimo',2,.FALSE.)

  !Se da la salida
  C=Valor(2,1)

END SUBROUTINE BuscarMinimo

!-------------------------------------------------------------------------------
! Calcula el desplazamiento en las coordenadas para el próximo paso por el
! método RFO
!  J. Baker; J. Comput. Chem. 7 (1968) 385-395
!-------------------------------------------------------------------------------
! Hess:      Matriz hessiana
! Grad:      Vector del gradiente
! Modo:      Modo que se maximizará, 0 para localizar un mínimo
! Newt:      0 -> método RFO, 1 -> método Newton-Raphson
! Paso:      Incremento calculado para el próximo paso
! Error:     Variable para controlar los errores
! PasoNewt:  Paso Newton-Raphson
! Mat:       Copia de la hessiana para diagonalizar
! Vec:       Matriz que contiene los vectores propios de la hessiana
! Val:       Vector que contiene los valores propios de la hessiana
! Desp:      Valores propios desplazados
! Fuer:      Proyección del gradiente sobre los vectores propios
! Preci:     Precisión para no dividir por cero
! MaxIt:     Número máximo de iteraciones para el cálculo de lambda
! Suma:      Variable para el cálculo de lambda
! Lambda:    Desplazamiento de los valores propios
! Tau:       Desplazamiento del valor propio a maximizar
! Lim1,Lim2: Variables auxiliares para el cálculo de lambda
! Num:       Dimensión de las matrices y vectores
! i,j,k:     Contadores
!-------------------------------------------------------------------------------
! Error=-1: Se tomó el paso Newton-Raphson
! Error=0:  Todo correcto
!-------------------------------------------------------------------------------
SUBROUTINE IncrementoRFO(Hess,Grad,Modo,Newt,Paso,Error)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Hess
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Grad
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Paso
  INTEGER, INTENT(IN) :: Modo,Newt
  INTEGER, INTENT(OUT) :: Error

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Mat,Vec
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Val,Desp,Fuer,PasoNewt
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-8
  DOUBLE PRECISION :: Suma,Lambda,Tau,Lim1,Lim2
  INTEGER :: Num,i,j,k,MaxIt

  MaxIt=10000
  Error=0
  Num=SIZE(Grad,1)
  ALLOCATE(Mat(Num,Num),Vec(Num,Num),Val(Num),Desp(Num),Fuer(Num),PasoNewt(Num))

  !Se diagonaliza la hessiana
  Mat(:,:)=Hess(:,:)
  CALL Diagonalizar(Mat,Vec,Val,1)
  DEALLOCATE(Mat)

  !Se calcula la proyección del gradiente sobre los modos propios
  Fuer(:)=MATMUL(TRANSPOSE(Vec(:,:)),Grad(:))

  !Se calcula el paso por el método Newton-Raphson
  Desp(:)=-ABS(Val(:))
  IF (Modo > 0) Desp(Modo)=-Desp(Modo)
  PasoNewt(:)=0.0D0
  DO i=1,Num
    IF (ABS(Desp(i)) < Preci) CYCLE
    PasoNewt(:)=PasoNewt(:)+Vec(:,i)*Fuer(i)/Desp(i)
  END DO

  IF (Newt == 0) THEN

    !Se calcula Lambda: Sum(Fuer[i]^2/(Lambda-Val[i]))=Lambda
    j=1
    IF (Modo == 1) j=2
    Tau=5.0D-2
    Lim1=Val(j)
    Lim2=-1.0D3
    IF (Val(j) > 0.0D0) THEN
      Lambda=0.0D0
     ELSE
      Lambda=Val(j)-Tau
    END IF
    DO k=1,MaxIt
      Suma=0.0D0
      DO i=1,Num
        IF (i == Modo) CYCLE
        Suma=Suma+Fuer(i)*Fuer(i)/(Lambda-Val(i))
      END DO
      IF (ABS(Suma-Lambda) < Preci) EXIT
      IF (Val(j) > 0.0D0) THEN
        Lambda=Suma
       ELSE
        IF (Suma < Lambda) THEN
          Lim1=Lambda
         ELSE
          Lim2=Lambda
        END IF
        IF (Lim2 > -1.0D3) THEN
          Lambda=0.5D0*(Lim1+Lim2)
         ELSE
          Lambda=Lim1-Tau
        END IF
      END IF
    END DO

    !Se comprueba que el resultado es correcto
    Suma=0.0D0
    DO i=1,Num
      IF (i == Modo) CYCLE
      Suma=Suma+Fuer(i)*Fuer(i)/(Lambda-Val(i))
    END DO
    IF (((ABS(Lambda-Suma) > Preci) .AND. (Val(j)-Lambda > Preci)) .OR. &
        (Lambda > Val(j))) THEN
      CALL Mensaje('IncrementoRFO',13,.TRUE.)
    END IF

    !Se calcula Tau si se va a optimizar un estado de transición
    IF (Modo > 0) THEN
      Tau=0.5D0*Val(Modo)
      Tau=Tau+SQRT(Tau*Tau+Fuer(Modo)*Fuer(Modo))
    END IF

    !Se calcula el paso por el método P-RFO
    Desp(:)=Lambda-Val(:)
    IF (Modo > 0) Desp(Modo)=Tau-Val(Modo)
    Paso(:)=0.0D0
    DO i=1,Num
      IF (ABS(Desp(i)) < Preci) CYCLE
      Paso(:)=Paso(:)+Vec(:,i)*Fuer(i)/Desp(i)
    END DO

    !Si el paso es demasiado grande, es mejor el Newton-Raphson
    IF (SUM(PasoNewt(:)*PasoNewt(:))/SUM(Paso(:)*Paso(:)) < Preci) THEN
      Error=-1
      CALL Mensaje('IncrementoRFO',14,.FALSE.)
    END IF
  END IF

  IF ((Newt /= 0) .OR. (Error == -1)) Paso(:)=PasoNewt(:)

  DEALLOCATE(Vec,Val,Desp,Fuer,PasoNewt)

END SUBROUTINE IncrementoRFO

!-------------------------------------------------------------------------------
! Actualiza la matriz Hessiana por distintos métodos
! 0 -> Ninguna
! 1 -> BFGS
! 2 -> MS
! 3 -> PSB
! 4 -> Bofill
!-------------------------------------------------------------------------------
! Inc:     Diferencia entre las dos últimas posiciones
! IncGrad: Diferencia entre los dos últimos gradientes
! Hess:    Matriz hessiana
! Tipo:    Método utilizado para la actualización
! Cor:     Matriz de corrección
! Aux:     Vector auxiliar
! a,b:     Escalares
! i,j:     Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ActualizarHessiana(Inc,IncGrad,Hess,Tipo)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Inc,IncGrad
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: Hess
  INTEGER, INTENT(IN) :: Tipo

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Cor
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Aux
  DOUBLE PRECISION :: a,b
  INTEGER :: i,j

  ALLOCATE(Cor(SIZE(Inc,1),SIZE(Inc,1)),Aux(SIZE(Inc,1)))

  SELECT CASE (Tipo)

   CASE (1) !Broyden-Fletcher-Goldfarb-Shanno
    Aux(:)=MATMUL(Hess(:,:),Inc(:))
    Cor(:,:)=ProductoTens(IncGrad(:),IncGrad(:))/ &
             DOT_PRODUCT(IncGrad(:),Inc(:)) - &
             ProductoTens(Aux(:),Aux(:))/DOT_PRODUCT(Aux(:),Inc(:))

   CASE (2) !Murtagh-Sargent
    Aux(:)=IncGrad(:)-MATMUL(Hess(:,:),Inc(:))
    Cor(:,:)=ProductoTens(Aux(:),Aux(:))/DOT_PRODUCT(Aux(:),Inc(:))

   CASE (3) !Powell-symmetric-Broyden
    b=DOT_PRODUCT(Inc(:),Inc(:))
    Aux(:)=IncGrad(:)-MATMUL(Hess(:,:),Inc(:))
    Cor(:,:)=(ProductoTens(Aux(:),Inc(:))+ProductoTens(Inc(:),Aux(:)))/b - &
             DOT_PRODUCT(Aux(:),Inc(:))*ProductoTens(Inc(:),Inc(:))/(b*b)

   CASE (4) !Bofill
    b=DOT_PRODUCT(Inc(:),Inc(:))
    Aux(:)=IncGrad(:)-MATMUL(Hess,Inc(:))
    a=DOT_PRODUCT(Aux(:),Inc(:))**2/(b*DOT_PRODUCT(Aux(:),Aux(:)))
    Cor(:,:)=a*ProductoTens(Aux(:),Aux(:))/DOT_PRODUCT(Aux(:),Inc(:)) + &
             (1.0D0-a)* &
             (ProductoTens(Aux(:),Inc(:))+ProductoTens(Inc(:),Aux(:)))/b - &
             (1.0D0-a)* &
             DOT_PRODUCT(Aux(:),Inc(:))*ProductoTens(Inc(:),Inc(:))/(b*b)

   CASE DEFAULT !Ninguna
    Cor(:,:)=0.0D0
  END SELECT

  Hess(:,:)=Hess(:,:)+Cor(:,:)

  !Asegura que la matriz es simétrica
  DO i=1,SIZE(Inc,1)
    DO j=i+1,SIZE(Inc,1)
      Hess(i,j)=0.5D0*(Hess(i,j)+Hess(j,i))
      Hess(j,i)=Hess(i,j)
    END DO
  END DO

  DEALLOCATE(Cor,Aux)

END SUBROUTINE ActualizarHessiana

!-------------------------------------------------------------------------------
! Calcula una hessiana inicial aproximada
! 1 -> Matriz unidad
! 2 -> Hessiana simple (matriz unidad escalada según el tipo de coordenada)
! 3 -> Hessiana modelo (R. Lindh et al.; Chem. Phys. Lett. 241 (1995) 423-428)
!-------------------------------------------------------------------------------
! Tipo:                   Tipo de hessiana inicial que se calculará
! KDis,KAng,KDie,Aij,Rij: Constantes para la hessiana modelo
! a,b,c,d:                Átomos implicados
! i:                      Contador
!-------------------------------------------------------------------------------
SUBROUTINE HessianaInicial(Tipo)
  USE Coordenadas
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Tipo

  INTEGER :: i,a,b,c,d
  DOUBLE PRECISION, PARAMETER :: KDis=0.45D0,KAng=0.15D0,KDie=0.005
  DOUBLE PRECISION, DIMENSION(3,3), PARAMETER :: &
    Aij=RESHAPE((/ 1.0000D0,0.3949D0,0.3949D0, &
                   0.3949D0,0.2800D0,0.2800D0, &
                   0.3949D0,0.2800D0,0.2800D0 /), (/3,3/)), &
    Rij=RESHAPE((/ 1.35D0,2.10D0,2.53D0, &
                   2.10D0,2.87D0,3.40D0, &
                   2.53D0,3.40D0,3.40D0 /), (/3,3/))

  !Se crea la matriz hessiana
  IF (ALLOCATED(Hessiana)) DEALLOCATE(Hessiana)
  ALLOCATE(Hessiana(SIZE(DefCoord,1),SIZE(DefCoord,1)))
  Hessiana(:,:)=0.0D0

  SELECT CASE(Tipo)

   CASE (2) !Hessiana simple (unidad escalada)
    DO i=1,SIZE(Hessiana,1)
      IF (DefCoord(i,4) > 0) THEN      !Diedros
        Hessiana(i,i)=0.1D0
      ELSE IF (DefCoord(i,3) > 0) THEN !Ángulos
        Hessiana(i,i)=0.2D0
      ELSE IF (DefCoord(i,2) > 0) THEN !Distancias
        Hessiana(i,i)=0.5D0
       ELSE
        Hessiana(i,i)=1.0D0 !Hessiana unidad para coordenadas cartesianas
      END IF
    END DO

   CASE (3) !Hessiana modelo (Chem. Phys. Lett. 241, 423)
            !(debería incluir todas las distancias y ángulos)
    DO i=1,SIZE(Hessiana,1)
      a=DefCoord(i,1)
      b=DefCoord(i,2)
      c=DefCoord(i,3)
      d=DefCoord(i,4)
      IF (d > 0) THEN      !Diedros
        Hessiana(i,i)=KDie*Rho(a,b)*Rho(b,c)*Rho(c,d)
      ELSE IF (c > 0) THEN !Ángulos
        Hessiana(i,i)=KAng*Rho(a,b)*Rho(b,c)
      ELSE IF (b > 0) THEN !Distancias
        Hessiana(i,i)=KDis*Rho(a,b)
       ELSE
        Hessiana(i,i)=1.0D0 !Hessiana unidad para coordenadas cartesianas
      END IF
    END DO

   CASE DEFAULT !Hessiana unidad
    DO i=1,SIZE(Hessiana,1)
      Hessiana(i,i)=1.0D0
    END DO
  END SELECT

  CONTAINS

  !Función que devuelve rho_ij para la hessiana modelo
  !At1,At2:   Átomos implicados
  !Dist:      Distancia entre los dos átomos
  !Alfa,Rref: Constantes para este par de átomos
  !Per1,Per2: Periodos a los que pertenecen los dos átomos
  FUNCTION Rho(At1,At2)
    USE Utilidades
    DOUBLE PRECISION :: Rho
    INTEGER, INTENT(IN) :: At1,At2

    DOUBLE PRECISION :: Dist,Alfa,Rref
    INTEGER :: Per1,Per2

    !Se obtienen los periodos de los dos átomos implicados
    !(z = 0 se supone del periodo 1, para periodo > 3 se usa 3)
    SELECT CASE (MolQM(At1)%z)
     CASE (0:2)
      Per1=1
     CASE (3:10)
      Per1=2
     CASE (11:18)
      Per1=3
     CASE DEFAULT
      Per1=3
    END SELECT
    SELECT CASE (MolQM(At2)%z)
     CASE (0:2)
      Per2=1
     CASE (3:10)
      Per2=2
     CASE (11:18)
      Per2=3
     CASE DEFAULT
      Per2=3
    END SELECT

    !Se calcula el parámetro rho_ij
    Alfa=Aij(Per1,Per2)
    Rref=Rij(Per1,Per2)
    Dist=Distancia(MolQM(At1)%pos(:),MolQM(At2)%pos(:))
    Rho=EXP(Alfa*(Rref*Rref-Dist*Dist))
  END FUNCTION Rho

END SUBROUTINE HessianaInicial

END MODULE Optimizacion
