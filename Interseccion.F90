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

MODULE Interseccion

#ifndef LANG
#define LANG Espanol
#endif

#ifdef __NAG__
USE F90_UNIX_IO
#endif

USE Ejecutar
USE Optimizacion
USE LANG

CONTAINS
!OptimizarInterseccion(Mol,Fich)

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
SUBROUTINE OptimizarInterseccion(Mol,Fich)
  USE TipoAtomo
  USE Parametros
  USE Coordenadas
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(INOUT) :: Mol
  INTEGER, INTENT(IN) :: Fich

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: CoordCart,GradCart,IncCart, &
                                                 IncPrev,GradPrev,CargasPrev, &
                                                 GradF,Dir
  DOUBLE PRECISION, DIMENSION(3) :: DipoloPrev
  DOUBLE PRECISION :: Energ,EnergPrev,Paso,Conv1,Conv2,Conv3,Alfa,Sigma
  LOGICAL :: ModCoord
  CHARACTER(LEN=LEN(Extension)) :: Ext
  INTEGER :: Num,Iter,Conv,Modo,Calc,Hess,UMod,i,j

  Alfa=0.02D0
  Sigma=3.5D0

  !Se inician los datos para los cálculos QM
  Num=SIZE(Mol,1)
  IF (ALLOCATED(MolQM)) DEALLOCATE(MolQM)
  ALLOCATE(MolQM(Num),CoordCart(3*Num),GradCart(3*Num),IncCart(3*Num), &
           CoordPrev(3*Num),CargasPrev(Num),GradF(3*Num),Dir(3*Num))
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
!IFG
      CALL CalcularFunciones(Energ,GradF,EnergiaQM,Energia2QM,GradQM,Grad2QM,Sigma,Alfa)
!      CALL FuncionGradiente(GradF,GradQM,Grad2QM,EnergiaQM,Energia2QM,Sigma,Alfa)
      CALL ConvertirGradHess(1,GradF(:))
      CALL Proyectar(Grad=Gradiente(:))
     ELSE
      !Se calcula o no la hessiana inicial según el valor de CalcHessiana
      IF ((CalcHessiana > 0) .OR. (HessInicial == 0)) THEN
        Hess=1
        CALL EjecutarQM(2)
!IFG
        CALL CalcularFunciones(Energ,GradF,EnergiaQM,Energia2QM,GradQM,Grad2QM,Sigma,Alfa)
!        CALL FuncionGradiente(GradF,GradQM,Grad2QM,EnergiaQM,Energia2QM,Sigma,Alfa)
        CALL ConvertirGradHess(1,GradF(:),HessQM(:,:))
       ELSE
        CALL EjecutarQM(1)
!IFG
        CALL CalcularFunciones(Energ,GradF,EnergiaQM,Energia2QM,GradQM,Grad2QM,Sigma,Alfa)
!        CALL FuncionGradiente(GradF,GradQM,Grad2QM,EnergiaQM,Energia2QM,Sigma,Alfa)
        CALL ConvertirGradHess(1,GradF(:))
        CALL HessianaInicial(HessInicial)
      END IF
      CALL Proyectar(Grad=Gradiente(:),Hess=Hessiana(:,:))
    END IF

    !Se convierte el gradiente a cartesianas
    CALL ConvertirGradHess(-1,GradCart(:))
!    Energ=0.0D0
    EnergPrev=Energ
    Inc(:)=0.0D0
  END IF

  CALL EscribirGeometria

  Conv=0
  DO
    Iter=Iter+1

    !Se actualizan energías y coordenadas
!    EnergPrev=Energ
!IFG
!    CALL FuncionEnergia(Energ,EnergiaQM,Energia2QM,Sigma,Alfa)
    CoordPrev(:)=CoordCart(:)
    IncPrev(:)=Inc(:)

    !Se comprueba la convergencia (gradiente y energía)
    Dir(:)=Normalizar(GradQM-Grad2QM)
    Conv1=DOT_PRODUCT(GradF,Dir)
    Conv2=Norma(GradF-Conv1*Dir)
    Conv3=Energ-EnergPrev
    Conv1=Conv1/Sigma
    WRITE(Fich,101) 'Paralelo',Conv1,'Eh/a0'
    WRITE(Fich,101) 'Perpendicular',Conv2,'Eh/a0'
    WRITE(Fich,101) 'Diferencia funcion',Conv3,'Eh'
    IF ((ABS(Conv1) < ConvGradOpt) .AND. &
        (ABS(Conv2) < ConvGradOpt) .AND. &
        (ABS(Conv3) < ConvEnerOpt)) THEN
      Sigma=4.0D0*Sigma
      MaxPasoOpt=0.5D0*MaxPasoOpt
      CALL HessianaInicial(HessInicial)
      WRITE(Fich,101) 'Nuevo sigma:',Sigma,''
    END IF
    CALL FLUSH(Fich)
    EnergPrev=Energ
!    Conv1=MAXVAL(ABS(GradCart(:)))
!    Conv3=Energ-EnergPrev
!    WRITE(Fich,101) TRIM(Textos(1)),Conv1,'Eh/a0'
!    WRITE(Fich,101) TRIM(Textos(2)),Conv3,'Eh'
!    CALL FLUSH(Fich)
!    IF ((Conv1 <= ConvGradOpt) .AND. &
!        ((ABS(Conv3) <= ConvEnerOpt) .OR. (Iter == 1))) THEN
!      Conv=1
!      EXIT
!    END IF

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
      GradCart(:)=GradQM(:)
      !Se hace un cálculo con el paso completo
      Paso=MIN(1.0D0,MaxPasoOpt/Norma(Inc))
      Inc(:)=Paso*Inc(:)
      CALL ConvertirIncremento(Inc(:),CoordPrev(:),IncCart(:))
      CoordCart(:)=CoordPrev(:)+IncCart(:)
      CALL CambiarCartesianas(MolQM,CoordCart(:))
      CALL ConvertirCoordenadas(CoordCart(:),1)
      CALL EjecutarQM(1)
      CALL ConvertirGradHess(1,GradQM(:))
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
!    Conv2=MAXVAL(ABS(IncCart))
!    WRITE(Fich,101) TRIM(Textos(3)),Conv2,'a0'
!    CALL FLUSH(Fich)
!    IF ((Conv1 <= ConvGradOpt) .AND. (Conv2 <= ConvPasoOpt)) THEN
!      !Si ha habido busqueda lineal, se recuperan los datos iniciales
!      IF (BusquedaLineal > 0) THEN
!        CALL CambiarCartesianas(MolQM,CoordPrev(:))
!        CALL ConvertirCoordenadas(CoordPrev(:),0)
!        MolQM(:)%q=CargasPrev(:)
!        EnergiaQM=Energ
!        GradQM(:)=GradCart(:)
!        Gradiente(:)=GradPrev(:)
!        DipoloQM(:)=DipoloPrev(:)
!      END IF
!      Conv=1
!      EXIT
!    END IF
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
!IFG
        CALL CalcularFunciones(Energ,GradF,EnergiaQM,Energia2QM,GradQM,Grad2QM,Sigma,Alfa)
!        CALL FuncionGradiente(GradF,GradQM,Grad2QM,EnergiaQM,Energia2QM,Sigma,Alfa)
        CALL ConvertirGradHess(1,GradF(:))
        CALL Proyectar(Grad=Gradiente(:))
      END IF
     ELSE
      !Métodos de segundo orden
      CALL ConvertirCoordenadas(CoordCart(:),2)
      IF (Hess == 1) THEN
        !Se calcula la hessiana
        CALL EjecutarQM(2)
!IFG
        CALL CalcularFunciones(Energ,GradF,EnergiaQM,Energia2QM,GradQM,Grad2QM,Sigma,Alfa)
!        CALL FuncionGradiente(GradF,GradQM,Grad2QM,EnergiaQM,Energia2QM,Sigma,Alfa)
        CALL ConvertirGradHess(1,GradF(:),HessQM(:,:))
        CALL Proyectar(Grad=Gradiente(:),Hess=Hessiana(:,:))
       ELSE
        !Se actualiza la hessiana
        IF (Calc == 1) THEN
          CALL EjecutarQM(1)
!IFG
          CALL CalcularFunciones(Energ,GradF,EnergiaQM,Energia2QM,GradQM,Grad2QM,Sigma,Alfa)
!          CALL FuncionGradiente(GradF,GradQM,Grad2QM,EnergiaQM,Energia2QM,Sigma,Alfa)
          CALL ConvertirGradHess(1,GradF(:))
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
      WRITE(Fich,102) ((HessQM(i,j),j=1,i),i=1,SIZE(HessQM,1))
     ELSE
      WRITE(Fich,100) TRIM(Textos(20))
      CALL ConvertirGradHess(-1,Hess=HessQM(:,:))
      WRITE(Fich,102) ((HessQM(i,j),j=1,i),i=1,SIZE(HessQM,1))
    END IF
  END IF

  Mol=MolQM
  Extension=Ext

  DEALLOCATE(MolQM,CoordCart,GradCart,IncCart,Inc,CoordPrev,IncPrev,GradPrev, &
             CargasPrev,GradF,Dir)

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
    WRITE(UCon) Energ,EnergiaQM,DipoloQM,Hess
    WRITE(UCon) SIZE(MolQM)
    WRITE(UCon) MolQM(:)%q
    WRITE(UCon) (MolQM(i)%pos(:),i=1,SIZE(MolQM))
    WRITE(UCon) GradQM(:)
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
    WRITE(Fich,107) 'Funcion energia:',Energ,'Eh'
    WRITE(Fich,107) 'Energia media:',0.5D0*(EnergiaQM+Energia2QM),'Eh'
    WRITE(Fich,107) 'Diferencia entre estados:',EnergiaQM-Energia2QM,'Eh'
    WRITE(Fich,107) 'Modulo de la funcion gradiente:',Norma(GradCart),'Eh/a0'
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
  107 FORMAT(A,T32,F18.10,1X,A)

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
    ALLOCATE(GradQM(N))

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

    READ(UCon) Energ,EnergiaQM,DipoloQM,Hess

    READ(UCon) N
    IF (N /= SIZE(MolQM)) CALL Mensaje('OptimizarGeometria',26,.TRUE.)
    READ(UCon) MolQM(:)%q
    READ(UCon) (MolQM(i)%pos(:),i=1,SIZE(MolQM))

    READ(UCon) GradQM(:)
    CLOSE(UCon)

    WRITE(Fich,100) TRIM(Textos(22))

  100 FORMAT(A)

  END SUBROUTINE LeerContinuacion

END SUBROUTINE OptimizarInterseccion

SUBROUTINE FuncionEnergia(E,E1,E2,S,A)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: E1,E2,S,A
  DOUBLE PRECISION, INTENT(OUT) :: E
  E=0.5D0*(E1+E2)+S*(E1-E2)**2.0D0/(A+E1-E2)
END SUBROUTINE FuncionEnergia

SUBROUTINE FuncionGradiente(G,G1,G2,E1,E2,S,A)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: G1,G2
  DOUBLE PRECISION, INTENT(IN) :: E1,E2,S,A
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: G
  G=0.5D0*(G1+G2)+S*((E1-E2)**2+2.0D0*A*(E1-E2))/(A+E1-E2)**2*(G1-G2)
END SUBROUTINE FuncionGradiente

SUBROUTINE CalcularFunciones(E,G,E1,E2,G1,G2,S,A)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: G1,G2
  DOUBLE PRECISION, INTENT(IN) :: E1,E2,S,A
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: G
  DOUBLE PRECISION, INTENT(OUT) :: E
  E=0.5D0*(E1+E2)+S*(E1-E2)**2.0D0/(A+E1-E2)
  G=0.5D0*(G1+G2)+S*((E1-E2)**2+2.0D0*A*(E1-E2))/(A+E1-E2)**2*(G1-G2)
END SUBROUTINE CalcularFunciones

END MODULE Interseccion
