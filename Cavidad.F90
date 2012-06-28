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

!-------------------------------------------------------------------------------
! Definición de la cavidad que contiene al soluto y a la primera capa de
! solvatación
!-------------------------------------------------------------------------------
MODULE Cavidad

!-------------------------------------------------------------------------------
! CavVert:  Posiciones de los vértices
! CavCent:  Posiciones de los centros de las caras (en la superficie)
! CavTrian: Índices de los vértices que forman cada cara
! Esferas:  Esferas que definen la cavidad
! Esferas2: Esferas auxiliares para las intersecciones de dos esferas
! Esferas3: Esferas auxiliares para las intersecciones de tres esferas
! Toros:    Toros para las intersecciones de dos esferas
! Int2:     Matriz con los índices de las intersecciones de dos esferas
! Int3:     Matriz con los índices de las intersecciones de tres esferas
! RDis:     Radio de la esfera que representa al disolvente (radio de exclusión)
! RMax:     Radio máximo de la cavidad
!-------------------------------------------------------------------------------
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CavVert,CavCent,Esferas, &
                                                 Esferas2,Esferas3,Toros
INTEGER, DIMENSION(:,:), ALLOCATABLE :: CavTrian,Int2
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: Int3
DOUBLE PRECISION :: RDis,RMax

CONTAINS
!ConstruirCavidad(Centros,Radios,R,NDiv)
!Interior(Pun)
!Desplazar(Pun,Vec)
!Distribuir
!Pentakisdodecaedro
!Dividir

!-------------------------------------------------------------------------------
! Define una cavidad adaptada a la forma del soluto, sin esquinas,
! y un conjunto de puntos en su superfice. Ver:
!  M.L. Connolly; J. Am. Chem. Soc. 107 (1985), 1118-1124
!  C.S. Pomelli, J. Tomasi; J. Comput. Chem. 19 (1998), 1758-1776
!-------------------------------------------------------------------------------
! Centros:           Posiciones de las esferas que definen la cavidad
! Radios:            Radios de las esferas que definen la cavidad
! R:                 Radio de la esfera para eliminar las esquinas
! NDiv:              Número de subdivisiones en el poliedro (opcional)
! N:                 Número de esferas que definen la cavidad
! Norm,Norm2,Vec:    Vectores auxiliares para los cálculos
! Rij:               Distancia entre dos esferas
! b,c,ae,be:         Variables auxiliares para los cálculos
! xf,xg,yg,xh,yh,zh: Variables auxiliares para los cálculos
! m:                 Número de intersecciones de dos o tres esferas
! UTmp:              Unidad temporal
! i,j,k:             Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ConstruirCavidad(Centros,Radios,R,NDiv)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Centros
  DOUBLE PRECISION, DIMENSION(SIZE(Centros,1)), INTENT(IN) :: Radios
  DOUBLE PRECISION, INTENT(IN) :: R
  INTEGER, INTENT(IN), OPTIONAL :: NDiv

  DOUBLE PRECISION, DIMENSION(3) :: Norm,Norm2,Vec
  DOUBLE PRECISION :: Rij,b,c,ae,be,xf,xg,yg,xh,yh,zh
  INTEGER :: N,m,UTmp,i,j,k

  !Prepara las matrices
  N=SIZE(Centros,1)
  IF (ALLOCATED(Esferas)) DEALLOCATE(Esferas)
  IF (ALLOCATED(Int2)) DEALLOCATE(Int2)
  IF (ALLOCATED(Int3)) DEALLOCATE(Int3)
  ALLOCATE(Esferas(N,4),Int2(N,N),Int3(N,N,N))

  !Define las esferas y RDis
  Esferas(:,1:3)=Centros(:,:)
  Esferas(:,4)=Radios(:)
  RDis=R

  !Calcula RMax
  RMax=0.0D0
  DO i=1,SIZE(Centros,1)
    RMax=MAX(RMax,Norma(Esferas(i,1:3))+Esferas(i,4))
  END DO

  !Se calculan las intersecciones de dos esferas
  UTmp=NuevaUnidad()
  OPEN(UTmp,STATUS='SCRATCH',ACTION='READWRITE')
  Int2(:,:)=0
  m=0
  IF (RDis > 0.0D0) THEN
    DO i=1,N
      DO j=i+1,N
        Norm(:)=Esferas(j,1:3)-Esferas(i,1:3)
        Rij=Norma(Norm(:))
        Norm(:)=Normalizar(Norm(:))
        !Si las esferas no están ni demasiado lejos ni demasiado cerca,
        !intersecan
        IF (Esferas(i,4)+Esferas(j,4)+2.0D0*RDis < Rij) CYCLE
        IF (ABS(Esferas(i,4)-Esferas(j,4)) > Rij) CYCLE
        !Se almacena el índice de la intersección
        m=m+1
        Int2(i,j)=m
        !Se calculan los parámetros geométricos
        ae=(Rij**2+(Esferas(i,4)+RDis)**2-(Esferas(j,4)+RDis)**2)/(2.0D0*Rij)
        be=SQRT((Esferas(i,4)+RDis)**2-ae**2)
        xg=Esferas(i,4)/(Esferas(i,4)+RDis)*ae
        yg=Esferas(i,4)/(Esferas(i,4)+RDis)*be
        xh=Rij-Esferas(j,4)/(Esferas(j,4)+RDis)*(Rij-ae)
        yh=Esferas(j,4)/(Esferas(j,4)+RDis)*be
        xf=(xg*xg+yg*yg-xh*xh-yh*yh)/(2.0D0*(xg-xh))
        !Se guarda el centro y el radio de la esfera auxiliar
        WRITE(UTmp,*) Esferas(i,1:3)+xf*Norm(:), &
                     SQRT(xg*xg+yg*yg-2.0D0*xg*xf+xf*xf)
        !Se guarda el centro y el radio mayor del toro auxiliar
        WRITE(UTmp,*) Esferas(i,1:3)+ae*Norm(:),be
      END DO
    END DO
  END IF

  !Se leen las esferas y toros auxiliares de las intersecciones de dos esferas
  IF (ALLOCATED(Esferas2)) DEALLOCATE(Esferas2)
  IF (ALLOCATED(Toros)) DEALLOCATE(Toros)
  ALLOCATE(Esferas2(m,4),Toros(m,4))
  REWIND(UTmp)
  DO i=1,m
    READ(UTmp,*) Esferas2(i,:)
    READ(UTmp,*) Toros(i,:)
  END DO
  CLOSE(UTmp)

  !Se calculan las intersecciones de tres esferas
  UTmp=NuevaUnidad()
  OPEN(UTmp,STATUS='SCRATCH',ACTION='READWRITE')
  Int3(:,:,:)=0
  m=0
  IF (RDis > 0.0D0) THEN
    DO i=1,N
      DO j=i+1,N
        IF (Int2(i,j) == 0) CYCLE
        Vec(:)=Esferas(j,1:3)-Esferas(i,1:3)
        Rij=Norma(Vec(:))
        Norm(:)=Vec(:)/Rij
        DO k=j+1,N
          !Sólo se considera la intersección si las esferas se cortan dos a dos
          IF ((Int2(i,k) == 0) .OR. (Int2(j,k) == 0)) CYCLE
          !Se calculan los parámetros geométricos
          Vec(:)=Esferas(k,1:3)-Esferas(i,1:3)
          b=DOT_PRODUCT(Vec(:),Norm(:))
          Norm2(:)=Vec(:)-b*Norm(:)
          c=Norma(Norm2(:))
          Norm2(:)=Norm2(:)/c
          xh=((Esferas(i,4)+RDis)**2-(Esferas(j,4)+RDis)**2+Rij**2)/(2.0D0*Rij)
          yh=((Esferas(i,4)+RDis)**2-(Esferas(k,4)+RDis)**2+ &
             (xh-b)**2)/(2.0D0*c)+0.5D0*c-xh**2/(2.0D0*c)
          zh=(Esferas(i,4)+RDis)**2-xh**2-yh**2
          IF ((zh <= 0.0D0) .OR. (c == 0.0D0)) CYCLE
          !Si existe intersección, se almacena el índice
          m=m+1
          Int3(i,j,k)=m
          zh=SQRT(zh)
          Vec(:)=ProductoVect(Norm(:),Norm2(:))
          !Se guarda el centro y el radio de las dos esferas auxiliares
          WRITE(UTmp,*) Esferas(i,1:3)+xh*Norm(:)+yh*Norm2(:)+zh*Vec(:)
          WRITE(UTmp,*) Esferas(i,1:3)+xh*Norm(:)+yh*Norm2(:)-zh*Vec(:)
        END DO
      END DO
    END DO
  END IF

  !Se leen las esferas auxiliares de las intersecciones de tres esferas
  IF (ALLOCATED(Esferas3)) DEALLOCATE(Esferas3)
  ALLOCATE(Esferas3(m,6))
  REWIND(UTmp)
  DO i=1,m
    READ(UTmp,*) Esferas3(i,1:3)
    READ(UTmp,*) Esferas3(i,4:6)
  END DO
  CLOSE(UTmp)

  !Si no hay que calcular los puntos en la superficie, se acaba
  IF (.NOT. PRESENT(NDiv)) RETURN

  !Se crea el poliedro
  CALL Pentakisdodecaedro
  !Se cambia de tamaño para que abarque toda la cavidad
  Rij=0.0D0
  DO i=1,N
    Rij=MAX(Rij,Norma(Esferas(i,1:3))+Esferas(i,4))
  END DO
  CavVert(:,:)=Rij*CavVert(:,:)

  !Se va subdividiendo
  N=0
  DO i=0,NDiv
    IF (i > 0) THEN
      N=SIZE(CavVert,1)
      CALL Dividir
    END IF
    !Se colocan los nuevos vértices en la superficie de la cavidad
!$OMP PARALLEL DO PRIVATE(j,Vec)
    DO j=N+1,SIZE(CavVert,1)
      Vec(:)=CavVert(j,:)
      CALL Desplazar(CavVert(j,:),Vec(:))
    END DO
!$OMP END PARALLEL DO
    !Se intentan redistribuir en la superficie
    CALL Distribuir
  END DO

  !Se colocan los centros de las caras curvas en la superficie
!$OMP PARALLEL DO PRIVATE(i,Vec)
  DO i=1,SIZE(CavCent,1)
    CavCent(i,:)=(CavVert(CavTrian(i,1),:)+CavVert(CavTrian(i,2),:)+ &
                  CavVert(CavTrian(i,3),:))/3.0D0
    Vec(:)=ProductoVect(CavVert(CavTrian(i,2),:)-CavVert(CavTrian(i,1),:), &
                        CavVert(CavTrian(i,3),:)-CavVert(CavTrian(i,1),:))
    CALL Desplazar(CavCent(i,:),Vec(:))
  END DO
!$OMP END PARALLEL DO

END SUBROUTINE ConstruirCavidad

!-------------------------------------------------------------------------------
! Función que devuelve T si un punto esta dentro de la cavidad y si no, F
!  M.L. Connolly; J. Am. Chem. Soc. 107 (1985), 1118-1124
!  C.S. Pomelli, J. Tomasi; J. Comput. Chem. 19 (1998), 1758-1776
!-------------------------------------------------------------------------------
! Pun:     Punto tridimensional que se comprueba
! N:       Número de esferas que definen la cavidad
! m:       Índice de la intersección correspondiente
! V:       Vector de posición del punto con respecto a un elemento geométrico
! CosG:    Coseno de gamma
! pepp:    Variable para comprobar la pertenencia a un toro
! a,b,c:   Variables auxiliares
! i,j,k:   Contadores
!-------------------------------------------------------------------------------
FUNCTION Interior(Pun)
  USE Utilidades
  IMPLICIT NONE
  LOGICAL :: Interior
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Pun

  DOUBLE PRECISION, DIMENSION(3) :: V
  DOUBLE PRECISION :: CosG,pepp,a,b,c
  INTEGER :: N,m,i,j,k

  !Si no hay cavidad definida, siempre está en el interior
  IF (.NOT. ALLOCATED(Esferas)) THEN
    Interior=.TRUE.
    RETURN
  END IF

  N=SIZE(Esferas,1)
  Interior=.FALSE.

  !Si está fuera de RMax, está fuera de la cavidad
  IF (Norma(Pun(:)) > RMax) RETURN

  !Si el punto está dentro de alguno de los elementos, sale del bucle
  Externo: DO i=1,N
    !Comprueba el interior de todas las esferas que componen la cavidad
    V(:)=Pun(:)-Esferas(i,1:3)
    IF (DOT_PRODUCT(V(:),V(:)) < Esferas(i,4)*Esferas(i,4)) THEN
      Interior=.TRUE.
      EXIT Externo
    END IF
    DO j=i+1,N
      m=Int2(i,j)
      IF (m == 0) CYCLE
      !Comprueba que esté dentro de la esfera auxiliar, pero fuera del toro
      V(:)=Pun(:)-Esferas2(m,1:3)
      IF (DOT_PRODUCT(V(:),V(:)) < Esferas2(m,4)*Esferas2(m,4)) THEN
        a=DOT_PRODUCT(Pun(:)-Toros(m,1:3),Pun(:)-Toros(m,1:3))
        CosG=DOT_PRODUCT(Normalizar(Esferas(j,1:3)-Esferas(i,1:3)), &
                         Normalizar(Pun(:)-Toros(m,1:3)))
        pepp=SQRT(Toros(m,4)*Toros(m,4)+a- &
                  2.0D0*Toros(m,4)*SQRT(a*(1.0D0-CosG*CosG)))
        IF (pepp >= RDis) THEN
          Interior=.TRUE.
          EXIT Externo
        END IF
      END IF
      DO k=j+1,N
        m=Int3(i,j,k)
        IF (m == 0) CYCLE
        !Comprueba a qué lado del plano ijk se encuentra el punto
        V(:)=ProductoVect(Esferas(j,1:3)-Esferas(i,1:3), &
                          Esferas(k,1:3)-Esferas(i,1:3))
        !En cada caso, comprueba que el punto está dentro de la pirámide,
        !pero fuera de la esfera
        IF (DOT_PRODUCT(Pun(:)-Esferas(i,1:3),V(:)) > 0.0D0) THEN
          V(:)=Pun(:)-Esferas3(m,1:3)
          a=DOT_PRODUCT(ProductoVect(Esferas(j,1:3)-Esferas3(m,1:3), &
                                     Esferas(i,1:3)-Esferas3(m,1:3)),V(:))
          b=DOT_PRODUCT(ProductoVect(Esferas(k,1:3)-Esferas3(m,1:3), &
                                     Esferas(j,1:3)-Esferas3(m,1:3)),V(:))
          c=DOT_PRODUCT(ProductoVect(Esferas(i,1:3)-Esferas3(m,1:3), &
                                     Esferas(k,1:3)-Esferas3(m,1:3)),V(:))
          IF ((a > 0.0D0) .AND. (b > 0.0D0) .AND. (c > 0.0D0) .AND. &
              (DOT_PRODUCT(V(:),V(:)) > RDis*RDis)) THEN
            Interior=.TRUE.
            EXIT Externo
          END IF
         ELSE
          V(:)=Pun(:)-Esferas3(m,4:6)
          a=DOT_PRODUCT(ProductoVect(Esferas(i,1:3)-Esferas3(m,4:6), &
                                     Esferas(j,1:3)-Esferas3(m,4:6)),V(:))
          b=DOT_PRODUCT(ProductoVect(Esferas(j,1:3)-Esferas3(m,4:6), &
                                     Esferas(k,1:3)-Esferas3(m,4:6)),V(:))
          c=DOT_PRODUCT(ProductoVect(Esferas(k,1:3)-Esferas3(m,4:6), &
                                     Esferas(i,1:3)-Esferas3(m,4:6)),V(:))
          IF ((a > 0.0D0) .AND. (b > 0.0D0) .AND. (c > 0.0D0) .AND. &
              (DOT_PRODUCT(V(:),V(:)) > RDis*RDis)) THEN
            Interior=.TRUE.
            EXIT Externo
          END IF
        END IF
      END DO
    END DO
  END DO Externo

END FUNCTION Interior

!-------------------------------------------------------------------------------
! Desplaza un punto en una dirección hasta que queda situado en la superficie
!-------------------------------------------------------------------------------
! Pun:     Punto original que se desplaza
! Vec:     Dirección en la que se desplaza el punto
! V:       Vector de desplazamiento del punto
! N:       Número por el que se divide el desplazamiento inicial
! R:       Nuevo punto
! Paso:    Magnitud del desplazamiento
! Delta:   Criterio de convergencia para el desplazamiento
!-------------------------------------------------------------------------------
SUBROUTINE Desplazar(Pun,Vec)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3), INTENT(INOUT) :: Pun
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Vec

  DOUBLE PRECISION, DIMENSION(3) :: V,R
  DOUBLE PRECISION :: Paso
  DOUBLE PRECISION, PARAMETER :: Delta=1.0D-4
  INTEGER, PARAMETER :: N=10

  !Se calcula el desplazamiento inicial
  V(:)=Vec(:)/DBLE(N)

  !Se desplaza el punto hasta que queda en el interior
  DO WHILE (.NOT. Interior(Pun(:)))
    Pun(:)=Pun(:)-V(:)
  END DO
  !Se desplaza el punto hasta que queda en el exterior
  DO WHILE (Interior(Pun(:)))
    Pun(:)=Pun(:)+V(:)
  END DO

  !Por bisección, se ajusta la posicion del punto hasa que queda en la
  !superficie
  Paso=Norma(Vec(:))
  DO WHILE (Paso > Delta)
    V(:)=0.5D0*V(:)
    Paso=0.5D0*Paso
    R(:)=Pun(:)-V(:)
    IF (.NOT. Interior(R(:))) Pun(:)=R(:)
  END DO
  Pun(:)=Pun(:)-0.5D0*V(:)

END SUBROUTINE Desplazar

!-------------------------------------------------------------------------------
! Distribuye más o menos regularmente los vértices por la superficie,
! se intenta que las longitudes de todas las aristas sean iguales
!-------------------------------------------------------------------------------
! Normales:  Vectores normales en los vértices
! Fuerzas:   Fuerzas calculadas en los vértices
! Vec:       Vector normal de cada cara
! MaxF:      Máxima fuerza en los vértices
! Lon:       Longitud promedio de las aristas
! MinL,MaxL: Longitudes mínima y máxima de las aristas
! F:         Factor para calcular desplazamientos a partir de las fuerzas
! Delta:     Criterio de convergencia para la fuerza máxima
! i,j,k,l:   Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Distribuir
  USE Utilidades
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Normales,Fuerzas
  DOUBLE PRECISION, DIMENSION(3) :: Vec
  DOUBLE PRECISION :: MaxF,A,Lon,MaxL,MinL
  INTEGER :: i,j,k,l
  DOUBLE PRECISION, PARAMETER :: F=0.2D0,Delta=2.0D-2

  ALLOCATE(Normales(SIZE(CavVert,1),3),Fuerzas(SIZE(CavVert,1),3))

  !Se repite el ciclo hasta que se alcanza la convergencia
  MaxF=1.0D0+Delta
  DO WHILE (MaxF > Delta)
    Normales(:,:)=0.0D0
    Fuerzas(:,:)=0.0D0
    Lon=0.0D0
    MaxL=TINY(MinL)
    MinL=HUGE(MaxL)
!$OMP PARALLEL PRIVATE(i,j,k,l,A,Vec)
!$OMP DO REDUCTION(+:Lon) REDUCTION(MAX:MaxL) REDUCTION(MIN:MinL)
    !Primero se calculan las normales de las caras
    !y las longitudes de las aristas
    DO i=1,SIZE(CavTrian,1)
      CavCent(i,:)=(CavVert(CavTrian(i,1),:)+CavVert(CavTrian(i,2),:)+ &
                    CavVert(CavTrian(i,3),:))/3.0D0
      Vec(:)=ProductoVect(CavVert(CavTrian(i,2),:)-CavVert(CavTrian(i,1),:), &
                          CavVert(CavTrian(i,3),:)-CavVert(CavTrian(i,1),:))
      DO j=1,3
        k=CavTrian(i,j)
        l=CavTrian(i,MOD(j,3)+1)
        A=Norma(CavVert(l,:)-CavVert(k,:))
        Lon=Lon+A
        MaxL=MAX(MaxL,A)
        MinL=MIN(MinL,A)
!$OMP CRITICAL
        !La normal en cada vértice sera el promedio de las de las caras
        Normales(k,:)=Normales(k,:)+Normalizar(Vec(:))
!$OMP END CRITICAL
      END DO
    END DO
!$OMP END DO
!$OMP SINGLE
    !Se calcula la longitud promedio (cada arista se ha contado dos veces)
    Lon=Lon/(3.0D0*SIZE(CavTrian,1))
    MaxF=0.0D0
!$OMP END SINGLE
!$OMP DO
    !Se calcula la fuerza en cada vértice
    DO i=1,SIZE(CavTrian,1)
      DO j=1,3
        k=CavTrian(i,j)
        l=CavTrian(i,MOD(j,3)+1)
        !La fuerza es proporcional a la desviación en la longitud de arista
        A=Norma(CavVert(l,:)-CavVert(k,:))-Lon
!$OMP CRITICAL
        Fuerzas(k,:)=Fuerzas(k,:)+A*Normalizar(CavVert(l,:)-CavVert(k,:))
!$OMP END CRITICAL
      END DO
    END DO
!$OMP END DO
!$OMP DO REDUCTION(MAX:MaxF)
    !Finalmente se desplazan los vértices
    DO i=1,SIZE(CavVert,1)
      Normales(i,:)=Normalizar(Normales(i,:))
      !Se elimina la componente normal de la fuerza
      Fuerzas(i,:)=Fuerzas(i,:)- &
                   DOT_PRODUCT(Fuerzas(i,:),Normales(i,:))*Normales(i,:)
      MaxF=MAX(MaxF,Norma(Fuerzas(i,:)))
      CavVert(i,:)=CavVert(i,:)+F*Fuerzas(i,:)
      CALL Desplazar(CavVert(i,:),Normales(i,:))
    END DO
!$OMP END DO
!$OMP END PARALLEL
    !La fuerza máxima se pone en relación con la longitud de arista
    MaxF=MaxF/Lon
  END DO

  DEALLOCATE(Normales,Fuerzas)

END SUBROUTINE Distribuir

!-------------------------------------------------------------------------------
! Genera los vértices de un pentakisdodecaedro (60 caras) de radio 1
!-------------------------------------------------------------------------------
! VZeta:   Ángulos zeta de las 6 "capas" de vértices
! CZeta:   Coseno de zeta
! SZeta:   Seno de zeta
! VFi:     Ángulo fi inicial de una "capa"
! DFi:     Diferencia de ángulos fi
! Fi:      Ángulo fi de un vértice
! i,j,k:   Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Pentakisdodecaedro
  USE Unidades
  IMPLICIT NONE

  DOUBLE PRECISION :: CZeta,SZeta,VFi,DFi,Fi
  DOUBLE PRECISION, DIMENSION(6) :: VZeta
  INTEGER :: i,j,k

  !Se crean las matrices necesarias
  IF (ALLOCATED(CavVert)) DEALLOCATE(CavVert)
  IF (ALLOCATED(CavCent)) DEALLOCATE(CavCent)
  IF (ALLOCATED(CavTrian)) DEALLOCATE(CavTrian)
  ALLOCATE(CavVert(32,3),CavCent(60,3),CavTrian(60,3))

  !Se definen los angulos zeta y fi
  VZeta(1)=ACOS(SQRT((5.0D0+2.0D0*SQRT(5.0D0))/15.0D0))
  VZeta(2)=ACOS(SQRT(0.2D0))
  VZeta(3)=ACOS(SQRT((5.0D0-2.0D0*SQRT(5.0D0))/15.0D0))
  VZeta(4)=Pi-VZeta(3)
  VZeta(5)=Pi-VZeta(2)
  VZeta(6)=Pi-VZeta(1)
  VFi=0.2D0*Pi
  DFi=0.4D0*Pi

  !Vértices superior e inferior
  CavVert(1,:) =(/0.0D0,0.0D0,1.0D0/)
  CavVert(32,:)=(/0.0D0,0.0D0,-1.0D0/)

  !Se calculan las posiciones de los vértices
  !k -> vértices
  !i -> capas
  !j -> posición en cada capa
  k=1
  DO i=1,6
    CZeta=COS(VZeta(i))
    SZeta=SIN(VZeta(i))
    DO j=1,5
      k=k+1
      Fi=MOD(i,2)*VFi+(j-1)*DFi
      CavVert(k,:)=(/SZeta*COS(Fi),SZeta*SIN(Fi),CZeta/)
    END DO
  END DO

  !Se definen los índices de los triángulos que forman cada cara
  CavTrian(1,:)= (/ 1, 6, 2/)
  CavTrian(2,:)= (/ 1, 2, 3/)
  CavTrian(3,:)= (/ 1, 3, 4/)
  CavTrian(4,:)= (/ 1, 4, 5/)
  CavTrian(5,:)= (/ 1, 5, 6/)
  CavTrian(6,:)= (/ 7, 2, 6/)
  CavTrian(7,:)= (/ 8, 3, 2/)
  CavTrian(8,:)= (/ 9, 4, 3/)
  CavTrian(9,:)= (/10, 5, 4/)
  CavTrian(10,:)=(/11, 6, 5/)
  CavTrian(11,:)=(/ 8, 2,12/)
  CavTrian(12,:)=(/ 9, 3,13/)
  CavTrian(13,:)=(/10, 4,14/)
  CavTrian(14,:)=(/11, 5,15/)
  CavTrian(15,:)=(/ 7, 6,16/)
  CavTrian(16,:)=(/ 7,12, 2/)
  CavTrian(17,:)=(/ 8,13, 3/)
  CavTrian(18,:)=(/ 9,14, 4/)
  CavTrian(19,:)=(/10,15, 5/)
  CavTrian(20,:)=(/11,16, 6/)
  CavTrian(21,:)=(/ 8,12,18/)
  CavTrian(22,:)=(/ 9,13,19/)
  CavTrian(23,:)=(/10,14,20/)
  CavTrian(24,:)=(/11,15,21/)
  CavTrian(25,:)=(/ 7,16,17/)
  CavTrian(26,:)=(/ 7,17,12/)
  CavTrian(27,:)=(/ 8,18,13/)
  CavTrian(28,:)=(/ 9,19,14/)
  CavTrian(29,:)=(/10,20,15/)
  CavTrian(30,:)=(/11,21,16/)
  CavTrian(31,:)=(/22,12,17/)
  CavTrian(32,:)=(/23,13,18/)
  CavTrian(33,:)=(/24,14,19/)
  CavTrian(34,:)=(/25,15,20/)
  CavTrian(35,:)=(/26,16,21/)
  CavTrian(36,:)=(/22,18,12/)
  CavTrian(37,:)=(/23,19,13/)
  CavTrian(38,:)=(/24,20,14/)
  CavTrian(39,:)=(/25,21,15/)
  CavTrian(40,:)=(/26,17,16/)
  CavTrian(41,:)=(/22,17,27/)
  CavTrian(42,:)=(/23,18,28/)
  CavTrian(43,:)=(/24,19,29/)
  CavTrian(44,:)=(/25,20,30/)
  CavTrian(45,:)=(/26,21,31/)
  CavTrian(46,:)=(/22,28,18/)
  CavTrian(47,:)=(/23,29,19/)
  CavTrian(48,:)=(/24,30,20/)
  CavTrian(49,:)=(/25,31,21/)
  CavTrian(50,:)=(/26,27,17/)
  CavTrian(51,:)=(/22,27,28/)
  CavTrian(52,:)=(/23,28,29/)
  CavTrian(53,:)=(/24,29,30/)
  CavTrian(54,:)=(/25,30,31/)
  CavTrian(55,:)=(/26,31,27/)
  CavTrian(56,:)=(/32,28,27/)
  CavTrian(57,:)=(/32,29,28/)
  CavTrian(58,:)=(/32,30,29/)
  CavTrian(59,:)=(/32,31,30/)
  CavTrian(60,:)=(/32,27,31/)

END SUBROUTINE Pentakisdodecaedro

!-------------------------------------------------------------------------------
! Subdivide un conjunto de triángulos, cada uno en cuatro
!-------------------------------------------------------------------------------
! CopiaVert:  Copia de la matriz de vértices para redimensionarla
! CopiaTrian: Copia de la matriz de triángulos para redimensionarla
! Lados:      Matriz con los índices de los nuevos puntos
! NV:         Número de vértices originales
! NT:         Número de triángulos originales
! n1,n2,n3:   Vértices de un triángulo
! p1,p2,p3:   Nuevos vértices de los triángulos
! i,j,k:      Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Dividir
  USE Utilidades
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: CopiaVert
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: CopiaTrian,Lados
  INTEGER :: NV,NT,i,j,k,n1,n2,n3,p1,p2,p3

  !Numero de vértices y triángulos
  NV=SIZE(CavVert,1)
  NT=SIZE(CavTrian,1)

  !Redimensiona las matrices
  ALLOCATE(CopiaVert(NV,3),CopiaTrian(NT,3),Lados(NV,NV))
  CopiaVert(:,:)=CavVert(:,:)
  CopiaTrian(:,:)=CavTrian(:,:)
  DEALLOCATE(CavVert,CavCent,CavTrian)
  ALLOCATE(CavVert(NV+3*NT/2,3),CavCent(4*NT,3),CavTrian(4*NT,3))
  CavVert(1:NV,:)=CopiaVert(:,:)
  CavTrian(1:NT,:)=CopiaTrian(:,:)
  DEALLOCATE(CopiaVert,CopiaTrian)

  Lados(:,:)=0
  j=NV
  k=NT
  !Para cada triángulo original
  DO i=1,NT
    !Se obtienen los tres vértices
    n1=CavTrian(i,1)
    n2=CavTrian(i,2)
    n3=CavTrian(i,3)
    !En cada uno de los tres lados, se halla el punto medio.
    !Si el punto medio ya ha sido asignado, se recupera su valor.
    p1=Lados(n2,n3)
    IF (p1 == 0) THEN
      j=j+1
      CavVert(j,:)=0.5D0*(CavVert(n2,:)+CavVert(n3,:))
      Lados(n3,n2)=j
      p1=j
    END IF
    p2=Lados(n3,n1)
    IF (p2 == 0) THEN
      j=j+1
      CavVert(j,:)=0.5D0*(CavVert(n3,:)+CavVert(n1,:))
      Lados(n1,n3)=j
      p2=j
    END IF
    p3=Lados(n1,n2)
    IF (p3 == 0) THEN
      j=j+1
      CavVert(j,:)=0.5D0*(CavVert(n1,:)+CavVert(n2,:))
      Lados(n2,n1)=j
      p3=j
    END IF
    !Se generan los nuevos 4 triángulos
    CavTrian(i,:)=(/p1,p2,p3/)
    CavTrian(k+1,:)=(/n1,p3,p2/)
    CavTrian(k+2,:)=(/n2,p1,p3/)
    CavTrian(k+3,:)=(/n3,p2,p1/)
    k=k+3
  END DO

  DEALLOCATE(Lados)

END SUBROUTINE Dividir

END MODULE Cavidad
