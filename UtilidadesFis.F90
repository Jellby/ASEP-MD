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

!-------------------------------------------------------------------------------
! Subrutinas útiles con mayor sentido fisico
!-------------------------------------------------------------------------------
MODULE UtilidadesFis

CONTAINS
!AjusteCargas(PosQ,Pot,QTot,ValQ,Errores,Rest,Dip,PotQ,Energia)
!OrientarMolecula(Mol,Tipo)
!IgualarCargas(Id,Cargas)
!PotencialCargas(Cargas,Puntos,Potencial)
!SuperponerMoleculas(Ref,Mol,Pesos)

!-------------------------------------------------------------------------------
! Ajusta los valores de una serie de cargas para reproducir un potencial dado
! por el método de los multiplicadores de Lagrange (restricción de carga total,
! de dipolo o de energía de interacción)
!-------------------------------------------------------------------------------
! PosQ:       Posiciones de las cargas cuyo valor se quiere ajustar
! Pot:        Posiciones de los puntos donde se conoce el potencial y su valor
! QTot:       Carga total del sistema a ajustar
! ValQ:       Valores de las cargas ajustadas
! Errores:    Errores del ajuste
! Rest:       Factor para restringir el valor absoluto de las cargas
! Dip:        Momento dipolar que deben tener las cargas ajustadas
! PotQ:       Potencial generado por un conjunto de cargas externas en las
!             posiciones de las cargas  a ajustar
! Energia:    Energía de interacción entre las cargas a ajustar y las cargas
!             externas
! Ecuaciones: Coeficientes del sistema de ecuaciones que hay que resolver
! TerminosI:  Términos independientes del sistema de ecuaciones
! Dist:       Vector con las distancias de todas las cargas a un punto dado
! DifPot:     Diferencia entre el potencial inicial y el calculado en cada punto
! Num:        Número de cargas a ajustar
! i,j:        Contadores
!-------------------------------------------------------------------------------
SUBROUTINE AjusteCargas(PosQ,Pot,QTot,ValQ,Errores,Rest,Dip,PotQ,Energia)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: PosQ,Pot
  DOUBLE PRECISION, INTENT(IN) :: QTot
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: ValQ
  DOUBLE PRECISION, DIMENSION(3), INTENT(OUT), OPTIONAL :: Errores
  DOUBLE PRECISION, OPTIONAL :: Rest
  DOUBLE PRECISION, DIMENSION(3), OPTIONAL :: Dip
  DOUBLE PRECISION, DIMENSION(:), OPTIONAL :: PotQ
  DOUBLE PRECISION, OPTIONAL :: Energia

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ecuaciones
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TerminosI,Dist,DifPot
  INTEGER :: Num,i,j

  !Se comprueba que el sistema no es indeterminado
  Num=SIZE(ValQ,1)
  IF (Num > SIZE(Pot,1)) CALL Mensaje('AjusteCargas',36,.TRUE.)

  !Se construyen las ecuaciones:
  !a(i,j)=SUMA[1/(r(ik)*r(jk))]  r(ik) -> distancia entre carga i y punto k
  !b(i)=SUMA[V(k)/r(ik)]         V(k)  -> potencial en el punto k
  i=1
  IF (PRESENT(Dip) .AND. (Num > 3)) i=i+3
  IF (PRESENT(PotQ) .AND. PRESENT(Energia)) i=i+1
  ALLOCATE(Ecuaciones(Num+i,Num+i),TerminosI(Num+i), &
           Dist(Num),DifPot(SIZE(Pot,1)))
  Ecuaciones(:,:)=0.0D0
  TerminosI(:)=0.0D0
  DO j=1,SIZE(Pot,1)
    DO i=1,Num
      Dist(i)=1.0D0/Distancia(PosQ(i,1:3),Pot(j,1:3))
    END DO
    Ecuaciones(1:Num,1:Num)=Ecuaciones(1:Num,1:Num)+ &
                            ProductoTens(Dist(:),Dist(:))
    TerminosI(1:Num)=TerminosI(1:Num)+Pot(j,4)*Dist(1:Num)
  END DO
  !Se añade la restricción para el valor de las cargas
  !a(i,j)=SUMA[1/(r(ik)*r(jk))]+R*delta(i,j)
  IF (PRESENT(Rest)) THEN
    DO i=1,Num
      Ecuaciones(i,i)=Ecuaciones(i,i)+Rest
    END DO
  END IF
  Ecuaciones(Num+1,1:Num)=1.0D0
  Ecuaciones(1:Num,Num+1)=1.0D0
  TerminosI(Num+1)=QTot
  j=Num+1
  IF (PRESENT(Dip) .AND. (Num > 3)) THEN
    DO i=1,3
      Ecuaciones(j+i,1:Num)=PosQ(:,i)
      Ecuaciones(1:Num,j+i)=PosQ(:,i)
    END DO
    TerminosI(j+1:j+3)=Dip(:)
    j=j+3
  END IF
  IF (PRESENT(PotQ) .AND. PRESENT(Energia)) THEN
    IF (.NOT. PRESENT(Energia)) Energia=0.0D0
    Ecuaciones(j+1,1:Num)=PotQ(:)
    Ecuaciones(1:Num,j+1)=PotQ(:)
    TerminosI(j+1)=Energia
  END IF

  !Se resuelve el sistema
  CALL ResolverLU(Ecuaciones,TerminosI)
  ValQ(:)=TerminosI(1:Num)

  !Se calculan los errores
  !1. RErr=SUMA[|(V-V')/V|)]/n           n -> número de puntos en el potencial
  !2. Rms=SQRT[SUMA[(V-V')**2]/n]        V -> potencial real
  !3. RRms=SQRT[SUMA[((V-V')/V)**2]/n]   V'-> potencial calculado
  IF (PRESENT(Errores)) THEN
    Errores(:)=0.0D0
    DifPot(:)=-Pot(:,4)
    DO j=1,SIZE(Pot,1)
      DO i=1,Num
        Dist(i)=1.0D0/Distancia(PosQ(i,1:3),Pot(j,1:3))
        DifPot(j)=DifPot(j)+ValQ(i)*Dist(i)
      END DO
    END DO
    Errores(1)=SUM(ABS(DifPot(:)/Pot(:,4)))/SIZE(Pot,1)
    Errores(2)=SQRT(SUM(DifPot(:)*DifPot(:))/SIZE(Pot,1))
    Errores(3)=SQRT(SUM(DifPot(:)*DifPot(:)/(Pot(:,4)*Pot(:,4)))/SIZE(Pot,1))
  END IF

  DEALLOCATE(Ecuaciones,TerminosI,Dist,DifPot)

END SUBROUTINE AjusteCargas

!-------------------------------------------------------------------------------
! Mueve y rota una molécula para centrarla en su centro de masas y orientarla
! según sus ejes principales de inercia
!-------------------------------------------------------------------------------
! Mol:      Molécula que se va a transformar
! Tipo:     1 -> se ordenan los ejes de inercia. 0 -> no se ordenan.
!           -1 -> sólo se centra la molécula (no se gira)
! Coord:    Matriz con las coordenadas cartesianas de la molécula
! Inercia:  Tensor de inercia de la molécula
! Ejes:     Ejes principales de inercia (vectores propios del tensor)
! Momentos: Momentos de inercia (valores propios del tensor)
! Masa:     Masa total de la molécula
! Aux:      Variable auxiliar
! i,j:      Contadores
!-------------------------------------------------------------------------------
SUBROUTINE OrientarMolecula(Mol,Tipo)
  USE TipoAtomo
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(INOUT) :: Mol
  INTEGER, INTENT(IN) :: Tipo

  DOUBLE PRECISION, DIMENSION(SIZE(Mol,1),3) :: Coord
  DOUBLE PRECISION, DIMENSION(3,3) :: Inercia,Ejes
  DOUBLE PRECISION, DIMENSION(3) :: Momentos
  DOUBLE PRECISION :: Masa,Aux
  INTEGER :: i,j

  !Si la molécula no existe, no se hace nada
  IF (SIZE(Mol,1) == 0) RETURN

  !Se copian las coordenadas en una matriz
  DO i=1,SIZE(Mol,1)
    Coord(i,:)=Mol(i)%pos(:)
  END DO

  !Se centra la molécula en el centro de masas
  Masa=SUM(Mol(:)%m)
  DO i=1,3
    Coord(:,i)=Coord(:,i)-SUM(Coord(:,i)*Mol(:)%m)/Masa
  END DO

  !Si Tipo < 0 sólo se centra la molécula
  IF (Tipo >= 0) THEN

    !Se calcula el tensor de inercia
    Aux=SUM((Coord(:,1)*Coord(:,1)+Coord(:,2)*Coord(:,2)+ &
            Coord(:,3)*Coord(:,3))*Mol(:)%m)
    DO i=1,3
      DO j=i,3
        Inercia(i,j)=-SUM(Coord(:,i)*Coord(:,j)*Mol(:)%m)
        Inercia(j,i)=Inercia(i,j)
      END DO
      Inercia(i,i)=Inercia(i,i)+Aux
    END DO

    !Se diagonaliza el tensor de inercia
    CALL Diagonalizar(Inercia,Ejes,Momentos,Tipo)

    IF (Tipo == 0) THEN
      !Intenta mantener una orientación parecida a la original
      i=SUM(MAXLOC(ABS(Ejes(1,:))))
      IF (i /= 1) CALL Intercambiar(Ejes(:,1),Ejes(:,i))
      i=SUM(MAXLOC(ABS(Ejes(2,2:))))+1
      IF (i /= 2) CALL Intercambiar(Ejes(:,2),Ejes(:,i))
    END IF

    !Mantiene la quiralidad de la molécula
    !(invierte uno de los vectores si el determinante es negativo)
    Aux=DOT_PRODUCT(Ejes(:,1),ProductoVect(Ejes(:,2),Ejes(:,3)))
    IF (Aux < 0.0D0) THEN
      i=SUM(MINLOC(MAXVAL(ABS(Ejes(:,:)),DIM=1)))
      Ejes(:,i)=-Ejes(:,i)
    END IF

    !Se gira la molécula
#ifdef __PORTLAND__
    Coord=TRANSPOSE(MATMUL(TRANSPOSE(Ejes),TRANSPOSE(Coord)))
#else
    Coord=MATMUL(Coord,Ejes)
#endif

  END IF

  !Se vuelven a copiar las coordenadas
  Mol(:)%pos(1)=Coord(:,1)
  Mol(:)%pos(2)=Coord(:,2)
  Mol(:)%pos(3)=Coord(:,3)

END SUBROUTINE

!-------------------------------------------------------------------------------
! Se igualan las cargas (haciendo la media) de los átomos con igual "id"
!-------------------------------------------------------------------------------
! Id:      Matriz con los "id" de los átomos
! Cargas:  Matriz con las cargas de los átomos
! Rep:     Variable para controlar los átomos con igual "id"
! i,j:     Contadores
!-------------------------------------------------------------------------------
SUBROUTINE IgualarCargas(Id,Cargas)
  IMPLICIT NONE
  INTEGER, DIMENSION(:), INTENT(IN) :: Id
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Cargas

  INTEGER :: Rep,i,j

  !Se igualan las cargas de los átomos equivalentes, haciendo la media
  DO i=1,SIZE(Id,1)
    Rep=0
    !Si es equivalente a algún átomo anterior, se igualan las cargas
    DO j=1,i-1
      IF (Id(j) == Id(i)) THEN
        Rep=j
        Cargas(i)=Cargas(Rep)
        EXIT
      END IF
    END DO
    !Si no, se halla la media de todos los equivalentes que vengan después
    IF (Rep == 0) THEN
      Rep=1
      DO j=i+1,SIZE(Id,1)
        IF (Id(j) == Id(i)) THEN
          Cargas(i)=Cargas(i)+Cargas(j)
          Rep=Rep+1
        END IF
      END DO
      Cargas(i)=Cargas(i)/DBLE(Rep)
    END IF
  END DO

END SUBROUTINE IgualarCargas

!-------------------------------------------------------------------------------
! Calcula el potencial electrostático generado por un conjunto de cargas en un
! conjunto de puntos
!-------------------------------------------------------------------------------
! Cargas:    Matriz con las cargas que generan el potencial
! Puntos:    Matriz con los puntos donde se calcula el potencial
! Potencial: Vector con los valores del potencial en los puntos
! Punto:     Localización de cada carga
! Q:         Valor de cada carga
! i,j:       Contadores
!-------------------------------------------------------------------------------
SUBROUTINE PotencialCargas(Cargas,Puntos,Potencial)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Cargas,Puntos
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Potencial

  DOUBLE PRECISION, DIMENSION(3) :: Punto
  DOUBLE PRECISION :: Q
  INTEGER :: i,j

  Potencial(:)=0.0D0
  DO i=1,SIZE(Cargas,1)
    Punto=Cargas(i,1:3)
    Q=Cargas(i,4)
!$OMP PARALLEL DO PRIVATE(j)
    DO j=1,SIZE(Puntos,1)
      Potencial(j)=Potencial(j)+Q/Distancia(Puntos(j,:),Punto(:))
    END DO
!$OMP END PARALLEL DO
  END DO

END SUBROUTINE PotencialCargas

!-------------------------------------------------------------------------------
! Superpone las geometrías de dos moléculas
!   Ver: Acta Chrystallogr. Sec. A 61 (2005), 478
!        J. Comput. Chem. 31 (2010), 1561
!-------------------------------------------------------------------------------
! Ref:             Molécula de referencia
! Mol:             Molécula que se superpone sobre la referencia
! Pesos:           Pesos relativos para considerar sólo parte de las moléculas
! Trans:           Transformación necesaria: cuaternión, desplazamiento
! MRef:            Molécula de referencia centrada
! W:               Pesos reales que se usan en el cálculo
! Coords:          Coordenadas de la molécula que se mueve
! MatrizM,MatrizK: Matrices clave para el cálculo
! Aux:             Matriz 3×3 auxiliar
! Centro1.Centro2: Centros de las dos moléculas
! Vector:          Vector (cuaternión) que da la rotación óptima
! P:               Coeficientes del polinomio característico
! G1,G2:           Productos internos de las dos moléculas
! Valor:           Valor que da la distancia RMSD de las dos estructuras
! ValorAnt:        Valor anterior para calcular la raíz del polinomio
! Preci:           Precisión para el cálculo de la raíz
! i,j:             Contadores
!-------------------------------------------------------------------------------
SUBROUTINE SuperponerMoleculas(Ref,Mol,Pesos,Trans)
  USE TipoAtomo
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(IN) :: Ref
  TYPE(Atomo), DIMENSION(SIZE(Ref,1)), INTENT(INOUT) :: Mol
  DOUBLE PRECISION, DIMENSION(SIZE(Ref,1)), INTENT(IN), OPTIONAL :: Pesos
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT), OPTIONAL :: Trans

  TYPE(Atomo), DIMENSION(SIZE(Ref,1)) :: MRef
  DOUBLE PRECISION, DIMENSION(SIZE(Ref,1)) :: W
  DOUBLE PRECISION, DIMENSION(SIZE(Ref,1),3) :: Coords
  DOUBLE PRECISION, DIMENSION(3,3) :: MatrizM,Aux
  DOUBLE PRECISION, DIMENSION(4,4) :: MatrizK
  DOUBLE PRECISION, DIMENSION(3) :: Centro1,Centro2
  DOUBLE PRECISION, DIMENSION(4) :: Vector
  DOUBLE PRECISION, DIMENSION(5) :: P
  DOUBLE PRECISION :: G1,G2,Valor,ValorAnt
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-6
  INTEGER :: i,j

  !Por defecto, todos los pesos son iguales
  IF (PRESENT(Pesos)) THEN
    W(:)=Pesos(:)
   ELSE
    W(:)=1.0D0
  END IF

  !Se centran las dos moléculas en el origen
  DO i=1,3
    Centro1(i)=SUM(W(:)*Ref(:)%pos(i))/SUM(W)
    MRef(:)%pos(i)=Ref(:)%pos(i)-Centro1(i)
    Centro2(i)=SUM(W(:)*Mol(:)%pos(i))/SUM(W)
    Mol(:)%pos(i)=Mol(:)%pos(i)-Centro2(i)
  END DO

  !Se calcula el producto interno de cada estructura
  G1=SUM(W(:)*(MRef(:)%pos(1)**2+MRef(:)%pos(2)**2+MRef(:)%pos(3)**2))
  G2=SUM(W(:)*(Mol(:)%pos(1)**2+Mol(:)%pos(2)**2+Mol(:)%pos(3)**2))

  !Se construye la matriz M, producto interno de las dos estructuras
  DO i=1,3
    DO j=1,3
      MatrizM(i,j)=SUM(W(:)*MRef(:)%pos(i)*Mol(:)%pos(j))
    END DO
  END DO

  !Se construye y diagonaliza la matriz K
  MatrizK(1,1)=MatrizM(1,1)+MatrizM(2,2)+MatrizM(3,3)
  MatrizK(1,2)=MatrizM(2,3)-MatrizM(3,2)
  MatrizK(1,3)=MatrizM(3,1)-MatrizM(1,3)
  MatrizK(1,4)=MatrizM(1,2)-MatrizM(2,1)
  MatrizK(2,2)=MatrizM(1,1)-MatrizM(2,2)-MatrizM(3,3)
  MatrizK(2,3)=MatrizM(1,2)+MatrizM(2,1)
  MatrizK(2,4)=MatrizM(1,3)+MatrizM(3,1)
  MatrizK(3,3)=MatrizM(2,2)-MatrizM(1,1)-MatrizM(3,3)
  MatrizK(3,4)=MatrizM(2,3)+MatrizM(3,2)
  MatrizK(4,4)=MatrizM(3,3)-MatrizM(1,1)-MatrizM(2,2)
  DO i=2,4
    DO j=1,i-1
      MatrizK(i,j)=MatrizK(j,i)
    END DO
  END DO

  !Se construye el polinomio característico de la matriz K
  P(1)=0.0D0
  !El determinante de una matriz 4×4
  j=1
  DO i=1,4
    Aux(:i-1,:j-1)=MatrizK(:i-1,:j-1)
    Aux(:i-1,j:)=MatrizK(:i-1,j+1:)
    Aux(i:,:j-1)=MatrizK(i+1:,:j-1)
    Aux(i:,j:)=MatrizK(i+1:,j+1:)
    P(1)=P(1)+(-1)**(i+j)*MatrizK(i,j)*Determinante(Aux)
  END DO
  P(2)=-8.0D0*Determinante(MatrizM)
  P(3)=-2.0D0*SUM(MatrizM*MatrizM)

  !Se halla la raíz mayor del polinomio por el método Newton-Raphson
  Valor=0.5D0*(G1+G2)
  ValorAnt=HUGE(ValorAnt)
  DO WHILE (ABS(Valor-ValorAnt) > Preci)
    ValorAnt=Valor
    !P(4) es el polinomio, P(5) la derivada
    P(4)=P(1)+P(2)*Valor+P(3)*Valor**2+Valor**4
    P(5)=P(2)+2.0D0*P(3)*Valor+4.0D0*Valor**3
    Valor=Valor-P(4)/P(5)
  END DO

  !Se obtiene el vector propio correspondiente
  MatrizK(1,1)=MatrizK(1,1)-Valor
  MatrizK(2,2)=MatrizK(2,2)-Valor
  MatrizK(3,3)=MatrizK(3,3)-Valor
  MatrizK(4,4)=MatrizK(4,4)-Valor
  !... como columna no nula de la matriz de adjuntos
  DO j=1,4
    DO i=1,4
      Aux(:i-1,:j-1)=MatrizK(:i-1,:j-1)
      Aux(:i-1,j:)=MatrizK(:i-1,j+1:)
      Aux(i:,:j-1)=MatrizK(i+1:,:j-1)
      Aux(i:,j:)=MatrizK(i+1:,j+1:)
      Vector(i)=(-1)**(i+j)*Determinante(Aux)
    END DO
    IF (Norma(Vector) > Preci) EXIT
  END DO
  !Parece que la convención de signos es distinta
  Vector(1)=-Vector(1)
  Vector=Normalizar(Vector)

  !Finalmente se rota y desplaza la segunda molécula
  DO i=1,SIZE(Ref,1)
    Coords(i,:)=Mol(i)%pos(:)
  END DO
  Coords(:,:)=RotarCuaternion(Coords(:,:),Vector)
  DO i=1,SIZE(Ref,1)
    Mol(i)%pos(:)=Coords(i,:)+Centro1(:)
  END DO

  !Se calcula la distancia RMSD
  Valor=SQRT(ABS(G1+G2-2.0D0*Valor)/SUM(W))

  !Se devuelve la transformación
  IF (PRESENT(Trans)) THEN
    !Primero el cuaternión
    Trans(1:4)=Vector(:)
    !Después el desplazamiento una vez rotado
    Aux(1,:)=Centro2(:)
    Aux(:,:)=RotarCuaternion(Aux(:,:),Vector)
    Trans(5:7)=Centro1(:)-Aux(1,:)
  END IF

END SUBROUTINE SuperponerMoleculas

END MODULE UtilidadesFis
