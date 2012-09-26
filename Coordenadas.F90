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

MODULE Coordenadas

#ifndef LLL
#define LLL 256
#endif

!-------------------------------------------------------------------------------
! TipoCoord:   Tipo de coordenadas de trabajo para utilizar en el cálculo
!              0 -> cartesianas, 1 -> cartesianas ponderadas, 2-> internas
! DefCoord:    Definición de las coordenadas internas (0 para cartesianas)
! Combinacion: Definición de las combinaciones lineales de coordenadas
! Cong:        Coordenadas congeladas, que no se optimizan
! Masa:        Raíz de la masa, para coordenadas ponderadas
! Geometria:   Valores de las coordenadas de trabajo
! Gradiente:   Vector gradiente en coordenadas de trabajo
! Hessiana:    Matriz hessiana en coordenadas de trabajo
! MatrizB:     Matriz B de conversión de cartesianas a internas
! InvB:        Pseudoinversa de la matriz B
! Der:         Matriz de segundas derivadas de las coordenadas internas
!-------------------------------------------------------------------------------
INTEGER :: TipoCoord
INTEGER, DIMENSION(:,:), ALLOCATABLE :: DefCoord
INTEGER, DIMENSION(:), ALLOCATABLE :: Cong
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Masa,Geometria,Gradiente
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MatrizB,InvB,Hessiana, &
                                                 Combinacion
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: Der

CONTAINS
!Cartesianas(Mol)
!CambiarCartesianas(Mol,Coord)
!DefinirCoordenadas(Mol)
!CoordenadasInternas(Mol)
!ConvertirCoordenadas(Coord,Deriv)
!ModificarCoordenadas(Fich)
!ConvertirGradHess(Tipo,Grad,Hess)
!ConvertirIncremento(Paso,CoordCart,PasoCart)
!Proyectar(Desp,Grad,Hess)

!-------------------------------------------------------------------------------
! Función que devuelve el vector de coordenadas cartesianas de una molécula
!-------------------------------------------------------------------------------
! Mol:     Molécula cuyas coordenadas se piden
! i:       Contador
!-------------------------------------------------------------------------------
FUNCTION Cartesianas(Mol)
  USE TipoAtomo
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(IN) :: Mol
  DOUBLE PRECISION, DIMENSION(3*SIZE(Mol,1)) :: Cartesianas

  INTEGER :: i

  DO i=1,SIZE(Mol,1)
    Cartesianas(3*i-2:3*i)=Mol(i)%pos(:)
  END DO

END FUNCTION Cartesianas

!-------------------------------------------------------------------------------
! Cambia las coordenadas de una molécula
!-------------------------------------------------------------------------------
! Mol:     Molécula cuyas coordenadas se cambian
! Coord:   Vector de coordenadas cartesianas
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE CambiarCartesianas(Mol,Coord)
  USE TipoAtomo
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(INOUT) :: Mol
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Coord

  INTEGER :: i

  DO i=1,SIZE(Mol,1)
    Mol(i)%pos(:)=Coord(3*i-2:3*i)
  END DO

END SUBROUTINE CambiarCartesianas

!-------------------------------------------------------------------------------
! Subrutina que asigna las coordenadas de trabajo (internas o cartesianas)
!-------------------------------------------------------------------------------
! Mol:     Molécula cuyas coordenadas se van a definir
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE DefinirCoordenadas(Mol)
  USE TipoAtomo
  USE Unidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(IN) :: Mol

  INTEGER :: i

  IF (ALLOCATED(Masa)) DEALLOCATE(Masa)
  ALLOCATE(Masa(3*SIZE(Mol,1)))
  IF (TipoCoord == 0) THEN
    !Coordenadas cartesianas normales
    Masa(:) = 1.0D0
   ELSE
    !Raíz cuadrada de la masa atómica para las coordenadas ponderadas
    DO i=1,SIZE(Mol,1)
      Masa(3*i-2:3*i)=SQRT(Mol(i)%m/AmuAtomica)
      !IF (Mol(i)%m == 0.0D0) Masa(3*i-2:3*i) = 1.0D0
    END DO
  END IF

  SELECT CASE (TipoCoord)
   CASE (0,1)
    !Se hace DefCoord igual a cero: coordenadas cartesianas
    IF (ALLOCATED(DefCoord)) DEALLOCATE(DefCoord)
    ALLOCATE(DefCoord(3*SIZE(Mol,1),4))
    DefCoord(:,:)=0
   CASE (2)
    !Se definen las coordenadas internas
    CALL CoordenadasInternas(Mol)
  END SELECT

  IF (ALLOCATED(Geometria)) DEALLOCATE(Geometria)
  IF (ALLOCATED(Cong)) DEALLOCATE(Cong)
  ALLOCATE(Geometria(SIZE(DefCoord,1)),Cong(SIZE(DefCoord,1)))
  Cong=0

END SUBROUTINE DefinirCoordenadas

!-------------------------------------------------------------------------------
! En esta subrutina se definen las coordenadas internas según los criterios de:
!  V. Bakken, T. Helgaker; J. Chem. Phys. 117 (2002), 9160-9174
!-------------------------------------------------------------------------------
! Mol:         Molécula cuyas coordenadas se van a definir
! Dists:       Matriz de distancias interatómicas
! Enlaces:     Matriz con el tipo de enlace entre cada dos átomos
! Frag:        Fragmento al que pertenece cada átomo
! Dist:        Distancia para hacer comparaciones
! Ang:         Ángulo para hacer comparaciones
! FactEnl,FactH,LimAngH,FactFrag,LimFrag,LimAng: Criterios para las coordenadas
! Num:         Número de átomos en la molécula
! At1,At2:     Átomos más cercanos entre dos fragmentos
! UTmp:        Unidad temporal
! i,j,k,l,m,n: Contadores
!-------------------------------------------------------------------------------
SUBROUTINE CoordenadasInternas(Mol)
  USE TipoAtomo
  USE Unidades
  USE Utilidades
  IMPLICIT NONE
  TYPE(Atomo), DIMENSION(:), INTENT(IN) :: Mol

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Dists
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Enlaces
  INTEGER, DIMENSION(:), ALLOCATABLE :: Frag
  DOUBLE PRECISION :: Dist,Ang
  DOUBLE PRECISION, PARAMETER :: FactEnl=1.3D0,FactH=0.9D0,LimAngH=Pi/2.0D0, &
                                 FactFrag=1.3D0,LimFrag=2.0D0*AngstromAtomica, &
                                 LimAng=0.99D0*Pi
  INTEGER :: Num,At1,At2,UTmp,i,j,k,l,m,n

  Num=SIZE(Mol,1)
  ALLOCATE(Dists(Num,Num),Enlaces(Num,Num),Frag(Num))

  !Se calculan todas las distancias y se localizan los enlaces
  Enlaces(:,:)=0
  DO i=1,Num
    DO j=i+1,Num
      Dists(i,j)=Distancia(Mol(i)%pos(:),Mol(j)%pos(:))
      Dists(j,i)=Dists(i,j)
      Dist=FactEnl*(RadiosCov(Mol(i)%z)+RadiosCov(Mol(j)%z))
      IF (Dists(i,j) > Dist) CYCLE
      !Enlaces normales
      Enlaces(i,j)=5
      Enlaces(j,i)=5
    END DO
  END DO
  
  !Se añaden los enlaces de hidrógeno
  DO i=1,Num
    !Sólo los H
    IF (Mol(i)%z /= 1) CYCLE
    DO j=1,Num
      !Todos los X unidos a H
      IF (Enlaces(i,j) /= 5) CYCLE
      IF ((Mol(j)%z /= 7) .AND. (Mol(j)%z /= 8) .AND. &
          (Mol(j)%z /= 9) .AND. (Mol(j)%z /=15) .AND. &
          (Mol(j)%z /=16) .AND. (Mol(j)%z /=17)) CYCLE
      DO k=1,Num
        !Todos los Y no unidos a H
        IF (Enlaces(i,k) /= 0) CYCLE
        IF ((k == i) .OR. (k == j)) CYCLE
        IF ((Mol(k)%z /= 7) .AND. (Mol(k)%z /= 8) .AND. &
            (Mol(k)%z /= 9) .AND. (Mol(k)%z /=15) .AND. &
            (Mol(k)%z /=16) .AND. (Mol(k)%z /=17)) CYCLE
        Dist=FactH*(RadiosVdW(Mol(i)%z)+RadiosVdW(Mol(k)%z))
        Ang=Angulo(Mol(j)%pos(:),Mol(i)%pos(:),Mol(k)%pos(:))
        IF ((Dists(i,k) > Dist) .OR. (Ang <= LimAngH)) CYCLE
        !Enlaces de hidrógeno
        Enlaces(i,k)=4
        Enlaces(k,i)=4
      END DO
    END DO
  END DO

  !Se detectan posibles fragmentos y se añaden los enlaces
  Frag(:)=0
  DO i=1,Num
    DO j=1,i-1
      IF (Enlaces(i,j) == 5) THEN
        IF ((Frag(i) > Frag(j)) .OR. (Frag(i) == 0)) THEN
          WHERE (Frag(1:i) == Frag(i)) Frag(1:i)=Frag(j)
        ELSE IF (Frag(i) < Frag(j)) THEN
          WHERE (Frag(1:i) == Frag(j)) Frag(1:i)=Frag(i)
        END IF
      END IF
    END DO
    IF (Frag(i) == 0) Frag(i)=MAXVAL(Frag(:))+1
  END DO
  At1=0
  At2=0
  DO i=1,MAXVAL(Frag(:))
    DO j=i+1,MAXVAL(Frag(:))
      !Mínima distancia entre cada dos fragmentos
      Dist=100.0D0
      DO k=1,Num
        IF (Frag(k) /= i) CYCLE
        DO l=1,Num
          IF (Frag(l) /= j) CYCLE
          IF (Dists(k,l) < Dist) THEN
            Dist=Dists(k,l)
            At1=k
            At2=l
          END IF
        END DO
      END DO
      !Enlaces interfragmento
      Enlaces(At1,At2)=3
      Enlaces(At2,At1)=3
      Dist=MAX(FactFrag*Dist,LimFrag)
      DO k=1,Num
        IF (Frag(k) /= i) CYCLE
        DO l=1,Num
          IF (Frag(l) /= j) CYCLE
          IF (Dists(k,l) > Dist) CYCLE
          IF (Enlaces(k,l) == 0) THEN
            !Enlaces auxiliares interfragmento
            Enlaces(k,l)=2
            Enlaces(l,k)=2
          END IF
        END DO
      END DO
    END DO
  END DO

  !Se cuentan los ángulos y diedros
  UTmp=NuevaUnidad()
  OPEN(UTmp,STATUS='SCRATCH',ACTION='READWRITE')
  m=0
  n=0
  DO i=1,Num
    DO j=i+1,Num
      DO k=1,Num
        IF ((k == i) .OR. (k == j)) CYCLE
        IF (Enlaces(k,i) < 3) CYCLE
        DO l=1,Num
          IF ((l == i) .OR. (l == j) .OR. (l == k)) CYCLE
          IF (Enlaces(l,j) < 3) CYCLE
          IF (Enlaces(l,k) < 3) CYCLE
          IF ((Angulo(Mol(i)%pos,Mol(k)%pos,Mol(l)%pos) > LimAng) .OR. &
              (Angulo(Mol(k)%pos,Mol(l)%pos,Mol(j)%pos) > LimAng)) CYCLE
          n=n+1
          WRITE(UTmp,*) 2
          WRITE(UTmp,*) i,k,l,j
        END DO
        IF (Enlaces(k,j) < 3) CYCLE
        m=m+1
        WRITE(UTmp,*) 1
        WRITE(UTmp,*) i,k,j
      END DO
    END DO
  END DO

  !Se asignan todas las coordenadas
  i=COUNT(Enlaces(:,:) /= 0)/2
  IF (ALLOCATED(DefCoord)) DEALLOCATE(DefCoord)
  ALLOCATE(DefCoord(6+i+m+n,4))
  DefCoord(:,:)=0
  DefCoord(1,4)=1
  DefCoord(2,4)=2
  DefCoord(3,4)=3
  DefCoord(4,3)=1
  DefCoord(5,3)=2
  DefCoord(6,3)=3
  k=7
  !Enlaces
  DO i=1,Num
    DO j=i+1,Num
      IF (Enlaces(i,j) == 0) CYCLE
      DefCoord(k,1:2)=(/i,j/)
      k=k+1
    END DO
  END DO
  REWIND(UTmp)
  !Ángulos
  DO i=1,m+n
    READ(UTmp,*) j
    IF (j == 1) THEN
      READ(UTmp,*) DefCoord(k,1:3)
      k=k+1
     ELSE
      READ(UTmp,*)
    END IF
  END DO
  REWIND(UTmp)
  !Diedros
  DO i=1,m+n
    READ(UTmp,*) j
    IF (j == 2) THEN
      READ(UTmp,*) DefCoord(k,1:4)
      k=k+1
     ELSE
      READ(UTmp,*)
    END IF
  END DO
  CLOSE(UTmp)

  DEALLOCATE(Dists,Enlaces,Frag)

END SUBROUTINE CoordenadasInternas

!-------------------------------------------------------------------------------
! Convierte las coordenadas cartesianas a coordenadas internas
! También crea las matrices de conversion para el gradiente y la hessiana
!  V. Bakken, T. Helgaker; J. Chem. Phys. 117 (2002) 9160-9174 (corregido)
!-------------------------------------------------------------------------------
! Coord:                   Coordenadas cartesianas
! Deriv:                   Orden de derivadas a calcular
! Num1:                    Número de coordenadas internas
! Num2:                    Número de coordenadas cartesianas
! Vec1-Vec4,Vec13,Vec23:   Vectores auxiliares
! Dist1,Dist2,Dist3:       Distancias auxiliares
! Cos1,Cos2:               Cosenos auxiliares
! Sen1,Sen2:               Senos auxiliares (a veces al cuadrado)
! Uni:                     Matriz unidad 3x3
! Aux1 - Aux8:             Matrices auxiliares 3x3
! Preci:                   Precisión para comparar con 1
! MatB,MatAux:             Matrices auxiliares para la inversión
! Diag:                    Valores singulares de B
! a,b,c,d,a1,a2,...,d1,d2: Variables auxiliares
! i,j,k,l                  Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ConvertirCoordenadas(Coord,Deriv)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Coord
  INTEGER, INTENT(IN) :: Deriv

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MatB,MatAux
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Diag
  DOUBLE PRECISION, DIMENSION(3,3) :: Aux1,Aux2,Aux3,Aux4,Aux5,Aux6,Aux7,Aux8, &
                                      Uni
  DOUBLE PRECISION, DIMENSION(3) :: Vec1,Vec2,Vec3,Vec4,Vec13,Vec23
  DOUBLE PRECISION :: Dist1,Dist2,Dist3,Cos1,Cos2,Sen1,Sen2
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D0-1.0D-4
  INTEGER :: Num1,Num2,i,j,k,l,a,b,c,d,a1,a2,b1,b2,c1,c2,d1,d2

  Num1=SIZE(DefCoord,1)
  Num2=SIZE(Coord,1)
  IF (ALLOCATED(MatrizB)) DEALLOCATE(MatrizB)
  IF (ALLOCATED(InvB)) DEALLOCATE(InvB)
  IF (ALLOCATED(Der)) DEALLOCATE(Der)
  ALLOCATE(MatrizB(Num1,Num2),InvB(Num2,Num1),Der(Num1,Num2,Num2))
  MatrizB(:,:)=0.0D0
  Der(:,:,:)=0.0D0
  Uni(:,:)=0.0D0
  DO i=1,3
    Uni(i,i)=1.0D0
  END DO

  !Se calculan las coordenadas internas
  ! y sus primeras derivadas (matriz B)
  ! y sus segundas derivadas (matriz B'=Der)
  DO i=1,Num1
    a=DefCoord(i,1); a1=3*a-2; a2=3*a
    b=DefCoord(i,2); b1=3*b-2; b2=3*b
    c=DefCoord(i,3); c1=3*c-2; c2=3*c
    d=DefCoord(i,4); d1=3*d-2; d2=3*d
    IF (a == 0) THEN
      IF (d /= 0) THEN
        !Traslación (centro geométrico, funciona mejor que el centro de masas)
        Geometria(i)=0.0D0
        DO j=1,SIZE(Coord,1),3
          Geometria(i)=Geometria(i)+Coord(j+d-1)
          MatrizB(i,j+d-1)=1.0D0
        END DO
        Geometria(i)=Geometria(i)*3.0D0/Num2
        MatrizB(i,:)=MatrizB(i,:)*3.0D0/Num2
        !Geometria(i)=Geometria(i)*3.0D0/SUM(Masa(:)*Masa(:))
        !MatrizB(i,:)=MatrizB(i,:)*3.0D0/SUM(Masa(:)*Masa(:))
       ELSE IF (c /= 0) THEN
        !Rotación (alrededor del centro geometrico)
        k=MOD(c,3)+1
        l=MOD(c+1,3)+1
        Vec1(:)=0.0D0
        DO j=1,SIZE(Coord,1),3
          Vec1(k)=Vec1(k)+Coord(j+k-1)
          Vec1(l)=Vec1(l)+Coord(j+l-1)
        END DO
        Vec1(:)=Vec1(:)*3.0D0/Num2
        DO j=1,SIZE(Coord,1),3
          Vec2(k)=Coord(j+k-1)-Vec1(k)
          Vec2(l)=Coord(j+l-1)-Vec1(l)
          Vec2(c)=Vec2(k)**2+Vec2(l)**2
          IF (Vec2(c) > 1.0D-4) THEN
            MatrizB(i,j+k-1)=Vec2(l)/Vec2(c)
            MatrizB(i,j+l-1)=-Vec2(k)/Vec2(c)
            IF (Deriv > 1) THEN
              Der(i,j+k-1,j+k-1)=-2.0D0*Vec2(k)*Vec2(l)/(Vec2(c)**2)
              Der(i,j+l-1,j+l-1)=2.0D0*Vec2(k)*Vec2(l)/(Vec2(c)**2)
              Der(i,j+k-1,j+l-1)=(Vec2(k)**2-Vec2(l)**2)/(Vec2(c)**2)
              Der(i,j+l-1,j+k-1)=Der(i,j+k-1,j+l-1)
            END IF
          END IF
        END DO
        Geometria(i)=0.0D0
        !MatrizB(i,:)=MatrizB(i,:)*3.0D0/Num2
        !Der(i,:,:)=Der(i,:,:)*3.0D0/Num2
        !MatrizB(i,:)=MatrizB(i,:)*(1.0D0-3.0D0/Num2-1.0D0)
        !Der(i,:,:)=Der(i,:,:)*(1.0D0-3.0D0/Num2-1.0D0)
       ELSE
        !Coordenadas cartesianas ponderadas
        Geometria(i)=Coord(i)*Masa(i)
        MatrizB(i,i)=Masa(i)
      END IF
    ELSE IF (a < 0) THEN
      !Combinación lineal de coordenadas internas
      a=-a
      Geometria(i)=0.0D0
      MatrizB(i,:)=0.0D0
      IF (Deriv > 1) Der(i,:,:)=0.0D0
      DO j=1,SIZE(Combinacion,2)
        Geometria(i)=Geometria(i)+Combinacion(a,j)*Geometria(j)
        MatrizB(i,:)=MatrizB(i,:)+Combinacion(a,j)*MatrizB(j,:)
        IF (Deriv > 1) Der(i,:,:)=Der(i,:,:)+Combinacion(a,j)*Der(j,:,:)
      END DO
    ELSE IF (b == 0) THEN
      !Coordenadas cartesianas de un átomo
      Geometria(i)=Coord(a1+d-1)
      MatrizB(i,a1+d-1)=1.0D0
    ELSE IF (c == 0) THEN
      !Conversión de distancias
      Dist1=Distancia(Coord(a1:a2),Coord(b1:b2))
      Geometria(i)=Dist1
      !Primeras derivadas de las distancias
      Vec1(:)=(Coord(a1:a2)-Coord(b1:b2))/Dist1
      MatrizB(i,a1:a2)=Vec1(:)
      MatrizB(i,b1:b2)=-Vec1(:)
      !Segundas derivadas de las distancias
      IF (Deriv > 1) THEN
        Aux1(:,:)=-(ProductoTens(Vec1(:),Vec1(:))-Uni(:,:))/Dist1
        Der(i,a1:a2,a1:a2)=Aux1(:,:)
        Der(i,b1:b2,b1:b2)=Aux1(:,:)
        Der(i,a1:a2,b1:b2)=-Aux1(:,:)
        Der(i,b1:b2,a1:a2)=-Aux1(:,:)
      END IF
    ELSE IF (d == 0) THEN
      !Conversión de ángulos
      Geometria(i)=Angulo(Coord(a1:a2),Coord(b1:b2),Coord(c1:c2))
      !Primeras derivadas de los ángulos
      Vec1(:)=Coord(a1:a2)-Coord(b1:b2)
      Dist1=SQRT(DOT_PRODUCT(Vec1(:),Vec1(:)))
      Vec1(:)=Vec1(:)/Dist1
      Vec2(:)=Coord(c1:c2)-Coord(b1:b2)
      Dist2=SQRT(DOT_PRODUCT(Vec2(:),Vec2(:)))
      Vec2(:)=Vec2(:)/Dist2
      IF (ABS(DOT_PRODUCT(Vec1(:),Vec2(:))) < Preci) THEN
        Vec3(:)=ProductoVect(Vec1(:),Vec2(:))
       ELSE
        Vec3(:)=(/1.0D0,-1.0D0,1.0D0/)
        Dist3=SQRT(3.0D0)
        IF ((ABS(DOT_PRODUCT(Vec1(:),Vec3(:)))/Dist3 < Preci) .AND. &
            (ABS(DOT_PRODUCT(Vec2(:),Vec3(:)))/Dist3 < Preci)) THEN
          Vec3(:)=ProductoVect(Vec1(:),Vec3(:))
         ELSE
          Vec3(:)=(/-1.0D0,1.0D0,1.0D0/)
          Vec3(:)=ProductoVect(Vec1(:),Vec3(:))
        END IF
      END IF
      Dist3=SQRT(DOT_PRODUCT(Vec3(:),Vec3(:)))
      Vec3(:)=Vec3(:)/Dist3
      Vec13(:)=ProductoVect(Vec1(:),Vec3(:))/Dist1
      Vec23(:)=ProductoVect(Vec3(:),Vec2(:))/Dist2
      MatrizB(i,a1:a2)=Vec13(:)
      MatrizB(i,b1:b2)=-Vec13(:)-Vec23(:)
      MatrizB(i,c1:c2)=Vec23(:)
      !Segundas derivadas de los ángulos
      IF (Deriv > 1) THEN
        Cos1=DOT_PRODUCT(Vec1(:),Vec2(:))
        Sen1=SQRT(1.0D0-Cos1*Cos1)
        Aux1(:,:)=(ProductoTens(Vec1(:),Vec2(:))+ProductoTens(Vec2(:),Vec1(:))-&
                  3.0D0*Cos1*ProductoTens(Vec1(:),Vec1(:))+Cos1*Uni(:,:))/ &
                  (Dist1*Dist1*Sen1)
        Aux2(:,:)=(ProductoTens(Vec2(:),Vec1(:))+ProductoTens(Vec1(:),Vec2(:))-&
                  3.0D0*Cos1*ProductoTens(Vec2(:),Vec2(:))+Cos1*Uni(:,:))/ &
                  (Dist2*Dist2*Sen1)
        Aux3(:,:)=(ProductoTens(Vec1(:),Vec1(:))+ProductoTens(Vec2(:),Vec2(:))-&
                  Cos1*ProductoTens(Vec1(:),Vec2(:))-Uni(:,:))/ &
                  (Dist1*Dist2*Sen1)
        Aux4(:,:)=(ProductoTens(Vec2(:),Vec2(:))+ProductoTens(Vec1(:),Vec1(:))-&
                  Cos1*ProductoTens(Vec2(:),Vec1(:))-Uni(:,:))/ &
                  (Dist1*Dist2*Sen1)
        Cos1=Cos1/Sen1
        Der(i,a1:a2,a1:a2)=Aux1(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,a1:a2),MatrizB(i,a1:a2))
        Der(i,b1:b2,b1:b2)=Aux1(:,:)+Aux2(:,:)+Aux3(:,:)+Aux4(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,b1:b2),MatrizB(i,b1:b2))
        Der(i,c1:c2,c1:c2)=Aux2(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,c1:c2),MatrizB(i,c1:c2))
        Der(i,a1:a2,b1:b2)=-Aux1(:,:)-Aux3(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,a1:a2),MatrizB(i,b1:b2))
        Der(i,a1:a2,c1:c2)=Aux3(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,a1:a2),MatrizB(i,c1:c2))
        Der(i,b1:b2,c1:c2)=-Aux2(:,:)-Aux3(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,b1:b2),MatrizB(i,c1:c2))
        Der(i,b1:b2,a1:a2)=-Aux1(:,:)-Aux4(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,b1:b2),MatrizB(i,a1:a2))
        Der(i,c1:c2,a1:a2)=Aux4(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,c1:c2),MatrizB(i,a1:a2))
        Der(i,c1:c2,b1:b2)=-Aux2(:,:)-Aux4(:,:)- &
                           Cos1*ProductoTens(MatrizB(i,c1:c2),MatrizB(i,b1:b2))
      END IF
     ELSE
      !Conversión de diedros
      Geometria(i)=Diedro(Coord(a1:a2),Coord(b1:b2),Coord(c1:c2),Coord(d1:d2))
      !Primeras derivadas de los diedros
      Vec1(:)=Coord(a1:a2)-Coord(b1:b2)
      Dist1=SQRT(DOT_PRODUCT(Vec1(:),Vec1(:)))
      Vec1(:)=Vec1(:)/Dist1
      Vec2(:)=Coord(d1:d2)-Coord(c1:c2)
      Dist2=SQRT(DOT_PRODUCT(Vec2(:),Vec2(:)))
      Vec2(:)=Vec2(:)/Dist2
      Vec3(:)=Coord(c1:c2)-Coord(b1:b2)
      Dist3=SQRT(DOT_PRODUCT(Vec3(:),Vec3(:)))
      Vec3(:)=Vec3(:)/Dist3
      Cos1=DOT_PRODUCT(Vec1(:),Vec3(:))
      Cos2=-DOT_PRODUCT(Vec2(:),Vec3(:))
      Sen1=1.0D0-Cos1*Cos1
      Sen2=1.0D0-Cos2*Cos2
      Vec13(:)=ProductoVect(Vec1(:),Vec3(:))/Sen1
      Vec23(:)=ProductoVect(Vec2(:),Vec3(:))/Sen2
      MatrizB(i,a1:a2)=Vec13(:)/Dist1
      !Error en el signo    --------------------------v
      !MatrizB(i,b1:b2)=-Vec13(:)/Dist1+(Cos1*Vec13(:)-Cos2*Vec23(:))/Dist3
      !MatrizB(i,c1:c2)=Vec23(:)/Dist2-(Cos1*Vec13(:)-Cos2*Vec23(:))/Dist3
      MatrizB(i,b1:b2)=-Vec13(:)/Dist1+(Cos1*Vec13(:)+Cos2*Vec23(:))/Dist3
      MatrizB(i,c1:c2)=Vec23(:)/Dist2-(Cos1*Vec13(:)+Cos2*Vec23(:))/Dist3
      MatrizB(i,d1:d2)=-Vec23(:)/Dist2
      !Segundas derivadas de los diedros
      IF (Deriv > 1) THEN
        Vec13(:)=ProductoVect(Vec1(:),Vec3(:))
        Vec23(:)=ProductoVect(Vec2(:),Vec3(:))
        Vec4(:)=Cos1*Vec3(:)-Vec1(:)
        Aux1(:,:)=ProductoTens(Vec13(:),Vec4(:))+ProductoTens(Vec4(:),Vec13(:))
        Aux1(:,:)=Aux1(:,:)/(Dist1*Dist1*Sen1*Sen1)
        !Vec4(:)=Cos2*Vec3(:)-Vec2(:)
        Vec4(:)=Cos2*Vec3(:)+Vec2(:)
        Aux2(:,:)=ProductoTens(Vec23(:),Vec4(:))+ProductoTens(Vec4(:),Vec23(:))
        Aux2(:,:)=Aux2(:,:)/(Dist2*Dist2*Sen2*Sen2)
        Vec4(:)=Vec3(:)-2.0D0*Cos1*Vec1(:)+Cos1*Cos1*Vec3(:)
        Aux3(:,:)=ProductoTens(Vec13(:),Vec4(:))+ProductoTens(Vec4(:),Vec13(:))
        Aux3(:,:)=Aux3(:,:)/(2.0D0*Dist1*Dist3*Sen1*Sen1)
        Vec4(:)=Vec3(:)+2.0D0*Cos2*Vec2(:)+Cos2*Cos2*Vec3(:)
        Aux4(:,:)=ProductoTens(Vec23(:),Vec4(:))+ProductoTens(Vec4(:),Vec23(:))
        Aux4(:,:)=Aux4(:,:)/(2.0D0*Dist2*Dist3*Sen2*Sen2)
        Vec4(:)=Vec1(:)+Cos1*Cos1*Vec1(:)-3.0D0*Cos1*Vec3(:)+ &
                Cos1*Cos1*Cos1*Vec3(:)
        Aux5(:,:)=ProductoTens(Vec13(:),Vec4(:))+ProductoTens(Vec4(:),Vec13(:))
        Aux5(:,:)=Aux5(:,:)/(2.0D0*Dist3*Dist3*Sen1*Sen1)
        Vec4(:)=Vec2(:)+Cos2*Cos2*Vec2(:)+3.0D0*Cos2*Vec3(:)- &
                Cos2*Cos2*Cos2*Vec3(:)
        Aux6(:,:)=ProductoTens(Vec23(:),Vec4(:))+ProductoTens(Vec4(:),Vec23(:))
        Aux6(:,:)=Aux6(:,:)/(2.0D0*Dist3*Dist3*Sen2*Sen2)
        !Vec4(:)=(Cos1*Vec3(:)-Vec1(:))/(Dist1*Dist3*SQRT(Sen1))
        Vec4(:)=(Cos1*Vec3(:)-Vec1(:))/(Dist1*Dist3*Sen1)
        Aux7(:,:)=0.0D0
        Aux7(1,2)=-0.5D0*Vec4(3); Aux7(2,1)=-Aux7(1,2)
        Aux7(1,3)= 0.5D0*Vec4(2); Aux7(3,1)=-Aux7(1,3)
        Aux7(2,3)=-0.5D0*Vec4(1); Aux7(3,2)=-Aux7(2,3)
        !Vec4(:)=(Cos2*Vec3(:)-Vec2(:))/(Dist2*Dist3*SQRT(Sen2))
        Vec4(:)=(Cos2*Vec3(:)+Vec2(:))/(Dist2*Dist3*Sen2)
        Aux8(:,:)=0.0D0
        Aux8(1,2)=-0.5D0*Vec4(3); Aux8(2,1)=-Aux8(1,2)
        Aux8(1,3)= 0.5D0*Vec4(2); Aux8(3,1)=-Aux8(1,3)
        Aux8(2,3)=-0.5D0*Vec4(1); Aux8(3,2)=-Aux8(2,3)
        Der(i,a1:a2,a1:a2)=Aux1(:,:)
        Der(i,b1:b2,b1:b2)=Aux1(:,:)-2.0D0*Aux3(:,:)-Aux5(:,:)+Aux6(:,:)
        Der(i,c1:c2,c1:c2)=Aux2(:,:)-2.0D0*Aux4(:,:)-Aux5(:,:)+Aux6(:,:)
        Der(i,d1:d2,d1:d2)=Aux2(:,:)
        Der(i,a1:a2,b1:b2)=-Aux1(:,:)+Aux3(:,:)+Aux7(:,:)
        Der(i,a1:a2,c1:c2)=-Aux3(:,:)-Aux7(:,:)
        Der(i,b1:b2,c1:c2)=Aux3(:,:)+Aux4(:,:)+ &
                           Aux5(:,:)-Aux6(:,:)+Aux7(:,:)+Aux8(:,:)
        !Der(i,b1:b2,d1:d2)=-Aux4(:,:)+Aux8(:,:)
        Der(i,b1:b2,d1:d2)=-Aux4(:,:)-Aux8(:,:)
        !Der(i,c1:c2,d1:d2)=-Aux2(:,:)+Aux4(:,:)-Aux8(:,:)
        Der(i,c1:c2,d1:d2)=-Aux2(:,:)+Aux4(:,:)+Aux8(:,:)
        Der(i,b1:b2,a1:a2)=TRANSPOSE(Der(i,a1:a2,b1:b2))
        Der(i,c1:c2,a1:a2)=TRANSPOSE(Der(i,a1:a2,c1:c2))
        Der(i,c1:c2,b1:b2)=TRANSPOSE(Der(i,b1:b2,c1:c2))
        Der(i,d1:d2,b1:b2)=TRANSPOSE(Der(i,b1:b2,d1:d2))
        Der(i,d1:d2,c1:c2)=TRANSPOSE(Der(i,c1:c2,d1:d2))
      END IF
    END IF
  END DO

  !Se calcula la matriz pseudoinversa de B
  ALLOCATE(MatB(Num1,Num2),MatAux(Num2,Num2),Diag(Num2))

  MatB(:,:)=MatrizB(:,:)
  CALL Pseudoinversa(MatB(:,:),Diag(:),MatAux(:,:),InvB(:,:))

  DEALLOCATE(MatB,MatAux,Diag)

END SUBROUTINE ConvertirCoordenadas

!-------------------------------------------------------------------------------
! Modifica la definición de las coordenadas internas añadiendo o eliminando
! determinadas coordenadas
!-------------------------------------------------------------------------------
! Fich:                  Fichero que contiene la modificación de coordenadas
! Linea:                 Línea que se va leyendo
! Signo:                 Signo inicial de la línea que se lee
! Copia:                 Copia de la definición de coordenadas modificada
! Aux1,Aux2:             Matrices auxiliares para añadir y eliminar coordenadas
! a,b,c,d:               Átomos que definen una coordenada
! Rep:                   Variable para determinar si la coordenada ya existe
! CoordMas:              Número de coordenadas a añadir
! CoordMen:              Número de coordenadas a eliminar
! CoordCong:             Número de coordenadas a congelar
! CoordComb:             Número de coordenadas como combinación lineal
! UMas,UMen,UCong,UComb: Unidades auxiliares
! Error:                 Variable para controlar los errores
! i,j,k:                 Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ModificarCoordenadas(Fich)
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich

  CHARACTER(LEN=LLL) :: Linea
  CHARACTER :: Signo
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Copia,Aux1,Aux2
  INTEGER :: Error,a,b,c,d,i,j,k,Rep,CoordMas,CoordMen,CoordCong,CoordComb, &
             UMas,UMen,UCong,UComb

  !Si se usan coordenadas cartesianas, no se modifican
  IF (TipoCoord <= 1) RETURN

  UMas=NuevaUnidad()
  OPEN(UMas,STATUS='SCRATCH',ACTION='READWRITE')
  UMen=NuevaUnidad()
  OPEN(UMen,STATUS='SCRATCH',ACTION='READWRITE')
  UCong=NuevaUnidad()
  OPEN(UCong,STATUS='SCRATCH',ACTION='READWRITE')
  UComb=NuevaUnidad()
  OPEN(UComb,STATUS='SCRATCH',ACTION='READWRITE')

  !Se identifican las coordenadas a añadir o eliminar
  CoordMas=0
  CoordMen=0
  CoordCong=0
  CoordComb=0
  DO
    CALL LeerSiguienteLinea(Fich,Linea,Error)
    IF (Error /= 0) EXIT

    !Si la línea comienza por '--', se eliminan todas las coordenadas
    IF (Linea(1:2) == '--') THEN
      DO i=1,SIZE(DefCoord,1)
        WRITE(UMen,*) DefCoord(i,:)
      END DO
      CoordMen=CoordMen+SIZE(DefCoord,1)
      CYCLE
    END IF

    !Si la línea tiene paréntesis se lee una combinación lineal
    IF (INDEX(Linea,'(') /= 0) THEN
      CALL LeerCombinacion
      CYCLE
    END IF

    !De cada línea, se lee un signo y hasta 4 números
    a=0; b=0; c=0; d=0
    Signo=Linea(1:1)
    IF ((Signo /= '+') .AND. (Signo /= '-') .AND. (Signo /= '*')) &
       CALL Mensaje('ModificarCoordenadas',17,.TRUE.)
    Linea=ADJUSTL(Linea(2:))

    SELECT CASE (Linea(1:1))
     !Traslación (X/Y/Z)
     CASE ('X')
      d=1
     CASE ('Y')
      d=2
     CASE ('Z')
      d=3
     !Rotación (RX/RY/RZ)
     CASE ('R')
      SELECT CASE (Linea(2:2))
       CASE ('X')
        c=1
       CASE ('Y')
        c=2
       CASE ('Z')
        c=3
      END SELECT
     !Cartesianas (CX/CY/CZ a)
     CASE ('C')
      SELECT CASE (Linea(2:2))
       CASE ('X')
        d=1
       CASE ('Y')
        d=2
       CASE ('Z')
        d=3
      END SELECT
      Linea=ADJUSTL(Linea(INDEX(Linea,' '):))
      IF (TRIM(Linea) /= '') READ(Linea,*) a
     !Distancias y ángulos (a b [c [d]])
     CASE DEFAULT
      IF (TRIM(Linea) /= '') READ(Linea,*) a
      Linea=ADJUSTL(Linea(INDEX(Linea,' '):))
      IF (TRIM(Linea) /= '') READ(Linea,*) b
      Linea=ADJUSTL(Linea(INDEX(Linea,' '):))
      IF (TRIM(Linea) /= '') READ(Linea,*) c
      Linea=ADJUSTL(Linea(INDEX(Linea,' '):))
      IF (TRIM(Linea) /= '') READ(Linea,*) d
      IF ((a == 0) .OR. (b == 0)) CALL Mensaje('ModificarCoordenadas',17,.TRUE.)
      IF (3*MAX(a,b,c,d) > SIZE(Masa,1)) &
        CALL Mensaje('ModificarCoordenadas',18,.TRUE.)
  
      !Se colocan los números en el orden adecuado
      IF (c == 0) THEN
        IF (b < a) CALL Intercambiar(a,b)
      ELSE IF (d == 0) THEN
        IF (c < a) CALL Intercambiar(a,c)
      ELSE IF (d < a) THEN
        CALL Intercambiar(a,d)
        CALL Intercambiar(b,c)
      END IF
    END SELECT

    !Se comprueba si la coordenada ya existe
    Rep=0
    DO i=1,SIZE(DefCoord,1)
      IF ((a == DefCoord(i,1)) .AND. (b == DefCoord(i,2)) .AND. &
          (c == DefCoord(i,3)) .AND. (d == DefCoord(i,4))) THEN
        Rep=1
        EXIT
      END IF
    END DO

    !Se almacena la coordenada para añadir o eliminar según corresponda
    IF (Signo == '+') THEN
      WRITE(UMas,*) a,b,c,d
      CoordMas=CoordMas+1
    END IF
    IF ((Signo == '-') .AND. (Rep == 1)) THEN
      WRITE(UMen,*) a,b,c,d
      CoordMen=CoordMen+1
    END IF
    IF (Signo == '*') THEN
      WRITE(UCong,*) a,b,c,d
      WRITE(UMas,*) a,b,c,d
      CoordCong=CoordCong+1
      CoordMas=CoordMas+1
    END IF
  END DO

  !Se leen las coordenadas a añadir
  ALLOCATE(Aux1(CoordMas,4))
  REWIND(UMas)
  DO i=1,CoordMas
    READ(UMas,*) Aux1(i,:)
  END DO
  CLOSE(UMas)

  !Se leen las coordenadas a eliminar, sin repetir
  ALLOCATE(Aux2(CoordMen,4))
  REWIND(UMen)
  k=1
  DO i=1,CoordMen
    READ(UMen,*) Aux2(k,:)
    Rep=0
    DO j=1,k-1
      IF ((Aux2(k,1) == Aux2(j,1)) .AND. (Aux2(k,2) == Aux2(j,2)) .AND. &
          (Aux2(k,3) == Aux2(j,3)) .AND. (Aux2(k,4) == Aux2(j,4))) THEN
        Rep=1
        EXIT
      END IF
    END DO
    IF (Rep == 0) k=k+1
  END DO
  CLOSE(UMen)
  CoordMen=k-1

  !Se crea una copia de la definición de coordenadas y se modifica
  ALLOCATE(Copia(SIZE(DefCoord,1)+CoordMas-CoordMen+CoordComb,4))
  !Se copian las coordenadas, excepto las que hay que eliminar
  k=1
  DO i=1,SIZE(DefCoord,1)
    Rep=0
    DO j=1,CoordMen
      IF ((Aux2(j,1) == DefCoord(i,1)) .AND. (Aux2(j,2) == DefCoord(i,2)) .AND.&
          (Aux2(j,3) == DefCoord(i,3)) .AND. (Aux2(j,4) == DefCoord(i,4))) THEN
        Rep=1
        EXIT
      END IF
    END DO
    IF (Rep == 1) CYCLE
    Copia(k,:)=DefCoord(i,:)
    k=k+1
  END DO
  !Se añaden las coordenadas que no existen ya
  DO i=1,CoordMas
    Rep=0
    DO j=1,k-1
      IF ((Aux1(i,1) == Copia(j,1)) .AND. (Aux1(i,2) == Copia(j,2)) .AND. &
          (Aux1(i,3) == Copia(j,3)) .AND. (Aux1(i,4) == Copia(j,4))) THEN
        Rep=1
        EXIT
      END IF
    END DO
    IF (Rep == 1) CYCLE
    Copia(k,:)=Aux1(i,:)
    k=k+1
  END DO

  !Se sustituye la matriz DefCoord por la modificada
  DEALLOCATE(DefCoord,Geometria,Cong)
  ALLOCATE(DefCoord(k-1,4),Geometria(k-1),Cong(k-1))
  DefCoord(:,:)=Copia(1:k-1,:)
  DEALLOCATE(Copia,Aux1,Aux2)

  !Se leen las coordenadas congeladas
  Cong=0
  REWIND(UCong)
  DO i=1,CoordCong
    READ(UCong,*) a,b,c,d
    DO j=1,SIZE(DefCoord,1)
      IF ((a == DefCoord(j,1)) .AND. (b == DefCoord(j,2)) .AND. &
          (c == DefCoord(j,3)) .AND. (d == DefCoord(j,4))) THEN
        Cong(j) = 1
        EXIT
      END IF
    END DO
  END DO
  CLOSE(UCong)

  !Se leen los factores de las combinaciones
  k=k-1-CoordComb
  IF (ALLOCATED(Combinacion)) DEALLOCATE(Combinacion)
  ALLOCATE(Combinacion(CoordComb,k))
  REWIND(UComb)
  Combinacion(:,:)=0.0D0
  DO i=1,CoordComb
    DO
      READ(UComb,*,IOSTAT=Error) a,b,c,d
      IF (Error /= 0) EXIT
      DO j=1,SIZE(DefCoord,1)
        IF ((a == DefCoord(j,1)) .AND. (b == DefCoord(j,2)) .AND. &
            (c == DefCoord(j,3)) .AND. (d == DefCoord(j,4))) THEN
          READ(UComb,*) Combinacion(i,j)
          EXIT
        END IF
      END DO
    END DO
  END DO
  CLOSE(UComb)

  CONTAINS

  !Subrutina que procesa una combinación lineal de coordenadas
  !Lin:       Línea que se va procesando
  !Factor:    Factor de la combinación lineal
  !Num1,Num2: Posiciones de los paréntesis
  SUBROUTINE LeerCombinacion
    CHARACTER(LEN=LLL) :: Lin
    DOUBLE PRECISION :: Factor
    INTEGER :: Num1,Num2
    
    !Lee el signo
    Signo=Linea(1:1)
    IF ((Signo /= '+') .AND. (Signo /= '*')) &
       CALL Mensaje('ModificarCoordenadas',17,.TRUE.)
    Linea=ADJUSTL(Linea(2:))

    Num1=INDEX(Linea,'(')
    Num2=INDEX(Linea,')')
    DO WHILE ((Num1 > 0) .AND. (Num2 > 2))
      !Lee el factor de la combinación lineal
      Lin=Linea(:Num1-1)
      READ(Lin,*) Factor
      !Lee la coordenada de la combinación
      Lin=Linea(Num1+1:Num2-1)
      DO i=1,LEN(Lin)
        IF (Lin(i:i) == '-') Lin(i:i)=' '
      END DO
      a=0; b=0; c=0; d=0
      IF (TRIM(Lin) /= '') READ(Lin,*) a
      Lin=ADJUSTL(Lin(INDEX(Lin,' '):))
      IF (TRIM(Lin) /= '') READ(Lin,*) b
      Lin=ADJUSTL(Lin(INDEX(Lin,' '):))
      IF (TRIM(Lin) /= '') READ(Lin,*) c
      Lin=ADJUSTL(Lin(INDEX(Lin,' '):))
      IF (TRIM(Lin) /= '') READ(Lin,*) d
      IF (3*MAX(a,b,c,d) > SIZE(Masa,1)) &
        CALL Mensaje('ModificarCoordenadas',18,.TRUE.)
      IF (c == 0) THEN
        IF (b < a) CALL Intercambiar(a,b)
      ELSE IF (d == 0) THEN
        IF (c < a) CALL Intercambiar(a,c)
      ELSE IF (d < a) THEN
        CALL Intercambiar(a,d)
        CALL Intercambiar(b,c)
      END IF
      !Añade la coordenada
      WRITE(UMas,*) a,b,c,d
      CoordMas=CoordMas+1
      Linea=TRIM(Linea(Num2+1:))
      !Añade la combinación
      WRITE(UComb,*) a,b,c,d
      WRITE(UComb,*) Factor
      Num1=INDEX(Linea,'(')
      Num2=INDEX(Linea,')')
    END DO
    WRITE(UComb,*) '----'

    !Se añade la combinación
    !"a" es el número de la combinación en negativo
    !"c" es -1 si la última coordenada es un ángulo o diedro
    CoordMas=CoordMas+1
    CoordComb=CoordComb+1
    a=-CoordComb
    IF (c /= 0) c=-1
    WRITE(UMas,*) a,0,c,0
    IF (Signo == '*') THEN
      CoordCong=CoordCong+1
      WRITE(UCong,*) a,0,c,0
    END IF
  END SUBROUTINE LeerCombinacion

END SUBROUTINE ModificarCoordenadas

!-------------------------------------------------------------------------------
! Convierte el gradiente y la hessiana de coordenadas cartesianas a coordenadas
! de trabajo, o viceversa
!-------------------------------------------------------------------------------
! Tipo:    -1 -> a cartesianas, 1 -> de cartesianas
! Grad:    Gradiente en coordenadas cartesianas
! Hess:    Hessiana en coordenadas cartesianas
! MatK:    Matriz K para la conversión de la hessiana
! Num:     Número de coordenadas cartesianas
! i,j:     Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ConvertirGradHess(Tipo,Grad,Hess)
  USE Utilidades
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Tipo
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT), OPTIONAL :: Grad
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: Hess

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MatK
  INTEGER :: Num,i,j

  IF (PRESENT(Grad)) THEN
    SELECT CASE (Tipo)
     CASE (-1)
      !Convierte el gradiente de las coordenadas de trabajo a cartesianas
      Grad(:)=MATMUL(TRANSPOSE(MatrizB(:,:)),Gradiente(:))
     CASE (1)
      !Convierte el gradiente de coordenadas cartesianas a las de trabajo
      IF (ALLOCATED(Gradiente)) DEALLOCATE(Gradiente)
      ALLOCATE(Gradiente(SIZE(DefCoord,1)))
      Gradiente(:)=MATMUL(TRANSPOSE(InvB(:,:)),Grad(:))
    END SELECT
  END IF

  IF (PRESENT(Hess)) THEN
    !La conversión de la hessiana a internas necesita el gradiente
    IF ((Tipo == 1) .AND. (.NOT. PRESENT(Grad))) &
       CALL Mensaje('ConvertirGradHess',0,.TRUE.)
    !Se calcula la matriz auxiliar K
    Num=SIZE(Hess,1)
    ALLOCATE(MatK(Num,Num))
    MatK(:,:)=0.0D0
    IF (TipoCoord > 1) THEN
      DO i=1,Num
        DO j=1,Num
          MatK(i,j)=DOT_PRODUCT(Gradiente(:),Der(:,i,j))
        END DO
      END DO
    END IF
    SELECT CASE (Tipo)
     CASE (-1)
      !Convierte la hessiana de las coordenadas de trabajo a cartesianas
      Hess(:,:)=MATMUL(MATMUL(TRANSPOSE(MatrizB(:,:)),Hessiana(:,:)), &
                       MatrizB(:,:)) + MatK(:,:)
     CASE (1)
      !Convierte la hessiana de coordenadas cartesianas a las de trabajo
      IF (ALLOCATED(Hessiana)) DEALLOCATE(Hessiana)
      ALLOCATE(Hessiana(SIZE(DefCoord,1),SIZE(DefCoord,1)))
      Hessiana(:,:)=MATMUL(MATMUL(TRANSPOSE(InvB(:,:)),Hess(:,:)-MatK(:,:)), &
                           InvB(:,:))
    END SELECT
    DEALLOCATE(MatK)

    !Asegura que la matriz es simétrica
    DO i=1,SIZE(Hessiana,1)
      DO j=i+1,SIZE(Hessiana,2)
        Hessiana(i,j)=0.5D0*(Hessiana(i,j)+Hessiana(j,i))
        Hessiana(j,i)=Hessiana(i,j)
      END DO
    END DO
  END IF

END SUBROUTINE ConvertirGradHess

!-------------------------------------------------------------------------------
! Convierte un incremento de coordenadas
! de coordenadas de trabajo a coordenadas cartesianas
!  V. Bakken, T. Helgaker; J. Chem. Phys. 117 (2002) 9160-9174
!-------------------------------------------------------------------------------
! Paso:       Incremento de coordenadas de trabajo
!             A la salida contiene el paso real
! CoordCart:  Coordenadas cartesianas iniciales
! PasoCart:   Incremento de coordenadas cartesianas
! Obj:        Coordenadas de trabajo finales
! Dif:        Diferencia entre las coordenadas finales y las conseguidas
! Pas:        Paso parcial calculado en coordenadas cartesianas
! Ini:        Diferencia inicial tras la primera conversión
! Dist:       Tamaño del paso parcial
! DistAnt:    Tamaño del paso anterior
! GeomCart:   Coordenadas cartesianas obtenidas
! Conv,Conv2: Criterios de convergencia
! MaxIt:      Número máximo de iteraciones para la conversión
! b,c,d:      Átomos implicados (para tratar los ángulos)
! Error:      Variable para controlar los errores
! i,j:        Contadores
!-------------------------------------------------------------------------------
SUBROUTINE ConvertirIncremento(Paso,CoordCart,PasoCart)
  USE Unidades
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: CoordCart
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Paso
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: PasoCart

  DOUBLE PRECISION, DIMENSION(SIZE(Paso,1)) :: Obj,Dif
  DOUBLE PRECISION, DIMENSION(SIZE(CoordCart,1)) :: Pas,GeomCart
  DOUBLE PRECISION, PARAMETER :: Conv=1.0D-6,Conv2=1.0D-12
  DOUBLE PRECISION :: Dist,DistAnt,Ini
  INTEGER, PARAMETER :: MaxIt=25
  INTEGER :: i,j,b,c,d,Error

  !Se calculan las coordenadas finales deseadas
  Obj(:)=Geometria(:)+Paso(:)
  DO i=1,SIZE(DefCoord,1)
    b=DefCoord(i,2)
    c=DefCoord(i,3)
    d=DefCoord(i,4)
    IF (Cong(i) == 1) Obj(i)=Geometria(i)
    IF (c == 0) CYCLE
    IF (b == 0) THEN
      !No hay valor asignado a las coordenadas de rotación
      Obj(i)=0.0D0
     ELSE IF (d == 0) THEN
      !Corrección de los ángulos
      Obj(i)=MOD(ABS(Obj(i)),2.0D0*Pi)
      IF (Obj(i) > Pi) Obj(i)=2.0D0*Pi-Obj(i)
     ELSE
      !Corrección de los diedros
      Obj(i)=MOD(Obj(i),2.0D0*Pi)
      IF (ABS(Obj(i)) > Pi) Obj(i)=MODULO(Obj(i),-SIGN(2.0D0*Pi,Obj(i)))
    END IF
  END DO

  !Se convierte el paso a coordenadas cartesianas iterativamente
  Error=0
  Dif(:)=Paso(:)
  Dist=0.0D0
  PasoCart(:)=0.0D0
  Ini=0.0D0
  DO i=1,MaxIt
    !Estimación del paso equivalente en coordenadas cartesianas
    Pas(:)=MATMUL(InvB(:,:),Dif(:))
    !Comprobación de convergencia
    DistAnt=Dist
    Dist=Norma(Pas)/SQRT(DBLE(SIZE(Paso,1)))
    IF ((Dist < Conv) .OR. (ABS(Dist-DistAnt) < Conv2)) EXIT
    !Actualización del paso en coordenadas cartesianas
    !y de las coordenadas cartesianas obtenidas
    PasoCart(:)=PasoCart(:)+Pas(:)
    GeomCart(:)=CoordCart(:)+PasoCart(:)
    !Se halla la correspondencia en coordenadas de trabajo
    CALL ConvertirCoordenadas(GeomCart,1)
    !Diferencia entre el objetivo y las coordenadas obtenidas
    Dif(:)=Obj(:)-Geometria(:)
    DO j=1,SIZE(DefCoord,1)
      c=DefCoord(j,3)
      IF (c == 0) CYCLE
      Dif(j)=MOD(Dif(j),2.0D0*Pi)
      IF (ABS(Dif(j)) > Pi) Dif(j)=MODULO(Dif(j),-SIGN(2.0D0*Pi,Dif(j)))
    END DO
    !Se vuelve a la primera estimación si hay divergencia
    IF (i == 1) THEN
      Ini=Norma(Dif)
    ELSE IF (Norma(Dif) > Ini) THEN
      PasoCart(:)=MATMUL(InvB(:,:),Paso(:))
      CALL ConvertirCoordenadas(CoordCart+PasoCart,1)
      CALL Mensaje('ConvertirIncremento',19,.FALSE.)
      EXIT
    END IF
  END DO
  IF (i > MaxIt) Error=1

  !Se restablecen las coordenadas congeladas
  IF (ANY(Cong == 1)) THEN
    DO i=1,MaxIt
      !Se calcula la diferencia sólo en las coordandas congeladas
      Dif(:)=Cong*(Obj(:)-Geometria(:))
      DO j=1,SIZE(DefCoord,1)
        c=DefCoord(j,3)
        IF (c == 0) CYCLE
        Dif(j)=MOD(Dif(j),2.0D0*Pi)
        IF (ABS(Dif(j)) > Pi) Dif(j)=MODULO(Dif(j),-SIGN(2.0D0*Pi,Dif(j)))
      END DO
      !Se comprueba la convergencia
      Dist=Norma(Dif)/SQRT(DBLE(COUNT(Cong == 1)))
      IF (Dist < Conv2) EXIT
      !Se calcula el paso en coordenadas cartesianas
      Pas(:)=MATMUL(InvB(:,:),Dif(:))
      PasoCart(:)=PasoCart(:)+Pas(:)
      GeomCart(:)=CoordCart(:)+PasoCart(:)
      CALL ConvertirCoordenadas(GeomCart,1)
    END DO
  END IF
  IF (i > MaxIt) Error=2

  IF (Error > 0) CALL Mensaje('ConvertirIncremento',2,.FALSE.)

END SUBROUTINE ConvertirIncremento

!-------------------------------------------------------------------------------
! Proyecta el gradiente y la hessiana para mantenerlos en las coordenadas
! de trabajo y para aplicar restricciones
!  C. Peng, P.Y. Ayala, M.J. Frisch; J. Comput. Chem. 17 (1996) 49-56
!-------------------------------------------------------------------------------
! Desp:     Desplazamiento en coordenadas de trabajo
! Grad:     Gradiente en coordenadas de trabajo
! Hess:     Hessiana en coordenadas de trabajo
! Proy:     Proyector en el espacio de las coordenadas de trabajo
! P1,P2,P3: Matrices auxiliares para la proyección
! D:        Vector diagonal para la inversión
! i,j:      Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Proyectar(Desp,Grad,Hess)
  USE Utilidades
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT), OPTIONAL :: Desp,Grad
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT), OPTIONAL :: Hess

  DOUBLE PRECISION, DIMENSION(SIZE(DefCoord,1),SIZE(DefCoord,1)) :: Proy
  DOUBLE PRECISION, DIMENSION(SIZE(DefCoord,1),SIZE(DefCoord,1)) :: P1,P2,P3
  DOUBLE PRECISION, DIMENSION(SIZE(DefCoord,1)) :: D
  INTEGER :: i,j

  !Cálculo de la matriz de proyección
  Proy(:,:)=MATMUL(MatrizB(:,:),InvB(:,:))

  !Se aplican las restricciones si hay coordenadas congeladas
  ! P' = P - PC (CPC)^-1 CP
  IF (ANY(Cong == 1)) THEN
    !P3 = (CPC)^-1
    P1(:,:)=Proy(:,:)
    DO i=1,SIZE(Proy,1)
      IF (Cong(i) == 0) THEN
        P1(:,i)=0.0D0
        P1(i,:)=0.0D0
      END IF
    END DO
    CALL Pseudoinversa(P1(:,:),D(:),P2(:,:),P3(:,:))
    !P1 = PC, P2 = CP
    P1(:,:)=0.0D0
    P2(:,:)=0.0D0
    DO i=1,SIZE(Proy,1)
      IF (Cong(i) == 1) THEN
        P1(:,i)=Proy(:,i)
        P2(i,:)=Proy(i,:)
      END IF
    END DO
    Proy(:,:)=Proy(:,:)-MATMUL(MATMUL(P1(:,:),P3(:,:)),P2(:,:))
  END IF

  !Proyección del desplazamiento
  IF (PRESENT(Desp)) THEN
    Desp(:)=MATMUL(Proy(:,:),Desp(:))
  END IF

  !Proyección del gradiente
  IF (PRESENT(Grad)) THEN
    Grad(:)=MATMUL(Proy(:,:),Grad(:))
  END IF

  !Proyección de la hessiana
  IF (PRESENT(Hess)) THEN
    Hess(:,:)=MATMUL(MATMUL(Proy(:,:),Hess(:,:)),Proy(:,:))
    !Asegura que la matriz es simétrica
    DO i=1,SIZE(Hess,1)
      DO j=i+1,SIZE(Hess,2)
        Hess(i,j)=0.5D0*(Hess(i,j)+Hess(j,i))
        Hess(j,i)=Hess(i,j)
      END DO
    END DO
  END IF

END SUBROUTINE Proyectar

END MODULE Coordenadas
