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
! Diversas subrutinas y funciones útiles
!-------------------------------------------------------------------------------
MODULE Utilidades

#ifndef LANG
#define LANG Espanol
#endif

!-------------------------------------------------------------------------------
! SG:      Fichero de salida general
!-------------------------------------------------------------------------------
INTEGER, PARAMETER :: SG=21

INTERFACE Intercambiar
  MODULE PROCEDURE IntercambiarEntero,IntercambiarEscalar,IntercambiarVector
END INTERFACE

CONTAINS
!Mensaje(Nombre,Codigo,Error)
!NuevaUnidad()
!IniciarAleatorio()
!LeerSiguienteLinea(Fich,Lin,Fin)
!PasarMinusculas(Lin)
!Intercambiar(A,B)
!Norma(Vec)
!Normalizar(Vec)
!Determinante(Mat)
!ProductoTens(Vec1,Vec2)
!ProductoVect(Vec1,Vec2)
!ProductoCuat(Cuat1,Cuat2)
!RotarCuaternion(Cuat,Puntos)
!Distancia(Pos1,Pos2)
!Angulo(Pos1,Pos2,Pos3)
!Diedro(Pos1,Pos2,Pos3,Pos4)
!Hipotenusa(A,B)
!Diagonalizar(Mat,Vec,Val,Ord)
!ResolverLU(Coef,Term)
!Pseudoinversa(Matriz,Diag,Vect,Inv)

!-------------------------------------------------------------------------------
! Escribe un mensaje y detiene el programa si es un error
!-------------------------------------------------------------------------------
! Nombre:  Nombre de la subrutina que ha generado el error
! Codigo:  Código del error
! Error:   0 -> Continúa el programa, 1 -> Se detiene el programa
!-------------------------------------------------------------------------------
SUBROUTINE Mensaje(Nombre,Codigo,Error)
  USE LANG
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: Nombre
  INTEGER, INTENT(IN) :: Codigo
  LOGICAL, INTENT(IN) :: Error

  IF (Error) THEN
    WRITE(6,101) Codigo,Nombre
    WRITE(6,100) TRIM(Errores(Codigo))
    STOP
   ELSE
    WRITE(6,102) Codigo,Nombre
    WRITE(6,100) TRIM(Errores(Codigo))
  END IF

100 FORMAT(A)
101 FORMAT('(EE) ',I3.3,' - ',A)
102 FORMAT('(II) ',I3.3,' - ',A)

END SUBROUTINE Mensaje

!-------------------------------------------------------------------------------
! Obtiene el siguiente número de unidad sin usar
!-------------------------------------------------------------------------------
! Uso:     Variable para saber si la unidad está en uso
! i:       Contador
!-------------------------------------------------------------------------------
FUNCTION NuevaUnidad()
  IMPLICIT NONE
  INTEGER :: NuevaUnidad

  INTEGER :: i
  LOGICAL :: Uso

  NuevaUnidad=0
  DO i=11,999
    INQUIRE(i,OPENED=Uso)
    IF (.NOT. Uso) THEN
      NuevaUnidad=i
      EXIT
    END IF
  END DO
  IF (NuevaUnidad == 0) CALL Mensaje('NuevaUnidad',31,.TRUE.)

END FUNCTION

!-------------------------------------------------------------------------------
! Inicia el generador de números aleatorios según la fecha y hora
!-------------------------------------------------------------------------------
! Semilla: Vector semilla
! Datos:   Vector con los datos de fecha y hora
! N:       Longitud del vector semilla
!-------------------------------------------------------------------------------
SUBROUTINE IniciarAleatorio(Inicio)
  IMPLICIT NONE
  INTEGER, INTENT(IN), OPTIONAL :: Inicio

  INTEGER, DIMENSION(:), ALLOCATABLE :: Semilla
  INTEGER, DIMENSION(8) :: Datos
  INTEGER :: N

  !Obtiene los datos de la fecha y hora, y el tamaño de la semilla
  IF (PRESENT(Inicio)) THEN
    Datos(:)=Inicio
   ELSE
    CALL DATE_AND_TIME(VALUES=Datos)
  END IF
  CALL RANDOM_SEED(SIZE=N)
  ALLOCATE(Semilla(N))

  !Inicia la semilla
  !(1): número de milisegundos en el día
  !(2): número de días desde 1900
  Semilla(1)=Datos(8)+1000*(Datos(7)+60*(Datos(6)+60*Datos(5)))
  Semilla(2:)=(Datos(3)-1)+30*(Datos(2)-1)+366*(Datos(1)-1900)
  CALL RANDOM_SEED(PUT=Semilla)

  DEALLOCATE(Semilla)
  
END SUBROUTINE IniciarAleatorio

!-------------------------------------------------------------------------------
! Lee la siguiente línea que no esté vacía ni sea un comentario
!-------------------------------------------------------------------------------
! Fich:    Número del fichero a leer
! Lin:     La línea que se va leyendo
! Fin:     Indica si se ha acabado el fichero (1) o no (0)
!-------------------------------------------------------------------------------
SUBROUTINE LeerSiguienteLinea(Fich,Lin,Fin)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Fich
  CHARACTER(LEN=*), INTENT(OUT) :: Lin
  INTEGER, INTENT(OUT) :: Fin

  Fin=0
  DO
    READ(Fich,'(A)',IOSTAT=Fin) Lin
    IF (Fin /= 0) EXIT
    Lin=ADJUSTL(Lin)
    IF ((TRIM(Lin) /= '') .AND. (Lin(1:1) /= '#')) EXIT
  END DO
  IF (INDEX(Lin,' #') > 0) THEN
    Lin=Lin(:INDEX(Lin,' #'))
  END IF

END SUBROUTINE LeerSiguienteLinea

!-------------------------------------------------------------------------------
! Pasa una cadena de caracteres a minúsculas
!-------------------------------------------------------------------------------
! Linea:   La cadena de caracteres
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE PasarMinusculas(Lin)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(INOUT) :: Lin

  INTEGER :: Letra,i

  DO i=1,LEN(TRIM(Lin))
    Letra=ICHAR(Lin(i:i))
    IF ((Letra > 64) .AND. (Letra < 91)) Lin(i:i)=CHAR(Letra+32)
  END DO

END SUBROUTINE PasarMinusculas

!-------------------------------------------------------------------------------
! Intercambia dos números o dos vectores
!-------------------------------------------------------------------------------
! A,B:     Los números o vectores a intercambiar
! Aux:     Variable auxiliar para el intercambio
!-------------------------------------------------------------------------------
SUBROUTINE IntercambiarEntero(A,B)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: A,B
  INTEGER :: Aux
  Aux=A; A=B; B=Aux
END SUBROUTINE IntercambiarEntero
SUBROUTINE IntercambiarEscalar(A,B)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT) :: A,B
  DOUBLE PRECISION :: Aux
  Aux=A; A=B; B=Aux
END SUBROUTINE IntercambiarEscalar
SUBROUTINE IntercambiarVector(A,B)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: A,B
  DOUBLE PRECISION, DIMENSION(SIZE(A,1)) :: Aux
  Aux(:)=A(:); A(:)=B(:); B(:)=Aux(:)
END SUBROUTINE IntercambiarVector

!-------------------------------------------------------------------------------
! Calcula la norma (módulo) de un vector
!-------------------------------------------------------------------------------
! Vec:     El vector
!-------------------------------------------------------------------------------
FUNCTION Norma(Vec)
  IMPLICIT NONE
  DOUBLE PRECISION :: Norma
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Vec

  Norma=SQRT(DOT_PRODUCT(Vec(:),Vec(:)))

END FUNCTION Norma

!-------------------------------------------------------------------------------
! Normaliza un vector
!-------------------------------------------------------------------------------
! Vec:     El vector
!-------------------------------------------------------------------------------
FUNCTION Normalizar(Vec)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Vec
  DOUBLE PRECISION, DIMENSION(SIZE(Vec,1)) :: Normalizar

  Normalizar(:)=Vec(:)/Norma(Vec(:))

END FUNCTION Normalizar

!-------------------------------------------------------------------------------
! Calcula el determinante de una matriz 3×3
!-------------------------------------------------------------------------------
! Mat:     La matriz
!-------------------------------------------------------------------------------
FUNCTION Determinante(Mat)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: Mat
  DOUBLE PRECISION :: Determinante

  Determinante=Mat(1,1)*Mat(2,2)*Mat(3,3) + &
               Mat(1,2)*Mat(2,3)*Mat(3,1) + &
               Mat(1,3)*Mat(2,1)*Mat(3,2) - &
               Mat(1,1)*Mat(2,3)*Mat(3,2) - &
               Mat(2,2)*Mat(3,1)*Mat(1,3) - &
               Mat(3,3)*Mat(1,2)*Mat(2,1)

END FUNCTION Determinante

!-------------------------------------------------------------------------------
! Calcula el producto tensorial de dos vectores
! P(i,j)=V1(i)*V2(j)
!-------------------------------------------------------------------------------
! Vec1,Vec2: Los dos vectores a multiplicar
!-------------------------------------------------------------------------------
FUNCTION ProductoTens(Vec1,Vec2)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Vec1,Vec2
  DOUBLE PRECISION, DIMENSION(SIZE(Vec1,1),SIZE(Vec2,1)) :: ProductoTens

  ProductoTens=SPREAD(Vec1(:),DIM=2,NCOPIES=SIZE(Vec2(:),1))* &
               SPREAD(Vec2(:),DIM=1,NCOPIES=SIZE(Vec1(:),1))

END FUNCTION ProductoTens

!-------------------------------------------------------------------------------
! Calcula el producto vectorial de dos vectores 3D (V1 x v2)
!-------------------------------------------------------------------------------
! Vec1,Vec2: Los dos vectores
!-------------------------------------------------------------------------------
FUNCTION ProductoVect(Vec1,Vec2)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Vec1,Vec2
  DOUBLE PRECISION, DIMENSION(3) :: ProductoVect

  ProductoVect(1)=Vec1(2)*Vec2(3)-Vec1(3)*Vec2(2)
  ProductoVect(2)=Vec1(3)*Vec2(1)-Vec1(1)*Vec2(3)
  ProductoVect(3)=Vec1(1)*Vec2(2)-Vec1(2)*Vec2(1)

END FUNCTION ProductoVect

!-------------------------------------------------------------------------------
! Calcula el producto de dos cuaterniones (rotación compuesta)
!-------------------------------------------------------------------------------
! Cuat1,Cuat2: Los dos cuaterniones
!-------------------------------------------------------------------------------
FUNCTION ProductoCuat(Cuat1,Cuat2)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(4), INTENT(IN) :: Cuat1,Cuat2
  DOUBLE PRECISION, DIMENSION(4) :: ProductoCuat

  ProductoCuat(1)=Cuat1(1)*Cuat2(1)- &
                  Cuat1(2)*Cuat2(2)-Cuat1(3)*Cuat2(3)-Cuat1(4)*Cuat2(4)
  ProductoCuat(2)=Cuat1(1)*Cuat2(2)+ &
                  Cuat1(2)*Cuat2(1)+Cuat1(3)*Cuat2(4)-Cuat1(4)*Cuat2(3)
  ProductoCuat(3)=Cuat1(1)*Cuat2(3)- &
                  Cuat1(2)*Cuat2(4)+Cuat1(3)*Cuat2(1)+Cuat1(4)*Cuat2(2)
  ProductoCuat(4)=Cuat1(1)*Cuat2(4)+ &
                  Cuat1(2)*Cuat2(3)-Cuat1(3)*Cuat2(2)+Cuat1(4)*Cuat2(1)

  !Se cambia el signo para que el primer elemento sea positivo
  IF (ProductoCuat(1) < 0.0D0) ProductoCuat(:)=-ProductoCuat(:)

END FUNCTION ProductoCuat

!-------------------------------------------------------------------------------
! Rota una serie de puntos según un cuaternión dado
!-------------------------------------------------------------------------------
! Puntos:  Matriz que contiene los puntos por filas
! Cuat:    Cuaternión con el que se rotan los puntos
! Matriz:  Matriz de rotación correspondiente al cuaternión
!-------------------------------------------------------------------------------
FUNCTION RotarCuaternion(Puntos,Cuat)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: Puntos
  DOUBLE PRECISION, DIMENSION(4), INTENT(IN) :: Cuat
  DOUBLE PRECISION, DIMENSION(SIZE(Puntos,1),3) :: RotarCuaternion
  DOUBLE PRECISION, DIMENSION(3,3) :: Matriz

  !Se calcula la matriz de rotación para el cuaternión
  ! (es la transpuesta, ver abajo)
  Matriz(1,1)=Cuat(1)**2+Cuat(2)**2-Cuat(3)**2-Cuat(4)**2
  Matriz(2,2)=Cuat(1)**2-Cuat(2)**2+Cuat(3)**2-Cuat(4)**2
  Matriz(3,3)=Cuat(1)**2-Cuat(2)**2-Cuat(3)**2+Cuat(4)**2
  Matriz(1,2)=2.0D0*(Cuat(2)*Cuat(3)+Cuat(1)*Cuat(4))
  Matriz(2,1)=2.0D0*(Cuat(2)*Cuat(3)-Cuat(1)*Cuat(4))
  Matriz(2,3)=2.0D0*(Cuat(3)*Cuat(4)+Cuat(1)*Cuat(2))
  Matriz(3,2)=2.0D0*(Cuat(3)*Cuat(4)-Cuat(1)*Cuat(2))
  Matriz(3,1)=2.0D0*(Cuat(2)*Cuat(4)+Cuat(1)*Cuat(3))
  Matriz(1,3)=2.0D0*(Cuat(2)*Cuat(4)-Cuat(1)*Cuat(3))

  !Se aplica la matriz de rotación a todos los puntos
  !al estar los puntos por filas, se invierte la rotación: (M·P^T)^T = P·M^T
  RotarCuaternion(:,:)=MATMUL(Puntos(:,1:3),Matriz(:,:))

END FUNCTION

!-------------------------------------------------------------------------------
! Calcula la distancia entre dos puntos 3D
!-------------------------------------------------------------------------------
! Pos1,Pos2: Los dos puntos
!-------------------------------------------------------------------------------
FUNCTION Distancia(Pos1,Pos2)
  IMPLICIT NONE
  DOUBLE PRECISION :: Distancia
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Pos1,Pos2

  Distancia=Norma(Pos1(:)-Pos2(:))

END FUNCTION Distancia

!-------------------------------------------------------------------------------
! Calcula el ángulo entre tres puntos 3D (en radianes)
!-------------------------------------------------------------------------------
! Pos1,Pos2,Pos3: Los tres puntos
! V1,V2:          Vectores que forman el ángulo
!-------------------------------------------------------------------------------
FUNCTION Angulo(Pos1,Pos2,Pos3)
  IMPLICIT NONE
  DOUBLE PRECISION :: Angulo
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Pos1,Pos2,Pos3

  DOUBLE PRECISION, DIMENSION(3) :: V1,V2

  V1(:)=Pos1(:)-Pos2(:)
  V1(:)=V1(:)/Norma(V1(:))
  V2(:)=Pos3(:)-Pos2(:)
  V2(:)=V2(:)/Norma(V2(:))
  Angulo=ACOS(DOT_PRODUCT(V1(:),V2(:)))

END FUNCTION Angulo

!-------------------------------------------------------------------------------
! Calcula el ángulo diedro entre cuatro puntos 3D (en radianes)
!-------------------------------------------------------------------------------
! Pos1,Pos2,Pos3,Pos4: Los cuatro puntos
! V1,V2,V3:            Vectores que forman el diedro
! N1,N2:               Normales de los planos
! Aux:                 Variable auxiliar
!-------------------------------------------------------------------------------
FUNCTION Diedro(Pos1,Pos2,Pos3,Pos4)
  IMPLICIT NONE
  DOUBLE PRECISION :: Diedro
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: Pos1,Pos2,Pos3,Pos4

  DOUBLE PRECISION :: Aux
  DOUBLE PRECISION, DIMENSION(3) :: V1,V2,V3,N1,N2

  V1(:)=Pos1(:)-Pos2(:)
  V2(:)=Pos3(:)-Pos2(:)
  V3(:)=Pos4(:)-Pos3(:)
  N1(:)=ProductoVect(V1(:),V2(:))
  N1(:)=N1(:)/Norma(N1(:))
  N2(:)=ProductoVect(V3(:),V2(:))
  N2(:)=N2(:)/Norma(N2(:))
  V1(:)=ProductoVect(N1(:),N2(:))
  !El coseno tiene que estar en [-1,+1]
  Aux=MIN(1.0D0,MAX(-1.0D0,DOT_PRODUCT(N1(:),N2(:))))
  Diedro=SIGN(ACOS(Aux),DOT_PRODUCT(V1(:),V2(:)))

END FUNCTION Diedro

!-------------------------------------------------------------------------------
! Calcula (A^2+B^2)^(1/2) minimizando los errores
!-------------------------------------------------------------------------------
! A,B:       Valores de los "catetos"
! AbsA,AbsB: Valores absolutos de A y B
!-------------------------------------------------------------------------------
FUNCTION Hipotenusa(A,B)
  IMPLICIT NONE
  DOUBLE PRECISION :: Hipotenusa
  DOUBLE PRECISION, INTENT(IN) :: A,B

  DOUBLE PRECISION :: AbsA,AbsB

  AbsA=ABS(A)
  AbsB=ABS(B)
  IF (AbsA > AbsB) THEN
    Hipotenusa=AbsA*SQRT(1.0D0+(AbsB/AbsA)**2)
   ELSE
    IF (AbsB == 0.0D0) THEN
      Hipotenusa=0.0D0
     ELSE
      Hipotenusa=AbsB*SQRT(1.0D0+(AbsA/AbsB)**2)
    END IF
  END IF

END FUNCTION Hipotenusa

!-------------------------------------------------------------------------------
! Esta subrutina diagonaliza una matriz simétrica por el metodo de Jacobi y
! devuelve sus vectores y valores propios, ordenados o no
!-------------------------------------------------------------------------------
! Mat:      Matriz a diagonalizar
! Vec:      Matriz que contiene los vectores propios en columnas
! Val:      Vector que contiene los valores propios
! Ord:      1 -> se ordenan los vectores propios, 0 -> no se ordenan
! Num:      Dimensión de las matrices y vectores
! Preci:    Precisión para comparar con cero
! MaxIt:    Número máximo de vueltas para diagonalizar
! Suma:     Suma de los elementos extradiagonales
! Tol:      Tolerancia para hacer comparaciones
! Zeta,Tau: Variables para la diagonalización
! t,c,s,h:  Variables auxiliares para la diagonalización
! e:        Vector auxiliar para la diagonalización
! i,j,k:    Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Diagonalizar(Mat,Vec,Val,Ord)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: Mat
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: Vec
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Val
  INTEGER, INTENT(IN) :: Ord

  DOUBLE PRECISION, DIMENSION(SIZE(Mat,1)) :: e
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-12
  DOUBLE PRECISION :: t,c,s,h,Suma,Tol,Zeta,Tau
  INTEGER :: i,j,k,Num,MaxIt

  !Se inician los datos
  MaxIt=5000
  Num=SIZE(Mat,1)
  Vec(:,:)=0.0D0
  DO i=1,Num
    Vec(i,i)=1.0D0
    Val(i)=Mat(i,i)
  END DO

  !Se comprueba si la matriz es simétrica
  Tol=0.0D0
  DO i=1,Num
    DO j=i+1,Num
      Tol=MAX(Tol,ABS(Mat(i,j)-Mat(j,i)))
    END DO
  END DO
  IF (Tol > Preci) CALL Mensaje('Diagonalizar',1,.TRUE.)

  !Se diagonaliza la matriz
  DO i=1,MaxIt
    Suma=0.0D0
    DO j=1,Num-1
      Suma=Suma+SUM(ABS(Mat(j,j+1:Num)))
    END DO
    IF (Suma == 0.0D0) EXIT
    DO j=1,Num-1
      DO k=j+1,Num
        h=Mat(j,k)
        Tol=100.0D0*ABS(h)
        IF ((i < 4) .AND. (Tol < 20.0D0*Suma/(Num*Num))) CYCLE
        IF ((Val(k)+Tol == Val(k)) .AND. (Val(j)+Tol == Val(j))) THEN
          Mat(j,k) = 0.0D0
          CYCLE
        END IF
        Zeta=0.5D0*(Val(k)-Val(j))/h
        IF (Val(k)-Val(j)+Tol == Val(k)-Val(j)) THEN
          t=h/(Val(k)-Val(j))
         ELSE
          t=1.0D0/(ABS(Zeta)+SQRT(Zeta*Zeta+1.0D0))
          IF (Zeta < 0.0D0) t=-t
        END IF
        c=1.0D0/SQRT(t*t+1.0D0)
        s=t*c
        Tau=s/(1.0D0+c)
        Val(j)=Val(j)-t*h
        Val(k)=Val(k)+t*h
        e(1:j-1)=Mat(1:j-1,j)
        e(j+1:Num)=Mat(j,j+1:Num)
        Mat(1:j-1,j)=e(1:j-1)-s*(Mat(1:j-1,k)+Tau*e(1:j-1))
        Mat(1:j-1,k)=Mat(1:j-1,k)+s*(e(1:j-1)-Tau*Mat(1:j-1,k))
        Mat(j,j+1:k-1)=e(j+1:k-1)-s*(Mat(j+1:k-1,k)+Tau*e(j+1:k-1))
        Mat(j+1:k-1,k)=Mat(j+1:k-1,k)+s*(e(j+1:k-1)-Tau*Mat(j+1:k-1,k))
        Mat(j,k+1:Num)=e(k+1:Num)-s*(Mat(k,k+1:Num)+Tau*e(k+1:Num))
        Mat(k,k+1:Num)=Mat(k,k+1:Num)+s*(e(k+1:Num)-Tau*Mat(k,k+1:Num))
        Mat(j,k)=0.0D0
        e(:)=Vec(:,j)
        Vec(:,j)=e(:)-s*(Vec(:,k)+Tau*e(:))
        Vec(:,k)=Vec(:,k)+s*(e(:)-Tau*Vec(:,k))
      END DO
    END DO
  END DO

  IF (Suma /= 0.0D0) CALL Mensaje('Diagonalizar',2,.TRUE.)

  !Se regenera la matriz inicial
  DO i=1,Num-1
    DO j=i+1,Num
      Mat(i,j)=Mat(j,i)
    END DO
  END DO

  !Se ordenan los vectores y valores propios
  IF (Ord == 1) THEN
    DO i=1,Num-1
      j=SUM(MINLOC(Val(i:Num)))+i-1
      IF (j /= i) THEN
        CALL Intercambiar(Vec(:,i),Vec(:,j))
        CALL Intercambiar(Val(i),Val(j))
      END IF
    END DO
  END IF

END SUBROUTINE Diagonalizar

!-------------------------------------------------------------------------------
! Resuelve un sistema de ecuaciones lineales por el método de descomposición LU
!-------------------------------------------------------------------------------
! Coef:    Matriz con los coeficientes de las ecuaciones, a la salida contiene
!          la descomposición LU
! Term:    Vector con los términos independientes, a la salida contiene las sol.
! Num:     Número de ecuaciones e incógnitas en el sistema
! Factor:  Vector con los factores para la normalización de las ecuaciones
! Preci:   Precisión para comparar con cero
! p:       Número de la ecuación que hará de pivote
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE ResolverLU(Coef,Term)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: Coef
  DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: Term

  DOUBLE PRECISION, DIMENSION(SIZE(Term,1)) :: Factor
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-15
  INTEGER :: Num,i,p

  Num=SIZE(Term,1)

  !Se calculan los factores de normalización, si hay alguna fila de ceros,
  !la matriz es singular
  Factor(:)=MAXVAL(ABS(Coef(:,:)),DIM=2)
  IF (ANY(Factor(:) < Preci)) CALL Mensaje('ResolverLU',3,.TRUE.)
  Factor(:)=1.0D0/Factor(:)

  !Se localiza la mejor ecuación para hacer de pivote y se pone en el lugar
  !adecuado, también se intercambian los términos independientes
  !Después se modifican los coeficientes
  DO i=1,Num
    p=i-1+SUM(MAXLOC(Factor(i:Num)*ABS(Coef(i:Num,i))))
    IF (p /= i) THEN
      CALL Intercambiar(Coef(i,:),Coef(p,:))
      CALL Intercambiar(Term(i),Term(p))
      Factor(p)=Factor(i)
    END IF
    IF (ABS(Coef(i,i)) < Preci) CALL Mensaje('ResolverLU',3,.TRUE.)
    Coef(i+1:Num,i)=Coef(i+1:Num,i)/Coef(i,i)
    Coef(i+1:Num,i+1:Num)=Coef(i+1:Num,i+1:Num)- &
                          ProductoTens(Coef(i+1:Num,i),Coef(i,i+1:Num))
  END DO

  !Se resuelve el sistema por doble sustitución
  !(primero directa, y luego inversa)
  DO i=2,Num
    Term(i)=Term(i)-DOT_PRODUCT(Term(1:i-1),Coef(i,1:i-1))
  END DO
  Term(Num)=Term(Num)/Coef(Num,Num)
  DO i=Num-1,1,-1
    Term(i)=(Term(i)-DOT_PRODUCT(Term(i+1:Num),Coef(i,i+1:Num)))/Coef(i,i)
  END DO

END SUBROUTINE ResolverLU

!-------------------------------------------------------------------------------
! Esta subrutina calcula la descomposición en valores singulares y la
! pseudoinversa de una matriz rectangular dada
!    A(m,n)=U(m,n)*W(n,n)*V'(n,n)    W->diagonal
!   A+(n,m)=V(n,n)*1/W(n,n)*U'(m,n)
!-------------------------------------------------------------------------------
! Matriz:                 Matriz a descomponer, en la salida contiene U
! Diag:                   Elementos diagonales de W (valores singulares)
! Vect:                   Matriz V
! Inv:                    Pseudoinversa de la matriz inicial
! Preci:                  Precisión para los valores singulares
! Aux1,Aux2,RVec:         Vectores auxiliares
! Norm,c,e,f,g,h,s,x,y,z: Variables auxiliares
! Dim1,Dim2:              Dimensiones de la matriz a descomponer
! Iter:                   Número de iteraciones del algoritmo
! i,j,k,l,m:              Contadores
!-------------------------------------------------------------------------------
SUBROUTINE Pseudoinversa(Matriz,Diag,Vect,Inv)
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: Matriz
  DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Diag
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: Vect
  DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: Inv

  DOUBLE PRECISION, DIMENSION(SIZE(Matriz,1)) :: Aux1
  DOUBLE PRECISION, DIMENSION(SIZE(Matriz,2)) :: Aux2,RVec
  DOUBLE PRECISION, PARAMETER :: Preci=1.0D-12
  DOUBLE PRECISION :: Norm,c,e,f,g,h,s,x,y,z
  INTEGER :: Dim1,Dim2,Iter,i,j,k,l,m

  Dim1=SIZE(Matriz,1)
  Dim2=SIZE(Matriz,2)
  g=0.0D0
  e=0.0D0
  DO i=1,Dim2
    l=i+1
    RVec(i)=e*g
    g=0.0D0
    e=0.0D0
    IF (i <= Dim1) THEN
      e=SUM(ABS(Matriz(i:Dim1,i)))
      IF (e /= 0.0D0) THEN
        Matriz(i:Dim1,i)=Matriz(i:Dim1,i)/e
        s=DOT_PRODUCT(Matriz(i:Dim1,i),Matriz(i:Dim1,i))
        f=Matriz(i,i)
        g=-SIGN(SQRT(s),f)
        h=f*g-s
        Matriz(i,i)=f-g
        Aux2(l:Dim2)=MATMUL(Matriz(i:Dim1,i),Matriz(i:Dim1,l:Dim2))/h
        Matriz(i:Dim1,l:Dim2)=Matriz(i:Dim1,l:Dim2)+ &
                              ProductoTens(Matriz(i:Dim1,i),Aux2(l:Dim2))
        Matriz(i:Dim1,i)=e*Matriz(i:Dim1,i)
      END IF
    END IF
    Diag(i)=e*g
    g=0.0D0
    e=0.0D0
    IF ((i <= Dim1) .AND. (i /= Dim2)) THEN
      e=SUM(ABS(Matriz(i,l:Dim2)))
      IF (e /= 0.0D0) THEN
        Matriz(i,l:Dim2)=Matriz(i,l:Dim2)/e
        s=DOT_PRODUCT(Matriz(i,l:Dim2),Matriz(i,l:Dim2))
        f=Matriz(i,l)
        g=-SIGN(SQRT(s),f)
        h=f*g-s
        Matriz(i,l)=f-g
        RVec(l:Dim2)=Matriz(i,l:Dim2)/h
        Aux1(l:Dim1)=MATMUL(Matriz(l:Dim1,l:Dim2),Matriz(i,l:Dim2))
        Matriz(l:Dim1,l:Dim2)=Matriz(l:Dim1,l:Dim2)+ &
                              ProductoTens(Aux1(l:Dim1),RVec(l:Dim2))
        Matriz(i,l:Dim2)=e*Matriz(i,l:Dim2)
      END IF
    END IF
  END DO
  Norm=MAXVAL(ABS(Diag(:))+ABS(RVec(:)))
  DO i=Dim2,1,-1
    IF (i < Dim2) THEN
      IF (g /= 0.0D0) THEN
        Vect(l:Dim2,i)=(Matriz(i,l:Dim2)/Matriz(i,l))/g
        Aux2(l:Dim2)=MATMUL(Matriz(i,l:Dim2),Vect(l:Dim2,l:Dim2))
        Vect(l:Dim2,l:Dim2)=Vect(l:Dim2,l:Dim2)+ &
                            ProductoTens(Vect(l:Dim2,i),Aux2(l:Dim2))
      END IF
      Vect(i,l:Dim2)=0.0D0
      Vect(l:Dim2,i)=0.0D0
    END IF
    Vect(i,i)=1.0D0
    g=RVec(i)
    l=i
  END DO
  DO i=MIN(Dim1,Dim2),1,-1
    l=i+1
    g=Diag(i)
    Matriz(i,l:Dim2)=0.0D0
    IF (g /= 0.0D0) THEN
      g=1.0D0/g
      Aux2(l:Dim2)=g*MATMUL(Matriz(l:Dim1,i),Matriz(l:Dim1,l:Dim2))/Matriz(i,i)
      Matriz(i:Dim1,l:Dim2)=Matriz(i:Dim1,l:Dim2)+ &
                            ProductoTens(Matriz(i:Dim1,i),Aux2(l:Dim2))
      Matriz(i:Dim1,i)=g*Matriz(i:Dim1,i)
     ELSE
      Matriz(i:Dim1,i)=0.0D0
    END IF
    Matriz(i,i)=Matriz(i,i)+1.0D0
  END DO
  DO k=Dim2,1,-1
    DO Iter=1,30
      DO l=k,1,-1
        m=l-1
        IF (ABS(RVec(l))+Norm == Norm) EXIT
        IF (ABS(Diag(m))+Norm == Norm) THEN
          c=0.0D0
          s=1.0D0
          DO i=l,k
            f=s*RVec(i)
            RVec(i)=c*RVec(i)
            IF (ABS(f)+Norm == Norm) EXIT
            g=Diag(i)
            h=Hipotenusa(f,g)
            Diag(i)=h
            h=1.0D0/h
            c=g*h
            s=-f*h
            Aux1(1:Dim1)=Matriz(1:Dim1,m)
            Matriz(1:Dim1,m)=c*Matriz(1:Dim1,m)+s*Matriz(1:Dim1,i)
            Matriz(1:Dim1,i)=c*Matriz(1:Dim1,i)-s*Aux1(1:Dim1)
          END DO
          EXIT
        END IF
      END DO
      z=Diag(k)
      IF (l == k) THEN
        IF (z < 0.0D0) THEN
          Diag(k)=-z
          Vect(1:Dim2,k)=-Vect(1:Dim2,k)
        END IF
        EXIT
      END IF
      IF (Iter == 30) CALL Mensaje('Pseudoinversa',2,.TRUE.)
      x=Diag(l)
      m=k-1
      y=Diag(m)
      g=RVec(m)
      h=RVec(k)
      f=((y+z)*(y-z)+(g+h)*(g-h))/(2.0D0*h*y)
      g=Hipotenusa(f,1.0D0)
      f=((x+z)*(x-z)+h*((y/(f+SIGN(g,f)))-h))/x
      c=1.0D0
      s=1.0D0
      DO j=l,m
        i=j+1
        g=RVec(i)
        y=Diag(i)
        h=s*g
        g=c*g
        z=Hipotenusa(f,h)
        RVec(j)=z
        c=f/z
        s=h/z
        f=c*x+s*g
        g=c*g-s*x
        h=s*y
        y=c*y
        Aux2(1:Dim2)=Vect(1:Dim2,j)
        Vect(1:Dim2,j)=c*Vect(1:Dim2,j)+s*Vect(1:Dim2,i)
        Vect(1:Dim2,i)=c*Vect(1:Dim2,i)-s*Aux2(1:Dim2)
        z=Hipotenusa(f,h)
        Diag(j)=z
        IF (z /= 0.0D0) THEN
          z=1.0D0/z
          c=f*z
          s=h*z
        END IF
        f=c*g+s*y
        x=c*y-s*g
        Aux1(1:Dim1)=Matriz(1:Dim1,j)
        Matriz(1:Dim1,j)=c*Matriz(1:Dim1,j)+s*Matriz(1:Dim1,i)
        Matriz(1:Dim1,i)=c*Matriz(1:Dim1,i)-s*Aux1(1:Dim1)
      END DO
      RVec(l)=0.0D0
      RVec(k)=f
      Diag(k)=x
    END DO
  END DO

  !Se reordenan los valores y los vectores singulares
  DO i=1,Dim2-1
    j=SUM(MAXLOC(Diag(i:Dim2)))+i-1
    IF (j /= i) THEN
      CALL Intercambiar(Matriz(:,i),Matriz(:,j))
      CALL Intercambiar(Vect(:,i),Vect(:,j))
      CALL Intercambiar(Diag(i),Diag(j))
    END IF
  END DO

  !Se anulan los valores pequeños y se calcula la matriz pseudoinversa
  x=MAXVAL(Diag(:))*Preci
  WHERE(Diag <= x) Diag=0.0D0
  DO i=1,Dim2
    IF (Diag(i) /= 0.0D0) THEN
      Inv(i,:)=Matriz(:,i)/Diag(i)
     ELSE
      Inv(i,:)=0.0D0
    END IF
  END DO
  Inv(:,:)=MATMUL(Vect(:,:),Inv(:,:))

END SUBROUTINE Pseudoinversa

END MODULE Utilidades
