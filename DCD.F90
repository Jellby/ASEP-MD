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

MODULE DCD

#ifndef LLL
#define LLL 256
#endif

!-------------------------------------------------------------------------------
! Tipo especial para una trayectoria DCD
!-------------------------------------------------------------------------------
! Nombre:   Nombre del fichero asociado a la trayectoria
! Unidad:   Número de la unidad del fichero asociado
! Abierto:  Variable que indica si el fichero está o no abierto
! Head:     Etiqueta de la cabecera del fichero
! Configs:  Número de configuraciones en la trayectoria
! Inicio:   Configuración inicial
! Fin:      Configuración final
! NAtomos:  Número total de átomos
! NFijos:   Número de átomos fijos
! NLibres:  Número de átomos móviles
! Version:  Versión del fichero DCD
! Posicion: Número de configuración en la que se encuentra el fichero
! PosConf:  Número de configuración a la que corresponden las coordenadas
! Delta:    Tamaño del paso temporal
! Cristal:  Variable que indica si hay información de la celda
! Titulo:   Líneas de título del fichero
! Libres:   Índices de los átomos móviles
! Coords:   Coordenadas cartesianas de todos los átomos, correspondientes a la
!           configuración indicada por PosConf
!-------------------------------------------------------------------------------
TYPE TipoDCD
  CHARACTER(LEN=LLL) :: Nombre
  LOGICAL :: Abierto=.FALSE.,Cristal
  INTEGER :: Unidad=0,Configs,Inicio,Fin,NAtomos,NFijos,NLibres,Version, &
             Posicion,PosConf
  DOUBLE PRECISION :: Delta
  CHARACTER(LEN=4) :: Head
  CHARACTER(LEN=80), DIMENSION(:), ALLOCATABLE :: Titulo
  INTEGER, DIMENSION(:), ALLOCATABLE :: Libres
  REAL, DIMENSION(:,:), ALLOCATABLE :: Coords
END TYPE TipoDCD

CONTAINS
!AbrirDCD(Tray)
!ContarConfigs(Tray)
!LeerDCD(Tray,Pos)
!CerrarDCD(Tray)

!-------------------------------------------------------------------------------
! Abre un fichero DCD asociado a una trayectoria y lee los datos iniciales
! (incluida la primera configuración) 
!-------------------------------------------------------------------------------
! Tray:    Trayectoria DCD que contiene los datos
! Control: Vector con las variables de la cabecera
! Delta2:  Tamaño de paso en precisión sencilla
! U:       Unidad del fichero asociado a la trayectoria
! NTit:    Número de líneas de título
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE AbrirDCD(Tray)
  USE Utilidades
  IMPLICIT NONE
  TYPE(TipoDCD), INTENT(INOUT) :: Tray

  INTEGER, DIMENSION(20) :: Control
  REAL :: Delta2
  INTEGER :: U,NTit,i

  IF (Tray%Abierto) CALL CerrarDCD(Tray)

  !Se abre el fichero binario
  !(puede que haya que cambiar entre big- y little-endian)
  U = NuevaUnidad()
  OPEN(U,FILE=TRIM(Tray%Nombre),FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')

  Tray%Unidad=U
  Tray%Abierto=.TRUE.

  IF (ALLOCATED(Tray%Titulo)) DEALLOCATE(Tray%Titulo)
  IF (ALLOCATED(Tray%Libres)) DEALLOCATE(Tray%Libres)
  IF (ALLOCATED(Tray%Coords)) DEALLOCATE(Tray%Coords)

  !Lee los datos de la cabecera y los almacena en su lugar
  READ(U) Tray%Head,Control(:9),Delta2,Control(11:)
  Tray%Configs=Control(1)
  Tray%Inicio=Control(2)
  Tray%Fin=Control(3)
  Tray%NFijos=Control(9)
  Tray%Delta=Delta2
  Tray%Cristal=(Control(11)==1)
  Tray%Version=Control(20)
  !Si la versión es 0, no hay información de la celda y Delta es
  !de doble precisión, hay que volver a leer la línea
  IF (Tray%Version==0) THEN
    BACKSPACE(U)
    READ(U) Tray%Head,Control(:9),Tray%Delta,Control(12:)
    Tray%Cristal=.FALSE.
  END IF

  !Si el fichero no dice cuántas configuraciones hay (moldy), se cuentan
  IF (Tray%Configs == 0) CALL ContarConfigs(Tray)

  !Se leen las líneas de título
  READ(U) NTit
  ALLOCATE(Tray%Titulo(NTit))
  BACKSPACE(U)
  READ(U) NTit,(Tray%Titulo(i),i=1,NTit)

  !Se lee el número total de átomos y, en el caso de que haya átomos fijos,
  !también los índices de los átomos móviles
  READ(U) Tray%NAtomos
  Tray%NLibres=Tray%NAtomos-Tray%NFijos
  ALLOCATE(Tray%Coords(Tray%NAtomos,3),Tray%Libres(Tray%NLibres))
  IF (Tray%NFijos /= 0) READ(U) Tray%Libres(:)

  !Se leen las coordenadas de todos los átomos correspondientes a la primera
  !configuración y se actualizan las posiciones
  IF (Tray%Cristal) READ(U)
  READ(U) Tray%Coords(:,1)
  READ(U) Tray%Coords(:,2)
  READ(U) Tray%Coords(:,3)
  Tray%Posicion=1
  Tray%PosConf=Tray%Posicion

END SUBROUTINE AbrirDCD

!-------------------------------------------------------------------------------
! Subrutina para contar el número de configuraciones en un fichero DCD cuando
! no está especificado en la cabecera
!-------------------------------------------------------------------------------
! Tray:    Trayectoria DCD que contiene los datos
! U:       Unidad del fichero asociado
! Num:     Número de configuraciones
! Error:   Variable para controlar los errores
!-------------------------------------------------------------------------------
SUBROUTINE ContarConfigs(Tray)
  IMPLICIT NONE
  TYPE(TipoDCD), INTENT(INOUT) :: Tray

  INTEGER :: U,Error,Num

  !Si la trayectoria no está abierta, no hay nada que hacer,
  IF (.NOT. Tray%Abierto) RETURN
  U=Tray%Unidad

  !Esta subrutina debe llamarse justo después de leer la primera línea del
  !fichero, por lo tanto hay que saltar dos o tres líneas
  READ(U,IOSTAT=Error)
  READ(U,IOSTAT=Error)
  IF (Tray%NFijos /= 0) READ(U,IOSTAT=Error)

  !Cada configuración ocupa tres o cuatro líneas
  Num=0
  DO WHILE (Error == 0)
    IF (Tray%Cristal) READ(U,IOSTAT=Error)
    READ(U,IOSTAT=Error)
    READ(U,IOSTAT=Error)
    READ(U,IOSTAT=Error)
    Num=Num+1
  END DO

  !Puesto que la última ha dado error, el número de configuraciones
  !es uno menos
  Tray%Configs=Num-1

  !Se rebobina el fichero y se vuelve a saltar la primera línea para dejarlo
  !como estaba
  REWIND(U)
  READ(U)

END SUBROUTINE

!-------------------------------------------------------------------------------
! Esta subrutina lee las coordenadas de una configuración del fichero DCD
! (por defecto, la siguiente)
!-------------------------------------------------------------------------------
! Tray:    Trayectoria DCD que contiene los datos
! Pos:     Número de configuración que se quiere leer
! U:       Unidad del fichero asociado
! X,Y,Z:   Coordenadas de los átomos móviles
! i:       Contador
!-------------------------------------------------------------------------------
SUBROUTINE LeerDCD(Tray,Pos)
  IMPLICIT NONE
  TYPE(TipoDCD), INTENT(INOUT) :: Tray
  INTEGER, INTENT(IN), OPTIONAL :: Pos

  INTEGER :: U,i
  REAL, DIMENSION(:), ALLOCATABLE :: X,Y,Z

  !Si el fichero no está abierto, se abre
  IF (.NOT. Tray%Abierto) THEN
    CALL AbrirDCD(Tray)
    !Si no se ha pedido una configuración concreta, o si se ha pedido la
    !primera, ya no hay más que hacer
    IF (.NOT. PRESENT(Pos)) THEN
      RETURN
     ELSE
      IF (Pos == 1) RETURN
    END IF
  END IF
  U=Tray%Unidad

  IF (PRESENT(Pos)) THEN
    !Si la configuración pedida es anterior en el fichero o mayor que el
    !número existente, no hace nada. Tampoco si ya es la configuración actual
    IF ((Pos < Tray%Posicion) .OR. (Pos > Tray%Configs)) RETURN
    IF (Pos == Tray%PosConf) RETURN
    !En otro caso, salta configuraciones hasta llegar a la anterior
    DO i=Tray%Posicion+1,Pos-1
      IF (Tray%Cristal) READ(U)
      READ(U)
      READ(U)
      READ(U)
      Tray%Posicion=Tray%Posicion+1
    END DO
   ELSE
    !Si no se ha pedido ninguna configuración, se lee la siguiente. Pero no se
    !hace nada si ya estamos en la última
    IF (Tray%Posicion >= Tray%Configs) RETURN
  END IF

  !A continuación se lee la configuración pedida del fichero DCD
  IF (Tray%Cristal) READ(U)
  IF (Tray%NFijos == 0) THEN
    !Si no hay átomos fijos, se leen todas las coordenadas
    READ(U) Tray%Coords(:,1)
    READ(U) Tray%Coords(:,3)
    READ(U) Tray%Coords(:,2)
   ELSE
    !Si hay átomos fijos, se leen sólo las coordenadas de los átomos móviles
    ALLOCATE(X(Tray%NLibres),Y(Tray%NLibres),Z(Tray%NLibres))
    READ(U) X(:)
    READ(U) Y(:)
    READ(U) Z(:)
    !Se asignan las coordenadas de los átomos móviles según sus índices
    DO i=1,SIZE(X)
      Tray%Coords(Tray%Libres(i),1)=X(i)
      Tray%Coords(Tray%Libres(i),2)=Y(i)
      Tray%Coords(Tray%Libres(i),3)=Z(i)
    END DO
  END IF

  !Finalmente se actualiza la posición del fichero y la configuración
  Tray%Posicion=Tray%Posicion+1
  Tray%PosConf=Tray%Posicion

END SUBROUTINE LeerDCD

!-------------------------------------------------------------------------------
! Cierra un fichero DCD
!-------------------------------------------------------------------------------
! Tray:    Trayectoria DCD que contiene los datos
!-------------------------------------------------------------------------------
SUBROUTINE CerrarDCD(Tray)
  IMPLICIT NONE
  TYPE(TipoDCD), INTENT(INOUT) :: Tray

  IF (.NOT. Tray%Abierto) RETURN

  !Cierra el fichero
  CLOSE(Tray%Unidad)
  
  !Establece las variables correspondientes
  Tray%Unidad=0
  Tray%Abierto=.FALSE.

  !Destruye las matrices
  DEALLOCATE(Tray%Titulo,Tray%Libres,Tray%Coords)

END SUBROUTINE CerrarDCD

END MODULE DCD
