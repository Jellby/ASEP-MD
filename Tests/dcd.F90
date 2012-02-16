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

PROGRAM prueba
  USE DCD
  IMPLICIT NONE
  TYPE(TipoDCD) :: Trayectoria
  INTEGER :: i

  READ(5,*) Trayectoria%Nombre

  CALL AbrirDCD(Trayectoria)

  WRITE(6,*) 'Head: ',Trayectoria%Head
  WRITE(6,*) 'Version: ',Trayectoria%Version
  WRITE(6,*) 'Configs: ',Trayectoria%Configs
  WRITE(6,*) 'Inicio: ',Trayectoria%Inicio
  WRITE(6,*) 'Fin: ',Trayectoria%Fin
  WRITE(6,*) 'Delta: ',Trayectoria%Delta
  WRITE(6,*) 'NAtomos: ',Trayectoria%NAtomos
  WRITE(6,*) 'NFijos: ',Trayectoria%NFijos
  WRITE(6,*) 'Cristal: ',Trayectoria%Cristal
  WRITE(6,*) 'Titulo: ',Trayectoria%Titulo(1)
  DO i=2,SIZE(Trayectoria%Titulo)
    WRITE(6,*) '        ',Trayectoria%Titulo(i)
  END DO
  WRITE(6,*) 'PosConf: ',Trayectoria%PosConf

!  DO i=1,Trayectoria%Natomos
!    WRITE(6,*) Trayectoria%Coords(i,:)
!  END DO

  DO i=2,Trayectoria%Configs,3
    CALL LeerDCD(Trayectoria,i)
    WRITE(6,*) 'PosConf: ',Trayectoria%PosConf
  END DO

  CALL CerrarDCD(Trayectoria)

END PROGRAM prueba
