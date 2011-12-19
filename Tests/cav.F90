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
  USE Parametros
  USE Entrada
  USE Sistema
  USE Moldy
  USE Cavidad
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Esf
  DOUBLE PRECISION, DIMENSION(3) :: Aux
  DOUBLE PRECISION :: R
  INTEGER :: NV,U,m,i,j,k

  CALL LeerEntrada(5)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaMM))
  CALL LeerControlMoldy(U)
  CLOSE(U)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(MoldyInput))
  CALL LeerSistemaMoldy(U)
  CLOSE(U)

  CALL OrientarMolecula(Soluto,1)

  IF (TipoCavidad == 0) THEN
    ALLOCATE(Esf(1,4))
    Esf(1,:)=(/0.0D0,0.0D0,0.0D0,RadioCavidad/)
   ELSE
    ALLOCATE(Esf(SIZE(Soluto,1),4))
    Esf(:,4)=1.0D0
    DO i=1,SIZE(Esf,1)
      Esf(i,1:3)=Soluto(i)%pos(:)
      IF (TipoCavidad == 1) Esf(i,4)=RadiosVdW(Soluto(i)%z)
    END DO
    Esf(:,4)=RadioCavidad*Esf(:,4)
  END IF

  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4),RadioDisolvente,Subdivisiones)

  WRITE(6,*) '/*'
  !Formato PDB
  DO i=1,SIZE(CavVert,1)
    WRITE(6,100) i,'HVER',1,CavVert(i,:),1.0D0,0.0D0,0
  END DO
  NV=SIZE(CavVert,1)
  DO i=1,SIZE(CavCent,1)
    WRITE(6,100) NV+i,'HCEN',1,CavCent(i,:),1.0D0,0.0D0,0
  END DO
  DO i=1,SIZE(CavTrian,1)
    WRITE(6,101) CavTrian(i,1),NV+i
    WRITE(6,101) CavTrian(i,2),NV+i
    WRITE(6,101) CavTrian(i,3),NV+i
    IF (CavTrian(i,1) < CavTrian(i,2)) WRITE(6,101) CavTrian(i,1),CavTrian(i,2)
    IF (CavTrian(i,2) < CavTrian(i,3)) WRITE(6,101) CavTrian(i,2),CavTrian(i,3)
    IF (CavTrian(i,3) < CavTrian(i,1)) WRITE(6,101) CavTrian(i,3),CavTrian(i,1)
  END DO
  WRITE(6,'(A)') 'END'
  WRITE(6,*) '*/'

  R=0.0D0
  DO i=1,SIZE(Esferas,1)
    R=MAX(R,Norma(Esferas(i,1:3))+Esferas(i,4))
  END DO

  !Formato POV-Ray
  WRITE(6,*) '//POV-Ray'
  WRITE(6,*) '#include "transforms.inc"'
  WRITE(6,*) 'camera { location -30*z look_at 0 }'
  WRITE(6,*) 'light_source { 100*<1,1,-1>, 1 }'
  WRITE(6,*) 'light_source { 100*<-1,1,-1>, 0.5 }'
  WRITE(6,*) 'merge {'
  DO i=1,SIZE(Esferas,1)
    WRITE(6,*)   '  sphere {'
    WRITE(6,102) '    ',Esferas(i,1:3),','
    WRITE(6,*)   '    ',Esferas(i,4)
    WRITE(6,*)   '    pigment { color rgbt <0,0,1,0> }'
    WRITE(6,*)   '  }'
    DO j=i+1,SIZE(Esferas,1)
      m=Int2(i,j)
      IF (m == 0) CYCLE
      Aux(:)=Normalizar(Esferas(i,1:3)-Esferas(j,1:3))
      WRITE(6,*)   '  difference {'
      WRITE(6,*)   '    sphere { '
      WRITE(6,102) '      ',Esferas2(m,1:3)-Toros(m,1:3),','
      WRITE(6,*)   '      ',Esferas2(m,4)
      WRITE(6,*)   '      pigment { color rgbt <1,1,1,0> }'
      WRITE(6,*)   '    }'
      WRITE(6,*)   '    torus{'
      WRITE(6,*)   '      ',Toros(m,4),',',RDis
      WRITE(6,*)   '      Reorient_Trans(y,'
      WRITE(6,102) '        ',Aux(:),')'
      WRITE(6,*)   '      pigment { color rgbt <0,1,0,0> }'
      WRITE(6,*)   '    }'
      IF (Toros(m,4) < RDis) THEN
        WRITE(6,*) '    sphere { 0,',SQRT(RDis*RDis-Toros(m,4)*Toros(m,4)),' }'
      END IF
      WRITE(6,*)   '    translate'
      WRITE(6,102) '      ',Toros(m,1:3)
      WRITE(6,*)   '  }'
      DO k=j+1,SIZE(Esferas,1)
        m=Int3(i,j,k)
        IF (m == 0) CYCLE
        Aux(:)=0.5D0*(0.8D0*Esferas3(m,1:3)+0.2D0*Esferas3(m,4:6))
        WRITE(6,*)   '  difference {'
        WRITE(6,*)   '    mesh {'
        WRITE(6,*)   '      triangle {'
        WRITE(6,102) '        ',Esferas3(m,1:3),','
        WRITE(6,102) '        ',Esferas(j,1:3),','
        WRITE(6,102) '        ',Esferas(k,1:3)
        WRITE(6,*)   '      }'
        WRITE(6,*)   '      triangle {'
        WRITE(6,102) '        ',Esferas(i,1:3),','
        WRITE(6,102) '        ',Esferas3(m,1:3),','
        WRITE(6,102) '        ',Esferas(k,1:3)
        WRITE(6,*)   '      }'
        WRITE(6,*)   '      triangle {'
        WRITE(6,102) '        ',Esferas(i,1:3),','
        WRITE(6,102) '        ',Esferas(j,1:3),','
        WRITE(6,102) '        ',Esferas3(m,1:3)
        WRITE(6,*)   '      }'
        WRITE(6,*)   '      triangle {'
        WRITE(6,102) '        ',Esferas(k,1:3),','
        WRITE(6,102) '        ',Esferas(j,1:3),','
        WRITE(6,102) '        ',Esferas3(m,4:6)
        WRITE(6,*)   '      }'
        WRITE(6,*)   '      triangle {'
        WRITE(6,102) '        ',Esferas(k,1:3),','
        WRITE(6,102) '        ',Esferas3(m,4:6),','
        WRITE(6,102) '        ',Esferas(i,1:3)
        WRITE(6,*)   '      }'
        WRITE(6,*)   '      triangle {'
        WRITE(6,102) '        ',Esferas3(m,4:6),','
        WRITE(6,102) '        ',Esferas(j,1:3),','
        WRITE(6,102) '        ',Esferas(i,1:3)
        WRITE(6,*)   '      }'
        WRITE(6,*)   '      inside_vector'
        WRITE(6,102) '        ',Aux(:)
        WRITE(6,*)   '      pigment { color rgbt <1,1,1,0> }'
        WRITE(6,*)   '    }'
        WRITE(6,*)   '    sphere {'
        WRITE(6,102) '      ',Esferas3(m,1:3),','
        WRITE(6,*)   '      ',RDis
        WRITE(6,*)   '      pigment { color rgbt <1,0,0,0> }'
        WRITE(6,*)   '    }'
        WRITE(6,*)   '    sphere {'
        WRITE(6,102) '      ',Esferas3(m,4:6),','
        WRITE(6,*)   '      ',RDis
        WRITE(6,*)   '      pigment { color rgbt <1,0,0,0> }'
        WRITE(6,*)   '    }'
        WRITE(6,*)   '  }'
      END DO
    END DO
  END DO
  WRITE(6,*) '  bounded_by { sphere { 0,',R,' } }'
  WRITE(6,*) '  no_shadow'
  WRITE(6,*) '}'
  WRITE(6,*) 'mesh{'
  DO i=1,SIZE(CavTrian,1)
    WRITE(6,*)   '  triangle{'
    WRITE(6,102) '    ',CavVert(CavTrian(i,1),1),CavVert(CavTrian(i,1),2), &
                        CavVert(CavTrian(i,1),3),','
    WRITE(6,102) '    ',CavVert(CavTrian(i,2),1),CavVert(CavTrian(i,2),2), &
                        CavVert(CavTrian(i,2),3),','
    WRITE(6,102) '    ',CavVert(CavTrian(i,3),1),CavVert(CavTrian(i,3),2), &
                        CavVert(CavTrian(i,3),3)
    WRITE(6,*)   '  }'
  END DO
  WRITE(6,*) '  pigment { color rgbt <1,1,1,0> }'
  WRITE(6,*) '  no_shadow'
  WRITE(6,*) '}'

100 FORMAT ('HETATM',I5,1X,A4,' ','CAV',1X,'X',I4,' ',3X,3F8.3,2F6.2,6X,'    ',' H',I2)
101 FORMAT ('CONECT',11I5)
102 FORMAT (' ',A,'<',F10.6,',',F10.6,',',F10.6,'>',A)

END PROGRAM prueba
