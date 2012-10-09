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

PROGRAM prueba
  USE Parametros
  USE Entrada
  USE Moldy
  USE Malla
  USE Cavidad
  USE Configuraciones
  USE Sistema
  USE Ejecutar
  USE Utilidades
  USE UtilidadesFis
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Potencial,ASEP
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Esf,Pot,CargasAjustadas
  DOUBLE PRECISION, DIMENSION(3) :: ErrAjuste
  INTEGER :: U,Utmp,i,Num,NumCargas,Error

  CALL LeerEntrada(5)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(EntradaMM))
  CALL LeerControlMoldy(U)
  CLOSE(U)

  U=NuevaUnidad()
  OPEN(U,FILE=TRIM(MoldyInput))
  CALL LeerSistemaMoldy(U)
  CLOSE(U)

  CALL OrientarMolecula(Soluto(:),0)
  CALL OrientarMolecula(Disolvente(:),0)
  CALL OrientarMolecula(Disolvente2(:),0)

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

  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4)+3.0D0,RadioDisolvente,Subdivisiones)
  CALL ConstruirCavidad(Esf(:,1:3),Esf(:,4),RadioDisolvente)

  Extension='.mod'
  CALL EjecutarMM(.TRUE.)

  CALL MallaMolecula(Soluto,MallaSoluto)

  U=NuevaUnidad()
  OPEN(U,FILE='malla.xyz')
  WRITE(U,*) SIZE(MallaSoluto,1)
  WRITE(U,*)
  DO i=1,SIZE(MallaSoluto,1)
    WRITE(U,*) 'X',MallaSoluto(i,:)/AngstromAtomica
  END DO
  CLOSE(U)

  ALLOCATE(Potencial(SIZE(MallaSoluto,1)),ASEP(SIZE(MallaSoluto,1)))

  Utmp=NuevaUnidad()
  OPEN(Utmp,STATUS='SCRATCH',ACTION='READWRITE',FORM='UNFORMATTED')

  U=NuevaUnidad()
  OPEN(U,FILE='configs.xyz')
  ASEP(:)=0.0D0
  Num=0
  NumCargas=0
  DO
    CALL LeerConfig(Error)
    IF (Error /= 0) EXIT
    CALL PotencialMM(MallaSoluto,Potencial)
    ASEP(:)=Asep(:)+Potencial(:)
    CALL SeleccionarMols(Utmp,NumCargas)
    Num=Num+1
    WRITE(6,101) Num
  END DO
  CLOSE(U)

  CLOSE(UConf)

  CALL ReducirCargas(Utmp,NumCargas,CargasDisolvente)
  CargasDisolvente(:,4)=CargasDisolvente(:,4)/DBLE(Num)
  CALL PotencialCargas(CargasDisolvente(:,:),MallaSoluto(:,:),Potencial(:))

  CLOSE(Utmp)

  ALLOCATE(CargasAjustadas(SIZE(CavVert,1),4),Pot(SIZE(ASEP),4))

  ASEP(:)=ASEP(:)/DBLE(Num)
  CargasAjustadas(:,1:3)=CavVert(:,1:3)
  Pot(:,1:3)=MallaSoluto(:,1:3)
  Pot(:,4)=ASEP(:)-Potencial(:)
  CALL AjusteCargas(CavVert(:,1:3),Pot(:,:),0.0D0,CargasAjustadas(:,4),ErrAjuste)

  Pot(:,4)=Potencial(:)
  CALL PotencialCargas(CargasAjustadas(:,:),MallaSoluto(:,:),Potencial(:))
  Potencial(:)=Potencial(:)+Pot(:,4)

  WRITE(6,100)
  WRITE(6,102) NumCargas,SIZE(CargasDisolvente,1)
  WRITE(6,103) ErrAjuste
  WRITE(6,104) MAXVAL(ABS(Potencial(:)-ASEP(:)))

  U=NuevaUnidad()
  OPEN(U,FILE='cargas.xyz')
  WRITE(U,*) SIZE(CargasAjustadas,1)+SIZE(CargasDisolvente,1)
  WRITE(U,*)
  DO i=1,SIZE(CargasAjustadas,1)
    WRITE(U,*) 'Bq',CargasAjustadas(i,1:3)/AngstromAtomica,CargasAjustadas(i,4)
  END DO
  DO i=1,SIZE(CargasDisolvente,1)
    WRITE(U,*) 'X',CargasDisolvente(i,1:3)/AngstromAtomica,CargasDisolvente(i,4)
  END DO
  CLOSE(U)

  DEALLOCATE(Potencial,ASEP,Pot)

100 FORMAT (A)
101 FORMAT ('Configuracion: ',I4)
102 FORMAT ('Cargas iniciales: ',I6,'; Despues de reducir: ',I6)
103 FORMAT ('RErr: ',F12.6,'; Rms: ',F12.6,'; RRms: ',F12.6)
104 FORMAT ('Error maximo en el potencial: ',F12.6)

END PROGRAM prueba
