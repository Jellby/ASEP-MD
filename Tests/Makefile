include ../Makefile.cfg
FFLAGS+= -I ..
COMP=$(FC) $(FFLAGS)
RMF=rm -f

CAV=../Cavidad.o
CONF=../Configuraciones.o
COMUN=../Comun.o ../Utilidades.o ../UtilidadesFis.o
COORD=../Coordenadas.o
GAUS=../Gaussian.o
GENMM=../GenericoMM.o
GENQM=../GenericoQM.o
MAL=../Malla.o
MOLCAS=../Molcas.o
EJEC=../Ejecutar.o $(GENMM) $(GENQM) $(GAUS) $(MOLCAS) $(MOLDY)
ENT=../Entrada.o
OPT=../Optimizacion.o $(COORD) $(ENT) $(EJEC)
INT=../Interseccion.o $(OPT)
MOLDY=../Moldy.o $(CONF) $(CAV) $(ENT)

TESTS= eq diag inv rfo busq cargas sist calcqm dernum coord opt int cav malla calcmm cargasatom asep asepmd interp elibre configs

.PHONY: all clean

all: $(TESTS)

eq: eq.F90 $(COMUN)
	$(COMP) $^ -o $@

diag: diag.F90 $(COMUN)
	$(COMP) $^ -o $@

inv: inv.F90 $(COMUN)
	$(COMP) $^ -o $@

rfo: rfo.F90 $(OPT) $(COMUN)
	$(COMP) $^ -o $@

busq: busq.F90 $(OPT) $(COMUN)
	$(COMP) $^ -o $@

cargas: cargas.F90 $(COMUN)
	$(COMP) $^ -o $@

sist: sist.F90 $(GENMM) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

calcqm: calcqm.F90 $(MOLDY) $(ENT) $(EJEC) $(COMUN)
	$(COMP) $^ -o $@

dernum: dernum.F90 $(MOLDY) $(ENT) $(EJEC) $(COMUN)
	$(COMP) $^ -o $@

coord: coord.F90 $(COORD) $(MOLDY) $(ENT) $(COMUN)
	$(COMP) $^ -o $@

opt: opt.F90 $(OPT) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

int: int.F90 $(INT) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

cav: cav.F90 $(CAV) $(ENT) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

malla: malla.F90 $(MAL) $(ENT) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

calcmm: calcmm.F90 $(ENT) $(MOLDY) $(MAL) $(CAV) $(EJEC) $(COMUN)
	$(COMP) $^ -o $@

cargasatom: cargasatom.F90 $(MOLDY) $(ENT) $(COMUN)
	$(COMP) $^ -o $@

asep: asep.F90 $(ENT) $(MOLDY) $(MAL) $(CAV) $(EJEC) $(COMUN)
	$(COMP) $^ -o $@

asepmd: asepmd.F90 $(ENT) $(MOLDY) $(MAL) $(CAV) $(EJEC) $(OPT) $(COMUN)
	$(COMP) $^ -o $@

interp: interp.F90 $(MOLDY) $(COORD) $(COMUN)
	$(COMP) $^ -o $@

elibre: elibre.F90 $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

configs: configs.F90 $(GENMM) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

clean:
	$(RMF) $(TESTS)
