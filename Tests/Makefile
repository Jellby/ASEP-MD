include ../Makefile.cfg
FFLAGS+= -I ..
COMP=$(FC) $(FFLAGS)
RMF=rm -f
BINDIR_TESTS=../bin/tests
BINDIR_SCRIPTS=../bin/scripts

CAV=../Cavidad.o
CONF=../Configuraciones.o
COMUN=../Comun.o ../Utilidades.o ../UtilidadesFis.o
COORD=../Coordenadas.o
GAUS=../Gaussian.o
GENMM=../GenericoMM.o $(DCD)
GENQM=../GenericoQM.o
MAL=../Malla.o
MOLCAS=../Molcas.o
EJEC=../Ejecutar.o $(GENMM) $(GENQM) $(GAUS) $(MOLCAS) $(MOLDY)
ENT=../Entrada.o
OPT=../Optimizacion.o $(COORD) $(ENT) $(EJEC)
INT=../Interseccion.o $(OPT)
MOLDY=../Moldy.o $(CONF) $(CAV) $(ENT)
DCD=../DCD.o

TESTS= eq diag inv rfo busq cargas sist calcqm dernum coord opt int cav malla calcmm cargasatom asep asepmd interp elibre configs dcd

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

coord: coord.F90 $(COORD) $(GENMM) $(MOLDY) $(ENT) $(COMUN)
	$(COMP) $^ -o $@

opt: opt.F90 $(OPT) $(GENMM) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

int: int.F90 $(INT) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

cav: cav.F90 $(CAV) $(ENT) $(GENMM) $(MOLDY) $(COMUN)
	$(COMP) $^ -o $@

malla: malla.F90 $(MAL) $(ENT) $(GENMM) $(MOLDY) $(COMUN)
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

dcd: dcd.F90 $(DCD) $(COMUN)
	$(COMP) $^ -o $@

install:
	install -d $(BINDIR_TESTS)
	install $(TESTS) $(BINDIR_TESTS)
	install -d $(BINDIR_SCRIPTS)
	install scripts/* $(BINDIR_SCRIPTS)

clean:
	$(RMF) $(TESTS) functions.mod
