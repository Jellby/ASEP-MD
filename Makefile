include Makefile.cfg

COMP   = $(FC) $(FFLAGS)
RMF    = rm -f
TAR    = tar -cjf
BINDIR = ./bin

TESTS = Tests

LANG       = espanol.mod english.mod

CAV        = Cavidad.o
CAV_MOD    = cavidad.mod

COMUN      = Comun.o
COMUN_MOD  = $(LANG) datosqm.mod datosqmmm.mod sistema.mod unidades.mod tipoatomo.mod

CONF       = Configuraciones.o
CONF_MOD   = configuraciones.mod

COORD      = Coordenadas.o
COORD_MOD  = coordenadas.mod

DCD        = DCD.o
DCD_MOD    = DCD.mod

EJEC       = Ejecutar.o
EJEC_MOD   = ejecutar.mod

ENT        = Entrada.o
ENT_MOD    = entrada.mod parametros.mod

GAUSS      = Gaussian.o
GAUSS_MOD  = gaussian.mod

GENMM      = GenericoMM.o
GENMM_MOD  = genericomm.mod

GENQM      = GenericoQM.o
GENQM_MOD  = genericoqm.mod

MALLA      = Malla.o
MALLA_MOD  = malla.mod

MOLCAS     = Molcas.o
MOLCAS_MOD = molcas.mod

MOLDY      = Moldy.o
MOLDY_MOD  = moldy.mod

OPTIM      = Optimizacion.o
OPTIM_MOD  = optimizacion.mod

INTER      = Interseccion.o
INTER_MOD  = Interseccion.mod

UTIL       = Utilidades.o
UTIL_MOD   = utilidades.mod

UTFIS      = UtilidadesFis.o
UTFIS_MOD  = utilidadesfis.mod

CUANT  = $(GENQM_MOD) $(GAUSS_MOD) $(MOLCAS_MOD)

OBJS  = $(CAV) $(COMUN) $(CONF) $(COORD) $(DCD) $(EJEC) $(ENT) $(GAUSS) $(GENMM) $(GENQM) $(MALLA) $(MOLCAS) $(MOLDY) $(OPTIM) $(INTER) $(UTIL) $(UTFIS)

PROGS = asepmd

#==============================================================================

.PHONY: all tests pack install clean

all: $(OBJS) asepmd tests

$(CAV):       Cavidad.F90 unidades.mod utilidades.mod
	$(COMP) -c $<
$(CAV_MOD): | Cavidad.F90 unidades.mod utilidades.mod
	$(COMP) -c $(firstword $|)

$(COMUN):       Comun.F90
	$(COMP) -c $<
$(COMUN_MOD): | Comun.F90
	$(COMP) -c $(firstword $|)

$(CONF):        Configuraciones.F90 cavidad.mod parametros.mod sistema.mod datosqm.mod datosqmmm.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $<
$(CONF_MOD):  | Configuraciones.F90 cavidad.mod parametros.mod sistema.mod datosqm.mod datosqmmm.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $(firstword $|)

$(COORD):       Coordenadas.F90 unidades.mod utilidades.mod tipoatomo.mod
	$(COMP) -c $<
$(COORD_MOD): | Coordenadas.F90 unidades.mod utilidades.mod tipoatomo.mod
	$(COMP) -c $(firstword $|)

$(DCD):       DCD.F90 utilidades.mod
	$(COMP) -c $<
$(DCD_MOD): | DCD.F90 utilidades.mod
	$(COMP) -c $(firstword $|)

$(EJEC):       Ejecutar.F90 $(CUANT) genericoqm.mod gaussian.mod molcas.mod genericomm.mod moldy.mod parametros.mod sistema.mod configuraciones.mod utilidades.mod
	$(COMP) -c $<
$(EJEC_MOD): | Ejecutar.F90 $(CUANT) genericoqm.mod gaussian.mod molcas.mod genericomm.mod moldy.mod parametros.mod sistema.mod configuraciones.mod utilidades.mod
	$(COMP) -c $(firstword $|)

$(ENT):     Entrada.F90 $(LANG) utilidades.mod
	$(COMP) -c $<
$(ENT_MOD): Entrada.F90 | $(LANG) utilidades.mod
	$(COMP) -c $<

$(GAUSS):       Gaussian.F90 datosqm.mod parametros.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $<
$(GAUSS_MOD): | Gaussian.F90 datosqm.mod parametros.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $(firstword $|)

$(GENMM):       GenericoMM.F90 parametros.mod sistema.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $<
$(GENMM_MOD): | GenericoMM.F90 parametros.mod sistema.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $(firstword $|)

$(GENQM):       GenericoQM.F90 datosqm.mod parametros.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $<
$(GENQM_MOD): | GenericoQM.F90 datosqm.mod parametros.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $(firstword $|)

$(MALLA):       Malla.F90 parametros.mod sistema.mod tipoatomo.mod unidades.mod utilidades.mod
	$(COMP) -c $<
$(MALLA_MOD): | Malla.F90 parametros.mod sistema.mod tipoatomo.mod unidades.mod utilidades.mod
	$(COMP) -c $(firstword $|)

$(MOLCAS):       Molcas.F90 datosqm.mod parametros.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $<
$(MOLCAS_MOD): | Molcas.F90 datosqm.mod parametros.mod unidades.mod utilidades.mod utilidadesfis.mod
	$(COMP) -c $(firstword $|)

$(MOLDY):       Moldy.F90 configuraciones.mod parametros.mod sistema.mod tipoatomo.mod unidades.mod utilidades.mod
	$(COMP) -c $<
$(MOLDY_MOD): | Moldy.F90 configuraciones.mod parametros.mod sistema.mod tipoatomo.mod unidades.mod utilidades.mod
	$(COMP) -c $(firstword $|)

$(OPTIM):       Optimizacion.F90 $(LANG) coordenadas.mod ejecutar.mod parametros.mod tipoatomo.mod datosqmmm.mod unidades.mod utilidades.mod
	$(COMP) -c $<
$(OPTIM_MOD): | Optimizacion.F90 $(LANG) coordenadas.mod ejecutar.mod parametros.mod tipoatomo.mod datosqmmm.mod unidades.mod utilidades.mod
	$(COMP) -c $(firstword $|)

$(INTER):       Interseccion.F90 $(LANG) coordenadas.mod ejecutar.mod optimizacion.mod parametros.mod tipoatomo.mod unidades.mod utilidades.mod
	$(COMP) -c $<
$(INTER_MOD): | Interseccion.F90 $(LANG) coordenadas.mod ejecutar.mod optimizacion.mod parametros.mod tipoatomo.mod unidades.mod utilidades.mod
	$(COMP) -c $(firstword $|)

$(UTIL):       Utilidades.F90 $(LANG)
	$(COMP) -c $<
$(UTIL_MOD): | Utilidades.F90 $(LANG)
	$(COMP) -c $(firstword $|)

$(UTFIS):       UtilidadesFis.F90 tipoatomo.mod utilidades.mod
	$(COMP) -c $<
$(UTFIS_MOD): | UtilidadesFis.F90 tipoatomo.mod utilidades.mod
	$(COMP) -c $(firstword $|)

asepmd:	asepmd.F90 $(LANG) $(OBJS) parametros.mod entrada.mod moldy.mod malla.mod cavidad.mod configuraciones.mod sistema.mod ejecutar.mod optimizacion.mod utilidades.mod utilidadesfis.mod
	$(COMP) asepmd.F90 $(OBJS) -o asepmd

tests: $(OBJS)
	@echo
	@$(MAKE) -C $(TESTS)

pack:
	$(TAR) asepmd2.tar.bz2 \
	pendiente \
        COPYING README \
	*.F90 \
	Tests/*.F90 \
	Makefile Makefile.cfg Tests/Makefile \
	Tests/*.in \
	Tests/scripts/*.sh Tests/scripts/*.py Tests/scripts/*.pl \
	Tests/coord.modif \
	Tests/test.ctr Tests/test.system Tests/test.g Tests/test.input

install:
	install -d $(BINDIR)
	install $(PROGS) $(BINDIR)
	$(MAKE) -C $(TESTS) install

clean:
	$(RMF) *.o *.mod $(PROGS)
	$(MAKE) -C $(TESTS) clean
