bin_PROGRAMS = $(top_builddir)/bin/h1jet

SOURCES = h1jet.f90\
	  scalar_integrals.f90\
	  ew_parameters.f90\
	  vboson.f90\
	  pdfs_tools.f90\
	  lcl_dec.f90\
	  hboson.f90\
	  input.f90\
	  cross_sections.f90\
	  banner.f90\
	  common_vars.f90\
	  io_utils.f90\
	  gauss_integrator.f90

if USERPATH_SET
## Use user_interface.f90 from USERPATH 
USER_INTERFACE_SOURCE = ${USERPATH}
else
## Use default user_interface.f90 
USER_INTERFACE_SOURCE = user_interface.f90
endif

__top_builddir__bin_h1jet_SOURCES = $(SOURCES)  
EXTRA___top_builddir__bin_h1jet_SOURCES = user_interface.f90 

## LHAPDF_CONFIG = "lhapdf-config"
## HOPPET_CONFIG = "hoppet-config"

## AM_FCFLAGS = -fPIC `$(HOPPET_CONFIG) --fflags`
## AM_FFLAGS = -fPIC `$(HOPPET_CONFIG) --fflags` 

AM_FCFLAGS = -fPIC ${HOPPET_FCFLAGS}
AM_FFLAGS = -fPIC ${HOPPET_FCFLAGS}

## h1jet_LDFLAGS = `$(HOPPET_CONFIG) --ldflags` `$(LHAPDF_CONFIG) --ldflags`  
## h1jet_LDADD = `$(HOPPET_CONFIG) --libs` `$(LHAPDF_CONFIG) --libs` -lchaplin

__top_builddir__bin_h1jet_LDADD = ${HOPPET_LIBS} ${LHAPDF_LIBS} -lchaplin 

__top_builddir__bin_h1jet_LDADD += $(USER_INTERFACE_SOURCE)
__top_builddir__bin_h1jet_DEPENDENCIES = $(USER_INTERFACE_SOURCE)

ew_parameters.o: 
common_vars.o: 
gauss_integrator.o: 
h1jet.o: pdfs_tools.o io_utils.o input.o cross_sections.o banner.o common_vars.o gauss_integrator.o 
hboson.o: scalar_integrals.o
io_utils.o: lcl_dec.o
pdfs_tools.o: ew_parameters.o hboson.o user_interface.o common_vars.o 
scalar_integrals.o:
vboson.o: ew_parameters.o
input.o: ew_parameters.o io_utils.o hboson.o pdfs_tools.o common_vars.o 
user_interface.o: ew_parameters.o common_vars.o  
cross_sections.o: vboson.o io_utils.o hboson.o input.o user_interface.o 
banner.o: user_interface.o common_vars.o 

