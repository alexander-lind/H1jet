# Edit the following two lines if lhapdf-config and hoppet-config are not in your path
LHAPDF_CONFIG="lhapdf-config"
HOPPET_CONFIG="hoppet-config"

# Location where modules will be stored
MODULEPATH = $(PWD)/modules

# Default fortran compiler
FC = gfortran
# store modules
FFLAGS = -O2 -fPIC 
# The option that causes modules to be stored in the $(MODULEPATH) directory
# depends on the compiler. Here it's provided for gfortran and ifort.  
ifeq ("$(FC)","gfortran")	
INCLUDE= -J$(MODULEPATH) -I$(MODULEPATH) `$(HOPPET_CONFIG) --fflags`
else 
ifeq ("$(FC)","ifort")
INCLUDE= -module $(MODULEPATH) -I$(MODULEPATH) `$(HOPPET_CONFIG) --fflags`
else
INCLUDE= -I. `$(HOPPET_CONFIG) --fflags`
endif
endif

LIBS= `$(HOPPET_CONFIG) --ldflags` `$(LHAPDF_CONFIG) --ldflags`

#CHAPLIN= `locate libchaplin.a | tail -1 | sed s/libchaplin.a//`
#CHAPLIN=/its/home/al629/scratch/heptools/chaplin-1.2/chaplin-install/lib
#LDFLAGS = -L$(CHAPLIN) -lchaplin
LDFLAGS = -lchaplin 

# Select between smallR version of the code and the svn one
SOURCEDIR  = $(PWD)/src

VPATH = $(USERPATH):$(SOURCEDIR):$(PWD)/obj

PROG = h1jet 

SRCS = scalar_integrals.f90 dgset.f dgquad.f d107d1.f frcgauss.f90 vboson.f90 pdfs_tools.f90 ew_parameters.f90 io_utils.f90 lcl_dec.f90 hboson.f90 input.f90 cross_sections.f90 user_interface.f90 banner.f90 common_vars.f90 

OBJS = scalar_integrals.o dgset.o dgquad.o d107d1.o frcgauss.o vboson.o pdfs_tools.o ew_parameters.o io_utils.o lcl_dec.o hboson.o input.o cross_sections.o user_interface.o banner.o common_vars.o 

# Trick to enable old 'make PROG=xxx' form to still work
ALL: $(PROG)__

$(PROG)__: $(PROG)

h1jet: h1jet.o $(OBJS) 
	$(FC)  -o $@ obj/$@.o $(patsubst %,obj/%,$(OBJS)) $(LIBS) $(LDFLAGS)

libclean:
	rm -f   obj/*.o  

clean:
	rm -f   obj/*.o  modules/*.mod *.d

distclean: clean
	rm -f  $(PROG)

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90 
	@ mkdir -p obj modules 
	$(FC) $(FFLAGS)  $(INCLUDE) -c -o obj/$@ $<

%.o: %.f
	@ mkdir -p obj modules 
	$(FC) $(FFLAGS) $(INCLUDE) -c -o obj/$@ $<

d107d1.o:
dgquad.o: d107d1.o
dgset.o: dgquad.o d107d1.o
ew_parameters.o: 
common_vars.o: 
frcgauss.o: dgset.o 
h1jet.o: pdfs_tools.o io_utils.o frcgauss.o input.o cross_sections.o banner.o common_vars.o 
hboson.o: scalar_integrals.o
io_utils.o: lcl_dec.o
pdfs_tools.o: ew_parameters.o hboson.o user_interface.o common_vars.o 
scalar_integrals.o:
vboson.o: ew_parameters.o
input.o: ew_parameters.o io_utils.o hboson.o pdfs_tools.o common_vars.o 
user_interface.o: ew_parameters.o common_vars.o  
cross_sections.o: vboson.o io_utils.o hboson.o input.o user_interface.o 
banner.o: user_interface.o common_vars.o 
