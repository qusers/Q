################################################################################
#  Q V.6.0.1 Makefile                                                          #
#  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      #
#  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          #
#  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          #
#  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer      #
#  Miha Purg, Irek Szeler, Mauricio Esguerra                                   #
#  Latest update: August 29, 2017                                              #
#  make [option] [COMP=compiler]                                               #
################################################################################

################################################################################
# Intel Fortran ifort ( tested in linux and osx)
################################################################################
ifeq ($(COMP),ifort)
	FC=             ifort
	FC_OPT=         -O3
	FFLAGS=         -ip -ipo -unroll -static-intel -fstack-protector \
	                ${FC_OPT}
	FPP_FLG=        -fpp

	MPIFC=          mpiifort
	MPIFLAGS=       -f90=${FC}
	OMPI_FC=        ${FC}
	MPI_FC=         ${FC}
	MPI_FCFLAGS=    ${FFLAGS}
	MPI=            ${MPI_FCFLAGS}

	DEBUG=          -debug -traceback -g  -I${VT_ROOT}/include \
	                -L${VT_LIB_DIR} ${VT_ADD_LIBS} ${FC_OPT}
	DEBUGMPI=       -debug -traceback -g  -I${VT_ROOT}/include \
	                -L${VT_LIB_DIR} ${VT_ADD_LIBS} ${FC_OPT}


################################################################################
# OSX Mojave using gfortran 8.1.0
################################################################################
else ifeq ($(COMP),osx)
	FC=             gfortran
	FC_OPT=         -O3
#	FFLAGS=         -ffree-line-length-none -fcray-pointer \
#	                -static-libgfortran -fall-intrinsics -std=legacy \
#	                -Wall -Wtabs -fstack-protector ${FC_OPT} -DG95=1
	FFLAGS=         -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 \
	                -ffast-math -funroll-loops ${FC_OPT}

	FPP_FLG=        -cpp

	MPIFC=          mpif90
	OMPI_FC=        ${FC}
	OMPI_FCFLAGS=   ${FFLAGS}
	MPI=            ${OMPI_FCFLAGS}

	DEBUG=          -fbacktrace -g ${FC_OPT}
	DEBUGMPI=       -fbacktrace -g ${FC_OPT}


################################################################################
# PGI Compiler -- compiles at csb
################################################################################
else ifeq ($(COMP),pgi)
	FC=             pgf90
	FC_OPT=         -fastsse -Mipa=fast,inline -O3 -Munroll -mcmodel=medium
	FFLAGS=         ${FC_OPT}
	FPP_FLG=        -Mpreprocess

	MPIFC=          mpif90
	MPIFLAGS=       -lmpi  #-f90=${FC}
	MPI=            ${MPIFLAGS} ${FFLAGS}
	MPI_FC=         ${FC}


################################################################################
# IBM Compiler -- untested
################################################################################
else ifeq ($(COMP),ibm)
	FC=             xlf90
	FC_OPT=         -qsmp=omp
	FFLAGS=


################################################################################
# GNU gcc-gfortran
################################################################################
else ifeq ($(COMP),gcc)
	FC=             gfortran
	FC_OPT=         -O3
	FFLAGS=         -funroll-loops -ffree-line-length-none -fcray-pointer \
	                -static-libgfortran -fall-intrinsics -std=legacy \
	                -fstack-protector ${FC_OPT} -DG95=1 -w

	FPP_FLG=        -cpp

	MPIFC=          mpif90
	OMPI_FC=        ${FC}
	OMPI_FCFLAGS=   ${FFLAGS}
	MPI=            ${OMPI_FCFLAGS}

	DEBUG=          -fbacktrace -g ${FC_OPT}	
	DEBUGMPI=       -fbacktrace -g ${FC_OPT}


################################################################################
# If no COMP option is given
################################################################################
else
	FC=             gfortran
	FC_OPT=         -O3
	FFLAGS=         -ffree-line-length-none -fcray-pointer \
	                -fall-intrinsics -std=legacy \
	                ${FC_OPT} -DG95=1 -w

	FPP_FLG=        -cpp

	MPIFC=          mpif90
	OMPI_FC=        ${FC}
	OMPI_FCFLAGS=   ${FFLAGS}
	MPI=            ${OMPI_FCFLAGS}

	DEBUG=          -debug -g ${FC_OPT}
	DEBUGMPI=       -debug -g ${FC_OPT}
endif


################################################################################
# Display options
################################################################################
default: what


################################################################################
# Targets
################################################################################
all:	qfep qprep qdyn qdum qcalc move1 move2

test:   qfep qprep qdyn_test qdum qcalc move1 move4

debug:
	@make FFLAGS="${DEBUG}" \
	qfep qprep qdyn qdum qcalc move1 move2

mpi:
	@make FFLAGS="${MPI}" \
	qdynp move1 move3

mpidebug:
	@make FFLAGS="${DEBUGMPI}" \
	qdynp move1 move3

clean:
	-rm -f *.o  *.mod *.M *.kmo *.il *.oo

nuke:
	-rm -rf *.o  *.mod *.M *.kmo *.il *.oo qfep qdynp qdyn qprep qcalc \
	qdum ../../bin/q6 ../../obj

qcalc qdyn qdyn_test qdum qdynp qprep qfep: misc.o mpiglob.o

qcalc qdyn qdyn_test qdum qdynp qprep: mask.o prmfile.o sizes.o topo.o trj.o index.o

qcalc qdyn qdyn_test qdum qdynp qfep: nrgy.o

qcalc qdyn qdyn_test qdum qdynp: qatom.o

qcalc qprep qfep: parse.o

qcalc qprep: maskmanip.o

qdyn qdum: qdyn.o

qdyn_test: qdyn_test.o

qprep: qprep.o prefs.o prep.o avetr.o nrgy.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qfep: qfep.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qcalc: calc_base.o calc_chemscore.o calc_fit.o calc_geom.o calc_pmfscore.o \
	calc_com_ke.o calc_com.o calc_rdf.o calc_rms.o calc_rmsf.o \
	calc_entropy.o calc_nb.o calc_xscore.o eigen.o qcalc.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn : md.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn_test : md_test.o #qdyn_test.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@	

qdum : md_dum.o
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdynp : md_mpi.o qdyn_mpi.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@


################################################################################
# Object modules
################################################################################
avetr.o: avetr.f90 prep.o
	${FC} ${FFLAGS} -c $<

calc_base.o:calc_base.f90 topo.o
	${FC} ${FFLAGS} -c $<

calc_chemscore.o: calc_chemscore.f90 maskmanip.o trj.o prmfile.o index.o qatom.o
	${FC} ${FFLAGS} -c $<

calc_entropy.o:calc_entropy.f90 calc_base.o maskmanip.o trj.o calc_fit.o
	${FC} ${FFLAGS} -c $<

calc_fit.o:calc_fit.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_geom.o:calc_geom.f90 calc_base.o
	${FC} ${FFLAGS} -c $<

calc_nb.o:calc_nb.f90 calc_base.o maskmanip.o parse.o qatom.o prmfile.o
	${FC} ${FFLAGS} -c $<

calc_pmfscore.o: calc_pmfscore.f90 calc_base.o maskmanip.o trj.o topo.o \
	prmfile.o index.o qatom.o misc.o
	${FC} ${FFLAGS} -c $<

calc_xscore.o: calc_xscore.f90 calc_base.o maskmanip.o trj.o topo.o prmfile.o \
	index.o qatom.o misc.o
	${FC} ${FFLAGS} -c $<

calc_rdf.o:calc_rdf.f90 calc_base.o parse.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_rms.o:calc_rms.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_com_ke.o:calc_com_ke.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_com.o:calc_com.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

calc_rmsf.o:calc_rmsf.f90 calc_base.o maskmanip.o
	${FC} ${FFLAGS} -c $<

eigen.o:eigen.f90
	${FC} ${FFLAGS} -c $<

index.o:index.f90
	${FC} ${FFLAGS} -c $<

mask.o:	mask.f90 topo.o
	${FC} ${FFLAGS} -c $<

maskmanip.o:maskmanip.f90 mask.o misc.o parse.o
	${FC} ${FFLAGS} -c $<

md.o:	md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o
	${FC} ${FFLAGS} ${FPP_FLG} -c md.f90

md_test.o:	md_test.f90 mpiglob.o qatom.o sizes.o trj.o topo.o
	${FC} ${FFLAGS} ${FPP_FLG} -c md_test.f90

md_dum.o:md.f90 mpiglob.o qatom.o sizes.o trj.o topo.o
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.f90 -o md_dum.o

md_mpi.o: md.f90 mpiglob.o qatom.o sizes.o topo.o trj.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c md.f90 -o md_mpi.o

misc.o: misc.f90 sizes.o
	${FC} ${FFLAGS} -c $<

mpiglob.o: mpiglob.f90 sizes.o nrgy.o
	${FC} ${FFLAGS} -c $<

nrgy.o: nrgy.f90 sizes.o
	${FC} ${FFLAGS} -c $<

parse.o: parse.f90 misc.o
	${FC} ${FFLAGS} -c $<

prefs.o: prefs.f90
	${FC} ${FFLAGS} -c $<

prep.o: prep.f90 maskmanip.o sizes.o parse.o prmfile.o trj.o index.o prefs.o
	${FC} ${FFLAGS} -c $<

prmfile.o: prmfile.f90 misc.o mpiglob.o
	${FC} ${FFLAGS} -c $<

qprep.o: qprep.f90 prep.o avetr.o
	${FC} ${FFLAGS} -c $<

qatom.o: qatom.f90 misc.o nrgy.o prmfile.o sizes.o index.o topo.o
	${FC} ${FFLAGS} -c $<

qcalc.o: qcalc.f90 calc_chemscore.o calc_pmfscore.o calc_xscore.o trj.o \
	calc_base.o calc_rms.o calc_fit.o calc_geom.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qcalc.f90

qdyn.o: qdyn.f90 md.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qdyn.f90

qdyn_test.o: qdyn_test.f90 md_test.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qdyn_test.f90	-o qdyn_test.o

qdyn_dum.o: qdyn.f90 md.o mpiglob.o
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c qdyn.f90 -o qdyn_dum.o

qdyn_mpi.o: qdyn.f90 md_mpi.o mpiglob.o
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
        -c qdyn.f90 -o qdyn_mpi.o

qfep.o: qfep.f90 nrgy.o parse.o
	${FC} ${FFLAGS} ${FPP_FLG} -c qfep.f90

sizes.o: sizes.f90
	${FC} ${FFLAGS} -c $<

topo.o:  topo.f90 misc.o mpiglob.o sizes.o
	${FC} ${FFLAGS} -c $<

trj.o:  trj.f90 mask.o misc.o
	${FC} ${FFLAGS} -c $<

move1:
	mkdir -p ../../obj ; mv *.o *.mod ../../obj/

move2:
	mkdir -p ../../bin/q6 ; mv qfep qprep qdyn qdum qcalc ../../bin/q6

move3:
	mkdir -p ../../bin/q6 ; mv qdynp ../../bin/q6

move4:
	mkdir -p ../../bin/q6 ; mv qfep qprep qdyn_test qdum qcalc ../../bin/q6	

################################################################################


################################################################################
# Build instructions
################################################################################
what:
	@echo "Use: make [target] [COMP=compiler]"
	@echo
	@echo "[target] is one of:"
	@echo "all      Everything except the MPI parallel version of Q (qdynp)"
	@echo "debug    All with debug information (stacktraces)."
	@echo "mpi      qdynp using a currently loaded MPI library"
	@echo "mpidebug qdynp with debug information (stacktraces)."
	@echo
	@echo "For compiler/os pick one of:"
	@echo "ifort    Intel Fortran compiler"
	@echo "osx      Gnu Compiler Collection (GCC) in mac"
	@echo "gcc      GCC in linux"
	@echo "pgi      Portland Group compiler"
	@echo
	@echo "Example: make all COMP=ifort"
	@echo
################################################################################
