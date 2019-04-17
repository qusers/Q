################################################################################
#  Q V.5.7 Makefile                                                            #
#  Code authors: Johan Aqvist, Martin Almlof, Martin Ander, Jens Carlson,      #
#  Isabella Feierberg, Peter Hanspers, Anders Kaplan, Karin Kolmodin,          #
#  Petra Wennerstrom, Kajsa Ljunjberg, John Marelius, Martin Nervall,          #
#  Johan Sund, Ake Sandgren, Alexandre Barrozo, Masoud Kazemi, Paul Bauer      #
#  Miha Purg, Irek Szeler, Mauricio Esguerra                                   #
#  Latest update: August 29, 2017                                              #
#  make [option] [COMP=compiler]                                               #
################################################################################

################################################################################
# Environment-specific variables
################################################################################

ifeq ($(OS),Windows_NT)
	OBJ=		obj
	EXE=		.exe
	MDFLAGS=	-p
else
 	OBJ=		o
	MDFLAGS=	-p
endif

################################################################################
# Intel Fortran ifort ( tested in linux and osx)
################################################################################

ifeq ($(COMP),ifort)
	FC=             ifort
	FC_OPT=		-O3
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
# OSX Sierra using gfortran 5.1
################################################################################
else ifeq ($(COMP),osx)
	FC=             gfortran
	FC_OPT=         -O3
	FFLAGS=         -ffree-line-length-none -fcray-pointer \
	                -static-libgfortran -fall-intrinsics -std=legacy \
	                -Wall -Wtabs -fstack-protector ${FC_OPT} -DG95=1

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
	FC=             pgf90$(EXE)
	FC_OPT=         -O0 -g -Mcuda
	FFLAGS=         ${FC_OPT}
	FPP_FLG=        -Mpreprocess
	PGI_TERM=       debug

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
	                -Wall -Wtabs -fstack-protector ${FC_OPT} -DG95=1

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
	                -fall-intrinsics -std=legacy -Wall \
	                ${FC_OPT} -DG95=1

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
all:	qfep qprep qdyn qcalc qdum move1 move2

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
	-rm -f *.$(OBJ) *.mod *.M *.kmo *.il *.oo *.dwf *.pdb

nuke:
	-rm -rf *.$(OBJ) *.mod *.M *.kmo *.il *.oo *.dwf *.pdb qfep$(EXE) qdynp$(EXE) qdyn$(EXE) qprep$(EXE) qcalc \
	qdum$(EXE) ../bin ../obj

print-%  : ; @echo $* = $($*)

qcalc qdyn qdum qdynp qprep qfep: misc.$(OBJ) mpiglob.$(OBJ)

qcalc qdyn qdum qdynp qprep: mask.$(OBJ) prmfile.$(OBJ) sizes.$(OBJ) topo.$(OBJ) trj.$(OBJ) index.$(OBJ)

qcalc qdyn qdum qdynp qfep: nrgy.$(OBJ)

qcalc qdyn qdum qdynp: qatom.$(OBJ)

qcalc qprep qfep: parse.$(OBJ)

qcalc qprep: maskmanip.$(OBJ)

qdyn qdum: qdyn.$(OBJ)

qprep: qprep.$(OBJ) prefs.$(OBJ) prep.$(OBJ) avetr.$(OBJ) nrgy.$(OBJ)
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qfep: qfep.$(OBJ)
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qcalc: calc_base.$(OBJ) calc_chemscore.$(OBJ) calc_fit.$(OBJ) calc_geom.$(OBJ) calc_pmfscore.$(OBJ) \
	calc_com_ke.$(OBJ) calc_com.$(OBJ) calc_rdf.$(OBJ) calc_rms.$(OBJ) calc_rmsf.$(OBJ) \
	calc_entropy.$(OBJ) calc_nb.$(OBJ) calc_xscore.$(OBJ) eigen.$(OBJ) qcalc.$(OBJ)
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn : md.$(OBJ)
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdum : md_dum.$(OBJ)
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdynp : md_mpi.$(OBJ) qdyn_mpi.$(OBJ)
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@


################################################################################
# $(OBJ)ect modules
################################################################################
avetr.$(OBJ): avetr.f90 prep.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_base.$(OBJ):calc_base.f90 topo.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_chemscore.$(OBJ): calc_chemscore.f90 maskmanip.$(OBJ) trj.$(OBJ) prmfile.$(OBJ) index.$(OBJ) qatom.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_entropy.$(OBJ):calc_entropy.f90 calc_base.$(OBJ) maskmanip.$(OBJ) trj.$(OBJ) calc_fit.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_fit.$(OBJ):calc_fit.f90 calc_base.$(OBJ) maskmanip.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_geom.$(OBJ):calc_geom.f90 calc_base.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_nb.$(OBJ):calc_nb.f90 calc_base.$(OBJ) maskmanip.$(OBJ) parse.$(OBJ) qatom.$(OBJ) prmfile.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_pmfscore.$(OBJ): calc_pmfscore.f90 calc_base.$(OBJ) maskmanip.$(OBJ) trj.$(OBJ) topo.$(OBJ) \
	prmfile.$(OBJ) index.$(OBJ) qatom.$(OBJ) misc.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_xscore.$(OBJ): calc_xscore.f90 calc_base.$(OBJ) maskmanip.$(OBJ) trj.$(OBJ) topo.$(OBJ) prmfile.$(OBJ) \
	index.$(OBJ) qatom.$(OBJ) misc.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_rdf.$(OBJ):calc_rdf.f90 calc_base.$(OBJ) parse.$(OBJ) maskmanip.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_rms.$(OBJ):calc_rms.f90 calc_base.$(OBJ) maskmanip.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_com_ke.$(OBJ):calc_com_ke.f90 calc_base.$(OBJ) maskmanip.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_com.$(OBJ):calc_com.f90 calc_base.$(OBJ) maskmanip.$(OBJ)
	${FC} ${FFLAGS} -c $<

calc_rmsf.$(OBJ):calc_rmsf.f90 calc_base.$(OBJ) maskmanip.$(OBJ)
	${FC} ${FFLAGS} -c $<

eigen.$(OBJ):eigen.f90
	${FC} ${FFLAGS} -c $<

index.$(OBJ):index.f90
	${FC} ${FFLAGS} -c $<

mask.$(OBJ):	mask.f90 topo.$(OBJ)
	${FC} ${FFLAGS} -c $<

maskmanip.$(OBJ):maskmanip.f90 mask.$(OBJ) misc.$(OBJ) parse.$(OBJ)
	${FC} ${FFLAGS} -c $<

md.$(OBJ):	md.f90 mpiglob.$(OBJ) qatom.$(OBJ) sizes.$(OBJ) trj.$(OBJ) topo.$(OBJ)
	${FC} ${FFLAGS} ${FPP_FLG} -c md.f90

md_dum.$(OBJ):md.f90 mpiglob.$(OBJ) qatom.$(OBJ) sizes.$(OBJ) trj.$(OBJ) topo.$(OBJ)
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.f90 -o md_dum.$(OBJ)

md_mpi.$(OBJ): md.f90 mpiglob.$(OBJ) qatom.$(OBJ) sizes.$(OBJ) topo.$(OBJ) trj.$(OBJ)
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c md.f90 -o md_mpi.$(OBJ)

misc.$(OBJ): misc.f90 sizes.$(OBJ)
	${FC} ${FFLAGS} -c $<

mpiglob.$(OBJ): mpiglob.f90 sizes.$(OBJ) nrgy.$(OBJ)
	${FC} ${FFLAGS} -c $<

nrgy.$(OBJ): nrgy.f90 sizes.$(OBJ)
	${FC} ${FFLAGS} -c $<

parse.$(OBJ): parse.f90 misc.$(OBJ)
	${FC} ${FFLAGS} -c $<

prefs.$(OBJ): prefs.f90
	${FC} ${FFLAGS} -c $<

prep.$(OBJ): prep.f90 maskmanip.$(OBJ) sizes.$(OBJ) parse.$(OBJ) prmfile.$(OBJ) trj.$(OBJ) index.$(OBJ) prefs.$(OBJ)
	${FC} ${FFLAGS} -c $<

prmfile.$(OBJ): prmfile.f90 misc.$(OBJ) mpiglob.$(OBJ)
	${FC} ${FFLAGS} -g -c $<

qprep.$(OBJ): qprep.f90 prep.$(OBJ) avetr.$(OBJ)
	${FC} ${FFLAGS} -c $<

qatom.$(OBJ): qatom.f90 misc.$(OBJ) nrgy.$(OBJ) prmfile.$(OBJ) sizes.$(OBJ) index.$(OBJ) topo.$(OBJ)
	${FC} ${FFLAGS} -c $<

qcalc.$(OBJ): qcalc.f90 calc_chemscore.$(OBJ) calc_pmfscore.$(OBJ) calc_xscore.$(OBJ) trj.$(OBJ) \
	calc_base.$(OBJ) calc_rms.$(OBJ) calc_fit.$(OBJ) calc_geom.$(OBJ)
	${FC} ${FFLAGS} ${FPP_FLG} -c qcalc.f90

qdyn.$(OBJ): qdyn.f90 md.$(OBJ) mpiglob.$(OBJ)
	${FC} ${FFLAGS} ${FPP_FLG} -c qdyn.f90

qdyn_dum.$(OBJ): qdyn.f90 md.$(OBJ) mpiglob.$(OBJ)
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c qdyn.f90 -o qdyn_dum.$(OBJ)

qdyn_mpi.$(OBJ): qdyn.f90 md_mpi.$(OBJ) mpiglob.$(OBJ)
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
        -c qdyn.f90 -o qdyn_mpi.$(OBJ)

qfep.$(OBJ): qfep.f90 nrgy.$(OBJ) parse.$(OBJ)
	${FC} ${FFLAGS} ${FPP_FLG} -c qfep.f90

sizes.$(OBJ): sizes.f90
	${FC} ${FFLAGS} -c $<

topo.$(OBJ):  topo.f90 misc.$(OBJ) mpiglob.$(OBJ) sizes.$(OBJ)
	${FC} ${FFLAGS} -c $<

trj.$(OBJ):  trj.f90 mask.$(OBJ) misc.$(OBJ)
	${FC} ${FFLAGS} -c $<

move1:
	mkdir $(MDFLAGS) ../obj ; mv *.$(OBJ) *.mod ../obj/

move2:
	mkdir $(MDFLAGS) ../bin ; mv qfep$(EXE) qprep$(EXE) qdyn$(EXE) qdum$(EXE) qcalc$(EXE) ../bin/

move3:
	mkdir $(MDFLAGS) ../bin ; mv qdynp$(EXE) ../bin/

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
