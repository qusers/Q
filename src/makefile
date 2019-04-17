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
	FC=             pgf90.exe
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
all:	qfep qprep qdyn qcalc move1 move2

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
	-rm -f *.obj  *.mod *.M *.kmo *.il *.oo

nuke:
	-rm -rf *.obj  *.mod *.M *.kmo *.il *.oo qfep qdynp qdyn qprep qcalc \
	qdum ../bin ../obj

qcalc qdyn qdum qdynp qprep qfep: misc.obj mpiglob.obj

qcalc qdyn qdum qdynp qprep: mask.obj prmfile.obj sizes.obj topo.obj trj.obj index.obj

qcalc qdyn qdum qdynp qfep: nrgy.obj

qcalc qdyn qdum qdynp: qatom.obj

qcalc qprep qfep: parse.obj

qcalc qprep: maskmanip.obj

qdyn qdum: qdyn.obj

qprep: qprep.obj prefs.obj prep.obj avetr.obj nrgy.obj
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qfep: qfep.obj
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qcalc: calc_base.obj calc_chemscore.obj calc_fit.obj calc_geom.obj calc_pmfscore.obj \
	calc_com_ke.obj calc_com.obj calc_rdf.obj calc_rms.obj calc_rmsf.obj \
	calc_entropy.obj calc_nb.obj calc_xscore.obj eigen.obj qcalc.obj
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdyn : md.obj
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdum : md_dum.obj
	${FC} ${FFLAGS} ${FLIBS} $+ -o $@

qdynp : md_mpi.obj qdyn_mpi.obj
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FLIBS} $+ -o $@


################################################################################
# Object modules
################################################################################
avetr.obj: avetr.f90 prep.obj
	${FC} ${FFLAGS} -c $<

calc_base.obj:calc_base.f90 topo.obj
	${FC} ${FFLAGS} -c $<

calc_chemscore.obj: calc_chemscore.f90 maskmanip.obj trj.obj prmfile.obj index.obj qatom.obj
	${FC} ${FFLAGS} -c $<

calc_entropy.obj:calc_entropy.f90 calc_base.obj maskmanip.obj trj.obj calc_fit.obj
	${FC} ${FFLAGS} -c $<

calc_fit.obj:calc_fit.f90 calc_base.obj maskmanip.obj
	${FC} ${FFLAGS} -c $<

calc_geom.obj:calc_geom.f90 calc_base.obj
	${FC} ${FFLAGS} -c $<

calc_nb.obj:calc_nb.f90 calc_base.obj maskmanip.obj parse.obj qatom.obj prmfile.obj
	${FC} ${FFLAGS} -c $<

calc_pmfscore.obj: calc_pmfscore.f90 calc_base.obj maskmanip.obj trj.obj topo.obj \
	prmfile.obj index.obj qatom.obj misc.obj
	${FC} ${FFLAGS} -c $<

calc_xscore.obj: calc_xscore.f90 calc_base.obj maskmanip.obj trj.obj topo.obj prmfile.obj \
	index.obj qatom.obj misc.obj
	${FC} ${FFLAGS} -c $<

calc_rdf.obj:calc_rdf.f90 calc_base.obj parse.obj maskmanip.obj
	${FC} ${FFLAGS} -c $<

calc_rms.obj:calc_rms.f90 calc_base.obj maskmanip.obj
	${FC} ${FFLAGS} -c $<

calc_com_ke.obj:calc_com_ke.f90 calc_base.obj maskmanip.obj
	${FC} ${FFLAGS} -c $<

calc_com.obj:calc_com.f90 calc_base.obj maskmanip.obj
	${FC} ${FFLAGS} -c $<

calc_rmsf.obj:calc_rmsf.f90 calc_base.obj maskmanip.obj
	${FC} ${FFLAGS} -c $<

eigen.obj:eigen.f90
	${FC} ${FFLAGS} -c $<

index.obj:index.f90
	${FC} ${FFLAGS} -c $<

mask.obj:	mask.f90 topo.obj
	${FC} ${FFLAGS} -c $<

maskmanip.obj:maskmanip.f90 mask.obj misc.obj parse.obj
	${FC} ${FFLAGS} -c $<

md.obj:	md.f90 mpiglob.obj qatom.obj sizes.obj trj.obj topo.obj
	${FC} ${FFLAGS} ${FPP_FLG} -c md.f90

md_dum.obj:md.f90 mpiglob.obj qatom.obj sizes.obj trj.obj topo.obj
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c md.f90 -o md_dum.obj

md_mpi.obj: md.f90 mpiglob.obj qatom.obj sizes.obj topo.obj trj.obj
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
	-c md.f90 -o md_mpi.obj

misc.obj: misc.f90 sizes.obj
	${FC} ${FFLAGS} -c $<

mpiglob.obj: mpiglob.f90 sizes.obj nrgy.obj
	${FC} ${FFLAGS} -c $<

nrgy.obj: nrgy.f90 sizes.obj
	${FC} ${FFLAGS} -c $<

parse.obj: parse.f90 misc.obj
	${FC} ${FFLAGS} -c $<

prefs.obj: prefs.f90
	${FC} ${FFLAGS} -c $<

prep.obj: prep.f90 maskmanip.obj sizes.obj parse.obj prmfile.obj trj.obj index.obj prefs.obj
	${FC} ${FFLAGS} -c $<

prmfile.obj: prmfile.f90 misc.obj mpiglob.obj
	${FC} ${FFLAGS} -g -c $<

qprep.obj: qprep.f90 prep.obj avetr.obj
	${FC} ${FFLAGS} -c $<

qatom.obj: qatom.f90 misc.obj nrgy.obj prmfile.obj sizes.obj index.obj topo.obj
	${FC} ${FFLAGS} -c $<

qcalc.obj: qcalc.f90 calc_chemscore.obj calc_pmfscore.obj calc_xscore.obj trj.obj \
	calc_base.obj calc_rms.obj calc_fit.obj calc_geom.obj
	${FC} ${FFLAGS} ${FPP_FLG} -c qcalc.f90

qdyn.obj: qdyn.f90 md.obj mpiglob.obj
	${FC} ${FFLAGS} ${FPP_FLG} -c qdyn.f90

qdyn_dum.obj: qdyn.f90 md.obj mpiglob.obj
	${FC} ${FFLAGS} ${FPP_FLG} -DDUM -c qdyn.f90 -o qdyn_dum.obj

qdyn_mpi.obj: qdyn.f90 md_mpi.obj mpiglob.obj
	${MPIFC} ${MPIFLAGS} ${FFLAGS} ${FPP_FLG} -DUSE_MPI \
        -c qdyn.f90 -o qdyn_mpi.obj

qfep.obj: qfep.f90 nrgy.obj parse.obj
	${FC} ${FFLAGS} ${FPP_FLG} -c qfep.f90

sizes.obj: sizes.f90
	${FC} ${FFLAGS} -c $<

topo.obj:  topo.f90 misc.obj mpiglob.obj sizes.obj
	${FC} ${FFLAGS} -c $<

trj.obj:  trj.f90 mask.obj misc.obj
	${FC} ${FFLAGS} -c $<

move1:
	mkdir -p ..\obj ; mv *.obj *.mod ..\obj\

move2:
	mkdir -p ..\bin ; mv qfep.exe qprep.exe qdyn.exe qdum.exe qcalc.exe ..\bin\

move3:
	mkdir -p ..\bin ; mv qdynp.exe ..\bin\

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
