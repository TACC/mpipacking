################################################################
################################################################
####
#### Performance of MPI packing modes
####
################################################################
################################################################

info ::
	@echo "make"

ifdef TACC_TAU_DIR
  CX = tau_cc.sh
  CXX = tau_cxx.sh
else
  CC  = mpicc
  CXX = mpicxx
endif

OPTLEVEL = 2
info ::
	@echo "    [ OPTLEVELS=nnn ] "
FFLAGS := -g -O${OPTLEVEL}
CFLAGS := ${FFLAGS} -std=c99 
CXXFLAGS := ${FFLAGS}

info ::
	@echo "    [ DEBUG=yes ] "
	@echo "    [ MULTI=yes ] "
	@echo "    [ EXTRAOPTIONS=... ]"
	@echo "      supported: -DFLUSH -DL3SIZE=nnn"
DEBUG = no
MULTI = no
YES = yes
NO = no
ifeq  "${DEBUG}" "${YES}"
  CFLAGS += -DDEBUG
  CXXFLAGS += -DDEBUG
  FFLAGS += -DDEBUG
endif
ifeq  "${MULTI}" "${YES}"
  CFLAGS += -DMULTI
  CXXFLAGS += -DMULTI
  FFLAGS += -DMULTI
endif

##
## we generate programs with the hostname attached
##
.PHONY: packperf
info ::
	@echo "make packperf (will create packperf.HOSTNAME)"
HOSTNAME := $(shell hostname -f | cut -d. -f2 )-${TACC_FAMILY_MPI}
packperf : packperf.c
	@echo ; echo "Using CC=`which ${CC}`" ; echo
	${CC} ${CFLAGS} ${EXTRAOPTIONS} \
	    $^ -o $@.${HOSTNAME} -lm
	@echo "made: $@.${HOSTNAME}" && echo
packperfmalloc : packperf.c
	@echo ; echo "Using CC=`which ${CC}`" ; echo
	${CC} ${CFLAGS} ${EXTRAOPTIONS} -DMALLOC \
	    $^ -o $@ -lm
clean ::
	@/bin/rm -f packperf

##
## batch file
##
.PHONY: sbatch
info ::
	@echo 
	@echo "make sbatch: generate batch file"
	@echo "    [ PARTITION=... (default=${PARTITION}) ]"
PARTITION = development
sbatch :
	@./sbatch.make -p ${PARTITION}

##
## generate plot
##
.PHONY: out plot plots
PYTHON = python3
info ::
	@echo "make out : copy slurm log to timings directory"
	@echo "make plot : convert slurm about to pdf plot"
	@echo "    [ HOSTNAME=... (default: ${HOSTNAME}) ]"
	@echo "    [ TAG=... (optional tag for pdf file) ]"
	@echo "    [ PLOT_OPTOINS=... (pass to packplot.py) ]"
	@echo "make plots : plot all in archs"
TAG =
out :
	@lastout=`ls packingtest-${HOSTNAME}.o* | tail -n 1` \
	&& cp $$lastout run_timings/out.${HOSTNAME}${TAG} \
	&& hg add run_timings/out.${HOSTNAME}${TAG}
PLOTHEADER =
plot :
	@cd run_timings \
	&& ${PYTHON} ../packplot.py -h ${HOSTNAME}${TAG} \
	             ${PLOT_OPTIONS} \
	             `if [ ! -z "${PLOTHEADER}" ] ; then \
	                 echo "-H ${PLOTHEADER}" ; fi ` \
	             out.${HOSTNAME}${TAG} \
	&& ${PYTHON} ../packplot.py -h ${HOSTNAME}${TAG} -s 0,1,3,4,7 \
	             ${PLOT_OPTIONS} \
	             -t 1 out.${HOSTNAME}${TAG} \
	&& ${PYTHON} ../packplot.py -h ${HOSTNAME}${TAG} -s 0,1,2,5,6\
	             ${PLOT_OPTIONS} \
	             -t 2 out.${HOSTNAME}${TAG} \
	&& hg add ${HOSTNAME}*.pdf
plots :
	@for a in `cat archs | awk '{print $$1}' ` ; do \
	  echo && echo "Plotting $$a" && echo \
	  && make plot HOSTNAME=$$a \
	               PLOTHEADER=" `awk '/'$$a'/ {print $$2}' archs ` "  \
	  ; done

####
#### report
####
.PHONY: report
report :
	@echo "Please go into the report directory and make there" \
	&& sleep 2 && echo "on second thought" && sleep 1 \
	&& cd report && make pdf FILE=ieee-report

info ::
	@echo "make clean"
	@echo "make arch_clean : architecture speficic cleanup"
.PHONY: arch_clean clean
clean ::
	@/bin/rm -rf *.o *~ *.gch a.out *.dSYM MULTI__* events.* *.runlog texput.* \
	    core.[0-9]* idev[0-9]*.o[0-9]* ddt.o[0-9]* \
	    jobtest.o* packingtest.o* tautrace_*
	@cd report && make --no-print-directory clean
arch_clean :: clean
	@/bin/rm -f packperf.${HOSTNAME} packingtest-${HOSTNAME}.o*
