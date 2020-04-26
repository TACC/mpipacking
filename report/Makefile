################################################################
################################################################
####
#### This makefile is part of the source of 
#### "Parallel Programming in MPI and OpenMP"
#### by Victor Eijkhout, copyright 2018
####
#### driver makefile for technical reports
####
################################################################
################################################################

.PHONY: info
info ::
	@echo "make format report [FILE=... (default: ${FILE})]"

FILE = packing
FILES = cpe-report ieee-report packing report-body scimacs wilmacs

.PHONY: onepass format pdf report
onepass format :
	TEXINPUTS=${TEXINPUTS}:../run_timings pdflatex ${FILE}
pdf report :
	make onepass FILE=${FILE}
	-bibtex ${FILE}
	make onepass FILE=${FILE}
	make onepass FILE=${FILE}
	@echo ".. Copying to dropbox"
	cp ${FILE}.pdf ${HOME}/Dropbox/Scicomp/performance/

clean ::
	@rm -f *.{a,aux,blg,idx,ilg,ind,log,mod,o,out,svn,toc,xwm} *~
	@for f in ${FILES} ; do \
	    rm -f $${f}.{pdf,bbl} ; \
	 done
	@rm -f acronyms.tex 
