################################################################
################################################################
####
#### Makefile for IPDPS SNACS presentation
####
################################################################
################################################################

info::

COURSE = mpi
FILE = snacs_eijkhout
FILES = ${FILE}
EXTRATEXINPUTS = :${HOME}/Software/TeX/LaTeX//

.PHONY: tacc-full online
VERSIONS = tacc-full tacc-advanced online
info ::
	@echo "make ${VERSIONS}"
${VERSIONS} :
	make course 
	cp ${FILE}.pdf $@.pdf
	make --no-print-directory upload FILE=$@

ENGINE = pdflatex -shell-escape
info::
	@echo "make presentation onepass pdf [FILE=... (default: ${FILE})]"
	@if [ ! -z "${FILES}" ] ; then \
	   echo "     available files: ${FILES}" ; \
	 fi

.PHONY: presentation onepass pdf 
onepass :
	TEXINPUTS=.:..//:../../istc//helix-ppse-private//:${TEXINPUTS}${EXTRATEXINPUTS} ${ENGINE} ${FILE}

presentation pdf : 
	@export presentation=${FILE} && export presentation=$${presentation%%.tex} \
	 && make onepass FILE=$${presentation} \
	 && if [ -z "${BIB}" ] ; then \
	        bib=0 ; else bib=1 ; fi \
	 && if [ -z "${INDEX}" ] ; then \
	        index=0 ; else index=1 ; fi \
	 && if [ $${bib} -eq 1 ] ; then \
	        ( bibtex          $${presentation} || echo true ) \
	    ; fi \
	 && if [ $${index} -eq 1 ] ; then \
	        ( makeindex       $${presentation} || echo true ) \
	    ; fi \
	 && if [ $${index} -eq 1 -o $${bib} -eq 1 ] ; then \
	        make onepass FILE=$${presentation} \
	    ; fi \
	 && make upload  FILE=$${presentation}

info ::
	@echo "make upload"
.PHONY: upload
upload ::
	@echo "uploading ${FILE} course" \
	&& if [ -d "${HOME}/Dropbox" ] ; then \
	      cp ${FILE}.pdf ${HOME}/Dropbox/Scicomp/Courses/ \
	   ; fi \
	&& cp ${FILE}.pdf ${HOME}/Current/istc/pdf/pcse

#
# Cleaning
#
info ::
	@echo "make clean"
.PHONY: clean
clean ::
	@/bin/rm -f *~ *.{aux,bbl,blg,idx,ilg,ind,log,nav,out,snm,toc,vrb} \
	    acronyms.tex *.cut
