EXTRA_DIST= craig_image.pdf header.html footer.html header.tex craig.doxygen

pdf: Makefile
	@echo "#"
	@echo "# pdf"
	@echo "#"
	cd latex &&\
	pdflatex refman.tex &&\
	pdflatex refman.tex
	latex_count=5 ; \
	while egrep -s 'Rerun (LaTeX|to get cross-references right)' latex/refman.log && [ $$latex_count -gt 0 ] ;\
		do \
		echo "Rerunning latex...." ;\
		cd latex && pdflatex refman.tex ;\
		latex_count=`expr $$latex_count - 1` ;\
	done
	mv latex/refman.pdf .
	mv latex/refman.log .
	$(MAKE) pdfclean

rundox: Makefile
	doxygen craig.doxygen 2> doxygen.out

html: Makefile
	tar cvf - html/ | gzip > refman.html.tar.gz
	$(MAKE) htmlclean

clean: pdfclean htmlclean
	rm -f refman.html.tar.gz
	rm -f refman.log
	rm -f refman.pdf.gz 
	rm -f doxygen.out
	rm -f refman.html.tar
	rm -f refman.pdf

pdfclean:
	rm -fr latex

htmlclean:
	rm -fr html
	rm -f doc

doc: rundox html pdf
	gzip -f refman.pdf



