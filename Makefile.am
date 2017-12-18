# adding doxygen support

SUBDIRS = models python/bin python/lib perl/bin perl/lib regex sparsevec lless src tests doc datasets

AUX_DIST =  @top_srcdir@/$(ac_aux_dir)/config.guess @top_srcdir@/$(ac_aux_dir)/config.sub @top_srcdir@/$(ac_aux_dir)/install-sh @top_srcdir@/$(ac_aux_dir)/ltmain.sh @top_srcdir@/$(ac_aux_dir)/missing

EXTRA_DIST= @top_srcdir@/ACKNOWLEDGEMENTS @top_srcdir@/autogen.sh @top_srcdir@/config @top_srcdir@/README.md @top_srcdir@/README.extended

MAINTAINERCLEANFILES  = Makefile.in aclocal.m4 configure config.in stamp-h1 $(AUX_DIST)

doc: doxygen-docs

doxygen-docs:
	cd doc && make doc
