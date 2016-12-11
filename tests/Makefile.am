## Make Autoconf tests.

# Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006
# Free Software Foundation, Inc.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

# We don't actually distribute the testsuite, since one only
# needs m4 to build it, m4 being required anyway to install Autoconf.

## ------------ ##
## Test suite.  ##
## ------------ ##

TESTSUITE = $(srcdir)/testsuite
TESTSUITE_AT = programs.at
AUTOM4TE = autom4te
EXTRA_DIST = $(TESTSUITE_AT) testsuite.at atlocal.in package.m4

# Running the uninstalled scripts.
DISTCLEANFILES = atconfig atlocal $(TESTSUITE)
MAINTAINERCLEANFILES = Makefile.in

## ------------ ##
## package.m4.  ##
## ------------ ##

$(srcdir)/package.m4: $(top_srcdir)/configure.ac
	{					\
	  echo '# Signature of the current package.'; \
	  echo 'm4_define([AT_PACKAGE_NAME],      [$(PACKAGE_NAME)])'; \
	  echo 'm4_define([AT_PACKAGE_TARNAME],   [$(PACKAGE_TARNAME)])'; \
	  echo 'm4_define([AT_PACKAGE_VERSION],   [$(PACKAGE_VERSION)])'; \
	  echo 'm4_define([AT_PACKAGE_STRING],    [$(PACKAGE_STRING)])'; \
	  echo 'm4_define([AT_PACKAGE_BUGREPORT], [$(PACKAGE_BUGREPORT)])'; \
	} >$(srcdir)/package.m4



# Run autom4te.

AUTOTEST = $(AUTOM4TE) --language=autotest
$(TESTSUITE): $(srcdir)/package.m4 $(TESTSUITE_AT)
	$(AUTOTEST) -I $(srcdir) $@.at -o $@.tmp
	mv $@.tmp $@

subdir = tests

atconfig: $(top_builddir)/config.status
	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

clean-local:
	test ! -f $(TESTSUITE) || $(SHELL) $(TESTSUITE) --clean
	rm -f *.tmp
	rm -f -r autom4te.cache

check-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE) AUTOTEST_PATH="$(CRAIG_HOME)/src/" $(TESTSUITEFLAGS)

# Run the test suite on the *installed* tree.
installcheck-local: atconfig atlocal $(TESTSUITE)
	$(SHELL) $(TESTSUITE) AUTOTEST_PATH="$(bindir)" $(TESTSUITEFLAGS)
