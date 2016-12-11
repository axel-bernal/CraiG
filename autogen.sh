#!/bin/sh
aclocal -I config && \
    case `uname` in
	Darwin*) glibtoolize --copy --force;;
	*) libtoolize --copy --force;;
    esac && \
    autoheader && \
    automake --gnu --add-missing --copy && \
    autoconf
