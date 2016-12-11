#!/bin/sh
aclocal -I config \
&& glibtoolize --force --copy \
&& autoheader \
&& automake --gnu --add-missing --copy \
&& autoconf
