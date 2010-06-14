##
##  Makefile for all PLATFORMs
##
##  Note:  the $PLATFORM environment variable must be set to the
##         platform that you want to compile for, e.g.,
##         setenv PLATFORM platform
##

default	:
	@if [ $${PLATFORM:-notset} = notset ] ; then \
	 echo "make: **ERROR - PLATFORM environment variable not set!" ; \
	else \
	 echo "make -f Makefile.$$PLATFORM" ; \
	 make -f Makefile.$$PLATFORM ; \
	fi

clean	:
	@if [ $${PLATFORM:-notset} = notset ] ; then \
	 echo "make: **ERROR - PLATFORM environment variable not set!" ; \
	else \
	 echo "make -f Makefile.$$PLATFORM clean" ; \
	 make -f Makefile.$$PLATFORM clean ; \
	fi
