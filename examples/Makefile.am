check_PROGRAMS = cmb_test gaussian_test rosenbrock_test himmelblau_test

TESTS=$(check_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir)/include

#  NOTE: for some versions of GNU compilers, errors like "undefined reference to ..." might happen,
# in that case, just un-comment " #gaussian_test_LDFLAGS = -Wl,--no-as-needed ", and run autoreconf to update configure file,
# after that run ./conifgure to generate new Makefiles.


# notice here that the compiler option: -Wl,--no-as-needed is for new versions of gnu compilers

gaussian_test_CPPFLAGS = -I$(top_srcdir)/include -g
gaussian_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
gaussian_test_SOURCES = gaussian.cpp
if Linux
	gaussian_test_LDFLAGS = -Wl,--no-as-needed
else
	gaussian_test_LDFLAGS = 
endif

rosenbrock_test_CPPFLAGS = -I$(top_srcdir)/include -g
rosenbrock_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
rosenbrock_test_SOURCES = rosenbrock.cpp
if Linux
	rosenbrock_test_LDFLAGS = -Wl,--no-as-needed
else
	rosenbrock_test_LDFLAGS =
endif
#rosenbrock_test_LDFLAGS = -Wl,--no-as-needed

himmelblau_test_CPPFLAGS = -I$(top_srcdir)/include -g
himmelblau_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
himmelblau_test_SOURCES = himmelblau.cpp
if Linux
	himmelblau_test_LDFLAGS = -Wl,--no-as-needed
else
	himmelblau_test_LDFLAGS =
endif
#himmelblau_test_LDFLAGS = -Wl,--no-as-needed

cmb_test_CPPFLAGS = -I$(top_srcdir)/include -g
cmb_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
cmb_test_SOURCES = cmb.cpp
if Linux
	cmb_test_LDFLAGS = -Wl,--no-as-needed
else
	cmb_test_LDFLAGS = 
endif
#cmb_test_LDFLAGS = -Wl,--no-as-needed
