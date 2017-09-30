check_PROGRAMS = cmb_test gaussian_test rosenbrock_test himmelblau_test

TESTS=$(check_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir)/include

#  NOTE: for some versions of GNU compilers, errors like "undefined reference to ..." might happen,
# in that case, just un-comment " #gaussian_test_LDFLAGS = -Wl,--no-as-needed ", and run autoreconf to update configure file,
# after that run ./conifgure to generate new Makefiles.


# notice here that the compiler option: -Wl,--no-as-needed is for new versions of gnu compilers

gaussian_test_CPPFLAGS = -I$(top_srcdir)/include -g
gaussian_test_LDFLAGS = -Wl,--no-as-needed
gaussian_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
gaussian_test_SOURCES = gaussian.cpp

rosenbrock_test_CPPFLAGS = -I$(top_srcdir)/include -g
rosenbrock_test_LDFLAGS = -Wl,--no-as-needed
rosenbrock_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
rosenbrock_test_SOURCES = rosenbrock.cpp

himmelblau_test_CPPFLAGS = -I$(top_srcdir)/include -g
himmelblau_test_LDFLAGS = -Wl,--no-as-needed
himmelblau_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
himmelblau_test_SOURCES = himmelblau.cpp

cmb_test_CPPFLAGS = -I$(top_srcdir)/include -g
cmb_test_LDFLAGS = -Wl,--no-as-needed
cmb_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
cmb_test_SOURCES = cmb.cpp
