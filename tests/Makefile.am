check_PROGRAMS = cmb_test gaussian_test rosenbrock_test himmelblau_test

TESTS=$(check_PROGRAMS)

AM_CPPFLAGS = -I$(top_srcdir)/include

gaussian_test_CPPFLAGS = -I$(top_srcdir)/include
gaussian_test_LDFLAGS = -Wl,--no-as-needed
gaussian_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
gaussian_test_SOURCES = gaussian.cpp

rosenbrock_test_CPPFLAGS = -I$(top_srcdir)/include
rosenbrock_test_LDFLAGS = -Wl,--no-as-needed
rosenbrock_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
rosenbrock_test_SOURCES = rosenbrock.cpp

himmelblau_test_CPPFLAGS = -I$(top_srcdir)/include
himmelblau_test_LDFLAGS = -Wl,--no-as-needed
himmelblau_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
himmelblau_test_SOURCES = himmelblau.cpp

cmb_test_CPPFLAGS = -I$(top_srcdir)/include
cmb_test_LDFLAGS = -Wl,--no-as-needed
cmb_test_LDADD = ../src/libimcmc.la -lgsl -lgslcblas
cmb_test_SOURCES = cmb.cpp

