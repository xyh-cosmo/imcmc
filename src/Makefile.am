EXTRA_DIST =

lib_LTLIBRARIES = libimcmc.la
libimcmc_la_CXXFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/include
libimcmc_la_LDFLAGS = $(OPENMP_CFLAGS)
libimcmc_la_SOURCES =   ensemble.cpp               \
                        ensemble_init.cpp          \
                        ensemble_io.cpp            \
                        ensemble_prior.cpp         \
                        ensemble_likelihood.cpp    \
                        ensemble_sampling.cpp      \
                        parser++.cpp

#   $(includedir)/@PACKAGE@ is set to /usr/local/include/imcmc by default
imcmcincdir = $(includedir)/@PACKAGE@
imcmcinc_HEADERS =  $(top_srcdir)/include/ensemble.hpp       \
                    $(top_srcdir)/include/imcmc.hpp         \
                    $(top_srcdir)/include/parser++.hpp


#AM_CPPFLAGS = -I$(top_srcdir)/include

#   pepp is used to post-process chains sampled by Paraensemble, since the parameter orders might not
#   the same as we want, so we need to move the corresponding colums to their right places
#bin_PROGRAMS = pepp
#pepp_SOURCES = Paraensemble_PostProcess.cpp
#pepp_CXXFLAGS= -I$(top_srcdir)/include
