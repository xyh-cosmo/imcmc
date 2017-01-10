EXTRA_DIST =

lib_LTLIBRARIES = libimcmc.la
libimcmc_la_CXXFLAGS = $(OPENMP_CFLAGS) -I$(top_srcdir)/include
libimcmc_la_LDFLAGS = $(OPENMP_CFLAGS)
libimcmc_la_SOURCES =   imcmc.cpp                   \
                        ensemble.cpp                \
                        ensemble_init.cpp           \
                        ensemble_io.cpp             \
                        ensemble_prior.cpp          \
                        ensemble_likelihood.cpp     \
                        ensemble_sampling.cpp       \
                        ensemble_state.cpp          \
                        hamiltonian_config.cpp      \
                        hamiltonian.cpp             \
                        metropolis_config.cpp       \
                        metropolis.cpp              \
                        parser++.cpp

#   $(includedir)/@PACKAGE@ is set to /usr/local/include/imcmc by default
imcmcincdir = $(includedir)/@PACKAGE@
imcmcinc_HEADERS =  $(top_srcdir)/include/ensemble.hpp      \
                    $(top_srcdir)/include/imcmc.hpp         \
                    $(top_srcdir)/include/hamiltonian.hpp   \
                    $(top_srcdir)/include/metropolis.hpp    \
                    $(top_srcdir)/include/parser++.hpp
