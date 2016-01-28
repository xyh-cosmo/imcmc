#ifndef __IMCMC__
#define __IMCMC__

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <map>
#include <cmath>
#include <ctime>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdexcept>

//  use macro __IMCMC_MPI__ to control source files.
#ifndef __IMCMC_MPI__
#define __IMCMC_MPI__
    #include "mpi.h"
#endif

extern "C"{
    #include <gsl/gsl_math.h>
    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>
}

//  ============================================================================
//  ============================================================================

//  ================
//   Version Number
//  ================

//  Version information
#ifndef __IMCMC_VERSION__
    #define __IMCMC_VERSION__   "0.1.1"
#endif

//  ==============
//  Author's Name
//  ==============
#ifndef __IMCMC_AUTHOR__
    #define __IMCMC_AUTHOR__    "Youhua Xu"
#endif

//  ==========
//  My E-mail
//  ==========
#ifndef __IMCMC_EMAIL__
    #define __IMCMC_EMAIL__     "yhxu@nao.cas.cn"
#endif

//  ============================================================================
//  ============================================================================

#ifndef _OUT_WIDTH_
    #define _OUT_WIDTH_ 18  //  the width of the outputs of the chains.
#endif

#ifndef ROOT_RANK
    #define ROOT_RANK    0
#endif

//  Maximum length of error information passed to imcmc_likelihood_state.
#ifndef _IMCMC_mesg_lenth_
    #define _IMCMC_MESG_LENGTH_ 1024
#endif

#ifndef _IMCMC_LNPOST_MAX_
    #define _IMCMC_LNPOST_MIN_ -1.0E50
#endif

#ifndef _IMCMC_CHISQ_MAX_
	#define _IMCMC_CHISQ_MAX_ 1.0E50
#endif


namespace imcmc{

//  imcmc key-value(s) date type
    typedef std::map<std::string, double>       imcmc_double;
    typedef std::map<std::string, double*>      imcmc_double_pointer;   //  used to construct walkers
    typedef std::map<std::string, std::string>  imcmc_string;

    typedef std::vector<double>         imcmc_vector_double;
    typedef std::vector<std::string>    imcmc_vector_string;

    typedef std::map<std::string, double>::iterator         imcmc_double_iterator;
    typedef std::map<std::string, std::string>::iterator    imcmc_string_iterator;

    typedef std::vector<double>::iterator       imcmc_vector_double_iterator;
    typedef std::vector<std::string>::iterator  imcmc_vector_string_iterator;

    struct imcmc_likelihood_state{  //  this will be used by user, so there is no suffix '_'
        bool this_like_is_ok;
        bool stop_on_error;
        bool prompt_warning;
        bool prompt_error;
        char errmesg[_IMCMC_MESG_LENGTH_];

        void store_mesg( std::string& why );
        void store_mesg( char *why );
        void what_happened();

        imcmc_likelihood_state();
		void init();
    };

    typedef struct imcmc_likelihood_state istate;

    struct likelihood_{

        void   *data;
        void   *model;

    //  use the name 'logpost' might be better
    //  double (*logpost)( imcmc_double&, double&, double&, void*, void* );
        // double (*loglike)( imcmc_double&, double&, double&, void*, void* );
        double (*loglike)( imcmc_double&, double&, double&, void*, void*, imcmc_likelihood_state& );

        ~likelihood_(){
            data    = NULL;
            model   = NULL;
            loglike = NULL;
        }
    };

    void init_random_seed();

    void imcmc_runtime_error( std::string err_info );
    void imcmc_runtime_warning( std::string warn_info );

//	new added verbose information displayer (@ Jan-28-2016)
	void imcmc_verbose( int rank, std::string& verbose_info );
}

#endif
