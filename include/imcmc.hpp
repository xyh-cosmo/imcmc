#ifndef __IMCMC__
#define __IMCMC__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>

//  Version information
#ifndef __IMCMC_VERSION__
    #define __IMCMC_VERSION__   "0.0.5"
#endif

//  Author
#ifndef __IMCMC_AUTHOR__
    #define __IMCMC_AUTHOR__    "Youhua Xu"
#endif

//  E-mail
#ifndef __IMCMC_EMAIL__
    #define __IMCMC_EMAIL__     "yhxu@nao.cas.cn"
#endif

namespace imcmc{

    #define _OUT_WIDTH_ 18  //  the width of the outputs of the chains.
    #define ROOT_RANK    0

    typedef std::map<std::string, double>       imcmc_double;
    typedef std::map<std::string, double*>      imcmc_double_pointer;   //  used to construct walkers
    typedef std::map<std::string, std::string>  imcmc_string;

    typedef std::vector<double>         imcmc_vector_double;
    typedef std::vector<std::string>    imcmc_vector_string;

    typedef std::map<std::string, double>::iterator         imcmc_double_iterator;
    typedef std::map<std::string, std::string>::iterator    imcmc_string_iterator;

    typedef std::vector<double>::iterator       imcmc_vector_double_iterator;
    typedef std::vector<std::string>::iterator  imcmc_vector_string_iterator;

    struct likelihood_{
        void   *data;
        void   *model;

    //  use the name 'logpost' might be better
    //  double (*logpost)( imcmc_double&, double&, double&, void*, void* );
        double (*loglike)( imcmc_double&, double&, double&, void*, void* );

        ~likelihood_(){
            data    = NULL;
            model   = NULL;
            loglike = NULL;
        }
    };

    void imcmc_runtime_error( std::string err_info );
    void imcmc_runtime_warning( std::string warn_info );

}

#endif
