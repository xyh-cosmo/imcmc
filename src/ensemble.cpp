#include "ensemble++.hpp"
#include "parser++.hpp"

using namespace imcmc::parser;

namespace imcmc{

    ensemble_workspace::ensemble_workspace(){
        rank_size           = 0;
        walker_num          = 0;
        burnin_step         = 0;
        chain_num           = 0;
        chain_size          = 0;
        total_accepts       = 0;
        total_rejects       = 0;
        sampling_param_num  = 0;
        full_param_num      = 0;
        chain_root          = "";
        chain_name          = "";
        config_file         = "";

        walker_initialized              = false;
        _searched_lndet_min_chisq_min_  = false;
    }

    ensemble_workspace::~ensemble_workspace(){

        if( walker_initialized ){
            gsl_rng_free(rand_seed);

            delete[] walker["LnPost"];
            delete[] walker["LnDet"];
            delete[] walker["Chisq"];

            imcmc_vector_string_iterator it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){
                delete[] walker[*it];
                ++it;
            }
        }

        if( likelihood.size() >= 1 ){
            std::vector<likelihood_>::size_type it_like;
            for( it_like = 0; it_like !=likelihood.size(); ++it_like )
                delete likelihood[it_like];
        }

        if( walker_initialized )
            std::cout << "\n** ~ensemble_workspace(): sampling is over, clearing ensemble workspace ...\n\n";
        else
            std::cout << "\n** ~ensemble_workspace(): no sampling, normal quit ...\n\n";
    }

}
