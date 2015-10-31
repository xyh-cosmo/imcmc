#include "ensemble.hpp"
#include "parser++.hpp"

using namespace imcmc::parser;

namespace imcmc{

    ensemble_workspace::ensemble_workspace(){

		if( MPI::COMM_WORLD.Get_rank() == ROOT_RANK ){
			std::cout << "\n #  =================================================================\n";
			std::cout << " #  =============  Creating Ensemble Sampler WorkSpace  =============\n";
			std::cout << " #  =================================================================\n";
		}

        rank_size           = 0;
        walker_num          = 0;
        burnin_step         = 0;
        chain_num           = 0;
        sample_step         = 0;
        total_accepts       = 0;
        total_rejects       = 0;
        total_errors        = 0;
        sampling_param_num  = 0;
        full_param_num      = 0;
        chain_root          = "";
        chain_name          = "";
        config_file         = "";

        efficient_a         = 2.0;  // default value, will be updated from ini file
        init_ball_radius    = 0.5;  // default value, will be updated from ini file, but if not found in ini file, 0.5 is still OK.

        walker_initialized              = false;
        _searched_lndet_min_chisq_min_  = false;

        use_cosmomc_format  = true; //  default use CosmoMC format
        write_chain_header  = true; //  default true
        save_burned_ashes   = true;
    }

    ensemble_workspace::~ensemble_workspace(){

        if( walker_initialized ){

            gsl_rng_free(rand_seed);

            delete[] walker["LnPost"];
            delete[] walker["LnDet"];
            delete[] walker["Chisq"];

            delete[] accept;
            delete[] error;

            if( use_cosmomc_format ){
                delete[] walker_io["Weight"];
                delete[] walker_io["LnPost"];
                delete[] walker_io["LnDet"];
                delete[] walker_io["Chisq"];
            }

            imcmc_vector_string_iterator it;

            it = sampling_param_name.begin();

            while( it != sampling_param_name.end() ){
                delete[] walker[*it];

                if( use_cosmomc_format )
                    delete[] walker_io[*it];

                ++it;
            }

            it = derived_param_name.begin();

            while( it != derived_param_name.end() ){
                delete[] walker[*it];

                if( use_cosmomc_format )
                    delete[] walker_io[*it];

                ++it;
            }
        }

        if( likelihood.size() >= 1 ){

            std::vector<likelihood_>::size_type it_like;

            for( it_like = 0; it_like !=likelihood.size(); ++it_like ){
                delete likelihood[it_like];
            }
        }

        if( rank == ROOT_RANK ){

            std::cout << "#  ========================================================================\n";

            if( walker_initialized )
                std::cout << "#  ~ensemble_workspace(): sampling is over, clearing ensemble workspace ...\n";
            else
                std::cout << "#  ~ensemble_workspace(): no sampling, normal quit ...\n";

            std::cout << "#  ========================================================================\n";
        }
    }

}
