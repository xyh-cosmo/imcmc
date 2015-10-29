#include "ensemble.hpp"
#include "parser++.hpp"

using namespace imcmc::parser;

namespace imcmc{

    void ensemble_workspace::add_likelihood( double (*loglike)( imcmc_double&, double&, double&, void*, void*, istate& ),
                                             imcmc_vector_string    modelparam,
                                             void*                  model,
                                             void*                  data ){

        likelihood_ *like   = new likelihood_;
        like->loglike       = loglike;
        like->model         = model;
        like->data          = data;

        likelihood.push_back( like );

        imcmc_vector_string_iterator it = modelparam.begin();

        while( it != modelparam.end() ){
            if( full_param.count(*it) == 0 ){    //    if not found
                full_param[*it]        = -9999;
                full_param_min[*it]    = -9999;
                full_param_max[*it]    = -9999;
            }
            ++it;
        }

        MPI::COMM_WORLD.Barrier();
    }

    void ensemble_workspace::add_likelihood( double (*loglike)( imcmc_double&, double&, double&, void*, void*, istate& ),
                                             imcmc_vector_string    modelparam,
                                             imcmc_vector_string    derivedparam,
                                             void*                  model,
                                             void*                  data ){

        likelihood_ *like   = new likelihood_;
        like->loglike       = loglike;
        like->model         = model;
        like->data          = data;

        likelihood.push_back( like );

        imcmc_vector_string_iterator it = modelparam.begin();

        while( it != modelparam.end() ){
            if( full_param.count(*it) == 0 ){    //    if not found
                full_param[*it]     = -9999;
                full_param_min[*it] = -9999;
                full_param_max[*it] = -9999;
            }
            ++it;
        }

        it = derivedparam.begin();

    //  explicitly adding derived parameters from the input argument derivedparam
        while( it != derivedparam.end() ){
            if( derived_param.count(*it) == 0 ){
                derived_param[*it]     = 0.0;  //  initialize to zeros
                derived_param_name.push_back(*it);
            }
            ++it;
        }

        MPI::COMM_WORLD.Barrier();
    }

    // //  return log(posterior) = -lndet - 0.5*chisq
    // double ensemble_workspace::likelihood_eval( imcmc_double& full_param,
    //                                             double& lndet,
    //                                             double& chisq   ){
    //
    //     double lndet_temp   = 0.0;
    //     double chisq_temp   = 0.0;
    //     double ln_post      = 0.0;
    //
    // //  DONOT ever forget to set these two to zero!!!!
    //     lndet = 0;
    //     chisq = 0;
    //
    //     std::vector<likelihood_>::size_type it_like;
    //
    //     for( it_like = 0; it_like !=likelihood.size(); ++it_like ){
    //         ln_post += likelihood[it_like]->loglike( full_param,
    //                                                  lndet_temp,     //    this will be first set to zero inside likelihood functions
    //                                                  chisq_temp,     //    this will be first set to zero inside likelihood functions
    //                                                  likelihood[it_like]->model,
    //                                                  likelihood[it_like]->data   );
    //
    //         lndet += lndet_temp;
    //         chisq += chisq_temp;
    //     }
    //
    //     return ln_post;
    // }

    double ensemble_workspace::likelihood_eval( imcmc_double& full_param,
                                                double& lndet,
                                                double& chisq   ){

        double lndet_temp   = 0.0;
        double chisq_temp   = 0.0;
        double ln_post      = 0.0;

    //  DONOT ever forget to set these two to zero!!!!
        lndet = 0;
        chisq = 0;

        std::vector<likelihood_>::size_type it_like;

        for( it_like = 0; it_like !=likelihood.size(); ++it_like ){
            
            ln_post += likelihood[it_like]->loglike( full_param,
                                                     lndet_temp,     //    this will be first set to zero inside likelihood functions
                                                     chisq_temp,     //    this will be first set to zero inside likelihood functions
                                                     likelihood[it_like]->model,
                                                     likelihood[it_like]->data,
                                                     likelihood_state );

        //  test whether there is any error happened
            if( likelihood_state.this_like_is_ok ){
                lndet += lndet_temp;
                chisq += chisq_temp;
            }
            else{
                if( likelihood_state.stop_on_error ){
                    likelihood_state.what_happened();
                    throw std::runtime_error("######## fatal error happened, see message printed above! ########");
                }
                else{// still going on the sampling, but will print the error information
                    likelihood_state.what_happened();
                    ln_post = _IMCMC_NEGINF_;
                    break;  //  jump out the loop
                }
            }
        }

        return ln_post;
    }

}
