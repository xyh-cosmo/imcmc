#include "ensemble.hpp"
#include "parser++.hpp"

using namespace imcmc::parser;

namespace imcmc{

    bool ensemble_workspace::prior( imcmc_double& full_param_temp ){

        bool status = true;

        imcmc_vector_string_iterator it = sampling_param_name.begin();

        while( it != sampling_param_name.end() ){
            if( (full_param_temp[*it] < full_param_min[*it]) || (full_param_temp[*it] > full_param_max[*it]) ){
                status = false;
                break;
            }
            ++it;
        }

        return status;
    }

}
