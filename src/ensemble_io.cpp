#include "ensemble.hpp"

namespace imcmc{

    void ensemble_workspace::write_walkers( std::ofstream& of ){

        imcmc_vector_string_iterator it;

        for( int i=0; i<walker_num; ++i ){

            it = output_param_name.begin();

            if( use_cosmomc_std_format ){    //    just set all the elements of the first column to 1

                of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << 1.0 << ""
                    << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << walker["Chisq"][i] << "";

            }
            else{

                of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << exp(walker["LnPost"][i] + _lndet_min_ + 0.5*_chisq_min_ ) << ""
                    << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << walker["Chisq"][i] << "";

            }

            while( it != output_param_name.end() ){
                of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << walker[*it][i] << "";
                ++it;
            }

            of << "\n";
        }
    }

    void ensemble_workspace::save_walker_state( std::ofstream& of ){
//      not implemented yet
    }

    void ensemble_workspace::read_walker_state( std::ofstream& of ){
//      not implemented yet
    }

}
