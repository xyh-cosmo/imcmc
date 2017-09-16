#include "ensemble.hpp"

namespace imcmc{

    void ensemble_workspace::write_walkers( std::ofstream& of ){

        std::cout << "# --> saving walkers into chain files ...\n\n";

        imcmc_vector_string_iterator it;

        for( int i=0; i<walker_num; ++i ){

            if( accept[i] == 0 ){
            //  this old walker stays where it was, so we just increase its weight by 1.0
                walker_io["Weight"][i] += 1.0;
            }
            else if( accept[i] == 1 ){ //  this old walker has been replaced by a new one, so we have to output it and update walker_io

                of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << walker_io["Weight"][i] << ""
                    << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << walker_io["Chisq"][i] << "";

                it = output_param_name.begin();

                while( it != output_param_name.end() ){

                    of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                        << walker_io[*it][i] << "";
                    walker_io[*it][i] = walker[*it][i]; //  update to new walker
                    ++it;
                }

                of << "\n";

            //  update weight, lnpost, lndet and chisq
                walker_io["Weight"][i] = 1.0;  //  reset to 1.0
                walker_io["LnPost"][i] = walker["LnPost"][i];
                walker_io["LnDet"][i]  = walker["LnDet"][i];
                walker_io["Chisq"][i]  = walker["Chisq"][i];
            }
            else{
            //  This should never happen.
                // imcmc_runtime_error("unknown accept value, must be 0 or 1!");
                MPI_IMCMC_ERROR("unknown accept value, must be 0 or 1!");
            }
        }

    }

}
