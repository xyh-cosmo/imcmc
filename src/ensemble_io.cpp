#include "ensemble.hpp"

namespace imcmc{

//    void ensemble_workspace::write_walkers( std::ofstream& of, bool last ){
    void ensemble_workspace::write_walkers( std::ofstream& of ){

        imcmc_vector_string_iterator it;

        if( use_cosmomc_format ){

            for( int i=0; i<walker_num; ++i ){

                if( accept[i] == 0 ){
                //  this old walker stays where it was, so we just increase its weight by 1.0
                    walker_io["Weight"][i] += 1.0;
                }
                else if( accept[i] == 1 ){
                //  this old walker has been replaced by a new one, so we have to output it and update walker_io

                    if( walker_io["Chisq"][i] < _IMCMC_CHISQ_MAX_ ){

                        of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                            << walker_io["Weight"][i] << ""
                            << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                            << walker_io["Chisq"][i] << "";

                        it = output_param_name.begin();

                        while( it != output_param_name.end() ){

                            of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                                << walker_io[*it][i] << "";

                        //  update to new walker
                            walker_io[*it][i] = walker[*it][i];
                            ++it;
                        }

                        of << "\n";
                    }
					else{
						of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
							<< 0.0 << "" //	set weight to zero.
							<< std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
							<< walker_io["Chisq"][i] << "";
					}


                //  update weight, lnpost, lndet and chisq
                    walker_io["Weight"][i] = 1.0;  //  reset to 1.0
                    walker_io["LnPost"][i] = walker["LnPost"][i];
                    walker_io["LnDet"][i]  = walker["LnDet"][i];
                    walker_io["Chisq"][i]  = walker["Chisq"][i];
                }
                else{
                    imcmc_runtime_error("unknown accept value, must be 0 or 1!");
                }
            }
        }
        else{   //  the following needs update

            for( int i=0; i<walker_num; ++i ){

                of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << exp(walker["LnPost"][i] + _lndet_min_ + 0.5*_chisq_min_ ) << ""
                    << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                    << walker["Chisq"][i] << "";

                it = output_param_name.begin();
                while( it != output_param_name.end() ){
                    of  << std::setw(_OUT_WIDTH_) << std::scientific << std::setprecision(10) << std::uppercase
                        << walker[*it][i] << "";
                    ++it;
                }

                of << "\n";
            }
        }

    }

    void ensemble_workspace::save_walker_state( std::ofstream& of ){
//      not implemented yet
    }

    void ensemble_workspace::read_walker_state( std::ofstream& of ){
//      not implemented yet
    }

}
