#ifndef __EMCEE__
#define __EMCEE__

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

#include "mpi.h"

extern "C"{
	#include <gsl/gsl_math.h>
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
}

#include "parser++.hpp"
#include "imcmc.hpp"

#define ROOT_RANK	0

namespace imcmc{

//	TODO: finish the walker_state structure
	struct walker_state{
//	...
	};

    class emcee_workspace{
        public:
	        emcee_workspace();
            ~emcee_workspace();

            std::string     config_file;

	        //	Parallel mode: 0, 1, 2
	        //	0: serial version
	        //	1: parallel version with small number of cores
	        //	2: parallel version with large number of cores
	        int				parallel_mode;
	        int				rank_size, rank;
	        int				burnin_step, skip_step, chain_num, chain_size;
	        int				walker_num;

	        //	name of current writing chain
	        std::string		chain_root;
	        std::string		chain_name;
	        std::ofstream	out_stream;

			//	whether use cosmomc standard form? "weight like param1 param2 param3 ..."
			//	NOTE: like is actually -2*ln(L) = chisq
			bool			use_cosmomc_std_format;	//	default: false

            //  filename of the parameters' limits
            std::string     param_limits;
            std::ofstream   param_limits_os;

	        //	random numbers
	        gsl_rng			*rand_seed;

	        //	efficient controling parameter
	        double			efficient_a;

	        //	state of acception or rejection, 1 or 0
	        int				total_accepts;
	        int				total_rejects;

	        //	**************************************
	        int				full_param_num;		//	number of full parameters
	        int				sampling_param_num;

	        imcmc_double	full_param;
	        imcmc_double	full_param_min;		//	max values
	        imcmc_double	full_param_max;		//	min values

            //  Walkers
	        imcmc_double_pointer	walker;     //  includes LnPost, LnDet and Chisq
            imcmc_vector_string     sampling_param_name;
            imcmc_vector_string     output_param_name;  //  if not set, then output_param_name = sampling_param_name

	        //	Likelihoof functions , contains MODEL and DATA
	        std::vector<likelihood_*>   likelihood;

	        void            set_efficient_a( double a );
	        inline	double  g( double z );
	        inline	double  gz();

	        void init( std::string infile );
	        void init_param();

            bool walker_initialized;
	        void init_walkers();
            void reset_walkers();

	        void add_likelihood( double (*like)( imcmc_double&, double&, double&, void*, void* ), 
			                    imcmc_vector_string 	modelparam,
			                    void 					*model,
			                    void 					*data );
        	//	check Samplingparams, if out of prior range, return false
	        bool            prior( imcmc_double& full_param );

            //  return log(posterior) = -lndet - 0.5*chisq
            //  det is the determinat in the denominator of the prefactor, plus some constant
	        double          likelihood_eval( imcmc_double& full_param, double& lndet, double& chisq );

            //  these two numbers will be used to re-scale the probability, in which case the chisq is too large
            //  so that exp(-lndet-0.5*chisq) --> 0
            double          _lndet_min_, _chisq_min_;
            bool            _searched_lndet_min_chisq_min_;
	
            int             update_a_walker( imcmc_double& full_param, int current_id, int rand_id );
	        void            update_walkers( bool do_sampling, int ith, int num ); //	Update walkers
	        void            write_walkers( std::ofstream& of ); //	Write walkers into text files
            void            do_sampling();

		//	TODO: add some new methods
		//	1) save walker state
			void save_walker_state( std::ofstream& of );	//	save current walkers' state into some file
			void read_walker_state( std::ofstream& of );	//	read the old walkers' state and then continue the MCMC sampling
    };

}

#endif // __PARAEMCEE__
