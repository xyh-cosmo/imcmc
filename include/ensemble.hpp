/*  ============================================================================
    -----
    NEWS:
    -----
1)  Oct-28-2015: add derived_param_names, and the derived parameters will also
    be strored in full_param, but these derived parameters, will be calculated
    only in user-provied likelihoods, so these functionality actually depends
    on the users.
    ============================================================================*/

#ifndef __IMCMC_ENSEMBLE__
#define __IMCMC_ENSEMBLE__

#include "parser++.hpp"
#include "imcmc.hpp"

namespace imcmc{

    class ensemble_workspace{
        public:
            ensemble_workspace();
            ~ensemble_workspace();

            std::string config_file;

            std::ofstream   used_settings; // to save used settings

        //  Parallel mode: 0, 1, 2
        //  0: serial version
        //  1: parallel version with small number of cores
        //  2: parallel version with large number of cores
            int    parallel_mode;
            int    rank_size, rank;
            int    burnin_step, skip_step, sample_step;
            int    walker_num, chain_num;

        //  name of current writing chain
            std::string     chain_root;
            std::string     chain_name;
            std::ofstream   out_stream;

        //  whether use cosmomc standard form? "weight like param1 param2 param3 ..."
            bool    use_cosmomc_format;  //  default: false. NOTE: like is actually -2*ln(L) = chisq
            bool    write_chain_header;  //  default: true
            bool    save_burned_ashes;   //  whether save the burn-in chains, default yes

        //  filename of the parameters' limits
            std::string     param_limits;
            std::ofstream   param_limits_os;

        //  random numbers
            gsl_rng *rand_seed;
            gsl_rng *rand_seed_walker_id;   // for preparing walker_id[]

        //  efficient controling parameter
            double  efficient_a;

        //  radius of the initialization "ball", default value: 0.5
            double  init_ball_radius;
            bool    start_from_fiducial;

        //  state of acception or rejection, 1 or 0
            int *accept;            // = new int[walker_num];
            int *error;
            int total_accepts;      //  total number of acceptance of each chain
            int total_rejects;      //  total number of rejections of each chain
            int total_errors;       //  record how many likelihood error happens

            int full_param_num;     //  number of full parameters
            int sampling_param_num; //  number of sampling parameters
            int derived_param_num;  //  number of derived parameters

            imcmc_double    full_param;
            imcmc_double    full_param_min; // max values
            imcmc_double    full_param_max; // min values

            imcmc_double    derived_param;  //  just for adding

        //  Walkers
            imcmc_double_pointer    walker;                 //  includes LnPost, LnDet and Chisq
            imcmc_double_pointer    walker_io;              //  this is acutally a backup of walker, and it will be used to write chains into files.
            imcmc_vector_string     sampling_param_name;    //  hold the names of parameters being sampled
            imcmc_vector_string     derived_param_name;     //  hold the names of derived parameters
            imcmc_vector_string     output_param_name;      //  if not set, then output_param_name = sampling_param_name

            imcmc_likelihood_state  likelihood_state;       //  save current likelihood state, i.e., possible error information

        //  Likelihoof functions , include both MODELs and DATA
            std::vector<likelihood_*>   likelihood;

            void set_efficient_a( double a );   //    control the range of Z, narrower Z range will increase the acceptance ratio.

        //  used to generate random variable 'Z'
            inline double  g( double z );
            inline double  gz();

            bool walker_initialized;            // inidicate whether the walkers are initialized.

            void init( std::string infile );    // read initialization & other settings from the input *ini file
            void init_param();                  // initialize relavant parameters, some might be set to default values.
            void init_walkers();                // initialize the walkers.
            void init_walkers_from_chains();    // initialize walkers from existing chains, can save a lot of time for re-burning.
            void reset_walkers();               // reset the walkers to a random state.

        //  =======================================================================================================================
        //  add likelihood functions.  If you have many likelihoods functions, which share some common parameters, it'd better to
        //  combine those likelihood functions into a BIGGER one, and adding some flags to control which likelihood functions will
        //  be used.  This is quite common when someone is trying to constraining model parameters using different combinations of
        //  data sets, which may share some common nuisance parameters.

            void add_likelihood( double (*like)( imcmc_double&, double&, double&, void*, void*, imcmc_likelihood_state& ),
                                 imcmc_vector_string modelparam,
                                 void                *model,
                                 void                *data );

            void add_likelihood( double (*like)( imcmc_double&, double&, double&, void*, void*, imcmc_likelihood_state& ),
                                 imcmc_vector_string modelparam,
                                 imcmc_vector_string derivedparam,
                                 void                *model,
                                 void                *data );

            bool    prior( imcmc_double& full_param );    //    check Samplingparams, if out of prior range, return false

        //  return log(posterior) = -lndet - 0.5*chisq
        //  det is the determinat in the denominator of the prefactor, plus some constant
            double  likelihood_eval( imcmc_double& full_param, double& lndet, double& chisq );

        //  these two numbers will be used to re-scale the probability, in case that the chisq might be too large
        //  so that exp(-lndet-0.5*chisq) --> 0
            double  _lndet_min_, _chisq_min_;
            bool    _searched_lndet_min_chisq_min_;

            int     update_a_walker( imcmc_double& full_param, int current_id, int rand_id );
            void    update_walkers( bool do_sampling, int ith, int num );   // Update walkers
            void    update_walkers_io();                                    // needed when use_cosmomc_format == true
            void    write_walkers( std::ofstream& of );                     // Write walkers into text files
            void    do_sampling();

        //  new added to start sampling from existing chains.
        //  very important: To enable re-starting sampling from existing chains, ALL
        //  sampling parameters MUST be written into chains_*.txt, otherwise the walkers
        //  will not be re-initialized correctly.
            bool    start_from_existing_chians;     // default: false
            int     num_of_existing_chains;         // number of existing chains. This is needed to assign correct post-fix for new chains.
    };

}

#endif // __IMCMC_ENSEMBLE__
