#include "imcmc.hpp"
#include "ensemble.hpp"
#include "parser++.hpp"

using namespace imcmc::parser;

namespace imcmc{

    void ensemble_workspace::set_efficient_a( double a ){

        if( a <= 0 )
            imcmc_runtime_error("\nensemble_workspace::set_efficient_a():\n a must be a positive number !");
        else if( a < 1.0 )
            efficient_a = 1.0 / a;
        else
            efficient_a = a;

        // MPI::COMM_WORLD.Barrier();
    }

    void ensemble_workspace::init( std::string paramfile ){

        rank        = MPI::COMM_WORLD.Get_rank();
        rank_size   = MPI::COMM_WORLD.Get_size();

        if( rank == ROOT_RANK ){
            std::cout   << "\n"
                        << "  #  =======================================================\n"
                        << "  #                    Welcome to imcmc \n"
                        << "  #\n"
                        << "  #                Version: " << __IMCMC_VERSION__ << "\n"
                        << "  #                Author : " << __IMCMC_AUTHOR__ << "\n"
                        << "  #                E-mail : " << __IMCMC_EMAIL__ << "\n"
                        << "  #  =======================================================\n";
        }

        likelihood_state.init();

        srand(time(NULL));

        config_file = paramfile;

        if( !Read::Has_File(paramfile) ){
            imcmc_runtime_error("Cannot open parameter file: " + paramfile + " !");
        }

        if( Read::Has_Key_in_File( paramfile, "walker_num" ) ){

            Read::Read_Value_from_File(paramfile, "walker_num", walker_num);

            if( (walker_num%2) != 0 ){
                imcmc_runtime_warning("We strongly suggest to use even number of walkers, while you set a odd number, so we increase it by 1.");
            }

            if( rank_size == 1 )
                parallel_mode   = 0;
            else if( rank_size < (walker_num/2) )
                parallel_mode   = 1;
            else if( rank_size >= (walker_num/2) ){
                parallel_mode   = 2;
                walker_num      = 2*rank_size;
                std::cout << "ensemble_workspace::Init() --> you have lots of cores, set the number of \
                             walkers to twice of number of ranks\n";
            }

            accept  = new int[walker_num];
            error   = new int[walker_num];

        }
        else
            imcmc_runtime_error("\'walker_num\' not found in:" + paramfile + " !");

        if( Read::Has_Key_in_File( paramfile, "burnin_step" ) ){
            Read::Read_Value_from_File(paramfile, "burnin_step", burnin_step);

            if( burnin_step <=0 ){
                burnin_step = 10;
                imcmc_runtime_warning("\'burnin_step\' not found in:" + paramfile + ", so we set it to default value 10.");
            }
        }
        else{
            burnin_step = 10;   //  this period is necessary, cannot be ignored.
            imcmc_runtime_warning("\'burnin_step\' not found in:" + paramfile + ", so we set it to default value 10.");
        }

        if( Read::Has_Key_in_File( paramfile, "skip_step" ) ){
            Read::Read_Value_from_File(paramfile, "skip_step", skip_step);
        }
        else{    //    set to default value 10
            skip_step = 10;
            imcmc_runtime_warning("\'skip_step\' not found in: " + paramfile + ", so we set it to default value 10.");
        }

        if( Read::Has_Key_in_File( paramfile, "chain_num" ) ){
            Read::Read_Value_from_File(paramfile, "chain_num", chain_num);

            if( chain_num <= 0 ){
                Info::ErrorInfo( "chain_num is <= 0, which means you do want any chains, so quit..." );
            }
        }

        if( Read::Has_Key_in_File( paramfile, "sample_step" ) ){
            Read::Read_Value_from_File(paramfile, "sample_step", sample_step);
        }
        else{
            imcmc_runtime_warning("\'sample_step\' not found, so the default value 50 will be used.");
            sample_step = 50;
        }

        if( Read::Has_Key_in_File( paramfile, "efficient_a" ) ){
            Read::Read_Value_from_File(paramfile, "efficient_a", efficient_a);
        }
        else{
            imcmc_runtime_warning("\'efficient_a\' not found, so the default value 2.0 will be used.");
        }

        if( Read::Has_Key_in_File( paramfile, "init_ball_radius" ) ){
            Read::Read_Value_from_File(paramfile, "init_ball_radius", init_ball_radius);

            if( init_ball_radius <=0.0 || init_ball_radius >= 0.999999 ){
                // imcmc_runtime_error("init_ball_radius should be greater than 0.0 and samller than 1.0");
                imcmc_runtime_warning("init_ball_radius should be greater than 0.0 and samller than 1.0, so I reset it to default value 0.5");
                init_ball_radius = 0.5;
            }
        }
        else{
            init_ball_radius    = 0.5;  //  set to default value 0.5
        }

        if( Read::Has_Key_in_File( paramfile, "start_from_fiducial" ) ){
            start_from_fiducial = Read::Read_Bool_from_File(paramfile, "start_from_fiducial");
        }

        if( Read::Has_Key_in_File( paramfile, "chain_root" ) ){
            Read::Read_Value_from_File(paramfile, "chain_root", chain_root);
            param_limits = chain_root + ".ranges";

            if( rank == ROOT_RANK ){
                std::string settings = chain_root + ".settings";
                used_settings.open( settings.c_str() );    //  will be closed at the end.
                std::cout << " ===> " << settings << " has been opened to save used settings\n";
            }
        }

        if( Read::Has_Key_in_File( paramfile, "use_cosmomc_format" ) ){
            use_cosmomc_format    = Read::Read_Bool_from_File(paramfile, "use_cosmomc_format");
        }

        if( Read::Has_Key_in_File( paramfile, "write_chain_header" ) ){
            write_chain_header = Read::Read_Bool_from_File(paramfile, "write_chain_header");
        }

        if( Read::Has_Key_in_File( paramfile, "save_burned_ashes" ) ){
            save_burned_ashes   = Read::Read_Bool_from_File(paramfile, "save_burned_ashes");

			if( rank == ROOT_RANK ){
				if( save_burned_ashes )
					std::cout << " ===> " << "burned chains (ashes) will be save to disk\n";
				else
					std::cout << " ===> " << "burned chains (ashes) will not be written into disk\n";
			}
        }

        if( Read::Has_Key_in_File( paramfile, "stop_on_error" ) ){
            likelihood_state.stop_on_error = Read::Read_Bool_from_File(paramfile, "stop_on_error");

            if( likelihood_state.stop_on_error )
                imcmc_runtime_warning("your sampling will stop when error(s) encountered inside likelihoods!");
        }

        if( Read::Has_Key_in_File( paramfile, "prompt_warning") ){
            likelihood_state.prompt_warning = Read::Read_Bool_from_File(paramfile, "prompt_warning");

            if( !likelihood_state.prompt_warning )
                imcmc_runtime_warning("you choose to ignore likelihood warning messages ...");
        }

    //  setup seeds for the random number generators
        unsigned long seed, rand_num;
        unsigned long *random_seeds = new unsigned long[rank_size];
        unsigned long *first_number = new unsigned long[rank_size];

//        rand_seed = gsl_rng_alloc(gsl_rng_taus2);
        rand_seed = gsl_rng_alloc(gsl_rng_mt19937);


        seed = rand();  //  inital seed

        for( int i=0; i<=rank; ++i ){ //  make sure differen ranks use differen random number seed
            gsl_rng_set( rand_seed, seed );
            seed = gsl_rng_get(rand_seed);
        }

        rand_num = gsl_rng_get(rand_seed);  // get the first random number for this rank

        MPI::COMM_WORLD.Gather( &seed, 1, MPI::UNSIGNED_LONG,
                                random_seeds, 1, MPI::UNSIGNED_LONG,
                                ROOT_RANK );

        MPI::COMM_WORLD.Gather( &rand_num, 1, MPI::UNSIGNED_LONG,
                                first_number, 1, MPI::UNSIGNED_LONG,
                                ROOT_RANK );

        if( rank == ROOT_RANK ){

            std::string     seedfile = chain_root + ".used_seeds";
            std::ofstream    outfile( seedfile.c_str() );
            outfile << std::setw(10) << "#    rank " << std::setw(12) << "rng_name "
                    << std::setw(15) << "ran_seeds " << std::setw(15) << "first_rand\n";

            for( int i = 0; i < rank_size; ++i ){
                outfile << std::setw(10) << i << " "
                        << std::setw(12) << gsl_rng_name(rand_seed) << " "
                        << std::setw(15) << random_seeds[i] << " "
                        << std::setw(15) << first_number[i] << "\n";
            }

            outfile.close();
        }

        delete[] random_seeds;
        delete[] first_number;

    //  save used settings:
        if( rank == ROOT_RANK ){

            int iterm_width = 40;

            used_settings << " # ========================= used settings ======================\n";
            used_settings << std::setw(iterm_width) << std::left << " chain_root" << " = " << chain_root << "\n";
            used_settings << std::setw(iterm_width) << std::left << " chain_num" << " = " << chain_num << "\n";
            used_settings << std::setw(iterm_width) << std::left << " walker_num" << " = " << walker_num << "\n";
            used_settings << std::setw(iterm_width) << std::left << " burnin_step" << " = " << burnin_step << "\n";
            used_settings << std::setw(iterm_width) << std::left << " sample_step" << " = " << sample_step << "\n";
            used_settings << std::setw(iterm_width) << std::left << " skip_step" << " = " << skip_step << "\n";
            used_settings << std::setw(iterm_width) << std::left << " efficient_a" << " = " << efficient_a << "\n";
            used_settings << std::setw(iterm_width) << std::left << " init_ball_radius" << " = " << init_ball_radius << "\n";

            if( start_from_fiducial )
                used_settings << std::setw(iterm_width) << std::left << " start_from_fiducial" << " = true\n";
            else
                used_settings << std::setw(iterm_width) << std::left << " start_from_fiducial" << " = false\n";

            if( use_cosmomc_format )
                used_settings << std::setw(iterm_width) << std::left << " used_cosmomc_format" << " = true\n";
            else
                used_settings << std::setw(iterm_width) << std::left << " used_cosmomc_format" << " = false\n";

            if( save_burned_ashes )
                used_settings << std::setw(iterm_width) << std::left << " save_burned_ashes" << " = true\n";
            else
                used_settings << std::setw(iterm_width) << std::left << " save_burned_ashes" << " = false\n";

            if( likelihood_state.stop_on_error )
                used_settings << std::setw(iterm_width) << std::left << " stop_on_error" << " = true\n";
            else
                used_settings << std::setw(iterm_width) << std::left << " stop_on_error" << " = false\n";

            if( likelihood_state.prompt_warning )
                used_settings << std::setw(iterm_width) << std::left << " prompt_warning" << " = true\n";
            else
                used_settings << std::setw(iterm_width) << std::left << " prompt_warning" << " = false\n";

            if( write_chain_header )
                used_settings << std::setw(iterm_width) << std::left << " write_chain_header" << " = true\n";
            else
                used_settings << std::setw(iterm_width) << std::left << " write_chain_header" << " = false\n";

            used_settings.close();
            std::cout << " ===> seetings have been saved\n\n";
        }

        // MPI::COMM_WORLD.Barrier();

        init_param();

        // MPI::COMM_WORLD.Barrier();

        init_walkers();

        walker_initialized = true;

        MPI::COMM_WORLD.Barrier();
    }


    void ensemble_workspace::init_param(){    //    loop over FullParams

        int params_width = 7;
        int derived_params_width = params_width;

        imcmc_double_iterator it = full_param.begin();

    //  obtain the max-width of parameters.
        while( it != full_param.end() ){
            if( it->first.size() > params_width )
                params_width = it->first.size();

            ++it;
        }

        if( rank == ROOT_RANK ){
            std::cout << "\n#  =============================================================\n"
                      << "#  ensemble_workspace::init_param():\n"
                      << "#  reading sampling parameters from: " + config_file << "\n";

            std::cout << std::setw(params_width) << std::left << " params:" << "  "
                      << std::right << std::setw(15) << "fid-value" << "  "
                      << std::right << std::setw(15) << "min-value" << "  "
                      << std::right << std::setw(15) << "max-value" << "\n";

            param_limits_os.open(param_limits.c_str());

            param_limits_os << "\n# This file contains parameters' limits, which might be useful when using\n"
                            << "# CosmoMC::getdist to process the sampled chains and make interesting plots.\n"
                            << "# Just copy the following into getdist parameter files (some *ini file)\n\n";

            param_limits_os << "# -------------------------------------------------------------------------\n"
                            << "# UPDATE (Nov-12-2015): GetDist now use new format of ranges: \n"
                            << "# 'parname min max',\n"
                            << "# so the old format limits[param] = min max is no longer used.\n"
                            << "# -------------------------------------------------------------------------\n\n";
        }

        // MPI::COMM_WORLD.Barrier();

        it = full_param.begin();

        while( it != full_param.end() ){

            ++full_param_num;

            int nvalue = Read::Num_of_Value_for_Key( config_file, it->first );

            if( nvalue > 0 ){

                double *par = new double[nvalue];

                Read::Read_Array_from_File(config_file, it->first, par, nvalue);

                if( nvalue == 1 ){   //  this parameter is fixed
                    full_param[it->first]       = par[0];
                    full_param_min[it->first]   = par[0];
                    full_param_max[it->first]   = par[0];
                }
                else if( nvalue == 3 ){
                    if( (par[1] < par[0]) && (par[0] < par[2]) ){       //  do MCMC
                        full_param[it->first]          = par[0];
                        full_param_min[it->first]      = par[1];
                        full_param_max[it->first]      = par[2];

                        sampling_param_name.push_back(it->first);
                        ++sampling_param_num;

                        if( rank == ROOT_RANK ){
                            // GetDist no loner use limits[param] = min max format

                            // param_limits_os << std::left << std::setw(params_width+8) << "limits[" + it->first + "]"
                            //                 << " = "
                            //                 << std::setw(10) << par[1] << " "
                            //                 << std::setw(10) << par[2] << "\n";

                            param_limits_os << std::left << std::setw(params_width+8) << it->first
                                            << std::setw(10) << par[1] << " "
                                            << std::setw(10) << par[2] << "\n";

                        }
                    }
                    else if( fabs(par[1]-par[2]) < 1E-15 ){ //  not do MCMC, be careful here
                        full_param[it->first]       = par[0];
                        full_param_min[it->first]   = par[1];
                        full_param_max[it->first]   = par[1];
                    }
                    else{
                        throw std::runtime_error("\nensemble_workspace::init_param( std::string paramfile\
                                             ) --> Check parameter value setting for " + it->first );
                    }
                }
                else{
                    throw std::runtime_error("\nensemble_workspace::init_param( std::string paramfile ) -\
                                            -> Check parameter value setting for " + it->first );
                }

                if( rank == ROOT_RANK ){
                    std::cout   << " " << std::setw(params_width) << std::left << it->first << " "
                                << std::right << std::setw(15) << full_param[it->first] << "  "
                                << std::right << std::setw(15) << full_param_min[it->first] << "  "
                                << std::right << std::setw(15) << full_param_max[it->first] << "\n";
                }

                delete[] par;
            }
            else{
                throw std::runtime_error( "\nensemble_workspace::InitParam( std::string paramfile ) --> \
                        input value(s) for " + it->first + " has not been found in " + config_file);
            }

            ++it;
        }

        if( rank == ROOT_RANK )
            param_limits_os.close();

    //  =========================================
    //  add derived parameters into full_param
    //  =========================================
        imcmc_double_iterator itd = derived_param.begin();

        derived_params_width = params_width;

        while( itd != derived_param.end() ){
            if( full_param.count(itd->first) == 0 ){    // make sure that the derived parameter is NOT in full_param

                full_param[itd->first]    = itd->second;

                if( itd->first.size() > derived_params_width )
                    derived_params_width = itd->first.size();
            }
            else
                imcmc_runtime_error("found duplicate parameter: " + itd->first + ", check your list of derived parameters!");

            ++itd;
        }

        if( Read::Has_Key_in_File( config_file, "output_params" ) ){

            int nvalue = Read::Num_of_Value_for_Key( config_file, "output_params" );

            if( nvalue == 0 ){

                output_param_name = sampling_param_name;

                if( rank == ROOT_RANK ){
                    std::cout << "\n#  =============================================================\n"
                              << "#  ensemble_workspace::init_param():\n"
                              << "#  keyword : output_params found in " + config_file + ", but no parameters\n"
                              << "#  were listed, so all sampling parameters will be output.\n";
                }
            }
            else{

                if( rank == ROOT_RANK ){
                    std::cout << "\n#  ==============================================\n"
                              << "#  ensemble_workspace::init_param():\n"
                              << "#  " << nvalue << " sampling parameters will be output:\n";
                }

                std::string *name = new std::string[nvalue];
                Read::Read_Array_from_File(config_file, "output_params", name, nvalue);

                int i=0;

                while( i < nvalue ){

                    //  check if there is any duplicate in name[]

                    for( int m=0; m<nvalue; ++m ){
                        for( int n=m+1; n<nvalue; ++n ){
                            if( name[m] == name[n] ){
                                std::string errmesg = "\nensemble_workspace::init_param()\n\tfound duplicate of sampling parameter: "
                                                    + name[m] + " in output_params, remove it.";
                                throw std::runtime_error(errmesg);
                            }
                        }
                    }

                    //  check if the readed param is in sampling_param_name:
                    int count = 0;
                    imcmc_vector_string_iterator it = sampling_param_name.begin();
                    while( it != sampling_param_name.end() ){
                        if( *it == name[i] )
                            ++count;
                        ++it;
                    }

                    if ( count == 1 ){
                        if( rank == ROOT_RANK )
                            std::cout << "*\tparam[" << std::setw(4) << i << "] : " << name[i] << "\t\t... sampling\n";

                        output_param_name.push_back( name[i] );
                        ++i;
                    }
                    else{
                        imcmc_runtime_error( name[i] + " is not in the sampling parameter list, check your input file.");
                    }
                }
            }
        }
        else{

            output_param_name = sampling_param_name;
            imcmc_vector_string_iterator it = derived_param_name.begin();

            if( rank == ROOT_RANK )
                std::cout << "\n # ensemble_workspace::init_param():\n"
                    << " # --> no keyword : output_params found in " + config_file + ", so all sampling\n"
                    << " # --> parameters will be written into chains.\n\n";
        }

        if( Read::Has_Key_in_File( config_file, "output_dparams" ) ){

            int nvalue = Read::Num_of_Value_for_Key( config_file, "output_dparams" );

            if( nvalue == 0 ){

            //  add all derived parameters to the end of output_param_name directly.
                imcmc_vector_string_iterator it = derived_param_name.begin();

                while( it != derived_param_name.end() ){
                    output_param_name.push_back(*it);
                    ++it;
                }

                if( rank == ROOT_RANK ){
                    std::cout << "\n#  =============================================================\n"
                              << "#  ensemble_workspace::init_param():\n"
                              << "#  keyword : output_dparams found in " + config_file + ", but no parameters\n"
                              << "#  were listed, so all derived parameters will be output.\n";
                }
            }
            else{

                if( rank == ROOT_RANK ){
                    std::cout << "\n#  ==============================================\n"
                              << "#  ensemble_workspace::init_param():\n"
                              << "#  " << derived_param.size() << " derived parameters will be output:\n";
                }

                std::string *name = new std::string[nvalue];
                Read::Read_Array_from_File(config_file, "output_dparams", name, nvalue);

                int i=0;

                while( i < nvalue ){

                    //  check if there is any duplicate in name[]

                    for( int m=0; m<nvalue; ++m ){
                        for( int n=m+1; n<nvalue; ++n ){
                            if( name[m] == name[n] ){
                                std::string errmesg = "\nensemble_workspace::init_param()\n\tfound duplicate of derived parameter: "
                                                    + name[m] + " in output_params, remove it.";
                                throw std::runtime_error(errmesg);
                            }
                        }
                    }

                    //  check if the readed derived parameter is in derived_param_name:
                    int count = 0;
                    imcmc_vector_string_iterator it = derived_param_name.begin();
                    while( it != derived_param_name.end() ){
                        if( *it == name[i] )
                            ++count;
                        ++it;
                    }

                    if ( count == 1 ){
                        if( rank == ROOT_RANK )
                            std::cout << "*\tparam[" << std::setw(4) << i << "] : " << name[i] << "\t\t... sampling\n";

                        output_param_name.push_back( name[i] );
                        ++i;
                    }
                    else{
                        imcmc_runtime_error( name[i] + " is not in the derived parameter list, check your input file.");
                    }
                }
            }
        }
        else{
            imcmc_vector_string_iterator it = derived_param_name.begin();

            while( it != derived_param_name.end() ){
                output_param_name.push_back(*it);
                ++it;
            }
        }

        if( rank == ROOT_RANK ){

            std::string ofile = chain_root + ".params";
            std::ofstream outfile( ofile.c_str() );

            if( !outfile.good() ){
                throw std::runtime_error("\nensemble_workspace::init_param(): failed to open " + ofile );
            }

            //  write the full parameters ...
            outfile << "\n# =========================================\n";
            outfile << "# === names and types of all parameters ===\n";
            outfile << "# =========================================\n";

        //    derived_params_width >= params_width !!
            outfile << " " << std::setw(derived_params_width) << std::left <<"params:" << " "
                    << std::setw(15) << std::right << "min-value" << " "
                    << std::setw(15) << std::right << "max-value" << " "
                    << std::setw(15) << std::right << "do sampling?" << "\n";

            it = full_param.begin();

            while( it != full_param.end() ){

                if( (full_param_min.count(it->first) == 1) && (full_param_max.count(it->first) == 1) ){

                    outfile << " " << std::setw(derived_params_width) << std::left << it->first << " "
                            << std::setw(15) << std::right << full_param_min[it->first] << " "
                            << std::setw(15) << std::right << full_param_max[it->first] << " ";

                    if( fabs(full_param_min[it->first] - full_param_max[it->first]) < 1.E-15 )
                        outfile << std::setw(15) << std::right << "no\n";
                    else
                        outfile << std::setw(15) << std::right << "yes\n";
                }
                else{

                    outfile << " " << std::setw(derived_params_width) << std::left << it->first << " "
                            << std::setw(15) << std::right << "--" << " "
                            << std::setw(15) << std::right << "--" << " ";

                    outfile << std::setw(15) << std::right << "derived\n";
                }

                ++it;
            }

            outfile << "\n# ==============================\n";
            outfile << "# == parameters being sampled ==\n";
            outfile << "# ==============================\n";

            imcmc_vector_string_iterator itx;

            itx = sampling_param_name.begin();

            while( itx != sampling_param_name.end() ){
                outfile << " " << std::setw(15) << std::left << *itx << "\n";
                ++itx;
            }

            outfile << "\n# ==============================\n";
            outfile << "# == parameters being derived ==\n";
            outfile << "# ==============================\n";

            itd = derived_param.begin();

            while( itd != derived_param.end() ){
                outfile << " " << std::setw(15) << std::left << itd->first << "\n";
                ++itd;
            }

            outfile << "\n# ==============================\n";
            outfile << "# ==== parameters to output ====\n";
            outfile << "# ==============================\n";

            itx = output_param_name.begin();

            while( itx != output_param_name.end() ){
                outfile << " " << std::setw(15) << std::left << *itx << "\n";
                ++itx;
            }

            outfile.close();
        }

        MPI::COMM_WORLD.Barrier();
    }


    void ensemble_workspace::init_walkers(){   //  NOTE: intialized walkers MUST lie in the valid prior!!!

        if( rank == ROOT_RANK ){

            std::cout << "\n#  ============================================\n"
                      << "#  imcmc::ensemble_workspace::init_walkers():\n"
                      << "#  initializing walkers ...\n"
                      << "#  searching _lndet_min_ & _chisq_min_ ...\n"
                      << "#  ============================================\n";

            //  _lndet_min_, _chisq_min_ are used only when writing probability into chains.

            std::cout << "\n#  ==============================    NOTE    ==============================\n"
                      << "#  this version of init_walkers() has been optimized to support parallel \n"
                      << "#  initialization, so that the time used to finish the initialization will \n"
                      << "#  be greatly reduced especially when the likelihoods need long time to \n"
                      << "#  compute.\n"
                      << "#  ========================================================================\n";
        }
        // else{
        //     std::cout << "rank: " << rank << " : waiting root rank ...." << std::endl;
        // }

        // MPI::COMM_WORLD.Barrier();

        imcmc_vector_string_iterator it         = sampling_param_name.begin();
        imcmc_vector_string_iterator it_derived = derived_param_name.begin();

        if( rank == ROOT_RANK ){
            std::cout << std::setw(60) << std::left << "# --> Creating Walkers for Ensemble_Workspace ... ";
        }

        while( it != sampling_param_name.end() ){

            walker[*it] = new double[walker_num];

            if( use_cosmomc_format )
                walker_io[*it] = new double[walker_num];

            ++it;
        }

        if( rank == ROOT_RANK ){
			std::cout << " [Done]\n";
            std::cout << std::setw(60) << std::left << "# --> Adding Derived Parameters ...";
        }

        // MPI::COMM_WORLD.Barrier();

        while( it_derived != derived_param_name.end() ){

            walker[*it_derived] = new double[walker_num];

            for( int i=0; i<walker_num; ++i )
                walker[*it_derived][i] = 0.0;

            if( use_cosmomc_format ){
                walker_io[*it_derived] = new double[walker_num];

                for( int i=0; i<walker_num; ++i )
                    walker_io[*it_derived][i] = 0.0;
            }

            ++it_derived;
        }

		if( rank == ROOT_RANK ){
			std::cout << " [Done]\n";
		}

        //  add LnPost, LnDet, and Chisq to walkers
        walker["LnPost"]    = new double[walker_num];
        walker["LnDet"]     = new double[walker_num];
        walker["Chisq"]     = new double[walker_num];

        if( use_cosmomc_format ){
            walker_io["Weight"] = new double[walker_num];
            walker_io["LnPost"] = new double[walker_num];
            walker_io["LnDet"]  = new double[walker_num];
            walker_io["Chisq"]  = new double[walker_num];
        }

        //  Now initialize
        imcmc_double    full_param_temp(full_param);

        int *sendcounts = new int[rank_size];
        int *recvcounts = new int[rank_size];
        int *displace   = new int[rank_size];

        if( rank == ROOT_RANK ){
            std::cout << std::setw(60) << std::left << "# --> start initializing walkers ...";
            std::cout << "\n";
        }

        //MPI::COMM_WORLD.Barrier();

        //  each rank will only calculate some of the likelihoods of the walkers.
        int i_start, i_end;

        if( (walker_num % rank_size) == 0 ){    //    each rank will evaluate the same number of likelihoods

            i_start = ( walker_num / rank_size ) * rank;
            i_end   = i_start + ( walker_num / rank_size ) - 1;

            for( int i=0; i<rank_size; ++i ){
                sendcounts[i]    = (i_end - i_start) + 1;
                recvcounts[i]    = sendcounts[i];
                displace[i]      = i * ( walker_num / rank_size );
            }
        }
        else{    //    rank_root will evaluate more likelihoods

            //  set for root rank
            i_start         = 0;
            i_end           = walker_num/rank_size + walker_num%rank_size - 1;

            sendcounts[0]   = (i_end - i_start) + 1;
            recvcounts[0]   = sendcounts[0];
            displace[0]     = 0;

            for( int i=1; i<rank_size; ++i ){   //  set for the remaining ranks

                i_start = walker_num/rank_size + walker_num%rank_size + (i-1)*(walker_num/rank_size);
                i_end   = i_start + (walker_num/rank_size - 1);

                sendcounts[i] = (i_end - i_start) + 1;
                recvcounts[i] = sendcounts[i];
                displace[i]   = sendcounts[0] + (i-1) * ( walker_num / rank_size );
            }
        }

        // MPI::COMM_WORLD.Barrier();

        for(int i=i_start; i<=i_end; ++i){

//            std::cout << " @ RANK: " << std::setw(4) << std::right << rank
//                      << " initializing " << std::setw(4) << std::right << i
//                      << " -th walker, [ i_start = " << std::setw(4) << std::right << i_start
//                      << ", i_end = " << std::setw(4) << std::right << i_end << "]\n";


            double lndet, chisq;
            double start_value, value_width;

            //  initialize full_param randomly. Note that full_param includes sampling parameters, fixed parameters and derived parameters
            it = sampling_param_name.begin();

            while( it != sampling_param_name.end() ){

                if( start_from_fiducial )
                    start_value = full_param[*it];
                else
                    start_value = 0.5*(full_param_min[*it] + full_param_max[*it]);

            //  initial guess of value_width
                value_width  = full_param_max[*it] - full_param_min[*it];

                if( start_from_fiducial ){  // fid-values are not necessarily the middle values.
                    double width_left  = fabs(start_value - full_param_min[*it]);
                    double width_right = fabs(start_value - full_param_max[*it]);

                    if( value_width > width_left )
                        value_width = width_left;

                    if( value_width > width_right )
                        value_width = width_right;
                }

            //  ==============================================================================
            //  just give the walkers some reasonable values...
            //  0.25 can be replaced by other values whoes absolute values are less than 0.5
            //  ==============================================================================
                walker[*it][i]      = gsl_ran_flat( rand_seed,
                                                    start_value - 0.5*init_ball_radius*value_width,
                                                    start_value + 0.5*init_ball_radius*value_width );

                full_param_temp[*it] = walker[*it][i];
                ++it;
            }

            //  if error happens during initialization, LnPost will be _IMCMC_LNPOST_MIN_

            walker["LnPost"][i] = likelihood_eval( full_param_temp, lndet, chisq );
            walker["LnDet"][i]  = lndet;
            walker["Chisq"][i]  = chisq;

            if( likelihood_state.this_like_is_ok ){
                error[i] = 0;
            }
            else{
                error[i] = 1;
                std::cout << " # ++++> RANK: " << rank << "  error happened when initializing walker["
                          << i << "],  [ i_start = " << i_start << ", i_end = " << i_end << "]\n";
            }
        }

        if( rank == ROOT_RANK ){
			std::cout << std::setw(60) << std::left << "# --> Initializing walkers ...";
            std::cout << " [Done]\n";
        }

        MPI::COMM_WORLD.Barrier();

        //  ===================================
        //  collecting all sampling parameters
        //  ===================================

        it = sampling_param_name.begin();

        while( it != sampling_param_name.end() ){

            MPI::COMM_WORLD.Gatherv(    &walker[*it][i_start],
                                        sendcounts[rank], MPI::DOUBLE,
                                        walker[*it],
                                        recvcounts, displace, MPI::DOUBLE,
                                        ROOT_RANK );

            // MPI::COMM_WORLD.Barrier();

            MPI::COMM_WORLD.Bcast(  walker[*it],
                                    walker_num,
                                    MPI::DOUBLE,
                                    ROOT_RANK    );

            ++it;
        }

        //  ===================================
        //  collecting all derived parameters
        //  ===================================

        it = derived_param_name.begin();

        while( it != derived_param_name.end() ){

            MPI::COMM_WORLD.Gatherv(    &walker[*it][i_start],
                                        sendcounts[rank], MPI::DOUBLE,
                                        walker[*it],
                                        recvcounts, displace, MPI::DOUBLE,
                                        ROOT_RANK );

            // MPI::COMM_WORLD.Barrier();

            MPI::COMM_WORLD.Bcast(  walker[*it],
                                    walker_num,
                                    MPI::DOUBLE,
                                    ROOT_RANK    );

            ++it;
        }

        MPI::COMM_WORLD.Gatherv(    &walker["LnPost"][i_start],
                                    sendcounts[rank], MPI::DOUBLE,
                                    walker["LnPost"],
                                    recvcounts, displace, MPI::DOUBLE,
                                    ROOT_RANK );

        MPI::COMM_WORLD.Gatherv(    &walker["LnDet"][i_start],
                                    sendcounts[rank], MPI::DOUBLE,
                                    walker["LnDet"],
                                    recvcounts, displace, MPI::DOUBLE,
                                    ROOT_RANK );

        MPI::COMM_WORLD.Gatherv(    &walker["Chisq"][i_start],
                                    sendcounts[rank], MPI::DOUBLE,
                                    walker["Chisq"],
                                    recvcounts, displace, MPI::DOUBLE,
                                    ROOT_RANK );

        MPI::COMM_WORLD.Gatherv(    &error[i_start],
                                    sendcounts[rank], MPI::INT,
                                    error,
                                    recvcounts, displace, MPI::INT,
                                    ROOT_RANK );

    //  broadcast walkers to each proc
        MPI::COMM_WORLD.Bcast(  walker["LnPost"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

        MPI::COMM_WORLD.Bcast(  walker["LnDet"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

        MPI::COMM_WORLD.Bcast(  walker["Chisq"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

        // copy walker into walker_io.
        if( use_cosmomc_format && (rank == ROOT_RANK) ){

            for( int i=0; i<walker_num; ++i ){  // set all initial weights to 1.0
                walker_io["Weight"][i] = 1.0;   //  because these walkers were born for then first time
                walker_io["LnPost"][i] = walker["LnPost"][i];
                walker_io["LnDet"][i]  = walker["LnDet"][i];
                walker_io["Chisq"][i]  = walker["Chisq"][i];
            }

        //  =====================
        //  initialize walker_io
        //  =====================

            it = sampling_param_name.begin();

            while( it != sampling_param_name.end() ){

                for( int i=0; i<walker_num; ++i )
                    walker_io[*it][i] = walker[*it][i];

                ++it;
            }

            it = derived_param_name.begin();

            while( it != derived_param_name.end() ){
                for( int i=0; i<walker_num; ++i )
                    walker_io[*it][i] = walker[*it][i]; //  all these initial values are zero

                ++it;
            }
        }

    //    search the _lndet_min_ & _chisq_min_
        _lndet_min_ = 1.0E99;
        _chisq_min_ = 1.0E99;

        int error_tot = 0;

        for( int i=0; i<walker_num; ++i ){
            if( walker["LnDet"][i] < _lndet_min_ )
                _lndet_min_ = walker["LnDet"][i];

            if( walker["Chisq"][i] < _chisq_min_ )
                _chisq_min_ = walker["Chisq"][i];

            error_tot += error[i];
        }

        if( rank == ROOT_RANK ){

            std::cout << "\n#  =============================================================\n"
                      << ">  searched _lndet_min_: " << _lndet_min_ << "\n"
                      << ">  searched _chisq_min_: " << _chisq_min_ << "\n"
                      << "#  they will be updated again in the burn-in process\n\n";

            std::cout << "#  ***  Note that this searching will be done only once  ***\n"
                      << "#  =============================================================\n";
            //  _lndet_min_, _chisq_min_ are used only when writing probability into chains.
            //
            std::cout << "#  " << error_tot << " likelihood errors encountered ...\n";
        }

        // MPI::COMM_WORLD.Barrier();

        delete[] sendcounts;
        delete[] recvcounts;
        delete[] displace;

        MPI::COMM_WORLD.Barrier();
    }


    void ensemble_workspace::reset_walkers(){  //  do similar stuff as ini_walkers(), no need to search _lndet_min_ & _chisq_min_

        if( rank == ROOT_RANK ){
            std::cout   << "imcmc::ensemble_workspace::reset_walkers():\n"
                        << "\tresetting walkers ...\n\n";
        }

        imcmc_vector_string_iterator it = sampling_param_name.begin();

        imcmc_double    full_param_temp(full_param);

        for(int i=0; i<walker_num; ++i){

            double lndet=0, chisq=0;

            it = sampling_param_name.begin();

            while( it != sampling_param_name.end() ){

                double mean_value   = 0.5*(full_param_min[*it]+full_param_max[*it]);
                double value_width  = full_param_max[*it] - full_param_min[*it];

                walker[*it][i]      = gsl_ran_flat( rand_seed,
                                                    mean_value - 0.5*init_ball_radius*value_width,
                                                    mean_value + 0.5*init_ball_radius*value_width );

                full_param_temp[*it] = walker[*it][i];
                ++it;
            }

            walker["LnPost"][i] = likelihood_eval( full_param_temp, lndet, chisq );
            walker["LnDet"][i]  = lndet;
            walker["Chisq"][i]  = chisq;
        }

        MPI::COMM_WORLD.Barrier();
    }

}
