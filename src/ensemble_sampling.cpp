#include "ensemble.hpp"
#include "parser++.hpp"

using namespace imcmc::parser;

namespace imcmc{

    double ensemble_workspace::g( double z ){
        if( z >= 1.0/efficient_a && z <= efficient_a )
            return 1.0 / sqrt(z);
        else
            return 0.0;
    }

    double ensemble_workspace::gz(){

        double z;
        double zmin = 1./efficient_a;
        double zmax = efficient_a;
        double gmax = g(zmin);
        double gx;

        do{
            z     = gsl_ran_flat(rand_seed, zmin, zmax);
            gx    = gsl_ran_flat(rand_seed,  0.0, gmax);
        } while( gx > g(z) );

        return z;
    }

    void ensemble_workspace::do_sampling(){

        int burnin_loops, sampling_loops;

        burnin_loops    = burnin_step;
        sampling_loops  = sample_step;

        // make sure searching lndet_min and chisq_min will be done at least one time during burn-in
        _searched_lndet_min_chisq_min_ = false;

        if( rank == ROOT_RANK ){
            if( start_from_existing_chains != true ){
                std::cout << "\n#  =====================================================\n";
                std::cout << "#  ensemble_workspace::do_sampling(): start burning\n";
                std::cout << "#  total evaluations: " << burnin_loops*walker_num << "\n";
                std::cout << "#  =====================================================\n\n";
            }
            else{
                std::cout << "\n#  =====================================================\n";
                std::cout << "#  --> start from existing chains !!!\n";
                std::cout << "#  ensemble_workspace::do_sampling(): start re-burning\n";
                std::cout << "#  total evaluations: " << burnin_loops*walker_num << "\n";
                std::cout << "#  =====================================================\n\n";
            }
        }

        if( save_burned_ashes && (rank == ROOT_RANK) ){

            if( start_from_existing_chains == false )
                chain_name = chain_root + "_ashes.txt";
            else
                chain_name = chain_root + "_ashes2.txt";    // improvement is needed !!!

            out_stream.open(chain_name.c_str(), std::ofstream::out);

            if( !out_stream.good() )
                imcmc_runtime_error( "Filed to open file: " + chain_name );

            //  write parameter names into the first line of the chain file
            imcmc_vector_string_iterator it = output_param_name.begin();

            //  Write the first line
            if( use_cosmomc_format && write_chain_header ){
                out_stream << "# standard cosmomc format.  This is the burned ashes ...\n";
                out_stream << "#";
                out_stream << std::setw(_OUT_WIDTH_-1) << "weight" << std::setw(_OUT_WIDTH_) << "-2log(L)";
            }
            else{
                out_stream << "# using MultiNest, the first column is exp(-0.5*chisq)\n";
                out_stream << "#";
                out_stream << std::setw(_OUT_WIDTH_-1) << "probability" << std::setw(_OUT_WIDTH_) << "chisq";
            }

            //  Write the second line
            if( write_chain_header ){
                while( it != output_param_name.end() ){   //  updates only the sampling parameters
                    out_stream << std::setw(_OUT_WIDTH_) << *it << "";
                    ++it;
                }
            }

            out_stream << "\n";
        }

        MPI::COMM_WORLD.Barrier();
        
        if( start_from_existing_chains == true ){
            burnin_loops = 50;
        }

        for( int j=0; j<burnin_loops; ++j ){    //  Burn in

            update_walkers( false, j, burnin_loops );

            if( use_cosmomc_format ){

                if( save_burned_ashes ){
                    if( rank == ROOT_RANK )
                        write_walkers(out_stream);
                }
                else
                    update_walkers_io();
            }

            // MPI::COMM_WORLD.Barrier();
        }

        if( rank == ROOT_RANK ) //  DONT forget to close out_stream before start REAL sampling
            out_stream.close();

        _searched_lndet_min_chisq_min_ = true;  //  once searched during the burning, change its state to TRUE.

        for( int i=1; i<=chain_num; ++i ){

            if( rank == ROOT_RANK ){

                if( chain_num == 1)
                    chain_name = chain_root + ".txt";
                else
                    chain_name = chain_root + "_" + Read::IntToString(i) + ".txt";

                out_stream.open(chain_name.c_str(), std::ofstream::out);

                if( !out_stream.good() )
                    imcmc_runtime_error("Filed to open file: " + chain_name);

                //  write parameter names into the first line of the chain file
                imcmc_vector_string_iterator it = output_param_name.begin();

                //  Write the first line
                if( use_cosmomc_format && write_chain_header ){
                    out_stream << "# standard cosmomc format\n";
                    out_stream << "#";
                    out_stream << std::setw(_OUT_WIDTH_-1) << "weight" << std::setw(_OUT_WIDTH_) << "-2log(L)";
                }
                else{
                    out_stream << "# ensemble format, exactly the same with MultiNest\n";
                    out_stream << "#";
                    out_stream << std::setw(_OUT_WIDTH_-1) << "probability" << std::setw(_OUT_WIDTH_) << "chisq";
                }

                //  Write the second line
                if( write_chain_header ){
                    while( it != output_param_name.end() ){   //  updates only the sampling parameters
                        out_stream << std::setw(_OUT_WIDTH_) << *it << "";
                        ++it;
                    }
                }

                out_stream << "\n";
            }

            if( rank == ROOT_RANK ){
                std::cout << "\n#  ===================================================\n";
                std::cout << "#  ensemble_workspace::do_sampling(): start sampling\n";
                std::cout << "#  ===================================================\n\n";
            }

            MPI::COMM_WORLD.Barrier();

            for( int j=0; j<sampling_loops; ++j ){

                update_walkers( true, j, sampling_loops );

                if( rank == ROOT_RANK )
                    write_walkers(out_stream);

                // MPI::COMM_WORLD.Barrier();
            }

            if( rank == ROOT_RANK ){

                out_stream.close();

                std::cout << "\n"
                    << "#  ===============================================\n"
                    << "#  summary of " << std::setw(3) << i << " -th chain:\n"
                    << "#  " << total_accepts << " of " << total_accepts + total_rejects << " walkers been accepted ...\n"
                    << "#  " << total_rejects << " of " << total_accepts + total_rejects << " walkers been rejected ...\n"
                    << "#  " << "acceptance ratio = " << double(total_accepts)/(total_accepts + total_rejects) << "\n"
                    << "#  " << "total number of likelihood erros: " << total_errors << "\n"
                    << "#  ===============================================\n\n";
            }

            MPI::COMM_WORLD.Barrier();

        //  after one chain is finished, you can choose to skip some steps
            if( (skip_step > 0)  && (i < chain_num) ){

                if( rank == ROOT_RANK )
                    std::cout << "\n ****** skipping some chains, can be viewed as extra burn-in ******\n";

                for( int j=0; j<skip_step; ++j ){
                    update_walkers( false, j, skip_step );
                    // MPI::COMM_WORLD.Barrier();
                }
            }

            MPI::COMM_WORLD.Barrier();
        }
    }


    int ensemble_workspace::update_a_walker( imcmc_double& full_param_temp, int current_id, int rand_id ){

        double z = ensemble_workspace::gz();

        imcmc_vector_string_iterator it = sampling_param_name.begin();
        while( it != sampling_param_name.end() ){   //  updates only the sampling parameters
            full_param_temp[*it] = walker[*it][rand_id] + z * (walker[*it][current_id] - walker[*it][rand_id]);
            ++it;
        }

        int state = 1;
        double lndet, chisq;

        likelihood_state.this_like_is_ok = true; //  MUST frist set to true

        if( !prior(full_param_temp) )    //  run out of prior
            state = 0;
        else{

            double lnpost_current   = walker["LnPost"][current_id];
            double lnpost_new       = likelihood_eval( full_param_temp, lndet, chisq );

            if( likelihood_state.this_like_is_ok ){

                double dlnpost      = lnpost_new - lnpost_current;
                double alpha        = GSL_MIN( 1., pow(z, sampling_param_num-1)*exp(dlnpost) );
                double ran          = gsl_ran_flat(rand_seed, 0, 1);

                if( ran <= alpha ){

                    state = 1;

                    walker["LnPost"][current_id]    = lnpost_new;
                    walker["LnDet"][current_id]     = lndet;
                    walker["Chisq"][current_id]     = chisq;

                //  =========================================================
                //  new position has been accepted, so update the parameters
                //  (both MCMC pars and derived ones)
                //  =========================================================
                    imcmc_vector_string_iterator it;

                    it = sampling_param_name.begin();
                    while( it != sampling_param_name.end() ){
                        walker[*it][current_id] = full_param_temp[*it];
                        ++it;
                    }

                    it = derived_param_name.begin();
                    while( it != derived_param_name.end() ){
                        walker[*it][current_id] = full_param_temp[*it];
                        ++it;
                    }
                }
                else
                    state = 0;  //  nothing to do ...
            }
            else
                state = 0;
        }

//      DEBUG
/*
         if( likelihood_state.this_like_is_ok )
             std::cout << "@rank: " << rank << " --- This Like is OK \n";
         else
             std::cout << "@rank: " << rank << " --- This Like is BAD \n";
*/

        return state;
    }

    void ensemble_workspace::update_walkers( bool do_sampling, int ith, int num ){   //  do_sampling = true means start to sample

    //  make a copy of full_param, and full_param will be updated inside update_a_walker()
        imcmc_double    full_param_temp(full_param);

        for( int i=0; i<walker_num; ++i )
            accept[i] = 0;

        int accept_tot = 0;
        int error_tot = 0;

        if( parallel_mode == 0 ){   //  serial mode

            int id;

            for( int i=0; i<walker_num; ++i ){
                do{    //    randomly select a walker in the complementary set of walkers
                    id    = gsl_rng_get(rand_seed) % walker_num;
                } while( id == i );

                accept[i] = update_a_walker( full_param_temp, i, id );

                if( likelihood_state.this_like_is_ok )
                    error[i] = 0;
                else
                    error[i] = 1;
            }
        }
        else{   //  opempi-parallel mode

            int num_each_proc;
            int *id_min     = new int[rank_size];
            int *id_max     = new int[rank_size];
            int *id_width   = new int[rank_size];

            int *sendcounts = new int[rank_size];
            int *recvcounts = new int[rank_size];
            int *displace   = new int[rank_size];

            int *walker_id = new int[walker_num];
            int walker_num_half = walker_num / 2;   //  walker_num should be an even number.

            for( int i=0; i<walker_num_half; ++i ){
                walker_id[i]                    = gsl_rng_get(rand_seed) % walker_num_half + walker_num_half;
                walker_id[walker_num_half + i]  = gsl_rng_get(rand_seed) % walker_num_half;
            }

            imcmc_vector_string_iterator it;

            if( parallel_mode == 1 ){

                if( walker_num_half%rank_size == 0 ){ //  each proc will handle same number of walkers

                    num_each_proc   = walker_num_half / rank_size;

                    for( int i=0; i<rank_size; ++i ){
                        id_width[i]   = num_each_proc;
                        id_min[i]     = i * num_each_proc;
                        id_max[i]     = id_min[i] + ( num_each_proc - 1 );
                        sendcounts[i] = id_width[i];
                        recvcounts[i] = sendcounts[i];
                        displace[i]   = i * num_each_proc;
                    }
                }
                else{

                    id_width[0]    = walker_num_half - (walker_num_half/rank_size + 1)*(rank_size-1);

                    if( id_width[0] <= 0 )
                        imcmc_runtime_error("id_width[0] must be positive integer, please adjust the number of walkers / processros");

                    id_min[0]     = 0;
                    id_max[0]     = id_width[0] - 1;
                    sendcounts[0] = id_width[0];
                    recvcounts[0] = sendcounts[0];
                    displace[0]   = 0;

                    num_each_proc    = walker_num_half/rank_size + 1;

                    for( int i=1; i<rank_size; ++i ){
                        id_width[i]   = num_each_proc;
                        id_min[i]     = id_width[0] + (i-1) * num_each_proc;
                        id_max[i]     = id_min[i] + ( num_each_proc - 1 );

                        sendcounts[i] = id_width[i] ;
                        recvcounts[i] = sendcounts[i];
                        displace[i]   = (i-1) * num_each_proc + id_width[0];
                    }
                }

                //  start to update walkers
                //  update the first half
                for( int i=id_min[rank]; i<=id_max[rank]; ++i ){

                    accept[i] = update_a_walker( full_param_temp, i, walker_id[i] );

                    if( likelihood_state.this_like_is_ok )
                        error[i] = 0;
                    else
                        error[i] = 1;
                }

                MPI::COMM_WORLD.Barrier();

            //  ==============================
            //  collecting sampling parameters
            //  ==============================
                it = sampling_param_name.begin();

                while( it != sampling_param_name.end() ){

                    MPI::COMM_WORLD.Gatherv(    &walker[*it][id_min[rank]],
                                                sendcounts[rank], MPI::DOUBLE,
                                                walker[*it],
                                                recvcounts, displace, MPI::DOUBLE,
                                                ROOT_RANK );

                    MPI::COMM_WORLD.Bcast(  walker[*it],
                                            walker_num,
                                            MPI::DOUBLE,
                                            ROOT_RANK    );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

            //  =============================
            //  collecting derived parameters
            //  =============================
                it = derived_param_name.begin();

                while( it != derived_param_name.end() ){

                    MPI::COMM_WORLD.Gatherv(    &walker[*it][id_min[rank]],
                                                sendcounts[rank], MPI::DOUBLE,
                                                walker[*it],
                                                recvcounts, displace, MPI::DOUBLE,
                                                ROOT_RANK );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

            //  count accepts
                MPI::COMM_WORLD.Gatherv(    &accept[id_min[rank]],
                                            sendcounts[rank], MPI::INT,
                                            accept,
                                            recvcounts, displace, MPI::INT,
                                            ROOT_RANK );

            //  count likelihood errors
                MPI::COMM_WORLD.Gatherv(    &error[id_min[rank]],
                                            sendcounts[rank], MPI::INT,
                                            error,
                                            recvcounts, displace, MPI::INT,
                                            ROOT_RANK );

                MPI::COMM_WORLD.Gatherv(    &walker["LnPost"][id_min[rank]],
                                            sendcounts[rank], MPI::DOUBLE,
                                            walker["LnPost"],
                                            recvcounts, displace, MPI::DOUBLE,
                                            ROOT_RANK );

                MPI::COMM_WORLD.Gatherv(    &walker["LnDet"][id_min[rank]],
                                            sendcounts[rank], MPI::DOUBLE,
                                            walker["LnDet"],
                                            recvcounts, displace, MPI::DOUBLE,
                                            ROOT_RANK );

                MPI::COMM_WORLD.Gatherv(    &walker["Chisq"][id_min[rank]],
                                            sendcounts[rank], MPI::DOUBLE,
                                            walker["Chisq"],
                                            recvcounts, displace, MPI::DOUBLE,
                                            ROOT_RANK );

//                MPI::COMM_WORLD.Barrier();

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

            //  =========================
            //  update the second half
            //  =========================
                id_min[rank] += walker_num_half;
                id_max[rank] += walker_num_half;

                for( int i=0; i<rank_size; ++i )
                    displace[i] += walker_num_half;

                for( int i=id_min[rank]; i<=id_max[rank]; ++i ){

                    accept[i] = update_a_walker( full_param_temp, i, walker_id[i] );

                    if( likelihood_state.this_like_is_ok )
                        error[i] = 0;
                    else
                        error[i] = 1;
                }

                MPI::COMM_WORLD.Barrier();

            //  ==============================
            //  collecting sampling parameters
            //  ==============================
                it = sampling_param_name.begin();

                while( it != sampling_param_name.end() ){

                    MPI::COMM_WORLD.Gatherv(    &walker[*it][id_min[rank]],
                                                sendcounts[rank], MPI::DOUBLE,
                                                walker[*it],
                                                recvcounts, displace, MPI::DOUBLE,
                                                ROOT_RANK );

                    MPI::COMM_WORLD.Bcast(  walker[*it],
                                            walker_num,
                                            MPI::DOUBLE,
                                            ROOT_RANK    );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

            //  =============================
            //  collecting derived parameters
            //  =============================

                it = derived_param_name.begin();

                while( it != derived_param_name.end() ){

                    MPI::COMM_WORLD.Gatherv(    &walker[*it][id_min[rank]],
                                                sendcounts[rank], MPI::DOUBLE,
                                                walker[*it],
                                                recvcounts, displace, MPI::DOUBLE,
                                                ROOT_RANK );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

            //  count accepts
                MPI::COMM_WORLD.Gatherv(    &accept[id_min[rank]],
                                            sendcounts[rank], MPI::INT,
                                            accept,
                                            recvcounts, displace, MPI::INT,
                                            ROOT_RANK );

            //  count likelihood errors
                MPI::COMM_WORLD.Gatherv(    &error[id_min[rank]],
                                            sendcounts[rank], MPI::INT,
                                            error,
                                            recvcounts, displace, MPI::INT,
                                            ROOT_RANK );

                MPI::COMM_WORLD.Gatherv(    &walker["LnPost"][id_min[rank]],
                                            sendcounts[rank], MPI::DOUBLE,
                                            walker["LnPost"],
                                            recvcounts, displace, MPI::DOUBLE,
                                            ROOT_RANK );

                MPI::COMM_WORLD.Gatherv(    &walker["LnDet"][id_min[rank]],
                                            sendcounts[rank], MPI::DOUBLE,
                                            walker["LnDet"],
                                            recvcounts, displace, MPI::DOUBLE,
                                            ROOT_RANK );

                MPI::COMM_WORLD.Gatherv(    &walker["Chisq"][id_min[rank]],
                                            sendcounts[rank], MPI::DOUBLE,
                                            walker["Chisq"],
                                            recvcounts, displace, MPI::DOUBLE,
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

                MPI::COMM_WORLD.Barrier();
            }
            else if( parallel_mode == 2 ){

            //  update first half
                accept[rank] = update_a_walker( full_param_temp, rank, walker_id[rank] );

                if( likelihood_state.this_like_is_ok )
                    error[rank] = 0;
                else
                    error[rank] = 1;

                MPI::COMM_WORLD.Barrier();

            //  ==============================
            //  collecting sampling parameters
            //  ==============================
                it = sampling_param_name.begin();

                while( it != sampling_param_name.end() ){

                    MPI::COMM_WORLD.Gather( &walker[*it][rank],
                                            1, MPI::DOUBLE,
                                            walker[*it],
                                            1, MPI::DOUBLE,
                                            ROOT_RANK );

                    MPI::COMM_WORLD.Bcast(  walker[*it],
                                            walker_num,
                                            MPI::DOUBLE,
                                            ROOT_RANK    );
                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

            //  ==============================
            //  collecting derived parameters
            //  ==============================
                it = derived_param_name.begin();

                while( it != derived_param_name.end() ){

                    MPI::COMM_WORLD.Gather( &walker[*it][rank],
                                            1, MPI::DOUBLE,
                                            walker[*it],
                                            1, MPI::DOUBLE,
                                            ROOT_RANK );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

                MPI::COMM_WORLD.Gather( &accept[rank],
                                        1, MPI::INT,
                                        accept,
                                        1, MPI::INT,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &error[rank],
                                        1, MPI::INT,
                                        error,
                                        1, MPI::INT,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &walker["LnPost"][rank],
                                        1, MPI::DOUBLE,
                                        walker["LnPost"],
                                        1, MPI::DOUBLE,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &walker["LnDet"][rank],
                                        1, MPI::DOUBLE,
                                        walker["LnDet"],
                                        1, MPI::DOUBLE,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &walker["Chisq"][rank],
                                        1, MPI::DOUBLE,
                                        walker["Chisq"],
                                        1, MPI::DOUBLE,
                                        ROOT_RANK );

                // MPI::COMM_WORLD.Barrier();

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

            //  update second half
                int rankx = rank + walker_num_half;
                accept[rankx] = update_a_walker( full_param_temp, rankx, walker_id[rankx] );

                if( likelihood_state.this_like_is_ok )
                    error[rankx] = 0;
                else
                    error[rankx] = 1;

                MPI::COMM_WORLD.Barrier();

            //  ==============================
            //  collecting sampling parameters
            //  ==============================
                it = sampling_param_name.begin();

                while( it != sampling_param_name.end() ){

                    MPI::COMM_WORLD.Gather( &walker[*it][rankx],
                                            1, MPI::DOUBLE,
                                            walker[*it],
                                            1, MPI::DOUBLE,
                                            ROOT_RANK );

                    MPI::COMM_WORLD.Bcast(  walker[*it],
                                            walker_num,
                                            MPI::DOUBLE,
                                            ROOT_RANK    );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

            //  ==============================
            //  collecting derived parameters
            //  ==============================
                it = derived_param_name.begin();

                while( it != derived_param_name.end() ){

                    MPI::COMM_WORLD.Gather( &walker[*it][rankx],
                                            1, MPI::DOUBLE,
                                            walker[*it],
                                            1, MPI::DOUBLE,
                                            ROOT_RANK );

                    ++it;
                    MPI::COMM_WORLD.Barrier();
                }

                MPI::COMM_WORLD.Gather( &accept[rankx],
                                        1, MPI::INT,
                                        accept,
                                        1, MPI::INT,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &error[rankx],
                                        1, MPI::INT,
                                        error,
                                        1, MPI::INT,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &walker["LnPost"][rankx],
                                        1, MPI::DOUBLE,
                                        walker["LnPost"],
                                        1, MPI::DOUBLE,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &walker["LnDet"][rankx],
                                        1, MPI::DOUBLE,
                                        walker["LnDet"],
                                        1, MPI::DOUBLE,
                                        ROOT_RANK );

                MPI::COMM_WORLD.Gather( &walker["Chisq"][rankx],
                                        1, MPI::DOUBLE,
                                        walker["Chisq"],
                                        1, MPI::DOUBLE,
                                        ROOT_RANK );

                // MPI::COMM_WORLD.Barrier();

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
            }
            else
                imcmc_runtime_error("Wrong parallel_mode !");

            delete[] walker_id;
        }

        //  =================================
        //  breifly summarize the samplings
        //  =================================

        if( rank == ROOT_RANK ){

            for( int i=0; i<walker_num; ++i ){

                accept_tot += accept[i];
                error_tot += error[i];

            //  continue to search for _lndet_min_ & _chisq_min_
                if( (do_sampling == false) && (_searched_lndet_min_chisq_min_==false) ){

                    if( walker["LnDet"][i] < _lndet_min_ )
                        _lndet_min_ = walker["LnDet"][i];

                    if( walker["Chisq"][i] < _chisq_min_ )
                        _chisq_min_ = walker["Chisq"][i];
                }
            }

            if( do_sampling == false ){
                std::cout   << "imcmc::ensemble " << std::setw(15) << "- burning - "
                            << "[" << std::setw(5) << ith+1 << " of " << std::setw(5) << num << "]\n"
                            << "# --> "
                            << std::setw(5) << accept_tot << " of " << std::setw(5) << walker_num
                            << " walkers updated ...\n"
                            << "# --> " << error_tot << " likelihood errors happened ..\n";
            }
            else{
                total_accepts += accept_tot;
                total_rejects += (walker_num-accept_tot);
                total_errors += error_tot;

                std::cout   << "imcmc::ensemble " << std::setw(15) << "- sampling - "
                            << "[" << std::setw(5) << ith+1 << " of " << std::setw(5) << num << "]\n"
                            << "# --> " 
                            << std::setw(5) << accept_tot << " of " << std::setw(5) << walker_num
                            << " walkers updated ...\n"
                            << "# --> " << error_tot << " likelihood errors happened ..\n";
            }
        }

        MPI::COMM_WORLD.Barrier();
    }

    void ensemble_workspace::update_walkers_io(){

        imcmc_vector_string_iterator it;

        if( rank == ROOT_RANK ){    //  only root rank need to do this

            for( int i=0; i<walker_num; ++i ){

                if( accept[i] == 0 ){   //  this old walker stays where it was, so we just increase its weight by 1.0
                    walker_io["Weight"][i] += 1.0;
                }
                else if( accept[i] == 1 ){  //  this old walker has been replaced by a new one, so we have to output it and update walker_io

                    //  update weight, lnpost, lndet and chisq
                    walker_io["Weight"][i] = 1.0;  //  reset to 1.0
                    walker_io["LnPost"][i] = walker["LnPost"][i];
                    walker_io["LnDet"][i]  = walker["LnDet"][i];
                    walker_io["Chisq"][i]  = walker["Chisq"][i];

                    it = output_param_name.begin();

                    while( it != output_param_name.end() ){ //  update to new walker
                        walker_io[*it][i] = walker[*it][i];
                        ++it;
                    }
                }
                else{
                    imcmc_runtime_error("unknown accept value, must be 0 or 1!");
                }
            }
        }
    }

}
