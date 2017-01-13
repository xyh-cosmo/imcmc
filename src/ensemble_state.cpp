#include "imcmc.hpp"
#include "ensemble.hpp"

using namespace imcmc::parser;

namespace imcmc{

    void ensemble_workspace::set_chkfile_width(int width){
        chkfile_width = width;
    }

    void ensemble_workspace::set_chkfile_precision(int precision){
        chkfile_precision = precision;
        if( precision >= chkfile_width ){
            chkfile_width = precision+2;
        }
    }

    void ensemble_workspace::save_state( int idx ){

        imcmc_vector_string_iterator it;
        std::string chkfile = chain_root + ".chk";
        std::string chkfile2 = chain_root + ".chk.temp";

        if( MPI::COMM_WORLD.Get_rank() == ROOT_RANK ){

            if( Read::Has_File(chkfile) ){ // save a backup
                std::cout << "==> making a backup of pre-existing check point file ...\n";
                backup_file(chkfile,chkfile2);
            }

            std::ofstream outfile(chkfile.c_str());

            std::cout << "==> saving new check point file ...\n\n";

            outfile << " walker_num = " << walker_num << std::endl;
            outfile << " chain_idx = " << idx << std::endl;

            outfile << "LnPost = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(chkfile_width)
                        << std::scientific
                        << std::setprecision(chkfile_precision)
                        << std::uppercase
                        << walker["LnPost"][i] << " ";
            }
            outfile << "\n";

            outfile << "LnDet = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(chkfile_width)
                        << std::scientific
                        << std::setprecision(chkfile_precision)
                        << std::uppercase
                        << walker["LnDet"][i] << " ";
            }
            outfile << "\n";

            outfile << "Chisq = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(chkfile_width)
                        << std::scientific
                        << std::setprecision(chkfile_precision)
                        << std::uppercase
                        << walker["Chisq"][i] << " ";
            }
            outfile << "\n";

            it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){
                outfile << *it << " = ";
                for( int i=0; i<walker_num; ++i ){
                    outfile << std::setw(chkfile_width)
                            << std::scientific
                            << std::setprecision(chkfile_precision)
                            << std::uppercase
                            << walker[*it][i] << " ";
                }
                ++it;
                outfile << "\n";
            }

            it = derived_param_name.begin();
            while( it != derived_param_name.end() ){
                outfile << *it << " = ";
                for( int i=0; i<walker_num; ++i ){
                    outfile << std::setw(chkfile_width)
                            << std::scientific
                            << std::setprecision(chkfile_precision)
                            << std::uppercase
                            << walker[*it][i] << " ";
                }
                ++it;
                outfile << "\n";
            }

            outfile << "\n\n";

        //  ====================================================================
        //  parameters for walker_io will be of the form par*, with an extra star

            outfile << "Weight* = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(5)
                        << std::setprecision(3)
                        << std::uppercase
                        << walker_io["Weight"][i] << " ";
            }
            outfile << "\n";

            outfile << "LnPost* = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(chkfile_width)
                        << std::scientific
                        << std::setprecision(chkfile_precision)
                        << std::uppercase
                        << walker_io["LnPost"][i] << " ";
            }
            outfile << "\n";

            outfile << "LnDet* = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(chkfile_width)
                        << std::scientific
                        << std::setprecision(chkfile_precision)
                        << std::uppercase
                        << walker_io["LnDet"][i] << " ";
            }
            outfile << "\n";

            outfile << "Chisq* = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(chkfile_width)
                        << std::scientific
                        << std::setprecision(chkfile_precision)
                        << std::uppercase
                        << walker_io["Chisq"][i] << " ";
            }
            outfile << "\n";

            it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){
                outfile << *it + "*"  << " = ";
                for( int i=0; i<walker_num; ++i ){
                    outfile << std::setw(chkfile_width)
                            << std::scientific
                            << std::setprecision(chkfile_precision)
                            << std::uppercase
                            << walker_io[*it][i] << " ";
                }
                ++it;
                outfile << "\n";
            }

            it = derived_param_name.begin();
            while( it != derived_param_name.end() ){
                outfile << *it + "*"  << " = ";
                for( int i=0; i<walker_num; ++i ){
                    outfile << std::setw(chkfile_width)
                            << std::scientific
                            << std::setprecision(chkfile_precision)
                            << std::uppercase
                            << walker_io[*it][i] << " ";
                }
                ++it;
                outfile << "\n";
            }

            outfile << "\n\n";

            outfile.close();
        }

//        MPI::COMM_WORLD.Barrier();
    }

    bool ensemble_workspace::read_state(){

        existed_chain_num = 0;
        bool read_success = true;

        std::string errmsg;
        imcmc_vector_string_iterator it;

        std::string chkfile = chain_root + ".chk";

        if( Read::Has_File(chkfile) == false ){
            std::cout << "==> failed to find check point file: " + chkfile;
            read_success = false;
            return read_success;
        }

        existed_chain_num = Read::Read_Int_from_File(chkfile,"chain_idx");

    //  only the root rank reads backup check point file.
        if( rank == ROOT_RANK ){

        //  make a simple check of the chkfile
            int walker_num_last_time = Read::Read_Int_from_File(chkfile,"walker_num");
            if( walker_num_last_time != walker_num ){
                std::cout << "==> *** Fatal Error in reading check point file ***\n"
                          << "==> Number of walkers used in the last MCMC sampling is: "
                          << walker_num_last_time << "\n"
                          << "==> But in this new run, you're using "
                          << walker_num << " walkers.\n"
                          << "==> please use the same number of walkers!\n";
                read_success = false;
            }

            if( !read_success ){
                std::string nw = Read::IntToString(walker_num_last_time);
                throw std::runtime_error("==> please reset a smaller walker number, no more than: "+nw);
            }

        //  now let's read the chkfile ....
            double *temp = new double[walker_num];

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "LnPost",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker["LnPost"][i] = temp[i];
            }

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "LnDet",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker["LnDet"][i] = temp[i];
            }

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "Chisq",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker["Chisq"][i] = temp[i];
            }

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "Weight*",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker_io["Weight"][i] = temp[i];
            }

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "LnPost*",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker_io["LnPost"][i] = temp[i];
            }

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "LnDet*",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker_io["LnDet"][i] = temp[i];
            }

            Read::Read_Array_of_Double_from_File(   chkfile,
                                                    "Chisq*",
                                                    temp,
                                                    walker_num );
            for( int i=0; i<walker_num; ++i ){
                walker_io["Chisq"][i] = temp[i];
            }

            it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){

                Read::Read_Array_of_Double_from_File(   chkfile,
                                                        *it,
                                                        temp,
                                                        walker_num );
                for( int i=0; i<walker_num; ++i ){
                    walker[*it][i] = temp[i];
                }

                Read::Read_Array_of_Double_from_File(   chkfile,
                                                        *it+"*",
                                                        temp,
                                                        walker_num );
                for( int i=0; i<walker_num; ++i ){
                    walker_io[*it][i] = temp[i];
                }

                ++it;
            }

            it = derived_param_name.begin();
            while( it != derived_param_name.end() ){

                Read::Read_Array_of_Double_from_File(   chkfile,
                                                        *it,
                                                        temp,
                                                        walker_num );
                for( int i=0; i<walker_num; ++i ){
                    walker[*it][i] = temp[i];
                }

                Read::Read_Array_of_Double_from_File(   chkfile,
                                                        *it+"*",
                                                        temp,
                                                        walker_num );
                for( int i=0; i<walker_num; ++i ){
                    walker_io[*it][i] = temp[i];
                }

                ++it;
            }

            delete[] temp;
        }

//        MPI::COMM_WORLD.Barrier();

        // std::cout << "==> Root rank finished re-loading check point file!\n";
        // exit(0);

    //  broadcast root rank's backup to all other ranks.

        it = sampling_param_name.end();
        while( it != sampling_param_name.end() ){

            if( rank == ROOT_RANK ){
                std::cout << "Bacsting: " << *it << std::endl;
            }

            MPI::COMM_WORLD.Bcast(  walker[*it],
                                    walker_num,
                                    MPI::DOUBLE,
                                    ROOT_RANK    );

            MPI::COMM_WORLD.Bcast(  walker_io[*it],
                                    walker_num,
                                    MPI::DOUBLE,
                                    ROOT_RANK    );
            ++it;
        }

        it = derived_param_name.begin();
        while( it != derived_param_name.end() ){

            if( rank == ROOT_RANK ){
                std::cout << "Bacsting: " << *it << std::endl;
            }

            MPI::COMM_WORLD.Bcast(  walker[*it],
                                    walker_num,
                                    MPI::DOUBLE,
                                    ROOT_RANK   );

            MPI::COMM_WORLD.Bcast(  walker_io[*it],
                                    walker_num,
                                    MPI::DOUBLE,
                                    ROOT_RANK   );

            ++it;
        }

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

        MPI::COMM_WORLD.Bcast(  walker_io["LnPost"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

        MPI::COMM_WORLD.Bcast(  walker_io["LnDet"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

        MPI::COMM_WORLD.Bcast(  walker_io["Chisq"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

        MPI::COMM_WORLD.Bcast(  walker_io["Weight"],
                                walker_num,
                                MPI::DOUBLE,
                                ROOT_RANK    );

// output the readed walkers ...
        // if( rank == ROOT_RANK ){
        //     imcmc_vector_string_iterator it = sampling_param_name.begin();
        //     while( it != sampling_param_name.end() ){
        //         std::cout << "param: " << *it << " = ";
        //         for( int i=0; i<walker_num; ++i ){
        //             std::cout << walker[*it][i] << " ";
        //         }
        //         std::cout << "\n";
        //         ++it;
        //     }
        // }
        //
        // exit(0);

//        MPI::COMM_WORLD.Barrier();

        return read_success;
    }

}
