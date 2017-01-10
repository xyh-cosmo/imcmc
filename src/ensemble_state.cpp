#include "imcmc.hpp"
#include "ensemble.hpp"

using namespace imcmc::parser;

namespace imcmc{

    ensemble_state::ensemble_state(){
        walker_num  = 0;
        sampling_param_num = 0;
        derived_param_num  = 0;
        use_cosmomc_format = true;
        walker.clear();
        walker_io.clear();
        sampling_param_name.clear();
        derived_param_name.clear();
    }

    ensemble_state::~ensemble_state(){

        delete[] walker["LnPost"];
        delete[] walker["LnDet"];
        delete[] walker["Chisq"];

        if( use_cosmomc_format ){
            delete[] walker_io["Weight"];
            delete[] walker_io["LnPost"];
            delete[] walker_io["LnDet"];
            delete[] walker_io["Chisq"];
        }

        imcmc_vector_string_iterator it;

        if( sampling_param_name.size() > 0 ){
            it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){
                delete[] walker[*it];
                if( use_cosmomc_format )
                    delete[] walker_io[*it];
                ++it;
            }
        }

        if( derived_param_name.size() > 0 ){
            it = derived_param_name.begin();
            while( it != derived_param_name.end() ){
                delete[] walker[*it];
                if( use_cosmomc_format )
                    delete[] walker_io[*it];
                ++it;
            }
        }
    }

    void ensemble_state::init( ensemble_workspace& ew ){

        chkfile_root= ew.chain_root;
        walker_num  = ew.walker_num;
        use_cosmomc_format = ew.use_cosmomc_format;
        
        setw();
        setp();

        imcmc_vector_string_iterator it;

    //  allocate memory for walkers and walker_io
        it = ew.sampling_param_name.begin();
        sampling_param_num = 0;
        while( it != ew.sampling_param_name.end() ){
            sampling_param_name.push_back(*it);
            walker[*it]     = new double[walker_num];

            if( use_cosmomc_format )
                walker_io[*it]  = new double[walker_num];

            ++it;
            ++sampling_param_num;
        }

        it = ew.derived_param_name.begin();
        while( it != ew.derived_param_name.end() ){
            derived_param_name.push_back(*it);
            walker[*it]     = new double[walker_num];

            if( use_cosmomc_format )
                walker_io[*it]  = new double[walker_num];

            ++it;
            ++derived_param_num;
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

    //  copy initialized walkers & walker_io to ensemble_state
        take_a_snapshot(ew);
    }
    
    void ensemble_state::setw(int width){
        this->width = width;
    }
    
    void ensemble_state::setp(int precision){
        this->precision = precision;

        if( this->width < precision )
            this->width = precision + 2;
    }

    void ensemble_state::take_a_snapshot( ensemble_workspace& ew ){

        if( ew.walker_initialized == false ){
            std::string errmsg = "";
            errmsg += "==> void ensemble_state::take_shapshot( ensemble_workspace& ew ):\n";
            errmsg += " ew.walker_initialized = false, stop!\n";
            throw std::runtime_error(errmsg);
        }

        imcmc_vector_string_iterator it;

        it = sampling_param_name.begin();
        while( it != sampling_param_name.end() ){

            for( int i=0; i<walker_num; ++i ){
                walker[*it][i]      = ew.walker[*it][i];
                walker_io[*it][i]   = ew.walker_io[*it][i];
            }

            ++it;
        }

        it = derived_param_name.begin();
        while( it != derived_param_name.end() ){

            for( int i=0; i<walker_num; ++i ){
                walker[*it][i]      = ew.walker[*it][i];
                walker_io[*it][i]   = ew.walker_io[*it][i];
            }

            ++it;
        }
        
        for( int i=0; i<walker_num; ++i ){

            walker["LnPost"][i] = ew.walker["LnPost"][i];
            walker["LnDet"][i] = ew.walker["LnDet"][i];
            walker["Chisq"][i] = ew.walker["Chisq"][i];

            if( use_cosmomc_format == true ){
                walker_io["Weight"][i] = ew.walker_io["Weight"][i];
                walker_io["LnPost"][i] = ew.walker_io["LnPost"][i];
                walker_io["LnDet"][i] = ew.walker_io["LnDet"][i];
                walker_io["Chisq"][i] = ew.walker_io["Chisq"][i];
            }
        }
        
        MPI::COMM_WORLD.Barrier();
    }

    bool ensemble_state::save_state( int idx ){

        imcmc_vector_string_iterator it;
        std::string chkfile = chkfile_root + ".chk";
        std::string chkfile2 = chkfile_root + ".chk.temp";

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
                outfile << std::setw(width)
                        << std::scientific
                        << std::setprecision(precision)
                        << std::uppercase
                        << walker["LnPost"][i] << " ";
            }
            outfile << "\n";

            outfile << "LnDet = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(width)
                        << std::scientific
                        << std::setprecision(precision)
                        << std::uppercase
                        << walker["LnDet"][i] << " ";
            }
            outfile << "\n";

            outfile << "Chisq = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(width)
                        << std::scientific
                        << std::setprecision(precision)
                        << std::uppercase
                        << walker["Chisq"][i] << " ";
            }
            outfile << "\n";

            it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){
                outfile << *it << " = ";
                for( int i=0; i<walker_num; ++i ){
                    outfile << std::setw(width)
                            << std::scientific
                            << std::setprecision(precision)
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
                    outfile << std::setw(width)
                            << std::scientific
                            << std::setprecision(precision)
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
                outfile << std::setw(width)
                        << std::scientific
                        << std::setprecision(precision)
                        << std::uppercase
                        << walker_io["LnPost"][i] << " ";
            }
            outfile << "\n";

            outfile << "LnDet* = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(width)
                        << std::scientific
                        << std::setprecision(precision)
                        << std::uppercase
                        << walker_io["LnDet"][i] << " ";
            }
            outfile << "\n";

            outfile << "Chisq* = ";
            for( int i=0; i<walker_num; ++i ){
                outfile << std::setw(width)
                        << std::scientific
                        << std::setprecision(precision)
                        << std::uppercase
                        << walker_io["Chisq"][i] << " ";
            }
            outfile << "\n";

            it = sampling_param_name.begin();
            while( it != sampling_param_name.end() ){
                outfile << *it + "*"  << " = ";
                for( int i=0; i<walker_num; ++i ){
                    outfile << std::setw(width)
                            << std::scientific
                            << std::setprecision(precision)
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
                    outfile << std::setw(width)
                            << std::scientific
                            << std::setprecision(precision)
                            << std::uppercase
                            << walker_io[*it][i] << " ";
                }
                ++it;
                outfile << "\n";
            }

            outfile << "\n\n";

            outfile.close();
        }
        
        MPI::COMM_WORLD.Barrier();

    }

    bool ensemble_state::read_state(){
    
        bool read_success = true;

        std::string errmsg;
        imcmc_vector_string_iterator it;

        std::string chkfile = chkfile_root + ".chk";

//  Update is needed here !
//  should only the root rank do the check !

        if( Read::Has_File(chkfile) == false ){
            errmsg = "";
            errmsg += "failed to detect check point file: " + chkfile;
            throw std::runtime_error(chkfile);
        }

        MPI::COMM_WORLD.Barrier();

        existed_chain_num = Read::Read_Int_from_File(chkfile,"chain_idx") + 1;

        //  only the root rank reads backup file.
        if( MPI::COMM_WORLD.Get_rank() == ROOT_RANK ){

        //  make a simple check of the chkfile
            int walker_num_last_time = Read::Read_Int_from_File(chkfile,"walker_num");
            if( walker_num_last_time < walker_num ){
                std::cout << "==> *** Fatal Error in reading check point file ***\n"
                          << "==> Number of walkers used in the last MCMC sampling is: "
                          << walker_num_last_time << "\n"
                          << "==> But in this new run, you're using "
                          << walker_num << " walkers.\n"
                          << "==> please reset a smaller walker number, no more than: "
                          << walker_num_last_time << "\n";
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

            delete[] temp;
        }

        MPI::COMM_WORLD.Barrier();

    //  broadcast root rank's backup to all other ranks.

        it = sampling_param_name.end();
        while( it != sampling_param_name.end() ){

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

        if( use_cosmomc_format == true ){

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
        }
    }

    void ensemble_state::reset_ensemble_workspace( ensemble_workspace& ew ){
        // each rank copies exactly the same ensemble_state !

        imcmc_vector_string_iterator it;

        it = ew.sampling_param_name.begin();
        while( it != ew.sampling_param_name.end() ){

            for( int i=0; i<walker_num; ++i ){
                ew.walker[*it][i] = walker[*it][i];

                if( use_cosmomc_format == true )
                    ew.walker_io[*it][i] = walker_io[*it][i];
            }

            ++it;
        }

        it = ew.derived_param_name.begin();
        while( it != ew.derived_param_name.end() ){

            for( int i=0; i<walker_num; ++i ){
                ew.walker[*it][i] = walker[*it][i];

                if( use_cosmomc_format == true )
                    ew.walker_io[*it][i] = walker_io[*it][i];
            }

            ++it;
        }

        for( int i=0; i<walker_num; ++i ){

            ew.walker["LnPost"][i] = walker["LnPost"][i];
            ew.walker["LnDet"][i] = walker["LnDet"][i];
            ew.walker["Chisq"][i] = walker["Chisq"][i];

            if( use_cosmomc_format == true ){
                ew.walker_io["Weight"][i] = walker_io["Weight"][i];
                ew.walker_io["LnPost"][i] = walker_io["LnPost"][i];
                ew.walker_io["LnDet"][i] = walker_io["LnDet"][i];
                ew.walker_io["Chisq"][i] = walker_io["Chisq"][i];
            }
        }

    }
}
