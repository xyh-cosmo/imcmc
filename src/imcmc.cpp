#include "imcmc.hpp"
#include "mpi.h"

namespace imcmc{

    imcmc_likelihood_state::imcmc_likelihood_state(){
        this_like_is_ok = true;
        stop_on_error   = false;
    }

    void imcmc_likelihood_state::store_mesg( std::string& why ){
        strcpy(errmesg, why.c_str());
    }

    void imcmc_likelihood_state::store_mesg( char *why ){
        strcpy(errmesg, why);
    }

    void imcmc_likelihood_state::what_happened(){
        int rank = MPI::COMM_WORLD.Get_rank();
        std::cout << "\n# ===========  Likelihood Error  ========== #\n"
                  << "# we captured the following error:\n"
                  << errmesg << "\n\n";
        std::cout << "# =========================================\n";
    }

    void imcmc_runtime_error( std::string err_info ){

        int rank = MPI::COMM_WORLD.Get_rank();

        if( rank == ROOT_RANK ){
            std::cout << "\n# ===========  Error Message  ========== #\n";
            std::cout << "---@ FileName: " << __FILE__ << "\n";
            std::cout << "---@ Line Num: " << __LINE__ << "\n";
            std::cout << "---@ Fun Name: " << __FUNCTION__ << "\n";
            std::cout << "+++> ";
        }

        throw std::runtime_error(err_info);
    }

    void imcmc_runtime_warning( std::string warn_info ){

        int rank = MPI::COMM_WORLD.Get_rank();

        if( rank == ROOT_RANK ){
            std::cout << "\n# ==========  Warning Message  ========= #\n";
            std::cout << "---@ FileName: " << __FILE__ << "\n";
            std::cout << "---@ Line Num: " << __LINE__ << "\n";
            std::cout << "---@ Fun Name: " << __FUNCTION__ << "\n";
            std::cout << "---> ";
            std::cout << warn_info << std::endl;
        }
    }

}
