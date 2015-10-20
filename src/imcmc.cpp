#include "imcmc.hpp"
#include "mpi.h"

namespace imcmc{

    void imcmc_runtime_error( std::string err_info ){

        int rank = MPI::COMM_WORLD.Get_rank();

        std::cout << "# ===========  Error Message from rank: " << rank << "  ========== #\n";
        std::cout << "FileName: " << __FILE__ << "\n";
        std::cout << "Line Num: " << __LINE__ << "\n";
        std::cout << "Fun Name: " << __FUNCTION__ << "\n";
        std::cout << "####> ";

        throw std::runtime_error(err_info);
    }

    void imcmc_runtime_warning( std::string warn_info ){

        int rank = MPI::COMM_WORLD.Get_rank();

        std::cout << "# ==========  Warning Message from rank: << " << rank << "  ========== #\n";
        std::cout << "FileName: " << __FILE__ << "\n";
        std::cout << "Line Num: " << __LINE__ << "\n";
        std::cout << "Fun Name: " << __FUNCTION__ << "\n";
        std::cout << "----> ";
        std::cout << warn_info << std::endl;
    }

}
