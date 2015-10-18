#include "imcmc.hpp"

namespace imcmc{

    void imcmc_runtime_error( std::string err_info ){

        std::cout << "# ===========  Error Message  ========== #\n";
        std::cout << "FileName: " << __FILE__ << "\n";
        std::cout << "Line Num: " << __LINE__ << "\n";
        std::cout << "Fun Name: " << __FUNCTION__ << "\n";
        std::cout << "####> ";

        throw std::runtime_error(err_info);
    }

    void imcmc_runtime_warning( std::string warn_info ){

        std::cout << "# ==========  Warning Message  ========== #\n";
        std::cout << "FileName: " << __FILE__ << "\n";
        std::cout << "Line Num: " << __LINE__ << "\n";
        std::cout << "Fun Name: " << __FUNCTION__ << "\n";
        std::cout << "----> ";
        std::cout << warn_info << std::endl;
    }

}
