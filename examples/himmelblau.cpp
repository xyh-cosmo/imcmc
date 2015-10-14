//  Testing example: Himmelblau function

#include "ensemble.hpp"
using namespace std;
using namespace imcmc;

struct Himmelblau{
    imcmc_vector_string p_name;
    imcmc_double        p;

    Himmelblau(){
        p_name.push_back("x");
        p_name.push_back("y");

        imcmc_vector_string_iterator it = p_name.begin();
        while( it != p_name.end() ){
            p[*it] = -9999;
            ++it;
        }
    }

    void Update( imcmc_double full_param ){
        imcmc_double_iterator it = p.begin();
        while( it != p.end() ){
            p[it->first] = full_param[it->first];
            ++it;
        }
    }

    double Chisq(){
        return pow(p["x"]*p["x"]+p["y"]-11,2) + pow(p["x"]+p["y"]*p["y"]-7,2);
    }
};

double TestLike( imcmc_double& full_param, double& lndet, double& chisq, void* model, void* data ){
    Himmelblau *h = static_cast<Himmelblau *>(model);
    h->Update(full_param);
    chisq = h->Chisq();
    lndet = 0;
    return -0.5*chisq;
}

int main( int argc, char **argv )
{
    MPI::Init(argc, argv);

    ensemble_workspace ew;

    Himmelblau H;

    ew.add_likelihood( TestLike, H.p_name, &H, NULL );
    ew.init("himmelblau.ini");
    ew.do_sampling();

    MPI::Finalize();
}
