//  Testing example: Rosenbrock

#include "ensemble.hpp"

using namespace std;
using namespace imcmc;

struct Rosenbrock{
    imcmc_vector_string p_name;
    imcmc_double        p;

    Rosenbrock(){
        p_name.push_back("x1");
        p_name.push_back("x2");

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
        double chisq = 0;
        chisq = (1-p["x1"])*(1-p["x1"]) + 100*pow(p["x2"]-p["x1"]*p["x1"], 2);
        chisq /= 20.;
        return chisq;
    }
};

double TestLike( imcmc_double&  full_param, 
                 double&        lndet, 
                 double&        chisq, 
                 void*          model, 
                 void*          data, 
                 istate&        state ){

    state.this_like_is_ok = true;
    state.store_mesg("nothing happened!");

    Rosenbrock *r = static_cast<Rosenbrock *>(model);
    r->Update(full_param);
    chisq = r->Chisq();
    lndet = 0;
    return -0.5*chisq;
}

int main( int argc, char **argv )
{
    MPI::Init(argc, argv);

    ensemble_workspace ew;

    Rosenbrock R;

    ew.add_likelihood( TestLike, R.p_name, &R, NULL );
    ew.init("rosenbrock.ini");
    ew.do_sampling();

    MPI::Finalize();
}
