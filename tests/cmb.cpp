#include <iostream>
#include <emcee++.hpp>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_integration.h>

using namespace std;
using namespace imcmc;

const double PI             = 3.14159265358979323846264338328;
const double C_Light_km_s   = 2.99782458e5;

double GSL_Integrator( double (*f)(double, void *), double x1, double x2, void *param ){
    double result, error;
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function  = f;
    F.params    = param;
    gsl_integration_qags( &F, x1, x2, 1e-10, 1e-10, 1000, w, &result, &error );
    gsl_integration_workspace_free(w);
    return result;
}

struct LCDM{
    imcmc_vector_string p_name;
    imcmc_double        p;

    double h;

    LCDM(){
//        p_name.push_back("H0");
//        p_name.push_back("Omegab");
//        p_name.push_back("Omegam");
//        p_name.push_back("Omegar");
//        p_name.push_back("Omegag");
//        p_name.push_back("Omegal");

//  add only parameters you want to sample ...
        p_name.push_back("H0");
        p_name.push_back("Omegab");
        p_name.push_back("Omegam");
        p_name.push_back("Neff");


        imcmc_vector_string_iterator it = p_name.begin();
        while( it != p_name.end() ){
            p[*it] = -9999;
            ++it;
        }

        p["H0"]     = 70;
        p["Omegab"] = 0.05;
        p["Omegam"] = 0.25;
        p["Neff"]   = 3.04;

        h           = p["H0"]/100.0;
		p["Omegag"]	= 2.469E-5/h/h;
        p["Omegar"]	= p["Omegag"]*(1+0.2271*p["Neff"]);
		p["Omegal"]	= 1.0 - p["Omegam"] - p["Omegar"];        
    }

    void Update( imcmc_double param ){
        imcmc_double_iterator it = p.begin();
        while( it != p.end() ){
            p[it->first] = param[it->first];
            ++it;
        }

        h           = p["H0"]/100.0;
		p["Omegag"]	= 2.469E-5/h/h;
        p["Omegar"]	= p["Omegag"]*(1+0.2271*p["Neff"]);
		p["Omegal"]	= 1.0 - p["Omegam"] - p["Omegar"];
    }

    //  inverse of the expansion rate, 1/E(z, param)
	static double iE(double z, void *param){
        LCDM *lcdm = static_cast<LCDM *>(param);
		double z2 = (1+z)*(1+z);
		double z3 = z2*(1+z);
		double z4 = z3*(1+z);

        double Omegam = lcdm->p["Omegam"];
        double Omegal = lcdm->p["Omegal"];
        double Omegar = lcdm->p["Omegar"];

		return  1 / sqrt(   Omegam*z3 
		                +   Omegal
		                +   Omegar*z4 );
	}

	double DC(double z){  //  Unit = Mpc
		double dc = GSL_Integrator(iE, 0, z, this);
        dc = dc/p["H0"];
        return dc * C_Light_km_s;
	}

	double DA(double z){  //  Unit = Mpc
		return this->DC(z) / (1+z);
	}

    double rs(){    //  fitting formula
        double zdec     = z_dec();
        double zeq      = p["Omegam"]/p["Omegar"];
        double R0       = 0.75*p["Omegab"]/p["Omegag"];
        double Rrec     = R0/(1+zdec);
        double Req      = R0/(1+zeq);

        return 2*C_Light_km_s*log((sqrt(1+Rrec)+sqrt(Rrec+Req))/(1+sqrt(Req)))/p["H0"]/sqrt(3*p["Omegam"]*zeq*Req);
    }

    double z_dec(){
        double omegab = p["Omegab"] * h * h;
        double omegam = p["Omegam"] * h * h;
        double g1 = 0.0783*pow(omegab, -0.238)/(1+39.5*pow(omegab, 0.763));
        double g2 = 0.560/(1+21.1*pow(omegab, 1.81));
        double zdec = 1048.*(1+0.00124*pow(omegab,-0.738))*(1+g1*pow(omegam,g2));
        return zdec;
    }

    double lA(){
        double zdec = z_dec();
        double dA = DA(zdec);
        return (1+zdec)*M_PI*dA/rs();
    }

	double R(){
        double zdec = z_dec();
		return p["H0"]*sqrt(p["Omegam"])*(1+zdec)*DA(zdec)/C_Light_km_s;
	}
};

double Like_CMB( imcmc_double& full_param, double& lndet, double& chisq, void* model, void* data){
    LCDM *lcdm = static_cast<LCDM *>(model);
    lcdm->Update(full_param);

	double lA_wmap7 = 302.09;
	double R_wmap7 = 1.725;
	double z_wmap7 = 1091.3;

	gsl_matrix *wmap7_cov_inv	= gsl_matrix_alloc(3,3);

	gsl_matrix_set( wmap7_cov_inv, 0, 0, 2.305 );
	gsl_matrix_set( wmap7_cov_inv, 0, 1, 29.698);
	gsl_matrix_set( wmap7_cov_inv, 0, 2, -1.333);

	gsl_matrix_set( wmap7_cov_inv, 1, 0, 29.698);
	gsl_matrix_set( wmap7_cov_inv, 1, 1, 6825.270);
	gsl_matrix_set( wmap7_cov_inv, 1, 2, -113.180);

	gsl_matrix_set( wmap7_cov_inv, 2, 0, -1.333);
	gsl_matrix_set( wmap7_cov_inv, 2, 1, -113.180);
	gsl_matrix_set( wmap7_cov_inv, 2, 2, 3.414);

	double wmap7_data[3] = {lA_wmap7, R_wmap7, z_wmap7};

    double zdec = lcdm->z_dec();
    double lA   = lcdm->lA();
    double R    = lcdm->R();

	double wmap7_model[3]  = { lA, R, zdec };

    lndet = 0;
    chisq = 0;

	for(int i=0; i<3; ++i){
		for( int j=0; j<3; ++j){
			chisq 	+= (wmap7_model[i] - wmap7_data[i]) 
					*  (wmap7_model[j] - wmap7_data[j]) 
					* gsl_matrix_get(wmap7_cov_inv, i, j);
		}
	}

	gsl_matrix_free(wmap7_cov_inv);
	return -0.5*chisq;
}


int main(int argc, char *argv[])
{
    MPI::Init(argc,argv);

    emcee_workspace ew;

    LCDM lcdm;

    ew.add_likelihood( Like_CMB, lcdm.p_name, &lcdm, NULL );
    ew.init("cmb.ini");
    ew.do_sampling();

    MPI::Finalize();        
}
