#include <Rcpp.h>
#include <math.h>
#include <cmath>
using namespace Rcpp;

double Tref = 15.;
double kb   = 8.62E-5;
double T0   = 273.15;

// [[Rcpp::export]]
NumericVector TK(SEXP x_){
  //DESCRIPTION:
  //
  //INPUT PARAMETERS:
  //x: 
  //Normalize temperature in celcius 
  NumericVector  x(clone(x_));
  NumericVector Y;

  Y= 1./((Tref+T0)*kb) - 1./((x + T0)*kb);
  return Y;
}

// [[Rcpp::export]]
//Reverse function of T.K
NumericVector TKrev(SEXP x_){
  NumericVector x(clone(x_));
  NumericVector Y;
  Y= 1./(1./(Tref+T0) - x*kb) - T0;
  return Y;
}

/*Average light within MLD */
// [[Rcpp::export]]
NumericVector par_avg(SEXP par0_, SEXP mld_, SEXP chl_){
      double kw  =0.04; 
      double kchl=0.025;
      NumericVector par0(clone(par0_));
      NumericVector mld(clone(mld_));
      NumericVector chl(clone(chl_));
      NumericVector kext;
      NumericVector par_bom;
      NumericVector par_avg;
      kext     = kw + kchl * chl;
      par_bom  = par0 * exp(-mld*kext);
      par_avg  = (par0 - par_bom)/log(par0/par_bom);
      return par_avg;
}

// [[Rcpp::export]]
// Calculate thermal traits based on Topt
NumericVector mu0(SEXP x_, SEXP Y_, SEXP E0_){
      NumericVector x(clone(x_));
      NumericVector Y(clone(Y_));
      NumericVector E0(clone(E0_));
      return Y*exp(TK(x) * E0);
}
//
//  int N = NPPm.size();
//  int F = L * M * NMo;
//  double npp;
//  NumericVector NPPI(F);
//  NumericVector dep(NZ);
//  NumericVector  hz(NZ);
//
//  //Integrated NPP to be calculated (unit: mgC m-2 d-1)
//  //for (int i = 0; i < N; i++) {
//   for (int j =0; j < M;  j++){
//     for (int i =0; i < L; i++){
//       for (int k =0; k < NZ; k++){
//         dep[k] = Depth[k*L*M + j*L + i] ;  //Depth profile
//          hz[k] =    Hz[k*L*M + j*L + i] ;  //Depth profile
//         for (int t = 0; t < NMo; t++){
//             //NPP depth profile
//             if (dep[k] > -260.){
//               npp = NPPm[t*L*M*NZ + k*L*M + j*L + i];  
//
//               //Integrate through the water column
//               NPPI[t*L*M + j*L + i]+= npp*hz[k]; //Unit: mgC m-2 d-1
//             }
//         }
//       }
//     }
//   }
//  return NPPI;
//} 
/*
NumericVector Johnson(SEXP tC_, SEXP mu0_, SEXP Ea_, SEXP dED_, SEXP Topt_){
     NumericVector tC(clone(tC_));
     NumericVector Ea(clone(Ea_));
     NumericVector Topt(clone(Topt_));
     NumericVector mu0(clone(mu0_));
     NumericVector dED(clone(dED_));
     NumericVector b;
     NumericVector Y;
     NumericVector h;
     NumericVector ED;
     NumericVector a;

     if (any(dED <= 0.0)) {         	// dED has to be positive
        throw std::range_error("dED has to be positive");
     }
      ED = dED + Ea;
       a = Ea/dED;
      b = exp(ED/kb*(1./(Topt+T0)-1./(tC+T0)));
      b = a*b;
      Y = Ea/kb*(1./(Tref+T0) - 1./(tC+T0));
      h = mu0*exp(Y)/(1.+b);
      return h;
}

*/
/*phytoplankton growth rate, Chl:C and N:C function
NumericVector r(SEXP t_, SEXP N_, SEXP par_, SEXP Topt_, SEXP KN_, SEXP aI0_, SEXP mu0p_)
{

    NumericVector Q0N = NumericVector::create(0.05);
    NumericVector Qmax;
    NumericVector Ec1 = NumericVector::create(-0.69);
    NumericVector Ec3 = NumericVector::create(0.17);
    NumericVector Ea  = NumericVector::create(0.81);
    NumericVector ED0 = NumericVector::create(2.46);
    NumericVector thetamin=NumericVector::create(0.02);
    NumericVector thetamax=NumericVector::create(0.63);
    NumericVector t(clone(t_));
    NumericVector N(clone(N_));
    NumericVector par(clone(par_));
    NumericVector KN(clone(KN_));
    NumericVector aI0(clone(aI0_));
    NumericVector Topt(clone(Topt_));
    NumericVector mu0p(clone(mu0p_));

    NumericVector um;  // Nutrient and light replete growth rate for one species under  one particular temperature
    NumericVector mu15;  // Maximal growth rate normalized to 15 degree (depending on Topt)
    NumericVector dED;  // (Eh - Ei) normalized to 15 degree (depending on Topt)

    NumericVector SI;  // Light limitation index
    NumericVector fN;  // Nutrient limitation index
    NumericVector gr;  // Realised growth rate
    NumericVector QN;  // N:C ratio (mol:mol)
    NumericVector theta;  // Chl:C ratio (gChl:molC)

    par = par/.4;  //Change unit to W m-2

    mu15= mu0(Topt, mu0p, Ec1); 
    dED = mu0(Topt, ED0,  Ec3);
    um  = Johnson(t,mu15, Ea, dED, Topt);
    SI  = 1. - exp(-aI0*par/um);
    fN  = N/(N+KN);
    gr  = um*SI*fN;  //Growth rate

    Qmax= 3.*Q0N;
    QN  = Q0N/(1.-(1.-Q0N/Qmax)*fN);
  theta = thetamin + gr/par/aI0 * (thetamax-thetamin);   //Unit: gChl/molC
  return List::create(_["r"]     = gr,
                      _["QN"]    = QN,
                      _["theta"] = theta);
} */
 
