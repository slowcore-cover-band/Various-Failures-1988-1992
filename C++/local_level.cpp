///Local Level Model in C++
///12/25/2020

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>

using namespace std;

//this program reads in one dimensional data and estimates a local level model on that data, then writes the draws from the MCMC to .csv files

void onedraw_states(double mu0, double sigm, double sig2eps, double sig2nu, double obs[], int nT, double mu[], int nT1, double sig[], int nT2){
  //setting inital values
  random_device rd;
  mt19937 mt(rd());
  normal_distribution<double> normal_dr(mu0,sigm);
  auto norm_draw = bind(normal_dr,mt);

  mu[0] = norm_draw();
  sig[0] = sig2nu+sig2eps;
  //compute filtered values
  for(int ii=0; ii<nT; ii++){
    mu[ii+1] = mu[ii] + (obs[ii]-mu[ii])*(sig[ii]+sig2nu)/(sig[ii]+sig2nu+sig2eps);
    sig[ii+1] = sig[ii] - (sig[ii]+sig2nu)/(sig[ii]+sig2nu+sig2eps);
  }
  //drawing states

  //draw time T
  normal_distribution<double> normal_dr1(mu[nT],pow(sig[nT],0.5));
  auto norm_draw1 = bind(normal_dr1,mt);

  mu[nT] = norm_draw1();

  //draw other states
  for(int jj=(nT-1);jj>0;jj--){
    double mumu = sig[jj]*mu[jj+1]/(sig[jj]+sig2nu)+sig2nu*mu[jj]/(sig[jj]+sig2nu);
    double sigP = sig[jj]*sig2nu/(sig[jj]+sig2nu);
    sigP = pow(sigP,0.5);
    normal_distribution<double> normal_dr2(mumu,sigP);
    auto norm_draw2 = bind(normal_dr2,mt);
    mu[jj] = norm_draw2();
  }
  return;
}

void onedraw_sigs(double a, double b, double c, double d, double* sig2eps, double* sig2nu, double obs[], int nT, double mu[], int nT1){
  double sumnu2 = 0;
  double sumeps2 = 0;
  for(int ii=0; ii<nT;ii++){
    sumnu2+=pow((mu[ii+1]-mu[ii]),2);
    sumeps2+=pow((mu[ii+1]-obs[ii]),2);
  }
  double shanu = a + nT/2;
  double shaeps = c + nT/2;
  double ratnu = 1/(b + 0.5*sumnu2);
  double rateps = 1/(d + 0.5*sumeps2);

  random_device rd;
  mt19937 mt(rd());


  //draw nu
  gamma_distribution<double> gamma_dr1(shanu,ratnu);
  auto gamma_draw1 = bind(gamma_dr1,mt);
  *sig2nu = 1/gamma_draw1();
  //draw eps
  gamma_distribution<double> gamma_dr2(shaeps,rateps);
  auto gamma_draw2 = bind(gamma_dr2,mt);
  *sig2eps = 1/gamma_draw2();
  return;
}

void make_draws_llm(double a, double b, double c, double d, double mu0, double sigm, double data[], int nT, int num_draws, int print_ev, string* out_loc){
  //set up files
  string sifilstring = *out_loc+"/sigdraws.csv";
  string sdfilstring = *out_loc+"/stadraws.csv";

  //allocate space
  double* sig_draws = new double[2*(num_draws+1)];
  double* sta_draws = new double[(nT+1)*num_draws];

  //initialize sigma values
  double cursignu = 5;
  double cursigeps = 5;
  sig_draws[0] = cursignu;
  sig_draws[1] = cursigeps;

  //make draws
  for(int ii = 0; ii<num_draws; ii++){
    double curmu[nT+1];
    double cursig[nT+1];

    onedraw_states(mu0,sigm,cursigeps,cursignu,data,nT,curmu,nT+1,cursig,nT+1);
    onedraw_sigs(a,b,c,d,&cursigeps,&cursignu,data,nT,curmu,nT+1);

    sig_draws[ii*2+2] = cursignu;
    sig_draws[ii*2+3] = cursigeps;

    for(int jj=0; jj<(nT+1); jj++){
      sta_draws[jj+ii*(nT+1)] = curmu[jj];
    }
    if(ii%print_ev==0){
      cout<<ii+1<<" draws done!"<<endl;
    }
  }

  //write draws to files
  ofstream sig_file;
  sig_file.open(sifilstring);
  sig_file << "sig_nu2,sig_eps2" <<endl;
  for(int nr =0; nr<(num_draws+1); nr++){
    sig_file << sig_draws[2*nr] << "," << sig_draws[2*nr+1] << endl;
  }
  sig_file.close();

  ofstream state_file;
  state_file.open(sdfilstring);
  for(int nc=0;nc<num_draws;nc++){
    for(int nr=0;nr<(nT+1);nr++){
      if(nr == nT){
        state_file << sta_draws[nc*(nT+1)+nr]<< endl;
      }else{
        state_file << sta_draws[nc*(nT+1)+nr]<< "," ;
      }
    }
  }
  state_file.close();
  return;
}

int main(){
  //read in data (fake data for example)
  double tempvals[10] = {50,55,60,40,45,50,60,55,50,70};
  int nobs = 10;
  //set sigma parameters
  double aval = 10;
  double bval = 10;
  double cval = 30;
  double dval = 10;

  //set mu0 parameters
  double mu0val = 50;
  double sigmval = 5;

  //set number of draws
  int do_draws = 10000;
  int printev = 1000;

  //set location to write to
  string file_loc = "Insert File";
  string* fl_point = &file_loc;

  make_draws_llm(aval,bval,cval,dval,mu0val,sigmval,tempvals, nobs, do_draws, printev, fl_point);

  return 0;
}
