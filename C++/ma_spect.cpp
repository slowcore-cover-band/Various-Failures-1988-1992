///Moving Average Spectra in C++
///Insert File
///08/14/2020

#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>

using namespace std;

void ma_spect(int order, int n_grid, string* file_loc){
  double domain[n_grid];
  double spectrum[n_grid];
  string filstring = *file_loc+"ma_spect"+to_string(order)+".csv";
  domain[0] = 0.0;
  for(int i=1;i<n_grid;i++){
    domain[i] = domain[i-1] + M_PI/(n_grid-1);
  }
  for(int j=0;j<n_grid;j++){
    double gain;
    complex<double> lp(0,0);
    complex<double> lpc(0,0);
    for(int k=1;k<(order+1);k++){
      complex<double> inexp(0,domain[j]*k);
      lp += exp(inexp);
      lpc += exp(-inexp);
      ///cout<<lp<<endl;
    }
    gain = sqrt(real(lp*lpc))/order;
    spectrum[j] = gain;
    gain = 0.0;
  }
  ofstream ma_file;
  ma_file.open(filstring);
  ma_file<<"domain,spectrum"<<endl;
  for(int l=0;l<n_grid;l++){
      ma_file<<domain[l]<<","<<spectrum[l]<<endl;
  }
  ma_file.close();
  return;
}

int main(){
  string file_loc = "Insert File";
  ma_spect(8,1000,&fl_loc);
  return 0;
}
