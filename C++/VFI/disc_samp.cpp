///Sampling from a discrete distribution C++
///08/17/20

#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>

using namespace std;

void disc_samp(double output[], int n_samp, double weights[], int n_weights){
  //sample from uniform distribution
  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<double> unif(0.0,1.0);
  auto unif_draw = bind(unif,mt);
  for(int ns = 0;ns<n_samp;ns++){
    output[ns] = unif_draw();
  }
  double* c_weights = new double[n_weights];
  //normalize weights
  double total;
  for(int i = 0; i<n_weights;i++){
    total+=weights[i];
  }
  for(int j = 0; j<n_weights; j++){
    if(j==0){
      c_weights[0] = weights[0]/total;
    }else{
      c_weights[j] = c_weights[j-1] + weights[j]/total;
    }
  }
  for(int k =0; k<n_samp;k++){
    for(int l = 0; l<n_weights;l++){
      if(l==0){
        if((output[k]>0)&(output[k]<c_weights[l])){
          output[k] = l;
        }
      }else{
        if((output[k]>c_weights[l-1])&(output[k]<c_weights[l])){
          output[k] = l;
        }
      }
    }
  }
  delete[] c_weights;
  return;
}

int main(){
  int n_s = 10;
  double rand_test[n_s];
  double samp_weights[4] = {0.1,0.2,0.3,0.4};
  disc_samp(rand_test,n_s,samp_weights,4);
  for(int run = 0; run<n_s;run++){
    cout<<rand_test[run]<<endl;
  }
  cout<<"--------------------------"<<endl;
  disc_samp(rand_test,n_s,samp_weights,4);
  for(int run = 0; run<n_s;run++){
    cout<<rand_test[run]<<endl;
  }
  return 0;
}
