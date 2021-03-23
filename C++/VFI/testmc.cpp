///Testing mchain_sim
///08/17/20


#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <functional>

using namespace std;

double normal_cdf(double value){
  double num = 0.5 + 0.5*erf((value/pow(2,0.5)));
  return num;
}


void tauchen(double rho, double sigma, double mult, double trans_mat[], int entries, double z_grid[], int n_z){
    z_grid[0] = -mult*sigma/pow((1-pow(rho,2)),0.5);
    double inc = -2*z_grid[0]/(n_z-1);
    for(int i = 1; i<n_z;i++){
      z_grid[i] = z_grid[i-1] + inc;
    }
    for(int j = 0; j<n_z;j++){
      for(int k = 0; k<n_z;k++){
        if(k==0){
          double val1 = (z_grid[k]+0.5*inc-rho*z_grid[j])/sigma;
          trans_mat[j+n_z*k] = normal_cdf(val1);
        }else if(k==(n_z-1)){
          double val2 = (z_grid[k]-0.5*inc-rho*z_grid[j])/sigma;
          trans_mat[j+n_z*k] = 1-normal_cdf(val2);
        }else{
          double val1 = (z_grid[k]+0.5*inc-rho*z_grid[j])/sigma;
          double val2 = (z_grid[k]-0.5*inc-rho*z_grid[j])/sigma;
          trans_mat[j+n_z*k] = normal_cdf(val1)-normal_cdf(val2);
        }
      }
    }
  for(int l = 0;l<n_z;l++){
    z_grid[l] = exp(z_grid[l]);
  }
  return;
}


void rouwenhorst(double rho, double sigma, double mu, double trans_mat[],int entries, double z_grid[], int n_z){
    double qu = (rho+1)*0.5;
    double nu = pow((n_z-1)/(1-pow(rho,2)),0.5) * sigma;
    //cout<<qu<<endl;
    double inc = 2*nu/(n_z-1);
    z_grid[0] = mu/(1-rho)-nu;
    for(int i=1;i<n_z;i++){
        z_grid[i] = z_grid[i-1] + inc;
        z_grid[i-1] = exp(z_grid[i-1]);
    }
    z_grid[(n_z-1)] = exp(z_grid[(n_z-1)]);
    double* theta_nm1 = new double[n_z*n_z];
    theta_nm1[0] = qu;
    theta_nm1[n_z] = 1-qu;
    theta_nm1[1] = 1-qu;
    theta_nm1[1+n_z] = qu;
    //cout<<theta_nm1[0]<<" "<<theta_nm1[1]<<endl;
    //cout<<theta_nm1[2]<<" "<<theta_nm1[3]<<endl;
    if(n_z>2){
        for(int j=3;j<(n_z+1);j++){
            double* theta_n = new double[j*j];
            for(int r=0;r<j;r++){
                for(int c=0;c<j;c++){
                    if(r==0){
                        if(c==0){
                            theta_n[r+j*c] = qu*theta_nm1[r+n_z*c];
                        }else if(c==j-1){
                            theta_n[r+j*c] = (1-qu)*theta_nm1[r+n_z*(c-1)];
                        }else{
                            theta_n[r+j*c] = qu*theta_nm1[r+n_z*c] + (1-qu)*theta_nm1[r+n_z*(c-1)];
                        }
                    }else if(r==j-1){
                        if(c==0){
                            theta_n[r+j*c] = (1-qu)*theta_nm1[r-1+n_z*c];
                        }else if(c==j-1){
                            theta_n[r+j*c] = qu*theta_nm1[r-1-n_z+n_z*c];
                        }else{
                            theta_n[r+j*c] = (1-qu)*theta_nm1[r-1+n_z*c] + qu*theta_nm1[r-1-n_z+n_z*c];
                        }
                    }else{
                        if(c==0){
                            theta_n[r+j*c] = 0.5*(qu*theta_nm1[r+n_z*c]+(1-qu)*theta_nm1[r-1+n_z*c]);
                        }else if(c==j-1){
                            theta_n[r+j*c] = 0.5*(qu*theta_nm1[r-1+n_z*(c-1)]+(1-qu)*theta_nm1[r+n_z*(c-1)]);
                        }else{
                            theta_n[r+j*c] = 0.5*(qu*theta_nm1[r+n_z*c]+(1-qu)*theta_nm1[r-1+n_z*c]+qu*theta_nm1[r-1+n_z*(c-1)]+(1-qu)*theta_nm1[r+n_z*(c-1)]);
                        }
                    }
                }
            }
            //cout<<theta_n<<endl;
            for(int tr =0; tr<j;tr++){
              for(int tc = 0; tc<j;tc++){
                  theta_nm1[tr+n_z*tc] = theta_n[tr+j*tc];
              }
            }

            delete[] theta_n;

            for(int l = 0; l<n_z;l++){
                for(int m = 0; m<n_z;m++){
                    if(m==(n_z-1)){
                        cout<<theta_nm1[l+n_z*m]<<endl;
                    }else{
                        cout<<theta_nm1[l+n_z*m]<<" ";
                    }
                }
            }
            cout<<"-----------------------------------------------"<<endl;
            //cout<<theta_nm1[j*j]<<endl;
        }
    }
    for(int ent = 0;ent<entries;ent++){
        trans_mat[ent]=theta_nm1[ent];
    }
    delete[] theta_nm1;
  return;
}

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


void mchain(double trans_mat[], int entries, int i_state, double states[], int n_st, double results[], int length){
  for(int tr = 0; tr<length;tr++){
    double* probs = new double[n_st];
    if(tr == 0){
      for(int i = 0;i<n_st;i++){
        probs[i] = trans_mat[i_state+n_st*i];
      }
    }else{
      for(int i = 0;i<n_st;i++){
        probs[i] = trans_mat[(int)results[tr-1]+n_st*i];
      }
    }
    double cur_res[1];
    disc_samp(cur_res,1,probs,n_st);
    results[tr] = cur_res[0];
    results[tr-1] = states[(int)results[tr-1]];
    delete[] probs;
  }
  results[length-1] = states[(int)results[length-1]];
  return;
}


int main() {
    int num_z = 20;
    double t_mat_t[num_z*num_z];
    double z_vals_t[num_z];
    double t_mat_r[num_z*num_z];
    double z_vals_r[num_z];
    rouwenhorst(0.9,0.1,0,t_mat_r,num_z*num_z,z_vals_r,num_z);
    tauchen(0.9,0.1,3,t_mat_t,num_z*num_z,z_vals_t,num_z);
    int n_sim = 10000;
    double z_chain_t[n_sim];
    double z_chain_r[n_sim];
    mchain(t_mat_t,num_z*num_z,3,z_vals_t,num_z,z_chain_t,n_sim);
    mchain(t_mat_r,num_z*num_z,3,z_vals_r,num_z,z_chain_r,n_sim);
    string fileout = "Insert File";
    ofstream outchain;
    outchain.open(fileout);
    outchain<<"Tauchen,Rouwenhorst"<<endl;
    for(int i=0;i<n_sim;i++){
        outchain<<z_chain_t[i]<<","<<z_chain_r[i]<<endl;
    }
    cout<<"Writing done."<<endl;
    outchain.close();
    return 0;
}
