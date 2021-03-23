///Stochastic Value Function Iteration (Neoclassical Growth Model)
///08/15/20

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

//normal_cdf function
double normal_cdf(double value){
  double num = 0.5 + 0.5*erf((value/pow(2,0.5)));
  return num;
}

//tauchen grid method
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

//rouwenhorst method
void rouwenhorst(double rho, double sigma, double mu, double trans_mat[],int entries, double z_grid[], int n_z){
    double qu = (rho+1)*0.5;
    double nu = pow((n_z-1)/(1-pow(rho,2)),0.5) * sigma;
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
        }
    }
    for(int ent = 0;ent<entries;ent++){
        trans_mat[ent]=theta_nm1[ent];
    }
    delete[] theta_nm1;
  return;
}

//utility function
double crra(double cons, double sigma){
  double val;
  if(sigma == 1){
    val = log(cons);
  }else{
    val = 1/((1-sigma)*pow(cons,sigma-1));
  }
  return val;
}

//computes steady state capital
double ss_cap(double alpha, double beta, double delta){
  double val = pow((alpha/((1/beta)+(delta-1))),(1/(1-alpha)));
  return val;
}

//computes maximum capital
double max_cap(double alpha, double delta, double Amax){
  double val =  pow((Amax*alpha/delta),(1/(1-alpha)));
  return val;
}

//equivalent to maximum(V[j,:,z]) in Julia or max(V[j,,z]) in R
void v_max(double v_mat[],int entries,int cols, int depth, int select_k_row, int select_z_col, double* max_val, int* max_loc){
  int rows = entries/(cols*depth);
  //cout<<rows<<endl;
  int cur_max_loc = 0;
  double cur_max_val = v_mat[select_k_row+cols*depth*select_z_col];
    for(int i=0; i<cols; i++){
      if(v_mat[select_k_row+i*cols+cols*depth*select_z_col]>cur_max_val){
        cur_max_val = v_mat[select_k_row+i*cols+cols*depth*select_z_col];
        cur_max_loc = i;
      }
    }
    *max_val = cur_max_val;
    *max_loc = cur_max_loc;
  return;
}

//implements VFI
void st_vfi_ncgm(double alpha, double beta, double delta, double sigma,int n_its,int n_grid, double rho, double mu, double sigma_z, int n_z, int print_every, string* out_loc){
  //setting the directory where we send results
  string kfilstring = *out_loc+"st_vfi_ncgm_kap_"+to_string(n_grid)+"_"+to_string(n_z)+".csv";
  string vfilstring = *out_loc+"st_vfi_ncgm_val_"+to_string(n_grid)+"_"+to_string(n_z)+".csv";

  //setting up the productivity grid and transition matrix
  //can use rouwenhorst or tauchen
  double* z_vals = new double[n_z];
  double* t_mat = new double[n_z*n_z];
  //tauchen(rho, sigma_z, mu, t_mat, n_z*n_z, z_vals, n_z);
  rouwenhorst(rho, sigma_z, mu, t_mat, n_z*n_z, z_vals, n_z);

  //setting up the capital grid
  double maxk = max_cap(alpha, delta, z_vals[n_z-1]);
  double lowerk = 0.1*maxk;
  double* kgrid = new double[n_grid];
  kgrid[0] = lowerk;
  for(int i=1; i<n_grid; i++){
    kgrid[i] = kgrid[i-1] + maxk/(n_grid-1);
  }

  //allocating value and policy functions
  double* Vikold = new double[n_grid*n_z];
  double* Viknew = new double[n_grid*n_z];
  double* kjknew = new double[n_grid*n_z];
  double* Vijk = new double[n_grid*n_grid*n_z];
  for(int nn=0;nn<n_its;nn++){
    for(int j = 0; j<n_grid;j++){
      for(int k=0; k<n_grid;k++){
        for(int z=0;z<n_z;z++){
          //computing consumption at each gridpoint
          double cons = z_vals[z]*pow((kgrid[j]),(alpha))+(1-delta)*kgrid[j]-kgrid[k];
          double cij;
          if(cons<0.0){
            cij = 0.0;
          }else{
            cij = cons;
          }
          //computing the expected value
          double cont_val_old;
          for(int cva = 0; cva<n_z;cva++){
            cont_val_old+=Vikold[k+cva*n_z]*t_mat[z+cva*n_z];
          }
          //computing the value function at each grid point
          Vijk[j+k*n_grid+z*n_grid*n_z] = crra(cij,sigma) + beta*cont_val_old;
          cont_val_old = 0.0;

          //finding the highest value and the capital choice associated with it
          int best_choice;
          double best_val;
          v_max(Vijk,n_grid*n_grid*n_z,n_grid,n_z,j,z,&best_val,&best_choice);
          Viknew[j+n_z*z] = best_val;
          kjknew[j+n_z*z] = kgrid[best_choice];
        }
      }
    }
    //computing the percentage change at each grid point
    double max_devs;
    max_devs = 0.0;
    for(int l=0; l<(n_grid*n_z);l++){
      double cur_devs = abs((Viknew[l]-Vikold[l]))/abs(Vikold[l]);
      if(cur_devs>max_devs){
        max_devs = cur_devs;
      }
    }
    if(nn < 2){
      max_devs = 100.0;
    }
    if(nn % print_every == 0){
      cout<<"Iterations = "<<nn+1<<endl;
      cout<<"Maximum Deviation = "<<max_devs<<endl;
    }
    //if iteration or convergence criteria are met, write value and policy
    //functions to .csv, delete allocated arrays, and break
    if((max_devs<1e-10)|(nn==n_its-1)){
      ///Write policy function to .csv
      ofstream vfi_file_k;
      vfi_file_k.open(kfilstring);
      ///Writing header
      //vfi_file_k<<"k,z,k'"<<endl;
      vfi_file_k<<"k,k'1,";
      for(int pz = 1;pz<n_z;pz++){
        if(pz == (n_z-1)){
          vfi_file_k<<"k'"<<n_z<<endl;
        }else{
          vfi_file_k<<"k'"<<(pz+1)<<",";
        }
      }
      ///Writing values
      //for(int mm=0;mm<n_grid;mm++){
        //for(int zz=0;zz<n_z;zz++){
            //vfi_file_k<<kgrid[mm]<<","<<z_vals[zz]<<","<<kjknew[mm+n_z*zz]<<endl;
        //}
      //}
      for(int mm = 0;mm<n_grid;mm++){
        vfi_file_k<<kgrid[mm]<<",";
        for(int zz= 0;zz<n_z;zz++){
          if(zz ==(n_z-1)){
            vfi_file_k<<kjknew[mm+n_z*zz]<<endl;
          }else{
            vfi_file_k<<kjknew[mm+n_z*zz]<<",";
          }
        }
      }
      vfi_file_k.close();
      ///Write value function to .csv
      ofstream vfi_file_v;
      vfi_file_v.open(vfilstring);
      vfi_file_v<<"k,v1,";
      for(int pz = 1;pz<n_z;pz++){
        if(pz == (n_z-1)){
          vfi_file_v<<"v"<<n_z<<endl;
        }else{
          vfi_file_v<<"v"<<(pz+1)<<",";
        }
      }
      //vfi_file_v<<"k,z,v"<<endl;
      //for(int mm=0;mm<n_grid;mm++){
        //for(int zz=0;zz<n_z;zz++){
          //vfi_file_v<<kgrid[mm]<<","<<z_vals[zz]<<","<<Viknew[mm+n_z*zz]<<endl;
        //}
      //}
      for(int mm = 0;mm<n_grid;mm++){
        vfi_file_v<<kgrid[mm]<<",";
        for(int zz= 0;zz<n_z;zz++){
          if(zz ==(n_z-1)){
            vfi_file_v<<Viknew[mm+n_z*zz]<<endl;
          }else{
            vfi_file_v<<Viknew[mm+n_z*zz]<<",";
          }
        }
      }
      vfi_file_v.close();
      ///clean heap
      delete[] kgrid;
      delete[] Vikold;
      delete[] Viknew;
      delete[] kjknew;
      delete[] Vijk;
      delete[] z_vals;
      delete[] t_mat;
      break;
    }
    for(int kk =0; kk<(n_grid*n_z);kk++){
        Vikold[kk] = Viknew[kk];
    }
  }
  return;
}


int main(){
  string file_loc = "Insert File";
  string* fl_point = &file_loc;
  double alpha = 0.3;
  double beta = 0.96;
  double delta = 0.08;
  double sigma = 2.0;
  int n_its = 1000;
  int n_grid = 100;
  double rho = 0.9;
  double mu = 0;
  double sigma_z = 0.1;
  int num_z = 4;
  double print_ev = 100;
  st_vfi_ncgm(alpha, beta, delta, sigma, n_its, n_grid, rho, mu, sigma_z, num_z, print_ev, fl_point);
  return 0;
}
