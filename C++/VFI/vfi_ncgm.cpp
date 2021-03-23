///Value Function Iteration in C++ (Neoclassical growth model)
///08/14/2020

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

double crra(double cons, double sigma){
  double val;
  if(sigma == 1){
    val = log(cons);
  }else{
    val = 1/((1-sigma)*pow(cons,sigma-1));
  }
  //cout<<val<<endl;
  return val;
}

double ss_cap(double alpha, double beta, double delta){
  double val = pow((alpha/((1/beta)+(delta-1))),(1/(1-alpha)));
  return val;
}

double max_cap(double alpha, double delta){
  double val =  pow((alpha/delta),(1/(1-alpha)));
  return val;
}

void v_max(double v_mat[],int entries,int cols, int select_row, double* max_val, int* max_loc){
  int rows = entries/cols;
  //cout<<rows<<endl;
  int cur_max_loc = 0;
  double cur_max_val;
  if((rows==1)|(cols==1)){
    cur_max_val = v_mat[0];
    for(int i = 0; i<max(rows,cols); i++){
      if(v_mat[i]>cur_max_val){
        cur_max_val = v_mat[i];
        cur_max_loc = i;
      }
    }
    *max_val = cur_max_val;
    *max_loc = cur_max_loc;
  }else{
    cur_max_val = v_mat[select_row];
    for(int i=0; i<cols; i++){
      if(v_mat[select_row+i*cols]>cur_max_val){
        cur_max_val = v_mat[select_row+i*cols];
        cur_max_loc = i;
      }
    }
    *max_val = cur_max_val;
    *max_loc = cur_max_loc;
  }
  return;
}


void vfi_ncgm(double alpha, double beta, double delta, double sigma,int n_its,int n_grid, bool ez_pref, int print_every, string* out_loc){
  string filstring;
  if(ez_pref==true){
    filstring = *out_loc+"vfi_ncgm_ez_"+to_string(n_grid)+".csv";
  }else{
    filstring = *out_loc+"vfi_ncgm_"+to_string(n_grid)+".csv";
  }
  double lowerk = 0.1*ss_cap(alpha,beta,delta);
  double maxk = max_cap(alpha, delta);
  double* kgrid = new double[n_grid];;
  kgrid[0] = lowerk;
  for(int i=1; i<n_grid; i++){
    kgrid[i] = kgrid[i-1] + maxk/(n_grid-1);
  }

  double* Viold = new double[n_grid];
  double* Vinew = new double[n_grid];
  double* kjnew = new double[n_grid];
  double* Vij = new double[n_grid*n_grid];
  for(int nn=0;nn<n_its;nn++){
    for(int j = 0; j<n_grid;j++){
      for(int k=0; k<n_grid;k++){
          double cons = pow((kgrid[j]),(alpha))+(1-delta)*kgrid[j]-kgrid[k];
          double cij;
          if(cons<0.0){
            cij = 0.0;
          }else{
            cij = cons;
          }
          if(ez_pref==false){
            Vij[j+k*n_grid] = crra(cij,sigma) + beta*Viold[k];
          }else{
            double pre;
            if(nn==0){
              pre = crra(cij,sigma);
            }else{
              pre = crra(cij,sigma) + beta*pow(Viold[k],(1-sigma));
            }
            //cout << pre << endl;
            Vij[j+k*n_grid] =pow(pre,(1/(1-sigma)));
          }


          int best_choice;
          double best_val;
          v_max(Vij,n_grid*n_grid,n_grid,j,&best_val,&best_choice);
          //cout<<j<<","<<k<<endl;
          Vinew[j] = best_val;
          kjnew[j] = kgrid[best_choice];
      }
    }
    double max_devs;
    max_devs = 0.0;
    for(int l=0; l<n_grid;l++){
      double cur_devs = abs((Vinew[l]-Viold[l]))/abs(Viold[l]);
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
    if((max_devs<1e-7)|(nn==n_its-1)){
      ///Write to .csv
      ofstream vfi_file;
      vfi_file.open(filstring);
      vfi_file<<"k,k',v(k)"<<endl;
      for(int mm=0;mm<n_grid;mm++){
          vfi_file<<kgrid[mm]<<","<<kjnew[mm]<<","<<Vinew[mm]<<endl;
      }
      vfi_file.close();
      ///clean heap
      delete[] kgrid;
      delete[] Viold;
      delete[] Vinew;
      delete[] kjnew;
      delete[] Vij;
      break;
    }
    for(int kk =0; kk<n_grid;kk++){
        Viold[kk] = Vinew[kk];
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
    int n_its = 10000;
    int n_grid = 1000;
    double print_ev = 10;
    bool ez = false;

    clock_t start, end;

    start = clock();
    ios_base::sync_with_stdio(false);
    vfi_ncgm(alpha, beta, delta, sigma, n_its, n_grid, ez, print_ev, fl_point);
    end = clock();
    double time_taken = double(end - start)/double(CLOCKS_PER_SEC);
    cout << fixed << time_taken << setprecision(5);
    cout << " sec elapsed " << endl;
    return 0;
}
