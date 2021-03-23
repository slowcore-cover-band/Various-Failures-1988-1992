///Rouwenhorst Method
///08/15/20

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

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

int main() {
    int num_z = 6;
    double* t_mat = new double[num_z*num_z];
    double* z_vals = new double[num_z];
    rouwenhorst(0.9,0.1,0,t_mat,num_z*num_z,z_vals,num_z);
    return 0;
}
