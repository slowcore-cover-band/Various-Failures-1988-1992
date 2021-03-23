///Sudoku Solver in C++
///12/24/20


#include <iostream>
#include <iomanip>
using namespace std;

bool sud_check(int grid[][9],int rows, int cols, int pro_val){
  int sr = 3;
  int sc = 3;
  if(rows>5){
    sr = 6;
  }else if(rows<3){
    sr = 0;
  }
  if(cols>5){
    sc = 6;
  }else if(cols<3){
    sc = 0;
  }
  //check values above and below
  for(int ii=0; ii<9; ii++){
    if(grid[ii][cols]==pro_val){
      return false;
    }
  }
  //check values to left and right
  for(int jj = 0; jj<9; jj++){
    if(grid[rows][jj]==pro_val){
      return false;
    }
  }
  //check box
  for(int kk = sr; kk<(sr+3); kk++){
    for(int ll = sc; ll<(sc+3); ll++){
      if(grid[kk][ll]== pro_val){
        return false;
      }
    }
  }
  return true;
}

void sud_solver(int grid[][9]){
  for(int rr=0; rr<9; rr++){
    for(int cc=0; cc<9; cc++){
      if(grid[rr][cc]==0){
        for(int vv = 1; vv<10; vv++){
          if(sud_check(grid,rr,cc,vv)){
            grid[rr][cc] = vv;
            //cout<< grid[rr][cc]<<endl;
            sud_solver(grid);
            grid[rr][cc] = 0;
          }
        }
        return;
      }
    }
  }
  cout<< grid[0][0] << grid[0][1] << grid[0][2] << "|" << grid[0][3] << grid[0][4] << grid[0][5] << "|" << grid[0][6] << grid[0][7] << grid[0][8]<< endl;
  cout<< grid[1][0] << grid[1][1] << grid[1][2] << "|" << grid[1][3] << grid[1][4] << grid[1][5] << "|" << grid[1][6] << grid[1][7] << grid[1][8]<< endl;
  cout<< grid[2][0] << grid[2][1] << grid[2][2] << "|" << grid[2][3] << grid[2][4] << grid[2][5] << "|" << grid[2][6] << grid[2][7] << grid[2][8]<< endl;
  cout<<"-----------"<<endl;
  cout<< grid[3][0] << grid[3][1] << grid[3][2] << "|" << grid[3][3] << grid[3][4] << grid[3][5] << "|" << grid[3][6] << grid[3][7] << grid[3][8]<< endl;
  cout<< grid[4][0] << grid[4][1] << grid[4][2] << "|" << grid[4][3] << grid[4][4] << grid[4][5] << "|" << grid[4][6] << grid[4][7] << grid[4][8]<< endl;
  cout<< grid[5][0] << grid[5][1] << grid[5][2] << "|" << grid[5][3] << grid[5][4] << grid[5][5] << "|" << grid[5][6] << grid[5][7] << grid[5][8]<< endl;
  cout<<"-----------"<<endl;
  cout<< grid[6][0] << grid[6][1] << grid[6][2] << "|" << grid[6][3] << grid[6][4] << grid[6][5] << "|" << grid[6][6] << grid[6][7] << grid[6][8]<< endl;
  cout<< grid[7][0] << grid[7][1] << grid[7][2] << "|" << grid[7][3] << grid[7][4] << grid[7][5] << "|" << grid[7][6] << grid[7][7] << grid[7][8]<< endl;
  cout<< grid[8][0] << grid[8][1] << grid[8][2] << "|" << grid[8][3] << grid[8][4] << grid[8][5] << "|" << grid[8][6] << grid[8][7] << grid[8][8]<< endl;
  cout<<""<<endl;
}

int main(){
  string rowvals;
  int grid[9][9];
  for(int er = 1; er<10; er++){
    cout<< "Enter Row #" << er<<":"<< endl;
    cin>>rowvals;
    for(unsigned ro = 0; ro < rowvals.length(); ro++){
      int curVal = rowvals[ro]-48;
      grid[er-1][ro] = curVal;
      //cout<<curVal;
    }
    cout<<endl;
  }
  clock_t start, end;

  start = clock();
  ios_base::sync_with_stdio(false);
  sud_solver(grid);
  end = clock();

  double time_taken = double(end - start)/double(CLOCKS_PER_SEC);
    cout << fixed << time_taken << setprecision(5);
    cout << " sec elapsed " << endl;
  return 0;
}
