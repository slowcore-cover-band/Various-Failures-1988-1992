###Sudoku Solver
###12/24/20
library(tictoc)
check_fn = function(grid,row,col,pro_val){
  if(length(which(grid[row,]==pro_val)) ==0){
    if(length(which(grid[,col]==pro_val)) ==0){
      sr = 0
      sc = 0
      if(row<4){
        sr = 1
      }else if(row>6){
        sr = 7
      }else{
        sr = 4
      }
      if(col<4){
        sc = 1
      }else if(col>6){
        sc = 7
      }else{
        sc = 4
      }
      if(length(which(grid[sr:(sr+2),sc:(sc+2)]==pro_val))==0){
        return(TRUE)
      }else{return(FALSE)}
    }else{return(FALSE)}
  }else{return(FALSE)}
}

pretty_grid = function(grid){
  print(paste(grid[1,1],grid[1,2],grid[1,3],"|",grid[1,4],grid[1,5],grid[1,6],"|",grid[1,7],grid[1,8],grid[1,9]))
  print(paste(grid[2,1],grid[2,2],grid[2,3],"|",grid[2,4],grid[2,5],grid[2,6],"|",grid[2,7],grid[2,8],grid[2,9]))
  print(paste(grid[3,1],grid[3,2],grid[3,3],"|",grid[3,4],grid[3,5],grid[3,6],"|",grid[3,7],grid[3,8],grid[3,9]))
  print("---------------------")
  print(paste(grid[4,1],grid[4,2],grid[4,3],"|",grid[4,4],grid[4,5],grid[4,6],"|",grid[4,7],grid[4,8],grid[4,9]))
  print(paste(grid[5,1],grid[5,2],grid[5,3],"|",grid[5,4],grid[5,5],grid[5,6],"|",grid[5,7],grid[5,8],grid[5,9]))
  print(paste(grid[6,1],grid[6,2],grid[6,3],"|",grid[6,4],grid[6,5],grid[6,6],"|",grid[6,7],grid[6,8],grid[6,9]))
  print("---------------------")
  print(paste(grid[7,1],grid[7,2],grid[7,3],"|",grid[7,4],grid[7,5],grid[7,6],"|",grid[7,7],grid[7,8],grid[7,9]))
  print(paste(grid[8,1],grid[8,2],grid[8,3],"|",grid[8,4],grid[8,5],grid[8,6],"|",grid[8,7],grid[8,8],grid[8,9]))
  print(paste(grid[9,1],grid[9,2],grid[9,3],"|",grid[9,4],grid[9,5],grid[9,6],"|",grid[9,7],grid[9,8],grid[9,9]))
  print("")
  print("")
}

solve_fn = function(){
  for(i in 1:9){
    for(j in 1:9){
      if(grid[i,j]==0){
        for(k in 1:9){
          if(check_fn(grid,i,j,k)){
            grid[i,j] <<- k
            solve_fn()
            grid[i,j] <<- 0
          }
        }
        return("DONE!")
      }
    }
  }
  pretty_grid(grid)
}

time_sf = function(){
  tic()
  solve_fn()
  toc()
}
