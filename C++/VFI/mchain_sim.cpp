///Sampling from a Markov Chain
///08/17/20

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
