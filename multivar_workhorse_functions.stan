





  {int counter_sizes;
  int counter_ks;
  counter_sizes = 1;
  counter_ks = 1;
  for (i in 1:n_code) {
    {int k_inds[code_ks[i]];
    matrix[code_ks[i], code_ks[i]] Yprec;
    matrix[code_ks[i], code_sizes[i]] tmp_cent_Y;
    vector[code_sizes[i]] LL_last_term;
    k_inds = segment(code_key, counter_ks, code_ks[i]);
    Yprec = inverse_spd(Ycov[k_inds, k_inds]);
    tmp_cent_Y = cent_Y[k_inds, segment(obs_code, counter_sizes, code_sizes[i])];
    for (j in 1:code_sizes[i]){
      LL_last_term[j] = tmp_cent_Y[,j]' * Yprec * tmp_cent_Y[,j];
}
target += -0.5 * (code_sizes[i] * code_ks[i] * (log(2 * pi())) + 
                    code_sizes[i] * -log_determinant(Yprec) + 
                    sum(LL_last_term));}
    counter_sizes = counter_sizes + code_sizes[i];
    counter_ks = counter_ks + code_ks[i];