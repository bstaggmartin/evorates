data {
  //basic data
  int obs; //number of observations
  int n; //number of tips
  int n_e; //number of edges
  int n_obs[n]; //number of observations per tip
  vector[obs] X_obs; //vector of observed trait values
  matrix[n_e, n_e] eV; //edge variance-covariance matrix
  
  //data for pruning algorithm
  vector[2 * n - 1] prune_T; //vector of edge lengths for pruning-compatible (i.e, bifurcating) tree (+1 for stem edge)
  int des_e[2 * n - 1, 2]; //integer matrix indicating the pruning-compatible tree descendant edge numbers of each edge (+1 for stem edge)
  int tip_e[n]; //integer vector indicating pruning-compatible tree edge numbers giving rise to nodes in order of node numbers (+1 for stem edge)
  int real_e[n_e]; //integer vector indicating which pruning-compatible tree edges are 'real' and have associated R values (+1 for stem edge)
  int prune_seq[n - 1]; //integer vector indicating the sequence of edge numbers to iterate over for pruning algorithm (+1 for stem edge)
  
  //prior specification
  real X_prior;
  real Xsig2_prior;
  real R0_prior;
  real Rsig2_prior;
  real X0_prior;
  real Rmu_prior;
}

transformed data {
  //use cholesky decomposition of edge variance-covariance matrix and edge midpts to tranform raw_R to R (see below)
  matrix[n_e, n_e] chol_eV;
  vector[n_e] T_midpts;
  chol_eV = cholesky_decompose(eV);
  T_midpts = diagonal(eV) + prune_T[real_e] / 6;
}

parameters {
  vector[n] X; //mean trait values for each tip
  real<lower=0> Xsig2; //tip variance
  real R0; //rate at root of tree (log scale)
  real Rmu; //trend in rate (log scale)
  real<lower=0> Rsig2; //accumulation in rate variance per unit time
  real X0; //trait value at root of tree
  vector[n_e] raw_R; //vector of untransformed rate values along edges (log scale)
}

transformed parameters {
  vector[n_e] R; //vector of transformed rate values along edges (log scale)
  R = R0 + Rmu * T_midpts + sqrt(Rsig2) * chol_eV * raw_R; //implies prior on R to be multinormal(R0, Rsig2 * eV)
}

model {
  vector[2 * n - 1] SS; //vector of pruning-compatible tree edge lengths multiplied by R values (+1 for stem edge)
  vector[2 * n - 1] XX; //vector of node-wise trait values in order of their ancestral edge numbers (+1 for stem edge)
  vector[2 * n - 1] VV; //vector of node-wise trait variances in order of their ancestral edge numbers (+1 for stem edge)
  vector[n - 1] LL; //vector of log-likelihoods associated with node-wise contrast, in order of their ancestral edge numbers as ordered in prune_seq
  int counter; //temporary integer to indicate position along LL for pruning loop
  vector[2] des_X; //temporary vector to indicate trait values of descendant nodes for pruning loop
  vector[2] des_V; //temporary vector to indicate trait variances of descendant nodes for pruning loop
  
  //priors
  X ~ cauchy(0, X_prior); //prior on X
  Xsig2 ~ cauchy(0, Xsig2_prior); //prior on Xsig2
  R0 ~ cauchy(0, R0_prior); //prior on R0
  Rmu ~ cauchy(0, Rmu_prior); //prior on Rmu
  Rsig2 ~ cauchy(0, Rsig2_prior); //prior on Rsig2
  X0 ~ cauchy(0, X0_prior); //prior on x0
  raw_R ~ std_normal(); //implies prior on R to be multinormal(R0, Rsig2 * eV)
  
  //likelihood of X
  counter = 1;
  for (i in 1:n) {
    segment(X_obs, counter, n_obs[i]) ~ normal(X[i], sqrt(Xsig2));
    counter = counter + n_obs[i];
  }
  
  //pruning algorithm = calculate likelihood of estimated mean trait values given rate values
  SS = rep_vector(0, 2 * n - 1);
  SS[real_e] = prune_T[real_e] .* exp(R);
  XX[tip_e] = X;
  VV[tip_e] = rep_vector(0, n);
  counter = 0;
  for(i in prune_seq){
    des_X = XX[des_e[i, ]];
    des_V = VV[des_e[i, ]] + SS[des_e[i, ]];
    counter = counter + 1;
    LL[counter] = - 0.5 * (log(2 * pi()) + log(sum(des_V)) + (des_X[1] - des_X[2])^2 / sum(des_V));
    XX[i] = des_V[2] / sum(des_V) * des_X[1] + des_V[1] / sum(des_V) * des_X[2];
    VV[i] = 1 / (1 / des_V[1] + 1 / des_V[2]);
  }
  target += sum(LL) - 0.5 * (log(2 * pi()) + log(VV[1]) + (X0 - XX[1])^2 / sum(VV[des_e[1, ]]));
}