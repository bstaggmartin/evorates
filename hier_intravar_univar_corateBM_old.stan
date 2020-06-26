functions {
  vector get_X (int n, real X0, vector prune_T, vector R, vector raw_X, int[] XX_seq, int[] real_e, int[,] des_e, int[] tip_e){
    vector[2 * n - 1] XX; //vector of node-wise trait values in order of their ancestral edge numbers (+1 for stem edge)
    vector[2 * n - 1] SS;
    XX[1] = X0;
    SS = rep_vector(0, 2 * n - 1);
    SS[real_e] = sqrt(prune_T[real_e] .* exp(R)) .* raw_X;
    for(i in XX_seq){
      XX[des_e[i, ]] = XX[i] + SS[des_e[i, ]];
    }
    return(XX[tip_e]);
  }
}
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
  
  //optional parameters
	int constr_Rsig2;
	int constr_Rmu;
}

transformed data {
  //use cholesky decomposition of edge variance-covariance matrix and edge midpts to tranform raw_R to R (see below)
  matrix[n_e, n_e] chol_eV;
  vector[n_e] T_midpts;
  int XX_seq[n - 1];
  int XX_counter;
  chol_eV = cholesky_decompose(eV);
  T_midpts = diagonal(eV) + prune_T[real_e] / 6;
  XX_counter = 0;
  for(i in 1:(2 * n - 1)){
    if(des_e[i, 1] == -1){
      continue;
    }
    XX_counter = XX_counter + 1;
    XX_seq[XX_counter] = i;
  }
}

parameters {
  vector[n_e] raw_X; //edge-wise deviations from ancestral values
  real<lower=0> Xsig2; //tip variance
  real R0; //rate at root of tree (log scale)
  real Rmu[constr_Rmu ? 0 : 1]; //trend in rate (log scale)
  real<lower=0> Rsig2[constr_Rsig2 ? 0 : 1]; //accumulation in rate variance per unit time
  real X0; //trait value at root of tree
  vector[constr_Rsig2 ? 0 : n_e] raw_R; //vector of untransformed rate values along edges (log scale)
}

transformed parameters {
  vector[n_e] R; //vector of transformed rate values along edges (log scale)
  vector[n] X; //mean trait values for each tip
  //implies prior on R to be multinormal(R0 + Rmu * T_midpts, Rsig2 * eV), with Rmu and Rsig2 set to 0 as appropriate
  R = rep_vector(R0, n_e);
  if(!constr_Rmu){
    R = R + Rmu[1] * T_midpts;
  }
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
  X = get_X(n, X0, prune_T, R, raw_X, XX_seq, real_e, des_e, tip_e);
}

model {
  int counter; //temporary integer to indicate position along LL for pruning loop
  
  //priors
  raw_X ~ std_normal(); //implies prior on X to be be multinormal...
  Xsig2 ~ cauchy(0, Xsig2_prior); //prior on Xsig2
  R0 ~ cauchy(0, R0_prior); //prior on R0
  if(!constr_Rmu){
	  Rmu ~ cauchy(0, Rmu_prior); //prior on Rmu
	}
	if(!constr_Rsig2){
	  Rsig2 ~ cauchy(0, Rsig2_prior); //prior on Rsig2
	  raw_R ~ std_normal(); //implies prior on R to be multinormal(R0 + Rmu * T_midpts, Rsig2 * eV), with Rmu and Rsig2 set to 0 as appropriate
	}
  X0 ~ cauchy(0, X0_prior); //prior on x0
  
  //likelihood of X
  counter = 1;
  for (i in 1:n) {
    if(!n_obs[i]){
      continue;
    }
    segment(X_obs, counter, n_obs[i]) ~ normal(X[i], sqrt(Xsig2));
    counter = counter + n_obs[i];
  }
}
