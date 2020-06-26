data {
  //basic data
  int obs; //number of observations
	int n; // number of tips
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
  int which_obs[obs];
  int which_unobs[n - obs];
  int obs_counter;
  int unobs_counter;
  chol_eV = cholesky_decompose(eV);
  T_midpts = diagonal(eV) + prune_T[real_e] / 6;
  obs_counter = 0;
  unobs_counter = 0;
  for(i in 1:n){
    if(!n_obs[i]){
      unobs_counter = unobs_counter + 1;
      which_unobs[unobs_counter] = i;
      continue;
    }
    obs_counter = obs_counter + 1;
    which_obs[obs_counter] = i;
  }
}

parameters {
  vector[n - obs] X_unobs;
	real R0; //rate at root of tree (log scale)
	real Rmu[constr_Rmu ? 0 : 1]; //trend in rate (log scale)
	real<lower=0> Rsig2[constr_Rsig2 ? 0 : 1]; //accumulation in rate variance per unit time
	real X0; //trait value at root of tree
	vector[constr_Rsig2 ? 0 : n_e] raw_R; //vector of untransformed rate values along edges (log scale)
}

transformed parameters {
  vector[n] X; //vector of trait values at tips
  vector[n_e] R; //vector of transformed rate values along edges (log scale)
  //implies prior on R to be multinormal(R0 + Rmu * T_midpts, Rsig2 * eV), with Rmu and Rsig2 set to 0 as appropriate
  X[which_obs] = X_obs;
  X[which_unobs] = X_unobs;
  R = rep_vector(R0, n_e);
  if(!constr_Rmu){
    R = R + Rmu[1] * T_midpts;
  }
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
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
	R0 ~ cauchy(0, R0_prior); //prior on R0
	if(!constr_Rmu){
	  Rmu ~ cauchy(0, Rmu_prior); //prior on Rmu
	}
	if(!constr_Rsig2){
	  Rsig2 ~ cauchy(0, Rsig2_prior); //prior on Rsig2
	  raw_R ~ std_normal(); //implies prior on R to be multinormal(R0 + Rmu * T_midpts, Rsig2 * eV), with Rmu and Rsig2 set to 0 as appropriate
	}
	X0 ~ cauchy(0, X0_prior); //prior on x0
	
	//pruning algorithm = calculate likelihood of trait data given rate values
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