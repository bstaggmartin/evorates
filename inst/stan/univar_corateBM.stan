data {
  //basic data
	int n; // number of tips
	int n_e; //number of edges
	vector[n] X; //vector of trait values at tips
	matrix[n_e, n_e] eV; //edge variance-covariance matrix
	
	//data for pruning algorithm
	vector[2 * n - 1] prune_T; //vector of edge lengths for pruning-compatible (i.e, bifurcating) tree (+1 for stem edge)
	int des_e[2 * n - 1, 2]; //integer matrix indicating the pruning-compatible tree descendant edge numbers of each edge (+1 for stem edge)
	int tip_e[n]; //integer vector indicating pruning-compatible tree edge numbers giving rise to nodes in order of node numbers (+1 for stem edge)
	int real_e[n_e]; //integer vector indicating which pruning-compatible tree edges are 'real' and have associated R values (+1 for stem edge)
	int prune_seq[n - 1]; //integer vector indicating the sequence of edge numbers to iterate over for pruning algorithm (+1 for stem edge)
}

transformed data {
  //use cholesky decomposition of edge variance-covariance matrix to tranform raw_R to R (see below)
  matrix[n_e, n_e] chol_eV;
  chol_eV = cholesky_decompose(eV);
}

parameters {
	real R0; //rate at root of tree (log scale)
	real<lower=0> Rsig2; //accumulation in rate variance per unit time
	vector[n_e] raw_R; //vector of untransformed rate values along edges (log scale)
	real X0; //trait value at root of tree
}

transformed parameters {
  vector[n_e] R; //vector of transformed rate values along edges (log scale)
  R = R0 + sqrt(Rsig2) * chol_eV * raw_R; //implies prior on R to be multinormal(R0, Rsig2 * eV)
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
	R0 ~ cauchy(0, 10); //prior on R0
	Rsig2 ~ cauchy(0, 20); //prior on Rsig2
	raw_R ~ std_normal(); //implies prior on R to be multinormal(R0, Rsig2 * eV)
	X0 ~ normal(0, 100); //prior on x0
	
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