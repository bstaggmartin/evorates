functions {
    //function to combine fixed, observed and sampled, unobserved tip means
    vector get_X (int n, vector Y, vector mis_Y, int[] which_mis){
    vector[n] X;           
    X = Y;
    X[which_mis] = mis_Y;
    return(X);
  }
}

data {
  //basic data
  int n; //number of tips
  int e; //number of edges
  vector[n] Y; //observed trait values
  matrix[e, e] eV; //edge variance-covariance matrix
  
  
  //missing data handling
  int mis;
  int which_mis[mis];
  
  
  //tip priors
  int n_tp; //number of tip priors
  int which_tp[n_tp]; //indicates to which tip each prior belongs
  vector[n_tp] tp_mu; //tip prior means
  vector[n_tp] tp_sig; //tip prior sds
  
  
  //data for pruning algorithm: note tree is coerced to be bifurcating and 1st edge is zero-length stem
  vector[2 * n - 1] prune_T; //edge lengths
  int des_e[2 * n - 1, 2]; //edge descendants (-1 for no descendants)
  int tip_e[n]; //edges with no descendants
  int real_e[e]; //non-zero length edges
  int postorder[n - 1]; //sequence of nodes (denoted by ancestral edge) to prune over
	
	
	//prior specification: see below for parameter definitions
  real R0_prior_mu;
  real R0_prior_sig;
  real Rsig2_prior;
  real X0_prior_mu;
  real X0_prior_sig;
  real Rmu_prior_mu;
  real Rmu_prior_sig;
	
	
	//parameter constraints
	int constr_Rsig2;
	int constr_Rmu;
	
	
	//sampling from prior/data cloning
	real lik_power;
}

transformed data {
  matrix[constr_Rsig2 ? 0:e, constr_Rsig2 ? 0:e] chol_eV; //cholesky decomp of edge variance-covariance matrix
  vector[constr_Rmu ? 0:e] T_midpts; //overall 'height' of edge mid-points
  
  
  //for sampling from R prior
  if(!constr_Rsig2){
    chol_eV = cholesky_decompose(eV);
  }
  if(!constr_Rmu){
    T_midpts = diagonal(eV) + prune_T[real_e] / 6;
  }
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> unif_R0;
  real<lower=-pi()/2, upper=pi()/2> unif_X0;
  real<lower=0, upper=pi()/2> unif_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> unif_Rmu[constr_Rmu ? 0:1];
  vector[constr_Rsig2 ? 0:e] raw_R;
  
  
  vector[mis] mis_Y; //unobserved tip means
}

transformed parameters {
  real R0; //(ln)rate at root
  real X0; //trait value at root of tree
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  vector[e] R; //edge-wise average (ln)rates
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * tan(unif_R0); //R0 prior: cauchy(R0_prior_mu, R0_prior_sig)
  X0 = X0_prior_mu + X0_prior_sig * tan(unif_X0); //X0 prior: cauchy(X0_prior_mu, X0_prior_sig)
	if(!constr_Rsig2){
	  Rsig2[1] = Rsig2_prior * tan(unif_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior)
	}
	if(!constr_Rmu){
	  Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * tan(unif_Rmu[1]); //Rmu prior: cauchy(Rmu_prior_mu, Rmu_prior_sig)
	}
  
  
  //R prior: multinormal(R0 + Rmu * T_midpts, Rsig2 * eV); Rsig2/Rmu = 0 when constrained
  R = rep_vector(R0, e);
  if(!constr_Rmu){
    R = R + Rmu[1] * T_midpts;
  }
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
}

model {
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//tip priors
	if(lik_power != 0){
	  if(n_tp != 0){
	    (mis_Y[which_tp] - tp_mu) .* tp_sig ~ std_normal();
	  }
	}
	
	
	//likelihood of X: pruning algorithm
	if(lik_power != 0){
	  {vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
    vector[2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] VV; //node-wise trait variances, indexed by ancestral edge
    vector[n - 1] LL; //log-likelihoods of node-wise contrasts, indexed by postorder sequence
    int counter; //position along LL
    vector[2] des_X; //temporary: descendant node trait values for given iteration in loop
    vector[2] des_V; //temporary: descendant node trait variances for given iteration in loop
	  SS = rep_vector(0, 2 * n - 1);
	  SS[real_e] = prune_T[real_e] .* exp(R);
    XX[tip_e] = get_X(n, Y, mis_Y, which_mis);
    VV[tip_e] = rep_vector(0, n);
    counter = 0;
    for(i in postorder){
      des_X = XX[des_e[i, ]];
      des_V = VV[des_e[i, ]] + SS[des_e[i, ]];
      counter = counter + 1;
      LL[counter] = - 0.5 * (log(sum(des_V)) + (des_X[1] - des_X[2])^2 / sum(des_V));
      XX[i] = des_V[2] / sum(des_V) * des_X[1] + des_V[1] / sum(des_V) * des_X[2];
      VV[i] = 1 / (1 / des_V[1] + 1 / des_V[2]);
    }
	  target += lik_power * 
	            (sum(LL) - 0.5 * (n * log(2 * pi()) + log(VV[1]) + (X0 - XX[1])^2 / sum(VV[des_e[1, ]])));}
	}

}
