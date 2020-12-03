functions {
  //get tip means given rates, root trait value, tree info, and 'seed' for trait change along each edge
  vector get_X (int n, real X0, vector prune_T, vector R, vector raw_X,
                int[] preorder, int[] real_e, int[,] des_e, int[] tip_e){
    vector[2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] SS; //trait change along each edge
    XX[1] = X0;
    SS = rep_vector(0, 2 * n - 1);
    SS[real_e] = sqrt(prune_T[real_e] .* exp(R)) .* raw_X;
    for(i in preorder){
      XX[des_e[i, ]] = XX[i] + SS[des_e[i, ]];
    }
    return(XX[tip_e]);
  }
}
data {
  //basic data
  int obs; //number of observations
  int n; //number of tips
  int e; //number of edges
  vector[obs] Y; //observed trait values
  int X_id[obs]; //indicates to which tip each observation belongs
  matrix[e, e] eV; //edge variance-covariance matrix
  
  
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
  
  
  //prior specification: see below for parameter definitions
  real Ysig2_prior;
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
  int preorder[n - 1]; //'cladewise' sequence of nodes excluding tips, numbered by ancestral edge
  int lik_pow_ind; //indicates if lik_power is 0
  int has_tp; //any tip priors?
  int has_obs; //any observations?
  
  
  //for sampling from R prior
  if(!constr_Rsig2){
    chol_eV = cholesky_decompose(eV);
  }
  if(!constr_Rmu){
    T_midpts = diagonal(eV) + prune_T[real_e] / 6;
  }
  
  
  //for sampling from X prior
  {int counter;
  counter = 0;
  for(i in 1:(2 * n - 1)){
    if(des_e[i, 1] == -1){
      continue;
    }
    counter = counter + 1;
    preorder[counter] = i;
  }}
  
  
  if(lik_power == 0){
    lik_pow_ind = 0;
  }else{
    lik_pow_ind = 1;
  }
  if(n_tp == 0){
    has_tp = 0;
  }else{
    has_tp = 1;
  }
  if(obs == 0){
    has_obs = 0;
  }else{
    has_obs = 1;
  }
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real std_R0;
  real std_X0;
  real<lower=0> std_Rsig2[constr_Rsig2 ? 0:1];
  real std_Rmu[constr_Rmu ? 0:1];
  real<lower=0> std_Ysig2;
  vector[constr_Rsig2 ? 0:e] raw_R;
  vector[e] raw_X;
}

transformed parameters {
  real R0; //(ln)rate at root
  real X0; //trait value at root of tree
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  real<lower=0> Ysig2; //tip variance
  vector[e] R; //edge-wise average (ln)rates
  vector[n] X; //tip means
  vector[obs] cent_Y; //centered data
  vector[n_tp] trans_tp; //tip means centered and scaled w/ respect to tip priors
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * std_R0; //R0 prior: normal(R0_prior_mu, R0_prior_sig)
  X0 = X0_prior_mu + X0_prior_sig * std_X0; //X0 prior: normal(X0_prior_mu, X0_prior_sig)
	if(!constr_Rsig2){
	  Rsig2[1] = Rsig2_prior * std_Rsig2[1]; //Rsig2 prior: half-normal(0, Rsig2_prior)
	}
	if(!constr_Rmu){
	  Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * std_Rmu[1]; //Rmu prior: normal(Rmu_prior_mu, Rmu_prior_sig)
	}
	Ysig2 = Ysig2_prior * std_Ysig2; //Ysig2 prior: half-normal(0, Ysig2_prior)
  
  
  //R prior: multinormal(R0 + Rmu * T_midpts, Rsig2 * eV); Rsig2/Rmu = 0 when constrained
  R = rep_vector(R0, e);
  if(!constr_Rmu){
    R = R + Rmu[1] * T_midpts;
  }
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
  
  
  //X prior: multinormal(X0, tree topology with branch lengths multiplied by R)
  X = get_X(n, X0, prune_T, R, raw_X, preorder, real_e, des_e, tip_e);
  
  
  if(lik_pow_ind){
    //center tip means with priors based on tp_mu and scale based on tp_sig
    if(has_tp){
      trans_tp = (X[which_tp] - tp_mu) ./ tp_sig;
    }
    //center observations based on X
    if(has_obs){
      cent_Y = Y - X[X_id];
    }
  }
}

model {
  //high level priors
  std_R0 ~ std_normal();
  std_X0 ~ std_normal();
  if(!constr_Rsig2){
    std_Rsig2[1] ~ std_normal();
  }
  if(!constr_Rmu){
    std_Rmu[1] ~ std_normal();
  }
  std_Ysig2 ~ std_normal();
  
  
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//'seed' for sampling from X prior: see above
	raw_X ~ std_normal();
	
	
	//trait data
	if(lik_pow_ind){
	  //tip priors
	  //transform adjust unneeded since tp_sig is fixed and therefore constant with respect to params
	  if(has_tp){
	    target += -0.5 * lik_power * dot_self(trans_tp);
	  }
	  //likelihood of Y
	  if(has_obs){
	    target += -0.5 * lik_power * (obs * log(Ysig2) + dot_self(cent_Y) / Ysig2);
	  }
	}
}
