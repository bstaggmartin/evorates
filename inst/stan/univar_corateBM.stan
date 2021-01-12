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
  int k_mis;
  int which_mis[k_mis];
  
  
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
  vector[constr_Rmu ? 0:e] T_l; //"real" edge lengths
  vector[constr_Rmu ? 0:e] T_1; //overall 'height' of edge ancestors
  vector[constr_Rmu ? 0:e] T_2; //overall 'height' of edge descendants
  int lik_pow_ind; //indicates if lik_power is 0
  int has_tp; //has tip priors?
  
  
  //for sampling from R prior
  if(!constr_Rsig2){
    chol_eV = cholesky_decompose(eV);
  }
  if(!constr_Rmu){
    T_l = prune_T[real_e];
    T_1 = diagonal(eV) - prune_T[real_e] / 3;
    T_2 = diagonal(eV) + 2 * prune_T[real_e] / 3;
  }
  
  
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
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> std_R0;
  real<lower=-pi()/2, upper=pi()/2> std_X0;
  real<lower=0, upper=pi()/2> std_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> std_Rmu[constr_Rmu ? 0:1];
  vector[constr_Rsig2 ? 0:e] raw_R;
  
  
  vector[k_mis] mis_Y; //unobserved tip means
}

transformed parameters {
  real R0; //(ln)rate at root
  real X0; //trait value at root of tree
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  vector[e] R; //edge-wise average (ln)rates
  vector[n_tp] trans_tp; //tip means centered and scaled w/ respect to tip priors
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * tan(std_R0); //R0 prior: cauchy(R0_prior_mu, R0_prior_sig)
  X0 = X0_prior_mu + X0_prior_sig * tan(std_X0); //X0 prior: cauchy(X0_prior_mu, X0_prior_sig)
	if(!constr_Rsig2){
	  Rsig2[1] = Rsig2_prior * tan(std_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior)
	}
	if(!constr_Rmu){
	  Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * tan(std_Rmu[1]); //Rmu prior: cauchy(Rmu_prior_mu, Rmu_prior_sig)
	}
  
  
  //R prior: multinormal(R0 + Rmu * T_midpts, Rsig2 * eV); Rsig2/Rmu = 0 when constrained
  //trend now calc'd differently to reflect AM, rather than GM; follows early burst model formula
  R = rep_vector(R0, e);
  if(!constr_Rmu){
    R = R - log(fabs(Rmu[1])) - log(T_l) + log(fabs(exp(Rmu[1] * T_2) - exp(Rmu[1] * T_1)));
  }
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
  
  
  //center tip means with priors based on tp_mu and scale based on tp_sig
  if(has_tp){
    trans_tp = (mis_Y[which_tp] - tp_mu) ./ tp_sig;
  }
}

model {
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//tip priors
  //transform adjust unneeded since tp_sig is fixed and therefore constant with respect to params
	if(has_tp != 0){
	 target += -0.5 * dot_self(trans_tp);
	}
	
	
	//trait data
	if(lik_pow_ind){
	  //likelihood of X: pruning algorithm
	  {vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
    vector[2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] VV; //node-wise trait variances, indexed by ancestral edge
    vector[n - 1] LL; //log-likelihoods of node-wise contrasts, indexed by postorder sequence
    int counter; //position along LL
    vector[2] des_X; //temporary: descendant node trait values for given iteration in loop
    vector[2] des_V; //temporary: descendant node trait variances for given iteration in loop
    real sum_des_V; //temporary: sum of descendant node trait variances for given iteration in loop
	  SS = rep_vector(0, 2 * n - 1);
	  SS[real_e] = prune_T[real_e] .* exp(R);
    XX[tip_e] = get_X(n, Y, mis_Y, which_mis);
    VV[tip_e] = rep_vector(0, n);
    counter = 0;
    for(i in postorder){
      des_X = XX[des_e[i, ]];
      des_V = VV[des_e[i, ]] + SS[des_e[i, ]];
      sum_des_V = sum(des_V);
      counter = counter + 1;
      LL[counter] = - 0.5 * (log(sum_des_V) + (des_X[1] - des_X[2])^2 / sum_des_V);
      XX[i] = des_V[2] / sum_des_V * des_X[1] + des_V[1] / sum_des_V * des_X[2];
      //make sure variance of an observed node stays 0!
      if(des_V[1] == 0 || des_V[2] == 0){
        VV[i] = 0;
      }else{
        VV[i] = 1 / (1 / des_V[1] + 1 / des_V[2]);
      }
    }
	  target += lik_power * (sum(LL) - 0.5 * (log(VV[1]) + (X0 - XX[1])^2 / VV[1]));}
	}
}
