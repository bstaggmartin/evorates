functions {
  //get tip means given rates, root trait value, tree info, and 'seed' for trait change along each edge
  vector get_X (int n, real X0, vector edge_scalars, vector raw_X,
                int[] preorder, int[] real_e, int[,] des_e, int[] tip_e, real lik_power){
    vector[2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] SS; //trait change along each edge
    XX[1] = X0;
    SS = rep_vector(0, 2 * n - 1);
    SS[real_e] = sqrt(edge_scalars / lik_power) .* raw_X;
    for(i in preorder){
      XX[des_e[i, ]] = XX[i] + SS[des_e[i, ]];
    }
    return(XX[tip_e]);
  }
  //prior probability for X: pruning algorithm
  real prune(int n, vector edge_scalars, int[,] des_e, int[] tip_e, int[] real_e, int[] postorder, 
             vector X, real X0){
	  vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
    vector[2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] VV; //node-wise trait variances, indexed by ancestral edge
    vector[n - 1] LL; //log-likelihoods of node-wise contrasts, indexed by postorder sequence
    int counter; //position along LL
    vector[2] des_X; //temporary: descendant node trait values for given iteration in loop
    vector[2] des_V; //temporary: descendant node trait variances for given iteration in loop
    real sum_des_V; //temporary: sum of descendant node trait variances for given iteration in loop
	  SS = rep_vector(0, 2 * n - 1);
	  SS[real_e] = edge_scalars;
    XX[tip_e] = X;
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
	  return(sum(LL) - 0.5 * (log(VV[1]) + (X0 - XX[1])^2 / VV[1]));
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
  
  
  //tree data: note tree is coerced to be bifurcating and 1st edge is zero-length stem
  vector[2 * n - 1] prune_T; //edge lengths
  int des_e[2 * n - 1, 2]; //edge descendants (-1 for no descendants)
  int tip_e[n]; //edges with no descendants
  int real_e[e]; //non-zero length edges
  int postorder[n - 1]; //sequence of nodes (denoted by ancestral edge) to prune over
  
  
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
  vector[constr_Rmu ? 0:e] T_l; //"real" edge lengths
  vector[constr_Rmu ? 0:e] T_1; //overall 'height' of edge ancestors
  vector[constr_Rmu ? 0:e] T_2; //overall 'height' of edge descendants
  int preorder[n - 1]; //'cladewise' sequence of nodes excluding tips, numbered by ancestral edge
  int lik_pow_ind; //indicates if lik_power is 0
  int c_lik_pow_ind; //indicates if lik_power is 1 
  int has_tp; //any tip priors?
  int has_obs; //any observations?
  real lik_const; //additional normalizing constants for log likelihood
  real prior_const; //additional normalizing constant for log prior probability
  
  
  //for sampling from R prior
  if(!constr_Rsig2){
    chol_eV = cholesky_decompose(eV);
  }
  if(!constr_Rmu){
    T_l = prune_T[real_e];
    T_1 = diagonal(eV) - prune_T[real_e] / 3;
    T_2 = diagonal(eV) + 2 * prune_T[real_e] / 3;
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
  if(lik_power ==1){
    c_lik_pow_ind = 1;
  }else{
    c_lik_pow_ind = 0;
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
  
  
  lik_const = -0.5 * (obs + n) * log(2 * pi());
  prior_const = 0;
  if(has_tp){
    prior_const += -0.5 * (n_tp * log(2 * pi()) + sum(log(tp_sig)));
  }
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> std_R0;
  real<lower=-pi()/2, upper=pi()/2> std_X0;
  real<lower=0, upper=pi()/2> std_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> std_Rmu[constr_Rmu ? 0:1];
  real<lower=0, upper=pi()/2> std_Ysig2;
  vector[constr_Rsig2 ? 0:e] raw_R;
  vector[e] raw_X;
  
  
  vector[lik_pow_ind ? 0:n] no_lik_X;
}

transformed parameters {
  real R0; //(ln)rate at root
  real X0; //trait value at root of tree
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  real<lower=0> Ysig2; //tip variance
  vector[e] R; //edge-wise average (ln)rates
  vector[e] edge_scalars; //branch lengths scaled by R for "real" edges
  vector[n] X; //tip means
  real tmp_tp[has_tp ? 1:0]; //probability of tip means given tip priors
  real tmp_lik; //log likelihood of trait data (set to 0 if lik_power is 0)
  vector[e] tmp_mu; //expected edge-wise average (ln)rates given R0 and Rmu
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * tan(std_R0); //R0 prior: cauchy(R0_prior_mu, R0_prior_sig)
  X0 = X0_prior_mu + X0_prior_sig * tan(std_X0); //X0 prior: cauchy(X0_prior_mu, X0_prior_sig)
	if(!constr_Rsig2){
	  Rsig2[1] = Rsig2_prior * tan(std_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior)
	}
	if(!constr_Rmu){
	  Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * tan(std_Rmu[1]); //Rmu prior: cauchy(Rmu_prior_mu, Rmu_prior_sig)
	}
	Ysig2 = Ysig2_prior * tan(std_Ysig2); //Ysig2 prior: half-cauchy(0, Ysig2_prior)
  
  
  //R prior: multinormal(f(R0, Rmu, tree topology), Rsig2 * eV); Rsig2/Rmu = 0 when constrained
  //trend formula identical to that used for calculating edgewise rates under early burst model
  //formula undefined when Rmu = 0, but should never happen when Rmu is unconstrained, so no need to worry...
  R = rep_vector(R0, e);
  if(!constr_Rmu){
    R = R - log(fabs(Rmu[1])) - log(T_l) + log(fabs(exp(Rmu[1] * T_2) - exp(Rmu[1] * T_1)));
  }
  tmp_mu = R;
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
  
  
  //scale edges of tree
  edge_scalars = prune_T[real_e] .* exp(R);
  
  
  //sampling X: multinormal(X0, tree topology with branch lengths multiplied by R)
  if(lik_pow_ind){
    X = get_X(n, X0, edge_scalars, raw_X, preorder, real_e, des_e, tip_e, lik_power);
  }else{
    X = no_lik_X;
  }
  
  
  //center tip means with priors based on tp_mu and scale based on tp_sig
  //square using dot_self and multiply to get associated normal density (up to constant)
  //transform adjust unneeded since tp_sig is fixed and therefore constant with respect to params
  if(has_tp){
    tmp_tp[1] = -0.5 * dot_self((X[which_tp] - tp_mu) ./ tp_sig);
  }
  
  
  //trait data likelihood calculations
  //center observations based on X
  //multinormal component unneeded; the X sampling algorithm takes care of it
  if(has_obs && lik_pow_ind){
    tmp_lik = -0.5 * (obs * log(Ysig2) + dot_self(Y - X[X_id]) / Ysig2);
  }else{
    tmp_lik = 0;
  }
}

model {
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//'seed' for sampling X: see above
	raw_X ~ std_normal();
	
	
	//tip priors
	if(has_tp){
	  target += tmp_tp[1];
	}
	
	
	//add trait data likelihood to posterior density
	target += lik_power * tmp_lik;
	//adjust for X sampling if lik_power isn't 0 or 1
	if(lik_pow_ind && !c_lik_pow_ind){
	  target += -0.5 * (lik_power - 1) * sum(log(edge_scalars));
	}
}

generated quantities{
  //should this be shifted to the multinormal in the no_obs case???
  real lik; //log likelihood of trait data
  real prior; //log prior probability of parameter estimates
  
  
  lik = tmp_lik + lik_const;
  if(has_obs && !lik_pow_ind){
    lik += -0.5 * (obs * log(Ysig2) + dot_self(Y - X[X_id]) / Ysig2);
  }
  lik += prune(n, edge_scalars, des_e, tip_e, real_e, postorder, 
               X, X0);
  
  
  prior = prior_const;
  prior += cauchy_lpdf(R0 | R0_prior_mu, R0_prior_sig);
  prior += cauchy_lpdf(X0 | X0_prior_mu, X0_prior_sig);
  if(!constr_Rsig2){
    prior += cauchy_lpdf(Rsig2[1] | 0, Rsig2_prior) + log(2);
    prior += multi_normal_cholesky_lpdf(R | tmp_mu, sqrt(Rsig2[1]) * chol_eV);
  }
  if(!constr_Rmu){
    prior += cauchy_lpdf(Rmu[1] | Rmu_prior_mu, Rmu_prior_sig);
  }
  prior += cauchy_lpdf(Ysig2 | 0, Ysig2_prior) + log(2);
  if(has_tp){
    prior += tmp_tp[1];
  }
}