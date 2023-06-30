functions {
  //function to combine fixed and unfixed tip SEs
  vector get_SE (int n, vector p_SE, vector inv_n_obs, real[] Ysig2, int[] which_mis_SE, int has_mis_SE){
    vector[n] SE;           
    SE = p_SE;
    if(has_mis_SE){
      SE[which_mis_SE] = inv_n_obs[which_mis_SE] * Ysig2[1];
    }
    return(SE);
  }
  //likelihood of X: pruning algorithm
  real prune(int n, vector prune_T, int[,] des_e, int[] tip_e, int[] real_e, int[] postorder, 
             int[] mis_code, int[] which_non_mis,
             vector X, vector SE, vector R){
	  vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
    vector[2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] VV; //node-wise trait variances, indexed by ancestral edge
    vector[n - 1] LL; //log-likelihoods of node-wise contrasts, indexed by postorder sequence
    int counter; //position along postorder sequence
    int counter_non_mis; //position along which_non_mis
    vector[2] des_X; //temporary: descendant node trait values for given iteration in loop
    vector[2] des_V; //temporary: descendant node trait variances for given iteration in loop
    real sum_des_V; //temporary: sum of descendant node trait variances for given iteration in loop
	  SS = rep_vector(0, 2 * n - 1);
	  SS[real_e] = prune_T[real_e] .* exp(R);
    XX[tip_e] = X;
    VV[tip_e] = SE;
    counter = 0;
    counter_non_mis = 0;
    for(i in postorder){
      des_X = XX[des_e[i, ]];
      des_V = VV[des_e[i, ]] + SS[des_e[i, ]];
      counter = counter + 1;
      if(mis_code[counter] == 2){
        sum_des_V = sum(des_V);
        LL[counter] = log(sum_des_V) + (des_X[1] - des_X[2])^2 / sum_des_V;
        XX[i] = des_V[2] / sum_des_V * des_X[1] + des_V[1] / sum_des_V * des_X[2];
        //Stan can't divide by 0, so we do this manually
        if(des_V[1] == 0 || des_V[2] == 0){
          VV[i] = 0;
        }else{
          VV[i] = 1 / (1 / des_V[1] + 1 / des_V[2]);
        }
      }else if(mis_code[counter] == 1){
        LL[counter] = 0;
        counter_non_mis = counter_non_mis + 1;
        XX[i] = des_X[which_non_mis[counter_non_mis]];
        VV[i] = des_V[which_non_mis[counter_non_mis]];
      }else if(mis_code[counter] == 0){
        LL[counter] = 0;
        XX[i] = 0;
        VV[i] = -2;
      }
    }
	  return(-0.5 * sum(LL));
  }
}

data {
  //tree data
  int n; //number of tips
  int e; //number of edges
  matrix[e, e] eV; //edge variance-covariance matrix
  vector[2 * n - 1] prune_T; //edge lengths
  int des_e[2 * n - 1, 2]; //edge descendants (-1 for no descendants)
  int tip_e[n]; //edges with no descendants
  int real_e[e]; //edges with estimated rate values (non-zero length or tips)
  int postorder[n - 1]; //sequence of nodes (denoted by ancestral edge) to prune over
  
  
  //missing data handling
  //basic structure: a "code" for each node in postorder sequence indicating how many descendants have non-missing
  //trait info, e.g., a 0 indicates both descendants, 1 indicates 1 descendant, and 2 indicates no descendants.
  //n_mis_2 indicates the number of nodes in postorder sequence with "code" 2
  //n_mis_1 indicates the number of nodes in postorder sequence with "code" 1
  //which_non_mis indicates which descendant has non-missing trait info for each node in postorder sequence with "code" 1
  int mis_code[n - 1];
  int n_mis_2;
  int n_mis_1;
  int which_non_mis[n_mis_1];
  
  
  //trait information
  vector[n] X; //fixed tip means (set to 0 if missing)
  vector[n] p_SE; //fixed tip mean standard errors (set to -1 if unfixed, -2 if missing (i.e., infinite))
  vector[n] inv_n_obs; //the reciprocoal of the number of observations per tip (set to 0 if p_SE is fixed)
  int n_contra; //number of contrasts to calc likelihood of intra-tip observations
  vector[n_contra] contra; //contrasts used to calc likelihood of intra-tip observations
  vector[n_contra] contra_var; //variance scalars of contrasts used to calc likelihood of intra-tip observations
  //explanation: for a tip with l intra-tip observations, Y, given some estimated intra-tip variance Y_sig2,
  //the likelihood for all l observation while marginalizing over uncertainty in the tip mean, X, is:
  //sum(dnorm(Y[2:l] - cumsum(Y[1:(l - 1)]) / (1:(l - 1)), 0, sqrt(Y_sig2 * (2:l)/(1:(l - 1)))))
  //In other words, the mean of the normal dists are the difference between each observaiton (excluding the first)
  //and the cumulative mean of observations up to it, and the variances are decreasing scalars of Y_sig2, starting
  //at 2 and asymptotically approaching Y_sig2
  //The estimated "ancestral state" at the tip mean is simply the average of Y with variance Y_sig2 / l
  //target += -0.5 * (contra^2 / (contra_var * Y_sig2) + sum(log(contra_var * Y_sig2)))
  //scale_sq_contra = dot_self(contra / sqrt(contra_var))
  //sum_contra_var = sum(log(contra_var))
  // = -0.5 * (scale_sq_contra / Y_sig2 + sum_contra_var + n_contra * log(Y_sig2))
  //sum_contra_var can be dropped because it's a constant
  //we're left with target += -0.5 * (scale_sq_contra / Y_sig2 + n_contra * log(Y_sig2))
  
  
  //unfixed tip data handling
  int n_mis_SE; //number of tips with unfixed standard errors
  int which_mis_SE[n_mis_SE]; //which tips have unfixed standard errors
  
	
	//prior specification: see below for parameter definitions
  real R0_prior_mu;
  real R0_prior_sig;
  real Ysig2_prior_sig;
  real Ysig2_prior_df;
  real Rsig2_prior_sig;
  real Rmu_prior_mu;
  real Rmu_prior_sig;
	
	
	//parameter constraints
	int constr_Rsig2;
	int constr_Rmu;
	
	
	//sampling from prior/data cloning
	real lik_power;
}

transformed data {
  int has_mis_SE;
  int Ysig2_prior_cauchy;
  int Ysig2_prior_t;
  int needs_std_Ysig2;
  real half_df;
  int has_contra;
  real sum_sq_scale_contra;
  matrix[constr_Rsig2 ? 0:e, constr_Rsig2 ? 0:e] chol_eV; //cholesky decomp of edge variance-covariance matrix
  vector[constr_Rmu ? 0:e] T_l; //"real" edge lengths
  vector[constr_Rmu ? 0:e] T_1; //overall 'height' of edge ancestors
  vector[constr_Rmu ? 0:e] T_2; //overall 'height' of edge descendants
  int lik_pow_ind; //indicates if lik_power is NOT 0
  real lik_const; //additional normalizing constants for log likelihood
  
  Ysig2_prior_cauchy = 0;
  Ysig2_prior_t = 0;
  needs_std_Ysig2 = 0;
  if(n_mis_SE == 0){
    has_mis_SE = 0;
  }else{
    has_mis_SE = 1;
    if(Ysig2_prior_df == 1){
      Ysig2_prior_cauchy = 1;
    }else if(Ysig2_prior_df > 0){
      Ysig2_prior_t = 1;
      needs_std_Ysig2 = 1;
      half_df = Ysig2_prior_df / 2;
    }else{
      needs_std_Ysig2 = 1;
    }
  }
  
  
  //for calculating likelihood of intra-tip observations
  if(n_contra == 0){
    has_contra = 0;
  }else{
    has_contra = 1;
    sum_sq_scale_contra = dot_self(contra ./ sqrt(contra_var));
  }
  
  
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
  
  
  lik_const = -0.5 * (n_mis_2) * log(2 * pi());
  if(has_contra){
    lik_const += -0.5 * (n_contra * log(2 * pi()) + sum(log(contra_var)));
  }
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real std_R0;
  real<lower=0> tau_Ysig2[Ysig2_prior_t ? 1:0];
  real<lower=0> std_Ysig2[needs_std_Ysig2 ? 1:0];
  real<lower=0, upper=pi()/2> std_Ysig2_cauchy[Ysig2_prior_cauchy ? 1:0];
  real<lower=0, upper=pi()/2> std_Rsig2[constr_Rsig2 ? 0:1];
  real std_Rmu[constr_Rmu ? 0:1];
  vector[constr_Rsig2 ? 0:e] raw_R;
}

transformed parameters {
  real R0; //(ln)rate at root
  real Ysig2[has_mis_SE ? 1:0]; //intra-tip variance
  vector[n] SE; //combined fixed and unfixed tip standard errors
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  vector[e] R; //edge-wise average (ln)rates
  real tmp_lik; //log likelihood of trait data (set to 0 if lik_power is 0)
  vector[e] tmp_mu; //expected edge-wise average (ln)rates given R0 and Rmu
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * std_R0; //R0 prior: normal(R0_prior_mu, R0_prior_sig)
  if(has_mis_SE){
    if(Ysig2_prior_cauchy){
      Ysig2[1] = Ysig2_prior_sig * tan(std_Ysig2_cauchy[1]); //Ysig2 prior: half-cauchy(0, Ysig2_prior_sig)
    }else if(Ysig2_prior_t){
      Ysig2[1] = Ysig2_prior_sig * std_Ysig2[1] / sqrt(tau_Ysig2[1]); //Ysig2 prior: half-t(Ysig2_prior_df, 0, Ysig2_prior_sig)
    }else{
      Ysig2[1] = Ysig2_prior_sig * std_Ysig2[1]; //Ysig2 prior: half-normal(0, Ysig2_prior_sig)
    }
  }
	if(!constr_Rsig2){
	  Rsig2[1] = Rsig2_prior_sig * tan(std_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior_sig)
	}
	if(!constr_Rmu){
	  Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * std_Rmu[1]; //Rmu prior: normal(Rmu_prior_mu, Rmu_prior_sig)
	}
  
  
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
  
  
  //trait data likelihood calculations
	SE = get_SE(n, p_SE, inv_n_obs, Ysig2, which_mis_SE, has_mis_SE);
	if(lik_pow_ind){
	  tmp_lik = prune(n, prune_T, des_e, tip_e, real_e, postorder,
	                  mis_code, which_non_mis,
                    X, SE, R);
    if(has_contra){
      tmp_lik += -0.5 * (sum_sq_scale_contra / Ysig2[1] + n_contra * log(Ysig2[1]));
    }
	}else{
	  tmp_lik = 0;
	}
}

model {
  //high level priors
  std_R0 ~ std_normal();
  if(has_mis_SE){
    if(needs_std_Ysig2){
      std_Ysig2[1] ~ std_normal();
    }
    if(Ysig2_prior_t){
      tau_Ysig2[1] ~ gamma(half_df, half_df);
    }
  }
  if(!constr_Rmu){
    std_Rmu[1] ~ std_normal();
  }
  
  
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//add trait data likelihood to posterior density
	target += lik_power * tmp_lik;
}

generated quantities{
  real lik; //log likelihood of trait data
  real prior; //log prior probability of parameter estimates
  
  
  if(lik_pow_ind){
    lik = tmp_lik + lik_const;
  }else{
	  lik = prune(n, prune_T, des_e, tip_e, real_e, postorder, 
	              mis_code, which_non_mis,
                X, SE, R) + lik_const;
    if(has_contra){
      lik += -0.5 * (sum_sq_scale_contra / Ysig2[1] + n_contra * log(Ysig2[1]));
    }
  }
  
  
  prior = normal_lpdf(R0 | R0_prior_mu, R0_prior_sig);
  if(has_mis_SE){
    if(Ysig2_prior_cauchy){
      prior += cauchy_lpdf(Ysig2[1] | 0, Ysig2_prior_sig) + log(2);
    }else if(Ysig2_prior_t){
      prior += student_t_lpdf(Ysig2[1] | Ysig2_prior_df, 0, Ysig2_prior_sig) + log(2);
    }else{
      prior += normal_lpdf(Ysig2[1] | 0, Ysig2_prior_sig) + log(2);
    }
  }
  if(!constr_Rsig2){
    prior += cauchy_lpdf(Rsig2[1] | 0, Rsig2_prior_sig) + log(2);
    prior += multi_normal_cholesky_lpdf(R | tmp_mu, sqrt(Rsig2[1]) * chol_eV);
  }
  if(!constr_Rmu){
    prior += normal_lpdf(Rmu[1] | Rmu_prior_mu, Rmu_prior_sig);
  }
}
