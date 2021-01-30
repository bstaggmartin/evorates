functions {
  //function to combine fixed, observed and sampled, unobserved tip means
  matrix get_X (int n, int k, matrix Y, vector mis_Y, int[] k_mis, int[] which_mis){
    matrix[k, n] X;
    int counter;
    X = Y;
    counter = 1;
    for(i in 1:k){
      if(!k_mis[i]){
        continue;
      }
      X[i, segment(which_mis, counter, k_mis[i])] = segment(mis_Y, counter, k_mis[i])';
      counter = counter + k_mis[i];
    }
    return(X);
  }
  //likelihood of X: pruning algorithm
  real prune(int n, int k, vector prune_T, int[,] des_e, int[] tip_e, int[] real_e, int[] postorder, 
             matrix Y, int[] k_mis, int[] which_mis, vector mis_Y, vector X0, matrix chol_Xcov, vector R){
	  matrix[k, k] Xprec; //inverse of Xcov
	  vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
    matrix[k, 2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    vector[2 * n - 1] VV; //node-wise trait variances, indexed by ancestral edge
    vector[n - 1] LL; //PARTIAL log-likelihoods of node-wise contrasts, indexed by postorder sequence
    int counter; //position along LL
    matrix[k, 2] des_X; //temporary: descendant node trait values for given iteration in loop
    vector[2] des_V; //temporary: descendant node trait variances for given iteration in loop
    real sum_des_V; //temporary: sum of descendant node trait variances for given iteration in loop
    vector[k] contr; //temporary: contrast between descendant node trait values
    Xprec = inverse_spd(chol_Xcov * chol_Xcov');
	  SS = rep_vector(0, 2 * n - 1);
	  SS[real_e] = prune_T[real_e] .* exp(R) ;
    XX[, tip_e] = get_X(n, k, Y, mis_Y, k_mis, which_mis);
    VV[tip_e] = rep_vector(0, n);
    counter = 0;
    for(i in postorder){
      des_X = XX[, des_e[i, ]];
      des_V = VV[des_e[i, ]] + SS[des_e[i, ]];
      sum_des_V = sum(des_V);
      contr = des_X[, 2] - des_X[, 1];
      counter = counter + 1;
      LL[counter] = -0.5 * (k * log(sum_des_V) + (contr' * Xprec * contr) / sum_des_V);
      //equivalent to eq'n A.4 in Caetano et al. 2019 when R matrices are the same b/t branches
      XX[, i] = des_V[2] / sum_des_V * des_X[, 1] + des_V[1] / sum_des_V * des_X[, 2];
      //make sure variance of an observed node stays 0!
      if(des_V[1] == 0 || des_V[2] == 0){
        VV[i] = 0;
      }else{
        VV[i] = 1 / (1 / des_V[1] + 1 / des_V[2]);
      }
    }
    return(sum(LL) - n * sum(log(diagonal(chol_Xcov))) -
           0.5 * (k * log(VV[1]) + ((X0 - XX[, 1])' * Xprec * (X0 - XX[, 1])) / VV[1]));
  }
}

data {
  //basic data
  int n; //number of tips
  int e; //number of edges
  int k; //number of traits
  matrix[k, n] Y; //observed trait values
  matrix[e, e] eV; //edge variance-covariance matrix
  
  
  //missing data handling
  int k_mis[k]; //number of missing observations per trait
  int which_mis[sum(k_mis)];  //indicates which tips are missing observations for each trait
  
  
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
	vector[k] Xsig2_prior;
	real Xcor_prior;
  real R0_prior_mu;
  real R0_prior_sig;
  real Rsig2_prior;
  vector[k] X0_prior_mu;
  vector[k] X0_prior_sig;
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
  int has_tp; //any tip priors?
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
  
  
  lik_const = -0.5 * n * k * log(2 * pi());
  prior_const = 0;
  if(has_tp){
    prior_const = -0.5 * (sum(k_mis) * log(2 * pi()) + sum(log(tp_sig)));
  }
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> std_R0;
  vector<lower=-pi()/2, upper=pi()/2>[k] std_X0;
  real<lower=0, upper=pi()/2> std_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> std_Rmu[constr_Rmu ? 0:1];
  simplex[k] std_Xsig2;
  cholesky_factor_corr[k] chol_Xcor; //evolutionary correlation matrix
  vector[constr_Rsig2 ? 0:e] raw_R;
  
  
  vector[sum(k_mis)] mis_Y; //unobserved tip values
}

transformed parameters {
  real R0; //(ln)rate at root
  vector[k] X0; //trait value at root of tree
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  vector<lower=0>[k] Xsig2; //relative rates of trait evolution
  cholesky_factor_cov[k] chol_Xcov; //cholesky decomp of evolutionary covariance matrix
  vector[e] R; //edge-wise average (ln)rates
  real tmp_tp[has_tp ? 1:0]; //probability of tip means given tip priors
  real tmp_lik; //log likelihood of trait data (set to 0 if lik_power is 0)
  vector[e] tmp_mu; //expected edge-wise average (ln)rates given R0 and Rmu
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * tan(std_R0); //R0 prior: cauchy(R0_prior_mu, R0_prior_sig)
  X0 = X0_prior_mu + X0_prior_sig .* tan(std_X0); //X0 prior: cauchy(X0_prior_mu, X0_prior_sig)
  if(!constr_Rsig2){
    Rsig2[1] = Rsig2_prior * tan(std_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior)
  }
  if(!constr_Rmu){
    Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * tan(std_Rmu[1]); //Rmu prior: cauchy(Rmu_prior_mu, Rmu_prior_sig)
  }
	Xsig2 = std_Xsig2 * k; //standardize Xsig2 to be mean 1
  
  
  //Xcov = sqrt(Xsig2') * Xcor * sqrt(Xsig2)
  chol_Xcov = diag_pre_multiply(sqrt(Xsig2), chol_Xcor);
  
  
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
  

  //center tip means with priors based on tp_mu and scale based on tp_sig
  //square using dot_self and multiply to get associated normal density (up to constant)
  //transform adjust unneeded since tp_sig is fixed and therefore constant with respect to params
  if(has_tp){
    tmp_tp[1] = -0.5 * dot_self((mis_Y[which_tp] - tp_mu) ./ tp_sig);
  }
  
  
  //trait data likelihood calculations
	if(lik_pow_ind){
	  tmp_lik = prune(n, k, prune_T, des_e, tip_e, real_e, postorder,
	                  Y, k_mis, which_mis, mis_Y, X0, chol_Xcov, R);
	}else{
	  tmp_lik = 0;
	}
}

model {
  //Xsig2 prior: dirichlet(Xsig2_prior)
  std_Xsig2 ~ dirichlet(Xsig2_prior);
  
  
  //Xcor prior: LKJcorr(Xcor_prior)
	chol_Xcor ~ lkj_corr_cholesky(Xcor_prior);
  
  
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//tip priors
	if(has_tp){
	  target += tmp_tp[1];
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
    lik = prune(n, k, prune_T, des_e, tip_e, real_e, postorder,
                Y, k_mis, which_mis, mis_Y, X0, chol_Xcov, R) +
          lik_const;
  }
  
  
  prior = prior_const;
  prior += cauchy_lpdf(R0 | R0_prior_mu, R0_prior_sig);
  for(i in 1:k){
    prior += cauchy_lpdf(X0[i] | X0_prior_mu[i], X0_prior_sig[i]);
  }
  if(!constr_Rsig2){
    prior += cauchy_lpdf(Rsig2[1] | 0, Rsig2_prior) + log(2);
    prior += multi_normal_cholesky_lpdf(R | tmp_mu, sqrt(Rsig2[1]) * chol_eV);
  }
  if(!constr_Rmu){
    prior += cauchy_lpdf(Rmu[1] | Rmu_prior_mu, Rmu_prior_sig);
  }
  prior += lkj_corr_lpdf(chol_Xcor * chol_Xcor' | Xcor_prior);
  prior += dirichlet_lpdf(std_Xsig2 | Xsig2_prior);
  if(has_tp){
    prior += tmp_tp[1];
  }
}
