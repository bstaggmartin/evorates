functions {
  //get tip means given rates, root trait value, tree data, 'seed' for trait change along each edge, and 
  //evolutionary covariance
  matrix get_X (int n, int k, vector X0, vector edge_scalars, matrix raw_X, matrix chol_Xcov,
                int[] preorder, int[] real_e, int[,] des_e, int[] tip_e, real lik_power){
    matrix[k, 2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    matrix[k, 2 * n - 1] SS; //trait change along each edge
    XX[, 1] = X0;
    SS = rep_matrix(0, k, 2 * n - 1);
    for(i in 1:k){
      SS[i, real_e] = (sqrt(edge_scalars / lik_power) .* to_vector(raw_X[i, ]))';
    }
    for(i in preorder){
      XX[, des_e[i, ]] = rep_matrix(XX[, i], 2) + chol_Xcov * SS[, des_e[i, ]];
    }
    return(XX[, tip_e]);
  }
  //likelihood of X: pruning algorithm
  real prune(int n, int k, vector edge_scalars, int[,] des_e, int[] tip_e, int[] real_e, int[] postorder, 
             matrix X, vector X0, matrix chol_Xcov){
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
	  SS[real_e] = edge_scalars;
    XX[, tip_e] = X;
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
  real get_lik(int obs, int k, matrix Y, matrix X, int[] X_id, matrix Ycov,
               int n_code, int[] code_ks, int[] code_sizes, int[] code_key, int[] obs_code){
    matrix[k, obs] cent_Y; //centered data;
    int counter_sizes;
    int counter_ks;
    real out;
    //center observations based on X
    cent_Y = Y - X[, X_id];
    //likelihood of Y
    counter_sizes = 1;
    counter_ks = 1;
    out = 0;
    for (i in 1:n_code){
      //for each observation code, get corresponding observations and tip covariance matrix
      int k_inds[code_ks[i]];
      matrix[code_ks[i], code_ks[i]] tmp_Ycov;
      matrix[code_ks[i], code_sizes[i]] tmp_cent_Y;
      k_inds = segment(code_key, counter_ks, code_ks[i]);
      tmp_Ycov = Ycov[k_inds, k_inds];
      tmp_cent_Y = cent_Y[k_inds, segment(obs_code, counter_sizes, code_sizes[i])];
      //centered obs.[!unobs. traits] ~ multinormal(0, Ycov[!unobs. traits, !unobs. traits])
      out += code_sizes[i] * log_determinant(tmp_Ycov) +
             sum(rows_dot_product(mdivide_right_spd(tmp_cent_Y', tmp_Ycov), tmp_cent_Y'));
      counter_sizes = counter_sizes + code_sizes[i];
      counter_ks = counter_ks + code_ks[i];
    }
    return(-0.5 * out);
  }
}
data {
  //basic data
  int obs; //number of observations
  int n; //number of tips
  int e; //number of edges
  int k; //number of traits
  matrix[k, obs] Y; //observed trait values
  int X_id[obs]; //indicates to which tip each observation belongs
  matrix[e, e] eV; //edge variance-covariance matrix
  
  
  //missing data: 'observation codes' denote which traits are known for each observation
  int n_code; //number of unique observation codes
  int code_sizes[n_code]; //number of observations per observation code
  int obs_code[obs]; //observations corresponding to each observation code
  int code_ks[n_code]; //number of trait values per observation code
  int code_key[sum(code_ks)]; //trait values corresponding to each observation code
  
  
  //tip priors
  int n_tp[k]; //number of tip priors for each trait
  int which_tp[sum(n_tp)]; //indicates to which tip each prior belongs
  vector[sum(n_tp)] tp_mu; //tip prior means
  vector[sum(n_tp)] tp_sig; //tip prior sds
  
  
  //tree data: note tree is coerced to be bifurcating and 1st edge is zero-length stem
  vector[2 * n - 1] prune_T; //edge lengths
  int des_e[2 * n - 1, 2]; //edge descendants (-1 for no descendants)
  int tip_e[n]; //edges with no descendants
  int real_e[e]; //non-zero length edges
  int postorder[n - 1]; //sequence of nodes (denoted by ancestral edge) to prune over
  
  
  //prior specification: see below for parameter definitions
  vector[k] Xsig2_prior;
  vector[k] Ysig2_prior;
  real Xcor_prior;
  real Ycor_prior;
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
  if(sum(n_tp) == 0){
    has_tp = 0;
  }else{
    has_tp = 1;
  }
  if(obs == 0){
    has_obs = 0;
  }else{
    has_obs = 1;
  }
  
  
  lik_const = -0.5 * (sum(to_vector(code_sizes) .* to_vector(code_ks)) + n * k) * log(2 * pi());
  prior_const = 0;
  if(has_tp){
    prior_const += -0.5 * (sum(n_tp) * log(2 * pi()) + sum(log(tp_sig)));
  }
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> std_R0;
  vector<lower=-pi()/2, upper=pi()/2>[k] std_X0;
  real<lower=0, upper=pi()/2> std_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> std_Rmu[constr_Rmu ? 0:1];
  simplex[k] std_Xsig2;
  vector<lower=0, upper=pi()/2>[k] std_Ysig2;
  cholesky_factor_corr[k] chol_Xcor; //evolutionary correlation matrix
  cholesky_factor_corr[k] chol_Ycor; //tip correlation matrix
  vector[constr_Rsig2 ? 0:e] raw_R;
  matrix[k, e] raw_X;
  
  
  matrix[lik_pow_ind ? 0:k, lik_pow_ind ? 0:n] no_lik_X;
}

transformed parameters {
  real R0; //(ln)rate at root
  vector[k] X0; //trait value at root of tree
  real<lower=0> Rsig2[constr_Rsig2 ? 0:1]; //rate of (ln)rate variance accumulation
  real Rmu[constr_Rmu ? 0:1]; //trend in (ln)rate
  vector<lower=0>[k] Xsig2; //relative rates of trait evolution
  vector<lower=0>[k] Ysig2; //tip variances
  cholesky_factor_cov[k] chol_Xcov; //cholesky decomp of evolutionary covariance matrix
  cholesky_factor_cov[k] chol_Ycov; //cholesky decomp of tip covariance matrix
  cov_matrix[k] Ycov; //tip covariance matrix
  vector[e] R; //edge-wise average (ln)rates
  vector[e] edge_scalars; //branch lengths scaled by R for "real" edges
  matrix[k, n] X; //tip means
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
  Ysig2 = Ysig2_prior .* tan(std_Ysig2); //Ysig2 prior: half-cauchy(0, Ysig2_prior)
  
  
  //Xcov | Ycov = sqrt(Xsig2' | Ysig2') * (Xcor | Ycor) * sqrt(Xsig2 | Ysig2)
  chol_Xcov = diag_pre_multiply(sqrt(Xsig2), chol_Xcor);
  chol_Ycov = diag_pre_multiply(sqrt(Ysig2), chol_Ycor);
  Ycov = chol_Ycov * chol_Ycov';
	
	
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
  //well, the kronecker product of the tree topology with Xcov, but you get the idea!
  if(lik_pow_ind){
    X = get_X(n, k, X0, edge_scalars, raw_X, chol_Xcov, preorder, real_e, des_e, tip_e, lik_power);
  }else{
    X = no_lik_X;
  }
  
  
  //center tip means with priors based on tp_mu and scale based on tp_sig
  //square using dot_self and multiply to get associated normal density (up to constant)
  //transform adjust unneeded since tp_sig is fixed and therefore constant with respect to params
  if(has_tp){
    {int counter;
    counter = 1;
    for(i in 1:k){
      if(!n_tp[i]){
        continue;
      }
		  tmp_tp[1] += dot_self((X[i, segment(which_tp, counter, n_tp[i])]' - segment(tp_mu, counter, n_tp[i])) ./
		                        segment(tp_sig, counter, n_tp[i]));
		  counter = counter + n_tp[i];
		}}
		tmp_tp[1] = -0.5 * tmp_tp[1];
  }
  
  
  //trait data likelihood calculations
  if(has_obs && lik_pow_ind){
	  tmp_lik = get_lik(obs, k, Y, X, X_id, Ycov,
	                    n_code, code_ks, code_sizes, code_key, obs_code);
	}else{
	  tmp_lik = 0;
	}
}

model {
  //Xsig2 prior: dirichlet(Xsig2_pior)
  std_Xsig2 ~ dirichlet(Xsig2_prior);
  
  
  //Xcor | Ycor prior: LKJcorr(Xcor_prior | Ycor_prior)
	chol_Xcor ~ lkj_corr_cholesky(Xcor_prior);
	chol_Ycor ~ lkj_corr_cholesky(Ycor_prior);
  
  
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//'seed' for sampling from X prior: see above
	for(i in 1:k){
	  raw_X[i, ] ~ std_normal();
	}
	
	
	//tip priors
	if(has_tp){
	  target += tmp_tp[1];
	}
	
	
	//add trait data likelihood to posterior density
	target += lik_power * tmp_lik;
	//adjust for X sampling if lik_power isn't 0 or 1
	if(lik_pow_ind && !c_lik_pow_ind){
	  target += (1 - lik_power) * (e * sum(log(diagonal(chol_Xcov))) + 0.5 * k * sum(log(edge_scalars)));
	}
}

generated quantities{
  //should this be shifted to the multinormal in the no_obs case???
  real lik; //log likelihood of trait data
  real prior; //log prior probability of parameter estimates
  
  
  lik = tmp_lik + lik_const;
  if(has_obs && !lik_pow_ind){
    lik += get_lik(obs, k, Y, X, X_id, Ycov,
	                 n_code, code_ks, code_sizes, code_key, obs_code);
  }
  lik += prune(n, k, edge_scalars, des_e, tip_e, real_e, postorder, 
               X, X0, chol_Xcov);
  
  
  prior = prior_const;
  prior += cauchy_lpdf(R0 | R0_prior_mu, R0_prior_sig);
  for(i in 1:k){
    prior += cauchy_lpdf(X0[i] | X0_prior_mu[i], X0_prior_sig[i]);
    prior += cauchy_lpdf(Ysig2[i] | 0, Ysig2_prior[i]) + log(2);
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
  prior += lkj_corr_lpdf(chol_Ycor * chol_Ycor' | Ycor_prior);
  if(has_tp){
    prior += tmp_tp[1];
  }
}