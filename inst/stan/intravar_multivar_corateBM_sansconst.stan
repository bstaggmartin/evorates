functions {
  //get tip means given rates, root trait value, tree data, 'seed' for trait change along each edge, and 
  //evolutionary covariance
  matrix get_X (int n, int k, vector X0, vector prune_T, vector R, matrix raw_X, matrix chol_Xcov,
                int[] preorder, int[] real_e, int[,] des_e, int[] tip_e){
    matrix[k, 2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
    matrix[k, 2 * n - 1] SS; //trait change along each edge
    XX[, 1] = X0;
    SS = rep_matrix(0, k, 2 * n - 1);
    for(i in 1:k){
      SS[i, real_e] = (sqrt(prune_T[real_e] .* exp(R)) .* to_vector(raw_X[i, ]))';
  }
  for(i in preorder){
    XX[, des_e[i, ]] = rep_matrix(XX[, i], 2) + chol_Xcov * SS[, des_e[i, ]];
  }
  return(XX[, tip_e]);
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
  vector[constr_Rmu ? 0:e] T_midpts; //overall 'height' of edge mid-points
  int preorder[n - 1]; //'cladewise' sequence of nodes excluding tips, numbered by ancestral edge
  
  
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
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> unif_R0;
  vector<lower=-pi()/2, upper=pi()/2>[k] unif_X0;
  real<lower=0, upper=pi()/2> unif_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> unif_Rmu[constr_Rmu ? 0:1];
  vector<lower=0, upper=pi()/2>[k] unif_Xsig2;
  vector<lower=0, upper=pi()/2>[k] unif_Ysig2;
  cholesky_factor_corr[k] chol_Xcor; //evolutionary correlation matrix
  cholesky_factor_corr[k] chol_Ycor; //tip correlation matrix
  vector[constr_Rsig2 ? 0:e] raw_R;
  matrix[k, e] raw_X;
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
  matrix[k, n] X; //tip means
  matrix[k, obs] cent_Y; //centered data
  vector[sum(n_tp)] trans_tp; //tip means centered and scaled w/ respect to tip priors
  
  
  //high level priors
  R0 = R0_prior_mu + R0_prior_sig * tan(unif_R0); //R0 prior: cauchy(R0_prior_mu, R0_prior_sig)
  X0 = X0_prior_mu + X0_prior_sig .* tan(unif_X0); //X0 prior: cauchy(X0_prior_mu, X0_prior_sig)
  if(!constr_Rsig2){
    Rsig2[1] = Rsig2_prior * tan(unif_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior)
  }
  if(!constr_Rmu){
    Rmu[1] = Rmu_prior_mu + Rmu_prior_sig * tan(unif_Rmu[1]); //Rmu prior: cauchy(Rmu_prior_mu, Rmu_prior_sig)
  }
  Xsig2 = Xsig2_prior .* tan(unif_Xsig2); //Xsig2 prior: half-cauchy(0, Xsig2_prior)
  Xsig2 = Xsig2 / mean(Xsig2); //standardize Xsig2 to be mean 1 (prevent rate unidentifiability)
  Ysig2 = Ysig2_prior .* tan(unif_Ysig2); //Ysig2 prior: half-cauchy(0, Ysig2_prior)
  
  
  //Xcov | Ycov = sqrt(Xsig2' | Ysig2') * (Xcor | Ycor) * sqrt(Xsig2 | Ysig2)
  chol_Xcov = diag_pre_multiply(sqrt(Xsig2), chol_Xcor);
  chol_Ycov = diag_pre_multiply(sqrt(Ysig2), chol_Ycor);
  Ycov = chol_Ycov * chol_Ycov';
	
	
  //R prior: multinormal(R0 + Rmu * T_midpts, Rsig2 * eV); Rsig2/Rmu = 0 when constrained
  R = rep_vector(R0, e);
  if(!constr_Rmu){
    R = R + Rmu[1] * T_midpts;
  }
  if(!constr_Rsig2){
    R = R + sqrt(Rsig2[1]) * chol_eV * raw_R;
  }
  
  
  //X prior: multinormal(X0, tree topology with branch lengths multiplied by R)
  X = get_X(n, k, X0, prune_T, R, raw_X, chol_Xcov, preorder, real_e, des_e, tip_e);
  
  
  //center observations based on X
  if(obs != 0){
    if(lik_power != 0){
      cent_Y = Y - X[, X_id];
    }
  }
  
  
  //center tip means with priors based on tp_mu and scale based on tp_sig
  if(sum(n_tp) != 0 && lik_power != 0){
    {int counter;
	  counter = 1;
	  for(i in 1:k){
		  if(!n_tp[i]){
		    continue;
		  }
		  trans_tp[1:n_tp[i] + (counter - 1)] =
		  (X[i, segment(which_tp, counter, n_tp[i])]' - segment(tp_mu, counter, n_tp[i])) ./
		  segment(tp_sig, counter, n_tp[i]);
	    counter = counter + n_tp[i];
	  }}
  }
}

model {
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
	//transform adjust unneeded since tp_sig is fixed and therefore constant with respect to params
	if(lik_power == 1){
	  trans_tp ~ std_normal();
	}else if(lik_power != 0){
	  target += -0.5 * lik_power * dot_self(trans_tp);
	}

	
	
	//likelihood of Y
	if(obs != 0){
	  if(lik_power != 0){
	    {int counter_sizes;
      int counter_ks;
      counter_sizes = 1;
      counter_ks = 1;
      for (i in 1:n_code) {
        //for each observation code, get corresponding observations and tip covariance matrix
        {int k_inds[code_ks[i]];
        matrix[code_ks[i], code_ks[i]] tmp_Ycov;
        matrix[code_ks[i], code_sizes[i]] tmp_cent_Y;
        k_inds = segment(code_key, counter_ks, code_ks[i]);
        tmp_Ycov = Ycov[k_inds, k_inds];
        tmp_cent_Y = cent_Y[k_inds, segment(obs_code, counter_sizes, code_sizes[i])];
        //centered obs.[!unobs. traits] ~ multinormal(0, Ycov[!unobs. traits, !unobs. traits])
        target += -0.5 * lik_power * (code_sizes[i] * log_determinant(tmp_Ycov) +
                  sum(rows_dot_product(mdivide_right_spd(tmp_cent_Y', tmp_Ycov), tmp_cent_Y')));
        counter_sizes = counter_sizes + code_sizes[i];
        counter_ks = counter_ks + code_ks[i];}
        }
      }
	  }
	}
}
