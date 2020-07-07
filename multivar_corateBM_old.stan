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
      counter = counter + 1;
    }
    return(X);
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
  int k_mis[k];
  int which_mis[sum(k_mis)]; 
  
  
  //data for pruning algorithm: note tree is coerced to be bifurcating and 1st edge is zero-length stem
  vector[2 * n - 1] prune_T; //edge lengths
  int des_e[2 * n - 1, 2]; //edge descendants (-1 for no descendants)
  int tip_e[n]; //edges with no descendants
  int real_e[e]; //non-zero length edges
  int postorder[n - 1]; //sequence of nodes (denoted by ancestral edge) to prune over
	
	
	//prior specification: see below for parameter definitions
	vector[k] Xsig2_prior;
	real Xcor_prior;
	real R0_prior;
	real Rsig2_prior;
	vector[k] X0_prior;
	real Rmu_prior;
	
	
	//parameter constraints
	int constr_Rsig2;
	int constr_Rmu;
}

transformed data {
  matrix[e, e] chol_eV; //cholesky decomp of edge variance-covariance matrix
  vector[e] T_midpts; //overall 'height' of edge mid-points
  
  //for sampling from R prior
  chol_eV = cholesky_decompose(eV);
  T_midpts = diagonal(eV) + prune_T[real_e] / 6;
}

parameters {
  //parameters on sampling scale: see below for parameter definitions
  real<lower=-pi()/2, upper=pi()/2> unif_R0;
  vector<lower=-pi()/2, upper=pi()/2>[k] unif_X0;
  real<lower=0, upper=pi()/2> unif_Rsig2[constr_Rsig2 ? 0:1];
  real<lower=-pi()/2, upper=pi()/2> unif_Rmu[constr_Rmu ? 0:1];
  vector<lower=0, upper=pi()/2>[k] unif_Xsig2;
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
  cov_matrix[k] Xcov;
  vector[e] R; //edge-wise average (ln)rates
  
  
  //high level priors
  R0 = R0_prior * tan(unif_R0); //R0 prior: cauchy(0, R0_prior)
  X0 = X0_prior .* tan(unif_X0); //X0 prior: cauchy(0, X0_prior)
	if(!constr_Rsig2){
	  Rsig2[1] = Rsig2_prior * tan(unif_Rsig2[1]); //Rsig2 prior: half-cauchy(0, Rsig2_prior)
	}
	if(!constr_Rmu){
	  Rmu[1] = Rmu_prior * tan(unif_Rmu[1]); //Rmu prior: cauchy(0, Rmu_prior)
	}
	Xsig2 = Xsig2_prior .* tan(unif_Xsig2); //Xsig2 prior: half-cauchy(0, Xsig2_prior)
  Xsig2 = Xsig2 / mean(Xsig2); //standardize Xsig2 to be mean 1 (prevent rate unidentifiability)
  
  
  //Xcov = sqrt(Xsig2') * Xcor * sqrt(Xsig2)
  chol_Xcov = diag_pre_multiply(sqrt(Xsig2), chol_Xcor);
  Xcov = chol_Xcor * chol_Xcor';
  
  
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
  //Xcor prior: LKJcorr(Xcor_prior)
	chol_Xcor ~ lkj_corr_cholesky(Xcor_prior);
  
  
  //'seed' for sampling from R prior: see above
	if(!constr_Rsig2){
	  raw_R ~ std_normal();
	}
	
	
	//likelihood of X: pruning algorithm
	{vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
  matrix[k, 2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
  matrix[k, k] VV[2 * n - 1]; //node-wise trait covariance matrices, indexed by ancestral edge
  vector[n - 1] LL; //log-likelihoods of node-wise contrasts, indexed by postorder sequence
  matrix[k, k] Xprec; //inverse of Xcov
  int counter; //position along LL
  matrix[k, 2] des_X; //temporary: descendant node trait values for given iteration in loop
  vector[k] contr;
  matrix[k, k] des_V[2]; //temporary: descendant node trait covariance matrices for given iteration in loop
  vector[2] des_S; //temporary: descendant node ancestral edge lengths * rates for given iteration in loop
	SS = rep_vector(0, 2 * n - 1);
	SS[real_e] = prune_T[real_e] .* exp(R);
  XX[, tip_e] = get_X(n, k, Y, mis_Y, k_mis, which_mis);
  VV[tip_e] = rep_array(rep_matrix(0, k, k), n);
  Xprec = inverse_spd(Xcov);
  counter = 0;
  for(i in postorder){
    des_X = XX[, des_e[i, ]];
    contr = des_X[, 2] - des_X[, 1];
    des_S = SS[des_e[i, ]];
    des_V[1] = VV[des_e[i, 1]] + des_S[1] * Xcov;
    des_V[2] = VV[des_e[i, 2]] + des_S[2] * Xcov;
    counter = counter + 1;
    LL[counter] = (k * (log(2 * pi())) + log_determinant(des_V[1] + des_V[2]) +
                  mdivide_right_spd(contr', des_V[1] + des_V[2]) * contr);
    XX[, i] = mdivide_right_spd(des_V[2], des_V[1] + des_V[2]) * des_X[, 1] + 
              mdivide_right_spd(des_V[1], des_V[1] + des_V[2]) * des_X[, 2];
    VV[i] = inverse_spd(inverse_spd(des_V[1]) + inverse_spd(des_V[2]));
  }
	target += sum(LL) - 
	          0.5 * (k* log(2 * pi()) + log_determinant(VV[1]) + 
	          (X0 - XX[, 1])' * inverse_spd(VV[des_e[1, 1]] + VV[des_e[1, 2]]) * (X0 - XX[, 1]));}
}