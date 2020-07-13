model.code<-
"functions {
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
  Xcov = chol_Xcov * chol_Xcov';
  
  
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
	{matrix[k, k] Xprec; //inverse of Xcov
	vector[2 * n - 1] SS; //edge lengths multiplied by respective average rate
  matrix[k, 2 * n - 1] XX; //node-wise trait values, indexed by ancestral edge
  vector[2 * n - 1] VV; //node-wise trait variances, indexed by ancestral edge
  vector[n - 1] LL; //PARTIAL log-likelihoods of node-wise contrasts, indexed by postorder sequence
  int counter; //position along LL
  matrix[k, 2] des_X; //temporary: descendant node trait values for given iteration in loop
  vector[2] des_V; //temporary: descendant node trait variances for given iteration in loop
  vector[k] contr; //temporary: contrast between descendant node trait values
  Xprec = inverse_spd(Xcov);
	SS = rep_vector(0, 2 * n - 1);
	SS[real_e] = prune_T[real_e] .* exp(R) ;
  XX[, tip_e] = get_X(n, k, Y, mis_Y, k_mis, which_mis);
  VV[tip_e] = rep_vector(0, n);
  counter = 0;
  for(i in postorder){
    des_X = XX[, des_e[i, ]];
    des_V = VV[des_e[i, ]] + SS[des_e[i, ]];
    contr = des_X[, 2] - des_X[, 1];
    counter = counter + 1;
    LL[counter] = -0.5 * (k * log(sum(des_V)) + (contr' * Xprec * contr) / sum(des_V));
    XX[, i] = des_V[2] / sum(des_V) * des_X[, 1] + des_V[1] / sum(des_V) * des_X[, 2];
    VV[i] = 1 / (1 / des_V[1] + 1 / des_V[2]);
  }
	target += sum(LL) - 0.5 * (k * n * log(2 * pi()) + n * log_determinant(Xcov) + 
	          k * log(VV[1]) + ((X0 - XX[, 1])' * Xprec * (X0 - XX[, 1])) / VV[1]);}
}

"
#key insight: you WERE looking for a dotproduct function t(x) * inv(sigma) does get you the intermediate
#product with points by rows variable by columns--you were doing a bunch of extraneous calculations
#doing matrix multiplication--all you were trying to do was take the sum of the products for each row
#in the intermediate stage!

#scrap dirichlet distribution? Instead just do cauchy's?
#did the above
#actually, I think dirichlet might be the way to go
#no, it isn't
library(rstan)
library(phytools)
library(contSimmap)
my.model<-stan_model(model_code = model.code)
set.seed(123)
tree<-pbtree(n=50)
plot(tree)
# test<-gen.corateBM(tree,X0=rep(0,3),intra.var=T,n_obs=rep(4,length.out=length(tree$tip.label)),
#                    Xsig2=rWishart(1,3,diag(3))[,,1])
test<-gen.corateBM(tree,X0=rep(0,3),Xsig2=rWishart(1,3,diag(3))[,,1],Rsig2=0)
test
#sometimes last column doesn't show up with alpha-->just a graphical error
plot(test,tree,trait=1:test$k,lwd=3)

##INTRAVAR##
Y<-test$Y
#Y[cbind(sample(nrow(Y),100,replace=T),sample(ncol(Y),100,replace=T))]<-NA
is.obs<-!t(apply(Y,1,is.na))
codes<-apply(is.obs,1,function(ii) sum(ii*10^((length(ii)-1):0)))
if(any(codes==0)){
  Y<-Y[-which(codes==0),]
  is.obs<-is.obs[-which(codes==0),]
  codes<-codes[-which(codes==0)]
}
n<-length(tree$tip.label)
e<-nrow(tree$edge)
eV<-edge.vcv(tree)
k<-ncol(Y)
X_id<-match(rownames(Y),tree$tip.label)
#obs_code stuff...
code_sizes<-tapply(codes,codes,length)
obs_code<-unlist(lapply(names(code_sizes),function(ii) which(codes==as.numeric(ii))))
names(code_sizes)<-paste(lapply(k-nchar(names(code_sizes)),function(ii) paste(rep(0,ii),collapse='')),
                         names(code_sizes),sep='')
n_code<-length(code_sizes)
parsed.codes<-gregexpr('1',names(code_sizes))
code_ks<-lengths(parsed.codes)
code_key<-unlist(parsed.codes)
####

##NO INTRAVAR##
Y<-test$X
Y[cbind(sample(nrow(Y),20,replace=T),sample(ncol(Y),20,replace=T))]<-NA
Y<-Y[-sample(nrow(Y),1),]
missing.X<-names(which(sapply(tree$tip.label,function(ii) sum(rownames(Y)==ii))==0))
Y<-do.call(rbind,c(list(Y),setNames(rep(list(NA),length(missing.X)),missing.X)))
Y<-Y[tree$tip.label,]
is.obs<-!t(apply(Y,1,is.na))
k<-ncol(Y)
which_mis<-apply(Y,2,function(ii) which(is.na(ii)))
if(length(which_mis)==0){
  k_mis<-rep(0,k)
}else{
  k_mis<-lengths(which_mis)
  which_mis<-unlist(which_mis)
}
n<-length(tree$tip.label)
e<-nrow(tree$edge)
eV<-edge.vcv(tree)
####


o.hgt<-max(eV)
tree$edge.length<-tree$edge.length/o.hgt
eV<-eV/o.hgt
o.Xsig2<-rep(NA,k)
o.X0<-rep(NA,k)
for(i in 1:k){
  tmp.Y<-Y[,i]
  missing.dat<-names(which(tapply(tmp.Y,names(tmp.Y),function(ii) sum(!is.na(ii)))==0))
  missing.dat<-c(missing.dat,tree$tip.label[!(tree$tip.label%in%names(tmp.Y))])
  if(length(missing.dat)>0){
    tmp.tree<-ape::drop.tip(tree,missing.dat)
  }else{
    tmp.tree<-tree
  }
  tmp.X<-tapply(tmp.Y,names(tmp.Y),mean,na.rm=T)[tmp.tree$tip.label]
  o.Xsig2[i]<-sum(ape::pic(tmp.X,ape::multi2di(tmp.tree))^2)/n
  xx<-c(tmp.X,rep(NA,tmp.tree$Nnode))
  for(ee in nrow(tmp.tree$edge):0){
    if(ee==0){
      des<-which(tmp.tree$edge[,1]==n+1)
      xx[n+1]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
      break
    }
    nn<-tmp.tree$edge[ee,2]
    if(nn<=n){
      next
    }
    des<-which(tmp.tree$edge[,1]==nn)
    xx[nn]<-sum((1/tmp.tree$edge.length[des])/sum(1/tmp.tree$edge.length[des])*xx[tmp.tree$edge[des,2]])
  }
  o.X0[i]<-xx[n+1]
  Y[,i]<-(Y[,i]-o.X0[i])/sqrt(o.Xsig2[i])
}
plot(density(test$Y[,1]))
plot(density(Y[,1][!is.na(Y[,1])]))
#checks out

if(!ape::is.binary(tree)){
  poly.nodes<-which(sapply(1:max(tree$edge),function(ii) length(which(ii==tree$edge[,1])))>2)
  d_poly<-lapply(poly.nodes,function(ii) which(ii==tree$edge[,1]))
  tmp<-sort(unlist(lapply(d_poly,function(ii) ii[seq(2,length(ii)-1)])))
  tmp<-tmp+0:(length(tmp)-1)
  real_e<-(1:(2*n-2))[-tmp]+1
  tree<-multi2di(tree,random=F)
}else{
  real_e<-1:e+1
}
des_e<-sapply(1:nrow(tree$edge),function(ii) which(tree$edge[,1]==tree$edge[ii,2])+1)
tip_e<-which(lengths(des_e)==0)
tip_e<-tip_e[order(tree$edge[tip_e,2])]
des_e[tip_e]<-list(rep(-1,2))
des_e<-matrix(unlist(des_e),ncol=2,byrow=T)
root.edges<-which(tree$edge[,1]==n+1)+1
des_e<-rbind(root.edges,des_e)
prune_T<-c(0,tree$edge.length)
postorder<-((2*n-2):1)[-((2*n-2)-tip_e)]
tip_e<-tip_e+1
Y[!is.obs]<- 0
##FOR INTRAVAR##
dat<-list('obs'=nrow(Y),'n'=n,'e'=e,'X_id'=X_id,'Y'=t(Y),'eV'=eV,'k'=k,
          'obs_code'=obs_code,'n_code'=n_code,'code_sizes'=code_sizes,'code_ks'=as.array(code_ks),'code_key'=code_key,
          'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,
          'R0_prior'=10,'Rsig2_prior'=20,'X0_prior'=rep(100,k),'Rmu_prior'=20,
          'constr_Rsig2'=F,'constr_Rmu'=T,
          'Ysig2_prior'=rep(50,k),'Xsig2_prior'=rep(1,k),'Xcor_prior'=1,'Ycor_prior'=1)
##ELSE##
dat<-list('n'=n,'e'=e,'Y'=t(Y),'eV'=eV,'k'=k,
          'k_mis'=k_mis,'which_mis'=which_mis,
          'prune_T'=prune_T,'des_e'=des_e,'tip_e'=tip_e,'real_e'=real_e,'postorder'=postorder,
          'R0_prior'=10,'Rsig2_prior'=20,'X0_prior'=rep(100,k),'Rmu_prior'=20,
          'constr_Rsig2'=T,'constr_Rmu'=T,
          'Xsig2_prior'=rep(1,k),'Xcor_prior'=1)
####
ret<-sampling(object=my.model,data=dat,chains=1,refresh=1,iter=1000)
plot(fit@sim$samples[[1]]$`cent_Y[1,1]`~fit@sim$samples[[1]]$`X[1,1]`)
hmmm<-as.array(fit)[,1,]
hmmm[,grep('Ycov',dimnames(hmmm)[[2]])]
i<-1;j<-1
plot(hmmm[,paste('Ycov[',i,',',j,']',sep='')]*sqrt(o.Xsig2[i]*o.Xsig2[j]),type='l')
plot(hmmm[,'Xsig2[3]'],type='l')


Xcov<-array(0,dim=c(k,k,nrow(hmmm)))
for(i in 1:k){
  for(j in 1:k){
    #Xcov[i,j,]<-hmmm[,paste('chol_Xcov[',i,',',j,']',sep='')]
    Xcov[i,j,]<-hmmm[,paste('Xcov[',i,',',j,']',sep='')]
  }
}
#Xcov[,,]<-apply(apply(Xcov,3,function(ii) ii%*%t(ii)),2,matrix,nrow=k,ncol=k)
Xcov<-lapply(asplit(Xcov,3),function(ii) diag(sqrt(o.Xsig2))%*%ii%*%diag(sqrt(o.Xsig2)))
new.fac<-sapply(1:length(Xcov),function(ii) 1/mean(diag(Xcov[[ii]])))
Xcov<-lapply(1:length(Xcov),function(ii) Xcov[[ii]]*new.fac[ii])
Xcov<-array(unlist(Xcov),dim=c(k,k,nrow(hmmm)))
i<-1;j<-3
plot(Xcov[i,j,],type='l')
abline(h=test$Xsig2[i,j])
#seems good! Can't figure out how to back-transform yet, though...
#might have to create new factor to standardize diagonal of Xcov to 1 and add log(factor) to R...?
#Okay, standardizing Xcov makes sense, but adding the factors to the R estimates seems to make things
#worse--why?
#Ohhh, have to subtract from R since Xsig2 is being MULTIPLIED by some factor--seems to work...
R<-extract(fit,'R',permute=F)[,1,]-log(o.hgt)
R<-apply(R,2,function(ii) ii-log(new.fac))
plot(0,type='n',xlim=range(test$R),ylim=range(R),
     xlab='true ln(rate)',ylab='ln(rate) posterior distribution')
for(i in 1:nrow(tree$edge)){
  dens<-density(R[,i],bw=0.2)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.1*c(dens$y,-1*rev(dens$y))+test$R[i]
  polygon(yy~xx,border=NA,col='gray')
}
points(apply(R,2,mean)~test$R,pch=16)
abline(0,1,col='red')

out.X<-aperm(array(hmmm[,grep('^X\\[\\d+,\\d+\\]',colnames(hmmm))],dim=c(nrow(hmmm),k,n)),c(3,2,1))
for(i in 1:k){
  out.X[,i,]<-out.X[,i,]*sqrt(o.Xsig2[i])+o.X0[i]
}
out.X[,,1]
plot(as.vector(aperm(out.X,c(3,1,2)))~rep(as.vector(test$X[1:n,]),each=500))
abline(0,1)
#perfect!

Ycov<-array(0,dim=c(k,k,nrow(hmmm)))
for(i in 1:k){
  for(j in 1:k){
    Ycov[i,j,]<-hmmm[,paste('Ycov[',i,',',j,']',sep='')]*sqrt(o.Xsig2[i]*o.Xsig2[j])
  }
}
plot(Ycov[2,2,],type='l')
plot(R[,2],type='l')
abline(h=test$R[2])
#it all seems to be working, despite the max treedepth warnings! Weird...
#still got the mostly max tree depth exceeding with tighter priors--who knows? :/

plot(hmmm[,'R0']-log(o.hgt)-log(new.fac),type='l')
abline(h=test$R0)

axes<-c(1,2)
old.par<-par(no.readonly=T)
par(pty='s')
ellipse<-lapply(asplit(Xcov[axes,axes,],3),eigen)
true.ellipse<-eigen(test$Xsig2[axes,axes])
xx<-seq(0,2*pi,length.out=100)
coords<-vector('list',length(ellipse))
for(i in 1:length(ellipse)){
  coords[[i]]$x<-sqrt(ellipse[[i]]$values[1])*cos(xx)*ellipse[[i]]$vectors[1,1]+
    sqrt(ellipse[[i]]$values[2])*sin(xx)*ellipse[[i]]$vectors[2,1]
  coords[[i]]$y<-sqrt(ellipse[[i]]$values[1])*cos(xx)*ellipse[[i]]$vectors[1,2]+
    sqrt(ellipse[[i]]$values[2])*sin(xx)*ellipse[[i]]$vectors[2,2]
}
true.coords<-list(
  x=sqrt(true.ellipse$values[1])*cos(xx)*true.ellipse$vectors[1,1]+
    sqrt(true.ellipse$values[2])*sin(xx)*true.ellipse$vectors[2,1],
  y=sqrt(true.ellipse$values[1])*cos(xx)*true.ellipse$vectors[1,2]+
    sqrt(true.ellipse$values[2])*sin(xx)*true.ellipse$vectors[2,2]
)
plot(0,col='white',xlim=range(sapply(coords,'[','x')),ylim=range(sapply(coords,'[','y')),asp=1)
for(i in 1:length(coords)){
  lines(coords[[i]]$y~coords[[i]]$x,col=rgb(0,0,0,0.1))
}
lines(true.coords$y~true.coords$x,lwd=4,col='red')
par(old.par)

i<-3
hist(hmmm[,paste0('X0[',i,']')]*sqrt(o.Xsig2[i])+o.X0[i])
#7/3 tests:

#2000 iter for 50 tip tree with Rsig2 turned on, finished in ~30 minutes with no divergences or maximized
#tree depths--posterior looks good! max_treedepth set to 12

#k=5 stress test with 50 tips; finished in just under 15 minutes with about half of transitions exceeding
#max tree depth--posterior looks good, though (500 iter)

#n=100 stress test with 3 traits;took just over 30 minutes with no divergences and all transitions exceeding
#max tree depth--posterior still looks good (500 iter)

#Rsig2 turned on, Rmu  turned off, and max_treedepth=12 for the above 2 cases

#So, both tips and traits seem to affect the max_treedepth issue, with the number of species increasing
#mcmc  time seemingly far more than increasing the number of traits. In all cases, though, the posterior
#looks pretty good. Model seems to work, even it's slow and doesn't seem to play nice by stan's standards...

#just ran a quick test to see how the model would do with  a completely missing tip--some divergent
#transitions, but honestly, the posterior looked SHOCKINGLY good given that I accidentally fed it a Y mat
#with only 90/147 entries! Some Xsig2's were off, but this tree only had 15 TIPS! Also took like 5 mins to
#run

#really think I'm onto something here...

#7/6: alright--you probs can use uniform priors on tips in the non-intra-tip var cases, but otherwise
#you get severe problems due to the way you use transformation to sample tip values. I think this is
#fine for now, but it's worth considering in the future if it would be worth to make versions with intra-tip
#var where Xs are directly sampled and their likelihood calculated, so you could use bounded probability
#dists for tips

#7/6: with a lot of missing data, treedepth issues are still common in fixed tip mean case (and it's kinda
#slow). However, with no missing data, sampling goes extremely fast with no issues!

#Wrong likelihood function--now fixed
#something still wrong--root credible intervals are stupid wide, R values seem under-estimated
#something odd is going on here
#found a missing matrix inverse--think I got it!
#seems prone to model less evolutionary correlation than there really is, oddly enough...
#okay, after some more exploration--changing priors did little to change the resulting inference.
#probs something to do with the confounding influences of R and evolutionary correlation

#running for more iterations helps; also seems to match phylopars max likelihood ests

#7/7: As you suspected, there were ways to simplify the pruning alg given that you're assuming a constant
#correlation matrix--updated the code per Freckleton 2012's formulae; about to see how it works out!
#Beautiful--works perfectly, it seems

#missing data and/or model misspecification generally seems to erode support for evolutionary correlation
#and tend towards correctly estimating variances while underestimating covariance!