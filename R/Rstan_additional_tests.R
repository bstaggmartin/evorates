stan.model <- 
"functions {
	matrix get_PVCV(int n, vector T, matrix V, vector R) {
	  matrix[n, n] VCV;
	  matrix[2 * n - 2, 2 * n - 2] rates;
	  rates = diag_matrix(exp(R)) * diag_matrix(T);
	  VCV = V * rates * V';
	  return VCV;
	}
	real precomp_multinorm(int k, vector dat, vector mu, matrix inv, real det, real c) {
	  vector[k] dat_centered;
	  real L;
	  dat_centered = dat - mu;
	  L = -0.5 * dat_centered' * (1/c * inv) * dat_centered - log(((2 * pi())^(k / 2)) * sqrt(c^k * det));
	  return L;
	}
}
data {
	int n;	    			  // number of tips
	vector[n] X;	    			  // vector of trait values at tips
	vector[2 * n - 2] T;	    			  // vector of branch lengths
	matrix[n, 2 * n - 2] V;	    			  // phylogeny in matrix form
	matrix[2 * n - 2, 2 * n - 2] eV;	    			  // edge variance-covariance matrix
}
transformed data {
  matrix[2 * n - 2, 2 * n - 2] inv_eV;
  real det_eV;
  inv_eV = inverse(eV);
  det_eV = determinant(eV);
}
parameters {
	real<lower=0> R0;	    			  // rate at root of tree
	real<lower=0> Rsig2;	    			  // accumulation in rate variance per unit time
	vector[2 * n - 2] R;	    			  // vector of rate values along edges (log scale)
	real X0;	    			  // trait value at root of tree
}
transformed parameters {

}
model {
	R0 ~ lognormal(0, 7);	    			  // prior on R0 *note that exp(7) is about 1000
	Rsig2 ~ exponential(sqrt(1e-3));	    			  // prior on Rsig2
	target += precomp_multinorm(2 * n - 2, R, rep_vector(log(R0), 2 * n - 2), inv_eV, det_eV, Rsig2); // prior on rate values along edges (on log scale)
	X0 ~ normal(0, 1e3);	    			  // prior on x0
	X ~ multi_normal(rep_vector(X0, n), get_PVCV(n, T, V, R));	    			  // likelihood function
}"
my.mod <- stan_model(model_code=stan.model)

library(phytools)
library(mvnfast)
#function to get phylogenetic 'indicator matrix' (can be multiplied with branch lengths to get phylogenetic variance-covariance matrix)
get.V<-function(tree){
  V<-matrix(nrow=length(tree$tip.label),ncol=nrow(tree$edge))
  for(i in nrow(tree$edge):1){
    d.edges<-which(tree$edge[,1]==tree$edge[i,2])
    if(length(d.edges)==0){
      V[,i]<-c(rep(0,tree$edge[i,2]-1),1,rep(0,length(tree$tip.label)-tree$edge[i,2]))
    }else{
      V[,i]<-as.logical(rowSums(V[,d.edges]))
    }
  }
  V
}
#function to get phylogenetic variance-covariance matrix for edge-wise means of trait evolving under BM
edge.vcv<-function(tree){
  ancs<-vector(mode='list',length=nrow(tree$edge))
  for(i in 1:nrow(tree$edge)){
    adj.edge<-which(tree$edge[,2]==tree$edge[i,1])
    if(length(adj.edge)>0){
      ancs[[i]]<-c(ancs[[adj.edge]],i)
    }else{
      ancs[[i]]<-i
    }
  }
  for(i in 1:length(ancs)){
    ancs[[i]]<-ancs[[i]][-which(ancs[[i]]==i)]
  }
  edge.hgts<-node.depth.edgelength(tree)[tree$edge[,2]]
  mrcas<-mrca(tree,full=T)
  mat<-matrix(NA,nrow=nrow(tree$edge),ncol=nrow(tree$edge))
  for(i in 1:nrow(tree$edge)){
    for(j in 1:nrow(tree$edge)){
      if(i==j){
        mat[i,j]<-edge.hgts[i]-2*tree$edge.length[i]/3
      }else{
        des1<-tree$edge[i,2];des2<-tree$edge[j,2]
        tmp.mrca<-mrcas[des1,des2]
        tmp.edge<-which(tree$edge[,2]==tmp.mrca)
        anc.des<-ifelse(j%in%ancs[[i]]|i%in%ancs[[j]],1,0)
        if(length(tmp.edge)>0){
          mat[i,j]<-edge.hgts[tmp.edge]-anc.des*tree$edge.length[tmp.edge]/2
        }else{
          mat[i,j]<-0
        }
      }
    }
  }
  mat
}
set.seed(321)
n<-50
tree<-pbtree(n=n,b=2,d=1.5,scale=1,extant.only=T,nsim=20)
tree<-tree[lengths(tree)>1][[1]]
TT<-tree$edge.length
H<-node.depth.edgelength(tree)[tree$edge[,2]]-tree$edge.length/2
V<-get.V(tree)
eV<-edge.vcv(tree)


####pick parameter values####

R0<-1
Rsig2<-0.5
if(is.infinite(Rsig2)){
  R<-rep(R0,2*n-2)
}else{
  R<-as.vector(rmvn(1,rep(log(R0),2*n-2),Rsig2*eV))
}
X0<-0


###simulate and view rest of data####

X<-as.vector(rmvn(1,rep(X0,n),V%*%diag(exp(R)*TT)%*%t(V)))
colramp<-colorRampPalette(c('blue','red'))(100)
colvec<-colramp[round((R-min(R))/(max(R)-min(R))*99)+1]
tmp<-fastAnc(tree,X)
plot(c(X,tmp)~node.depth.edgelength(tree),col='white')
segments(x0=node.depth.edgelength(tree)[tree$edge[,1]],x1=node.depth.edgelength(tree)[tree$edge[,2]],
         y0=c(X,tmp)[tree$edge[,1]],y1=c(X,tmp)[tree$edge[,2]],col=colvec,lwd=2)
exp(max(R)-min(R))
#~20-fold difference in rate, tree seems to display signals of rate variation


####run model####

my.data <- list("n" = n,
                "X" = X,
                "T" = TT,
                "V" = V,
                "eV" = eV)
fit <- sampling(object = my.mod,
                data = my.data,
                iter = 2000,chains = 1,
                control=list(adapt_delta=0.8,max_treedepth=10))


#looking at the branch-wise rates, which are the main parameter of interest for me
rates<-extract(fit,"R",permute=FALSE,inc_warmup=FALSE)

#sample a random branch and plot the trace
branch<-sample(2*n-2,1)
plot(rates[,1,branch],type='l',main=paste('ln(rate) for branch',branch))
abline(h=R[branch],col='red')
##actually looks pretty good, minus a few brief intervals where the effective step size was reduced (look ~iteration 2500ish)

#look at overall correlation
plot(0,type='n',xlim=range(R),ylim=range(rates),xlab='true ln(rate)',ylab='ln(rate) posterior distribution')
for(i in 1:(dim(rates)[3])){
  dens<-density(rates[,1,i],bw=0.2)
  yy<-c(dens$x,rev(dens$x))
  xx<-0.1*c(dens$y,-1*rev(dens$y))+R[i]
  polygon(yy~xx,border=NA,col='gray')
}
points(as.vector(apply(rates,c(2,3),median))~R,pch=16)
abline(0,1,col='red')
#points(R~H,col='red')
##again, this time I finally seemed to get a decent result (on the one I'm showing to you, go figure)
##credible intervals are wide, but that's to be expected with such an over-parameterized model in my experience

#even when credible intervals overlap, the differences between branch-wise rates in any given sample may look decent
hist(rates[,1,which.max(R)]-rates[,1,which.min(R)],col='black',breaks=50)
abline(v=R[which.max(R)]-R[which.min(R)],col='red')
branches<-sample(2*n-2,2)
hist(rates[,1,branches[1]]-rates[,1,branches[2]],col='black',breaks=50)
abline(v=R[branches[1]]-R[branches[2]],col='red')
##it's alright generally speaking, not great, but not terrible given what I've seen before out of this model

#looking at trace of other parameters
est.Rsig2<-extract(fit,"Rsig2",permute=FALSE,inc_warmup=FALSE)
plot(est.Rsig2,type='l');abline(h=log(Rsig2),col='red')

est.R0<-extract(fit,"R0",permute=FALSE,inc_warmup=FALSE)
plot(est.R0,type='l');abline(h=log(R0),col='red')

est.X0<-extract(fit,"X0",permute=FALSE,inc_warmup=FALSE)
plot(est.X0,type='l');abline(h=X0,col='red')
##overall, these look pretty good, only two notes:
##1) Rsig2 displays higher autocorrelation than the other parameters
##2) X0 displays preculiarly long, thin tails compared to other parameters

summary(lm(as.vector(apply(rates,c(2,3),median))~R))

plot(tree,edge.color=colvec,edge.width=3)
med.rates<-as.vector(apply(rates,c(2,3),mean))
tmp<-round((med.rates-min(R))/(max(R)-min(R))*99)+1
tmp[tmp<1]<-1;tmp[tmp>100]<-100
emp.colvec<-colramp[tmp]
plot(tree,edge.color=emp.colvec,edge.width=3)
