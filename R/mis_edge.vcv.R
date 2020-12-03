#get the edge-wise expected variance covariance evolving under BM with branch lengths equal to evolutionary rate
#' @export
edge.vcv<-function(tree){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  edge.hgts<-node.depth.edgelength(tree)[tree$edge[,2]]
  #base.hgts<-node.depth.edgelength(tree)[tree$edge[,1]]
  e<-nrow(tree$edge)
  mat<-matrix(0,e,e)
  end.pts<-vector(mode='integer',length=e)
  for(i in 1:e){
    sister.edge<-which(tree$edge[,1]==tree$edge[i,1])
    sister.edge<-suppressWarnings(min(sister.edge[sister.edge>i]))
    if(!is.infinite(sister.edge)){
      end.pts[i]<-sister.edge-1
    }else{
      anc.edge<-which(tree$edge[,2]==tree$edge[i,1])
      if(length(anc.edge)==0){
        end.pts[i]<-e
      }else{
        end.pts[i]<-end.pts[anc.edge]
      }
    }
    mat[i:end.pts[i],i:end.pts[i]]<-edge.hgts[i]-tree$edge.length[i]/2
    mat[(i:end.pts[i])[-1],(i:end.pts[i])[-1]]<-edge.hgts[i]
    mat[i,i]<-edge.hgts[i]-2*tree$edge.length[i]/3
  }
  mat
}

#now works with polytomies!