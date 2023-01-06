#get the edge-wise expected variance covariance evolving under BM with branch lengths equal to evolutionary rate
#now works with polytomies!
#way faster now using efficient tree-walking functions
#would break in case of singleton edges, but this should never really happen
#' @export
edge.vcv<-function(tree){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  edge.hgts<-node.depth.edgelength(tree)[tree$edge[,2]]
  e<-nrow(tree$edge)
  mat<-matrix(0,e,e)
  end.pts<-vector(mode='integer',length=e)
  anc<-anc.edges(tree)
  sis<-sis.edges(tree)
  for(i in seq_len(e)){
    tmp<-sis[[i]]>i
    if(any(tmp)){
      end.pts[i]<-min(sis[[i]][tmp])-1
    }else{
      if(length(anc[[i]])){
        end.pts[i]<-end.pts[anc[[i]]]
      }else{
        end.pts[i]<-e
      }
    }
    inds<-seq.int(i,end.pts[i])
    tmp<-tree$edge.length[i]/2
    mat[inds,inds]<-edge.hgts[i]-tmp
    if(length(inds)>1){
      mat[inds[-1],inds[-1]]<-mat[inds[-1],inds[-1]]+tmp
    }
    mat[i,i]<-mat[i,i]-tmp/3
  }
  mat
}

# edge.vcv2<-function(tree,zlim){
#   attr(tree,'order')<-NULL
#   tree<-reorder(tree,'cladewise')
#   tree<-di2multi(tree,tol=0)
#   edge.hgts<-node.depth.edgelength(tree)[tree$edge[,2]]
#   e<-nrow(tree$edge)
#   mat<-matrix(0,e,e)
#   end.pts<-vector(mode='integer',length=e)
#   anc<-anc.edges(tree)
#   sis<-sis.edges(tree)
#   for(i in seq_len(e)){
#     if(sis[[i]]>i){
#       end.pts[i]<-min(sis[[i]])-1
#     }else{
#       if(length(anc[[i]])){
#         end.pts[i]<-end.pts[anc[[i]]]
#       }else{
#         end.pts[i]<-e
#       }
#     }
#     inds<-seq.int(i,end.pts[i])
#     tmp<-tree$edge.length[i]/2
#     mat[inds,inds]<-edge.hgts[i]-tmp
#     mat[inds[-1],inds[-1]]<-mat[inds[-1],inds[-1]]+tmp
#     mat[i,i]<-mat[i,i]-tmp/3
#     png(paste0(i,'.png'),width=98*4,height=98*4)
#     image(mat,zlim=zlim)
#     dev.off()
#   }
#   mat
# }
# #quick  demonstration of algorithm
# set.seed(123)
# setwd('~/../Desktop/alg_demonstration')
# par(oma=c(0,0,0,0),mar=c(0,0,0,0))
# tree<-pbtree(n=50)
# zlim<-range(edge.vcv(tree))
# edge.vcv2(tree,zlim)
