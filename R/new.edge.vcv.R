new.edge.vcv<-function(tree){
  attr(tree,'order')<-NULL
  tree<-reorder(tree,'cladewise')
  tree<-di2multi(tree,tol=0)
  edge.hgts<-node.depth.edgelength(tree)[tree$edge[,2]]
  e<-nrow(tree$edge)
  mat<-matrix(0,e,e)
  end.pts<-vector(mode='integer',length=e)
  for(i in 1:e){
    sister.edge<-which(tree$edge[,1]==tree$edge[i,1])
    sister.edge<-min(sister.edge[sister.edge!=i])
    if(sister.edge>i){
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

# #I feel like there's some way to speed this up more, but I can't figure it out...
# e<-nrow(tree$edge)
# mat2<-matrix(NA,e,e)
# end.pts<-vector(mode='integer',length=e)
# for(i in 1:e){
#   sister.edge<-which(tree$edge[,1]==tree$edge[i,1])
#   sister.edge<-min(sister.edge[sister.edge!=i])
#   if(sister.edge>i){
#     end.pts[i]<-sister.edge-1
#   }else{
#     anc.edge<-which(tree$edge[,2]==tree$edge[i,1])
#     if(length(anc.edge)==0){
#       end.pts[i]<-e
#     }else{
#       end.pts[i]<-end.pts[anc.edge]
#     }
#   }
#   mat2[i:end.pts[i],i:end.pts[i]]<-i
# }
# mat2[upper.tri(mat2)]<-NA
# get.combos<-function(vec){
#   cbind(unlist(lapply((1:length(vec))[-length(vec)],function(ii) vec[-(1:ii)])),rep(vec,(length(vec)-1):0))
# }
# 
# ancs<-vector(mode='list',length=nrow(tree$edge))
# for(i in 1:nrow(tree$edge)){
#   adj.edge<-which(tree$edge[,2]==tree$edge[i,1])
#   if(length(adj.edge)>0){
#     ancs[[i]]<-c(ancs[[adj.edge]],i)
#   }else{
#     ancs[[i]]<-i
#   }
# }
# mat<-matrix(NA,nrow=nrow(tree$edge),ncol=nrow(tree$edge))
# # get.mrca<-function(other.edge,ref.edge,ancs){
# #   ind<-which(ancs[[ref.edge]]%in%ancs[[other.edge]])
# #   if(length(ind)==0){
# #     NA
# #   }else{
# #     ancs[[ref.edge]][max(ind)]
# #   }
# # }
# for(i in 1:nrow(tree$edge)){
#   tmp.edge<-unlist(lapply((1:i)[-i],get.mrca,ref.edge=i,ancs=ancs))
#   mat[i,1:i]<-c(tmp.edge,i)
# }
# 
# #wooo! new one is much, much faster!