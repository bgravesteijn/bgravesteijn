#' Bootstrap with clustered data
#'
#' @param data data.frame or data.table 
#' @param cluster.v the grouping variable of interest
#' @param statistic the function to use on the bootstrap sample
#' @param R the number of bootstrap replicates
#' @description  Bootstrap with clustered data
cluster.boot<-function(data=NULL, cluster.v=NULL,statistic=NULL, R=NULL){
  #the function for making a clustered bootstrap sample
  cluster.boot.sample<-function(data=NULL, cluster.v=NULL){ 
    #cluster.v  is the vector in the dataframe which is the cluster variable
    n.cluster<-length(unique(cluster.v)) #how many clusters are there
    names.cluster<-unique(cluster.v)     #which clusters are present
    
    #Which clusters will be sampled from:
    samp.clust<-sample(x = names.cluster, size = n.cluster)
    
    #which patients will be sampled in each cluster
    samp.clust.pat<-vector(mode = "list", length = n.cluster)
    for(i in 1:n.cluster){
      cluster.pat<-data[cluster.v==samp.clust[i],]
      cluster.pat$index<-1:nrow(cluster.pat)
      samp.clust.pat[[i]]<-sample(x = cluster.pat$index, nrow(cluster.pat))
    }
    
    #From each cluster, sample patients from that cluster
    bootsamp<-vector(mode = "list", length = n.cluster)
    for (i in 1:n.cluster){
      bootsamp[[i]]<-data[cluster.v==samp.clust[i],][samp.clust.pat[[i]],]
    }
    bootsample<-do.call(what = "rbind", args = bootsamp)
  }
  
#the function for applying the statistic to the sample
  boot.res<-rep(NA, R)
  for (i in 1:R){
    cluster.boot.sample(data = data, cluster.v=cluster.v)
    boot.res[i]<-statistic(data)
  }
  return(boot.res)
}
